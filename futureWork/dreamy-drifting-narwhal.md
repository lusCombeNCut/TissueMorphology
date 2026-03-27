# Parallelisation Opportunities for CryptBuddingApp ECM Calculations

## Context

ECM (ghost node) calculations consume ~90% of wall time in CryptBuddingApp. Currently the entire codebase is **single-threaded** -- no OpenMP, no MPI, no GPU. The user wants to identify where parallelisation could help, referencing Yalla's GPU approach as inspiration. The focus is on analysis rather than implementation (unless trivial).

---

## Bottleneck Anatomy

The hot path per timestep (at ghost spacing=1.0, ~9,500 active nodes, ~18-200 cells):

| Phase | Location | Work | % of ECM time (est.) |
|-------|----------|------|----------------------|
| Ghost-ghost forces | `GhostNodeEcmField.hpp:421-478` | ~28,500 pair evaluations (sqrt, log/exp, vector ops) | **60-70%** |
| Cell-ghost forces | `GhostNodeEcmForce.hpp:162-290` | ~200 cells x ~15 nearby ghosts = 3,000 pair evals | 15-25% |
| Spatial hash rebuild | `GhostNodeEcmField.hpp:1174-1202` | O(N_ghost) insert | 5% |
| Degradation + remodeling | `GhostNodeEcmField.hpp:508-520, 812-839` | O(N_cells x K) | 5% |

---

## Tier 1: Trivial Wins (hours of work)

### T1. Eliminate `FindReciprocalIndex` in viscoelastic variant
- **File**: `ViscoelasticGhostNodeEcmField.hpp` ~line 433, definition ~line 1181
- **Problem**: Linear search O(K) through neighbour list every pair, every timestep
- **Fix**: Precompute reciprocal indices at connectivity build time: `reciprocal_indices[nb] = index of i in j's neighbour list`
- **Impact**: 5-15% speedup on the viscoelastic ghost-ghost loop
- **Effort**: ~10 lines in `BuildNeighbourConnectivity` + 1 line in force loop

### T2. Pass precomputed distances to Degrade/Remodel
- **Files**: `GhostNodeEcmForce.hpp`, `GhostNodeEcmField.hpp`
- **Problem**: `DegradeNearby` and `RemodelNearby` recompute `norm_2(gn.position - cellPos)` for each nearby ghost -- already computed in the force loop
- **Fix**: Add `std::vector<double>& rDistances` output to `GetNearbyGhostNodes`, pass to degradation/remodeling
- **Impact**: Eliminates ~8,000 redundant sqrt/timestep at 200 cells. Modest but free
- **Effort**: ~15 lines across 3 functions

### T3. Granular profiling of sub-phases
- **Files**: Both force/field files
- **Fix**: Add `ScopedTimer` calls to separate ghost-ghost, cell-ghost, spatial hash, degradation, remodeling phases
- **Impact**: 0% performance, but essential baseline for all further work

---

## Tier 2: Medium-Effort CPU Parallelisation (days of work)

### M1. OpenMP-parallel ghost-ghost force loop (highest impact)
- **File**: `GhostNodeEcmField.hpp:421-478`
- **Current**: Newton's 3rd law loop -- each pair (i,j) processed once, both `node_i.force` and `node_j.force` modified. Creates write races if parallelised naively.
- **Fix**: **Drop Newton's 3rd law**. Each node independently computes forces from all its neighbours:
  ```cpp
  #pragma omp parallel for schedule(dynamic)
  for (unsigned i = 0; i < mNodes.size(); i++) {
      if (!mNodes[i].is_active) continue;
      for (unsigned n_idx : mNodes[i].neighbours) {
          // compute force from j on i (no j > i filter, no write to j)
          mNodes[i].force += compute_spring_force(i, n_idx);
      }
  }
  ```
  Doubles FLOPs (~60K evals vs ~30K), but on 4 cores: net ~2x speedup. On 8 cores: ~3-4x.
- **Build change**: Add `-fopenmp` to project's CMakeLists.txt via `target_compile_options`
- **Risk**: Chaste's `RandomNumberGenerator` is a singleton (not thread-safe), but the force loop doesn't use it. Safe.
- **Impact**: ~3-6x on dominant loop. **~2-4x overall** since this is 60-70% of ECM time.
- **Effort**: 1-2 days including testing

### M2. OpenMP-parallel cell-ghost force loop
- **File**: `GhostNodeEcmForce.hpp:162-290`
- **Challenge**: Cell iterator is `std::list`-based (not random-access). Also, reaction forces on ghost nodes create write races.
- **Fix**:
  1. Collect cell pointers into `std::vector` first
  2. Each thread uses a private `std::vector<c_vector<double,DIM>>` for ghost forces
  3. After parallel region, reduce (merge) thread-local ghost forces into actual ghost forces
- **Impact**: 2-3x on cell-ghost phase (~15-25% of ECM time)
- **Effort**: 1 day

### M3. Active-node index list
- **File**: `GhostNodeEcmField.hpp` (throughout)
- **Problem**: All loops iterate over full `mNodes` array checking `is_active`. As nodes degrade, an increasing fraction of iterations are wasted.
- **Fix**: Maintain `std::vector<unsigned> mActiveIndices`; update on removal
- **Impact**: Up to 30% if many nodes are removed during simulation
- **Effort**: ~0.5 day

### M4. Same as M1/M2 for Viscoelastic variant
- **Files**: `ViscoelasticGhostNodeEcmField.hpp`, `ViscoelasticGhostNodeEcmForce.hpp`
- **Additional complexity**: Per-pair state (dashpot rest lengths) requires read-only access from both i's and j's perspective when dropping N3L
- **Effort**: 1 day (after M1/M2 are validated)

---

## Tier 3: GPU Port (major overhaul, weeks)

### Yalla-Inspired GPU Architecture for Ghost Node ECM

**What Yalla does well:**
- SoA data layout for coalesced GPU memory access
- Sort-based spatial hashing: `compute_cube_id` kernel -> `thrust::sort_by_key` -> contiguous bucket iteration
- Stateless pairwise force functions as `__device__` functors
- `atomicAdd` for parallel force accumulation
- `thrust::fill` for zeroing, `thrust::reduce` for global constraints

**What would need to change:**

| Component | Current | GPU Target |
|-----------|---------|------------|
| Data layout | AoS (`std::vector<GhostNode>`) | SoA (separate arrays per field) |
| Spatial hash | Bucket-of-vectors, CPU | Sort-based, GPU (thrust) |
| Force compute | Sequential neighbour-list loop | CUDA kernel, 1 thread per node |
| N3L optimisation | Process each pair once | Drop N3L; each thread handles own node |
| Viscoelastic state | Per-pair rest_lengths in neighbour list | CSR-format edge storage on GPU |
| Cell-ghost interop | Same memory space | GPU->CPU position copy, CPU->GPU reaction force copy |

**Phase breakdown:**
1. **SoA migration** (1 week) -- benefits CPU too, prerequisite for GPU
2. **CUDA spatial hash** (1 week) -- Yalla's `solvers.cuh:346-503` as template
3. **CUDA force kernels** (1 week) -- ghost-ghost forces, position update
4. **Host-device interop** (3-5 days) -- cell-ghost stays on CPU, ghost data transfers
5. **Validation** (3-5 days) -- regression tests, benchmarking

**Break-even analysis:** GPU port is worthwhile only if ghost node count > ~50,000. Current 2D simulations (~10K nodes) would see modest GPU gains due to PCIe transfer overhead and kernel launch latency. **The CPU OpenMP approach (M1) gives comparable speedup with far less effort for current problem sizes.** GPU becomes compelling for 3D simulations with fine ghost node grids.

**Viscoelastic complication:** The SLS constitutive model stores per-pair state (dashpot rest lengths). Yalla's design is fundamentally stateless per pair. Options:
- Store edge state in CSR format on GPU (adds complexity)
- Redesign as node-centric Kelvin-Voigt model (changes physics)
- Keep viscoelastic ghost-ghost on CPU, only port elastic variant to GPU

---

## Tier 4: Algorithmic Improvements (no parallelism)

### A1. SoA data layout (CPU benefit)
- Convert `GhostNode` struct-of-arrays enables compiler auto-vectorisation (SIMD)
- Potential 20-50% speedup on force loops from better cache utilisation
- Effort: 2-3 days. Large restructuring but benefits all tiers.

### A2. Reduce spatial hash rebuild frequency
- Ghost nodes move slowly (overdamped). Track max displacement; only rebuild when displacement > cell_size/2
- Potential 5-10x fewer rebuilds, but rebuild is already cheap (~5% of ECM time)

---

## Recommended Order

| Priority | Item | Expected Gain | Effort |
|----------|------|---------------|--------|
| 1 | T3: Granular profiling | Diagnostic baseline | 1 hour |
| 2 | T1: Eliminate FindReciprocalIndex | 5-15% on VE loop | 30 min |
| 3 | T2: Pass precomputed distances | 2-5% on cell-ghost | 1 hour |
| 4 | **M1: OpenMP ghost-ghost loop** | **2-4x overall** | 1-2 days |
| 5 | M2: OpenMP cell-ghost loop | Additional 1.3-1.5x | 1 day |
| 6 | M3: Active-node index list | Up to 1.3x | 0.5 day |
| 7 | M4: Viscoelastic variants | Same as M1/M2 | 1 day |
| 8 | A1: SoA layout | 1.2-1.5x (enables SIMD) | 2-3 days |
| 9 | GPU port | 5-20x (only at >50K nodes) | 3-5 weeks |

**Realistic cumulative gain from Tiers 1+2 (M1 being the main driver): ~3-5x overall speedup on ECM, achievable in ~1 week of work.**

---

## Key Files

- `src/forces/GhostNodeEcmField.hpp` -- ghost node data, spatial hash, ghost-ghost forces
- `src/forces/GhostNodeEcmForce.hpp` -- cell-ghost force computation
- `src/forces/ViscoelasticGhostNodeEcmField.hpp` -- viscoelastic ghost variant
- `src/forces/ViscoelasticGhostNodeEcmForce.hpp` -- viscoelastic cell-ghost variant
- `src/forces/ECMConfinementForce.hpp` -- grid-based ECM (similar but less critical)
- `CMakeLists.txt` -- needs `-fopenmp` for M1/M2
