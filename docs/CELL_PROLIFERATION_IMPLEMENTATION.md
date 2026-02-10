# Cell Proliferation and Death Implementation for Intestinal Organoids

## Current State

### Cell-Cell Adhesion (Factor 1)
- **Implementation**: `GeneralisedLinearSpringForce`
- **Strength**: 0.4 (from Metzcar 2025)
- **Cutoff**: 50 µm
- **Model**: Linear spring force between neighboring cells
- **Limitation**: Uniform adhesion - doesn't distinguish between apical/basal cell types

### Cell Cycle Model
- **Current**: `UniformG1GenerationalCellCycleModel`
- **Proliferation**: DISABLED (`SetMaxTransitGenerations(0)`)
- **Birth time**: All cells born at t=0
- **Issue**: No cell division, population remains at 30 cells

## Built-in Chaste Components for Intestinal Organoids

### From Crypt Examples (Intestinal Epithelium)

Chaste has extensive crypt (intestinal) simulation support in `Chaste/crypt/`:

#### 1. **Cell Cycle Models Available**
```cpp
// Simple models with generation counting:
FixedG1GenerationalCellCycleModel         // Fixed G1 duration
UniformG1GenerationalCellCycleModel       // Random G1 from uniform distribution
ExponentialG1GenerationalCellCycleModel   // Random G1 from exponential distribution

// Wnt-dependent models (for stem cell niche):
WntCellCycleModel                         // Simple Wnt-based
SimpleWntCellCycleModel                   // Simplified version
StochasticWntCellCycleModel              // Stochastic Wnt
VanLeeuwen2009WntSwatCellCycleModel       // Full ODE-based Wnt

// Contact inhibition:
ContactInhibitionCellCycleModel           // Quiescence when compressed

// ODE-based:
TysonNovakCellCycleModel                  // Full molecular cell cycle
```

#### 2. **Cell Killers (for Apoptosis/Anoikis)**
```cpp
RandomCellKiller                          // Random probability of death
ApoptoticCellKiller                       // Kills apoptotic cells
SloughingCellKiller                       // Kills cells above height threshold
PlaneBasedCellKiller                      // Kills cells past a plane
TargetedCellKiller                        // Kills specific cells
IsolatedLabelledCellKiller               // Kills isolated labeled cells
```

#### 3. **Cell Generators**
```cpp
CryptCellsGenerator<CELL_CYCLE_MODEL>    // Generates cells by height (generations)
```

#### 4. **Simulation Modifiers**
```cpp
VolumeTrackingModifier                    // Track cell volumes for contact inhibition
AbstractTargetAreaModifier                // Target area/volume for cells
```

## Recommended Implementation for Organoids

### Phase 1: Basic Proliferation with Contact Inhibition

Use `ContactInhibitionCellCycleModel` - this is PERFECT for organoids!

```cpp
// In Test2dDynamicECMInvasion.hpp

// Replace UniformG1GenerationalCellCycleModel with:
#include "ContactInhibitionCellCycleModel.hpp"

for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
{
    ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
    p_cycle_model->SetDimension(2);
    
    // Cell cycle phase durations (from typical epithelial cells)
    p_cycle_model->SetStemCellG1Duration(8.0);   // 8 hours
    p_cycle_model->SetTransitCellG1Duration(12.0); // 12 hours  
    p_cycle_model->SetSDuration(5.0);            // 5 hours
    p_cycle_model->SetG2Duration(4.0);           // 4 hours  
    p_cycle_model->SetMDuration(1.0);            // 1 hour
    
    // Contact inhibition: quiescence at 70% of equilibrium volume
    p_cycle_model->SetQuiescentVolumeFraction(0.7);
    p_cycle_model->SetEquilibriumVolume(M_PI * 5.0 * 5.0); // ~78.5 µm² for radius=5µm
    
    CellPtr p_cell(new Cell(p_state, p_cycle_model));
    p_cell->SetCellProliferativeType(p_transit_type);
    p_cell->SetBirthTime(-p_gen->ranf() * 22.0); // Random age 0-22 hours
    
    cells.push_back(p_cell);
}
```

**Benefits**:
- ✅ Automatic contact inhibition (cells stop dividing when compressed)
- ✅ Realistic cell cycle: ~22 hours total (8+5+4+1=18h minimum, +random G1)
- ✅ No manual generation counting needed

### Phase 2: Add Cell Death (Anoikis)

```cpp
// Add after creating simulator
#include "RandomCellKiller.hpp"

// Low basal death rate (homeostasis)
MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.0005)); 
// 0.05% per hour = ~1% per day
simulator.AddCellKiller(p_cell_killer);
```

### Phase 3: Volume Tracking for Contact Inhibition

```cpp
#include "VolumeTrackingModifier.hpp"

// Enable volume tracking (needed by ContactInhibitionCellCycleModel)
MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
simulator.AddSimulationModifier(p_volume_modifier);
```

### Phase 4: ECM-Dependent Death (Anoikis)

Create custom cell killer that kills cells in low-ECM regions:

```cpp
// New file: ECMAnoikisCellKiller.hpp
template<unsigned DIM>
class ECMAnoikisCellKiller : public AbstractCellKiller<DIM>
{
private:
    boost::shared_ptr<DynamicECMField> mpECMField;
    double mCriticalECMDensity;  // Below this, cells die
    double mDeathProbability;     // Probability per hour
    
public:
    ECMAnoikisCellKiller(AbstractCellPopulation<DIM>* pCellPopulation,
                         boost::shared_ptr<DynamicECMField> pECMField,
                         double criticalDensity = 0.2,
                         double deathProb = 0.01)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mpECMField(pECMField),
          mCriticalECMDensity(criticalDensity),
          mDeathProbability(deathProb)
    {}
    
    void CheckAndLabelCellsForApoptosisOrDeath()
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            
            // Query ECM density at cell location
            c_vector<double, 2> pos;
            pos[0] = cell_location[0];
            pos[1] = cell_location[1];
            
            double ecm_density = mpECMField->GetDensityAt(pos);
            
            // Anoikis: detachment-induced death in low ECM
            if (ecm_density < mCriticalECMDensity)
            {
                double dt = SimulationTime::Instance()->GetTimeStep();
                double death_this_timestep = 1.0 - exp(-mDeathProbability * dt);
                
                if (RandomNumberGenerator::Instance()->ranf() < death_this_timestep)
                {
                    cell_iter->Kill();
                }
            }
        }
    }
};
```

## Expected Organoid Behavior

### With ContactInhibitionCellCycleModel + VolumeTracker:
1. **Early stage (t=0-24h)**: Cells divide freely, low density
2. **Growth phase (t=1-3 days)**: Exponential growth until crowding
3. **Homeostasis (t=3-5 days)**: Contact inhibition kicks in
   - Cells in center: Compressed → Quiescent (G0-like state)
   - Cells at edge: More space → Continue dividing
4. **Lumen formation**: Central cells die (anoikis in low ECM)

### Growth Curve Prediction:
- **t=0**: 30 cells
- **t=1 day**: ~60-90 cells (doubling time ~22h)
- **t=2 days**: ~120-180 cells
- **t=3 days**: ~200-300 cells (contact inhibition starts)
- **t=5 days**: ~400-600 cells (homeostatic balance)

Unlike Metzcar's invasion front (900 cells), organoids reach homeostasis earlier due to:
- ✅ Contact inhibition (not in Metzcar model)
- ✅ 3D-like packing constraints
- ✅ Cell death balancing proliferation

## Advanced: Spatial Heterogeneity

### Option A: Wnt Gradient (Stem Cell Niche)
```cpp
// Use SimpleWntCellCycleModel instead
#include "SimpleWntCellCycleModel.hpp"
#include "WntConcentration.hpp"

// Set up Wnt gradient (high at base)
WntConcentration<2>::Instance()->SetType(LINEAR);
WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
WntConcentration<2>::Instance()->SetCryptLength(1000.0); // Domain height

// Cells in high Wnt (bottom): Stem-like, fast proliferation
// Cells in low Wnt (top): Differentiated, slow/no proliferation
```

### Option B: ECM-Modulated Proliferation
Modify `ContactInhibitionCellCycleModel` to also check ECM:

```cpp
// In custom ECMContactInhibitionCellCycleModel:
double GetG1Duration()
{
    double base_duration = ContactInhibitionCellCycleModel::GetG1Duration();
    
    // Query ECM density
    double ecm_density = QueryECMAtCellLocation();
    
    // High ECM → Faster division (structural support)
    // Low ECM → Slower division  
    double ecm_factor = 0.5 + ecm_density; // Range: 0.5x to 1.5x
    
    return base_duration / ecm_factor;
}
```

## Implementation Checklist

- [ ] Replace `UniformG1GenerationalCellCycleModel` with `ContactInhibitionCellCycleModel`
- [ ] Add `VolumeTrackingModifier` 
- [ ] Set realistic cell cycle durations (G1=8-12h, S=5h, G2=4h, M=1h)
- [ ] Configure contact inhibition (quiescent fraction = 0.7)
- [ ] Add `RandomCellKiller` for basal death rate (~0.05%/hour)
- [ ] Optional: Create `ECMAnoikisCellKiller` for ECM-dependent death
- [ ] Optional: Add Wnt gradient for stem cell zonation
- [ ] Increase simulation time to 5 days (7200 minutes)
- [ ] Update visualization to show cell count over time

## Key Advantages Over Custom Implementation

1. **Built-in Contact Inhibition**: `ContactInhibitionCellCycleModel` handles density-dependent growth automatically
2. **Volume Tracking**: `VolumeTrackingModifier` computes cell volumes from Voronoi tessellation
3. **Validated Models**: Crypt models extensively tested for intestinal epithelium
4. **Wnt Signaling**: Ready-to-use stem cell niche gradients
5. **Multiple Cell Killers**: Combine different death mechanisms easily

## Files to Modify

1. `Test2dDynamicECMInvasion.hpp`:
   - Include `ContactInhibitionCellCycleModel.hpp`
   - Include `VolumeTrackingModifier.hpp`  
   - Include `RandomCellKiller.hpp`
   - Replace cell cycle model creation
   - Add volume modifier
   - Add cell killer
   - Increase EndTime to 7200 minutes (5 days)

2. **Optional** `ECMAnoikisCellKiller.hpp` (new file)

3. `visualize_dynamic_ecm.py`:
   - Plot cell count vs time
   - Show proliferation/death rates

## References

- Chaste Crypt Tutorials: `Chaste/crypt/test/`
- Contact Inhibition: `ContactInhibitionCellCycleModel.hpp`
- Volume Tracking: `VolumeTrackingModifier.hpp`  
- Cell Death: `AbstractCellKiller.hpp` and subclasses
- Wnt Signaling: `WntCellCycleModel.hpp`, `WntConcentration.hpp`
