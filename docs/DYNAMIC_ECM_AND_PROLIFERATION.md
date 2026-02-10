# Dynamic ECM Remodeling and Cell Proliferation

## Issues with Current Implementation

Based on your observations, the current model has several limitations:

### 1. Static ECM (Fixed) ❌
- **Current**: ECM orientation is fixed throughout simulation
- **Reality**: Cells remodel ECM through:
  - Mechanical pulling (fiber alignment)
  - Enzymatic degradation (MMPs)
  - ECM secretion and deposition
  - Contact guidance feedback loops

### 2. Unidirectional Fiber Movement (Fixed) ✅
- **Current**: Cells moved only in arrow direction
- **Fixed**: Now bidirectional + perpendicular diffusion + random walk
- **Reality**: Fibers are non-polar - cells can move either direction along them

### 3. Missing Lateral Diffusion (Fixed) ✅
- **Current**: Pure alignment with no cross-fiber motion
- **Fixed**: Added perpendicular mobility (reduced) and random walk
- **Reality**: Cells have intrinsic motility even in aligned ECM

### 4. No Cell Proliferation/Death ❌
- **Current**: Fixed cell count (30 cells, no division)
- **Reality**: Cells proliferate, differentiate, undergo apoptosis

## Dynamic ECM Remodeling (From Metzcar 2025)

### Model from Paper

Cells mechanically remodel ECM through:

```cpp
// ECM fiber alignment from cell traction
θ_ECM(t+1) = θ_ECM(t) + α * (θ_cell - θ_ECM(t))

where:
  α = remodeling rate (0.01 - 0.1 per timestep)
  θ_cell = cell's traction direction
  θ_ECM = current ECM orientation
```

### Implementation: Dynamic ECM Field

To implement this, we need to:

1. **Store ECM as a field** (not compute on-the-fly)
2. **Update ECM based on cell forces**
3. **Track ECM density/degradation**

#### Option A: Grid-Based ECM Field

Create a spatial grid storing ECM properties:

```cpp
class ECMField
{
private:
    // Grid parameters
    double mGridSpacing;  // e.g., 10 µm
    std::map<std::pair<int,int>, ECMProperties> mGrid;
    
    struct ECMProperties
    {
        double fiber_angle;      // Current orientation
        double density;          // ECM amount (0-1)
        double anisotropy;       // How aligned (0-1)
    };
    
public:
    void UpdateFromCellTraction(c_vector<double,2> position,
                               c_vector<double,2> traction_force)
    {
        // Get grid cell
        auto grid_idx = GetGridIndex(position);
        
        // Calculate cell's preferred direction
        double cell_angle = atan2(traction_force[1], traction_force[0]);
        
        // Update ECM orientation (weighted average)
        double alpha = 0.05;  // Remodeling rate
        double current_angle = mGrid[grid_idx].fiber_angle;
        
        // Handle angle wrapping
        double angle_diff = AngleDifference(cell_angle, current_angle);
        mGrid[grid_idx].fiber_angle += alpha * angle_diff;
        
        // Update density (cells deposit ECM)
        mGrid[grid_idx].density += 0.001;  // Slow deposition
        mGrid[grid_idx].density = std::min(1.0, mGrid[grid_idx].density);
    }
    
    double GetOrientationAt(c_vector<double,2> position)
    {
        auto grid_idx = GetGridIndex(position);
        return mGrid[grid_idx].fiber_angle;
    }
};
```

#### Option B: Cell-Based ECM Remodeling (Simpler)

Store ECM state with cells, interpolate between them:

```cpp
// In ECMContactGuidanceForce::AddForceContribution()

// Get cell's migration direction
c_vector<double, DIM> cell_velocity = p_node->rGetAppliedForce();

// Store as "ECM deformation" - cells pull ECM in their direction
double cell_direction_angle = atan2(cell_velocity[1], cell_velocity[0]);
cell_iter->GetCellData()->SetItem("ecm_deformation_angle", cell_direction_angle);
cell_iter->GetCellData()->SetItem("ecm_remodeling_strength", norm_2(cell_velocity));

// In next timestep, blend ECM field with cell deformations
// This creates feedback: cells align ECM → ECM guides cells → ...
```

### ECM Degradation

Add MMP (matrix metalloproteinase) activity:

```cpp
// Cells degrade ECM as they migrate
double degradation_rate = 0.01;  // per minute
double current_density = GetECMDensityAt(position);
new_density = current_density * (1.0 - degradation_rate * dt);

// ECM guides less when degraded
effective_anisotropy = base_anisotropy * current_density;
```

## Cell Proliferation and Death

### From Metzcar 2025 Paper

The paper mentions:
- **Contact inhibition**: Cells stop dividing when crowded
- **Nutrient/oxygen gradients**: Affect proliferation rate
- **Cell cycle**: G1/S/G2/M phases with stochastic transitions

### Implementation Options

#### Option 1: Simple Contact-Dependent Proliferation

```cpp
// In Test2dInvasiveFront, change cell cycle model:

// Replace:
p_cycle_model->SetMaxTransitGenerations(0);  // No division

// With:
p_cycle_model->SetTransitCellG1Duration(12.0);  // 12 hours
p_cycle_model->SetMaxTransitGenerations(5);     // Allow divisions

// Add contact inhibition force
MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
simulator.AddSimulationModifier(p_modifier);
```

#### Option 2: Boundary Birth Modifier (Like Paper)

Create custom modifier for cell influx at boundary:

```cpp
template<unsigned DIM>
class BoundaryBirthModifier : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    double mBirthRate;        // Cells per time
    double mBoundaryY;        // Y position for birth
    double mLastBirthTime;
    unsigned mCellsPerBirth;
    
public:
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        double current_time = SimulationTime::Instance()->GetTime();
        
        // Birth every 180 minutes (from paper)
        if (current_time - mLastBirthTime >= 180.0)
        {
            // Add 30 new cells at bottom boundary
            for (unsigned i = 0; i < 30; i++)
            {
                double x = RandomNumberGenerator::Instance()->ranf() * 600.0;
                double y = mBoundaryY + RandomNumberGenerator::Instance()->ranf() * 10.0;
                
                // Create new cell
                Node<DIM>* p_node = new Node<DIM>(GetNextNodeIndex(), 
                                                  false, x, y);
                rCellPopulation.rGetMesh().AddNode(p_node);
                
                // Create cell with cycle model
                // ... (similar to initial setup)
            }
            
            mLastBirthTime = current_time;
        }
    }
};
```

#### Option 3: Nutrient-Dependent Growth

Add diffusible nutrient field:

```cpp
// Add PDE for oxygen/nutrient diffusion
// Cells consume nutrients, nutrients diffuse from boundary
// Proliferation rate depends on local nutrient concentration

MAKE_PTR(CellwiseSourceEllipticPde<2>, p_pde);
p_pde->SetCoefficient(0.1);  // Diffusion
p_pde->SetLinearInUTerm(-0.01);  // Consumption

MAKE_PTR(AveragedSourceEllipticPde<2>, p_nutrient_pde);
// ... setup PDE solver
```

## Biochemical Signaling (Advanced)

### Growth Factors and Chemotaxis

The paper can include:

1. **Chemoattractants**: Cells migrate up/down gradients
2. **Growth factors**: Control proliferation spatially
3. **Morphogens**: Pattern formation (Turing patterns, etc.)

### Implementation with Subcellular Reaction Networks (SRN)

Chaste supports SRNs for biochemical signaling:

```cpp
// Add Srn model to cells
#include "DeltaNotchSrnModel.hpp"

DeltaNotchSrnModel* p_srn = new DeltaNotchSrnModel();
CellPtr p_cell(new Cell(p_state, p_cycle_model, p_srn));

// Or custom ODE system
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"

class GrowthFactorOdeSystem : public AbstractOdeSystem
{
    // Implement signal transduction pathway
    // E.g., MAPK/ERK, Wnt, Notch, etc.
};
```

## Recommended Next Steps

### Priority 1: Improve Migration Realism ✅
- [x] Bidirectional fiber movement
- [x] Perpendicular diffusion
- [x] Random walk component

### Priority 2: Add Cell Proliferation
- [ ] Implement BoundaryBirthModifier
- [ ] Add to Test2dInvasiveFront
- [ ] Verify cell count increases (30 → ~900 over 5 days)

### Priority 3: Dynamic ECM (Optional)
- [ ] Create ECMField class with grid storage
- [ ] Implement ECM-cell feedback loop
- [ ] Test ECM alignment from cell traction

### Priority 4: Biochemical Signals (Advanced)
- [ ] Add nutrient PDE
- [ ] Implement chemotaxis force
- [ ] Couple to proliferation rate

## Testing Strategy

### Phase 1: Verify Improved Migration
Run with new ECMContactGuidanceForce:
```bash
make Test2dInvasiveFront && ./projects/TissueMorphology/test/Test2dInvasiveFront
python3 create_invasion_animations.py
```

**Expected changes**:
- Cells spread laterally (not just vertically in perpendicular case)
- Some cells migrate "backwards" along fibers
- Invasion front is less sharp, more diffusive

### Phase 2: Add Proliferation
- Cell count should grow from 30 → ~900 in 5 days
- Denser packing in rear, sparse at front
- Matches paper's Figure 3

### Phase 3: Dynamic ECM
- ECM should align perpendicular to invasion front in "random" case
- "Tracks" form behind leading cells
- Feedback creates finger-like protrusions

## Code Examples

See:
- **ECMContactGuidanceForce.hpp** (updated with realistic migration)
- **Test2dInvasiveFront.hpp** (cell influx TODO comment)
- **Chaste tutorials**: `TestRunningDifferentialAdhesionSimulationsTutorial.hpp`
- **Proliferation example**: `TestOffLatticeSimulationWithPdes.hpp`

## References

Metzcar et al. 2025:
- Section 3.1: Invasive front model
- Supplementary Section 5.1: Full parameter list
- Supplementary Table 19: Invasion depth statistics
- Note: "instantiate 30 new cells every 180 minutes along domain edge"

Painter 2009:
- Original ECM contact guidance model
- Showed perpendicular > random > parallel invasion speed
