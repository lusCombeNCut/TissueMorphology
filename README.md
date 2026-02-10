# TissueMorphology Project

A Chaste-based computational biology project for studying organoid formation and tissue morphology with basement membrane effects.

## Project Structure

```
TissueMorphology/
├── src/                          # Source code
│   ├── OrganoidCellFactory.hpp   # Cell factory with basement membrane properties
│   └── BasementMembraneForce.hpp # Basement membrane force implementation
├── test/                         # Test suite
│   └── TestOrganoidFormation.hpp # Comprehensive organoid formation tests
├── CMakeLists.txt               # Build configuration
├── visualize_chaste_results.py  # Main visualization script
├── create_visualizations.py     # Comparison visualization script
├── VISUALIZATION_GUIDE.md       # Detailed visualization documentation
└── README.md                    # This file
```

## Features

### Computational Models
- **OrganoidCellFactory**: Creates cell populations with basement membrane stiffness properties
- **BasementMembraneForce**: Applies radial constraints simulating basement membrane effects
- **Cell-based simulations**: Uses Chaste's MeshBasedCellPopulation framework

#### Mathematical Framework

The basement membrane force is modeled as a radial constraint:

$$\vec{F}_{BM}(i) = -k_{BM} \cdot \max(0, |\vec{r}_i - \vec{r}_{center}| - R_0) \cdot \hat{r}_i$$

where:
- $k_{BM}$ is the basement membrane stiffness parameter
- $\vec{r}_i$ is the position of cell $i$
- $\vec{r}_{center}$ is the organoid center (updated dynamically)
- $R_0$ is the baseline organoid radius
- $\hat{r}_i$ is the unit radial vector from center to cell $i$

Cell dynamics follow the overdamped equation:

$$\eta {d\vec{r}_i \over dt} = \vec{F}_{spring}(i) + \vec{F}_{BM}(i) + \vec{F}_{random}(i)$$

#### Key Assumptions
- Cells behave as point masses in 2D
- Basement membrane provides purely radial constraints
- Cell-cell interactions via generalized linear springs
- Overdamped dynamics (no inertial effects)
- Stochastic cell division and aging processes

### Test Suite  
- **Basic organoid formation**: Fundamental cell population dynamics
- **Basement membrane stiffness effects**: Comparison of low vs. high stiffness
- **Cell factory properties**: Validation of cell creation and property assignment
- **Integrated force testing**: Complete simulation with basement membrane forces

### Visualization System
- **Static plots**: Cell positions colored by age at specific time points
- **Animations**: Dynamic visualization of organoid formation process
- **Comparison figures**: Side-by-side analysis of different scenarios
- **Time series**: Evolution of morphology over simulation time

## Quick Start

### Building the Project
```bash
cd /home/orlando/Thesis
make TestOrganoidFormation
```

### Running Tests

**2D Organoid Tests (Basic validation):**
```bash
./projects/TissueMorphology/test/TestOrganoidFormation
```

**3D Organoid Tests (Full simulations):**
```bash
./projects/TissueMorphology/test/Test3dOrganoidFormation
```

### Creating Visualizations

**3D Organoid Visualization:**
```bash
cd /home/orlando/Thesis
python3 organoid_viewer_clean.py LongTermDevelopment
```

**2D Organoid Analysis:**
```bash
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology
python3 create_visualizations.py
```

## Generated Output

### Test Results

**2D Tests (`TestOrganoidFormation`):**
- Output: `/home/orlando/Thesis/testoutput/OrganoidFormation/`
- 4 scenarios: BasicTest, LowStiffness, HighStiffness, WithBasementMembrane
- Quick validation tests (0.2 time units, ~3 timesteps)

**3D Tests (`Test3dOrganoidFormation`):**
- Output: `/home/orlando/Thesis/testoutput/Organoid3d/`
- 3 scenarios: SphericalFormation, StiffnessTest, LongTermDevelopment
- Long simulations (20 time units, ~21 timesteps with VTU files)
- Rich 3D visualization data with basement membrane effects

### Visualization Files
- `organoid_comparison.png`: Multi-scenario comparison
- `organoid_time_series.png`: Temporal evolution analysis  
- `organoid_formation.gif`: Animated formation process
- Individual scenario plots for detailed analysis

## Scientific Applications

This project provides a foundation for studying:
- **Organoid morphogenesis**: How basement membrane properties affect tissue organization
- **Cell migration patterns**: Dynamic analysis of cell movement during development
- **Mechanical constraints**: Role of basement membrane forces in tissue architecture
- **Parameter sensitivity**: Effects of varying basement membrane stiffness on outcomes

## Dependencies

### Chaste Framework
- Chaste 2024.2 or compatible version
- cell_based components for cell population modeling
- VTK support for visualization output

### Python Visualization
- Python 3 with matplotlib and numpy
- PIL for animation generation
- Standard libraries for data processing

## Documentation

See `VISUALIZATION_GUIDE.md` for detailed visualization system documentation including:
- Complete usage instructions
- Scientific interpretation of results  
- Customization options
- Technical implementation details

## Contact

This project is part of ongoing computational biology research using the Chaste framework for tissue morphology studies.
