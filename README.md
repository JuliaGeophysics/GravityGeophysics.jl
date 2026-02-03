# GravityGeophysics.jl

[![Julia](https://img.shields.io/badge/Julia-1.10+-blue.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Julia package for 3D gravity forward modeling and analysis with UBC-GIF format support.

## Features

- **3D Mesh Creation**: Uniform and variable cell mesh support
- **UBC Format I/O**: Read/write mesh, model, and data files in UBC-GIF format
- **Forward Modeling Methods**:
  - **Prism (Analytical)**: Closed-form solution based on Plouff (1976) / Blakely (1995)
  - **Finite Difference**: Poisson equation solver with CG and Jacobi preconditioner
- **Model Building**: Block and sphere anomaly utilities
- **Visualization**: Comparison plots and model slices (requires Plots.jl)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/GravityGeophysics.jl")
```

Or for development:
```julia
] dev path/to/GravityGeophysics.jl
```

## Quick Start

```julia
using GravityGeophysics

# Create mesh (160 km × 160 km × 80 km domain)
mesh = create_mesh(
    nx = 64, ny = 64, nz = 40,
    lx = 160000.0, ly = 160000.0, lz = 80000.0
)

# Create model with density anomaly
model = create_model(mesh)
add_block_anomaly!(model;
    xmin = -4000.0, xmax = 4000.0,
    ymin = -4000.0, ymax = 4000.0,
    zmin = 8000.0, zmax = 16000.0,
    drho = 300.0  # kg/m³
)

# Create receivers
receivers = create_receivers(mesh; nrx=51, nry=51, z=50.0)

# Forward modeling - Prism method
data_prism = forward(model, receivers, Prism())

# Forward modeling - Finite Difference method
data_fd = forward(model, receivers, FiniteDifference())

# Compare methods
stats = compute_stats(data_prism, data_fd)

# Save in UBC format
save_data_ubc("gravity.obs", data_prism)
save_model_ubc("density.den", model)
```

## Package Structure

```
GravityGeophysics.jl/
├── Project.toml
├── README.md
├── LICENSE.md
├── src/
│   ├── GravityGeophysics.jl   # Main module
│   ├── Mesh.jl                # Mesh and model types
│   ├── UBCFormat.jl           # UBC-GIF I/O
│   ├── ForwardPrism.jl        # Analytical prism method
│   ├── ForwardFD.jl           # Finite difference method
│   ├── Forward.jl             # Unified interface
│   └── Visualization.jl       # Plotting utilities
├── examples/
│   └── 01_gravity_forward.jl  # Basic forward modeling example
└── test/
    └── runtests.jl
```

## API Reference

### Mesh Creation

```julia
# Create uniform mesh
mesh = create_mesh(nx=64, ny=64, nz=32, lx=10000.0, ly=10000.0, lz=5000.0)

# Get cell centers and edges
xc, yc, zc = get_cell_centers(mesh)
xe, ye, ze = get_cell_edges(mesh)

# Print mesh info
mesh_info(mesh)
```

### Model Creation

```julia
# Create empty model
model = create_model(mesh; background=0.0)

# Add anomalies
add_block_anomaly!(model; xmin, xmax, ymin, ymax, zmin, zmax, drho)
add_sphere_anomaly!(model; cx, cy, cz, radius, drho)
```

### Forward Modeling

```julia
# Prism method (exact, slower for large meshes)
data = forward(model, receivers, Prism())

# Finite Difference method (faster for large meshes)
fd_opts = FDOptions(tol=1e-10, maxiter=5000, use_prism_bc=false)
data = forward(model, receivers, FiniteDifference(fd_opts))

# Compare methods
data_prism, data_fd, stats = compare_methods(model, receivers)
```

### UBC Format I/O

```julia
# Mesh
save_mesh_ubc("mesh.msh", mesh)
mesh = load_mesh_ubc("mesh.msh")

# Model
save_model_ubc("density.den", model)
model = load_model_ubc("density.den", mesh)

# Data
save_data_ubc("gravity.obs", data)
data = load_data_ubc("gravity.obs")
```

## Forward Modeling Methods

### Prism Method (Analytical)

Computes the vertical gravity component (gz) using the closed-form rectangular prism formula:

$$g_z = G \rho \sum_{i=1}^{2} \sum_{j=1}^{2} \sum_{k=1}^{2} \mu_{ijk} \left[ x \ln(y+r) + y \ln(x+r) - z \arctan\frac{xy}{zr} \right]$$

where $r = \sqrt{x^2 + y^2 + z^2}$ and $\mu_{ijk} = (-1)^{i+j+k}$.

**Characteristics:**
- Exact solution for homogeneous rectangular prisms
- O(N_cells × N_receivers) complexity
- Best for small to medium-sized problems

### Finite Difference Method

Solves the 3D Poisson equation:

$$\nabla^2 \Phi = 4\pi G \rho$$

Using:
- 7-point Laplacian stencil
- Dirichlet boundary conditions (monopole approximation or exact prism)
- CG solver with Jacobi preconditioner

Gravity computed via: $g_z = -\partial\Phi/\partial z$

**Characteristics:**
- Better scaling for large meshes
- Approximate solution
- 2nd or 4th order derivative schemes available

## References

1. Plouff, D., 1976, Gravity and magnetic fields of polygonal prisms and application to magnetic terrain corrections: Geophysics, 41, 727-741.

2. Blakely, R.J., 1995, Potential Theory in Gravity and Magnetic Applications, Cambridge University Press.

3. UBC-GIF, GRAV3D Documentation, https://grav3d.readthedocs.io/

## License

MIT License - see [LICENSE.md](LICENSE.md)

## Contributing

Contributions welcome! Please open an issue or submit a pull request.
