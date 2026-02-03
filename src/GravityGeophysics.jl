# Author: @pankajkmishra
"""
# GravityGeophysics.jl

Julia package for 3D gravity forward modeling and analysis.

## Features
- 3D mesh creation with UBC format support
- Forward modeling using:
  - Prism (Analytical) method - Plouff (1976) / Blakely (1995)
  - Finite Difference (FD) method - Poisson equation solver
- UBC format I/O for models and data
- Visualization and comparison tools

## Quick Start
```julia
using GravityGeophysics

# Create mesh
mesh = create_mesh(nx=64, ny=64, nz=32, 
                   lx=160000.0, ly=160000.0, lz=80000.0)

# Load model in UBC format
model = load_model_ubc("model.ubc", mesh)

# Create receivers
receivers = GravityReceivers(mesh, nrx=51, nry=51, z=50.0)

# Forward modeling (Prism method)
data = forward(model, receivers, Prism())

# Or use FD method
data_fd = forward(model, receivers, FiniteDifference())

# Save data in UBC format
save_data_ubc("gravity_data.ubc", data, receivers)
```
"""
module GravityGeophysics

using LinearAlgebra
using Statistics
using Printf
using DelimitedFiles

# Physical constants
const G_NEWTON = 6.67430e-11  # Gravitational constant [m³/(kg·s²)]
const SI_TO_MGAL = 1e5        # Conversion factor SI to mGal

# Include submodules
include("Mesh.jl")
include("UBCFormat.jl")
include("ForwardPrism.jl")
include("ForwardFD.jl")
include("Forward.jl")

# Conditionally include visualization
has_visualization = false
try
    using Plots
    include("Visualization.jl")
    global has_visualization = true
    export plot_gravity_comparison, plot_model_slice
catch e
    @warn "Plots not available, visualization functionality disabled"
end

# Export mesh types and functions
export GravityMesh, create_mesh, mesh_info
export get_cell_centers, get_cell_edges, get_nodes

# Export model types and functions
export GravityModel, create_model, set_density!
export add_block_anomaly!, add_sphere_anomaly!

# Export receiver types
export GravityReceivers, create_receivers

# Export UBC I/O functions
export load_model_ubc, save_model_ubc
export load_data_ubc, save_data_ubc
export load_mesh_ubc, save_mesh_ubc

# Export forward modeling
export ForwardMethod, Prism, FiniteDifference, FDOptions
export forward, forward_prism, forward_fd, forward_with_info

# Export data types
export GravityData

# Export utilities
export G_NEWTON, SI_TO_MGAL
export compute_stats, compare_methods

end # module
