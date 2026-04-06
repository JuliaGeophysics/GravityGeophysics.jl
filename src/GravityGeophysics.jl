# Author: @pankajkmishra
"""
Julia package for 3D gravity forward modeling, inversion, and visualization.

The active exchange workflow uses voxel models (`.vox`) and gravity or receiver tables
in xyz text format (`.xyz`).
"""
module GravityGeophysics

using CairoMakie
using Flux
using GLMakie
using GeoInterface
using LinearAlgebra
using Random
using Proj
using Shapefile
using SparseArrays
using Statistics
using Printf

const G_NEWTON = 6.67430e-11
const SI_TO_MGAL = 1e5
const has_visualization = true

include("Mesh.jl")
include("Geometry.jl")
include("VoxelIO.jl")
include("ForwardPrism.jl")
include("ForwardFD.jl")
include("Forward.jl")
include("Synthetic.jl")
include("Inversion.jl")
include("INRInversion.jl")
include("GeoOverlay.jl")
include("PlotData.jl")
include("PlotModel3D.jl")

export GravityMesh, create_mesh, mesh_info
export get_cell_centers, get_cell_edges, get_nodes
export cell_volumes, mesh_cell_count

export GravityModel, create_model, set_density!
export add_block_anomaly!, add_sphere_anomaly!

export GravityReceivers, create_receivers
export load_model_vox, save_model_vox
export load_data_xyz, save_data_xyz, save_receivers_xyz

export ForwardMethod, Prism, FiniteDifference, FDOptions
export forward, forward_prism, forward_fd, forward_with_info
export GravityData

export edges_from_centers, sample_polyline, build_section_surface_polyline

export create_terrascope_demo_mesh
export create_terrascope_synthetic_model
export create_demo_receivers
export generate_demo_bundle
export create_block_benchmark_mesh
export create_block_benchmark_model
export create_block_benchmark_receivers
export generate_block_benchmark_bundle
export create_pyhasalmi_mesh
export create_pyhasalmi_synthetic_model
export create_pyhasalmi_receivers
export generate_pyhasalmi_bundle

export sensitivity_matrix_prism
export gradient_operators
export depth_weighting
export invert_tikhonov
export invert_inr
export compare_inversion_methods

export plot_gravity_comparison
export plot_gravity_data
export plot_model_slice
export plot_inversion_comparison
export plot_orthogonal_inversion_comparison
export plot_convergence
export plot_receiver_locations

export plot_shapefile_on_3d!
export plot_shapefile_on_axis!
export launch_model_viewer

export G_NEWTON, SI_TO_MGAL
export compute_stats, compare_methods

end # module
