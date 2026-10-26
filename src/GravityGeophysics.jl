# GravityGeophysics package module entry.
# Author: @pankajkmishra
# This file includes all core source scripts and defines the exported public API.

module GravityGeophysics

using LinearAlgebra
using Statistics
using Dates
using Printf
using SparseArrays
using Proj

#----- Data and model I/O (.xyz data, block-layout models) -----------------#

include("Data.jl")
include("Model.jl")

#----- Headless core/padding utilities (always available) ------------------#

include("CoreUtils3D.jl")

#----- Control files -------------------------------------------------------#

include("Control3D.jl")

#----- Forward solvers, synthetic models and misfit -------------------------#

include("Fwd3D.jl")
include("Fwd3D_FD.jl")
include("Synthetic3D.jl")

#----- Gauss-Newton inversion ----------------------------------------------#

include("Inv3D.jl")
include("Inv3D_GN.jl")

#----- Implicit-neural-representation inversion (Lux.jl) -------------------#

include("Inv3D_INR.jl")

#----- Mesh design ---------------------------------------------------------#

include("MakeMesh3D.jl")

#----- Headless plots (CairoMakie) -----------------------------------------#

include("PlotModel.jl")
include("PlotData3D.jl")

#----- Interactive model viewers (optional, requires GLMakie) --------------#

try
    using GLMakie
    include("PlotModel3D.jl")
    export depth_slice_viewer, xz_slice_viewer, yz_slice_viewer, grav_3d_viewer, draw_north_and_scale!
    export PlotModelXY, PlotModelXZ, PlotModelYZ, PlotModelXYZ
catch err
    @warn "GLMakie not available, the interactive model viewers are disabled" err
end

#----- Public API ----------------------------------------------------------#

export GravData3D, ReadData3D, WriteData3D, data_like, ndata, has_geographic
export GravModel3D, ReadModel3D, WriteModel3D, model_like, same_mesh, air_mask, cell_volumes, surface_depth, sensor_depth
export edges_from_centers, core_indices, z_indices_for_max_depth, lateral_core_ranges
export FwdCtrl3D, ReadFwdCtrl3D, WriteFwdCtrl3D
export InvCtrl3D, ReadInvCtrl3D, WriteInvCtrl3D
export InvCtrlINR3D, ReadInvCtrlINR3D, WriteInvCtrlINR3D, ReadInvCtrl
export Cov3D, ReadCov3D, WriteCov3D
export GravKernel3D, sensitivity_row, forward3d, Sensitivity3D, SaveSensitivity3D, LoadSensitivity3D, check_receivers, ForwardSolve3D
export FDCtrl3D, forward3d_fd
export sphere_anomaly, add_prism!, add_dyke!, add_noise, chi2_rms3d
export AbstractInversion3D, Inv3DOptions, Inv3DResult, Invert3D, model_weights
export GaussNewton3DConfig, GaussNewton3DResult, GaussNewton3D
export INR3DConfig, INR3DResult, INR3D
export MakeMesh3D, design_mesh, build_mesh_bundle, data_spacing
export compute_colorrange, prepare_model_arrays, model_slice
export plot_grav3d_data, plot_inv3d_convergence, PlotData3D
export GAMMA, SI_TO_MGAL

end # module
