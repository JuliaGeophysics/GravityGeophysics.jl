# GravityGeophysics package module entry.
# Author: @pankajkmishra
# 3D gravity: density-contrast models on the same tensor mesh and model file as the resistivity
# and susceptibility models, so a density model opens in the same viewers. x north, y east, z DOWN.
module GravityGeophysics

using CairoMakie
using Dates
using GeoInterface
using LinearAlgebra
using Printf
using Random
using Shapefile
using SparseArrays
using Statistics: mean, median, quantile, std
import Proj

#----- Model files and masks ---------------------------------------------#
include("Model.jl")
include("Mask3D.jl")
include("MeshToMesh.jl")
#----- Data files (with the survey frame) and rotation -------------------#
include("Data.jl")
include("Rotate.jl")
#----- Core/padding utilities, shapefile overlays, georeferencing ---------#
include("CoreUtils3D.jl")
include("ShapefileOverlay.jl")
include("GeoRef3D.jl")
#----- Forward modelling --------------------------------------------------#
include("Fwd3D.jl")
#----- Inversion ----------------------------------------------------------#
include("Control3D.jl")
include("Inv3D.jl")
include("Inv3D_LSQR.jl")
#----- Mesh design ---------------------------------------------------------#
include("Lakes3D.jl")
include("MakeMesh3D.jl")

#----- Visualization (optional, requires GLMakie) -------------------------#

try
    using GLMakie
    include("PlotModel.jl")
    include("PlotModel3D.jl")
    include("EditModel3D.jl")
    include("MakeMesh3DGUI.jl")
    include("PlotData3D.jl")
    export compute_colorrange, prepare_model_arrays, save_figure
    export PlotModelXY, PlotModelXZ, PlotModelYZ, PlotModelXYZ
    export EditModelByLayers, EditModelByDrawing
    export plot_data, plot_data_fit, PlotDataMap, write_data_gis
catch LoadError
    @warn "GLMakie not available, interactive visualization functionality disabled"
end

export Model
export load_model, load_mesh, write_model
export air_cells, active_cells, same_mesh
export load_mask, write_mask, water_mask
export extract_topography, write_topography, read_topography, air_mask_from_topography, air_mask_from_model
export Frame, load_frame, to_frame, north_azimuth
export Data, ndata, load_data_grav, write_data_grav, COMPONENTS
export rotate_data, drape!, check_frame
export GravForward, forward, jvec, jtvec, model_vector, model_array, ForwardSolve3D
export FwdCtrl3D, InvCtrl3D
export ReadFwdCtrl3D, WriteFwdCtrl3D, ReadInvCtrl3D, WriteInvCtrl3D
export Invert3D, Inv3DOptions, Inv3DResult, Inv3DProblem, AbstractInversion3D
export inv3d_tag, inv3d_init, inv3d_prepare!, inv3d_direction, inv3d_reject!, inv3d_accept!, inv3d_info, inv3d_validate
export LSQR3D, LSQR3DConfig, LSQR3DResult
export objective, gradient, data_hessian, reg_hessian, phi_d, phi_m, rms, predicted
export edges_from_centers, core_indices, z_indices_for_max_depth, lateral_core_ranges, core_view
export detect_shapefile_crs, shapefile_coord_transform, load_shapefile_geometries, prepare_shapefiles, pick_shapefile
export MakeMesh3D, read_lakes
export MeshToMesh

export AIR_VALUE, COLORMAP, GAMMA, SI_TO_MGAL

end # module
