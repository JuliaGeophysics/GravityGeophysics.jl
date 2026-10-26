# Design an inversion mesh around a survey.
# Author: @pankajkmishra
# Writes the start model and the covariance next to the data. Edit the settings block, then:
#
#   julia --project=. examples/make_mesh3D.jl examples/demo1/data.xyz

using GravityGeophysics

data_file      = length(ARGS) >= 1 ? ARGS[1] : "examples/demo1/data.xyz"
out_model      = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(abspath(data_file)), "mesh_start_model.den")
out_covariance = length(ARGS) >= 3 ? ARGS[3] : joinpath(dirname(out_model), "C3.dat")

background_density = 0.0      # kg/m3 contrast the mesh starts at
cell_width_frac    = 0.5      # core cell as a fraction of the station spacing
n_pad              = 6        # padding cells per side
pad_factor         = 1.4      # how fast the padding grows
first_layer_div    = 5.0      # top layer thickness = core cell / this
vertical_factor    = 1.15     # layer growth with depth
depth_mult         = 0.5      # model depth as a multiple of the narrower core width
cov_value          = 0.3      # smoothing written into the covariance file
cov_apply          = 2

# a callable (north, east) -> elevation in metres, or `nothing` for flat ground
topography = nothing

MakeMesh3D(data_file;
    out_model          = out_model,
    out_covariance     = out_covariance,
    topography         = topography,
    background_density = background_density,
    cell_width_frac    = cell_width_frac,
    n_pad              = n_pad,
    pad_factor         = pad_factor,
    first_layer_div    = first_layer_div,
    vertical_factor    = vertical_factor,
    depth_mult         = depth_mult,
    cov_value          = cov_value,
    cov_apply          = cov_apply,
)
