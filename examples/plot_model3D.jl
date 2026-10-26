# Open a density model in an interactive viewer
# Author: @pankajkmishra
# `view` picks the window: xy for depth slices, xz and yz for sections, xyz for the three planes in 3D.
# Each has a slider with Prev/Next/First/Last, a core/padding toggle and an Export Figure button that
# writes a timestamped PNG. Passing the survey draws the stations, passing a true model outlines its bodies
#
#   julia --project=. examples/plot_model3D.jl examples/demo1/model.true xy

using GravityGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : "examples/demo1/model.true"
view       = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "xy"
data_file  = length(ARGS) >= 3 ? ARGS[3] : joinpath(dirname(abspath(model_file)), "data.xyz")

outline_file  = ""       # a true model on the same mesh, to outline its bodies on every slice
with_padding  = false    # false shows the core only, which is what the survey constrains
max_depth     = nothing  # metres, or nothing for the whole model
density_range = nothing  # (lo, hi) kg/m3, or nothing to scale from the model
show_grid     = false    # draw the cell edges over the slice
output_dir    = dirname(abspath(model_file))

data = isfile(data_file) ? data_file : ""

if view == "xyz"
    PlotModelXYZ(model_file, data; with_padding, max_depth, density_range, output_dir)
else
    plot = view == "xy" ? PlotModelXY : view == "xz" ? PlotModelXZ :
           view == "yz" ? PlotModelYZ : error("view must be xy, xz, yz or xyz")
    plot(model_file, data; outline_file, with_padding, max_depth, density_range, show_grid, output_dir)
end
