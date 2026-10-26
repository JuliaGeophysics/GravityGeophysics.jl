# Map a gravity survey, and the fit when a prediction is given.
# Author: @pankajkmishra
# With one file you get the observed map; with a .pred file as well you get observed, predicted and the
# residual in mGal.
#
#   julia --project=. examples/plot_data3D.jl examples/demo1/data.xyz

using GravityGeophysics
using CairoMakie: Reverse

data_file      = length(ARGS) >= 1 ? ARGS[1] : "examples/demo1/data.xyz"
predicted_file = length(ARGS) >= 2 ? ARGS[2] : ""

data_range = nothing   # (lo, hi) mGal on the observed/predicted panels, or nothing to autoscale
markersize = 5
colormap   = Reverse(:Spectral)

path = PlotData3D(data_file;
    predicted_file = predicted_file,
    data_range     = data_range,
    markersize     = markersize,
    colormap       = colormap,
)

println("wrote ", path)
