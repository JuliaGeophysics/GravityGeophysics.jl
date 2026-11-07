# 3D gravity forward run: model + data file (stations, data type, errors) + FwdCtrl -> data.pred
# Author: @pankajkmishra
# data.pred is written next to the data file, with a data-fit plot beside it
# FwdCtrl is shipped in examples/ctrl/3D; generate the models and data once: julia --project=. helpers/benchmarks_3D.jl
# Usage: julia --project=. --threads=8 examples/run_fwd3D.jl [model.den data.dat FwdCtrl]

using GravityGeophysics

const USER_ARGS = abspath.(ARGS)                 # file arguments, taken before the cd below
cd(@__DIR__)                                      # every path from here on is relative to this script
const CASE_DIR = "synthetic"
const CTRL_DIR = joinpath("ctrl", "3D")

inputs = if isempty(ARGS)
    [joinpath(CASE_DIR, "model.true"), joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "FwdCtrl")]
elseif length(ARGS) == 3
    relpath.(USER_ARGS)
else
    error("usage: julia --project=. examples/run_fwd3D.jl model.den data.dat FwdCtrl")
end
all(isfile, inputs) || error("missing inputs $(filter(!isfile, inputs)); run julia --project=. helpers/benchmarks_3D.jl first")

pred_path = ForwardSolve3D(inputs...)
d = load_data_grav(inputs[2])
plot_path = save_figure(joinpath(dirname(pred_path), "DataFit.png"), plot_data_fit(d, load_data_grav(pred_path).value))

println("Predicted : ", relpath(pred_path))
println("Plot      : ", relpath(plot_path))
