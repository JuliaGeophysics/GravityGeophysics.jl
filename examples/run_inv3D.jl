# 3D gravity deterministic inversion (LSQR) from six files
# Author: @pankajkmishra
# The algorithm and its settings come from the inversion control; results go to run_YYYYmmdd_HHMMSS/ next to the data
# Controls are shipped in examples/ctrl/3D; generate the models and data once: julia --project=. helpers/benchmarks_3D.jl
# Usage: julia --project=. --threads=8 examples/run_inv3D.jl [LSQR]
#        julia --project=. --threads=8 examples/run_inv3D.jl model.start data.dat FwdCtrl InvCtrl model.mask model.prior

using GravityGeophysics
using Printf

const USER_ARGS = abspath.(ARGS)                 # file arguments, taken before the cd below
cd(@__DIR__)                                      # every path from here on, and in the run log, is relative to this script
const CASE_DIR = "synthetic"
const CTRL_DIR = joinpath("ctrl", "3D")

inputs = if length(ARGS) <= 1
    algorithm = uppercase(isempty(ARGS) ? "LSQR" : ARGS[1])
    [joinpath(CASE_DIR, "model.start"), joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "FwdCtrl"),
     joinpath(CTRL_DIR, "InvCtrl.$algorithm"), joinpath(CASE_DIR, "model.mask"), joinpath(CASE_DIR, "model.prior")]
elseif length(ARGS) == 6
    relpath.(USER_ARGS)
else
    error("usage: julia --project=. examples/run_inv3D.jl [LSQR] or model.start data.dat FwdCtrl InvCtrl model.mask model.prior")
end
all(isfile, inputs) || error("missing inputs $(filter(!isfile, inputs)); run julia --project=. helpers/benchmarks_3D.jl first")

elapsed = @elapsed run = Invert3D(inputs...)

println()
@printf("Algorithm   : %s\n", uppercase(string(run.algorithm)))
@printf("Termination : %s (converged = %s), %d iterations\n", run.reason, run.converged, length(run.history) - 1)
@printf("Final RMS   : %.3f, %.1f s\n", run.rms, elapsed)
println("Outputs     : ", relpath(run.run_dir), "  (model.den, data.pred, run log)")
