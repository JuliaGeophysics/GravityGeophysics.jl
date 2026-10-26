# demo2: Gauss-Newton inversion of the synthetic survey.
# Run from the repository root with the project activated:
#
#   julia --project=. -t auto examples/demo2/00_run_inv_GN.jl
#
# Writes run_gn_<timestamp>/ next to this script: model.den, data.pred, History.log and
# Summary.txt. Gitignored; only the input files in this directory are tracked. Plot the
# result with 01_plot_data.jl, 02_plot_model_XY.jl, 03_plot_model_XZ.jl or 04_plot_model_YZ.jl.

using GravityGeophysics

const HERE = @__DIR__
here(name) = joinpath(HERE, name)

# GN inverts the padding; its sensitivity weighting pays for the volume of a padding cell
gn = Invert3D(here("model.start"), here("data.xyz"), here("FwdCtrl"),
              here("InvCtrl.GN"), here("C3.dat"), here("model.prior");
              run_dir = here("run_gn"))

println("GN: rms $(round(gn.rms; digits = 2)) in $(gn.run_dir)")
