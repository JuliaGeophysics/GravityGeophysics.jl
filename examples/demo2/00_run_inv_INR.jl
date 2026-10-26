# demo2: neural (INR) inversion of the synthetic survey.
# Run from the repository root with the project activated:
#
#   julia --project=. -t auto examples/demo2/00_run_inv_INR.jl
#
# Writes run_inr_<timestamp>/ next to this script: model.den, data.pred, History.log and
# Summary.txt. Gitignored; only the input files in this directory are tracked. Plot the
# result with 01_plot_data.jl, 02_plot_model_XY.jl, 03_plot_model_XZ.jl or 04_plot_model_YZ.jl.

using GravityGeophysics

const HERE = @__DIR__
here(name) = joinpath(HERE, name)

# the INR has no volume term, so InvCtrl.INR holds the padding fixed with "Invert padding : no"
inr = Invert3D(here("model.start"), here("data.xyz"), here("FwdCtrl"),
               here("InvCtrl.INR"), here("C3.dat"), here("model.prior");
               run_dir = here("run_inr"))

println("INR: rms $(round(inr.rms; digits = 2)) in $(inr.run_dir)")
