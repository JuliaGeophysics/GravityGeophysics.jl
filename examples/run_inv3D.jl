# Run a 3D inversion from the six input files.
# Author: @pankajkmishra
# The control file decides the algorithm: pass InvCtrl.GN for Gauss-Newton, InvCtrl.INR for the neural
# inversion. Everything else is identical. The run goes to a timestamped folder next to the data, and the
# standard figures (the model, the data fit and the convergence) go to run_dir/plots
#
#   julia --project=. -t auto examples/run_inv3D.jl examples/demo1 InvCtrl.GN

using GravityGeophysics
using Printf

demo_dir  = length(ARGS) >= 1 ? ARGS[1] : "examples/demo1"
inv_ctrl  = length(ARGS) >= 2 ? ARGS[2] : "InvCtrl.GN"
run_label = length(ARGS) >= 3 ? ARGS[3] : ""

at(name) = joinpath(demo_dir, name)

run = Invert3D(at("model.start"), at("data.xyz"), at("FwdCtrl"),
               at(inv_ctrl), at("C3.dat"), at("model.prior");
               run_dir = isempty(run_label) ? nothing : at(run_label))

plots_dir = joinpath(run.run_dir, "plots")

# the true model is only there in a synthetic test; it outlines the real bodies on every panel
true_path  = isfile(at("model.true")) ? at("model.true") : nothing
model_name = uppercase(string(run.algorithm))
models     = true_path === nothing ? [model_name => run.final] : ["True" => true_path, model_name => run.final]
for (view, plot) in (("XY", PlotModelXY), ("XZ", PlotModelXZ), ("YZ", PlotModelYZ))
    fig, parts = plot(models, at("data.xyz"); outline_file = something(true_path, ""), interactive = false)
    parts.export_figure(joinpath(plots_dir, "Model$(view).png"))
end
PlotData3D(at("data.xyz"); predicted_file = joinpath(run.run_dir, "data.pred"),
          output_path = joinpath(plots_dir, "DataFit.png"))
plot_inv3d_convergence(run.history; output_path = joinpath(plots_dir, "Convergence.png"),
                       target_rms = run.ctrl.target_rms)

@printf("%s finished: rms %.3f after %d iterations (%s)\n",
        model_name, run.rms, length(run.history) - 1, run.reason)
@printf("results and figures in %s\n", run.run_dir)
