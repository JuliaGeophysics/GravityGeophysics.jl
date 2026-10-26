# Predict the gravity of a model at the stations of a survey.
# Author: @pankajkmishra
# Writes a .pred file beside the data, with the same positions and errors. Run with -t auto:
#
#   julia --project=. -t auto examples/run_fwd3D.jl examples/demo1/model.true examples/demo1/data.xyz

using GravityGeophysics
using Printf

model_file = length(ARGS) >= 1 ? ARGS[1] : "examples/demo1/model.true"
data_file  = length(ARGS) >= 2 ? ARGS[2] : "examples/demo1/data.xyz"
fwd_file   = length(ARGS) >= 3 ? ARGS[3] : "examples/demo1/FwdCtrl"

out_path = joinpath(dirname(abspath(data_file)), "data.pred")

pred_path = ForwardSolve3D(model_file, data_file, fwd_file; output_path = out_path)

observed  = ReadData3D(data_file)
predicted = ReadData3D(pred_path)
fit       = chi2_rms3d(observed, predicted)

@printf("wrote %s\n", pred_path)
@printf("predicted gz %.4f .. %.4f mGal\n", minimum(predicted.value), maximum(predicted.value))
@printf("rms against the observed data: %.3f\n", fit.rms)
