# demo2: compare True, GN and INR on an XY (map) slice.
# Run from the repository root with the project activated:
#
#   julia --project=. examples/demo2/02_plot_model_XY.jl
#
# Run 00_run_inv_GN.jl / 00_run_inv_INR.jl first; without them this plots model.true alone.
# Figure goes to plots/ next to this script.

using GravityGeophysics

const HERE = @__DIR__
here(name) = joinpath(HERE, name)
plots_dir  = here("plots")
mkpath(plots_dir)

# the most recently started run_<label>_* folder next to this script, or nothing
function latest_run(label)
    cands = filter(d -> isdir(joinpath(HERE, d)) && startswith(d, label * "_"), readdir(HERE))
    isempty(cands) ? nothing : joinpath(HERE, maximum(cands))
end

models = ["True" => here("model.true")]
for (label, tag) in (("run_gn", "GN"), ("run_inr", "INR"))
    dir = latest_run(label)
    dir === nothing || push!(models, tag => joinpath(dir, "model.den"))
end

fig, parts = PlotModelXY(models, here("data.xyz"); outline_file = here("model.true"), interactive = false)
out = joinpath(plots_dir, "ModelXY.png")
parts.export_figure(out)
println("plot in $out")
