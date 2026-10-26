# demo1: map the survey, and the data fit and convergence of the latest GN/INR runs if present.
# Run from the repository root with the project activated:
#
#   julia --project=. examples/demo1/01_plot_data.jl
#
# Run 00_run_inv_GN.jl / 00_run_inv_INR.jl first for the fit and convergence figures; without
# them this only maps the observed data. Figures go to plots/ next to this script.

using GravityGeophysics
using Printf

const HERE = @__DIR__
here(name) = joinpath(HERE, name)
plots_dir  = here("plots")
mkpath(plots_dir)

# the most recently started run_<label>_* folder next to this script, or nothing
function latest_run(label)
    cands = filter(d -> isdir(joinpath(HERE, d)) && startswith(d, label * "_"), readdir(HERE))
    isempty(cands) ? nothing : joinpath(HERE, maximum(cands))
end

# a History.log written by Invert3D, parsed back into plot_inv3d_convergence's NamedTuples
function read_history(path)
    lines = readlines(path)
    rule  = findfirst(startswith("-"), lines)
    names = Symbol.(split(lines[rule + 1]))
    first_row = findnext(startswith("-"), lines, rule + 1) + 1
    rows  = filter(!startswith("-"), lines[first_row:end])
    [NamedTuple{Tuple(names)}(Tuple(parse(Float64, x) for x in split(r))) for r in rows]
end

PlotData3D(here("data.xyz"); output_path = joinpath(plots_dir, "Data.png"))

for (label, tag) in (("run_gn", "GN"), ("run_inr", "INR"))
    dir = latest_run(label)
    dir === nothing && continue
    PlotData3D(here("data.xyz"); predicted_file = joinpath(dir, "data.pred"),
              output_path = joinpath(plots_dir, "DataFit$(tag).png"))
    ctrl = ReadInvCtrl(here("InvCtrl.$tag"))
    plot_inv3d_convergence(read_history(joinpath(dir, "History.log"));
                           output_path = joinpath(plots_dir, "Convergence$(tag).png"),
                           target_rms = ctrl.target_rms)
    @printf("%s fit and convergence from %s\n", tag, dir)
end

println("plots in $plots_dir")
