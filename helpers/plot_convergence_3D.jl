# Case figures: LSQR rms per iteration
# Author: @pankajkmishra
# Runs: the newest run_*/ of each algorithm in the case directory (from its 0<algorithm>.log).
# Writes figures/<case>_4_convergence.png in the case directory.
# Usage: julia --project=. helpers/plot_convergence_3D.jl [case_dir]   (default examples/synthetic)

using GravityGeophysics
using CairoMakie
using Printf
using Statistics
CairoMakie.activate!()   # figures go to files with Cairo: the package also loads GLMakie, which needs a display

#----- style -----------------------------------------------------------------#

const AX = (spinewidth = 2.5, xtickwidth = 2.5, ytickwidth = 2.5, xticksize = 8, yticksize = 8,
            xgridvisible = false, ygridvisible = false, titlefont = :regular, titlesize = 22)
const CMAP = COLORMAP
labels!(fig, pos) = for (n, p) in enumerate(pos)
    Label(fig[p..., TopLeft()], "($('a' + n - 1))"; fontsize = 24, padding = (0, 14, 18, 0), halign = :left)
end
pct(v, q) = (s = sort(filter(isfinite, v)); s[clamp(ceil(Int, q * length(s)), 1, length(s))])


#----- case ------------------------------------------------------------------#

const DIR = abspath(isempty(ARGS) ? joinpath(dirname(@__DIR__), "examples", "synthetic") : ARGS[1])
const NAME = basename(rstrip(DIR, '/'))
const MAX_DEPTH = 2000.0                                                # shown below the highest ground (m)
out(k, what) = joinpath(mkpath(joinpath(DIR, "figures")), "$(NAME)_$(k)_$(what).png")

# finished runs: a run log 0lsqr.log with its result ("# key = value" lines)
logof(r) = (f = filter(f -> occursin(r"^0lsqr\.log$", f), readdir(r)); isempty(f) ? "" : joinpath(r, only(f)))
result(r) = Dict(strip(a) => strip(b) for (a, b) in (split(l[3:end], " = "; limit = 2) for l in readlines(logof(r)) if startswith(l, "# ") && occursin(" = ", l)))
runs = sort(filter(r -> isfile(logof(r)) && haskey(result(r), "rms"), [joinpath(DIR, r) for r in readdir(DIR) if startswith(r, "run_") && isdir(joinpath(DIR, r))]))
latest = Dict(uppercase(basename(logof(r))[2:end-4]) => r for r in runs)   # sorted by time: the newest wins
algs = filter(a -> haskey(latest, a), ["LSQR"])
isempty(algs) && error("no runs in $DIR")
rmsof(a) = parse(Float64, result(latest[a])["rms"])


const RTICKS = [0.5, 1, 2, 3, 5, 10, 20, 50, 100]
rows(f) = [parse.(Float64, split(l)) for l in readlines(f) if !isempty(l) && tryparse(Int, first(split(l))) !== nothing]   # log rows
fig = Figure(size = (900, 650), fontsize = 22)
det = algs                                                              # LSQR: rms per iteration
if !isempty(det)
    ax = Axis(fig[1, 1]; title = join(det, " and "), xlabel = "Iteration", ylabel = "rms", yscale = log10, yticks = RTICKS, AX...)
    for (k, a) in enumerate(det)
        t = rows(logof(latest[a]))                                      # Iter Objective Chi2 RMS …
        scatterlines!(ax, first.(t), getindex.(t, 4); color = (:black, :gray55, :steelblue, :darkorange)[k], linewidth = 2.5,
                      markersize = 12, label = a)
    end
    hlines!(ax, 1.0; color = :gray40, linestyle = :dash, linewidth = 2)
    length(det) > 1 && axislegend(ax; framevisible = false)
end
labels!(fig, [(1, 1)])
save(out(4, "convergence"), fig; px_per_unit = 3)

println("figure → ", relpath(out(4, "convergence"), DIR))
