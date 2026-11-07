# Case figures: Observed and predicted gz of the newest run of each algorithm, normalised residuals and their histograms
# Author: @pankajkmishra
# Map view: easting horizontal, northing vertical, in km about the frame centre.
# Runs: the newest run_*/ of each algorithm in the case directory (from its 0<algorithm>.log).
# Writes figures/<case>_2_fit_predicted.png and <case>_2_fit_residuals.png in the case directory.
# Usage: julia --project=. helpers/plot_fit_3D.jl [case_dir]   (default examples/synthetic)

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

d = load_data_grav(joinpath(DIR, "data.dat"))
start = load_model(joinpath(DIR, "model.start"))

# finished runs: a run log 0lsqr.log with its result ("# key = value" lines)
logof(r) = (f = filter(f -> occursin(r"^0lsqr\.log$", f), readdir(r)); isempty(f) ? "" : joinpath(r, only(f)))
result(r) = Dict(strip(a) => strip(b) for (a, b) in (split(l[3:end], " = "; limit = 2) for l in readlines(logof(r)) if startswith(l, "# ") && occursin(" = ", l)))
runs = sort(filter(r -> isfile(logof(r)) && haskey(result(r), "rms"), [joinpath(DIR, r) for r in readdir(DIR) if startswith(r, "run_") && isdir(joinpath(DIR, r))]))
latest = Dict(uppercase(basename(logof(r))[2:end-4]) => r for r in runs)   # sorted by time: the newest wins
algs = filter(a -> haskey(latest, a), ["LSQR"])
isempty(algs) && error("no runs in $DIR")
rmsof(a) = parse(Float64, result(latest[a])["rms"])

pred = Dict(a => load_data_grav(joinpath(latest[a], "data.pred")).value for a in algs)


east, north = d.y ./ 1000, d.x ./ 1000
# data maps: inverse-distance grid (6 nearest stations) on a grid of half the along-line spacing, blanked farther
# than 0.75 line spacings from any station; small black dots mark the stations
line = first.(split.(d.code, ':'))
near_any(i, same) = minimum((hypot(east[i] - east[j], north[i] - north[j]) for j in eachindex(east) if j != i && (line[j] == line[i]) == same); init = Inf)
along = median(filter(isfinite, [near_any(i, true) for i in eachindex(east)]))
across = length(unique(line)) > 1 ? median(filter(isfinite, [near_any(i, false) for i in eachindex(east)])) : along
ge = range(extrema(east)...; step = along / 2); gn = range(extrema(north)...; step = along / 2)
gw = [begin
          r2 = (east .- e) .^ 2 .+ (north .- n) .^ 2
          k = partialsortperm(r2, 1:6)
          (k, 1 ./ max.(r2[k], 1e-6), sqrt(r2[k[1]]) <= 0.75across)
      end for e in ge, n in gn]
grid(v) = [w[3] ? sum(w[2] .* v[w[1]]) / sum(w[2]) : NaN for w in gw]
function station_map!(pos, v, title; lim, label, ylab = true, bar = true)
    ax = Axis(pos[1, 1]; title, xlabel = "Easting (km)", ylabel = ylab ? "Northing (km)" : "", aspect = DataAspect(), AX...)
    hm = heatmap!(ax, ge, gn, grid(v); colormap = CMAP, colorrange = lim, nan_color = :white)
    scatter!(ax, east, north; color = :black, markersize = 1.5)
    bar && Colorbar(pos[1, 2], hm; label, width = 18)
    ax, hm
end
dlim = (pct(d.value, 0.005), pct(d.value, 0.995))



nc = length(algs) + 1
fig = Figure(size = (560nc + 160, 560), fontsize = 22)
station_map!(fig[1, 1], d.value, "Observed"; lim = dlim, label = "", bar = false)
for (c, a) in enumerate(algs)
    station_map!(fig[1, c+1], pred[a], @sprintf("%s predicted (rms %.2f)", a, rmsof(a)); lim = dlim, label = "",
                 ylab = false, bar = false)
end
Colorbar(fig[1, nc+1]; colormap = CMAP, limits = dlim, label = "gz (mGal)", width = 18)
labels!(fig, [(1, c) for c in 1:nc])
save(out(2, "fit_predicted"), fig; px_per_unit = 3)

fig = Figure(size = (560nc + 160, 560), fontsize = 22)
for (c, a) in enumerate(algs)
    station_map!(fig[1, c], (d.value .- pred[a]) ./ d.error, "$a normalised residual"; lim = (-3.0, 3.0), label = "",
                 ylab = c == 1, bar = false)
end
Colorbar(fig[1, nc]; colormap = CMAP, limits = (-3.0, 3.0), label = "Residual / error", width = 18)
ax = Axis(fig[1, nc+1]; title = "Normalised residuals", xlabel = "(observed − predicted) / error", ylabel = "Stations", AX...)
for (c, a) in enumerate(algs)
    hist!(ax, clamp.((d.value .- pred[a]) ./ d.error, -6, 6); bins = range(-6, 6; length = 49),
          color = (Makie.wong_colors()[c], 0.55), strokewidth = 1.5, strokecolor = Makie.wong_colors()[c], label = a)
end
axislegend(ax; framevisible = false, position = :rt)
labels!(fig, [[(1, c) for c in 1:length(algs)]; [(1, nc + 1)]])
save(out(2, "fit_residuals"), fig; px_per_unit = 3)

foreach(w -> println("figure → ", relpath(out(2, w), DIR)), ("fit_predicted", "fit_residuals"))
