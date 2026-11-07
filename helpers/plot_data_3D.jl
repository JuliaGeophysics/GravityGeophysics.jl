# Case figures: Observed Bouguer anomaly (and station elevation for surveys not on one level)
# Author: @pankajkmishra
# Map view: easting horizontal, northing vertical, in km about the frame centre.
# Runs: the newest run_*/ of each algorithm in the case directory (from its 0<algorithm>.log).
# Writes figures/<case>_1_data.png in the case directory.
# Usage: julia --project=. helpers/plot_data_3D.jl [case_dir]   (default examples/synthetic)

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



# top of the active cells in every column of the core (NaN = air)
ix, iy = (start.npad[1] + 1):(start.nx - start.npad[1]), (start.npad[2] + 1):(start.ny - start.npad[2])
topk = [something(findfirst(k -> !isnan(start.A[i, j, k]), 1:start.nz), start.nz) for i in 1:start.nx, j in 1:start.ny]
ground = start.datum .- start.z[topk]                                   # elevation of the top active cell (m)
field = maximum(-d.z) - minimum(-d.z) > 1                              # stations not on one level

fig = Figure(size = (field ? 1700 : 900, 800), fontsize = 22)
station_map!(fig[1, 1], d.value, "Observed Bouguer anomaly"; lim = dlim, label = "gz (mGal)")
if field
    station_map!(fig[1, 2], -d.z, "Station elevation"; lim = extrema(-d.z), label = "Elevation (m)", ylab = false)
    labels!(fig, [(1, 1), (1, 2)])
end
save(out(1, "data"), fig; px_per_unit = 3)

println("figure → ", relpath(out(1, "data"), DIR))
