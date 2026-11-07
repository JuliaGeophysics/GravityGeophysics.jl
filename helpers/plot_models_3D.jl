# Case figures: depth slice and east-west / north-south sections through the centre of the body (the true model's for
# a synthetic, else the first recovered model's); a map and a sections figure per model: true and LSQR
# Author: @pankajkmishra
# Map view: easting horizontal, northing vertical, in km about the frame centre. Sections at true scale, depth (km)
# positive down from 0 at the top of the topography (the highest active cell of the core); thick black line: top of
# the active cells as the mesh discretises it, air in white, water (mask 9: lakes) in blue; thin black outlines: true bodies (synthetic only),
# along cell edges. Recovered models share one colour range. Only the core (no padding) is shown.
# Runs: the newest run_*/ of each algorithm in the case directory (from its 0<algorithm>.log).
# Writes figures/<case>_3_models_<n>_<model>_map.png and _sections.png (two sections stacked) in the case directory.
# The body: the connected cells (faces shared) of at least half the largest Δρ in the core that hold the largest; the
# slice and the sections pass through its centre of mass. Sections reach 1.5 × the depth of its base, or max_depth_m.
# Usage: julia --project=. helpers/plot_models_3D.jl [case_dir [max_depth_m]]   (default examples/synthetic)

using GravityGeophysics
using CairoMakie
using Printf
CairoMakie.activate!()   # figures go to files with Cairo: the package also loads GLMakie, which needs a display

#----- style -----------------------------------------------------------------#

const AX = (spinewidth = 2.5, xtickwidth = 2.5, ytickwidth = 2.5, xticksize = 8, yticksize = 8,
            xgridvisible = false, ygridvisible = false, titlefont = :regular, titlesize = 22)
const CMAP = COLORMAP
labels!(pos) = for (n, (g, r, c)) in enumerate(pos)
    Label(g[r, c, TopLeft()], "($('a' + n - 1))"; fontsize = 24, padding = (0, 14, 18, 0), halign = :left)
end

#----- case ------------------------------------------------------------------#

const DIR = abspath(isempty(ARGS) ? joinpath(dirname(@__DIR__), "examples", "synthetic") : ARGS[1])
const NAME = basename(rstrip(DIR, '/'))
const MAX_DEPTH = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : nothing # shown below the lowest ground of the core (m)
const W = 460                                                           # width of a map panel (px)
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

#----- models and geometry ---------------------------------------------------#

panels = Tuple{String,Model}[]
isfile(joinpath(DIR, "model.true")) && push!(panels, ("True", load_model(joinpath(DIR, "model.true"))))
for a in algs
    push!(panels, (@sprintf("%s (rms %.2f)", a, rmsof(a)), load_model(joinpath(latest[a], "model.den"))))
end
truth = first(panels)[1] == "True"
rec = truth ? panels[2:end] : panels

ix, iy = (start.npad[1] + 1):(start.nx - start.npad[1]), (start.npad[2] + 1):(start.ny - start.npad[2])
topk = [something(findfirst(k -> !isnan(start.A[i, j, k]), 1:start.nz), start.nz) for i in 1:start.nx, j in 1:start.ny]
z0 = minimum(start.z[topk[ix, iy]])                                    # top of the topography: depth 0

# the body of the reference model: cells of the core of at least half its largest Δρ, connected to the largest
ref = last(first(panels))
C = fill(-Inf, size(ref.A)); C[ix, iy, :] .= replace(ref.A[ix, iy, :], NaN => -Inf)
seed = argmax(C)
inbody = C .>= C[seed] / 2
body = falses(size(C)); body[seed] = true
stack = [seed]
while !isempty(stack)
    c = pop!(stack)
    for s in (CartesianIndex(1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, 0, 1))
        for n in (c + s, c - s)
            checkbounds(Bool, C, n) && inbody[n] && !body[n] && (body[n] = true; push!(stack, n))
        end
    end
end
cells = findall(body)
vol(c) = start.dx[c[1]] * start.dy[c[2]] * start.dz[c[3]]
centre(axis, c) = (start.cx, start.cy, start.cz)[axis][c[axis]]
com = [sum(vol(c) * centre(a, c) for c in cells) / sum(vol, cells) for a in 1:3]
I = (argmin(abs.(start.cx .- com[1])), argmin(abs.(start.cy .- com[2])), argmin(abs.(start.cz .- com[3])))
zbase = maximum(start.z[c[3] + 1] for c in cells) - z0                  # depth of the base of the body
reach = something(MAX_DEPTH, 1.5 * zbase)
kz = searchsortedlast(start.z, z0):min(start.nz, searchsortedlast(start.z, maximum(start.z[topk[ix, iy]]) + reach))
edges(e, r) = e[first(r):last(r)+1]
he, hn, ve = edges(start.y, iy) ./ 1000, edges(start.x, ix) ./ 1000, (edges(start.z, kz) .- z0) ./ 1000
@printf("body: %d cells, centre northing %.0f m, easting %.0f m, depth %.0f m; base %.0f m\n", length(cells),
        com[1], com[2], com[3] - z0, zbase)
depth(k) = (start.z[k] - z0, start.z[k+1] - z0)
@printf("slice point: northing %.0f m, easting %.0f m, depth slice %.0f–%.0f m\n", start.cx[I[1]], start.cy[I[2]], depth(I[3])...)

staircase(h, t) = (vcat([[h[j], h[j+1]] for j in eachindex(t)]...), vcat([[t[j], t[j]] for j in eachindex(t)]...))
topo_ew = staircase(he, (start.z[topk[I[1], iy]] .- z0) ./ 1000)
topo_ns = staircase(hn, (start.z[topk[ix, I[2]]] .- z0) ./ 1000)
planes(m) = (permutedims(m.A[ix, iy, I[3]]), m.A[I[1], iy, kz], m.A[ix, I[2], kz])   # (east × north), (east × depth), (north × depth)
maskfile = joinpath(DIR, "model.mask")
wet = isfile(maskfile) ? water_mask(maskfile, start) : falses(size(start.A))       # lakes, drawn in blue
wetplanes = (wet[I[1], iy, kz], wet[ix, I[2], kz])

# outline of a mask (horizontal × vertical cells) along cell edges
function outline(mask, h, v)
    seg = Point2f[]
    inside(a, b) = checkbounds(Bool, mask, a, b) && mask[a, b]
    for a in axes(mask, 1), b in axes(mask, 2)
        mask[a, b] || continue
        inside(a - 1, b) || push!(seg, Point2f(h[a], v[b]), Point2f(h[a], v[b+1]))
        inside(a + 1, b) || push!(seg, Point2f(h[a+1], v[b]), Point2f(h[a+1], v[b+1]))
        inside(a, b - 1) || push!(seg, Point2f(h[a], v[b]), Point2f(h[a+1], v[b]))
        inside(a, b + 1) || push!(seg, Point2f(h[a], v[b+1]), Point2f(h[a+1], v[b+1]))
    end
    seg
end
bodies = nothing
if truth
    T = planes(last(first(panels)))
    bodies = (outline(T[1] .> 0, he, hn), outline(replace(T[2], NaN => 0.0) .> 0, he, ve),
              outline(replace(T[3], NaN => 0.0) .> 0, hn, ve))
end
fin(A) = filter(isfinite, vec(A))
tmax = truth ? maximum(maximum(fin(A)) for A in planes(last(first(panels)))) : 0.0
cmax = maximum(maximum(fin(A)) for (_, m) in rec for A in planes(m))
# density contrasts may be negative (water, light rock): the range reaches down to the lowest value
tmin = truth ? min(0.0, minimum(minimum(fin(A)) for A in planes(last(first(panels))))) : 0.0
cmin = min(0.0, minimum(minimum(fin(A)) for (_, m) in rec for A in planes(m)))

#----- figures: a map and a sections figure per model -------------------------#

px = W / (he[end] - he[1])                                              # px per km, every panel
tag(t) = lowercase(replace(first(split(t, " (")), r"[^A-Za-z0-9]+" => "_"))
files = String[]
for (r, (title, m)) in enumerate(panels)
    istrue = truth && r == 1
    lim = istrue ? (tmin, tmax) : (cmin, cmax)
    clabel = (istrue ? "true" : "recovered") * " Δρ (kg/m³)"
    P = planes(m)
    name = "models_$(r)_$(tag(title))"

    # map: depth slice with the stations and the section lines
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1, 1]; title = title * @sprintf(", depth %.0f–%.0f m", depth(I[3])...), xlabel = "Easting (km)",
              ylabel = "Northing (km)", width = W, height = px * (hn[end] - hn[1]), AX...)
    hm = heatmap!(ax, he, hn, P[1]; colormap = CMAP, colorrange = lim, nan_color = :white)
    scatter!(ax, d.y ./ 1000, d.x ./ 1000; color = :black, markersize = 2)
    hlines!(ax, start.cx[I[1]] / 1000; color = :black, linestyle = :dash, linewidth = 1.5)
    vlines!(ax, start.cy[I[2]] / 1000; color = :black, linestyle = :dash, linewidth = 1.5)
    truth && linesegments!(ax, bodies[1]; color = :black, linewidth = 1.5)
    xlims!(ax, he[1], he[end]); ylims!(ax, hn[1], hn[end])
    Colorbar(fig[1, 2], hm; label = clabel, width = 18, height = Relative(0.8))
    resize_to_layout!(fig)
    save(out(3, name * "_map"), fig; px_per_unit = 3)
    push!(files, out(3, name * "_map"))

    # sections: east-west above north-south, through the slice point
    fig = Figure(fontsize = 22)
    local hs
    for (k, (h, A, t, lab, xl)) in enumerate(((he, P[2], topo_ew, @sprintf("%s, northing %.2f km", title, start.cx[I[1]] / 1000), "Easting (km)"),
                                              (hn, P[3], topo_ns, @sprintf("%s, easting %.2f km", title, start.cy[I[2]] / 1000), "Northing (km)")))
        sx = Axis(fig[k, 1]; title = lab, xlabel = xl, ylabel = "Depth (km)", yreversed = true,
                  width = px * (h[end] - h[1]), height = px * (ve[end] - ve[1]), AX...)
        hs = heatmap!(sx, h, ve, A; colormap = CMAP, colorrange = lim, nan_color = :white)
        any(wetplanes[k]) && heatmap!(sx, h, ve, map(w -> w ? 1.0 : NaN, wetplanes[k]); colormap = [:dodgerblue, :dodgerblue],
                                      colorrange = (0, 2), nan_color = :transparent)
        truth && linesegments!(sx, bodies[k+1]; color = :black, linewidth = 1.5)
        lines!(sx, t...; color = :black, linewidth = 3.5)
        xlims!(sx, h[1], h[end]); ylims!(sx, ve[end], ve[1])
    end
    Colorbar(fig[1:2, 2], hs; label = clabel, width = 18, height = Relative(0.8))
    labels!([(fig.layout, 1, 1), (fig.layout, 2, 1)])
    resize_to_layout!(fig)
    save(out(3, name * "_sections"), fig; px_per_unit = 3)
    push!(files, out(3, name * "_sections"))
end

foreach(f -> println("figure → ", relpath(f, DIR)), files)
