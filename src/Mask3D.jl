# topography extraction, air masks and inversion masks for density models
# Author: @pankajkmishra
# air cells are the NaN tag the loader applies above 1e15; the topography file stores the base of the
# tagged column run at physical cell centers in the model frame (x, y and z carry the model origin), so
# it can be re-applied to any grid on the same coordinate frame and plotted directly on a model section.
# An inversion mask is a LINEAR model file on the model mesh: 1 = variable, 0 = fixed at the starting
# value, 9 = water (lakes, as MTGeophysics; 1000 − ρ_ref in the model, fixed, outside the smoothing); air cells are
# fixed whatever the mask says. `load_mesh` still makes the sea of a resistivity model air

using Dates
using Printf

#---------- topography and air masks ----------

# base of the contiguous run of tagged cells at the top of each column, as depth
# below the model top and as z in the model frame
function _tagged_column_base(m::Model, istagged)
    zbot = cumsum(m.dz)
    xs = Float64[]; ys = Float64[]; ds = Float64[]; zs = Float64[]
    n_stray = 0
    for j in 1:m.ny, i in 1:m.nx
        klast = 0
        for k in 1:m.nz
            istagged(m.A[i, j, k]) || break
            klast = k
        end
        n_stray += count(k -> istagged(m.A[i, j, k]), klast+1:m.nz)
        if klast > 0
            push!(xs, m.cx[i]); push!(ys, m.cy[j])
            push!(ds, zbot[klast]); push!(zs, m.z[klast+1])
        end
    end
    return (x=xs, y=ys, depth=ds, z=zs), n_stray
end

function _write_surface_table(path::AbstractString, cols, title::AbstractString,
                              base_label::AbstractString, meta::AbstractString,
                              origin::AbstractVector{<:Real})
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# GravityGeophysics ", title, " — ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "# ", meta, "  n_columns=", length(cols.depth))
        @printf(io, "# origin_x_m=%.3f  origin_y_m=%.3f  origin_z_m=%.3f\n",
                origin[1], origin[2], origin[3])
        println(io, "# x_center_m  y_center_m  ", base_label, "  base_z_m")
        hasz = hasproperty(cols, :z)
        for t in eachindex(cols.depth)
            @printf(io, "%15.3f %15.3f %12.3f %12.3f\n",
                    cols.x[t], cols.y[t], cols.depth[t], hasz ? cols.z[t] : NaN)
        end
    end
    return path
end

function _read_surface_table(path::AbstractString)
    xs = Float64[]; ys = Float64[]; ds = Float64[]; zs = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        p = split(s)
        push!(xs, parse(Float64, p[1]))
        push!(ys, parse(Float64, p[2]))
        push!(ds, parse(Float64, p[3]))
        push!(zs, length(p) >= 4 ? parse(Float64, p[4]) : NaN)
    end
    return (x=xs, y=ys, depth=ds, z=zs)
end

# depth is measured from the model top while m.cz carries origin[3], so the cell
# test uses cumsum(dz); comparing centers leaves half a cell of slack, which the
# 3-decimal file format needs
function _mask_above_base(m::Model, cols)
    mask = falses(m.nx, m.ny, m.nz)
    zcen = cumsum(m.dz) .- m.dz ./ 2
    for t in eachindex(cols.depth)
        i = argmin(abs.(m.cx .- cols.x[t]))
        j = argmin(abs.(m.cy .- cols.y[t]))
        for k in 1:m.nz
            zcen[k] < cols.depth[t] || break
            mask[i, j, k] = true
        end
    end
    return mask
end

"""
    extract_topography(m) -> (x, y, depth, z)

Ground-surface depth (m, positive down) per column of a model that carries topographic
air. Air is the `NaN` tag applied by [`load_model`](@ref) above 1e15; depth is the bottom
edge of the contiguous air run at the top of the column, `z` the same interface in the
model frame. Air below the ground surface triggers a warning and is ignored.
"""
function extract_topography(m::Model)
    cols, n_stray = _tagged_column_base(m, isnan)
    n_stray > 0 && @warn "extract_topography: $n_stray air cells below the ground surface ignored"
    return cols
end

"""
    write_topography(path, topo; origin=[0,0,0])

Write a topography table (`x_center_m  y_center_m  ground_depth_m  base_z_m`, one row per
column with air) as produced by [`extract_topography`](@ref). `x`/`y` are cell centers and
`base_z_m` the ground surface in the model frame, both carrying `origin`, so the surface
plots directly on a model section.
"""
write_topography(path::AbstractString, topo;
                 origin::AbstractVector{<:Real}=[0.0, 0.0, 0.0]) =
    _write_surface_table(path, topo, "topography", "ground_depth_m",
                         "air_tag=NaN(value>1e15)", origin)

"""
    read_topography(path) -> (x, y, depth, z)

Read a topography table written by [`write_topography`](@ref).
"""
read_topography(path::AbstractString) = _read_surface_table(path)

"""
    air_mask_from_topography(m, topo) -> BitArray{3}

Full-grid mask of air cells: for each column with air (matched to the nearest
cell-center column of `m`), cells whose centers lie above the ground surface.
"""
air_mask_from_topography(m::Model, topo) = _mask_above_base(m, topo)

"""
    air_mask_from_model(m) -> BitArray{3}

Full-grid mask of topographic air, i.e. the cells the loader tagged `NaN`.
"""
air_mask_from_model(m::Model) = BitArray(isnan.(m.A))

#---------- inversion masks ----------

"""
    load_mask(path, m) -> BitArray{3}

Variable cells of `m`: mask value 1 and not air. The mask must be on the same mesh as `m`; water (9) and
fixed (0) cells are not variable.
"""
load_mask(path::AbstractString, m::Model) = (_mask_values(path, m) .== 1) .& active_cells(m)

"""
    water_mask(path, m) -> BitArray{3}

Water cells of `m` (lakes): mask value 9 in the mask file at `path`.
"""
water_mask(path::AbstractString, m::Model) = (_mask_values(path, m) .== 9) .& active_cells(m)

function _mask_values(path::AbstractString, m::Model)
    r = _model_read(path)
    r.scale == "LINEAR" || throw(ArgumentError("mask files are LINEAR; $path is $(r.scale)"))
    mask = Model(r.dx, r.dy, r.dz, r.A; origin = r.origin, rotation = r.rotation, datum = r.datum)
    same_mesh(m, mask) || throw(ArgumentError("mask and model meshes differ"))
    all(v -> isnan(v) || v == 0 || v == 1 || v == 9, mask.A) ||
        throw(ArgumentError("mask values must be 0 (fixed), 1 (variable) or 9 (water)"))
    mask.A
end

"""
    write_mask(path, m, variable; water=falses)

Write `variable` (true = inverted) as a mask on the mesh of `m`, with `water` cells as 9; air cells are
written as air.
"""
function write_mask(path::AbstractString, m::Model, variable::AbstractArray{Bool,3};
                    water::AbstractArray{Bool,3} = falses(size(m.A)))
    size(variable) == size(water) == size(m.A) || throw(ArgumentError("mask size differs from mesh"))
    any(variable .& water) && throw(ArgumentError("a water cell cannot be variable"))
    A = [isnan(m.A[I]) ? NaN : water[I] ? 9.0 : Float64(variable[I]) for I in CartesianIndices(m.A)]
    mask = Model(m.dx, m.dy, m.dz, A; origin = m.origin, rotation = m.rotation, datum = m.datum)
    write_model(path, mask; comment = "mask: 1 variable, 0 fixed, 9 water")
end
