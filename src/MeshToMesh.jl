# projection of a density model onto another mesh, with its inversion mask and data file
# Author: @pankajkmishra
# The target mesh may be a density model or any other model file (e.g. a resistivity model with
# topography); it keeps its own air (and sea water with `water_log10`), and the source fills the rest by
# nearest cell or trilinear interpolation. Station heights are kept as measured unless `snap_sites`
# drapes them `clearance` m above the target ground. Meshes may differ in rotation (frames turned about the
# frame origin, as rotate_data) and datum: target cell centres are taken into the source frame, and the data
# file is turned to the target's rotation

using Printf
using Statistics
using Dates

#----- cell lookup ---------------------------------------------------------#

# index of the cell containing coordinate q, clamped to the grid: target cells
# outside the source extent take the edge value rather than a fill constant
function _cell_index(edges::AbstractVector{<:Real}, q::Real)
    n = length(edges) - 1
    q <= edges[1] && return 1
    q >= edges[end] && return n
    return clamp(searchsortedlast(edges, q), 1, n)
end

# trilinear weights over the 8 surrounding source cell centers, dropping any
# neighbour that is not usable (air) and renormalising what is left
function _interp_linear(A, cx, cy, cz, x, y, z, usable)
    i = clamp(searchsortedlast(cx, x), 1, length(cx) - 1)
    j = clamp(searchsortedlast(cy, y), 1, length(cy) - 1)
    k = clamp(searchsortedlast(cz, z), 1, length(cz) - 1)
    tx = clamp((x - cx[i]) / (cx[i+1] - cx[i]), 0.0, 1.0)
    ty = clamp((y - cy[j]) / (cy[j+1] - cy[j]), 0.0, 1.0)
    tz = clamp((z - cz[k]) / (cz[k+1] - cz[k]), 0.0, 1.0)
    num = 0.0; den = 0.0
    for (di, wi) in ((0, 1 - tx), (1, tx)),
        (dj, wj) in ((0, 1 - ty), (1, ty)),
        (dk, wk) in ((0, 1 - tz), (1, tz))
        ii, jj, kk = i + di, j + dj, k + dk
        usable[ii, jj, kk] || continue
        w = wi * wj * wk
        num += w * A[ii, jj, kk]; den += w
    end
    return den > 0 ? num / den : NaN
end

# the target mesh: a density model as it is, any other model file through load_mesh
function _target_mesh(path::AbstractString, water_log10)
    _model_read(path).scale == "LINEAR" && return load_model(path)
    load_mesh(path; water_log10 = water_log10)
end

#----- the projection ------------------------------------------------------#

"""
    project_model_data(src_model, src_data, dst_mesh; kwargs...)

Resample the density model `src_model` onto the grid of `dst_mesh`, write the
inversion mask of the target mesh and a copy of the data file `src_data`.

Works in either direction (coarse to fine or fine to coarse) and does not require the two
meshes to share extent, origin, rotation or datum. Target cells outside the source
extent take the nearest source value. The data file is written in the target's frame
(`rotate_data` to the target rotation).

The target mesh keeps **its own** topography: cells that are air in `dst_mesh` (and sea water
below `water_log10` in a resistivity mesh) stay air and are never filled from the source.

By default the three files land in a fresh timestamped directory beside the target mesh and
are named after their inputs, so nothing is ever overwritten:

    <dst dir>/projected/20260726_150322/modelC_to_modelF.den
                                        /dataC-W.dat
                                        /modelC_to_modelF.mask

Keywords
- `outdir`      override the timestamped directory
- `stem`        override the `<src>_to_<dst>` name shared by the model and mask
- `method`      `:nearest` (default, preserves contrasts) or `:linear`
- `water_log10` log10 ρ below which a resistivity target cell is sea water (air here), `NaN` for none
- `snap_sites`  drape the stations `clearance` m above the target ground (default false: airborne
  heights are measured and kept)
- `clearance`   height above the ground of draped stations (m)

Returns a NamedTuple of the three output paths plus diagnostics.
"""
function project_model_data(src_model::AbstractString, src_data::AbstractString,
                            dst_mesh::AbstractString;
                            outdir::AbstractString = joinpath(dirname(dst_mesh),
                                "projected", Dates.format(now(), "yyyymmdd_HHMMSS")),
                            stem::AbstractString = string(
                                first(splitext(basename(src_model))), "_to_",
                                first(splitext(basename(dst_mesh)))),
                            method::Symbol = :nearest,
                            water_log10::Real = NaN,
                            snap_sites::Bool = false,
                            clearance::Real = 1.0)

    method in (:nearest, :linear) || error("method must be :nearest or :linear")

    src = load_model(src_model)
    dst = _target_mesh(dst_mesh, water_log10)

    src_air = isnan.(src.A)
    dst_air = isnan.(dst.A)
    usable = .!src_air                       # never pull a value out of the air

    # source z in the target's frame: the same elevation, z = datum − elevation
    shift = dst.datum - src.datum
    sx, sy, sz = src.x, src.y, src.z .+ shift
    scz = src.cz .+ shift
    # target frame → source frame: turned by the difference of the rotations
    φ = src.rotation - dst.rotation
    c, s = cosd(φ), sind(φ)

    A = similar(dst.A)
    n_outside = 0
    n_kept = 0
    @inbounds for k in 1:dst.nz, j in 1:dst.ny, i in 1:dst.nx
        # the target mesh's own surface wins
        if dst_air[i, j, k]
            A[i, j, k] = NaN; n_kept += 1; continue
        end
        x, y, z = c * dst.cx[i] + s * dst.cy[j], -s * dst.cx[i] + c * dst.cy[j], dst.cz[k]
        outside = x < sx[1] || x > sx[end] || y < sy[1] || y > sy[end] ||
                  z < sz[1] || z > sz[end]
        outside && (n_outside += 1)
        v = if method === :linear && !outside
            _interp_linear(src.A, src.cx, src.cy, scz, x, y, z, usable)
        else
            NaN
        end
        if !isfinite(v)
            # nearest cell, walking up the column if we land in source air
            ii = _cell_index(sx, x); jj = _cell_index(sy, y); kk = _cell_index(sz, z)
            while kk >= 1 && src_air[ii, jj, kk]; kk -= 1; end
            kk < 1 && (kk = findfirst(!, view(src_air, ii, jj, :)))
            v = kk === nothing ? median(filter(isfinite, src.A)) : src.A[ii, jj, kk]
        end
        A[i, j, k] = v
    end

    mkpath(outdir)
    # -W marks the data file as written for the target mesh
    dstem, dext = splitext(basename(src_data))
    model_out = joinpath(outdir, stem * ".den")
    data_out  = joinpath(outdir, dstem * "-W" * (isempty(dext) ? ".dat" : dext))
    mask_out  = joinpath(outdir, stem * ".mask")

    projected = Model(dst.dx, dst.dy, dst.dz, A; origin = dst.origin, rotation = dst.rotation, datum = dst.datum)
    write_model(model_out, projected; comment = "$(basename(src_model)) projected onto $(basename(dst_mesh))")
    write_mask(mask_out, projected, .!dst_air)

    #---------- stations onto the target mesh ----------
    d = load_data_grav(src_data)
    isapprox(rem(d.azimuth - dst.rotation, 360, RoundNearest), 0; atol = 1e-9) || (d = rotate_data(d, dst.rotation - d.azimuth))
    z_old = copy(d.z)
    off = count(n -> d.x[n] < dst.x[1] || d.x[n] > dst.x[end] || d.y[n] < dst.y[1] || d.y[n] > dst.y[end], 1:ndata(d))
    snap_sites && drape!(d, projected, clearance)
    write_data_grav(data_out, d; comment = snap_sites ? @sprintf("stations draped %.3g m above %s", clearance, basename(dst_mesh)) :
                                         "for $(basename(dst_mesh))")
    shifts = d.z .- z_old

    return (model = model_out, data = data_out, mask = mask_out,
            src_grid = (src.nx, src.ny, src.nz),
            dst_grid = (dst.nx, dst.ny, dst.nz),
            cells_outside_source = n_outside,
            surface_cells_kept = n_kept,
            n_air = count(dst_air),
            sites = ndata(d), sites_off_grid = off,
            elev_shift = isempty(shifts) ? (0.0, 0.0, 0.0) :
                (minimum(shifts), median(shifts), maximum(shifts)))
end

"""
    report(r)

Print the NamedTuple returned by [`project_model_data`](@ref).
"""
function report(r::NamedTuple)
    @printf("  source grid %s  ->  target grid %s\n", r.src_grid, r.dst_grid)
    @printf("  target cells outside the source extent : %d (filled from nearest)\n",
            r.cells_outside_source)
    @printf("  target air cells kept as air           : %d\n", r.surface_cells_kept)
    @printf("  stations %d (%d outside the target grid)\n", r.sites, r.sites_off_grid)
    @printf("  station height shift min/median/max    : %.1f / %.1f / %.1f m\n",
            r.elev_shift...)
    println("  wrote ", r.model)
    println("        ", r.data)
    println("        ", r.mask)
    return nothing
end

"""
    MeshToMesh(src_model, src_data, dst_mesh; outdir="", stem="", method=:nearest, water_log10=NaN,
               snap_sites=false, clearance=1.0)

Project a density model and its data file onto another mesh and report it; see
[`project_model_data`](@ref).
"""
function MeshToMesh(src_model::AbstractString, src_data::AbstractString,
                    dst_mesh::AbstractString;
                    outdir::AbstractString = "",
                    stem::AbstractString = "",
                    method::Symbol = :nearest,
                    water_log10::Real = NaN,
                    snap_sites::Bool = false,
                    clearance::Real = 1.0)

    isempty(src_model) && error("source model path is required")
    isfile(src_model) || error("Source model not found: $src_model")
    isempty(src_data) && error("source data path is required")
    isfile(src_data) || error("Source data not found: $src_data")
    isempty(dst_mesh) && error("target mesh path is required")
    isfile(dst_mesh) || error("Target mesh not found: $dst_mesh")

    kwargs = Dict{Symbol, Any}(
        :method => method,
        :water_log10 => water_log10,
        :snap_sites => snap_sites,
        :clearance => clearance,
    )
    isempty(outdir) || (kwargs[:outdir] = outdir)
    isempty(stem) || (kwargs[:stem] = stem)

    r = project_model_data(src_model, src_data, dst_mesh; kwargs...)
    report(r)
    return r
end
