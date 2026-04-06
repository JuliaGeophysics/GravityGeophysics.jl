function edges_from_centers(c::AbstractVector{<:Real})
    n = length(c)
    n == 0 && return Float64[]
    n == 1 && return [Float64(c[1]) - 0.5, Float64(c[1]) + 0.5]
    mids = (Float64.(c[1:end-1]) .+ Float64.(c[2:end])) ./ 2
    first_edge = Float64(c[1]) - (Float64(c[2]) - Float64(c[1])) / 2
    last_edge = Float64(c[end]) + (Float64(c[end]) - Float64(c[end - 1])) / 2
    return vcat(first_edge, mids, last_edge)
end

function compute_colorrange(values::AbstractArray)
    finite_values = Float64.(values[isfinite.(values)])
    isempty(finite_values) && return (0.0, 1.0)
    qlo, qhi = quantile(vec(finite_values), (0.02, 0.98))
    lo = min(qlo, qhi)
    hi = max(qlo, qhi)
    if lo == hi
        delta = max(1e-12, 1e-6 * abs(lo))
        lo -= delta
        hi += delta
    end
    return lo, hi
end

function bracket_index_and_weight(vals::AbstractVector{<:Real}, query::Real)
    n = length(vals)
    n <= 1 && return 1, 1, 0.0
    q = Float64(query)
    if q <= vals[1]
        return 1, 2, 0.0
    elseif q >= vals[end]
        return n - 1, n, 1.0
    end
    i1 = searchsortedfirst(vals, q)
    i0 = i1 - 1
    v0 = Float64(vals[i0])
    v1 = Float64(vals[i1])
    if v0 == v1
        return i0, i1, 0.0
    end
    w = (q - v0) / (v1 - v0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

function sample_polyline(points::Vector{Tuple{Float64, Float64}}, nsamp::Int)
    length(points) < 2 && return Float64[], Float64[], Float64[]

    seglen = Float64[0.0]
    for i in 2:length(points)
        push!(seglen, seglen[end] + hypot(points[i][1] - points[i - 1][1], points[i][2] - points[i - 1][2]))
    end

    total = seglen[end]
    if total <= 0
        return fill(points[1][1], nsamp), fill(points[1][2], nsamp), zeros(nsamp)
    end

    s_query = collect(range(0.0, total; length = max(8, nsamp)))
    xs = Vector{Float64}(undef, length(s_query))
    ys = similar(xs)
    j = 2
    for (k, s) in enumerate(s_query)
        while j < length(points) && seglen[j] < s
            j += 1
        end
        j0 = max(1, j - 1)
        j1 = min(length(points), j)
        s0 = seglen[j0]
        s1 = seglen[j1]
        t = s1 > s0 ? (s - s0) / (s1 - s0) : 0.0
        x0, y0 = points[j0]
        x1, y1 = points[j1]
        xs[k] = (1.0 - t) * x0 + t * x1
        ys[k] = (1.0 - t) * y0 + t * y1
    end
    return xs, ys, s_query
end

function build_section_surface_polyline(xv, yv, zv, values, path::Vector{Tuple{Float64, Float64}}; nsamp::Int = 360)
    xs, ys, sdist = sample_polyline(path, nsamp)
    ns = length(xs)
    nz = length(zv)
    X = Matrix{Float64}(undef, ns, nz)
    Y = Matrix{Float64}(undef, ns, nz)
    Z = Matrix{Float64}(undef, ns, nz)
    C = Matrix{Float64}(undef, ns, nz)

    for i in 1:ns
        ix0, ix1, wx = bracket_index_and_weight(xv, xs[i])
        iy0, iy1, wy = bracket_index_and_weight(yv, ys[i])
        for k in 1:nz
            X[i, k] = xs[i]
            Y[i, k] = ys[i]
            Z[i, k] = zv[k]
            v00 = values[ix0, iy0, k]
            v10 = values[ix1, iy0, k]
            v01 = values[ix0, iy1, k]
            v11 = values[ix1, iy1, k]
            v0 = (1.0 - wx) * v00 + wx * v10
            v1 = (1.0 - wx) * v01 + wx * v11
            C[i, k] = (1.0 - wy) * v0 + wy * v1
        end
    end
    return X, Y, Z, C, sdist
end

function _sanitize_iso_range(vmin::Real, vmax::Real, cmin::Real, cmax::Real)
    lo = Float64(vmin)
    hi = Float64(vmax)
    if !isfinite(lo) || !isfinite(hi)
        lo, hi = Float64(cmin), Float64(cmax)
    end
    if lo > hi
        lo, hi = hi, lo
    end
    if lo == hi
        delta = max(1e-8, 1e-4 * max(abs(lo), abs(Float64(cmax - cmin)), 1.0))
        lo -= delta
        hi += delta
    end
    return lo, hi
end

function _cube_triangles(x0::Float64, x1::Float64, y0::Float64, y1::Float64, z0::Float64, z1::Float64)
    p000 = (x0, y0, z0)
    p100 = (x1, y0, z0)
    p010 = (x0, y1, z0)
    p110 = (x1, y1, z0)
    p001 = (x0, y0, z1)
    p101 = (x1, y0, z1)
    p011 = (x0, y1, z1)
    p111 = (x1, y1, z1)
    return NTuple{3, NTuple{3, Float64}}[
        (p000, p100, p110), (p000, p110, p010),
        (p001, p111, p101), (p001, p011, p111),
        (p000, p001, p101), (p000, p101, p100),
        (p010, p110, p111), (p010, p111, p011),
        (p000, p010, p011), (p000, p011, p001),
        (p100, p101, p111), (p100, p111, p110),
    ]
end

function build_range_volume_triangles(xv, yv, zv, values, vmin::Real, vmax::Real, dstart_m::Real, dend_m::Real; stride::Int = 1)
    nx, ny, nz = size(values)
    if nx == 0 || ny == 0 || nz == 0
        return NTuple{3, NTuple{3, Float64}}[], 0
    end

    x_edges = edges_from_centers(xv)
    y_edges = edges_from_centers(yv)
    z_edges = edges_from_centers(zv)
    lo, hi = _sanitize_iso_range(vmin, vmax, minimum(values), maximum(values))
    d0 = max(0.0, min(Float64(dstart_m), Float64(dend_m)))
    d1 = max(0.0, max(Float64(dstart_m), Float64(dend_m)))

    step = max(1, stride)
    selected = 0
    triangles = NTuple{3, NTuple{3, Float64}}[]
    for i in 1:step:nx, j in 1:step:ny, k in 1:step:nz
        val = Float64(values[i, j, k])
        depth_m = Float64(zv[k])
        if !(val >= lo && val <= hi && depth_m >= d0 && depth_m <= d1)
            continue
        end
        append!(triangles, _cube_triangles(Float64(x_edges[i]), Float64(x_edges[i + 1]), Float64(y_edges[j]), Float64(y_edges[j + 1]), -Float64(z_edges[k]), -Float64(z_edges[k + 1])))
        selected += 1
    end
    return triangles, selected
end

function triangles_to_vertices_faces(triangles::Vector{NTuple{3, NTuple{3, Float64}}})
    vertices = GLMakie.Point3f[]
    faces = GLMakie.TriangleFace{Int32}[]
    for triangle in triangles
        base = length(vertices) + 1
        push!(vertices, GLMakie.Point3f(triangle[1][1], triangle[1][2], triangle[1][3]))
        push!(vertices, GLMakie.Point3f(triangle[2][1], triangle[2][2], triangle[2][3]))
        push!(vertices, GLMakie.Point3f(triangle[3][1], triangle[3][2], triangle[3][3]))
        push!(faces, GLMakie.TriangleFace(Int32(base), Int32(base + 1), Int32(base + 2)))
    end
    return vertices, faces
end