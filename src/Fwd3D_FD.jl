# 3D gravity forward solver, finite-difference alternative to Fwd3D.jl
# Author: @pankajkmishra
# Solves the Poisson equation ∇²Φ = 4πγρ (Blakely 1995 notation, γ the gravitational constant, ρ the density
# contrast) by a node-centred finite-volume discretisation of the
# model's tensor grid, extended outward by geometrically growing air padding so that every sensor and its
# differencing stencil lie strictly inside the domain: gravity stations sit in the air above the model, which
# the model's own grid does not cover. The finite-volume operator is symmetric positive definite on any
# tensor grid, so the padded and stretched meshes MakeMesh3D writes are handled as they are, and air cells
# (NaN under topography) enter as zero density. Dirichlet boundary values come from a monopole approximation
# of the model's mass, which the padding pushes far from it. gz = -∂Φ/∂z with z depth positive down, as in
# Fwd3D.jl; GAMMA and SI_TO_MGAL come from Fwd3D.jl.

"""
    FDCtrl3D(; tol=1e-10, maxiter=5000, n_pad=8, pad_factor=1.4, use_4th_order=false)

Input of `forward3d_fd()`: how hard to solve the Poisson equation and how far to carry the boundary from the mass.
"""
Base.@kwdef struct FDCtrl3D
    tol::Float64        = 1e-10  # conjugate-gradient relative residual to stop at
    maxiter::Int        = 5000   # maximum conjugate-gradient iterations
    n_pad::Int          = 8      # air cells added beyond each face of the model
    pad_factor::Float64 = 1.4    # how fast those padding cells grow outwards
    use_4th_order::Bool = false  # 5-point instead of 3-point central difference for gz
end

#---------- extended grid ----------

# padding widths beyond one face: at least `n_pad` cells growing by `factor`, continued to span `reach`
function _pad_side(w0::Real, n_pad::Integer, factor::Real, reach::Real)
    out   = Float64[]
    h     = Float64(w0)
    total = 0.0
    while length(out) < n_pad || total < reach
        push!(out, h)
        total += h
        h *= factor
        length(out) >= 256 && break
    end
    out
end

# the model's tensor grid padded on all six faces, with the offsets of the model block inside it
function _extended_grid(m::GravModel3D, depth::AbstractVector{<:Real}, fd::FDCtrl3D)
    zlo = max(0.0, m.z[1]   - minimum(depth)) + 3 * m.dz[1]
    zhi = max(0.0, maximum(depth) - m.z[end]) + 3 * m.dz[end]
    pxl = _pad_side(m.dx[1],   fd.n_pad, fd.pad_factor, 0.0)
    pxr = _pad_side(m.dx[end], fd.n_pad, fd.pad_factor, 0.0)
    pyl = _pad_side(m.dy[1],   fd.n_pad, fd.pad_factor, 0.0)
    pyr = _pad_side(m.dy[end], fd.n_pad, fd.pad_factor, 0.0)
    pzl = _pad_side(m.dz[1],   fd.n_pad, fd.pad_factor, zlo)
    pzr = _pad_side(m.dz[end], fd.n_pad, fd.pad_factor, zhi)
    dx  = vcat(reverse(pxl), m.dx, pxr)
    dy  = vcat(reverse(pyl), m.dy, pyr)
    dz  = vcat(reverse(pzl), m.dz, pzr)
    xn  = cumsum(vcat(m.x[1] - sum(pxl), dx))
    yn  = cumsum(vcat(m.y[1] - sum(pyl), dy))
    zn  = cumsum(vcat(m.z[1] - sum(pzl), dz))
    (; dx, dy, dz, xn, yn, zn, off = (length(pxl), length(pyl), length(pzl)))
end

# density contrast on the extended cells, model air (NaN) and padding as zero
function _extended_density(m::GravModel3D, g)
    A = zeros(Float64, length(g.dx), length(g.dy), length(g.dz))
    ox, oy, oz = g.off
    @inbounds for k in axes(m.A, 3), j in axes(m.A, 2), i in axes(m.A, 1)
        v = m.A[i, j, k]
        isnan(v) || (A[i+ox, j+oy, k+oz] = v)
    end
    A
end

# width of a node's dual cell along one axis, half a cell at each end
_dual_widths(d::Vector{Float64}) =
    [((i > 1 ? d[i-1] : 0.0) + (i <= length(d) ? d[i] : 0.0)) / 2 for i in 1:(length(d) + 1)]

#---------- boundary condition ----------

# total mass and centroid of the density contrast, in the model's (north, east, depth) frame
function _mass_centroid(A::Array{Float64, 3}, xc, yc, zc, V)
    M  = 0.0
    sx = sy = sz = 0.0
    @inbounds for k in eachindex(zc), j in eachindex(yc), i in eachindex(xc)
        ρ = A[i, j, k]
        ρ == 0.0 && continue
        v   = ρ * V[i, j, k]
        M  += v
        sx += v * xc[i]
        sy += v * yc[j]
        sz += v * zc[k]
    end
    M == 0.0 ? (0.0, (0.0, 0.0, 0.0)) : (M, (sx / M, sy / M, sz / M))
end

# Dirichlet potential Φ = -γM/r on the faces of the node grid, interior left zero
function _boundary_potential(xn, yn, zn, M::Real, centroid)
    nx, ny, nz = length(xn), length(yn), length(zn)
    Φb = zeros(Float64, nx, ny, nz)
    cx, cy, cz = centroid
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        (i == 1 || i == nx || j == 1 || j == ny || k == 1 || k == nz) || continue
        r = sqrt((xn[i] - cx)^2 + (yn[j] - cy)^2 + (zn[k] - cz)^2)
        r > 0 && (Φb[i, j, k] = -GAMMA * M / r)
    end
    Φb
end

#---------- finite-volume Poisson system ----------

# conductances of the faces meeting each interior node, and the diagonal they sum to
function _coefficients(g)
    hx, hy, hz = _dual_widths(g.dx), _dual_widths(g.dy), _dual_widths(g.dz)
    nxi, nyi, nzi = length(g.xn) - 2, length(g.yn) - 2, length(g.zn) - 2
    cxm = [1 / g.dx[ii]     for ii in 1:nxi]
    cxp = [1 / g.dx[ii+1]   for ii in 1:nxi]
    cym = [1 / g.dy[jj]     for jj in 1:nyi]
    cyp = [1 / g.dy[jj+1]   for jj in 1:nyi]
    czm = [1 / g.dz[kk]     for kk in 1:nzi]
    czp = [1 / g.dz[kk+1]   for kk in 1:nzi]
    (; hx = hx[2:end-1], hy = hy[2:end-1], hz = hz[2:end-1], cxm, cxp, cym, cyp, czm, czp, nxi, nyi, nzi)
end

function _diagonal(c)
    D = Array{Float64, 3}(undef, c.nxi, c.nyi, c.nzi)
    @inbounds for kk in 1:c.nzi, jj in 1:c.nyi, ii in 1:c.nxi
        D[ii, jj, kk] = (c.cxm[ii] + c.cxp[ii]) * c.hy[jj] * c.hz[kk] +
                        (c.cym[jj] + c.cyp[jj]) * c.hx[ii] * c.hz[kk] +
                        (c.czm[kk] + c.czp[kk]) * c.hx[ii] * c.hy[jj]
    end
    D
end

# right-hand side -4πγ times the mass in each node's dual cell, with the Dirichlet faces folded in
function _assemble_rhs(A, g, c, Φb)
    rhs = zeros(Float64, c.nxi, c.nyi, c.nzi)
    @inbounds for kk in 1:c.nzi, jj in 1:c.nyi, ii in 1:c.nxi
        i, j, k = ii + 1, jj + 1, kk + 1
        mass = 0.0
        for dk in 0:1, dj in 0:1, di in 0:1
            ci, cj, ck = i - 1 + di, j - 1 + dj, k - 1 + dk
            mass += A[ci, cj, ck] * g.dx[ci] * g.dy[cj] * g.dz[ck] / 8
        end
        r = -4π * GAMMA * mass
        ii == 1      && (r += c.cxm[ii] * c.hy[jj] * c.hz[kk] * Φb[1, j, k])
        ii == c.nxi  && (r += c.cxp[ii] * c.hy[jj] * c.hz[kk] * Φb[end, j, k])
        jj == 1      && (r += c.cym[jj] * c.hx[ii] * c.hz[kk] * Φb[i, 1, k])
        jj == c.nyi  && (r += c.cyp[jj] * c.hx[ii] * c.hz[kk] * Φb[i, end, k])
        kk == 1      && (r += c.czm[kk] * c.hx[ii] * c.hy[jj] * Φb[i, j, 1])
        kk == c.nzi  && (r += c.czp[kk] * c.hx[ii] * c.hy[jj] * Φb[i, j, end])
        rhs[ii, jj, kk] = r
    end
    rhs
end

# matrix-free symmetric finite-volume Laplacian, Dirichlet (zero) outside the interior block
function _laplacian!(out::Array{Float64, 3}, u::Array{Float64, 3}, c, D)
    @inbounds for kk in 1:c.nzi, jj in 1:c.nyi, ii in 1:c.nxi
        v = D[ii, jj, kk] * u[ii, jj, kk]
        ii > 1     && (v -= c.cxm[ii] * c.hy[jj] * c.hz[kk] * u[ii-1, jj, kk])
        ii < c.nxi && (v -= c.cxp[ii] * c.hy[jj] * c.hz[kk] * u[ii+1, jj, kk])
        jj > 1     && (v -= c.cym[jj] * c.hx[ii] * c.hz[kk] * u[ii, jj-1, kk])
        jj < c.nyi && (v -= c.cyp[jj] * c.hx[ii] * c.hz[kk] * u[ii, jj+1, kk])
        kk > 1     && (v -= c.czm[kk] * c.hx[ii] * c.hy[jj] * u[ii, jj, kk-1])
        kk < c.nzi && (v -= c.czp[kk] * c.hx[ii] * c.hy[jj] * u[ii, jj, kk+1])
        out[ii, jj, kk] = v
    end
    out
end

# conjugate gradients with Jacobi preconditioning, returning the interior solution and iteration count
function _cg_poisson(rhs::Array{Float64, 3}, c, D; tol::Real = 1e-10, maxiter::Integer = 5000)
    x   = zeros(size(rhs))
    r   = copy(rhs)
    z   = r ./ D
    p   = copy(z)
    Ap  = similar(r)
    rs0 = sum(r .* z)
    rs0 == 0.0 && return x, 0
    rs   = rs0
    nrm0 = sqrt(sum(abs2, r))
    for it in 1:maxiter
        _laplacian!(Ap, p, c, D)
        denom = sum(p .* Ap)
        denom == 0.0 && return x, it
        α = rs / denom
        @. x += α * p
        @. r -= α * Ap
        sqrt(sum(abs2, r)) <= tol * nrm0 && return x, it
        @. z = r / D
        rsnew = sum(r .* z)
        rs == 0.0 && return x, it
        β  = rsnew / rs
        @. p = z + β * p
        rs = rsnew
    end
    x, maxiter
end

#---------- interpolation and the vertical derivative ----------

function _find_cell(nodes::Vector{Float64}, v::Real)
    n = length(nodes)
    v <= nodes[1]   && return 1, 0.0
    v >= nodes[end] && return n - 1, 1.0
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        nodes[mid] <= v ? (lo = mid) : (hi = mid)
    end
    lo, (v - nodes[lo]) / (nodes[lo+1] - nodes[lo])
end

function _trilinear(Φ::Array{Float64, 3}, xn, yn, zn, x::Real, y::Real, z::Real)
    ix, tx = _find_cell(xn, x)
    iy, ty = _find_cell(yn, y)
    iz, tz = _find_cell(zn, z)
    c00 = Φ[ix, iy,   iz  ] * (1 - tx) + Φ[ix+1, iy,   iz  ] * tx
    c10 = Φ[ix, iy+1, iz  ] * (1 - tx) + Φ[ix+1, iy+1, iz  ] * tx
    c01 = Φ[ix, iy,   iz+1] * (1 - tx) + Φ[ix+1, iy,   iz+1] * tx
    c11 = Φ[ix, iy+1, iz+1] * (1 - tx) + Φ[ix+1, iy+1, iz+1] * tx
    c0  = c00 * (1 - ty) + c10 * ty
    c1  = c01 * (1 - ty) + c11 * ty
    c0 * (1 - tz) + c1 * tz
end

# gz = -dPhi/dz by central differences over the local cell height, 4th order where the wider stencil fits
function _gz(Φ, g, x::Real, y::Real, z::Real, fourth::Bool)
    iz, _ = _find_cell(g.zn, z)
    h     = g.dz[iz]
    lo, hi = g.zn[1], g.zn[end]
    Φat(zi) = _trilinear(Φ, g.xn, g.yn, g.zn, x, y, zi)
    if fourth && z - 2h > lo && z + 2h < hi
        return -(-Φat(z + 2h) + 8Φat(z + h) - 8Φat(z - h) + Φat(z - 2h)) / (12h)
    end
    z1, z2 = max(z - h, lo), min(z + h, hi)
    z2 > z1 || error("the finite-difference grid is too thin to difference the potential at depth $z")
    -(Φat(z2) - Φat(z1)) / (z2 - z1)
end

#---------- forward ----------

"""
    forward3d_fd(model, data, fd=FDCtrl3D()) -> Vector{Float64}

Predicted gz (mGal) of a model by finite differences: an independent check of `forward3d` (Fwd3D.jl), or a base for meshes too large for a sensitivity matrix.
"""
function forward3d_fd(m::GravModel3D, d::GravData3D, fd::FDCtrl3D = FDCtrl3D())
    check_receivers(m, d)
    depth = sensor_depth(m, d)
    g     = _extended_grid(m, depth, fd)
    A     = _extended_density(m, g)
    xc    = g.xn[1:end-1] .+ g.dx ./ 2
    yc    = g.yn[1:end-1] .+ g.dy ./ 2
    zc    = g.zn[1:end-1] .+ g.dz ./ 2
    V     = [a * b * c for a in g.dx, b in g.dy, c in g.dz]
    M, centroid = _mass_centroid(A, xc, yc, zc, V)
    Φb    = _boundary_potential(g.xn, g.yn, g.zn, M, centroid)
    c     = _coefficients(g)
    D     = _diagonal(c)
    rhs   = _assemble_rhs(A, g, c, Φb)
    uint, iters = _cg_poisson(rhs, c, D; tol = fd.tol, maxiter = fd.maxiter)
    Φ = copy(Φb)
    Φ[2:end-1, 2:end-1, 2:end-1] .= uint
    pred = Vector{Float64}(undef, ndata(d))
    for i in eachindex(pred)
        pred[i] = _gz(Φ, g, d.x[i], d.y[i], depth[i], fd.use_4th_order) * SI_TO_MGAL
    end
    iters >= fd.maxiter &&
        @warn "forward3d_fd: conjugate gradients did not reach tol = $(fd.tol) in $(fd.maxiter) iterations"
    pred
end
