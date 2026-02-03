# Author: @pankajkmishra
"""
    ForwardFD.jl - Finite Difference Method for Gravity Forward Modeling

Solves the 3D Poisson equation: ∇²Φ = 4πGρ
Uses 7-point Laplacian stencil with Dirichlet boundary conditions.
Gravity computed via: gz = -∂Φ/∂z

Features:
- CG solver with Jacobi preconditioner
- Boundary conditions: monopole approximation or exact prism formula
- 2nd and 4th order derivative schemes for gz computation
"""

"""
    FDOptions

Options for finite difference forward modeling.

# Fields
- `tol`: CG solver tolerance (default: 1e-10)
- `maxiter`: Maximum CG iterations (default: 5000)
- `use_prism_bc`: Use exact prism formula for boundary conditions (default: false)
- `use_4th_order`: Use 4th-order derivative for gz (default: false)
"""
Base.@kwdef struct FDOptions
    tol::Float64 = 1e-10
    maxiter::Int = 5000
    use_prism_bc::Bool = false
    use_4th_order::Bool = false
end

"""
    mass_centroid(model::GravityModel)

Compute total mass and centroid location of density distribution.
Used for monopole boundary condition approximation.

# Returns
- `M`: Total mass [kg]
- `centroid`: (cx, cy, cz) centroid coordinates [m]
"""
function mass_centroid(model::GravityModel)
    mesh = model.mesh
    xc, yc, zc = get_cell_centers(mesh)
    
    M = 0.0
    sx = 0.0
    sy = 0.0
    sz = 0.0
    
    for k in 1:mesh.nz
        for j in 1:mesh.ny
            for i in 1:mesh.nx
                ρ = model.density[i, j, k]
                if ρ != 0.0
                    dv = mesh.dx[i] * mesh.dy[j] * mesh.dz[k]
                    m = ρ * dv
                    M += m
                    sx += m * xc[i]
                    sy += m * yc[j]
                    sz += m * zc[k]
                end
            end
        end
    end
    
    cx = M > 0 ? sx / M : 0.0
    cy = M > 0 ? sy / M : 0.0
    cz = M > 0 ? sz / M : 0.0
    
    return M, (cx, cy, cz)
end

"""
    boundary_potential_monopole(xnodes, ynodes, znodes, M, centroid)

Compute boundary potential using monopole approximation.
Φ = -G*M/r where r is distance from mass centroid.
"""
function boundary_potential_monopole(xnodes::Vector{Float64},
                                     ynodes::Vector{Float64},
                                     znodes::Vector{Float64},
                                     M::Float64,
                                     centroid::NTuple{3,Float64})
    nx1 = length(xnodes)
    ny1 = length(ynodes)
    nz1 = length(znodes)
    
    Φb = zeros(Float64, nx1, ny1, nz1)
    cx, cy, cz = centroid
    
    for k in 1:nz1, j in 1:ny1, i in 1:nx1
        isb = (i == 1 || i == nx1 || j == 1 || j == ny1 || k == 1 || k == nz1)
        if isb
            dx = xnodes[i] - cx
            dy = ynodes[j] - cy
            dz = znodes[k] - cz
            r = sqrt(dx*dx + dy*dy + dz*dz)
            if r > 0.0
                Φb[i, j, k] = -G_NEWTON * M / r
            end
        end
    end
    
    return Φb
end

"""
    cells_to_nodes(density::Array{Float64,3})

Average cell-centered density values to nodes (vertices).
Uses simple 8-cell averaging for each interior node.
"""
function cells_to_nodes(density::Array{Float64,3})
    nx, ny, nz = size(density)
    ρn = zeros(Float64, nx+1, ny+1, nz+1)
    cnt = zeros(Int, nx+1, ny+1, nz+1)
    
    for k in 1:nz, j in 1:ny, i in 1:nx
        v = density[i, j, k]
        for dk in 0:1, dj in 0:1, di in 0:1
            ii = i + di
            jj = j + dj
            kk = k + dk
            ρn[ii, jj, kk] += v
            cnt[ii, jj, kk] += 1
        end
    end
    
    for k in 1:(nz+1), j in 1:(ny+1), i in 1:(nx+1)
        c = cnt[i, j, k]
        if c > 0
            ρn[i, j, k] /= c
        end
    end
    
    return ρn
end

"""
    assemble_rhs(ρnodes, Φb, dx, dy, dz)

Assemble right-hand side vector for Poisson equation.
Includes boundary condition contributions.
"""
function assemble_rhs(ρnodes::Array{Float64,3},
                      Φb::Array{Float64,3},
                      dx::Float64, dy::Float64, dz::Float64)
    nx1, ny1, nz1 = size(ρnodes)
    nxi = nx1 - 2  # Interior nodes
    nyi = ny1 - 2
    nzi = nz1 - 2
    
    ax = 1.0 / (dx * dx)
    ay = 1.0 / (dy * dy)
    az = 1.0 / (dz * dz)
    
    rhs = zeros(Float64, nxi, nyi, nzi)
    c = 4.0 * π * G_NEWTON
    
    for kk in 1:nzi, jj in 1:nyi, ii in 1:nxi
        i = ii + 1
        j = jj + 1
        k = kk + 1
        
        # Source term
        r = -c * ρnodes[i, j, k]
        
        # Boundary contributions
        if ii == 1
            r += ax * Φb[1, j, k]
        end
        if ii == nxi
            r += ax * Φb[nx1, j, k]
        end
        if jj == 1
            r += ay * Φb[i, 1, k]
        end
        if jj == nyi
            r += ay * Φb[i, ny1, k]
        end
        if kk == 1
            r += az * Φb[i, j, 1]
        end
        if kk == nzi
            r += az * Φb[i, j, nz1]
        end
        
        rhs[ii, jj, kk] = r
    end
    
    return rhs
end

"""
    apply_laplacian!(out, u, ax, ay, az)

Apply 7-point Laplacian operator (matrix-free).
"""
function apply_laplacian!(out::Array{Float64,3}, u::Array{Float64,3},
                          ax::Float64, ay::Float64, az::Float64)
    nxi, nyi, nzi = size(u)
    fill!(out, 0.0)
    diag = 2.0 * (ax + ay + az)
    
    @inbounds for k in 1:nzi, j in 1:nyi, i in 1:nxi
        v = diag * u[i, j, k]
        
        if i > 1
            v -= ax * u[i-1, j, k]
        end
        if i < nxi
            v -= ax * u[i+1, j, k]
        end
        if j > 1
            v -= ay * u[i, j-1, k]
        end
        if j < nyi
            v -= ay * u[i, j+1, k]
        end
        if k > 1
            v -= az * u[i, j, k-1]
        end
        if k < nzi
            v -= az * u[i, j, k+1]
        end
        
        out[i, j, k] = v
    end
    
    return out
end

"""
    cg_solve(rhs, dx, dy, dz; tol=1e-10, maxiter=5000)

Solve Poisson equation using Conjugate Gradient with Jacobi preconditioner.

# Returns
- `x`: Solution array
- `iters`: Number of iterations
"""
function cg_solve(rhs::Array{Float64,3},
                  dx::Float64, dy::Float64, dz::Float64;
                  tol::Float64=1e-10, maxiter::Int=5000)
    nxi, nyi, nzi = size(rhs)
    
    ax = 1.0 / (dx * dx)
    ay = 1.0 / (dy * dy)
    az = 1.0 / (dz * dz)
    diag = 2.0 * (ax + ay + az)
    
    # Initialize
    x = zeros(Float64, nxi, nyi, nzi)
    r = copy(rhs)
    z = similar(r)
    @. z = r / diag  # Jacobi preconditioner
    p = copy(z)
    Ap = similar(r)
    
    rs0 = sum(r .* z)
    rs = rs0
    
    if rs0 == 0.0
        return x, 0
    end
    
    for it in 1:maxiter
        apply_laplacian!(Ap, p, ax, ay, az)
        denom = sum(p .* Ap)
        
        if denom == 0.0
            return x, it
        end
        
        α = rs / denom
        @. x = x + α * p
        @. r = r - α * Ap
        @. z = r / diag
        
        rsnew = sum(r .* z)
        
        if sqrt(abs(rsnew)) <= tol * sqrt(abs(rs0))
            return x, it
        end
        
        β = rsnew / rs
        @. p = z + β * p
        rs = rsnew
    end
    
    return x, maxiter
end

"""
    build_full_potential(Φb, uint)

Combine boundary and interior potential values.
"""
function build_full_potential(Φb::Array{Float64,3}, uint::Array{Float64,3})
    Φ = copy(Φb)
    nxi, nyi, nzi = size(uint)
    
    for k in 1:nzi, j in 1:nyi, i in 1:nxi
        Φ[i+1, j+1, k+1] = uint[i, j, k]
    end
    
    return Φ
end

"""
    find_cell(nodes, x)

Find cell index and local coordinate for interpolation.
"""
function find_cell(nodes::Vector{Float64}, x::Float64)
    n = length(nodes)
    
    if x <= nodes[1]
        return 1, 0.0
    end
    if x >= nodes[end]
        return n-1, 1.0
    end
    
    # Binary search
    lo = 1
    hi = n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if nodes[mid] <= x
            lo = mid
        else
            hi = mid
        end
    end
    
    t = (x - nodes[lo]) / (nodes[lo+1] - nodes[lo])
    return lo, t
end

"""
    trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z)

Trilinear interpolation of potential field.
"""
function trilinear_interp(Φ::Array{Float64,3},
                          xnodes::Vector{Float64},
                          ynodes::Vector{Float64},
                          znodes::Vector{Float64},
                          x::Float64, y::Float64, z::Float64)
    ix, tx = find_cell(xnodes, x)
    iy, ty = find_cell(ynodes, y)
    iz, tz = find_cell(znodes, z)
    
    x0, x1 = ix, ix + 1
    y0, y1 = iy, iy + 1
    z0, z1 = iz, iz + 1
    
    c000 = Φ[x0, y0, z0]
    c100 = Φ[x1, y0, z0]
    c010 = Φ[x0, y1, z0]
    c110 = Φ[x1, y1, z0]
    c001 = Φ[x0, y0, z1]
    c101 = Φ[x1, y0, z1]
    c011 = Φ[x0, y1, z1]
    c111 = Φ[x1, y1, z1]
    
    c00 = c000 * (1-tx) + c100 * tx
    c10 = c010 * (1-tx) + c110 * tx
    c01 = c001 * (1-tx) + c101 * tx
    c11 = c011 * (1-tx) + c111 * tx
    
    c0 = c00 * (1-ty) + c10 * ty
    c1 = c01 * (1-ty) + c11 * ty
    
    return c0 * (1-tz) + c1 * tz
end

"""
    gz_from_potential(Φ, xnodes, ynodes, znodes, x, y, z, h)

Compute gz from potential using 2nd-order central difference.
gz = -∂Φ/∂z
"""
function gz_from_potential(Φ::Array{Float64,3},
                           xnodes::Vector{Float64},
                           ynodes::Vector{Float64},
                           znodes::Vector{Float64},
                           x::Float64, y::Float64, z::Float64,
                           h::Float64)
    zmin = znodes[1] + 1e-9
    zmax = znodes[end] - 1e-9
    
    z1 = clamp(z - h, zmin, zmax)
    z2 = clamp(z + h, zmin, zmax)
    
    ϕ1 = trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z1)
    ϕ2 = trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z2)
    
    return -(ϕ2 - ϕ1) / (z2 - z1)
end

"""
    gz_from_potential_4th(Φ, xnodes, ynodes, znodes, x, y, z, h)

Compute gz using 4th-order accurate 5-point stencil.
(-f₋₂ + 8f₋₁ - 8f₁ + f₂) / (12h)
"""
function gz_from_potential_4th(Φ::Array{Float64,3},
                               xnodes::Vector{Float64},
                               ynodes::Vector{Float64},
                               znodes::Vector{Float64},
                               x::Float64, y::Float64, z::Float64,
                               h::Float64)
    zmin = znodes[1] + 1e-9
    zmax = znodes[end] - 1e-9
    
    z_m2 = clamp(z - 2h, zmin, zmax)
    z_m1 = clamp(z - h, zmin, zmax)
    z_p1 = clamp(z + h, zmin, zmax)
    z_p2 = clamp(z + 2h, zmin, zmax)
    
    ϕ_m2 = trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z_m2)
    ϕ_m1 = trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z_m1)
    ϕ_p1 = trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z_p1)
    ϕ_p2 = trilinear_interp(Φ, xnodes, ynodes, znodes, x, y, z_p2)
    
    dϕdz = (-ϕ_p2 + 8ϕ_p1 - 8ϕ_m1 + ϕ_m2) / (12h)
    return -dϕdz
end

"""
    forward_fd(model::GravityModel, receivers::GravityReceivers; 
               options::FDOptions=FDOptions())

Compute gravity forward response using Finite Difference method.

# Arguments
- `model`: GravityModel with density distribution
- `receivers`: Observation locations
- `options`: FDOptions for solver parameters

# Returns
- `GravityData`: Computed gravity data
- `info::Dict`: Solver information (iterations, etc.)
"""
function forward_fd(model::GravityModel, receivers::GravityReceivers;
                    options::FDOptions=FDOptions())
    mesh = model.mesh
    
    # Get uniform cell sizes (required for current FD implementation)
    dx = mesh.dx[1]
    dy = mesh.dy[1]
    dz = mesh.dz[1]
    
    # Build node grid
    xnodes = collect(range(mesh.x0, mesh.x0 + mesh.lx; length=mesh.nx+1))
    ynodes = collect(range(mesh.y0, mesh.y0 + mesh.ly; length=mesh.ny+1))
    znodes = collect(range(mesh.z0, mesh.z0 + mesh.lz; length=mesh.nz+1))
    
    # Interpolate density to nodes
    ρnodes = cells_to_nodes(model.density)
    
    # Compute boundary potential
    if options.use_prism_bc
        Φb = compute_boundary_potential_prism(xnodes, ynodes, znodes, model)
    else
        M, centroid = mass_centroid(model)
        Φb = boundary_potential_monopole(xnodes, ynodes, znodes, M, centroid)
    end
    
    # Assemble and solve
    rhs = assemble_rhs(ρnodes, Φb, dx, dy, dz)
    uint, iters = cg_solve(rhs, dx, dy, dz; tol=options.tol, maxiter=options.maxiter)
    
    # Build full potential field
    Φ = build_full_potential(Φb, uint)
    
    # Compute gz at receivers
    gz = Vector{Float64}(undef, receivers.n)
    step = dz
    
    for ir in 1:receivers.n
        rx = receivers.x[ir]
        ry = receivers.y[ir]
        rz = receivers.z[ir]
        
        if options.use_4th_order && rz > 2.5 * dz
            gz[ir] = gz_from_potential_4th(Φ, xnodes, ynodes, znodes, rx, ry, rz, step)
        else
            gz[ir] = gz_from_potential(Φ, xnodes, ynodes, znodes, rx, ry, rz, step)
        end
    end
    
    info = Dict("iterations" => iters, "method" => "FD")
    
    return GravityData(receivers, gz), info
end

"""
    forward_fd_2d(model::GravityModel, receivers::GravityReceivers;
                  options::FDOptions=FDOptions())

Compute gravity using FD and return as 2D array (for regular receiver grids).

# Returns
- `gz::Array{Float64,2}`: Gravity values reshaped to (nrx, nry)
- `info::Dict`: Solver information
"""
function forward_fd_2d(model::GravityModel, receivers::GravityReceivers;
                       options::FDOptions=FDOptions())
    data, info = forward_fd(model, receivers; options=options)
    return reshape(data.gz, receivers.nx, receivers.ny), info
end
