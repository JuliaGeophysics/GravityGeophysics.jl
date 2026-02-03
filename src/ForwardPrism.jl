# Author: @pankajkmishra
"""
    ForwardPrism.jl - Analytical Prism Method for Gravity Forward Modeling

Computes vertical gravity anomaly (gz) using closed-form rectangular prism formula.
Based on Plouff (1976) and Blakely (1995).

Reference:
- Plouff, D., 1976, Gravity and magnetic fields of polygonal prisms and 
  application to magnetic terrain corrections: Geophysics, 41, 727-741.
- Blakely, R.J., 1995, Potential Theory in Gravity and Magnetic Applications,
  Cambridge University Press.
"""

"""
    gz_prism(rx, ry, rz, x1, x2, y1, y2, z1, z2, drho)

Compute vertical gravity component from a single rectangular prism.

# Arguments
- `rx, ry, rz`: Receiver/observation point coordinates [m]
- `x1, x2`: Prism x bounds [m]
- `y1, y2`: Prism y bounds [m]
- `z1, z2`: Prism z bounds [m] (z positive downward)
- `drho`: Density contrast [kg/m³]

# Returns
- `gz`: Vertical gravity component [m/s²]
"""
function gz_prism(rx::Float64, ry::Float64, rz::Float64,
                  x1::Float64, x2::Float64,
                  y1::Float64, y2::Float64,
                  z1::Float64, z2::Float64,
                  drho::Float64)
    xs = (x1 - rx, x2 - rx)
    ys = (y1 - ry, y2 - ry)
    zs = (z1 - rz, z2 - rz)
    
    s = 0.0
    eps = 1e-30  # Small number to avoid log(0)
    
    @inbounds for ii in 1:2, jj in 1:2, kk in 1:2
        x = xs[ii]
        y = ys[jj]
        z = zs[kk]
        r = sqrt(x*x + y*y + z*z)
        
        # Prism formula terms
        t1 = x * log(max(y + r, eps))
        t2 = y * log(max(x + r, eps))
        t3 = -z * atan(x*y, z*r + eps)
        
        term = t1 + t2 + t3
        sign = ((ii + jj + kk) % 2 == 0) ? 1.0 : -1.0
        s += sign * term
    end
    
    return -G_NEWTON * drho * s
end

"""
    phi_prism(px, py, pz, x1, x2, y1, y2, z1, z2, drho)

Compute gravitational potential from a rectangular prism.

Used for computing exact boundary conditions in FD method.

# Arguments
- `px, py, pz`: Point coordinates [m]
- `x1, x2, y1, y2, z1, z2`: Prism bounds [m]
- `drho`: Density contrast [kg/m³]

# Returns
- `Φ`: Gravitational potential [m²/s²]
"""
function phi_prism(px::Float64, py::Float64, pz::Float64,
                   x1::Float64, x2::Float64,
                   y1::Float64, y2::Float64,
                   z1::Float64, z2::Float64,
                   drho::Float64)
    xs = (x1 - px, x2 - px)
    ys = (y1 - py, y2 - py)
    zs = (z1 - pz, z2 - pz)
    
    s = 0.0
    eps = 1e-30
    
    @inbounds for ii in 1:2, jj in 1:2, kk in 1:2
        x = xs[ii]
        y = ys[jj]
        z = zs[kk]
        r = sqrt(x*x + y*y + z*z)
        
        # Standard potential kernel
        t1 = x * y * log(z + r + eps)
        t2 = y * z * log(x + r + eps)
        t3 = x * z * log(y + r + eps)
        t4 = -0.5 * x * x * atan(y*z, x*r + eps)
        t5 = -0.5 * y * y * atan(x*z, y*r + eps)
        t6 = -0.5 * z * z * atan(x*y, z*r + eps)
        
        term = t1 + t2 + t3 + t4 + t5 + t6
        sign = ((ii + jj + kk) % 2 == 0) ? 1.0 : -1.0
        s += sign * term
    end
    
    return -G_NEWTON * drho * s
end

"""
    forward_prism(model::GravityModel, receivers::GravityReceivers)

Compute gravity forward response using the analytical prism method.

# Arguments
- `model`: GravityModel with density distribution
- `receivers`: Observation locations

# Returns
- `GravityData`: Computed gravity data

# Notes
- Complexity: O(N_cells × N_receivers)
- Exact solution for homogeneous rectangular prisms
- Can be slow for large meshes; consider FD method for large problems
"""
function forward_prism(model::GravityModel, receivers::GravityReceivers)
    mesh = model.mesh
    xe, ye, ze = get_cell_edges(mesh)
    
    # Find non-zero cells for efficiency
    nonzero_idx = findall(!iszero, model.density)
    
    gz = Vector{Float64}(undef, receivers.n)
    
    # Parallel-friendly loop (can add @threads later)
    for ir in 1:receivers.n
        rx = receivers.x[ir]
        ry = receivers.y[ir]
        rz = receivers.z[ir]
        
        acc = 0.0
        @inbounds for I in nonzero_idx
            i, j, k = I.I
            drho = model.density[i, j, k]
            acc += gz_prism(rx, ry, rz,
                           xe[i], xe[i+1],
                           ye[j], ye[j+1],
                           ze[k], ze[k+1],
                           drho)
        end
        gz[ir] = acc
    end
    
    return GravityData(receivers, gz)
end

"""
    forward_prism_2d(model::GravityModel, receivers::GravityReceivers)

Compute gravity and return as 2D array (for regular receiver grids).

# Returns
- `gz::Array{Float64,2}`: Gravity values reshaped to (nrx, nry)
"""
function forward_prism_2d(model::GravityModel, receivers::GravityReceivers)
    data = forward_prism(model, receivers)
    return reshape(data.gz, receivers.nx, receivers.ny)
end

"""
    compute_boundary_potential_prism(xnodes, ynodes, znodes, model, mesh)

Compute exact boundary potential using prism formula for FD boundary conditions.
"""
function compute_boundary_potential_prism(xnodes::Vector{Float64},
                                          ynodes::Vector{Float64},
                                          znodes::Vector{Float64},
                                          model::GravityModel)
    mesh = model.mesh
    nx1 = length(xnodes)
    ny1 = length(ynodes)
    nz1 = length(znodes)
    
    Φb = zeros(Float64, nx1, ny1, nz1)
    
    # Get cell edges
    xe, ye, ze = get_cell_edges(mesh)
    
    # Find non-zero cells
    nonzero_idx = findall(!iszero, model.density)
    
    for k in 1:nz1, j in 1:ny1, i in 1:nx1
        # Only boundary nodes
        isb = (i == 1 || i == nx1 || j == 1 || j == ny1 || k == 1 || k == nz1)
        if isb
            px, py, pz = xnodes[i], ynodes[j], znodes[k]
            acc = 0.0
            @inbounds for I in nonzero_idx
                ci, cj, ck = I.I
                drho = model.density[ci, cj, ck]
                acc += phi_prism(px, py, pz,
                                xe[ci], xe[ci+1],
                                ye[cj], ye[cj+1],
                                ze[ck], ze[ck+1],
                                drho)
            end
            Φb[i, j, k] = acc
        end
    end
    
    return Φb
end
