# Author: @pankajkmishra
"""
    Mesh.jl - 3D Mesh/Grid Creation for Gravity Forward Modeling

Provides mesh creation, cell center/edge computation, and mesh utilities.
Supports both uniform and variable cell sizes in UBC format.
"""

"""
    GravityMesh

3D mesh structure for gravity forward modeling.

# Fields
- `nx, ny, nz::Int`: Number of cells in each direction
- `dx, dy, dz::Vector{Float64}`: Cell widths (can be uniform or variable)
- `x0, y0, z0::Float64`: Origin coordinates (top-SW corner in UBC convention)
- `lx, ly, lz::Float64`: Total domain lengths
"""
struct GravityMesh
    # Number of cells
    nx::Int
    ny::Int
    nz::Int
    
    # Cell widths (for variable mesh)
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    
    # Origin (top-SW corner in UBC convention, z positive down)
    x0::Float64
    y0::Float64
    z0::Float64
    
    # Total domain size
    lx::Float64
    ly::Float64
    lz::Float64
end

"""
    create_mesh(; nx=64, ny=64, nz=32, lx=10000.0, ly=10000.0, lz=5000.0,
                x0=0.0, y0=0.0, z0=0.0, dx=nothing, dy=nothing, dz=nothing)

Create a 3D mesh for gravity forward modeling.

# Keyword Arguments
- `nx, ny, nz`: Number of cells in x, y, z directions
- `lx, ly, lz`: Total domain lengths [m]
- `x0, y0, z0`: Origin coordinates [m] (default: centered for x,y; 0 for z)
- `dx, dy, dz`: Optional cell width vectors for variable mesh

# Returns
- `GravityMesh`: Mesh structure

# Example
```julia
# Uniform mesh centered at origin
mesh = create_mesh(nx=64, ny=64, nz=32, lx=160000.0, ly=160000.0, lz=80000.0)

# Variable mesh (provide cell widths)
dx = [100.0, 200.0, 500.0, ...]  # cell widths
mesh = create_mesh(nx=length(dx), dx=dx, ...)
```
"""
function create_mesh(; nx::Int=64, ny::Int=64, nz::Int=32,
                     lx::Float64=10000.0, ly::Float64=10000.0, lz::Float64=5000.0,
                     x0::Union{Float64,Nothing}=nothing, 
                     y0::Union{Float64,Nothing}=nothing, 
                     z0::Float64=0.0,
                     dx::Union{Vector{Float64},Nothing}=nothing,
                     dy::Union{Vector{Float64},Nothing}=nothing,
                     dz::Union{Vector{Float64},Nothing}=nothing)
    
    # Create cell width vectors
    if isnothing(dx)
        dx_vec = fill(lx / nx, nx)
    else
        @assert length(dx) == nx "dx length must match nx"
        dx_vec = dx
        lx = sum(dx)
    end
    
    if isnothing(dy)
        dy_vec = fill(ly / ny, ny)
    else
        @assert length(dy) == ny "dy length must match ny"
        dy_vec = dy
        ly = sum(dy)
    end
    
    if isnothing(dz)
        dz_vec = fill(lz / nz, nz)
    else
        @assert length(dz) == nz "dz length must match nz"
        dz_vec = dz
        lz = sum(dz)
    end
    
    # Default origin: centered in x,y, surface at z=0
    x0_final = isnothing(x0) ? -0.5 * lx : x0
    y0_final = isnothing(y0) ? -0.5 * ly : y0
    
    return GravityMesh(nx, ny, nz, dx_vec, dy_vec, dz_vec, 
                       x0_final, y0_final, z0, lx, ly, lz)
end

"""
    get_cell_centers(mesh::GravityMesh)

Get cell center coordinates.

# Returns
- `xc, yc, zc::Vector{Float64}`: Cell center coordinates
"""
function get_cell_centers(mesh::GravityMesh)
    # X centers
    xc = Vector{Float64}(undef, mesh.nx)
    x = mesh.x0
    for i in 1:mesh.nx
        xc[i] = x + 0.5 * mesh.dx[i]
        x += mesh.dx[i]
    end
    
    # Y centers
    yc = Vector{Float64}(undef, mesh.ny)
    y = mesh.y0
    for j in 1:mesh.ny
        yc[j] = y + 0.5 * mesh.dy[j]
        y += mesh.dy[j]
    end
    
    # Z centers (positive down)
    zc = Vector{Float64}(undef, mesh.nz)
    z = mesh.z0
    for k in 1:mesh.nz
        zc[k] = z + 0.5 * mesh.dz[k]
        z += mesh.dz[k]
    end
    
    return xc, yc, zc
end

"""
    get_cell_edges(mesh::GravityMesh)

Get cell edge/face coordinates.

# Returns
- `xe, ye, ze::Vector{Float64}`: Cell edge coordinates (n+1 points each)
"""
function get_cell_edges(mesh::GravityMesh)
    # X edges
    xe = Vector{Float64}(undef, mesh.nx + 1)
    xe[1] = mesh.x0
    for i in 1:mesh.nx
        xe[i+1] = xe[i] + mesh.dx[i]
    end
    
    # Y edges
    ye = Vector{Float64}(undef, mesh.ny + 1)
    ye[1] = mesh.y0
    for j in 1:mesh.ny
        ye[j+1] = ye[j] + mesh.dy[j]
    end
    
    # Z edges
    ze = Vector{Float64}(undef, mesh.nz + 1)
    ze[1] = mesh.z0
    for k in 1:mesh.nz
        ze[k+1] = ze[k] + mesh.dz[k]
    end
    
    return xe, ye, ze
end

"""
    get_nodes(mesh::GravityMesh)

Alias for get_cell_edges. Get node coordinates (vertices of cells).
"""
get_nodes(mesh::GravityMesh) = get_cell_edges(mesh)

mesh_cell_count(mesh::GravityMesh) = mesh.nx * mesh.ny * mesh.nz

function cell_volumes(mesh::GravityMesh)
    volumes = Array{Float64}(undef, mesh.nx, mesh.ny, mesh.nz)
    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        volumes[i, j, k] = mesh.dx[i] * mesh.dy[j] * mesh.dz[k]
    end
    return volumes
end

"""
    mesh_info(mesh::GravityMesh)

Print mesh information.
"""
function mesh_info(mesh::GravityMesh)
    println("═"^60)
    println("GravityMesh Information")
    println("═"^60)
    @printf("  Cells:       %d × %d × %d = %d total\n", 
            mesh.nx, mesh.ny, mesh.nz, mesh.nx * mesh.ny * mesh.nz)
    @printf("  Domain:      %.2f × %.2f × %.2f m\n", mesh.lx, mesh.ly, mesh.lz)
    @printf("  Domain (km): %.2f × %.2f × %.2f km\n", 
            mesh.lx/1000, mesh.ly/1000, mesh.lz/1000)
    @printf("  Origin:      (%.2f, %.2f, %.2f) m\n", mesh.x0, mesh.y0, mesh.z0)
    
    # Check if uniform
    is_uniform_x = all(d ≈ mesh.dx[1] for d in mesh.dx)
    is_uniform_y = all(d ≈ mesh.dy[1] for d in mesh.dy)
    is_uniform_z = all(d ≈ mesh.dz[1] for d in mesh.dz)
    
    if is_uniform_x && is_uniform_y && is_uniform_z
        @printf("  Cell size:   %.2f × %.2f × %.2f m (uniform)\n",
                mesh.dx[1], mesh.dy[1], mesh.dz[1])
    else
        @printf("  Cell dx:     [%.2f ... %.2f] m\n", minimum(mesh.dx), maximum(mesh.dx))
        @printf("  Cell dy:     [%.2f ... %.2f] m\n", minimum(mesh.dy), maximum(mesh.dy))
        @printf("  Cell dz:     [%.2f ... %.2f] m\n", minimum(mesh.dz), maximum(mesh.dz))
    end
    println("═"^60)
end

"""
    GravityModel

3D density model on a mesh.

# Fields
- `mesh::GravityMesh`: Associated mesh
- `density::Array{Float64,3}`: Density values [kg/m³]
"""
struct GravityModel
    mesh::GravityMesh
    density::Array{Float64,3}  # (nx, ny, nz)
end

"""
    create_model(mesh::GravityMesh; background=0.0)

Create a gravity model with optional background density.

# Arguments
- `mesh`: GravityMesh structure
- `background`: Background density [kg/m³]
"""
function create_model(mesh::GravityMesh; background::Float64=0.0)
    density = fill(background, mesh.nx, mesh.ny, mesh.nz)
    return GravityModel(mesh, density)
end

"""
    set_density!(model::GravityModel, density::Array{Float64,3})

Set density values for the entire model.
"""
function set_density!(model::GravityModel, density::Array{Float64,3})
    @assert size(density) == (model.mesh.nx, model.mesh.ny, model.mesh.nz)
    model.density .= density
    return model
end

"""
    add_block_anomaly!(model::GravityModel; 
                       xmin, xmax, ymin, ymax, zmin, zmax, drho)

Add a rectangular block density anomaly to the model.

# Arguments
- `xmin, xmax`: X bounds [m]
- `ymin, ymax`: Y bounds [m]
- `zmin, zmax`: Z bounds [m] (positive down)
- `drho`: Density contrast [kg/m³]
"""
function add_block_anomaly!(model::GravityModel;
                            xmin::Float64, xmax::Float64,
                            ymin::Float64, ymax::Float64,
                            zmin::Float64, zmax::Float64,
                            drho::Float64)
    xc, yc, zc = get_cell_centers(model.mesh)
    
    for k in 1:model.mesh.nz
        for j in 1:model.mesh.ny
            for i in 1:model.mesh.nx
                if xmin <= xc[i] <= xmax && 
                   ymin <= yc[j] <= ymax && 
                   zmin <= zc[k] <= zmax
                    model.density[i,j,k] += drho
                end
            end
        end
    end
    return model
end

"""
    add_sphere_anomaly!(model::GravityModel;
                        cx, cy, cz, radius, drho)

Add a spherical density anomaly to the model.

# Arguments
- `cx, cy, cz`: Center coordinates [m]
- `radius`: Sphere radius [m]
- `drho`: Density contrast [kg/m³]
"""
function add_sphere_anomaly!(model::GravityModel;
                             cx::Float64, cy::Float64, cz::Float64,
                             radius::Float64, drho::Float64)
    xc, yc, zc = get_cell_centers(model.mesh)
    r2 = radius^2
    
    for k in 1:model.mesh.nz
        for j in 1:model.mesh.ny
            for i in 1:model.mesh.nx
                dist2 = (xc[i] - cx)^2 + (yc[j] - cy)^2 + (zc[k] - cz)^2
                if dist2 <= r2
                    model.density[i,j,k] += drho
                end
            end
        end
    end
    return model
end

"""
    GravityReceivers

Receiver/observation locations for gravity data.

# Fields
- `x, y, z::Vector{Float64}`: Receiver coordinates
- `n::Int`: Number of receivers
"""
struct GravityReceivers
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    n::Int
    # Optional: for regular grid
    nx::Int
    ny::Int
end

"""
    create_receivers(mesh::GravityMesh; nrx=51, nry=51, z=50.0, 
                     xfrac=0.9, yfrac=0.9)

Create a regular grid of receivers.

# Arguments
- `mesh`: Reference mesh for domain extent
- `nrx, nry`: Number of receivers in x, y
- `z`: Observation height (negative = above surface) [m]
- `xfrac, yfrac`: Fraction of domain to cover (default 90%)
"""
function create_receivers(mesh::GravityMesh; 
                          nrx::Int=51, nry::Int=51, z::Float64=50.0,
                          xfrac::Float64=0.9, yfrac::Float64=0.9)
    # Receiver extent
    rx_min = mesh.x0 + (1 - xfrac) / 2 * mesh.lx
    rx_max = mesh.x0 + (1 + xfrac) / 2 * mesh.lx
    ry_min = mesh.y0 + (1 - yfrac) / 2 * mesh.ly
    ry_max = mesh.y0 + (1 + yfrac) / 2 * mesh.ly
    
    rx = range(rx_min, rx_max; length=nrx) |> collect
    ry = range(ry_min, ry_max; length=nry) |> collect
    
    # Create full coordinate vectors
    n = nrx * nry
    x = Vector{Float64}(undef, n)
    y = Vector{Float64}(undef, n)
    zv = fill(z, n)
    
    idx = 1
    for iy in 1:nry
        for ix in 1:nrx
            x[idx] = rx[ix]
            y[idx] = ry[iy]
            idx += 1
        end
    end
    
    return GravityReceivers(x, y, zv, n, nrx, nry)
end

"""
    create_receivers(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})

Create receivers from arbitrary coordinate vectors.
"""
function create_receivers(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    @assert length(x) == length(y) == length(z)
    return GravityReceivers(x, y, z, length(x), length(x), 1)
end

"""
    GravityData

Gravity data container.

# Fields
- `receivers::GravityReceivers`: Observation locations
- `gz::Vector{Float64}`: Vertical gravity component [m/s²]
- `error::Vector{Float64}`: Data uncertainty [m/s²]
"""
struct GravityData
    receivers::GravityReceivers
    gz::Vector{Float64}
    error::Vector{Float64}
end

"""
    GravityData(receivers, gz; error=nothing)

Create gravity data with optional error estimates.
"""
function GravityData(receivers::GravityReceivers, gz::Vector{Float64}; 
                     error::Union{Vector{Float64},Nothing}=nothing)
    if isnothing(error)
        error = zeros(length(gz))
    end
    return GravityData(receivers, gz, error)
end
