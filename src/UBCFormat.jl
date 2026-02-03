# Author: @pankajkmishra
"""
    UBCFormat.jl - UBC-GIF Format I/O for Gravity Models and Data

Supports reading and writing:
- Mesh files (.msh)
- Model files (density)
- Data files (gravity observations)
- Locations files

Reference: UBC-GIF GRAV3D documentation
"""

"""
    load_mesh_ubc(filename::String)

Load a mesh from UBC format mesh file.

# UBC Mesh Format:
```
nE nN nZ
E0 N0 Z0
dE1 dE2 ... dEnE
dN1 dN2 ... dNnN
dZ1 dZ2 ... dZnZ
```

Where:
- nE, nN, nZ: Number of cells in East, North, Vertical
- E0, N0, Z0: Origin (SW-Top corner)
- dE, dN, dZ: Cell widths
"""
function load_mesh_ubc(filename::String)
    lines = readlines(filename)
    idx = 1
    
    # Skip comments
    while startswith(lines[idx], "#") || startswith(lines[idx], "!")
        idx += 1
    end
    
    # Read mesh dimensions
    dims = parse.(Int, split(strip(lines[idx])))
    nx, ny, nz = dims[1], dims[2], dims[3]
    idx += 1
    
    # Read origin
    origin = parse.(Float64, split(strip(lines[idx])))
    x0, y0, z0 = origin[1], origin[2], origin[3]
    idx += 1
    
    # Read cell widths
    # Can be on single line or multiple lines
    dx = _read_cell_widths(lines, idx, nx)
    idx += _count_value_lines(lines, idx, nx)
    
    dy = _read_cell_widths(lines, idx, ny)
    idx += _count_value_lines(lines, idx, ny)
    
    dz = _read_cell_widths(lines, idx, nz)
    
    return GravityMesh(nx, ny, nz, dx, dy, dz, x0, y0, z0, sum(dx), sum(dy), sum(dz))
end

"""Helper to read cell widths (may span multiple lines)"""
function _read_cell_widths(lines::Vector{String}, start_idx::Int, n::Int)
    values = Float64[]
    idx = start_idx
    while length(values) < n && idx <= length(lines)
        line = strip(lines[idx])
        if !isempty(line) && !startswith(line, "#")
            # Check for repeated value format: n*value
            for token in split(line)
                if occursin("*", token)
                    parts = split(token, "*")
                    count = parse(Int, parts[1])
                    value = parse(Float64, parts[2])
                    append!(values, fill(value, count))
                else
                    push!(values, parse(Float64, token))
                end
            end
        end
        idx += 1
    end
    return values[1:n]
end

"""Count how many lines are needed to read n values"""
function _count_value_lines(lines::Vector{String}, start_idx::Int, n::Int)
    count = 0
    values = 0
    idx = start_idx
    while values < n && idx <= length(lines)
        line = strip(lines[idx])
        if !isempty(line) && !startswith(line, "#")
            for token in split(line)
                if occursin("*", token)
                    parts = split(token, "*")
                    values += parse(Int, parts[1])
                else
                    values += 1
                end
            end
        end
        count += 1
        idx += 1
    end
    return count
end

"""
    save_mesh_ubc(filename::String, mesh::GravityMesh)

Save mesh to UBC format.
"""
function save_mesh_ubc(filename::String, mesh::GravityMesh)
    open(filename, "w") do io
        # Write dimensions
        @printf(io, "%d %d %d\n", mesh.nx, mesh.ny, mesh.nz)
        
        # Write origin
        @printf(io, "%.6f %.6f %.6f\n", mesh.x0, mesh.y0, mesh.z0)
        
        # Write cell widths (check if uniform for compact format)
        _write_cell_widths(io, mesh.dx)
        _write_cell_widths(io, mesh.dy)
        _write_cell_widths(io, mesh.dz)
    end
    return nothing
end

"""Helper to write cell widths with optional compression"""
function _write_cell_widths(io::IO, widths::Vector{Float64}; max_per_line::Int=10)
    # Check if uniform
    if all(w ≈ widths[1] for w in widths)
        @printf(io, "%d*%.6f\n", length(widths), widths[1])
    else
        # Write individual values
        for i in 1:length(widths)
            @printf(io, "%.6f", widths[i])
            if i < length(widths)
                print(io, (i % max_per_line == 0) ? "\n" : " ")
            end
        end
        println(io)
    end
end

"""
    load_model_ubc(filename::String, mesh::GravityMesh)

Load density model from UBC format.

# UBC Model Format:
Values listed in order (fastest to slowest): E, N, Z
One value per line or space-separated.
"""
function load_model_ubc(filename::String, mesh::GravityMesh)
    lines = readlines(filename)
    idx = 1
    
    # Skip header comments
    while idx <= length(lines) && (startswith(lines[idx], "#") || startswith(lines[idx], "!"))
        idx += 1
    end
    
    # Read all values
    values = Float64[]
    ntotal = mesh.nx * mesh.ny * mesh.nz
    
    while length(values) < ntotal && idx <= length(lines)
        line = strip(lines[idx])
        if !isempty(line) && !startswith(line, "#")
            for token in split(line)
                if !isempty(token)
                    push!(values, parse(Float64, token))
                end
            end
        end
        idx += 1
    end
    
    @assert length(values) >= ntotal "Model file has insufficient values: $(length(values)) < $ntotal"
    
    # Reshape to 3D array (UBC order: E fastest, then N, then Z)
    # UBC stores top-to-bottom for Z, we may need to flip
    density = zeros(Float64, mesh.nx, mesh.ny, mesh.nz)
    vidx = 1
    for k in 1:mesh.nz
        for j in 1:mesh.ny
            for i in 1:mesh.nx
                density[i, j, k] = values[vidx]
                vidx += 1
            end
        end
    end
    
    return GravityModel(mesh, density)
end

"""
    save_model_ubc(filename::String, model::GravityModel; header::String="")

Save density model to UBC format.
"""
function save_model_ubc(filename::String, model::GravityModel; 
                        header::String="# UBC Gravity Model")
    mesh = model.mesh
    
    open(filename, "w") do io
        # Write header
        if !isempty(header)
            println(io, header)
        end
        
        # Write dimensions
        @printf(io, "%d %d %d\n", mesh.nx, mesh.ny, mesh.nz)
        
        # Write cell widths as single representative value (uniform) or first value
        @printf(io, "%.1f %.1f %.1f\n", mesh.dx[1], mesh.dy[1], mesh.dz[1])
        
        # Write density values (UBC order)
        for k in 1:mesh.nz
            for j in 1:mesh.ny
                for i in 1:mesh.nx
                    @printf(io, "%.4f\n", model.density[i, j, k])
                end
            end
        end
    end
    return nothing
end

"""
    load_data_ubc(filename::String)

Load gravity data from UBC format.

# UBC Data Format:
```
! Comment
N_data
E1 N1 Z1 gz1 [error1]
E2 N2 Z2 gz2 [error2]
...
```
"""
function load_data_ubc(filename::String)
    lines = readlines(filename)
    idx = 1
    
    # Skip comments
    while idx <= length(lines) && (startswith(lines[idx], "!") || startswith(lines[idx], "#"))
        idx += 1
    end
    
    # Read number of data points (optional in some formats)
    first_line = split(strip(lines[idx]))
    if length(first_line) == 1
        ndata = parse(Int, first_line[1])
        idx += 1
    else
        # Count non-comment, non-empty lines
        ndata = count(i -> !isempty(strip(lines[i])) && 
                          !startswith(lines[i], "!") && 
                          !startswith(lines[i], "#"), idx:length(lines))
    end
    
    # Read data
    x = Vector{Float64}(undef, ndata)
    y = Vector{Float64}(undef, ndata)
    z = Vector{Float64}(undef, ndata)
    gz = Vector{Float64}(undef, ndata)
    err = zeros(Float64, ndata)
    
    didx = 1
    while didx <= ndata && idx <= length(lines)
        line = strip(lines[idx])
        if !isempty(line) && !startswith(line, "!") && !startswith(line, "#")
            parts = parse.(Float64, split(line))
            x[didx] = parts[1]
            y[didx] = parts[2]
            z[didx] = parts[3]
            gz[didx] = parts[4]
            if length(parts) >= 5
                err[didx] = parts[5]
            end
            didx += 1
        end
        idx += 1
    end
    
    receivers = GravityReceivers(x, y, z, ndata, ndata, 1)
    return GravityData(receivers, gz, err)
end

"""
    save_data_ubc(filename::String, data::GravityData; 
                  header::String="! Gravity Data (UBC Format)")

Save gravity data to UBC format.
"""
function save_data_ubc(filename::String, data::GravityData;
                       header::String="! Gravity Data (UBC Format)")
    open(filename, "w") do io
        println(io, header)
        @printf(io, "%d\n", data.receivers.n)
        
        for i in 1:data.receivers.n
            @printf(io, "%.6f %.6f %.6f %.6e %.6e\n",
                    data.receivers.x[i],
                    data.receivers.y[i],
                    data.receivers.z[i],
                    data.gz[i],
                    data.error[i])
        end
    end
    return nothing
end

"""
    save_data_ubc(filename::String, receivers::GravityReceivers, 
                  gz::Array{Float64,2}; header::String="")

Save 2D array gravity data to UBC format (for regular grids).
"""
function save_data_ubc(filename::String, receivers::GravityReceivers, 
                       gz::Array{Float64,2}; 
                       header::String="! Gravity Data (UBC Format)")
    gz_vec = vec(gz)
    data = GravityData(receivers, gz_vec)
    save_data_ubc(filename, data; header=header)
end

"""
    load_locations_ubc(filename::String)

Load observation locations from UBC format.
"""
function load_locations_ubc(filename::String)
    lines = readlines(filename)
    
    # Skip comments and find data
    locs = Vector{NTuple{3,Float64}}()
    
    for line in lines
        sline = strip(line)
        if !isempty(sline) && !startswith(sline, "!") && !startswith(sline, "#")
            parts = parse.(Float64, split(sline))
            if length(parts) >= 3
                push!(locs, (parts[1], parts[2], parts[3]))
            end
        end
    end
    
    n = length(locs)
    x = [l[1] for l in locs]
    y = [l[2] for l in locs]
    z = [l[3] for l in locs]
    
    return GravityReceivers(x, y, z, n, n, 1)
end

"""
    save_locations_ubc(filename::String, receivers::GravityReceivers)

Save observation locations to UBC format.
"""
function save_locations_ubc(filename::String, receivers::GravityReceivers)
    open(filename, "w") do io
        println(io, "! Observation Locations")
        @printf(io, "%d\n", receivers.n)
        for i in 1:receivers.n
            @printf(io, "%.6f %.6f %.6f\n",
                    receivers.x[i], receivers.y[i], receivers.z[i])
        end
    end
    return nothing
end
