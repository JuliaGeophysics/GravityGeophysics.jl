function _metadata_line(io::IO, key::AbstractString, value)
    println(io, key * "=" * string(value))
end

function _mesh_from_centers(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    x_edges = edges_from_centers(x)
    y_edges = edges_from_centers(y)
    z_edges = edges_from_centers(z)
    return GravityMesh(
        length(x),
        length(y),
        length(z),
        collect(diff(x_edges)),
        collect(diff(y_edges)),
        collect(diff(z_edges)),
        x_edges[1],
        y_edges[1],
        z_edges[1],
        x_edges[end] - x_edges[1],
        y_edges[end] - y_edges[1],
        z_edges[end] - z_edges[1],
    )
end

function save_model_vox(path::AbstractString, model::GravityModel; name::AbstractString = "Gravity Density Model", kind::Symbol = :density, units::AbstractString = "kg/m^3")
    mkpath(dirname(path))
    xc, yc, zc = get_cell_centers(model.mesh)
    open(path, "w") do io
        println(io, "# GravityGeophysics voxel volume")
        _metadata_line(io, "name", name)
        _metadata_line(io, "kind", String(kind))
        _metadata_line(io, "units", units)
        _metadata_line(io, "dims", "$(length(xc)) $(length(yc)) $(length(zc))")
        _metadata_line(io, "x", join(xc, ' '))
        _metadata_line(io, "y", join(yc, ' '))
        _metadata_line(io, "z", join(zc, ' '))
        for k in 1:model.mesh.nz
            println(io, "slice=$k")
            for j in 1:model.mesh.ny
                println(io, join(model.density[:, j, k], ' '))
            end
        end
    end
    return path
end

function load_model_vox(path::AbstractString)
    lines = readlines(path)
    meta = Dict{String, String}()
    data_start = 1
    for (idx, line) in enumerate(lines)
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "slice=")
            data_start = idx
            break
        elseif startswith(stripped, "#")
            continue
        elseif occursin('=', stripped)
            key, value = split(stripped, '='; limit = 2)
            meta[strip(key)] = strip(value)
        end
    end

    dims = parse.(Int, split(meta["dims"]))
    x = parse.(Float64, split(meta["x"]))
    y = parse.(Float64, split(meta["y"]))
    z = parse.(Float64, split(meta["z"]))
    density = zeros(Float64, dims[1], dims[2], dims[3])

    current_k = 0
    current_j = 0
    for line in lines[data_start:end]
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "slice=")
            current_k = parse(Int, split(stripped, '='; limit = 2)[2])
            current_j = 0
            continue
        end
        current_j += 1
        density[:, current_j, current_k] = parse.(Float64, split(stripped))
    end

    mesh = _mesh_from_centers(x, y, z)
    return GravityModel(mesh, density)
end

function save_receivers_xyz(path::AbstractString, receivers::GravityReceivers)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# x y z")
        for i in 1:receivers.n
            @printf(io, "%.6f %.6f %.6f\n", receivers.x[i], receivers.y[i], receivers.z[i])
        end
    end
    return path
end

function save_data_xyz(path::AbstractString, data::GravityData; include_error::Bool = true)
    mkpath(dirname(path))
    open(path, "w") do io
        if include_error
            println(io, "# x y z gz_si error_si")
        else
            println(io, "# x y z gz_si")
        end
        for i in 1:data.receivers.n
            if include_error
                @printf(io, "%.6f %.6f %.6f %.12e %.12e\n", data.receivers.x[i], data.receivers.y[i], data.receivers.z[i], data.gz[i], data.error[i])
            else
                @printf(io, "%.6f %.6f %.6f %.12e\n", data.receivers.x[i], data.receivers.y[i], data.receivers.z[i], data.gz[i])
            end
        end
    end
    return path
end

function load_data_xyz(path::AbstractString)
    x = Float64[]
    y = Float64[]
    z = Float64[]
    gz = Float64[]
    err = Float64[]

    for line in readlines(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, '#') && continue
        parts = split(stripped)
        length(parts) < 4 && continue
        push!(x, parse(Float64, parts[1]))
        push!(y, parse(Float64, parts[2]))
        push!(z, parse(Float64, parts[3]))
        push!(gz, parse(Float64, parts[4]))
        push!(err, length(parts) >= 5 ? parse(Float64, parts[5]) : 0.0)
    end

    receivers = create_receivers(x, y, z)
    return GravityData(receivers, gz, err)
end