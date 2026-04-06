function _linear_index(i::Int, j::Int, k::Int, nx::Int, ny::Int)
    return i + nx * (j - 1) + nx * ny * (k - 1)
end

function _reshape_density(mesh::GravityMesh, density_vec::AbstractVector{<:Real})
    @assert length(density_vec) == mesh_cell_count(mesh)
    return reshape(Float64.(density_vec), mesh.nx, mesh.ny, mesh.nz)
end

function sensitivity_matrix_prism(mesh::GravityMesh, receivers::GravityReceivers)
    xe, ye, ze = get_cell_edges(mesh)
    ndata = receivers.n
    nmodel = mesh_cell_count(mesh)
    G = Matrix{Float64}(undef, ndata, nmodel)

    col = 1
    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        for ir in 1:ndata
            G[ir, col] = gz_prism(
                receivers.x[ir], receivers.y[ir], receivers.z[ir],
                xe[i], xe[i + 1],
                ye[j], ye[j + 1],
                ze[k], ze[k + 1],
                1.0,
            )
        end
        col += 1
    end
    return G
end

function sensitivity_matrix_prism(model::GravityModel, receivers::GravityReceivers)
    return sensitivity_matrix_prism(model.mesh, receivers)
end

function gradient_operators(mesh::GravityMesh)
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    nmodel = mesh_cell_count(mesh)

    rows_x = Int[]
    cols_x = Int[]
    vals_x = Float64[]
    row_x = 0
    for k in 1:nz, j in 1:ny, i in 1:(nx - 1)
        row_x += 1
        c1 = _linear_index(i, j, k, nx, ny)
        c2 = _linear_index(i + 1, j, k, nx, ny)
        h = 0.5 * (mesh.dx[i] + mesh.dx[i + 1])
        append!(rows_x, (row_x, row_x))
        append!(cols_x, (c1, c2))
        append!(vals_x, (-1.0 / h, 1.0 / h))
    end

    rows_y = Int[]
    cols_y = Int[]
    vals_y = Float64[]
    row_y = 0
    for k in 1:nz, j in 1:(ny - 1), i in 1:nx
        row_y += 1
        c1 = _linear_index(i, j, k, nx, ny)
        c2 = _linear_index(i, j + 1, k, nx, ny)
        h = 0.5 * (mesh.dy[j] + mesh.dy[j + 1])
        append!(rows_y, (row_y, row_y))
        append!(cols_y, (c1, c2))
        append!(vals_y, (-1.0 / h, 1.0 / h))
    end

    rows_z = Int[]
    cols_z = Int[]
    vals_z = Float64[]
    row_z = 0
    for k in 1:(nz - 1), j in 1:ny, i in 1:nx
        row_z += 1
        c1 = _linear_index(i, j, k, nx, ny)
        c2 = _linear_index(i, j, k + 1, nx, ny)
        h = 0.5 * (mesh.dz[k] + mesh.dz[k + 1])
        append!(rows_z, (row_z, row_z))
        append!(cols_z, (c1, c2))
        append!(vals_z, (-1.0 / h, 1.0 / h))
    end

    Dx = sparse(rows_x, cols_x, vals_x, row_x, nmodel)
    Dy = sparse(rows_y, cols_y, vals_y, row_y, nmodel)
    Dz = sparse(rows_z, cols_z, vals_z, row_z, nmodel)
    return Dx, Dy, Dz
end

function depth_weighting(mesh::GravityMesh; beta::Float64 = 1.0, z0::Union{Nothing, Float64} = nothing)
    _, _, zc = get_cell_centers(mesh)
    shift = isnothing(z0) ? 0.5 * minimum(mesh.dz) : z0
    weights = Vector{Float64}(undef, mesh_cell_count(mesh))
    idx = 1
    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        weights[idx] = (zc[k] + shift)^(-beta / 2)
        idx += 1
    end
    return weights
end

function _weighted_system(G::AbstractMatrix, data::GravityData)
    weights = if isempty(data.error) || all(iszero, data.error)
        ones(length(data.gz))
    else
        1.0 ./ max.(Float64.(data.error), eps(Float64))
    end
    Wd = Diagonal(weights)
    return Wd * G, Wd * Float64.(data.gz), weights
end

function invert_tikhonov(data::GravityData, mesh::GravityMesh;
    G::Union{Nothing, AbstractMatrix} = nothing,
    lambda::Float64 = 0.1,
    alpha_s::Float64 = 0.05,
    alpha_xyz::NTuple{3, Float64} = (1.0, 1.0, 1.0),
    depth_beta::Float64 = 0.0,
    reference_model::Union{Nothing, AbstractVector{<:Real}} = nothing,
    positivity::Bool = false)

    Gmat = isnothing(G) ? sensitivity_matrix_prism(mesh, data.receivers) : Matrix{Float64}(G)
    Gw, dw, weights = _weighted_system(Gmat, data)
    nmodel = mesh_cell_count(mesh)
    reference = isnothing(reference_model) ? zeros(nmodel) : Float64.(reference_model)
    @assert length(reference) == nmodel

    Dx, Dy, Dz = gradient_operators(mesh)
    depth_w = depth_beta > 0 ? depth_weighting(mesh; beta = depth_beta) : ones(nmodel)
    Wm = Diagonal(depth_w)

    A = vcat(
        Gw,
        lambda * sqrt(alpha_s) * Wm,
        lambda * sqrt(alpha_xyz[1]) * Matrix(Dx),
        lambda * sqrt(alpha_xyz[2]) * Matrix(Dy),
        lambda * sqrt(alpha_xyz[3]) * Matrix(Dz),
    )
    b = vcat(
        dw,
        lambda * sqrt(alpha_s) * (Wm * reference),
        zeros(size(Dx, 1)),
        zeros(size(Dy, 1)),
        zeros(size(Dz, 1)),
    )

    density_vec = A \ b
    if positivity
        density_vec = max.(density_vec, 0.0)
    end

    model = GravityModel(mesh, _reshape_density(mesh, density_vec))
    predicted = Gmat * density_vec
    info = Dict(
        "predicted_gz" => predicted,
        "rmse_mgal" => sqrt(mean((predicted .- data.gz).^2)) * SI_TO_MGAL,
        "lambda" => lambda,
        "alpha_s" => alpha_s,
        "alpha_xyz" => alpha_xyz,
        "data_weights" => weights,
        "sensitivity" => Gmat,
    )
    return model, info
end

function invert_tikhonov(data::GravityData, model::GravityModel; kwargs...)
    return invert_tikhonov(data, model.mesh; kwargs...)
end