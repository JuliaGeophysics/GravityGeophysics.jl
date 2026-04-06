function _normalize_coordinates(mesh::GravityMesh)
    xc, yc, zc = get_cell_centers(mesh)
    xmean = mean(xc)
    ymean = mean(yc)
    zmean = mean(zc)
    xscale = max(std(xc), eps(Float64))
    yscale = max(std(yc), eps(Float64))
    zscale = max(std(zc), eps(Float64))

    coords = Matrix{Float32}(undef, 3, mesh_cell_count(mesh))
    idx = 1
    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        coords[1, idx] = Float32((xc[i] - xmean) / xscale)
        coords[2, idx] = Float32((yc[j] - ymean) / yscale)
        coords[3, idx] = Float32((zc[k] - zmean) / zscale)
        idx += 1
    end
    return coords
end

function positional_encoding(coords::AbstractMatrix{<:Real}; nfreq::Int = 2)
    x = Float32.(coords)
    features = [x]
    for k in 0:(nfreq - 1)
        scale = Float32(2.0^k)
        push!(features, sin.(scale .* x))
        push!(features, cos.(scale .* x))
    end
    return vcat(features...)
end

leakyrelu01(x) = Flux.leakyrelu(x, 0.01f0)

function build_inr_network(input_dim::Int; hidden_dim::Int = 256, depth::Int = 4)
    layers = Any[Flux.Dense(input_dim, hidden_dim, leakyrelu01)]
    for _ in 2:depth
        push!(layers, Flux.Dense(hidden_dim, hidden_dim, leakyrelu01))
    end
    final_layer = Flux.Dense(hidden_dim, 1)
    final_layer.weight .*= 0.01f0
    final_layer.bias .= 0.0f0
    push!(layers, final_layer)
    return Flux.Chain(layers...)
end

function invert_inr(data::GravityData, mesh::GravityMesh;
    G::Union{Nothing, AbstractMatrix} = nothing,
    epochs::Int = 500,
    learning_rate::Float64 = 1e-2,
    hidden_dim::Int = 256,
    depth::Int = 4,
    nfreq::Int = 2,
    max_density::Float64 = 600.0,
    seed::Int = 42,
    verbose::Bool = true)

    Random.seed!(seed)
    Gmat = isnothing(G) ? sensitivity_matrix_prism(mesh, data.receivers) : Matrix{Float64}(G)
    _, _, weights = _weighted_system(Gmat, data)
    observed = Float32.(data.gz)
    w = Float32.(weights)

    coords = _normalize_coordinates(mesh)
    features = positional_encoding(coords; nfreq = nfreq)
    network = build_inr_network(size(features, 1); hidden_dim = hidden_dim, depth = depth)
    state = Flux.setup(Flux.Adam(learning_rate), network)

    loss_history = Float64[]
    misfit_history = Float64[]

    function density_prediction(net)
        response = vec(net(features))
        return Float32(max_density) .* tanh.(response)
    end

    for epoch in 1:epochs
        loss, grads = Flux.withgradient(network) do net
            density = density_prediction(net)
            predicted = Float32.(Gmat) * density
            residual = w .* (predicted .- observed)
            mean(abs2, residual)
        end
        Flux.update!(state, network, grads[1])

        density = density_prediction(network)
        predicted = Float32.(Gmat) * density
        residual = w .* (predicted .- observed)
        misfit = mean(abs2, residual)

        push!(loss_history, Float64(loss))
        push!(misfit_history, Float64(misfit))

        if verbose && (epoch == 1 || epoch % max(10, epochs ÷ 10) == 0 || epoch == epochs)
            @printf("INR epoch %d/%d | loss=%.6e | misfit=%.6e\n", epoch, epochs, Float64(loss), Float64(misfit))
        end
    end

    density_vec = Float64.(density_prediction(network))
    predicted = Float64.(Gmat * density_vec)
    model = GravityModel(mesh, _reshape_density(mesh, density_vec))
    info = Dict(
        "predicted_gz" => predicted,
        "loss_history" => loss_history,
        "misfit_history" => misfit_history,
        "rmse_mgal" => sqrt(mean((predicted .- data.gz).^2)) * SI_TO_MGAL,
        "network" => network,
        "features" => features,
        "sensitivity" => Gmat,
    )
    return model, info
end

function compare_inversion_methods(observed::GravityData, true_model::GravityModel;
    G::Union{Nothing, AbstractMatrix} = nothing,
    tikhonov_kwargs::NamedTuple = (;),
    inr_kwargs::NamedTuple = (;))

    Gmat = isnothing(G) ? sensitivity_matrix_prism(true_model.mesh, observed.receivers) : Matrix{Float64}(G)
    standard_model, standard_info = invert_tikhonov(observed, true_model.mesh; G = Gmat, pairs(tikhonov_kwargs)...)
    inr_model, inr_info = invert_inr(observed, true_model.mesh; G = Gmat, pairs(inr_kwargs)...)

    true_vec = vec(true_model.density)
    standard_vec = vec(standard_model.density)
    inr_vec = vec(inr_model.density)

    comparison = Dict(
        "standard_rmse_model" => sqrt(mean((standard_vec .- true_vec).^2)),
        "inr_rmse_model" => sqrt(mean((inr_vec .- true_vec).^2)),
        "standard_rmse_data_mgal" => standard_info["rmse_mgal"],
        "inr_rmse_data_mgal" => inr_info["rmse_mgal"],
        "standard" => standard_info,
        "inr" => inr_info,
    )
    return standard_model, inr_model, comparison
end