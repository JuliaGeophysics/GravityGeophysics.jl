function _map_plot_axes(xvals::AbstractVector{<:Real}, yvals::AbstractVector{<:Real})
    return Float64.(xvals), Float64.(yvals), "Easting (m)", "Northing (m)"
end

function _density_colorbar_label(use_bulk_density::Bool)
    return use_bulk_density ? "Density (kg/m^3)" : "Density contrast (kg/m^3)"
end

function _density_colormap()
    return :Spectral
end

function _save_figure(fig, save_path::String)
    if !isempty(save_path)
        CairoMakie.activate!()
        CairoMakie.save(save_path, fig)
    end
    return fig
end

function plot_gravity_data(data::GravityData; title::String = "Gravity Data", unit::Symbol = :mgal, save_path::String = "")
    values = unit == :mgal ? data.gz .* SI_TO_MGAL : data.gz
    unit_label = unit == :mgal ? "mGal" : "m/s^2"
    fig = CairoMakie.Figure(size = (1100, 850), backgroundcolor = :white)
    rx_plot, ry_plot, xlabel, ylabel = _map_plot_axes(data.receivers.x, data.receivers.y)

    if data.receivers.nx > 1 && data.receivers.ny > 1 && data.receivers.n == data.receivers.nx * data.receivers.ny
        rx = rx_plot[1:data.receivers.nx]
        ry = [ry_plot[(iy - 1) * data.receivers.nx + 1] for iy in 1:data.receivers.ny]
        grid = reshape(values, data.receivers.nx, data.receivers.ny)
        ax = CairoMakie.Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, title = title, aspect = CairoMakie.DataAspect())
        hm = CairoMakie.heatmap!(ax, rx, ry, grid'; colormap = :viridis)
        CairoMakie.Colorbar(fig[1, 2], hm, label = unit_label, width = 18)
    else
        ax = CairoMakie.Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, title = title, aspect = CairoMakie.DataAspect())
        sc = CairoMakie.scatter!(ax, rx_plot, ry_plot; color = values, markersize = 10, colormap = :viridis)
        CairoMakie.Colorbar(fig[1, 2], sc, label = unit_label, width = 18)
    end
    return _save_figure(fig, save_path)
end

function plot_receiver_locations(receivers::GravityReceivers; title::String = "Receiver Locations", save_path::String = "")
    xplot, yplot, xlabel, ylabel = _map_plot_axes(receivers.x, receivers.y)
    fig = CairoMakie.Figure(size = (1100, 850), backgroundcolor = :white)
    ax = CairoMakie.Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, title = title, aspect = CairoMakie.DataAspect())
    CairoMakie.scatter!(ax, xplot, yplot; color = :black, markersize = 8)
    return _save_figure(fig, save_path)
end

function plot_gravity_comparison(receivers::GravityReceivers, gz_prism::Vector{Float64}, gz_fd::Vector{Float64}; title_prefix::String = "Gravity", save_path::String = "")
    nrx, nry = receivers.nx, receivers.ny
    gz_pr_2d = reshape(gz_prism, nrx, nry)
    gz_fd_2d = reshape(gz_fd, nrx, nry)
    xplot, yplot, xlabel, ylabel = _map_plot_axes(receivers.x, receivers.y)
    rx = xplot[1:nrx]
    ry = [yplot[(i - 1) * nrx + 1] for i in 1:nry]
    gz_pr_mgal = gz_pr_2d .* SI_TO_MGAL
    gz_fd_mgal = gz_fd_2d .* SI_TO_MGAL
    residual = gz_fd_mgal .- gz_pr_mgal
    vmin = min(minimum(gz_pr_mgal), minimum(gz_fd_mgal))
    vmax = max(maximum(gz_pr_mgal), maximum(gz_fd_mgal))
    res_max = maximum(abs.(residual))
    mid_y = div(nry, 2) + 1

    fig = CairoMakie.Figure(size = (1600, 1200), backgroundcolor = :white)
    CairoMakie.Label(fig[0, 1:2], "$title_prefix: Prism vs FD", fontsize = 28)
    ax1 = CairoMakie.Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, title = "Prism", aspect = CairoMakie.DataAspect())
    ax2 = CairoMakie.Axis(fig[1, 2], xlabel = xlabel, ylabel = ylabel, title = "FD", aspect = CairoMakie.DataAspect())
    ax3 = CairoMakie.Axis(fig[2, 1], xlabel = xlabel, ylabel = ylabel, title = "FD - Prism", aspect = CairoMakie.DataAspect())
    ax4 = CairoMakie.Axis(fig[2, 2], xlabel = xlabel, ylabel = "gz (mGal)", title = "Center Profile")
    hm1 = CairoMakie.heatmap!(ax1, rx, ry, gz_pr_mgal'; colormap = :viridis, colorrange = (vmin, vmax))
    hm2 = CairoMakie.heatmap!(ax2, rx, ry, gz_fd_mgal'; colormap = :viridis, colorrange = (vmin, vmax))
    hm3 = CairoMakie.heatmap!(ax3, rx, ry, residual'; colormap = :balance, colorrange = (-res_max, res_max))
    CairoMakie.lines!(ax4, rx, gz_pr_mgal[:, mid_y]; color = :black, linewidth = 3, label = "Prism")
    CairoMakie.lines!(ax4, rx, gz_fd_mgal[:, mid_y]; color = :firebrick3, linewidth = 3, linestyle = :dash, label = "FD")
    CairoMakie.axislegend(ax4, position = :rb)
    CairoMakie.Colorbar(fig[1, 3], hm1, label = "mGal", width = 16)
    CairoMakie.Colorbar(fig[2, 3], hm3, label = "mGal", width = 16)
    return _save_figure(fig, save_path)
end

function plot_gravity_comparison(data_prism::GravityData, data_fd::GravityData; kwargs...)
    return plot_gravity_comparison(data_prism.receivers, data_prism.gz, data_fd.gz; kwargs...)
end

function _extract_model_slice(model::GravityModel, slice_axis::Symbol, slice_index::Int)
    mesh = model.mesh
    xc, yc, zc = get_cell_centers(mesh)
    xc_plot, yc_plot, map_xlabel, map_ylabel = _map_plot_axes(xc, yc)

    if slice_axis == :z
        idx = slice_index == 0 ? div(mesh.nz, 2) + 1 : clamp(slice_index, 1, mesh.nz)
        return slice_data = model.density[:, :, idx], xc_plot, yc_plot, idx, map_xlabel, map_ylabel, "Density at depth $(round(zc[idx], digits = 1)) m", false
    elseif slice_axis == :y
        idx = slice_index == 0 ? div(mesh.ny, 2) + 1 : clamp(slice_index, 1, mesh.ny)
        return slice_data = model.density[:, idx, :], xc_plot, zc, idx, map_xlabel, "Depth (m)", "Density at N=$(round(yc_plot[idx], digits = 1)) m", true
    else
        idx = slice_index == 0 ? div(mesh.nx, 2) + 1 : clamp(slice_index, 1, mesh.nx)
        return slice_data = model.density[idx, :, :], yc_plot, zc, idx, map_ylabel, "Depth (m)", "Density at E=$(round(xc_plot[idx], digits = 1)) m", true
    end
end

function plot_model_slice(model::GravityModel; slice_axis::Symbol = :z, slice_index::Int = 0, save_path::String = "", title::String = "", clims = nothing, use_bulk_density::Bool = false)
    slice_data, axis1, axis2, _, xlabel, ylabel, default_title, yflip = _extract_model_slice(model, slice_axis, slice_index)
    plot_title = isempty(title) ? default_title : title
    fig = CairoMakie.Figure(size = (1100, 850), backgroundcolor = :white)
    ax = CairoMakie.Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, title = plot_title, aspect = slice_axis == :z ? CairoMakie.DataAspect() : nothing)
    ax.yreversed = yflip
    hm = CairoMakie.heatmap!(ax, axis1, axis2, slice_data'; colormap = _density_colormap(), colorrange = clims)
    CairoMakie.Colorbar(fig[1, 2], hm, label = _density_colorbar_label(use_bulk_density), width = 18)
    return _save_figure(fig, save_path)
end

function plot_inversion_comparison(true_model::GravityModel, standard_model::GravityModel, inr_model::GravityModel; slice_axis::Symbol = :z, slice_index::Int = 0, save_path::String = "", use_bulk_density::Bool = false)
    true_slice, _, _, _, _, _, _, _ = _extract_model_slice(true_model, slice_axis, slice_index)
    standard_slice, _, _, _, _, _, _, _ = _extract_model_slice(standard_model, slice_axis, slice_index)
    inr_slice, _, _, _, _, _, _, _ = _extract_model_slice(inr_model, slice_axis, slice_index)
    density_limits = extrema(vcat(vec(true_slice), vec(standard_slice), vec(inr_slice)))

    fig = CairoMakie.Figure(size = (1800, 650), backgroundcolor = :white)
    CairoMakie.Label(fig[0, 1:3], "True vs Standard vs INR Inversion", fontsize = 32)

    function add_panel!(gridpos, model, panel_title, limits)
        slice_data, axis1, axis2, _, xlabel, ylabel, _, yflip = _extract_model_slice(model, slice_axis, slice_index)
        ax = CairoMakie.Axis(gridpos, xlabel = xlabel, ylabel = ylabel, title = panel_title, aspect = slice_axis == :z ? CairoMakie.DataAspect() : nothing)
        ax.yreversed = yflip
        hm = CairoMakie.heatmap!(ax, axis1, axis2, slice_data'; colormap = _density_colormap(), colorrange = limits)
        return hm
    end

    hm1 = add_panel!(fig[1, 1], true_model, "True Model", density_limits)
    hm2 = add_panel!(fig[1, 2], standard_model, "Standard Tikhonov", density_limits)
    hm3 = add_panel!(fig[1, 3], inr_model, "INR Reconstruction", density_limits)
    density_label = _density_colorbar_label(use_bulk_density)
    CairoMakie.Colorbar(fig[1, 4], hm3, label = density_label, width = 18)
    return _save_figure(fig, save_path)
end

function plot_orthogonal_inversion_comparison(true_model::GravityModel, standard_model::GravityModel, inr_model::GravityModel;
    x_slice_index::Int = 0,
    y_slice_index::Int = 0,
    z_slice_index::Int = 0,
    save_path::String = "",
    use_bulk_density::Bool = false,
    clims = nothing,
    residual_clims = nothing)

    mesh = true_model.mesh
    x_idx = x_slice_index == 0 ? div(mesh.nx, 2) + 1 : clamp(x_slice_index, 1, mesh.nx)
    y_idx = y_slice_index == 0 ? div(mesh.ny, 2) + 1 : clamp(y_slice_index, 1, mesh.ny)
    z_idx = z_slice_index == 0 ? div(mesh.nz, 2) + 1 : clamp(z_slice_index, 1, mesh.nz)

    sections = ((:z, z_idx, "XY"), (:y, y_idx, "XZ"), (:x, x_idx, "YZ"))
    models = (("True Model", true_model), ("Standard Tikhonov", standard_model), ("INR Reconstruction", inr_model))
    density_label = _density_colorbar_label(use_bulk_density)

    all_density_values = Float64[]
    for (_, model) in models
        for (axis, idx, _) in sections
            slice_data, _, _, _, _, _, _, _ = _extract_model_slice(model, axis, idx)
            append!(all_density_values, vec(Float64.(slice_data)))
        end
    end
    density_limits = isnothing(clims) ? extrema(all_density_values) : clims

    fig = CairoMakie.Figure(size = (2100, 1350), backgroundcolor = :white)
    CairoMakie.Label(fig[0, 1:4], "Orthogonal Inversion Comparison", fontsize = 32)

    for (row, (row_title, model)) in enumerate(models)
        for (col, (axis, idx, plane_label)) in enumerate(sections)
            slice_data, axis1, axis2, _, xlabel, ylabel, _, yflip = _extract_model_slice(model, axis, idx)
            ax = CairoMakie.Axis(fig[row, col], xlabel = xlabel, ylabel = ylabel, title = "$row_title: $plane_label", aspect = axis == :z ? CairoMakie.DataAspect() : nothing)
            ax.yreversed = yflip
            hm = CairoMakie.heatmap!(ax, axis1, axis2, slice_data'; colormap = _density_colormap(), colorrange = density_limits)
            if col == 3
                CairoMakie.Colorbar(fig[row, 4], hm, label = density_label, width = 18)
            end
        end
    end

    return _save_figure(fig, save_path)
end

function plot_convergence(history; save_path::String = "")
    loss = get(history, "loss_history", Float64[])
    misfit = get(history, "misfit_history", Float64[])
    smooth = get(history, "smoothness_history", Float64[])
    fig = CairoMakie.Figure(size = (1100, 800), backgroundcolor = :white)
    ax = CairoMakie.Axis(fig[1, 1], xlabel = "Epoch", ylabel = "Objective", title = "INR Convergence", yscale = log10)
    if !isempty(loss)
        CairoMakie.lines!(ax, 1:length(loss), loss; color = :dodgerblue3, linewidth = 3, label = "Loss")
    end
    if !isempty(misfit)
        CairoMakie.lines!(ax, 1:length(misfit), misfit; color = :orangered3, linewidth = 3, linestyle = :dash, label = "Misfit")
    end
    if !isempty(smooth)
        CairoMakie.lines!(ax, 1:length(smooth), smooth; color = :darkgreen, linewidth = 2, linestyle = :dot, label = "Smoothness")
    end
    CairoMakie.axislegend(ax, position = :rt)
    return _save_figure(fig, save_path)
end