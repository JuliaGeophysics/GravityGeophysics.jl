# Author: @pankajkmishra
"""
    Visualization.jl - Plotting Functions for Gravity Forward Modeling

Provides plotting utilities for gravity data comparison and model visualization.
Requires Plots.jl to be loaded.
"""

"""
    plot_gravity_comparison(receivers::GravityReceivers,
                            gz_prism::Vector{Float64},
                            gz_fd::Vector{Float64};
                            title_prefix::String="Gravity",
                            save_path::String="")

Create a 2x2 comparison plot of Prism vs FD gravity results.

# Arguments
- `receivers`: Receiver locations (must be regular grid)
- `gz_prism`: Prism method results [m/s²]
- `gz_fd`: FD method results [m/s²]
- `title_prefix`: Prefix for plot titles
- `save_path`: If non-empty, save plot to this path

# Returns
- `fig`: Plot figure object
"""
function plot_gravity_comparison(receivers::GravityReceivers,
                                 gz_prism::Vector{Float64},
                                 gz_fd::Vector{Float64};
                                 title_prefix::String="Gravity",
                                 save_path::String="")
    # Reshape to 2D
    nrx, nry = receivers.nx, receivers.ny
    gz_pr_2d = reshape(gz_prism, nrx, nry)
    gz_fd_2d = reshape(gz_fd, nrx, nry)
    
    # Get unique x and y coordinates
    rx = receivers.x[1:nrx]
    ry = [receivers.y[(i-1)*nrx + 1] for i in 1:nry]
    
    # Convert to km and mGal
    rx_km = rx ./ 1000
    ry_km = ry ./ 1000
    gz_pr_mGal = gz_pr_2d .* SI_TO_MGAL
    gz_fd_mGal = gz_fd_2d .* SI_TO_MGAL
    residual = gz_fd_mGal .- gz_pr_mGal
    
    # Color limits
    vmin = min(minimum(gz_pr_mGal), minimum(gz_fd_mGal))
    vmax = max(maximum(gz_pr_mGal), maximum(gz_fd_mGal))
    
    # Create plots
    p1 = heatmap(rx_km, ry_km, gz_pr_mGal',
                 xlabel="X (km)", ylabel="Y (km)",
                 title="Prism (Analytical)",
                 c=:viridis, clims=(vmin, vmax),
                 aspect_ratio=:equal,
                 colorbar_title="mGal")
    
    p2 = heatmap(rx_km, ry_km, gz_fd_mGal',
                 xlabel="X (km)", ylabel="Y (km)",
                 title="FD (Numerical)",
                 c=:viridis, clims=(vmin, vmax),
                 aspect_ratio=:equal,
                 colorbar_title="mGal")
    
    res_max = maximum(abs.(residual))
    p3 = heatmap(rx_km, ry_km, residual',
                 xlabel="X (km)", ylabel="Y (km)",
                 title="Residual (FD - Prism)",
                 c=:RdBu, clims=(-res_max, res_max),
                 aspect_ratio=:equal,
                 colorbar_title="mGal")
    
    # Profile at Y=0 (center row)
    mid_y = div(nry, 2) + 1
    rmse_mGal = sqrt(mean((gz_fd_mGal .- gz_pr_mGal).^2))
    
    p4 = plot(rx_km, gz_pr_mGal[:, mid_y], lw=2, label="Prism",
              xlabel="X (km)", ylabel="gz (mGal)",
              title="Profile at Y=0 (RMSE=$(round(rmse_mGal, digits=4)) mGal)")
    plot!(p4, rx_km, gz_fd_mGal[:, mid_y], lw=2, ls=:dash, label="FD")
    
    # Combine
    fig = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 900),
               plot_title="$title_prefix: Prism vs FD Comparison")
    
    if !isempty(save_path)
        savefig(fig, save_path)
        @printf("Plot saved to %s\n", save_path)
    end
    
    return fig
end

"""
    plot_gravity_comparison(data_prism::GravityData, data_fd::GravityData; kwargs...)

Convenience method using GravityData objects.
"""
function plot_gravity_comparison(data_prism::GravityData, data_fd::GravityData; kwargs...)
    return plot_gravity_comparison(data_prism.receivers, data_prism.gz, data_fd.gz; kwargs...)
end

"""
    plot_model_slice(model::GravityModel; 
                     slice_axis::Symbol=:z,
                     slice_index::Int=0,
                     save_path::String="")

Plot a 2D slice of the density model.

# Arguments
- `model`: GravityModel to plot
- `slice_axis`: Axis perpendicular to slice (:x, :y, or :z)
- `slice_index`: Index along slice axis (0 = center)
- `save_path`: If non-empty, save plot

# Returns
- `fig`: Plot figure object
"""
function plot_model_slice(model::GravityModel;
                          slice_axis::Symbol=:z,
                          slice_index::Int=0,
                          save_path::String="")
    mesh = model.mesh
    xc, yc, zc = get_cell_centers(mesh)
    
    # Convert to km
    xc_km = xc ./ 1000
    yc_km = yc ./ 1000
    zc_km = zc ./ 1000
    
    if slice_axis == :z
        idx = slice_index == 0 ? div(mesh.nz, 2) + 1 : slice_index
        slice_data = model.density[:, :, idx]
        fig = heatmap(xc_km, yc_km, slice_data',
                      xlabel="X (km)", ylabel="Y (km)",
                      title="Density at z = $(round(zc_km[idx], digits=2)) km",
                      c=:viridis,
                      colorbar_title="kg/m³",
                      aspect_ratio=:equal)
    elseif slice_axis == :y
        idx = slice_index == 0 ? div(mesh.ny, 2) + 1 : slice_index
        slice_data = model.density[:, idx, :]
        fig = heatmap(xc_km, zc_km, slice_data',
                      xlabel="X (km)", ylabel="Z (km)",
                      title="Density at y = $(round(yc_km[idx], digits=2)) km",
                      c=:viridis,
                      colorbar_title="kg/m³",
                      yflip=true)  # Z positive down
    else  # :x
        idx = slice_index == 0 ? div(mesh.nx, 2) + 1 : slice_index
        slice_data = model.density[idx, :, :]
        fig = heatmap(yc_km, zc_km, slice_data',
                      xlabel="Y (km)", ylabel="Z (km)",
                      title="Density at x = $(round(xc_km[idx], digits=2)) km",
                      c=:viridis,
                      colorbar_title="kg/m³",
                      yflip=true)
    end
    
    if !isempty(save_path)
        savefig(fig, save_path)
        @printf("Plot saved to %s\n", save_path)
    end
    
    return fig
end

"""
    plot_residual_histogram(residual::Vector{Float64}; 
                            unit::String="mGal",
                            save_path::String="")

Plot histogram of residuals.
"""
function plot_residual_histogram(residual::Vector{Float64};
                                 unit::String="mGal",
                                 save_path::String="")
    fig = histogram(residual, bins=50,
                    xlabel="Residual ($unit)",
                    ylabel="Count",
                    title="Residual Distribution",
                    legend=false,
                    alpha=0.7)
    vline!(fig, [mean(residual)], lw=2, ls=:dash, color=:red)
    
    if !isempty(save_path)
        savefig(fig, save_path)
        @printf("Histogram saved to %s\n", save_path)
    end
    
    return fig
end
