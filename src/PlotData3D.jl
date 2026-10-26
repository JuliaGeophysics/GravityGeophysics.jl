# 3D gravity data plots
# Author: @pankajkmishra
# Maps of observed and predicted gz and their residual in mGal, inversion convergence, and PlotData3D,
# the file entry used by the examples. Gridded surveys are drawn as heatmaps with the stations on top,
# scattered surveys as coloured points. The axes carry longitude and latitude whenever the survey's CRS
# can be transformed, and the stations are thin black outlines so the field underneath stays readable.

using CairoMakie

#---------- maps ----------

# values on the regular grid of a gridded survey, or nothing when the stations are scattered
function _grid_values(d::GravData3D, values)
    es, ns = sort(unique(d.y)), sort(unique(d.x))
    length(es) * length(ns) == ndata(d) || return nothing
    Z = fill(NaN, length(es), length(ns))
    for i in eachindex(values)
        Z[searchsortedfirst(es, d.y[i]), searchsortedfirst(ns, d.x[i])] = values[i]
    end
    any(isnan, Z) ? nothing : (es, ns, Z)
end

# the stations as open circles with a thin black edge
_station_markers!(ax, d; markersize = 5) =
    scatter!(ax, d.y, d.x; color = :transparent, strokecolor = :black, strokewidth = _GRAV_STATION_STROKE,
             markersize = markersize)

"""
    plot_grav3d_data(observed, predicted=nothing; output_path, colormap=Reverse(:Spectral),
                     data_range=nothing, markersize=5) -> path

Map of the observed data; with predictions three maps: observed and predicted on one colour scale, and the residual (predicted − observed) in mGal on its own symmetric scale.
"""
function plot_grav3d_data(observed::GravData3D, predicted::Union{Nothing, GravData3D} = nothing;
                          output_path::AbstractString, colormap = _GRAV_COLORMAP,
                          data_range = nothing, markersize::Real = 5)
    east, north = _data_centre(observed)
    ha, va      = _slice_axes(:xy, observed.crs, east, north)
    vals   = predicted === nothing ? observed.value : vcat(observed.value, predicted.value)
    crange = data_range === nothing ? (a = maximum(abs, vals); (-a, a)) : Tuple(Float64.(data_range))
    crange[1] == crange[2] && (crange = (crange[1] - 1, crange[2] + 1))
    panels = predicted === nothing ? 1 : 3
    fig    = Figure(size = (540 * panels + 120, 560), fontsize = 19)
    function panel!(col, values, title; cmap = colormap, cr = crange, label = rich("g", subscript("z"), " (mGal)"))
        ax   = _grav_axis(fig[1, 2col-1]; _axis_kw(ha, va)..., ylabel = col == 1 ? va.label : "",
                          title, aspect = DataAspect())
        grid = _grid_values(observed, values)
        if grid === nothing
            scatter!(ax, observed.y, observed.x; color = values, colormap = cmap, colorrange = cr,
                     markersize = 4markersize, strokecolor = :black, strokewidth = _GRAV_STATION_STROKE)
        else
            heatmap!(ax, grid[1], grid[2], grid[3]; colormap = cmap, colorrange = cr)
            _station_markers!(ax, observed; markersize)
        end
        Colorbar(fig[1, 2col], colormap = cmap, colorrange = cr, label = label, width = 16, height = Relative(0.8))
    end
    panel!(1, observed.value, "Observed")
    if predicted !== nothing
        ndata(predicted) == ndata(observed) || throw(DimensionMismatch("observed and predicted data differ in length"))
        fit = chi2_rms3d(observed, predicted)
        panel!(2, predicted.value, @sprintf("Predicted, RMS %.2f", fit.rms))
        # the residual is a small fraction of the anomaly, so it takes its own symmetric scale
        res = predicted.value .- observed.value
        a   = maximum(abs, res)
        a > 0 || (a = 1.0)
        panel!(3, res, @sprintf("Residual, max |Δ| %.3g mGal", a);
               cmap = :RdBu, cr = (-a, a), label = "predicted − observed (mGal)")
    end
    _grav_labels!(fig, [(1, 2c - 1) for c in 1:panels])
    _grav_save(output_path, fig)
end

#---------- convergence ----------

"""
    plot_inv3d_convergence(history; output_path, target_rms=nothing) -> path

RMS against the iterations, with the regularization weight β beside it when the algorithm cools one (β is zero throughout an INR run).
"""
function plot_inv3d_convergence(history::AbstractVector; output_path::AbstractString, target_rms = nothing)
    it     = [h.iteration for h in history]
    beta   = [h.beta for h in history]
    cooled = any(>(0), beta)
    fig    = Figure(size = (cooled ? 980 : 560, 420), fontsize = 19)
    ax1    = _grav_axis(fig[1, 1]; xlabel = "Iteration", ylabel = "RMS", yscale = log10)
    many   = length(it) > 60                                # an INR run has one entry per epoch
    many ? lines!(ax1, it, [h.rms for h in history]; color = :black, linewidth = 2.2) :
           scatterlines!(ax1, it, [h.rms for h in history]; color = :black, linewidth = 2.2, markersize = 9)
    target_rms === nothing || target_rms <= 0 ||
        hlines!(ax1, [target_rms]; color = :firebrick, linestyle = :dash, linewidth = 1.8)
    if cooled
        ax2 = _grav_axis(fig[1, 2]; xlabel = "Iteration", ylabel = "β (regularization weight)", yscale = log10)
        scatterlines!(ax2, it, beta; color = RGBf(0.12, 0.38, 0.72), linewidth = 2.2, markersize = 9)
    end
    _grav_save(output_path, fig)
end

#---------- file entry ----------

"""
    PlotData3D(data_file; predicted_file="", output_path=<data stem>.png, kwargs...) -> path

Map of a data file, and with `predicted_file` the observed, predicted and residual maps.
"""
function PlotData3D(data_file::AbstractString; predicted_file::AbstractString = "",
                    output_path::AbstractString = splitext(abspath(isempty(predicted_file) ? data_file : predicted_file))[1] *
                                                  (isempty(predicted_file) ? ".png" : "_fit.png"),
                    kwargs...)
    observed  = ReadData3D(data_file)
    predicted = isempty(predicted_file) ? nothing : ReadData3D(predicted_file)
    plot_grav3d_data(observed, predicted; output_path, kwargs...)
end
