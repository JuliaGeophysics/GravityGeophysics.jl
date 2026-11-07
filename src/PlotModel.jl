# Shared plotting utilities for the model and data viewers.
# Author: @pankajkmishra
# Property label and colour range of density-contrast plots, model arrays for the viewers, and the window:
# the viewers are Makie blocks, drawn headless with CairoMakie (`interactive = false`) or opened with GLMakie.

# colour-bar label of density contrast, linear or log10
_value_label(log10scale::Bool) = log10scale ? "log₁₀ Δρ (kg/m³)" : "Δρ (kg/m³)"

"""
    compute_colorrange(R; value_range=nothing) -> (cmin, cmax)

`value_range` when given, otherwise the 2-98 % range of the finite values of `R` (air is NaN).
"""
function compute_colorrange(R::AbstractArray;
        value_range::Union{Nothing, Tuple{<:Real,<:Real}, AbstractVector{<:Real}} = nothing)
    if isnothing(value_range)
        vals = R[isfinite.(R)]
        isempty(vals) && (vals = [0.0, 1.0])
        qlo, qhi = quantile(vec(vals), (0.02, 0.98))
        cmin, cmax = min(qlo, qhi), max(qlo, qhi)
        if cmin == cmax
            ϵ = max(1e-12, 1e-6 * abs(cmin))
            cmin -= ϵ; cmax += ϵ
        end
    else
        lo, hi = Float64(value_range[1]), Float64(value_range[2])
        (!isfinite(lo) || !isfinite(hi)) && error("value_range must contain finite numbers")
        lo > hi && ((lo, hi) = (hi, lo))
        if lo == hi
            ϵ = max(1e-12, 1e-6 * abs(lo))
            lo -= ϵ; hi += ϵ
        end
        cmin, cmax = lo, hi
    end
    return cmin, cmax
end

# Depth below the top of the topography (the top of the highest earth cell), positive down: cell centres,
# layer edges and the z of that top in the model frame. Air above it has negative depth.
function _depths(M)
    k0 = something(findfirst(k -> any(isfinite, view(M.A, :, :, k)), axes(M.A, 3)), 1)
    z0 = M.z[k0]
    M.cz .- z0, M.z .- z0, z0
end

# layers whose top lies above `max_depth` (m below the top of the topography); all when nothing
_depth_layers(ze, max_depth) = isnothing(max_depth) ? (1:length(ze)-1) :
    (1:clamp(searchsortedfirst(ze, float(max_depth)) - 1, 1, length(ze) - 1))

"""
    prepare_model_arrays(M; log10scale=false, withPadding=true, max_depth=nothing, pad_tol=0.2)

Cell centres and values of the part of `M` shown: `(x, y, z, R, ix, iy, kz)`, z the depth below the top of
the topography (positive down); `max_depth` counts from there.
"""
function prepare_model_arrays(M;
        log10scale::Bool = false,
        withPadding::Bool = true,
        max_depth::Union{Nothing, Real} = nothing,
        pad_tol::Real = 0.2)
    x_all = M.cx
    y_all = M.cy
    z_all, ze_all, _ = _depths(M)
    A_all = log10scale ? log10.(M.A) : copy(M.A)

    ix = withPadding ? (1:length(x_all)) : core_indices(x_all; tol = pad_tol)
    iy = withPadding ? (1:length(y_all)) : core_indices(y_all; tol = pad_tol)
    kz = _depth_layers(ze_all, max_depth)

    return x_all[ix], y_all[iy], z_all[kz], A_all[ix, iy, kz], ix, iy, kz
end

#---------- window ----------

# open `fig` in a GLMakie window and wait until it is closed
function _open_window(fig)
    GLMakie.activate!()
    Makie.update_state_before_display!(fig)
    screen = GLMakie.Screen(fig.scene)
    if screen isa GLMakie.Screen && GLMakie.requires_update(screen)
        GLMakie.render_frame(screen)
    end
    println("\nClose the figure window to exit...")
    wait(screen)
    return nothing
end

"""
    save_figure(path, fig; px_per_unit=2)

Write a figure (PNG, PDF or SVG by extension) with CairoMakie.
"""
save_figure(path::AbstractString, fig; px_per_unit::Real = 2) =
    (save(path, fig; px_per_unit, backend = CairoMakie); String(path))
