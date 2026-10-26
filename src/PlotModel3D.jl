# Interactive 3D density model viewers
# Author: @pankajkmishra
# GLMakie viewers for a density model: depth slices (XY), north-south and east-west sections (XZ and YZ) and
# the three orthogonal planes together (XYZ). Each has a slider with Prev/Next/First/Last, a core/padding
# toggle, a zoom reset and an Export Figure button that writes a PNG through CairoMakie
# PlotModelXY, PlotModelXZ, PlotModelYZ and PlotModelXYZ are the file entries used by the examples
# Every entry takes one model or a list of label => model pairs, which are drawn side by side on one slider
# and one colour scale, so a true model and the models recovered by two algorithms share every image
# Depth is drawn as -z so the mesh top is at the top of every section, with positive depths on the ticks
# Air cells are NaN in a model and are left transparent, so the topography shows as the edge of the colour

using CairoMakie
using Dates

const _GRAV_DENSITY_LABEL = "Δρ (kg/m³)"
const _GRAV_FULL_INFO     = "View: full mesh, padding included"
const _GRAV_CORE_INFO     = "View: core mesh, padding dropped"

#---------- shared preparation ----------

# centres, edges and values of the full and the core model over the depths to plot
function _viewer_arrays(M::GravModel3D; max_depth = nothing, pad_tol::Real = 0.2)
    ix_full, iy_full = 1:length(M.cx), 1:length(M.cy)
    ix_core, iy_core = lateral_core_ranges(M; tol = pad_tol)
    kz = max_depth === nothing ? (1:length(M.cz)) : z_indices_for_max_depth(M.cz, float(max_depth))
    edges(e, r) = e[first(r):last(r)+1]
    (; kz,
       z       = M.cz[kz],
       z_edges = edges(M.z, kz),
       x_full  = M.cx[ix_full], y_full  = M.cy[iy_full],
       x_core  = M.cx[ix_core], y_core  = M.cy[iy_core],
       xe_full = edges(M.x, ix_full), ye_full = edges(M.y, iy_full),
       xe_core = edges(M.x, ix_core), ye_core = edges(M.y, iy_core),
       R_full  = M.A[ix_full, iy_full, kz],
       R_core  = M.A[ix_core, iy_core, kz],
       ix_full, iy_full, ix_core, iy_core)
end

# one label, model and set of prepared arrays per panel, every model checked against the first mesh
function _viewer_panels(models::AbstractVector{<:Pair{String, GravModel3D}}; max_depth = nothing,
                        pad_tol::Real = 0.2)
    isempty(models) && throw(ArgumentError("no models to plot"))
    ref = last(first(models))
    for (label, M) in models
        same_mesh(M, ref) || @warn "$label is on another mesh; its panel is not cell for cell with the first"
    end
    [(; label = first(p), M = last(p), a = _viewer_arrays(last(p); max_depth, pad_tol)) for p in models]
end

# the colour range shared by every panel, from the strongest contrasts in the core of all of them
function _panel_colorrange(panels, density_range)
    pooled = reduce(vcat, (vec(p.a.R_core) for p in panels))
    compute_colorrange(pooled;
                       density_range = density_range === nothing ? _body_density_range(pooled) : density_range)
end

# true where a reference model on the same mesh holds a body, over the given index ranges
_outline_mask(outline, M, dims...) = outline === nothing || !same_mesh(outline, M) ? nothing :
                                     replace(outline.A[dims...], NaN => 0.0) .!= 0

# positive depths on an axis drawn with negated z
_grav_depth_tickformat(values) = [string(round(Int, abs(v))) for v in values]

# a slider position as a label; map coordinates stay in metres, however large they are
_format_position(val::Real) = "$(round(val, digits = 1)) m"

# a bar length as a label, in km past a thousand metres
_format_length(len::Real) = len >= 1000 ? "$(round(len / 1000, digits = 2)) km" : "$(round(len, digits = 1)) m"

# the staircase of a section's topography as one polyline
_topography_points(values, he, ve) = (p = _topography_steps(values, he, ve); Point2f.(p[1], p[2]))

# the stations of a survey as open circles with a thin black edge
_viewer_stations!(ax, station_x, station_y) =
    isempty(station_x) ? nothing :
    scatter!(ax, station_y, station_x; color = :transparent, strokecolor = :black,
             strokewidth = _GRAV_STATION_STROKE, markersize = 5)

# axis keywords of a viewer panel, leaving the ticks automatic where the CRS gives none
_viewer_axis_kw(xlabel, ylabel, xticks, yticks) =
    merge((; xlabel, ylabel), xticks === nothing ? NamedTuple() : (; xticks),
          yticks === nothing ? NamedTuple() : (; yticks))

# where an export goes: the asked-for path, or a timestamped name so repeated exports accumulate
_export_path(given, dir, model_name, view, tag, full) =
    given === nothing ?
    joinpath(dir, "$(model_name)_$(view)_$(tag)_$(full ? "full" : "core")_$(Dates.format(now(), "yyyymmdd_HHMMSS")).png") :
    String(given)

# a round number near the asked-for length, for a scale bar
function _nice_scale_length(target::Real)
    target <= 0 && return 1.0
    p = 10.0^floor(log10(target))
    for m in (1.0, 2.0, 5.0)
        m * p >= target && return m * p
    end
    10p
end

#---------- depth slices ----------

"""
    depth_slice_viewer(models; model_name="model", cmap=Reverse(:Spectral), figsize=nothing,
                       with_padding=false, max_depth=nothing, pad_tol=0.2, density_range=nothing,
                       show_grid=false, outline=nothing, station_x=Float64[], station_y=Float64[],
                       axis_xlabel="Y East (m)", axis_ylabel="X North (m)", axis_xticks=nothing,
                       axis_yticks=nothing, export_dpi=3, export_figsize=nothing,
                       output_dir=pwd()) -> (fig, parts)

Slider-driven depth slices of one model or of several side by side, east horizontal and north vertical, one layer at a time.
"""
function depth_slice_viewer(models::AbstractVector{<:Pair{String, GravModel3D}};
                            model_name::AbstractString = "model",
                            cmap = _GRAV_COLORMAP,
                            figsize = nothing,
                            with_padding::Bool = false,
                            max_depth = nothing,
                            pad_tol::Real = 0.2,
                            density_range = nothing,
                            show_grid::Bool = false,
                            grid_color = :black,
                            grid_linewidth::Real = 0.5,
                            grid_alpha::Real = 0.15,
                            outline::Union{Nothing, GravModel3D} = nothing,
                            station_x::AbstractVector = Float64[],
                            station_y::AbstractVector = Float64[],
                            axis_xlabel::AbstractString = "Y East (m)",
                            axis_ylabel::AbstractString = "X North (m)",
                            axis_xticks = nothing,
                            axis_yticks = nothing,
                            export_dpi::Int = 3,
                            export_figsize = nothing,
                            output_dir::AbstractString = pwd())
    panels = _viewer_panels(models; max_depth, pad_tol)
    ncol   = length(panels)
    nz     = length(first(panels).a.z)
    edges  = first(panels).a.z_edges
    cmin, cmax = _panel_colorrange(panels, density_range)
    start  = cld(nz, 2)                                    # opens on the mesh centre, not an edge layer

    show_full = Observable(with_padding)
    layer     = Observable(start)
    state     = [(; p...,
                    xe       = Ref(with_padding ? p.a.xe_full : p.a.xe_core),
                    ye       = Ref(with_padding ? p.a.ye_full : p.a.ye_core),
                    values   = Ref(with_padding ? p.a.R_full  : p.a.R_core),
                    ranges   = Ref(with_padding ? (p.a.ix_full, p.a.iy_full) : (p.a.ix_core, p.a.iy_core)),
                    slice    = Observable(permutedims((with_padding ? p.a.R_full : p.a.R_core)[:, :, start])),
                    segments = Observable(Point2f[]),
                    ax       = Ref{Any}(nothing)) for p in panels]

    fig   = Figure(size = figsize === nothing ? (520 * ncol + 180, 950) : figsize, fontsize = 16)
    title = Observable(" ")
    Label(fig[0, 1:ncol+1], title, fontsize = 18, font = :bold)
    for (c, s) in enumerate(state)
        s.ax[] = _grav_axis(fig[1, c];
                            merge(_viewer_axis_kw(axis_xlabel, axis_ylabel, axis_xticks, axis_yticks),
                                  (; ylabel = c == 1 ? axis_ylabel : ""))...,
                            aspect = DataAspect(), title = ncol > 1 ? s.label : "")
    end
    Colorbar(fig[1, ncol+1]; colormap = cmap, colorrange = (cmin, cmax), label = _GRAV_DENSITY_LABEL,
             width = 16, height = Relative(0.8))

    reset_zoom!() = for s in state
        limits!(s.ax[], s.ye[][1], s.ye[][end], s.xe[][1], s.xe[][end])
    end

    # the heatmap and the outline follow their observables, so only a view change redraws the axes
    draw_axes! = function ()
        for s in state
            ax = s.ax[]
            empty!(ax)
            heatmap!(ax, s.ye[], s.xe[], s.slice; colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            linesegments!(ax, s.segments; color = :black, linewidth = _GRAV_OUTLINE_WIDTH)
            _viewer_stations!(ax, station_x, station_y)
            if show_grid
                for v in s.xe[]
                    lines!(ax, [s.ye[][1], s.ye[][end]], [v, v]; color = (grid_color, grid_alpha),
                           linewidth = grid_linewidth)
                end
                for v in s.ye[]
                    lines!(ax, [v, v], [s.xe[][1], s.xe[][end]]; color = (grid_color, grid_alpha),
                           linewidth = grid_linewidth)
                end
            end
        end
        reset_zoom!()
    end

    update_slice! = function (k)
        k       = clamp(k, 1, nz)
        layer[] = k
        for s in state
            s.slice[]    = permutedims(s.values[][:, :, k])
            mask         = _outline_mask(outline, s.M, s.ranges[][1], s.ranges[][2], s.a.kz[k])
            s.segments[] = mask === nothing ? Point2f[] : _cell_outline(permutedims(mask), s.ye[], s.xe[])
        end
        top, bot = edges[k], edges[k+1]
        title[]  = @sprintf("Depth layer %d / %d  |  %.0f to %.0f m  |  cell %.0f m", k, nz, top, bot, bot - top)
    end

    update_slice!(start)
    draw_axes!()

    bar      = fig[2, 1:ncol+1] = GridLayout()
    btn_prev = Button(bar[1, 1], label = "<< Prev")
    Label(bar[1, 2], "Depth layer:")
    slider   = Slider(bar[1, 3], range = 1:nz, startvalue = start, width = 400)
    count    = Observable("$start / $nz")
    Label(bar[1, 4], count)
    btn_next = Button(bar[1, 5], label = "Next >>")

    row        = fig[3, 1:ncol+1] = GridLayout()
    btn_first  = Button(row[1, 1], label = "|<< First")
    btn_last   = Button(row[1, 2], label = "Last >>|")
    btn_label  = Observable(with_padding ? "Show core" : "Show full")
    btn_toggle = Button(row[1, 3], label = btn_label)
    btn_reset  = Button(row[1, 4], label = "Reset zoom")
    btn_export = Button(row[1, 5], label = "Export figure")
    view_info  = Observable(with_padding ? _GRAV_FULL_INFO : _GRAV_CORE_INFO)
    Label(fig[4, 1:ncol+1], view_info, fontsize = 12, color = :gray35)

    on(slider.value) do k
        count[] = "$k / $nz"
        update_slice!(k)
    end
    on(btn_prev.clicks)  do _; set_close_to!(slider, max(1, slider.value[] - 1)) end
    on(btn_next.clicks)  do _; set_close_to!(slider, min(nz, slider.value[] + 1)) end
    on(btn_first.clicks) do _; set_close_to!(slider, 1) end
    on(btn_last.clicks)  do _; set_close_to!(slider, nz) end
    on(btn_reset.clicks) do _; reset_zoom!() end
    on(btn_toggle.clicks) do _
        show_full[] = !show_full[]
        for s in state
            s.xe[]     = show_full[] ? s.a.xe_full : s.a.xe_core
            s.ye[]     = show_full[] ? s.a.ye_full : s.a.ye_core
            s.values[] = show_full[] ? s.a.R_full  : s.a.R_core
            s.ranges[] = show_full[] ? (s.a.ix_full, s.a.iy_full) : (s.a.ix_core, s.a.iy_core)
        end
        btn_label[] = show_full[] ? "Show core" : "Show full"
        view_info[] = show_full[] ? _GRAV_FULL_INFO : _GRAV_CORE_INFO
        update_slice!(layer[])
        draw_axes!()
    end

    export_figure = function (output_path = nothing)
        f = Figure(size = export_figsize === nothing ? (520 * ncol + 180, 620) : export_figsize, fontsize = 18)
        Label(f[1, 1:ncol+1], title[], fontsize = 20, font = :bold)
        for (c, s) in enumerate(state)
            exax = _grav_axis(f[2, c];
                              merge(_viewer_axis_kw(axis_xlabel, axis_ylabel, axis_xticks, axis_yticks),
                                    (; ylabel = c == 1 ? axis_ylabel : ""))...,
                              aspect = DataAspect(), title = ncol > 1 ? s.label : "")
            heatmap!(exax, s.ye[], s.xe[], s.slice[]; colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            isempty(s.segments[]) || linesegments!(exax, s.segments[]; color = :black,
                                                   linewidth = _GRAV_OUTLINE_WIDTH)
            _viewer_stations!(exax, station_x, station_y)
            limits!(exax, s.ye[][1], s.ye[][end], s.xe[][1], s.xe[][end])
        end
        Colorbar(f[2, ncol+1]; colormap = cmap, colorrange = (cmin, cmax), label = _GRAV_DENSITY_LABEL,
                 width = 16, height = Relative(0.8))
        path = _export_path(output_path, output_dir, model_name, "XY", "layer$(layer[])", show_full[])
        mkpath(dirname(abspath(path)))
        save(path, f; px_per_unit = export_dpi, backend = CairoMakie)
        println("Figure exported: $path")
        path
    end
    on(btn_export.clicks) do _; export_figure() end

    fig, (; axes = [s.ax[] for s in state], ax = first(state).ax[], slider, layer,
            slice = first(state).slice, slices = [s.slice for s in state],
            segments = first(state).segments, show_full, z = first(panels).a.z, z_edges = edges,
            colorrange = (; cmin, cmax), reset_zoom!, export_figure)
end

depth_slice_viewer(M::GravModel3D; model_name::AbstractString = "model", kwargs...) =
    depth_slice_viewer([String(model_name) => M]; model_name, kwargs...)

#---------- sections ----------

# one section viewer, for a north-south (:xz) or an east-west (:yz) cut
function _section_viewer(models::AbstractVector{<:Pair{String, GravModel3D}}, view::Symbol;
                         model_name::AbstractString = "model",
                         cmap = _GRAV_COLORMAP,
                         figsize = nothing,
                         with_padding::Bool = false,
                         max_depth = nothing,
                         pad_tol::Real = 0.2,
                         density_range = nothing,
                         show_grid::Bool = false,
                         grid_color = :black,
                         grid_linewidth::Real = 0.5,
                         grid_alpha::Real = 0.15,
                         outline::Union{Nothing, GravModel3D} = nothing,
                         show_topography::Bool = true,
                         horiz_label::AbstractString = view == :xz ? "X North (m)" : "Y East (m)",
                         horiz_ticks = nothing,
                         slider_label::AbstractString = view == :xz ? "East" : "North",
                         export_dpi::Int = 3,
                         export_figsize = nothing,
                         output_dir::AbstractString = pwd())
    view in (:xz, :yz) || throw(ArgumentError("view must be :xz or :yz"))
    panels = _viewer_panels(models; max_depth, pad_tol)
    ncol   = length(panels)
    neg    = -first(panels).a.z_edges                      # depth drawn downward from the mesh top
    cmin, cmax = _panel_colorrange(panels, density_range)

    # the slider runs across the section: east for :xz, north for :yz
    across_of(a, full) = view == :xz ? (full ? a.y_full : a.y_core) : (full ? a.x_full : a.x_core)
    along_of(a, full)  = view == :xz ? (full ? a.xe_full : a.xe_core) : (full ? a.ye_full : a.ye_core)
    cut(R, k)          = view == :xz ? R[:, k, :] : R[k, :, :]

    start = cld(length(across_of(first(panels).a, with_padding)), 2)   # opens on the mesh centre, not an edge

    show_full = Observable(with_padding)
    index     = Observable(start)
    state     = [(; p...,
                    he       = Ref(along_of(p.a, with_padding)),
                    values   = Ref(with_padding ? p.a.R_full : p.a.R_core),
                    across   = Ref(across_of(p.a, with_padding)),
                    ranges   = Ref(with_padding ? (p.a.ix_full, p.a.iy_full) : (p.a.ix_core, p.a.iy_core)),
                    slice    = Observable(cut(with_padding ? p.a.R_full : p.a.R_core, start)),
                    segments = Observable(Point2f[]),
                    topo     = Observable(Point2f[]),
                    ax       = Ref{Any}(nothing)) for p in panels]
    nacross() = length(first(state).across[])

    fig   = Figure(size = figsize === nothing ? (560 * ncol + 180, 950) : figsize, fontsize = 16)
    title = Observable(" ")
    Label(fig[0, 1:ncol+1], title, fontsize = 18, font = :bold)
    for (c, s) in enumerate(state)
        s.ax[] = _grav_axis(fig[1, c];
                            merge(_viewer_axis_kw(horiz_label, "Depth (m)", horiz_ticks, nothing),
                                  (; ylabel = c == 1 ? "Depth (m)" : ""))...,
                            ytickformat = _grav_depth_tickformat, title = ncol > 1 ? s.label : "")
    end
    Colorbar(fig[1, ncol+1]; colormap = cmap, colorrange = (cmin, cmax), label = _GRAV_DENSITY_LABEL,
             width = 16, height = Relative(0.8))

    reset_zoom!() = for s in state
        limits!(s.ax[], s.he[][1], s.he[][end], neg[end], neg[1])
    end

    draw_axes! = function ()
        for s in state
            ax = s.ax[]
            empty!(ax)
            heatmap!(ax, s.he[], neg, s.slice; colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            linesegments!(ax, s.segments; color = :black, linewidth = _GRAV_OUTLINE_WIDTH)
            show_topography && lines!(ax, s.topo; color = :black, linewidth = _GRAV_TOPO_WIDTH)
            if show_grid
                for v in neg
                    lines!(ax, [s.he[][1], s.he[][end]], [v, v]; color = (grid_color, grid_alpha),
                           linewidth = grid_linewidth)
                end
                for v in s.he[]
                    lines!(ax, [v, v], [neg[end], neg[1]]; color = (grid_color, grid_alpha),
                           linewidth = grid_linewidth)
                end
            end
        end
        reset_zoom!()
    end

    update_slice! = function (k)
        k       = clamp(k, 1, nacross())
        index[] = k
        for s in state
            s.slice[] = cut(s.values[], k)
            mask      = view == :xz ? _outline_mask(outline, s.M, s.ranges[][1], s.ranges[][2][k], s.a.kz) :
                                      _outline_mask(outline, s.M, s.ranges[][1][k], s.ranges[][2], s.a.kz)
            s.segments[] = mask === nothing ? Point2f[] : _cell_outline(mask, s.he[], neg)
            s.topo[]     = show_topography ? _topography_points(s.slice[], s.he[], neg) : Point2f[]
        end
        title[] = @sprintf("%s section %d / %d  |  %s = %s", slider_label, k, nacross(), slider_label,
                           _format_position(first(state).across[][k]))
    end

    update_slice!(start)
    draw_axes!()

    bar      = fig[2, 1:ncol+1] = GridLayout()
    btn_prev = Button(bar[1, 1], label = "<< Prev")
    Label(bar[1, 2], "$(slider_label) slice:")
    slider   = Slider(bar[1, 3], range = 1:nacross(), startvalue = start, width = 400)
    count    = Observable("$start / $(nacross())")
    Label(bar[1, 4], count)
    btn_next = Button(bar[1, 5], label = "Next >>")

    row        = fig[3, 1:ncol+1] = GridLayout()
    btn_first  = Button(row[1, 1], label = "|<< First")
    btn_last   = Button(row[1, 2], label = "Last >>|")
    btn_label  = Observable(with_padding ? "Show core" : "Show full")
    btn_toggle = Button(row[1, 3], label = btn_label)
    btn_reset  = Button(row[1, 4], label = "Reset zoom")
    btn_export = Button(row[1, 5], label = "Export figure")
    view_info  = Observable(with_padding ? _GRAV_FULL_INFO : _GRAV_CORE_INFO)
    Label(fig[4, 1:ncol+1], view_info, fontsize = 12, color = :gray35)

    on(slider.value) do k
        count[] = "$k / $(nacross())"
        update_slice!(k)
    end
    on(btn_prev.clicks)  do _; set_close_to!(slider, max(1, slider.value[] - 1)) end
    on(btn_next.clicks)  do _; set_close_to!(slider, min(nacross(), slider.value[] + 1)) end
    on(btn_first.clicks) do _; set_close_to!(slider, 1) end
    on(btn_last.clicks)  do _; set_close_to!(slider, nacross()) end
    on(btn_reset.clicks) do _; reset_zoom!() end
    on(btn_toggle.clicks) do _
        show_full[] = !show_full[]
        for s in state
            s.he[]     = along_of(s.a, show_full[])
            s.values[] = show_full[] ? s.a.R_full : s.a.R_core
            s.across[] = across_of(s.a, show_full[])
            s.ranges[] = show_full[] ? (s.a.ix_full, s.a.iy_full) : (s.a.ix_core, s.a.iy_core)
        end
        btn_label[]    = show_full[] ? "Show core" : "Show full"
        view_info[]    = show_full[] ? _GRAV_FULL_INFO : _GRAV_CORE_INFO
        slider.range[] = 1:nacross()
        update_slice!(clamp(index[], 1, nacross()))
        draw_axes!()
    end

    export_figure = function (output_path = nothing)
        f = Figure(size = export_figsize === nothing ? (560 * ncol + 180, 560) : export_figsize, fontsize = 18)
        Label(f[1, 1:ncol+1], title[], fontsize = 20, font = :bold)
        for (c, s) in enumerate(state)
            exax = _grav_axis(f[2, c];
                              merge(_viewer_axis_kw(horiz_label, "Depth (m)", horiz_ticks, nothing),
                                    (; ylabel = c == 1 ? "Depth (m)" : ""))...,
                              ytickformat = _grav_depth_tickformat, title = ncol > 1 ? s.label : "")
            heatmap!(exax, s.he[], neg, s.slice[]; colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            isempty(s.segments[]) || linesegments!(exax, s.segments[]; color = :black,
                                                   linewidth = _GRAV_OUTLINE_WIDTH)
            isempty(s.topo[]) || lines!(exax, s.topo[]; color = :black, linewidth = _GRAV_TOPO_WIDTH)
            limits!(exax, s.he[][1], s.he[][end], neg[end], neg[1])
        end
        Colorbar(f[2, ncol+1]; colormap = cmap, colorrange = (cmin, cmax), label = _GRAV_DENSITY_LABEL,
                 width = 16, height = Relative(0.8))
        path = _export_path(output_path, output_dir, model_name, uppercase(string(view)),
                            "slice$(index[])", show_full[])
        mkpath(dirname(abspath(path)))
        save(path, f; px_per_unit = export_dpi, backend = CairoMakie)
        println("Figure exported: $path")
        path
    end
    on(btn_export.clicks) do _; export_figure() end

    fig, (; axes = [s.ax[] for s in state], ax = first(state).ax[], slider, index,
            slice = first(state).slice, slices = [s.slice for s in state],
            segments = first(state).segments, topo = first(state).topo, show_full,
            z = first(panels).a.z, z_edges = first(panels).a.z_edges,
            colorrange = (; cmin, cmax), reset_zoom!, export_figure)
end

"""
    xz_slice_viewer(models; kwargs...) -> (fig, parts)

Slider-driven north-south sections of one or several density models, the slider stepping east through the mesh.
"""
xz_slice_viewer(models::AbstractVector{<:Pair{String, GravModel3D}}; kwargs...) =
    _section_viewer(models, :xz; kwargs...)

xz_slice_viewer(M::GravModel3D; model_name::AbstractString = "model", kwargs...) =
    _section_viewer([String(model_name) => M], :xz; model_name, kwargs...)

"""
    yz_slice_viewer(models; kwargs...) -> (fig, parts)

Slider-driven east-west sections of one or several density models, the slider stepping north through the mesh.
"""
yz_slice_viewer(models::AbstractVector{<:Pair{String, GravModel3D}}; kwargs...) =
    _section_viewer(models, :yz; kwargs...)

yz_slice_viewer(M::GravModel3D; model_name::AbstractString = "model", kwargs...) =
    _section_viewer([String(model_name) => M], :yz; model_name, kwargs...)

#---------- three planes together ----------

"""
    draw_north_and_scale!(scene; xv, yv, z_fixed, show_north=true, show_scale=true, color=:black,
                          line_width=2.0)

A north arrow and a scale bar on a 3D scene, drawn on the plane `z_fixed`.
"""
function draw_north_and_scale!(scene; xv, yv, z_fixed::Real, show_north::Bool = true, show_scale::Bool = true,
                               color = :black, line_width::Real = 2.0)
    (isempty(xv) || isempty(yv)) && return nothing
    xmin, xmax = extrema(xv)
    ymin, ymax = extrema(yv)
    dx, dy = xmax - xmin, ymax - ymin
    (dx > 0 && dy > 0) || return nothing
    z = Float64(z_fixed)
    if show_north                                          # x is north, so the arrow points along +x
        x0, y0 = xmin + 0.78 * dx, ymin + 0.10 * dy
        len    = 0.09 * min(dx, dy)
        head   = 0.3 * len
        lines!(scene, [x0, x0 + len], [y0, y0], [z, z]; color, linewidth = line_width, overdraw = true)
        lines!(scene, [x0 + len - head, x0 + len], [y0 - head, y0], [z, z]; color, linewidth = line_width,
               overdraw = true)
        lines!(scene, [x0 + len - head, x0 + len], [y0 + head, y0], [z, z]; color, linewidth = line_width,
               overdraw = true)
        text!(scene, [x0 + len + 0.03 * dx], [y0], [z]; text = ["N"], color, fontsize = 16, overdraw = true)
    end
    if show_scale
        len  = _nice_scale_length(0.25 * dy)
        x0   = xmin + 0.08 * dx
        y0   = ymin + 0.10 * dy
        tick = 0.012 * dx
        lines!(scene, [x0, x0], [y0, y0 + len], [z, z]; color, linewidth = line_width, overdraw = true)
        lines!(scene, [x0 - tick, x0 + tick], [y0, y0], [z, z]; color, linewidth = line_width, overdraw = true)
        lines!(scene, [x0 - tick, x0 + tick], [y0 + len, y0 + len], [z, z]; color, linewidth = line_width,
               overdraw = true)
        text!(scene, [x0 + 0.02 * dx], [y0 + len / 2], [z]; text = [_format_length(len)], color, fontsize = 12,
              overdraw = true)
    end
    nothing
end

"""
    grav_3d_viewer(models; model_name="model", cmap=Reverse(:Spectral), figsize=nothing,
                   with_padding=false, max_depth=nothing, pad_tol=0.2, density_range=nothing,
                   station_x=Float64[], station_y=Float64[], station_z=Float64[],
                   show_north_arrow=true, show_scale_bar=true, export_dpi=3, output_dir=pwd()) -> (fig, parts)

The three orthogonal planes of one or several density models, a 3D scene each, on shared sliders and visibility toggles.
"""
function grav_3d_viewer(models::AbstractVector{<:Pair{String, GravModel3D}};
                        model_name::AbstractString = "model",
                        cmap = _GRAV_COLORMAP,
                        figsize = nothing,
                        with_padding::Bool = false,
                        max_depth = nothing,
                        pad_tol::Real = 0.2,
                        density_range = nothing,
                        station_x::AbstractVector = Float64[],
                        station_y::AbstractVector = Float64[],
                        station_z::AbstractVector = Float64[],
                        show_north_arrow::Bool = true,
                        show_scale_bar::Bool = true,
                        annotation_color = :black,
                        annotation_line_width::Real = 2.0,
                        default_view_direction = Vec3f(-1.05f0, -0.80f0, 0.72f0),
                        default_view_scale = 1.12f0,
                        export_dpi::Int = 3,
                        output_dir::AbstractString = pwd())
    panels = _viewer_panels(models; max_depth, pad_tol)
    ncol   = length(panels)
    cmin, cmax = _panel_colorrange(panels, density_range)

    show_full = Observable(with_padding)
    state     = [(; p...,
                    x      = Ref(with_padding ? p.a.x_full : p.a.x_core),
                    y      = Ref(with_padding ? p.a.y_full : p.a.y_core),
                    z      = Ref(-p.a.z),
                    R      = Ref(with_padding ? p.a.R_full : p.a.R_core),
                    scene  = Ref{Any}(nothing),
                    plt    = Ref{Any}(nothing),
                    planes = Ref{NamedTuple}((yz = nothing, xz = nothing, xy = nothing))) for p in panels]
    head = first(state)

    fig = Figure(size = figsize === nothing ? (940 * ncol + 100, 920) : figsize, fontsize = 14)
    for (c, s) in enumerate(state)
        ncol > 1 && Label(fig[0, c], s.label; fontsize = 18, font = :bold)
        s.scene[] = LScene(fig[1, c], show_axis = false)
    end
    controls = fig[2, 1:ncol+1] = GridLayout()
    info     = Observable(" ")
    Label(fig[3, 1:ncol+1], info, fontsize = 11, color = :gray35, halign = :left, tellwidth = false)
    rowsize!(fig.layout, 1, Relative(0.84))
    Colorbar(fig[1, ncol+1]; colormap = cmap, colorrange = (cmin, cmax), label = _GRAV_DENSITY_LABEL, width = 16)

    # plots are drawn in order of their z translation, so the annotations are built low and lifted back up,
    # which puts them in place but drawn after the slice planes
    lift_by   = 2 * maximum(abs, head.a.z_edges) + 1.0
    overlays! = function (s)
        scene = s.scene[].scene
        n0    = length(scene.plots)
        draw_north_and_scale!(scene; xv = s.x[], yv = s.y[], z_fixed = -first(s.a.z_edges) - lift_by,
                              show_north = show_north_arrow, show_scale = show_scale_bar,
                              color = annotation_color, line_width = annotation_line_width)
        if !isempty(station_x)
            zs = length(station_z) == length(station_x) ? collect(float.(station_z)) :
                                                          fill(-first(s.a.z_edges), length(station_x))
            scatter!(scene, station_x, station_y, zs .- lift_by; color = annotation_color,
                     markersize = 7, overdraw = true)
        end
        for p in scene.plots[n0+1:end]
            translate!(p, 0, 0, lift_by)
        end
    end

    build! = function ()
        for s in state
            scene = s.scene[].scene
            for p in copy(scene.plots)
                delete!(scene, p)
            end
            s.plt[]    = volumeslices!(s.scene[], s.x[], s.y[], s.z[], s.R[]; colormap = cmap,
                                       colorrange = (cmin, cmax), nan_color = :transparent, bbox_visible = true)
            s.planes[] = (yz = s.plt[][:heatmap_yz][], xz = s.plt[][:heatmap_xz][], xy = s.plt[][:heatmap_xy][])
            overlays!(s)
        end
    end

    build!()
    for s in state
        cam3d!(s.scene[].scene, projectiontype = Makie.Perspective)
    end

    fit_camera! = function ()
        for s in state
            xmn, xmx = extrema(s.x[]); ymn, ymx = extrema(s.y[]); zmn, zmx = extrema(s.z[])
            centre = Vec3f((xmn + xmx) / 2, (ymn + ymx) / 2, (zmn + zmx) / 2)
            horiz  = Float32(max(xmx - xmn, ymx - ymn))
            zspan  = max(Float32(zmx - zmn), 0.45f0 * horiz)
            d      = default_view_direction ./ sqrt(sum(abs2, default_view_direction))
            eye    = centre + Vec3f(default_view_scale * d[1] * horiz, default_view_scale * d[2] * horiz,
                                    default_view_scale * d[3] * zspan)
            update_cam!(s.scene[].scene, eye, centre, Vec3f(0, 0, 1))
        end
    end
    fit_camera!()

    nx()  = length(head.x[])
    ny()  = length(head.y[])
    ndep() = length(head.z[])

    Label(controls[1, 1], "North:", halign = :right)
    btn_prev_yz = Button(controls[1, 2], label = "Prev")
    sl_yz       = Slider(controls[1, 3], range = 1:nx(), startvalue = cld(nx(), 2), width = 180)
    btn_next_yz = Button(controls[1, 4], label = "Next")
    tog_yz      = Toggle(controls[1, 5], active = true)
    lbl_yz      = Label(controls[1, 6], _format_position(head.x[][sl_yz.value[]]), color = :gray35)

    Label(controls[1, 7], "East:", halign = :right)
    btn_prev_xz = Button(controls[1, 8], label = "Prev")
    sl_xz       = Slider(controls[1, 9], range = 1:ny(), startvalue = cld(ny(), 2), width = 180)
    btn_next_xz = Button(controls[1, 10], label = "Next")
    tog_xz      = Toggle(controls[1, 11], active = true)
    lbl_xz      = Label(controls[1, 12], _format_position(head.y[][sl_xz.value[]]), color = :gray35)

    Label(controls[1, 13], "Depth:", halign = :right)
    btn_prev_xy = Button(controls[1, 14], label = "Prev")
    sl_xy       = Slider(controls[1, 15], range = 1:ndep(), startvalue = cld(ndep(), 2), width = 180)
    btn_next_xy = Button(controls[1, 16], label = "Next")
    tog_xy      = Toggle(controls[1, 17], active = true)
    lbl_xy      = Label(controls[1, 18], _format_position(-head.z[][sl_xy.value[]]), color = :gray35)

    btn_toggle = Button(controls[1, 19], label = with_padding ? "Show core" : "Show full")
    btn_reset  = Button(controls[1, 20], label = "Reset view")
    btn_export = Button(controls[1, 21], label = "Export")
    colgap!(controls, 6)

    update_info! = function ()
        x, y, z = head.x[], head.y[], head.z[]
        info[] = @sprintf("North [%.0f, %.0f] m   East [%.0f, %.0f] m   Depth [0, %.0f] m   |   slice north %s, east %s, depth %s",
                          minimum(x), maximum(x), minimum(y), maximum(y), maximum(-z),
                          _format_position(x[clamp(sl_yz.value[], 1, length(x))]),
                          _format_position(y[clamp(sl_xz.value[], 1, length(y))]),
                          _format_position(-z[clamp(sl_xy.value[], 1, length(z))]))
    end

    plane_visibility! = function ()
        on_yz, on_xz, on_xy = tog_yz.active[], tog_xz.active[], tog_xy.active[]
        if !(on_yz || on_xz || on_xy)
            tog_xy.active[] = true
            on_xy = true
        end
        for s in state
            s.planes[].yz === nothing || (s.planes[].yz.visible[] = on_yz)
            s.planes[].xz === nothing || (s.planes[].xz.visible[] = on_xz)
            s.planes[].xy === nothing || (s.planes[].xy.visible[] = on_xy)
        end
    end

    on(sl_yz.value) do v
        for s in state
            s.plt[] === nothing || s.plt[][:update_yz][](v)
        end
        lbl_yz.text[] = _format_position(head.x[][clamp(v, 1, nx())])
        update_info!()
    end
    on(sl_xz.value) do v
        for s in state
            s.plt[] === nothing || s.plt[][:update_xz][](v)
        end
        lbl_xz.text[] = _format_position(head.y[][clamp(v, 1, ny())])
        update_info!()
    end
    on(sl_xy.value) do v
        for s in state
            s.plt[] === nothing || s.plt[][:update_xy][](v)
        end
        lbl_xy.text[] = _format_position(-head.z[][clamp(v, 1, ndep())])
        update_info!()
    end
    on(btn_prev_yz.clicks) do _; set_close_to!(sl_yz, max(1, sl_yz.value[] - 1)) end
    on(btn_next_yz.clicks) do _; set_close_to!(sl_yz, min(nx(), sl_yz.value[] + 1)) end
    on(btn_prev_xz.clicks) do _; set_close_to!(sl_xz, max(1, sl_xz.value[] - 1)) end
    on(btn_next_xz.clicks) do _; set_close_to!(sl_xz, min(ny(), sl_xz.value[] + 1)) end
    on(btn_prev_xy.clicks) do _; set_close_to!(sl_xy, max(1, sl_xy.value[] - 1)) end
    on(btn_next_xy.clicks) do _; set_close_to!(sl_xy, min(ndep(), sl_xy.value[] + 1)) end
    on(tog_yz.active) do _; plane_visibility!() end
    on(tog_xz.active) do _; plane_visibility!() end
    on(tog_xy.active) do _; plane_visibility!() end
    on(btn_reset.clicks) do _; fit_camera!() end
    on(btn_toggle.clicks) do _
        show_full[] = !show_full[]
        for s in state
            s.x[] = show_full[] ? s.a.x_full : s.a.x_core
            s.y[] = show_full[] ? s.a.y_full : s.a.y_core
            s.R[] = show_full[] ? s.a.R_full : s.a.R_core
        end
        btn_toggle.label[] = show_full[] ? "Show core" : "Show full"
        build!()
        sl_yz.range[] = 1:nx()
        sl_xz.range[] = 1:ny()
        set_close_to!(sl_yz, cld(nx(), 2))
        set_close_to!(sl_xz, cld(ny(), 2))
        plane_visibility!()
        fit_camera!()
        update_info!()
    end
    update_info!()

    # the software backend cannot render a 3D scene, so the export is the same planes drawn flat,
    # one row per model
    export_figure = function (output_path = nothing)
        i, j, k = sl_yz.value[], sl_xz.value[], sl_xy.value[]
        neg = -head.a.z_edges
        f   = Figure(size = (1620, 520 * ncol), fontsize = 18)
        for (r, s) in enumerate(state)
            xe     = show_full[] ? s.a.xe_full : s.a.xe_core
            ye     = show_full[] ? s.a.ye_full : s.a.ye_core
            prefix = ncol > 1 ? "$(s.label): " : ""
            ax1 = _grav_axis(f[r, 1]; xlabel = "Y East (m)", ylabel = "X North (m)", aspect = DataAspect(),
                             title = prefix * @sprintf("Depth %.0f m", s.a.z[k]))
            heatmap!(ax1, ye, xe, permutedims(s.R[][:, :, k]); colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            ax2 = _grav_axis(f[r, 2]; xlabel = "X North (m)", ylabel = "Depth (m)",
                             ytickformat = _grav_depth_tickformat,
                             title = prefix * @sprintf("East %.0f m", s.y[][j]))
            heatmap!(ax2, xe, neg, s.R[][:, j, :]; colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            limits!(ax2, xe[1], xe[end], neg[end], neg[1])
            ax3 = _grav_axis(f[r, 3]; xlabel = "Y East (m)", ylabel = "Depth (m)",
                             ytickformat = _grav_depth_tickformat,
                             title = prefix * @sprintf("North %.0f m", s.x[][i]))
            heatmap!(ax3, ye, neg, s.R[][i, :, :]; colormap = cmap, colorrange = (cmin, cmax),
                     nan_color = :transparent)
            limits!(ax3, ye[1], ye[end], neg[end], neg[1])
        end
        Colorbar(f[1:ncol, 4]; colormap = cmap, colorrange = (cmin, cmax), label = _GRAV_DENSITY_LABEL,
                 width = 16, height = Relative(0.8))
        _grav_labels!(f, [(r, c) for r in 1:ncol for c in 1:3])
        path = _export_path(output_path, output_dir, model_name, "XYZ", "n$(i)_e$(j)_d$(k)", show_full[])
        mkpath(dirname(abspath(path)))
        save(path, f; px_per_unit = export_dpi, backend = CairoMakie)
        println("Figure exported: $path")
        path
    end
    on(btn_export.clicks) do _; export_figure() end

    fig, (; scenes = [s.scene[] for s in state], scene = head.scene[],
            sliders = (; yz = sl_yz, xz = sl_xz, xy = sl_xy),
            toggles = (; yz = tog_yz, xz = tog_xz, xy = tog_xy),
            planes = head.planes, show_full,
            data = (; x = head.x, y = head.y, z = head.z, R = head.R),
            colorrange = (; cmin, cmax), fit_camera!, export_figure)
end

grav_3d_viewer(M::GravModel3D; model_name::AbstractString = "model", kwargs...) =
    grav_3d_viewer([String(model_name) => M]; model_name, kwargs...)

#---------- file entries ----------

# the models of a file entry, which takes one path, a model, or a list of label => model pairs
_entry_model(m::GravModel3D) = m
_entry_model(m::AbstractString) = (isfile(m) || error("model file not found: $m"); ReadModel3D(m))

_entry_label(m::GravModel3D) = isempty(m.name) ? "model" : splitext(basename(m.name))[1]
_entry_label(m::AbstractString) = splitext(basename(String(m)))[1]

function _plot_model_files(model_file)
    list = model_file isa AbstractString || model_file isa GravModel3D ? [model_file] : collect(model_file)
    isempty(list) && error("model_file is required")
    out = Pair{String, GravModel3D}[]
    for item in list
        label = item isa Pair ? String(first(item)) : _entry_label(item)
        M     = _entry_model(item isa Pair ? last(item) : item)
        @printf("Model %s: %d x %d x %d cells, %d padding cells per lateral side\n", label, size(M.A)...,
                M.npad[1])
        push!(out, label => M)
    end
    out
end

# models, optional survey and the name the exports are built from
function _plot_model_inputs(model_file, data_file::AbstractString)
    models = _plot_model_files(model_file)
    d = if isempty(data_file)
        nothing
    else
        isfile(data_file) || error("data file not found: $data_file")
        ReadData3D(data_file)
    end
    d === nothing || @printf("Survey %s: %d stations, %s\n", basename(data_file), ndata(d), d.crs)
    models, d, length(models) > 1 ? "compare" : first(first(models))
end

# the axes of a model panel follow the survey's CRS, the same choice the data plots make
function _model_axes(models, d, view::Symbol)
    east, north = _model_centre(last(first(models)))
    _slice_axes(view, d === nothing ? "local" : d.crs, east, north)
end

# open the window and block until it is closed, the way the examples are run
function _show_viewer(fig, interactive::Bool)
    interactive || return fig
    Makie.update_state_before_display!(fig)
    screen = GLMakie.Screen(fig.scene)
    println("Close the figure window to exit")
    wait(screen)
    fig
end

"""
    PlotModelXY(model_file, data_file=""; outline_file="", with_padding=false, max_depth=nothing,
                density_range=nothing, colormap=Reverse(:Spectral), show_grid=false, pad_tol=0.2,
                viewer_figsize=nothing, export_dpi=3, export_figsize=nothing,
                output_dir=pwd(), interactive=true) -> (fig, parts)

Depth slices of a model file, with the stations of `data_file` and the bodies of `outline_file` drawn on every layer; `model_file` may also be a list of `label => model` pairs, drawn side by side on one slider and one colour scale.
"""
function PlotModelXY(model_file, data_file::AbstractString = "";
                     outline_file::AbstractString = "",
                     with_padding::Bool = false,
                     max_depth = nothing,
                     density_range = nothing,
                     colormap = _GRAV_COLORMAP,
                     show_grid::Bool = false,
                     pad_tol::Real = 0.2,
                     viewer_figsize = nothing,
                     export_dpi::Int = 3,
                     export_figsize = nothing,
                     output_dir::AbstractString = pwd(),
                     interactive::Bool = true)
    models, d, name = _plot_model_inputs(model_file, data_file)
    outline = isempty(outline_file) ? nothing : ReadModel3D(outline_file)
    ha, va  = _model_axes(models, d, :xy)
    fig, parts = depth_slice_viewer(models; model_name = name, cmap = colormap, figsize = viewer_figsize,
                                    with_padding, max_depth, pad_tol, density_range, show_grid, outline,
                                    station_x = d === nothing ? Float64[] : d.x,
                                    station_y = d === nothing ? Float64[] : d.y,
                                    axis_xlabel = ha.label, axis_ylabel = va.label,
                                    axis_xticks = ha.ticks, axis_yticks = va.ticks,
                                    export_dpi, export_figsize, output_dir)
    _show_viewer(fig, interactive)
    fig, parts
end

"""
    PlotModelXZ(model_file, data_file=""; outline_file="", show_topography=true, kwargs...) -> (fig, parts)

North-south sections of a model file, the slider stepping east through the mesh.
"""
function PlotModelXZ(model_file, data_file::AbstractString = "";
                     outline_file::AbstractString = "",
                     with_padding::Bool = false,
                     max_depth = nothing,
                     density_range = nothing,
                     colormap = _GRAV_COLORMAP,
                     show_grid::Bool = false,
                     pad_tol::Real = 0.2,
                     show_topography::Bool = true,
                     viewer_figsize = nothing,
                     export_dpi::Int = 3,
                     export_figsize = nothing,
                     output_dir::AbstractString = pwd(),
                     interactive::Bool = true)
    models, d, name = _plot_model_inputs(model_file, data_file)
    outline = isempty(outline_file) ? nothing : ReadModel3D(outline_file)
    ha, _   = _model_axes(models, d, :xz)
    fig, parts = _section_viewer(models, :xz; model_name = name, cmap = colormap, figsize = viewer_figsize,
                                 with_padding, max_depth, pad_tol, density_range, show_grid, outline,
                                 show_topography, horiz_label = ha.label, horiz_ticks = ha.ticks,
                                 export_dpi, export_figsize, output_dir)
    _show_viewer(fig, interactive)
    fig, parts
end

"""
    PlotModelYZ(model_file, data_file=""; outline_file="", show_topography=true, kwargs...) -> (fig, parts)

East-west sections of a model file, the slider stepping north through the mesh.
"""
function PlotModelYZ(model_file, data_file::AbstractString = "";
                     outline_file::AbstractString = "",
                     with_padding::Bool = false,
                     max_depth = nothing,
                     density_range = nothing,
                     colormap = _GRAV_COLORMAP,
                     show_grid::Bool = false,
                     pad_tol::Real = 0.2,
                     show_topography::Bool = true,
                     viewer_figsize = nothing,
                     export_dpi::Int = 3,
                     export_figsize = nothing,
                     output_dir::AbstractString = pwd(),
                     interactive::Bool = true)
    models, d, name = _plot_model_inputs(model_file, data_file)
    outline = isempty(outline_file) ? nothing : ReadModel3D(outline_file)
    ha, _   = _model_axes(models, d, :yz)
    fig, parts = _section_viewer(models, :yz; model_name = name, cmap = colormap, figsize = viewer_figsize,
                                 with_padding, max_depth, pad_tol, density_range, show_grid, outline,
                                 show_topography, horiz_label = ha.label, horiz_ticks = ha.ticks,
                                 export_dpi, export_figsize, output_dir)
    _show_viewer(fig, interactive)
    fig, parts
end

"""
    PlotModelXYZ(model_file, data_file=""; with_padding=false, max_depth=nothing, density_range=nothing,
                 colormap=Reverse(:Spectral), pad_tol=0.2, viewer_figsize=nothing,
                 show_north_arrow=true, show_scale_bar=true, export_dpi=3, output_dir=pwd(),
                 interactive=true) -> (fig, parts)

The three orthogonal planes of a model file in one 3D scene, with the stations of `data_file` at their elevations.
"""
function PlotModelXYZ(model_file, data_file::AbstractString = "";
                      with_padding::Bool = false,
                      max_depth = nothing,
                      density_range = nothing,
                      colormap = _GRAV_COLORMAP,
                      pad_tol::Real = 0.2,
                      viewer_figsize = nothing,
                      show_north_arrow::Bool = true,
                      show_scale_bar::Bool = true,
                      export_dpi::Int = 3,
                      output_dir::AbstractString = pwd(),
                      interactive::Bool = true)
    models, d, name = _plot_model_inputs(model_file, data_file)
    top = last(first(models)).top_elevation
    fig, parts = grav_3d_viewer(models; model_name = name, cmap = colormap, figsize = viewer_figsize,
                                with_padding, max_depth, pad_tol, density_range,
                                station_x = d === nothing ? Float64[] : d.x,
                                station_y = d === nothing ? Float64[] : d.y,
                                station_z = d === nothing ? Float64[] : d.elevation .- top,
                                show_north_arrow, show_scale_bar, export_dpi, output_dir)
    _show_viewer(fig, interactive)
    fig, parts
end
