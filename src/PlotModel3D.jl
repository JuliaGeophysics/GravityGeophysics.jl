function _draw_box!(ax, model::GravityModel)
    xc, yc, zc = get_cell_centers(model.mesh)
    x_edges, y_edges, z_edges = get_cell_edges(model.mesh)
    xmin, xmax = first(x_edges), last(x_edges)
    ymin, ymax = first(y_edges), last(y_edges)
    zmin, zmax = -last(z_edges), -first(z_edges)
    GLMakie.lines!(ax, [xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], fill(zmax, 5); color = :gray65, linewidth = 1.2)
    GLMakie.lines!(ax, [xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], fill(zmin, 5); color = :gray65, linewidth = 1.2)
    for (x0, y0) in ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax))
        GLMakie.lines!(ax, [x0, x0], [y0, y0], [zmin, zmax]; color = :gray65, linewidth = 1.2)
    end
end

function _fit_camera!(scene, model::GravityModel)
    x_edges, y_edges, z_edges = get_cell_edges(model.mesh)
    xmin, xmax = first(x_edges), last(x_edges)
    ymin, ymax = first(y_edges), last(y_edges)
    zmin, zmax = -last(z_edges), -first(z_edges)
    center = GLMakie.Vec3f((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
    span = Float32(max(xmax - xmin, ymax - ymin))
    zspan = Float32(zmax - zmin)
    eye = center + GLMakie.Vec3f(-1.05f0 * span, -0.82f0 * span, max(0.65f0 * span, 1.25f0 * zspan))
    GLMakie.update_cam!(scene, eye, center, GLMakie.Vec3f(0, 0, 1))
end

function _display_figure(fig; fullscreen::Bool = false)
    screen = fullscreen ? GLMakie.Screen(; fullscreen = true, float = false, focus_on_show = true) : GLMakie.Screen(; focus_on_show = true)
    display(screen, fig)
    return screen
end

function launch_model_viewer(model::GravityModel;
    data::Union{Nothing, GravityData} = nothing,
    shapefile_paths::Vector{String} = String[],
    open_fullscreen::Bool = false,
    colormap = :balance,
    export_png_scale::Int = 3,
    block::Bool = !isinteractive())

    xc, yc, zc = get_cell_centers(model.mesh)
    x_edges = edges_from_centers(xc)
    y_edges = edges_from_centers(yc)
    cmin, cmax = compute_colorrange(model.density)

    fig = GLMakie.Figure(size = (1600, 950))
    ax3 = GLMakie.LScene(fig[1, 1], show_axis = false)
    controls = fig[2, 1] = GLMakie.GridLayout()
    selector_ax = GLMakie.Axis(fig[3, 1], title = "Map slice", aspect = GLMakie.DataAspect(), xlabel = "X (m)", ylabel = "Y (m)")
        selector_ax = GLMakie.Axis(fig[3, 1], title = "Map slice", aspect = GLMakie.DataAspect(), xlabel = "Easting (m)", ylabel = "Northing (m)")
    colorbar = GLMakie.Colorbar(fig[1, 2], colormap = colormap, limits = (cmin, cmax), label = "Density contrast (kg/m^3)", width = 18)

    GLMakie.rowsize!(fig.layout, 1, GLMakie.Relative(0.62))
    GLMakie.rowsize!(fig.layout, 2, GLMakie.Auto(90))
    GLMakie.rowsize!(fig.layout, 3, GLMakie.Relative(0.38))

    depth_idx = GLMakie.Observable(clamp(div(model.mesh.nz, 2) + 1, 1, model.mesh.nz))
    show_slice = GLMakie.Observable(true)
    show_receivers = GLMakie.Observable(data !== nothing)
    show_isovolume = GLMakie.Observable(false)
    export_status = GLMakie.Observable("Viewer ready")

    GLMakie.Label(controls[1, 1], "Depth", halign = :right)
    btn_prev = GLMakie.Button(controls[1, 2], label = "Prev")
    slider = GLMakie.Slider(controls[1, 3], range = 1:model.mesh.nz, startvalue = depth_idx[], width = 260)
    btn_next = GLMakie.Button(controls[1, 4], label = "Next")
    depth_label = GLMakie.Label(controls[1, 5], "Depth: $(round(zc[depth_idx[]], digits = 0)) m")
    btn_slice = GLMakie.Button(controls[1, 6], label = "Hide Slice")
    btn_receivers = GLMakie.Button(controls[1, 7], label = data === nothing ? "Receivers N/A" : "Hide Receivers")
    btn_iso = GLMakie.Button(controls[1, 8], label = "Show Iso")
    btn_export = GLMakie.Button(controls[1, 9], label = "Export PNG")
    GLMakie.Label(controls[2, 1:9], export_status, halign = :left, color = :gray30)

    selector_slice = GLMakie.Observable(copy(model.density[:, :, depth_idx[]]))
    GLMakie.heatmap!(selector_ax, x_edges, y_edges, selector_slice; colormap = colormap, colorrange = (cmin, cmax))
    for shp in shapefile_paths
        plot_shapefile_on_axis!(selector_ax, shp; xlim = extrema(xc), ylim = extrema(yc))
    end

    dynamic_plots = Any[]
    function clear_dynamic!()
        for plot in reverse(dynamic_plots)
            try
                delete!(ax3.scene, plot)
            catch
                try
                    delete!(ax3, plot)
                catch
                end
            end
        end
        empty!(dynamic_plots)
    end

    function redraw!()
        clear_dynamic!()
        idx = clamp(depth_idx[], 1, model.mesh.nz)
        selector_slice[] = copy(model.density[:, :, idx])
        depth_label.text[] = "Depth: $(round(zc[idx], digits = 0)) m"
        _draw_box!(ax3, model)

        if show_slice[]
            push!(dynamic_plots, GLMakie.surface!(ax3, xc, yc, fill(-zc[idx], length(xc), length(yc)); color = model.density[:, :, idx], colormap = colormap, colorrange = (cmin, cmax), shading = GLMakie.NoShading, transparency = true, alpha = 0.92))
        end

        if show_isovolume[]
            triangles, selected = build_range_volume_triangles(xc, yc, zc, model.density, cmin + 0.65 * (cmax - cmin), cmax, minimum(zc), maximum(zc))
            if !isempty(triangles)
                vertices, faces = triangles_to_vertices_faces(triangles)
                push!(dynamic_plots, GLMakie.mesh!(ax3, vertices, faces; color = fill(Float32(cmax), length(vertices)), colormap = colormap, colorrange = (cmin, cmax), transparency = true, alpha = 0.35))
                export_status[] = "Iso cells shown: $selected"
            end
        end

        if data !== nothing && show_receivers[]
            push!(dynamic_plots, GLMakie.scatter!(ax3, data.receivers.x, data.receivers.y, -data.receivers.z; color = data.gz .* SI_TO_MGAL, markersize = 7, colormap = :viridis, colorrange = extrema(data.gz .* SI_TO_MGAL)))
        end

        for shp in shapefile_paths
            plot_shapefile_on_3d!(ax3, shp; z_fixed = 0.0, xlim = extrema(xc), ylim = extrema(yc))
        end
    end

    redraw!()
    _fit_camera!(ax3.scene, model)

    on(slider.value) do value
        depth_idx[] = clamp(round(Int, value), 1, model.mesh.nz)
        redraw!()
    end
    on(btn_prev.clicks) do _
        depth_idx[] = max(1, depth_idx[] - 1)
        slider.value[] = depth_idx[]
        redraw!()
    end
    on(btn_next.clicks) do _
        depth_idx[] = min(model.mesh.nz, depth_idx[] + 1)
        slider.value[] = depth_idx[]
        redraw!()
    end
    on(btn_slice.clicks) do _
        show_slice[] = !show_slice[]
        btn_slice.label[] = show_slice[] ? "Hide Slice" : "Show Slice"
        redraw!()
    end
    on(btn_receivers.clicks) do _
        if data === nothing
            return
        end
        show_receivers[] = !show_receivers[]
        btn_receivers.label[] = show_receivers[] ? "Hide Receivers" : "Show Receivers"
        redraw!()
    end
    on(btn_iso.clicks) do _
        show_isovolume[] = !show_isovolume[]
        btn_iso.label[] = show_isovolume[] ? "Hide Iso" : "Show Iso"
        redraw!()
    end
    on(btn_export.clicks) do _
        outdir = joinpath(@__DIR__, "..", "examples", "demo_output")
        mkpath(outdir)
        outpath = joinpath(outdir, "gravity_viewer_$(round(Int, time())).png")
        GLMakie.save(outpath, fig; px_per_unit = export_png_scale)
        export_status[] = "Saved viewer snapshot: $(basename(outpath))"
    end

    screen = _display_figure(fig; fullscreen = open_fullscreen)
    if block
        wait(screen)
    end
    return fig, screen
end