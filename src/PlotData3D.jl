# Data plots: station maps, observed / predicted / residual panels, the data viewer and the GIS export.
# Author: @pankajkmishra
# x north on the vertical axis, y east on the horizontal, in the data frame (crs = "model"), or easting and
# northing in WGS84 or any EPSG through the data origin (as the model viewers). Observed and predicted
# share one colour range; residuals are normalised by the errors on a symmetric diverging scale.

"""
    plot_data(d; value=d.value, colormap=COLORMAP, title="", size=(800, 700)) -> Figure

Map of the stations of `d` coloured by `value` (default the data).
"""
function plot_data(d::Data; value::AbstractVector = d.value, colormap = COLORMAP, title::AbstractString = "",
                   size = (800, 700))
    length(value) == ndata(d) || throw(ArgumentError("one value per station"))
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = "Y (m)", ylabel = "X (m)", aspect = DataAspect(),
              title = isempty(title) ? _COMPONENT_NAME[d.component] : title)
    sc = scatter!(ax, d.y, d.x; color = value, colormap, markersize = 9)
    Colorbar(fig[1, 2], sc; label = "$(_COMPONENT_NAME[d.component]) ($(units(d.component)))")
    fig
end

"""
    plot_data_fit(d, predicted; colormap=COLORMAP, size=(1500, 520)) -> Figure

Observed, predicted and normalised residual (d − p)/σ maps side by side, with the rms in the title.
"""
function plot_data_fit(d::Data, predicted::AbstractVector; colormap = COLORMAP, size = (1500, 520))
    length(predicted) == ndata(d) || throw(ArgumentError("one predicted value per station"))
    r = (d.value .- predicted) ./ d.error
    lo, hi = extrema(vcat(d.value, predicted))
    rmax = max(maximum(abs, r), eps())
    u = "$(_COMPONENT_NAME[d.component]) ($(units(d.component)))"
    fig = Figure(; size)
    for (k, (v, t, cm, cr, lab)) in enumerate(((d.value, "Observed", colormap, (lo, hi), u),
                                               (predicted, "Predicted", colormap, (lo, hi), u),
                                               (r, @sprintf("Residual, rms %.3f", sqrt(sum(abs2, r) / length(r))), colormap,
                                                (-rmax, rmax), "(d − p) / σ")))
        ax = Axis(fig[1, 2k-1]; xlabel = "Y (m)", ylabel = k == 1 ? "X (m)" : "", aspect = DataAspect(), title = t)
        sc = scatter!(ax, d.y, d.x; color = v, colormap = cm, colorrange = cr, markersize = 8)
        Colorbar(fig[1, 2k], sc; label = lab)
    end
    fig
end

#---------- data viewer ----------

"""
    PlotDataMap(data_file, predicted_file=""; crs, shapefiles, ..., interactive)

Interactive map of a data file: observed values, and with a predicted file (e.g. data.pred of a
run) the predicted values and normalised residuals (d − p)/σ, switched by buttons, with shapefile
overlays (given, or added from a file chooser), PNG export and GIS export of the stations.

In:  path to a native data file, optionally a predicted file on the same stations, and the display
     options; see `examples/plot_data_map.jl` for the full set with defaults.
Out: the Figure and its parts, after the window closes (`interactive = false`: at once, no window).

`write_data_gis` writes the same shapefile headless. `crs = "model"` (the default) plots frame
coordinates, so data without latitude/longitude can be shown; any EPSG goes through the data origin.
"""
function PlotDataMap(data_file::AbstractString, predicted_file::AbstractString = "";
    crs::AbstractString = "model",
    shapefiles = [],
    shapefile_color = :grey30,
    shapefile_alpha::Real = 0.9,
    shapefile_line_width::Real = 1.2,
    shapefile_point_size::Real = 6,
    colormap = COLORMAP,
    value_range = nothing,
    residual_range = (-3.0, 3.0),
    site_markersize::Real = 9,
    map_pad::Real = 0.05,
    viewer_figsize = (1100, 950),
    export_dpi::Int = 3,
    export_figsize = (1100, 900),
    gis_output_dir::AbstractString = "",
    interactive::Bool = true)

    isempty(data_file) && error("data_file is required")
    isfile(data_file) || error("Data file not found: $data_file")
    println("Loading data: $data_file")
    d = load_data_grav(data_file)
    p = nothing
    if !isempty(predicted_file)
        isfile(predicted_file) || error("Predicted data file not found: $predicted_file")
        println("Loading predicted data: $predicted_file")
        p = load_data_grav(predicted_file)
        p.x == d.x && p.y == d.y || error("$predicted_file holds other stations than $data_file")
    end
    name = _COMPONENT_NAME[d.component]
    unit = units(d.component)
    println("  $(ndata(d)) stations, $name ($unit)")

    crs_up = uppercase(strip(crs))
    sx, sy = _stations_in_crs(d, crs)
    xlabel, ylabel = crs_up == "MODEL" ? ("Y (m)", "X (m)") : crs_up == "EPSG:4326" ? ("Longitude (°)", "Latitude (°)") :
                     ("Easting (m)", "Northing (m)")
    aspect = crs_up == "EPSG:4326" ? AxisAspect(Float64(first(_distance_based_aspect(sy, sx)))) : DataAspect()
    loaded_shapefiles = Any[]
    if !isempty(shapefiles)
        println("\nShapefile overlay:")
        append!(loaded_shapefiles, prepare_shapefiles(shapefiles, crs))
    end
    shp_def(path) = (enabled = true, path = String(path), color = shapefile_color,
                     alpha = shapefile_alpha, point_size = shapefile_point_size,
                     line_width = shapefile_line_width)

    panels = Dict(:observed => (d.value, "Observed $name", "$name ($unit)"))
    if p !== nothing
        panels[:predicted] = (p.value, "Predicted $name", "$name ($unit)")
        panels[:residual] = ((d.value .- p.value) ./ d.error, "Normalised residual, rms " *
                             @sprintf("%.3f", sqrt(sum(abs2, (d.value .- p.value) ./ d.error) / ndata(d))), "(d − p) / σ")
    end
    vrange = compute_colorrange(p === nothing ? d.value : vcat(d.value, p.value); value_range)
    range_of(k) = k == :residual ? Tuple(Float64.(residual_range)) : vrange

    current = Observable(:observed)
    values = Observable(panels[:observed][1])
    crange = Observable(vrange)
    title = Observable(panels[:observed][2])
    cblabel = Observable(panels[:observed][3])

    fig = Figure(size = viewer_figsize)
    Label(fig[0, 1:2], title, fontsize = 18, font = :bold)
    ax = Axis(fig[1, 1]; xlabel, ylabel, aspect, xtickformat = _plain_tickformat, ytickformat = _plain_tickformat)
    sc = scatter!(ax, sx, sy; color = values, colormap, colorrange = crange, markersize = site_markersize)
    shp_plots = Any[]
    function redraw_overlays!()
        foreach(p -> delete!(ax, p), shp_plots); empty!(shp_plots)
        n0 = length(ax.scene.plots)
        _draw_all_shapefiles!(ax, loaded_shapefiles)
        append!(shp_plots, ax.scene.plots[n0+1:end])
    end
    redraw_overlays!()
    pad(v) = (lo = minimum(v); hi = maximum(v); w = max(hi - lo, eps()); (lo - map_pad * w, hi + map_pad * w))
    limits!(ax, pad(sx)..., pad(sy)...)
    Colorbar(fig[1, 2], sc; label = cblabel)

    buttons = fig[2, 1:2] = GridLayout()
    show = Dict(:observed => Button(buttons[1, 1], label = "Observed"))
    p === nothing || (show[:predicted] = Button(buttons[1, 2], label = "Predicted");
                      show[:residual] = Button(buttons[1, 3], label = "Residual"))
    btn_export = Button(buttons[1, 4], label = "Export Figure")
    btn_gis = Button(buttons[1, 5], label = "Export GIS")
    btn_browse = Button(buttons[1, 6], label = "Add Shapefile")
    btn_clear = Button(buttons[1, 7], label = "Clear Shapefiles")
    info = Observable("CRS: $crs | $(ndata(d)) stations")
    Label(fig[3, 1:2], info, fontsize = 12, color = :blue)

    function select!(k)
        current[] = k
        values[], title[], cblabel[] = panels[k]
        crange[] = range_of(k)
    end
    for (k, b) in show
        on(_ -> select!(k), b.clicks)
    end

    stem = splitext(basename(data_file))[1]
    function export_figure()
        k = current[]
        v, t, lab = panels[k]
        ef = Figure(size = export_figsize, fontsize = 16)
        Label(ef[1, 1:2], t, fontsize = 20, font = :bold)
        eax = Axis(ef[2, 1]; xlabel, ylabel, aspect, xtickformat = _plain_tickformat, ytickformat = _plain_tickformat)
        esc = scatter!(eax, sx, sy; color = v, colormap, colorrange = range_of(k), markersize = site_markersize)
        _draw_all_shapefiles!(eax, loaded_shapefiles)
        limits!(eax, pad(sx)..., pad(sy)...)
        Colorbar(ef[2, 2], esc; label = lab)
        filename = joinpath(dirname(abspath(data_file)), "$(stem)_$(k)_$(Dates.format(now(), "yyyymmdd_HHMMSS")).png")
        save(filename, ef; px_per_unit = export_dpi, backend = CairoMakie)
        println("Figure exported: $filename")
        info[] = "Exported $filename"
        filename
    end
    on(_ -> export_figure(), btn_export.clicks)
    on(btn_gis.clicks) do _
        try
            out = write_data_gis(data_file, predicted_file; crs = crs_up == "MODEL" ? "model" : crs, output_dir = gis_output_dir)
            info[] = "GIS exported: $out"
        catch err
            info[] = "GIS export failed: $(err)"
            println("GIS export error: ", err)
        end
    end

    on(btn_browse.clicks) do _
        path = String(strip(pick_shapefile()))
        if isempty(path)
            info[] = "No file chosen (or no file chooser on this system)."
            return
        elseif !isfile(path)
            info[] = "Not a file: $path"
            return
        end
        try
            added = prepare_shapefiles([shp_def(path)], crs)
            if isempty(added)
                info[] = "Nothing loaded from $(basename(path))"
            else
                append!(loaded_shapefiles, added)
                redraw_overlays!()
                info[] = "Added $(basename(path)) ($(length(added[1].geoms)) features)"
            end
        catch e
            showerror(stderr, e, catch_backtrace()); println(stderr)
            info[] = "Import failed: $(sprint(showerror, e))"
        end
    end
    on(btn_clear.clicks) do _
        empty!(loaded_shapefiles)
        redraw_overlays!()
        info[] = "Shapefiles cleared"
    end

    println("\nViewer ready: Observed$(p === nothing ? "" : " / Predicted / Residual") buttons, Export Figure, Export GIS, Add/Clear Shapefiles")
    isempty(loaded_shapefiles) || println("  - $(length(loaded_shapefiles)) shapefile(s) overlaid")
    interactive && _open_window(fig)
    return fig, (ax = ax, scatter = sc, select! = select!, export_figure = export_figure, x = sx, y = sy,
                 shapefiles = loaded_shapefiles, add_shapefile! = path -> (append!(loaded_shapefiles, prepare_shapefiles([shp_def(path)], crs)); redraw_overlays!()),
                 clear_shapefiles! = () -> notify(btn_clear.clicks))
end

"""
    write_data_gis(data_file, predicted_file=""; crs="EPSG:4326", output_dir="") -> path

Stations of a data file as a point shapefile in `crs` (WGS84, any EPSG, or "model" for frame
coordinates without a .prj), with the fields code, elev_m, value, error and, with a predicted file,
pred and resid ((d − p)/σ). Default directory: <data>-GIS next to the data file.
"""
function write_data_gis(data_file::AbstractString, predicted_file::AbstractString = "";
                        crs::AbstractString = "EPSG:4326", output_dir::AbstractString = "")
    d = load_data_grav(data_file)
    sx, sy = _stations_in_crs(d, crs)
    dir = mkpath(isempty(output_dir) ? joinpath(dirname(abspath(data_file)), splitext(basename(data_file))[1] * "-GIS") : output_dir)
    feats = (code = d.code, elev_m = -d.z, value = d.value, error = d.error)
    if !isempty(predicted_file)
        p = load_data_grav(predicted_file)
        p.x == d.x && p.y == d.y || error("$predicted_file holds other stations than $data_file")
        feats = merge(feats, (pred = p.value, resid = (d.value .- p.value) ./ d.error))
    end
    path = joinpath(dir, splitext(basename(data_file))[1] * "_stations.shp")
    _write_shapefile_with_sidecars(path, [Shapefile.Point(x, y) for (x, y) in zip(sx, sy)], feats,
                                   _resolve_prj_wkt_for_crs(crs))
    println("Stations exported: $path ($(ndata(d)) points, CRS $crs)")
    path
end
