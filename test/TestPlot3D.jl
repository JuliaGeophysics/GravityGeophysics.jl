# Plots and viewers, headless (interactive = false): figures, CRS axes and GIS export.


# block problem written to files; stations get lat/lon around (60° N, 25° E) through a local TM frame
function plot_case()
    m, d, true_A = block_problem()
    mt = Model(m.dx, m.dy, m.dz, true_A; origin = m.origin)
    mt.A[:, :, 1] .= NaN                                        # an air layer
    mt = Model(mt.dx, mt.dy, mt.dz, mt.A; origin = mt.origin)
    f = Frame(60.0, 25.0)
    lat = 60.0 .+ d.x ./ 111_320; lon = 25.0 .+ d.y ./ (111_320 * cosd(60))
    x, y = to_frame(f, lat, lon)
    d.lat .= lat; d.lon .= lon; d.x .= x; d.y .= y; d.origin = [60.0, 25.0]
    op = GravForward(mt, d)
    dir = mktempdir()
    (dir = dir, m = mt, d = d,
     model = write_model(joinpath(dir, "model.den"), mt),
     data = write_data_grav(joinpath(dir, "data.dat"), d),
     pred = write_data_grav(joinpath(dir, "data.pred"), d; value = d.value .+ 0.5 .* d.error))
end

@testset "model viewers" begin
    c = plot_case()
    fig, parts = PlotModelXY(c.model, c.data; interactive = false)
    @test parts.ax.xlabel[] == "Y (m)" && parts.ax.ylabel[] == "X (m)"
    @test all(isnan, parts.R_full[:, :, 1])                              # air stays NaN
    @test (parts.colorrange.cmin, parts.colorrange.cmax) == compute_colorrange(parts.R_full)
    @test filesize(save_figure(joinpath(c.dir, "xy.png"), fig)) > 10_000
    _, p4326 = PlotModelXY(c.model, c.data; crs = "EPSG:4326", interactive = false)
    @test p4326.ax.xlabel[] == "Longitude (°)" && 24.9 < mean(p4326.y_full) < 25.1 && 59.9 < mean(p4326.x_full) < 60.1
    _, putm = PlotModelXY(c.model, c.data; crs = "EPSG:32635", value_range = (0.0, 300.0), interactive = false)
    @test putm.ax.xlabel[] == "Easting (m)" && putm.colorrange == (cmin = 0.0, cmax = 300.0)
    fig, parts = PlotModelXZ(c.model; max_depth = 200, interactive = false)
    @test parts.ax.ylabel[] == "Depth (m)" && filesize(save_figure(joinpath(c.dir, "xz.png"), fig)) > 10_000
    # depth positive down from 0 at the top of the topography: the air layer above it is negative
    @test parts.ax.yreversed[] && parts.z[1] < 0 < parts.z[2] && parts.z[2] ≈ c.m.dz[2] / 2
    fig, parts = PlotModelYZ(c.model; log10_scale = true, interactive = false)
    @test filesize(save_figure(joinpath(c.dir, "yz.png"), fig)) > 10_000
    fig, parts = PlotModelXYZ(c.model, c.data; max_depth = 300, interactive = false)
    @test size(parts.data.R[]) == (length(parts.data.x[]), length(parts.data.y[]), length(parts.data.z[]))
    # frame coordinates have no lat/lon: geographic axes are refused, model axes work
    m, d, _ = block_problem()
    nolatlon = write_data_grav(joinpath(c.dir, "frame.dat"), d)
    @test_throws ErrorException PlotModelXY(c.model, nolatlon; crs = "EPSG:4326", interactive = false)
end

@testset "GIS export of depth slices" begin
    c = plot_case()
    _, parts = PlotModelXY(c.model, c.data; crs = "EPSG:32635", interactive = false)
    lo, hi = parts.colorrange
    out = joinpath(c.dir, "gis")
    n = GG._export_all_depth_slices_gis(output_dir = out, model_name = "model", x_edges = parts.x_edges_core,
                                        y_edges = parts.y_edges_core, R = parts.R_core, z_edges = parts.z_edges,
                                        crs = "EPSG:32635", log10scale = false, vis_range = (lo, hi))
    shp = sort(filter(f -> endswith(f, ".shp"), readdir(out)))
    @test n == length(shp) == count(k -> any(isfinite, parts.R_core[:, :, k]), axes(parts.R_core, 3)) == size(parts.R_core, 3) - 1
    @test all(isfile, joinpath.(out, replace.(shp, ".shp" => ".prj")))
    @test any(startswith(f, "model_depth_0m_") for f in shp)                 # depths from the top of the topography
    tables = [GG.Shapefile.Table(joinpath(out, f)) for f in shp]
    v = reduce(vcat, [Float64.(t.value) for t in tables]); vis = reduce(vcat, [Float64.(t.value_vis) for t in tables])
    @test count(isfinite, parts.R_core) == length(v)                   # one polygon per earth cell
    @test vis == clamp.(v, lo, hi) && maximum(v) > 0
    @test !any(f -> endswith(f, ".qml") || endswith(f, ".lyrx"), readdir(out))   # no style files shipped
end

@testset "data viewer and plots" begin
    c = plot_case()
    fig, parts = PlotDataMap(c.data, c.pred; interactive = false)
    @test maximum(abs, parts.x .- c.d.y) < 1e-3 && maximum(abs, parts.y .- c.d.x) < 1e-3
    parts.select!(:residual)
    @test collect(parts.scatter.colorrange[]) ≈ [-3.0, 3.0]
    path = parts.export_figure()
    @test dirname(path) == c.dir && filesize(path) > 10_000
    _, pll = PlotDataMap(c.data; crs = "EPSG:4326", interactive = false)
    @test all(isapprox.(pll.y, c.d.lat))
    shp = write_data_gis(c.data, c.pred; output_dir = joinpath(c.dir, "gis"))
    t = GG.Shapefile.Table(shp)
    @test length(t) == ndata(c.d) && isapprox(Float64.(t.resid), fill(-0.5, ndata(c.d)); atol = 1e-3)
    @test isfile(replace(shp, ".shp" => ".prj"))
    fig, pm = PlotDataMap(c.data; crs = "model", interactive = false)
    model_shp = write_data_gis(c.data; crs = "model", output_dir = joinpath(c.dir, "gis-model"))
    pm.add_shapefile!(model_shp)                                   # as the viewer's Add Shapefile button
    @test length(pm.shapefiles) == 1
    pm.clear_shapefiles!()
    @test isempty(pm.shapefiles)
    m, d, true_A = block_problem()
    mt = Model(m.dx, m.dy, m.dz, true_A; origin = m.origin)
    op = GravForward(m, d)
    files = [save_figure(joinpath(c.dir, "d.png"), plot_data(d)),
             save_figure(joinpath(c.dir, "fit.png"), plot_data_fit(d, forward(op, model_vector(op, mt))))]
    @test all(f -> filesize(f) > 10_000, files)
    @test_throws ArgumentError plot_data(d; value = [1.0])
end

@testset "exports go next to the inputs" begin
    c = plot_case()
    @test dirname(write_data_gis(c.data)) == joinpath(c.dir, "data-GIS")
end
