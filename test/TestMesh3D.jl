# Mesh design: core, padding, layers, topography (function and GeoTIFF), station checks, MakeMesh3D and the window.

grid_data(; s = 40.0, n = 9, z = -30.0) = (x = range(-(n - 1) / 2 * s, (n - 1) / 2 * s; length = n);
    station_data([a for a in x for _ in x], [b for _ in x for b in x], fill(z, n^2)))

bundle(d; kw...) = GG.build_mesh_bundle(d, d.x, d.y; dx_core = 20.0, dy_core = 20.0,
    nx_core = GG.design_mesh(d.x, d.y; dx_core = 20.0, dy_core = 20.0, pad_factor = 1.4, z_first = 20.0, z_factor = 1.4, depth_mult = 0.5).nx_core,
    ny_core = GG.design_mesh(d.x, d.y; dx_core = 20.0, dy_core = 20.0, pad_factor = 1.4, z_first = 20.0, z_factor = 1.4, depth_mult = 0.5).ny_core,
    pad_factor = 1.4, z_first = 20.0, z_factor = 1.4, depth_mult = 0.5, kw...)

@testset "design_mesh" begin
    d = grid_data()
    @test GG.nearest_neighbour_spacing(d.x, d.y).median ≈ 40
    m = GG.design_mesh(d.x, d.y; dx_core = 20.0, dy_core = 20.0, pad_factor = 1.4, z_first = 20.0, z_factor = 1.4, depth_mult = 0.5)
    ix = (m.nx_pad + 1):(m.nx_pad + m.nx_core); iy = (m.ny_pad + 1):(m.ny_pad + m.ny_core)
    @test all(m.dx[ix] .== 20) && all(m.dy[iy] .== 20) && isodd(m.nx_core)
    @test !any(x -> x in m.x_edges, d.x) && !any(y -> y in m.y_edges, d.y)       # stations off the mesh planes
    @test m.pad_ok && m.pad_extent_km >= m.core_width_km                          # padding ≥ core
    @test m.maxper == 1 && m.outside == 0 && m.occupied > 0
    # regular 20 m layers down to the depth of investigation (half the core width), then padding growing
    # by the vertical factor until it is as thick as the regular block
    doi = 1000m.doi_km
    @test doi ≈ (m.x_edges[last(ix)+1] - m.x_edges[first(ix)]) / 2
    nb = m.nz_regular
    @test all(m.dz[1:nb] .== 20) && sum(m.dz[1:nb]) >= doi && sum(m.dz[1:nb-1]) < doi
    @test m.dz[nb+1:end] ≈ 20 .* 1.4 .^ (1:m.nz-nb) && sum(m.dz[nb+1:end]) >= sum(m.dz[1:nb])
    @test GG.design_mesh(d.x, d.y; dx_core = 20.0, dy_core = 20.0, nx_pad = 3, ny_pad = 2, pad_factor = 1.5, z_first = 5.0,
                         z_factor = 1.4, depth_mult = 0.5).dx[1:3] == 20 .* 1.5 .^ (3:-1:1)
    coarse = GG.design_mesh(d.x, d.y; dx_core = 100.0, dy_core = 100.0, nx_pad = 2, ny_pad = 2, pad_factor = 1.4,
                            z_first = 20.0, z_factor = 1.4, depth_mult = 0.5)
    @test coarse.maxper > 1 && !coarse.pad_ok
    s = GG.suggestions((; coarse..., site_surface = (impossible = [],), site_ground = (ok = true,)))
    @test !s.ok && any(occursin("Lateral padding", l) for l in s.lines) && any(occursin("More than one station", l) for l in s.lines)
end

@testset "topography, air layers and station columns" begin
    d = grid_data(z = -200.0)                                                   # 200 m above sea level
    slope(x, y) = 50 + 0.5y                                                     # ground 50 ± 100 m
    m = bundle(d; surface = slope, n_air = 3)
    ix = (m.nx_pad + 1):(m.nx_pad + m.nx_core); iy = (m.ny_pad + 1):(m.ny_pad + m.ny_core)
    @test m.nz_air == 3 && all(isnan, m.A[ix, iy, 1:3])                         # air layers over the core
    @test m.dz[1:3] == reverse(10 .* 1.35 .^ (0:2))                             # thinnest at the ground (layer / 2)
    @test m.z_edges[4] ≈ -slope(0, m.y_centers[last(iy)])                       # earth starts at the core's highest ground
    relief = 0.5 * (m.y_centers[last(iy)] - m.y_centers[first(iy)])
    @test m.nz_regular == ceil(Int, (relief + 1000m.doi_km) / 20 - 1e-9)       # regular across the relief and the depth
    for i in (1, m.nx), j in (first(iy), last(iy))                              # earth starts below the local ground
        k = findfirst(!isnan, m.A[i, j, :])
        @test m.z_centers[k] >= -slope(m.x_centers[i], m.y_centers[j]) && m.z_centers[k-1] < -slope(m.x_centers[i], m.y_centers[j])
    end
    @test m.site_ground.ok && m.site_surface.nadjusted == 0 && GG.suggestions(m).ok
    @test m.mask == .!isnan.(m.A)
    # stations at sea level stand in the eastern ground: their columns are lowered until they are on or above it
    low = grid_data(z = 0.0)
    l = bundle(low; surface = slope)
    @test l.site_surface.nadjusted > 0 && l.site_ground.ok
    r = l.site_surface.adjustments[1]
    @test r.surface_to > r.surface_from && l.z_edges[findfirst(!isnan, l.A[r.ix, r.iy, :])] >= r.site_z
    model = GG.mesh_model(m)                                                    # usable as is
    op = GravForward(model, d)
    @test length(forward(op, fill(0.01, length(op.cells)))) == ndata(d)
    @test check_frame(d, model) === nothing
end

@testset "GeoTIFF topography" begin
    d = grid_data(z = -400.0)
    d.origin .= [60.0, 25.0]
    tr = GG.Proj.Transformation("+proj=tmerc +lat_0=60 +lon_0=25 +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs",
                                "EPSG:4326"; always_xy = true)
    for i in eachindex(d.x)
        lon, lat = tr((d.y[i] + 500000.0, d.x[i])); d.lat[i], d.lon[i] = lat, lon
    end
    dir = mktempdir()
    dp = write_data_grav(joinpath(dir, "obs.dat"), d)
    # DEM in geographic coordinates: elevation rising 1 m per 0.0001° of longitude
    lon = range(24.97, 25.03; length = 241); lat = range(59.97, 60.03; length = 241)
    elev(lo, la) = 100 + 10_000 * (lo - 25)
    Z = Float32[elev(lo, la) for lo in lon, la in reverse(lat)]                 # (x, y), y from the north
    tif = joinpath(dir, "dem.tif")
    GG.ArchGDAL.create(tif; driver = GG.ArchGDAL.getdriver("GTiff"), width = length(lon), height = length(lat), nbands = 1, dtype = Float32) do ds
        Δ = step(lon)
        GG.ArchGDAL.setgeotransform!(ds, [first(lon) - Δ / 2, Δ, 0, last(lat) + Δ / 2, 0, -Δ])
        GG.ArchGDAL.setproj!(ds, GG.ArchGDAL.toWKT(GG.ArchGDAL.importEPSG(4326)))
        GG.ArchGDAL.write!(ds, Z, 1)
    end
    m = MakeMesh3D(dp; topo_file = tif, mode = :nogui, out_model = joinpath(dir, "start.den"))
    @test m.topo_active && m.topo_name == "dem.tif"
    for (i, j) in ((m.nx ÷ 2, 1), (m.nx ÷ 2, m.ny), (1, m.ny ÷ 2))
        lo, la = tr((m.y_centers[j] + 500000.0, m.x_centers[i]))
        @test isapprox(m.surface_local[i, j], -elev(clamp(lo, 24.97, 25.03), la); atol = 1.0)   # z down = −elevation
    end
    @test isfile(joinpath(dir, "start.den")) && isfile(joinpath(dir, "mesh_start_model.mask"))
    nolatlon = write_data_grav(joinpath(dir, "frame.dat"), grid_data())
    @test_throws ErrorException MakeMesh3D(nolatlon; topo_file = tif, mode = :nogui)
end

@testset "frame to CRS on a national-grid frame" begin
    # stations on the EPSG:3067 grid near Kuusamo, 2° east of its central meridian: there the survey's local TM is
    # rotated and scaled against the frame, which a constant shift cannot absorb
    tr = GG.Proj.Transformation("EPSG:3067", "EPSG:4326"; always_xy = true)
    e = [590_000.0 + a for a in 0:1000:10_000 for _ in 0:1000:10_000]; n = [7_330_000.0 + b for _ in 0:1000:10_000 for b in 0:1000:10_000]
    d = station_data(n, e, fill(-300.0, length(e)))
    ll = [tr((a, b)) for (a, b) in zip(e, n)]
    d.lon .= first.(ll); d.lat .= last.(ll); d.origin .= [mean(d.lat), mean(d.lon)]
    ctx = GG._local_tm_to_target_context(d, "EPSG:3067")
    err = maximum(hypot((GG._local_xy_to_target_xy(ctx, d.x[k], d.y[k]) .- (e[k], n[k]))...) for k in eachindex(e))
    report("frame → EPSG:3067 at the stations (m)", err, 0.5)
    back = GG.Proj.inv(ctx.trans)
    p = GG._target_xy_to_local_xy(ctx, back, 610_000.0, 7_350_000.0)            # 10 km beyond the stations
    report("EPSG:3067 → frame, 10 km outside the stations (m)", hypot(p[1] - 7_350_000, p[2] - 610_000), 2.0)
end

@testset "lakes and uniform layers" begin
    # geometry: the even-odd rule makes an island dry; distance to the nearest shore
    sq(a) = [(-a, -a), (a, -a), (a, a), (-a, a)]
    @test GG._in_rings(0.0, 0.0, [sq(10.0)]) && !GG._in_rings(0.0, 0.0, [sq(10.0), sq(2.0)]) && GG._in_rings(5.0, 0.0, [sq(10.0), sq(2.0)])
    @test GG._ring_distance(0.0, 0.0, [sq(10.0)]) ≈ 10 && GG._ring_distance(12.0, 0.0, [sq(10.0)]) ≈ 2
    # flat ground at 200 m a.s.l., a lake at its level with a modelled bed 30 m deep, reached 40 m from the shore;
    # stations 30 m above the water; 5 m layers through the relief (the bed) and 4 more, then growing
    d = grid_data(z = -230.0)
    dir = mktempdir(); dp = write_data_grav(joinpath(dir, "obs.dat"), d)
    lake = (polygon = sq(70.0), level = 200.0, depth = 30.0, taper = 40.0)
    m = MakeMesh3D(dp; mode = :nogui, surface = (x, y) -> 200.0, lakes = [lake], lakes_crs = "frame", z_first = 5.0,
                   uniform_layers = 4)
    @test m.nlakes == 1 && count(m.water) > 0 && all(==(1000.0 - 2670.0), m.A[m.water]) && !any(m.mask .& m.water)
    inlake = [abs(x) < 70 && abs(y) < 70 for x in m.x_centers, y in m.y_centers]
    @test all(inlake[I[1], I[2]] for I in findall(m.water))                   # water only under the lake
    i, j = argmin(abs.(m.x_centers)), argmin(abs.(m.y_centers))               # the middle: 30 m of water, 6 layers
    @test count(m.water[i, j, :]) == 6 && all(-200 .< m.z_centers[m.water[i, j, :]] .< -170)
    @test m.z_edges[findfirst(m.mask[i, j, :])] ≈ -170                         # the earth starts at the bed
    shore = [GG._ring_distance(m.x_centers[I[1]], m.y_centers[I[2]], [sq(70.0)]) for I in findall(m.water)]
    @test minimum(shore) < 40                                                  # the taper reaches the shallows
    @test m.nz_regular == 6 + 4 && all(m.dz[1:10] .== 5)                       # 30 m of relief + 4 layers
    @test m.dz[11:end] ≈ 5 .* 1.4 .^ (1:m.nz-10) && sum(m.dz) >= 2 * (30 + 1000m.doi_km)
    # the written model: the ground at 200 m a.s.l. is z = 0
    w = load_model(joinpath(dir, "mesh_start_model.den"))
    @test w.datum == 200 && w.Z[1, 1] == 0 && check_frame(d, w) === nothing
    # the mask file: 9 for water, not variable; water_mask reads it back
    r = load_model(joinpath(dir, "mesh_start_model.den"))
    @test all(==(-1670.0), r.A[m.water])
    # the water contrast follows the reduction density
    @test all(==(1000.0 - 2200.0), MakeMesh3D(dp; mode = :nogui, surface = (x, y) -> 200.0, lakes = [lake], lakes_crs = "frame",
                                              z_first = 5.0, uniform_layers = 4, ρ_ref = 2200.0,
                                              out_model = joinpath(dir, "r2200.den")).A[m.water])
    mp = joinpath(dir, "mesh_start_model.mask")
    @test water_mask(mp, r) == m.water && load_mask(mp, r) == m.mask
    @test count(==(9.0), load_model(mp).A) == count(m.water)
    @test_throws ArgumentError write_mask(joinpath(dir, "bad.mask"), r, m.water; water = m.water)
    # without a depth the ground is the bed: a DEM at the water level gives no water, a lowered ground gives some
    flat = MakeMesh3D(dp; mode = :nogui, surface = (x, y) -> 200.0, lakes = [(polygon = sq(70.0), level = 200.0)],
                      lakes_crs = "frame", z_first = 5.0, uniform_layers = 4)
    @test count(flat.water) == 0
    dip = MakeMesh3D(dp; mode = :nogui, surface = (x, y) -> abs(x) < 70 && abs(y) < 70 ? 180.0 : 200.0,
                     lakes = [(polygon = sq(70.0), level = 200.0)], lakes_crs = "frame", z_first = 5.0, uniform_layers = 4)
    @test count(dip.water[i, j, :]) == 4
    # a station below the lake level is an error; the default layers are unchanged by the option
    low = write_data_grav(joinpath(dir, "low.dat"), grid_data(z = -190.0))
    @test_throws ErrorException MakeMesh3D(low; mode = :nogui, surface = (x, y) -> 150.0, lakes = [lake], lakes_crs = "frame")
    @test MakeMesh3D(dp; mode = :nogui, surface = (x, y) -> 200.0, z_first = 5.0).nz_regular ==
          ceil(Int, 1000GG.design_mesh(d.x, d.y; dx_core = 20.0, dy_core = 20.0, pad_factor = 1.4, z_first = 5.0,
                                       z_factor = 1.4, depth_mult = 0.5).doi_km / 5 - 1e-9)
    # read_lakes from a shapefile with an island, its level from an attribute, the bed from a constant
    shp = joinpath(dir, "lakes.shp")
    poly = GG.Shapefile.Polygon(GG.Shapefile.Rect(-1.0, -1.0, 1.0, 1.0), Int32[0, 5],
                                [GG.Shapefile.Point(x, y) for (x, y) in vcat(sq(1.0), [(-1.0, -1.0)], reverse(sq(0.2)), [(-0.2, -0.2)])])
    GG.Shapefile.write(shp, [(geometry = poly, level = 123.0, name = "A")])
    l = read_lakes(shp; level = "level", depth = 12.0, taper = 50.0, crs = "EPSG:4326")
    @test length(l) == 1 && l[1].level == 123 && l[1].depth == 12 && l[1].taper == 50 && length(l[1].polygon) == 2
    @test GG._in_rings(0.5, 0.5, l[1].polygon) && !GG._in_rings(0.0, 0.0, l[1].polygon)    # the island is dry
    @test_throws ErrorException read_lakes(shp; level = 1.0, crs = "EPSG:4326", select = r -> r.name == "B")
end

@testset "MakeMesh3D and the designer" begin
    d = grid_data()
    dir = mktempdir()
    dp = write_data_grav(joinpath(dir, "obs.dat"), d)
    m = MakeMesh3D(dp; mode = :nogui)
    r = load_model(joinpath(dir, "mesh_start_model.den"))
    @test same_mesh(r, GG.mesh_model(m)) && all(r.A .== 0) && m.dx_core == 20
    # z = 0 at the highest ground: the top of the highest earth cell, its elevation the datum
    k = findfirst(k -> any(!isnan, r.A[:, :, k]), 1:r.nz)
    @test r.z[k] == 0 && r.datum == -m.z_edges[k] && check_frame(d, r) === nothing
    @test r.z ≈ m.z_edges .+ r.datum && minimum(r.Z) == 0
    # the same earth as the bundle at sea level: the same gz
    sl = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
    report("MakeMesh model vs the sea-level bundle", maximum(abs, forward(GravForward(r, d), fill(100.0, count(active_cells(r)))) .-
           forward(GravForward(sl, d), fill(100.0, count(active_cells(sl))))), 1e-10)
    # the highest earth cell in the padding: z = 0 is there, not over the core
    hill = MakeMesh3D(dp; mode = :nogui, surface = (x, y) -> abs(x) > 400 ? 120.0 : 0.0, n_pad = 4,
                      out_model = joinpath(dir, "hill.den"), out_mask = joinpath(dir, "hill.mask"))
    h = load_model(joinpath(dir, "hill.den"))
    kh = findfirst(k -> any(!isnan, h.A[:, :, k]), 1:h.nz)
    @test h.z[kh] == 0 && minimum(h.Z) == 0 && h.datum == -hill.z_edges[kh]
    # data in a rotated frame: the mesh carries the data's azimuth, and the pair is accepted
    rp = write_data_grav(joinpath(dir, "rot.dat"), rotate_data(d, 30.0))
    MakeMesh3D(rp; mode = :nogui, out_model = joinpath(dir, "rot.den"), out_mask = joinpath(dir, "rot.mask"))
    rm_ = load_model(joinpath(dir, "rot.den"))
    @test rm_.rotation == 30.0 && check_frame(load_data_grav(rp), rm_) === nothing
    @test m.site_surface.nadjusted == 0 && m.site_ground.ok                   # stations 30 m above flat ground
    mask = load_mask(joinpath(dir, "mesh_start_model.mask"), r)
    @test mask == active_cells(r)
    fig, mesh = MakeMesh3D(dp; out_model = joinpath(dir, "gui.den"), interactive = false)
    @test same_mesh(GG.mesh_model(mesh[]), GG.mesh_model(m))
    @test filesize(save_figure(joinpath(dir, "designer.png"), fig)) > 10_000
    @test_throws ArgumentError MakeMesh3D(dp; mode = :web)
end
