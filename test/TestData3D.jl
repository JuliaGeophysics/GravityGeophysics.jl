# Data: frame, CSV / Geosoft XYZ import, native file round trip, rotation, draping, frame checks.

import Proj

const LAT0, LON0 = 63.671883, 26.065008

# MT data file with sites on a lat/lon grid whose x/y are linear in lat/lon, turned by θ (as in
# real MT files, which are not made with one projection)
function write_mt_sites(path; θ = 0.5, shift = (120.0, -30.0))
    lat = [LAT0 + a for a in -0.4:0.2:0.4 for _ in 1:5]; lon = [LON0 + b for _ in 1:5 for b in -0.8:0.4:0.8]
    xt, yt = 111_195.0 .* (lat .- LAT0), 111_195.0 * cosd(LAT0) .* (lon .- LON0)
    x = cosd(θ) .* xt .+ sind(θ) .* yt .+ shift[1]
    y = -sind(θ) .* xt .+ cosd(θ) .* yt .+ shift[2]
    open(path, "w") do io
        println(io, "# test sites\n# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error")
        println(io, "> Full_Impedance\n> exp(+i\\omega t)\n> [mV/km]/[nT]\n> 0.00\n> $LAT0 $LON0\n> 1 $(length(lat))")
        for k in eachindex(lat), c in ("ZXY", "ZYX")
            @printf(io, "1.0 S%02d %.6f %.6f %.3f %.3f 0.0 %s 1.0 1.0 0.1\n", k, lat[k], lon[k], x[k], y[k], c)
        end
    end
    path
end

flat_model(; rotation = 0.0) = Model(fill(1000.0, 20), fill(1000.0, 20), fill(100.0, 10), zeros(20, 20, 10);
                                     origin = [-10_000.0, -10_000.0, 0.0], rotation = rotation)

@testset "frame" begin
    f = Frame(LAT0, LON0)
    x, y = to_frame(f, [LAT0], [LON0])
    @test abs(x[1]) < 1e-6 && abs(y[1]) < 1e-6
    @test abs(north_azimuth(f, LAT0, LON0)) < 1e-6
    # east of the central meridian true north leans west: azimuth ≈ −Δλ sin φ
    @test north_azimuth(f, LAT0, LON0 + 1) ≈ -sind(LAT0) rtol = 0.01
    path = write_mt_sites(joinpath(mktempdir(), "mt.dat"))
    g = load_frame(path)
    @test g.rms < 1e-3 && g.lat0 == LAT0
    @test north_azimuth(g, LAT0, LON0) ≈ -0.5 atol = 1e-6
    nohdr = joinpath(dirname(path), "nohdr.dat"); write(nohdr, "# no header\n1.0 S01 63 26 0 0 0 ZXY 1 1 1\n")
    @test_throws ArgumentError load_frame(nohdr)
    # a projected CRS as the frame: x, y are its northing, easting (ETRS-TM35FIN at the Koillismaa Deep Hole)
    tm = Frame("EPSG:3067", (65.6, 65.86), (28.47, 29.01))
    @test tm.rms < 0.1 && occursin("EPSG:3067", tm.source)
    p = Proj.Transformation("EPSG:3067", "EPSG:4326"; always_xy = true)((579_780.0, 7_290_550.0))
    x, y = to_frame(tm, [p[2]], [p[1]])
    @test abs(x[1] - 7_290_550.0) < 0.1 && abs(y[1] - 579_780.0) < 0.1
    # grid north at 28.74° E, east of the 27° E central meridian, leans west of true north: ≈ −Δλ sin φ
    @test north_azimuth(tm, p[2], p[1]) ≈ -(p[1] - 27) * sind(p[2]) rtol = 0.01
end

@testset "CSV import" begin
    dir = mktempdir()
    p = joinpath(dir, "survey.csv")
    write(p, """
    # airborne test
    LINE,LAT,LON,ALT,GZ,ERR
    L10,63.60,26.00,150.0,12.5,1.0
    L10,63.61,26.00,151.0,*,1.0
    L20,63.62,26.10,152.0,-3.25,1.5
    "L20",63.63,26.10,153.0,4.0,2.0
    """)
    d = load_data_grav(p; columns = (lat = "LAT", lon = "lon", elevation = "ALT", value = "GZ", error = "ERR", line = "LINE"))
    @test ndata(d) == 3 && d.component == :gz && d.code == ["L10:1", "L20:3", "L20:4"]
    @test d.z == [-150.0, -152.0, -153.0] && d.value == [12.5, -3.25, 4.0]
    f = Frame(d.origin...)
    x, y = to_frame(f, d.lat, d.lon)
    @test d.x == x && d.y == y
    @test d.origin ≈ [sum(d.lat) / 3, sum(d.lon) / 3]
    @test_throws ArgumentError load_data_grav(p; component = :tmi,
                                         columns = (lat = "LAT", lon = "LON", elevation = "ALT", value = "GZ", error = "ERR"))
    @test_throws ArgumentError load_data_grav(p; columns = (lat = "LAT", lon = "LON", elevation = "ALT", value = "GRAV", error = "ERR"))
    # projected coordinates give the same stations
    tr = Proj.Transformation("EPSG:4326", "EPSG:32635"; always_xy = true)
    q = joinpath(dir, "utm.csv")
    open(q, "w") do io
        println(io, "E,N,Z,V,S")
        for k in 1:ndata(d)
            e, n = tr((d.lon[k], d.lat[k])); @printf(io, "%.4f,%.4f,%.1f,%.3f,%.3f\n", e, n, -d.z[k], d.value[k], d.error[k])
        end
    end
    u = load_data_grav(q; crs = "EPSG:32635", frame = Frame(d.origin...),
                  columns = (easting = "E", northing = "N", elevation = "Z", value = "V", error = "S"))
    @test maximum(abs, u.x .- d.x) < 1e-3 && maximum(abs, u.y .- d.y) < 1e-3
end

@testset "Geosoft XYZ import" begin
    p = joinpath(mktempdir(), "survey.xyz")
    write(p, """
    / Geosoft export
    /     X        Y      Z     GRAV   ERR
    Line 1010
     100.0    200.0  -80.0   5.25   0.5
     110.0    200.0  -80.0      *   0.5
    Tie 20
     120.0    250.0  -81.0  -1.75   0.5
    """)
    d = load_data_grav(p; columns = (north = "Y", east = "X", elevation = "Z", value = "GRAV", error = "ERR"))
    @test d.code == ["L1010:1", "T20:3"] && d.x == [200.0, 250.0] && d.y == [100.0, 120.0] && d.z == [80.0, 81.0]
    @test d.north == 0 && isnan(d.lat[1])
    e = load_data_grav(p; columns = (north = 2, east = 1, elevation = 3, value = 4, error = 5))
    @test e.value == d.value
end

@testset "native file round trip" begin
    path = write_mt_sites(joinpath(mktempdir(), "mt.dat"))
    p = joinpath(dirname(path), "s.csv")
    write(p, "lat,lon,h,v,e\n63.6,26.0,100,1.5,0.1\n63.7,26.1,120,-2.5,0.2\n63.8,26.2,140,3.5,0.3\n")
    d = load_data_grav(p; frame = load_frame(path),
                  columns = (lat = "lat", lon = "lon", elevation = "h", value = "v", error = "e"))
    d = rotate_data(rotate_data(d, 20.0), -5.0; kind = :strike)
    out = write_data_grav(joinpath(dirname(path), "d.dat"), d; comment = "test")
    r = load_data_grav(out)
    @test r.component == :gz
    @test r.rotations == [(:mesh, 20.0), (:strike, -5.0)] && r.azimuth == 15.0
    @test r.north ≈ d.north atol = 1e-6
    @test r.frame_rms ≈ d.frame_rms atol = 1e-3
    @test maximum(abs, r.x .- d.x) < 1e-3 && r.z == d.z && r.code == d.code
    @test isapprox(r.value, d.value; rtol = 1e-6) && isapprox(r.error, d.error; rtol = 1e-6)
    @test r.lat ≈ d.lat atol = 1e-6
    bad = joinpath(dirname(path), "bad.dat")
    write(bad, replace(read(out, String), "> 3" => "> 4"))
    @test_throws ArgumentError load_data_grav(bad)
    # a magnetic file (one header line more) is not a gravity file
    mag = joinpath(dirname(path), "mag.dat")
    write(mag, replace(read(out, String), "> gz\n" => "> TMI\n> 51000.000 75.500000 9.250000\n"))
    @test_throws ArgumentError load_data_grav(mag)
end

@testset "rotation" begin
    p = joinpath(mktempdir(), "s.csv")
    write(p, "n,e,z,v,s\n100,0,-50,1,1\n0,200,-50,2,1\n-300,-400,-50,3,1\n")
    d = load_data_grav(p; columns = (north = "n", east = "e", elevation = "z", value = "v", error = "s"))
    r = rotate_data(d, 90.0)
    @test isapprox(r.x, d.y; atol = 1e-12) && isapprox(r.y, -d.x; atol = 1e-12) && r.value == d.value
    back = rotate_data(r, -90.0)
    @test isapprox(back.x, d.x; atol = 1e-9) && isapprox(back.y, d.y; atol = 1e-9) && back.azimuth == 0
    @test_throws ArgumentError rotate_data(d, 1.0; kind = :declination)
end

@testset "draping and frame checks" begin
    m = flat_model()
    m.A[:, :, 1] .= NaN; m.A[5, 5, 2] = NaN; m = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
    p = joinpath(mktempdir(), "s.csv")
    cols = (north = "n", east = "e", elevation = "z", value = "v", error = "s")
    write(p, "n,e,z,v,s\n-5500,-5500,0,1,1\n500,500,0,2,1\n")
    d = drape!(load_data_grav(p; columns = cols), m, 30.0)
    @test d.z == [m.Z[5, 5] - 30, m.Z[11, 11] - 30] == [170.0, 70.0]
    @test check_frame(d, m) === nothing
    # ground stations: on the top of the earth cells is allowed, below is not
    g = drape!(deepcopy(d), m, 0.0)
    @test g.z == [m.Z[5, 5], m.Z[11, 11]] && check_frame(g, m) === nothing
    @test_throws ArgumentError drape!(d, m, -1.0)
    write(p, "n,e,z,v,s\n500,500,-150,1,1\n")
    @test_throws ArgumentError check_frame(load_data_grav(p; columns = cols), m)   # below ground
    write(p, "n,e,z,v,s\n50000,500,500,1,1\n")
    @test_throws ArgumentError check_frame(load_data_grav(p; columns = cols), m)   # outside
    # a rotated mesh needs data turned to its frame
    write(p, "n,e,z,v,s\n1000,0,500,1,1\n")
    mr = flat_model(rotation = 30.0)
    d = load_data_grav(p; columns = cols)
    @test_throws ArgumentError check_frame(d, mr)
    r = rotate_data(d, mr.rotation)
    @test check_frame(r, mr) === nothing
    @test r.x[1] ≈ 1000cosd(30) && r.y[1] ≈ -1000sind(30)
end

@testset "stations on a mesh with a datum" begin
    # flat ground at 250 m a.s.l. as z = 0 of the mesh; the data keep z = −elevation
    m = Model(fill(1000.0, 20), fill(1000.0, 20), fill(100.0, 10), zeros(20, 20, 10); origin = [-10_000.0, -10_000.0, 0.0],
              datum = 250.0)
    p = joinpath(mktempdir(), "s.csv")
    cols = (north = "n", east = "e", elevation = "z", value = "v", error = "s")
    write(p, "n,e,z,v,s\n500,500,250,1,1\n-500,500,260,1,1\n")
    d = load_data_grav(p; columns = cols)
    @test d.z == [-250.0, -260.0] && check_frame(d, m) === nothing        # on the ground and 10 m above it
    write(p, "n,e,z,v,s\n500,500,249,1,1\n")
    @test_throws ArgumentError check_frame(load_data_grav(p; columns = cols), m)   # 1 m under the ground
    g = drape!(deepcopy(d), m, 30.0)
    @test g.z == [-280.0, -280.0]                                          # −elevation, 30 m above 250 m a.s.l.
end
