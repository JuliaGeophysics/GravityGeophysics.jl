# MeshToMesh: a density model onto a finer mesh with its own topography; mask and data files.

@testset "MeshToMesh" begin
    dir = mktempdir()
    src = Model(fill(100.0, 8), fill(100.0, 8), fill(50.0, 6); origin = [-400.0, -400.0, 0.0])
    src.A[5:8, :, :] .= 250.0                                     # northern half dense
    sp = write_model(joinpath(dir, "coarse.den"), src)
    d = station_data([-100.0, 100.0], [0.0, 0.0], [-30.0, -30.0])
    dp = write_data_grav(joinpath(dir, "obs.dat"), d)
    # target: a resistivity mesh, 50 m cells, top layer air over the west half and sea water east of y = 200 m
    ρ = fill(100.0, 16, 16, 12); ρ[:, 1:8, 1] .= 1e17; ρ[:, 13:16, 1] .= 0.3
    open(joinpath(dir, "fine.rho"), "w") do io
        println(io, "16 16 12 0 LOGE"); println(io, join(fill(50.0, 16), " ")); println(io, join(fill(50.0, 16), " "))
        println(io, join(fill(25.0, 12), " "))
        for k in 1:12, j in 1:16
            println(io, join((log(ρ[i, j, k]) for i in 16:-1:1), " "))
        end
        println(io, "-400 -400 0"); println(io, "0")
    end
    r = MeshToMesh(sp, dp, joinpath(dir, "fine.rho"); outdir = joinpath(dir, "out"), water_log10 = 0.0)
    m = load_model(r.model)
    @test (m.nx, m.ny, m.nz) == (16, 16, 12) && r.dst_grid == (16, 16, 12)
    @test all(isnan, m.A[:, 1:8, 1]) && all(isnan, m.A[:, 13:16, 1])       # target air and sea water stay out
    @test all(m.A[9:16, :, 2:end] .== 250.0) && all(m.A[1:8, :, 2:end] .== 0.0)
    @test r.n_air == 16 * 12 && r.cells_outside_source == 0
    @test load_mask(r.mask, m) == .!isnan.(m.A)
    @test load_data_grav(r.data).z == d.z                                     # heights kept
    lin = MeshToMesh(sp, dp, joinpath(dir, "fine.rho"); outdir = joinpath(dir, "lin"), method = :linear, snap_sites = true,
                     clearance = 5.0)
    ml = load_model(lin.model)
    @test 0 < ml.A[8, 9, 6] < 250.0 && ml.A[12, 9, 6] == 250.0                  # interpolated across the contact
    dl = load_data_grav(lin.data)
    @test all(dl.z .≈ 25.0 - 5.0) || all(isapprox.(dl.z, ml.Z[GG._column_index(ml, d.x[1], d.y[1])...] - 5.0))
    @test_throws ErrorException MeshToMesh(sp, dp, joinpath(dir, "missing.rho"))
    # a source with its z = 0 at 100 m a.s.l. lands 100 m higher on a sea-level target
    hi = Model(src.dx, src.dy, src.dz, src.A; origin = src.origin, datum = 100.0)
    hi.A[:, :, 1:2] .= 0.0; hi.A[:, :, 3:end] .= 250.0          # contact 100 m below its z = 0: at sea level
    hp = write_model(joinpath(dir, "high.den"), hi)
    target = write_model(joinpath(dir, "target.den"), Model(fill(100.0, 8), fill(100.0, 8), fill(25.0, 12); origin = [-400.0, -400.0, -100.0]))
    t = load_model(MeshToMesh(hp, dp, target; outdir = joinpath(dir, "datum")).model)
    @test t.datum == 0 && all(t.A[:, :, 1:4] .== 0.0) && all(t.A[:, :, 5:end] .== 250.0)
    # a target turned 90° clockwise: north is −y′ there, so a body 200–400 m north lands at y′ −400…−200
    body = Model(fill(100.0, 8), fill(100.0, 8), fill(50.0, 4); origin = [-400.0, -400.0, 0.0])
    body.A[7:8, 4:5, :] .= 250.0
    bp = write_model(joinpath(dir, "body.den"), body)
    tq = write_model(joinpath(dir, "turned.den"), Model(fill(100.0, 8), fill(100.0, 8), fill(50.0, 4);
                                                        origin = [-400.0, -400.0, 0.0], rotation = 90.0))
    r90 = MeshToMesh(bp, dp, tq; outdir = joinpath(dir, "turned"))
    t90 = load_model(r90.model)
    I = findall(t90.A .> 0)
    @test t90.rotation == 90 && extrema(t90.cx[i[1]] for i in I) == (-50.0, 50.0) && extrema(t90.cy[i[2]] for i in I) == (-350.0, -250.0)
    d90 = load_data_grav(r90.data)
    @test d90.azimuth == 90 && isapprox(d90.x, d.y; atol = 1e-9) && isapprox(d90.y, -d.x; atol = 1e-9)
    @test check_frame(d90, t90) === nothing
end
