# Model files: layout, round trip, resistivity meshes, errors.
# data/resistivity.rho is a LOGE resistivity model, 5 × 4 × 6, one air layer plus two
# topography air cells at [1:2, 1, 2], a 0.3 Ω·m water cell at [4, 3, 2], 10 Ω·m at [3, 2, 4],
# and ρ = 1000 + 100i + 10j + k in layers 5 and 6.

function demo_model()
    dx = [3000.0, 1000, 1000, 1000, 3000]; dy = [2000.0, 500, 500, 2000]; dz = [400.0, 100, 100, 200, 400, 800]
    A = [1e-3 * i + 1e-4 * j + 1e-5 * k for i in 1:5, j in 1:4, k in 1:6]
    A[:, :, 1] .= NaN
    A[1:2, 1, 2] .= NaN
    A[3, 2, 4] = -2.5e-4                 # negative contrasts are allowed
    A[2, 3, 5] = 0.0
    Model(dx, dy, dz, A; origin = [-4500.0, -2500.0, -400.0], rotation = 12.5)
end

@testset "x stored backwards within each y row" begin
    path = joinpath(mktempdir(), "tiny.den")
    write(path, "2 2 1 0 LINEAR\n10 20\n5 5\n7\n 1 2\n 3 4\n0 0 0\n0\n")
    m = load_model(path)
    @test m.A[:, :, 1] == [2.0 4.0; 1.0 3.0]
    @test m.x == [0, 10, 30] && m.y == [0, 5, 10] && m.z == [0, 7]
    @test m.cx == [5, 20]
end

@testset "round trip" begin
    m = demo_model()
    path = write_model(joinpath(mktempdir(), "demo.den"), m)
    header = first(l for l in eachline(path) if !startswith(l, "#"))
    @test header == "5 4 6 1 LINEAR"
    r = load_model(path)
    @test same_mesh(m, r)
    @test isequal(isnan.(r.A), isnan.(m.A))
    act = active_cells(m)
    @test maximum(abs, r.A[act] .- m.A[act]) <= 1e-9 * maximum(abs, m.A[act])
    @test r.A[3, 2, 4] < 0 && r.A[2, 3, 5] == 0
    @test r.rotation == 12.5 && r.origin == [-4500.0, -2500.0, -400.0]
    @test count(air_cells(r)) == 5 * 4 + 2
    # surface: z of the first earth cell under the air
    @test r.Z[1, 1] == -400 + 400 + 100 && r.Z[3, 3] == -400 + 400
    @test r.npad == (1, 1)
    @test r.cz[1] == -200
end

@testset "resistivity file supplies the mesh" begin
    path = joinpath(DATA, "resistivity.rho")
    r = GravityGeophysics._model_read(path)
    @test r.scale == "LOGE"
    @test isnan(r.A[1, 1, 2]) && isnan(r.A[2, 1, 2]) && !isnan(r.A[3, 1, 2])
    @test all(isnan, r.A[:, :, 1])
    @test r.A[4, 3, 2] ≈ 0.3 rtol = 1e-4
    @test r.A[3, 2, 4] ≈ 10 rtol = 1e-4
    @test all(isapprox(r.A[i, j, k], 1000 + 100i + 10j + k; rtol = 1e-4) for i in 1:5, j in 1:4, k in 5:6)
    m = load_mesh(path; value = 0.01)
    @test count(air_cells(m)) == 22 && all(m.A[active_cells(m)] .== 0.01)
    @test m.dx == [3000.0, 1000, 1000, 1000, 3000] && m.origin == [-4500.0, -2500.0, -400.0]
    w = load_mesh(path; water_log10 = 0.0)
    @test count(air_cells(w)) == 23 && isnan(w.A[4, 3, 2])
    @test_throws ArgumentError load_model(path)
end

@testset "errors" begin
    dir = mktempdir()
    bad(s) = (p = joinpath(dir, "bad.den"); write(p, s); p)
    @test_throws ArgumentError load_model(bad("2 1 1 0 LINEAR\n1 1\n1\n1\n 1\n"))           # too few values
    @test_throws ArgumentError load_model(bad("2 1 1 0\n1 1\n1\n1\n 1 2\n"))                # no scale
    @test_throws ArgumentError load_model(bad("2 1 1 0 LINEAR\n1 -1\n1\n1\n 1 2\n"))        # width ≤ 0
    @test_throws ArgumentError load_model(bad("2 1 1 1 LINEAR\n1 1\n1\n1\n 1 2\n"))         # nzAir ≥ nz
    m = demo_model(); m.A[2, 2, 3] = Inf
    @test_throws ArgumentError write_model(joinpath(dir, "inf.den"), m)
end

@testset "datum" begin
    dir = mktempdir()
    m = demo_model()
    @test m.datum == 0
    s = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin, rotation = m.rotation, datum = 344.125)
    p = write_model(joinpath(dir, "s.den"), s)
    @test any(==("# Datum elevation (m) : 344.125"), readlines(p))
    r = load_model(p)
    @test r.datum == 344.125 && same_mesh(r, s) && !same_mesh(r, m)
    # a file without the line (MT, magnetic, older) is at sea level
    write(joinpath(dir, "old.den"), "2 2 1 0 LINEAR\n10 20\n5 5\n7\n 1 2\n 3 4\n0 0 0\n0\n")
    @test load_model(joinpath(dir, "old.den")).datum == 0
    @test load_mesh(p).datum == 344.125
end
