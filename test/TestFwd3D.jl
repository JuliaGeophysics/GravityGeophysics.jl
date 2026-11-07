# Forward operator: dense = matrix-free, adjoint, model vectors, frame check, file workflow.

# 6 × 5 × 4 mesh with topography (two air cells) and a random density model
function small_problem(; seed = 1, height = 30.0)
    m = Model([200.0, 50, 50, 50, 50, 200], [200.0, 50, 50, 50, 200], [25.0, 25, 50, 100])
    m.A[1:2, 1, 1] .= NaN
    m = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
    x = range(-150, 150; length = 7); y = range(-100, 100; length = 5)
    d = station_data([a for a in x for _ in y], [b for _ in x for b in y], fill(-height, length(x) * length(y)))
    ρ = -200 .+ 400 .* abs.(sin.(seed .* (1:count(active_cells(m)))))
    m, d, ρ
end

@testset "dense = matrix-free, adjoint" begin
    m, d, ρ = small_problem()
    D, F = GravForward(m, d), GravForward(m, d; solver = :matrixfree)
    v = cos.(1:length(ρ)); w = sin.(1:ndata(d))
    @test isapprox(forward(D, ρ), forward(F, ρ); rtol = 1e-12)
    @test isapprox(jvec(D, ρ, v), jvec(F, ρ, v); rtol = 1e-12)
    @test isapprox(jtvec(D, ρ, w), jtvec(F, ρ, w); rtol = 1e-12)
    for op in (D, F)
        Jv = jvec(op, ρ, v)
        @test abs(dot(Jv, w) - dot(v, jtvec(op, ρ, w))) / (norm(Jv) * norm(w)) < 1e-12
    end
end

@testset "linear data and model vectors" begin
    m, d, ρ = small_problem()
    op = GravForward(m, d)
    @test forward(op, 2ρ) ≈ 2forward(op, ρ)
    @test jvec(op, ρ, ρ) == forward(op, ρ)
    @test isequal(isnan.(model_array(op, ρ)), air_cells(m))
    mm = Model(m.dx, m.dy, m.dz, model_array(op, ρ); origin = m.origin)
    @test model_vector(op, mm) == ρ
    @test_throws ArgumentError forward(op, ρ[1:end-1])
    @test_throws ArgumentError GravForward(m, rotate_data(d, 5.0))
    # a dense body below gives positive (downward) gz
    @test all(>(0), forward(op, fill(100.0, length(ρ))))
end

@testset "ground stations" begin
    # stations on the top of the earth cells (z = 0) and on the top of the topography column
    m, d, ρ = small_problem(; height = 0.0)
    @test check_frame(d, m) === nothing
    g = forward(GravForward(m, d), ρ)
    @test all(isfinite, g)
    # continuity: 1 mm above the ground gives the same gz to 1e-4
    _, d1, _ = small_problem(; height = 1e-3)
    report("gz on the ground vs 1 mm above", maximum(abs, g .- forward(GravForward(m, d1), ρ)) / maximum(abs, g), 1e-4)
    # below the ground is rejected
    _, below, _ = small_problem(; height = -1.0)
    @test_throws ArgumentError GravForward(m, below)
end

@testset "file workflow" begin
    m, d, ρ = small_problem()
    dir = mktempdir()
    op = GravForward(m, d)
    mp = write_model(joinpath(dir, "model.den"), Model(m.dx, m.dy, m.dz, model_array(op, ρ); origin = m.origin))
    dp = write_data_grav(joinpath(dir, "data.dat"), d)
    fp = WriteFwdCtrl3D(joinpath(dir, "FwdCtrl"), FwdCtrl3D(solver = :matrixfree))
    out = ForwardSolve3D(mp, dp, fp)
    @test out == joinpath(dir, "data.pred")
    p = load_data_grav(out)
    @test isapprox(p.value, forward(op, ρ); rtol = 1e-6) && p.error == d.error
    noisy = load_data_grav(ForwardSolve3D(mp, dp, fp; output_path = joinpath(dir, "noisy.pred"), add_noise = true))
    @test !isapprox(noisy.value, p.value; rtol = 1e-6)
end

@testset "datum: the same earth at another z = 0" begin
    # the model moved down by s with its datum raised by s: the same elevations, the same gz
    m, d, ρ = small_problem()
    s = 312.5
    t = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin .+ [0, 0, s], datum = s)
    report("gz with the datum shifted", maximum(abs, forward(GravForward(t, d), ρ) .- forward(GravForward(m, d), ρ)) /
                                        maximum(abs, forward(GravForward(m, d), ρ)), 1e-12)
    @test GravForward(t, d).zs == d.z .+ s
    # depth weights see the same depths below the stations
    for op in (GravForward(m, d), GravForward(t, d))
        @test GG._weights(op, ones(length(ρ)), collect(1:length(ρ)), Inv3DOptions(weighting = :depth)) ≈
              GG._weights(GravForward(m, d), ones(length(ρ)), collect(1:length(ρ)), Inv3DOptions(weighting = :depth))
    end
    # draping and the frame check work in the model's datum, the data keep z = −elevation
    g = drape!(deepcopy(d), t, 0.0)
    @test g.z == drape!(deepcopy(d), m, 0.0).z && check_frame(g, t) === nothing
    @test_throws ArgumentError check_frame(g, Model(t.dx, t.dy, t.dz, t.A; origin = t.origin, datum = s + 1))
end
