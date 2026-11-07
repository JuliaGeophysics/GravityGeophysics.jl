# Masks: round trip, air cells always fixed, mesh and value checks.

@testset "mask" begin
    m = demo_model()
    dir = mktempdir()
    variable = falses(size(m.A)); variable[2:4, 2:3, :] .= true
    path = write_mask(joinpath(dir, "mask.den"), m, variable)
    v = load_mask(path, m)
    @test v == variable .& active_cells(m)
    @test !any(v[air_cells(m)])
    other = Model(m.dx, m.dy, m.dz; origin = m.origin .+ [1.0, 0, 0], rotation = m.rotation)
    @test_throws ArgumentError load_mask(path, other)
    m2 = deepcopy(m); m2.A[active_cells(m)] .= 0.5
    bad = write_model(joinpath(dir, "bad.den"), m2)
    @test_throws ArgumentError load_mask(bad, m)
end

@testset "topography and air masks" begin
    m = demo_model()
    dir = mktempdir()
    air = air_mask_from_model(m)
    @test air == air_cells(m) && any(air)
    topo = extract_topography(m)
    cols = count(any(air[i, j, :]) for i in 1:m.nx, j in 1:m.ny)
    @test length(topo.depth) == cols
    # depth below the model top and z in the model frame describe the same interface
    @test topo.z ≈ topo.depth .+ m.origin[3]
    path = write_topography(joinpath(dir, "topo.txt"), topo; origin = m.origin)
    back = read_topography(path)
    @test isapprox(back.x, topo.x; atol = 1e-3) && isapprox(back.depth, topo.depth; atol = 1e-3)
    @test air_mask_from_topography(m, back) == air
    # stray air below the ground is reported and left out of the surface
    m2 = deepcopy(m); I = findlast(!isnan, m2.A); m2.A[I] = NaN
    t2 = @test_logs (:warn, r"air cells below the ground surface") extract_topography(m2)
    @test t2.depth == topo.depth
end
