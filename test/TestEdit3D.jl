# Model editors, driven headless: layer replacement with a blend, polygon zones, air kept, files written.

@testset "EditModelByLayers" begin
    m = demo_model()                                         # air in layer 1 and two cells of layer 2
    dir = mktempdir()
    p = write_model(joinpath(dir, "start.den"), m); m = load_model(p)          # values as stored
    fig, parts, st = EditModelByLayers(p; target_density = 30.0, blend_previous_percent = 50, interactive = false)
    set_close_to!(parts.slider, 4)                           # cut at layer 4, full view: the whole layer
    st.apply!()
    A = st.A_modified
    @test all(v -> v == 30.0, filter(!isnan, A[:, :, 4:end]))
    @test filter(!isnan, A[:, :, 3]) ≈ filter(!isnan, 0.5 .* m.A[:, :, 3] .+ 0.5 * 30.0)     # linear blend
    @test isnan.(A) == isnan.(m.A)                           # air stays air
    @test A[:, :, 1:2] == m.A[:, :, 1:2] || isequal(A[:, :, 1:2], m.A[:, :, 1:2])
    st.save!()
    out = only(filter(f -> occursin("_modified_layer4_rho30.0_blendprev50_withpadding", f), readdir(dir)))
    r = load_model(joinpath(dir, out))
    @test same_mesh(r, m) && isnan.(r.A) == isnan.(A) && isapprox(filter(!isnan, r.A), filter(!isnan, A); rtol = 1e-6)
end

@testset "EditModelByDrawing" begin
    m = demo_model()
    dir = mktempdir()
    p = write_model(joinpath(dir, "start.den"), m); m = load_model(p)          # values as stored
    fig, parts = EditModelByDrawing(p; replacement_density = 70.0, interactive = false)
    # polygon around the column of cells (3, 2) and (3, 3) in (y, x) plot coordinates
    yl, yh = m.cy[2] - 1, m.cy[3] + 1; xl, xh = m.cx[3] - 1, m.cx[3] + 1
    parts.polygon_points[] = GG.Point2f[(yl, xl), (yh, xl), (yh, xh), (yl, xh)]
    parts.depth_toggle.active[] = true                       # all depths
    notify(parts.apply_button.clicks)
    A = parts.A_working
    @test all(v -> v == 70.0, filter(!isnan, A[3, 2:3, :]))
    @test isequal(A[1:2, :, :], m.A[1:2, :, :]) && isequal(A[4:5, :, :], m.A[4:5, :, :])
    @test isnan.(A) == isnan.(m.A)
    notify(parts.save_button.clicks)
    out = only(filter(f -> occursin("_edited_rho70.0", f), readdir(dir)))
    @test isequal(load_model(joinpath(dir, out)).A, A)
end

@testset "editors keep lake water" begin
    m = demo_model()
    dir = mktempdir()
    W = falses(size(m.A)); W[3, 2:3, 3] .= true                       # two water cells in the edited column
    m.A[W] .= 1000.0 - 2670.0
    p = write_model(joinpath(dir, "start.den"), m); m = load_model(p)
    write_mask(joinpath(dir, "start.mask"), m, active_cells(m) .& .!W; water = W)     # found beside the model
    fig, parts, st = EditModelByLayers(p; target_density = 30.0, interactive = false)
    set_close_to!(parts.slider, 2)
    st.apply!()
    deep = falses(size(m.A)); deep[:, :, 2:end] .= true                # layers 2 and below are replaced
    @test all(st.A_modified[W] .== -1670.0) && all(v -> v == 30.0, filter(!isnan, st.A_modified[deep .& .!W]))
    fig, parts = EditModelByDrawing(p; replacement_density = 70.0, interactive = false)
    yl, yh = m.cy[2] - 1, m.cy[3] + 1; xl, xh = m.cx[3] - 1, m.cx[3] + 1
    parts.polygon_points[] = GG.Point2f[(yl, xl), (yh, xl), (yh, xh), (yl, xh)]
    parts.depth_toggle.active[] = true
    notify(parts.apply_button.clicks)
    A = parts.A_working
    @test all(A[W] .== -1670.0) && all(v -> v == 70.0, filter(!isnan, A[3, 2:3, :][.!W[3, 2:3, :]]))
    # an explicit mask of another mesh is an error; one found beside the model is skipped
    other = Model(m.dx[2:end], m.dy, m.dz; origin = m.origin)
    q = write_mask(joinpath(dir, "other.mask"), other, active_cells(other))
    @test_throws ArgumentError EditModelByLayers(p; mask_file = q, interactive = false)
end
