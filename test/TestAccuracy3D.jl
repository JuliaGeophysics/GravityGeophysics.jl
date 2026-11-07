# Solver accuracy. Each check states its tolerance and prints the measured error.

using HCubature
include(joinpath(@__DIR__, "FwdFD3D.jl"))

relmax(a, b) = maximum(abs, a .- b) / maximum(abs, b)

# gz (mGal, positive down) of a point mass M at R = source − station (z down)
point_gz(M, R) = GAMMA * SI_TO_MGAL * M * R[3] / norm(R)^3

# 10 m cube, Δρ = 300 kg/m³, at 100–110 m depth, in a mesh padded with empty cells
cube_mesh() = Model([1000.0, 10, 1000], [1000.0, 10, 1000], [100.0, 10, 100]; origin = [-1005.0, -1005.0, 0.0])
cube_ρ() = (ρ = zeros(27); ρ[LinearIndices((3, 3, 3))[2, 2, 2]] = 300.0; ρ)

@testset "4.1 cube vs point mass (far field)" begin
    m, ρ = cube_mesh(), cube_ρ()
    x = [400.0, -300, 50, 0]; y = [100.0, 250, -500, 0]; z = [-50.0, -20, -80, -350]
    d = station_data(x, y, z)
    ref = [point_gz(1000 * 300.0, [0 - x[k], 0 - y[k], 105 - z[k]]) for k in eachindex(x)]
    # a cube has no quadrupole: the difference is the l = 4 term, (a/r)⁴ ≈ 1e-7 here
    report("cube vs point mass", relmax(forward(GravForward(m, d), ρ), ref), 1e-6)
end

@testset "4.2 near field vs numerical integration" begin
    # prism column sticks up out of the topography: air down to 110 m beside it
    m = Model([100.0, 10, 100], [100.0, 10, 100], [100.0, 10, 100]; origin = [-100.0, -100.0, 0.0])
    m.A[:, :, 1:2] .= NaN; m.A[2, 2, 2] = 0.0
    m = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
    cells = findall(active_cells(m)); ρ = [I == CartesianIndex(2, 2, 2) ? 300.0 : 0.0 for I in cells]
    stations = [(5.0, 5.0, 99.0, "1 m above the centre"), (0.0, 5.0, 99.0, "1 m above an edge"),
                (0.0, 0.0, 99.0, "1 m above a corner"), (-1.0, 5.0, 102.0, "1 m beside, below the top"),
                (12.0, 13.0, 95.0, "off a corner, above")]
    for (x, y, z, label) in stations
        d = station_data([x], [y], [z])
        ref = hcubature(p -> point_gz(300.0, p .- [x, y, z]), [0.0, 0, 100], [10.0, 10, 110]; rtol = 1e-11, maxevals = 10^7)[1]
        got = forward(GravForward(m, d), ρ)[1]
        report("gz, $label", abs(got - ref) / abs(ref), 1e-8)
    end
end

@testset "4.3 superposition, Bouguer slab and limits" begin
    one = Model([1000.0, 10, 1000], [1000.0, 10, 1000], [100.0, 10, 100]; origin = [-1005.0, -1005.0, 0.0])
    split = Model([1000.0, 5, 5, 1000], [1000.0, 5, 5, 1000], [100.0, 5, 5, 100]; origin = [-1005.0, -1005.0, 0.0])
    ρs = [I[1] in 2:3 && I[2] in 2:3 && I[3] in 2:3 ? 300.0 : 0.0 for I in CartesianIndices((4, 4, 4))]
    x = [0.0, 30, -7]; y = [0.0, -12, 9]; z = [-1.0, -40, -3]
    d = station_data(x, y, z)
    report("one prism = 8 sub-prisms", relmax(forward(GravForward(split, d), vec(ρs)), forward(GravForward(one, d), cube_ρ())), 1e-12)
    report("zero model gives zero data", maximum(abs, forward(GravForward(one, d), zeros(27))), 0.0)
    # a 10 m slab 200 km wide: 2πγΔρt, less the edge deficit of about 2t/(πa) ≈ 6e-5
    slab = Model([1e5, 1e5], [1e5, 1e5], [10.0]; origin = [-1e5, -1e5, 0.0])
    g = forward(GravForward(slab, station_data([0.0], [0.0], [-1.0])), fill(300.0, 4))[1]
    report("Bouguer slab, 1 m above", abs(g - 2π * GAMMA * 300 * 10 * SI_TO_MGAL) / g, 2e-4)
    # on a vertex, an edge and a face of the cube: finite, and continuous with a station 1 mm above
    for (x, y, label) in ((-5.0, -5.0, "vertex"), (0.0, -5.0, "edge"), (0.0, 0.0, "face"))
        m = Model([1000.0, 10, 1000], [1000.0, 10, 1000], [100.0, 10, 100]; origin = [-1005.0, -1005.0, -100.0])
        m.A[:, :, 1] .= NaN; m = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
        ρ = [I == CartesianIndex(2, 2, 2) ? 300.0 : 0.0 for I in findall(active_cells(m))]
        on = forward(GravForward(m, station_data([x], [y], [0.0])), ρ)[1]
        above = forward(GravForward(m, station_data([x], [y], [-1e-3])), ρ)[1]
        report("station on a cube $label vs 1 mm above", abs(on - above) / abs(on), 1e-3)
    end
end

@testset "4.4 round-off: far stations and thin cells vs BigFloat" begin
    setprecision(BigFloat, 256) do
        function prism(T, st, lo, hi)
            x0, y0, z0 = T.(st); s = zero(T)
            for a in 0:1, b in 0:1, c in 0:1
                X = T(a == 1 ? hi[1] : lo[1]) - x0; Y = T(b == 1 ? hi[2] : lo[2]) - y0; Z = T(c == 1 ? hi[3] : lo[3]) - z0
                r = sqrt(X^2 + Y^2 + Z^2)
                f = X * GG._log_plus(Y, X^2 + Z^2, r) + Y * GG._log_plus(X, Y^2 + Z^2, r) - Z * GG._atan_ratio(X * Y, Z * r)
                s += (-1)^(a + b + c + 1) * f
            end
            s
        end
        # the vertex sum cancels as the station moves away and the cells thin; the error is measured
        # against the response right above the cell, the scale of the data
        for (lo, hi, label, tol) in (((-5.0, -5.0, 100.0), (5.0, 5.0, 110.0), "10 m cell", 1e-7),
                                     ((-5.0, -5.0, 100.0), (5.0, 5.0, 101.0), "1 m thick cell", 1e-6),
                                     ((-5.0, -5.0, 100.0), (5.0, 5.0, 100.01), "1 cm thick cell (stress)", 1e-4))
            near = abs(prism(BigFloat, (0.0, 0.0, -100.0), lo, hi))
            for R in (5_000.0, 50_000.0, 200_000.0)
                st = (0.6R, 0.8R, -100.0)
                err = Float64(abs(prism(Float64, st, lo, hi) - prism(BigFloat, st, lo, hi)) / near)
                report("$label at $(Int(R ÷ 1000)) km, err / response above", err, tol)
            end
        end
    end
end

@testset "4.6 finite-difference Poisson solver" begin
    # a sphere of 400 kg/m³ on a 20 m mesh, stations 40 m above it; the FD solver discretises the same
    # cells, so the two solvers differ only by the FD truncation and boundary errors
    h, n, nz = 20.0, 24, 12
    A = zeros(n, n, nz)
    m = Model(fill(h, n), fill(h, n), fill(h, nz), A; origin = [-n * h / 2, -n * h / 2, 0.0])
    for k in 1:nz, j in 1:n, i in 1:n
        hypot(m.cx[i], m.cy[j], m.cz[k] - 120) <= 80 && (m.A[i, j, k] = 400.0)
    end
    c = ((1:7) .- 4) .* 60.0
    d = station_data(repeat(c; inner = 7), repeat(c; outer = 7), fill(-40.0, 49))
    op = GravForward(m, d)
    prism = forward(op, model_vector(op, m))
    fd = forward3d_fd(m, d, FDCtrl3D(tol = 1e-11, maxiter = 20_000))
    relrms(a, b) = sqrt(mean(abs2, a .- b)) / sqrt(mean(abs2, b))
    report("FD Poisson vs prisms, relative rms", relrms(fd, prism), 0.05)
    @test cor(fd, prism) > 0.99
    # the staircase sphere holds about the sphere's mass: point-mass field to a few per cent
    M = 400.0 * 4π / 3 * 80^3
    exact = [point_gz(M, [-d.x[k], -d.y[k], 120 - d.z[k]]) for k in 1:ndata(d)]
    report("prisms vs the sphere (staircase)", relmax(prism, exact), 0.05)
end

@testset "4.9 rotation: the same body on a 90°-turned frame" begin
    # mesh B is mesh A seen from a frame turned 90° clockwise: x' = y, y' = −x
    A = Model([100.0, 20, 30, 20, 100], [100.0, 10, 20, 40, 100], [50.0, 20, 40]; origin = [-135.0, -135.0, 0.0])
    ρA = [10.0 * (I[1] + 2I[2] + 3I[3]) for I in CartesianIndices(size(A.A))]
    B = Model(A.dy, reverse(A.dx), A.dz, [ρA[A.nx + 1 - j, i, k] for i in 1:A.ny, j in 1:A.nx, k in 1:A.nz];
              origin = [A.origin[2], -(A.origin[1] + sum(A.dx)), A.origin[3]], rotation = 90.0)
    d = station_data([10.0, -40, 33], [5.0, 20, -27], [-5.0, -15, -30])
    a = forward(GravForward(A, d), vec(ρA))
    b = forward(GravForward(B, rotate_data(d, 90.0)), vec(B.A))
    report("rotated frame", relmax(b, a), 1e-12)
end

@testset "4.10 air cells above the earth" begin
    flat = Model([200.0, 50, 50, 200], [200.0, 50, 50, 200], [30.0, 60]; origin = [-250.0, -250.0, 0.0])
    ρ = [10.0 * (I[1] + I[2] * I[3]) for I in CartesianIndices((4, 4, 2))]
    air = Model(flat.dx, flat.dy, [40.0, 30, 60], cat(fill(NaN, 4, 4, 1), ρ; dims = 3); origin = [-250.0, -250.0, -40.0])
    d = station_data([10.0, -60], [5.0, 70], [-45.0, -50])
    report("air layer on top changes nothing", relmax(forward(GravForward(air, d), vec(ρ)), forward(GravForward(flat, d), vec(ρ))), 1e-13)
end
