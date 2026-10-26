# Model and data files
# Author: @pankajkmishra
# Round trips of the density model (air as 1e17, the top-elevation header) and of the .xyz
# survey, the sensitivity of readers to malformed files, and the geographic-to-UTM projection on reading

using Test

@testset "model and data I/O" begin
    mktempdir() do dir
        m = small_mesh(n = 5, nz = 4, h = 30.0)
        m.A[2:3, 2:4, 2] .= 275.5
        m.A[:, :, 1] .= NaN                      # a whole air layer
        m.top_elevation = 143.0
        path = WriteModel3D(joinpath(dir, "model.den"), m; title = "round trip")
        back = ReadModel3D(path)
        @test isequal(back.A, m.A) && back.dx == m.dx && back.dy == m.dy && back.dz == m.dz
        @test back.origin == m.origin && back.top_elevation == 143.0
        @test same_mesh(back, m) && air_mask(back) == air_mask(m)
        @test count(isnan, back.A) == 25
        @test cell_volumes(back)[1, 1, 1] ≈ 30.0^3
        @test all(surface_depth(back)[:] .≈ 30.0)      # the first earth layer starts one cell down

        # air is written as the file sentinel, not as NaN text
        text = read(path, String)
        @test occursin("1e+17", lowercase(text)) || occursin("1.000000e+17", lowercase(text))
        @test occursin("Top elevation", text)

        # a file in another block layout is refused
        open(joinpath(dir, "other.den"), "w") do io
            println(io, "5 5 4")
            println(io, "0.0")
        end
        @test_throws Exception ReadModel3D(joinpath(dir, "other.den"))
        @test_throws ErrorException ReadModel3D(joinpath(dir, "absent.den"))

        # data round trip, errors included
        d = small_survey(n = 4, h = 25.0, height = 12.0)
        d.value .= range(-1.5, 2.5; length = ndata(d))
        d.error .= 0.05
        dpath = WriteData3D(joinpath(dir, "data.xyz"), d; title = "survey")
        bd = ReadData3D(dpath)
        @test bd.x ≈ d.x && bd.y ≈ d.y && bd.elevation ≈ d.elevation
        @test bd.value ≈ d.value rtol = 1e-5         # values are written with six significant digits
        @test bd.error ≈ d.error && bd.title == "survey"
        @test ndata(bd) == 16 && sensor_depth(m, bd) ≈ m.top_elevation .- d.elevation

        # data_like keeps the geometry and replaces the values
        p = data_like(d, 2 .* d.value; title = "predicted")
        @test p.x == d.x && p.value ≈ 2 .* d.value && p.title == "predicted"
        @test chi2_rms3d(d, d).rms == 0.0
        @test chi2_rms3d(d, data_like(d, d.value .+ d.error)).rms ≈ 1.0

        # geographic coordinates are projected to UTM on reading
        lines = ["# Title : geographic", "# CRS : EPSG:4326",
                 "# longitude latitude elevation value error"]
        for (lon, lat) in ((140.50, -21.75), (140.51, -21.75), (140.50, -21.74), (140.51, -21.74))
            push!(lines, "$lon $lat 0.0 1.0 0.1")
        end
        gpath = joinpath(dir, "geo.xyz")
        write(gpath, join(lines, "\n") * "\n")
        g = ReadData3D(gpath)
        @test ndata(g) == 4 && occursin("327", g.crs)          # UTM zone 54 south
        @test all(>(1e5), g.y) && all(<(0), g.x) == false      # eastings in metres, northings positive
        @test maximum(g.y) - minimum(g.y) > 500                # about 1 km across
    end
end
