# 3D mesh design
# Author: @pankajkmishra
# MakeMesh3D on a flat survey and on one over a hill: the core covers the data at the station spacing, the
# padding grows outward, the air mask follows the ground, every sensor ends up in air (check_receivers is the
# judge) and the covariance mask matches the air cells

using Test

@testset "3D mesh design" begin
    mktempdir() do dir
        d = small_survey(n = 9, h = 50.0, height = 20.0)
        data = WriteData3D(joinpath(dir, "data.xyz"), d)

        @test data_spacing(d.x, d.y).median ≈ 50.0

        m = MakeMesh3D(data; out_model = joinpath(dir, "model.start"), out_covariance = joinpath(dir, "C3.dat"),
                       n_pad = 4, verbose = false)
        model = ReadModel3D(joinpath(dir, "model.start"))
        @test m.dx_core == 50.0 && m.nx_core == 8 && m.ny_core == 8
        @test size(model.A) == (16, 16, m.nz) && model.npad == (4, 4)
        @test model.x[1] < minimum(d.x) && model.x[end] > maximum(d.x)
        @test issorted(m.dx[1:4]) == false && m.dx[1] > m.dx[4] > m.dx[5]     # padding grows outward
        @test m.pad_ok && m.depth_ok

        # flat ground at the lowest sensor elevation: no air, every sensor above the mesh
        @test !m.topo_active && m.n_air == 0 && model.top_elevation ≈ 20.0
        @test check_receivers(model, d) === nothing

        cov = ReadCov3D(joinpath(dir, "C3.dat"))
        @test size(cov.mask) == size(model.A) && all(cov.mask .== 1)

        # a hill: air follows the ground, the mesh top is the highest ground, sensors stay in air
        hill(x, y) = 60.0 * exp(-(x^2 + y^2) / (2 * 150.0^2))
        draped = data_like(d, d.value; error = d.error)
        draped.elevation .= [hill(x, y) + 30.0 for (x, y) in zip(d.x, d.y)]
        dpath = WriteData3D(joinpath(dir, "draped.xyz"), draped)
        mt = MakeMesh3D(dpath; out_model = joinpath(dir, "topo.start"), out_covariance = joinpath(dir, "C3topo.dat"),
                        topography = hill, n_pad = 4, verbose = false)
        tmodel = ReadModel3D(joinpath(dir, "topo.start"))
        @test mt.topo_active && mt.n_air > 0
        @test tmodel.top_elevation ≈ mt.ground_max ≈ maximum(hill(x, y) for x in mt.x_centers, y in mt.y_centers)
        @test check_receivers(tmodel, draped) === nothing
        @test count(isnan, tmodel.A) == mt.n_air
        @test ReadCov3D(joinpath(dir, "C3topo.dat")).mask == Int.(.!isnan.(tmodel.A))

        # the deepest air column is the one under the lowest ground
        surface = surface_depth(tmodel)
        @test argmax(filter(isfinite, surface)) > 0 && maximum(filter(isfinite, surface)) > 0
        @test surface[cld(size(surface, 1), 2), cld(size(surface, 2), 2)] ≈ 0.0 atol = mt.z_first
    end
end
