# Synthetic models, analytic responses and noise
# Author: @pankajkmishra
# The analytic sphere (a uniform sphere attracts like a point mass outside itself) checks the prism solver;
# prisms and dipping dykes build benchmark models on a mesh
# Notation follows Blakely (1995): γ is the gravitational constant (GAMMA, from Fwd3D.jl), ρ the density contrast

using LinearAlgebra
using Random

#---------- analytic ----------

"""
    sphere_anomaly(data, centre, radius, density_contrast; top_elevation=0.0) -> Vector

Exact gz (mGal) of a uniform sphere at the data points; `centre` is `(x, y, z)` with z the depth below the datum whose elevation is `top_elevation`, as in a model.
"""
function sphere_anomaly(d::GravData3D, centre, radius::Real, Δρ::Real; top_elevation::Real = 0.0)
    M = Δρ * 4π / 3 * radius^3
    map(eachindex(d.x)) do i
        dz = centre[3] - (top_elevation - d.elevation[i])
        R  = sqrt((centre[1] - d.x[i])^2 + (centre[2] - d.y[i])^2 + dz^2)
        R > radius || throw(ArgumentError("datum $i lies inside the sphere"))
        GAMMA * M * dz / R^3 * SI_TO_MGAL
    end
end

#---------- model builders ----------

"""
    add_prism!(model, xrange, yrange, zrange, density_contrast) -> model

Set the earth cells whose centres fall in the box, z as depth below the mesh top, to the density contrast.
"""
function add_prism!(m::GravModel3D, xr, yr, zr, Δρ::Real)
    for k in eachindex(m.cz), j in eachindex(m.cy), i in eachindex(m.cx)
        isnan(m.A[i, j, k]) && continue
        (xr[1] <= m.cx[i] <= xr[2] && yr[1] <= m.cy[j] <= yr[2] && zr[1] <= m.cz[k] <= zr[2]) && (m.A[i, j, k] = Δρ)
    end
    m
end

"""
    add_dyke!(model; top_centre, strike, dip, length, width, depth_extent, density_contrast) -> model

Set the earth cells inside a thin dipping slab: `top_centre` `(x, y, z)` is the middle of its upper edge (z depth below the mesh top), `strike` clockwise from north (deg), `dip` from horizontal (deg) to the right of strike, `width` the true thickness, `depth_extent` the vertical extent.
"""
function add_dyke!(m::GravModel3D; top_centre, strike::Real, dip::Real, length::Real, width::Real,
                   depth_extent::Real, density_contrast::Real)
    s = (cosd(strike), sind(strike))                       # along strike (north, east)
    n = (-sind(strike), cosd(strike))                      # horizontal normal, dip side
    for k in eachindex(m.cz), j in eachindex(m.cy), i in eachindex(m.cx)
        isnan(m.A[i, j, k]) && continue
        dz = m.cz[k] - top_centre[3]
        0 <= dz <= depth_extent || continue
        dxh, dyh = m.cx[i] - top_centre[1], m.cy[j] - top_centre[2]
        along  = dxh * s[1] + dyh * s[2]
        across = dxh * n[1] + dyh * n[2] - dz / tand(dip)  # offset from the dipping mid-plane
        abs(along) <= length / 2 && abs(across * sind(dip)) <= width / 2 && (m.A[i, j, k] = density_contrast)
    end
    m
end

#---------- noise ----------

"""
    add_noise(data; floor=0.01, percent=2.0, seed=20261025) -> GravData3D

Data plus Gaussian noise of standard deviation `floor + percent/100 |datum|` (mGal), which is also written as the error.
"""
function add_noise(d::GravData3D; floor::Real = 0.01, percent::Real = 2.0, seed::Integer = 20261025)
    rng = MersenneTwister(seed)
    σ = floor .+ percent / 100 .* abs.(d.value)
    data_like(d, d.value .+ σ .* randn(rng, length(σ)); error = σ)
end

#---------- misfit ----------

"""
    chi2_rms3d(observed, predicted) -> (chi2, rms, count)

Weighted misfit with the observed errors, `rms = sqrt(chi2 / count)`.
"""
function chi2_rms3d(o::GravData3D, p::GravData3D)
    length(o.value) == length(p.value) || throw(DimensionMismatch("observed and predicted data differ in length"))
    r    = (p.value .- o.value) ./ o.error
    chi2 = sum(abs2, r)
    (; chi2, rms = sqrt(chi2 / length(r)), count = length(r))
end
