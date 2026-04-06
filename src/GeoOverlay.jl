function _prj_path_from_shp(shp::AbstractString)
    root, _ = splitext(shp)
    return root * ".prj"
end

function _detect_crs_type(shp::AbstractString)
    prj_path = _prj_path_from_shp(shp)
    if !isfile(prj_path)
        return :unknown, ""
    end
    wkt = read(prj_path, String)
    upper = uppercase(wkt)
    if occursin("PROJCS", upper) || occursin("PROJCRS", upper)
        return :projected, wkt
    elseif occursin("GEOGCS", upper) || occursin("GEOGRAPHICCRS", upper)
        return :geographic, wkt
    end
    return :unknown, wkt
end

function _make_coord_transform(shp_path::AbstractString, crs_type::Symbol; auto_reproject_to_wgs84::Bool = false)
    if !(auto_reproject_to_wgs84 && crs_type == :projected)
        return (x, y) -> (Float64(x), Float64(y))
    end
    prj_path = _prj_path_from_shp(shp_path)
    !isfile(prj_path) && return (x, y) -> (Float64(x), Float64(y))
    try
        trans = Proj.Transformation(read(prj_path, String), "EPSG:4326"; always_xy = true)
        return (x, y) -> begin
            ll = trans((Float64(x), Float64(y)))
            return Float64(ll[1]), Float64(ll[2])
        end
    catch
        return (x, y) -> (Float64(x), Float64(y))
    end
end

function shapefile_segments(shapefile_path::AbstractString; auto_reproject_to_wgs84::Bool = false, post_transform = (x, y) -> (x, y), xlim = nothing, ylim = nothing)
    !isfile(shapefile_path) && return Vector{Vector{Tuple{Float64, Float64}}}()
    table = Shapefile.Table(shapefile_path)
    crs_type, _ = _detect_crs_type(shapefile_path)
    coord_transform = _make_coord_transform(shapefile_path, crs_type; auto_reproject_to_wgs84 = auto_reproject_to_wgs84)
    segments = Vector{Vector{Tuple{Float64, Float64}}}()
    xlo, xhi = isnothing(xlim) ? (-Inf, Inf) : extrema(Float64.(collect(xlim)))
    ylo, yhi = isnothing(ylim) ? (-Inf, Inf) : extrema(Float64.(collect(ylim)))

    is_xy(c) = (c isa Tuple || c isa AbstractVector) && length(c) >= 2 && c[1] isa Real && c[2] isa Real
    inside(x, y) = (x >= xlo && x <= xhi && y >= ylo && y <= yhi)

    function collect_recursive!(coords)
        if coords isa AbstractVector && !isempty(coords)
            first_item = first(coords)
            if is_xy(first_item)
                segment = Tuple{Float64, Float64}[]
                for point in coords
                    x0, y0 = coord_transform(Float64(point[1]), Float64(point[2]))
                    x, y = post_transform(x0, y0)
                    inside(x, y) && push!(segment, (x, y))
                end
                length(segment) >= 2 && push!(segments, segment)
            else
                for child in coords
                    collect_recursive!(child)
                end
            end
        end
    end

    for row in table
        geom = GeoInterface.geometry(row)
        coords = try
            GeoInterface.coordinates(geom)
        catch
            nothing
        end
        isnothing(coords) || collect_recursive!(coords)
    end
    return segments
end

function plot_shapefile_on_3d!(ax, shapefile_path::AbstractString; z_fixed::Real = 0.0, line_color = :black, line_width::Real = 1.3, auto_reproject_to_wgs84::Bool = false, post_transform = (x, y) -> (x, y), xlim = nothing, ylim = nothing)
    count = 0
    for segment in shapefile_segments(shapefile_path; auto_reproject_to_wgs84 = auto_reproject_to_wgs84, post_transform = post_transform, xlim = xlim, ylim = ylim)
        GLMakie.lines!(ax, first.(segment), last.(segment), fill(Float64(z_fixed), length(segment)); color = line_color, linewidth = line_width)
        count += 1
    end
    return count
end

function plot_shapefile_on_axis!(ax, shapefile_path::AbstractString; line_color = :black, line_width::Real = 1.3, auto_reproject_to_wgs84::Bool = false, post_transform = (x, y) -> (x, y), xlim = nothing, ylim = nothing)
    count = 0
    for segment in shapefile_segments(shapefile_path; auto_reproject_to_wgs84 = auto_reproject_to_wgs84, post_transform = post_transform, xlim = xlim, ylim = ylim)
        GLMakie.lines!(ax, first.(segment), last.(segment); color = line_color, linewidth = line_width)
        count += 1
    end
    return count
end