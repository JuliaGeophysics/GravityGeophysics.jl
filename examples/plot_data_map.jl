using GravityGeophysics

flags      = filter(a -> startswith(a, "--"), ARGS)
positional = filter(a -> !startswith(a, "--"), ARGS)

data_file      = length(positional) >= 1 ? positional[1] : ""
predicted_file = length(positional) >= 2 ? positional[2] : ""
coordinate_system = length(positional) >= 3 ? positional[3] : "model"

# extra overlays: (path, enabled, color, alpha, point_size, line_width); more can be added from the viewer
shapefiles = []
shapefile_color      = :grey30    # style of the shapefiles added from the viewer
shapefile_alpha      = 0.9
shapefile_line_width = 1.2
shapefile_point_size = 6

colormap       = COLORMAP
value_range    = nothing          # observed / predicted colour range; nothing: 2-98 %
residual_range = (-3.0, 3.0)
site_markersize = 9
map_pad        = 0.05

viewer_figsize = (1100, 950)
export_dpi     = 3
export_figsize = (1100, 900)

gis_output_dir = ""               # default <data>-GIS next to the data file

if "--gis-only" in flags
    write_data_gis(data_file, predicted_file;
        crs        = coordinate_system == "model" ? "model" : coordinate_system,
        output_dir = gis_output_dir)
    exit(0)
end

PlotDataMap(data_file, predicted_file;
    crs            = coordinate_system,
    shapefiles     = shapefiles,
    shapefile_color      = shapefile_color,
    shapefile_alpha      = shapefile_alpha,
    shapefile_line_width = shapefile_line_width,
    shapefile_point_size = shapefile_point_size,
    colormap       = colormap,
    value_range    = value_range,
    residual_range = residual_range,
    site_markersize = site_markersize,
    map_pad        = map_pad,
    viewer_figsize = viewer_figsize,
    export_dpi     = export_dpi,
    export_figsize = export_figsize,
    gis_output_dir = gis_output_dir,
)
