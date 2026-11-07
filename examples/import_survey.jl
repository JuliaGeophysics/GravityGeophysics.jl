using GravityGeophysics

length(ARGS) >= 1 || error("usage: julia --project=. examples/import_survey.jl survey.csv [mt_data.dat [mt_model.rho]]")
survey_file  = ARGS[1]
mt_data_file = length(ARGS) >= 2 ? ARGS[2] : ""
mt_model     = length(ARGS) >= 3 ? ARGS[3] : ""
out_file     = splitext(survey_file)[1] * ".dat"

crs = "EPSG:4326"                  # of easting/northing columns
# edit to the column names of the survey file (lat/lon, easting/northing or north/east); the values are the
# Bouguer anomaly (mGal) at the reduction density the model's contrasts are taken against
columns = (lat = "LAT", lon = "LON", elevation = "ELEV", value = "BOUGUER", error = "ERR", line = "LINE")

frame = isempty(mt_data_file) ? nothing : load_frame(mt_data_file)
d = load_data_grav(survey_file; columns, crs, frame)
if !isempty(mt_model)
    m = load_mesh(mt_model)
    m.rotation == 0 || (d = rotate_data(d, m.rotation))
    check_frame(d, m)
end
println(write_data_grav(out_file, d; comment = "imported from $(basename(survey_file))"))
