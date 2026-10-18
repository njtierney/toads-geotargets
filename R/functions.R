# Read one data source, with a `source` column named from its file:
# data-raw/sources/cane-toad-awc.parquet becomes "awc".
read_toads <- function(path) {
  source_name <- path |>
    path_file() |>
    path_ext_remove() |>
    str_remove("cane-toad-")

  read_parquet(path) |>
    mutate(source = source_name, .before = 1)
}

tidy_toads <- function(toads_raw) {
  toads_raw |>
    clean_names() |>
    rename(
      lat = decimal_latitude,
      lon = decimal_longitude,
      date = event_date,
      coord_var_m = coordinate_uncertainty_in_meters,
      resource_name = data_resource_name
    ) |>
    # Year is computed from date, decade computed from year
    mutate(
      year = year(date),
      decade = floor(year / 10) * 10,
      .after = date
    ) |>
    # Drop records with no date/lat/lon info
    drop_na(date, lon, lat) |>
    # Toads were introduced in 1935, so we filter that out.
    filter(year >= 1935)
}

# Which state each record falls inside, as a new `point_in_state` column.
find_state <- function(toads, map) {
  # `map` is ozmap_states, or anything with a NAME column for each polygon.
  map_states <- map |>
    st_transform(crs = 4326) |>
    select(point_in_state = NAME)

  toads |>
    # Transformed to WGS84 to match ALA. Records outside every state get NA.
    # remove = FALSE keeps the lon/lat
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    st_join(map_states) |>
    st_drop_geometry()
}

# The front in each decade: the mean position of the `n` westernmost records.
# One record can be an error; the mean of ten is steadier. summarise() keeps
# only decade, lon and lat, because an average of ten records has no single
# state.
filter_most_west <- function(toads, n = 10) {
  toads |>
    group_by(decade) |>
    slice_min(lon, n = n) |>
    summarise(
      lon = mean(lon),
      lat = mean(lat)
    ) |>
    ungroup()
}

# The front for one source, labelled with that source, so the fronts can be
# compared once the branches are bound together.
source_front <- function(toads, n = 10) {
  filter_most_west(toads, n = n) |>
    mutate(source = first(toads$source), .before = 1)
}

# A straight line through the front: how many degrees of longitude it moves
# each year. The front has a decade, not a year, so the slope is per year of
# decade.
fit_front <- function(front) {
  lm(lon ~ decade, data = front)
}

# How far each row is from the row before it, in kilometres.
add_distance <- function(data, lon = "lon", lat = "lat") {
  # geodist() calculates distance between lon/lat pairs
  dist_mat <- cbind(lon = data[[lon]], lat = data[[lat]])

  distances_m <- dist_mat |>
    geodist(
      measure = "geodesic",
      # `sequential = TRUE` makes it row-to-row rather than pairwise.
      sequential = TRUE,
      # No step into the first row, so it is padded with NA.
      pad = TRUE
    )

  data |>
    mutate(
      distance_km = distances_m / 1000,
      .after = decade
    )
}

# Every record on a map, one panel per decade. `map` is a ggplot of the states.
gg_toads_decade <- function(map, toads) {
  map +
    geom_point(
      data = toads,
      aes(x = lon, y = lat),
      alpha = 0.2
    ) +
    facet_wrap(
      facets = vars(decade),
      nrow = 2
    )
}

# "Front" is the westernmost record in each decade, in orange
# `map` is a ggplot of the states.
gg_toads_west_path <- function(map, toads, toads_west) {
  map +
    geom_point(
      data = toads,
      aes(x = lon, y = lat),
      alpha = 0.2
    ) +
    geom_point(
      data = toads_west,
      aes(x = lon, y = lat),
      colour = "orange"
    ) +
    geom_path(
      data = toads_west,
      aes(x = lon, y = lat),
      colour = "orange"
    )
}

# Write a data frame to CSV and return the path, so the file can be a target.
# readr::write_csv() returns the data, not the path, which is why the last line
# is needed.
write_csv_path <- function(x, path) {
  dir_create(path_dir(path))
  write_csv(x, path)
  path
}

# ---- geospatial ----

# Two WorldClim layers, cropped to Australia: bio_1 is annual mean temperature
# (degrees C), bio_12 is annual precipitation (mm).
read_bioclim <- function(
  dir,
  layers = c(1, 12),
  extent = c(112, 154, -44, -9)
) {
  paths <- path(dir, glue("wc2.1_10m_bio_{layers}.tif"))

  climate <- rast(paths) |>
    crop(ext(extent))

  names(climate) <- glue("bio_{layers}")
  climate
}

# The records as terra points. Records outside every state are at sea, or
# overseas, so they go.
as_toad_points <- function(toads_states) {
  toads_states |>
    drop_na(point_in_state) |>
    vect(geom = c("lon", "lat"), crs = "EPSG:4326")
}

# Random points on land, to compare the records against: where toads were
# recorded, against everywhere they could have been.
sample_background <- function(climate, n = 5000) {
  spatSample(
    climate,
    size = n,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE
  )
}

# The climate at each record (presence = 1) and each background point
# (presence = 0), as one data frame.
extract_climate <- function(points, background, climate) {
  presences <- extract(climate, points, ID = FALSE) |>
    mutate(presence = 1)

  absences <- extract(climate, background, ID = FALSE) |>
    mutate(presence = 0)

  bind_rows(presences, absences) |>
    drop_na()
}

# Your model goes here. A logistic regression of records against background,
# on temperature and rainfall.
fit_sdm <- function(model_data) {
  glm(
    presence ~ bio_1 + bio_12,
    family = binomial,
    data = model_data
  )
}

# The model's prediction for every cell of the climate raster.
predict_suitability <- function(sdm, climate) {
  predict(climate, sdm, type = "response")
}

# The prediction as a map. The raster becomes a data frame first, so the plot
# doesn't hold on to a terra object.
gg_suitability <- function(suitability) {
  suitability_df <- as.data.frame(suitability, xy = TRUE)
  names(suitability_df) <- c("lon", "lat", "suitability")

  ggplot(suitability_df, aes(x = lon, y = lat, fill = suitability)) +
    geom_raster() +
    scale_fill_viridis_c() +
    coord_sf(crs = 4326)
}
