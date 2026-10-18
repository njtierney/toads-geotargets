# Climate raster, simple model of toad location, and map
source("packages.R")
tar_source("R")

# Demo 2 workers.
tar_option_set(
  controller = crew_controller_local(workers = 2),
  workspace_on_error = TRUE
)

tar_assign({
  # ---- read, one branch per source ----
  # Every file in the folder. A new file is a new branch.
  toad_files <- dir_ls("data-raw/sources") |> tar_files_input()

  toads <- toad_files |>
    read_toads() |>
    tidy_toads() |>
    tar_target(pattern = map(toad_files))

  # Each source's own front, labelled with its source
  source_fronts <- source_front(toads) |>
    tar_target(pattern = map(toads))

  # ---- analysis, every source together ----
  # Without a pattern, `toads` is every branch bound together
  toads_states <- find_state(toads, map = ozmap_states) |> tar_target()
  toads_west_front <- filter_most_west(toads_states) |> tar_target()
  toad_distances <- add_distance(toads_west_front) |> tar_target()
  front_model <- fit_front(toads_west_front) |> tar_target()

  # ---- results ----
  toads_source_count <- toads_states |>
    count(source, point_in_state) |>
    tar_target()

  toads_decade_count <- toads |> count(decade) |> tar_target()

  gg_decade_count <- tar_target(
    ggplot(toads, aes(x = decade)) + geom_bar()
  )

  gg_distances <- tar_target(
    ggplot(toad_distances, aes(x = decade, y = distance_km)) +
      geom_col() +
      scale_x_continuous(breaks = scales::breaks_width(width = 10))
  )

  # ---- maps ----
  gg_oz_map <- tar_target(
    ozmap_states |>
      filter(NAME %in% unique(toads_states$point_in_state)) |>
      ggplot() +
      geom_sf()
  )

  gg_decades <- gg_toads_decade(gg_oz_map, toads_states) |> tar_target()
  gg_front <- gg_toads_west_path(gg_oz_map, toads_states, toads_west_front) |>
    tar_target()

  # ---- where could they go? ----
  # WorldClim, downloaded once with
  # geodata::worldclim_global(var = "bio", res = 10, path = "data-raw/worldclim")
  bioclim_dir <- "data-raw/worldclim/climate/wc2.1_10m" |> tar_file()

  # terra objects need tar_terra_rast() and tar_terra_vect(), not tar_target()
  bioclim <- read_bioclim(bioclim_dir) |> tar_terra_rast()
  toad_points <- as_toad_points(toads_states) |> tar_terra_vect()
  background_points <- sample_background(bioclim) |> tar_terra_vect()

  model_data <- extract_climate(toad_points, background_points, bioclim) |>
    tar_target()
  sdm <- fit_sdm(model_data) |> tar_target()

  suitability <- predict_suitability(sdm, bioclim) |> tar_terra_rast()
  gg_suitability_map <- gg_suitability(suitability) |> tar_target()

  # ---- outputs ----
  toad_distances_csv <- toad_distances |>
    write_csv_path("output/toad-distances.csv") |>
    tar_file()

  report <- tar_quarto(path = "docs/toad-report.qmd")
}) |>
  # Each {crew} worker is a new R session, so the conflicts_prefer() calls in
  # packages.R don't reach it. This runs them before every target.
  tar_hook_before(
    hook = conflicts_prefer(
      dplyr::filter,
      dplyr::select,
      terra::extract,
      .quiet = TRUE
    )
  )
