library(sf)
library(terra)
library(stars)
library(mapview)

aoi <- st_read("data/polygons.gpkg", layer = "study_area", quiet = TRUE)
zones <- st_read("data/polygons.gpkg", layer = "zones", quiet = TRUE)
terra_aoi <- terra::vect("data/polygons.gpkg", layer = "study_area")
terra_zones <- terra::vect("data/polygons.gpkg", layer = "zones")

points_tbl <- read.csv("data/points.csv")
points_sf <- st_as_sf(points_tbl, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
terra_points <- terra::vect(points_tbl, geom = c("lon", "lat"), crs = "EPSG:4326")

tracks_sf <- points_sf |>
  dplyr::arrange(track_id, timestamp) |>
  dplyr::group_by(track_id) |>
  dplyr::summarise(do_union = FALSE) |>
  st_cast("LINESTRING")

tracks_tbl <- points_tbl |>
  dplyr::arrange(track_id, timestamp) |>
  dplyr::summarise(
    wkt = paste0(
      "LINESTRING (",
      paste(paste(lon, lat), collapse = ", "),
      ")"
    ),
    .by = track_id
  )

terra_tracks <- terra::vect(tracks_tbl, geom = "wkt", crs = "EPSG:4326")

aoi_qc <- st_transform(aoi, 32198)
zones_qc <- st_transform(zones, 32198)
points_sf_qc <- st_transform(points_sf, 32198)
tracks_sf_qc <- st_transform(tracks_sf, 32198)

terra_aoi_qc <- terra::project(terra_aoi, "EPSG:32198")
terra_zones_qc <- terra::project(terra_zones, "EPSG:32198")
terra_points_qc <- terra::project(terra_points, "EPSG:32198")
terra_tracks_qc <- terra::project(terra_tracks, "EPSG:32198")

focus_area <- zones |>
  dplyr::filter(zone_id %in% c("NE", "SE")) |>
  dplyr::summarise()

points_focus_sf <- st_filter(points_sf, focus_area)
zones_focus_sf <- st_crop(zones, st_bbox(focus_area))

terra_focus_area <- terra_zones[terra_zones$zone_id %in% c("NE", "SE")]
terra_points_focus <- crop(terra_points, terra_focus_area)
terra_zones_focus <- crop(terra_zones, ext(terra_focus_area))

zone_area_sf <- data.frame(
  zone_id = zones_qc$zone_id,
  area_km2 = round(as.numeric(units::set_units(st_area(zones_qc), "km^2")), 1)
)

track_length_sf <- data.frame(
  track_id = tracks_sf_qc$track_id,
  length_km = round(as.numeric(units::set_units(st_length(tracks_sf_qc), "km")), 1)
)

point_boundary_sf <- data.frame(
  obs_id = points_sf_qc$obs_id,
  boundary_km = round(as.numeric(units::set_units(st_distance(points_sf_qc, st_boundary(aoi_qc)), "km")), 1)
)

zone_area_terra <- data.frame(
  zone_id = terra_zones_qc$zone_id,
  area_km2 = round(expanse(terra_zones_qc, unit = "km"), 1)
)

track_length_terra <- data.frame(
  track_id = terra_tracks_qc$track_id,
  length_km = round(perim(terra_tracks_qc) / 1000, 1)
)

point_boundary_terra <- data.frame(
  obs_id = terra_points_qc$obs_id,
  boundary_km = round(distance(terra_points_qc, as.lines(terra_aoi_qc))[, 1] / 1000, 1)
)

points_joined_sf <- st_join(points_sf, zones[, c("zone_id", "zone_type")])
terra_point_attrs <- terra::extract(terra_zones, terra_points)
terra_points_joined <- cbind(terra_points, terra_point_attrs[, c("zone_id", "zone_type")])

track_segments_sf <- suppressWarnings(st_intersection(tracks_sf_qc, zones_qc))
terra_track_segments <- intersect(terra_tracks_qc, terra_zones_qc)

point_buffers_sf <- st_buffer(points_sf_qc, dist = 10000)
terra_point_buffers <- buffer(terra_points_qc, width = 10000)

surface_stars <- read_stars("data/surface.tif")
surface_terra <- rast("data/surface.tif")
surface_stars_qc <- st_warp(surface_stars, crs = st_crs(aoi_qc))
surface_stars_crop <- st_crop(surface_stars_qc, st_bbox(aoi_qc))
surface_stars_mask <- surface_stars_crop[aoi_qc]
surface_terra_qc <- project(surface_terra, "EPSG:32198")
surface_terra_crop <- crop(surface_terra_qc, ext(terra_aoi_qc))
surface_terra_mask <- mask(surface_terra_crop, terra_aoi_qc)

crop_zone_sf <- zones[zones$zone_id == "NE", ]
crop_zone_terra <- terra_zones[terra_zones$zone_id == "NE"]
single_buffer_sf <- st_transform(point_buffers_sf[2, ], st_crs(surface_stars))
single_buffer_terra <- project(terra_point_buffers[2], crs(surface_terra))

surface_stars_crop_demo <- st_crop(surface_stars, st_bbox(crop_zone_sf))
surface_stars_mask_demo <- surface_stars[single_buffer_sf]

surface_terra_crop_demo <- crop(surface_terra, crop_zone_terra)
surface_terra_mask_demo <- mask(surface_terra, single_buffer_terra)

surface_stars_template <- st_as_stars(st_bbox(surface_stars), dx = 0.05, dy = 0.05)
st_crs(surface_stars_template) <- st_crs(surface_stars)
surface_stars_near <- suppressWarnings(
  st_warp(surface_stars, dest = surface_stars_template, method = "near", use_gdal = TRUE)
)
surface_stars_bilinear <- suppressWarnings(
  st_warp(surface_stars, dest = surface_stars_template, method = "bilinear", use_gdal = TRUE)
)

surface_terra_template <- rast(ext(surface_terra), resolution = 0.05, crs = crs(surface_terra))
surface_terra_near <- resample(surface_terra, surface_terra_template, method = "near")
surface_terra_bilinear <- resample(surface_terra, surface_terra_template, method = "bilinear")

surface_stars_alt <- surface_stars * 1.1
names(surface_stars_alt) <- "surface_alt"
surface_stars_stack <- c(surface_stars, surface_stars_alt)
surface_stars_mean <- st_apply(surface_stars_stack, MARGIN = c("x", "y"), FUN = mean)

surface_terra_alt <- surface_terra * 1.1
names(surface_terra_alt) <- "surface_alt"
surface_terra_stack <- c(surface_terra, surface_terra_alt)
surface_terra_mean <- app(surface_terra_stack, mean)

points_surface_vals_sf <- st_extract(surface_stars, points_joined_sf)
points_surface_sf <- points_joined_sf
points_surface_sf$surface_value <- points_surface_vals_sf[[1]]

tracks_surface_vals_sf <- st_extract(
  surface_stars,
  tracks_sf,
  FUN = function(x) mean(x, na.rm = TRUE)
)
tracks_surface_sf <- tracks_sf
tracks_surface_sf$surface_value <- as.numeric(tracks_surface_vals_sf[[1]])

zones_surface_vals_sf <- st_extract(
  surface_stars,
  zones,
  FUN = function(x) mean(x, na.rm = TRUE)
)
zones_surface_sf <- zones
zones_surface_sf$surface_value <- as.numeric(zones_surface_vals_sf[[1]])

point_zone_summary_sf <- points_surface_sf |>
  st_drop_geometry() |>
  dplyr::group_by(zone_id, zone_type) |>
  dplyr::summarise(
    point_n = dplyr::n(),
    mean_surface = round(mean(surface_value, na.rm = TRUE), 2),
    mean_value = round(mean(value, na.rm = TRUE), 2),
    .groups = "drop"
  )

analysis_table_sf <- point_zone_summary_sf |>
  dplyr::left_join(zone_area_sf, by = "zone_id") |>
  dplyr::select(zone_id, zone_type, point_n, mean_surface, mean_value, area_km2)

terra_points_surface_vals <- terra::extract(surface_terra, terra_points_joined)
terra_points_surface <- cbind(
  terra_points_joined,
  terra_points_surface_vals[, "surface_value", drop = FALSE]
)

terra_tracks_surface_vals <- terra::extract(surface_terra, terra_tracks, fun = mean, na.rm = TRUE)
terra_tracks_surface <- cbind(
  terra_tracks,
  terra_tracks_surface_vals[, "surface_value", drop = FALSE]
)

terra_zones_surface_vals <- terra::extract(surface_terra, terra_zones, fun = mean, na.rm = TRUE)
terra_zones_surface <- cbind(
  terra_zones,
  terra_zones_surface_vals[, "surface_value", drop = FALSE]
)

point_zone_summary_terra <- as.data.frame(terra_points_surface) |>
  dplyr::group_by(zone_id, zone_type) |>
  dplyr::summarise(
    point_n = dplyr::n(),
    mean_surface = round(mean(surface_value, na.rm = TRUE), 2),
    mean_value = round(mean(value, na.rm = TRUE), 2),
    .groups = "drop"
  )

analysis_table_terra <- point_zone_summary_terra |>
  dplyr::left_join(zone_area_terra, by = "zone_id") |>
  dplyr::select(zone_id, zone_type, point_n, mean_surface, mean_value, area_km2)

analysis_data_sf <- points_surface_sf |>
  st_drop_geometry() |>
  dplyr::select(obs_id, track_id, zone_id, zone_type, category, value, surface_value) |>
  dplyr::filter(!is.na(zone_id), !is.na(surface_value)) |>
  dplyr::mutate(zone_type = factor(zone_type))

analysis_data_terra <- as.data.frame(terra_points_surface) |>
  dplyr::select(obs_id, track_id, zone_id, zone_type, category, value, surface_value) |>
  dplyr::filter(!is.na(zone_id), !is.na(surface_value)) |>
  dplyr::mutate(zone_type = factor(zone_type))

glm_sf <- glm(value ~ surface_value + zone_type, family = poisson(), data = analysis_data_sf)
glm_terra <- glm(value ~ surface_value + zone_type, family = poisson(), data = analysis_data_terra)

glm_coef_sf <- data.frame(
  term = rownames(summary(glm_sf)$coefficients),
  round(as.data.frame(summary(glm_sf)$coefficients), 3),
  row.names = NULL
)

glm_coef_terra <- data.frame(
  term = rownames(summary(glm_terra)$coefficients),
  round(as.data.frame(summary(glm_terra)$coefficients), 3),
  row.names = NULL
)

set.seed(42)
pseudo_points_sf <- st_as_sf(st_sample(aoi_qc, size = nrow(points_sf_qc), type = "random", exact = TRUE))
pseudo_points_sf$pseudo_id <- paste0("PA", sprintf("%02d", seq_len(nrow(pseudo_points_sf))))
pseudo_points_sf <- st_join(pseudo_points_sf, zones_qc[, c("zone_id", "zone_type")])
pseudo_surface_vals_sf <- st_extract(surface_stars_qc, pseudo_points_sf)
pseudo_points_sf$surface_value <- pseudo_surface_vals_sf[[1]]
pseudo_points_sf$presence <- 0L

presence_points_sf <- st_transform(points_surface_sf, 32198)
presence_points_sf$presence <- 1L

pa_analysis_sf <- dplyr::bind_rows(
  presence_points_sf |>
    dplyr::select(zone_id, zone_type, surface_value, presence),
  pseudo_points_sf |>
    dplyr::select(zone_id, zone_type, surface_value, presence)
) |>
  dplyr::filter(!is.na(zone_id), !is.na(surface_value)) |>
  dplyr::mutate(zone_type = factor(zone_type))

glm_pa_sf <- glm(presence ~ surface_value + zone_type, family = binomial(), data = st_drop_geometry(pa_analysis_sf))
glm_pa_coef_sf <- data.frame(
  term = rownames(summary(glm_pa_sf)$coefficients),
  round(as.data.frame(summary(glm_pa_sf)$coefficients), 3),
  row.names = NULL
)

set.seed(42)
pseudo_points_terra <- spatSample(terra_aoi_qc, size = nrow(terra_points_qc), method = "random")
pseudo_points_terra$pseudo_id <- paste0("PA", sprintf("%02d", seq_len(nrow(pseudo_points_terra))))
pseudo_zone_vals_terra <- terra::extract(terra_zones_qc, pseudo_points_terra)
pseudo_surface_vals_terra <- terra::extract(surface_terra_qc, pseudo_points_terra)
pseudo_points_terra <- cbind(
  pseudo_points_terra,
  pseudo_zone_vals_terra[, c("zone_id", "zone_type"), drop = FALSE],
  pseudo_surface_vals_terra[, "surface_value", drop = FALSE]
)
pseudo_points_terra$presence <- 0L

presence_points_terra <- terra::project(terra_points_surface, "EPSG:32198")
presence_points_terra$presence <- 1L

pa_analysis_terra <- dplyr::bind_rows(
  as.data.frame(presence_points_terra)[, c("zone_id", "zone_type", "surface_value", "presence")],
  as.data.frame(pseudo_points_terra)[, c("zone_id", "zone_type", "surface_value", "presence")]
) |>
  dplyr::filter(!is.na(zone_id), !is.na(surface_value)) |>
  dplyr::mutate(zone_type = factor(zone_type))

glm_pa_terra <- glm(presence ~ surface_value + zone_type, family = binomial(), data = pa_analysis_terra)
glm_pa_coef_terra <- data.frame(
  term = rownames(summary(glm_pa_terra)$coefficients),
  round(as.data.frame(summary(glm_pa_terra)$coefficients), 3),
  row.names = NULL
)

points_xy_sf <- st_coordinates(points_sf_qc)
aoi_bbox_qc <- st_bbox(aoi_qc)
points_kde_sf_raw <- MASS::kde2d(
  points_xy_sf[, 1],
  points_xy_sf[, 2],
  n = 80,
  h = c(10000, 10000),
  lims = c(aoi_bbox_qc["xmin"], aoi_bbox_qc["xmax"], aoi_bbox_qc["ymin"], aoi_bbox_qc["ymax"])
)
points_kde_sf <- st_as_stars(
  list(kde = points_kde_sf_raw$z),
  dimensions = st_dimensions(x = points_kde_sf_raw$x, y = points_kde_sf_raw$y)
) |>
  st_set_crs(32198)
points_kde_sf <- points_kde_sf[aoi_qc]
points_kde_sf_plot <- points_kde_sf
points_kde_sf_plot[[1]] <- points_kde_sf_plot[[1]] / max(points_kde_sf_plot[[1]], na.rm = TRUE)

points_xy_terra <- crds(terra_points_qc)
aoi_ext_qc <- ext(terra_aoi_qc)
points_kde_terra_raw <- MASS::kde2d(
  points_xy_terra[, 1],
  points_xy_terra[, 2],
  n = 80,
  h = c(10000, 10000),
  lims = c(aoi_ext_qc$xmin, aoi_ext_qc$xmax, aoi_ext_qc$ymin, aoi_ext_qc$ymax)
)
points_kde_terra <- rast(
  t(points_kde_terra_raw$z)[nrow(t(points_kde_terra_raw$z)):1, ],
  extent = ext(
    min(points_kde_terra_raw$x),
    max(points_kde_terra_raw$x),
    min(points_kde_terra_raw$y),
    max(points_kde_terra_raw$y)
  ),
  crs = "EPSG:32198"
)
points_kde_terra <- mask(points_kde_terra, terra_aoi_qc)
points_kde_terra_plot <- points_kde_terra / global(points_kde_terra, "max", na.rm = TRUE)[1, 1]
