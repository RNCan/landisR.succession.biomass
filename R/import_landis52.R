# Import a LANDIS-II Biomass Succession input directory into the canonical
# data-only format consumed by the autonomous R engine.

bs52_import_stop <- function(...) stop(..., call. = FALSE)

bs52_import_active <- function(ecoregions) {
  tolower(as.character(ecoregions$active)) %in% c("yes", "true", "active")
}

bs52_import_validate_grid <- function(grid, config) {
  required <- c("cell_id", "row", "col", "ecoregion_map_code",
                "ecoregion", "active", "community")
  if (!identical(names(grid), required) || !nrow(grid))
    bs52_import_stop("Grille canonique invalide.")
  nr <- max(grid$row); nc <- max(grid$col)
  if (anyNA(grid) || nrow(grid) != nr * nc ||
      !identical(grid$cell_id, seq_len(nr * nc)) ||
      any(grid$cell_id != grid$row + (grid$col - 1L) * nr))
    bs52_import_stop("La grille ne respecte pas l'ordre LANDIS colonne-major.")
  key <- as.character(config$ecoregions$map_code)
  i <- match(as.character(grid$ecoregion_map_code), key)
  if (anyNA(i) || any(as.character(grid$ecoregion) !=
                       as.character(config$ecoregions$ecoregion[i])) ||
      any(grid$active != bs52_import_active(config$ecoregions)[i]))
    bs52_import_stop("La grille est incompatible avec le fichier d'ecoregions.")
  if (any(grid$active & !grid$community %in% config$community_codes) ||
      any(!grid$active & grid$community != 0L))
    bs52_import_stop("Les communautes initiales sont incompatibles avec la grille.")
  grid
}

bs52_import_grid_from_cells <- function(cells, config) {
  required <- c("cell_id", "row", "col", "ecoregion", "community")
  if (!identical(names(cells), required))
    bs52_import_stop("cells.csv doit contenir exactement : ",
                     paste(required, collapse = ", "))
  # cells.csv may hold either LANDIS ecoregion labels (for example "216")
  # or the numeric map codes.  The latter is required for inactive cells:
  # map code 0 corresponds to the label "none" in landtypes.
  i <- match(as.character(cells$ecoregion),
             as.character(config$ecoregions$ecoregion))
  unmatched <- is.na(i)
  if (any(unmatched)) {
    i[unmatched] <- match(as.character(cells$ecoregion[unmatched]),
                          as.character(config$ecoregions$map_code))
  }
  if (anyNA(i))
    bs52_import_stop("cells.csv contient une ecoregion ou un code absent de landtypes.")
  grid <- data.frame(
    cell_id = as.integer(cells$cell_id), row = as.integer(cells$row),
    col = as.integer(cells$col),
    ecoregion_map_code = as.integer(config$ecoregions$map_code[i]),
    ecoregion = as.character(config$ecoregions$ecoregion[i]),
    active = as.logical(bs52_import_active(config$ecoregions)[i]),
    community = as.integer(cells$community), stringsAsFactors = FALSE
  )
  bs52_import_validate_grid(grid, config)
}

bs52_import_landis_matrix <- function(path, nr = NULL, nc = NULL) {
  if (!requireNamespace("terra", quietly = TRUE))
    bs52_import_stop("Importer les cartes LANDIS requiert terra : install.packages('terra').")
  r <- terra::rast(path)
  if (terra::nlyr(r) != 1L || (!is.null(nr) &&
      (terra::nrow(r) != nr || terra::ncol(r) != nc)))
    bs52_import_stop("Dimensions raster LANDIS incompatibles : ", path)
  raw <- terra::values(r, mat = FALSE)
  # terra::values() reads the first raster row first.  This is already the
  # LANDIS row order in the original TIFF files, as verified against cells.csv.
  # Do not reverse Y here: plotting an unreferenced LANDIS raster can look
  # mirrored, but reversing the stored values changes the landscape itself.
  matrix(raw, nrow = terra::nrow(r), ncol = terra::ncol(r), byrow = TRUE)
}

bs52_import_grid_from_maps <- function(config) {
  eco <- bs52_import_landis_matrix(config$ecoregions_map)
  community <- bs52_import_landis_matrix(config$initial_map, nrow(eco), ncol(eco))
  if (anyNA(eco) || anyNA(community) || any(!is.finite(eco)) ||
      any(!is.finite(community)) || any(eco != trunc(eco)) ||
      any(community != trunc(community)))
    bs52_import_stop("Les cartes LANDIS doivent contenir seulement des codes entiers finis.")
  code <- as.integer(as.vector(eco))
  i <- match(as.character(code), as.character(config$ecoregions$map_code))
  if (anyNA(i))
    bs52_import_stop("La carte d'ecoregions contient un code absent de landtypes.")
  nr <- nrow(eco); nc <- ncol(eco)
  position <- expand.grid(row = seq_len(nr), col = seq_len(nc))
  grid <- data.frame(
    cell_id = seq_len(nr * nc), row = as.integer(position$row),
    col = as.integer(position$col), ecoregion_map_code = code,
    ecoregion = as.character(config$ecoregions$ecoregion[i]),
    active = as.logical(bs52_import_active(config$ecoregions)[i]),
    community = as.integer(as.vector(community)), stringsAsFactors = FALSE
  )
  bs52_import_validate_grid(grid, config)
}

bs52_import_source_fingerprint <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[!file.info(files)$isdir]
  normalized <- normalizePath(files, winslash = "/", mustWork = TRUE)
  data.frame(
    source_relative_path = substring(normalized, nchar(normalizePath(root,
      winslash = "/", mustWork = TRUE)) + 2L),
    md5 = unname(tools::md5sum(normalized)), stringsAsFactors = FALSE
  )
}

bs52_import_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

bs52_import_landis_to_r <- function(source_directory = "inputs",
                                    destination = "inputs_R/landis_import",
                                    scenario = "scenario.txt",
                                    prefer_cells_csv = TRUE) {
  source_directory <- normalizePath(source_directory, winslash = "/", mustWork = TRUE)
  destination <- normalizePath(destination, winslash = "/", mustWork = FALSE)
  if (dir.exists(destination) && length(list.files(destination, all.files = TRUE,
                                                    no.. = TRUE)))
    bs52_import_stop("Le dossier destination existe deja et n'est pas vide : ",
                     destination, ".")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)

  config <- bs52_read_mauricie(source_directory, scenario)
  cells_path <- file.path(source_directory, "cells.csv")
  grid <- if (prefer_cells_csv && file.exists(cells_path)) {
    bs52_import_grid_from_cells(utils::read.csv(cells_path,
                                                stringsAsFactors = FALSE), config)
  } else bs52_import_grid_from_maps(config)

  initial_mode <- if (is.null(config$explicit_initial)) "age_spinup" else "explicit_biomass"
  initial_cohorts <- if (identical(initial_mode, "explicit_biomass")) {
    config$explicit_initial$cohorts
  } else transform(config$initial_cohorts, biomass = NA_real_)
  community_codes <- data.frame(map_code = as.integer(config$community_codes),
                                stringsAsFactors = FALSE)

  for (folder in c("landscape", "ecology", "succession", "initial", "provenance"))
    dir.create(file.path(destination, folder), recursive = TRUE, showWarnings = FALSE)
  bs52_import_write_csv(grid, file.path(destination, "landscape", "grid.csv"))
  bs52_import_write_csv(config$species, file.path(destination, "ecology", "species.csv"))
  bs52_import_write_csv(config$ecoregions, file.path(destination, "ecology", "ecoregions.csv"))
  bs52_import_write_csv(config$species_parameters,
                        file.path(destination, "succession", "species_parameters.csv"))
  bs52_import_write_csv(config$ecoregion_parameters,
                        file.path(destination, "succession", "ecoregion_parameters.csv"))
  bs52_import_write_csv(config$dynamic,
                        file.path(destination, "succession", "dynamic.csv"))
  bs52_import_write_csv(config$shade,
                        file.path(destination, "succession", "shade.csv"))
  bs52_import_write_csv(config$light,
                        file.path(destination, "succession", "sufficient_light.csv"))
  bs52_import_write_csv(config$fire_reductions,
                        file.path(destination, "succession", "fire_reductions.csv"))
  bs52_import_write_csv(community_codes,
                        file.path(destination, "initial", "community_codes.csv"))
  bs52_import_write_csv(initial_cohorts,
                        file.path(destination, "initial", "cohorts.csv"))
  settings <- data.frame(
    key = c("format", "source_duration", "cell_length", "timestep",
            "seeding_algorithm", "spinup_cohorts", "spinup_mortality_fraction",
            "initial_mode"),
    value = c("landisR52/canonical-inputs-v1", as.character(config$duration),
              as.character(config$cell_length), as.character(config$timestep),
              as.character(config$seeding), as.character(config$spinup_cohorts),
              as.character(config$spinup_mortality), initial_mode),
    stringsAsFactors = FALSE
  )
  bs52_import_write_csv(settings, file.path(destination, "model_settings.csv"))
  source_fingerprint <- bs52_import_source_fingerprint(source_directory)
  bs52_import_write_csv(source_fingerprint,
                        file.path(destination, "provenance", "source_files.csv"))
  manifest <- list(
    format = "landisR52/canonical-inputs-v1", imported_at = Sys.time(),
    source_directory = source_directory, source_scenario = scenario,
    import_grid = if (prefer_cells_csv && file.exists(cells_path)) "cells_csv" else "landis_maps",
    source_fingerprint = source_fingerprint
  )
  saveRDS(manifest, file.path(destination, "provenance", "import_manifest.rds"),
          compress = "gzip")

  destination
}
