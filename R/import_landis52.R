# Import a LANDIS-II Biomass Succession input directory into the canonical
# data-only format consumed by the autonomous R engine.

bs52_import_stop <- function(...) stop(..., call. = FALSE)

bs52_import_active <- function(ecoregions) {
  tolower(as.character(ecoregions$active)) %in% c("yes", "true", "active")
}

bs52_import_validate_grid <- function(grid, ecoregions, community_codes) {
  required <- c("cell_id", "row", "col", "ecoregion_map_code",
                "ecoregion", "active", "community")
  if (!identical(names(grid), required) || !nrow(grid))
    bs52_import_stop("Grille canonique invalide.")
  nr <- max(grid$row); nc <- max(grid$col)
  if (anyNA(grid) || nrow(grid) != nr * nc ||
      !identical(grid$cell_id, seq_len(nr * nc)) ||
      any(grid$cell_id != grid$row + (grid$col - 1L) * nr))
    bs52_import_stop("La grille ne respecte pas l'ordre LANDIS colonne-major.")
  key <- as.character(ecoregions$map_code)
  i <- match(as.character(grid$ecoregion_map_code), key)
  if (anyNA(i) || any(as.character(grid$ecoregion) !=
                       as.character(ecoregions$ecoregion[i])) ||
      any(grid$active != bs52_import_active(ecoregions)[i]))
    bs52_import_stop("La grille est incompatible avec le fichier d'ecoregions.")
  if (any(grid$active & !grid$community %in% community_codes) ||
      any(!grid$active & grid$community != 0L))
    bs52_import_stop("Les communautes initiales sont incompatibles avec la grille.")
  grid
}

bs52_import_landis_matrix <- function(path) {
  if (!requireNamespace("terra", quietly = TRUE))
    bs52_import_stop("Importer les cartes LANDIS requiert terra : install.packages('terra').")
  r <- terra::rast(path)
  if (terra::nlyr(r) != 1L)
    bs52_import_stop("Dimensions raster LANDIS incompatibles : ", path)
  raw <- terra::values(r, mat = FALSE)
  # terra::values() reads the first raster row first.  This is already the
  # LANDIS row order in the original TIFF files, as verified against cells.csv.
  # Do not reverse Y here: plotting an unreferenced LANDIS raster can look
  # mirrored, but reversing the stored values changes the landscape itself.
  matrix(raw, nrow = terra::nrow(r), ncol = terra::ncol(r), byrow = TRUE)
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

bs52_assemble_landis_inputs <- function(project, config) {
  common <- project$canonical
  species <- common$species$species
  if (!all(species %in% config$species_parameters$species) ||
      !all(config$dynamic$species %in% species) ||
      !all(config$initial_cohorts$species %in% species))
    bs52_import_stop("The Biomass Succession species must match the landscape species catalog.")
  columns <- c("cell_id", "row", "col", "ecoregion_map_code", "ecoregion", "active")
  if (!all(columns %in% names(common$landscape$cells)))
    bs52_import_stop("Biomass Succession requires ecoregion_map_code in the landscape cells.")
  grid <- common$landscape$cells[columns]
  if (is.null(common$ecoregions)) {
    common$ecoregions <- unique(grid[c("ecoregion_map_code", "ecoregion", "active")])
    names(common$ecoregions)[1L] <- "map_code"
  }
  community_codes <- as.integer(config$community_codes)
  if (!is.null(config$initial_map)) {
    if (!identical(dim(config$initial_map), c(common$landscape$nrow, common$landscape$ncol)))
      bs52_import_stop("The initial communities map must match the landscape dimensions.")
    community <- as.vector(config$initial_map)
  } else {
    cells <- config$community_by_cell
    if (anyDuplicated(cells$cell_id) || !setequal(cells$cell_id, grid$cell_id))
      bs52_import_stop("The initial communities in cells.csv must match the landscape cell IDs.")
    community <- cells$community[match(grid$cell_id, cells$cell_id)]
  }
  if (!is.numeric(community) || any(!is.finite(community) | community != trunc(community) |
                                    community < 0 | community > .Machine$integer.max))
    bs52_import_stop("Initial communities must have finite, nonnegative integer codes.")
  grid$community <- as.integer(community)
  grid <- bs52_import_validate_grid(grid, common$ecoregions, community_codes)
  explicit <- !is.null(config$explicit_initial)
  cohorts <- if (explicit) config$explicit_initial$cohorts else config$initial_cohorts
  # The imported ages and codes are integral, as in the equivalent canonical CSV.
  cohorts$age <- as.integer(cohorts$age)
  cohorts$map_code <- as.integer(cohorts$map_code)
  assembled <- landisR::landis_build(
    landscape = grid,
    species = common$species, species_parameters = config$species_parameters,
    ecoregions = common$ecoregions, ecoregion_parameters = config$ecoregion_parameters,
    dynamic = config$dynamic, initial_cohorts = cohorts,
    shade = config$shade, light = config$light, fire_reductions = config$fire_reductions,
    duration = common$duration, cell_length = common$cell_length,
    timestep = config$timestep, seeding_algorithm = config$seeding,
    spinup_cohorts = config$spinup_cohorts, spinup_mortality_fraction = config$spinup_mortality,
    initial_mode = if (explicit) "explicit_biomass" else "age_spinup",
    community_codes = community_codes, name = project$name
  )
  assembled$source_type <- "LANDIS-II"
  assembled$source <- list(landscape = project$source, succession = config$scenario_file)
  assembled
}

bs52_import_landis_to_r <- function(source_directory = "inputs",
                                    destination = "inputs_R/landis_import",
                                    scenario = "scenario.txt",
                                    prefer_cells_csv = TRUE) {
  landscape <- landisR::landis_read(source_directory, scenario,
                                     prefer_cells_csv = prefer_cells_csv)
  succession <- biomass_succession_read(source_directory, scenario, prefer_cells_csv)
  project <- bs52_assemble_landis_inputs(landscape, succession)
  config <- project$canonical
  grid <- config$landscape$cells
  source_directory <- landscape$source$directory
  destination <- normalizePath(destination, winslash = "/", mustWork = FALSE)
  if (dir.exists(destination) && length(list.files(destination, all.files = TRUE,
                                                    no.. = TRUE)))
    bs52_import_stop("Le dossier destination existe deja et n'est pas vide : ",
                     destination, ".")
  initial_mode <- if (is.null(config$explicit_initial)) "age_spinup" else "explicit_biomass"
  initial_cohorts <- if (identical(initial_mode, "explicit_biomass")) {
    config$explicit_initial$cohorts
  } else transform(config$initial_cohorts, biomass = NA_real_)
  community_codes <- data.frame(map_code = config$community_codes)

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
    import_grid = if (is.null(succession$community_by_cell)) "landis_maps" else "cells_csv",
    source_fingerprint = source_fingerprint
  )
  saveRDS(manifest, file.path(destination, "provenance", "import_manifest.rds"),
          compress = "gzip")

  destination
}
