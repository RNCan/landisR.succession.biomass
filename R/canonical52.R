# Reader and validator for canonical autonomous-R Biomass Succession inputs.

bs52_canonical_stop <- function(...) stop(..., call. = FALSE)

bs52_canonical_path <- function(root, ...) {
  path <- file.path(root, ...)
  if (!file.exists(path)) bs52_canonical_stop("Entree canonique absente : ", path)
  path
}

bs52_canonical_csv <- function(root, ...) {
  utils::read.csv(bs52_canonical_path(root, ...), stringsAsFactors = FALSE,
                  check.names = FALSE)
}

bs52_canonical_required_columns <- function(x, expected, label) {
  if (!identical(names(x), expected))
    bs52_canonical_stop(label, " doit contenir exactement : ",
                        paste(expected, collapse = ", "))
  x
}

bs52_canonical_settings <- function(root) {
  x <- bs52_canonical_csv(root, "model_settings.csv")
  bs52_canonical_required_columns(x, c("key", "value"), "model_settings.csv")
  if (anyDuplicated(x$key)) bs52_canonical_stop("Cles dupliquees dans model_settings.csv.")
  value <- setNames(as.character(x$value), x$key)
  required <- c("format", "source_duration", "cell_length", "timestep",
                "seeding_algorithm", "spinup_cohorts", "spinup_mortality_fraction",
                "initial_mode")
  if (!all(required %in% names(value)) ||
      !identical(value[["format"]], "landisR52/canonical-inputs-v1"))
    bs52_canonical_stop("model_settings.csv ne correspond pas au format attendu.")
  value
}

bs52_read_canonical52 <- function(directory) {
  root <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  settings <- bs52_canonical_settings(root)
  species <- bs52_canonical_csv(root, "ecology", "species.csv")
  species <- landisR_normalize_serotiny_traits(species)
  ecoregions <- bs52_canonical_csv(root, "ecology", "ecoregions.csv")
  species_parameters <- bs52_canonical_csv(root, "succession", "species_parameters.csv")
  ecoregion_parameters <- bs52_canonical_csv(root, "succession", "ecoregion_parameters.csv")
  dynamic <- bs52_canonical_csv(root, "succession", "dynamic.csv")
  shade <- bs52_canonical_csv(root, "succession", "shade.csv")
  light <- bs52_canonical_csv(root, "succession", "sufficient_light.csv")
  fire_reductions <- bs52_canonical_csv(root, "succession", "fire_reductions.csv")
  community_codes <- bs52_canonical_csv(root, "initial", "community_codes.csv")
  cohorts <- bs52_canonical_csv(root, "initial", "cohorts.csv")
  grid <- bs52_canonical_csv(root, "landscape", "grid.csv")

  bs52_canonical_required_columns(community_codes, c("map_code"),
                                  "initial/community_codes.csv")
  bs52_canonical_required_columns(cohorts, c("map_code", "species", "age", "biomass"),
                                  "initial/cohorts.csv")
  bs52_canonical_required_columns(grid,
    c("cell_id", "row", "col", "ecoregion_map_code", "ecoregion", "active", "community"),
    "landscape/grid.csv")
  if (anyDuplicated(community_codes$map_code) || anyDuplicated(cohorts[c("map_code", "species", "age")]))
    bs52_canonical_stop("Communautes ou cohortes initiales dupliquees.")
  if (any(!cohorts$map_code %in% community_codes$map_code) ||
      any(!cohorts$species %in% species$species) || any(cohorts$age < 1L) ||
      any(!is.finite(cohorts$biomass)) || any(cohorts$biomass < 0))
    bs52_canonical_stop("Cohortes initiales canoniques invalides.")
  if (identical(settings[["initial_mode"]], "explicit_biomass") && anyNA(cohorts$biomass))
    bs52_canonical_stop("La biomasse explicite manque dans initial/cohorts.csv.")
  if (!identical(settings[["initial_mode"]], "explicit_biomass") &&
      !identical(settings[["initial_mode"]], "age_spinup"))
    bs52_canonical_stop("Mode d'initialisation canonique inconnu.")

  nr <- max(grid$row); nc <- max(grid$col)
  active <- as.logical(grid$active)
  if (anyNA(active) || nrow(grid) != nr * nc ||
      !identical(grid$cell_id, seq_len(nr * nc)) ||
      any(grid$cell_id != grid$row + (grid$col - 1L) * nr) || !any(active))
    bs52_canonical_stop("Grille canonique invalide ou vide.")
  eco_index <- match(as.character(grid$ecoregion_map_code),
                     as.character(ecoregions$map_code))
  expected_active <- tolower(as.character(ecoregions$active[eco_index])) %in%
    c("yes", "true", "active")
  if (anyNA(eco_index) || any(as.character(grid$ecoregion) !=
      as.character(ecoregions$ecoregion[eco_index])) || any(active != expected_active) ||
      any(active & !grid$community %in% community_codes$map_code) ||
      any(!active & grid$community != 0L))
    bs52_canonical_stop("Grille et tables ecologiques canoniques incompatibles.")

  explicit_initial <- if (identical(settings[["initial_mode"]], "explicit_biomass")) {
    list(community_codes = as.integer(community_codes$map_code), cohorts = cohorts)
  } else NULL
  list(
    root = root, duration = as.integer(settings[["source_duration"]]),
    cell_length = as.numeric(settings[["cell_length"]]),
    timestep = as.integer(settings[["timestep"]]),
    spinup_cohorts = identical(tolower(settings[["spinup_cohorts"]]), "true"),
    seeding = settings[["seeding_algorithm"]],
    spinup_mortality = as.numeric(settings[["spinup_mortality_fraction"]]),
    species = species, species_parameters = species_parameters,
    ecoregions = ecoregions, ecoregion_parameters = ecoregion_parameters,
    dynamic = dynamic, initial_cohorts = cohorts[c("map_code", "species", "age")],
    community_codes = as.integer(community_codes$map_code), shade = shade,
    light = light, fire_reductions = fire_reductions,
    explicit_initial = explicit_initial,
    landscape = list(cells = grid, active = active,
                     ids = grid$cell_id[active], nrow = nr, ncol = nc,
                     order = "cell_id = row + (col - 1) * nrow; LANDIS rows")
  )
}
