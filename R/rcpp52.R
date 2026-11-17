bs52_rcpp_loaded <- function() {
  exists("bs52_grow_live_block_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_grow_live_block_dead_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_grow_live_block_dead_pools_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_ward_failure_log_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_recruit_scan_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_shade_light_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_append_recruits_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_recruit_records_cpp", mode = "function", inherits = TRUE) &&
    exists("bs52_record_landscape_cpp", mode = "function", inherits = TRUE)
}

bs52_record_landscape_rcpp <- function(cohorts, year, ctx) {
  out <- bs52_record_landscape_cpp(cohorts, as.integer(ctx$ids), ctx$config$species$species, as.numeric(year))
  # record_landscape() builds max_age with max(), therefore this column is a
  # double in the R reference.  Rcpp emits an IntegerVector for efficiency;
  # restore the public table type before the strict identical() comparison.
  out$max_age <- as.numeric(out$max_age)
  rownames(out) <- NULL
  out
}

bs52_ward_rcpp_loaded <- function() {
  exists("bs52_ward_failure_log_cpp", mode = "function", inherits = TRUE)
}

# sourceCpp creates several nested folders below cacheDir.  Keeping that
# directory beside a long simulation output path can exceed the Windows path
# limit before compilation starts.  A short per-session temporary cache also
# prevents result folders from accumulating compiler artefacts.
bs52_rcpp_cache_dir <- function() {
  cache_dir <- file.path(tempdir(), "bs52")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
}

bs52_init_rcpp <- function(cache_dir = bs52_rcpp_cache_dir()) {
  if (bs52_rcpp_loaded()) return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE))
    stop("Le moteur Rcpp est selectionne, mais le paquet Rcpp est absent. ",
         "Executer install.packages('Rcpp'), puis relancer le script.")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  source_file <- normalizePath("src/growth52.cpp", mustWork = TRUE)
  tryCatch(
    Rcpp::sourceCpp(
      file = source_file, env = .GlobalEnv, rebuild = FALSE,
      cacheDir = normalizePath(cache_dir), showOutput = FALSE,
      verbose = FALSE
    ),
    error = function(e) stop(
      "Compilation du moteur Rcpp impossible. Verifier Rtools sous Windows. ",
      "Erreur originale : ", conditionMessage(e), call. = FALSE
    )
  )
  if (!bs52_rcpp_loaded())
    stop("La compilation Rcpp n'a pas cree la fonction attendue.")
  invisible(TRUE)
}

bs52_grow_live_block_rcpp <- function(cohorts, params, timestep,
                                      combine_last = TRUE,
                                      spinup_fraction = 0,
                                      spinup = FALSE) {
  if (!bs52_rcpp_loaded())
    stop("Le moteur Rcpp doit etre initialise avec bs52_init_rcpp().")
  bs52_grow_live_block_cpp(
    cohorts, params, timestep, combine_last, spinup_fraction, spinup
  )
}

# Diagnostic path: emits annual dead-wood and litter inputs in cohort-processing
# order.  Production simulations use bs52_grow_live_block_dead_pools_rcpp()
# below to retain the pool vectors in C++ throughout the annual block.
bs52_grow_live_block_dead_rcpp <- function(cohorts, params,
                                           static_parameters, years, timestep,
                                           combine_last = TRUE,
                                           spinup_fraction = 0,
                                           spinup = FALSE) {
  if (!bs52_rcpp_loaded())
    stop("Le moteur Rcpp doit etre initialise avec bs52_init_rcpp().")
  if (!is.data.frame(static_parameters) ||
      !all(c("species", "leaf_longevity", "woody_decay", "leaf_lignin") %in%
           names(static_parameters)))
    stop("Parametres statiques invalides pour les pools morts.")
  years <- as.integer(years)
  if (length(years) != length(params) || anyNA(years))
    stop("Un an absolu est requis pour chaque table de parametres.")
  out <- bs52_grow_live_block_dead_cpp(
    cohorts, params, static_parameters, years, timestep, combine_last,
    spinup_fraction, spinup
  )
  rownames(out$cohorts) <- NULL
  rownames(out$events) <- NULL
  out
}

# Production dead-pool path.  A living-cohort state is shared by a group of
# cells, but its woody and litter pools are still individual vectors because
# stochastic recruitment may have made their prior stocks differ.  Updating
# those vectors next to cohort processing avoids serializing event data.frames
# back to R for every annual cohort change.
bs52_grow_live_block_dead_pools_rcpp <- function(
    cohorts, params, static_parameters, aet, years, timestep,
    wood, wood_decay, litter, litter_decay, combine_last = TRUE,
    spinup_fraction = 0, spinup = FALSE) {
  if (!bs52_rcpp_loaded())
    stop("Le moteur Rcpp doit etre initialise avec bs52_init_rcpp().")
  if (!is.data.frame(static_parameters) ||
      !all(c("species", "leaf_longevity", "woody_decay", "leaf_lignin") %in%
           names(static_parameters)))
    stop("Parametres statiques invalides pour les pools morts.")
  years <- as.integer(years)
  vectors <- list(wood = wood, wood_decay = wood_decay,
                  litter = litter, litter_decay = litter_decay)
  if (length(years) != length(params) || anyNA(years) || !is.finite(aet) ||
      !length(wood) || any(vapply(vectors, length, 0L) != length(wood)) ||
      any(!vapply(vectors, is.numeric, FALSE)) ||
      any(!is.finite(unlist(vectors))) || any(wood < 0 | litter < 0))
    stop("Etat ou sequence invalide pour les pools morts.")
  out <- bs52_grow_live_block_dead_pools_cpp(
    cohorts, params, static_parameters, as.numeric(aet), years, timestep,
    as.numeric(wood), as.numeric(wood_decay), as.numeric(litter),
    as.numeric(litter_decay), combine_last, spinup_fraction, spinup
  )
  rownames(out$cohorts) <- NULL
  if (!is.numeric(out$pre_final_site_biomass) ||
      length(out$pre_final_site_biomass) != 1L ||
      !is.finite(out$pre_final_site_biomass) ||
      out$pre_final_site_biomass < 0)
    stop("Le moteur Rcpp a retourne un contexte de feu invalide.")
  out
}

bs52_ward_probability_rcpp <- function(mature, active, targets, effective,
                                       maximum, cell_length, neighborhood,
                                       variant = c("ward2019", "ward2022")) {
  variant <- match.arg(variant)
  if (!bs52_ward_rcpp_loaded())
    stop("Le module Ward Rcpp doit etre initialise avec bs52_init_rcpp().")
  if (!is.matrix(mature) || !is.logical(mature) || !is.matrix(active) ||
      !is.logical(active) || !identical(dim(mature), dim(active)) ||
      anyNA(mature) || anyNA(active)) stop("Invalid source/active matrices")
  targets <- as.matrix(targets)
  if (ncol(targets) != 2 || anyNA(targets) || any(targets != trunc(targets)) ||
      any(targets[, 1] < 1 | targets[, 1] > nrow(active)) ||
      any(targets[, 2] < 1 | targets[, 2] > ncol(active)))
    stop("Invalid targets")
  destination_active <- active[targets]
  if (effective == -1) return(as.numeric(destination_active))
  probability <- bs52_ward_kernel(
    neighborhood$distance, effective, maximum, cell_length)
  if (variant == "ward2019") {
    support_limit <- bs52_ward2019_support_limit(
      maximum,
      cell_length
    )
    probability[neighborhood$distance > support_limit] <- 0
  }
  local <- mature[targets] & destination_active
  pending <- which(destination_active & !local)
  if (!length(pending)) return(as.numeric(local))
  pending_targets <- targets[pending, , drop = FALSE]
  storage.mode(pending_targets) <- "integer"
  active_failure <- probability > 0
  failure_log <- bs52_ward_failure_log_cpp(
    mature, active, pending_targets,
    as.integer(neighborhood$row[active_failure]),
    as.integer(neighborhood$col[active_failure]),
    log1p(-probability[active_failure]),
    identical(variant, "ward2019")
  )
  result <- as.numeric(local)
  result[pending] <- -expm1(failure_log)
  result
}

bs52_establish_landscape_rcpp <- function(cohorts, ids, targets, ecoregion,
                                          config, params_by_eco, year,
                                          scan, arrival, stochastic = TRUE) {
  species_names <- config$species$species
  n_sites <- length(ids)
  n_species <- length(species_names)
  mature <- scan$mature
  eco_active <- as.integer(ecoregion[ids])
  threshold_count <- nrow(config$shade)
  thresholds <- matrix(NA_real_, n_sites, threshold_count)
  sep <- matrix(NA_real_, n_sites, n_species)
  max_anpp <- matrix(NA_real_, n_sites, n_species)
  eco_max <- numeric(n_sites)
  for (eco in unique(eco_active)) {
    i <- which(eco_active == eco)
    p <- params_by_eco[[as.character(eco)]]
    if (!identical(as.character(p$species), species_names))
      stop("Ordre des especes incompatible dans les parametres dynamiques")
    threshold <- config$shade[[as.character(eco)]][order(config$shade$shade)]
    if (length(threshold) != threshold_count)
      stop("Nombre de seuils d'ombrage incompatible")
    for (h in seq_len(threshold_count)) thresholds[i, h] <- threshold[h]
    sep[i, ] <- matrix(rep(p$sep, each = length(i)), nrow = length(i))
    max_anpp[i, ] <- matrix(rep(p$max_anpp, each = length(i)), nrow = length(i))
    eco_max[i] <- max(p$max_biomass)
  }
  max_tolerance <- max(config$light$shade_tolerance)
  light <- matrix(NA_real_, max_tolerance, 6L)
  light[config$light$shade_tolerance, ] <- as.matrix(
    config$light[, paste0("shade_", 0:5), drop = FALSE]
  )
  shade_tolerance <- as.integer(
    config$species$shade_tolerance[match(species_names, config$species$species)]
  )
  shade_light <- bs52_shade_light_cpp(
    scan$shade_biomass, eco_max, thresholds, sep, shade_tolerance, light
  )
  probability <- shade_light$probability * arrival
  recruited <- matrix(FALSE, n_sites, n_species)
  if (stochastic) {
    # Conserver le meme ordre cellulaire des tirages que le chemin R valide.
    for (k in seq_len(n_sites))
      recruited[k, ] <- runif(n_species) < probability[k, ]
  } else recruited <- probability > 0
  cohorts <- bs52_append_recruits_cpp(
    cohorts, as.integer(ids), recruited, species_names, max_anpp, eco_max,
    scan$establishment_biomass
  )
  records <- bs52_recruit_records_cpp(
    species_names, mature, arrival, probability, recruited, as.integer(ids),
    as.integer(targets[, 1]), as.integer(targets[, 2]), as.integer(year)
  )
  # data.frame() conserve ici les noms d'especes comme noms de lignes; les
  # restaurer pour une equivalence exacte avec le chemin R historique.
  records$year <- as.numeric(records$year)
  # do.call(rbind, records) in the reference path disambiguates repeated
  # species names without a separator (ABIE.BAL1, ...).
  rownames(records) <- make.unique(rep(species_names, n_sites), sep = "")
  list(cohorts = cohorts, establishment = records)
}
