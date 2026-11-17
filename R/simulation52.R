# Experimental forward simulator for the living-biomass portion of Biomass
# Succession 5.2. It is intentionally explicit about unsupported processes.

bs52_dynamic_at <- function(config, year, ecoregion) {
  d <- config$dynamic[config$dynamic$ecoregion == ecoregion &
                        config$dynamic$year <= year, , drop = FALSE]
  if (!nrow(d)) stop("No dynamic parameters for year ", year,
                     " and ecoregion ", ecoregion)
  d <- d[order(d$species, d$year), , drop = FALSE]
  d <- d[!duplicated(d$species, fromLast = TRUE), , drop = FALSE]
  p <- merge(config$species, config$species_parameters,
             by = "species", sort = FALSE)
  p <- merge(p, d, by = "species", sort = FALSE)
  if (nrow(p) != nrow(config$species))
    stop("Dynamic parameters do not cover every species")
  p
}

bs52_light_probability <- function(config, shade_tolerance, actual_shade) {
  row <- config$light[config$light$shade_tolerance == shade_tolerance,
                      , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing sufficient-light row")
  as.numeric(row[[paste0("shade_", actual_shade)]])
}

bs52_establish <- function(cohorts, params, config, ecoregion, timestep,
                           stochastic = TRUE, seed_probability = NULL) {
  eco_col <- as.character(ecoregion)
  if (!eco_col %in% names(config$shade)) stop("Missing shade thresholds")
  thresholds <- config$shade[[eco_col]][order(config$shade$shade)]
  eco_max <- max(params$max_biomass)
  shade <- bs52_shade(cohorts, eco_max, thresholds)
  mature <- vapply(params$species, function(s) {
    x <- cohorts[cohorts$species == s, , drop = FALSE]
    any(x$age >= params$maturity[params$species == s])
  }, logical(1))
  # Local mature-seed-source approximation, NOT a port of Ward dispersal.
  # One stochastic site does not represent the distribution across a landscape.
  light <- mapply(bs52_light_probability,
                  MoreArgs = list(config = config, actual_shade = shade),
                  shade_tolerance = params$shade_tolerance)
  probability <- pmin(1, pmax(0, params$sep * light))
  arrival <- if (is.null(seed_probability)) as.numeric(mature) else {
    if (is.null(names(seed_probability)) || anyDuplicated(names(seed_probability)))
      stop("Seed probabilities require unique species names")
    as.numeric(seed_probability[params$species])
  }
  if (length(arrival) != nrow(params) || anyNA(arrival) ||
      any(!is.finite(arrival)) || any(arrival < 0 | arrival > 1))
    stop("Invalid seed probabilities")

  # Strictly serotinous species cannot use ordinary sexual recruitment.
  strict <- landisR_serotiny_strict_index(config$species)
  if (length(strict)) {
    strict_species <- as.character(config$species$species[strict])
    strict_params <- match(strict_species, as.character(params$species))
    strict_params <- strict_params[!is.na(strict_params)]
    arrival[strict_params] <- 0
  }

  probability <- probability * arrival
  recruited <- if (stochastic) runif(nrow(params)) < probability else probability > 0
  for (i in which(recruited)) {
    s <- params$species[i]
    b <- bs52_recruit_biomass(params$max_anpp[i], eco_max, cohorts, timestep)
    cohorts <- rbind(cohorts, data.frame(species = s, age = 1, biomass = b))
  }
  list(cohorts = cohorts, shade = shade,
       establishment = data.frame(species = params$species,
                                  mature = mature, seed_probability = arrival, probability = probability,
                                  recruited = recruited))
}

bs52_simulate_site <- function(config, community_code = 1, ecoregion = 314,
                               duration = config$duration, seed = 1,
                               stochastic_establishment = TRUE,
                               record_cohorts = TRUE, track_dead = FALSE,
                               explicit_initial = NULL, dispersal_provider = NULL) {
  if (duration < 0 || duration != trunc(duration)) stop("Invalid duration")
  bs52_scalar(duration, "duration")
  bs52_scalar(config$timestep, "timestep", 1)
  if (config$timestep != trunc(config$timestep)) stop("Timestep must be integer")
  set.seed(seed)
  dead <- if (track_dead) bs52_dead_state(config, ecoregion) else NULL
  if (!is.null(explicit_initial) && track_dead)
    stop("Explicit initial dead pools must be specified before tracking dead biomass")
  cohorts <- if (is.null(explicit_initial))
    bs52_spinup_live(config, community_code, ecoregion, dead) else
    bs52_initialize_explicit(explicit_initial, community_code)
  initial_dead <- if (track_dead) bs52_dead_snapshot(dead, 0) else NULL
  if (track_dead) dead$history <- list()
  totals <- list(); details <- list(); recruitment <- list()
  record <- function(year) {
    x <- if (nrow(cohorts)) aggregate(biomass ~ species, cohorts, sum) else
      data.frame(species = character(), biomass = numeric())
    out <- merge(data.frame(species = config$species$species), x,
                 by = "species", all.x = TRUE, sort = FALSE)
    out$biomass[is.na(out$biomass)] <- 0
    out$year <- year
    out[, c("year", "species", "biomass")]
  }
  totals[[1]] <- record(0)
  if (record_cohorts) details[[1]] <- transform(cohorts, year = rep(0, nrow(cohorts)))
  if (duration > 0) for (year in seq_len(duration)) {
    params <- bs52_dynamic_at(config, year, ecoregion)
    cohorts <- bs52_grow_live_year(
      cohorts, params, config$timestep,
      combine_young = year %% config$timestep == 0,
      spinup_fraction = 0, spinup = FALSE, dead = dead
    )
    if (track_dead) bs52_dead_end_year(dead, year)
    if (year %% config$timestep == 0) {
      arrival <- if (is.null(dispersal_provider)) NULL else
        dispersal_provider(year = year, cohorts = cohorts, params = params)
      est <- bs52_establish(cohorts, params, config, ecoregion,
                            config$timestep, stochastic_establishment, arrival)
      cohorts <- est$cohorts
      recruitment[[length(recruitment) + 1L]] <-
        transform(est$establishment, year = year, shade = est$shade)
    }
    if (year %% config$timestep == 0 || year == duration) {
      totals[[length(totals) + 1L]] <- record(year)
      if (record_cohorts)
        details[[length(details) + 1L]] <- transform(cohorts, year = rep(year, nrow(cohorts)))
    }
  }
  species_biomass <- do.call(rbind, totals)
  total_biomass <- aggregate(biomass ~ year, species_biomass, sum)
  list(dead_initial = initial_dead,
       dead_history = if (track_dead && length(dead$history)) do.call(rbind, dead$history) else data.frame(),
       final_cohorts = cohorts, species_biomass = species_biomass,
       total_biomass = total_biomass,
       cohort_history = if (record_cohorts) do.call(rbind, details) else NULL,
       recruitment = if (length(recruitment)) do.call(rbind, recruitment) else
         data.frame())
}
