# Experimental living-biomass spinup. Sources: InitialBiomass.cs at 490f4f2;
# Library-Biomass-Cohort at e9ab242d9740212e6717d63dfe539b95a45b49b6.
# Original authors: Robert M. Scheller, James B. Domingo. Apache-2.0.
# Optional dead-pool collector receives annual growth and longevity deaths.
# Preserve each species collection storage order. Validated against 330 cohort
# growth calculations from five LANDIS-II traces supplied by the user.
bs52_grow_live_year <- function(cohorts, params, timestep,
                                combine_young = FALSE, spinup_fraction = 0,
                                spinup = TRUE, dead = NULL) {
  site_biomass <- sum(cohorts$biomass[cohorts$age >= timestep])
  spp_order <- unique(cohorts$species)
  if (combine_young && timestep > 1) {
    for (s in spp_order) {
      young <- cohorts$species == s & cohorts$age <= timestep
      if (any(young)) {
        b <- sum(cohorts$biomass[young])
        cohorts <- cohorts[!young, , drop = FALSE]
        cohorts <- rbind(cohorts, data.frame(species = s, age = timestep - 1,
                                            biomass = b))
      }
    }
  }
  if (!nrow(cohorts)) return(cohorts)
  cohorts <- cohorts[order(match(cohorts$species, spp_order), seq_len(nrow(cohorts))), ]
  rownames(cohorts) <- NULL
  cohorts$id <- seq_len(nrow(cohorts))
  queues <- lapply(spp_order, function(s) cohorts$id[cohorts$species == s])
  names(queues) <- spp_order
  age_of <- function(id) cohorts$age[match(id, cohorts$id)]
  insert <- function(s, iters) {
    age <- age_of(queues[[s]][1])
    before <- which(vapply(iters, function(z) age_of(queues[[z]][1]) < age, TRUE))
    if (!length(before)) c(iters, s) else append(iters, s, after = before[1] - 1)
  }
  iters <- character()
  for (s in spp_order) iters <- insert(s, iters)
  while (length(iters)) {
    s <- iters[1]; id <- queues[[s]][1]; i <- match(id, cohorts$id)
    p <- params[match(s, params$species), ]
    if (is.na(p$species)) stop("Missing parameters for ", s)
    if (cohorts$age[i] >= p$longevity) {
      if (!is.null(dead)) {
        comp <- bs52_competition(cohorts$biomass[i], cohorts$age[i], s, cohorts)
        r <- bs52_cohort_change(cohorts$biomass[i], cohorts$age[i], p$longevity,
          p$mortality_shape, p$growth_shape, p$max_anpp, p$max_biomass,
          site_biomass, comp, spinup = spinup,
          spinup_mortality_fraction = spinup_fraction)
        bs52_dead_death(dead, cohorts[i, ], p, max(0, r$raw_anpp - r$mortality_age))
      }
      cohorts <- cohorts[-i, , drop = FALSE]
    } else {
      new_age <- cohorts$age[i] + 1
      comp <- bs52_competition(cohorts$biomass[i], new_age, s, cohorts)
      r <- bs52_cohort_change(cohorts$biomass[i], new_age, p$longevity,
        p$mortality_shape, p$growth_shape, p$max_anpp, p$max_biomass,
        site_biomass, comp, spinup = spinup,
        spinup_mortality_fraction = spinup_fraction)
      if (!is.null(dead)) bs52_dead_growth(dead, r, p)
      cohorts$age[i] <- new_age
      cohorts$biomass[i] <- max(0, r$new_biomass)
      if (cohorts$biomass[i] == 0) cohorts <- cohorts[-i, , drop = FALSE]
    }
    queues[[s]] <- queues[[s]][-1]
    if (!length(queues[[s]])) iters <- iters[-1] else if (length(iters) > 1 &&
        age_of(queues[[s]][1]) < age_of(queues[[iters[2]]][1]))
      iters <- insert(s, iters[-1])
  }
  cohorts$id <- NULL
  cohorts
}

bs52_spinup_live <- function(config, community_code = 1, ecoregion = 314, dead = NULL) {
  initial <- config$initial_cohorts[config$initial_cohorts$map_code == community_code, ]
  if (!community_code %in% config$community_codes) stop("Unknown community")
  cohorts <- data.frame(species = character(), age = numeric(), biomass = numeric())
  if (!nrow(initial)) return(cohorts)
  step <- config$timestep
  if (any(initial$age %% step != 0)) stop("Initial ages must be timestep multiples")
  initial <- initial[order(-initial$age, seq_len(nrow(initial))), ]
  dy <- config$dynamic[config$dynamic$year == 0 & config$dynamic$ecoregion == ecoregion, ]
  params <- merge(config$species, config$species_parameters, by = "species", sort = FALSE)
  params <- merge(params, dy, by = "species", sort = FALSE)
  if (any(!initial$species %in% params$species)) stop("Missing spinup parameters")
  eco_max <- max(dy$max_biomass)
  end <- if (step == 1) -1 else 0
  for (time in seq(-max(initial$age), end, by = step)) {
    for (year in seq_len(step)) {
      cohorts <- bs52_grow_live_year(cohorts, params,
        step, combine_young = year == step,
        spinup_fraction = config$spinup_mortality, dead = dead)
      if (!is.null(dead)) bs52_dead_end_year(dead, time + year - step)
    }
    births <- initial[initial$age == -time, ]
    for (j in seq_len(nrow(births))) {
      s <- births$species[j]
      b <- bs52_recruit_biomass(params$max_anpp[match(s, params$species)],
                                eco_max, cohorts, step)
      cohorts <- rbind(cohorts, data.frame(species = s, age = 1, biomass = b))
    }
  }
  rownames(cohorts) <- NULL
  cohorts
}
