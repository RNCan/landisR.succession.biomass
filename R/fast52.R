# Sequential performance path for living-cohort growth.
#
# The public cohort state remains a data.frame.  During an annual update this
# function operates on atomic vectors and stable cohort identifiers, avoiding
# repeated data.frame copies.  It preserves the physical storage order within
# each species, which was validated against the supplied LANDIS-II traces.
bs52_new_biomass_unchecked <- function(biomass, age, longevity,
                                       mortality_shape, growth_shape,
                                       max_anpp, max_biomass, site_biomass,
                                       competition, spinup,
                                       spinup_mortality_fraction) {
  mortality_age <- min(biomass * exp(age / longevity * mortality_shape) /
                         exp(mortality_shape), biomass)
  if (spinup && spinup_mortality_fraction > 0)
    mortality_age <- min(mortality_age + biomass * spinup_mortality_fraction,
                         biomass)
  potential_biomass <- max(1, max_biomass - site_biomass + biomass)
  # Previous-year mortality is zero in the current living-biomass simulator.
  potential_biomass <- max(potential_biomass, 0)
  b_ap <- biomass / potential_biomass
  power <- b_ap^growth_shape
  raw_anpp <- max_anpp * exp(1) * power * exp(-power) * competition
  raw_anpp <- min(max_anpp * competition, raw_anpp)
  actual_anpp <- max(1, raw_anpp - mortality_age)
  mortality_growth <- if (b_ap > 1) max_anpp * competition else
    max_anpp * (2 * b_ap) / (1 + b_ap) * competition
  mortality_growth <- min(biomass, mortality_growth,
                          max_anpp * competition)
  mortality_growth <- min(max(0, mortality_growth - mortality_age),
                          actual_anpp)
  biomass + trunc(actual_anpp - mortality_age - mortality_growth)
}

bs52_grow_live_vectors_fast <- function(species, age, biomass, params,
                                        timestep, combine_young = FALSE,
                                        spinup_fraction = 0,
                                        spinup = TRUE) {
  if (!length(age))
    return(list(species = species, age = age, biomass = biomass))
  spp_order <- unique(species)
  site_biomass <- sum(biomass[age >= timestep])

  if (combine_young && timestep > 1) {
    keep <- rep(TRUE, length(age))
    add_species <- character()
    add_age <- numeric()
    add_biomass <- numeric()
    for (s in spp_order) {
      young <- species == s & age <= timestep
      if (any(young)) {
        keep[young] <- FALSE
        add_species <- c(add_species, s)
        add_age <- c(add_age, timestep - 1)
        add_biomass <- c(add_biomass, sum(biomass[young]))
      }
    }
    species <- c(species[keep], add_species)
    age <- c(age[keep], add_age)
    biomass <- c(biomass[keep], add_biomass)
  }
  if (!length(age))
    return(list(species = character(), age = numeric(), biomass = numeric()))

  # Stable grouping by species.  Do not sort cohorts by age: LANDIS-II keeps
  # their physical order inside SpeciesCohorts for these explicit communities.
  ord <- order(match(species, spp_order), seq_along(species))
  species <- species[ord]
  age <- age[ord]
  biomass <- biomass[ord]
  cohort_id <- seq_along(age)
  alive <- rep(TRUE, length(age))

  queues <- lapply(spp_order, function(s) cohort_id[species == s])
  names(queues) <- spp_order
  head_age <- function(s) age[queues[[s]][1L]]
  insert <- function(s, iters) {
    a <- head_age(s)
    before <- which(vapply(iters, function(z) head_age(z) < a, TRUE))
    if (!length(before)) c(iters, s) else append(iters, s, after = before[1L] - 1L)
  }
  iters <- character()
  for (s in spp_order) iters <- insert(s, iters)

  param_index <- match(spp_order, params$species)
  if (anyNA(param_index)) stop("Missing species parameters")
  names(param_index) <- spp_order

  while (length(iters)) {
    s <- iters[1L]
    id <- queues[[s]][1L]
    pi <- param_index[[s]]
    if (age[id] >= params$longevity[pi]) {
      alive[id] <- FALSE
    } else {
      new_age <- age[id] + 1
      live <- which(alive)
      focal_weight <- max(biomass[id]^0.95, 1)
      include <- age[live] + 1 != new_age | species[live] != s
      competition <- focal_weight /
        (focal_weight + sum(pmax(biomass[live[include]]^0.95, 1)))
      new_biomass <- bs52_new_biomass_unchecked(
        biomass[id], new_age, params$longevity[pi],
        params$mortality_shape[pi], params$growth_shape[pi],
        params$max_anpp[pi], params$max_biomass[pi], site_biomass,
        competition, spinup, spinup_fraction
      )
      age[id] <- new_age
      biomass[id] <- max(0, new_biomass)
      if (biomass[id] == 0) alive[id] <- FALSE
    }
    queues[[s]] <- queues[[s]][-1L]
    if (!length(queues[[s]])) {
      iters <- iters[-1L]
    } else if (length(iters) > 1L && head_age(s) < head_age(iters[2L])) {
      iters <- insert(s, iters[-1L])
    }
  }
  list(species = species[alive], age = age[alive], biomass = biomass[alive])
}

bs52_vectors_to_cohorts <- function(x) {
  data.frame(species = x$species, age = x$age, biomass = x$biomass,
             row.names = NULL)
}

bs52_grow_live_year_fast <- function(cohorts, params, timestep,
                                     combine_young = FALSE,
                                     spinup_fraction = 0,
                                     spinup = TRUE, dead = NULL) {
  # Dead-pool bookkeeping still uses the fully instrumented reference path.
  if (!is.null(dead))
    return(bs52_grow_live_year(cohorts, params, timestep, combine_young,
                               spinup_fraction, spinup, dead))
  if (!nrow(cohorts)) return(cohorts)
  bs52_vectors_to_cohorts(bs52_grow_live_vectors_fast(
    as.character(cohorts$species), cohorts$age, cohorts$biomass, params,
    timestep, combine_young, spinup_fraction, spinup
  ))
}

# Keeps one state in vector form for a complete deterministic block. This
# avoids rebuilding a data.frame after every annual update.
bs52_grow_live_block_fast <- function(cohorts, params, timestep,
                                      combine_last = TRUE,
                                      spinup_fraction = 0,
                                      spinup = FALSE) {
  if (!nrow(cohorts)) return(cohorts)
  state <- list(species = as.character(cohorts$species), age = cohorts$age,
                biomass = cohorts$biomass)
  for (i in seq_along(params)) {
    state <- bs52_grow_live_vectors_fast(
      state$species, state$age, state$biomass, params[[i]], timestep,
      combine_young = combine_last && i == length(params),
      spinup_fraction = spinup_fraction, spinup = spinup
    )
  }
  bs52_vectors_to_cohorts(state)
}

# Sequential Ward recruitment with the same cell-major RNG order as the
# reference implementation. The validated Ward kernel is unchanged.
bs52_ward_recruit_landscape_fast <- function(cohorts, ecoregion, config, year,
                                              variant = "ward2019",
                                              mode = c("ward", "local"),
                                              stochastic = TRUE,
                                              use_rcpp = FALSE,
                                              target_ids = NULL) {
  mode <- match.arg(mode)
  if (!is.matrix(ecoregion) || length(cohorts) != length(ecoregion) ||
      anyNA(ecoregion))
    stop("Need ecoregion matrix and one cohort table per cell")
  active <- ecoregion > 0
  all_ids <- which(active)
  ids <- if (is.null(target_ids)) all_ids else as.integer(target_ids)
  if (anyNA(ids) || anyDuplicated(ids) || any(ids < 1L | ids > length(cohorts)) ||
      any(!active[ids]))
    stop("Invalid recruitment target cells")
  targets <- arrayInd(ids, dim(ecoregion))
  species_names <- config$species$species
  n_species <- length(species_names)
  if (!length(all_ids) || !length(ids))
    return(list(cohorts = cohorts, establishment = data.frame()))

  maturity <- setNames(config$species$maturity, species_names)
  if (use_rcpp) {
    recruitment_scan <- bs52_recruit_scan_cpp(
      cohorts, as.integer(all_ids), species_names,
      as.integer(maturity[species_names]), as.integer(config$timestep)
    )
    mature_active <- recruitment_scan$mature
  } else {
    mature_active <- matrix(FALSE, length(all_ids), n_species,
                            dimnames = list(NULL, species_names))
    for (k in seq_along(all_ids)) {
      x <- cohorts[[all_ids[k]]]
      if (nrow(x)) {
        good <- x$biomass > 0 & x$age >= maturity[x$species]
        if (any(good))
          mature_active[k, match(unique(x$species[good]), species_names)] <- TRUE
      }
    }
  }

  arrival <- matrix(0, length(ids), n_species,
                    dimnames = list(NULL, species_names))
  target_row <- match(ids, all_ids)
  neighborhood <- bs52_ward_neighborhood(
    config$cell_length, max(config$species$max_seed_distance))
  for (j in seq_len(n_species)) {
    if (mode == "local") {
      arrival[, j] <- as.numeric(mature_active[target_row, j])
    } else {
      mature <- matrix(FALSE, nrow(ecoregion), ncol(ecoregion))
      mature[all_ids] <- mature_active[, j]
      sp <- config$species[j, ]
      probability_function <- if (use_rcpp) bs52_ward_probability_rcpp else
        bs52_ward_probability
      arrival[, j] <- probability_function(
        mature, active, targets, sp$effective_seed_distance,
        sp$max_seed_distance, config$cell_length, neighborhood, variant)
    }
  }

  # Strictly serotinous species are excluded from the ordinary seeding pathway.
  strict <- landisR_serotiny_strict_index(config$species)
  if (length(strict)) {
    arrival[, strict] <- 0
  }

  eco_values <- unique(ecoregion[ids])
  params_by_eco <- lapply(eco_values, function(e)
    bs52_dynamic_at(config, year, e))
  names(params_by_eco) <- as.character(eco_values)
  if (use_rcpp && identical(ids, all_ids))
    return(bs52_establish_landscape_rcpp(
      cohorts, ids, targets, ecoregion, config, params_by_eco, year,
      recruitment_scan, arrival, stochastic
    ))
  records <- vector("list", length(ids))
  for (k in seq_along(ids)) {
    id <- ids[k]
    eco <- ecoregion[id]
    p <- params_by_eco[[as.character(eco)]]
    x <- cohorts[[id]]
    eco_max <- max(p$max_biomass)
    thresholds <- config$shade[[as.character(eco)]][order(config$shade$shade)]
    actual <- sum(x$biomass[x$age > 5])
    actual <- min(eco_max, actual)
    shade_classes <- which(actual / eco_max >= thresholds)
    shade <- if (length(shade_classes)) max(shade_classes) else 0L
    light_rows <- match(p$shade_tolerance, config$light$shade_tolerance)
    if (anyNA(light_rows)) stop("Missing sufficient-light row")
    light <- as.numeric(config$light[light_rows, paste0("shade_", shade)])
    probability <- pmin(1, pmax(0, p$sep * light)) * arrival[k, ]
    recruited <- if (stochastic) runif(n_species) < probability else
      probability > 0
    if (any(recruited)) {
      selected <- x$age >= config$timestep & x$age > 1
      site_biomass <- sum(x$biomass[selected])
      b <- vapply(p$max_anpp[recruited], function(max_anpp)
        trunc(max(2, min(max_anpp,
                         max_anpp * exp(-1.6 * site_biomass / eco_max)))), 0)
      cohorts[[id]] <- rbind(
        x, data.frame(species = p$species[recruited], age = 1,
                      biomass = b, row.names = NULL))
    }
    records[[k]] <- data.frame(
      species = p$species, mature = mature_active[target_row[k], ],
      seed_probability = as.numeric(arrival[k, ]), probability = probability,
      recruited = recruited, cell_id = id, row = targets[k, 1],
      col = targets[k, 2], year = year
    )
  }
  list(cohorts = cohorts, establishment = do.call(rbind, records))
}
