# Experimental R translation of selected Biomass Succession 5.2 calculations.
# Original C# authors: Robert M. Scheller and James B. Domingo.
# Source: LANDIS-II-Foundation/Extension-Biomass-Succession
# Revision: 490f4f260825c252fca3e3d5a24213af171ffed6. Apache-2.0.
# Modified: functional R interface, explicit state inputs and input validation.
# This is a calculation kernel, NOT a complete LANDIS-II simulator.

bs52_scalar <- function(x, name, lower = 0, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper) stop("Invalid ", name, call. = FALSE)
  invisible(x)
}

# CalculateCompetition: caller supplies the live collection at the exact
# moment of evaluation. The focal cohort has already been aged by one year
# relative to its entry in that collection. Do not pre-age the collection.
bs52_competition <- function(biomass, age, species, live_cohorts) {
  bs52_scalar(biomass, "biomass")
  bs52_scalar(age, "age")
  if (!is.data.frame(live_cohorts) ||
      !all(c("species", "age", "biomass") %in% names(live_cohorts)))
    stop("live_cohorts requires species, age, biomass")
  if (length(species) != 1L || is.na(species) ||
      anyNA(live_cohorts) || any(!is.finite(live_cohorts$biomass)) ||
      any(!is.finite(live_cohorts$age)) ||
      any(live_cohorts$biomass < 0) || any(live_cohorts$age < 0))
    stop("Invalid cohort state")
  focal_weight <- max(biomass^0.95, 1)
  include <- live_cohorts$age + 1 != age | live_cohorts$species != species
  focal_weight / (focal_weight + sum(pmax(live_cohorts$biomass[include]^0.95, 1)))
}

# ComputeChange's living-biomass calculations. All biomass/productivity units
# must match the original input units. Mortality/litter bookkeeping is separate.
# competition must come from the original live-state semantics above.
bs52_cohort_change <- function(biomass, age, longevity, mortality_shape,
                              growth_shape, max_anpp, max_biomass,
                              site_biomass, competition,
                              previous_year_mortality = 0,
                              capacity_reduction = 0, growth_reduction = 0,
                              defoliation = 0, spinup = FALSE,
                              spinup_mortality_fraction = 0) {
  vals <- list(biomass = biomass, age = age, mortality_shape = mortality_shape,
               growth_shape = growth_shape, max_anpp = max_anpp,
               max_biomass = max_biomass, site_biomass = site_biomass,
               previous_year_mortality = previous_year_mortality,
               spinup_mortality_fraction = spinup_mortality_fraction)
  for (nm in names(vals)) bs52_scalar(vals[[nm]], nm)
  bs52_scalar(longevity, "longevity", .Machine$double.eps)
  bs52_scalar(competition, "competition", .Machine$double.eps, 1)
  for (nm in c("capacity_reduction", "growth_reduction", "defoliation"))
    bs52_scalar(get(nm), nm, 0, 1)
  if (!is.logical(spinup) || length(spinup) != 1L || is.na(spinup))
    stop("spinup must be TRUE or FALSE")
  if (biomass != trunc(biomass) || age != trunc(age))
    stop("Original cohort biomass and age must be integers")

  mortality_age <- min(biomass * exp(age / longevity * mortality_shape) /
                         exp(mortality_shape), biomass)
  if (spinup && spinup_mortality_fraction > 0)
    mortality_age <- mortality_age + biomass * spinup_mortality_fraction
  mortality_age <- min(mortality_age, biomass)
  capacity_remaining <- 1 - capacity_reduction
  potential_biomass <- max(1, max_biomass * capacity_remaining -
                             site_biomass + biomass)
  if (capacity_remaining >= 1)
    potential_biomass <- max(potential_biomass, previous_year_mortality)
  b_ap <- biomass / potential_biomass
  b_pm <- competition
  power <- b_ap^growth_shape
  raw_anpp <- max_anpp * exp(1) * power * exp(-power) * b_pm
  raw_anpp <- min(max_anpp * b_pm, raw_anpp)
  if (growth_reduction > 0) raw_anpp <- raw_anpp * (1 - growth_reduction)
  # Preserve the actual source minimum of 1, including with zero maxANPP.
  actual_anpp <- max(1, raw_anpp - mortality_age)
  mortality_growth <- if (b_ap > 1) max_anpp * b_pm else
    max_anpp * (2 * b_ap) / (1 + b_ap) * b_pm
  mortality_growth <- min(biomass, mortality_growth, max_anpp * b_pm)
  if (growth_reduction > 0)
    mortality_growth <- mortality_growth * (1 - growth_reduction)
  mortality_growth <- min(max(0, mortality_growth - mortality_age), actual_anpp)
  total_mortality <- mortality_age + mortality_growth
  if (!is.finite(total_mortality) || total_mortality > biomass)
    stop("Mortality exceeds cohort biomass or is nonfinite")
  defoliation_loss <- 0.35 * actual_anpp * defoliation
  # C# integer cast truncates toward zero, including negative changes.
  delta <- trunc(actual_anpp - total_mortality - defoliation_loss)
  if (!is.finite(delta) || delta < -2147483648 || delta > 2147483647)
    stop("Delta is outside supported C# Int32 range")
  list(delta_biomass = delta, new_biomass = biomass + delta,
       raw_anpp = raw_anpp, actual_anpp = actual_anpp, mortality_age = mortality_age,
       mortality_growth = mortality_growth, total_mortality = total_mortality,
       defoliation_loss = defoliation_loss, potential_biomass = potential_biomass,
       b_ap = b_ap, competition = b_pm)
}

# InitialBiomass for a newly recruited cohort, NOT the initial-community spinup.
bs52_recruit_biomass <- function(max_anpp, ecoregion_max_biomass,
                                 live_cohorts, succession_timestep) {
  bs52_scalar(max_anpp, "max_anpp")
  bs52_scalar(ecoregion_max_biomass, "ecoregion_max_biomass", .Machine$double.eps)
  bs52_scalar(succession_timestep, "succession_timestep", 1)
  if (!is.data.frame(live_cohorts) ||
      !all(c("age", "biomass") %in% names(live_cohorts)) ||
      anyNA(live_cohorts[c("age", "biomass")]) ||
      any(!is.finite(live_cohorts$age)) || any(!is.finite(live_cohorts$biomass)) ||
      any(live_cohorts$age < 0) || any(live_cohorts$biomass < 0))
    stop("Invalid live_cohorts")
  selected <- live_cohorts$age >= succession_timestep & live_cohorts$age > 1
  actual <- sum(live_cohorts$biomass[selected])
  trunc(max(2, min(max_anpp, max_anpp * exp(-1.6 * actual / ecoregion_max_biomass))))
}

# ComputeShade: class thresholds must be in order 1, 2, 3, 4, 5.
bs52_shade <- function(live_cohorts, ecoregion_max_biomass,
                       minimum_relative_biomass, previous_year_mortality = 0) {
  bs52_scalar(ecoregion_max_biomass, "ecoregion_max_biomass", .Machine$double.eps)
  bs52_scalar(previous_year_mortality, "previous_year_mortality")
  if (!is.data.frame(live_cohorts) ||
      !all(c("age", "biomass") %in% names(live_cohorts)) ||
      anyNA(live_cohorts[c("age", "biomass")]) ||
      any(!is.finite(live_cohorts$age)) || any(!is.finite(live_cohorts$biomass)) ||
      any(live_cohorts$age < 0) || any(live_cohorts$biomass < 0))
    stop("Invalid live_cohorts")
  if (length(minimum_relative_biomass) != 5L ||
      any(!is.finite(minimum_relative_biomass)) ||
      any(minimum_relative_biomass < 0) || any(minimum_relative_biomass > 1) ||
      is.unsorted(minimum_relative_biomass)) stop("Invalid shade thresholds")
  actual <- sum(live_cohorts$biomass[live_cohorts$age > 5])
  actual <- min(ecoregion_max_biomass - previous_year_mortality, actual)
  classes <- which(actual / ecoregion_max_biomass >= minimum_relative_biomass)
  if (length(classes)) max(classes) else 0L
}
