# Ward source references: Library-Succession 8c8617d (2019), 008c84f (2022).
# See reference/*.cs. Analytical probabilities reproduce marginal event laws,
# not LANDIS random-number streams. Input matrices use [row, column].
bs52_ward_neighborhood <- function(cell_length, max_seed_distance) {
  bs52_scalar(cell_length, "cell length", .Machine$double.eps)
  bs52_scalar(max_seed_distance, "maximum seed distance")
  # LANDIS-II coarse-grid safeguard: Ward integrates dispersal across the
  # destination-cell width. Always retain at least the first cardinal ring
  # even when a species MaxSeedDistance is shorter than one raster cell.
  radius <- max(
    cell_length,
    max_seed_distance + cell_length / 2
  )
  n <- trunc(radius / cell_length) + 1L
  x <- expand.grid(row = seq_len(n), col = 0:n)
  x$distance <- cell_length * sqrt(x$row^2 + x$col^2)
  x <- x[x$distance <= radius, ]
  x[order(x$distance, x$row, x$col), ]
}

# Maximum centre-to-centre distance considered by the landisR Ward-2019
# implementation. LANDIS-II evaluates the first adjacent cardinal cell even
# when MaxSeedDistance is shorter than cell size; retain that one-cell
# safeguard while preserving the existing species-specific outer support.
bs52_ward2019_support_limit <- function(maximum, cell_length) {
  bs52_scalar(maximum, "maximum distance")
  bs52_scalar(cell_length, "cell length", .Machine$double.eps)
  max(
    cell_length,
    maximum + cell_length / 2 * 1.414
  )
}

bs52_ward_kernel <- function(distance, effective, maximum, cell_length) {
  bs52_scalar(effective, "effective distance", .Machine$double.eps)
  bs52_scalar(maximum, "maximum distance", effective)
  bs52_scalar(cell_length, "cell length", .Machine$double.eps)
  if (anyNA(distance) || any(!is.finite(distance)) || any(distance < 0))
    stop("Invalid distances")
  lo <- pmax(0, distance - cell_length)
  l1 <- log(.05) / effective; l2 <- log(.01) / maximum
  tail <- .05 * (exp(l2 * (lo - effective)) - exp(l2 * (distance - effective)))
  if (cell_length <= effective)
    p <- ifelse(distance <= effective, exp(l1 * lo) - exp(l1 * distance), tail)
  else p <- ifelse(distance <= cell_length,
                  exp(l1 * lo) - .05 * exp(l2 * (distance - effective)), tail)
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid Ward kernel")
  p
}
bs52_ward_probability <- function(mature, active, targets, effective, maximum,
                                  cell_length, neighborhood,
                                  variant = c("ward2019", "ward2022")) {
  variant <- match.arg(variant)
  if (!is.matrix(mature) || !is.logical(mature) || !is.matrix(active) ||
      !is.logical(active) || !identical(dim(mature), dim(active)) ||
      anyNA(mature) || anyNA(active)) stop("Invalid source/active matrices")
  targets <- as.matrix(targets)
  if (ncol(targets) != 2 || anyNA(targets) || any(targets != trunc(targets)) ||
      any(targets[,1] < 1 | targets[,1] > nrow(active)) ||
      any(targets[,2] < 1 | targets[,2] > ncol(active))) stop("Invalid targets")
  destination_active <- active[targets]
  if (effective == -1) return(as.numeric(destination_active))
  prob <- bs52_ward_kernel(neighborhood$distance, effective, maximum, cell_length)
  if (variant == "ward2019") {
    support_limit <- bs52_ward2019_support_limit(
      maximum,
      cell_length
    )
    prob[neighborhood$distance > support_limit] <- 0
  }
  local <- mature[targets] & destination_active
  pending <- which(destination_active & !local)
  failure_log <- numeric(length(pending))
  rr <- targets[pending, 1]; cc <- targets[pending, 2]
  for (j in which(prob > 0)) {
    dr <- neighborhood$row[j]; dc <- neighborhood$col[j]
    offsets <- if (dc == 0) rbind(c(dr,0),c(0,dr),c(-dr,0),c(0,-dr)) else
      rbind(c(dr,dc),c(-dr,dc),c(-dr,-dc),c(dr,-dc))
    count <- integer(length(pending))
    for (k in seq_len(4L)) {
      r <- rr + offsets[k,1]; c <- cc + offsets[k,2]
      good <- which(r >= 1 & r <= nrow(active) & c >= 1 & c <= ncol(active))
      idx <- cbind(r[good], c[good])
      count[good] <- count[good] + as.integer(mature[idx] & active[idx])
    }
    power <- if (variant == "ward2019") count else as.integer(count > 0)
    # Only positive counts: avoid 0 * -Inf for a kernel of exactly 1.
    good <- which(power > 0)
    failure_log[good] <- failure_log[good] + power[good] * log1p(-prob[j])
  }
  result <- as.numeric(local)
  result[pending] <- -expm1(failure_log)
  result
}

# One synchronous sexual-reproduction step for a small/medium R landscape.
# cohorts is a list of length length(ecoregion), in R column-major cell order.
# No growth or disturbance is performed here. Freeze source maturity for this step.
bs52_ward_recruit_landscape <- function(cohorts, ecoregion, config, year,
                                        variant = "ward2019", mode = c("ward", "local"),
                                        stochastic = TRUE) {
  mode <- match.arg(mode)
  if (!is.matrix(ecoregion) || length(cohorts) != length(ecoregion) || anyNA(ecoregion))
    stop("Need ecoregion matrix and one cohort table per cell")
  active <- ecoregion > 0
  ids <- which(active); targets <- arrayInd(ids, dim(ecoregion))
  species <- config$species$species
  if (!length(ids)) return(list(cohorts = cohorts, establishment = data.frame()))
  probabilities <- matrix(0, length(ids), length(species), dimnames = list(NULL,species))
  neighborhood <- bs52_ward_neighborhood(config$cell_length,
                                         max(config$species$max_seed_distance))
  for (j in seq_along(species)) {
    sp <- config$species[j, ]
    mature <- matrix(FALSE, nrow(ecoregion), ncol(ecoregion))
    mature[ids] <- vapply(ids, function(id) {
      x <- cohorts[[id]]
      any(x$species == sp$species & x$age >= sp$maturity & x$biomass > 0)
    }, logical(1))
    probabilities[,j] <- if (mode == "local") as.numeric(mature[ids]) else
      bs52_ward_probability(mature, active, targets, sp$effective_seed_distance,
        sp$max_seed_distance, config$cell_length, neighborhood, variant)
  }
  params_by_eco <- lapply(unique(ecoregion[ids]), function(e)
    bs52_dynamic_at(config, year, e))
  names(params_by_eco) <- as.character(unique(ecoregion[ids]))
  records <- vector("list", length(ids))
  for (k in seq_along(ids)) {
    id <- ids[k]; eco <- ecoregion[id]
    result <- bs52_establish(cohorts[[id]], params_by_eco[[as.character(eco)]],
      config, eco, config$timestep, stochastic, probabilities[k, ])
    cohorts[[id]] <- result$cohorts
    records[[k]] <- transform(result$establishment, cell_id = id,
                              row = targets[k,1], col = targets[k,2], year = year)
  }
  list(cohorts = cohorts, establishment = do.call(rbind, records))
}
