# The serial v7 release tests one direction at a time after the light and SEP gates.
bs52_ward531_groups <- function(config, native_order = FALSE) {
  species <- config$species
  width <- config$cell_length
  radius <- max(0, species$max_seed_distance) + width / 2
  reach <- floor(radius / width) + 1L
  quarter <- expand.grid(dc = 0:reach, dr = seq_len(reach))
  quarter$distance <- width * sqrt(quarter$dr^2 + quarter$dc^2)
  quarter <- quarter[quarter$distance <= radius, ]
  by_distance <- if (native_order) bs52_v7_distance_order_cpp(quarter$distance) else
    order(quarter$distance)
  quarter <- quarter[by_distance, ]
  lapply(seq_len(nrow(species)), function(j) {
    effective <- species$effective_seed_distance[j]
    maximum <- species$max_seed_distance[j]
    if (effective == -1) return(NULL)
    d <- quarter[quarter$distance <= maximum + width / 2 * 1.414, ]
    list(dr = d$dr, dc = d$dc, p = bs52_ward_kernel(d$distance, effective, maximum, width))
  })
}

bs52_ward531_context <- function(ctx) {
  ctx$ward531_groups <- bs52_ward531_groups(ctx$config)
  coordinates <- ctx$targets
  raster_index <- (coordinates[, 1L] - 1L) * ctx$ncol + coordinates[, 2L]
  ctx$ward531_positions <- integer(ctx$nrow * ctx$ncol)
  ctx$ward531_positions[raster_index] <- seq_along(ctx$ids)
  ctx$ward531_order <- order(raster_index, method = "radix")
  ctx
}

bs52_ward531_recruit <- function(features, ecology, ctx, cell_state, state_eco,
                                 establishment, rng, deferred = FALSE, pending = NULL) {
  species <- as.character(ctx$config$species$species)
  mask <- integer(length(ctx$ids))
  if (NROW(establishment$planting)) {
    mask <- landisR_establishment_mask(mask, establishment, ctx$ids,
      cell_state, state_eco, species, ecology$sep)
  }
  strict <- if (deferred) rep(FALSE, length(species)) else
    seq_along(species) %in% landisR_serotiny_strict_index(ctx$config$species)
  groups <- if (deferred) ctx$ward_v7_groups else ctx$ward531_groups
  bs52_ward531_draw_cpp(features$mature, cell_state, state_eco, features$shade_class,
    ctx$targets, ctx$nrow, ctx$ncol, ctx$ward531_positions, ctx$ward531_order,
    ecology$sep, ecology$light, ecology$shade_tolerance,
    groups, strict, mask, match(establishment$excluded_cells, ctx$ids), rng$state, rng$offset,
    deferred, match(pending$cell_id, ctx$ids), as.integer(pending$species),
    as.logical(pending$serotiny), as.numeric(ctx$config$species$vegetative_probability))
}
