# ==============================================================================
# Serotiny API for landisR beta.3-dev
#
# New functionality relative to LANDIS-II Biomass Succession:
#   - facultative serotiny (LANDIS-II-compatible default);
#   - strict serotiny (no natural sexual recruitment without fire);
#   - post-fire serotinous seed dispersal through either local or Ward kernel
#     dispersal.
#
# Fire is represented as a structured event through landis_fire(). The function
# captures PRE-FIRE mature serotinous sources, applies an optional post-fire
# cohort update, then evaluates post-fire establishment and recruitment.
# ==============================================================================

landisR_normalize_serotiny_traits <- function(species) {
  if (!is.data.frame(species) || !"species" %in% names(species)) {
    stop("species doit etre un data.frame contenant la colonne species.",
         call. = FALSE)
  }

  n <- nrow(species)
  if (!n) return(species)

  has_post_fire <- "post_fire_regeneration" %in% names(species)
  post_fire <- if (has_post_fire) {
    tolower(trimws(as.character(species$post_fire_regeneration)))
  } else {
    rep("none", n)
  }
  post_fire[is.na(post_fire) | !nzchar(post_fire)] <- "none"

  default_type <- ifelse(post_fire == "serotiny", "facultative", "none")

  if (!"serotiny_type" %in% names(species)) {
    species$serotiny_type <- default_type
  } else {
    x <- tolower(trimws(as.character(species$serotiny_type)))
    missing <- is.na(x) | !nzchar(x)
    x[missing] <- default_type[missing]
    species$serotiny_type <- x
  }

  valid_type <- c("none", "facultative", "strict")
  if (any(!species$serotiny_type %in% valid_type)) {
    bad <- unique(species$serotiny_type[!species$serotiny_type %in% valid_type])
    stop(
      "serotiny_type invalide : ", paste(bad, collapse = ", "),
      ". Valeurs permises : none, facultative, strict.",
      call. = FALSE
    )
  }

  default_dispersal <- ifelse(species$serotiny_type == "none", "none", "local")

  if (!"serotiny_dispersal" %in% names(species)) {
    species$serotiny_dispersal <- default_dispersal
  } else {
    x <- tolower(trimws(as.character(species$serotiny_dispersal)))
    missing <- is.na(x) | !nzchar(x)
    x[missing] <- default_dispersal[missing]
    species$serotiny_dispersal <- x
  }

  valid_disp <- c("none", "local", "kernel")
  if (any(!species$serotiny_dispersal %in% valid_disp)) {
    bad <- unique(species$serotiny_dispersal[
      !species$serotiny_dispersal %in% valid_disp
    ])
    stop(
      "serotiny_dispersal invalide : ", paste(bad, collapse = ", "),
      ". Valeurs permises : none, local, kernel.",
      call. = FALSE
    )
  }

  if (any(species$serotiny_type == "none" &
          species$serotiny_dispersal != "none")) {
    stop(
      "Une espece avec serotiny_type='none' doit avoir ",
      "serotiny_dispersal='none'.",
      call. = FALSE
    )
  }

  if (any(species$serotiny_type != "none" &
          species$serotiny_dispersal == "none")) {
    stop(
      "Une espece serotineuse doit utiliser serotiny_dispersal='local' ",
      "ou 'kernel'.",
      call. = FALSE
    )
  }

  species
}

landisR_serotiny_strict_index <- function(species) {
  species <- landisR_normalize_serotiny_traits(species)
  which(species$serotiny_type == "strict")
}

landisR_serotiny_index <- function(species) {
  species <- landisR_normalize_serotiny_traits(species)
  which(species$serotiny_type != "none")
}

landisR_serotiny_rng_draw <- function(simulation, nrow, ncol) {
  if (!is.numeric(simulation$rng_state) || !length(simulation$rng_state)) {
    stop("Etat RNG absent de la simulation.", call. = FALSE)
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, .GlobalEnv)
    } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  assign(".Random.seed", simulation$rng_state, .GlobalEnv)

  draw <- matrix(
    stats::runif(nrow * ncol),
    nrow = nrow,
    ncol = ncol,
    byrow = TRUE
  )

  simulation$rng_state <- get(".Random.seed", .GlobalEnv, inherits = FALSE)

  list(simulation = simulation, draw = draw)
}

.landisR_legacy_prefire_270F <- function(simulation, cells) {
  config <- simulation$ctx$config
  species <- landisR_normalize_serotiny_traits(config$species)
  sero <- which(species$serotiny_type != "none")

  out <- setNames(vector("list", length(sero)), species$species[sero])
  if (!length(sero) || !length(cells)) return(out)

  active_ids <- as.integer(simulation$ctx$ids)
  pos <- match(cells, active_ids)
  if (anyNA(pos)) {
    stop("cells contient une cellule absente ou inactive.", call. = FALSE)
  }

  states <- .biomass_core$services$state_list(simulation)
  maturity <- setNames(as.numeric(species$maturity), species$species)

  for (j in seq_along(sero)) {
    sp_index <- sero[j]
    sp <- species$species[sp_index]
    mature_cells <- logical(length(cells))

    for (k in seq_along(cells)) {
      state <- states[[simulation$cell_state[pos[k]]]]
      if (nrow(state)) {
        mature_cells[k] <- any(
          state$species == sp &
            state$age >= maturity[[sp]] &
            state$biomass > 0
        )
      }
    }

    out[[j]] <- as.integer(cells[mature_cells])
  }

  out
}

landisR_serotiny_arrival <- function(simulation, source_cells, species_row) {
  active_ids <- as.integer(simulation$ctx$ids)
  arrival <- numeric(length(active_ids))

  if (!length(source_cells)) return(arrival)

  mode <- as.character(species_row$serotiny_dispersal)

  if (identical(mode, "local")) {
    pos <- match(source_cells, active_ids)
    pos <- pos[!is.na(pos)]
    arrival[pos] <- 1
    return(arrival)
  }

  if (!identical(mode, "kernel")) {
    stop("Mode de dispersion serotineuse inconnu : ", mode, call. = FALSE)
  }

  source_mask <- matrix(
    FALSE,
    nrow = simulation$ctx$nrow,
    ncol = simulation$ctx$ncol
  )
  source_mask[source_cells] <- TRUE

  probability_function <- if (exists(
    "bs52_ward_probability_rcpp",
    mode = "function",
    inherits = TRUE
  )) {
    bs52_ward_probability_rcpp
  } else {
    bs52_ward_probability
  }

  probability_function(
    source_mask,
    simulation$ctx$active,
    simulation$ctx$targets,
    species_row$effective_seed_distance,
    species_row$max_seed_distance,
    simulation$ctx$config$cell_length,
    simulation$ctx$neighborhood,
    simulation$ward_variant
  )
}

.landisR_legacy_recruit_270F <- function(simulation,
                                                prefire_sources,
                                                stochastic = TRUE) {
  config <- simulation$ctx$config
  species <- landisR_normalize_serotiny_traits(config$species)
  sero <- which(species$serotiny_type != "none")

  source_count <- integer(length(sero))
  names(source_count) <- species$species[sero]

  if (!length(sero)) {
    return(list(
      simulation = simulation,
      source_count = source_count,
      recruited_cells = 0L,
      recruited_cohorts = 0L,
      recruitment = data.frame()
    ))
  }

  source_count[] <- vapply(
    species$species[sero],
    function(sp) length(prefire_sources[[sp]]),
    integer(1L)
  )

  active_sero <- sero[source_count > 0L]
  if (!length(active_sero)) {
    return(list(
      simulation = simulation,
      source_count = source_count,
      recruited_cells = 0L,
      recruited_cohorts = 0L,
      recruitment = data.frame()
    ))
  }

  states <- .biomass_core$services$state_list(simulation)

  # This helper can return unmasked SEP x light probabilities for a fire event.
  features <- bs52_large_state_features52(
    states,
    simulation$state_eco,
    config,
    simulation$year,
    regular_recruitment = FALSE
  )

  n_sites <- length(simulation$ctx$ids)
  n_sero <- length(active_sero)

  arrival <- matrix(
    0,
    nrow = n_sites,
    ncol = n_sero,
    dimnames = list(NULL, species$species[active_sero])
  )

  for (j in seq_along(active_sero)) {
    sp_index <- active_sero[j]
    sp <- species$species[sp_index]
    arrival[, j] <- landisR_serotiny_arrival(
      simulation,
      prefire_sources[[sp]],
      species[sp_index, , drop = FALSE]
    )
  }

  base <- features$base_probability[
    simulation$cell_state,
    active_sero,
    drop = FALSE
  ]

  # Keep the probability object explicitly two-dimensional.
  # pmin()/pmax() with scalar first arguments can drop matrix dimensions,
  # which breaks arr.ind indexing when only one serotinous species is active.
  if (!is.matrix(base) || !identical(dim(base), dim(arrival))) {
    stop(
      "Dimensions incompatibles entre probabilite d'etablissement et arrivee serotineuse.",
      call. = FALSE
    )
  }

  probability <- base * arrival
  probability[!is.finite(probability)] <- 0
  probability[probability < 0] <- 0
  probability[probability > 1] <- 1

  # Defensive guarantee for the single-species case.
  dim(probability) <- dim(arrival)
  dimnames(probability) <- dimnames(arrival)

  if (isTRUE(stochastic)) {
    rng <- landisR_serotiny_rng_draw(simulation, n_sites, n_sero)
    simulation <- rng$simulation
    recruited <- rng$draw < probability
  } else {
    recruited <- probability > 0
  }

  pairs <- which(recruited, arr.ind = TRUE)

  if (!nrow(pairs)) {
    return(list(
      simulation = simulation,
      source_count = source_count,
      recruited_cells = 0L,
      recruited_cohorts = 0L,
      recruitment = data.frame()
    ))
  }

  active_ids <- as.integer(simulation$ctx$ids)
  cell_groups <- split(pairs[, 2L], pairs[, 1L])
  replacement <- vector("list", length(cell_groups))
  record <- vector("list", nrow(pairs))
  record_i <- 0L

  params_cache <- new.env(parent = emptyenv())

  get_params <- function(eco) {
    key <- as.character(eco)
    if (!exists(key, params_cache, inherits = FALSE)) {
      assign(
        key,
        bs52_dynamic_at(config, simulation$year, eco),
        envir = params_cache
      )
    }
    get(key, params_cache, inherits = FALSE)
  }

  group_names <- names(cell_groups)
  for (g in seq_along(cell_groups)) {
    active_pos <- as.integer(group_names[g])
    cell_id <- active_ids[active_pos]
    state_id <- simulation$cell_state[active_pos]
    x <- states[[state_id]]
    eco <- as.integer(simulation$ctx$active_eco[active_pos])
    p <- get_params(eco)
    eco_max <- max(p$max_biomass)

    selected_existing <- x$age >= config$timestep & x$age > 1
    site_biomass <- sum(x$biomass[selected_existing])

    recruit_columns <- cell_groups[[g]]
    recruit_species_index <- active_sero[recruit_columns]

    for (m in seq_along(recruit_species_index)) {
      sp_index <- recruit_species_index[m]
      sp <- species$species[sp_index]
      pi <- match(sp, p$species)
      max_anpp <- p$max_anpp[pi]
      biomass <- trunc(max(
        2,
        min(
          max_anpp,
          max_anpp * exp(-1.6 * site_biomass / eco_max)
        )
      ))

      existing_age1 <- which(x$species == sp & x$age == 1)
      if (length(existing_age1)) {
        x$biomass[existing_age1[1L]] <- x$biomass[existing_age1[1L]] + biomass
      } else {
        x <- rbind(
          x,
          data.frame(
            species = sp,
            age = 1,
            biomass = biomass,
            stringsAsFactors = FALSE
          )
        )
      }

      record_i <- record_i + 1L
      record[[record_i]] <- data.frame(
        cell_id = cell_id,
        species = sp,
        seed_probability = arrival[active_pos, recruit_columns[m]],
        establishment_probability = probability[active_pos, recruit_columns[m]],
        biomass = biomass,
        stringsAsFactors = FALSE
      )
    }

    x <- x[order(match(x$species, species$species), x$age), , drop = FALSE]
    rownames(x) <- NULL
    replacement[[g]] <- transform(x, cell_id = cell_id)[
      c("cell_id", "species", "age", "biomass")
    ]
  }

  replacement <- do.call(rbind, replacement)
  rownames(replacement) <- NULL

  simulation <- .biomass_core$services$apply_cohorts(
    simulation,
    cohorts = replacement,
    clear_cells = integer()
  )

  recruitment <- do.call(rbind, record[seq_len(record_i)])
  rownames(recruitment) <- NULL

  list(
    simulation = simulation,
    source_count = source_count,
    recruited_cells = length(unique(recruitment$cell_id)),
    recruited_cohorts = nrow(recruitment),
    recruitment = recruitment
  )
}

