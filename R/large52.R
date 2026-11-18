# ==============================================================================
# large52.R
# LANDIS-R Biomass Succession v5.2
#
# Large-landscape production engine, v3.6 CSR + bitset Ward + exact age-mortality lookup.
#
# This file is an OVERLAY. It does not replace the validated reference engine.
# Source the usual modules first, including production52.R, fast52.R, rcpp52.R,
# and growth52.cpp through bs52_init_rcpp().
#
# Design:
#   - one cohort table per UNIQUE ecological/cohort state, not per raster cell;
#   - active cells keep only an integer state_id;
#   - deterministic growth is evaluated once per unique state;
#   - Ward dispersal is still evaluated at the active-cell level using the
#     already validated Rcpp Ward probability function;
#   - stochastic establishment draws remain cell-major, species-major;
#   - compact state-compressed snapshots are written one at a time;
#   - checkpoints save state tables + cell_state + RNG state;
#   - no recruitment table is built when save_recruitment = FALSE.
#
# v1 scope:
#   - living biomass only;
#   - track_dead_pools must be FALSE;
#   - save_recruitment must be FALSE.
#
# The biological equations are unchanged. Exact equivalence must be verified
# with 33_valider_moteur_large_5km.R before full-landscape production.
# ==============================================================================

bs52_large_stop <- function(...) stop(..., call. = FALSE)

bs52_large_progress_enabled52 <- function() {
  isTRUE(
    getOption(
      "landisR52.progress",
      FALSE
    )
  )
}



bs52_large_n_cohorts52 <- function(states) {
  if (!is.list(states) ||
      is.null(states$species_id)) {
    bs52_large_stop(
      "Store d'etats CSR v3 invalide."
    )
  }

  as.integer(
    length(states$species_id)
  )
}


bs52_large_integer <- function(x, label, minimum = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != trunc(x) || x < minimum) {
    bs52_large_stop(label, " doit etre un entier >= ", minimum, ".")
  }
  as.integer(x)
}

bs52_large_elapsed <- function(start) {
  proc.time()[["elapsed"]] - start
}


bs52_large_input_signature52 <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !file.info(files)$isdir]
  files <- sort(normalizePath(files, winslash = "/", mustWork = TRUE))
  info <- file.info(files)

  data.frame(
    path = substring(files, nchar(root) + 2L),
    size = as.numeric(info$size),
    mtime = as.numeric(info$mtime),
    stringsAsFactors = FALSE
  )
}

bs52_large_read_canonical_cached52 <- function(inputs, cache_file = NULL) {
  inputs <- normalizePath(inputs, winslash = "/", mustWork = TRUE)
  if (is.null(cache_file)) return(bs52_read_canonical52(inputs))

  if (!is.character(cache_file) || length(cache_file) != 1L ||
      is.na(cache_file) || !nzchar(cache_file))
    bs52_large_stop("Le chemin du cache canonique doit etre une chaine non vide.")

  cache_file <- path.expand(cache_file)
  cache_file <- file.path(
    normalizePath(dirname(cache_file), winslash = "/", mustWork = FALSE),
    basename(cache_file)
  )
  if (dir.exists(cache_file))
    bs52_large_stop("Le cache canonique doit etre un fichier : ", cache_file)

  input_signature <- function() {
    files <- sort(list.files(inputs, pattern = "[.]csv$", ignore.case = TRUE,
                             recursive = TRUE, full.names = TRUE), method = "radix")
    paths <- normalizePath(files, winslash = "/", mustWork = TRUE)
    target <- normalizePath(cache_file, winslash = "/", mustWork = FALSE)
    if (.Platform$OS.type == "windows") {
      paths <- tolower(paths)
      target <- tolower(target)
    }
    if (target %in% paths)
      bs52_large_stop("Le cache canonique ne peut pas remplacer un CSV d'entree.")
    checksums <- unname(tools::md5sum(files))
    if (anyNA(checksums))
      bs52_large_stop("Impossible de lire les CSV pour verifier le cache canonique.")
    setNames(checksums, substring(files, nchar(inputs) + 2L))
  }

  signature <- input_signature()
  format <- "landisR52/canonical-cache-v2"
  reader_version <- as.character(getNamespaceVersion("landisR.succession.biomass"))

  if (file.exists(cache_file)) {
    cache <- tryCatch(readRDS(cache_file), error = function(e) NULL,
                      warning = function(w) NULL)
    if (is.list(cache) && identical(cache$format, format) &&
        identical(cache$source, inputs) &&
        identical(cache$reader_version, reader_version) &&
        identical(cache$signature, signature) && is.list(cache$canonical)) {
      message("Cache canonique reutilise : ", cache_file)
      return(cache$canonical)
    }
    rm(cache)
  }

  message("Lecture CSV canonique et creation du cache binaire...")
  canonical <- bs52_read_canonical52(inputs)
  if (!identical(signature, input_signature()))
    bs52_large_stop("Les CSV canoniques ont change pendant la lecture. Relancer l'import.")

  .biomass_core$services$write_rds(
    list(format = format, source = inputs, reader_version = reader_version,
         signature = signature, canonical = canonical),
    cache_file, compress = "none", overwrite = TRUE
  )

  message("Cache canonique cree : ", cache_file)
  canonical
}


bs52_large_function_body_signature52 <- function(fun) {
  paste(
    deparse(
      body(
        fun
      ),
      width.cutoff = 500L
    ),
    collapse = "\n"
  )
}

bs52_large_context_code_signature52 <- function() {
  list(
    large_context =
      bs52_large_function_body_signature52(
        bs52_large_context52
      ),
    internalize_ecoregions =
      bs52_large_function_body_signature52(
        bs52_production_internalize_ecoregions
      ),
    initialize_explicit =
      bs52_large_function_body_signature52(
        bs52_initialize_explicit
      ),
    spinup_live =
      bs52_large_function_body_signature52(
        bs52_spinup_live
      ),
    ward_neighborhood =
      bs52_large_function_body_signature52(
        bs52_ward_neighborhood
      ),
    ward_table =
      bs52_large_function_body_signature52(
        bs52_large_ward_table52
      )
  )
}

bs52_large_climate_signature52 <- function(climate) {
  dynamic_file <- if (is.null(
        climate$dynamic_file
      )) {
    NULL
  } else {
    normalizePath(
      climate$dynamic_file,
      winslash = "/",
      mustWork = TRUE
    )
  }

  list(
    name = if (is.null(
          climate$name
        )) {
      NA_character_
    } else {
      as.character(
        climate$name
      )
    },
    dynamic_file =
      dynamic_file,
    dynamic_md5 =
      if (is.null(
            dynamic_file
          )) {
        NA_character_
      } else {
        unname(
          tools::md5sum(
            dynamic_file
          )
        )
      }
  )
}

bs52_large_context_signature52 <- function(
    canonical,
    inputs,
    inputs_signature,
    climate) {

  explicit_initial <-
    !is.null(
      canonical$explicit_initial
    )

  list(
    format =
      "landisR52/context-cache-signature-v1",
    source =
      normalizePath(
        inputs,
        winslash = "/",
        mustWork = TRUE
      ),
    inputs_signature =
      inputs_signature,
    initialization =
      if (explicit_initial) {
        "explicit-initial-static-across-climate"
      } else {
        "spinup-climate-specific"
      },
    climate =
      if (explicit_initial) {
        NULL
      } else {
        bs52_large_climate_signature52(
          climate
        )
      },
    code =
      bs52_large_context_code_signature52()
  )
}

bs52_large_refresh_dynamic52 <- function(
    ctx,
    canonical) {

  if (!is.list(ctx) ||
      is.null(ctx$config) ||
      is.null(ctx$ecoregion_code)) {
    bs52_large_stop(
      "Contexte invalide pour le rafraichissement climatique."
    )
  }

  if (is.null(
        canonical$dynamic
      )) {
    bs52_large_stop(
      "Table dynamique absente des entrees canoniques."
    )
  }

  code <-
    ctx$ecoregion_code

  labels <-
    names(
      code
    )

  dynamic <-
    canonical$dynamic

  dynamic_labels <-
    as.character(
      dynamic$ecoregion
    )

  if (!all(
        labels %in%
          dynamic_labels
      )) {
    bs52_large_stop(
      "Parametres dynamiques incomplets pour une ecoregion active."
    )
  }

  keep <-
    dynamic_labels %in%
      labels

  dynamic <-
    dynamic[
      keep,
      ,
      drop = FALSE
    ]

  mapped <-
    unname(
      code[
        as.character(
          dynamic$ecoregion
        )
      ]
    )

  if (anyNA(
        mapped
      )) {
    bs52_large_stop(
      "Codage climatique interne invalide."
    )
  }

  dynamic$ecoregion <-
    mapped

  ctx$config$dynamic <-
    dynamic

  ctx
}

bs52_large_derive_context_cache_file52 <- function(
    canonical_cache_file) {

  if (is.null(
        canonical_cache_file
      )) {
    return(
      NULL
    )
  }

  path <-
    normalizePath(
      canonical_cache_file,
      winslash = "/",
      mustWork = FALSE
    )

  name <-
    basename(
      path
    )

  if (grepl(
        "_canonical\\.rds$",
        name
      )) {
    name <-
      sub(
        "_canonical\\.rds$",
        "_context.rds",
        name
      )
  } else if (grepl(
               "\\.rds$",
               name
             )) {
    name <-
      sub(
        "\\.rds$",
        "_context.rds",
        name
      )
  } else {
    name <-
      paste0(
        name,
        "_context.rds"
      )
  }

  file.path(
    dirname(
      path
    ),
    name
  )
}

bs52_large_read_context_cached52 <- function(
    canonical,
    inputs,
    inputs_signature,
    climate,
    cache_file = NULL) {

  if (is.null(
        cache_file
      )) {
    return(
      bs52_large_context52(
        canonical
      )
    )
  }

  cache_file <-
    normalizePath(
      cache_file,
      winslash = "/",
      mustWork = FALSE
    )

  signature <-
    bs52_large_context_signature52(
      canonical = canonical,
      inputs = inputs,
      inputs_signature = inputs_signature,
      climate = climate
    )

  if (file.exists(
        cache_file
      )) {

    cache <-
      tryCatch(
        readRDS(
          cache_file
        ),
        error = function(e) {
          NULL
        }
      )

    if (is.list(cache) &&
        identical(
          cache$format,
          "landisR52/context-cache-v1"
        ) &&
        identical(
          cache$signature,
          signature
        ) &&
        is.list(
          cache$context
        )) {

      ctx <-
        bs52_large_refresh_dynamic52(
          cache$context,
          canonical
        )

      message(
        "Cache contexte reutilise : ",
        cache_file
      )

      return(
        ctx
      )
    }

    message(
      "Cache contexte obsolete; reconstruction."
    )
  }

  message(
    "Construction du contexte compact et creation du cache..."
  )

  ctx <-
    bs52_large_context52(
      canonical
    )

  ctx <-
    bs52_large_refresh_dynamic52(
      ctx,
      canonical
    )

  .biomass_core$services$write_rds(
    list(
      format =
        "landisR52/context-cache-v1",
      signature =
        signature,
      context =
        ctx
    ),
    cache_file,
    compress = "gzip1",
    overwrite = TRUE
  )

  message(
    "Cache contexte cree : ",
    cache_file
  )

  ctx
}


bs52_large_ward_table52 <- function(config, neighborhood,
                                    variant = c("ward2019", "ward2022")) {
  variant <- match.arg(variant)

  species <- config$species
  n_species <- nrow(species)
  n_neighbors <- nrow(neighborhood)

  log_failure <- matrix(
    0,
    nrow = n_neighbors,
    ncol = n_species
  )

  omnidisperse <- logical(
    n_species
  )

  for (j in seq_len(n_species)) {
    effective <- species$effective_seed_distance[j]
    maximum <- species$max_seed_distance[j]

    if (identical(
          as.numeric(effective),
          -1
        )) {
      omnidisperse[j] <- TRUE
      next
    }

    probability <- bs52_ward_kernel(
      neighborhood$distance,
      effective,
      maximum,
      config$cell_length
    )

    if (identical(
          variant,
          "ward2019"
        )) {
      support_limit <- bs52_ward2019_support_limit(
        maximum,
        config$cell_length
      )
      probability[
        neighborhood$distance > support_limit
      ] <- 0
    }

    active <- probability > 0

    log_failure[active, j] <- log1p(
      -probability[active]
    )
  }

  list(
    log_failure = log_failure,
    omnidisperse = omnidisperse
  )
}


bs52_large_context52 <- function(canonical) {
  if (!is.list(canonical) ||
      is.null(canonical$landscape) ||
      is.null(canonical$landscape$cells)) {
    bs52_large_stop("Entrees canoniques invalides.")
  }

  cells <- canonical$landscape$cells
  active <- as.logical(canonical$landscape$active)

  if (length(active) != nrow(cells) || anyNA(active)) {
    bs52_large_stop("Etat actif canonique invalide.")
  }

  internal <- bs52_production_internalize_ecoregions(
    canonical,
    cells,
    active
  )

  config <- internal$config

  eco_code <- integer(nrow(cells))
  eco_code[active] <- unname(
    internal$code[as.character(cells$ecoregion[active])]
  )

  if (anyNA(eco_code[active])) {
    bs52_large_stop("Codage interne des ecoregions invalide.")
  }

  eco <- matrix(
    eco_code,
    nrow = canonical$landscape$nrow,
    ncol = canonical$landscape$ncol
  )

  ids <- as.integer(cells$cell_id[active])
  active_eco <- as.integer(eco[ids])
  active_community <- as.integer(cells$community[active])

  allowed_communities <- if (!is.null(config$explicit_initial)) {
    as.integer(config$explicit_initial$community_codes)
  } else {
    as.integer(config$community_codes)
  }

  if (!length(ids) ||
      anyNA(active_eco) ||
      anyNA(active_community) ||
      any(active_eco <= 0L) ||
      any(!active_community %in% allowed_communities)) {
    bs52_large_stop("Cellules actives invalides.")
  }

  # LANDIS initial-community MapCode 0 is legitimate when it is explicitly
  # declared in the initial-community table. The canonical reader already
  # validates active cells against community_codes, so do not impose > 0 here.
  #
  # Build a collision-free numeric key that also supports MapCode 0.
  community_min <- min(active_community)
  community_span <- max(active_community) - community_min + 1

  initial_key <- as.double(active_eco - 1L) * community_span +
    as.double(active_community - community_min)

  message(
    "Large52 contexte : ",
    format(length(ids), big.mark = " "),
    " cellules actives; MapCode initial min=",
    community_min,
    ", max=",
    max(active_community),
    ", n(MapCode 0)=",
    sum(active_community == 0L),
    "."
  )

  first <- !duplicated(initial_key)
  unique_key <- initial_key[first]
  cell_state <- match(initial_key, unique_key)

  state_eco <- active_eco[first]
  state_community <- active_community[first]

  states <- if (!is.null(config$explicit_initial)) {
    bs52_initialize_explicit_states(config$explicit_initial, state_community)
  } else {
    lapply(seq_along(state_community), function(s) {
      state <- bs52_spinup_live(config, state_community[s], state_eco[s])
      rownames(state) <- NULL
      state
    })
  }

  bs52_large_init_bulk52()

  states <- .biomass_core$services$encode_csr(
    states,
    as.character(config$species$species)
  )

  active_matrix <- eco > 0L
  targets <- arrayInd(ids, dim(eco))
  storage.mode(targets) <- "integer"

  neighborhood <- bs52_ward_neighborhood(
    config$cell_length,
    max(config$species$max_seed_distance)
  )

  ward_tables <- list(
    ward2019 = bs52_large_ward_table52(
      config,
      neighborhood,
      "ward2019"
    ),
    ward2022 = bs52_large_ward_table52(
      config,
      neighborhood,
      "ward2022"
    )
  )

  list(
    config = config,
    eco = eco,
    active = active_matrix,
    ids = ids,
    targets = targets,
    active_eco = active_eco,
    active_community = active_community,
    cell_state = as.integer(cell_state),
    states = states,
    state_eco = as.integer(state_eco),
    state_community = as.integer(state_community),
    ecoregion_code = internal$code,
    nrow = canonical$landscape$nrow,
    ncol = canonical$landscape$ncol,
    neighborhood = neighborhood,
    ward_tables = ward_tables
  )
}

bs52_large_params52 <- function(ctx, duration) {
  ecoregions <- as.character(sort(unique(ctx$active_eco)))

  lapply(seq_len(duration), function(year) {
    values <- lapply(
      ecoregions,
      function(eco) {
        bs52_dynamic_at(
          ctx$config,
          year,
          as.integer(eco)
        )
      }
    )
    setNames(values, ecoregions)
  })
}

bs52_large_block_params_by_eco52 <- function(params, block_years,
                                               state_eco) {
  eco_values <- sort(unique(as.integer(state_eco)))

  if (!length(eco_values) || anyNA(eco_values) || any(eco_values < 1L)) {
    bs52_large_stop("Ecoregions d'etat invalides.")
  }

  out <- vector(
    "list",
    max(eco_values)
  )

  for (eco in eco_values) {
    key <- as.character(eco)
    out[[eco]] <- lapply(
      block_years,
      function(year) {
        value <- params[[year]][[key]]
        if (is.null(value)) {
          bs52_large_stop(
            "Parametres climatiques absents pour ecoregion ",
            eco,
            " annee ",
            year,
            "."
          )
        }
        value
      }
    )
  }

  out
}

bs52_large_grow_states52 <- function(states, state_eco, block_years,
                                     params, ctx,
                                     engine = c("Rcpp", "R"),
                                     progress = bs52_large_progress_enabled52(),
                                     progress_label = "Croissance CSR + features",
                                     feature_year = tail(block_years, 1L), report = NULL) {
  engine <- match.arg(engine)

  if (!identical(engine, "Rcpp")) {
    bs52_large_stop(
      "large52 v3 compact exige engine='Rcpp'."
    )
  }

  bs52_large_init_bulk52()

  params_by_eco <- bs52_large_block_params_by_eco52(
    params,
    block_years,
    state_eco
  )

  ecology <- bs52_large_ecology_tables52(
    ctx$config,
    feature_year,
    state_eco
  )

  out <- landisR::landis_progress_native(report, 0.05, 0.34, bs52_grow_csr_indexed_features_cpp(
    states,
    as.integer(state_eco),
    params_by_eco,
    as.character(ctx$config$species$species),
    ecology$maturity,
    as.integer(ctx$config$timestep),
    ecology$thresholds,
    ecology$sep,
    ecology$eco_max,
    ecology$shade_tolerance,
    ecology$light,
    combine_last = tail(block_years, 1L) %% ctx$config$timestep == 0L,
    spinup_fraction = 0,
    spinup = FALSE,
    show_progress = progress,
    progress_label = progress_label
  ))

  base_probability <- out$base_probability

  # Strictly serotinous species cannot recruit through the ordinary
  # succession-step sexual reproduction pathway.
  strict <- landisR_serotiny_strict_index(ctx$config$species)
  if (length(strict)) {
    base_probability[, strict] <- 0
  }

  list(
    states = out$states,
    features = list(
      mature = out$mature,
      establishment_biomass = as.numeric(
        out$establishment_biomass
      ),
      base_probability = base_probability,
      shade_class = out$shade_class
    ),
    ecology = ecology,
    fused = TRUE,
    weight_lookup_fallbacks = as.numeric(
      out$weight_lookup_fallbacks
    )
  )
}

bs52_large_ecology_tables52 <- function(config, year, state_eco) {
  species_names <- as.character(config$species$species)
  n_species <- length(species_names)
  eco_values <- sort(unique(as.integer(state_eco)))

  if (!length(eco_values) || anyNA(eco_values) || any(eco_values < 1L)) {
    bs52_large_stop("Ecoregions d'etat invalides.")
  }

  n_eco <- max(eco_values)
  threshold_count <- nrow(config$shade)

  thresholds <- matrix(
    NA_real_,
    nrow = n_eco,
    ncol = threshold_count
  )

  sep <- matrix(
    NA_real_,
    nrow = n_eco,
    ncol = n_species
  )

  max_anpp <- matrix(
    NA_real_,
    nrow = n_eco,
    ncol = n_species
  )

  eco_max <- rep(
    NA_real_,
    n_eco
  )

  for (eco in eco_values) {
    p <- bs52_dynamic_at(
      config,
      year,
      eco
    )

    if (!identical(as.character(p$species), species_names)) {
      bs52_large_stop(
        "Ordre des especes incompatible dans les parametres dynamiques."
      )
    }

    threshold <- config$shade[[as.character(eco)]][
      order(config$shade$shade)
    ]

    if (length(threshold) != threshold_count) {
      bs52_large_stop(
        "Nombre de seuils d'ombrage incompatible."
      )
    }

    thresholds[eco, ] <- as.numeric(threshold)
    sep[eco, ] <- as.numeric(p$sep)
    max_anpp[eco, ] <- as.numeric(p$max_anpp)
    eco_max[eco] <- max(p$max_biomass)
  }

  max_tolerance <- max(config$light$shade_tolerance)

  light <- matrix(
    NA_real_,
    nrow = max_tolerance,
    ncol = 6L
  )

  light[config$light$shade_tolerance, ] <- as.matrix(
    config$light[
      ,
      paste0("shade_", 0:5),
      drop = FALSE
    ]
  )

  list(
    thresholds = thresholds,
    sep = sep,
    max_anpp = max_anpp,
    eco_max = eco_max,
    shade_tolerance = as.integer(
      config$species$shade_tolerance[
        match(species_names, config$species$species)
      ]
    ),
    maturity = as.numeric(
      config$species$maturity[
        match(species_names, config$species$species)
      ]
    ),
    light = light
  )
}

bs52_large_state_features52 <- function(states, state_eco, config, year,
                                        need_output_biomass = FALSE,
                                        regular_recruitment = TRUE) {
  species_names <- as.character(config$species$species)
  n_species <- length(species_names)
  n_states <- length(states)

  maturity <- setNames(
    as.numeric(config$species$maturity),
    species_names
  )

  mature <- matrix(FALSE, n_states, n_species)
  shade_biomass <- numeric(n_states)
  establishment_biomass <- numeric(n_states)
  max_age <- numeric(n_states)

  state_biomass <- if (isTRUE(need_output_biomass)) {
    matrix(
      0,
      nrow = n_states,
      ncol = n_species,
      dimnames = list(NULL, species_names)
    )
  } else {
    NULL
  }

  for (s in seq_len(n_states)) {
    x <- states[[s]]

    if (!nrow(x)) {
      next
    }

    sp <- as.character(x$species)
    age <- as.numeric(x$age)
    biomass <- as.numeric(x$biomass)

    sp_index <- match(sp, species_names)

    if (anyNA(sp_index)) {
      bs52_large_stop("Espece de cohorte absente du registre.")
    }

    shade_biomass[s] <- sum(biomass[age > 5])
    establishment_biomass[s] <- sum(
      biomass[
        age >= config$timestep &
        age > 1
      ]
    )
    max_age[s] <- max(age)

    good <- biomass > 0 &
      age >= unname(maturity[sp])

    if (any(good)) {
      mature[
        s,
        unique(sp_index[good])
      ] <- TRUE
    }

    if (isTRUE(need_output_biomass)) {
      totals <- rowsum(
        biomass,
        group = sp_index,
        reorder = FALSE
      )

      state_biomass[
        s,
        as.integer(rownames(totals))
      ] <- totals[, 1L]
    }
  }

  params_by_eco <- lapply(
    sort(unique(state_eco)),
    function(eco) {
      bs52_dynamic_at(
        config,
        year,
        eco
      )
    }
  )
  names(params_by_eco) <- as.character(sort(unique(state_eco)))

  base_probability <- matrix(
    0,
    nrow = n_states,
    ncol = n_species
  )

  eco_max <- numeric(n_states)
  shade_class <- integer(n_states)

  for (s in seq_len(n_states)) {
    eco <- state_eco[s]
    p <- params_by_eco[[as.character(eco)]]

    if (!identical(as.character(p$species), species_names)) {
      bs52_large_stop(
        "Ordre des especes incompatible dans les parametres dynamiques."
      )
    }

    threshold <- config$shade[[as.character(eco)]][
      order(config$shade$shade)
    ]

    current_eco_max <- max(p$max_biomass)
    eco_max[s] <- current_eco_max

    actual <- min(
      current_eco_max,
      shade_biomass[s]
    )

    relative <- actual / current_eco_max
    shade_classes <- which(relative >= threshold)
    shade <- if (length(shade_classes)) {
      max(shade_classes)
    } else {
      0L
    }

    shade_class[s] <- shade

    light_rows <- match(
      p$shade_tolerance,
      config$light$shade_tolerance
    )

    if (anyNA(light_rows)) {
      bs52_large_stop("Ligne de lumiere suffisante absente.")
    }

    light <- as.numeric(
      config$light[
        light_rows,
        paste0("shade_", shade)
      ]
    )

    base_probability[s, ] <- pmin(
      1,
      pmax(
        0,
        p$sep * light
      )
    )
  }

  if (isTRUE(regular_recruitment)) {
    strict <- landisR_serotiny_strict_index(config$species)
    if (length(strict)) {
      base_probability[, strict] <- 0
    }
  }

  list(
    mature = mature,
    shade_biomass = shade_biomass,
    establishment_biomass = establishment_biomass,
    max_age = max_age,
    state_biomass = state_biomass,
    params_by_eco = params_by_eco,
    base_probability = base_probability,
    eco_max = eco_max,
    shade_class = shade_class
  )
}

bs52_large_masks_from_probability52 <- function(probability,
                                                stochastic = TRUE,
                                                chunk_cells = 100000L, report = NULL) {
  if (!is.matrix(probability) || !is.numeric(probability)) {
    bs52_large_stop("Matrice de probabilite invalide.")
  }

  n_sites <- nrow(probability)
  n_species <- ncol(probability)

  if (n_species < 1L || n_species > 30L) {
    bs52_large_stop(
      "Le codage compact des recrutements prend en charge 1 a 30 especes."
    )
  }

  powers <- 2^(seq_len(n_species) - 1L)
  masks <- integer(n_sites)

  if (!isTRUE(stochastic)) {
    recruited <- probability > 0
    masks[] <- as.integer(
      drop(recruited %*% powers)
    )
    return(masks)
  }

  chunk_cells <- bs52_large_integer(
    chunk_cells,
    "chunk_cells",
    1L
  )

  starts <- seq.int(
    1L,
    n_sites,
    by = chunk_cells
  )

  for (start in starts) {
    if (!is.null(report))
      report(0.64 + 0.07 * (start - 1) / n_sites, "Recruitment draws")
    end <- min(
      n_sites,
      start + chunk_cells - 1L
    )

    idx <- seq.int(start, end)
    n <- length(idx)

    # One contiguous RNG draw per cell chunk. With byrow=TRUE, the resulting
    # order is cell-major then species-major, exactly like repeated
    # runif(n_species) calls for successive cells.
    draw <- matrix(
      stats::runif(n * n_species),
      nrow = n,
      ncol = n_species,
      byrow = TRUE
    )

    recruited <- draw <
      probability[idx, , drop = FALSE]

    masks[idx] <- as.integer(
      drop(recruited %*% powers)
    )
  }

  masks
}

bs52_large_append_pattern52 <- function(state, mask, species_names,
                                        params, establishment_biomass,
                                        eco_max) {
  mask <- as.integer(mask)

  if (mask == 0L) {
    rownames(state) <- NULL
    return(state)
  }

  bit <- bitwShiftL(
    1L,
    seq_along(species_names) - 1L
  )

  selected <- which(
    bitwAnd(mask, bit) != 0L
  )

  if (!length(selected)) {
    rownames(state) <- NULL
    return(state)
  }

  anpp <- as.numeric(params$max_anpp[selected])

  candidate <- anpp * exp(
    -1.6 * establishment_biomass / eco_max
  )

  biomass <- trunc(
    pmax(
      2,
      pmin(
        anpp,
        candidate
      )
    )
  )

  recruits <- data.frame(
    species = as.character(params$species[selected]),
    age = rep(1, length(selected)),
    biomass = biomass,
    stringsAsFactors = FALSE
  )

  out <- rbind(
    state,
    recruits
  )

  rownames(out) <- NULL
  out
}

bs52_large_recruit52 <- function(states, state_eco, cell_state,
                                 ctx, year,
                                 variant = c("ward2019", "ward2022"),
                                 mode = c("ward", "local"),
                                 stochastic = TRUE,
                                 rng_chunk_cells = 100000L,
                                 progress = bs52_large_progress_enabled52(),
                                 precomputed_features = NULL,
                                 precomputed_ecology = NULL,
                                 establishment = NULL,
                                 landis_rng = NULL,
                                 deferred = FALSE, pending_regeneration = NULL, report = NULL) {
  variant <- match.arg(variant)
  mode <- match.arg(mode)

  bs52_large_init_bulk52()

  config <- ctx$config
  species_names <- as.character(config$species$species)
  n_species <- length(species_names)
  n_sites <- length(ctx$ids)

  if (length(cell_state) != n_sites ||
      anyNA(cell_state) ||
      any(cell_state < 1L | cell_state > .biomass_core$services$n_states(states))) {
    bs52_large_stop("cell_state invalide.")
  }

  if (is.null(precomputed_features) ||
      is.null(precomputed_ecology)) {
    bs52_large_stop(
      "large52 v3 exige les features pre-calcules par le kernel de croissance."
    )
  }

  features <- precomputed_features
  ecology <- precomputed_ecology
  feature_seconds <- 0

  ward_start <- proc.time()[["elapsed"]]
  if (!is.null(report)) report(0.35, "Preparing regeneration")

  if (is.null(landis_rng)) {
    ward_table <- ctx$ward_tables[[variant]]

    if (is.null(ward_table)) {
      bs52_large_stop(
        "Table Ward precompilee absente pour ",
        variant,
        "."
      )
    }

    probability <- landisR::landis_progress_native(report, 0.39, 0.63,
      bs52_ward_probability_fused_bitset_cpp(
      features$mature,
      as.integer(cell_state),
      ctx$targets,
      as.integer(ctx$nrow),
      as.integer(ctx$ncol),
      as.integer(ctx$neighborhood$row),
      as.integer(ctx$neighborhood$col),
      ward_table$log_failure,
      as.logical(ward_table$omnidisperse),
      features$base_probability,
      local_only = identical(
        mode,
        "local"
      ),
      ward2019 = identical(
        variant,
        "ward2019"
      ),
      show_progress = progress,
      progress_label = paste0(
        "Annee ",
        year,
        " - Ward bitset"
      )
    ))

    ward_seconds <- bs52_large_elapsed(
      ward_start
    )

    rng_start <- proc.time()[["elapsed"]]

    recruit_mask <- bs52_large_masks_from_probability52(
      probability,
      stochastic = stochastic,
      chunk_cells = rng_chunk_cells, report = report
    )

    if (!is.null(establishment)) {
      recruit_mask <- landisR_establishment_mask(
        recruit_mask, establishment, ctx$ids, cell_state, state_eco,
        species_names, ecology$sep
      )
    }

    rng_seconds <- bs52_large_elapsed(rng_start)

    rm(probability)
  } else {
    draw <- landisR::landis_progress_native(report, 0.39, 0.72,
      bs52_ward531_recruit(features, ecology, ctx, cell_state, state_eco,
                           establishment, landis_rng, deferred, pending_regeneration))
    recruit_mask <- draw$mask
    landis_rng <- list(state = draw$rng_state, offset = draw$rng_offset)
    ward_seconds <- bs52_large_elapsed(ward_start)
    rng_seconds <- 0
  }

  invisible(gc(FALSE))

  rebuild_start <- proc.time()[["elapsed"]]

  rebuilt <- landisR::landis_progress_native(report, 0.73, 0.91, bs52_rebuild_csr_cpp(
    states,
    as.integer(state_eco),
    as.integer(cell_state),
    as.integer(recruit_mask),
    as.integer(n_species),
    ecology$max_anpp,
    ecology$eco_max,
    as.numeric(features$establishment_biomass),
    show_progress = progress,
    progress_label = paste0(
      "Annee ",
      year,
      " - rebuild CSR"
    )
  ))

  rebuild_seconds <- bs52_large_elapsed(rebuild_start)

  list(
    states = rebuilt$states,
    state_eco = as.integer(rebuilt$state_eco),
    cell_state = as.integer(rebuilt$cell_state),
    n_states_before = .biomass_core$services$n_states(states),
    n_states_after = .biomass_core$services$n_states(rebuilt$states),
    n_recruited_cells = as.integer(rebuilt$n_recruited_cells),
    landis_rng = landis_rng,
    timing = data.frame(
      feature_seconds = feature_seconds,
      ward_seconds = ward_seconds,
      rng_seconds = rng_seconds,
      rebuild_seconds = rebuild_seconds
    )
  )
}


bs52_large_snapshot_summary52 <- function(snapshot) {
  if (!is.list(snapshot) ||
      !identical(snapshot$format, "landisR52/state-snapshot-v1")) {
    bs52_large_stop("Snapshot compact invalide.")
  }

  n_states <- nrow(snapshot$state_biomass)
  counts <- tabulate(
    snapshot$state_id,
    nbins = n_states
  )

  used <- counts > 0L
  n_cells <- length(snapshot$cell_id)

  biomass <- cbind(
    snapshot$state_biomass,
    TotalBiomass = snapshot$state_total_biomass
  )

  species <- c(
    snapshot$species,
    "TotalBiomass"
  )

  totals <- as.numeric(
    crossprod(
      counts,
      biomass
    )
  )

  means <- totals / n_cells

  mins <- vapply(
    seq_len(ncol(biomass)),
    function(j) min(biomass[used, j]),
    numeric(1)
  )

  maxs <- vapply(
    seq_len(ncol(biomass)),
    function(j) max(biomass[used, j]),
    numeric(1)
  )

  data.frame(
    year = rep(as.integer(snapshot$year), length(species)),
    species = species,
    n_cells = rep(n_cells, length(species)),
    biomass_total = totals,
    biomass_mean = means,
    biomass_min = mins,
    biomass_max = maxs,
    stringsAsFactors = FALSE
  )
}

bs52_large_snapshot_to_long52 <- function(snapshot) {
  if (!is.list(snapshot) ||
      !identical(snapshot$format, "landisR52/state-snapshot-v1")) {
    bs52_large_stop("Snapshot compact invalide.")
  }

  biomass <- snapshot$state_biomass[
    snapshot$state_id,
    ,
    drop = FALSE
  ]

  total <- snapshot$state_total_biomass[
    snapshot$state_id
  ]

  max_age <- snapshot$state_max_age[
    snapshot$state_id
  ]

  wide <- cbind(
    biomass,
    TotalBiomass = total
  )

  species <- c(
    snapshot$species,
    "TotalBiomass"
  )

  n_species <- length(species)

  out <- data.frame(
    year = rep(
      as.numeric(snapshot$year),
      length(snapshot$cell_id) * n_species
    ),
    cell_id = rep(
      as.integer(snapshot$cell_id),
      each = n_species
    ),
    species = rep(
      species,
      times = length(snapshot$cell_id)
    ),
    biomass = as.numeric(t(wide)),
    max_age = rep(
      as.numeric(max_age),
      each = n_species
    ),
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL
  out
}


bs52_large_profile_row52 <- function(year, n_states_before, n_states_after,
                                     growth_seconds, feature_seconds,
                                     ward_seconds, rng_seconds,
                                     rebuild_seconds, snapshot_seconds,
                                     checkpoint_seconds,
                                     n_recruited_cells) {
  data.frame(
    year = as.integer(year),
    n_states_before = as.integer(n_states_before),
    n_states_after = as.integer(n_states_after),
    growth_seconds = as.numeric(growth_seconds),
    feature_seconds = as.numeric(feature_seconds),
    ward_seconds = as.numeric(ward_seconds),
    rng_seconds = as.numeric(rng_seconds),
    rebuild_seconds = as.numeric(rebuild_seconds),
    snapshot_seconds = as.numeric(snapshot_seconds),
    checkpoint_seconds = as.numeric(checkpoint_seconds),
    n_recruited_cells = as.integer(n_recruited_cells),
    stringsAsFactors = FALSE
  )
}

bs52_run_large_landscape52 <- function(
    ctx,
    duration,
    seed,
    engine = c("Rcpp", "R"),
    dispersal = c("ward", "local"),
    ward_variant = c("ward2019", "ward2022"),
    output_years = NULL,
    checkpoint_every = 50L,
    output_directory = NULL,
    run_signature = NULL,
    resume_checkpoint = NULL,
    collect_snapshots = FALSE,
    rng_chunk_cells = 100000L,
    snapshot_compression = c(
      "gzip1",
      "gzip",
      "none"
    ),
    checkpoint_compression = c(
      "gzip1",
      "gzip",
      "none"
    ),
    establishment = NULL,
    landis_rng = NULL,
    deferred = FALSE, pending_regeneration = NULL,
    report_progress = NULL, detailed_progress = FALSE) {

  engine <- match.arg(engine)
  dispersal <- match.arg(dispersal)
  ward_variant <- match.arg(ward_variant)

  snapshot_compression <- match.arg(
    snapshot_compression
  )

  checkpoint_compression <- match.arg(
    checkpoint_compression
  )

  duration <- bs52_large_integer(
    duration,
    "duration",
    1L
  )

  seed <- bs52_large_integer(
    seed,
    "seed",
    0L
  )

  timestep <- bs52_large_integer(
    ctx$config$timestep,
    "Timestep canonique",
    1L
  )

  if (is.null(output_years)) {
    output_years <- seq.int(
      0L,
      duration,
      by = timestep
    )
  }

  years <- as.integer(output_years)

  if (duration %% timestep != 0L ||
      !length(years) ||
      !0L %in% years ||
      any(years < 0L | years > duration) ||
      any(years %% timestep != 0L) ||
      anyDuplicated(years) ||
      is.unsorted(years)) {
    bs52_large_stop(
      "Duree et annees de sortie invalides."
    )
  }

  checkpoint_every <- bs52_large_integer(
    checkpoint_every,
    "checkpoint_every",
    timestep
  )

  if (checkpoint_every %% timestep != 0L) {
    bs52_large_stop(
      "checkpoint_every doit etre un multiple du timestep."
    )
  }

  if (identical(engine, "Rcpp")) {
    bs52_large_init_bulk52()
  }

  params_start <- proc.time()[["elapsed"]]
  params <- bs52_large_params52(
    ctx,
    duration
  )
  params_seconds <- bs52_large_elapsed(params_start)

  old_kind <- RNGkind()
  had_seed <- exists(
    ".Random.seed",
    .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit({
    do.call(
      RNGkind,
      as.list(old_kind)
    )

    if (had_seed) {
      assign(
        ".Random.seed",
        old_seed,
        .GlobalEnv
      )
    } else if (exists(
      ".Random.seed",
      .GlobalEnv,
      inherits = FALSE
    )) {
      rm(
        ".Random.seed",
        envir = .GlobalEnv
      )
    }
  }, add = TRUE)

  RNGkind(
    "Mersenne-Twister",
    "Inversion",
    "Rejection"
  )

  if (is.null(resume_checkpoint)) {
    set.seed(seed)

    states <- ctx$states
    state_eco <- ctx$state_eco
    cell_state <- ctx$cell_state
    start_year <- 0L
  } else {
    if (!is.list(resume_checkpoint) ||
        !identical(
          resume_checkpoint$format,
          "landisR52/state-checkpoint-v3-csr"
        ) ||
        !identical(
          resume_checkpoint$run_signature,
          run_signature
        )) {
      bs52_large_stop(
        "Checkpoint incompatible avec le run actuel."
      )
    }

    states <- resume_checkpoint$states
    state_eco <- resume_checkpoint$state_eco
    cell_state <- resume_checkpoint$cell_state
    start_year <- as.integer(resume_checkpoint$year)

    assign(
      ".Random.seed",
      resume_checkpoint$rng_state,
      .GlobalEnv
    )

    if (!is.null(output_directory) && is.null(report_progress)) message(
      "Reprise depuis le checkpoint de l'annee ",
      start_year,
      "."
    )
  }

  collected <- list()
  summary_rows <- list()
  profile_rows <- list()

  write_snapshot <- function(year) {
    snapshot_start <- proc.time()[["elapsed"]]

    snapshot <- .biomass_core$services$snapshot(
      states,
      state_eco,
      cell_state,
      ctx,
      year
    )

    summary_rows[[as.character(year)]] <<-
      bs52_large_snapshot_summary52(snapshot)

    if (isTRUE(collect_snapshots)) {
      collected[[as.character(year)]] <<- snapshot
    }

    if (!is.null(output_directory)) {
      path <- file.path(
        output_directory,
        "snapshots",
        sprintf("snapshot_%03d.rds", year)
      )

      .biomass_core$services$write_rds(
        snapshot,
        path,
        compress = snapshot_compression,
        overwrite = TRUE
      )
    }

    bs52_large_elapsed(snapshot_start)
  }

  if (start_year == 0L &&
      0L %in% years) {
    snapshot_seconds <- write_snapshot(0L)

    profile_rows[["0"]] <- bs52_large_profile_row52(
      0L,
      .biomass_core$services$n_states(states),
      .biomass_core$services$n_states(states),
      0,
      0,
      0,
      0,
      0,
      snapshot_seconds,
      0,
      0L
    )
  }

  next_years <- if (start_year < duration) {
    seq.int(
      start_year + timestep,
      duration,
      by = timestep
    )
  } else {
    integer()
  }

  report <- getOption("landisR.progress.report")
  if (!is.null(report_progress))
    report <- function(fraction, label) report_progress(end_year, fraction, label)
  for (end_year in next_years) {
    n_states_before <- .biomass_core$services$n_states(states)

    if (detailed_progress) {
      message("")
      message(
        "------------------------------------------------------------"
      )
      message(
        "Time step ",
        end_year - timestep,
        " -> ",
        end_year,
        " / ",
        format(
          n_states_before,
          big.mark = " "
        ),
        " etats en entree"
      )
    }

    if (!is.null(report)) report(0.01, "Preparing growth")
    growth_start <- proc.time()[["elapsed"]]

    block_years <- seq.int(
      end_year - timestep + 1L,
      end_year
    )

    growth_bundle <- bs52_large_grow_states52(
      states,
      state_eco,
      block_years,
      params,
      ctx,
      engine = engine,
      progress = detailed_progress,
      progress_label = paste0(
        "Annee ",
        end_year,
        " - croissance + features"
      ),
      feature_year = end_year, report = report
    )

    states <- growth_bundle$states

    growth_seconds <- bs52_large_elapsed(
      growth_start
    )

    recruited <- bs52_large_recruit52(
      states,
      state_eco,
      cell_state,
      ctx,
      end_year,
      variant = ward_variant,
      mode = dispersal,
      stochastic = TRUE,
      rng_chunk_cells = rng_chunk_cells,
      progress = detailed_progress,
      precomputed_features = growth_bundle$features,
      precomputed_ecology = growth_bundle$ecology,
      establishment = establishment,
      landis_rng = landis_rng,
      deferred = deferred,
      pending_regeneration = pending_regeneration, report = report
    )

    pending_regeneration <- NULL
    landis_rng <- recruited$landis_rng
    states <- recruited$states
    state_eco <- recruited$state_eco
    cell_state <- recruited$cell_state

    if (!is.null(report) &&
        (end_year %in% years ||
         (!is.null(output_directory) && end_year %% checkpoint_every == 0L)))
      report(0.94, "Saving forest state")
    snapshot_seconds <- 0

    if (end_year %in% years) {
      snapshot_seconds <- write_snapshot(
        end_year
      )
    }

    checkpoint_seconds <- 0

    if (!is.null(output_directory) &&
        end_year %% checkpoint_every == 0L) {
      checkpoint_start <- proc.time()[["elapsed"]]

      checkpoint <- .biomass_core$services$checkpoint(
        end_year,
        states,
        state_eco,
        cell_state,
        get(
          ".Random.seed",
          .GlobalEnv,
          inherits = FALSE
        ),
        run_signature
      )

      .biomass_core$services$write_rds(
        checkpoint,
        file.path(
          output_directory,
          "checkpoints",
          sprintf(
            "checkpoint_%03d.rds",
            end_year
          )
        ),
        compress = checkpoint_compression,
        overwrite = FALSE
      )

      checkpoint_seconds <- bs52_large_elapsed(
        checkpoint_start
      )
    }

    profile_rows[[as.character(end_year)]] <-
      bs52_large_profile_row52(
        end_year,
        n_states_before,
        .biomass_core$services$n_states(states),
        growth_seconds,
        recruited$timing$feature_seconds,
        recruited$timing$ward_seconds,
        recruited$timing$rng_seconds,
        recruited$timing$rebuild_seconds,
        snapshot_seconds,
        checkpoint_seconds,
        recruited$n_recruited_cells
      )

    if (detailed_progress) message(
      "Large52 / annee ",
      end_year,
      " / ",
      format(length(ctx$ids), big.mark = " "),
      " cellules actives / ",
      format(.biomass_core$services$n_states(states), big.mark = " "),
      " etats / growth=",
      round(growth_seconds, 2),
      " s / Ward=",
      round(recruited$timing$ward_seconds, 2),
      " s / rebuild=",
      round(recruited$timing$rebuild_seconds, 2),
      " s"
    )
  }

  profile <- do.call(
    rbind,
    profile_rows
  )

  attr(profile, "params_seconds") <- params_seconds

  list(
    states = states,
    state_eco = state_eco,
    cell_state = cell_state,
    snapshots = collected,
    summary = if (length(summary_rows)) {
      do.call(rbind, summary_rows)
    } else {
      data.frame()
    },
    profile = profile,
    landis_rng = landis_rng,
    final_rng_state = get(
      ".Random.seed",
      .GlobalEnv,
      inherits = FALSE
    )
  )
}

bs52_large_run_signature52 <- function(scenario, inputs_signature) {
  climate <- scenario$climate

  climate_md5 <- if (is.null(climate$dynamic_file)) {
    NA_character_
  } else {
    unname(
      tools::md5sum(
        climate$dynamic_file
      )
    )
  }

  list(
    format = "landisR52/large-run-v1",
    scenario = scenario[names(scenario) %in% c(
      "name",
      "duration",
      "timestep",
      "seed",
      "engine",
      "dispersal",
      "ward_variant",
      "output_years",
      "checkpoint_every",
      "snapshot_compression",
      "checkpoint_compression"
    )],
    climate = list(
      name = climate$name,
      dynamic_file = climate$dynamic_file,
      dynamic_md5 = climate_md5
    ),
    inputs_signature = inputs_signature,
    code_fingerprint = bs52_large_code_fingerprint52(),
    r_version = R.version.string,
    rcpp_version = as.character(
      utils::packageVersion("Rcpp")
    )
  )
}

# The interface migration preserves the production kernel and checkpoint payload.
# Accept known earlier engine signatures only if the rest of the manifest matches.
bs52_compatible_run_signature52 <- function(saved, current) {
  if (identical(saved, current)) return(TRUE)

  kernel <- "v3.11-ward-interior-fastpath"
  if (!is.list(saved) ||
      !identical(current$code_fingerprint[c("succession_version", "kernel")],
                 c(succession_version = "0.1.0.9001", kernel = kernel)))
    return(FALSE)
  supported <- list(c(package_version = "0.1.0.9008", kernel = kernel),
                    c(package_version = "0.1.0.9009", kernel = kernel),
                    c(succession_version = "0.1.0.9000", core_version = "0.1.0.9010",
                      kernel = kernel))
  if (!any(vapply(supported, identical, logical(1L), y = saved$code_fingerprint)))
    return(FALSE)

  saved$code_fingerprint <- current$code_fingerprint
  identical(saved, current)
}

bs52_large_find_latest_checkpoint52 <- function(output_directory,
                                                run_signature) {
  directory <- file.path(
    output_directory,
    "checkpoints"
  )

  if (!dir.exists(directory)) {
    return(NULL)
  }

  files <- list.files(
    directory,
    pattern = "^checkpoint_[0-9]+\\.rds$",
    full.names = TRUE
  )

  if (!length(files)) {
    return(NULL)
  }

  years <- suppressWarnings(
    as.integer(
      sub(
        "^checkpoint_([0-9]+)\\.rds$",
        "\\1",
        basename(files)
      )
    )
  )

  files <- files[order(years, decreasing = TRUE)]

  for (path in files) {
    checkpoint <- tryCatch(
      readRDS(path),
      error = function(e) NULL
    )

    if (is.list(checkpoint) &&
        identical(
          checkpoint$format,
          "landisR52/state-checkpoint-v3-csr"
        ) &&
        identical(
          checkpoint$run_signature,
          run_signature
        )) {
      return(checkpoint)
    }
  }

  NULL
}

bs52_run_large_production52 <- function(scenario) {
  required <- c(
    "name",
    "inputs",
    "duration",
    "timestep",
    "seed",
    "engine",
    "dispersal",
    "ward_variant",
    "climate",
    "track_dead_pools",
    "save_recruitment",
    "output_years",
    "output_directory",
    "checkpoint_every"
  )

  if (!is.list(scenario) ||
      !all(required %in% names(scenario))) {
    bs52_large_stop(
      "Scenario large52 incomplet."
    )
  }

  progress <- scenario[["progress"]]
  if (is.null(progress)) progress <- getOption("landisR.progress", interactive())
  detailed_progress <- isTRUE(progress) && bs52_large_progress_enabled52()
  bar <- .biomass_core$services$progress_bar(progress && !detailed_progress,
                                            "Biomass Succession")
  if (!is.null(bar)) {
    on.exit(bar$close(), add = TRUE)
    bar$update(0, "Preparing inputs")
  }
  read_inputs <- if (is.null(bar)) force else suppressMessages

  if (is.null(
        scenario$snapshot_compression
      )) {
    scenario$snapshot_compression <-
      "gzip1"
  }

  if (is.null(
        scenario$checkpoint_compression
      )) {
    scenario$checkpoint_compression <-
      "gzip1"
  }

  scenario$snapshot_compression <- match.arg(
    scenario$snapshot_compression,
    c(
      "gzip1",
      "gzip",
      "none"
    )
  )

  scenario$checkpoint_compression <- match.arg(
    scenario$checkpoint_compression,
    c(
      "gzip1",
      "gzip",
      "none"
    )
  )

  if (isTRUE(scenario$track_dead_pools)) {
    bs52_large_stop(
      "large52 v1 ne prend pas encore en charge track_dead_pools=TRUE."
    )
  }

  if (isTRUE(scenario$save_recruitment)) {
    bs52_large_stop(
      "large52 v1 n'ecrit pas les tables de recrutement massives."
    )
  }

  if (!identical(scenario$engine, "Rcpp") &&
      !identical(scenario$engine, "R")) {
    bs52_large_stop(
      "engine doit etre Rcpp ou R."
    )
  }

  scenario$inputs <- normalizePath(
    scenario$inputs,
    winslash = "/",
    mustWork = TRUE
  )

  scenario$output_directory <- normalizePath(
    scenario$output_directory,
    winslash = "/",
    mustWork = FALSE
  )

  if (!is.null(scenario$climate$dynamic_file)) {
    scenario$climate$dynamic_file <- normalizePath(
      scenario$climate$dynamic_file,
      winslash = "/",
      mustWork = TRUE
    )
  }

  cache_file <- if (!is.null(scenario$canonical_cache_file)) {
    scenario$canonical_cache_file
  } else {
    NULL
  }

  context_cache_file <- if (!is.null(
        scenario$context_cache_file
      )) {
    scenario$context_cache_file
  } else {
    bs52_large_derive_context_cache_file52(
      cache_file
    )
  }

  total_start <- proc.time()[["elapsed"]]

  load_start <- proc.time()[["elapsed"]]

  inputs_signature <- bs52_large_input_signature52(
    scenario$inputs
  )

  canonical <- read_inputs(bs52_large_read_canonical_cached52(
    scenario$inputs,
    cache_file
  ))

  load_seconds <- bs52_large_elapsed(load_start)

  climate_start <- proc.time()[["elapsed"]]

  canonical <- bs52_apply_climate52(
    canonical,
    scenario$climate
  )

  climate_seconds <- bs52_large_elapsed(
    climate_start
  )

  if (!identical(
        as.integer(scenario$timestep),
        canonical$timestep
      ) ||
      scenario$duration %% canonical$timestep != 0L ||
      scenario$duration > canonical$duration) {
    bs52_large_stop(
      "Le scenario doit respecter la duree et le timestep canoniques."
    )
  }

  if (identical(scenario$engine, "Rcpp")) {
    bs52_large_init_bulk52()
  }

  context_start <- proc.time()[["elapsed"]]

  ctx <- read_inputs(bs52_large_read_context_cached52(
    canonical = canonical,
    inputs = scenario$inputs,
    inputs_signature = inputs_signature,
    climate = scenario$climate,
    cache_file = context_cache_file
  ))

  context_seconds <- bs52_large_elapsed(
    context_start
  )

  # The full 5.6M-row canonical grid is no longer needed once the compact
  # context has been built.
  rm(canonical)
  invisible(gc())

  run_signature <- bs52_large_run_signature52(
    scenario,
    inputs_signature
  )

  output <- scenario$output_directory
  manifest_file <- file.path(
    output,
    "run_manifest.rds"
  )
  complete_file <- file.path(
    output,
    "complete.rds"
  )

  if (dir.exists(output) &&
      file.exists(manifest_file) &&
      file.exists(complete_file)) {
    saved_manifest <- readRDS(
      manifest_file
    )

    if (bs52_compatible_run_signature52(saved_manifest, run_signature)) {
      if (is.null(bar)) message(
        "Simulation large52 deja terminee : ",
        output
      )
      if (!is.null(bar)) bar$update(1, "Already complete")

      return(
        invisible(
          list(
            output = output,
            resumed = TRUE
          )
        )
      )
    }

    bs52_large_stop(
      "Le dossier contient une autre simulation : ",
      output
    )
  }

  if (!dir.exists(output)) {
    dir.create(
      output,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  if (!dir.exists(output)) {
    bs52_large_stop(
      "Impossible de creer le dossier : ",
      output
    )
  }

  if (file.exists(manifest_file)) {
    saved_manifest <- readRDS(
      manifest_file
    )

    if (!bs52_compatible_run_signature52(saved_manifest, run_signature)) {
      bs52_large_stop(
        "Manifest existant incompatible : ",
        output
      )
    }
    # Keep the original signature on appended checkpoints to preserve the
    # consistency of a run that was started before engine extraction.
    run_signature <- saved_manifest
  } else {
    .biomass_core$services$write_rds(
      run_signature,
      manifest_file,
      compress = "gzip",
      overwrite = FALSE
    )
  }

  resume_checkpoint <- bs52_large_find_latest_checkpoint52(
    output,
    run_signature
  )

  simulation_start <- proc.time()[["elapsed"]]

  report_progress <- NULL
  if (!is.null(bar)) {
    step_count <- scenario$duration / scenario$timestep
    total_phases <- 3 * step_count + 2
    report_progress <- function(end_year, fraction, label) {
      completed <- 1 + 3 * (end_year / scenario$timestep - 1 + fraction)
      bar$update(completed / total_phases, paste0("Year ", end_year, " | ", label))
    }
    if (is.null(resume_checkpoint)) {
      bar$update(1 / total_phases, "Preparing initial state")
    } else {
      completed <- 1 + 3 * resume_checkpoint$year / scenario$timestep
      bar$update(completed / total_phases, paste("Resuming from year", resume_checkpoint$year))
    }
  }
  result <- bs52_run_large_landscape52(
    ctx = ctx,
    duration = scenario$duration,
    seed = scenario$seed,
    engine = scenario$engine,
    dispersal = scenario$dispersal,
    ward_variant = scenario$ward_variant,
    output_years = scenario$output_years,
    checkpoint_every = scenario$checkpoint_every,
    output_directory = output,
    run_signature = run_signature,
    resume_checkpoint = resume_checkpoint,
    collect_snapshots = FALSE,
    rng_chunk_cells = if (!is.null(scenario$rng_chunk_cells)) {
      scenario$rng_chunk_cells
    } else {
      100000L
    },
    snapshot_compression =
      scenario$snapshot_compression,
    checkpoint_compression =
      scenario$checkpoint_compression,
    report_progress = report_progress,
    detailed_progress = detailed_progress
  )

  simulation_seconds <- bs52_large_elapsed(
    simulation_start
  )

  if (!is.null(bar)) bar$update((total_phases - 1) / total_phases, "Writing outputs")
  snapshot_files <- list.files(
    file.path(output, "snapshots"),
    pattern = "^snapshot_[0-9]+\\.rds$",
    full.names = TRUE
  )

  if (!length(snapshot_files)) {
    bs52_large_stop("Aucun snapshot compact n'a ete produit.")
  }

  snapshot_years <- suppressWarnings(
    as.integer(
      sub(
        "^snapshot_([0-9]+)\\.rds$",
        "\\1",
        basename(snapshot_files)
      )
    )
  )

  snapshot_files <- snapshot_files[
    order(snapshot_years)
  ]

  all_summary <- do.call(
    rbind,
    lapply(
      snapshot_files,
      function(path) {
        bs52_large_snapshot_summary52(
          readRDS(path)
        )
      }
    )
  )

  utils::write.csv(
    all_summary,
    file.path(
      output,
      "biomass_summary.csv"
    ),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )

  profile <- result$profile

  profile$load_canonical_seconds <- load_seconds
  profile$apply_climate_seconds <- climate_seconds
  profile$build_context_seconds <- context_seconds
  profile$params_seconds <- attr(
    result$profile,
    "params_seconds"
  )

  utils::write.csv(
    profile,
    file.path(
      output,
      "performance_profile.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  total_seconds <- bs52_large_elapsed(
    total_start
  )

  complete <- list(
    format = "landisR52/large-complete-v1",
    finished_at = Sys.time(),
    duration = scenario$duration,
    output_years = scenario$output_years,
    checkpoint_every = scenario$checkpoint_every,
    snapshot_compression =
      scenario$snapshot_compression,
    checkpoint_compression =
      scenario$checkpoint_compression,
    context_cache_file =
      context_cache_file,
    active_cells = length(ctx$ids),
    final_states = .biomass_core$services$n_states(result$states),
    load_canonical_seconds = load_seconds,
    apply_climate_seconds = climate_seconds,
    build_context_seconds = context_seconds,
    simulation_seconds = simulation_seconds,
    total_seconds = total_seconds
  )

  .biomass_core$services$write_rds(
    complete,
    complete_file,
    compress = "gzip",
    overwrite = FALSE
  )

  if (!is.null(bar)) {
    bar$update(1, paste("Complete at year", scenario$duration))
    bar$close()
  }

  message("")
  message("============================================================")
  message("LARGE52 TERMINE")
  message("============================================================")
  message(
    "Cellules actives : ",
    format(length(ctx$ids), big.mark = " ")
  )
  message(
    "Etats finaux     : ",
    format(
      .biomass_core$services$n_states(result$states),
      big.mark = " "
    )
  )
  message(
    "Lecture canonique: ",
    round(load_seconds / 60, 2),
    " min"
  )
  message(
    "Construction ctx : ",
    round(context_seconds / 60, 2),
    " min"
  )
  message(
    "Simulation       : ",
    round(simulation_seconds / 60, 2),
    " min"
  )
  message(
    "Total            : ",
    round(total_seconds / 60, 2),
    " min"
  )
  message(
    "Profil           : ",
    normalizePath(
      file.path(
        output,
        "performance_profile.csv"
      ),
      winslash = "/",
      mustWork = TRUE
    )
  )

  invisible(
    list(
      output = output,
      resumed = !is.null(resume_checkpoint),
      summary = all_summary,
      profile = profile,
      complete = complete
    )
  )
}
