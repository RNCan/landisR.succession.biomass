landisR_serotiny_prefire_sources <- function(simulation,cells) {
  traits <- landisR_normalize_serotiny_traits(simulation$ctx$config$species)
  sero <- which(traits$serotiny_type != 'none')
  labels <- as.character(traits$species[sero])
  if (!length(sero) || !length(cells))
    return(stats::setNames(vector('list',length(sero)),labels))
  positions <- match(cells,as.integer(simulation$ctx$ids))
  if (anyNA(positions))
    return(.landisR_legacy_prefire_270F(simulation,cells))
  ans <- basefire_prefire_csr_252(simulation$states,
    as.integer(simulation$cell_state),as.integer(positions),
    as.integer(cells),as.integer(sero),as.numeric(traits$maturity[sero]),0L)
  names(ans) <- labels
  ans
}

.landisR_baseprob_270F <- function(simulation,config,year) {
  eco <- as.integer(simulation$state_eco)
  codes <- sort(unique(eco))
  sp <- as.character(config$species$species)
  nsp <- length(sp)
  thresholds <- vector('list',length(codes))
  probability <- vector('list',length(codes))
  maxbio <- numeric(length(codes))
  for (k in seq_along(codes)) {
    key <- as.character(codes[[k]])
    p <- bs52_dynamic_at(config,year,codes[[k]])
    if (!identical(as.character(p$species),sp))
      stop('270F: ordre des especes incompatible avec la reference 259')
    thresholds[[k]] <- as.numeric(config$shade[[key]][order(config$shade$shade)])
    if (!length(thresholds[[k]]) || anyNA(thresholds[[k]]))
      stop('270F: seuils de classe d ombre incomplets')
    maxbio[[k]] <- max(p$max_biomass)
    light_rows <- match(p$shade_tolerance,config$light$shade_tolerance)
    if (anyNA(light_rows)) stop('270F: shade_tolerance introuvable')
    class_prob <- matrix(0,nrow=length(thresholds[[k]])+1L,ncol=nsp)
    for (shade in seq.int(0L,length(thresholds[[k]]))) {
      col <- paste0('shade_',shade)
      if (!col %in% names(config$light)) stop('270F: lumiere absente: ',col)
      light <- as.numeric(config$light[light_rows,col])
      class_prob[shade+1L,] <- pmin(1,pmax(0,p$sep*light))
    }
    probability[[k]] <- class_prob
  }
  basefire_serotiny_baseprob_csr_259(simulation$states,eco,
    as.integer(codes),thresholds,as.numeric(maxbio),probability)
}

.landisR_native_site_build_270F <- function(simulation,pairs,active_sero,
    species_names,arrival,probability,get_params,timestep) {
  ecos <- as.integer(simulation$ctx$active_eco)
  eco_keys <- sort(unique(ecos[unique(pairs[,1L])]))
  anpp <- matrix(NA_real_,nrow=length(eco_keys),ncol=length(species_names))
  eco_max <- numeric(length(eco_keys))
  for (i in seq_along(eco_keys)) {
    p <- get_params(eco_keys[[i]])
    eco_max[[i]] <- max(p$max_biomass)
    anpp[i,] <- as.numeric(p$max_anpp[match(species_names,p$species)])
  }
  basefire_serotiny_build_256(simulation$states,as.integer(simulation$cell_state),
    as.integer(simulation$ctx$ids),ecos,pairs,as.integer(active_sero),
    as.character(species_names),arrival,probability,as.integer(eco_keys),
    anpp,eco_max,as.integer(timestep))
}

# landisR_serotiny_recruit_after_fire est genere AU BUILD, depuis
# la fonction source 0.1.0.9003 et deux ancres AST strictes.

landisR_serotiny_recruit_after_fire <- function (simulation, prefire_sources, stochastic = TRUE, year = simulation$year)
{
    config <- simulation$ctx$config
    species <- landisR_normalize_serotiny_traits(config$species)
    sero <- which(species$serotiny_type != "none")
    source_count <- integer(length(sero))
    names(source_count) <- species$species[sero]
    if (!length(sero)) {
        return(list(simulation = simulation, source_count = source_count, recruited_cells = 0L, recruited_cohorts = 0L, recruitment = data.frame()))
    }
    source_count[] <- vapply(species$species[sero], function(sp) length(prefire_sources[[sp]]), integer(1L))
    active_sero <- sero[source_count > 0L]
    if (!length(active_sero)) {
        return(list(simulation = simulation, source_count = source_count, recruited_cells = 0L, recruited_cohorts = 0L, recruitment = data.frame()))
    }
    states <- NULL
    features <- list(base_probability = .landisR_baseprob_270F(simulation, config, year))
    n_sites <- length(simulation$ctx$ids)
    n_sero <- length(active_sero)
    arrival <- matrix(0, nrow = n_sites, ncol = n_sero, dimnames = list(NULL, species$species[active_sero]))
    for (j in seq_along(active_sero)) {
        sp_index <- active_sero[j]
        sp <- species$species[sp_index]
        arrival[, j] <- landisR_serotiny_arrival(simulation, prefire_sources[[sp]], species[sp_index, , drop = FALSE])
    }
    base <- features$base_probability[simulation$cell_state, active_sero, drop = FALSE]
    if (!is.matrix(base) || !identical(dim(base), dim(arrival))) {
        stop("Dimensions incompatibles entre probabilite d'etablissement et arrivee serotineuse.", call. = FALSE)
    }
    probability <- base * arrival
    probability[!is.finite(probability)] <- 0
    probability[probability < 0] <- 0
    probability[probability > 1] <- 1
    dim(probability) <- dim(arrival)
    dimnames(probability) <- dimnames(arrival)
    if (isTRUE(stochastic)) {
        rng <- landisR_serotiny_rng_draw(simulation, n_sites, n_sero)
        simulation <- rng$simulation
        recruited <- rng$draw < probability
    }
    else {
        recruited <- probability > 0
    }
    excluded <- simulation$establishment$excluded_cells
    if (length(excluded))
        recruited[match(excluded, simulation$ctx$ids), ] <- FALSE
    pairs <- which(recruited, arr.ind = TRUE)
    if (!nrow(pairs)) {
        return(list(simulation = simulation, source_count = source_count, recruited_cells = 0L, recruited_cohorts = 0L, recruitment = data.frame()))
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
            assign(key, bs52_dynamic_at(config, year, eco), envir = params_cache)
        }
        get(key, params_cache, inherits = FALSE)
    }
    group_names <- names(cell_groups)
    {
        built_270F <- .landisR_native_site_build_270F(simulation, pairs, active_sero, as.character(species$species), arrival, probability, get_params, config$timestep)
        replacement <- list(built_270F$replacement)
        record <- list(built_270F$recruitment)
        record_i <- 1L
    }
    replacement <- do.call(rbind, replacement)
    rownames(replacement) <- NULL
    simulation <- .biomass_core$services$apply_cohorts(simulation, cohorts = replacement, clear_cells = integer())
    recruitment <- do.call(rbind, record[seq_len(record_i)])
    rownames(recruitment) <- NULL
    list(simulation = simulation, source_count = source_count, recruited_cells = length(unique(recruitment$cell_id)), recruited_cohorts = nrow(recruitment), recruitment = recruitment)
}
