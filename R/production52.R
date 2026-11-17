# Production runner for canonical autonomous-R Biomass Succession inputs.
# It deliberately has no dependency on LANDIS scenario files or rasters.

bs52_production_stop <- function(...) stop(..., call. = FALSE)

bs52_production_path <- function(root, path, label, must_work = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    bs52_production_stop("Chemin invalide pour ", label, ".")
  candidate <- if (grepl("^(?:[A-Za-z]:)?[/\\\\]", path)) path else
    file.path(root, path)
  normalizePath(candidate, winslash = "/", mustWork = must_work)
}

bs52_production_integer <- function(x, label, minimum = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != trunc(x) || x < minimum)
    bs52_production_stop(label, " doit etre un entier >= ", minimum, ".")
  as.integer(x)
}


# Create and validate a production output directory before any file is written.
# This is intentionally stricter than a bare dir.create(): on mapped/network
# drives, directory creation can fail or become visible with a short delay.
bs52_production_prepare_output_directory52 <- function(path,
                                                       retries = 10L,
                                                       delay = 0.1) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    bs52_production_stop("Chemin de sortie invalide.")

  retries <- bs52_production_integer(retries, "retries", 1L)
  if (!is.numeric(delay) || length(delay) != 1L || is.na(delay) ||
      !is.finite(delay) || delay < 0)
    bs52_production_stop("delay doit etre un nombre >= 0.")

  if (!dir.exists(path)) {
    for (attempt in seq_len(retries)) {
      suppressWarnings(
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
      )
      if (dir.exists(path)) break
      if (attempt < retries && delay > 0) Sys.sleep(delay)
    }
  }

  if (!dir.exists(path))
    bs52_production_stop(
      "Impossible de creer le dossier de sortie apres ", retries,
      " tentative(s) : ", path
    )

  # A directory can exist but still be unusable on a mapped/network drive.
  probe <- file.path(
    path,
    paste0(
      ".landisR52_write_test_",
      Sys.getpid(), "_",
      format(Sys.time(), "%Y%m%d%H%M%OS6")
    )
  )

  writable <- tryCatch({
    writeLines("ok", probe, useBytes = TRUE)
    file.exists(probe)
  }, error = function(e) FALSE)

  if (file.exists(probe))
    unlink(probe)

  if (!isTRUE(writable))
    bs52_production_stop(
      "Le dossier de sortie existe mais n'est pas inscriptible : ", path
    )

  invisible(path)
}

bs52_read_production_scenario52 <- function(path = "scenario_production_R.R") {
  scenario_file <- normalizePath(path, winslash = "/", mustWork = TRUE)
  e <- new.env(parent = baseenv())
  sys.source(scenario_file, envir = e)
  if (!exists("scenario", envir = e, inherits = FALSE) || !is.list(e$scenario))
    bs52_production_stop("Le fichier scenario doit definir une liste nommee scenario.")
  s <- e$scenario
  root <- dirname(scenario_file)

  if (!identical(s$format, "landisR52/production-scenario-v1"))
    bs52_production_stop("Format de scenario de production non pris en charge.")
  if (!is.character(s$name) || length(s$name) != 1L || !nzchar(s$name) ||
      grepl("[\\\\/:*?\"<>|]", s$name))
    bs52_production_stop("Nom de scenario invalide.")
  s$inputs <- bs52_production_path(root, s$inputs, "inputs")
  s$duration <- bs52_production_integer(s$duration, "duration", 1L)
  s$timestep <- bs52_production_integer(s$timestep, "timestep", 1L)
  s$seed <- bs52_production_integer(s$seed, "seed", 0L)
  if (!identical(s$engine, "Rcpp") && !identical(s$engine, "R"))
    bs52_production_stop("engine doit etre Rcpp ou R.")
  if (is.null(s$dispersal)) s$dispersal <- "ward"
  if (!s$dispersal %in% c("ward", "local"))
    bs52_production_stop("dispersal doit etre ward ou local.")
  if (is.null(s$ward_variant)) s$ward_variant <- "ward2019"
  if (!s$ward_variant %in% c("ward2019", "ward2022"))
    bs52_production_stop("ward_variant doit etre ward2019 ou ward2022.")
  if (is.null(s$climate)) s$climate <- list(name = "baseline", dynamic_file = NULL)
  if (!is.list(s$climate))
    bs52_production_stop("climate doit etre une liste avec name et dynamic_file.")
  climate_name <- s$climate$name
  if (!is.character(climate_name) || length(climate_name) != 1L ||
      is.na(climate_name) || !nzchar(climate_name) ||
      grepl("[\\\\/:*?\"<>|]", climate_name))
    bs52_production_stop("Nom de climat invalide.")
  climate_file <- s$climate$dynamic_file
  if (!is.null(climate_file))
    climate_file <- bs52_production_path(root, climate_file,
                                         "climate$dynamic_file")
  s$climate <- list(name = climate_name, dynamic_file = climate_file)
  if (is.null(s$track_dead_pools)) s$track_dead_pools <- FALSE
  if (!isTRUE(s$track_dead_pools) && !identical(s$track_dead_pools, FALSE))
    bs52_production_stop("track_dead_pools doit etre TRUE ou FALSE.")
  if (is.null(s$dead_pool_initialization)) s$dead_pool_initialization <- "auto"
  if (!s$dead_pool_initialization %in% c("auto", "zero", "spinup"))
    bs52_production_stop("dead_pool_initialization doit etre auto, zero ou spinup.")
  if (is.null(s$save_recruitment)) s$save_recruitment <- FALSE
  if (!isTRUE(s$save_recruitment) && !identical(s$save_recruitment, FALSE))
    bs52_production_stop("save_recruitment doit etre TRUE ou FALSE.")
  if (is.null(s$resume)) s$resume <- TRUE
  if (!isTRUE(s$resume) && !identical(s$resume, FALSE))
    bs52_production_stop("resume doit etre TRUE ou FALSE.")

  years <- as.integer(s$output_years)
  if (!length(years) || anyNA(years) || any(years < 0L | years > s$duration) ||
      anyDuplicated(years) || is.unsorted(years) || !0L %in% years ||
      any(years %% s$timestep != 0L))
    bs52_production_stop("output_years doit inclure 0 et seulement des pas de temps valides.")
  s$output_years <- years
  s$output_directory <- bs52_production_path(root, s$output_directory,
                                              "output_directory", must_work = FALSE)
  s$root <- root
  s$file <- scenario_file
  s
}

bs52_production_internalize_ecoregions <- function(config, cells, active) {
  labels <- unique(as.character(cells$ecoregion[active]))
  if (!length(labels) || anyNA(labels))
    bs52_production_stop("Aucune ecoregion active valide dans les entrees canoniques.")
  code <- setNames(seq_along(labels), labels)
  if (!all(labels %in% as.character(config$dynamic$ecoregion)) ||
      !all(labels %in% as.character(config$ecoregion_parameters$ecoregion)) ||
      !all(labels %in% names(config$shade)))
    bs52_production_stop("Parametres ecologiques incomplets pour une ecoregion active.")

  keep_dynamic <- as.character(config$dynamic$ecoregion) %in% labels
  config$dynamic <- config$dynamic[keep_dynamic, , drop = FALSE]
  config$dynamic$ecoregion <- unname(code[as.character(config$dynamic$ecoregion)])
  keep_parameters <- as.character(config$ecoregion_parameters$ecoregion) %in% labels
  config$ecoregion_parameters <- config$ecoregion_parameters[keep_parameters, , drop = FALSE]
  config$ecoregion_parameters$ecoregion <-
    unname(code[as.character(config$ecoregion_parameters$ecoregion)])
  shade_names <- names(config$shade)
  shade_names[match(labels, shade_names)] <- as.character(unname(code[labels]))
  names(config$shade) <- shade_names
  list(config = config, code = code)
}

bs52_production_context52 <- function(canonical) {
  if (!is.list(canonical) || is.null(canonical$landscape) ||
      is.null(canonical$landscape$cells))
    bs52_production_stop("Entrees canoniques invalides.")
  cells <- canonical$landscape$cells
  active <- as.logical(canonical$landscape$active)
  if (length(active) != nrow(cells) || anyNA(active))
    bs52_production_stop("Etat actif canonique invalide.")
  internal <- bs52_production_internalize_ecoregions(canonical, cells, active)
  config <- internal$config
  codes <- integer(nrow(cells))
  codes[active] <- unname(internal$code[as.character(cells$ecoregion[active])])
  if (anyNA(codes[active])) bs52_production_stop("Codage interne des ecoregions invalide.")
  ecoregion <- matrix(codes, nrow = canonical$landscape$nrow,
                      ncol = canonical$landscape$ncol)
  ids <- as.integer(cells$cell_id[active])
  empty <- data.frame(species = character(), age = numeric(), biomass = numeric(),
                      stringsAsFactors = FALSE)
  initial <- lapply(seq_len(nrow(cells)), function(id) {
    if (!active[id]) return(empty)
    if (!is.null(config$explicit_initial))
      bs52_initialize_explicit(config$explicit_initial, cells$community[id])
    else bs52_spinup_live(config, cells$community[id], ecoregion[id])
  })
  list(config = config, cells = cells, eco = ecoregion, ids = ids,
       initial = initial, active = active, ecoregion_code = internal$code)
}

bs52_production_params52 <- function(ctx, duration) {
  ecoregions <- as.character(sort(unique(ctx$eco[ctx$ids])))
  lapply(seq_len(duration), function(year) {
    values <- lapply(ecoregions, function(eco)
      bs52_dynamic_at(ctx$config, year, as.integer(eco)))
    setNames(values, ecoregions)
  })
}

bs52_production_state_groups52 <- function(cohorts, ctx) {
  keys <- vapply(ctx$ids, function(id) {
    x <- cohorts[[id]]
    paste(ctx$eco[id], paste(x$species, x$age, x$biomass, collapse = ";"),
          sep = "|")
  }, "")
  unname(split(ctx$ids, match(keys, unique(keys))))
}

bs52_production_grow_block52 <- function(cohorts, groups, years, params, ctx,
                                          engine = c("Rcpp", "R"), dead = NULL,
                                          timestep = ctx$config$timestep) {
  engine <- match.arg(engine)
  for (group in groups) {
    id <- group[1L]
    eco <- as.character(ctx$eco[id])
    block_params <- lapply(years, function(year) params[[year]][[eco]])
    if (is.null(dead)) {
      value <- if (identical(engine, "Rcpp"))
        bs52_grow_live_block_rcpp(cohorts[[id]], block_params, timestep,
                                  combine_last = tail(years, 1L) %% timestep == 0L,
                                  spinup = FALSE) else
        bs52_grow_live_block_fast(cohorts[[id]], block_params, timestep,
                                  combine_last = tail(years, 1L) %% timestep == 0L,
                                  spinup = FALSE)
      cohorts[group] <- rep(list(value), length(group))
      next
    }
    aet <- ctx$config$ecoregion_parameters$aet[
      match(ctx$eco[id], ctx$config$ecoregion_parameters$ecoregion)]
    z <- bs52_grow_live_block_dead_pools_rcpp(
      cohorts[[id]], block_params, ctx$config$species_parameters, aet, years, timestep,
      dead$wood[group], dead$wood_decay[group], dead$litter[group],
      dead$litter_decay[group],
      combine_last = tail(years, 1L) %% timestep == 0L, spinup = FALSE
    )
    cohorts[group] <- rep(list(z$cohorts), length(group))
    dead$wood[group] <- z$wood
    dead$wood_decay[group] <- z$wood_decay
    dead$litter[group] <- z$litter
    dead$litter_decay[group] <- z$litter_decay
    dead$wood_input[group] <- 0
    dead$litter_input[group] <- 0
  }
  if (is.null(dead)) cohorts else list(cohorts = cohorts, dead = dead)
}

bs52_production_snapshot52 <- function(cohorts, year, ctx,
                                        engine = c("Rcpp", "R")) {
  engine <- match.arg(engine)
  if (identical(engine, "Rcpp"))
    return(bs52_record_landscape_rcpp(cohorts, year, ctx))
  species <- ctx$config$species$species
  do.call(rbind, lapply(ctx$ids, function(id) {
    x <- cohorts[[id]]
    biomass <- setNames(numeric(length(species)), species)
    if (nrow(x)) {
      total <- tapply(x$biomass, x$species, sum)
      biomass[names(total)] <- total
    }
    data.frame(year = as.numeric(year), cell_id = id,
               species = c(names(biomass), "TotalBiomass"),
               biomass = as.numeric(c(biomass, sum(biomass))),
               max_age = if (nrow(x)) max(x$age) else 0)
  }))
}

bs52_production_dead_initial52 <- function(ctx, initialization) {
  mode <- initialization
  if (identical(mode, "auto"))
    mode <- if (is.null(ctx$config$explicit_initial)) "spinup" else "zero"
  if (identical(mode, "zero")) return(bs52_dead_landscape_zero(ctx$eco))
  bs52_dead_landscape_spinup(ctx)
}

bs52_run_landscape52 <- function(ctx, duration, seed, engine = c("Rcpp", "R"),
                                 dispersal = c("ward", "local"),
                                 ward_variant = c("ward2019", "ward2022"),
                                 output_years = NULL,
                                 track_dead_pools = FALSE,
                                 dead_pool_initialization = "auto") {
  engine <- match.arg(engine)
  dispersal <- match.arg(dispersal)
  ward_variant <- match.arg(ward_variant)
  duration <- bs52_production_integer(duration, "duration", 1L)
  seed <- bs52_production_integer(seed, "seed", 0L)
  timestep <- bs52_production_integer(ctx$config$timestep, "Timestep canonique", 1L)
  if (is.null(output_years)) output_years <- seq.int(0L, duration, by = timestep)
  years <- as.integer(output_years)
  if (duration %% timestep != 0L || !length(years) || !0L %in% years ||
      any(years < 0L | years > duration) || any(years %% timestep != 0L) ||
      anyDuplicated(years) || is.unsorted(years))
    bs52_production_stop("Duree et annees de sortie doivent etre des multiples du Timestep, avec l'annee 0.")
  if (isTRUE(track_dead_pools) && !identical(engine, "Rcpp"))
    bs52_production_stop("Les pools morts de production requierent le moteur Rcpp.")

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)

  params <- bs52_production_params52(ctx, duration)
  cohorts <- ctx$initial
  snapshots <- list()
  recruitment <- list()
  dead <- if (isTRUE(track_dead_pools))
    bs52_production_dead_initial52(ctx, dead_pool_initialization) else NULL
  dead_snapshots <- list()
  if (0L %in% years) {
    snapshots[["0"]] <- bs52_production_snapshot52(cohorts, 0L, ctx, engine)
    if (!is.null(dead))
      dead_snapshots[["0"]] <- bs52_dead_landscape_snapshot(dead, ctx$ids, 0L)
  }

  for (end_year in seq.int(timestep, duration, by = timestep)) {
    block_years <- seq.int(end_year - timestep + 1L, end_year)
    groups <- bs52_production_state_groups52(cohorts, ctx)
    grown <- bs52_production_grow_block52(cohorts, groups, block_years, params,
                                           ctx, engine, dead, timestep)
    if (is.null(dead)) cohorts <- grown else {
      cohorts <- grown$cohorts
      dead <- grown$dead
    }
    recruited <- bs52_ward_recruit_landscape_fast(
      cohorts, ctx$eco, ctx$config, end_year, variant = ward_variant,
      mode = dispersal, use_rcpp = identical(engine, "Rcpp")
    )
    cohorts <- recruited$cohorts
    recruitment[[as.character(end_year)]] <- recruited$establishment
    if (end_year %in% years) {
      snapshots[[as.character(end_year)]] <-
        bs52_production_snapshot52(cohorts, end_year, ctx, engine)
      if (!is.null(dead))
        dead_snapshots[[as.character(end_year)]] <-
          bs52_dead_landscape_snapshot(dead, ctx$ids, end_year)
    }
    message("Production / annee ", end_year, " / ", length(ctx$ids),
            " cellules actives")
  }
  list(snapshots = snapshots, recruitment = recruitment, final_cohorts = cohorts,
       dead_snapshots = dead_snapshots, final_dead = dead)
}

bs52_production_biomass_summary52 <- function(snapshots) {
  if (!length(snapshots)) return(data.frame())
  out <- lapply(snapshots, function(x) {
    species <- unique(as.character(x$species))
    n_cells <- length(unique(x$cell_id))
    expected <- n_cells * length(species)
    if (nrow(x) != expected)
      bs52_production_stop("Instantane de biomasse incomplet.")
    biomass <- matrix(x$biomass, nrow = n_cells, ncol = length(species), byrow = TRUE)
    data.frame(year = as.integer(x$year[1L]), species = species, n_cells = n_cells,
               biomass_total = colSums(biomass), biomass_mean = colMeans(biomass),
               biomass_min = apply(biomass, 2L, min),
               biomass_max = apply(biomass, 2L, max), row.names = NULL)
  })
  do.call(rbind, out)
}

bs52_production_fingerprint52 <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[!file.info(files)$isdir]
  files <- sort(normalizePath(files, winslash = "/", mustWork = TRUE))
  data.frame(path = substring(files, nchar(normalizePath(root, winslash = "/")) + 2L),
             md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
}

bs52_production_manifest52 <- function(scenario, canonical) {
  code_files <- c("R/biomass52.R", "R/dead52.R", "R/explicit_initial52.R",
                  "R/fast52.R", "R/rcpp52.R", "R/simulation52.R",
                  "R/spinup52.R", "R/ward52.R", "R/canonical52.R",
                  "R/climate52.R", "R/production52.R", "src/growth52.cpp")
  code_files <- code_files[file.exists(code_files)]
  climate <- scenario$climate
  climate$dynamic_md5 <- if (is.null(climate$dynamic_file)) NA_character_ else
    unname(tools::md5sum(climate$dynamic_file))
  list(
    format = "landisR52/production-run-v1",
    scenario = scenario[names(scenario) %in% c(
      "format", "name", "inputs", "duration", "timestep", "seed", "engine",
      "dispersal", "ward_variant", "track_dead_pools", "dead_pool_initialization",
      "save_recruitment", "output_years", "output_directory")],
    climate = climate,
    canonical_input_fingerprint = bs52_production_fingerprint52(canonical$root),
    code_fingerprint = unname(tools::md5sum(code_files)),
    r_version = R.version.string,
    rcpp_version = if (identical(scenario$engine, "Rcpp"))
      as.character(utils::packageVersion("Rcpp")) else NA_character_
  )
}

bs52_run_production52 <- function(scenario) {
  if (!is.list(scenario) || !identical(scenario$format,
                                       "landisR52/production-scenario-v1"))
    bs52_production_stop("Scenario de production invalide.")
  if (!exists("bs52_apply_climate52", mode = "function"))
    bs52_production_stop("Charger R/climate52.R avant le lanceur de production.")
  if (is.null(scenario$climate))
    scenario$climate <- list(name = "baseline", dynamic_file = NULL)
  canonical <- bs52_read_canonical52(scenario$inputs)
  canonical <- bs52_apply_climate52(canonical, scenario$climate)
  if (!identical(scenario$timestep, canonical$timestep) ||
      scenario$duration %% canonical$timestep != 0L ||
      scenario$duration > canonical$duration)
    bs52_production_stop("Le scenario doit respecter la duree et le Timestep des entrees canoniques.")
  if (identical(scenario$engine, "Rcpp")) bs52_init_rcpp()
  ctx <- bs52_production_context52(canonical)
  manifest <- bs52_production_manifest52(scenario, canonical)
  output <- scenario$output_directory
  manifest_file <- file.path(output, "run_manifest.rds")
  complete_file <- file.path(output, "complete.rds")
  if (dir.exists(output) && file.exists(manifest_file) && file.exists(complete_file)) {
    if (isTRUE(scenario$resume) && identical(readRDS(manifest_file), manifest)) {
      message("Simulation deja terminee : ", output)
      return(invisible(list(output = output, resumed = TRUE)))
    }
    bs52_production_stop("Le dossier de sortie contient deja une simulation differente : ", output)
  }
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)))
    bs52_production_stop("Le dossier de sortie est incomplet ou non vide : ", output,
                         ". Choisir un nouveau dossier de sortie.")

  # The output directory must exist and be writable before run_manifest.rds is
  # created. This avoids failures observed on mapped/network drives where a
  # bare dir.create() did not guarantee that the path was immediately usable.
  bs52_production_prepare_output_directory52(output)

  saveRDS(manifest, manifest_file, compress = "gzip")

  result <- bs52_run_landscape52(
    ctx, scenario$duration, scenario$seed, engine = scenario$engine,
    dispersal = scenario$dispersal, ward_variant = scenario$ward_variant,
    output_years = scenario$output_years,
    track_dead_pools = scenario$track_dead_pools,
    dead_pool_initialization = scenario$dead_pool_initialization
  )
  summary <- bs52_production_biomass_summary52(result$snapshots)
  saveRDS(result$snapshots, file.path(output, "biomass_snapshots.rds"), compress = "gzip")
  saveRDS(result$final_cohorts, file.path(output, "final_cohorts.rds"),
          compress = "gzip")
  if (isTRUE(scenario$save_recruitment))
    saveRDS(result$recruitment, file.path(output, "recruitment.rds"), compress = "gzip")
  if (isTRUE(scenario$track_dead_pools))
    saveRDS(result$dead_snapshots, file.path(output, "dead_pool_snapshots.rds"),
            compress = "gzip")
  utils::write.csv(summary, file.path(output, "biomass_summary.csv"),
                   row.names = FALSE, na = "", fileEncoding = "UTF-8")
  saveRDS(list(finished_at = Sys.time(), duration = scenario$duration,
               output_years = scenario$output_years), complete_file,
          compress = "gzip")
  invisible(list(output = output, resumed = FALSE, summary = summary))
}

bs52_run_production_scenario52 <- function(path = "scenario_production_R.R") {
  bs52_run_production52(bs52_read_production_scenario52(path))
}
