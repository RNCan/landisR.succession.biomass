# Reader for the supplied LANDIS text inputs, with age-only or explicit
# age (biomass) communities. Explicit biomass is retained in explicit_initial.
# Reject unsupported headers and missing inputs rather than guessing.
bs52_lines <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- sub("^\ufeff", "", x)
  x <- trimws(sub("<<.*$", "", x))
  x[nzchar(x) & !grepl("^>>", x)]
}
bs52_tokens <- function(x) scan(text = x, what = character(), quiet = TRUE,
                               quote = '"', comment.char = "")
bs52_value <- function(lines, key, required = TRUE) {
  rows <- lapply(lines, bs52_tokens)
  hits <- Filter(function(x) length(x) > 0 && x[1] == key, rows)
  if (!length(hits) && !required) return(NULL)
  if (length(hits) != 1L || length(hits[[1]]) != 2L)
    stop("Expected exactly one value for ", key)
  hits[[1]][2]
}
bs52_section <- function(lines, start, end = NULL) {
  i <- which(lines == start)
  j <- if (is.null(end)) length(lines) + 1L else which(vapply(lapply(lines, bs52_tokens), function(x) x[1] == end, logical(1)))
  if (length(i) != 1L || length(j) != 1L || j <= i)
    stop("Invalid section: ", start)
  if (j == i + 1L) return(character())
  lines[seq.int(i + 1L, j - 1L)]
}
bs52_table <- function(lines, columns) {
  if (!length(lines)) return(setNames(as.data.frame(replicate(length(columns),
    character(), simplify = FALSE)), columns))
  z <- lapply(lines, bs52_tokens)
  if (any(lengths(z) != length(columns))) stop("Unexpected table width")
  out <- as.data.frame(do.call(rbind, z), stringsAsFactors = FALSE)
  names(out) <- columns
  for (nm in names(out)) out[[nm]] <- type.convert(out[[nm]], as.is = TRUE)
  out
}
bs52_read_mauricie <- function(scenario_file, prefer_cells_csv) {
  root <- dirname(scenario_file)
  path <- function(x) {
    p <- chartr("\\", "/", x)
    if (!grepl("^(/|[A-Za-z]:/)", p)) p <- file.path(root, p)
    if (!file.exists(p)) stop("Missing Biomass Succession input: ", p)
    normalizePath(p, winslash = "/", mustWork = TRUE)
  }
  sc <- bs52_lines(scenario_file)
  main <- bs52_lines(path(bs52_value(sc, "Biomass Succession")))
  sp <- bs52_table(bs52_section(main, "SpeciesParameters", "EcoregionParameters"),
    c("species", "leaf_longevity", "woody_decay", "mortality_shape",
      "growth_shape", "leaf_lignin"))
  dy <- bs52_table(bs52_lines(path(bs52_value(main, "DynamicInputFile")))[-1],
    c("year", "ecoregion", "species", "sep", "max_anpp", "max_biomass"))
  if (anyDuplicated(sp$species) ||
      any(!dy$species %in% sp$species) ||
      anyDuplicated(dy[c("year", "ecoregion", "species")]))
    stop("Duplicate or unmatched species/dynamic keys")
  initial_file <- bs52_value(main, "InitialCommunities")
  initial <- bs52_lines(path(initial_file))
  if (!identical(bs52_tokens(initial[1]), c("LandisData", "Initial Communities")))
    stop("Unsupported initial-community format")
  code <- NULL; cohort_rows <- list(); community_codes <- numeric()
  for (line in initial[-1]) {
    t <- bs52_tokens(line)
    if (t[1] == "MapCode") {
      code <- as.numeric(t[2]); community_codes <- c(community_codes, code)
    } else {
      if (is.null(code) || !t[1] %in% sp$species || length(t) < 2)
        stop("Invalid initial-community row")
      age_tokens <- bs52_tokens(gsub("\\([^)]*\\)", "", line))[-1]
      ages <- suppressWarnings(as.numeric(age_tokens))
      if (any(!is.finite(ages)) || any(ages < 1) || any(ages != trunc(ages)))
        stop("Only explicit integer cohort ages are supported")
      cohort_rows[[length(cohort_rows) + 1L]] <- data.frame(
        map_code = code, species = t[1], age = ages)
    }
  }
  if (anyDuplicated(community_codes)) stop("Duplicate community code")
  cohorts <- if (length(cohort_rows)) do.call(rbind, cohort_rows) else
    data.frame(map_code = numeric(), species = character(), age = numeric())
  if (anyDuplicated(cohorts)) stop("Duplicate initial cohort")
  thresholds <- bs52_section(main, "MinRelativeBiomass", "SufficientLight")
  shade_ecos <- bs52_tokens(thresholds[1])
  shade <- bs52_table(gsub("%", "", thresholds[-1]), c("shade", shade_ecos))
  shade[-1] <- lapply(shade[-1], function(x) x / 100)
  light <- bs52_table(bs52_section(main, "SufficientLight", "SpeciesParameters"),
    c("shade_tolerance", paste0("shade_", 0:5)))
  fire <- bs52_table(
    bs52_section(main, "FireReductionParameters", "HarvestReductionParameters"),
    c("severity", "wood_reduction", "litter_reduction")
  )
  if (!nrow(fire) || anyNA(fire) || any(fire$severity != trunc(fire$severity)) ||
      any(fire$severity < 1L) || anyDuplicated(fire$severity) ||
      any(fire$wood_reduction < 0 | fire$wood_reduction > 1) ||
      any(fire$litter_reduction < 0 | fire$litter_reduction > 1))
    stop("Invalid FireReductionParameters")
  cells_path <- file.path(root, "cells.csv")
  community_by_cell <- NULL
  initial_map <- NULL
  if (isTRUE(prefer_cells_csv) && file.exists(cells_path)) {
    header <- utils::read.csv(cells_path, nrows = 0L)
    columns <- c("cell_id", "community")
    if (!all(columns %in% names(header)))
      stop("cells.csv must include cell_id and community for Biomass Succession.", call. = FALSE)
    community_by_cell <- utils::read.csv(
      cells_path, colClasses = ifelse(names(header) %in% columns, NA_character_, "NULL")
    )
  } else {
    initial_map <- bs52_import_landis_matrix(path(bs52_value(main, "InitialCommunitiesMap")))
  }

  list(timestep = as.numeric(bs52_value(main, "Timestep")),
       spinup_cohorts = isTRUE(tolower(bs52_value(main, "SpinupCohorts", required = FALSE)) %in% c("yes", "true")),
       seeding = bs52_value(main, "SeedingAlgorithm"),
       spinup_mortality = as.numeric(bs52_value(main, "SpinupMortalityFraction")),
       species_parameters = sp, ecoregion_parameters = bs52_table(
         bs52_section(main, "EcoregionParameters", "DynamicInputFile"),
         c("ecoregion", "aet")), dynamic = dy, initial_cohorts = cohorts,
       community_codes = community_codes, shade = shade, light = light,
       fire_reductions = fire,
       initial_map = initial_map, community_by_cell = community_by_cell,
       explicit_initial = if (any(grepl("\\(", initial)))
         bs52_read_explicit_initial_communities(path(initial_file), sp$species) else NULL)
}
