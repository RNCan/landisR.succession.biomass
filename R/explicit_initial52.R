# Reader for Biomass Succession explicit-initial-biomass communities.
# Syntax supported: SPECIES age (biomass) [age (biomass) ...].
# It corresponds to InitialBiomass.CreateInitialCommunitiesCSV: no spinup.
bs52_read_explicit_initial_communities <- function(path, allowed_species = NULL) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- sub("^\\ufeff", "", x)
  if (!any(grepl('^LandisData\\s+"Initial Communities"', trimws(x))))
    stop("Not an explicit Initial Communities file")
  code <- NA_integer_; codes <- integer(); out_code <- list(); out_species <- list(); out_age <- list(); out_biomass <- list()
  for (line in x) {
    line <- trimws(line)
    if (!nzchar(line) || grepl('^LandisData', line)) next
    m <- regexec('^MapCode\\s+([0-9]+)\\s*$', line); hit <- regmatches(line, m)[[1]]
    if (length(hit)) { code <- as.integer(hit[2]); codes <- c(codes, code); next }
    if (is.na(code)) stop("Cohort before MapCode")
    species <- sub('\\s+.*$', '', line)
    if (!is.null(allowed_species) && !species %in% allowed_species)
      stop("Unknown species: ", species)
    pairs <- regmatches(line, gregexpr('([0-9]+)\\s*\\(\\s*([0-9]+)\\s*\\)', line, perl = TRUE))[[1]]
    if (!length(pairs)) stop("Expected age (biomass) pairs: ", line)
    z <- regmatches(pairs, regexec('([0-9]+)\\s*\\(\\s*([0-9]+)\\s*\\)', pairs, perl = TRUE))
    n <- length(z)
    out_code[[length(out_code) + 1L]] <- rep.int(code, n)
    out_species[[length(out_species) + 1L]] <- rep.int(species, n)
    out_age[[length(out_age) + 1L]] <- as.integer(vapply(z, `[[`, "", 2L))
    out_biomass[[length(out_biomass) + 1L]] <- as.integer(vapply(z, `[[`, "", 3L))
  }
  if (anyDuplicated(codes)) stop("Duplicate MapCode")
  ans <- if (length(out_code)) data.frame(map_code = unlist(out_code, use.names = FALSE),
    species = unlist(out_species, use.names = FALSE), age = unlist(out_age, use.names = FALSE),
    biomass = unlist(out_biomass, use.names = FALSE), stringsAsFactors = FALSE) else
    data.frame(map_code=integer(),species=character(),age=integer(),biomass=integer())
  if (any(ans$age < 1L) || any(ans$biomass < 0L) || anyDuplicated(ans[c('map_code','species','age')]))
    stop("Invalid or duplicate explicit cohort")
  list(community_codes = codes, cohorts = ans)
}

bs52_initialize_explicit <- function(initial, map_code) {
  if (!is.list(initial) || is.null(initial$community_codes) || is.null(initial$cohorts))
    stop("Invalid explicit-initial object")
  if (length(map_code) != 1L || is.na(map_code) || !map_code %in% initial$community_codes)
    stop("Unknown explicit initial-community map code")
  out <- initial$cohorts[initial$cohorts$map_code == map_code, c('species','age','biomass'), drop = FALSE]
  rownames(out) <- NULL
  out
}
