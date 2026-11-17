# Climate-specific dynamic inputs for autonomous-R Biomass Succession.

bs52_climate_stop <- function(...) stop(..., call. = FALSE)

bs52_climate_columns52 <- c("year", "ecoregion", "species", "sep",
                            "max_anpp", "max_biomass")

bs52_climate_normalize_dynamic52 <- function(dynamic, label = "Table dynamique",
                                              allow_sep_above_one = FALSE) {
  if (!is.data.frame(dynamic) || !identical(names(dynamic), bs52_climate_columns52))
    bs52_climate_stop(label, " doit contenir exactement : ",
                      paste(bs52_climate_columns52, collapse = ", "))
  out <- data.frame(
    year = suppressWarnings(as.numeric(dynamic$year)),
    ecoregion = as.character(dynamic$ecoregion),
    species = as.character(dynamic$species),
    sep = suppressWarnings(as.numeric(dynamic$sep)),
    max_anpp = suppressWarnings(as.numeric(dynamic$max_anpp)),
    max_biomass = suppressWarnings(as.numeric(dynamic$max_biomass)),
    stringsAsFactors = FALSE
  )
  if (!nrow(out) || anyNA(out) || any(!is.finite(out$year)) ||
      any(out$year != trunc(out$year)) || any(out$year < 0L) ||
      any(!nzchar(out$ecoregion)) || any(!nzchar(out$species)) ||
      any(!is.finite(out$sep)) || any(!is.finite(out$max_anpp)) ||
      any(!is.finite(out$max_biomass)) || any(out$sep < 0) ||
      (!isTRUE(allow_sep_above_one) && any(out$sep > 1)) ||
      any(out$max_anpp < 0 | out$max_biomass < 0) ||
      anyDuplicated(out[c("year", "ecoregion", "species")]))
    bs52_climate_stop(label, " contient des valeurs dynamiques invalides.")
  out$year <- as.integer(out$year)
  out
}

bs52_read_landis_dynamic52 <- function(path, allow_sep_above_one = FALSE) {
  source_file <- normalizePath(path, winslash = "/", mustWork = TRUE)
  lines <- readLines(source_file, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^\\ufeff", "", lines)
  lines <- trimws(sub("<<.*$", "", lines))
  lines <- lines[nzchar(lines)]
  if (!length(lines) || !grepl('^LandisData +"?Dynamic Input', lines[1L]))
    bs52_climate_stop("Fichier dynamique LANDIS non reconnu : ", source_file)
  body <- lines[-1L]
  body <- body[!grepl("^>>", body)]
  tokens <- strsplit(body, "\\s+")
  if (!length(tokens) || any(lengths(tokens) != 6L))
    bs52_climate_stop("Largeur invalide dans le fichier dynamique LANDIS : ", source_file)
  dynamic <- data.frame(
    year = vapply(tokens, `[[`, "", 1L),
    ecoregion = vapply(tokens, `[[`, "", 2L),
    species = vapply(tokens, `[[`, "", 3L),
    sep = vapply(tokens, `[[`, "", 4L),
    max_anpp = vapply(tokens, `[[`, "", 5L),
    max_biomass = vapply(tokens, `[[`, "", 6L),
    stringsAsFactors = FALSE
  )
  bs52_climate_normalize_dynamic52(
    dynamic, "Fichier dynamique LANDIS", allow_sep_above_one
  )
}

bs52_write_landis_dynamic52 <- function(dynamic, path) {
  dynamic <- bs52_climate_normalize_dynamic52(dynamic,
                                               "Fichier dynamique LANDIS corrige")
  format_number <- function(x) format(x, scientific = FALSE, trim = TRUE,
                                      digits = 15L)
  rows <- vapply(seq_len(nrow(dynamic)), function(i) paste(
    dynamic$year[i], dynamic$ecoregion[i], dynamic$species[i],
    format_number(dynamic$sep[i]), format_number(dynamic$max_anpp[i]),
    format_number(dynamic$max_biomass[i]), sep = "\t"
  ), "")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    'LandisData "Dynamic Input Data"', "",
    ">>\tyear\tlandtype\tspecies\tprobEst\tmaxANPP\tmaxB", rows
  ), path, useBytes = TRUE)
  invisible(path)
}

bs52_cap_landis_sep52 <- function(source_file, destination) {
  source_file <- normalizePath(source_file, winslash = "/", mustWork = TRUE)
  destination <- normalizePath(destination, winslash = "/", mustWork = FALSE)
  if (file.exists(destination))
    bs52_climate_stop("Le fichier climatique corrige existe deja : ", destination)
  dynamic <- bs52_read_landis_dynamic52(source_file, allow_sep_above_one = TRUE)
  corrected <- dynamic$sep > 1
  original_sep <- dynamic$sep[corrected]
  dynamic$sep[corrected] <- 1
  bs52_write_landis_dynamic52(dynamic, destination)
  list(
    source_file = source_file,
    destination = destination,
    n_corrected = sum(corrected),
    corrections = data.frame(
      year = dynamic$year[corrected], ecoregion = dynamic$ecoregion[corrected],
      species = dynamic$species[corrected], sep_original = original_sep,
      sep_corrected = dynamic$sep[corrected], stringsAsFactors = FALSE
    )
  )
}

bs52_validate_climate_dynamic52 <- function(dynamic, canonical) {
  dynamic <- bs52_climate_normalize_dynamic52(dynamic, "Table dynamique climatique")
  if (!is.list(canonical) || is.null(canonical$species) ||
      is.null(canonical$ecoregions) || is.null(canonical$landscape))
    bs52_climate_stop("Entrees canoniques invalides pour la validation climatique.")
  valid_ecoregions <- as.character(canonical$ecoregions$ecoregion)
  if (any(!dynamic$ecoregion %in% valid_ecoregions) ||
      any(!dynamic$species %in% canonical$species$species))
    bs52_climate_stop("La table dynamique climatique contient une ecoregion ou une espece inconnue.")
  active <- as.logical(canonical$landscape$active)
  active_ecoregions <- unique(as.character(
    canonical$landscape$cells$ecoregion[active]
  ))
  expected <- expand.grid(ecoregion = active_ecoregions,
                          species = as.character(canonical$species$species),
                          stringsAsFactors = FALSE)
  key <- paste(dynamic$ecoregion, dynamic$species, sep = "\r")
  first_year <- tapply(dynamic$year, key, min)
  required_key <- paste(expected$ecoregion, expected$species, sep = "\r")
  missing <- required_key[is.na(first_year[required_key]) |
                            first_year[required_key] > 1L]
  if (length(missing))
    bs52_climate_stop("La table dynamique climatique ne couvre pas toutes les especes et ecoregions actives des la premiere annee.")
  dynamic[order(dynamic$ecoregion, dynamic$species, dynamic$year), , drop = FALSE]
}

bs52_read_climate_dynamic52 <- function(path, canonical) {
  dynamic_file <- normalizePath(path, winslash = "/", mustWork = TRUE)
  dynamic <- utils::read.csv(dynamic_file, stringsAsFactors = FALSE,
                             check.names = FALSE)
  bs52_validate_climate_dynamic52(dynamic, canonical)
}

bs52_apply_climate52 <- function(canonical, climate) {
  if (is.null(climate$dynamic_file)) return(canonical)
  canonical$dynamic <- bs52_read_climate_dynamic52(climate$dynamic_file, canonical)
  canonical
}

bs52_import_landis_climate52 <- function(source_file, destination,
                                          canonical_directory) {
  destination <- normalizePath(destination, winslash = "/", mustWork = FALSE)
  if (file.exists(destination))
    bs52_climate_stop("Le fichier climatique existe deja : ", destination)
  source_file <- normalizePath(source_file, winslash = "/", mustWork = TRUE)
  canonical <- bs52_read_canonical52(canonical_directory)
  dynamic <- bs52_validate_climate_dynamic52(
    bs52_read_landis_dynamic52(source_file), canonical
  )
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(dynamic, destination, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
  manifest <- list(
    format = "landisR52/climate-dynamic-v1",
    source_file = source_file,
    source_md5 = unname(tools::md5sum(source_file)),
    created_at = Sys.time()
  )
  saveRDS(manifest, paste0(destination, ".manifest.rds"), compress = "gzip")
  destination
}
