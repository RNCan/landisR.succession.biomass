biomass_succession_read <- function(source_directory, scenario = "scenario.txt",
                                    prefer_cells_csv = TRUE) {
  source_directory <- normalizePath(source_directory, winslash = "/", mustWork = TRUE)
  scenario_file <- chartr("\\", "/", scenario)
  if (!grepl("^(/|[A-Za-z]:/)", scenario_file))
    scenario_file <- file.path(source_directory, scenario_file)
  scenario_file <- normalizePath(scenario_file, winslash = "/", mustWork = TRUE)
  inputs <- bs52_read_mauricie(scenario_file, prefer_cells_csv)
  inputs$scenario_file <- scenario_file
  structure(inputs, class = "landisR_biomass_succession")
}

landis_initialize_succession.landisR_biomass_succession <- function(
    succession, project, seed, ...) {
  assembled <- bs52_assemble_landis_inputs(project, succession)
  landisR::landis_initialize(assembled, seed = seed, ...)
}
