test_that("LANDIS-II readers are independent and initialization assembles their inputs", {
  skip_if_not_installed("terra")
  inputs <- tempfile("succession-")
  dir.create(inputs)
  on.exit(unlink(inputs, recursive = TRUE), add = TRUE)

  writeLines(c(
    'LandisData "Scenario"', "Duration 20", "CellLength 100",
    "Species species.txt", "Ecoregions ecoregions.txt", "EcoregionsMap eco.tif",
    '"Biomass Succession" biomass.txt'
  ), file.path(inputs, "scenario.txt"))
  writeLines(c(
    'LandisData "Species"', "PINE 120 15 1 3 50 100 0 0 0 none"
  ), file.path(inputs, "species.txt"))
  writeLines(c(
    'LandisData "Ecoregions"', 'yes 1 north "Northern landtype"'
  ), file.path(inputs, "ecoregions.txt"))

  raster <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 200,
                       ymin = 0, ymax = 100)
  terra::writeRaster(terra::setValues(raster, c(1L, 1L)), file.path(inputs, "eco.tif"))
  landscape <- landisR::landis_read(inputs)
  expect_identical(landisR::landis_initialize(landscape)$year, 0L)
  expect_null(landscape$canonical$landscape$cells$community)

  writeLines(c(
    'LandisData "Biomass Succession"', "Timestep 10",
    "SeedingAlgorithm WardSeedDispersal", "SpinupMortalityFraction 0",
    "InitialCommunities initial.txt", "InitialCommunitiesMap initial.tif",
    "MinRelativeBiomass", "north", "1 0%", "2 20%", "3 40%", "4 60%", "5 80%",
    "SufficientLight", vapply(1:5, paste, "", "1 1 1 1 1 1"),
    "SpeciesParameters", "PINE 3 0.05 15 0.5 0.2",
    "EcoregionParameters", "north 600", "DynamicInputFile dynamic.txt",
    "FireReductionParameters", "1 0 0.5", "HarvestReductionParameters"
  ), file.path(inputs, "biomass.txt"))
  writeLines(c('LandisData "Initial Communities"', "MapCode 1", "PINE 30 (100)",
               "MapCode 2", "PINE 60 (500)"), file.path(inputs, "initial.txt"))
  writeLines(c('LandisData "Dynamic Input"', "0 north PINE 1 100 5000"),
             file.path(inputs, "dynamic.txt"))
  terra::writeRaster(terra::setValues(raster, c(1L, 2L)), file.path(inputs, "initial.tif"))

  unlink(file.path(inputs, c("species.txt", "ecoregions.txt", "eco.tif")))
  config <- biomass_succession_read(inputs)
  sim <- landisR::landis_initialize(landscape, seed = 1001L, succession = config)
  expect_identical(sim$seed, 1001L)
  expect_identical(sim$year, 0L)
  expect_identical(sim$project$canonical$timestep, 10L)
  expect_equal(landisR::landis_state(sim, representation = "cohorts")$biomass, c(100, 500))
  expect_identical(landisR::landis_sim(landscape, succession = config,
    steps = 1L, progress = FALSE)$simulation$year, 10L)

  wrong_map <- config
  wrong_map$initial_map <- matrix(1L, nrow = 1L, ncol = 1L)
  expect_error(landisR::landis_initialize(landscape, succession = wrong_map),
               "must match the landscape dimensions")
})
