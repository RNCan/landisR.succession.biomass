test_that("strict serotiny uses corrected Ward on coarse grids", {

  species <- data.frame(
    species = "PINU.BAN",
    longevity = 120L,
    maturity = 15L,
    shade_tolerance = 1L,
    fire_tolerance = 3L,
    effective_seed_distance = 10,
    max_seed_distance = 30,
    vegetative_probability = 0,
    min_sprout_age = 0L,
    max_sprout_age = 0L,
    post_fire_regeneration = "serotiny",
    serotiny_type = "strict",
    serotiny_dispersal = "kernel",
    stringsAsFactors = FALSE
  )

  species_parameters <- data.frame(
    species = "PINU.BAN",
    leaf_longevity = 3,
    woody_decay = 0.05,
    mortality_shape = 15,
    growth_shape = 0.5,
    leaf_lignin = 0.2,
    stringsAsFactors = FALSE
  )

  ecoregions <- data.frame(
    map_code = 1L,
    ecoregion = "eco1",
    active = TRUE,
    description = "test",
    stringsAsFactors = FALSE
  )

  ecoregion_parameters <- data.frame(
    ecoregion = "eco1",
    aet = 600,
    stringsAsFactors = FALSE
  )

  dynamic <- data.frame(
    year = 0L,
    ecoregion = "eco1",
    species = "PINU.BAN",
    sep = 1,
    max_anpp = 100,
    max_biomass = 5000,
    stringsAsFactors = FALSE
  )

  initial_cohorts <- data.frame(
    map_code = 1L,
    species = "PINU.BAN",
    age = 30L,
    biomass = 100,
    stringsAsFactors = FALSE
  )

  shade <- data.frame(
    shade = 1:5,
    eco1 = c(0, 0.2, 0.4, 0.6, 0.8)
  )

  light <- data.frame(
    shade_tolerance = 1:5,
    shade_0 = 1, shade_1 = 1, shade_2 = 1,
    shade_3 = 1, shade_4 = 1, shade_5 = 1
  )

  fire_reductions <- data.frame(
    severity = 1L,
    wood_reduction = 0,
    litter_reduction = 0.5
  )

  nr <- 9L
  nc <- 9L

  grid <- expand.grid(
    row = seq_len(nr),
    col = seq_len(nc)
  )

  grid$cell_id <- with(
    grid,
    row + (col - 1L) * nr
  )

  grid$ecoregion_map_code <- 1L
  grid$ecoregion <- "eco1"
  grid$active <- TRUE
  grid$community <- 2L

  centre_row <- 5L
  centre_col <- 5L
  centre_index <- which(
    grid$row == centre_row & grid$col == centre_col
  )
  centre_cell <- grid$cell_id[centre_index]
  grid$community[centre_index] <- 1L

  project <- landis_build(
    landscape = grid[
      c(
        "cell_id", "row", "col", "ecoregion_map_code",
        "ecoregion", "active", "community"
      )
    ],
    species = species,
    species_parameters = species_parameters,
    ecoregions = ecoregions,
    ecoregion_parameters = ecoregion_parameters,
    dynamic = dynamic,
    initial_cohorts = initial_cohorts,
    shade = shade,
    light = light,
    fire_reductions = fire_reductions,
    duration = 10L,
    cell_length = 250,
    timestep = 10L,
    seeding_algorithm = "WardSeedDispersal",
    spinup_cohorts = FALSE,
    spinup_mortality_fraction = 0,
    initial_mode = "explicit_biomass",
    community_codes = c(1L, 2L),
    name = "serotiny-ward-integration"
  )

  # Strict serotiny must not recruit through ordinary seeding.
  sim_no_fire <- landis_initialize(project, seed = 4357L)
  sim_no_fire <- landis_step(sim_no_fire)
  no_fire <- landis_state(sim_no_fire, representation = "cohorts")

  expect_false(any(
    no_fire$species == "PINU.BAN" &
      no_fire$age <= 10L
  ))

  # Fire-triggered recruitment must use the corrected Ward kernel.
  sim <- landis_initialize(project, seed = 4357L)

  species_row <- sim$ctx$config$species[
    sim$ctx$config$species$species == "PINU.BAN",
    ,
    drop = FALSE
  ]

  arrival <- getFromNamespace(
    "landisR_serotiny_arrival",
    "landisR.succession.biomass"
  )(
    sim,
    source_cells = centre_cell,
    species_row = species_row
  )

  sim <- landis_fire(
    sim,
    cells = centre_cell,
    update = landis_update(clear_cells = centre_cell),
    severity = 5,
    stochastic = FALSE
  )

  state <- landis_state(sim, representation = "cohorts")
  recruits <- unique(
    state$cell_id[
      state$species == "PINU.BAN" & state$age == 1L
    ]
  )

  grid$arrival <- as.numeric(arrival)
  grid$recruited <- grid$cell_id %in% recruits

  cardinal <- (
    (abs(grid$row - centre_row) == 1L & grid$col == centre_col) |
      (abs(grid$col - centre_col) == 1L & grid$row == centre_row)
  )

  diagonal <- (
    abs(grid$row - centre_row) == 1L &
      abs(grid$col - centre_col) == 1L
  )

  expect_equal(sum(cardinal), 4L)
  expect_true(all(grid$arrival[cardinal] > 0))
  expect_true(all(grid$recruited[cardinal]))

  expect_equal(sum(diagonal), 4L)
  expect_true(all(grid$arrival[diagonal] == 0))
  expect_false(any(grid$recruited[diagonal]))
})
