test_that("Ward neighborhood always includes first cardinal ring", {
  f <- getFromNamespace("bs52_ward_neighborhood", "landisR.succession.biomass")
  n <- f(cell_length = 250, max_seed_distance = 30)
  expect_true(any(n$row == 1L & n$col == 0L & n$distance == 250))
})

test_that("Ward 2019 allows adjacent cardinal dispersal when MaxD < cell size", {
  neighborhood <- getFromNamespace("bs52_ward_neighborhood", "landisR.succession.biomass")(
    250, 30
  )
  prob_fun <- getFromNamespace("bs52_ward_probability", "landisR.succession.biomass")

  active <- matrix(TRUE, 3, 3)
  mature <- matrix(FALSE, 3, 3)
  mature[2, 3] <- TRUE
  target <- matrix(c(2L, 2L), nrow = 1L)

  p <- prob_fun(
    mature, active, target,
    effective = 10, maximum = 30, cell_length = 250,
    neighborhood = neighborhood, variant = "ward2019"
  )

  expect_length(p, 1L)
  expect_gt(p, 0)
  expect_lte(p, 1)
})

test_that("short-distance species does not inherit another species long extent", {
  # Deliberately build a neighborhood large enough for another species.
  neighborhood <- getFromNamespace("bs52_ward_neighborhood", "landisR.succession.biomass")(
    250, 1000
  )
  prob_fun <- getFromNamespace("bs52_ward_probability", "landisR.succession.biomass")

  active <- matrix(TRUE, 9, 9)
  target <- matrix(c(5L, 5L), nrow = 1L)

  # Cardinal adjacent source: allowed despite MaxD = 30 m.
  mature_cardinal <- matrix(FALSE, 9, 9)
  mature_cardinal[5, 6] <- TRUE
  p_cardinal <- prob_fun(
    mature_cardinal, active, target,
    effective = 10, maximum = 30, cell_length = 250,
    neighborhood = neighborhood, variant = "ward2019"
  )
  expect_gt(p_cardinal, 0)

  # Diagonal adjacent source is 353.6 m away and remains outside the
  # one-cell cardinal safeguard for this very short-distance species.
  mature_diagonal <- matrix(FALSE, 9, 9)
  mature_diagonal[6, 6] <- TRUE
  p_diagonal <- prob_fun(
    mature_diagonal, active, target,
    effective = 10, maximum = 30, cell_length = 250,
    neighborhood = neighborhood, variant = "ward2019"
  )
  expect_equal(p_diagonal, 0)

  # A farther source also remains impossible.
  mature_far <- matrix(FALSE, 9, 9)
  mature_far[5, 7] <- TRUE
  p_far <- prob_fun(
    mature_far, active, target,
    effective = 10, maximum = 30, cell_length = 250,
    neighborhood = neighborhood, variant = "ward2019"
  )
  expect_equal(p_far, 0)
})

test_that("standard Ward neighborhoods are unchanged when MaxD is large", {
  f <- getFromNamespace("bs52_ward_neighborhood", "landisR.succession.biomass")
  new <- f(100, 700)

  radius_old <- 700 + 100 / 2
  n_old <- trunc(radius_old / 100) + 1L
  old <- expand.grid(row = seq_len(n_old), col = 0:n_old)
  old$distance <- 100 * sqrt(old$row^2 + old$col^2)
  old <- old[old$distance <= radius_old, ]
  old <- old[order(old$distance, old$row, old$col), ]
  rownames(old) <- NULL
  rownames(new) <- NULL

  expect_identical(new, old)
})

test_that("large52 Ward table carries the adjacent-cell safeguard", {
  neighborhood <- getFromNamespace("bs52_ward_neighborhood", "landisR.succession.biomass")(
    250, 30
  )
  table_fun <- getFromNamespace("bs52_large_ward_table52", "landisR.succession.biomass")

  config <- list(
    cell_length = 250,
    species = data.frame(
      species = "TEST",
      effective_seed_distance = 10,
      max_seed_distance = 30,
      stringsAsFactors = FALSE
    )
  )

  tab <- table_fun(config, neighborhood, "ward2019")
  cardinal <- which(neighborhood$row == 1L & neighborhood$col == 0L)
  expect_length(cardinal, 1L)
  expect_lt(tab$log_failure[cardinal, 1L], 0)
})

test_that("R and Rcpp Ward agree for the coarse-grid edge case", {
  skip_if_not(
    exists("bs52_ward_failure_log_cpp", mode = "function", inherits = TRUE),
    "Rcpp Ward engine not loaded in this test session"
  )

  neighborhood <- getFromNamespace("bs52_ward_neighborhood", "landisR.succession.biomass")(
    250, 30
  )
  r_fun <- getFromNamespace("bs52_ward_probability", "landisR.succession.biomass")
  cpp_fun <- getFromNamespace("bs52_ward_probability_rcpp", "landisR.succession.biomass")

  active <- matrix(TRUE, 3, 3)
  mature <- matrix(FALSE, 3, 3)
  mature[2, 3] <- TRUE
  target <- matrix(c(2L, 2L), nrow = 1L)

  args <- list(
    mature = mature, active = active, targets = target,
    effective = 10, maximum = 30, cell_length = 250,
    neighborhood = neighborhood, variant = "ward2019"
  )

  p_r <- do.call(r_fun, args)
  p_cpp <- do.call(cpp_fun, args)
  expect_equal(p_cpp, p_r, tolerance = 1e-12)
  expect_gt(p_cpp, 0)
})
