source(file.path("..", "..", "R", "10_dgp.R"))

test_that("simulate_series outputs expected structure and mutually exclusive labels", {
  sim <- simulate_series(
    n = 1000,
    phi = 0.6,
    noise = "gaussian",
    k0 = 5,
    jump = 2,
    outlier_regime = "isolated",
    p_outlier = 0.02,
    collective = 2,
    seed = 42
  )

  expect_named(sim, c("y", "truth"))
  expect_equal(length(sim$y), 1000)
  expect_equal(length(sim$truth$class), 1000)
  expect_equal(levels(sim$truth$class), c("nominal", "shift", "collective", "point"))

  # Check that labels never overlap
  # Every index in cp must be "shift"
  expect_true(all(sim$truth$class[sim$truth$cp] == "shift"))
  # Every index in point must be "point"
  expect_true(all(sim$truth$class[sim$truth$point] == "point"))
  # Every index in collective must be "collective"
  if (nrow(sim$truth$collective) > 0) {
    for (r in seq_len(nrow(sim$truth$collective))) {
      seg <- sim$truth$collective$start[r]:sim$truth$collective$end[r]
      expect_true(all(sim$truth$class[seg] == "collective"))
    }
  }

  # Mutual exclusivity of index sets
  coll_indices <- integer(0)
  if (nrow(sim$truth$collective) > 0) {
    coll_indices <- unlist(apply(sim$truth$collective, 1, function(row) row[1]:row[2]))
  }
  expect_equal(length(intersect(sim$truth$cp, sim$truth$point)), 0)
  expect_equal(length(intersect(sim$truth$cp, coll_indices)), 0)
  expect_equal(length(intersect(sim$truth$point, coll_indices)), 0)

  # Check decomposition: y == mu + psi + eps
  expect_equal(sim$y, sim$truth$mu + sim$truth$psi + sim$truth$eps, tolerance = 1e-12)
})

test_that("clean AR(1) noise sample ACF matches theoretical phi within MC tolerance", {
  phi_target <- 0.7
  n_pts <- 10000
  sim <- simulate_series(
    n = n_pts,
    phi = phi_target,
    noise = "gaussian",
    k0 = 0,
    outlier_regime = "none",
    collective = 0,
    seed = 12345
  )

  # Sample lag-1 autocorrelation
  acf_lag1 <- cor(sim$truth$eps[-1], sim$truth$eps[-n_pts])
  # SE of sample ACF is approx sqrt((1 - phi^2) / n) for lag 1
  mc_se <- sqrt((1 - phi_target^2) / n_pts)

  expect_equal(acf_lag1, phi_target, tolerance = 4 * mc_se)
})

test_that("simulate_series supports t3 innovations and clustered outliers", {
  sim_t3 <- simulate_series(
    n = 600,
    phi = 0.4,
    noise = "t3",
    k0 = 3,
    outlier_regime = "clustered",
    collective = 1,
    seed = 999
  )

  expect_equal(length(sim_t3$y), 600)
  expect_equal(length(sim_t3$truth$cp), 3)
  expect_true(length(sim_t3$truth$point) >= 1)
  expect_true(nrow(sim_t3$truth$collective) == 1)
})
