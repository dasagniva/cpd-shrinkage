source(file.path("..", "..", "R", "10_dgp.R"))
source(file.path("..", "..", "R", "20_frontend_nsp.R"))
source(file.path("..", "..", "R", "21_frontend_capa.R"))
source(file.path("..", "..", "R", "22_prior_map.R"))

test_that("NSP front-end detects structural step change", {
  set.seed(101)
  # Series with one distinct jump
  y <- c(rnorm(100, 0, 1), rnorm(100, 4, 1))
  nsp_res <- run_frontend_nsp(y, alpha = 0.1, M = 500, min_size = 20)

  expect_named(nsp_res, c("intervals", "threshold", "method", "raw"))
  expect_true(is.data.frame(nsp_res$intervals))
  expect_true(nrow(nsp_res$intervals) >= 1)
  # Interval should cover the true jump at t = 101
  first_int <- nsp_res$intervals[1, ]
  expect_true(first_int$start <= 101 && first_int$end >= 101)
})

test_that("CAPA front-end identifies point and collective anomalies", {
  set.seed(202)
  n <- 200
  y <- rnorm(n)
  y[50] <- 8.0 # Point anomaly
  y[120:135] <- y[120:135] + 5.0 # Collective anomaly

  capa_res <- run_frontend_capa(y, prewhiten = FALSE, type = "mean", min_seg_len = 10)
  expect_true(50L %in% capa_res$point)
  expect_true(nrow(capa_res$collective) >= 1)
  expect_true(capa_res$collective$start[1] <= 125 && capa_res$collective$end[1] >= 130)
})

test_that("CAPA prewhitening robustly estimates phi and reduces false flags under AR(1)", {
  set.seed(303)
  n <- 300
  sim <- simulate_series(n = n, phi = 0.7, noise = "gaussian", k0 = 0, outlier_regime = "none", collective = 0, seed = 303)

  capa_pw <- run_frontend_capa(sim$y, prewhiten = TRUE, type = "mean")
  expect_true(abs(capa_pw$phi_hat - 0.7) < 0.25)
  expect_true(is.numeric(capa_pw$point))
})

test_that("build_prior_maps correctly maps intervals and anomalies to kappa and a", {
  n <- 100
  nsp_ints <- data.frame(start = c(20L), end = c(30L))
  capa_res <- list(
    point = c(75L),
    collective = data.frame(start = c(45L), end = c(60L))
  )

  pmap <- build_prior_maps(n = n, nsp_intervals = nsp_ints, capa_results = capa_res, kappa_in = 3.0, endpoint_buffer = 2L)

  expect_equal(length(pmap$kappa), n)
  expect_equal(length(pmap$a), n)

  # Check kappa_t values
  # Inside NSP interval [20, 30]
  expect_true(all(pmap$kappa[20:30] == 3.0))
  # Outside candidates
  expect_true(all(pmap$kappa[1:15] == 1.0))

  # Around collective endpoints (45 ± 2 and 60 ± 2)
  expect_true(all(pmap$kappa[43:47] == 3.0))
  expect_true(all(pmap$kappa[58:62] == 3.0))

  # Check point anomaly indicator a_t
  expect_equal(pmap$a[75], 1.0)
  expect_true(all(pmap$a[-75] == 0.0))

  # Ablation check: kappa_in = 1 switches NSP information off
  pmap_ablation <- build_prior_maps(n = n, nsp_intervals = nsp_ints, capa_results = capa_res, kappa_in = 1.0)
  expect_true(all(pmap_ablation$kappa == 1.0))
})
