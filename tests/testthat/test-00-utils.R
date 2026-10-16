source(file.path("..", "..", "R", "00_utils.R"))

test_that("new_cpd_fit constructs expected default object", {
  fit <- new_cpd_fit(method = "test_method")
  expect_s3_class(fit, "cpd_fit")
  expect_true(is.list(fit))
  expect_equal(fit$method, "test_method")
  expect_equal(fit$cp, integer(0))
  expect_equal(nrow(fit$cp_int), 0)
  expect_equal(fit$point, integer(0))
  expect_equal(nrow(fit$collective), 0)
  expect_equal(length(fit$class), 0)
  expect_equal(levels(fit$class), CPD_CLASS_LEVELS)
  expect_equal(nrow(fit$scores), 0)
  expect_equal(fit$runtime_sec, 0)
  expect_equal(fit$diag, list())
})

test_that("validate_cpd_fit validates correct fit without error", {
  n <- 100
  class_vec <- factor(rep("nominal", n), levels = CPD_CLASS_LEVELS)
  class_vec[25] <- "shift"
  class_vec[50:60] <- "collective"
  class_vec[80] <- "point"

  scores_df <- data.frame(
    t = 1:n,
    s_cp = runif(n),
    s_out = runif(n)
  )

  fit <- new_cpd_fit(
    method = "proposed",
    cp = 25L,
    cp_int = data.frame(start = 24L, end = 26L),
    point = 80L,
    collective = data.frame(start = 50L, end = 60L),
    class = class_vec,
    scores = scores_df,
    runtime_sec = 1.23,
    diag = list(rhat = 1.002, ess = 450)
  )

  expect_silent(validate_cpd_fit(fit, n = n))
})

test_that("validate_cpd_fit catches invalid structure or out-of-bounds indices", {
  # Not a cpd_fit
  expect_error(validate_cpd_fit(list(a = 1)), "Object must inherit from 'cpd_fit'")

  # Missing field
  bad_fit <- new_cpd_fit("bad")
  bad_fit$runtime_sec <- NULL
  expect_error(validate_cpd_fit(bad_fit), "Missing required fields")

  # Wrong class levels
  bad_fit <- new_cpd_fit("bad")
  bad_fit$class <- factor(c("a", "b"))
  expect_error(validate_cpd_fit(bad_fit), "must be a factor with levels")

  # Inverted interval
  bad_fit <- new_cpd_fit("bad", collective = data.frame(start = 50L, end = 40L))
  expect_error(validate_cpd_fit(bad_fit), "start > end")

  # Out of bounds for given n
  fit <- new_cpd_fit(
    method = "test",
    cp = 150L,
    class = factor(rep("nominal", 100), levels = CPD_CLASS_LEVELS)
  )
  expect_error(validate_cpd_fit(fit, n = 100), "outside \\[1, n\\]")
})

test_that("init_rng_stream sets L'Ecuyer-CMRG and generates reproducible distinct streams", {
  init_rng_stream(cell_id = 1, rep_id = 1, base_seed = 42)
  expect_equal(RNGkind()[1], "L'Ecuyer-CMRG")
  draw_c1_r1_a <- runif(5)

  # Re-run same stream
  init_rng_stream(cell_id = 1, rep_id = 1, base_seed = 42)
  draw_c1_r1_b <- runif(5)
  expect_equal(draw_c1_r1_a, draw_c1_r1_b)

  # Advance replicate
  init_rng_stream(cell_id = 1, rep_id = 2, base_seed = 42)
  draw_c1_r2 <- runif(5)
  expect_false(identical(draw_c1_r1_a, draw_c1_r2))

  # Advance cell
  init_rng_stream(cell_id = 2, rep_id = 1, base_seed = 42)
  draw_c2_r1 <- runif(5)
  expect_false(identical(draw_c1_r1_a, draw_c2_r1))
  expect_false(identical(draw_c1_r2, draw_c2_r1))
})
