# ==============================================================================
# 50_copula.R
# Copula Markov Chain Comparator (Phase P6).
# Implements Section 5 (Phase P6) and benchmarks tail-dependence vs sparse contamination.
# ==============================================================================

suppressPackageStartupMessages({
  library(rvinecopulib)
})

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
}

#' Tukey's Biweight Loss Function and Location M-Estimator
#'
#' @param e Numeric vector of residuals
#' @param c Tuning constant (default 4.685 for 95% Gaussian efficiency)
#'
#' @return Numeric vector of loss values
#' @export
biweight_loss <- function(e, c = 4.685) {
  u <- abs(e) / c
  val <- ifelse(u <= 1.0, (c^2 / 6.0) * (1.0 - (1.0 - u^2)^3), c^2 / 6.0)
  val
}

#' Fit Robust Piecewise-Constant Location
#'
#' Estimates piecewise location using robust loss (Tukey's biweight) to ensure
#' residuals are robust to outliers and mean shifts, independently of the proposed model.
#'
#' @param y Numeric vector of observations (length n)
#' @param min_seg_len Minimum segment length between changes (default 10)
#' @param penalty Penalty for adding a change point (default 2 * log(n))
#'
#' @return A list with elements:
#'   \item{mu_hat}{Numeric vector of length n of fitted piecewise constants}
#'   \item{residuals}{Numeric vector of residuals r_t = y_t - mu_hat_t}
#'   \item{cp}{Integer vector of estimated change point locations}
#' @export
fit_robust_location <- function(y, min_seg_len = 10, penalty = NULL) {
  n <- length(y)
  stopifnot(is.numeric(y), n >= min_seg_len)

  if (is.null(penalty)) {
    penalty <- 2.5 * log(n)
  }

  # Estimate robust scale of innovations
  sigma_rob <- stats::mad(diff(y)) / sqrt(2)
  if (sigma_rob <= 1e-8) sigma_rob <- stats::mad(y)
  if (sigma_rob <= 1e-8) sigma_rob <- 1.0

  # If robseg package is available, use Rob_FPOP
  if (requireNamespace("robseg", quietly = TRUE)) {
    try_rob <- tryCatch({
      fit_rfpop <- robseg::Rob_FPOP(y, beta = penalty * sigma_rob^2)
      cp_locs <- fit_rfpop$t.est
      mu_hat <- fit_rfpop$step_levels
      list(mu_hat = mu_hat, residuals = y - mu_hat, cp = cp_locs)
    }, error = function(e) NULL)
    if (!is.null(try_rob)) return(try_rob)
  }

  # Standalone robust binary segmentation with biweight loss
  # Helper: biweight location M-estimator for a segment
  eval_biweight_loc <- function(x) {
    med_x <- median(x)
    # 2-step iteratively reweighted biweight
    loc <- med_x
    c_val <- 4.685 * sigma_rob
    for (iter in 1:5) {
      r <- x - loc
      u <- abs(r) / c_val
      w <- ifelse(u < 1.0, (1.0 - u^2)^2, 0.0)
      if (sum(w) > 1e-6) {
        loc <- sum(w * x) / sum(w)
      } else {
        loc <- med_x
        break
      }
    }
    loc
  }

  eval_seg_cost <- function(s, e) {
    sub_y <- y[s:e]
    loc <- eval_biweight_loc(sub_y)
    sum(biweight_loss(sub_y - loc, c = 4.685 * sigma_rob))
  }

  # Recursive binary segmentation
  cps <- integer(0)

  find_cp <- function(s, e) {
    len <- e - s + 1L
    if (len < 2L * min_seg_len) return()

    base_cost <- eval_seg_cost(s, e)
    best_gain <- 0.0
    best_split <- NA_integer_

    # Candidate split points
    split_cands <- seq(s + min_seg_len, e - min_seg_len)
    for (sp in split_cands) {
      split_cost <- eval_seg_cost(s, sp - 1L) + eval_seg_cost(sp, e)
      gain <- base_cost - split_cost
      if (gain > best_gain) {
        best_gain <- gain
        best_split <- sp
      }
    }

    if (!is.na(best_split) && best_gain > penalty * (sigma_rob^2)) {
      cps <<- c(cps, best_split)
      find_cp(s, best_split - 1L)
      find_cp(best_split, e)
    }
  }

  find_cp(1L, n)
  cps <- sort(cps)

  # Compute piecewise levels
  mu_hat <- numeric(n)
  seg_bounds <- c(1L, cps, n + 1L)
  for (i in seq_len(length(seg_bounds) - 1L)) {
    st <- seg_bounds[i]
    en <- seg_bounds[i + 1L] - 1L
    mu_hat[st:en] <- eval_biweight_loc(y[st:en])
  }

  residuals <- y - mu_hat

  list(
    mu_hat    = mu_hat,
    residuals = residuals,
    cp        = cps
  )
}

#' Simulate Stationary First-Order Copula Markov Chain
#'
#' Simulates a sequence of uniform pseudo-observations U_1, ..., U_n
#' from a bivariate copula Markov chain using sequential inverse h-functions:
#'   U_1 ~ Unif(0, 1)
#'   U_t = h_1^{-1}(U_{t-1}, W_t) where W_t ~ Unif(0, 1) i.i.d.
#'
#' @param cop A `bicop_dist` or fitted `bicop` object from `rvinecopulib`
#' @param n Integer length of simulated sequence
#' @param seed Optional random seed
#'
#' @return Numeric vector of length n with values in (0, 1)
#' @export
simulate_copula_markov <- function(cop, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n >= 2)

  u <- numeric(n)
  u[1] <- stats::runif(1)
  w <- stats::runif(n - 1L)

  for (t in 2:n) {
    # Evaluate inverse h-function: h_1^{-1}(u_{t-1}, w_{t-1})
    u[t] <- as.numeric(rvinecopulib::hbicop(
      c(u[t - 1L], w[t - 1L]),
      cond_var = 1,
      cop,
      inverse  = TRUE
    ))
  }

  # Bound strictly away from 0 and 1
  pmax(1e-7, pmin(1.0 - 1e-7, u))
}

#' Calibrate Copula Tail Threshold via Parametric Bootstrap
#'
#' Simulates B bootstrap paths from the fitted copula Markov chain to find
#' the empirical threshold tau_alpha such that the false flagging rate is alpha:
#'   P(q_t <= tau_alpha) = alpha
#'
#' @param fit_copula Fitted `bicop` object
#' @param n Series length
#' @param alpha Target false flagging rate (default 0.05)
#' @param B_boot Number of bootstrap replications (default 200)
#' @param seed Optional random seed
#'
#' @return Numeric scalar threshold tau_alpha in (0, 1)
#' @export
calibrate_copula_threshold <- function(fit_copula, n, alpha = 0.05, B_boot = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(alpha > 0, alpha < 1, B_boot >= 10, n >= 10)

  # Generate bootstrap q-values
  all_q <- vector("list", B_boot)

  for (b in seq_len(B_boot)) {
    u_boot <- simulate_copula_markov(fit_copula, n)
    d_boot <- cbind(u_boot[1:(n - 1L)], u_boot[2:n])
    h_boot <- as.numeric(rvinecopulib::hbicop(d_boot, cond_var = 1, fit_copula))
    # Two-sided conditional tail probability
    q_boot <- 2.0 * pmin(h_boot, 1.0 - h_boot)
    # Add t=1 marginal tail probability
    q1 <- 2.0 * min(u_boot[1], 1.0 - u_boot[1])
    all_q[[b]] <- c(q1, q_boot)
  }

  flat_q <- unlist(all_q, use.names = FALSE)
  # Empirical quantile at level alpha
  tau_alpha <- as.numeric(stats::quantile(flat_q, probs = alpha, na.rm = TRUE, names = FALSE))
  tau_alpha
}

#' Fit Copula Markov Chain Comparator
#'
#' Implements the 6-step comparator specified in Section 5 (Phase P6):
#' 1. Robust location fit (biweight FPOP) -> residuals r_t
#' 2. Pseudo-observations U_t = rank(r_t) / (n + 1)
#' 3. First-order copula Markov fit via `rvinecopulib::bicop` (candidate families:
#'    Gaussian, t, Clayton, Gumbel, Frank, BB1, selected by AIC)
#' 4. Conditional tail scores: h_t = C_{2|1}(U_t | U_{t-1}), q_t = 2 min(h_t, 1 - h_t)
#' 5. Calibrate threshold tau_alpha via parametric bootstrap at level alpha
#' 6. Extract lower and upper tail dependence coefficients lambda_L, lambda_U
#'
#' @param y Numeric vector of time series observations
#' @param alpha Target nominal false-flagging rate (default 0.05)
#' @param family_set Candidate copula families (default c("gaussian", "t", "clayton", "gumbel", "frank", "bb1"))
#' @param B_boot Number of bootstrap iterations for threshold calibration (default 200)
#' @param seed Optional random seed
#'
#' @return A validated `cpd_fit` object of class `cpd_fit`
#' @export
fit_copula_comparator <- function(
    y,
    alpha = 0.05,
    family_set = c("gaussian", "t", "clayton", "gumbel", "frank", "bb1"),
    B_boot = 200,
    seed = NULL
) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 15)

  if (!is.null(seed)) set.seed(seed)

  # Step 1: Location fit
  loc_fit <- fit_robust_location(y)
  mu_hat <- loc_fit$mu_hat
  r_t <- loc_fit$residuals
  cp_locs <- loc_fit$cp

  # Step 2: Pseudo-observations
  # U_t = rank(r_t) / (n + 1)
  U_t <- rank(r_t) / (n + 1.0)

  # Step 3: Copula fit
  # Pairs (U_{t-1}, U_t)
  D_pairs <- cbind(U_t[1:(n - 1L)], U_t[2:n])

  fit_cop <- rvinecopulib::bicop(
    data        = D_pairs,
    family_set  = family_set,
    selcrit     = "aic",
    par_method  = "mle"
  )

  # Step 4: Conditional scores
  # h_t = C_{2|1}(U_t | U_{t-1})
  h_pairs <- as.numeric(rvinecopulib::hbicop(D_pairs, cond_var = 1, fit_cop))
  # Two-sided conditional tail probabilities
  q_pairs <- 2.0 * pmin(h_pairs, 1.0 - h_pairs)

  # Score vector for t = 1..n
  q_vec <- numeric(n)
  q_vec[1] <- 2.0 * min(U_t[1], 1.0 - U_t[1])
  q_vec[2:n] <- q_pairs

  # Step 5: Bootstrap threshold calibration
  tau_alpha <- calibrate_copula_threshold(
    fit_copula = fit_cop,
    n          = n,
    alpha      = alpha,
    B_boot     = B_boot,
    seed       = seed
  )

  # Flagged point anomalies: observations with extreme conditional tail probability
  flagged_pts <- which(q_vec <= tau_alpha)

  # Step 6: Tail-dependence coefficients
  td_mat <- rvinecopulib::tail_dep(fit_cop)
  lambda_L <- as.numeric(td_mat["lower", "lower"])
  lambda_U <- as.numeric(td_mat["upper", "upper"])

  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  # Construct four-class classification vector
  class_vec <- rep("nominal", n)
  if (length(cp_locs) > 0) {
    class_vec[cp_locs] <- "shift"
  }
  if (length(flagged_pts) > 0) {
    pt_only <- setdiff(flagged_pts, cp_locs)
    if (length(pt_only) > 0) {
      class_vec[pt_only] <- "point"
    }
  }

  scores_df <- data.frame(
    t     = seq_len(n),
    s_cp  = rep(NA_real_, n),
    s_out = pmax(0.0, pmin(1.0, 1.0 - q_vec)) # Outlier score in [0, 1]: higher = more anomalous
  )

  new_cpd_fit(
    method      = "copula_markov",
    cp          = as.integer(cp_locs),
    cp_int      = NULL,
    point       = as.integer(setdiff(flagged_pts, cp_locs)),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = scores_df,
    runtime_sec = runtime_sec,
    diag        = list(
      copula_fit         = fit_cop,
      family             = fit_cop$family,
      rotation           = fit_cop$rotation,
      parameters         = fit_cop$parameters,
      lambda_L           = lambda_L,
      lambda_U           = lambda_U,
      threshold          = tau_alpha,
      q_t                = q_vec,
      pseudo_obs         = U_t,
      mu_hat             = mu_hat,
      residuals          = r_t
    )
  )
}
