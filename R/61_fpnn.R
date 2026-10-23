# ==============================================================================
# 61_fpnn.R
# JISPS two-phase EV-DeCAFS + Fourier Probabilistic Neural Network (FPNN), ported from
# ev-decafs-serve (src/evdecafs_serve: features/ar1.py, features/evt.py,
# models/decafs.py, features/extract.py, training/bocpd.py, training/labelling.py,
# training/smote.py, models/fpnn.py, training/pipeline.py), itself distilled from
# changepoint-evdecafs.
#
#   Phase I : AR(1) estimates -> flat BIC penalty C log n -> EV-DeCAFS grid DP
#   Phase II: features (delta_mu, S, phi_local, V, xi_local) at each Phase-I change
#             -> BOCPD-oracle labels refined by the Algorithm-3 rule -> SMOTE
#             -> FourierPNN fit -> classify each change as sustained (1) / recoiled (0)
#
# Indices are 1-based throughout: a change at t means mu_t != mu_{t-1}. Two deliberate
# differences from the Python code: SMOTE uses R's RNG (so it is not bit-identical; the
# serving repo makes the same choice, DECISIONS.md #6), and the grid DP is in C++
# (src/ev_decafs_dp.cpp) with the same recursion and tie-breaking.
# ==============================================================================

local({
  root <- if (exists("find_project_root")) find_project_root() else getwd()
  cpp <- file.path(root, "src", "ev_decafs_dp.cpp")
  if (!exists("ev_decafs_dp_cpp", mode = "function")) {
    Rcpp::sourceCpp(cpp, env = globalenv(), cacheDir = file.path(root, ".rcpp_cache"))
  }
})

np_var <- function(x) mean((x - mean(x))^2) # numpy.var (denominator n)

#' AR(1) parameters for EV-DeCAFS Phase I (features/ar1.py: estimate_ar1_params)
#'
#' phi by Yule-Walker (MLE autocovariances) on the mean-centred series; sigma_v^2 from the
#' AR(1) residuals; sigma_eta^2 from the excess variance of first differences, floored at
#' 1e-3 sigma_v^2.
#'
#' @param y Numeric series
#' @return list(phi, sigma_v_sq, sigma_eta_sq)
#' @keywords internal
fpnn_ar1_params <- function(y) {
  stopifnot(length(y) >= 3)
  yc <- y - mean(y)
  n <- length(yc)
  phi <- sum(yc[-n] * yc[-1]) / sum(yc^2)
  phi <- min(0.999, max(-0.999, phi))
  sigma_v_sq <- stats::var(y[-1] - phi * y[-n])
  raw <- stats::var(diff(y)) - 2.0 * sigma_v_sq / (1.0 + abs(phi) + 1e-12)
  list(phi = phi, sigma_v_sq = sigma_v_sq, sigma_eta_sq = max(raw, 1e-3 * sigma_v_sq))
}

#' Local extreme value index field (features/evt.py: compute_evi_field)
#'
#' Method-of-moments GPD shape of the exceedances of |y_s - mean(W_t)| over their q0
#' quantile, in a reflect-padded window of half-width w; 0 when there are fewer than
#' `min_exceedances` exceedances. Clipped to [-1, 2].
#'
#' @keywords internal
fpnn_evi_field <- function(y, w = 50L, q0 = 0.90, min_exceedances = 5L) {
  n <- length(y)
  stopifnot(n > w)
  # numpy.pad(mode = "reflect"): the edge value is not repeated
  yp <- c(y[(w + 1L):2L], y, y[(n - 1L):(n - w)])
  xi <- numeric(n)
  for (t in seq_len(n)) {
    win <- yp[t:(t + 2L * w)]
    dev <- abs(win - mean(win))
    thr <- stats::quantile(dev, q0, names = FALSE, type = 7)
    exc <- dev[dev > thr] - thr
    if (length(exc) < min_exceedances) next
    m <- mean(exc)
    if (m <= 0) next
    ratio <- np_var(exc) / (m * m)
    xi[t] <- min(2.0, max(-1.0, 0.5 * (1.0 - 1.0 / max(ratio, 1e-10))))
  }
  xi
}

#' Change-point features (features/extract.py: extract_features, Algorithm 2)
#'
#' For a change at tau: D- = y[tau - L, ..., tau - 1], D+ = y[tau + 1, ..., tau + L].
#' delta_mu = mean(D+) - mean(D-); S = share of D+ within epsilon of mean(D+); phi_local =
#' lag-1 coefficient of y - means over D+; V = var(D-) / var(D+); xi_local = mean |xi| over
#' tau - L, ..., tau + L.
#'
#' @keywords internal
fpnn_extract_features <- function(y, cps, means, L = 5L, epsilon = NULL, xi_field = NULL) {
  n <- length(y)
  if (is.null(epsilon)) epsilon <- max(1e-10, stats::median(abs(y - means)))
  feat_names <- c("delta_mu", "S", "phi_local", "V", if (!is.null(xi_field)) "xi_local")
  X <- matrix(NA_real_, nrow = length(cps), ncol = length(feat_names), dimnames = list(NULL, feat_names))

  for (i in seq_along(cps)) {
    tau <- cps[i]
    D_minus <- y[max(1L, tau - L):(tau - 1L)]
    hi <- min(n, tau + L)
    idx_plus <- if (tau + 1L <= hi) (tau + 1L):hi else integer(0)
    D_plus <- y[idx_plus]

    mu_minus <- if (length(D_minus) > 0) mean(D_minus) else 0.0
    mu_plus <- if (length(D_plus) > 0) mean(D_plus) else 0.0
    S <- if (length(D_plus) > 0) mean(abs(D_plus - mu_plus) <= epsilon) else 0.0

    r <- y[idx_plus] - means[idx_plus]
    phi_local <- 0.0
    if (length(r) >= 2) {
      den <- sum(r[-length(r)]^2)
      if (abs(den) > 1e-10) phi_local <- min(0.999, max(-0.999, sum(r[-length(r)] * r[-1]) / den))
    }

    var_minus <- if (length(D_minus) > 1) np_var(D_minus) else 0.0
    var_plus <- if (length(D_plus) > 1) np_var(D_plus) else 0.0
    row <- c(mu_plus - mu_minus, S, phi_local, var_minus / (var_plus + 1e-10))
    if (!is.null(xi_field)) row <- c(row, mean(abs(xi_field[max(1L, tau - L):min(n, tau + L)])))
    X[i, ] <- row
  }
  X
}

#' BOCPD labelling oracle (training/bocpd.py: run_bocpd)
#'
#' Adams-MacKay BOCPD with a Normal-Gamma prior on the AR(1)-prewhitened series, constant
#' hazard 1/n, run lengths truncated to the 500 most probable. Flags t when
#' P(r_t = 0 | y_1:t) > threshold.
#'
#' @keywords internal
fpnn_bocpd <- function(y, phi, threshold = 0.3, kappa0 = 1.0, alpha0 = 1.0) {
  n <- length(y)
  if (n < 2) return(logical(n))
  z <- c(y[1], y[-1] - phi * y[-n])
  mu0 <- mean(z)
  beta0 <- max(np_var(z) / 2.0, 1e-6)
  H <- 1.0 / n

  log_t <- function(x, df, mu, s2) {
    q <- (x - mu)^2 / (s2 + 1e-12)
    lgamma(0.5 * (df + 1)) - lgamma(0.5 * df) - 0.5 * log(pi * df * s2 + 1e-300) -
      0.5 * (df + 1) * log1p(q / (df + 1e-12))
  }
  df0 <- 2.0 * alpha0
  s2_0 <- max(beta0 * (kappa0 + 1) / (alpha0 * kappa0 + 1e-12), 1e-12)

  R <- 1.0
  kap <- kappa0
  mu <- mu0
  al <- alpha0
  be <- beta0
  cp_prob <- numeric(n)

  for (t in 2:n) {
    x <- z[t]
    s2 <- pmax(be * (kap + 1) / (al * kap + 1e-12), 1e-12)
    lp <- log_t(x, 2.0 * al, mu, s2)
    lp_max <- max(lp)
    pred_prior <- exp(log_t(x, df0, mu0, s2_0) - lp_max)
    R_new <- c(H * pred_prior, R * exp(lp - lp_max) * (1 - H))
    Z <- sum(R_new)
    R_new <- if (Z > 1e-300) R_new / Z else c(1, rep(0, length(R_new) - 1L))
    cp_prob[t] <- R_new[1]

    be <- c(beta0, be + 0.5 * kap / (kap + 1) * (x - mu)^2)
    mu <- c(mu0, (kap * mu + x) / (kap + 1))
    kap <- c(kappa0, kap + 1)
    al <- c(alpha0, al + 0.5)
    R <- R_new

    if (length(R) > 500L) {
      keep <- sort(order(R)[(length(R) - 499L):length(R)])
      R <- R[keep] / sum(R[keep])
      mu <- mu[keep]
      kap <- kap[keep]
      al <- al[keep]
      be <- be[keep]
    }
  }
  cp_prob > threshold
}

#' Phase-II training labels without ground truth (training/pipeline.py: _make_labels)
#'
#' A Phase-I change within `tol` of a BOCPD flag is pending and resolved by Algorithm 3
#' (sustained iff |delta_mu| > kappa_mu and S > kappa_S, kappa_mu the given percentile of
#' |delta_mu|); other changes are recoiled. Falls back to the Algorithm-3 labels when this
#' gives a single class.
#'
#' @keywords internal
fpnn_make_labels <- function(X, cps, bocpd_cps, tol, kappa_mu_percentile = 75, kappa_S = 0.5) {
  kappa_mu <- stats::quantile(abs(X[, "delta_mu"]), kappa_mu_percentile / 100, names = FALSE, type = 7)
  feature_labels <- as.integer(abs(X[, "delta_mu"]) > kappa_mu & X[, "S"] > kappa_S)
  near_bocpd <- vapply(cps, function(tau) length(bocpd_cps) > 0 && min(abs(bocpd_cps - tau)) <= tol, logical(1))
  labels <- ifelse(near_bocpd, feature_labels, 0L)
  if (length(unique(labels)) < 2) labels <- feature_labels
  as.integer(labels)
}

#' SMOTE oversampling of every non-majority class (training/smote.py)
#'
#' @keywords internal
fpnn_smote <- function(X, labels, k = 5L) {
  classes <- sort(unique(labels))
  if (length(classes) < 2) return(list(X = X, y = labels))
  counts <- table(factor(labels, levels = classes))
  X_out <- X
  y_out <- labels
  for (cls in classes) {
    n_needed <- max(counts) - counts[[as.character(cls)]]
    if (n_needed <= 0) next
    X_cls <- X[labels == cls, , drop = FALSE]
    k_eff <- min(k, nrow(X_cls) - 1L)
    if (k_eff < 1) next
    dist <- as.matrix(stats::dist(X_cls))
    diag(dist) <- Inf
    nbrs <- t(apply(dist, 1, function(d) order(d)[seq_len(k_eff)]))
    if (k_eff == 1) nbrs <- matrix(nbrs, ncol = 1)
    base <- sample.int(nrow(X_cls), n_needed, replace = TRUE)
    pick <- sample.int(k_eff, n_needed, replace = TRUE)
    gap <- stats::runif(n_needed)
    anchors <- X_cls[base, , drop = FALSE]
    partners <- X_cls[nbrs[cbind(base, pick)], , drop = FALSE]
    X_out <- rbind(X_out, anchors + gap * (partners - anchors))
    y_out <- c(y_out, rep(cls, n_needed))
  }
  list(X = X_out, y = y_out)
}

#' Fourier PNN training (models/fpnn.py: FourierPNN.fit, Algorithm 4)
#'
#' Features are min-max scaled to `scaling_range`; for each class and feature the
#' Fejer-weighted coefficients are A_j = w_j sum cos(pi j z), B_j = w_j sum sin(pi j z),
#' w_j = (J + 1 - j) / (N_c (J + 1)).
#'
#' @export
fourier_pnn_fit <- function(X, y, J = 10L, scaling_range = c(-0.5, 0.5)) {
  X <- as.matrix(X)
  data_min <- apply(X, 2, min)
  data_range <- apply(X, 2, max) - data_min
  data_range[data_range < 10 * .Machine$double.eps] <- 1.0 # sklearn _handle_zeros_in_scale
  model <- list(J = J, scaling_range = scaling_range, data_min = data_min, data_range = data_range,
                classes = sort(unique(y)), n_samples = length(y))
  Z <- fourier_pnn_scale(model, X)
  j <- seq_len(J)
  fejer <- (J + 1 - j) / (J + 1)
  model$coef_cos <- model$coef_sin <- model$counts <- list()
  for (cls in model$classes) {
    Zc <- Z[y == cls, , drop = FALSE]
    key <- as.character(cls)
    model$counts[[key]] <- nrow(Zc)
    w <- fejer / nrow(Zc)
    model$coef_cos[[key]] <- t(vapply(seq_len(ncol(Zc)), function(d) colSums(cos(pi * outer(Zc[, d], j))) * w, numeric(J)))
    model$coef_sin[[key]] <- t(vapply(seq_len(ncol(Zc)), function(d) colSums(sin(pi * outer(Zc[, d], j))) * w, numeric(J)))
  }
  model
}

fourier_pnn_scale <- function(model, X) {
  Z <- sweep(sweep(as.matrix(X), 2, model$data_min, "-"), 2, model$data_range, "/")
  Z * diff(model$scaling_range) + model$scaling_range[1]
}

#' Fourier PNN class probabilities (models/fpnn.py: FourierPNN.predict_proba, Algorithm 5)
#'
#' f_d = max(1e-10, 0.5 + sum_j A_j cos(pi j z_d) + B_j sin(pi j z_d)) per feature, and
#' P(c | x) proportional to prod_d f_d * N / (2 N_c).
#'
#' @return Matrix of probabilities, one column per class in `model$classes`
#' @export
fourier_pnn_predict_proba <- function(model, X) {
  Z <- fourier_pnn_scale(model, X)
  j <- seq_len(model$J)
  logp <- vapply(model$classes, function(cls) {
    key <- as.character(cls)
    f <- vapply(seq_len(ncol(Z)), function(d) {
      ang <- pi * outer(Z[, d], j)
      0.5 + cos(ang) %*% model$coef_cos[[key]][d, ] + sin(ang) %*% model$coef_sin[[key]][d, ]
    }, numeric(nrow(Z)))
    f <- matrix(pmax(f, 1e-10), nrow = nrow(Z))
    rowSums(log(f)) + log(model$n_samples / (2.0 * model$counts[[key]] + 1e-12))
  }, numeric(nrow(Z)))
  logp <- matrix(logp, nrow = nrow(Z))
  logp <- logp - apply(logp, 1, max)
  p <- exp(logp)
  colnames(p) <- as.character(model$classes)
  p / rowSums(p)
}
