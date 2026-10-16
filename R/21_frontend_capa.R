# ==============================================================================
# 21_frontend_capa.R
# Anomaly Front-End: Prewhitened CAPA wrapper for point and collective anomalies.
# Implements Section 3.3 and Section 5 (Phase P2).
# ==============================================================================

#' Estimate AR(1) phi coefficient robustly in the presence of jumps and outliers
#'
#' Removes local level using a running median, then fits a robust M-estimator
#' regression of lag-1 residuals.
#'
#' @param y Numeric vector of time series observations
#' @param max_phi Maximum absolute value for phi (default 0.95)
#'
#' @return Numeric scalar estimate phi_hat in [-max_phi, max_phi]
#' @export
estimate_phi_robust <- function(y, max_phi = 0.95) {
  n <- length(y)
  stopifnot(n >= 10)

  # Running median bandwidth (odd integer proportional to n)
  k_med <- min(51L, 2L * as.integer(floor(n / 25)) + 1L)
  if (k_med < 11L) k_med <- 11L
  if (k_med %% 2L == 0L) k_med <- k_med + 1L

  # Detrend locally
  r <- y - stats::runmed(y, k = k_med)

  r1 <- r[-n]
  r2 <- r[-1]

  # Robust linear regression without intercept
  fit <- tryCatch({
    MASS::rlm(r2 ~ r1 - 1, psi = MASS::psi.bisquare)
  }, error = function(e) {
    MASS::rlm(r2 ~ r1 - 1, psi = MASS::psi.huber)
  })

  phi_est <- as.numeric(coef(fit)[1])
  if (is.na(phi_est)) phi_est <- 0.0

  # Truncate to keep away from non-stationary boundary
  pmax(-max_phi, pmin(max_phi, phi_est))
}

#' Run CAPA with optional prewhitening
#'
#' @param y Numeric vector of time series observations
#' @param prewhiten Logical indicating whether to prewhiten before CAPA (default TRUE)
#' @param type Type of change in CAPA: "mean" or "meanvar" (default "mean")
#' @param min_seg_len Minimum segment length for collective anomalies (default 10)
#' @param ... Additional arguments passed to `anomaly::capa`
#'
#' @return A list with elements:
#'   \item{point}{Integer vector of point anomaly indices}
#'   \item{collective}{Data frame with columns `start` and `end` of collective anomalies}
#'   \item{phi_hat}{Estimated AR(1) persistence (0 if prewhiten = FALSE)}
#'   \item{y_tilde}{Prewhitened / scaled series fed to CAPA}
#'   \item{raw}{Raw capa.class object}
#' @export
run_frontend_capa <- function(
    y,
    prewhiten = TRUE,
    type = "mean",
    min_seg_len = 10,
    ...
) {
  n <- length(y)
  stopifnot(is.numeric(y), n >= min_seg_len)

  med_y <- median(y, na.rm = TRUE)

  if (isTRUE(prewhiten)) {
    phi_hat <- estimate_phi_robust(y)

    # Prewhitened series P_phi * (y - med_y)
    y_centered <- y - med_y
    y_tilde <- numeric(n)
    y_tilde[1] <- sqrt(1 - phi_hat^2) * y_centered[1]
    y_tilde[2:n] <- y_centered[2:n] - phi_hat * y_centered[1:(n - 1L)]

    # Robust scaling to standard normal scale
    scale_factor <- stats::mad(y_tilde)
    if (scale_factor <= 1e-8) scale_factor <- 1.0
    y_tilde <- y_tilde / scale_factor
  } else {
    phi_hat <- 0.0
    scale_factor <- stats::mad(y - med_y)
    if (scale_factor <= 1e-8) scale_factor <- 1.0
    y_tilde <- (y - med_y) / scale_factor
  }

  res <- anomaly::capa(
    x = y_tilde,
    type = type,
    min_seg_len = min_seg_len,
    ...
  )

  # Extract point anomalies
  pts <- anomaly::point_anomalies(res)
  if (is.null(pts) || nrow(pts) == 0) {
    point_locs <- integer(0)
  } else {
    point_locs <- sort(as.integer(pts$location))
  }

  # Extract collective anomalies
  colls <- anomaly::collective_anomalies(res)
  if (is.null(colls) || nrow(colls) == 0) {
    collective_df <- data.frame(start = integer(0), end = integer(0))
  } else {
    collective_df <- data.frame(
      start = as.integer(colls$start),
      end   = as.integer(colls$end)
    )
    collective_df <- collective_df[order(collective_df$start), , drop = FALSE]
    rownames(collective_df) <- NULL
  }

  list(
    point      = point_locs,
    collective = collective_df,
    phi_hat    = phi_hat,
    y_tilde    = y_tilde,
    raw        = res
  )
}
