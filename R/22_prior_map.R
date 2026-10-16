# ==============================================================================
# 22_prior_map.R
# Maps front-end outputs (NSP intervals, CAPA anomalies) to prior hyperparameters.
# Implements Section 3.2 (Eq. 5) and Section 3.3 (Eq. 7).
# ==============================================================================

#' Build significance-informed prior modifiers kappa_t and a_t
#'
#' @param n Integer length of the time series
#' @param nsp_intervals Data frame with columns `start` and `end` from NSP
#' @param capa_results List returned by `run_frontend_capa`, containing `point` and `collective`
#' @param kappa_in Local shrinkage inflation factor on candidate set I (default 2.0; sensitivity in {1, 2, 5})
#' @param endpoint_buffer Buffer around CAPA collective endpoints (default 2L, as in Section 3.2)
#'
#' @return A list with elements:
#'   \item{kappa}{Numeric vector of length n for local Horseshoe+ scales (Eq. 5)}
#'   \item{a}{Numeric binary vector of length n for dynamic outlier offset (Eq. 7)}
#'   \item{I_mask}{Logical vector of length n indicating membership in candidate set I}
#'   \item{I_indices}{Integer vector of indices in candidate set I}
#' @export
build_prior_maps <- function(
    n,
    nsp_intervals,
    capa_results,
    kappa_in = 2.0,
    endpoint_buffer = 2L
) {
  stopifnot(n >= 1, kappa_in >= 1.0, endpoint_buffer >= 0)

  # 1. Construct candidate set I (Section 3.2)
  # Union of NSP intervals at level alpha_NSP + endpoints (±2 points) of CAPA collective segments
  in_I <- rep(FALSE, n)

  # Add NSP intervals
  if (is.data.frame(nsp_intervals) && nrow(nsp_intervals) > 0) {
    for (r in seq_len(nrow(nsp_intervals))) {
      s <- max(1L, as.integer(nsp_intervals$start[r]))
      e <- min(n, as.integer(nsp_intervals$end[r]))
      if (s <= e) {
        in_I[s:e] <- TRUE
      }
    }
  }

  # Add endpoints (± endpoint_buffer) of CAPA collective segments
  if (is.list(capa_results) && is.data.frame(capa_results$collective) && nrow(capa_results$collective) > 0) {
    colls <- capa_results$collective
    for (r in seq_len(nrow(colls))) {
      c_start <- as.integer(colls$start[r])
      c_end   <- as.integer(colls$end[r])

      # Start endpoint buffer
      s_range <- max(1L, c_start - endpoint_buffer):min(n, c_start + endpoint_buffer)
      in_I[s_range] <- TRUE

      # End endpoint buffer
      e_range <- max(1L, c_end - endpoint_buffer):min(n, c_end + endpoint_buffer)
      in_I[e_range] <- TRUE
    }
  }

  # 2. Compute kappa_t (Eq. 5)
  # kappa_t = kappa_in if t in I, else 1
  kappa <- ifelse(in_I, as.numeric(kappa_in), 1.0)

  # 3. Compute a_t (Eq. 7)
  # a_t = 1 if CAPA flags t as point anomaly, else 0
  a <- rep(0.0, n)
  if (is.list(capa_results) && length(capa_results$point) > 0) {
    valid_pts <- capa_results$point[capa_results$point >= 1 & capa_results$point <= n]
    if (length(valid_pts) > 0) {
      a[valid_pts] <- 1.0
    }
  }

  list(
    kappa     = kappa,
    a         = a,
    I_mask    = in_I,
    I_indices = which(in_I)
  )
}
