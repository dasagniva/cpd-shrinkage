# ==============================================================================
# 42_classify.R
# Four-Class Labelling Rule: Collective Anomalies, Sustained Shifts,
# Point Anomalies, and Nominal Observations.
# Implements Section 3.6.
# ==============================================================================

#' Four-Class Labelling Rule (Section 3.6)
#'
#' Applies the four-class labelling rules in strict sequential priority
#' ("the first match wins"):
#'
#' 1. Collective anomaly:
#'    A pair of flagged increments at t1 < t2 with opposite signs, where
#'    2 <= t2 - t1 <= L_max and |mu_hat_{t2+1} - mu_hat_{t1-1}| <= r * |omega_hat_{t1}|
#'    (default r = 0.3).
#'    Alternatively, any pair with 2 <= t2 - t1 <= L_max and opposite signs having
#'    at least 50% overlap with a CAPA collective segment supported by flagged increment(s).
#'    Points t1, ..., t2 are labelled "collective".
#'
#' 2. Sustained shift:
#'    A flagged increment that is not paired under rule 1. Point t is labelled "shift".
#'
#' 3. Point anomaly:
#'    A flagged psi_t outside any collective segment. Point t is labelled "point".
#'
#' 4. Nominal:
#'    Everything else. Labelled "nominal".
#'
#' @param n Integer length of the series.
#' @param flagged_omega Integer vector of flagged increment locations in time (t in 2..n).
#' @param flagged_psi Integer vector of flagged outlier locations (t in 1..n).
#' @param mu_hat Numeric vector of posterior mean level of length n.
#' @param omega_hat Optional numeric vector of increment magnitudes. If NULL, computed as diff(mu_hat).
#' @param capa_collective Optional data frame with columns `start` and `end` from CAPA.
#' @param L_max Maximum length of a collective anomaly segment (default 50, the longest
#'   collective anomaly in the simulation design; spec change C4 from 0.1 n).
#' @param r Tolerance fraction for return to baseline (default 0.3).
#' @param capa_overlap_min Minimum overlap fraction with CAPA collective segment (default 0.5).
#' @param cp_window Half-width h used for windowed shift probabilities. Jump sizes and
#'   return levels are then measured across the window: mu_hat[t + h] - mu_hat[t - h - 1].
#'   With h = 0 this is exactly omega_hat_t = mu_hat_t - mu_hat_{t-1}.
#'
#' @return A list with elements:
#'   \item{class}{Factor of length n with levels c("nominal", "shift", "collective", "point")}
#'   \item{cp}{Integer vector of sustained shift change points}
#'   \item{collective}{Data frame with columns `start` and `end` of collective anomaly segments}
#'   \item{point}{Integer vector of point anomaly locations}
#' @export
classify_four_class <- function(
    n,
    flagged_omega,
    flagged_psi,
    mu_hat,
    omega_hat = NULL,
    capa_collective = NULL,
    L_max = NULL,
    r = 0.3,
    capa_overlap_min = 0.5,
    cp_window = 0L
) {
  stopifnot(is.numeric(n), n >= 2, length(mu_hat) == n)
  h <- as.integer(cp_window)

  if (is.null(L_max)) {
    L_max <- 50L
  }
  stopifnot(is.numeric(L_max), L_max >= 2)

  # Filter flagged_omega to valid range [2, n]
  flagged_w <- sort(unique(as.integer(flagged_omega)))
  flagged_w <- flagged_w[flagged_w >= 2L & flagged_w <= n]

  # Filter flagged_psi to valid range [1, n]
  flagged_p <- sort(unique(as.integer(flagged_psi)))
  flagged_p <- flagged_p[flagged_p >= 1L & flagged_p <= n]

  # Compute increment magnitudes if not provided
  if (is.null(omega_hat)) {
    # omega_hat[t] corresponds to mu_hat[t + h] - mu_hat[t - h - 1] (h = 0: mu_hat[t] - mu_hat[t-1])
    t_idx <- 2:n
    omega_hat_vec <- numeric(n)
    omega_hat_vec[t_idx] <- mu_hat[pmin(n, t_idx + h)] - mu_hat[pmax(1L, t_idx - h - 1L)]
  } else {
    stopifnot(length(omega_hat) %in% c(n - 1L, n))
    omega_hat_vec <- if (length(omega_hat) == n - 1L) c(0.0, omega_hat) else omega_hat
  }

  # Normalize CAPA collective segments
  has_capa <- !is.null(capa_collective) && is.data.frame(capa_collective) && nrow(capa_collective) > 0

  # Helper: compute overlap between [t1, t2] and CAPA segments
  calc_capa_overlap <- function(t1, t2) {
    if (!has_capa) return(0.0)
    max_ov <- 0.0
    for (m in seq_len(nrow(capa_collective))) {
      s_c <- capa_collective$start[m]
      e_c <- capa_collective$end[m]
      len_c <- e_c - s_c + 1L
      len_cand <- t2 - t1 + 1L
      inter <- max(0L, min(t2, e_c) - max(t1, s_c) + 1L)
      if (inter > 0) {
        # Overlap fraction relative to smaller of the two segments
        ov <- inter / min(len_cand, len_c)
        if (ov > max_ov) max_ov <- ov
      }
    }
    max_ov
  }

  # Initialize classification vector
  class_vec <- rep("nominal", n)

  # ----------------------------------------------------------------------------
  # Rule 1: Collective Anomalies
  # ----------------------------------------------------------------------------
  collective_starts <- integer(0)
  collective_ends   <- integer(0)
  paired_omega      <- integer(0)

  if (length(flagged_w) >= 2) {
    available_w <- flagged_w

    i <- 1L
    while (i <= length(available_w)) {
      t1 <- available_w[i]
      w1 <- omega_hat_vec[t1]
      sign1 <- sign(w1)

      if (sign1 != 0) {
        candidates <- available_w[available_w > t1]
        best_t2 <- NA_integer_
        best_diff <- Inf

        for (t2 in candidates) {
          sep <- t2 - t1
          if (sep < 2L) next
          if (sep > L_max) break # since available_w is sorted

          w2 <- omega_hat_vec[t2]
          sign2 <- sign(w2)

          # Opposite signs requirement
          if (sign1 * sign2 < 0) {
            # Baseline before jump: mu_{t1 - h - 1}
            b_level <- mu_hat[max(1L, t1 - h - 1L)]
            # Level after return: mu_{t2 + h} or mu_{t2 + h + 1}
            a_level1 <- mu_hat[min(n, t2 + h + 1L)]
            a_level2 <- mu_hat[min(n, t2 + h)]

            ret_diff <- min(abs(a_level1 - b_level), abs(a_level2 - b_level))
            tol <- r * abs(w1)

            # Criteria 1: Return condition satisfied
            cond_return <- ret_diff <= tol

            # Criteria 2: Overlap with CAPA collective segment >= capa_overlap_min
            cond_capa <- if (has_capa) {
              calc_capa_overlap(t1, t2) >= capa_overlap_min
            } else {
              FALSE
            }

            if (cond_return || cond_capa) {
              if (ret_diff < best_diff) {
                best_diff <- ret_diff
                best_t2 <- t2
              }
            }
          }
        }

        if (!is.na(best_t2)) {
          # omega_{t2} = mu_{t2} - mu_{t2-1} is the return step, so the excursion
          # occupies t1, ..., t2 - 1 (t2 - 1 >= t1 + 1 because t2 - t1 >= 2)
          collective_starts <- c(collective_starts, t1)
          collective_ends   <- c(collective_ends, best_t2 - 1L)
          paired_omega      <- c(paired_omega, t1, best_t2)

          class_vec[t1:(best_t2 - 1L)] <- "collective"

          # Remove paired points from available
          available_w <- setdiff(available_w, c(t1, best_t2))
          # Restart or continue from next available
          i <- 1L
          next
        }
      }
      i <- i + 1L
    }
  }

  # Build collective data frame; interleaved pairs can give overlapping segments,
  # which are merged into one collective anomaly
  if (length(collective_starts) > 0) {
    collective_df <- data.frame(
      start = as.integer(collective_starts),
      end   = as.integer(collective_ends)
    )
    collective_df <- collective_df[order(collective_df$start), , drop = FALSE]
    merged <- collective_df[1, , drop = FALSE]
    for (k in seq_len(nrow(collective_df))[-1]) {
      last <- nrow(merged)
      if (collective_df$start[k] <= merged$end[last]) {
        merged$end[last] <- max(merged$end[last], collective_df$end[k])
      } else {
        merged <- rbind(merged, collective_df[k, , drop = FALSE])
      }
    }
    collective_df <- merged
    rownames(collective_df) <- NULL
  } else {
    collective_df <- data.frame(start = integer(0), end = integer(0))
  }

  # ----------------------------------------------------------------------------
  # Rule 2: Sustained Shifts (Flagged increments not paired under Rule 1)
  # ----------------------------------------------------------------------------
  unpaired_w <- setdiff(flagged_w, paired_omega)
  shift_locs <- integer(0)

  if (length(unpaired_w) > 0) {
    for (sw in unpaired_w) {
      # Only label as shift if not already consumed by collective
      if (class_vec[sw] != "collective") {
        class_vec[sw] <- "shift"
        shift_locs <- c(shift_locs, sw)
      }
    }
    shift_locs <- sort(shift_locs)
  }

  # ----------------------------------------------------------------------------
  # Rule 3: Point Anomalies (Flagged psi outside any collective segment)
  # ----------------------------------------------------------------------------
  point_locs <- integer(0)

  if (length(flagged_p) > 0) {
    for (pt in flagged_p) {
      # Only label as point if outside any collective segment and not a shift
      if (class_vec[pt] == "nominal") {
        class_vec[pt] <- "point"
        point_locs <- c(point_locs, pt)
      }
    }
    point_locs <- sort(point_locs)
  }

  # ----------------------------------------------------------------------------
  # Rule 4: Nominal (Remaining points remain "nominal")
  # ----------------------------------------------------------------------------
  class_factor <- factor(class_vec, levels = CPD_CLASS_LEVELS)

  list(
    class      = class_factor,
    cp         = as.integer(shift_locs),
    collective = collective_df,
    point      = as.integer(point_locs)
  )
}

#' Classify and Construct cpd_fit Object from Posterior Scores
#'
#' Unified end-to-end wrapper: applies Bayesian FDR to p_omega and p_psi,
#' runs the four-class labelling rule, and packages the result into a validated cpd_fit.
#'
#' @param y Numeric vector of observed series (length n).
#' @param mu_hat Numeric vector of estimated posterior mean level (length n).
#' @param s_out Numeric vector of outlier scores o_t (length n).
#' @param p_omega Numeric vector of increment probabilities (length n - 1 or length n).
#' @param p_psi Numeric vector of outlier probabilities (length n).
#' @param alpha_fdr Target FDR level (default 0.05).
#' @param capa_collective Optional data frame of CAPA collective segments.
#' @param method Identifier string (default "proposed").
#' @param runtime_sec Execution time in seconds.
#' @param diag Arbitrary diagnostic list.
#' @param L_max Maximum length for collective anomalies.
#' @param r Return tolerance fraction (default 0.3).
#' @param cp_window Half-width h of windowed shift probabilities (0 = Eq. 9 per time point).
#'   When h > 0, each run of BFDR-flagged times (gaps of at most h) is collapsed to the
#'   single time with the largest `omega_strength`.
#' @param omega_strength Optional length n - 1 vector (e.g. posterior mean |omega_t|)
#'   used to pick one change location per flagged run.
#'
#' @return A validated `cpd_fit` object.
#' @export
classify_and_build_fit <- function(
    y,
    mu_hat,
    s_out,
    p_omega,
    p_psi,
    alpha_fdr = 0.05,
    capa_collective = NULL,
    method = "proposed",
    runtime_sec = 0.0,
    diag = list(),
    L_max = NULL,
    r = 0.3,
    cp_window = 0L,
    omega_strength = NULL
) {
  n <- length(y)
  stopifnot(n >= 2, length(mu_hat) == n, length(s_out) == n, length(p_psi) == n)

  # Check p_omega dimension
  if (length(p_omega) == n - 1L) {
    # Increment j in 1..(n-1) maps to time point j + 1
    bfdr_w <- bayesian_fdr(p_omega, alpha = alpha_fdr)
    flagged_omega <- if (bfdr_w$k > 0) bfdr_w$flagged + 1L else integer(0)
    s_cp <- c(0.0, p_omega)
  } else if (length(p_omega) == n) {
    bfdr_w <- bayesian_fdr(p_omega[-1], alpha = alpha_fdr)
    flagged_omega <- if (bfdr_w$k > 0) bfdr_w$flagged + 1L else integer(0)
    s_cp <- p_omega
  } else {
    stop("p_omega must have length n - 1 or n")
  }

  # Windowed probabilities flag a run of ~2h+1 times around each change: keep one
  # representative per run. Flagged times at most h apart belong to the same run,
  # since a window of half-width h cannot separate changes that close.
  if (cp_window > 0 && length(flagged_omega) > 1) {
    strength <- if (!is.null(omega_strength)) c(0.0, omega_strength) else s_cp
    runs <- split(flagged_omega, cumsum(c(1L, diff(flagged_omega) > cp_window)))
    flagged_omega <- vapply(runs, function(rn) rn[which.max(strength[rn])], integer(1))
    names(flagged_omega) <- NULL
  }

  bfdr_psi <- bayesian_fdr(p_psi, alpha = alpha_fdr)
  flagged_psi <- bfdr_psi$flagged

  # Run four-class labelling
  clf <- classify_four_class(
    n                = n,
    flagged_omega    = flagged_omega,
    flagged_psi      = flagged_psi,
    mu_hat           = mu_hat,
    capa_collective  = capa_collective,
    L_max            = L_max,
    r                = r,
    cp_window        = cp_window
  )

  scores_df <- data.frame(
    t     = seq_len(n),
    s_cp  = s_cp,
    s_out = s_out
  )

  # Attach FDR details to diag
  diag$bfdr_omega <- bfdr_w
  diag$bfdr_psi   <- bfdr_psi

  new_cpd_fit(
    method      = method,
    cp          = clf$cp,
    cp_int      = NULL,
    point       = clf$point,
    collective  = clf$collective,
    class       = clf$class,
    scores      = scores_df,
    runtime_sec = runtime_sec,
    diag        = diag
  )
}
