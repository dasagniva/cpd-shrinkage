# ==============================================================================
# 10_dgp.R
# Data Generating Processes (DGPs) for simulation studies.
# Implements Section 3 (Equations 1-2) and Section 5 (Phase P1, Phase P8).
# ==============================================================================

#' Simulate time series with change points, collective anomalies, outliers and AR(1) noise
#'
#' @param n Integer length of the series (default 1000)
#' @param phi Autoregressive parameter |phi| < 1 (default 0.6)
#' @param noise Innovation distribution: "gaussian" or "t3" (default "gaussian")
#' @param sigma Innovation standard deviation (default 1.0)
#' @param k0 Number of sustained structural change points (default 5; 0 for null series)
#' @param jump Jump size multiple of marginal standard deviation (default 2)
#' @param outlier_regime Outlier regime: "none", "isolated", "clustered", "mixed" (default "isolated")
#' @param p_outlier Marginal probability of isolated outliers (default 0.02)
#' @param collective Number of collective anomaly segments (default 2; 0 for none)
#' @param min_coll_len Minimum length of collective anomaly (default 10)
#' @param max_coll_len Maximum length of collective anomaly (default 50)
#' @param seed Optional random seed
#'
#' @return A list with elements:
#'   \item{y}{Observed numeric vector of length n}
#'   \item{truth}{List containing:
#'     \itemize{
#'       \item{mu}{Ground truth piecewise constant mean vector}
#'       \item{psi}{Ground truth point outlier vector}
#'       \item{eps}{Ground truth AR(1) noise vector}
#'       \item{class}{Factor of length n with levels nominal, shift, collective, point}
#'       \item{cp}{Integer vector of true change point indices}
#'       \item{collective}{Data frame with columns start, end of collective segments}
#'       \item{point}{Integer vector of true outlier indices}
#'       \item{phi}{True phi}
#'       \item{sigma}{True innovation sigma}
#'       \item{sigma_marg}{True marginal standard deviation}
#'     }}
#' @export
simulate_series <- function(
    n = 1000,
    phi = 0.6,
    noise = c("gaussian", "t3"),
    sigma = 1.0,
    k0 = 5,
    jump = 2,
    outlier_regime = c("isolated", "none", "clustered", "mixed"),
    p_outlier = 0.02,
    collective = 2,
    min_coll_len = 10,
    max_coll_len = 50,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  stopifnot(n >= 50, abs(phi) < 1, sigma > 0)
  noise <- match.arg(noise)
  outlier_regime <- match.arg(outlier_regime)

  sigma_marg <- sigma / sqrt(1 - phi^2)

  # 1. Generate stationary AR(1) noise: eps_t = phi * eps_{t-1} + u_t
  if (noise == "gaussian") {
    u <- rnorm(n, mean = 0, sd = sigma)
    init_eps <- rnorm(1, mean = 0, sd = sigma_marg)
  } else {
    # Scaled t3: Var(t_3) = 3, so scale by sqrt(1/3) to have Var(u) = sigma^2
    u <- sigma * sqrt(1 / 3) * rt(n, df = 3)
    init_eps <- sigma_marg * sqrt(1 / 3) * rt(1, df = 3)
  }

  eps <- numeric(n)
  eps[1] <- init_eps
  for (t in 2:n) {
    eps[t] <- phi * eps[t - 1] + u[t]
  }

  # 2. Allocation mask to prevent structural overlaps
  # Mask: 0 = free, 1 = shift, 2 = collective, 3 = point
  occupied <- rep(FALSE, n)

  # 3. Generate Sustained Change Points
  cp_locs <- integer(0)
  mu <- rep(0.0, n)

  if (k0 > 0) {
    # Place change points at fixed fractions as in thesis DGP
    fractions <- seq(0.12, 0.88, length.out = k0)
    raw_cps <- round(fractions * n)
    raw_cps <- pmax(5L, pmin(n - 5L, raw_cps))
    cp_locs <- sort(unique(raw_cps))
    k0_actual <- length(cp_locs)

    # Alternate jump directions
    jump_signs <- rep(c(1, -1), length.out = k0_actual)
    jump_mags <- rep(jump, length.out = k0_actual) * sigma_marg * jump_signs

    # Build piecewise constant mu
    current_mu <- 0.0
    last_cp <- 1L
    for (i in seq_along(cp_locs)) {
      cp_i <- cp_locs[i]
      mu[last_cp:(cp_i - 1L)] <- current_mu
      current_mu <- current_mu + jump_mags[i]
      last_cp <- cp_i
    }
    mu[last_cp:n] <- current_mu

    # Reserve change points and buffer (±2) to avoid collision with anomalies
    for (cp_i in cp_locs) {
      idx_range <- max(1L, cp_i - 2L):min(n, cp_i + 2L)
      occupied[idx_range] <- TRUE
    }
  }

  # 4. Generate Collective Anomalies (Shift-and-return excursion segments)
  collective_df <- data.frame(start = integer(0), end = integer(0))
  if (collective > 0) {
    coll_starts <- integer(0)
    coll_ends <- integer(0)

    # Search for available segments
    attempts <- 0
    while (length(coll_starts) < collective && attempts < 200) {
      attempts <- attempts + 1
      L <- sample(min(min_coll_len, max(2L, floor(n / 10L))):min(max_coll_len, max(3L, floor(n / 4L))), 1)
      max_start <- n - L - 3L
      if (max_start < 3L) next
      candidate_start <- sample(3L:max_start, 1)
      candidate_end <- candidate_start + L - 1L
      candidate_range <- max(1L, candidate_start - 2L):min(n, candidate_end + 2L)

      if (!any(occupied[candidate_range])) {
        coll_starts <- c(coll_starts, candidate_start)
        coll_ends <- c(coll_ends, candidate_end)
        occupied[candidate_range] <- TRUE

        # Shift during collective anomaly
        coll_mag <- sample(c(-1, 1), 1) * jump * sigma_marg
        mu[candidate_start:candidate_end] <- mu[candidate_start:candidate_end] + coll_mag
      }
    }

    if (length(coll_starts) > 0) {
      collective_df <- data.frame(
        start = as.integer(coll_starts),
        end = as.integer(coll_ends)
      )
      collective_df <- collective_df[order(collective_df$start), , drop = FALSE]
      rownames(collective_df) <- NULL
    }
  }

  # 5. Generate Outliers (psi)
  psi <- rep(0.0, n)
  point_indices <- integer(0)

  if (outlier_regime != "none") {
    free_indices <- which(!occupied)

    if (outlier_regime %in% c("isolated", "mixed")) {
      p_iso <- if (outlier_regime == "mixed") p_outlier / 2 else p_outlier
      n_iso <- rbinom(1, length(free_indices), p_iso)
      if (n_iso > 0) {
        iso_locs <- sample(free_indices, n_iso)
        signs <- sample(c(-1, 1), n_iso, replace = TRUE)
        mags <- runif(n_iso, 4, 8) * sigma_marg
        psi[iso_locs] <- signs * mags
        point_indices <- c(point_indices, iso_locs)
        occupied[iso_locs] <- TRUE
      }
    }

    if (outlier_regime %in% c("clustered", "mixed")) {
      # 2-state Markov chain with mean patch length 3: P(1|1) = 2/3
      # We insert 1 to 2 outlier clusters of average length 3 into free regions
      n_clusters <- if (outlier_regime == "mixed") 1 else max(1, round(n * p_outlier / 3))
      free_indices <- which(!occupied)

      for (cl in seq_len(n_clusters)) {
        patch_len <- rgeom(1, prob = 1 / 3) + 1L # Geometric with mean 3
        patch_len <- min(patch_len, 8L) # Cap max cluster to prevent becoming collective
        # Find continuous block in free_indices
        avail_starts <- free_indices[free_indices <= (n - patch_len)]
        valid_starts <- avail_starts[sapply(avail_starts, function(st) {
          !any(occupied[st:(st + patch_len - 1L)])
        })]

        if (length(valid_starts) > 0) {
          c_start <- sample(valid_starts, 1)
          c_range <- c_start:(c_start + patch_len - 1L)
          signs <- sample(c(-1, 1), patch_len, replace = TRUE)
          mags <- runif(patch_len, 4, 8) * sigma_marg
          psi[c_range] <- signs * mags
          point_indices <- c(point_indices, c_range)
          occupied[c_range] <- TRUE
        }
      }
    }
  }

  point_indices <- sort(unique(point_indices))

  # 6. Observed series y_t = mu_t + psi_t + eps_t (Eq. 1)
  y <- mu + psi + eps

  # 7. Ground truth labels (mutually exclusive)
  class_vec <- rep("nominal", n)

  if (nrow(collective_df) > 0) {
    for (r in seq_len(nrow(collective_df))) {
      class_vec[collective_df$start[r]:collective_df$end[r]] <- "collective"
    }
  }

  if (length(cp_locs) > 0) {
    class_vec[cp_locs] <- "shift"
  }

  if (length(point_indices) > 0) {
    class_vec[point_indices] <- "point"
  }

  class_factor <- factor(class_vec, levels = c("nominal", "shift", "collective", "point"))

  truth <- list(
    mu          = mu,
    psi         = psi,
    eps         = eps,
    class       = class_factor,
    cp          = cp_locs,
    collective  = collective_df,
    point       = point_indices,
    phi         = phi,
    sigma       = sigma,
    sigma_marg  = sigma_marg
  )

  list(y = y, truth = truth)
}
