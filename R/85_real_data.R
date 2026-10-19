# ==============================================================================
# 85_real_data.R
# Real-Data Application Benchmark (Phase P9).
# Runs proposed model, copula comparator, and all baselines on the four
# empirical benchmark series: Well-log, Oilwell, US Industrial Production,
# and Bitcoin / USD.
# ==============================================================================

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
  source("R/20_frontend_nsp.R")
  source("R/21_frontend_capa.R")
  source("R/22_prior_map.R")
  source("R/32_gibbs_updates.R")
  source("R/40_scoring.R")
  source("R/41_bfdr.R")
  source("R/42_classify.R")
  source("R/31_gibbs_sampler.R")
  source("R/50_copula.R")
  source("R/60_baselines.R")
  source("R/70_metrics.R")
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
  source("../../R/20_frontend_nsp.R")
  source("../../R/21_frontend_capa.R")
  source("../../R/22_prior_map.R")
  source("../../R/32_gibbs_updates.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
  source("../../R/31_gibbs_sampler.R")
  source("../../R/50_copula.R")
  source("../../R/60_baselines.R")
  source("../../R/70_metrics.R")
}

#' Load and Preprocess Empirical Benchmark Series
#'
#' @param name Character name of the dataset: "welllog", "oilwell", "indpro", or "btc"
#' @param raw_dir Directory containing raw CSV files (default "data/raw")
#'
#' @return A list with elements:
#'   \item{name}{Series identifier}
#'   \item{description}{Brief description}
#'   \item{y}{Numeric vector of time series observations}
#'   \item{n}{Integer length of the series}
#'   \item{dates}{Optional date vector if available}
#' @export
load_real_series <- function(name = c("welllog", "oilwell", "indpro", "btc"),
                             raw_dir = "data/raw") {
  name <- match.arg(name)

  if (name == "welllog") {
    csv_path <- file.path(raw_dir, "welllog.csv")
    df <- read.csv(csv_path)
    y <- as.numeric(df$welllog)
    stopifnot(length(y) == 4050, !anyNA(y))
    return(list(
      name        = "welllog",
      description = "Nuclear magnetic resonance well-log (raw, 4,050 obs)",
      y           = y,
      n           = length(y),
      dates       = NULL
    ))
  } else if (name == "oilwell") {
    csv_path <- file.path(raw_dir, "oilwell.csv")
    df <- read.csv(csv_path)
    y <- as.numeric(df$oilwell)
    stopifnot(length(y) == 4000, !anyNA(y))
    return(list(
      name        = "oilwell",
      description = "Oilwell pressure sensor series (4,000 obs)",
      y           = y,
      n           = length(y),
      dates       = NULL
    ))
  } else if (name == "indpro") {
    csv_path <- file.path(raw_dir, "us_indpro.csv")
    df <- read.csv(csv_path)
    # Monthly log growth rate: diff(log(INDPRO))
    indpro_raw <- as.numeric(df$INDPRO)
    y <- diff(log(indpro_raw))
    dates <- as.character(df$observation_date)[-1]
    stopifnot(length(y) == 312, !anyNA(y))
    return(list(
      name        = "indpro",
      description = "US Industrial Production Index monthly log growth (312 obs, 2000-01 to 2026)",
      y           = y,
      n           = length(y),
      dates       = dates
    ))
  } else if (name == "btc") {
    csv_path <- file.path(raw_dir, "btc_usd.csv")
    df <- read.csv(csv_path)
    # Filter 2018-01-01 to 2023-12-31 as in thesis/workflow benchmark (2,191 points)
    df_sub <- subset(df, date >= "2018-01-01" & date <= "2023-12-31")
    y <- as.numeric(df_sub$log_price)
    dates <- as.character(df_sub$date)
    stopifnot(length(y) == 2191, !anyNA(y))
    return(list(
      name        = "btc",
      description = "Bitcoin / USD daily log price (2,191 obs, 2018-01-01 to 2023-12-31)",
      y           = y,
      n           = length(y),
      dates       = dates
    ))
  }
}

#' Run a Method on Empirical Series
#'
#' @param y Numeric vector of observations
#' @param method Character: "proposed", "copula", or baseline name
#' @param m_args Optional list of arguments for method
#'
#' @return A validated `cpd_fit` object
#' @export
run_real_data_method <- function(y, method = "proposed", m_args = list()) {
  n <- length(y)
  stopifnot(n >= 10, is.numeric(y))

  t0 <- proc.time()[3]

  fit <- tryCatch({
    if (method == "proposed") {
      alpha_nsp <- if (!is.null(m_args$alpha_nsp)) m_args$alpha_nsp else 0.05
      kappa_in  <- if (!is.null(m_args$kappa_in)) m_args$kappa_in else 2.0
      delta     <- if (!is.null(m_args$delta)) m_args$delta else 1.0
      n_iter    <- if (!is.null(m_args$n_iter)) m_args$n_iter else 2500
      burn_in   <- if (!is.null(m_args$burn_in)) m_args$burn_in else 1000
      n_chains  <- if (!is.null(m_args$n_chains)) m_args$n_chains else 2
      alpha_fdr <- if (!is.null(m_args$alpha_fdr)) m_args$alpha_fdr else 0.05
      dynamic   <- if (!is.null(m_args$dynamic)) m_args$dynamic else TRUE

      # Front-ends
      intervals  <- run_frontend_nsp(y, alpha = alpha_nsp)
      capa_res   <- run_frontend_capa(y, prewhiten = TRUE)
      prior_maps <- build_prior_maps(n = n, nsp_intervals = intervals$intervals,
                                     capa_results = capa_res, kappa_in = kappa_in)

      gibbs_sampler(
        y               = y,
        kappa           = prior_maps$kappa,
        a_vec           = prior_maps$a,
        n_iter          = n_iter,
        burn_in         = burn_in,
        thin            = 1,
        n_chains        = n_chains,
        dynamic         = dynamic,
        delta           = delta,
        capa_collective = capa_res$collective,
        alpha_fdr       = alpha_fdr
      )
    } else if (method == "copula") {
      b_boot <- if (!is.null(m_args$B_boot)) m_args$B_boot else 50
      fit_copula_comparator(y, alpha = 0.05, B_boot = b_boot)
    } else {
      do.call(fit_baseline, c(list(method = method, y = y), m_args))
    }
  }, error = function(e) {
    warning(sprintf("Error running %s on real series: %s", method, e$message))
    NULL
  })

  if (is.null(fit)) {
    return(new_cpd_fit(
      method      = method,
      runtime_sec = as.numeric(proc.time()[3] - t0),
      diag        = list(failed = TRUE)
    ))
  }

  fit
}

#' Summarize Detection Results for a Real-Data Fit
#'
#' @param fit A `cpd_fit` object
#' @param series_name Character name of series
#'
#' @return A single-row data frame
#' @export
summarize_real_data_fit <- function(fit, series_name = "") {
  stopifnot(inherits(fit, "cpd_fit"))

  n_cp <- length(fit$cp)
  n_pt <- length(fit$point)
  n_coll_segs <- if (is.data.frame(fit$collective)) nrow(fit$collective) else 0L
  n_coll_pts  <- if (n_coll_segs > 0) sum(fit$collective$end - fit$collective$start + 1L) else 0L

  phi_est <- if (!is.null(fit$diag$phi_hat)) fit$diag$phi_hat else NA_real_
  sigma_est <- if (!is.null(fit$diag$sigma_hat)) fit$diag$sigma_hat else NA_real_

  data.frame(
    series       = series_name,
    method       = fit$method,
    n_cp         = n_cp,
    n_point      = n_pt,
    n_coll_segs  = n_coll_segs,
    n_coll_pts   = n_coll_pts,
    phi_est      = phi_est,
    sigma_est    = sigma_est,
    runtime_sec  = fit$runtime_sec,
    stringsAsFactors = FALSE
  )
}

#' Run Comprehensive Real-Data Benchmark Across Datasets and Methods
#'
#' @param series_names Vector of series names to process
#' @param methods Vector of methods to evaluate
#' @param m_args_list Named list of arguments per method
#' @param raw_dir Path to raw data directory
#'
#' @return A list containing:
#'   \item{fits}{Nested list of all fits: `fits[[series]][[method]]`}
#'   \item{summary_table}{Consolidated detection and timing table}
#'   \item{timestamp}{Execution timestamp}
#' @export
run_real_data_benchmark <- function(
    series_names = c("welllog", "oilwell", "indpro", "btc"),
    methods      = c("proposed", "copula", "decafs", "pelt", "robseg", "capa", "fpnn", "ch4"),
    m_args_list  = list(),
    raw_dir      = "data/raw"
) {
  cat("\n====================================================================\n")
  cat("     EMPIRICAL REAL-DATA BENCHMARK (Section 5, Phase P9)\n")
  cat("====================================================================\n\n")

  fits <- list()
  summary_rows <- list()
  row_idx <- 1L

  for (s_name in series_names) {
    cat(sprintf("\n--- Loading Dataset: %s ---\n", s_name))
    dat <- load_real_series(s_name, raw_dir = raw_dir)
    cat(sprintf("  %s: %d observations\n", dat$description, dat$n))
    fits[[s_name]] <- list()

    for (m in methods) {
      cat(sprintf("  Fitting method '%s' on %s (n = %d)... ", m, s_name, dat$n))
      m_args <- if (!is.null(m_args_list[[m]])) m_args_list[[m]] else list()
      t0 <- proc.time()[3]
      fit_res <- run_real_data_method(dat$y, method = m, m_args = m_args)
      t_elapsed <- as.numeric(proc.time()[3] - t0)

      fits[[s_name]][[m]] <- fit_res

      sum_row <- summarize_real_data_fit(fit_res, series_name = s_name)
      summary_rows[[row_idx]] <- sum_row
      row_idx <- row_idx + 1L

      cat(sprintf("done (%.2fs) | CPs: %d, Points: %d, Coll: %d\n",
                  t_elapsed, sum_row$n_cp, sum_row$n_point, sum_row$n_coll_segs))
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  cat("\n====================================================================\n")
  cat("     REAL-DATA BENCHMARK COMPLETE\n")
  cat("====================================================================\n\n")
  print(summary_df)

  list(
    fits          = fits,
    summary_table = summary_df,
    timestamp     = Sys.time()
  )
}
