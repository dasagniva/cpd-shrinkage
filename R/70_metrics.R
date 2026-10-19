# ==============================================================================
# 70_metrics.R
# Comprehensive Evaluation Metrics for Change Points, Anomalies, and Classification.
# Implements Section 6 of CPD2_Code_Workflow.md.
# ==============================================================================

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
}

# ------------------------------------------------------------------------------
# 1. Change Point Evaluation (Shift Class)
# ------------------------------------------------------------------------------
#' Evaluate Change Point Detection Metrics
#'
#' Computes:
#' - Precision, Recall, and F1 within a margin of +-5 points (scaled for n)
#' - Hausdorff distance between true and estimated change sets
#' - Partition covering metric (van den Burg & Williams, 2020)
#' - Absolute count error |K_hat - K|
#'
#' @param truth_cp Integer vector of true change point locations
#' @param est_cp Integer vector of estimated change point locations
#' @param n Integer length of the series
#' @param base_margin Base tolerance margin (default 5 for n = 1000)
#'
#' @return A list of change point evaluation metrics
#' @export
evaluate_change_points <- function(truth_cp, est_cp, n, base_margin = 5) {
  truth_cp <- sort(unique(as.integer(truth_cp)))
  est_cp   <- sort(unique(as.integer(est_cp)))

  # Scale margin with n: margin = max(2, round(base_margin * n / 1000))
  margin <- max(2L, as.integer(round(base_margin * n / 1000.0)))

  K_true <- length(truth_cp)
  K_est  <- length(est_cp)
  count_error <- abs(K_est - K_true)

  if (K_true == 0 && K_est == 0) {
    return(list(
      precision   = 1.0,
      recall      = 1.0,
      f1          = 1.0,
      hausdorff   = 0.0,
      covering    = 1.0,
      count_error = 0L,
      margin      = margin
    ))
  }

  if (K_true == 0 && K_est > 0) {
    return(list(
      precision   = 0.0,
      recall      = 1.0,
      f1          = 0.0,
      hausdorff   = as.numeric(n),
      covering    = 0.0,
      count_error = as.integer(K_est),
      margin      = margin
    ))
  }

  if (K_true > 0 && K_est == 0) {
    return(list(
      precision   = 1.0,
      recall      = 0.0,
      f1          = 0.0,
      hausdorff   = as.numeric(n),
      covering    = 0.0,
      count_error = as.integer(K_true),
      margin      = margin
    ))
  }

  # Precision: fraction of estimated CPs matched to a true CP within margin
  matched_est <- sapply(est_cp, function(e) any(abs(e - truth_cp) <= margin))
  precision <- mean(matched_est)

  # Recall: fraction of true CPs matched to an estimated CP within margin
  matched_true <- sapply(truth_cp, function(t) any(abs(t - est_cp) <= margin))
  recall <- mean(matched_true)

  f1 <- if (precision + recall > 0) 2.0 * (precision * recall) / (precision + recall) else 0.0

  # Hausdorff distance: max( max_t min_e |t - e|, max_e min_t |e - t| )
  dist_t_to_e <- sapply(truth_cp, function(t) min(abs(t - est_cp)))
  dist_e_to_t <- sapply(est_cp, function(e) min(abs(e - truth_cp)))
  hausdorff <- max(max(dist_t_to_e), max(dist_e_to_t))

  # Partition Covering Metric (van den Burg & Williams, 2020)
  covering <- compute_partition_covering(truth_cp, est_cp, n)

  list(
    precision   = as.numeric(precision),
    recall      = as.numeric(recall),
    f1          = as.numeric(f1),
    hausdorff   = as.numeric(hausdorff),
    covering    = as.numeric(covering),
    count_error = as.integer(count_error),
    margin      = margin
  )
}

#' Compute Partition Covering Metric (van den Burg & Williams, 2020)
#'
#' @param cp1 Change points defining partition 1
#' @param cp2 Change points defining partition 2
#' @param n Series length
#'
#' @return Symmetric covering score in [0, 1]
#' @export
compute_partition_covering <- function(cp1, cp2, n) {
  # Build segments from change points: [s_i, e_i]
  make_segments <- function(cps) {
    bounds <- c(1L, sort(unique(cps)), n + 1L)
    bounds <- unique(bounds)
    seg_list <- vector("list", length(bounds) - 1L)
    for (i in seq_len(length(bounds) - 1L)) {
      seg_list[[i]] <- bounds[i]:(bounds[i + 1L] - 1L)
    }
    seg_list
  }

  segs1 <- make_segments(cp1)
  segs2 <- make_segments(cp2)

  # One-sided covering C(P1, P2) = (1/n) * sum_{A in P1} |A| * max_{A' in P2} Jaccard(A, A')
  cover_one_sided <- function(P1, P2) {
    tot <- 0.0
    for (A in P1) {
      len_A <- length(A)
      max_jacc <- 0.0
      for (B in P2) {
        inter <- length(intersect(A, B))
        if (inter > 0) {
          union_len <- length(union(A, B))
          jacc <- inter / union_len
          if (jacc > max_jacc) max_jacc <- jacc
        }
      }
      tot <- tot + len_A * max_jacc
    }
    tot / n
  }

  c12 <- cover_one_sided(segs1, segs2)
  c21 <- cover_one_sided(segs2, segs1)
  min(c12, c21)
}

# ------------------------------------------------------------------------------
# 2. NSP Intervals Evaluation
# ------------------------------------------------------------------------------
#' Evaluate NSP Intervals
#'
#' @param truth_cp Integer vector of true change points
#' @param intervals Data frame with `start` and `end` columns
#' @param alpha_nsp Nominal significance level (default 0.05)
#'
#' @return A list with coverage, empirical nominal gap, and median width
#' @export
evaluate_nsp_intervals <- function(truth_cp, intervals, alpha_nsp = 0.05) {
  if (is.null(intervals) || nrow(intervals) == 0) {
    return(list(
      n_intervals  = 0L,
      coverage     = NA_real_,
      median_width = NA_real_
    ))
  }

  stopifnot(c("start", "end") %in% names(intervals))
  n_ints <- nrow(intervals)

  # Check if each interval contains at least one true change point
  covered <- sapply(seq_len(n_ints), function(i) {
    st <- intervals$start[i]
    en <- intervals$end[i]
    any(truth_cp >= st & truth_cp <= en)
  })

  coverage <- mean(covered)
  widths <- intervals$end - intervals$start + 1L
  med_width <- as.numeric(median(widths))

  list(
    n_intervals  = as.integer(n_ints),
    coverage     = as.numeric(coverage),
    target_cov   = 1.0 - alpha_nsp,
    median_width = med_width
  )
}

# ------------------------------------------------------------------------------
# 3. Point Anomalies Evaluation
# ------------------------------------------------------------------------------
#' Evaluate Point Anomaly Detection Metrics
#'
#' @param truth_pt Integer vector of true point anomaly locations
#' @param est_pt Integer vector of estimated point anomaly locations
#' @param o_t Optional numeric vector of outlier scores of length n
#' @param n Series length (required if o_t is provided)
#'
#' @return A list with precision, recall, f1, realised FDR, and optional PR-AUC
#' @export
evaluate_point_anomalies <- function(truth_pt, est_pt, o_t = NULL, n = NULL) {
  truth_pt <- sort(unique(as.integer(truth_pt)))
  est_pt   <- sort(unique(as.integer(est_pt)))

  K_true <- length(truth_pt)
  K_est  <- length(est_pt)

  if (K_true == 0 && K_est == 0) {
    prec <- 1.0
    rec  <- 1.0
    f1   <- 1.0
    fdr  <- 0.0
  } else if (K_true == 0 && K_est > 0) {
    prec <- 0.0
    rec  <- 1.0
    f1   <- 0.0
    fdr  <- 1.0
  } else if (K_true > 0 && K_est == 0) {
    prec <- 1.0
    rec  <- 0.0
    f1   <- 0.0
    fdr  <- 0.0
  } else {
    n_tp <- sum(est_pt %in% truth_pt)
    prec <- n_tp / K_est
    rec  <- n_tp / K_true
    f1   <- if (prec + rec > 0) 2.0 * (prec * rec) / (prec + rec) else 0.0
    fdr  <- (K_est - n_tp) / K_est
  }

  pr_auc <- NA_real_
  if (!is.null(o_t) && !is.null(n) && K_true > 0 && length(o_t) == n) {
    # Compute threshold-free PR-AUC
    # Binary label vector: 1 for outlier, 0 for clean
    y_true <- integer(n)
    y_true[truth_pt] <- 1L

    # Order scores decreasingly
    ord <- order(o_t, decreasing = TRUE)
    y_sorted <- y_true[ord]
    cum_tp <- cumsum(y_sorted)
    cum_fp <- cumsum(1L - y_sorted)

    p_curve <- cum_tp / seq_len(n)
    r_curve <- cum_tp / K_true

    # Trapezoidal integration for PR-AUC
    # Prepend (r=0, p=p_curve[1])
    r_pts <- c(0.0, r_curve)
    p_pts <- c(p_curve[1], p_curve)
    dr <- diff(r_pts)
    avg_p <- (p_pts[-1] + p_pts[-length(p_pts)]) / 2.0
    pr_auc <- sum(dr * avg_p)
  }

  list(
    precision    = as.numeric(prec),
    recall       = as.numeric(rec),
    f1           = as.numeric(f1),
    realised_fdr = as.numeric(fdr),
    pr_auc       = as.numeric(pr_auc)
  )
}

# ------------------------------------------------------------------------------
# 4. Collective Anomalies Evaluation
# ------------------------------------------------------------------------------
#' Evaluate Collective Anomalies Detection Metrics
#'
#' @param truth_coll Data frame of true collective anomaly segments (start, end)
#' @param est_coll Data frame of estimated collective anomaly segments (start, end)
#' @param min_overlap Minimum Jaccard overlap fraction (default 0.5)
#'
#' @return A list with segment-level precision, recall, F1, start error, end error
#' @export
evaluate_collective_anomalies <- function(truth_coll, est_coll, min_overlap = 0.5) {
  has_truth <- !is.null(truth_coll) && is.data.frame(truth_coll) && nrow(truth_coll) > 0
  has_est   <- !is.null(est_coll) && is.data.frame(est_coll) && nrow(est_coll) > 0

  if (!has_truth && !has_est) {
    return(list(
      precision   = 1.0,
      recall      = 1.0,
      f1          = 1.0,
      start_error = 0.0,
      end_error   = 0.0
    ))
  }

  if (!has_truth && has_est) {
    return(list(
      precision   = 0.0,
      recall      = 1.0,
      f1          = 0.0,
      start_error = NA_real_,
      end_error   = NA_real_
    ))
  }

  if (has_truth && !has_est) {
    return(list(
      precision   = 1.0,
      recall      = 0.0,
      f1          = 0.0,
      start_error = NA_real_,
      end_error   = NA_real_
    ))
  }

  N_tr <- nrow(truth_coll)
  N_es <- nrow(est_coll)

  # Match segments based on Jaccard overlap >= min_overlap
  matched_es <- logical(N_es)
  matched_tr <- logical(N_tr)
  start_diffs <- numeric(0)
  end_diffs   <- numeric(0)

  for (i in seq_len(N_es)) {
    s_e <- est_coll$start[i]
    e_e <- est_coll$end[i]
    for (j in seq_len(N_tr)) {
      s_t <- truth_coll$start[j]
      e_t <- truth_coll$end[j]
      inter <- max(0L, min(e_e, e_t) - max(s_e, s_t) + 1L)
      if (inter > 0) {
        union_len <- max(e_e, e_t) - min(s_e, s_t) + 1L
        jacc <- inter / union_len
        if (jacc >= min_overlap) {
          matched_es[i] <- TRUE
          matched_tr[j] <- TRUE
          start_diffs <- c(start_diffs, abs(s_e - s_t))
          end_diffs   <- c(end_diffs, abs(e_e - e_t))
        }
      }
    }
  }

  prec <- mean(matched_es)
  rec  <- mean(matched_tr)
  f1   <- if (prec + rec > 0) 2.0 * (prec * rec) / (prec + rec) else 0.0

  list(
    precision   = as.numeric(prec),
    recall      = as.numeric(rec),
    f1          = as.numeric(f1),
    start_error = if (length(start_diffs) > 0) mean(start_diffs) else NA_real_,
    end_error   = if (length(end_diffs) > 0) mean(end_diffs) else NA_real_
  )
}

# ------------------------------------------------------------------------------
# 5. Four-Class Labelling Evaluation
# ------------------------------------------------------------------------------
#' Evaluate Four-Class Classification
#'
#' @param truth_class Factor of length n with levels c("nominal", "shift", "collective", "point")
#' @param est_class Factor of length n with levels c("nominal", "shift", "collective", "point")
#'
#' @return A list with confusion matrix, per-class F1, and Macro-F1
#' @export
evaluate_four_class <- function(truth_class, est_class) {
  truth_factor <- factor(as.character(truth_class), levels = CPD_CLASS_LEVELS)
  est_factor   <- factor(as.character(est_class), levels = CPD_CLASS_LEVELS)

  conf_mat <- table(Truth = truth_factor, Predicted = est_factor)

  # Compute per-class Precision, Recall, F1
  classes <- CPD_CLASS_LEVELS
  per_class <- list()
  f1_vec <- numeric(length(classes))

  for (i in seq_along(classes)) {
    c_name <- classes[i]
    tp <- conf_mat[c_name, c_name]
    fp <- sum(conf_mat[, c_name]) - tp
    fn <- sum(conf_mat[c_name, ]) - tp

    prec <- if (tp + fp > 0) tp / (tp + fp) else 1.0
    rec  <- if (tp + fn > 0) tp / (tp + fn) else 1.0
    f1   <- if (prec + rec > 0) 2.0 * (prec * rec) / (prec + rec) else 0.0

    per_class[[c_name]] <- list(precision = prec, recall = rec, f1 = f1)
    f1_vec[i] <- f1
  }

  macro_f1 <- mean(f1_vec)

  list(
    confusion_matrix = conf_mat,
    per_class        = per_class,
    macro_f1         = as.numeric(macro_f1)
  )
}

# ------------------------------------------------------------------------------
# 6. Comprehensive Metrics Master Aggregator
# ------------------------------------------------------------------------------
#' Master Metrics Evaluation Function
#'
#' Evaluates a `cpd_fit` object against ground truth from `simulate_series`.
#'
#' @param truth The `truth` list returned by `simulate_series`
#' @param fit The `cpd_fit` object returned by proposed method or baseline
#'
#' @return A comprehensive named list of evaluation metrics
#' @export
evaluate_all_metrics <- function(truth, fit) {
  if (inherits(truth, "cpd_fit") && is.list(fit)) {
    tmp <- truth
    truth <- fit
    fit <- tmp
  }
  stopifnot(is.list(truth), inherits(fit, "cpd_fit"))
  n <- length(truth$class)

  # 1. Change points (shift class)
  cp_eval <- evaluate_change_points(truth$cp, fit$cp, n = n)

  # 2. Point anomalies
  pt_eval <- evaluate_point_anomalies(
    truth_pt = truth$point,
    est_pt   = fit$point,
    o_t      = fit$scores$s_out,
    n        = n
  )

  # 3. Collective anomalies
  coll_eval <- evaluate_collective_anomalies(truth$collective, fit$collective)

  # 4. Four-class labelling
  clf_eval <- evaluate_four_class(truth$class, fit$class)

  list(
    method      = fit$method,
    runtime_sec = fit$runtime_sec,
    cp          = cp_eval,
    point       = pt_eval,
    collective  = coll_eval,
    four_class  = clf_eval
  )
}
