# ==============================================================================
# 41_bfdr.R
# Bayesian False Discovery Rate (BFDR) Selection Rule.
# Implements Section 3.5.
# ==============================================================================

#' Bayesian False Discovery Rate (BFDR) Selection Rule
#'
#' Given posterior inclusion / exceedance probabilities p_t in [0, 1],
#' sorts p_t in descending order: p_{(1)} >= p_{(2)} >= ... >= p_{(n)}.
#'
#' Selects the largest k in {1, ..., n} such that the posterior expected
#' false discovery rate does not exceed nominal alpha:
#'   (1 / k) * sum_{j=1}^k (1 - p_{(j)}) <= alpha
#'
#' If no k satisfies the condition (for example, if 1 - p_{(1)} > alpha),
#' no points are flagged (k = 0).
#'
#' @param p_vec Numeric vector of posterior probabilities in [0, 1].
#' @param alpha Target Bayesian FDR level (default 0.05). Must be in (0, 1).
#' @param preserve_order Logical; if TRUE (default), flagged indices are returned
#'   in increasing order of original position (e.g. time order). If FALSE,
#'   they are returned in descending order of posterior probability.
#'
#' @return An object of class `bfdr_result` with fields:
#'   \item{flagged}{Integer vector of flagged indices}
#'   \item{k}{Number of discoveries (integer)}
#'   \item{threshold}{Cutoff probability p_{(k)} (or 1.0 if k = 0)}
#'   \item{est_fdr}{Estimated posterior FDR (or 0.0 if k = 0)}
#'   \item{alpha}{Target nominal alpha}
#'   \item{cum_fdr}{Full cumulative FDR vector along sorted probabilities}
#'   \item{ranked_indices}{Indices of p_vec sorted in descending order of p_t}
#' @export
bayesian_fdr <- function(p_vec, alpha = 0.05, preserve_order = TRUE) {
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)

  if (length(p_vec) == 0) {
    res <- list(
      flagged        = integer(0),
      k              = 0L,
      threshold      = 1.0,
      est_fdr        = 0.0,
      alpha          = alpha,
      cum_fdr        = numeric(0),
      ranked_indices = integer(0)
    )
    class(res) <- "bfdr_result"
    return(res)
  }

  stopifnot(is.numeric(p_vec))
  if (any(is.na(p_vec))) {
    warning("NA values detected in p_vec; setting them to 0.0")
    p_vec[is.na(p_vec)] <- 0.0
  }
  if (any(p_vec < 0.0 | p_vec > 1.0)) {
    stop("All probabilities in p_vec must lie in [0, 1]")
  }

  n <- length(p_vec)
  ord <- order(p_vec, decreasing = TRUE)
  sorted_p <- p_vec[ord]

  # Posterior expected false discovery rate for flagging top j elements:
  # FDR(j) = (1 / j) * sum_{m=1}^j (1 - sorted_p[m])
  cum_fdr <- cumsum(1.0 - sorted_p) / seq_len(n)

  # Find largest k such that cum_fdr[k] <= alpha
  valid_k <- which(cum_fdr <= alpha)

  if (length(valid_k) == 0) {
    k_star <- 0L
    flagged <- integer(0)
    threshold <- 1.0
    est_fdr <- 0.0
  } else {
    k_star <- as.integer(max(valid_k))
    flagged <- ord[seq_len(k_star)]
    threshold <- as.numeric(sorted_p[k_star])
    est_fdr <- as.numeric(cum_fdr[k_star])
    if (isTRUE(preserve_order)) {
      flagged <- sort(flagged)
    }
  }

  res <- list(
    flagged        = as.integer(flagged),
    k              = k_star,
    threshold      = threshold,
    est_fdr        = est_fdr,
    alpha          = alpha,
    cum_fdr        = cum_fdr,
    ranked_indices = as.integer(ord)
  )
  class(res) <- "bfdr_result"
  res
}

#' Extract Flagged Indices from a BFDR Result
#'
#' @param x An object of class `bfdr_result` or a numeric vector of probabilities
#' @param alpha Target alpha if x is a numeric vector
#'
#' @return Integer vector of flagged indices
#' @export
bfdr_flagged <- function(x, alpha = 0.05) {
  if (inherits(x, "bfdr_result")) {
    return(x$flagged)
  }
  bayesian_fdr(x, alpha = alpha)$flagged
}

#' Print method for bfdr_result
#'
#' @param x A `bfdr_result` object
#' @param ... Ignored
#' @export
print.bfdr_result <- function(x, ...) {
  cat(sprintf(
    "Bayesian FDR (alpha = %.3f): %d discoveries (threshold >= %.4f, est. FDR = %.4f)\n",
    x$alpha, x$k, x$threshold, x$est_fdr
  ))
  if (x$k > 0) {
    preview <- if (x$k <= 10) paste(x$flagged, collapse = ", ") else paste(c(x$flagged[1:10], "..."), collapse = ", ")
    cat(sprintf("Flagged indices: %s\n", preview))
  }
  invisible(x)
}
