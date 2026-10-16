# ==============================================================================
# 20_frontend_nsp.R
# Significance Front-End: Narrowest Significance Pursuit (NSP) wrapper.
# Implements Section 3.2 and Section 5 (Phase P2).
# ==============================================================================

#' Run Narrowest Significance Pursuit (NSP) to find significance intervals
#'
#' Implements Section 3.2 and Phase P2. Supports both self-normalised NSP
#' (`nsp_poly_selfnorm` with `deg = 0`) as specified in Section 2, and the
#' autoregressive extension (`nsp_poly_ar` with `ord = 1, deg = 0`).
#'
#' @param y Numeric vector of observed time series
#' @param alpha Significance level alpha_NSP (default 0.1; sensitivity in {0.05, 0.10})
#' @param method Character: "selfnorm" (default, per Section 2) or "ar" (for AR noise)
#' @param prewhiten Logical; if TRUE and method == "selfnorm", prewhitens y using estimated phi
#' @param M Number of random intervals to draw (default 1000)
#' @param min_size Minimum interval size (default 20)
#' @param power Power parameter for self-normalisation (default 0.5)
#' @param ... Additional arguments passed to the underlying NSP function
#'
#' @return A list with elements:
#'   \item{intervals}{Data frame with integer columns `start` and `end`}
#'   \item{threshold}{The threshold used by NSP}
#'   \item{method}{Method used: "selfnorm" or "ar"}
#'   \item{raw}{Raw NSP return object}
#' @export
run_frontend_nsp <- function(
    y,
    alpha = 0.1,
    method = c("selfnorm", "ar"),
    prewhiten = FALSE,
    M = 1000,
    min_size = 20,
    power = 0.5,
    ...
) {
  stopifnot(is.numeric(y), length(y) >= min_size, alpha > 0, alpha < 1)
  method <- match.arg(method)

  if (method == "selfnorm") {
    y_input <- y
    if (isTRUE(prewhiten)) {
      source("R/21_frontend_capa.R", local = TRUE)
      phi_hat <- estimate_phi_robust(y)
      n <- length(y)
      y_input <- numeric(n)
      y_input[1] <- sqrt(1 - phi_hat^2) * (y[1] - median(y))
      y_input[2:n] <- (y[2:n] - median(y)) - phi_hat * (y[1:(n - 1L)] - median(y))
    }

    res <- nsp::nsp_poly_selfnorm(
      y = y_input,
      M = M,
      power = power,
      min.size = min_size,
      alpha = alpha,
      deg = 0,
      ...
    )
  } else {
    # Autoregressive NSP (Fryzlewicz 2021)
    res <- nsp::nsp_poly_ar(
      y = y,
      ord = 1,
      M = M,
      power = power,
      min.size = min_size,
      alpha = alpha,
      deg = 0,
      ...
    )
  }

  raw_ints <- res$intervals
  if (is.null(raw_ints) || nrow(raw_ints) == 0) {
    intervals_df <- data.frame(
      start = integer(0),
      end   = integer(0)
    )
  } else {
    intervals_df <- data.frame(
      start = as.integer(raw_ints$starts),
      end   = as.integer(raw_ints$ends)
    )
    intervals_df <- intervals_df[order(intervals_df$start), , drop = FALSE]
    rownames(intervals_df) <- NULL
  }

  list(
    intervals = intervals_df,
    threshold = res$threshold.used,
    method    = method,
    raw       = res
  )
}
