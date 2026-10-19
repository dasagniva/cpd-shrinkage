# ==============================================================================
# 00_utils.R
# Core utilities: cpd_fit output schema constructor, validator, RNG stream
# helper (L'Ecuyer-CMRG), and lightweight structured logging.
# ==============================================================================

#' Valid classification categories for each time index
CPD_CLASS_LEVELS <- c("nominal", "shift", "collective", "point")

#' Construct a new `cpd_fit` object (Section 4.7)
#'
#' @param method Character identifier for the algorithm (e.g. "proposed", "decafs", "abco")
#' @param cp Integer vector of estimated change point locations (sustained shift class)
#' @param cp_int Data frame with columns `start` and `end` for interval estimates (optional)
#' @param point Integer vector of estimated isolated point anomaly locations
#' @param collective Data frame with columns `start` and `end` for collective anomaly segments
#' @param class Factor of length n with levels c("nominal", "shift", "collective", "point")
#' @param scores Data frame with columns `t`, `s_cp`, `s_out` (NA if method provides no score)
#' @param runtime_sec Execution time in seconds
#' @param diag Arbitrary list of diagnostics (e.g. R-hat, ESS, MCMC chains, convergence flags)
#'
#' @return An object of class `cpd_fit` and `list`
#' @export
new_cpd_fit <- function(
    method,
    cp = integer(0),
    cp_int = NULL,
    point = integer(0),
    collective = NULL,
    class = factor(character(0), levels = CPD_CLASS_LEVELS),
    scores = NULL,
    runtime_sec = 0.0,
    diag = list()
) {
  # Normalize cp
  cp <- as.integer(cp)

  # Normalize cp_int
  if (is.null(cp_int) || nrow(cp_int) == 0) {
    cp_int <- data.frame(start = integer(0), end = integer(0))
  } else {
    stopifnot(is.data.frame(cp_int), c("start", "end") %in% names(cp_int))
    cp_int$start <- as.integer(cp_int$start)
    cp_int$end <- as.integer(cp_int$end)
  }

  # Normalize point
  point <- as.integer(point)

  # Normalize collective
  if (is.null(collective) || nrow(collective) == 0) {
    collective <- data.frame(start = integer(0), end = integer(0))
  } else {
    stopifnot(is.data.frame(collective), c("start", "end") %in% names(collective))
    collective$start <- as.integer(collective$start)
    collective$end <- as.integer(collective$end)
  }

  # Normalize class
  if (!is.factor(class)) {
    class <- factor(as.character(class), levels = CPD_CLASS_LEVELS)
  } else {
    levels(class) <- CPD_CLASS_LEVELS
  }

  # Normalize scores
  if (is.null(scores) || nrow(scores) == 0) {
    scores <- data.frame(
      t = integer(0),
      s_cp = numeric(0),
      s_out = numeric(0)
    )
  } else {
    stopifnot(is.data.frame(scores), c("t", "s_cp", "s_out") %in% names(scores))
    scores$t <- as.integer(scores$t)
    scores$s_cp <- as.numeric(scores$s_cp)
    scores$s_out <- as.numeric(scores$s_out)
  }

  out <- list(
    method      = as.character(method)[1],
    cp          = cp,
    cp_int      = cp_int,
    point       = point,
    collective  = collective,
    class       = class,
    scores      = scores,
    runtime_sec = as.numeric(runtime_sec)[1],
    diag        = as.list(diag)
  )

  class(out) <- c("cpd_fit", "list")
  out
}

#' Validate a `cpd_fit` object
#'
#' @param fit An object expected to be of class `cpd_fit`
#' @param n Optional integer specifying expected series length. If provided, checks length of class and range of indices.
#'
#' @return The validated `fit` object invisibly, or errors with a descriptive message.
#' @export
validate_cpd_fit <- function(fit, n = NULL) {
  if (!inherits(fit, "cpd_fit") || !is.list(fit)) {
    stop("Object must inherit from 'cpd_fit' and be a list", call. = FALSE)
  }

  required_fields <- c(
    "method", "cp", "cp_int", "point", "collective",
    "class", "scores", "runtime_sec", "diag"
  )
  missing_fields <- setdiff(required_fields, names(fit))
  if (length(missing_fields) > 0) {
    stop("Missing required fields in cpd_fit: ", paste(missing_fields, collapse = ", "), call. = FALSE)
  }

  if (!is.character(fit$method) || length(fit$method) != 1 || is.na(fit$method)) {
    stop("'method' must be a single non-NA character string", call. = FALSE)
  }

  if (!is.integer(fit$cp)) {
    stop("'cp' must be an integer vector", call. = FALSE)
  }

  if (!is.data.frame(fit$cp_int) || !all(c("start", "end") %in% names(fit$cp_int))) {
    stop("'cp_int' must be a data.frame with columns 'start' and 'end'", call. = FALSE)
  }

  if (!is.integer(fit$point)) {
    stop("'point' must be an integer vector", call. = FALSE)
  }

  if (!is.data.frame(fit$collective) || !all(c("start", "end") %in% names(fit$collective))) {
    stop("'collective' must be a data.frame with columns 'start' and 'end'", call. = FALSE)
  }

  if (!is.factor(fit$class) || !identical(levels(fit$class), CPD_CLASS_LEVELS)) {
    stop("'class' must be a factor with levels: ", paste(CPD_CLASS_LEVELS, collapse = ", "), call. = FALSE)
  }

  if (!is.data.frame(fit$scores) || !all(c("t", "s_cp", "s_out") %in% names(fit$scores))) {
    stop("'scores' must be a data.frame with columns 't', 's_cp', 's_out'", call. = FALSE)
  }

  if (!is.numeric(fit$runtime_sec) || length(fit$runtime_sec) != 1 || fit$runtime_sec < 0) {
    stop("'runtime_sec' must be a single non-negative numeric value", call. = FALSE)
  }

  if (!is.list(fit$diag)) {
    stop("'diag' must be a list", call. = FALSE)
  }

  # Interval ordering checks
  if (nrow(fit$cp_int) > 0) {
    if (any(fit$cp_int$start > fit$cp_int$end)) {
      stop("cp_int contains rows where start > end", call. = FALSE)
    }
  }

  if (nrow(fit$collective) > 0) {
    if (any(fit$collective$start > fit$collective$end)) {
      stop("collective contains rows where start > end", call. = FALSE)
    }
  }

  # Length and index validation when series length n is given
  if (!is.null(n)) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 1)
    n <- as.integer(n)

    if (length(fit$class) != n) {
      stop(sprintf("Length of class (%d) does not match expected length n (%d)", length(fit$class), n), call. = FALSE)
    }

    if (length(fit$cp) > 0 && any(fit$cp < 1 | fit$cp > n)) {
      stop("One or more 'cp' indices are outside [1, n]", call. = FALSE)
    }

    if (length(fit$point) > 0 && any(fit$point < 1 | fit$point > n)) {
      stop("One or more 'point' indices are outside [1, n]", call. = FALSE)
    }

    if (nrow(fit$collective) > 0) {
      if (any(fit$collective$start < 1 | fit$collective$end > n)) {
        stop("One or more 'collective' intervals are outside [1, n]", call. = FALSE)
      }
    }

    if (nrow(fit$cp_int) > 0) {
      if (any(fit$cp_int$start < 1 | fit$cp_int$end > n)) {
        stop("One or more 'cp_int' intervals are outside [1, n]", call. = FALSE)
      }
    }

    if (nrow(fit$scores) > 0) {
      if (nrow(fit$scores) != n) {
        stop(sprintf("Number of score rows (%d) does not match n (%d)", nrow(fit$scores), n), call. = FALSE)
      }
      if (any(fit$scores$t != seq_len(n))) {
        stop("Scores 't' column must match 1:n", call. = FALSE)
      }
    }
  }

  invisible(fit)
}

#' Initialize reproducible parallel RNG stream using L'Ecuyer-CMRG (Section 1.5, P0)
#'
#' Ensures non-overlapping independent RNG streams for each simulation cell and replicate.
#'
#' @param cell_id Integer cell identifier in the simulation grid (>= 1)
#' @param rep_id Integer replicate identifier within the cell (>= 1)
#' @param base_seed Master integer seed
#'
#' @return The stream state invisibly
#' @export
init_rng_stream <- function(cell_id, rep_id, base_seed = 20261016L) {
  stopifnot(cell_id >= 1, rep_id >= 1)

  # Switch to L'Ecuyer-CMRG
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(base_seed))

  # Retrieve initial stream
  current_stream <- .Random.seed

  # Advance by cell_id streams
  for (c_idx in seq_len(cell_id)) {
    current_stream <- parallel::nextRNGStream(current_stream)
  }

  # Advance within cell by rep_id substreams
  for (r_idx in seq_len(rep_id)) {
    current_stream <- parallel::nextRNGSubStream(current_stream)
  }

  assign(".Random.seed", current_stream, envir = .GlobalEnv)
  invisible(current_stream)
}

#' Robust location/scale standardisation
#'
#' Centres at the median and scales by a robust innovation-scale estimate,
#' MAD(diff(y)) / sqrt(2), which is insensitive to level shifts and isolated
#' outliers. Every model and baseline that assumes unit-scale noise is fitted
#' on `z` and back-transformed with `center` and `scale`.
#'
#' @param y Numeric vector
#'
#' @return A list with `z` (standardised series), `center` and `scale`
#' @export
robust_standardize <- function(y) {
  stopifnot(is.numeric(y), length(y) >= 3)
  center <- stats::median(y)
  scale <- stats::mad(diff(y)) / sqrt(2)
  if (!is.finite(scale) || scale <= 1e-12) scale <- stats::mad(y)
  if (!is.finite(scale) || scale <= 1e-12) scale <- stats::sd(y)
  if (!is.finite(scale) || scale <= 1e-12) scale <- 1.0
  list(z = (y - center) / scale, center = center, scale = scale)
}

#' Simple structured logging
#'
#' @param ... Message components passed to cat / paste
#' @param level Character log level: "INFO", "WARN", "ERROR"
#' @export
log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(...)
  cat(sprintf("[%s] [%-5s] %s\n", ts, toupper(level), msg))
}
