# ==============================================================================
# test_geweke_getting_it_right.R
# Geweke (2004) "Getting it right" joint-distribution test of the production Gibbs
# sampler (Phase P4 gate), with stationary starts.
#
# If every conditional update leaves p(theta | y) invariant, the successive-conditional
# kernel -- one sweep of run_gibbs_chain() given y, then y ~ p(y | theta) -- leaves the
# joint p(theta, y) invariant. Each of N_CHAINS chains is therefore started from an exact
# joint draw (theta_0, y_0), including the Makalic-Schmidt and Polya-Gamma auxiliaries
# drawn from their conditionals, and run for K sweeps; under a correct sampler the final
# states are independent exact draws from p(theta, y), whatever the mixing. Their test
# functions are compared with an independent forward (marginal-conditional) sample by
# two-sample z-tests, Bonferroni-adjusted. A wrong update makes the distribution drift
# away from p(theta, y) within the K sweeps. (A single long chain, as in Geweke's
# original design, cannot be used here: under the horseshoe+ prior the joint chain
# random-walks over levels and scales spanning many orders of magnitude.)
#
# The successive-conditional kernel is the production chain itself (run_gibbs_chain with
# geweke_fun and init_state), so the test covers the update order, the marginal-scale
# Metropolis steps, dsp's dynamic horseshoe, the non-centred rho step, kappa and the
# CAPA offset. The model matches the code: mu_1 ~ N(0, MU1_PRIOR_VAR); omega_t ~ N(0,
# sigma^2 v_t), v_t = clip(tau^2 kappa_t^2 lambda_t^2 eta_t^2) with lambda_t, eta_t ~
# C+(0, 1), tau ~ C+(0, 1/sqrt(n)); (phi + 1)/2 ~ Beta(10, 2) truncated to |phi| <= 0.95;
# dynamic outliers psi_t ~ N(0, clip(exp(m_h + h~_t + delta a_t))), h~_1 = e_1,
# h~_t = rho h~_{t-1} + e_t, e_t = log(C+(0, 1)^2), exp(m_h / 2) ~ C+(0, 1/sqrt(n)),
# (rho + 1)/2 ~ Beta(10, 2); clip() is the sampler's numerical guard [1e-12, 1e12]. The
# spec's Jeffreys prior on sigma^2 is improper, so the test uses sigma^2 ~ IG(3, 2)
# (passed to the production update_sigma2).
#
# The horseshoe(+) priors give mu, omega, psi and y no finite mean, so the test functions
# are bounded (atan, indicators, correlations) or logs of scales. A planted error (the
# sampler given a different sigma^2 prior) checks that the test can fail.
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(posterior)
})
source("R/00_utils.R")
source("R/32_gibbs_updates.R")
source("R/40_scoring.R")
source("R/41_bfdr.R")
source("R/42_classify.R")
source("R/31_gibbs_sampler.R")

n <- 40L
N_FWD <- 100000L # independent marginal-conditional draws
N_CHAINS <- 4000L # chains, each started from an exact joint draw
K <- 100L # successive-conditional sweeps per chain
SIGMA2_PRIOR <- c(3, 2)
PHI_PRIOR <- c(10, 2)
ALPHA <- 0.01 # family-wise level (Bonferroni over test functions)
N_WORKERS <- cpd_workers(if (nzchar(Sys.getenv("GEWEKE_WORKERS"))) as.integer(Sys.getenv("GEWEKE_WORKERS")) else NULL)

kappa_vec <- rep(1.0, n)
kappa_vec[15:22] <- 2.0 # NSP-style inflated region
a_vec <- rep(0.0, n)
a_vec[c(10, 30)] <- 1.0 # CAPA point flags
T_IN <- 20L  # increment inside the inflated region
T_OUT <- 6L  # increment outside it

configs <- list(
  production = list(dynamic = TRUE, delta = 1.0, kappa = kappa_vec, a_vec = a_vec, tau_fixed = NULL,
                    sampler_sigma2_prior = SIGMA2_PRIOR),
  plugin_indep = list(dynamic = FALSE, delta = 0.0, kappa = rep(1.0, n), a_vec = rep(0.0, n),
                      tau_fixed = (0.02 / 0.98) / sqrt(n), sampler_sigma2_prior = SIGMA2_PRIOR),
  planted_error = list(dynamic = TRUE, delta = 1.0, kappa = kappa_vec, a_vec = a_vec, tau_fixed = NULL,
                       sampler_sigma2_prior = c(3, 3))
)

test_functions <- function(s, cfg) {
  w <- diff(s$mu)
  dy <- diff(s$y)
  g <- c(
    phi           = s$phi,
    phi_sq        = s$phi^2,
    log_sigma     = log(s$sigma),
    log_sigma_sq  = log(s$sigma)^2,
    atan_w_in     = atan(w[T_IN - 1L] / s$sigma),
    atan_w_out    = atan(w[T_OUT - 1L] / s$sigma),
    big_w_in      = as.numeric(abs(w[T_IN - 1L]) > s$sigma),
    big_w_out     = as.numeric(abs(w[T_OUT - 1L]) > s$sigma),
    log_le_in     = log(s$lambda_omega[T_IN - 1L] * s$eta_omega[T_IN - 1L]),
    log_le_out    = log(s$lambda_omega[T_OUT - 1L] * s$eta_omega[T_OUT - 1L]),
    log_nu2_flag  = log(s$nu2[10]),
    log_nu2_plain = log(s$nu2[20]),
    atan_psi_flag = atan(s$psi[10] / s$sigma),
    atan_psi      = atan(s$psi[20] / s$sigma),
    acf1_dy       = sum(dy[-1] * dy[-length(dy)]) / sum(dy^2),
    log_mad_dy    = log(stats::median(abs(dy))),
    atan_dy_in    = atan(dy[T_IN - 1L] / s$sigma)
  )
  if (is.null(cfg$tau_fixed)) g <- c(g, log_tau_omega = log(s$tau_omega))
  if (cfg$dynamic) g <- c(g, rho = s$rho, m_h = s$m_h) else g <- c(g, log_tau_psi = log(s$tau_psi))
  g
}

rhalfcauchy <- function(k, s = 1) abs(stats::rcauchy(k, 0, s))
clip <- function(v) pmin(1e12, pmax(1e-12, v))
# Makalic-Schmidt auxiliary of a half-Cauchy scale l: nu | l ~ IG(1, 1 + 1/l^2)
aux_of <- function(l) 1 / stats::rgamma(length(l), shape = 1, rate = 1 + 1 / l^2)

forward_draw <- function(cfg) {
  repeat {
    phi <- 2 * stats::rbeta(1, PHI_PRIOR[1], PHI_PRIOR[2]) - 1
    if (abs(phi) <= 0.95) break # truncation by rejection, not clipping
  }
  sigma <- sqrt(1 / stats::rgamma(1, shape = SIGMA2_PRIOR[1], rate = SIGMA2_PRIOR[2]))
  tau <- if (is.null(cfg$tau_fixed)) rhalfcauchy(1, 1 / sqrt(n)) else cfg$tau_fixed
  lam <- rhalfcauchy(n - 1)
  eta <- rhalfcauchy(n - 1)
  w <- stats::rnorm(n - 1, 0, sigma * sqrt(clip((tau * cfg$kappa[-1] * lam * eta)^2)))
  mu <- cumsum(c(stats::rnorm(1, 0, sqrt(MU1_PRIOR_VAR)), w))

  st <- list(tau_omega = tau, lambda_omega = lam, eta_omega = eta,
             aux_lambda_omega = aux_of(lam), aux_eta_omega = aux_of(eta))
  rho <- m_h <- tau_psi <- NA_real_
  if (cfg$dynamic) {
    m_h <- 2 * log(rhalfcauchy(1, 1 / sqrt(n)))
    rho <- 2 * stats::rbeta(1, 10, 2) - 1
    e <- log(rhalfcauchy(n)^2)
    h_til <- as.numeric(stats::filter(e, rho, method = "recursive"))
    nu2 <- clip(exp(m_h + h_til + cfg$delta * cfg$a_vec))
    # dsp state, with the Polya-Gamma auxiliaries drawn given the innovations
    st$dsp_state <- list(ht = matrix(m_h + h_til, ncol = 1), dhs_mean = m_h, dhs_phi = rho,
                         sigma_eta_t = matrix(1 / sqrt(pgdraw::pgdraw(1, e[-1])), ncol = 1),
                         sigma_eta_0 = 1 / sqrt(pgdraw::pgdraw(1, e[1])), dhs_mean0 = m_h)
  } else {
    tau_psi <- rhalfcauchy(1, 1 / sqrt(n))
    lp <- rhalfcauchy(n)
    ep <- rhalfcauchy(n)
    nu2 <- clip((tau_psi * lp * ep)^2)
    st <- c(st, list(tau_psi = tau_psi, lambda_psi = lp, eta_psi = ep, aux_lambda_psi = aux_of(lp),
                     aux_eta_psi = aux_of(ep)))
  }
  psi <- stats::rnorm(n, 0, sqrt(nu2))
  u <- c(stats::rnorm(1, 0, sigma / sqrt(1 - phi^2)), sigma * stats::rnorm(n - 1))
  y <- mu + psi + as.numeric(stats::filter(u, phi, method = "recursive"))
  c(st, list(y = y, mu = mu, psi = psi, phi = phi, sigma = sigma, nu2 = nu2, tau_psi = tau_psi,
             rho = rho, m_h = m_h))
}

run_config <- function(name, cfg, seed) {
  set.seed(seed)
  n_g <- length(test_functions(forward_draw(cfg), cfg))
  fwd <- t(vapply(seq_len(N_FWD), function(i) test_functions(forward_draw(cfg), cfg), numeric(n_g)))

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed + 1L)
  final <- mclapply(seq_len(N_CHAINS), function(k) {
    s0 <- forward_draw(cfg)
    out <- tryCatch(run_gibbs_chain(
      chain_id = k, y = s0$y, kappa = cfg$kappa, a_vec = cfg$a_vec, n_iter = K, burn_in = K - 1L,
      dynamic = cfg$dynamic, delta = cfg$delta, tau_omega_fixed = cfg$tau_fixed,
      phi_prior = PHI_PRIOR, sigma2_prior = cfg$sampler_sigma2_prior, init_state = s0,
      geweke_fun = function(s) test_functions(s, cfg)
    )$geweke_draws, error = function(e) NULL)
    if (is.null(out)) rep(NA_real_, n_g) else out[nrow(out), ]
  }, mc.cores = N_WORKERS, mc.set.seed = TRUE)
  RNGkind("default")
  fin <- do.call(rbind, final)
  n_fail <- sum(!stats::complete.cases(fin))
  if (n_fail > 0.001 * N_CHAINS) stop(sprintf("[%s] %d of %d chains failed", name, n_fail, N_CHAINS))
  fin <- fin[stats::complete.cases(fin), , drop = FALSE]

  tab <- do.call(rbind, lapply(colnames(fwd), function(v) {
    x1 <- fwd[, v]
    x2 <- fin[, v]
    se1 <- stats::sd(x1) / sqrt(length(x1))
    se2 <- stats::sd(x2) / sqrt(length(x2))
    z <- (mean(x1) - mean(x2)) / sqrt(se1^2 + se2^2)
    data.frame(config = name, variable = v, mean_forward = mean(x1), se_forward = se1, mean_chain = mean(x2),
               se_chain = se2, z = z, p_value = 2 * stats::pnorm(-abs(z)))
  }))
  tab$p_bonferroni <- pmin(1, tab$p_value * nrow(tab))
  attr(tab, "n_failed_chains") <- n_fail
  tab
}

t0 <- proc.time()[3]
results <- list()
for (k in seq_along(configs)) {
  nm <- names(configs)[k]
  cat(sprintf("\n=== Geweke test, configuration '%s' (%d forward draws; %d chains of %d sweeps from exact joint draws) ===\n",
              nm, N_FWD, N_CHAINS, K))
  results[[nm]] <- run_config(nm, configs[[nm]], seed = 20261021L + 100L * k)
  print(results[[nm]][, c("variable", "mean_forward", "mean_chain", "se_chain", "z", "p_bonferroni")],
        digits = 3, row.names = FALSE)
}
runtime_min <- (proc.time()[3] - t0) / 60

verdict <- function(tab) {
  list(min_p_bonferroni = min(tab$p_bonferroni), n_failed_chains = attr(tab, "n_failed_chains"),
       passed = min(tab$p_bonferroni) >= ALPHA)
}
v <- lapply(results, verdict)
passed <- v$production$passed && v$plugin_indep$passed
power_ok <- v$planted_error$min_p_bonferroni < ALPHA

cat(sprintf("\nproduction: %s | plugin_indep: %s | planted error detected: %s | %.1f min\n",
            v$production$passed, v$plugin_indep$passed, power_ok, runtime_min))

dir.create("output", showWarnings = FALSE)
saveRDS(list(
  test_table   = do.call(rbind, results),
  verdicts     = v,
  passed       = passed,
  power_ok     = power_ok,
  settings     = list(n = n, N_FWD = N_FWD, N_CHAINS = N_CHAINS, K = K, alpha = ALPHA,
                      sigma2_prior = SIGMA2_PRIOR, phi_prior = PHI_PRIOR),
  runtime_min  = runtime_min,
  timestamp    = Sys.time()
), "output/geweke_test_results.rds")

if (!power_ok) stop("Geweke test did not detect the planted error: it has too little power")
if (!passed) stop("Geweke joint distribution test FAILED")
cat(">>> GEWEKE GETTING-IT-RIGHT TEST: PASSED (planted error detected) <<<\n")
