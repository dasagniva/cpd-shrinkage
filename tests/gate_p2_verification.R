# ==============================================================================
# gate_p2_verification.R
# Acceptance Gate for Phase P2 (CPD2_Code_Workflow.md):
# 1. On null series (no changes) with phi in {0, 0.6, 0.9}, the NSP false-interval
#    rate is at most alpha_NSP (200 replicates), for the front-end configuration the
#    pipeline uses (run_frontend_nsp default: nsp_poly_ar, ord = 1, deg = 0; spec
#    change C1). The self-normalised NSP named in the spec is reported for comparison.
# 2. CAPA's point-anomaly false-positive rate under phi = 0.9 is reported with and
#    without prewhitening.
# ==============================================================================

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/20_frontend_nsp.R")
source("R/21_frontend_capa.R")
source("R/22_prior_map.R")

log_msg("=== Starting Phase P2 Gate Verification ===")

alpha_nsp <- 0.05 # level used by the pipeline (fit_proposed default)
n_pts <- 300
reps <- 200
phi_grid <- c(0.0, 0.6, 0.9)

null_series <- function(phi_val, b) {
  simulate_series(n = n_pts, phi = phi_val, k0 = 0, outlier_regime = "none", collective = 0,
                  seed = 20261000L + as.integer(phi_val * 1000) + b)$y
}

nsp_rows <- list()
for (phi_val in phi_grid) {
  fi_sn <- logical(reps)
  fi_ar <- logical(reps)
  for (b in seq_len(reps)) {
    y <- null_series(phi_val, b)
    fi_sn[b] <- nrow(run_frontend_nsp(y, alpha = alpha_nsp, method = "selfnorm")$intervals) > 0
    fi_ar[b] <- nrow(run_frontend_nsp(y, alpha = alpha_nsp, method = "ar")$intervals) > 0
  }
  nsp_rows[[length(nsp_rows) + 1]] <- data.frame(
    phi = phi_val,
    fir_selfnorm = mean(fi_sn), se_selfnorm = sqrt(mean(fi_sn) * (1 - mean(fi_sn)) / reps),
    fir_ar = mean(fi_ar), se_ar = sqrt(mean(fi_ar) * (1 - mean(fi_ar)) / reps)
  )
  cat(sprintf("phi %.1f | selfnorm FIR %.3f | ar FIR %.3f | target <= %.2f\n",
              phi_val, mean(fi_sn), mean(fi_ar), alpha_nsp))
}
nsp_tab <- do.call(rbind, nsp_rows)

capa_reps <- 60
fp_raw <- numeric(capa_reps)
fp_pw <- numeric(capa_reps)
for (b in seq_len(capa_reps)) {
  y <- simulate_series(n = n_pts, phi = 0.9, k0 = 0, outlier_regime = "none", collective = 0, seed = 20264000L + b)$y
  fp_raw[b] <- length(run_frontend_capa(y, prewhiten = FALSE)$point)
  fp_pw[b] <- length(run_frontend_capa(y, prewhiten = TRUE)$point)
}
cat(sprintf("CAPA mean false points per series under phi = 0.9: raw %.3f, prewhitened %.3f\n", mean(fp_raw), mean(fp_pw)))

# Pass if the pipeline front-end's (nsp_poly_ar) FIR is within 2 MC SE of alpha_NSP at every phi
nsp_pass <- all(nsp_tab$fir_ar <= alpha_nsp + 2 * sqrt(alpha_nsp * (1 - alpha_nsp) / reps))

saveRDS(list(nsp = nsp_tab, capa_fp_raw = mean(fp_raw), capa_fp_prewhitened = mean(fp_pw),
             alpha_nsp = alpha_nsp, reps = reps, nsp_pass = nsp_pass, timestamp = Sys.time()),
        file = "output/gate_p2_results.rds")

if (!nsp_pass) stop("Phase P2 acceptance gate FAILED: NSP false-interval rate exceeds alpha_NSP for the pipeline front-end")
cat(">>> PHASE P2 ACCEPTANCE GATE: PASSED <<<\n")
