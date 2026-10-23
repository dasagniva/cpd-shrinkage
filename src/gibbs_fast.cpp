// =============================================================================
// gibbs_fast.cpp
// O(n) linear algebra for the production Gibbs sampler (spec P4: move hot loops to Rcpp).
//
// 1. mu_psi_factor_cpp: the interleaved precision Q of x = (mu_1, psi_1, ..., mu_n, psi_n)
//    (bandwidth 3; entries exactly as build_mu_psi_system in R/32_gibbs_updates.R) and
//    its banded Cholesky Q = L L', the forward solve w = L^{-1} b and the marginal
//    log-likelihood used by the global-scale Metropolis steps.
// 2. band_backsolve_cpp: x = L^{-T} r, giving a draw x = L^{-T}(w + z) ~ N(Q^{-1} b, Q^{-1}).
// 3. tridiag_sample_cpp: draw from N(Q^{-1} b, Q^{-1}) for tridiagonal Q (dsp's
//    log-volatility step, replacing Matrix::bandSparse + Matrix::chol).
// 4. ar1_filter_cpp: h_1 = e_1, h_t = e_t + rho h_{t-1} (stats::filter, recursive).
// =============================================================================
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

static const int KD = 3; // bandwidth of Q

// Banded Cholesky of the symmetric matrix with upper band Qu (m x 4, Qu(i, d) = Q[i][i+d]).
// Returns the lower band Lb (m x 4, Lb(i, d) = L[i][i-d]); false if Q is not positive definite.
static bool band_chol(const std::vector<double>& Qu, std::vector<double>& Lb, int m) {
  for (int i = 0; i < m; ++i) {
    const int k0 = std::max(0, i - KD);
    for (int k = k0; k <= i; ++k) {
      double s = Qu[k * 4 + (i - k)];
      for (int j = std::max(k0, k - KD); j < k; ++j) s -= Lb[i * 4 + (i - j)] * Lb[k * 4 + (k - j)];
      if (k == i) {
        if (!(s > 0.0)) return false;
        Lb[i * 4] = std::sqrt(s);
      } else {
        Lb[i * 4 + (i - k)] = s / Lb[k * 4];
      }
    }
  }
  return true;
}

// [[Rcpp::export]]
List mu_psi_factor_cpp(NumericVector y, double phi, double sigma, NumericVector v_omega,
                       NumericVector nu2, double mu1_prior_var) {
  const int n = y.size();
  const int m = 2 * n;
  const double inv_s2 = 1.0 / (sigma * sigma);
  const double off = -phi * inv_s2;

  std::vector<double> Qu((size_t)m * 4, 0.0);
  auto q = [&](int i, int d) -> double& { return Qu[(size_t)i * 4 + d]; };
  for (int t = 0; t < n; ++t) {
    const double c = (t == 0 || t == n - 1) ? 1.0 : 1.0 + phi * phi;
    q(2 * t, 0) += c * inv_s2;
    q(2 * t + 1, 0) += c * inv_s2 + 1.0 / nu2[t];
    q(2 * t, 1) = c * inv_s2; // (mu_t, psi_t)
  }
  q(0, 0) += 1.0 / mu1_prior_var;
  for (int t = 0; t < n - 1; ++t) {
    const double w = 1.0 / v_omega[t];
    q(2 * t, 0) += w;
    q(2 * t + 2, 0) += w;
    q(2 * t + 1, 1) = off;     // (psi_t, mu_{t+1})
    q(2 * t, 2) = off - w;     // (mu_t, mu_{t+1})
    q(2 * t + 1, 2) = off;     // (psi_t, psi_{t+1})
    q(2 * t, 3) = off;         // (mu_t, psi_{t+1})
  }
  for (int i = 0; i < m; ++i) q(i, 0) += 1e-8;

  // Equilibrate to unit diagonal (Qs = S Q S, S = diag(Q)^{-1/2}), factorise Qs = Ls Ls',
  // and unscale (L = S^{-1} Ls). Under the heavy-tailed scale priors the diagonal can span
  // 24 orders of magnitude; if Qs is still numerically indefinite, a relative jitter
  // 1e-12, 1e-11, ..., 1e-6 is added to its unit diagonal.
  std::vector<double> sc(m);
  for (int i = 0; i < m; ++i) sc[i] = 1.0 / std::sqrt(q(i, 0));
  std::vector<double> Qs((size_t)m * 4, 0.0);
  for (int i = 0; i < m; ++i)
    for (int d = 0; d <= KD && i + d < m; ++d) Qs[(size_t)i * 4 + d] = q(i, d) * sc[i] * sc[i + d];
  std::vector<double> Lb((size_t)m * 4, 0.0);
  bool ok = band_chol(Qs, Lb, m);
  for (double jit = 1e-12; !ok && jit <= 1.0001e-6; jit *= 10.0) {
    std::vector<double> Qj(Qs);
    for (int i = 0; i < m; ++i) Qj[(size_t)i * 4] += jit;
    ok = band_chol(Qj, Lb, m);
  }
  if (!ok) stop("mu_psi_factor_cpp: precision matrix is not positive definite");
  for (int i = 0; i < m; ++i)
    for (int d = 0; d <= KD && i - d >= 0; ++d) Lb[(size_t)i * 4 + d] /= sc[i];

  // b = A' P' P y / sigma^2
  std::vector<double> u(n), Ptu(n);
  const double r1 = std::sqrt(1.0 - phi * phi);
  u[0] = r1 * y[0];
  for (int t = 1; t < n; ++t) u[t] = y[t] - phi * y[t - 1];
  for (int t = 0; t < n; ++t) Ptu[t] = u[t];
  Ptu[0] = r1 * u[0];
  for (int t = 0; t < n - 1; ++t) Ptu[t] -= phi * u[t + 1];

  NumericVector w(m);
  double sum_w2 = 0.0, sum_logdiag = 0.0;
  for (int i = 0; i < m; ++i) {
    double s = Ptu[i / 2] * inv_s2;
    for (int d = 1; d <= KD && i - d >= 0; ++d) s -= Lb[(size_t)i * 4 + d] * w[i - d];
    w[i] = s / Lb[(size_t)i * 4];
    sum_w2 += w[i] * w[i];
    sum_logdiag += std::log(Lb[(size_t)i * 4]);
  }
  double sum_logv = 0.0, sum_lognu = 0.0;
  for (int t = 0; t < n - 1; ++t) sum_logv += std::log(v_omega[t]);
  for (int t = 0; t < n; ++t) sum_lognu += std::log(nu2[t]);

  NumericMatrix L(m, 4);
  for (int i = 0; i < m; ++i)
    for (int d = 0; d < 4; ++d) L(i, d) = Lb[(size_t)i * 4 + d];

  return List::create(_["L"] = L, _["w"] = w,
                      _["logml"] = -0.5 * sum_logv - 0.5 * sum_lognu - sum_logdiag + 0.5 * sum_w2);
}

// [[Rcpp::export]]
NumericVector band_backsolve_cpp(NumericMatrix L, NumericVector r) {
  const int m = L.nrow();
  NumericVector x(m);
  for (int i = m - 1; i >= 0; --i) {
    double s = r[i];
    for (int d = 1; d <= KD && i + d < m; ++d) s -= L(i + d, d) * x[i + d];
    x[i] = s / L(i, 0);
  }
  return x;
}

// [[Rcpp::export]]
NumericVector tridiag_sample_cpp(NumericVector q_diag, NumericVector q_off, NumericVector b,
                                 NumericVector z) {
  const int n = q_diag.size();
  std::vector<double> l(n), e(n, 0.0), v(n);
  l[0] = std::sqrt(q_diag[0]);
  v[0] = b[0] / l[0];
  for (int i = 1; i < n; ++i) {
    e[i] = q_off[i - 1] / l[i - 1];
    const double s = q_diag[i] - e[i] * e[i];
    if (!(s > 0.0)) stop("tridiag_sample_cpp: matrix is not positive definite");
    l[i] = std::sqrt(s);
    v[i] = (b[i] - e[i] * v[i - 1]) / l[i];
  }
  NumericVector x(n);
  x[n - 1] = (v[n - 1] + z[n - 1]) / l[n - 1];
  for (int i = n - 2; i >= 0; --i) x[i] = (v[i] + z[i] - e[i + 1] * x[i + 1]) / l[i];
  return x;
}

// [[Rcpp::export]]
NumericVector ar1_filter_cpp(NumericVector e, double rho) {
  const int n = e.size();
  NumericVector h(n);
  if (n == 0) return h;
  h[0] = e[0];
  for (int t = 1; t < n; ++t) h[t] = e[t] + rho * h[t - 1];
  return h;
}
