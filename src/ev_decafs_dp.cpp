// =============================================================================
// ev_decafs_dp.cpp
// EV-DeCAFS grid dynamic programme (JISPS Algorithm 1), ported line for line from
// ev-decafs-serve/src/evdecafs_serve/models/decafs.py (ev_decafs).
//
//   Q_t(mu) = min_u { Q_{t-1}(u) + min(lambda (mu - u)^2, alpha_t)
//                     + gamma ((y_t - mu) - phi (y_{t-1} - u))^2 }
//
// on a uniform grid of n_grid levels spanning [min(y) - 2 sd(y), max(y) + 2 sd(y)]
// (sd with denominator n, as numpy.std). Ties in the minimisation go to the lowest
// grid index, as numpy.argmin.
// =============================================================================
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List ev_decafs_dp_cpp(NumericVector y, NumericVector alpha_t, double lambda_param,
                      double gamma, double phi, int n_grid) {
  const int n = y.size();
  if (n < 2) {
    return List::create(_["changepoints"] = IntegerVector(0), _["means"] = clone(y),
                        _["cost"] = 0.0);
  }

  double mean_y = 0.0;
  for (int t = 0; t < n; ++t) mean_y += y[t];
  mean_y /= n;
  double ss = 0.0, y_min = y[0], y_max = y[0];
  for (int t = 0; t < n; ++t) {
    ss += (y[t] - mean_y) * (y[t] - mean_y);
    if (y[t] < y_min) y_min = y[t];
    if (y[t] > y_max) y_max = y[t];
  }
  double std_y = std::sqrt(ss / n);
  if (std_y < 1e-10) std_y = 1.0;
  const double mu_min = y_min - 2.0 * std_y;
  const double mu_max = y_max + 2.0 * std_y;

  std::vector<double> grid(n_grid);
  for (int i = 0; i < n_grid; ++i) {
    grid[i] = (n_grid == 1) ? mu_min : mu_min + (mu_max - mu_min) * i / (n_grid - 1.0);
  }

  // Time-invariant quadratic level-change penalty lambda (mu_i - u_j)^2
  std::vector<double> quad((size_t)n_grid * n_grid);
  for (int i = 0; i < n_grid; ++i)
    for (int j = 0; j < n_grid; ++j) {
      const double d = grid[i] - grid[j];
      quad[(size_t)i * n_grid + j] = lambda_param * d * d;
    }

  std::vector<double> Q(n_grid), Q_new(n_grid);
  for (int i = 0; i < n_grid; ++i) Q[i] = gamma * (y[0] - grid[i]) * (y[0] - grid[i]);

  std::vector<int> back((size_t)n * n_grid, 0);

  for (int t = 1; t < n; ++t) {
    const double a_t = alpha_t[t];
    for (int i = 0; i < n_grid; ++i) {
      const double c_i = y[t] - grid[i] - phi * y[t - 1];
      const double* q_row = &quad[(size_t)i * n_grid];
      double best = R_PosInf;
      int best_j = 0;
      for (int j = 0; j < n_grid; ++j) {
        const double r = c_i + phi * grid[j];
        const double pen = q_row[j] < a_t ? q_row[j] : a_t;
        const double tot = Q[j] + pen + gamma * r * r;
        if (tot < best) {
          best = tot;
          best_j = j;
        }
      }
      Q_new[i] = best;
      back[(size_t)t * n_grid + i] = best_j;
    }
    Q.swap(Q_new);
  }

  int best_final = 0;
  for (int i = 1; i < n_grid; ++i)
    if (Q[i] < Q[best_final]) best_final = i;
  const double min_cost = Q[best_final];

  std::vector<int> idx(n);
  idx[n - 1] = best_final;
  for (int t = n - 2; t >= 0; --t) idx[t] = back[(size_t)(t + 1) * n_grid + idx[t + 1]];

  NumericVector mu_hat(n);
  for (int t = 0; t < n; ++t) mu_hat[t] = grid[idx[t]];

  // Change at t (0-based) when (mu_t - mu_{t-1})^2 > alpha_t / lambda; returned 1-based,
  // so a change at R index t means mu_t != mu_{t-1} (the cpd_fit convention)
  std::vector<int> cps;
  for (int t = 1; t < n; ++t) {
    const double jump = mu_hat[t] - mu_hat[t - 1];
    if (jump * jump > alpha_t[t] / (lambda_param + 1e-12)) cps.push_back(t + 1);
  }

  return List::create(_["changepoints"] = wrap(cps), _["means"] = mu_hat,
                      _["cost"] = min_cost);
}
