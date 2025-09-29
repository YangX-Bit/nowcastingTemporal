data {
  // -------------------- (A) Time scale --------------------
  int<lower=0>      H;                   // NEW: Length of historical period (totals only)
  int<lower=1>      T;                   // Current period length (with delay triangle)
  int<lower=0>      D;                   // Maximum delay

  // Historical totals (length H): corresponds to global time 1..H
  array[H] int<lower=0> Y_hist;           // NEW

  // Delay triangle for current period (size T × (D+1)): row t corresponds to global time H+t
  array[T, D+1] int<lower=0> Y;

  // -------------------- (B) Q-model hyperparameters --------------------
  real               mean_logit_phi;
  real<lower=0>      sd_logit_phi;
  real               mean_log_b;
  real<lower=0>      sd_log_b;

  // -------------------- (C) Covariates & seasonality (covering full period H+T) --------------------
  int<lower=1>       K;
  matrix[H + T, K]   x;                   // CHANGED: dimension from T to H+T
  int<lower=1>       S;
  matrix[H + T, S]   bases_s;              // CHANGED: dimension from T to H+T

  // -------------------- (D) Exposure (scalar form) --------------------
  int<lower=0>       E;                   // If needed, can be extended to vector form
}

transformed data {
  int TT = H + T;                         // NEW: total length (historical + current)
  real<lower=0> factor = (TT - 1) * (2*TT - 1) / (6.0 * TT); // kept from original (not used)
}

parameters {
  // -------------------- (1) Delay-model parameters: only for current period --------------------
  vector[T]          log_b;               // Corresponds to rows t=1..T (global index H+t)
  vector[T]          logit_phi;           // Same as above
  real               mu_log_b;
  real               mu_logit_phi;
  real<lower=0>      theta_log_b;
  real<lower=0>      theta_logit_phi;
  real<lower=0>      sigma_log_b;
  real<lower=0>      sigma_logit_phi;

  // -------------------- (2) Linear predictor --------------------
  real                beta0;
  vector[K]           beta_x;
  vector[S]           beta_s;

  // -------------------- (3) Correlated noise (covering full period) --------------------
  vector[H + T]       w;                  // CHANGED: length H+T
  real<lower=0>       sigma_w;
  real<lower=0,upper=1> rho_w;

  // -------------------- (4) White noise (covering full period) --------------------
  vector[H + T]       epsilon;            // CHANGED: length H+T
  real<lower=0>       sigma_epsilon;
}

transformed parameters {
  // -------------------- (A) Delay q: only for rows t=1..T --------------------
  vector<lower=0>[T]            b      = exp(log_b);
  vector<lower=0,upper=1>[T]    phi    = inv_logit(logit_phi);

  // Corresponding to global index g = H + t
  matrix[T, D+1]       q;
  matrix[T, D+1]       log_q;
  for (t in 1:T) {
    int g = H + t;
    // max_d: maximum usable delay for this row
    int max_d = (T - t < D) ? (T - t) : D;
    for (d in 0:D) {
      if (d <= max_d) {
        q[t, d+1]     = 1 - (1 - phi[t]) * exp(-b[t] * d);
        log_q[t, d+1] = log(q[t, d+1] + 1e-12);
      } else {
        q[t, d+1]     = 0;
        log_q[t, d+1] = -1e9;
      }
    }
  }

  // -------------------- (B) Seasonality: covering full period --------------------
  vector[H + T] season = bases_s * beta_s;

  // -------------------- (C) Intensity lambda: covering full period --------------------
  real intercept = beta0 + log(E);
  vector[H + T] linear_pred = rep_vector(intercept, H + T)
                              + x * beta_x
                              + season
                              + w
                              + epsilon;
  vector[H + T] log_lambda = linear_pred;
  vector[H + T] lambda     = exp(log_lambda);
}

model {
  // -------------------- (1) Delay OU-process priors (only for current period) --------------------
  mu_log_b        ~ normal(mean_log_b,   sd_log_b);
  mu_logit_phi    ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b     ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b     ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);

  log_b[1]        ~ normal(mu_log_b,     sigma_log_b);
  logit_phi[1]    ~ normal(mu_logit_phi, sigma_logit_phi);
  for (t in 2:T) {
    log_b[t]      ~ normal(log_b[t-1]     + theta_log_b     * (mu_log_b     - log_b[t-1]),     sigma_log_b);
    logit_phi[t]  ~ normal(logit_phi[t-1] + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]), sigma_logit_phi);
  }

  // -------------------- (2) Predictor/season/noise priors --------------------
  beta0         ~ normal(0, 1.5);
  beta_x        ~ normal(0, 0.5);
  for (i in 1:S) beta_s[i] ~ normal(0, 0.5);

  sigma_w ~ lognormal(-1.5, 1);
  rho_w   ~ lognormal(-0.8, 0.5);
  w[1]    ~ normal(0, sigma_w);
  for (g in 2:(H+T))
    w[g] ~ normal(rho_w * w[g-1], sigma_w);

  sigma_epsilon ~ lognormal(-1.5, 1);
  epsilon       ~ normal(0, sigma_epsilon);

  // -------------------- (3) Likelihood A: historical totals (no delays), g = 1..H --------------------
  if (H > 0) {
    for (h in 1:H)
      target += poisson_log_lpmf(Y_hist[h] | log_lambda[h]);
  }

  // -------------------- (4) Likelihood B: current period triangle (with delays),
  // row t=1..T mapped to global index g = H+t --------------------
  for (t in 1:T) {
    int g = H + t;
    int max_d = (T - t < D) ? (T - t) : D;
    for (d in 0:max_d) {
      target += poisson_log_lpmf(Y[t, d+1] | log_lambda[g] + log_q[t, d+1]);
    }
  }
}

generated quantities {
  // (Optional) Posterior predictive counts for full period
  vector[H + T] N;
  for (g in 1:(H+T))
    N[g] = poisson_rng(exp(log_lambda[g]));
}
