data {
  // (1) Nowcasting / delay-model data
  int<lower=1>      T;                   // # of time points
  int<lower=0>      D;                   // max delay
  array[T, D+1] int<lower=0> Y;          // reported counts

  // (2) Q-model hyperparameters
  real               mean_logit_phi;
  real<lower=0>      sd_logit_phi;
  real               mean_log_b;
  real<lower=0>      sd_log_b;
  
  // Hyperparameters for alpha
  real mean_logit_alpha;
  real<lower=0> sd_logit_alpha;

  // (3) Covariates (time-varying effects)
  int<lower=1>       K;                   // number of covariates
  matrix[T, K]       x;                   // covariate matrix

  // (4) Pre-scaled Fourier basis for seasonality
  int<lower=1>       S;                   // S = 2 * H_max
  matrix[T, S]       bases_s;                 // each column sd ≈ 1
}

transformed data {
  real<lower=0> factor = (T-1) * (2*T-1) / (6.0 * T);
}

parameters {
  // (A) Delay-model (OU) parameters
  vector[T]          log_b;
  vector[T]          logit_phi;
  vector[T]          logit_alpha;  
  real               mu_log_b;
  real               mu_logit_phi;
  real               mu_logit_alpha;   
  real<lower=0>      theta_log_b;
  real<lower=0>      theta_logit_phi;
  real<lower=0>      sigma_log_b;
  real<lower=0>      sigma_logit_phi;
  real<lower=0>      sigma_logit_alpha;

  // (B) Covariate effect
  real                beta0;
  matrix[K, T]        beta_x;             // time-varying coefficients for x
  real<lower=0>       sigma_beta_x;       // prior sd for beta_x

  // (C) Seasonality coefficients
  vector[S]           beta_s;
  real<lower=0>       sigma_s;

  // (D) White noise
  vector[T]           eta;
  real<lower=0>       sigma_eta;
}

transformed parameters {
  // (1) Delay-model
  vector<lower=0>[T]            b      = exp(log_b);
  vector<lower=0,upper=1>[T]    phi    = inv_logit(logit_phi);
  vector<lower=0,upper=1>[T]   alpha   = inv_logit(logit_alpha);
  matrix[T, D+1]       q;
  matrix[T, D+1]       log_q;
  for (t in 1:T) {
    for (d in 0:D) {
      if (t + d <= T) {
        q[t, d+1]     = alpha[t] - (alpha[t] - phi[t]) * exp(-b[t] * d);
        log_q[t, d+1] = log(q[t, d+1]);
      } else {
        q[t, d+1]     = 0;
        log_q[t, d+1] = negative_infinity();
      }
    }
  }

  // (2) Seasonal component
  vector[T]            season = bases_s * beta_s;

  // (3) Intensity & likelihood pivot
  vector[T]            log_lambda;
  vector[T]            lambda;
  for (t in 1:T) {
    log_lambda[t] = beta0
                  + dot_product(beta_x[, t], x[t])
                  + season[t]
                  + eta[t];
    lambda[t]     = exp(log_lambda[t]);
  }
}

model {
  // --- (A) Delay-model priors (original) ---
  mu_log_b        ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi    ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b     ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b     ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);
  sigma_logit_alpha ~ lognormal(-3, 1); // less var

  log_b[1]        ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1]    ~ normal(mu_logit_phi, sigma_logit_phi);
  logit_alpha[1] ~ normal(mean_logit_alpha, sqrt(sd_logit_alpha^2 + sigma_logit_alpha^2 * factor));
  for (t in 2:T) {
    log_b[t]      ~ normal(
                       log_b[t-1]
                     + theta_log_b * (mu_log_b - log_b[t-1]),
                     sigma_log_b
                   );
    logit_phi[t]  ~ normal(
                       logit_phi[t-1]
                     + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]),
                     sigma_logit_phi
                   );
    logit_alpha[t] ~ normal(logit_alpha[t-1], sigma_logit_alpha);
  }

  // --- (B) Covariate effect priors (adjusted) ---
  beta0         ~ normal(0, 1);
  for (k in 1:K)
    for (t in 1:T)
      beta_x[k, t] ~ normal(0, sigma_beta_x);
  sigma_beta_x  ~ normal(0, 0.2);

  // --- (C) Seasonality priors (match original α/γ prior) ---
  for (i in 1:S)
    beta_s[i]    ~ normal(0, 0.2);
  sigma_s      ~ normal(0, 1);

  // --- (D) Noise prior (original) ---
  eta           ~ normal(0, sigma_eta);
  sigma_eta     ~ normal(0, 0.2);

  // --- Likelihood over delays ---
  for (t in 1:T)
    for (d in 0:D)
      if (t + d <= T)
        target += poisson_log_lpmf(Y[t, d+1] | log_lambda[t] + log_q[t, d+1]);
}

generated quantities {
  vector[T] N;
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}