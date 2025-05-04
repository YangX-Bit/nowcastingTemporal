data {
  // Data
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
  
  int<lower=1> K;                       // number of covariates
  matrix[T, K] X;                       // covariates matrix (time points x covariates)
  
  // Hyperparameters for phi and b processes
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
}

parameters {
  real beta_0;                          // intercept
  matrix[K, T] beta_t;                  // time-varying coefficients β_{k,t}
  real<lower=0> sigma_beta;             // RW-1 sd for β

  vector[T] log_b;                      // OU process state for log b
  vector[T] logit_phi;                  // OU process state for logit φ
  real mu_log_b;                        // OU long-term mean for log b
  real mu_logit_phi;                    // OU long-term mean for logit φ
  real<lower=0> theta_log_b;            // OU mean-reversion rate for log b
  real<lower=0> theta_logit_phi;        // OU mean-reversion rate for logit φ
  real<lower=0> sigma_log_b;            // OU diffusion sd for log b
  real<lower=0> sigma_logit_phi;        // OU diffusion sd for logit φ

  vector[T] eta;                        // AR(1) residual
  real<lower=-1,upper=1> rho;           // AR(1) autocorrelation
  real<lower=0> sigma_eta;              // AR(1) innovation sd
}

transformed parameters {
  vector[T] lambda;                     // Poisson mean
  vector[T] log_lambda;                 // log-linear predictor
  vector[T] b = exp(log_b);             // rate of accumulated reporting prob
  vector[T] phi = inv_logit(logit_phi); // delayed reporting prob
  matrix[T, D+1] q;                     // accumulated reporting prob
  matrix[T, D+1] log_q;                 // log(q), used as offset

  for (t in 1:T) {
    // linear predictor including eta[t]
    log_lambda[t] = beta_0
                  + dot_product(beta_t[, t], X[t, ])
                  + eta[t];
    lambda[t] = exp(log_lambda[t]);
  }
  
  for (d in 0:D) {
    for (t in 1:(T-d)) {
      q[t, d+1]     = 1 - (1 - phi[t]) * exp(-b[t] * d);
      log_q[t, d+1] = log(q[t, d+1] + 1e-12); // add tiny to avoid log(0)
    }
  }
}

model {
  // Priors
  beta_0 ~ normal(3, 1);
  mu_log_b ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);

  // RW-1 for β_{k,t}
  for (k in 1:K) {
    beta_t[k,1] ~ normal(0.03, 0.01);
    for (t in 2:T)
      beta_t[k,t] ~ normal(beta_t[k,t-1], sigma_beta);
  }
  sigma_beta ~ normal(0, 0.005);

  // OU processes for log_b and logit_phi
  log_b[1]      ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1]  ~ normal(mu_logit_phi, sigma_logit_phi);
  for (t in 2:T) {
    log_b[t]     ~ normal(log_b[t-1] + theta_log_b * (mu_log_b - log_b[t-1]),
                           sigma_log_b);
    logit_phi[t] ~ normal(logit_phi[t-1] + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]),
                           sigma_logit_phi);
  }

  // AR(1) residual for eta
  eta[1] ~ normal(0, sigma_eta);
  for (t in 2:T)
    eta[t] ~ normal(rho * eta[t-1], sigma_eta);
  rho       ~ uniform(-1, 1);
  sigma_eta ~ normal(0, 0.2);

  // Likelihood with offset log_q
  for (d in 0:D)
    for (t in 1:(T-d))
      target += poisson_log_lpmf(Y[t, d+1] | log_lambda[t] + log_q[t, d+1]);
}

generated quantities {
  vector<lower=0>[T] N;
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}
