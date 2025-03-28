data {
  // Data
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
  
  int<lower=1> K;                       // number of covariates
  matrix[T, K] X;                       // covariates matrix (time points x covariates)
  
  // Hyperparameters
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
}

parameters {
  vector[K] beta;                       // regression coefficients for covariates
  vector[K] beta_quad;    
  real beta_0;                          // intercept
  
  vector[T] log_b;                      // log rate of accumulated reporting probability
  vector[T] logit_phi;                  // logit of delayed reporting probability
  real mu_log_b;                        // long-term mean for log b
  real mu_logit_phi;                    // long-term mean for logit phi
  real<lower=0> theta_log_b;            // mean-reversion rate for log b
  real<lower=0> theta_logit_phi;        // mean-reversion rate for logit phi
  real<lower=0> sigma_log_b;            // diffusion coefficient for log b
  real<lower=0> sigma_logit_phi;        // diffusion coefficient for logit phi
}

transformed parameters {
  vector<lower=0>[T] lambda;   
  vector<lower=0>[T] b = exp(log_b);                     // rate of accumulated reporting probability
  vector<lower=0,upper=1>[T] phi = inv_logit(logit_phi); // delayed reporting probability
  matrix[T, D+1] q;                                      // accumulated reporting probability
  
  for (t in 1:T)
    lambda[t] = exp(beta_0 + X[t] * beta + square(X[t]) * beta_quad);     // regression for lambda
  
  for (d in 0:D)
    for (t in 1:(T-d))
      q[t, d+1] = 1 - (1 - phi[t]) * exp(-b[t] * d);
}

model {
  // Priors
  // Priors for regression coefficients
  beta_0 ~ normal(0, 1);                // intercept prior
  beta ~ normal(0, 0.2);                  // covariates prior (can adjust as needed)
  beta_quad ~ normal(0, 0.2); 
  
  mu_log_b ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);

  // Ornstein-Uhlenbeck processes
  log_b[1] ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1] ~ normal(mu_logit_phi, sigma_logit_phi);
  for (t in 2:T) {
    log_b[t] ~ normal(log_b[t-1] + theta_log_b * (mu_log_b- log_b[t-1]),
      sigma_log_b);
    logit_phi[t] ~ normal(logit_phi[t-1] + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]),
      sigma_logit_phi);
  }

  // Likelihood
  for (d in 0:D)
    for (t in 1:(T-d))
      Y[t, d+1] ~ poisson(lambda[t] * q[t, d+1]);
}

generated quantities {
  vector<lower=0>[T] N;                 // number of cases
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}
