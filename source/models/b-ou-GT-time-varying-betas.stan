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
  real beta_0;                          // intercept
  //vector[K] beta;                       // regression coefficients for covariates
  // time-varying coefficients for ALL K covariates
  matrix[K, T]  beta_t;               // β_{k,t}
  
  // RW-1 
  real<lower=0> sigma_beta;
  
  vector[T] log_b;                      // log rate of accumulated reporting probability
  vector[T] logit_phi;                  // logit of delayed reporting probability
  real mu_log_b;                        // long-term mean for log b
  real mu_logit_phi;                    // long-term mean for logit phi
  real<lower=0> theta_log_b;            // mean-reversion rate for log b
  real<lower=0> theta_logit_phi;        // mean-reversion rate for logit phi
  real<lower=0> sigma_log_b;            // diffusion coefficient for log b
  real<lower=0> sigma_logit_phi;        // diffusion coefficient for logit phi
  
  //AR(1) residual
  vector[T]   eta;    
  real<lower=-1,upper=1> rho;
  real<lower=0> sigma_eta;
}

transformed parameters {
  vector<lower=0>[T] lambda;   
  vector<lower=0>[T] b = exp(log_b);                     // rate of accumulated reporting probability
  vector<lower=0,upper=1>[T] phi = inv_logit(logit_phi); // delayed reporting probability
  matrix[T, D+1] q;                                      // accumulated reporting probability
  
  for (t in 1:T){
    //lambda[t] = exp(beta_0 + X[t] * beta);     // regression for lambda
    real log_lambda = beta_0 + dot_product( beta_t[,t] , X[t,])+ eta[t];
    lambda[t] = exp(log_lambda);
  }

  
  for (d in 0:D)
    for (t in 1:(T-d))
      q[t, d+1] = 1 - (1 - phi[t]) * exp(-b[t] * d);
}

model {
  // Priors
  // Priors for regression coefficients
  beta_0 ~ normal(0, 1);                // intercept prior
  //beta ~ normal(0, 0.2);                  // covariates prior (can adjust as needed)
  
  mu_log_b ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);
  
  //  Random-walk for every β_k(t)
  for (k in 1:K) {
    beta_t[k,1] ~ normal(0, 1);
    for (t in 2:T)
      beta_t[k,t] ~ normal( beta_t[k,t-1] , sigma_beta );
  }


  // Ornstein-Uhlenbeck processes
  log_b[1] ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1] ~ normal(mu_logit_phi, sigma_logit_phi);
  for (t in 2:T) {
    log_b[t] ~ normal(log_b[t-1] + theta_log_b * (mu_log_b- log_b[t-1]),
      sigma_log_b);
    logit_phi[t] ~ normal(logit_phi[t-1] + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]),
      sigma_logit_phi);
  }
  
  // AR(1) residual
  eta[1]  ~ normal(0, sigma_eta);
  for (t in 2:T)
    eta[t] ~ normal(rho * eta[t-1], sigma_eta);
  sigma_eta ~ normal(0,0.2);

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
