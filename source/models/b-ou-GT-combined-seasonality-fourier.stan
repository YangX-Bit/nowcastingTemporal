data {
  // Data
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
  // Hyperparameters for qtd
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
  
  // For lambda
  // covariates x
  int<lower=1> K;                      
  matrix[T, K] x;    
  
  // seasonality
  int<lower=1> P; // number of time‐units in one full seasonal cycle (e.g., weeks per year)
  int<lower=1> H; //  number of harmonics (Fourier components) 
  array[T] int<lower=1, upper=P> week; 
}

parameters {
  //vector<lower=0>[T] lambda;            // expected number of cases
  vector[T] log_b;                      // log rate of accumulated reporting probability
  vector[T] logit_phi;                  // logit of delayed reporting probability
  real mu_log_b;                        // long-term mean for log b
  real mu_logit_phi;                    // long-term mean for logit phi
  real<lower=0> theta_log_b;            // mean-reversion rate for log b
  real<lower=0> theta_logit_phi;        // mean-reversion rate for logit phi
  real<lower=0> sigma_log_b;            // diffusion coefficient for log b
  real<lower=0> sigma_logit_phi;        // diffusion coefficient for logit phi
  
  // betas
  real beta_0; 
  matrix[K, T] beta_t; 
  real<lower=0> sigma_beta;
  
  // white noise
  vector[T] eta;                       
  real<lower=0> sigma_eta;  
  
  // seasonality
  vector[H] alpha;
  vector[H] gamma;
}

transformed parameters {
  vector<lower=0>[T] b = exp(log_b);                     // rate of accumulated reporting probability
  vector<lower=0,upper=1>[T] phi = inv_logit(logit_phi); // delayed reporting probability
  matrix[T, D+1] q;                                     // accumulated reporting probability
  matrix[T, D+1] log_q;                                      
  for (d in 0:D){
    for (t in 1:(T-d)){
      q[t, d+1] = 1 - (1 - phi[t]) * exp(-b[t] * d);
      log_q[t, d+1] = log( q[t, d+1] );
    }
  }
  
  vector[T] s;
  vector[T] log_lambda;
  vector[T] lambda;   
  for (t in 1:T) {
    // seasonality
    real angle = 2*pi() * (week[t] - 1) / P;
    real sum_h = 0;
    for( h in 1:H){
      sum_h += alpha[h] * cos(h * angle)
      + gamma[h] * sin(h * angle);
    }
    s[t] = sum_h;
    // final lambda
    log_lambda[t] = beta_0
                  + dot_product(beta_t[, t], x[t])
                  + s[t]
                  + eta[t];
    lambda[t] = exp(log_lambda[t]);
  }
}

model {
  // Priors
  lambda ~ lognormal(0, 2.5);
  mu_log_b ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);
  
  // Independent priors for β_{k,t}
  for (k in 1:K)
    for (t in 1:T)
      beta_t[k, t] ~ normal(0.5, sigma_beta);
  sigma_beta ~ normal(0, 0.2);
  
  // Seasonality
  alpha ~ normal(0,0.2);
  gamma ~ normal(0,0.2);
  
  // noise
  sigma_eta ~ normal(0,0.2);

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
      target += poisson_log_lpmf(Y[t, d+1] | log_lambda[t] + log_q[t, d+1]);
}

generated quantities {
  vector<lower=0>[T] N;                 // number of cases
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}

