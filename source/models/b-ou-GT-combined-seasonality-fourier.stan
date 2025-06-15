data {
  int<lower=0> T;                      
  int<lower=0> D;                      
  array[T, D+1] int<lower=0> Y;        
  
  int<lower=1> K;                      
  matrix[T, K] X;                      

  // seasonality via Fourier basis
  int<lower=1> P;                       // seasonal cycle length
  array[T] int<lower=1,upper=P> w;      // week‐of‐year index
  int<lower=1> H;                       // number of harmonics

  // hyperpriors for reporting‐delay process
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;

  real<lower=0> sigma_hyper;           // prior SD for Fourier coeffs
}

parameters {
  real beta_0;                         
  matrix[K, T] beta_t;                 
  real<lower=0> sigma_beta;            

  vector[T] log_b;                     
  vector[T] logit_phi;                 
  real mu_log_b;                       
  real mu_logit_phi;                   
  real<lower=0> theta_log_b;           
  real<lower=0> theta_logit_phi;       
  real<lower=0> sigma_log_b;           
  real<lower=0> sigma_logit_phi;       

  // independent residual noise
  vector[T] eta;                       
  real<lower=0> sigma_eta;             

  // Fourier seasonal coefficients
  vector[H] alpha;                     // cosines
  vector[H] gamma;                     // sines
}

transformed parameters {
  vector[T] log_lambda;
  vector[T] lambda;
  vector[T] b = exp(log_b);
  vector[T] phi = inv_logit(logit_phi);

  // build seasonal effect via Fourier series
  vector[T] s_t;
  for (t in 1:T) {
    real θ = 2 * pi() * (w[t] - 1) / P;
    real season = 0;
    for (h in 1:H)
      season += alpha[h] * cos(h * θ)
              + gamma[h] * sin(h * θ);
    s_t[t] = season;
  }

  // linear predictor + seasonality
  for (t in 1:T) {
    log_lambda[t] = beta_0
                  + dot_product(beta_t[, t], X[t])
                  + eta[t]
                  + s_t[t];
    lambda[t]     = exp(log_lambda[t]);
  }
}

model {
  // intercept & regression priors
  beta_0 ~ normal(3, 1);
  for (k in 1:K)
    for (t in 1:T)
      beta_t[k, t] ~ normal(0, sigma_beta);
  sigma_beta ~ normal(0, 0.005);

  // OU priors for reporting parameters
  mu_log_b     ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b     ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b     ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);

  log_b[1]     ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1] ~ normal(mu_logit_phi, sigma_logit_phi);
  for (t in 2:T) {
    log_b[t]     ~ normal(log_b[t-1]
                           + theta_log_b * (mu_log_b     - log_b[t-1]),
                           sigma_log_b);
    logit_phi[t] ~ normal(logit_phi[t-1]
                           + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]),
                           sigma_logit_phi);
  }

  // independent residuals
  eta ~ normal(0, sigma_eta);
  sigma_eta ~ normal(0, 0.2);

  // identifiability: zero‐mean priors on Fourier coeffs
  for (h in 1:H) {
    alpha[h] ~ normal(0, sigma_hyper);
    gamma[h] ~ normal(0, sigma_hyper);
  }

  // likelihood
  for (d in 0:D)
    for (t in 1:(T-d))
      target += poisson_log_lpmf(Y[t, d+1]
                               | log_lambda[t] + 
                                 log1m_exp(-b[t] * d) +         // delay offset
                                 log(inv_logit(phi[t])));
}

generated quantities {
  vector<lower=0>[T] N;
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}
