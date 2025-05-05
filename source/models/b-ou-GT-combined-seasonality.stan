data {
  int<lower=0> T;                      
  int<lower=0> D;                      
  array[T, D+1] int<lower=0> Y;        
  
  int<lower=1> K;                      
  matrix[T, K] X;                      

  // seasonality
  int<lower=1> P;                      // length of seasonal cycle (e.g. 52 weeks)
  array[T] int<lower=1,upper=P> w;           // week-of-year index for each t
  
  // Hyperparameters for phi and b
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
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

  // AR(1) residual
  vector[T] eta;                       
  real<lower=-1,upper=1> rho;          
  real<lower=0> sigma_eta;             

  // seasonal effect
  vector[P] s;                         // one seasonal effect per week
  real<lower=0> sigma_s;               // innovation sd for second-order RW
}


transformed parameters {
  vector[T] log_lambda;               
  vector[T] lambda;                   
  vector[T] b = exp(log_b);           
  vector[T] phi = inv_logit(logit_phi);
  matrix[T, D+1] q;                   
  matrix[T, D+1] log_q;               

  // map seasonal s to each t
  vector[T] s_t;
  for (t in 1:T)
    s_t[t] = s[w[t]];

  for (t in 1:T) {
    log_lambda[t] = beta_0
                  + dot_product(beta_t[, t], X[t])
                  + eta[t]
                  + s_t[t];          // add seasonality here
    lambda[t] = exp(log_lambda[t]);
  }
  
  for (d in 0:D) {
    for (t in 1:(T-d)) {
      q[t, d+1]     = 1 - (1 - phi[t]) * exp(-b[t] * d);
      log_q[t, d+1] = log(q[t, d+1] + 1e-12);
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
  
    // Independent priors for β_{k,t}
  for (k in 1:K)
    for (t in 1:T)
      beta_t[k, t] ~ normal(0, sigma_beta);
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
  
  // cyclic second-order RW prior for s[1:P]:
  //    η_w − 2 η_{w−1} + η_{w−2} ~ Normal(0, σ_s)
  // with η₀ ≡ η_P and η_{−1} ≡ η_{P−1}
  s[1] ~ normal(2*s[P] - s[P-1], sigma_s);
  s[2] ~ normal(2*s[1] - s[P],  sigma_s);
  for (p in 3:P)
    s[p] ~ normal(2*s[p-1] - s[p-2], sigma_s);
  sigma_s ~ normal(0, 1);  // half-normal on [0,∞)

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
