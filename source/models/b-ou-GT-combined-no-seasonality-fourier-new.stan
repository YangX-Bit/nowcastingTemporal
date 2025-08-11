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
  // real mean_logit_alpha;
  // real<lower=0> sd_logit_alpha;

  // (3) Predictor
  int<lower=1>       K;                   // number of covariates
  matrix[T, K]       x;                   // covariate matrix

  // (4) Pre-scaled Fourier basis for seasonality
  //int<lower=1>       S;                   // S = 2 * H_max
  //matrix[T, S]       bases_s;                 // each column sd ≈ 1
  
  // (5) Population at risk
  int<lower=0> E;
}

transformed data {
  real<lower=0> factor = (T-1) * (2*T-1) / (6.0 * T);
}

parameters {
  // (1) Delay-model (OU) parameters
  vector[T]          log_b;
  vector[T]          logit_phi;
  // vector[T]          logit_alpha;  
  real               mu_log_b;
  real               mu_logit_phi;
  // real               mu_logit_alpha;   
  real<lower=0>      theta_log_b;
  real<lower=0>      theta_logit_phi;
  real<lower=0>      sigma_log_b;
  real<lower=0>      sigma_logit_phi;
  // real<lower=0>      sigma_logit_alpha;

  // (2) Predictor
  real                beta0;
  vector[K]           beta_x;             // coefficients for x

  // (3) Seasonality coefficients
  //vector[S]           beta_s;
  
  // (4) Dependent noise
  vector[T]           w;
  real<lower=0>       sigma_w;
  real<lower=0,upper=1> rho_w;

  // (5) White noise
  vector[T]           epsilon;
  real<lower=0>       sigma_epsilon;
}

transformed parameters {
  // (1) Delay-model
  vector<lower=0>[T]            b      = exp(log_b);
  vector<lower=0,upper=1>[T]    phi    = inv_logit(logit_phi);
  matrix[T, D+1]       q;
  matrix[T, D+1]       log_q;
  for (d in 0:D) {
    for (t in 1:(T-d)) {
      q[t, d+1]     = 1 - (1 - phi[t]) * exp(-b[t] * d);
      log_q[t, d+1] = log(q[t, d+1] + 1e-12); // add tiny to avoid log(0)
    }
  }

  // (2) Seasonal component
  //vector[T]            season = bases_s * beta_s;

  // (3) Risk & Intensity
  // vector[T]            log_risk;
  // vector[T]            risk;
  // vector[T] log_lambda;
  // vector[T] lambda;
  // for (t in 1:T) {
  //   log_risk[t] = beta0
  //                 + dot_product(beta_x, x[t])
  //                 //+ beta_x * x[t]
  //                 + season[t]
  //                 + w[t]
  //                 + epsilon[t];
  //   log_lambda[t] = log(E) + log_risk[t];
  //   lambda[t]     = exp(log_lambda[t]);
  // }
  
  real intercept = beta0 + log(E);
  vector[T] linear_pred = rep_vector(intercept, T)
                   + x * beta_x 
                   //+ season
                   + w
                   + epsilon;

  vector[T] log_lambda = linear_pred;
  vector[T] lambda     = exp(log_lambda);
}


model {
  // --- (1) Delay-model priors (original) ---
  mu_log_b        ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi    ~ normal(mean_logit_phi, sd_logit_phi);
  theta_log_b     ~ lognormal(0, 1);
  theta_logit_phi ~ lognormal(0, 1);
  sigma_log_b     ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);
  // sigma_logit_alpha ~ lognormal(-3, 1); // less var

  log_b[1]        ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1]    ~ normal(mu_logit_phi, sigma_logit_phi);
  // logit_alpha[1] ~ normal(mean_logit_alpha, sqrt(sd_logit_alpha^2 + sigma_logit_alpha^2 * factor));
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
    // logit_alpha[t] ~ normal(logit_alpha[t-1], sigma_logit_alpha);
  }

  // --- (2) Predictor coefficient ---
  beta0         ~ normal(0, 0.1);
  
  beta_x        ~ normal(0, 0.5);
      
  // --- (3) Seasonality priors (match original α/γ prior) ---
  // for (i in 1:S)
  //   beta_s[i]    ~ normal(0, 0.5);

  // --- (4) Dependent noise prior ---
  sigma_w ~ lognormal(-1.5, 1);
  rho_w ~ lognormal(-0.3, 0.5);
  w[1] ~ normal(0, sigma_w);
  for (j in 2:T)
    w[j] ~ normal(rho_w * w[j-1],sigma_w);
  
  // --- (5) Independent noise prior ---
  sigma_epsilon     ~ lognormal(-1.5, 1);
  epsilon           ~ normal(0, sigma_epsilon);


  // --- Likelihood over delays ---
  for (d in 0:D)
    for (t in 1:(T-d))
        target += poisson_log_lpmf(Y[t, d+1] | log_lambda[t] + log_q[t, d+1]);
}

generated quantities {
  vector[T] N;
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}