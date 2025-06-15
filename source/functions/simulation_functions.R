library(dplyr)
library(tidyr)

simulateData <- function(
    params = list(
      data = list(
        lambda      = rep(1, 30),                 # ── now a numeric vector of length T
        T           = 30,
        date_start  = as.Date("2024-01-01"),
        D           = 20
      ),
      q_model = list(
        method        = "b_constant",
        method_params = list(b = 0.5, phi = 0.2)
      )
    )
){
  # (A) data
  lambda     <- params$data$lambda
  T          <- params$data$T
  date_start <- params$data$date_start
  D          <- params$data$D
  
  if (length(lambda) != T) {
    stop("`lambda` must be a numeric vector of length T.")
  }
  
  # (B) model for reporting proportions
  method        <- params$q_model$method
  method_params <- params$q_model$method_params
  
  simsQ_out <- generateQ(method = method,
                         params = method_params,
                         T      = T,
                         D      = D)
  
  simulation_result <- runSimulation(
    lambda     = lambda,
    T          = T,
    date_start = date_start,
    simsQ_out  = simsQ_out,
    D          = D
  )
  
  return(simulation_result)
}

runSimulation <- function(
    lambda,
    T,
    date_start,
    simsQ_out,
    D
) {
  # Generate the date sequence
  date_seq <- seq.Date(from = date_start, by = "day", length.out = T)
  
  # Initialize variables
  case_true     <- integer(T)
  case_reported <- matrix(0, nrow = T, ncol = D + 1)
  rownames(case_reported) <- as.character(date_seq)
  
  # Simulation loop
  for (tt in seq_len(T)) {
    # 1) True number of cases from given λ
    case_true[tt] <- rpois(1, lambda = lambda[tt])
    
    # 2) Get reporting proportions
    if (is.vector(simsQ_out$q)) {
      prob_temp <- simsQ_out$q
    } else {
      prob_temp <- simsQ_out$q[tt, ]
    }
    
    # 3) Single‐day reporting p’s
    p_temp <- c(prob_temp[1],
                diff(prob_temp),
                1 - prob_temp[D + 1])
    
    # 4) Multinomial split
    case_reported[tt, ] <- rmultinom(n = 1,
                                     size = case_true[tt],
                                     prob = p_temp)[1:(D + 1)]
  }
  
  # Prepare outputs
  qd_out <- if (is.vector(simsQ_out$q)) {
    simsQ_out$q[1:(D + 1)]
  } else {
    simsQ_out$q[, 1:(D + 1)]
  }
  
  case_true_mat <- matrix(case_true, ncol = 1)
  rownames(case_true_mat) <- as.character(date_seq)
  case_reported_cumulated <- t(apply(case_reported, 1, cumsum))
  
  return(list(
    case_reported           = case_reported,
    case_reported_cumulated = case_reported_cumulated,
    case_true               = case_true_mat,
    lambda                  = round(lambda, 4),
    b                       = round(simsQ_out$b, 4),
    phi                     = round(simsQ_out$phi, 4),
    q                       = round(qd_out, 4),
    date_seq                = date_seq,
    D                       = D
  ))
}

