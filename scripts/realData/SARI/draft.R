nowcasting_moving_window_with_temporal <- function(data, scoreRange, case_true = NULL,
                                     start_date = NULL, predict_length = NULL,
                                     D = 20,
                                     methods = c("baseline", "proposed"),
                                     compiled_models,
                                     S, bases_s,
                                     E,
                                     iter_sampling = 2500, iter_warmup = 1000, refresh = 500,
                                     num_chains = 4, thin = 2,suppress_output = TRUE,
                                     posterior_draws_path = file.path(path_proj, "source", "models",
                                                                      "posterior_draws"), hypers = NULL
){
  if(is.null(case_true)){
    stop("You must input true cases.")
  }
  
  if (is.null(compiled_models) || !all(methods %in% names(compiled_models))) {
    stop("You must provide compiled models matching 'methods'.")
  }
  
  # get the date
  if(is.null(start_date)){ start_date = rownames(data)[1] 
  }else {
    data <- data[rownames(data) >= start_date,]
  } 
  
  # prepare data
  data <- as.matrix(data)
  scoreRange <- as.Date(scoreRange)
  data_list <- slice_data(data, scoreRange, 
                          start_date = start_date, window_day_length = predict_length)
  scoreRange <- tail(scoreRange, length(data_list)) #remove invalid scoring date
  # result list
  model_fits <- list()
  for (i in 1:length(scoreRange)) {
    #What's "today"
    now <- scoreRange[i]
    # show the status
    cat(paste("====================\nnow=",now,
              " (",i,"/",length(scoreRange),")\n====================\n",sep=""))
    
    # prepare the data for Stan
    data_use <- data_list[[i]]
    data_trunc <- create_triangular_data(data_use, if_zero = F)
    
    # information for plot
    model_fits[["case_true"]][[i]] <- case_true[rownames(case_true) 
                                                %in% rownames(data_use), , drop = FALSE]
    model_fits[["case_reported"]][[i]] <- extract_last_valid(data_trunc)
    model_fits[["dates"]][[i]] <- as.Date(rownames(data_use))
    
    N_obs_local <- nrow(data_trunc) # num of obs
    indices_data_trunc <- find_non_na_coords(data_trunc) # coordinates for non-NAs
    data_trunc[is.na(data_trunc)] <- 0 # to avoid NAs in data
    #X_spline <- create_basis(N_obs_local, n_knots = 5) # functions to create basis
    if(nrow(data_trunc) <= D + 1){
      warning("The number of rows of the input data is smaller than number of max delay D, which might cause inaccuracy." )
    }
    # return(stan_data_trunc)
    # Fit models based on what is selected 
    for (model_name in methods) {
      compiled_model <- compiled_models[[model_name]]
      if (is.null(compiled_model)) {
        stop(paste("Model path for", model_name, "is not specified in model_paths."))
      }
      
      # Fit the Stan model
      if(model_name == "baseline"){
        stan_data_trunc <- c(list(T = N_obs_local, D = D, Y = data_trunc), hypers)
        
        sampling_code <- function() {
          compiled_model$sample(
            data = stan_data_trunc,
            iter_sampling = iter_sampling,
            iter_warmup = iter_warmup,
            chains = num_chains,
            refresh = refresh,
            thin = thin,
            output_dir = posterior_draws_path
          )
        }
      } else{
        stan_data_trunc <- c(list(T = N_obs_local, D = D, Y = data_trunc), hypers,
                             list(S = S, bases_s = bases_s[1:N_obs_local,]), E = E)
        
        sampling_code <- function(){
          compiled_model$sample(
            data = stan_data_trunc,
            iter_sampling = iter_sampling,
            iter_warmup = iter_warmup,
            chains = num_chains,
            refresh = refresh,
            thin = thin,
            output_dir = posterior_draws_path
          )
        }
      }

      
      if (suppress_output) {
        fit <- suppressWarnings(suppressMessages(sampling_code()))
      } else {
        fit <- sampling_code()
      }
      
      # Store the result
      model_fits[[model_name]][[i]] <- fit
    }
  }
  return(model_fits)
}

i = 1
model_name<- "proposed"
model_name<- "baseline"

data <- data_41002_input
scoreRange = scoreRange_SARI
case_true = case_true
start_date = first_date
D = D
methods = models_to_use
compiled_models = compiled_models
S = H_used * 2
bases_s = fourier_basis[, 1:(H_used * 2)]
E = 500
iter_sampling = 2500
iter_warmup = 1000
refresh = 0
num_chains = 4
suppress_output = T
posterior_draws_path = path_to_draws_BR
hypers = hypers

nowcasting_moving_window_with_temporal(data_41002_input,  scoreRange = scoreRange_SARI,
                                       case_true = case_true,
                                       start_date = first_date,
                                       D = D,
                                       methods = models_to_use,
                                       compiled_models = compiled_models,
                                       S = H_used * 2, bases_s = fourier_basis[, 1:(H_used * 2)],
                                       E = 500,
                                       iter_sampling = 2500, iter_warmup = 1000, refresh = 0,
                                       num_chains = 4, suppress_output = T,
                                       posterior_draws_path = path_to_draws_BR, hypers = hypers)
