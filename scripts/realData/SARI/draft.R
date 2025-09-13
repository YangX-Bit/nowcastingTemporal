nowcasting_moving_window_with_temporal <- function(data, scoreRange, case_true = NULL,
                                     start_date = NULL, start_date_triangular = NULL, predict_length = NULL,
                                     D = 20,
                                     methods = c("baseline", "proposed"),
                                     compiled_models,
                                     S, bases_s,
                                     E,
                                     iter_sampling = 2500, iter_warmup = 1000, refresh = 500,
                                     num_chains = 4, thin = 1,suppress_output = TRUE,
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
    
    # for historical
    nrow_recent <- as.numeric((now - start_date_triangular)/7)
    data_trunc_recent <- tail(data_trunc, nrow_recent)
    N_obs_local_recent <- nrow(data_trunc_recent)
    
    data_trunc_hist <- head(data_trunc, N_obs_local - nrow_recent)
    Y_hist <- data_trunc_hist[, c(ncol(data_trunc_hist))]
    
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
            parallel_chains = num_chains,
            refresh = refresh,
            thin = thin,
            output_dir = posterior_draws_path
          )
        }
      }else if(model_name == "baseline2"){
        stan_data_trunc <- c(list(T = N_obs_local_recent, D = D, Y = data_trunc_recent), hypers)
        
        sampling_code <- function(){
          compiled_model$sample(
            data = stan_data_trunc,
            iter_sampling = iter_sampling,
            iter_warmup = iter_warmup,
            chains = num_chains,
            parallel_chains = num_chains,
            refresh = refresh,
            thin = thin,
            output_dir = posterior_draws_path
          )
        }
      }else{
        stan_data_trunc <- c(list(H = N_obs_local - N_obs_local_recent, 
                                  T = N_obs_local_recent, D = D, 
                                  Y = data_trunc_recent, Y_hist = Y_hist), hypers,
                             list(S = S, bases_s = bases_s[1:N_obs_local,]), E = E)
        
        sampling_code <- function(){
          compiled_model$sample(
            data = stan_data_trunc,
            iter_sampling = iter_sampling,
            iter_warmup = iter_warmup,
            chains = num_chains,
            parallel_chains = num_chains,
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

data <- SARI_collapsed_input
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
thin = 1
refresh = 0
num_chains = 4
suppress_output = T
posterior_draws_path = path_to_draws_BR
hypers = hypers

nowcasts_table <- function(results_list, 
                           D = NULL,
                           report_unit = "week",
                           methods = c("fixed_q", "fixed_b", "linear_b", "ou_b"),
                           N_num = NULL,
                           replicate_id = NA_integer_, alpha = 0.05) {
  # Basic checks
  if (is.null(D)) stop("Parameter 'D' must be provided.")
  if (!report_unit %in% c("week", "day")) {
    stop("report_unit must be 'week' or 'day'.")
  }
  
  library(lubridate)
  library(dplyr)
  
  # Decide factor for date shifting
  factor_loc <- if (report_unit == "week") 7 else 1
  nowcasts_out <- list()
  n_runs <- length(N_num)  # how many sets of data we have
  for (i in 1:n_runs) {
    # Extract data
    case_true     <- results_list[["case_true"]][[i]]
    case_reported <- results_list[["case_reported"]][[i]]
    dates         <- results_list[["dates"]][[i]]
    
    # Basic date references
    now      <- as.Date(dplyr::last(dates))
    earliest <- as.Date(dplyr::first(dates))
    last_date_for_delay <- now - days(D * factor_loc)
    
    # Initialize a data frame for storing nowcasts
    nowcasts_df <- data.frame(
      date          = dates,
      case_true     = case_true,
      case_reported = case_reported,
      row.names     = seq_along(dates)
    )
    
    # Store these references as columns (the same for all rows in this i)
    nowcasts_df$now                <- now
    nowcasts_df$earliest           <- earliest
    nowcasts_df$last_date_for_delay <- last_date_for_delay
    N_all <- nrow(nowcasts_df)
    # cut

    nowcasts_df <- tail(nowcasts_df, N_num[i])
    
    if (!is.na(replicate_id)) {
      nowcasts_df$replicate_id       <- replicate_id # for replicate
    }
    
    lower_p <- alpha / 2
    upper_p <- 1 - alpha / 2
    
    last_vars <- paste0("N[", (N_all - N_num[i] + 1):N_all, "]")

    # Dynamically add model results
    for (model_name in methods) {
      # model_name might be "fixed_q", "fixed_b", ...
      if(model_name != "baseline2"){
        samples <- results_list[[model_name]][[i]]$draws(variables = "N", format = "draws_matrix")[, last_vars]
      }else {
        samples <- results_list[[model_name]][[i]]$draws(variables = "N", format = "draws_matrix")
      }
      nowcasts_df[[paste0("mean_", model_name)]]  <- apply(samples, 2, mean)
      nowcasts_df[[paste0("lower_", model_name)]] <- apply(samples, 2, quantile, probs = lower_p, na.rm = TRUE)# if this na.rm = TRUE is okay?
      nowcasts_df[[paste0("upper_", model_name)]] <- apply(samples, 2, quantile, probs = upper_p, na.rm = TRUE)
    }
    # Save to output list
    nowcasts_out[[i]] <- nowcasts_df
  }
  
  return(nowcasts_out)
}

model_name <- "proposed"

results_list <- out_SARI

nowcasts_table(out_SARI, D = D, report_unit = "week", 
               methods = models_to_use)
