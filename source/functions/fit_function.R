library(dplyr)
library(zoo)
library(lubridate)
library(readr)

# run_moving_window <- function(root_dir_infodengue, root_dir_GT,
#                               state, full_start, full_end,
#                               window_length = 8, step_size = 1,
#                               D, hypers, model_file,
#                               model_name = c("ARABM", "ARABMwGT"),
#                               E , bases_s, S,
#                               mode = c("fixed", "expanding"),
#                               iter_sampling = 10000, iter_warmup = 5000,
#                               chains = 3, thin = 1,
#                               posterior_draws_path
#                               ) {
#   mode <- match.arg(mode)
# 
#   # Convert full_start and full_end to Dates. Assume full_start and full_end are given in "YYYYWW"
#   start_date <- get_date(
#     year = substr(full_start, 1, 4),
#     week = substr(full_start, 5, 6)
#   )
#   end_date <- get_date(
#     year = substr(full_end, 1, 4),
#     week = substr(full_end, 5, 6)
#   )
# 
#   if (end_date < start_date) {
#     stop("End date cannot be earlier than start date.")
#   }
# 
#   # Generate a sequence of dates (weekly)
#   all_eweeks <- seq(from = full_start, to = full_end, by = 1)
#   all_dates <- seq(from = start_date, to = end_date, by = "week")
#   total_weeks <- length(all_dates)
# 
#   window_list <- list()
#   if (mode == "fixed") {
#     # Sliding windows of fixed length
#     for (i in seq(1+step_size, total_weeks - window_length + 1, by = step_size)) {
#       window_list[[length(window_list) + 1]] <- list(
#         start_date = all_dates[i],
#         end_date = all_dates[i + window_length - 1]
#       )
#     }
#   } else if (mode == "expanding") {
#     # Expanding window: start is fixed; end increases in steps
#     for (i in seq(1+step_size, total_weeks, by = step_size)) {
#       window_list[[length(window_list) + 1]] <- list(
#         start_date = start_date,
#         end_date = all_dates[i]
#       )
#     }
#   }
# 
#   # Compile Stan model
#   model <- cmdstan_model(model_file)
#   results <- list()
# 
#   # Loop over each window
#   for (win in seq_along(window_list)) {
#     win_start <- window_list[[win]]$start_date
#     win_end <- window_list[[win]]$end_date
# 
#     # Convert window start and end to EW codes (for get_infodengue_data)
#     win_start_ew <- date_to_ew(win_start)
#     win_end_ew <- date_to_ew(win_end)
# 
#     # Get delay matrix via your existing function (assumes it returns a list with element named by state)
#     window_data <- get_infodengue_data(root_dir_infodengue = root_dir_infodengue,
#                                        states = state,
#                                        start = win_start_ew,
#                                        end = win_end_ew,
#                                        D = D,
#                                        fill_missing = TRUE)
#     delay_matrix <- window_data[[state]]
#     delay_matrix[is.na(delay_matrix)] <- 0
#     T_window <- nrow(delay_matrix)
# 
# 
# 
#     # Arrange data for Stan
#     if(model_name == "ARABMwGT"){
#       gt <- read.csv(file.path(root_dir_GT,win_end_ew, paste0(state,"_trends.csv")),
#                stringsAsFactors = FALSE, skip = 2)[,c(1:3)] %>%
#         setNames(c("date", "dengue", "sintomas_dengue")) %>%
#         filter(date <= win_end & date >= win_start) %>%
#         mutate(
#           across(
#             .cols = -date,
#             .fns = ~ as.numeric(ifelse(. == "<1", "0.5", .))
#           )
#         ) %>%
#         select("dengue") %>%
#         mutate(dengue = 2*(dengue - min(dengue)) /
#                                      (max(dengue) - min(dengue)) - 1) %>%
#         as.matrix()
# 
# 
#       stan_data <- c(list(T = T_window, D = D, Y = delay_matrix),
#                      list(K = 1, x = gt, E = E,
#                           S=S, bases_s = bases_s[1:T_window, ]),
#                      hypers)
#     }else{
#       stan_data <- c(list(T = T_window, D = D, Y = delay_matrix), hypers)
#     }
# 
# 
#     # Fit the model
#     fit <- model$sample(
#       data = stan_data,
#       iter_sampling = iter_sampling,
#       iter_warmup = iter_warmup,
#       chains = chains,
#       parallel_chains = chains,
#       thin = thin,
#       refresh = 0,
#       output_dir = posterior_draws_path
#     )
# 
#     # Extract a summary (or any diagnostics you need)
#     fit_summary <- fit$summary()
# 
#     results[[win]] <- list(
#       window_start = win_start,
#       window_end = win_end,
#       fit_summary = fit_summary,
#       fit_object = fit
#     )
# 
#     cat("Finished window", win, "from", format(win_start, "%Y-%m-%d"),
#         "to", format(win_end, "%Y-%m-%d"), "\n")
#   }
# 
#   return(results)
# }

run_moving_window <- function(root_dir_infodengue, root_dir_GT,
                              state, full_start, full_end,
                              window_length = 8, step_size = 1,
                              D, hypers, model_file,
                              model_name = c("ARABM", "ARABMwGT"),
                              E , bases_s, S,
                              historical_year = 2,
                              smooth_para_gt = 3,
                              mode = c("fixed", "expanding"),
                              iter_sampling = 10000, iter_warmup = 5000,
                              chains = 3, thin = 1, refresh = 500,
                              if_show_messages = FALSE,
                              posterior_draws_path
) {
  mode <- match.arg(mode)
  
  # Convert full_start and full_end to Dates. Assume full_start and full_end are given in "YYYYWW"
  start_date <- get_date(
    year = substr(full_start, 1, 4),
    week = substr(full_start, 5, 6)
  )
  end_date <- get_date(
    year = substr(full_end, 1, 4),
    week = substr(full_end, 5, 6)
  )
  
  if (end_date < start_date) {
    stop("End date cannot be earlier than start date.")
  }
  
  # Generate a sequence of dates (weekly)
  all_eweeks <- seq(from = full_start, to = full_end, by = 1)
  all_dates <- seq(from = start_date, to = end_date, by = "week")
  total_weeks <- length(all_dates)
  
  window_list <- list()
  if (mode == "fixed") {
    # Sliding windows of fixed length
    for (i in seq(1+step_size, total_weeks - window_length + 1, by = step_size)) {
      window_list[[length(window_list) + 1]] <- list(
        start_date = all_dates[i],
        end_date = all_dates[i + window_length - 1]
      )
    }
  } else if (mode == "expanding") {
    # Expanding window: start is fixed; end increases in steps
    for (i in seq(1+step_size, total_weeks, by = step_size)) {
      window_list[[length(window_list) + 1]] <- list(
        start_date = start_date,
        end_date = all_dates[i]
      )
    }
  }
  
  # Compile Stan model
  model <- cmdstan_model(model_file)
  results <- list()
  
  weeks_per_year <- 52                        
  hist_len <- historical_year * weeks_per_year       # historial length
  # date to decide the historical data
  ealiest_date <- window_list[[1]]$start_date - weeks(historical_year * 52)
  ealiest_date_for_smooth <- ealiest_date - weeks(smooth_para_gt-1)

  # Loop over each window
  for (win in seq_along(window_list)) {
    win_start <- window_list[[win]]$start_date
    win_end <- window_list[[win]]$end_date
    
    if(win_end == "2024-06-09") {next}
    
    # Convert window start and end to EW codes (for get_infodengue_data)
    win_start_ew <- date_to_ew(win_start)
    win_end_ew <- date_to_ew(win_end)
    
    # Get delay matrix via your existing function (assumes it returns a list with element named by state)
    window_data <- get_infodengue_data(root_dir_infodengue = root_dir_infodengue,
                                       states = state,
                                       start = win_start_ew,
                                       end = win_end_ew,
                                       D = D,
                                       fill_missing = TRUE)
    delay_matrix <- window_data[[state]]
    delay_matrix[is.na(delay_matrix)] <- 0
    T_window <- nrow(delay_matrix)
    
    
    # Arrange data for Stan
    if(model_name == "ARABMwGT"){
      hist_start <- win_start - weeks(hist_len)       # hist start
      hist_end   <- win_start - weeks(1)   
      
      hist_case <- readr::read_csv(file.path(path_source_denguetracker_infodengue,
                                             win_end_ew,paste0(state,"_",as.character( win_end + weeks(1)),"_infodengue.csv")), 
                                   show_col_types = FALSE) %>%
        filter(ew_start >= hist_start & ew_start <= hist_end) %>%
        select(sum_of_cases)
      hist_case <- as.numeric(hist_case$sum_of_cases)
      
      gt_temp <- read.csv(file.path(root_dir_GT, win_end_ew, paste0(state, "_trends.csv")),
                          stringsAsFactors = FALSE, skip = 2)[, c(1:2)] %>%
        setNames(c("date", "dengue")) %>%
        filter(date >= ealiest_date_for_smooth & date <= win_end) %>%
        arrange(date) %>%
        mutate(dengue = as.numeric(ifelse(dengue == "<1", "0.5", dengue)))
      
      gt <- zoo::rollmean(as.numeric(gt_temp$dengue),
                          k = smooth_para_gt, align = "right", fill = NA)
      
      gt <- gt[-(1:(smooth_para_gt-1))] # remove NAs
      
      gt_dates <- gt_temp$date
      gt_dates <- gt_dates[-(1:(smooth_para_gt-1))]
      
      dates_needed <- seq.Date(from = hist_start, to = win_end, by = "week")     
      idx <- match(dates_needed, gt_dates)
      gt_vec <- gt[idx]
      
      k <- smooth_para_gt
      if (k > 1 && any(is.na(gt_vec[1:(k-1)]))) {                               
        gt_vec <- gt_vec[-(1:(k-1))]
        dates_needed <- dates_needed[-(1:(k-1))]
      }
      
      expected_len <- hist_len + T_window
      if (length(gt_vec) != expected_len) {                                       
        stop(sprintf("Length of gt=%d, which is different from %d (historical length) + %d (length of run-off triangle)!", 
                     length(gt_vec), hist_len, T_window))
      }
      
      stan_data <- c(list(H = length(hist_case),T = T_window, D = D, 
                          Y = delay_matrix, Y_hist = hist_case),
                     list(K = 1, x = as.matrix(gt_vec), E = E,
                          S=S, bases_s = bases_s[1:(length(hist_case)+T_window),]),
                     hypers)
    }else{
      stan_data <- c(list(T = T_window, D = D, Y = delay_matrix), hypers)
    }
    
    
    # Fit the model
    fit <- model$sample(
      data = stan_data,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      chains = chains,
      parallel_chains = chains,
      thin = thin,
      refresh = refresh,
      show_messages = if_show_messages,
      output_dir = posterior_draws_path
    )
    
    # Extract a summary (or any diagnostics you need)
    fit_summary <- fit$summary()
    
    results[[win]] <- list(
      window_start = win_start,
      window_end = win_end,
      fit_summary = fit_summary,
      fit_object = fit
    )
    
    cat("Finished window", win, "from", format(win_start, "%Y-%m-%d"),
        "to", format(win_end, "%Y-%m-%d"), "\n")
  }
  
  return(results)
}


extract_N_summary <- function(results, num_N = 1, model_name = "ARABM") {
  # num_N: number of N values to extract from the end (latest num_N elements)
  
  out_list <- list()
  for (i in seq_along(results)) {
    res <- results[[i]]
    window_start <- res$window_start  # window start date (Date object)
    window_end <- res$window_end      # window end date (Date object)
    fit <- res$fit_object        # the cmdstanr fit object
    
    # Extract summary for generated quantity "N" with desired quantiles.
    # This returns a data.frame with one row per element of N.
    n_sum <- fit$summary(variables = "N")
    total_N <- nrow(n_sum)
    if (total_N < num_N) {
      selected <- n_sum  # if there are fewer rows than requested, take all rows
    } else {
      selected <- n_sum[(total_N - num_N + 1):total_N, ]
    }
    
    ew = seq(as.numeric(date_to_ew(window_end)) - num_N + 1,
             as.numeric(date_to_ew(window_end)) , by = 1)
    # Add window information and a relative index for the selected N values.
    selected <- selected %>%
      mutate(method = model_name,
             ew_now = date_to_ew(window_end),
             ew = ew,
             !!model_name := mean,
             !!paste0(model_name, "_lwr") := q5,
             !!paste0(model_name, "_upr") := q95) %>%
      select(ew_now, ew,
             !!model_name, !!paste0(model_name, "_lwr"), !!paste0(model_name, "_upr"))
    out_list[[i]] <- selected
  }
  
  # Combine all window summaries into one data frame.
  out_df <- bind_rows(out_list)
  return(out_df)
}

#######################################################################


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
