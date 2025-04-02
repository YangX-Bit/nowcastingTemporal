
run_moving_window <- function(root_dir_infodengue, root_dir_GT,
                              state, full_start, full_end,
                              window_length = 8, step_size = 1, 
                              D, hypers, model_file,
                              model_name = c("ARABM", "ARABMwGT"),
                              mode = c("fixed", "expanding"),
                              iter_sampling = 10000, iter_warmup = 5000,
                              chains = 3, thin = 1,
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
    for (i in seq(1, total_weeks - window_length + 1, by = step_size)) {
      window_list[[length(window_list) + 1]] <- list(
        start_date = all_dates[i],
        end_date = all_dates[i + window_length - 1]
        # start_ew = all_eweeks[i],
        # end_ew = all_eweeks[i + window_length - 1]
      )
    }
  } else if (mode == "expanding") {
    # Expanding window: start is fixed; end increases in steps
    for (i in seq(1, total_weeks, by = step_size)) {
      window_list[[length(window_list) + 1]] <- list(
        start_date = start_date,
        end_date = all_dates[i]
        # start_ew = start_ew,
        # end_ew = all_eweeks[i]
      )
    }
  }
  
  # Compile Stan model
  model <- cmdstan_model(model_file)
  results <- list()

  # Loop over each window
  for (win in seq_along(window_list)) {
    win_start <- window_list[[win]]$start_date
    win_end <- window_list[[win]]$end_date
    
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
      gt <- read.csv(file.path(root_dir_GT,win_end_ew, paste0(state,"_trends.csv")), 
               stringsAsFactors = FALSE, skip = 2)[,c(1:3)] %>%
        setNames(c("date", "dengue", "sintomas_dengue")) %>%
        filter(date <= win_end & date >= win_start) %>%
        mutate(
          across(
            .cols = -date,
            .fns = ~ as.numeric(ifelse(. == "<1", "0.5", .))
          )
        ) %>% 
        select("dengue", "sintomas_dengue") %>%
        as.matrix()
      stan_data <- c(list(T = T_window, D = D, Y = delay_matrix),
                     list(K = 2, X =log(gt + 0.001)), hypers)
    }else{
      stan_data <- c(list(T = T_window, D = D, Y = delay_matrix), hypers)
    }
    
    
    # Fit the model
    fit <- model$sample(
      data = stan_data,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      chains = chains,
      thin = thin,
      refresh = 0,
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


extract_N_summary <- function(results, num_N = 8, model_name = "ARABM") {
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



