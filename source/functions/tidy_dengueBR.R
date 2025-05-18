library(aweek)
library(dplyr)
library(readr)

brazil_ufs <- c(
  "AC", "AL", "AP",
  "AM",
  "BA", "CE", "DF", "ES",
  "GO",
  "MA", "MT", "MS", "MG", "PA",
  "PB", "PR", "PE",
  "PI",
  "RJ", "RN", "RS", "RO", "RR",
  "SC",
  "SP", "SE", "TO"
)

# Helper functions (for reference file processing)
is_ew_format <- function(x) {
  grepl("^\\d{6,7}$", as.character(x))
}

is_date_format <- function(x) {
  d <- try(as.Date(x), silent = TRUE)
  !inherits(d, "try-error") && !is.na(d)
}

date_to_ew <- function(date_obj) {
  aw_obj <- date2week(date_obj, week_start = "Sunday")
  year <- substr(aw_obj, 1, 4)
  weeknum <- substr(aw_obj, 7, 8)
  paste0(year, weeknum)
}

normalize_ew <- function(x) {
  if (is_ew_format(x)) {
    return(as.character(x))
  } else if (is_date_format(x)) {
    return(date_to_ew(as.Date(x)))
  } else {
    stop("Invalid `start` or `end` format. Use YYYYWW or YYYY-MM-DD.")
  }
}

get_infodengue_data <- function(
    root_dir_infodengue,
    start,
    end,
    states,
    D = 0,
    week_start = "Sunday",
    fill_missing = TRUE,
    if_last_D_cols_NA = TRUE
) {
  # Validate function inputs
  validate_inputs <- function() {
    if (!dir.exists(root_dir_infodengue)) {
      stop("Root directory does not exist: ", root_dir_infodengue)
    }
    
    if (!is.character(states) || length(states) == 0) {
      stop("`states` must be a non-empty character vector.")
    }
    
    if (!is.numeric(D) || D < 0 || D != as.integer(D)) {
      stop("`D` must be a non-negative integer.")
    }
    
    if (!is.logical(fill_missing) || length(fill_missing) != 1) {
      stop("`fill_missing` must be a single logical value.")
    }
  }
  
  validate_inputs()
  set_week_start(week_start)
  
  start_ew <- normalize_ew(start)
  end_ew <- normalize_ew(end)
  
  start_date <- get_date(
    year = substr(start_ew, 1, 4),
    week = substr(start_ew, 5, 6)
  )
  end_date <- get_date(
    year = substr(end_ew, 1, 4),
    week = substr(end_ew, 5, 6)
  )
  
  if (end_date < start_date) {
    stop("End date cannot be earlier than start date.")
  }
  
  # Generate week sequence from start to end
  all_dates <- seq(from = start_date, to = end_date, by = "week")
  all_ews <- sapply(all_dates, date_to_ew, USE.NAMES = FALSE)
  
  if (D > 0) {
    extra_dates <- seq(from = end_date + 7, by = "week", length.out = D)
    extra_ews <- sapply(extra_dates, date_to_ew, USE.NAMES = FALSE)
    all_ews_extended <- c(all_ews, extra_ews)
    all_dates_extended <- c(all_dates, extra_dates)
  } else {
    all_ews_extended <- all_ews
    all_dates_extended <- all_dates
  }
  
  N <- length(all_ews)
  total_weeks <- length(all_ews_extended)
  data_log <- list()
  result_list <- vector("list", length(states))
  names(result_list) <- states
  
  for (st_idx in seq_along(states)) {
    st <- states[st_idx]
    M <- matrix(NA, nrow = N, ncol = D + 1)
    
    for (i in seq_len(total_weeks)) {
      ew_code <- all_ews_extended[i]
      current_ew_start_date <- format(all_dates_extended[i], "%Y-%m-%d")
      delayed_date <- as.Date(current_ew_start_date) + 7
      file_name <- paste0(st, "_", format(delayed_date, "%Y-%m-%d"), "_infodengue.csv")
      folder_path <- file.path(root_dir_infodengue, ew_code)
      file_path <- file.path(folder_path, file_name)
      
      if (file.exists(file_path)) {
        df <- readr::read_csv(file_path, show_col_types = FALSE)
        
        if ("ew" %in% names(df) && "sum_of_cases" %in% names(df)) {
          for (row_i in seq_len(nrow(df))) {
            W_code <- df$ew[row_i]
            cases_val <- df$sum_of_cases[row_i]
            W_idx <- match(W_code, all_ews)
            
            if (!is.na(W_idx)) {
              col <- i - W_idx + 1
              if (col >= 1 && col <= (D + 1)) {
                M[W_idx, col] <- cases_val
              }
            }
          }
        }
        
        data_log[[length(data_log) + 1]] <- list(
          state = st,
          ew = ew_code,
          file = file_name,
          status = "Found"
        )
      } else {
        data_log[[length(data_log) + 1]] <- list(
          state = st,
          ew = ew_code,
          file = file_name,
          status = "Missing"
        )
      }
    }
    
    # Fill missing values from left to right
    if (fill_missing) {
      for (r in seq_len(N)) {
        for (c in 2:(D + 1)) {
          if (is.na(M[r, c])) {
            M[r, c] <- M[r, c - 1]
          }
        }
      }
    }
    
    # Set last X rows of delayX columns to NA
    if(if_last_D_cols_NA){
      for (X in seq_len(D)) {
        start_index <- max(1, N - X + 1)
        M[start_index:N, X + 1] <- NA
      }
    }

    rownames(M) <- format(all_dates, "%Y-%m-%d")
    colnames(M) <- paste0("delay", 0:D)
    result_list[[st]] <- M
  }
  
  structure(result_list, data_log = data_log)
}

# get the baseline dataframe
get_baseline_data <- function(root_dir, states, start, end, last_n = 8, 
                              time_unit = c("week","day"),
                              all_dates = NULL) {
  
  # Process start and end dates
  start_ew <- normalize_ew(start)
  end_ew <- normalize_ew(end)
  
  start_date <- get_date(year = substr(start_ew, 1, 4), week = substr(start_ew, 5, 6))
  end_date <- get_date(year = substr(end_ew, 1, 4), week = substr(end_ew, 5, 6))
  
  if (end_date < start_date) {
    stop("End date cannot be earlier than start date.")
  }
  
  if(time_unit == "week"){
    time_gap = 7
  }else if(time_unit == "day"){
    time_gap = 1
  }else{
    stop("time_unit has to be 'week' or 'day'!")
  }
  
  # Generate weekly sequence of dates
  if (is.null(all_dates)) {
    all_dates <- seq(from = start_date, to = end_date, by = time_gap)
  }
  all_ews <- sapply(all_dates, date_to_ew, USE.NAMES = FALSE)
  
  results_list <- list()
  
  # Loop over each state and week
  for (st in states) {
    for (i_week in seq_along(all_ews)) {
      ew_now <- all_ews[i_week]
      current_ew_start_date <- all_dates[i_week]
      delayed_date <- as.Date(current_ew_start_date) + 7
      file_name <- paste0(st, "_", format(delayed_date, "%Y-%m-%d"), "_infodengue.csv")
      folder_path <- file.path(root_dir, "infodengue", ew_now)
      file_path <- file.path(folder_path, file_name)
      
      if (file.exists(file_path)) {
        df <- read_csv(file_path, show_col_types = FALSE)
        if(ew_now == "202409"){
          df$ew <- str_replace(substr(as.aweek(df$ew_start, week_start = "Sunday"), 1, 8),
                               "-W", "")
        }
        required_cols <- c("ew_start", "ew", "sum_of_cases", 
                           "cases_est_id", "cases_est_id_min", "cases_est_id_max")
        if (all(required_cols %in% names(df))) {
          n_rows <- nrow(df)
          n_to_take <- ifelse(n_rows >= last_n, last_n, n_rows)
          df_last <- tail(df, n_to_take)
          # prediction from linear model with GT
          df_GT_pred <- tail(generate_data(states, last_ew_start = delayed_date, ew = ew_now, save=F,
                                           root_dir = root_dir), last_n)
          
          df_last <- df_last %>%
            mutate(state = st, ew_now = ew_now,
                   file_date = format(current_ew_start_date, "%Y-%m-%d"),
                   GT = df_GT_pred$prediction, 
                   GT_lwr = df_GT_pred$lwr, GT_upr = df_GT_pred$upr) %>%
            select(state, file_date, ew_now, all_of(required_cols), 
                   GT, GT_lwr, GT_upr)
          
          results_list[[length(results_list) + 1]] <- df_last
        }
      }
    }
  }
  
  results_df <- bind_rows(results_list)
  return(results_df)
}

# choose ref data and add the "actual cases" to compare
add_actual_cases <- function(df, root_dir_infodengue, states, reference) {
  # Process the reference week
  ref_ew <- normalize_ew(reference)
  ref_start_date <- get_date(
    year = substr(ref_ew, 1, 4),
    week = substr(ref_ew, 5, 6)
  )
  ref_delayed_date <- as.Date(ref_start_date) + 7
  
  ref_list <- list()
  # Loop over each state to read its reference file
  for (st in states) {
    ref_file_name <- paste0(st, "_", format(ref_delayed_date, "%Y-%m-%d"), "_infodengue.csv")
    ref_folder <- file.path(root_dir_infodengue, ref_ew)
    ref_file_path <- file.path(ref_folder, ref_file_name)
    
    if (file.exists(ref_file_path)) {
      ref_df <- read_csv(ref_file_path, show_col_types = FALSE)
      if (all(c("ew", "sum_of_cases") %in% names(ref_df))) {
        ref_df <- ref_df %>%
          select(ew, sum_of_cases) %>%
          rename(sum_of_cases_final = sum_of_cases) %>%
          mutate(state = st)
        ref_list[[length(ref_list) + 1]] <- ref_df
      }
    }
  }
  
  if (length(ref_list) > 0) {
    ref_all <- bind_rows(ref_list)
    # Merge the reference data with the original data frame by state and ew
    df <- df %>% left_join(ref_all, by = c("state", "ew"))
  } else {
    df <- df %>% mutate(sum_of_cases_final = NA_real_)
  }
  
  return(df)
}

