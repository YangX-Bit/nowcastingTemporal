library(dplyr)
library(tidyr)
library(lubridate)
library(data.table)
library(zoo)

if (!dir.exists(data_dir)) {
  stop("Specified data directory does not exist: ", data_dir)
}

# Validate date/EW format
validate_date_input <- function(start, end) {
  # Check both date and EW format validity
  date_pattern <- "\\d{4}/\\d{2}/\\d{2}"
  ew_pattern <- "\\d{4}W\\d{2}"
  
  if (!grepl(date_pattern, start) && !grepl(ew_pattern, start)) {
    stop("Invalid start date format. Use YYYY/MM/DD or YYYYW##")
  }
  
  if (!grepl(date_pattern, end) && !grepl(ew_pattern, end)) {
    stop("Invalid end date format. Use YYYY/MM/DD or YYYYW##")
  }
}

# Add state validation
validate_states <- function(data, states) {
  if (!is.null(states)) {
    available_states <- unique(data$state)
    invalid_states <- setdiff(states, available_states)
    if (length(invalid_states) > 0) {
      warning("Some specified states not found in data: ", 
              paste(invalid_states, collapse = ", "))
    }
  }
}

# Add error handling for file reading
safe_read_csv <- function(file) {
  tryCatch({
    data <- read.csv(file)
    if (nrow(data) == 0) {
      warning("Empty file: ", file)
      return(NULL)
    }
    return(data)
  }, error = function(e) {
    warning("Error reading file ", file, ": ", e$message)
    return(NULL)
  })
}

# Modify file reading step
valid_files <- lapply(valid_files, safe_read_csv)
valid_files <- Filter(Negate(is.null), valid_files)

if (length(valid_files) == 0) {
  stop("No valid files found in the specified date range")
}

process_dengue_data <- function(data_dir, start, end, resolution = "week", states = NULL) {
  # Convert to data.table for faster processing
  dengue_data <- rbindlist(lapply(valid_files, function(file) {
    data <- fread(file)
    data[, state := basename(file) %>% sub("_.*", "", .)]
    data[ew_start >= start_date & ew_start <= end_date, 
         .(state, ew_start, ew, sum_of_cases)]
  }))
}


# Allow customizable delay matrix generation
generate_delay_matrix <- function(df, max_delays = 10) {
  N <- nrow(df)
  D <- min(N, max_delays)
  
  # More robust matrix creation
  mat <- sapply(0:(D - 1), function(d) {
    if (d == 0) return(df$sum_of_cases)
    c(rep(NA, d), df$sum_of_cases[1:(N - d)])
  })
  
  colnames(mat) <- paste0("Delay_", 0:(D - 1))
  return(as.data.frame(mat))
}

# Replace existing delay matrix generation
delay_matrices <- lapply(split(filled_data, filled_data$state), 
                         generate_delay_matrix)

# Add optional logging
process_dengue_data <- function(data_dir, start, end, 
                                resolution = "week", 
                                states = NULL, 
                                verbose = FALSE) {
  if (verbose) {
    cat("Processing data from:", data_dir, "\n")
    cat("Date range:", start, "to", end, "\n")
    cat("Resolution:", resolution, "\n")
  }
  
  # Rest of the function remains the same
  
  if (verbose) {
    cat("Processed", length(delay_matrices), "state matrices\n")
  }
  
  return(delay_matrices)
}