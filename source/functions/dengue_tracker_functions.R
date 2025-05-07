library(sf)
library(geobr)
library(ggplot2)
library(geofacet)
library(tidyverse)
library(leaflet)
library(apexcharter)
library(highcharter)
library(leaflet.extras)
#library(leafgl)
library(vroom)
library(plotly)
library(viridis)
library(quantreg)
library(lubridate)
library(rvest)
library(denguetracker)

# library(glmnet)
library(aweek)
# library(forecast)

library(xtable)
file.path("..", "Repo_B")

process_data <- function(uf, last_ew_start, ew = NULL,
                         root_dir) {
  # Set directory paths relative to the root directory
  infodengue_dir <- file.path(root_dir, "infodengue")
  gtrends_dir <- file.path(root_dir, "gtrends")
  
  # List directories in the infodengue folder
  dirs <- list.dirs(infodengue_dir, full.names = TRUE)
  
  # If ew is not provided, extract the maximum ew from directory names (ignoring "city" subdirectories)
  if (is.null(ew)) {
    ew <- max(gsub(".*/(\\d+)$", "\\1", gsub("/city", "", dirs))[-1])
  }
  
  # Build the file paths using the root directories
  gt_filename <- sprintf("%s/%s/%s_trends.csv", gtrends_dir, ew, uf)
  cases_filename <- sprintf("%s/%s/%s_%s_infodengue.csv", infodengue_dir, ew, uf, last_ew_start)
  
  # Read in the cases and trends data
  cases <- read.csv(cases_filename, stringsAsFactors = FALSE)
  trends <- read.csv(gt_filename, stringsAsFactors = FALSE, skip = 2)
  
  # Clean up trends column names
  colnames(trends) <- gsub("\\.{3}.*$", "", colnames(trends))
  colnames(trends)[colnames(trends) == "Semana"] <- "Week"
  
  # Get the minimum week from trends
  min_week <- min(trends$Week)
  
  # Filter cases data based on the minimum week from trends
  cases <- cases |>
    filter(ew_start >= min_week)
  
  # Convert date columns to Date objects
  cases$ew_start <- as.Date(cases$ew_start)
  trends$Week <- as.Date(trends$Week)
  
  # Remove unwanted column from trends
  trends <- trends[, !colnames(trends) %in% "tratamento.dengue"]
  
  # Store the topics (all columns except the first) from trends
  topics <- colnames(trends)[-1]
  
  # Merge cases and trends data by matching dates
  merged_data <- merge(cases, trends, by.x = "ew_start", by.y = "Week", all = TRUE)
  
  # Convert all "<1" values to 0 and ensure numeric conversion (skip the ew_start column)
  for (col in names(merged_data)) {
    if (col != "ew_start") {
      merged_data[[col]] <- ifelse(merged_data[[col]] == "<1", 0, as.numeric(merged_data[[col]]))
    }
  }
  
  # Add the uf column to the merged data
  merged_data$uf <- uf
  
  # Return a tibble version of the merged data and the topics
  return(list(as_tibble(merged_data), topics))
}


run_model <- function(merged_data, topics, gamma, K = 4) {
  if (unique(merged_data$uf == "RR")) topics <- c("dengue")
  formula_str <- paste("sum_of_cases ~ ", paste(topics, collapse = " + "))
  best_linear_transform <- lm(
    as.formula(formula_str),
    merged_data[1:(nrow(merged_data) - K), ]
  )
  prediction <- predict(best_linear_transform, merged_data)
  
  best_linear_transform_lower <- rq(as.formula(formula_str),
                                    merged_data[1:(nrow(merged_data) - K), ],
                                    tau = (1 - gamma) / 2
  )
  
  prediction_lower <- predict(best_linear_transform_lower, merged_data)
  
  best_linear_transform_upper <- rq(as.formula(formula_str),
                                    merged_data[1:(nrow(merged_data) - K), ],
                                    tau = 1 - (1 - gamma) / 2
  )
  prediction_upper <- predict(best_linear_transform_upper, merged_data)
  
  error <-
    apply(cbind(
      prediction_lower[1:(nrow(merged_data) - K)] - merged_data$sum_of_cases[1:(nrow(merged_data) - K)],
      merged_data$sum_of_cases[1:(nrow(merged_data) - K)] - prediction_upper[1:(nrow(merged_data) - K)]
    ), 1, max)
  
  quantile_error <- quantile(error, probs = gamma, na.rm = T)
  merged_data$lwr <- pmax(prediction_lower - quantile_error, 0)
  merged_data$upr <- pmax(prediction_upper + quantile_error, 0)
  merged_data$prediction <- pmax(prediction, 0)
  
  return(merged_data)
}

generate_data <- function(ufs,
                          last_ew_start = Sys.Date() - wday(Sys.Date()) + 1,
                          ew = NULL,
                          index_of_queries = c(1,2),
                          gamma = 0.95,
                          save = F,
                          root_dir) {
  final_df <- data.frame()
  for (uf in ufs) {
    #cat(last_ew_start, ew)
    out <- process_data(uf, last_ew_start, ew = ew, 
                        root_dir = root_dir)
    data <- out[[1]]
    topics <- out[[2]][index_of_queries]
    
    K <- 4
    if(uf == "ES") K <- 15
    
    # Filter the data and use last three years to train
    date_fil <- last_ew_start %m-% years(3)
    date_fil <- date_fil %m-% weeks(K+1)
    data <- data %>% filter(ew_start >= date_fil)
    
    merged_data <- run_model(data, topics, gamma, K = K)
    if (is.null(merged_data[nrow(merged_data), "sum_of_cases"])) {
      merged_data[nrow(merged_data), "ew"] <- max(merged_data$ew, na.rm=T) + 1
    }
    final_df <- rbind(final_df, merged_data)
  }
  
  final_df <- final_df |>
    select("ew_start", "ew", "sum_of_cases", "cases_est_id", "cases_est_id_min",
           "cases_est_id_max","dengue", "sintomas.dengue", "uf", "lwr", "upr",
           "prediction")
  if (save) {
    write.csv(final_df,
              sprintf("data/model_results/model_%s.csv", last_ew_start),
              row.names = F)
  }
  final_df
}
