calculate_metrics <- function(result_list, metric = "RMSE", digits_len = 2) {
  # Allowed metrics: RMSE, MAE, RMSPE, MAPE, Mean Credible Interval Width (MCIW), and Mean Coverage Rate (MCR)
  allowed_metrics <- c("RMSE", "MAE", "RMSPE", "MAPE", "MCIW", "MCR")
  if (!metric %in% allowed_metrics) {
    stop("Invalid metric. Please choose from: ", paste(allowed_metrics, collapse = ", "))
  }
  
  # Initialize an empty data frame for the output
  out_df <- data.frame(state = character(0), 
                       InfoDengue = numeric(0), 
                       GT = numeric(0), 
                       ARABM = numeric(0), 
                       ARABMwGT = numeric(0),
                       stringsAsFactors = FALSE)
  
  # Helper function to compute error metrics for RMSE, MAE, RMSPE, and MAPE
  compute_error <- function(pred, actual, metric) {
    if (metric == "RMSE") {
      return(sqrt(mean((pred - actual)^2, na.rm = TRUE)))
    } else if (metric == "MAE") {
      return(mean(abs(pred - actual), na.rm = TRUE))
    } else if (metric == "RMSPE") {
      valid <- actual != 0
      return(sqrt(mean(((pred[valid] - actual[valid]) / actual[valid])^2, na.rm = TRUE)))
    } else if (metric == "MAPE") {
      valid <- actual != 0
      return(mean(abs((pred[valid] - actual[valid]) / actual[valid]), na.rm = TRUE))
    }
  }
  
  # Function to compute Mean Credible Interval Width (MCIW)
  compute_MCIW <- function(lower, upper) {
    return(mean(upper - lower, na.rm = TRUE))
  }
  
  # Function to compute Mean Coverage Rate (MCR)
  compute_MCR <- function(lower, upper, actual) {
    covered <- (actual >= lower) & (actual <= upper)
    return(mean(covered, na.rm = TRUE))
  }
  
  # Loop through each state's data in the list
  for (df in result_list) {
    # Extract the state value; use unname() to remove any names
    state_val <- unname(unique(df$state))[1]
    
    # 'actual' holds the final reported cases
    actual <- df$sum_of_cases_final
    
    # Extract predictions and intervals for each model
    # InfoDengue model
    InfoDengue_pred  <- df$cases_est_id
    InfoDengue_lower <- df$cases_est_id_min
    InfoDengue_upper <- df$cases_est_id_max
    
    # GT model
    gt_pred  <- df$GT
    gt_lower <- df$GT_lwr
    gt_upper <- df$GT_upr
    
    # ARABM model
    arabm_pred  <- df$ARABM
    arabm_lower <- df$ARABM_lwr
    arabm_upper <- df$ARABM_upr
    
    # ARABMwGT model
    arabm_wgt_pred  <- df$ARABMwGT
    arabm_wgt_lower <- df$ARABMwGT_lwr
    arabm_wgt_upper <- df$ARABMwGT_upr
    
    # Initialize a vector to store the metrics for the current state
    model_metrics <- numeric(4)
    names(model_metrics) <- c("InfoDengue", "GT", "ARABM", "ARABMwGT")
    
    # Calculate metrics based on the selected metric type
    if (metric %in% c("RMSE", "MAE", "RMSPE", "MAPE")) {
      model_metrics["InfoDengue"]   <- compute_error(InfoDengue_pred, actual, metric)
      model_metrics["GT"]       <- compute_error(gt_pred, actual, metric)
      model_metrics["ARABM"]    <- compute_error(arabm_pred, actual, metric)
      model_metrics["ARABMwGT"] <- compute_error(arabm_wgt_pred, actual, metric)
    } else if (metric == "MCIW") {
      model_metrics["InfoDengue"]   <- compute_MCIW(InfoDengue_lower, InfoDengue_upper)
      model_metrics["GT"]       <- compute_MCIW(gt_lower, gt_upper)
      model_metrics["ARABM"]    <- compute_MCIW(arabm_lower, arabm_upper)
      model_metrics["ARABMwGT"] <- compute_MCIW(arabm_wgt_lower, arabm_wgt_upper)
    } else if (metric == "MCR") {
      model_metrics["InfoDengue"]   <- compute_MCR(InfoDengue_lower, InfoDengue_upper, actual)
      model_metrics["GT"]       <- compute_MCR(gt_lower, gt_upper, actual)
      model_metrics["ARABM"]    <- compute_MCR(arabm_lower, arabm_upper, actual)
      model_metrics["ARABMwGT"] <- compute_MCR(arabm_wgt_lower, arabm_wgt_upper, actual)
    }
    
    # Append the results for the current state to the output data frame
    out_df <- rbind(out_df, data.frame(state = state_val,
                                       InfoDengue = round(model_metrics["InfoDengue"], digits_len),
                                       GT = round(model_metrics["GT"], digits_len),
                                       ARABM = round(model_metrics["ARABM"], digits_len),
                                       ARABMwGT = round(model_metrics["ARABMwGT"], digits_len),
                                       stringsAsFactors = FALSE))
  }
  
  # Reset row names to sequential numbers from 1 to nrow(out_df)
  rownames(out_df) <- 1:nrow(out_df)
  
  return(out_df)
}

df_to_latex <- function(df, 
                        caption = "Results Table", 
                        label = "tab:results", 
                        highlight_mode = "max",   # "max" or "min"
                        highlight_n = 1           # either 1 or 2
) {
  # Check inputs
  if (!highlight_mode %in% c("max", "min")) {
    stop("highlight_mode must be either 'max' or 'min'.")
  }
  if (!(highlight_n %in% c(1, 2))) {
    stop("highlight_n must be either 1 or 2.")
  }
  
  # Define column specification: first column left-aligned, the rest right-aligned
  ncols <- ncol(df)
  col_spec <- paste0("l", paste(rep("r", ncols - 1), collapse = ""))
  
  # Build the header of the LaTeX table using booktabs style (without label here)
  header <- paste0(
    "\\begin{table}[ht]\n",
    "\\centering\n",
    "\\resizebox{\\linewidth}{!}{%\n",
    "\\begin{tabular}{", col_spec, "}\n",
    "\\toprule"
  )
  
  
  # Create the column header row
  header_row <- paste(colnames(df), collapse = " & ")
  header_row <- paste0(header_row, " \\\\ \\midrule")
  
  # Process each row to build the table body.
  # The function highlights numeric columns (assumed to be columns 2:ncols) for each row.
  body_rows <- c()
  for (i in 1:nrow(df)) {
    # Extract row (first column is assumed to be non-numeric: state)
    row_data <- df[i, ]
    state_val <- as.character(row_data[[1]])
    
    # Convert remaining columns to numeric (if not already)
    numeric_vals <- as.numeric(row_data[ , -1])
    
    # Determine ranking indices depending on highlight_mode
    if (highlight_mode == "max") {
      order_idx <- order(numeric_vals, decreasing = TRUE)
    } else {  # "min"
      order_idx <- order(numeric_vals, decreasing = FALSE)
    }
    
    # Choose which indices to highlight based on highlight_n
    highlight_idx <- order_idx[1:min(highlight_n, length(numeric_vals))]
    
    # Format numeric cells as strings; if they are highlighted, wrap them in textcolor.
    formatted_vals <- as.character(row_data[ , -1])
    
    if (length(highlight_idx) >= 1) {
      # First (best) value: highlight in red
      formatted_vals[highlight_idx[1]] <- paste0("\\textcolor{red}{", formatted_vals[highlight_idx[1]], "}")
    }
    if (highlight_n == 2 && length(highlight_idx) >= 2) {
      # Second best value: highlight in blue
      formatted_vals[highlight_idx[2]] <- paste0("\\textcolor{blue}{", formatted_vals[highlight_idx[2]], "}")
    }
    
    # Combine the state value with the formatted numeric values
    row_string <- paste(c(state_val, formatted_vals), collapse = " & ")
    row_string <- paste0(row_string, " \\\\")
    
    body_rows <- c(body_rows, row_string)
  }
  
  # Build the footer of the LaTeX table. Place the label at the bottom.
  footer <- paste0(
    "\\bottomrule\n",
    "\\end{tabular}%\n",
    "}% <-- end of resizebox\n",
    "\\label{", label, "}\n",
    "\\caption{", caption, "}\n",
    "\\end{table}"
  )
  
  # Combine header, header_row, body, and footer with newline separators.
  latex_code <- paste(c(header, header_row, body_rows, footer), collapse = "\n")
  
  return(latex_code)
}
