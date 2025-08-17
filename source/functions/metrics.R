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

################ helper function to calc results ###########
# only read last Lth N
# for RMSP, RMSPE, etc.
get_N_lastL_mean <- function(files, wks, L) {
  fit2 <- as_cmdstan_fit(files)
  vars <- paste0("N[", (wks - L + 1):wks, "]")
  s    <- fit2$summary(variables = vars)$mean  # 
  rm(fit2); gc()
  s
}

# for CRPS
get_N_lastL_draws <- function(files, wks, L) {
  fit2 <- as_cmdstan_fit(files)
  vars <- paste0("N[", (wks - L + 1):wks, "]")
  dm   <- fit2$draws(variables = vars, format = "draws_matrix")  # draws x L
  mat  <- t(dm)  # L x draws
  rm(fit2, dm); gc()
  mat
}

coverage_rate_sample <- function(y, dat, level = 0.95) {
  stopifnot(is.numeric(y), length(y) > 0)
  if (is.null(dim(dat))) stop("`dat` must be a matrix/data.frame of draws.")
  M <- as.matrix(dat)
  # Align shape so that columns correspond to time points
  if (ncol(M) == length(y)) {
    # ok
  } else if (nrow(M) == length(y)) {
    M <- t(M)
  } else {
    stop(sprintf("`dat` dims (%d x %d) not compatible with length(y) = %d",
                 nrow(M), ncol(M), length(y)))
  }
  alpha  <- (1 - level)/2
  lower  <- apply(M, 2, stats::quantile, probs = alpha,    na.rm = TRUE)
  upper  <- apply(M, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  mean(y >= lower & y <= upper)
}



################# result table for sims #####################

library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
############### tidy the table first ############### 
# ---  Filter & normalize to long-format ------------------------------------
# Returns a long table with columns: n_weeks, metric, Baseline, Proposed
filter_metric_pairs <- function(
    data,
    last_len,                  # expected to filter to a single last_len
    scenario,                  # expected to filter to a single scenario
    metrics,                   # e.g., c("RMSE","RMSPE","CRPS")
    baseline_col = "bsl",      # tolerate "bsp" -> "bsl"
    proposed_col = "prps",
    baseline_label = "Baseline",
    proposed_label = "Proposed"
){
  if (!baseline_col %in% names(data) && "bsp" %in% names(data)) {
    data <- dplyr::rename(data, !!baseline_col := .data[["bsp"]])
  }
  
  needed <- c("last_len","n_weeks","metric","scenario", baseline_col, proposed_col)
  miss   <- setdiff(needed, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  df <- data |>
    dplyr::filter(.data$last_len %in% !!last_len,
                  .data$scenario %in% !!scenario,
                  .data$metric   %in% !!metrics)
  
  if (nrow(df) == 0) stop("No rows after filtering.")
  if (dplyr::n_distinct(df$scenario) != 1)
    stop("This layout expects exactly one 'scenario' after filtering.")
  if (dplyr::n_distinct(df$last_len) != 1)
    stop("This layout expects exactly one 'last_len' after filtering.")
  
  df <- df |>
    dplyr::select(.data$n_weeks, .data$metric,
                  !!baseline_label := .data[[baseline_col]],
                  !!proposed_label := .data[[proposed_col]])
  
  # Respect the user's order of metrics in output
  df$metric <- factor(df$metric, levels = metrics)
  df
}


# --- Long -> wide (metric pairs as side-by-side columns) -------------------
# Returns a list with:
#   - wide_raw   : column keys like "RMSE__Baseline", "RMSE__Proposed", ...
#   - wide_print : human-friendly column names "RMSE: Baseline", ...
to_wide_metric_pairs <- function(
    df_long,
    metrics,
    baseline_label = "Baseline",
    proposed_label = "Proposed",
    digits = 2
){
  num_cols <- intersect(c(baseline_label, proposed_label), names(df_long))
  df_long  <- df_long |>
    dplyr::mutate(dplyr::across(dplyr::all_of(num_cols),
                                ~ round(.x, digits = digits)))
  
  long2 <- df_long |>
    tidyr::pivot_longer(cols = dplyr::all_of(c(baseline_label, proposed_label)),
                        names_to = "which_model", values_to = "value") |>
    dplyr::mutate(colkey = paste0(as.character(.data$metric), "__", .data$which_model))
  
  wide <- long2 |>
    dplyr::select(.data$n_weeks, .data$colkey, .data$value) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = "colkey", values_from = "value")
  
  pair_names <- unlist(lapply(metrics, function(m) {
    c(paste0(m, "__", baseline_label),
      paste0(m, "__", proposed_label))
  }))
  keep_cols <- c("n_weeks", intersect(pair_names, names(wide)))
  wide <- wide |>
    dplyr::select(dplyr::all_of(keep_cols)) |>
    dplyr::arrange(.data$n_weeks)
  
  pretty_names <- c("n_weeks",
                    unlist(lapply(metrics, function(m) c(
                      paste0(m, ":", baseline_label),
                      paste0(m, ":", proposed_label)
                    ))))
  wide_print <- wide
  names(wide_print) <- pretty_names
  
  list(wide_raw = wide, wide_print = wide_print)
}


# 1) Row-level highlighter: return two strings for Baseline/Proposed
highlight_pair_row <- function(b_val, p_val,
                               direction = c("lower","higher","closer"),
                               target = NA_real_,
                               color = "NavyBlue",
                               bold  = TRUE,
                               tie_highlight = TRUE) {
  direction <- match.arg(direction)
  
  # numeric for comparison; keep original string for printing
  b_num <- suppressWarnings(as.numeric(b_val))
  p_num <- suppressWarnings(as.numeric(p_val))
  b_chr <- as.character(b_val)
  p_chr <- as.character(p_val)
  
  # define loss: smaller loss = better
  loss_b <- switch(direction,
                   lower  = b_num,
                   higher = -b_num,
                   closer = abs(b_num - target))
  loss_p <- switch(direction,
                   lower  = p_num,
                   higher = -p_num,
                   closer = abs(p_num - target))
  
  # comparisons (NA-safe)
  better_b <- !is.na(loss_b) & (is.na(loss_p) | loss_b <  loss_p)
  better_p <- !is.na(loss_p) & (is.na(loss_b) | loss_p <  loss_b)
  tie_case <- !is.na(loss_b) & !is.na(loss_p) & (loss_b == loss_p)
  
  wrap <- function(x, on) {
    if (!on || is.na(x) || x == "") return(as.character(x))
    if (bold) paste0("\\textcolor{", color, "}{\\textbf{", x, "}}")
    else      paste0("\\textcolor{", color, "}{", x, "}")
  }
  
  if (tie_highlight && tie_case) {
    return(c(wrap(b_chr, TRUE), wrap(p_chr, TRUE)))
  } else {
    return(c(wrap(b_chr, better_b), wrap(p_chr, better_p)))
  }
}


# 2) Table-level highlighter: mutate your wide_print table in-place
# df_wide: columns like "n_weeks", "RMSE: Baseline", "RMSE: Proposed", ...
# Robust highlighter: works with "RMSE:Baseline" OR "RMSE: Baseline"
# Works with any number of metrics; metrics can be NULL to auto-detect.
highlight_metric_pairs_simple <- function(df_wide,
                                          metrics = NULL,          # <- allow NULL (auto-detect)
                                          directions = NULL,       # named: "lower"/"higher"/"closer"
                                          targets = NULL,          # named only for "closer" metrics
                                          baseline_label = "Baseline",
                                          proposed_label = "Proposed",
                                          color = "NavyBlue",
                                          bold  = TRUE,
                                          tie_highlight = FALSE,   # default: do NOT highlight ties
                                          sep_candidates = c(": ", ":")) {
  
  out <- df_wide
  text_cols <- setdiff(names(out), "n_weeks")
  out[text_cols] <- lapply(out[text_cols], as.character)
  
  # --- auto-detect metrics from column names if not provided ---
  if (is.null(metrics)) {
    # find columns that end with ": Baseline" (or ":Baseline")
    rx <- paste0("\\s*:\\s*", gsub("([.^$|()*+?{\\[\\]\\\\])","\\\\\\1", baseline_label), "$")
    base_cols <- grep(rx, names(out), value = TRUE)
    metric_from_col <- function(x) sub("\\s*:\\s*.*$", "", x)
    metrics <- metric_from_col(base_cols)
  }
  
  # directions: default "lower" for any metric not specified
  if (is.null(directions)) directions <- character(0)
  dir_full <- setNames(rep("lower", length(metrics)), metrics)
  if (length(directions)) {
    if (is.null(names(directions))) stop("`directions` must be named, e.g. c(RMSE='lower').")
    dir_full[names(directions)] <- directions
  }
  
  # helper: pick the first existing header among "Metric: Label" and "Metric:Label"
  resolve_col <- function(base, label) {
    for (sep in sep_candidates) {
      cand <- paste0(base, sep, label)
      if (cand %in% names(out)) return(cand)
    }
    rx2 <- paste0("^", gsub("([.^$|()*+?{\\[\\]\\\\])","\\\\\\1", base), "\\s*:\\s*",
                  gsub("([.^$|()*+?{\\[\\]\\\\])","\\\\\\1", label), "$")
    hits <- grep(rx2, names(out), value = TRUE)
    if (length(hits) > 0) return(hits[1])
    NA_character_
  }
  
  # infer target for "closer" when not provided
  infer_target <- function(metric_name, b_col, p_col) {
    tgt <- if (grepl("95", metric_name)) 0.95 else if (grepl("50", metric_name)) 0.50 else NA_real_
    b_num <- suppressWarnings(as.numeric(out[[b_col]]))
    p_num <- suppressWarnings(as.numeric(out[[p_col]]))
    medv  <- stats::median(c(b_num, p_num), na.rm = TRUE)
    if (!is.na(tgt) && is.finite(medv) && medv > 1.5) tgt <- tgt * 100  # percent scale
    tgt
  }
  
  wrap_latex <- function(x, on) {
    if (!on || is.na(x) || x == "") return(as.character(x))
    if (bold) paste0("\\textcolor{", color, "}{\\textbf{", x, "}}")
    else      paste0("\\textcolor{", color, "}{", x, "}")
  }
  
  for (m in metrics) {
    col_b <- resolve_col(m, baseline_label)
    col_p <- resolve_col(m, proposed_label)
    if (is.na(col_b) || is.na(col_p)) {
      warning(sprintf("Columns for metric '%s' not found; skipping.", m))
      next
    }
    
    dir_m <- dir_full[[m]]
    
    # target only when needed
    tgt_m <- NA_real_
    if (identical(dir_m, "closer")) {
      if (!is.null(targets) && m %in% names(targets)) {
        tgt_m <- targets[[m]]
      } else {
        tgt_m <- infer_target(m, col_b, col_p)
        if (is.na(tgt_m)) tgt_m <- 0.5
      }
    }
    
    res <- mapply(function(b, p) {
      b_num <- suppressWarnings(as.numeric(b))
      p_num <- suppressWarnings(as.numeric(p))
      loss_b <- switch(dir_m,
                       lower  = b_num,
                       higher = -b_num,
                       closer = abs(b_num - tgt_m)
      )
      loss_p <- switch(dir_m,
                       lower  = p_num,
                       higher = -p_num,
                       closer = abs(p_num - tgt_m)
      )
      better_b <- !is.na(loss_b) & (is.na(loss_p) | loss_b <  loss_p)
      better_p <- !is.na(loss_p) & (is.na(loss_b) | loss_p <  loss_b)
      tie_case <- !is.na(loss_b) & !is.na(loss_p) & (loss_b == loss_p)
      
      if (!tie_highlight && tie_case) {
        c(as.character(b), as.character(p))  # no highlight on ties
      } else if (tie_highlight && tie_case) {
        c(wrap_latex(b, TRUE), wrap_latex(p, TRUE))
      } else {
        c(wrap_latex(b, better_b), wrap_latex(p, better_p))
      }
    }, out[[col_b]], out[[col_p]])
    
    out[[col_b]] <- as.vector(res[1, ])
    out[[col_p]] <- as.vector(res[2, ])
  }
  
  out
}

# Derives header groups from df_print, so it never mismatches col count
render_latex_metric_pairs <- function(df_print,
                                      baseline_label = "Baseline",
                                      proposed_label = "Proposed",
                                      caption = NULL) {
  stopifnot(ncol(df_print) >= 3)  # n_weeks + at least one metric pair
  
  cols <- names(df_print)
  # Expect first column to be n_weeks; if not, we still proceed
  pair_cols <- cols[-1]
  if (length(pair_cols) %% 2 != 0) {
    stop("Expected an even number of columns after the first (pairs of Baseline/Proposed).")
  }
  
  # Parse metric names from headers like "RMSE:Baseline" or "RMSE: Baseline"
  metric_from_col <- function(x) sub("\\s*:\\s*.*$", "", x)
  metrics_present <- metric_from_col(pair_cols)[seq(1, length(pair_cols), by = 2)]
  
  # Build subheaders: "" for the first column, then Baseline/Proposed repeated
  subheads <- c("", rep(c(baseline_label, proposed_label), length(metrics_present)))
  if (length(subheads) != ncol(df_print)) {
    stop(sprintf("Header length mismatch: have %d columns but %d headers.",
                 ncol(df_print), length(subheads)))
  }
  
  knitr::kable(
    df_print,
    format   = "latex",
    booktabs = TRUE,
    caption  = caption,
    escape   = FALSE,         # allow \textcolor
    align    = "c",
    col.names = subheads
  ) |>
    kableExtra::add_header_above(c(" " = 1, stats::setNames(rep(2, length(metrics_present)), metrics_present))) |>
    kableExtra::kable_styling(latex_options = c("hold_position"))
}

################# combine those ######################
build_multi_scenario_table_latex <- function(
    data,
    last_len,
    scenarios = c("FR","NFRM","NFRS"),
    metrics   = c("RMSE","RMSPE","CR95","CR50","CRPS"),
    directions = NULL,                 # e.g., c(RMSE="lower", RMSPE="higher", CR95="closer", CR50="closer", CRPS="lower")
    targets    = NULL,                 # e.g., c(CR95=0.95, CR50=0.50) ; 若你的 coverage 是百分比则用 95/50
    baseline_col = "bsl",
    proposed_col = "prps",
    baseline_label = "Baseline",
    proposed_label = "Proposed",
    digits = 2,
    color  = "NavyBlue",
    caption = NULL,
    scenario_colname = "Scenario",
    n_weeks_header   = "$n_{\\text{weeks}}$"
){
  # default: unspecified metrics默认 "lower"
  if (is.null(directions)) {
    directions <- setNames(rep("lower", length(metrics)), metrics)
  } else if (is.null(names(directions))) {
    stop("`directions` must be named, e.g. c(RMSE='lower', CR95='closer').")
  }
  
  blocks <- lapply(scenarios, function(sc) {
    # 1) filter to long
    df_long <- filter_metric_pairs(
      data           = data,
      last_len       = last_len,
      scenario       = sc,
      metrics        = metrics,
      baseline_col   = baseline_col,
      proposed_col   = proposed_col,
      baseline_label = baseline_label,
      proposed_label = proposed_label
    )
    
    # 2) long -> wide (pretty names)
    ww <- to_wide_metric_pairs(
      df_long,
      metrics        = metrics,
      baseline_label = baseline_label,
      proposed_label = proposed_label,
      digits         = digits
    )
    
    # 3) highlight winners (no tie highlight)
    df_hi <- highlight_metric_pairs_simple(
      df_wide        = ww$wide_print,
      metrics        = metrics,
      directions     = directions,
      targets        = targets,
      baseline_label = baseline_label,
      proposed_label = proposed_label,
      color          = color,
      bold           = TRUE,
      tie_highlight  = FALSE
    )
    
    # 4) add scenario column and order columns: Scenario | n_weeks | pairs...
    df_hi[[scenario_colname]] <- sc
    df_hi <- df_hi[, c(ncol(df_hi), 1:(ncol(df_hi)-1)), drop = FALSE]
    df_hi[order(df_hi$n_weeks), , drop = FALSE]
  })
  
  df_all <- do.call(rbind, blocks)
  rownames(df_all) <- NULL
  
  # ---- Render LaTeX (derive headers from df_all) ----
  cols <- names(df_all)
  stopifnot(length(cols) >= 3)
  pair_cols <- cols[-(1:2)]
  if (length(pair_cols) %% 2 != 0) stop("After Scenario & n_weeks, columns must come in pairs.")
  
  metric_from_col <- function(x) sub("\\s*:\\s*.*$", "", x)
  metrics_present <- metric_from_col(pair_cols)[seq(1, length(pair_cols), by = 2)]
  
  subheads <- c(scenario_colname, n_weeks_header,
                rep(c(baseline_label, proposed_label), length(metrics_present)))
  if (length(subheads) != ncol(df_all)) {
    stop(sprintf("Header mismatch: %d columns vs %d headers.", ncol(df_all), length(subheads)))
  }
  
  knitr::kable(
    df_all,
    format   = "latex",
    booktabs = TRUE,
    caption  = caption,
    escape   = FALSE,                    
    align    = "c",
    col.names = subheads
  ) |>
    kableExtra::add_header_above(
      c(" " = 2, stats::setNames(rep(2, length(metrics_present)), metrics_present))
    ) |>
    kableExtra::kable_styling(latex_options = c("hold_position")) |>
    # Merge the Scenario column (multirow style) and draw \midrule between groups
    kableExtra::collapse_rows(columns = 1, latex_hline = "major", valign = "top")
}