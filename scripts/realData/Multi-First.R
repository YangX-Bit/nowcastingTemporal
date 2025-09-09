## =========================
## (A) 公共工具函数
## =========================
library(dplyr)
library(purrr)
library(tidyr)
library(readr)

# 更稳健的误差函数（避免 actual = 0 时除零）
safe_rmspe <- function(pred, actual, eps = 1e-8) {
  sqrt(mean(((pred - actual) / pmax(actual, eps))^2, na.rm = TRUE)) * 100
}
safe_mape <- function(pred, actual, eps = 1e-8) {
  mean(abs((pred - actual) / pmax(actual, eps)), na.rm = TRUE) * 100
}
rmse_fun <- function(pred, actual) sqrt(mean((pred - actual)^2, na.rm = TRUE))
mae_fun  <- function(pred, actual) mean(abs(pred - actual), na.rm = TRUE)

# 覆盖率
coverage_rate <- function(truth, lwr, upr) {
  mean(truth >= lwr & truth <= upr, na.rm = TRUE)
}

# 基于正态近似的 CRPS（用 95% 区间反推 σ，再计算 CRPS(N(μ,σ), y)）
crps_normal_vec <- function(y, mu, lwr95, upr95) {
  z975 <- qnorm(0.975)
  sigma <- (upr95 - lwr95) / (2 * z975)
  # 避免 sigma = 0
  sigma <- pmax(sigma, 1e-8)
  z <- (y - mu) / sigma
  # φ 和 Φ
  phi <- dnorm(z)
  Phi <- pnorm(z)
  # CRPS 公式：σ [ 1/√π - 2 φ(z) - z(2Φ(z)-1) ]
  crps <- sigma * (1 / sqrt(pi) - 2 * phi - z * (2 * Phi - 1))
  mean(crps, na.rm = TRUE)
}

# 由 95% 区间近似出 50% 区间（正态：μ ± z0.75 * σ）
make_50_ci_from_95 <- function(mu, lwr95, upr95) {
  z975 <- qnorm(0.975)
  z075 <- qnorm(0.75)  # 0.6744898
  sigma <- (upr95 - lwr95) / (2 * z975)
  sigma <- pmax(sigma, 1e-8)
  lwr50 <- mu - z075 * sigma
  upr50 <- mu + z075 * sigma
  list(lwr50 = lwr50, upr50 = upr50)
}

generate_fourier_basis <- function(weeks, H_max = 3, period = 52) {
  if (inherits(weeks, "Date")) {
    weeks <- seq_along(weeks)
  }
  
  angle <- 2 * pi * (weeks - 1) / period
  
  fourier_list <- lapply(1:H_max, function(h) {
    cbind(
      cos = cos(h * angle),
      sin = sin(h * angle)
    )
  })
  
  basis <- do.call(cbind, fourier_list)
  
  # colnames
  colnames(basis) <- unlist(lapply(1:H_max, function(h) {
    c(paste0("cos_h", h), paste0("sin_h", h))
  }))
  
  return(basis)
}

weeks_per_yr <- 52
t <- 1:(weeks_per_yr * 3)
fourier_basis <- generate_fourier_basis(t, H_max = 10, period = weeks_per_yr)
H_used <- 3

## =========================
## (B) 路径 & 全局设置（保持你现在的变量）
## =========================
# 这些路径/参数你已在上文定义：path_proj, path_source_nowcasting, path_source, 
# path_source_denguetracker_data, path_source_denguetracker_infodengue, path_source_denguetracker_GT,
# ew_start, ew_end, delay_D, step_size, iter_sampling, iter_warmup, chain_num, H_used, fourier_basis, etc.

# 拆分保存目录
fit_dir <- file.path(path_proj, "outputs", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

# 你希望批量跑的州 & 对应的 E（与州一一对应）
all_states <- c(
  "AC", "AL", "AP", "AM",
  "BA", "CE", "DF", "ES",
  "GO", "MA", "MT", "MS",
  "MG", "PA", "PB", "PR",
  "PE", "PI", "RJ", "RN",
  "RS", "RO", "RR", "SC",
  "SP", "SE", "TO"
)       # 示例；请改成你的 20+ 个州
E_vec      <- c(
  "AC", "AL", "AP", "AM",
  "BA", "CE", "DF", "ES",
  "GO", "MA", "MT", "MS",
  "MG", "PA", "PB", "PR",
  "PE", "PI", "RJ", "RN",
  "RS", "RO", "RR", "SC",
  "SP", "SE", "TO"
)      # 与州一一对应
include_GT <- FALSE                      # 是否在 metrics 表中包含 GT 指标

delay_D <- 15
iter_sampling <- 2500
iter_warmup <- 1000
window_length <- 10
thin <- 1
step_size <- 5
chain_num <- 4

ew_start <- "202410"
ew_end <- "202430"


## =========================
## (C) Loop 1: 仅运行 MCMC 并保存（不做指标）
## =========================
for (i in seq_along(all_states)) {
  state_i <- all_states[i]
  E_val   <- E_vec[i]
  cat(sprintf(">>> [MCMC %d/%d] State = %s ...\n", i, length(all_states), state_i))
  
  # baseline（含真实值）
  baseline_i <- get_baseline_data(
    root_dir = path_source_denguetracker_data,
    states   = state_i,
    start    = ew_start,
    end      = ew_end,
    time_unit = "week",
    last_n   = 10
  ) %>%
    add_actual_cases(root_dir_infodengue = path_source_denguetracker_infodengue,
                     state = state_i,
                     202511)
  saveRDS(baseline_i, file.path(fit_dir, sprintf("%s_baseline.rds", state_i)))
  
  # 模型 1：ARABM（不含 GT）
  results_i <- run_moving_window(
    root_dir_infodengue = path_source_denguetracker_infodengue,
    root_dir_GT         = path_source_denguetracker_GT,
    state               = state_i,
    full_start          = ew_start,
    full_end            = ew_end,
    step_size           = step_size,
    D                   = delay_D,
    hypers              = hypers_q(phi_ref = 0.2, D_ref = 15, type = "exponential",
                                   alpha_phi = 1.4, sd_log_b = 1, delay_seq = 0:15),
    model_file          = file.path(path_source_nowcasting, "models", "b-ou.stan"),
    model_name          = "ARABM",
    mode                = "expanding",
    iter_sampling       = iter_sampling,
    iter_warmup         = iter_warmup,
    chains              = chain_num,
    posterior_draws_path = file.path(path_proj, "draws", "length8")
  )
  saveRDS(results_i, file.path(fit_dir, sprintf("%s_ARABM_fit.rds", state_i)))
  
  # 模型 2：ARABMwGT（含 GT + seasonality + E）
  results_wGT_i <- run_moving_window(
    root_dir_infodengue = path_source_denguetracker_infodengue,
    root_dir_GT         = path_source_denguetracker_GT,
    state               = state_i,
    full_start          = ew_start,
    full_end            = ew_end,
    step_size           = step_size,
    D                   = delay_D,
    hypers              = hypers_q(phi_ref = 0.2, D_ref = 15, type = "exponential",
                                   alpha_phi = 1.4, sd_log_b = 1, delay_seq = 0:15),
    E                   = E_val,
    S                   = H_used * 2,
    bases_s             = fourier_basis[, 1:(H_used * 2)],
    model_file          = file.path(path_source, "models", "nowcast_hist_triangle.stan"),
    model_name          = "ARABMwGT",
    mode                = "expanding",
    iter_sampling       = iter_sampling,
    iter_warmup         = iter_warmup,
    chains              = chain_num,
    posterior_draws_path = file.path(path_proj, "draws", "length8")
  )
  saveRDS(results_wGT_i, file.path(fit_dir, sprintf("%s_ARABMwGT_fit.rds", state_i)))
}

## =========================
## (D) Loop 2: 读取结果并计算指标（RMSE/RMSPE/MAE/MAPE/CRPS/CR95/CR50）
## =========================
all_metrics <- list()
include_GT <- T

for (i in seq_along(all_states)) {
  state_i <- all_states[i]
  cat(sprintf(">>> [EVAL %d/%d] State = %s ...\n", i, length(all_states), state_i))
  
  baseline_i   <- readRDS(file.path(fit_dir, sprintf("%s_baseline.rds", state_i)))
  results_i    <- readRDS(file.path(fit_dir, sprintf("%s_ARABM_fit.rds", state_i)))
  results_wGT_i<- readRDS(file.path(fit_dir, sprintf("%s_ARABMwGT_fit.rds", state_i)))
  
  # 提取摘要（与你当前代码保持一致）
  arabm_i     <- extract_N_summary(results_i,    num_N = 6) # 默认名：ARABM, ARABM_lwr, ARABM_upr
  arabm_wGT_i <- extract_N_summary(results_wGT_i,num_N = 6, model_name = "ARABMwGT")
  
  result_out_i <- baseline_i %>%
    left_join(arabm_i %>% select(ew_now, ew, ARABM, ARABM_lwr, ARABM_upr),
              by = c("ew_now", "ew")) %>%
    left_join(arabm_wGT_i %>% select(ew_now, ew, ARABMwGT, ARABMwGT_lwr, ARABMwGT_upr),
              by = c("ew_now", "ew"))
  
  # 每个 ew_now 取“当前可见的最新一周”（与你原逻辑一致）
  n_last <- 1
  result_out_cut <- result_out_i %>%
    filter(!is.na(ARABM)) %>%
    group_by(ew_now) %>%
    slice_max(order_by = ew, n = n_last) %>%
    ungroup()
  
  truth <- result_out_cut$sum_of_cases_final
  
  #— 覆盖率（95% 用已有 lwr/upr；50% 用正态近似）
  cr95_ARABM   <- coverage_rate(truth, result_out_cut$ARABM_lwr,   result_out_cut$ARABM_upr)
  cr95_wGT     <- coverage_rate(truth, result_out_cut$ARABMwGT_lwr,result_out_cut$ARABMwGT_upr)
  # 50% 区间近似
  ci50_A <- make_50_ci_from_95(result_out_cut$ARABM,
                               result_out_cut$ARABM_lwr, result_out_cut$ARABM_upr)
  ci50_W <- make_50_ci_from_95(result_out_cut$ARABMwGT,
                               result_out_cut$ARABMwGT_lwr, result_out_cut$ARABMwGT_upr)
  cr50_ARABM <- coverage_rate(truth, ci50_A$lwr50, ci50_A$upr50)
  cr50_wGT   <- coverage_rate(truth, ci50_W$lwr50, ci50_W$upr50)
  
  #— CRPS（用正态近似）
  crps_ARABM <- crps_normal_vec(truth,
                                mu    = result_out_cut$ARABM,
                                lwr95 = result_out_cut$ARABM_lwr,
                                upr95 = result_out_cut$ARABM_upr)
  crps_wGT   <- crps_normal_vec(truth,
                                mu    = result_out_cut$ARABMwGT,
                                lwr95 = result_out_cut$ARABMwGT_lwr,
                                upr95 = result_out_cut$ARABMwGT_upr)
  
  #— 误差指标
  metric_tbl <- tibble(
    state  = state_i,
    model  = c("Bastos", "ARABM", "ARABMwGT"),
    RMSE   = c(
      rmse_fun(result_out_cut$cases_est_id, truth),
      rmse_fun(result_out_cut$ARABM,        truth),
      rmse_fun(result_out_cut$ARABMwGT,     truth)
    ),
    RMSPE  = c(
      safe_rmspe(result_out_cut$cases_est_id, truth),
      safe_rmspe(result_out_cut$ARABM,        truth),
      safe_rmspe(result_out_cut$ARABMwGT,     truth)
    ),
    MAE    = c(
      mae_fun(result_out_cut$cases_est_id, truth),
      mae_fun(result_out_cut$ARABM,        truth),
      mae_fun(result_out_cut$ARABMwGT,     truth)
    ),
    MAPE   = c(
      safe_mape(result_out_cut$cases_est_id, truth),
      safe_mape(result_out_cut$ARABM,        truth),
      safe_mape(result_out_cut$ARABMwGT,     truth)
    ),
    CRPS   = c(NA_real_, crps_ARABM, crps_wGT),  # Bastos 如无区间/σ，这里设 NA（或你也可加近似）
    CR95   = c(NA_real_, cr95_ARABM, cr95_wGT),
    CR50   = c(NA_real_, cr50_ARABM, cr50_wGT)
  )
  
  # 可选：加入 GT（若你要）
  if (include_GT && "GT" %in% names(result_out_cut)) {
    metric_tbl <- bind_rows(
      tibble(
        state = state_i,
        model = "GT",
        RMSE  = rmse_fun(result_out_cut$GT, truth),
        RMSPE = safe_rmspe(result_out_cut$GT, truth),
        MAE   = mae_fun(result_out_cut$GT, truth),
        MAPE  = safe_mape(result_out_cut$GT, truth),
        CRPS  = NA_real_,  # 若你能给出 GT 的区间，同样可用正态近似
        CR95  = NA_real_,
        CR50  = NA_real_
      ),
      metric_tbl
    )
  }
  
  all_metrics[[state_i]] <- metric_tbl
}

# 合并所有州指标
final_metrics <- bind_rows(all_metrics)

# 保存一份总表（可选）
out_dir <- file.path(path_proj, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(final_metrics, file.path(out_dir, "multi_state_metrics_with_crps_coverage.csv"))

cat(">>> Done. Combined metrics saved to: ",
    file.path(out_dir, "multi_state_metrics_with_crps_coverage.csv"), "\n")
