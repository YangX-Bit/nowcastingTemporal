test_states <- c("RJ","AC")

list_test <- list()
result_list_test <- list()
for (i in 1:length(test_states)) {
  list_test[[i]] <- run_moving_window(
    root_dir_infodengue = path_source_denguetracker_infodengue,
    root_dir_GT = path_source_denguetracker_GT,
    state = test_states[i],
    full_start = ew_start,
    full_end = "202420",
    window_length = 5,
    step_size = step_size,
    D = delay_D,
    hypers = hypers_q(phi_ref = 0.2, D_ref = 15, type = "exponential",
                      alpha_phi = 1.4, sd_log_b = 1, delay_seq = 0:15),
    model_file = file.path(path_source_nowcasting, "models", "b-ou.stan"),
    model_name = "ARABM",
    mode = "fixed",
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    chains = 3,
    posterior_draws_path = file.path(path_proj, "draws", "length5")
  )
  result_list_test[[i]] <- extract_N_summary(list_test[[i]], num_N = 1)
}

rbind(result_list_test)

baseline_list <- list()
result_out_list <- list()
for (i in 1:length(test_states)) {
  baseline_list[[i]] <- get_baseline_data(
    root_dir = path_source_denguetracker_data, # root file for data,
    states = test_states[i],
    start = ew_start,
    end = ew_end,
    last_n = 8
  )
  result_out_list[[i]] <- baseline_list[[i]] %>%
    left_join(
      result_list_test[[i]] %>% select(ew_now, ew, ARABM, ARABM_lwr, ARABM_upr),
      by = c("ew_now", "ew")
    ) 
    # left_join(
    #   arabm_wGT_length8 %>% select(ew_now, ew, ARABMwGT, ARABMwGT_lwr, ARABMwGT_upr),
    #   by = c("ew_now", "ew")
    # )
  
  result_out_list[[i]] <- add_actual_cases(result_out_list[[i]], 
                                         root_dir_infodengue = path_source_denguetracker_infodengue,
                                         test_states[i],
                                         202511)
  result_out_list[[i]] <- result_out_list[[i]] %>%
    filter(!is.na(ARABM))
}

result <- do.call(rbind, result_out_list)

