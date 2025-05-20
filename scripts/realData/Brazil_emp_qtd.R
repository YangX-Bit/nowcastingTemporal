# 1. Load required packages
library(tidyverse)

# 2. Define the list of federative units
brazil_ufs <- c(
  "AC", "AL", "AP", "AM", "BA", "CE", "DF", "GO",
  "MA", "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI",
  "RJ", "RN", "RS", "RO", "RR", "SC", "SP", "SE", "TO"
)

D <- 30
# 3. Loop over each UF, retrieve & standardize the matrix, pivot to long, and bind
all_dfs <- map_dfr(brazil_ufs, function(uf) {
  print(uf)
  temp_mat <- get_infodengue_data(
    path_source_denguetracker_infodengue,
    "2024-03-03", "2024-12-29",
    uf, D = 65, if_last_D_cols_NA = FALSE
  )[[1]]
  used_mat <- temp_mat[ , c(0:D+1)]
  used_mat_std <- temp_mat[,ncol(temp_mat)]
  used_mat_std <- round(used_mat / used_mat_std, 2)
  df <- as.data.frame(used_mat_std)
  colnames(df) <- paste0("Delay_", 0:D)
  df$Time <- seq_len(nrow(df))
  df_long <- df %>%
    pivot_longer(
      cols      = starts_with("Delay_"),
      names_to  = "Delay",
      values_to = "Value"
    ) %>%
    mutate(
      Delay = as.integer(str_remove(Delay, "Delay_")),
      UF    = uf
    )
  return(df_long)
})


# 1) For each week, each time, the week reach 95%
daily_q95 <- all_dfs %>%
  group_by(UF, Time) %>%
  summarise(
    week95 = {
      valid_idx <- !is.na(Value) & Value >= 0.95
      if (any(valid_idx)) {
        min(Delay[valid_idx], na.rm = TRUE)
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# 2) Average
state_summary <- daily_q95 %>%
  group_by(UF) %>%
  summarise(
    mean_week95   = mean(week95, na.rm = TRUE),
    median_week95 = median(week95, na.rm = TRUE),
    n_obs         = sum(!is.na(week95)),     
    .groups       = "drop"
  )
  


# 4a. Plot SP as a standalone panel using its average week95
sp_avg <- state_summary %>% filter(UF == "SP") %>% pull(mean_week95)

sp_df <- filter(all_dfs, UF == "SP")
p_sp <- ggplot(sp_df, aes(x = Delay, y = Value, group = Time)) +
  geom_line(color = "grey40", linewidth = 0.4, alpha = 0.5) +
  geom_vline(xintercept = sp_avg, color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0.95, color = "red", linewidth = 0.6) +
  labs(
    title    = "São Paulo (SP): Normalized Delay Distributions",
    subtitle = paste0("Average 95% quantile delay: ", round(sp_avg, 1), " weeks"),
    x        = "Delay (weeks)",
    y        = "Standardized Reporting Proportion"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    panel.grid    = element_line(color = "grey90", linewidth = 0.3)
  )

print(p_sp)

# 4b. Plot the remaining 25 states in a 5x5 grid, each with its own average week95 line
others_df <- filter(all_dfs, !UF %in% c("SP", "ES"))
# Plot with annotated red line labels
p_others <- ggplot(others_df, aes(x = Delay, y = Value, group = Time)) +
  geom_line(color = "grey40", linewidth = 0.3, alpha = 0.4) +
  # Vertical lines for each state's mean_week95
  geom_vline(
    data = state_summary %>% filter(!UF %in% c("SP", "ES")),
    aes(xintercept = mean_week95),
    color = "red", linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  # Annotate each red line with its numeric value
  geom_text(
    data = state_summary %>% filter(!UF %in% c("SP", "ES")),
    aes(
      x     = mean_week95+2,
      y     = 0.02,                   # slightly above the maximum of Value (assuming max ≤ 1)
      label = round(mean_week95, 1)
    ),
    color    = "red",
    size     = 3,
    fontface = "bold",
    inherit.aes = FALSE,
    vjust    = 0
  ) +
  geom_hline(yintercept = 0.95, color = "red", linewidth = 0.6) +
  facet_wrap(~ UF, ncol = 5) +
  labs(
    title    = "Brazilian States (excluding SP, ES): Normalized Delay Distributions",
    subtitle = "Red vertical line: state-specific average 95% quantile delay",
    x        = "Delay (weeks)",
    y        = "Standardized Reporting Proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    strip.text    = element_text(face = "bold"),
    panel.grid    = element_line(color = "grey90", linewidth = 0.2)
  )


print(p_others)


# After creating your ggplot (e.g. p_sp or p_others), save it like this:

# Save São Paulo plot
ggsave(
  filename = "SP_delay_plot_D30.png",  # output file
  plot     = p_sp,                 # which plot object
  device   = "png",                # file format
  width    = 10,                    # width in inches
  height   = 8,                    # height in inches
  units    = "in",                 # units: "in", "cm", or "mm"
  dpi      = 300                   # resolution
)

# Save the multi‐panel figure
ggsave(
  filename = "Brazil_states_delay_D30.png",
  plot     = p_others,
  device   = "png",
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)

