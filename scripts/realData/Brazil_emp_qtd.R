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
    uf, D = 51, if_last_D_cols_NA = FALSE
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


# 4a. Plot SP as a standalone panel
sp_df <- filter(all_dfs, UF == "SP")
p_sp <- ggplot(sp_df, aes(x = Delay, y = Value, group = Time)) +
  geom_line(color = "grey40", linewidth = 0.4, alpha = 0.5) +
  geom_vline(xintercept = 15, color = "red", linewidth = 0.8) +
  geom_vline(xintercept = 20, color = "red", linewidth = 0.6) +
  geom_vline(xintercept = 25, color = "red", linewidth = 0.6) +
  geom_hline(yintercept = 0.95, color = "red", linewidth = 0.6) +
  geom_hline(yintercept = 0.90, color = "red", linewidth = 0.6) +
  labs(
    title    = "São Paulo (SP): Normalized Delay Distributions",
    subtitle = "Red line at 15,20,25-week delay (D = 15,20,25)",
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

# 4b. Plot the remaining 25 states in a 5x5 grid
others_df <- filter(all_dfs, UF != "SP")
p_others <- ggplot(others_df, aes(x = Delay, y = Value, group = Time)) +
  geom_line(color = "grey40", linewidth = 0.3, alpha = 0.4) +
  geom_vline(xintercept = 15, color = "red", linewidth = 0.6) +
  geom_vline(xintercept = 20, color = "red", linewidth = 0.6) +
  geom_vline(xintercept = 25, color = "red", linewidth = 0.6) +
  geom_hline(yintercept = 0.95, color = "red", linewidth = 0.6) +
  geom_hline(yintercept = 0.90, color = "red", linewidth = 0.6) +
  facet_wrap(~ UF, ncol = 5) +
  labs(
    title    = "Brazilian States (excluding SP): Normalized Delay Distributions",
    subtitle = "Each panel shows one state; red line at 15,20,25-week delay",
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

