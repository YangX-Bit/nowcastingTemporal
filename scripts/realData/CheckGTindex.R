library(tidyverse)
path_proj = here::here()
path_source_denguetracker_data <- file.path(dirname(path_proj), "dengue-tracker/data/weekly_data/")
path_source_denguetracker_GT <- file.path(path_source_denguetracker_data, "gtrends")

win_end_ew <- "202452"
brazil_ufs <- c(
  "AC", "AL", "AP", "AM", "BA", "CE", "DF", "GO",
  "MA", "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI",
  "RJ", "RN", "RS", "RO", "RR", "SC", "SP", "SE", "TO" # "BR"
)

# to check the full year 2024
start_date <- "2024-01-01"
end_date <- "2024-12-31"

#google index
GT_index <- list() 
for (i in c(1:length(brazil_ufs))) {
  GT_index[[i]] <- read.csv(file.path(path_source_denguetracker_GT,win_end_ew, paste0(brazil_ufs[i],"_trends.csv")), 
           stringsAsFactors = FALSE, skip = 2)[,c(1:3)] %>%
    setNames(c("date", "dengue", "sintomas_dengue")) %>%
    filter(date <=  end_date & date >= start_date) %>%
    mutate(
      across(
        .cols = -date,
        .fns = ~ as.numeric(ifelse(. == "<1", "0.5", .))
      )
    ) %>% 
    select("dengue", "sintomas_dengue") %>%
    as.matrix()
}

# latest number of cases to compare
latest_week <- "202518"
latest_date <- "2025-05-04"

InfoDengue_root_dir <- file.path(path_source_denguetracker_data, "infodengue")

InfoDengue <- list()
for (i in c(1:length(brazil_ufs))) {
  InfoDengue[[i]] <- read.csv(file.path(InfoDengue_root_dir, latest_week, paste0(brazil_ufs[i], "_",latest_date,"_infodengue.csv"))) %>%
    filter(ew_start <=  end_date & ew_start >= start_date) %>%
    select(ew_start, sum_of_cases) %>%
    mutate(
      scaled_cases = round(sum_of_cases / max(sum_of_cases, na.rm = TRUE) * 100, 2),
      log_cases = round(log(sum_of_cases),2)
    )
}

p <- list()
cor.index <- numeric(length(InfoDengue))
cor.index.scale <- numeric(length(InfoDengue))
cor.index.log <- numeric(length(InfoDengue))

for (i in seq_along(InfoDengue)) {
  plot_data <- data.frame(
    date = as.Date(InfoDengue[[i]]$ew_start),
    sum_of_cases = InfoDengue[[i]]$sum_of_cases,
    scaled_cases = InfoDengue[[i]]$scaled_cases,
    log_cases = InfoDengue[[i]]$log_cases,
    gt_index = GT_index[[i]][, 1]
  ) 
  cor.index[i] <- cor(plot_data$sum_of_cases, plot_data$gt_index)
  cor.index.scale[i] <- cor(plot_data$scaled_cases, plot_data$gt_index)
  cor.index.log[i] <- cor(plot_data$log_cases, plot_data$gt_index)

  p[[i]] <- ggplot(plot_data, aes(x = date)) +
    geom_line(aes(y = scaled_cases), color = "blue", size = 1) +
    geom_line(aes(y = gt_index), color = "red", size = 1) +
    labs(
      title = paste0("Comparison Plot of ", brazil_ufs[i]),
      x = "Date",
      y = "Index Value",
      caption = "Blue: InfoDengue | Red: GT_index"
    ) +
    theme_minimal()
  print(p[[i]])
}

######################################
# Test the corr between GT and cases
cor.index_dengue <- numeric(length(InfoDengue))
cor.index_sintoma_dengue <- numeric(length(InfoDengue))
for (i in seq_along(InfoDengue)) {
  plot_data <- data.frame(
    date = as.Date(InfoDengue[[i]]$ew_start),
    sum_of_cases = InfoDengue[[i]]$sum_of_cases,
    gt_index_dengue = GT_index[[i]][, 1],
    gt_index_sintoma_dengue = GT_index[[i]][, 2]
  ) 
  cor.index_dengue[i] <- cor(plot_data$sum_of_cases, plot_data$gt_index_dengue)
  cor.index_sintoma_dengue[i] <- cor(plot_data$sum_of_cases, plot_data$gt_index_sintoma_dengue)
}
cor.index_dengue <- round(cor.index_dengue, 3)
cor.index_sintoma_dengue <- round(cor.index_sintoma_dengue, 3)

cor.index_dengue
cor.index_sintoma_dengue

t.test(cor.index_dengue, cor.index_sintoma_dengue, paired = T)

# Assume cor.index_dengue and cor.index_sintoma_dengue are two numeric vectors of equal length

#  Fisher z–transform and paired t‐test on z
z1 <- atanh(cor.index_dengue)
z2 <- atanh(cor.index_sintoma_dengue)
t_test_z <- t.test(z1, z2, paired = TRUE)
print(t_test_z)

#  Wilcoxon signed‐rank test on z (nonparametric)
wilcox_z <- wilcox.test(z1, z2, paired = TRUE)
print(wilcox_z)


#  Scatter plot and identity line
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
df <- data.frame(
  dengue = cor.index_dengue,
  sintoma = cor.index_sintoma_dengue
)

#  Bland–Altman plot (mean vs. difference)
df$mean_val <- rowMeans(df)
df$diff    <- df$dengue - df$sintoma
ba_plot <- ggplot(df, aes(x = mean_val, y = diff)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  geom_hline(yintercept = mean(df$diff), color = "blue") +
  geom_hline(yintercept = mean(df$diff) + 1.96*sd(df$diff),
             linetype = "dashed") +
  geom_hline(yintercept = mean(df$diff) - 1.96*sd(df$diff),
             linetype = "dashed") +
  labs(
    x = "Mean of Paired Correlations",
    y = "Difference (dengue – sintoma.dengue)",
    title = NULL
  ) +
  theme_classic()
print(ba_plot)


