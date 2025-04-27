path_source_denguetracker_GT <- file.path(path_source_denguetracker_data, "gtrends")

win_end_ew <- "202452"
brazil_ufs <- c(
  "AC", "AL", "AP", "AM", "BA", "CE", "DF", "GO",
  "MA", "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI",
  "RJ", "RN", "RS", "RO", "RR", "SC", "SP", "SE", "TO" # "BR"
)
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

latest_week <- "202512"
latest_date <- "2025-03-23"

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
cor.index.log <- numeric(length(InfoDengue))

for (i in seq_along(InfoDengue)) {
  plot_data <- data.frame(
    date = as.Date(InfoDengue[[i]]$ew_start),
    scaled_cases = InfoDengue[[i]]$scaled_cases,
    log_cases = InfoDengue[[i]]$log_cases,
    gt_index = GT_index[[i]][, 1]
  ) 
  cor.index[i] <- cor(plot_data$scaled_cases, plot_data$gt_index)
  cor.index.log[i] <- cor(plot_data$log_cases, plot_data$gt_index)

  p[[i]] <- ggplot(plot_data, aes(x = date)) +
    geom_line(aes(y = log_cases), color = "blue", size = 1) +
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




