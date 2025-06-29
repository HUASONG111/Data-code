# install.packages(c("dplyr", "ggplot2", "rstatix"))
library(dplyr)
library(ggplot2)
library(rstatix)

read_and_clean <- function(values, group_name) {
  df <- data.frame(value = values)
  df <- df %>%
    filter(!is.na(value), value != -999) %>%
    mutate(group = group_name)
  return(df)
}

## Since some groups in the small randomly selected subset of the real data have zero variance,
## which prevents the Games-Howell test from running properly, simulated data is used here 
## to ensure the code executes correctly.

##############Lag Onset time##########
set.seed(3)  # 
df_po4 <- read_and_clean(rnorm(22, mean = 2, sd = 1), 'PO4')
df_no3 <- read_and_clean(rnorm(25, mean = 1, sd = 1), 'NO3')
df_sst <- read_and_clean(rnorm(25, mean = 2.5, sd = 1.2), 'SST')
df_ssw <- read_and_clean(rnorm(25, mean = 3.5, sd = 1), 'SSW')

data <- bind_rows(df_po4, df_no3, df_sst, df_ssw)

group_means <- data %>%
  group_by(group) %>%
  summarise(mean_val = mean(value), .groups = "drop")

gh_res <- data %>%
  games_howell_test(value ~ group) %>%
  mutate(p.adj.signif = case_when(
    p.adj < 0.001  ~ "***",
    p.adj < 0.01   ~ "**",
    p.adj < 0.05   ~ "*",
    TRUE           ~ ""
  )) %>%

  left_join(group_means, by = c("group1" = "group")) %>%
  rename(mean1 = mean_val) %>%
  left_join(group_means, by = c("group2" = "group")) %>%
  rename(mean2 = mean_val) %>%
  rowwise() %>%
  mutate(
    large_group = if_else(mean1 >= mean2, group1, group2),
    small_group = if_else(mean1 < mean2, group1, group2),
    large_mean = max(mean1, mean2),
    small_mean = min(mean1, mean2),
    corrected_estimate = large_mean - small_mean,
    corrected_conf_low = if_else(mean1 >= mean2, conf.low, -conf.high),
    corrected_conf_high = if_else(mean1 >= mean2, conf.high, -conf.low),
    corrected_comparison = paste(large_group, "vs", small_group)
  ) %>%
  ungroup()


ggplot(gh_res, aes(x = reorder(corrected_comparison, corrected_estimate),
                    y = corrected_estimate,
                    label = p.adj.signif)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymax = -corrected_conf_low, ymin = -corrected_conf_high), width = 0.2) +
  geom_text(nudge_y = 0.05, size = 5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(x = "Pairwise Comparison",
       y = "Mean Difference") +
  theme_minimal()


#################Lag Duration###################

read_and_clean <- function(values, group_name) {
  df <- data.frame(value = values)
  df <- df %>%
    filter(!is.na(value), value != 0) %>%
    mutate(group = group_name)
  return(df)
}

set.seed(3)  # 固定随机种子，确保结果可复现

df_po4 <- read_and_clean(rnorm(22, mean = 6, sd = 1), 'PO4')
df_no3 <- read_and_clean(rnorm(25, mean = 7, sd = 1), 'NO3')
df_sst <- read_and_clean(rnorm(25, mean = 6.5, sd = 1.2), 'SST')
df_ssw <- read_and_clean(rnorm(25, mean = 8, sd = 1), 'SSW')

data <- bind_rows(df_po4, df_no3, df_sst, df_ssw)

group_means <- data %>%
  group_by(group) %>%
  summarise(mean_val = mean(value), .groups = "drop")

gh_res <- data %>%
  games_howell_test(value ~ group) %>%
  mutate(p.adj.signif = case_when(
    p.adj < 0.001  ~ "***",
    p.adj < 0.01   ~ "**",
    p.adj < 0.05   ~ "*",
    TRUE           ~ ""
  )) %>%
  left_join(group_means, by = c("group1" = "group")) %>%
  rename(mean1 = mean_val) %>%
  left_join(group_means, by = c("group2" = "group")) %>%
  rename(mean2 = mean_val) %>%
  rowwise() %>%
  mutate(
    large_group = if_else(mean1 >= mean2, group1, group2),
    small_group = if_else(mean1 < mean2, group1, group2),
    large_mean = max(mean1, mean2),
    small_mean = min(mean1, mean2),
    corrected_estimate = large_mean - small_mean,
    corrected_conf_low = if_else(mean1 >= mean2, conf.low, -conf.high),
    corrected_conf_high = if_else(mean1 >= mean2, conf.high, -conf.low),
    corrected_comparison = paste(large_group, "vs", small_group)
  ) %>%
  ungroup()

ggplot(gh_res, aes(x = reorder(corrected_comparison, corrected_estimate),
                   y = corrected_estimate,
                   label = p.adj.signif)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymax = -corrected_conf_low, ymin = -corrected_conf_high), width = 0.2) +
  geom_text(nudge_y = 0.05, size = 5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(x = "Pairwise Comparison",
       y = "Mean Difference") +
  theme_minimal()


