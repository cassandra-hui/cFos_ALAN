# 1/9/2026 Getting data for Jenny

library(tidyverse)
library(dplyr)

cfos_df <- read.csv('cfos_df.csv')

head(cfos_df)

table(cfos_df$area)

# NCL = PLN

# Filter for NCL area and calculate mean and SE for ZENK expression
ncl_zenk <- cfos_df %>%
  filter(area == "PLN") %>%
  summarise(
    n = length(ZENK_percent),
    mean_ZENK = mean(ZENK_percent, na.rm = TRUE),
    se_ZENK = sd(ZENK_percent, na.rm = TRUE) / sqrt(n)
  )

print(ncl_zenk)

# Mean and SE by group
ncl_zenk_by_group <- cfos_df %>%
  filter(area == "PLN") %>%
  group_by(group) %>%
  summarise(
    mean_ZENK = mean(ZENK_percent, na.rm = TRUE),
    se_ZENK = sd(ZENK_percent, na.rm = TRUE) / sqrt(sum(!is.na(ZENK_percent)))
  )

ncl_zenk_by_group


 