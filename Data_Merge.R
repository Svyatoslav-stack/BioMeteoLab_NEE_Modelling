# Iisaku Collars + Logger Merge (see README.md for data and usage)

library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)


# ---- 1. Paths ----

path_collars <- "data/Iisaku collars_2018-19.xlsx"
path_logger  <- "data/Iisaku loger 2018-19.xlsx"


# ---- 2. Load logger ----

logger <- read_excel(path_logger, sheet = 1)

logger <- logger %>%
  mutate(TIMESTAMP = round_date(TIMESTAMP, "10 minutes"))


# ---- 3. Load all collar sheets ----

sheet_names <- excel_sheets(path_collars)

collars_all <- map_dfr(seq_along(sheet_names), function(i) {
  read_excel(path_collars, sheet = sheet_names[i]) %>%
    rename(TIMESTAMP = date) %>%
    mutate(collar_id = paste0("collars_", i))
})


# ---- 4. Align timestamps and reshape to wide ----

# Collar times are irregular; round to 10-min grid to match logger
collars_long <- collars_all %>%
  mutate(TIMESTAMP = round_date(as.POSIXct(TIMESTAMP), "10 minutes"))

value_columns <- setdiff(names(collars_long), c("collar_id", "TIMESTAMP"))

collars_wide <- collars_long %>%
  pivot_wider(
    id_cols     = TIMESTAMP,
    names_from  = collar_id,
    values_from = all_of(value_columns)
  ) %>%
  arrange(TIMESTAMP)


# ---- 5. Merge with logger and export ----

merged_data <- collars_wide %>%
  left_join(logger, by = "TIMESTAMP")

write.csv(merged_data, "data/Merged_data.csv", row.names = FALSE)
