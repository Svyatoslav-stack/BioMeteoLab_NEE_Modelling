# =============================================================================
# NEE Random Forest — Separate Model per Collar
# =============================================================================
# For each collar (1–16), fits a separate Random Forest model using that
# collar's NEE, local predictors (PAR, temp, moisture), and shared logger
# predictors. No CSV/RDS written; summaries to console, plots only to plots/.
# If merged_data is not in the environment, it is loaded from data/Merged_data.csv
# (create it first by running Data_Merge.R).
# =============================================================================

# --- 1. Load packages ---
library(tidyverse)
library(lubridate)
library(ranger)
library(ggplot2)

# Load data if not already in environment
if (!exists("merged_data")) {
  merged_data <- read.csv("data/Merged_data.csv", stringsAsFactors = FALSE)
}

# --- 2. Ensure output directory for plots ---
if (!dir.exists("plots")) dir.create("plots", showWarnings = FALSE)

# --- 3. Constants and predictor definitions ---
n_collars <- 16L
min_obs_per_collar <- 30L
set.seed(42)

# Shared (logger) predictors — only these global columns
shared_predictors <- c(
  "Ts_air_Avg",
  "Ts_1_3_Avg",
  "Ts_1_5_Avg",
  "Ts_1_8_Avg",
  "VWC_1_3_Avg",
  "VWC_1_5_Avg",
  "VWC_1_8_Avg",
  "Ts_2_3_Avg",
  "Ts_2_5_Avg",
  "Ts_2_8_Avg",
  "VWC_2_3_Avg",
  "VWC_2_5_Avg",
  "VWC_2_8_Avg",
  "Rain_mm_Tot",
  "hour",
  "doy"
)

# --- 4. Prepare base data: TIMESTAMP as POSIXct, add hour and doy if missing ---
if (!inherits(merged_data$TIMESTAMP, "POSIXct")) {
  merged_data$TIMESTAMP <- as.POSIXct(merged_data$TIMESTAMP, tz = "UTC")
}
if (!"hour" %in% names(merged_data)) {
  merged_data$hour <- as.integer(lubridate::hour(merged_data$TIMESTAMP))
}
if (!"doy" %in% names(merged_data)) {
  merged_data$doy <- as.integer(lubridate::yday(merged_data$TIMESTAMP))
}

# --- 5. Helper: build modeling dataframe for one collar (safe column names) ---
# Returns a tibble with columns: NEE, PAR_collar, temp_collar, moisture_collar,
# plus all shared_predictors. Returns NULL if required columns are missing or
# too few complete rows.
build_collar_df <- function(data, i) {
  nee_col    <- paste0("NEE.CO2_collars_", i)
  par_col    <- paste0("PAR_collars_", i)
  temp_col   <- paste0("temp_collars_", i)
  moist_col  <- paste0("moisture_collars_", i)

  needed <- c(nee_col, par_col, temp_col, moist_col, shared_predictors)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0L) {
    warning("Collar ", i, ": missing columns: ", paste(missing, collapse = ", "))
    return(NULL)
  }

  # Select by name; then rename to syntactic names (no NSE with dynamic names)
  cols <- c(nee_col, par_col, temp_col, moist_col, shared_predictors)
  new_names <- c("NEE", "PAR_collar", "temp_collar", "moisture_collar", shared_predictors)

  df <- data %>% select(all_of(cols))
  names(df) <- new_names

  df <- df %>% filter(!is.na(.data$NEE)) %>% filter(complete.cases(.))
  if (nrow(df) < min_obs_per_collar) {
    message("Collar ", i, ": too few complete rows (n = ", nrow(df), "), skipping.")
    return(NULL)
  }

  df
}

# --- 6. Helper: compute test-set metrics ---
compute_metrics <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2     <- if (ss_tot > 0) 1 - (ss_res / ss_tot) else NA_real_
  rmse   <- sqrt(mean((obs - pred)^2))
  mae    <- mean(abs(obs - pred))
  list(r2 = r2, rmse = rmse, mae = mae)
}

# --- 7. Loop over collars: fit RF, evaluate, store results ---
model_summary_list <- list()
importance_list   <- list()

for (i in seq_len(n_collars)) {

  df <- build_collar_df(merged_data, i)
  if (is.null(df)) next

  n_obs <- nrow(df)
  train_idx <- sample.int(n_obs, size = floor(0.5 * n_obs))
  train_data <- df[train_idx, ]
  test_data  <- df[-train_idx, ]

  # Formula: local predictors + all shared_predictors (built from vector)
  rf_formula <- as.formula(paste(
    "NEE ~ PAR_collar + temp_collar + moisture_collar +",
    paste(shared_predictors, collapse = " + ")
  ))
  rf_fit <- ranger(
    formula   = rf_formula,
    data     = train_data,
    num.trees = 500L,
    importance = "permutation",
    replace  = TRUE,
    seed     = 42
  )

  pred_test <- predict(rf_fit, data = test_data)$predictions
  met <- compute_metrics(test_data$NEE, pred_test)

  model_summary_list[[i]] <- tibble(
    collar = i,
    n_obs  = n_obs,
    r2     = met$r2,
    rmse   = met$rmse,
    mae    = met$mae
  )

  imp <- importance(rf_fit)
  importance_list[[i]] <- tibble(
    collar     = i,
    variable   = names(imp),
    importance = as.numeric(imp)
  )
}

# --- 8. Combine results (drop NULL from skipped collars) ---
model_summary <- bind_rows(model_summary_list)
var_importance <- bind_rows(importance_list)

# --- 9. Print summary tables to console ---
cat("\n========== NEE Random Forest by Collar — Metrics ==========\n\n")
print(model_summary)
cat("\n========== Variable importance (long format) ==========\n\n")
print(var_importance)

# --- 10. Best / worst collar and average R2 ---
if (nrow(model_summary) > 0L) {
  best  <- model_summary %>% slice_max(r2, n = 1)
  worst <- model_summary %>% slice_min(r2, n = 1)
  avg_r2 <- mean(model_summary$r2, na.rm = TRUE)
  cat("\n---------- Summary stats ----------\n")
  cat("Best collar by R2:  ", best$collar, " (R2 = ", round(best$r2, 4), ")\n", sep = "")
  cat("Worst collar by R2: ", worst$collar, " (R2 = ", round(worst$r2, 4), ")\n", sep = "")
  cat("Average R2 (all modeled collars): ", round(avg_r2, 4), "\n", sep = "")

  top5 <- var_importance %>%
    group_by(variable) %>%
    summarise(mean_importance = mean(importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_importance)) %>%
    slice_head(n = 5)
  cat("\nTop 5 most important variables (mean across collars):\n")
  print(top5)
}
cat("\n========== End of summary ==========\n")

# --- 11. Plots (only into plots/) ---
if (nrow(model_summary) > 0L) {

  # A. R2 by collar
  p_r2 <- ggplot(model_summary, aes(x = factor(collar), y = r2, fill = r2)) +
    geom_col() +
    scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 0.9) +
    labs(title = "R² by collar (Random Forest, test set)", x = "Collar", y = "R²") +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.major.x = element_blank())
  ggsave("plots/rf_r2_by_collar.png", p_r2, width = 8, height = 5, dpi = 150)

  # B. Heatmap of variable importance by collar (красный → белый (0) → синий)
  # Топ-3 фактора по важности для каждого коллара помечаем рангом 1, 2, 3
  top3_per_collar <- var_importance %>%
    group_by(collar) %>%
    slice_max(importance, n = 3, with_ties = FALSE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  p_heat <- ggplot(var_importance, aes(x = factor(collar), y = variable, fill = importance)) +
    geom_tile(colour = "grey85", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      limits = c(-500, 7500),
      name = "Importance"
    ) +
    geom_text(
      data = top3_per_collar,
      aes(x = factor(collar), y = variable, label = rank),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold",
      colour = "black"
    ) +
    labs(title = "Variable importance by collar (permutation). 1–3 = top 3 factors per collar", x = "Collar", y = "Variable") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9), panel.grid = element_blank())
  ggsave("plots/rf_importance_heatmap_by_collar.png", p_heat, width = 10, height = 6, dpi = 150)

  # C. RMSE by collar
  p_rmse <- ggplot(model_summary, aes(x = factor(collar), y = rmse, fill = rmse)) +
    geom_col() +
    scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 0.9) +
    labs(title = "RMSE by collar (Random Forest, test set)", x = "Collar", y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.major.x = element_blank())
  ggsave("plots/rf_rmse_by_collar.png", p_rmse, width = 8, height = 5, dpi = 150)
}

# =============================================================================
# End of script
# =============================================================================
