# ============================================================
# Quantile Regression Forest simulation analysis
# Standard R / RStudio script 
# ============================================================

# Block 1 — Set up the R environment and load the simulation data

# Install missing packages from CRAN, then load them.
required_packages <- c(
  "quantregForest", "ggplot2", "tidyr", "dplyr",
  "ggridges", "patchwork", "ggnewscale", "viridis"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

cat("R packages loaded successfully.\n")
cat("R version:", R.version.string, "\n")

# Place the CSV in the current working directory, or update the path below.
# In RStudio, open the project folder before running this script.
# To inspect the current working directory, run getwd().
file_name <- "grid_for_QRF_analysis_overlapping_regions.csv"

if (!file.exists(file_name)) {
  stop(
    "Input CSV not found: ", file_name,
    "\nCurrent working directory: ", getwd(),
    "\nPlace the CSV in this directory or update 'file_name'."
  )
}

# Read the simulation data.
df_qrf <- read.csv(file_name, stringsAsFactors = FALSE)

# Convert model variables to the required data types.
df_qrf$DOY <- as.numeric(df_qrf$DOY)
df_qrf$mean_temp <- as.numeric(df_qrf$mean_temp)
df_qrf$std_temp <- as.numeric(df_qrf$std_temp)
if ("region" %in% names(df_qrf)) {
  df_qrf$region <- as.factor(df_qrf$region)
}

cat("Simulation data loaded successfully.\n")
cat("Number of observations:", nrow(df_qrf), "\n")
cat("Columns:", paste(names(df_qrf), collapse = ", "), "\n\n")

# Preview the dataset and summarize the response variable.
print(head(df_qrf))
cat("\nDOY summary statistics:\n")
print(summary(df_qrf$DOY))

# Display the number of observations in each sampled region.
if ("region" %in% names(df_qrf)) {
  cat("\nObservations by region:\n")
  print(table(df_qrf$region))
}

# =========================================
# Block 2 — Train QRF models using five sampled regions
# =========================================

# ---- Features and response ----
features <- c("mean_temp", "std_temp")
response <- "DOY"

# ---- Quantiles ----
quantiles_to_predict <- c(0.25, 0.50, 0.75)

# ---- Regions used for training ----
training_regions <- c(
  "Lower_left",
  "Upper_left",
  "Lower_right",
  "Upper_right",
  "Center"
)

# ---- Containers ----
qrf_models <- list()
qrf_train_results <- list()

set.seed(42)

# =========================================
# Train one QRF model for each region
# =========================================

 for (reg in training_regions) {

  cat("\n=========================================\n")
  cat("Training QRF model using region:", reg, "\n")
  cat("=========================================\n")

  # ---- Training data ----
  df_train <- df_qrf %>%
    filter(region == reg) %>%
    drop_na(all_of(c(features, response)))

  cat("Training sample size:", nrow(df_train), "\n")

  # ---- Fit QRF ----
  qrf_model <- quantregForest(
    x = df_train[, features],
    y = df_train[[response]],
    ntree = 1000
  )

  # ---- In-sample prediction ----
  qrf_pred <- predict(
    qrf_model,
    newdata = df_train[, features],
    what = quantiles_to_predict
  )

  # ---- Save predictions ----
  df_result <- df_train %>%
    mutate(
      pred_25 = qrf_pred[, 1],
      pred_50 = qrf_pred[, 2],
      pred_75 = qrf_pred[, 3],
      train_region = reg
    )

  # ---- Store model and result ----
  qrf_models[[reg]] <- qrf_model
  qrf_train_results[[reg]] <- df_result

  cat("✅ Finished:", reg, "\n")
}

# ---- Combine all in-sample results ----
df_qrf_train_results <- bind_rows(qrf_train_results)

cat("\n✅ All five QRF models trained successfully!\n")
cat("Combined result rows:", nrow(df_qrf_train_results), "\n")

# ---- Preview ----
print(head(df_qrf_train_results))

# ---- Check model names ----
print(names(qrf_models))

# =========================================
# Block 3 — Center-trained QRF:
# In-sample metrics first, then cross-region validation
# =========================================

# ---- Metric functions ----
pinball_loss <- function(y, y_pred, tau) {
  err <- y - y_pred
  mean(ifelse(err >= 0, tau * err, (tau - 1) * err))
}

calculate_qrf_metrics <- function(y, p25, p50, p75) {
  data.frame(
    ICP = mean(y >= p25 & y <= p75),
    MPIW = mean(p75 - p25),
    MAE_25 = mean(abs(y - p25)),
    MAE_50 = mean(abs(y - p50)),
    MAE_75 = mean(abs(y - p75)),
    Pinball_25 = pinball_loss(y, p25, 0.25),
    Pinball_50 = pinball_loss(y, p50, 0.50),
    Pinball_75 = pinball_loss(y, p75, 0.75)
  )
}

# ---- Training model: Center ----
train_region <- "Center"
center_model <- qrf_models[[train_region]]

# =========================================
# Part 1 — In-sample metrics for Center
# =========================================

df_center_train <- df_qrf %>%
  filter(region == train_region) %>%
  drop_na(all_of(c(features, response)))

center_pred_in <- predict(
  center_model,
  newdata = df_center_train[, features],
  what = c(0.25, 0.50, 0.75)
)

df_center_in_sample <- df_center_train %>%
  mutate(
    train_region = train_region,
    test_region = "Center_in_sample",
    pred_25 = center_pred_in[, 1],
    pred_50 = center_pred_in[, 2],
    pred_75 = center_pred_in[, 3],
    interval_width = pred_75 - pred_25,
    in_interval = DOY >= pred_25 & DOY <= pred_75
  )

df_center_in_sample_metrics <- calculate_qrf_metrics(
  y = df_center_in_sample$DOY,
  p25 = df_center_in_sample$pred_25,
  p50 = df_center_in_sample$pred_50,
  p75 = df_center_in_sample$pred_75
) %>%
  mutate(
    train_region = train_region,
    test_region = "Center_in_sample",
    n_test = nrow(df_center_in_sample),
    .before = 1
  )

cat("\n=========================================\n")
cat("In-sample metrics: Center-trained model on Center data\n")
cat("=========================================\n")
print(df_center_in_sample_metrics)

# =========================================
# Part 2 — Cross-region validation
# =========================================

test_regions <- c(
  "Lower_left",
  "Upper_left",
  "Lower_right",
  "Upper_right"
)

center_validation_results <- list()
center_prediction_results <- list()

for (test_reg in test_regions) {

  cat("\n=========================================\n")
  cat("Training region:", train_region, "\n")
  cat("Testing region :", test_reg, "\n")
  cat("=========================================\n")

  df_test <- df_qrf %>%
    filter(region == test_reg) %>%
    drop_na(all_of(c(features, response)))

  pred <- predict(
    center_model,
    newdata = df_test[, features],
    what = c(0.25, 0.50, 0.75)
  )

  df_pred <- df_test %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      pred_25 = pred[, 1],
      pred_50 = pred[, 2],
      pred_75 = pred[, 3],
      interval_width = pred_75 - pred_25,
      in_interval = DOY >= pred_25 & DOY <= pred_75
    )

  metric_df <- calculate_qrf_metrics(
    y = df_pred$DOY,
    p25 = df_pred$pred_25,
    p50 = df_pred$pred_50,
    p75 = df_pred$pred_75
  ) %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      n_test = nrow(df_pred),
      .before = 1
    )

  center_validation_results[[test_reg]] <- metric_df
  center_prediction_results[[test_reg]] <- df_pred

  print(metric_df)
}

df_center_validation_metrics <- bind_rows(
  df_center_in_sample_metrics,
  bind_rows(center_validation_results)
)

df_center_validation_predictions <- bind_rows(
  df_center_in_sample,
  bind_rows(center_prediction_results)
)

cat("\n✅ Center-trained in-sample + cross-region validation completed!\n\n")
print(df_center_validation_metrics)

# =========================================
# Block 4 — Upper-right-trained QRF:
# In-sample metrics first, then cross-region validation
# =========================================

# ---- Training model: Upper_right ----
train_region <- "Upper_right"
upper_right_model <- qrf_models[[train_region]]

# =========================================
# Part 1 — In-sample metrics for Upper_right
# =========================================

df_upper_right_train <- df_qrf %>%
  filter(region == train_region) %>%
  drop_na(all_of(c(features, response)))

upper_right_pred_in <- predict(
  upper_right_model,
  newdata = df_upper_right_train[, features],
  what = c(0.25, 0.50, 0.75)
)

df_upper_right_in_sample <- df_upper_right_train %>%
  mutate(
    train_region = train_region,
    test_region = "Upper_right_in_sample",
    pred_25 = upper_right_pred_in[, 1],
    pred_50 = upper_right_pred_in[, 2],
    pred_75 = upper_right_pred_in[, 3],
    interval_width = pred_75 - pred_25,
    in_interval = DOY >= pred_25 & DOY <= pred_75
  )

df_upper_right_in_sample_metrics <- calculate_qrf_metrics(
  y = df_upper_right_in_sample$DOY,
  p25 = df_upper_right_in_sample$pred_25,
  p50 = df_upper_right_in_sample$pred_50,
  p75 = df_upper_right_in_sample$pred_75
) %>%
  mutate(
    train_region = train_region,
    test_region = "Upper_right_in_sample",
    n_test = nrow(df_upper_right_in_sample),
    .before = 1
  )

cat("\n=========================================\n")
cat("In-sample metrics: Upper-right-trained model on Upper_right data\n")
cat("=========================================\n")
print(df_upper_right_in_sample_metrics)

# =========================================
# Part 2 — Cross-region validation
# =========================================

test_regions <- c(
  "Lower_left",
  "Upper_left",
  "Lower_right",
  "Center"
)

upper_right_validation_results <- list()
upper_right_prediction_results <- list()

for (test_reg in test_regions) {

  cat("\n=========================================\n")
  cat("Training region:", train_region, "\n")
  cat("Testing region :", test_reg, "\n")
  cat("=========================================\n")

  df_test <- df_qrf %>%
    filter(region == test_reg) %>%
    drop_na(all_of(c(features, response)))

  pred <- predict(
    upper_right_model,
    newdata = df_test[, features],
    what = c(0.25, 0.50, 0.75)
  )

  df_pred <- df_test %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      pred_25 = pred[, 1],
      pred_50 = pred[, 2],
      pred_75 = pred[, 3],
      interval_width = pred_75 - pred_25,
      in_interval = DOY >= pred_25 & DOY <= pred_75
    )

  metric_df <- calculate_qrf_metrics(
    y = df_pred$DOY,
    p25 = df_pred$pred_25,
    p50 = df_pred$pred_50,
    p75 = df_pred$pred_75
  ) %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      n_test = nrow(df_pred),
      .before = 1
    )

  upper_right_validation_results[[test_reg]] <- metric_df
  upper_right_prediction_results[[test_reg]] <- df_pred

  print(metric_df)
}

df_upper_right_validation_metrics <- bind_rows(
  df_upper_right_in_sample_metrics,
  bind_rows(upper_right_validation_results)
)

df_upper_right_validation_predictions <- bind_rows(
  df_upper_right_in_sample,
  bind_rows(upper_right_prediction_results)
)

cat("\n✅ Upper-right-trained in-sample + cross-region validation completed!\n\n")
print(df_upper_right_validation_metrics)

# =========================================
# Block 5 — Upper-left-trained QRF:
# In-sample metrics first, then cross-region validation
# =========================================

# ---- Training model: Upper_left ----
train_region <- "Upper_left"
upper_left_model <- qrf_models[[train_region]]

# =========================================
# Part 1 — In-sample metrics for Upper_left
# =========================================

df_upper_left_train <- df_qrf %>%
  filter(region == train_region) %>%
  drop_na(all_of(c(features, response)))

upper_left_pred_in <- predict(
  upper_left_model,
  newdata = df_upper_left_train[, features],
  what = c(0.25, 0.50, 0.75)
)

df_upper_left_in_sample <- df_upper_left_train %>%
  mutate(
    train_region = train_region,
    test_region = "Upper_left_in_sample",
    pred_25 = upper_left_pred_in[, 1],
    pred_50 = upper_left_pred_in[, 2],
    pred_75 = upper_left_pred_in[, 3],
    interval_width = pred_75 - pred_25,
    in_interval = DOY >= pred_25 & DOY <= pred_75
  )

df_upper_left_in_sample_metrics <- calculate_qrf_metrics(
  y = df_upper_left_in_sample$DOY,
  p25 = df_upper_left_in_sample$pred_25,
  p50 = df_upper_left_in_sample$pred_50,
  p75 = df_upper_left_in_sample$pred_75
) %>%
  mutate(
    train_region = train_region,
    test_region = "Upper_left_in_sample",
    n_test = nrow(df_upper_left_in_sample),
    .before = 1
  )

cat("\n=========================================\n")
cat("In-sample metrics: Upper-left-trained model on Upper_left data\n")
cat("=========================================\n")
print(df_upper_left_in_sample_metrics)

# =========================================
# Part 2 — Cross-region validation
# =========================================

test_regions <- c(
  "Lower_left",
  "Lower_right",
  "Upper_right",
  "Center"
)

upper_left_validation_results <- list()
upper_left_prediction_results <- list()

for (test_reg in test_regions) {

  cat("\n=========================================\n")
  cat("Training region:", train_region, "\n")
  cat("Testing region :", test_reg, "\n")
  cat("=========================================\n")

  df_test <- df_qrf %>%
    filter(region == test_reg) %>%
    drop_na(all_of(c(features, response)))

  pred <- predict(
    upper_left_model,
    newdata = df_test[, features],
    what = c(0.25, 0.50, 0.75)
  )

  df_pred <- df_test %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      pred_25 = pred[, 1],
      pred_50 = pred[, 2],
      pred_75 = pred[, 3],
      interval_width = pred_75 - pred_25,
      in_interval = DOY >= pred_25 & DOY <= pred_75
    )

  metric_df <- calculate_qrf_metrics(
    y = df_pred$DOY,
    p25 = df_pred$pred_25,
    p50 = df_pred$pred_50,
    p75 = df_pred$pred_75
  ) %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      n_test = nrow(df_pred),
      .before = 1
    )

  upper_left_validation_results[[test_reg]] <- metric_df
  upper_left_prediction_results[[test_reg]] <- df_pred

  print(metric_df)
}

df_upper_left_validation_metrics <- bind_rows(
  df_upper_left_in_sample_metrics,
  bind_rows(upper_left_validation_results)
)

df_upper_left_validation_predictions <- bind_rows(
  df_upper_left_in_sample,
  bind_rows(upper_left_prediction_results)
)

cat("\n✅ Upper-left-trained in-sample + cross-region validation completed!\n\n")
print(df_upper_left_validation_metrics)

# =========================================
# Block 6 — Lower-right-trained QRF:
# In-sample metrics first, then cross-region validation
# =========================================

# ---- Training model: Lower_right ----
train_region <- "Lower_right"
lower_right_model <- qrf_models[[train_region]]

# =========================================
# Part 1 — In-sample metrics for Lower_right
# =========================================

df_lower_right_train <- df_qrf %>%
  filter(region == train_region) %>%
  drop_na(all_of(c(features, response)))

lower_right_pred_in <- predict(
  lower_right_model,
  newdata = df_lower_right_train[, features],
  what = c(0.25, 0.50, 0.75)
)

df_lower_right_in_sample <- df_lower_right_train %>%
  mutate(
    train_region = train_region,
    test_region = "Lower_right_in_sample",
    pred_25 = lower_right_pred_in[, 1],
    pred_50 = lower_right_pred_in[, 2],
    pred_75 = lower_right_pred_in[, 3],
    interval_width = pred_75 - pred_25,
    in_interval = DOY >= pred_25 & DOY <= pred_75
  )

df_lower_right_in_sample_metrics <- calculate_qrf_metrics(
  y = df_lower_right_in_sample$DOY,
  p25 = df_lower_right_in_sample$pred_25,
  p50 = df_lower_right_in_sample$pred_50,
  p75 = df_lower_right_in_sample$pred_75
) %>%
  mutate(
    train_region = train_region,
    test_region = "Lower_right_in_sample",
    n_test = nrow(df_lower_right_in_sample),
    .before = 1
  )

cat("\n=========================================\n")
cat("In-sample metrics: Lower-right-trained model on Lower_right data\n")
cat("=========================================\n")
print(df_lower_right_in_sample_metrics)

# =========================================
# Part 2 — Cross-region validation
# =========================================

test_regions <- c(
  "Lower_left",
  "Upper_left",
  "Upper_right",
  "Center"
)

lower_right_validation_results <- list()
lower_right_prediction_results <- list()

for (test_reg in test_regions) {

  cat("\n=========================================\n")
  cat("Training region:", train_region, "\n")
  cat("Testing region :", test_reg, "\n")
  cat("=========================================\n")

  df_test <- df_qrf %>%
    filter(region == test_reg) %>%
    drop_na(all_of(c(features, response)))

  pred <- predict(
    lower_right_model,
    newdata = df_test[, features],
    what = c(0.25, 0.50, 0.75)
  )

  df_pred <- df_test %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      pred_25 = pred[, 1],
      pred_50 = pred[, 2],
      pred_75 = pred[, 3],
      interval_width = pred_75 - pred_25,
      in_interval = DOY >= pred_25 & DOY <= pred_75
    )

  metric_df <- calculate_qrf_metrics(
    y = df_pred$DOY,
    p25 = df_pred$pred_25,
    p50 = df_pred$pred_50,
    p75 = df_pred$pred_75
  ) %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      n_test = nrow(df_pred),
      .before = 1
    )

  lower_right_validation_results[[test_reg]] <- metric_df
  lower_right_prediction_results[[test_reg]] <- df_pred

  print(metric_df)
}

df_lower_right_validation_metrics <- bind_rows(
  df_lower_right_in_sample_metrics,
  bind_rows(lower_right_validation_results)
)

df_lower_right_validation_predictions <- bind_rows(
  df_lower_right_in_sample,
  bind_rows(lower_right_prediction_results)
)

cat("\n✅ Lower-right-trained in-sample + cross-region validation completed!\n\n")
print(df_lower_right_validation_metrics)

# =========================================
# Block 7 — Lower-left-trained QRF:
# In-sample metrics first, then cross-region validation
# =========================================

# ---- Training model: Lower_left ----
train_region <- "Lower_left"
lower_left_model <- qrf_models[[train_region]]

# =========================================
# Part 1 — In-sample metrics for Lower_left
# =========================================

df_lower_left_train <- df_qrf %>%
  filter(region == train_region) %>%
  drop_na(all_of(c(features, response)))

lower_left_pred_in <- predict(
  lower_left_model,
  newdata = df_lower_left_train[, features],
  what = c(0.25, 0.50, 0.75)
)

df_lower_left_in_sample <- df_lower_left_train %>%
  mutate(
    train_region = train_region,
    test_region = "Lower_left_in_sample",
    pred_25 = lower_left_pred_in[, 1],
    pred_50 = lower_left_pred_in[, 2],
    pred_75 = lower_left_pred_in[, 3],
    interval_width = pred_75 - pred_25,
    in_interval = DOY >= pred_25 & DOY <= pred_75
  )

df_lower_left_in_sample_metrics <- calculate_qrf_metrics(
  y = df_lower_left_in_sample$DOY,
  p25 = df_lower_left_in_sample$pred_25,
  p50 = df_lower_left_in_sample$pred_50,
  p75 = df_lower_left_in_sample$pred_75
) %>%
  mutate(
    train_region = train_region,
    test_region = "Lower_left_in_sample",
    n_test = nrow(df_lower_left_in_sample),
    .before = 1
  )

cat("\n=========================================\n")
cat("In-sample metrics: Lower-left-trained model on Lower_left data\n")
cat("=========================================\n")
print(df_lower_left_in_sample_metrics)

# =========================================
# Part 2 — Cross-region validation
# =========================================

test_regions <- c(
  "Upper_left",
  "Lower_right",
  "Upper_right",
  "Center"
)

lower_left_validation_results <- list()
lower_left_prediction_results <- list()

for (test_reg in test_regions) {

  cat("\n=========================================\n")
  cat("Training region:", train_region, "\n")
  cat("Testing region :", test_reg, "\n")
  cat("=========================================\n")

  df_test <- df_qrf %>%
    filter(region == test_reg) %>%
    drop_na(all_of(c(features, response)))

  pred <- predict(
    lower_left_model,
    newdata = df_test[, features],
    what = c(0.25, 0.50, 0.75)
  )

  df_pred <- df_test %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      pred_25 = pred[, 1],
      pred_50 = pred[, 2],
      pred_75 = pred[, 3],
      interval_width = pred_75 - pred_25,
      in_interval = DOY >= pred_25 & DOY <= pred_75
    )

  metric_df <- calculate_qrf_metrics(
    y = df_pred$DOY,
    p25 = df_pred$pred_25,
    p50 = df_pred$pred_50,
    p75 = df_pred$pred_75
  ) %>%
    mutate(
      train_region = train_region,
      test_region = test_reg,
      n_test = nrow(df_pred),
      .before = 1
    )

  lower_left_validation_results[[test_reg]] <- metric_df
  lower_left_prediction_results[[test_reg]] <- df_pred

  print(metric_df)
}

df_lower_left_validation_metrics <- bind_rows(
  df_lower_left_in_sample_metrics,
  bind_rows(lower_left_validation_results)
)

df_lower_left_validation_predictions <- bind_rows(
  df_lower_left_in_sample,
  bind_rows(lower_left_prediction_results)
)

cat("\n✅ Lower-left-trained in-sample + cross-region validation completed!\n\n")
print(df_lower_left_validation_metrics)

# =========================================
# Block 8 — Full-sample QRF using four quadrants
# 10-fold cross-validation
# Metrics: ICP, MPIW, MAE, Pinball Loss
# =========================================

# ---- Use only four quadrant regions ----
quadrant_regions <- c(
  "Upper_right",
  "Upper_left",
  "Lower_right",
  "Lower_left"
)

df_full_10000 <- df_qrf %>%
  filter(region %in% quadrant_regions) %>%
  drop_na(all_of(c(features, response)))

cat("Full-sample size:", nrow(df_full_10000), "\n")
print(table(df_full_10000$region))

# ---- Metric functions ----
pinball_loss <- function(y, y_pred, tau) {
  err <- y - y_pred
  mean(ifelse(err >= 0, tau * err, (tau - 1) * err))
}

calculate_qrf_metrics <- function(y, p25, p50, p75) {
  data.frame(
    ICP = mean(y >= p25 & y <= p75),
    MPIW = mean(p75 - p25),
    MAE_25 = mean(abs(y - p25)),
    MAE_50 = mean(abs(y - p50)),
    MAE_75 = mean(abs(y - p75)),
    Pinball_25 = pinball_loss(y, p25, 0.25),
    Pinball_50 = pinball_loss(y, p50, 0.50),
    Pinball_75 = pinball_loss(y, p75, 0.75)
  )
}

# =========================================
# Part 1 — Full-sample in-sample QRF metrics
# =========================================

set.seed(42)

qrf_full_model <- quantregForest(
  x = df_full_10000[, features],
  y = df_full_10000[[response]],
  ntree = 1000
)

pred_full_in <- predict(
  qrf_full_model,
  newdata = df_full_10000[, features],
  what = c(0.25, 0.50, 0.75)
)

df_full_in_sample <- df_full_10000 %>%
  mutate(
    pred_25 = pred_full_in[, 1],
    pred_50 = pred_full_in[, 2],
    pred_75 = pred_full_in[, 3],
    interval_width = pred_75 - pred_25,
    in_interval = DOY >= pred_25 & DOY <= pred_75
  )

df_full_in_sample_metrics <- calculate_qrf_metrics(
  y = df_full_in_sample$DOY,
  p25 = df_full_in_sample$pred_25,
  p50 = df_full_in_sample$pred_50,
  p75 = df_full_in_sample$pred_75
) %>%
  mutate(
    stage = "Full_sample_in_sample",
    n = nrow(df_full_in_sample),
    .before = 1
  )

cat("\n=========================================\n")
cat("Full-sample in-sample metrics\n")
cat("=========================================\n")
print(df_full_in_sample_metrics)

# =========================================
# Part 2 — 10-fold cross-validation
# =========================================

set.seed(42)

n <- nrow(df_full_10000)

df_full_10000 <- df_full_10000 %>%
  mutate(
    fold = sample(rep(1:10, length.out = n))
  )

cv_metrics_list <- list()
cv_predictions_list <- list()

for (k in 1:10) {

  cat("\n=========================================\n")
  cat("10-fold CV — Fold", k, "\n")
  cat("=========================================\n")

  df_train <- df_full_10000 %>%
    filter(fold != k)

  df_test <- df_full_10000 %>%
    filter(fold == k)

  cat("Train n:", nrow(df_train), "\n")
  cat("Test n :", nrow(df_test), "\n")

  qrf_cv_model <- quantregForest(
    x = df_train[, features],
    y = df_train[[response]],
    ntree = 1000
  )

  pred_cv <- predict(
    qrf_cv_model,
    newdata = df_test[, features],
    what = c(0.25, 0.50, 0.75)
  )

  df_cv_pred <- df_test %>%
    mutate(
      pred_25 = pred_cv[, 1],
      pred_50 = pred_cv[, 2],
      pred_75 = pred_cv[, 3],
      interval_width = pred_75 - pred_25,
      in_interval = DOY >= pred_25 & DOY <= pred_75
    )

  fold_metrics <- calculate_qrf_metrics(
    y = df_cv_pred$DOY,
    p25 = df_cv_pred$pred_25,
    p50 = df_cv_pred$pred_50,
    p75 = df_cv_pred$pred_75
  ) %>%
    mutate(
      stage = "10_fold_CV",
      fold = k,
      n = nrow(df_cv_pred),
      .before = 1
    )

  cv_metrics_list[[k]] <- fold_metrics
  cv_predictions_list[[k]] <- df_cv_pred

  print(fold_metrics)
}

df_full_cv_metrics_by_fold <- bind_rows(cv_metrics_list)
df_full_cv_predictions <- bind_rows(cv_predictions_list)

# ---- CV mean ± SD ----
df_full_cv_summary <- df_full_cv_metrics_by_fold %>%
  summarise(
    stage = "10_fold_CV_mean_sd",
    n_total = sum(n),

    ICP_mean = mean(ICP),
    ICP_sd = sd(ICP),

    MPIW_mean = mean(MPIW),
    MPIW_sd = sd(MPIW),

    MAE_25_mean = mean(MAE_25),
    MAE_25_sd = sd(MAE_25),
    MAE_50_mean = mean(MAE_50),
    MAE_50_sd = sd(MAE_50),
    MAE_75_mean = mean(MAE_75),
    MAE_75_sd = sd(MAE_75),

    Pinball_25_mean = mean(Pinball_25),
    Pinball_25_sd = sd(Pinball_25),
    Pinball_50_mean = mean(Pinball_50),
    Pinball_50_sd = sd(Pinball_50),
    Pinball_75_mean = mean(Pinball_75),
    Pinball_75_sd = sd(Pinball_75)
  )

cat("\n=========================================\n")
cat("10-fold CV summary: mean ± SD across folds\n")
cat("=========================================\n")
print(df_full_cv_summary)

# ---- Also calculate pooled CV metrics across all held-out predictions ----
df_full_cv_pooled_metrics <- calculate_qrf_metrics(
  y = df_full_cv_predictions$DOY,
  p25 = df_full_cv_predictions$pred_25,
  p50 = df_full_cv_predictions$pred_50,
  p75 = df_full_cv_predictions$pred_75
) %>%
  mutate(
    stage = "10_fold_CV_pooled",
    n = nrow(df_full_cv_predictions),
    .before = 1
  )

cat("\n=========================================\n")
cat("10-fold CV pooled metrics across all held-out predictions\n")
cat("=========================================\n")
print(df_full_cv_pooled_metrics)

# ---- Final compact comparison table ----
df_full_model_comparison <- bind_rows(
  df_full_in_sample_metrics,
  df_full_cv_pooled_metrics
)

cat("\n=========================================\n")
cat("Final comparison: in-sample vs 10-fold CV pooled\n")
cat("=========================================\n")
print(df_full_model_comparison)

# =========================================
# Block 9 — Histogram of DOY distributions across five sampled regions
# Bin width = 1 day; overlapping transparent histograms
# Center plotted first so it stays in the back layer
# =========================================

library(ggplot2)
library(dplyr)

# ---- Use five sampled regions ----
plot_regions <- c(
  "Center",
  "Lower_left",
  "Upper_left",
  "Lower_right",
  "Upper_right"
)

df_hist_5region <- df_qrf %>%
  filter(region %in% plot_regions) %>%
  mutate(
    region = factor(
      region,
      levels = c("Center", "Lower_left", "Upper_left", "Lower_right", "Upper_right")
    )
  )

# ---- Plot ----
p_hist_5region <- ggplot(df_hist_5region, aes(x = DOY, fill = region)) +
  geom_histogram(
    binwidth = 1,
    position = "identity",
    alpha = 0.35,
    color = "black",
    linewidth = 0.2
  ) +
  labs(
    x = "DOY",
    y = "Frequency",
    fill = "Sampling region"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_hist_5region)

# ---- Export ----
ggsave(
  "DOY_histogram_five_sampling_regions.png",
  plot = p_hist_5region,
  width = 10,
  height = 7,
  dpi = 500
)

ggsave(
  "DOY_histogram_five_sampling_regions.pdf",
  plot = p_hist_5region,
  width = 10,
  height = 7
)

cat("✅ Five-region DOY histogram saved as PNG and PDF.\n")

# ======================================================
# Figure: Simulated climatic-state grid +
# conditioned phenological distributions
# Output: PDF + 500 dpi PNG
# ======================================================

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("ggridges")) install.packages("ggridges", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("patchwork")) install.packages("patchwork", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("ggnewscale")) install.packages("ggnewscale", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("viridis")) install.packages("viridis", repos = "https://cran.csie.ntu.edu.tw/")

library(ggplot2)
library(dplyr)
library(ggridges)
library(patchwork)
library(ggnewscale)
library(viridis)

base_font <- "Arial"

# ======================================================
# Shared colours
# ======================================================

region_cols <- c(
  "Center" = "#D55E00",
  "Lower_left" = "#E69F00",
  "Upper_left" = "#009E73",
  "Lower_right" = "#56B4E9",
  "Upper_right" = "#CC79A7"
)

# ======================================================
# Panel a: Simulated climatic-state grid
# ======================================================

region_overlay <- data.frame(
  region = c(
    "Lower_left",
    "Upper_left",
    "Lower_right",
    "Upper_right",
    "Center"
  ),
  xmin = c(0.5, 0.5, 50.5, 50.5, 24.5),
  xmax = c(50.5, 50.5, 100.5, 100.5, 75.5),
  ymin = c(0.5, 50.5, 0.5, 50.5, 24.5),
  ymax = c(50.5, 100.5, 50.5, 100.5, 75.5)
)

p_heatmap <- ggplot(
  df_qrf,
  aes(x = x, y = y)
) +
  geom_tile(
    aes(fill = DOY)
  ) +
  scale_fill_gradient(
  low = "white",
  high = "grey40",
  name = "DOY"
  ) +

  ggnewscale::new_scale_fill() +

  geom_rect(
    data = region_overlay,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = region
    ),
    inherit.aes = FALSE,
    alpha = 0.22,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = region_cols,
    name = "Sampling region",
    breaks = c(
      "Center",
      "Lower_left",
      "Upper_left",
      "Lower_right",
      "Upper_right"
    ),
    labels = c(
      "Center",
      "Lower left",
      "Upper left",
      "Lower right",
      "Upper right"
    )
  ) +

  coord_equal(
    xlim = c(0.5, 100.5),
    ylim = c(0.5, 100.5),
    expand = FALSE
  ) +

  labs(
    x = "X grid position",
    y = "Y grid position"
  ) +

  theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),

    panel.grid = element_blank()
  )

# ======================================================
# Panel b: Conditioned phenological distributions
# ======================================================

ridge_order <- c(
  "Upper_right",
  "Upper_left",
  "Center",
  "Lower_right",
  "Lower_left"
)

df_ridge <- df_qrf %>%
  filter(region %in% ridge_order) %>%
  mutate(
    region = factor(
      region,
      levels = ridge_order
    )
  )

p_ridge <- ggplot(
  df_ridge,
  aes(
    x = DOY,
    y = region,
    fill = region
  )
) +
  geom_density_ridges(
    alpha = 0.85,
    scale = 2.4,
    colour = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = region_cols,
    guide = "none"
  ) +
  scale_y_discrete(
    labels = c(
      "Upper_right" = "Upper right",
      "Upper_left" = "Upper left",
      "Center" = "Center",
      "Lower_right" = "Lower right",
      "Lower_left" = "Lower left"
    )
  ) +
  labs(
    x = "DOY",
    y = NULL
  ) +
  theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    legend.position = "none",

    panel.grid.minor = element_blank()
  )

# ======================================================
# Combine panels
# ======================================================

p_final <- (
  p_heatmap /
  p_ridge
) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(
      family = base_font,
      face = "bold",
      size = 24
    )
  )

print(p_final)

# ======================================================
# Export
# ======================================================

ggsave(
  filename = "Figure_Simulation_Grid_Ridge.png",
  plot = p_final,
  width = 10,
  height = 12,
  dpi = 500
)

ggsave(
  filename = "Figure_Simulation_Grid_Ridge.pdf",
  plot = p_final,
  width = 10,
  height = 12
)

# ======================================================
# Figure: Simulated climatic-state grid +
# conditioned phenological distributions
# Output: PDF + 500 dpi PNG
# ======================================================

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("ggridges")) install.packages("ggridges", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("patchwork")) install.packages("patchwork", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("ggnewscale")) install.packages("ggnewscale", repos = "https://cran.csie.ntu.edu.tw/")

library(ggplot2)
library(dplyr)
library(ggridges)
library(patchwork)
library(ggnewscale)

base_font <- "Arial"

region_cols <- c(
  "Center" = "#D55E00",
  "Lower_left" = "#E69F00",
  "Upper_left" = "#009E73",
  "Lower_right" = "#56B4E9",
  "Upper_right" = "#CC79A7"
)

# ======================================================
# Panel a: Simulated climatic-state grid
# ======================================================

region_overlay <- data.frame(
  region = c(
    "Lower_left",
    "Upper_left",
    "Lower_right",
    "Upper_right",
    "Center"
  ),
  xmin = c(0.5, 0.5, 50.5, 50.5, 24.5),
  xmax = c(50.5, 50.5, 100.5, 100.5, 75.5),
  ymin = c(0.5, 50.5, 0.5, 50.5, 24.5),
  ymax = c(50.5, 100.5, 50.5, 100.5, 75.5)
)

p_heatmap <- ggplot(
  df_qrf,
  aes(x = x, y = y)
) +
  geom_tile(
    aes(fill = DOY)
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "grey35",
    name = "DOY"
  ) +

  ggnewscale::new_scale_fill() +

  geom_rect(
    data = region_overlay,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = region
    ),
    inherit.aes = FALSE,
    alpha = 0.35,
    colour = NA
  ) +
  scale_fill_manual(
    values = region_cols,
    name = "Sampling region",
    breaks = c(
      "Center",
      "Lower_left",
      "Upper_left",
      "Lower_right",
      "Upper_right"
    ),
    labels = c(
      "Center",
      "Lower left",
      "Upper left",
      "Lower right",
      "Upper right"
    )
  ) +
  coord_equal(
    xlim = c(0.5, 100.5),
    ylim = c(0.5, 100.5),
    expand = FALSE
  ) +
  labs(
    x = "X grid position",
    y = "Y grid position"
  ) +
  theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),

    panel.grid = element_blank()
  )

# ======================================================
# Panel b: Conditioned phenological distributions
# ======================================================

ridge_order <- c(
  "Upper_right",
  "Upper_left",
  "Center",
  "Lower_right",
  "Lower_left"
)

df_ridge <- df_qrf %>%
  filter(region %in% ridge_order) %>%
  mutate(
    region = factor(region, levels = ridge_order)
  )

p_ridge <- ggplot(
  df_ridge,
  aes(
    x = DOY,
    y = region,
    fill = region
  )
) +
  geom_density_ridges(
    alpha = 0.85,
    scale = 2.4,
    colour = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = region_cols,
    guide = "none"
  ) +
  scale_y_discrete(
    labels = c(
      "Upper_right" = "Upper right",
      "Upper_left" = "Upper left",
      "Center" = "Center",
      "Lower_right" = "Lower right",
      "Lower_left" = "Lower left"
    )
  ) +
  labs(
    x = "DOY",
    y = NULL
  ) +
  theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ======================================================
# Combine panels
# ======================================================

p_final <- (
  p_heatmap /
  p_ridge
) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = base_font,
      face = "bold",
      size = 24
    )
  )

print(p_final)

# ======================================================
# Export
# ======================================================

ggsave(
  filename = "Figure_Simulation_Grid_Ridge.png",
  plot = p_final,
  width = 10,
  height = 12,
  dpi = 500
)

ggsave(
  filename = "Figure_Simulation_Grid_Ridge.pdf",
  plot = p_final,
  width = 10,
  height = 12
)

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("ggridges")) install.packages("ggridges", repos = "https://cran.csie.ntu.edu.tw/")
if (!require("patchwork")) install.packages("patchwork", repos = "https://cran.csie.ntu.edu.tw/")

library(ggplot2)
library(dplyr)
library(ggridges)
library(patchwork)

base_font <- "Arial"

region_cols <- c(
  "Center" = "#D55E00",
  "Lower_left" = "#E69F00",
  "Upper_left" = "#009E73",
  "Lower_right" = "#56B4E9",
  "Upper_right" = "#CC79A7"
)

region_overlay_climate <- data.frame(
  region = c(
    "Lower_left",
    "Upper_left",
    "Lower_right",
    "Upper_right",
    "Center"
  ),
  xmin = c(15.0, 15.0, 25.0, 25.0, 19.8),
  xmax = c(24.8, 24.8, 35.0, 35.0, 29.8),
  ymin = c(1.0, 3.0, 1.0, 3.0, 1.96),
  ymax = c(2.96, 5.0, 2.96, 5.0, 3.96)
)

# ======================================================
# Panel a
# ======================================================

p_state <- ggplot() +
  geom_rect(
    data = region_overlay_climate,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = region
    ),
    alpha = 0.38,
    colour = NA
  ) +
  scale_fill_manual(
    values = region_cols,
    name = "Sampling region",
    breaks = c(
      "Center",
      "Lower_left",
      "Upper_left",
      "Lower_right",
      "Upper_right"
    ),
    labels = c(
      "Center",
      "Lower left",
      "Upper left",
      "Lower right",
      "Upper right"
    )
  ) +
  scale_x_continuous(
    limits = c(15, 35),
    breaks = seq(15, 35, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(1, 5),
    breaks = seq(1, 5, by = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Mean temperature (°C)",
    y = "Standard deviation"
  ) +
  theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),

    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
  )

# ======================================================
# Panel b
# ======================================================

ridge_order <- c(
  "Upper_right",
  "Upper_left",
  "Center",
  "Lower_right",
  "Lower_left"
)

df_ridge <- df_qrf %>%
  filter(region %in% ridge_order) %>%
  mutate(
    region = factor(region, levels = ridge_order)
  )

p_ridge <- ggplot(
  df_ridge,
  aes(
    x = DOY,
    y = region,
    fill = region
  )
) +
  geom_density_ridges(
    alpha = 0.85,
    scale = 2.4,
    colour = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = region_cols,
    guide = "none"
  ) +
  scale_y_discrete(
    labels = c(
      "Upper_right" = "Upper right",
      "Upper_left" = "Upper left",
      "Center" = "Center",
      "Lower_right" = "Lower right",
      "Lower_left" = "Lower left"
    )
  ) +
  labs(
    x = "DOY",
    y = NULL
  ) +
  theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
  )

# ======================================================
# Combine panels
# free() prevents panel b's long y-axis labels from pulling panel a left
# ======================================================

p_final <- (
  patchwork::free(p_state, side = "l") /
  p_ridge
) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = base_font,
      face = "bold",
      size = 24
    )
  )

print(p_final)

# ======================================================
# Export
# ======================================================

ggsave(
  filename = "Figure_Climatic_State_Ridge.png",
  plot = p_final,
  width = 10,
  height = 12,
  dpi = 500
)

ggsave(
  filename = "Figure_Climatic_State_Ridge.pdf",
  plot = p_final,
  width = 10,
  height = 12
)
