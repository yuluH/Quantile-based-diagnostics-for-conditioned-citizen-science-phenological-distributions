# Quantile Random Forest analysis of citizen-science phenological observations
# Run from the project root in a local R session or RStudio project.
# Required input: CSCI_Features_for_R_QRF.csv
# Additional figure inputs: Hc.csv, Zm.csv, Fl.csv, Mf.csv
# The cross-species comparison figures below contain manually compiled values.

# ============================================================
# Package setup
# ============================================================

required_packages <- c(
  "quantregForest",
  "ggplot2",
  "tidyr",
  "dplyr",
  "sf",
  "blockCV",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}

cat("All required R packages are available.\n")

required_packages <- c("quantregForest", "ggplot2", "tidyr", "dplyr", "sf", "blockCV", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install the required packages before running: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

input_file <- "CSCI_Features_for_R_QRF.csv"
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, "\nWorking directory: ", getwd())
}
df_qrf <- read.csv(input_file, stringsAsFactors = FALSE)
required_columns <- c("DOY", "time_observed_at", "latitude", "longitude")
missing_columns <- setdiff(required_columns, names(df_qrf))
if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))
df_qrf$time_observed_at <- as.Date(df_qrf$time_observed_at)
df_qrf$DOY <- as.numeric(df_qrf$DOY)
cat("Input rows:", nrow(df_qrf), "\n")

# Descriptive statistics for the observed day of year.
doy_stats <- data.frame(
  Statistic = c("Maximum DOY", "75th-percentile DOY", "Mean DOY", "Median DOY", "25th-percentile DOY", "Minimum DOY"),
  Value = c(max(df_qrf$DOY, na.rm = TRUE), quantile(df_qrf$DOY, .75, na.rm = TRUE),
            mean(df_qrf$DOY, na.rm = TRUE), median(df_qrf$DOY, na.rm = TRUE),
            quantile(df_qrf$DOY, .25, na.rm = TRUE), min(df_qrf$DOY, na.rm = TRUE))
)
print(doy_stats)

# Fit the full-sample quantile regression forest.
features <- setdiff(names(df_qrf), c("ID", "latitude", "longitude", "time_observed_at", "DOY"))
if (!length(features)) stop("No predictor columns found.")
df_modeling <- tidyr::drop_na(df_qrf, DOY, dplyr::all_of(features))
set.seed(42)
qrf_model <- quantregForest(x = df_modeling[, features, drop = FALSE], y = df_modeling$DOY, ntree = 1000)
quantiles_to_predict <- c(.25, .50, .75)
qrf_predictions <- predict(qrf_model, df_modeling[, features, drop = FALSE], what = quantiles_to_predict)
df_results <- cbind(df_modeling, qrf_predictions)
names(df_results)[(ncol(df_results) - 2):ncol(df_results)] <- c("pred_25", "pred_50", "pred_75")
print(summary(df_results[, c("DOY", "pred_25", "pred_50", "pred_75")]))
res_50 <- df_results$DOY - df_results$pred_50
r2_median <- 1 - sum(res_50^2) / sum((df_results$DOY - mean(df_results$DOY))^2)
cat("In-sample median-based R-squared:", round(r2_median, 3), "\n")

# Evaluate interval coverage, interval width, absolute error, and pinball loss.
pinball_loss <- function(y, y_pred, tau) {
  err <- y - y_pred
  return(mean(ifelse(err >= 0, tau * err, (tau - 1) * err)))
}

icp <- sum(df_results$DOY >= df_results$pred_25 & df_results$DOY <= df_results$pred_75) / nrow(df_results)

mpiw <- mean(df_results$pred_75 - df_results$pred_25)

loss_25 <- pinball_loss(df_results$DOY, df_results$pred_25, 0.25)
loss_50 <- pinball_loss(df_results$DOY, df_results$pred_50, 0.50)
loss_75 <- pinball_loss(df_results$DOY, df_results$pred_75, 0.75)

cat(" Interval Coverage Probability (ICP, target 0.5):", round(icp, 3), "\n")
cat(" Mean Prediction Interval Width (MPIW):", round(mpiw, 3), " days\n")
cat(" Pinball Loss (25%):", round(loss_25, 3), "\n")
cat(" Pinball Loss (50%):", round(loss_50, 3), "\n")
cat(" Pinball Loss (75%):", round(loss_75, 3), "\n\n")

calculate_metrics_set <- function(y, p25, p50, p75) {
  # 1. ICP (target 0.5 for [25, 75] interval)
  icp_val <- sum(y >= p25 & y <= p75) / length(y)
  # 2. MPIW (Resolution)
  mpiw_val <- mean(p75 - p25)
  # 3. MAE (Accuracy for each quantile)
  mae25 <- mean(abs(y - p25))
  mae50 <- mean(abs(y - p50))
  mae75 <- mean(abs(y - p75))
  # 4. Pinball Loss (Quantile Loss)
  pl25 <- pinball_loss(y, p25, 0.25)
  pl50 <- pinball_loss(y, p50, 0.50)
  pl75 <- pinball_loss(y, p75, 0.75)

  return(c(icp = icp_val, mpiw = mpiw_val,
           mae25 = mae25, mae50 = mae50, mae75 = mae75,
           pl25 = pl25, pl50 = pl50, pl75 = pl75))
}

# ----------------------------------------------------------------------
in_sample_m <- calculate_metrics_set(df_results$DOY,
                                     df_results$pred_25,
                                     df_results$pred_50,
                                     df_results$pred_75)

cat(" ICP (target 0.5):", round(in_sample_m["icp"], 3), "\n")
cat(" MPIW:", round(in_sample_m["mpiw"], 3), " days\n")
cat(" MAE (25%, 50%, 75%):",
    round(in_sample_m["mae25"], 2), ",",
    round(in_sample_m["mae50"], 2), ",",
    round(in_sample_m["mae75"], 2), "\n")
cat(" Pinball Loss (25%, 50%, 75%):",
    round(in_sample_m["pl25"], 3), ",",
    round(in_sample_m["pl50"], 3), ",",
    round(in_sample_m["pl75"], 3), "\n\n")

# ----------------------------------------------------------------------

set.seed(42)
folds <- sample(rep(1:5, length.out = nrow(df_modeling)))
cv_metrics_list <- list()

for(i in 1:5) {
  train_data <- df_modeling[folds != i, ]
  test_data  <- df_modeling[folds == i, ]

  tmp_qrf <- quantregForest(x = train_data[, features], y = train_data$DOY, ntree = 500)

  tmp_pred <- predict(tmp_qrf, test_data[, features], what = c(0.25, 0.50, 0.75))

  cv_metrics_list[[i]] <- calculate_metrics_set(test_data$DOY,
                                                tmp_pred[,1],
                                                tmp_pred[,2],
                                                tmp_pred[,3])
}

cv_results_df <- as.data.frame(do.call(rbind, cv_metrics_list))
cv_mean <- colMeans(cv_results_df)
cv_sd   <- apply(cv_results_df, 2, sd)

cat(" ICP:", round(cv_mean["icp"], 3), "±", round(cv_sd["icp"], 3), "\n")
cat(" MPIW:", round(cv_mean["mpiw"], 2), "±", round(cv_sd["mpiw"], 2), " days\n")
cat(" MAE (25%):", round(cv_mean["mae25"], 2), "±", round(cv_sd["mae25"], 2), "\n")
cat(" MAE (50%):", round(cv_mean["mae50"], 2), "±", round(cv_sd["mae50"], 2), "\n")
cat(" MAE (75%):", round(cv_mean["mae75"], 2), "±", round(cv_sd["mae75"], 2), "\n")
cat(" Pinball Loss (25%):", round(cv_mean["pl25"], 3), "±", round(cv_sd["pl25"], 3), "\n")
cat(" Pinball Loss (50%):", round(cv_mean["pl50"], 3), "±", round(cv_sd["pl50"], 3), "\n")
cat(" Pinball Loss (75%):", round(cv_mean["pl75"], 3), "±", round(cv_sd["pl75"], 3), "\n")




# ======================================================
# Part 2B: Spatial block cross-validation (SCV)
# ======================================================


df_sf <- st_as_sf(
  df_modeling,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

df_sf <- st_transform(df_sf, 3826)

set.seed(42)

scv_obj <- cv_spatial(
  x = df_sf,
  column = NULL,
  k = 5,
  size = 50000,
  selection = "random",
  iteration = 100,
  progress = TRUE
)

scv_metrics_list <- list()

for(i in seq_along(scv_obj$folds_list)) {

  train_id <- scv_obj$folds_list[[i]][[1]]
  test_id  <- scv_obj$folds_list[[i]][[2]]

  train_data <- df_modeling[train_id, ]
  test_data  <- df_modeling[test_id, ]

  tmp_qrf <- quantregForest(
    x = train_data[, features],
    y = train_data$DOY,
    ntree = 500
  )

  tmp_pred <- predict(
    tmp_qrf,
    test_data[, features],
    what = c(0.25, 0.50, 0.75)
  )

  scv_metrics_list[[i]] <- calculate_metrics_set(
    test_data$DOY,
    tmp_pred[, 1],
    tmp_pred[, 2],
    tmp_pred[, 3]
  )
}

scv_results_df <- as.data.frame(do.call(rbind, scv_metrics_list))
scv_mean <- colMeans(scv_results_df)
scv_sd   <- apply(scv_results_df, 2, sd)

cat(" ICP:", round(scv_mean["icp"], 3), "±", round(scv_sd["icp"], 3), "\n")
cat(" MPIW:", round(scv_mean["mpiw"], 2), "±", round(scv_sd["mpiw"], 2), " days\n")
cat(" MAE (25%):", round(scv_mean["mae25"], 2), "±", round(scv_sd["mae25"], 2), "\n")
cat(" MAE (50%):", round(scv_mean["mae50"], 2), "±", round(scv_sd["mae50"], 2), "\n")
cat(" MAE (75%):", round(scv_mean["mae75"], 2), "±", round(scv_sd["mae75"], 2), "\n")
cat(" Pinball Loss (25%):", round(scv_mean["pl25"], 3), "±", round(scv_sd["pl25"], 3), "\n")
cat(" Pinball Loss (50%):", round(scv_mean["pl50"], 3), "±", round(scv_sd["pl50"], 3), "\n")
cat(" Pinball Loss (75%):", round(scv_mean["pl75"], 3), "±", round(scv_sd["pl75"], 3), "\n")



