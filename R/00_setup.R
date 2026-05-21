# 00_setup.R — Package installation, library loading, seed, helper functions
cat("\n=== SADRA Model Analysis: Setup ===\n")

# ── Package installation ──────────────────────────────────────────────
required_pkgs <- c(
  "mice", "glmnet", "Boruta", "JMbayes2", "survival", "survminer",
  "ggplot2", "dplyr", "tidyr", "tibble", "readr", "purrr", "stringr",
  "timeROC", "riskRegression", "pec", "boot", "rms", "pROC",
  "randomForestSRC", "CoxBoost", "glmnet", "DoubleML", "mlr3", "mlr3learners",
  "shiny", "shinythemes", "plotly", "DT", "rmarkdown", "knitr",
  "ggalluvial", "UpSetR", "scales", "cowplot",
  "future", "future.apply", "doParallel", "foreach",
  "Cairo", "showtext", "sysfonts", "splines"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing package: %s", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in required_pkgs) {
  tryCatch(
    install_if_missing(pkg),
    error = function(e) warning(sprintf("Could not load '%s': %s", pkg, e$message))
  )
}

# ── Global settings ───────────────────────────────────────────────────
set.seed(2025)
options(scipen = 999)
options(future.globals.maxSize = 2 * 1024^3)
plan(sequential)  # safer default; switch to multisession for MCMC chains only when needed

# ── Helper functions ──────────────────────────────────────────────────

#' Compute integrated AUC via trapezoidal rule
#' @param auc_df data.frame with columns "time" and "AUC"
compute_auc_jm <- function(auc_df) {
  stopifnot(all(c("time", "AUC") %in% colnames(auc_df)))
  n <- nrow(auc_df)
  if (n < 2) return(NA_real_)
  sum(diff(auc_df$time) * (auc_df$AUC[-1] + auc_df$AUC[-n]) / 2)
}

#' Calibration plot for predicted vs observed risk
#' @param predicted numeric vector of predicted probabilities
#' @param observed   numeric vector of observed events (0/1)
#' @param n_groups   number of risk groups (default 10)
#' @param title      plot title
plot_calibration <- function(predicted, observed, n_groups = 10,
                             title = "Calibration Plot") {
  df <- data.frame(pred = predicted, obs = observed)
  df <- df[order(df$pred), ]
  df$group <- cut(seq_len(nrow(df)), breaks = n_groups, labels = FALSE)
  cal <- aggregate(cbind(pred, obs) ~ group, data = df, mean)

  p <- ggplot(cal, aes(x = pred, y = obs)) +
    geom_point(size = 3, color = "#2166AC") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "loess", se = TRUE, color = "#B2182B", fill = "#B2182B", alpha = 0.15) +
    labs(x = "Predicted Risk", y = "Observed Proportion",
         title = title) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 13)
  return(p)
}

#' Cluster bootstrap AUC by patient ID
#' Resamples patient IDs with replacement; for each selected patient,
#' all longitudinal observations are kept.
#' @param predictions numeric vector of predicted risks (one per patient)
#' @param outcomes    numeric vector of observed events (0/1, one per patient)
#' @param patient_ids patient ID vector
#' @param n_boot      number of bootstrap resamples
#' @param times       evaluation time points (for naming only)
cluster_bootstrap_auc <- function(predictions, outcomes, patient_ids,
                                  n_boot = 1000) {
  unique_ids <- unique(patient_ids)
  n_patients <- length(unique_ids)
  auc_boot <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    sampled_ids <- sample(unique_ids, n_patients, replace = TRUE)
    boot_idx <- which(patient_ids %in% sampled_ids)
    pred_b <- predictions[boot_idx]
    obs_b  <- outcomes[boot_idx]
    auc_boot[b] <- tryCatch(
      as.numeric(pROC::roc(obs_b, pred_b, quiet = TRUE)$auc),
      error = function(e) NA_real_
    )
  }

  auc_boot <- auc_boot[!is.na(auc_boot)]
  list(
    AUC       = mean(auc_boot, na.rm = TRUE),
    AUC_SE    = sd(auc_boot, na.rm = TRUE) / sqrt(length(auc_boot)),
    AUC_lower = quantile(auc_boot, 0.025, na.rm = TRUE),
    AUC_upper = quantile(auc_boot, 0.975, na.rm = TRUE),
    auc_dist  = auc_boot
  )
}

#' Standard error formatting
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

cat("Setup complete. Packages loaded, seed = 2025.\n")
