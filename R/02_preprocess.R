# 02_preprocess.R — Exclude high-missing variables, MICE imputation, long-format SOFA
cat("\n=== Step 2: Preprocessing & Imputation ===\n")

# ── Load raw data ─────────────────────────────────────────────────────
train <- readRDS("data/training_cohort.rds")
prosp <- readRDS("data/prospective_cohort.rds")
trace <- readRDS("data/trace_cohort.rds")

# ── Detect daily SOFA column pattern ──────────────────────────────────
sofa_cols_all <- grep("^SOFA_day\\d+$", colnames(train), value = TRUE)
if (length(sofa_cols_all) < 3) {
  stop("FATAL: Fewer than 3 daily SOFA columns found. Need SOFA_day1 through at least SOFA_day7.")
}
sofa_days <- as.numeric(gsub("SOFA_day", "", sofa_cols_all))
sofa_days <- sort(sofa_days)
message(sprintf("Daily SOFA columns: SOFA_day%d through SOFA_day%d (%d time points)",
                min(sofa_days), max(sofa_days), length(sofa_days)))

# ── Exclude variables with >=20% missing ───────────────────────────────
# Baseline variables only (not SOFA daily scores — those are kept as NA
# and handled by the joint model's longitudinal submodel)
baseline_cols <- setdiff(colnames(train), sofa_cols_all)

exclude_high_missing <- function(df, threshold = 0.2) {
  # Only evaluate baseline columns (not daily SOFA, not ID/outcomes)
  eval_cols <- setdiff(colnames(df), c(sofa_cols_all, "ID", "time_to_IPN", "IPN",
                                       "IPN_strict", "treatment", "Antibiotic_actual"))
  eval_cols <- intersect(eval_cols, colnames(df))
  miss_pct <- colMeans(is.na(df[, eval_cols, drop = FALSE]))
  drop_cols <- names(miss_pct[miss_pct >= threshold])
  if (length(drop_cols) > 0) {
    message(sprintf("  Excluding %d variables with >=%.0f%% missing:", length(drop_cols), threshold * 100))
    for (col in drop_cols) {
      message(sprintf("    - %s (%.1f%%)", col, miss_pct[col] * 100))
    }
    df <- df[, setdiff(colnames(df), drop_cols), drop = FALSE]
  } else {
    message("  No baseline variables exceed missing threshold.")
  }
  return(df)
}

message("\nTraining cohort:")
train <- exclude_high_missing(train)
message("\nProspective cohort:")
prosp <- exclude_high_missing(prosp)
message("\nTRACE cohort:")
trace <- exclude_high_missing(trace)

# ── Align columns across datasets ──────────────────────────────────────
common_cols <- Reduce(intersect, list(colnames(train), colnames(prosp), colnames(trace)))
# Always include daily SOFA columns even if not perfectly aligned
common_cols <- unique(c(common_cols, intersect(sofa_cols_all, colnames(train)),
                        intersect(sofa_cols_all, colnames(prosp)),
                        intersect(sofa_cols_all, colnames(trace))))
train <- train[, common_cols, drop = FALSE]
prosp <- prosp[, common_cols, drop = FALSE]
trace <- trace[, common_cols, drop = FALSE]
message(sprintf("\nCommon columns retained: %d", length(common_cols)))

# Update sofa_cols_all after column alignment
sofa_cols_all <- grep("^SOFA_day\\d+$", common_cols, value = TRUE)

# ── Ensure binary outcomes ─────────────────────────────────────────────
for (df_name in c("train", "prosp", "trace")) {
  df <- get(df_name)
  if ("IPN" %in% colnames(df)) {
    df$IPN <- as.integer(as.numeric(df$IPN) > 0)
  }
  if ("IPN_strict" %in% colnames(df)) {
    df$IPN_strict <- as.integer(as.numeric(df$IPN_strict) > 0)
  }
  # Create event indicator for survival analysis
  # Composite definition: IPN = 1 as event, 0 as censored
  df$event_type <- df$IPN
  assign(df_name, df)
}

# ── Time to event: censored at 90 days ────────────────────────────────
max_time <- 90
for (df_name in c("train", "prosp", "trace")) {
  df <- get(df_name)
  if ("time_to_IPN" %in% colnames(df)) {
    df$time_to_IPN <- pmin(as.numeric(df$time_to_IPN), max_time)
    # Ensure minimum positive time for survival analysis
    df$time_to_IPN <- pmax(df$time_to_IPN, 0.1)
  }
  # Alias for convenience in joint model
  df$event_time <- df$time_to_IPN
  assign(df_name, df)
}

# ── Multiple imputation with mice (m = 5) on training set ─────────────
message("\nRunning MICE imputation (m=5) on training set ...")

# Identify baseline columns to impute (exclude ID, outcomes, daily SOFA)
skip_cols <- c("ID", "IPN", "IPN_strict", "event_type",
               "time_to_IPN", "event_time",
               "Antibiotic_actual", "treatment",
               sofa_cols_all)

impute_cols <- setdiff(colnames(train), skip_cols)
impute_cols <- impute_cols[sapply(train[impute_cols, drop = FALSE], is.numeric)]

if (length(impute_cols) > 0) {
  imp_data <- train[, impute_cols, drop = FALSE]

  # Only run MICE if there are actually missing values
  has_missing <- any(colMeans(is.na(imp_data)) > 0)
  if (has_missing) {
    imp_obj <- tryCatch({
      mice::mice(imp_data, m = 5, maxit = 5, seed = 2025, printFlag = FALSE)
    }, error = function(e) {
      warning("MICE failed: ", e$message, ". Using median imputation instead.")
      NULL
    })

    if (!is.null(imp_obj)) {
      train_imp <- mice::complete(imp_obj, 1)
      train[, impute_cols] <- train_imp[, impute_cols]
      message("MICE completed successfully (m=5, maxit=5).")
    } else {
      for (col in impute_cols) {
        train[[col]][is.na(train[[col]])] <- median(train[[col]], na.rm = TRUE)
      }
      message("Median imputation used as fallback.")
    }
  } else {
    message("No missing values in baseline variables — skipping MICE.")
  }
}

# MICE imputation for prospective cohort
message("\nRunning MICE imputation (m=5) on prospective cohort ...")
prosp_impute_cols <- intersect(impute_cols, colnames(prosp))
prosp_impute_cols <- prosp_impute_cols[sapply(prosp[prosp_impute_cols, drop = FALSE], is.numeric)]
if (length(prosp_impute_cols) > 0) {
  prosp_imp_data <- prosp[, prosp_impute_cols, drop = FALSE]
  if (any(colMeans(is.na(prosp_imp_data)) > 0)) {
    prosp_imp <- tryCatch({
      mice::mice(prosp_imp_data, m = 5, maxit = 5, seed = 2025, printFlag = FALSE)
    }, error = function(e) {
      warning("MICE failed for prospective cohort: ", e$message, ". Using median imputation.")
      NULL
    })
    if (!is.null(prosp_imp)) {
      prosp_imp_complete <- mice::complete(prosp_imp, 1)
      prosp[, prosp_impute_cols] <- prosp_imp_complete[, prosp_impute_cols]
      message("Prospective cohort MICE completed.")
    } else {
      for (col in prosp_impute_cols) {
        prosp[[col]][is.na(prosp[[col]])] <- median(prosp[[col]], na.rm = TRUE)
      }
    }
  }
}

# MICE imputation for TRACE cohort
message("\nRunning MICE imputation (m=5) on TRACE cohort ...")
trace_impute_cols <- intersect(impute_cols, colnames(trace))
trace_impute_cols <- trace_impute_cols[sapply(trace[trace_impute_cols, drop = FALSE], is.numeric)]
if (length(trace_impute_cols) > 0) {
  trace_imp_data <- trace[, trace_impute_cols, drop = FALSE]
  if (any(colMeans(is.na(trace_imp_data)) > 0)) {
    trace_imp <- tryCatch({
      mice::mice(trace_imp_data, m = 5, maxit = 5, seed = 2025, printFlag = FALSE)
    }, error = function(e) {
      warning("MICE failed for TRACE cohort: ", e$message, ". Using median imputation.")
      NULL
    })
    if (!is.null(trace_imp)) {
      trace_imp_complete <- mice::complete(trace_imp, 1)
      trace[, trace_impute_cols] <- trace_imp_complete[, trace_impute_cols]
      message("TRACE cohort MICE completed.")
    } else {
      for (col in trace_impute_cols) {
        trace[[col]][is.na(trace[[col]])] <- median(trace[[col]], na.rm = TRUE)
      }
    }
  }
}
message("Imputation complete for all datasets.")

# ── Daily SOFA: keep missing as NA (handled by joint model) ───────────
message("Daily SOFA missing values retained as NA (joint model handles missing longitudinal data).")

# ── Create long format for SOFA scores ─────────────────────────────────
create_sofa_long <- function(df) {
  id_col <- "ID"
  if (!(id_col %in% colnames(df))) {
    stop("ID column missing — cannot create long format data.")
  }

  # Time-constant covariates (everything except daily SOFA)
  const_vars <- setdiff(colnames(df), sofa_cols_all)

  # Standardize SOFA column names to sequential order
  present_sofa <- intersect(sofa_cols_all, colnames(df))
  if (length(present_sofa) < 2) {
    stop("Insufficient daily SOFA columns to create long format.")
  }

  # Reshape to long
  long <- reshape(
    df,
    varying   = list(present_sofa),
    v.names   = "SOFA",
    timevar   = "day",
    times     = as.numeric(gsub("SOFA_day", "", present_sofa)),
    idvar     = id_col,
    direction = "long"
  )

  # Merge time-constant covariates
  const_df <- df[, const_vars, drop = FALSE]
  long <- merge(long, const_df, by = id_col, all.x = TRUE)
  long <- long[order(long[[id_col]], long$day), ]
  rownames(long) <- NULL

  message(sprintf("  Long format: %d observations, %d subjects, days %d-%d",
                  nrow(long), length(unique(long[[id_col]])),
                  min(long$day, na.rm = TRUE), max(long$day, na.rm = TRUE)))
  return(long)
}

message("\nCreating long-format SOFA data ...")
train_long <- create_sofa_long(train)
prosp_long <- create_sofa_long(prosp)
trace_long <- create_sofa_long(trace)

# ── Save processed data ───────────────────────────────────────────────
saveRDS(train,      "data/train_processed.rds")
saveRDS(prosp,      "data/prosp_processed.rds")
saveRDS(trace,      "data/trace_processed.rds")
saveRDS(train_long, "data/train_long.rds")
saveRDS(prosp_long, "data/prosp_long.rds")
saveRDS(trace_long, "data/trace_long.rds")
cat("Processed data saved.\n")
cat("=== Preprocessing complete ===\n")
