# 01_load_data.R — Read CSV files; STOP with error if missing
cat("\n=== Step 1: Data Loading ===\n")

# ── Required column names (exact specification per paper protocol) ──────
# Baseline variables (selected predictors from feature selection)
required_baseline <- c(
  "ID", "Age", "Sex",
  "Cr", "HR", "RR", "CTSI", "ALB"
)

# Daily SOFA scores (days 1 through 7 minimum; extended range supported)
required_daily_sofa <- paste0("SOFA_day", 1:7)

# Outcome variables
required_outcome <- c(
  "IPN",              # composite IPN outcome (0/1)
  "time_to_IPN"       # time to IPN or censoring (days, max 90)
)

# Treatment variable for TRACE cohort
required_treatment <- c("treatment")  # thymosin_alpha1 (0/1)

# Additional useful variables (warn if missing but do not stop)
optional_vars <- c(
  "IPN_strict",       # strict IPN definition (microbiological + imaging confirmed)
  "Antibiotic_actual", # actual antibiotic use (0/1) for simulation
  "SIRS_day1", "SIRS_day2", "SIRS_day3", "SIRS_day4",
  "SIRS_day5", "SIRS_day6", "SIRS_day7",
  "APACHEII_day1", "APACHEII_day2", "APACHEII_day3", "APACHEII_day4",
  "APACHEII_day5", "APACHEII_day6", "APACHEII_day7",
  "BISAP_day1", "BISAP_day2", "BISAP_day3", "BISAP_day4",
  "BISAP_day5", "BISAP_day6", "BISAP_day7",
  "CRP", "PCT", "WBC", "HCT", "PLT", "TBIL", "BUN",
  "Na", "K", "Ca", "Glu", "SBP", "DBP", "MAP", "Temp",
  "SpO2", "FiO2", "pH", "PaO2", "PaCO2", "Lac", "HCO3", "BE",
  "Comorbidity_DM", "Comorbidity_HTN", "Comorbidity_CVD",
  "Comorbidity_COPD", "Comorbidity_CKD",
  "Etiology_Biliary", "Etiology_Alcohol",
  "Etiology_Hypertriglyceridemia", "Etiology_Other",
  "LOS_ICU", "LOS_hospital", "Ventilation", "CRRT", "Center"
)

all_expected <- unique(c(required_baseline, required_daily_sofa,
                         required_outcome, required_treatment, optional_vars))

# ── Column name aliases (maps alternate names to canonical names) ───────
# Supports both SOFA_day1 and SOFA_0/1 style naming
col_aliases <- list(
  "SOFA_day0" = "SOFA_day0",
  "time_to_IPN" = "time_to_IPN",
  "Time_to_IPN" = "time_to_IPN",
  "IPN_strict" = "IPN_strict",
  "Status" = "IPN",
  "event_time" = "time_to_IPN",
  "Antibiotic" = "Antibiotic_actual"
)
# Also support SOFA_0, SOFA_1, ... as aliases for SOFA_day0, SOFA_day1, ...
for (d in 0:14) {
  col_aliases[[paste0("SOFA_", d)]] <- paste0("SOFA_day", d)
}

# ── Helper: load and validate a single cohort ───────────────────────────
load_cohort <- function(filepath, cohort_name, extra_required = NULL) {
  if (!file.exists(filepath)) {
    sample_size <- if (grepl("training", filepath)) "n=812"
                   else if (grepl("prospective", filepath)) "n=465"
                   else "n=437"
    extra_str <- if (!is.null(extra_required))
      paste0("Additional: ", paste(extra_required, collapse = ", ")) else ""
    msg <- paste0(
      "\n===================================================================\n",
      "  ERROR: ", filepath, " not found.\n",
      "  Please place your CSV file at: ", filepath, "\n",
      "  Expected sample size: ", sample_size, "\n",
      "===================================================================\n",
      "  Required columns (exact spelling):\n",
      "    Baseline: ", paste(required_baseline, collapse = ", "), "\n",
      "    Daily SOFA: SOFA_day1, SOFA_day2, SOFA_day3, ..., SOFA_day14",
      " (at minimum SOFA_day1 through SOFA_day7)\n",
      "    Outcome: ", paste(required_outcome, collapse = ", "), "\n",
      "    ", extra_str, "\n",
      "===================================================================\n"
    )
    stop(msg)
  }

  df <- tryCatch(
    read.csv(filepath, stringsAsFactors = FALSE),
    error = function(e) stop("Error reading ", filepath, ": ", e$message)
  )

  message(sprintf("  Loaded %s: %d rows x %d columns", cohort_name, nrow(df), ncol(df)))

  # Apply column name aliases
  for (old_name in names(col_aliases)) {
    if (old_name %in% colnames(df) && !(col_aliases[[old_name]] %in% colnames(df))) {
      colnames(df)[colnames(df) == old_name] <- col_aliases[[old_name]]
    }
  }

  # Check required columns
  required_all <- c(required_baseline, required_outcome)
  if (grepl("trace", filepath)) {
    required_all <- c(required_all, "treatment")
  }
  if (!is.null(extra_required)) {
    required_all <- c(required_all, extra_required)
  }

  missing_req <- setdiff(required_all, colnames(df))
  if (length(missing_req) > 0) {
    stop(sprintf(
      "  FATAL: Missing required columns in %s:\n    %s\n",
      filepath, paste(missing_req, collapse = ", ")
    ))
  }

  # Check for at least some daily SOFA columns
  sofa_cols_present <- grep("^SOFA_day\\d+$", colnames(df), value = TRUE)
  if (length(sofa_cols_present) < 3) {
    stop(sprintf(
      "  FATAL: %s has only %d daily SOFA columns. Need at least SOFA_day1 through SOFA_day7.\n",
      filepath, length(sofa_cols_present)
    ))
  }
  message(sprintf("  Daily SOFA columns found: %d (days %s to %s)",
                  length(sofa_cols_present),
                  min(as.numeric(gsub("SOFA_day", "", sofa_cols_present))),
                  max(as.numeric(gsub("SOFA_day", "", sofa_cols_present)))))

  # Warn about missing optional columns
  missing_opt <- setdiff(optional_vars, colnames(df))
  if (length(missing_opt) > 0 && length(missing_opt) < length(optional_vars)) {
    warning("  Optional columns missing in ", filepath, ":\n    ",
            paste(head(missing_opt, 10), collapse = ", "),
            if (length(missing_opt) > 10) sprintf(" ... and %d more", length(missing_opt) - 10))
  }

  return(df)
}

# ── Load three cohorts ──────────────────────────────────────────────────
message("\n[1/3] Loading training cohort (retrospective, expected n=812) ...")
train <- load_cohort("data/training_cohort.csv", "Training")

message("\n[2/3] Loading prospective cohort (expected n=465) ...")
prosp <- load_cohort("data/prospective_cohort.csv", "Prospective")

message("\n[3/3] Loading TRACE cohort (expected n=437) ...")
trace <- load_cohort("data/trace_cohort.csv", "TRACE",
                     extra_required = "treatment")

# ── Save RDS ──────────────────────────────────────────────────────────
saveRDS(train, "data/training_cohort.rds")
saveRDS(prosp, "data/prospective_cohort.rds")
saveRDS(trace, "data/trace_cohort.rds")
cat("\nRaw data saved to data/*.rds\n")
cat("=== Data loading complete ===\n")
