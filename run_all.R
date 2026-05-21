# run_all.R — Master script: runs the complete SADRA model analysis pipeline
#
# Usage: Rscript run_all.R
#        or source("run_all.R") in R/RStudio
#
# Environment variables:
#   SADRA_TESTING=TRUE   → Use reduced MCMC iterations for faster testing
#   GITHUB_PAT=<token>   → Personal access token for GitHub push
#
# Expected input files (in data/):
#   training_cohort.csv
#   prospective_cohort.csv
#   trace_cohort.csv

cat("
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║          SADRA Model — Full Analysis Pipeline                 ║
║  Sequential Assessment of Dynamic Risk in Acute Pancreatitis  ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
")

start_time <- Sys.time()

# ── Helper: source with error handling ────────────────────────────────
source_safe <- function(script, description) {
  cat(sprintf("\n%s\n  Executing: %s ...\n", paste(rep("─", 60), collapse = ""), script))
  tryCatch({
    source(script)
    cat(sprintf("  ✓ %s completed successfully.\n", description))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("  ✗ %s FAILED: %s\n", description, e$message))
    return(FALSE)
  })
}

# ── Check for data files (strict: STOP if missing) ────────────────────
cat("\nChecking for input data files ...\n")
missing_data <- c()
for (f in c("data/training_cohort.csv", "data/prospective_cohort.csv", "data/trace_cohort.csv")) {
  if (file.exists(f)) {
    cat(sprintf("  ✓ %s found\n", f))
  } else {
    cat(sprintf("  ✗ %s MISSING\n", f))
    missing_data <- c(missing_data, f)
  }
}
if (length(missing_data) > 0) {
  msg <- paste0(
    "\n===================================================================\n",
    "  FATAL: Required data files not found:\n",
    "    ", paste(missing_data, collapse = "\n    "), "\n",
    "  Please place your CSV files in the data/ folder:\n",
    "    data/training_cohort.csv      (retrospective, n=812)\n",
    "    data/prospective_cohort.csv   (prospective, n=465)\n",
    "    data/trace_cohort.csv         (TRACE trial, n=437)\n",
    "\n  Required columns (exact spelling):\n",
    "    ID, Age, Sex, Cr, HR, RR, CTSI, ALB\n",
    "    SOFA_day1, SOFA_day2, ..., SOFA_day7 (at minimum)\n",
    "    IPN, time_to_IPN\n",
    "    TRACE cohort additionally: treatment\n",
    "===================================================================\n"
  )
  stop(msg)
}

# ── Run pipeline ──────────────────────────────────────────────────────
scripts <- list(
  c("R/00_setup.R",               "Package setup"),
  c("R/01_load_data.R",           "Data loading"),
  c("R/02_preprocess.R",          "Preprocessing & imputation"),
  c("R/03_feature_selection.R",   "Feature selection (LASSO + Boruta)"),
  c("R/04_joint_model.R",         "Bayesian joint model"),
  c("R/05_performance.R",         "Performance evaluation"),
  c("R/06_sensitivity.R",         "Sensitivity analyses"),
  c("R/07_antibiotic_simulation.R","Antibiotic simulation"),
  c("R/08_treatment_heterogeneity.R", "Treatment heterogeneity"),
  c("R/09_github_upload.R",       "GitHub upload")
)

pipeline_ok <- TRUE
for (s in scripts) {
  ok <- source_safe(s[1], s[2])
  if (!ok) {
    cat(sprintf("\n⚠ Pipeline step '%s' failed. Continuing with remaining steps ...\n", s[2]))
    pipeline_ok <- FALSE
  }
}

# ── Summary ────────────────────────────────────────────────────────────
elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\n%s\n", paste(rep("═", 60), collapse = "")))

if (pipeline_ok) {
  cat("\n✓ Analysis completed successfully!\n")
  cat(sprintf("  Total time: %.1f minutes\n", as.numeric(elapsed)))
  cat("\nOutput files:\n")
  cat("  - output/selected_features.rds\n")
  cat("  - output/joint_model.rds\n")
  cat("  - output/AUC_results.csv\n")
  cat("  - output/risk_cutoff.rds\n")
  cat("  - output/antibiotic_simulation.csv\n")
  cat("  - output/TRACE_validation_AUC.csv\n")
  cat("  - output/Table_S2_longitudinal_params.csv/.txt\n")
  cat("  - output/Table_S3_survival_params.txt\n")
  cat("  - output/Table_S4_association_params.txt\n")
  cat("  - output/Table_S5_ML_comparison.csv\n")
  cat("  - output/Table_S6_CATE_by_risk.csv\n")
  cat("\nFigures:\n")
  cat("  - figures/Figure_2_AUC_curve.pdf\n")
  cat("  - figures/Figure_2_calibration.pdf\n")
  cat("  - figures/Figure_3_KM_curves.pdf\n")
  cat("  - figures/Figure_4_antibiotic_simulation.pdf\n")
  cat("  - figures/Figure_5_CATE_forest.pdf\n")
  cat("  - figures/Figure_S1_feature_selection.pdf\n")
  cat("  - figures/Figure_S3_ML_comparison.pdf\n")
  cat("  - figures/Figure_S4_SIRS.pdf\n")
  cat("  - figures/Figure_S5_APACHEII.pdf\n")
  cat("  - figures/Figure_S6_BISAP.pdf\n")
  cat("  - figures/Figure_S7_strict_IPN.pdf\n")
  cat("  - figures/Figure_S8_TRACE_validation.pdf\n")
  cat("\nShiny App:\n")
  cat("  Run: shiny::runApp('shiny_app')\n")
  cat("\nReport:\n")
  cat("  Knit: reports/report.Rmd\n")
  cat("\nGitHub:\n")
  cat("  Repository: https://github.com/xiongsx123-art/SADRA-model.git\n")
} else {
  cat("\n⚠ Analysis completed with some errors. Check messages above.\n")
  cat(sprintf("  Total time: %.1f minutes\n", as.numeric(elapsed)))
}

cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
cat("║         Analysis pipeline finished                            ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")
