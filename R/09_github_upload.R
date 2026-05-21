# 09_github_upload.R — Git operations and GitHub push
cat("\n=== Step 9: GitHub Upload ===\n")

repo_dir <- getwd()
repo_url <- "https://github.com/xiongsx123-art/SADRA-model.git"

# ── Check for required files ──────────────────────────────────────────
required_files <- c(
  "R/00_setup.R", "R/01_load_data.R", "R/02_preprocess.R",
  "R/03_feature_selection.R", "R/04_joint_model.R", "R/05_performance.R",
  "R/06_sensitivity.R", "R/07_antibiotic_simulation.R",
  "R/08_treatment_heterogeneity.R", "R/09_github_upload.R",
  "run_all.R", "reports/report.Rmd", "shiny_app/app.R",
  ".gitignore"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  warning("Some project files are missing:\n  ",
          paste(missing_files, collapse = "\n  "))
}

# ── Initialize git if needed ──────────────────────────────────────────
if (!dir.exists(".git")) {
  message("Initializing git repository ...")
  system2("git", c("init"))
}

# ── Configure remote ──────────────────────────────────────────────────
existing_remote <- tryCatch({
  system2("git", c("remote", "get-url", "origin"), stdout = TRUE, stderr = FALSE)
}, error = function(e) character(0))

if (length(existing_remote) == 0) {
  message("Adding remote origin ...")
  system2("git", c("remote", "add", "origin", repo_url))
} else {
  message(sprintf("Remote origin already configured: %s", existing_remote[1]))
  if (existing_remote[1] != repo_url) {
    message(sprintf("Updating remote origin to: %s", repo_url))
    system2("git", c("remote", "set-url", "origin", repo_url))
  }
}

# ── Stage files ──────────────────────────────────────────────────────
message("\nStaging files ...")
# Stage new and modified files (not data CSVs)
stage_files <- c(
  "R/00_setup.R", "R/01_load_data.R", "R/02_preprocess.R",
  "R/03_feature_selection.R", "R/04_joint_model.R", "R/05_performance.R",
  "R/06_sensitivity.R", "R/07_antibiotic_simulation.R",
  "R/08_treatment_heterogeneity.R", "R/09_github_upload.R",
  "run_all.R", "reports/report.Rmd", "shiny_app/app.R",
  ".gitignore", "README.md",
  "data/.gitkeep", "output/.gitkeep", "figures/.gitkeep"
)

for (f in stage_files) {
  if (file.exists(f)) {
    system2("git", c("add", f))
  }
}

# ── Commit ────────────────────────────────────────────────────────────
commit_msg <- "Complete SADRA model analysis pipeline

- 00_setup.R: Package loading and helper functions
- 01_load_data.R: CSV validation with required column checks
- 02_preprocess.R: Missing data handling (MICE) and long-format SOFA
- 03_feature_selection.R: LASSO Cox + Boruta with RSF importance
- 04_joint_model.R: Bayesian joint model (JMbayes2) with explicit priors
- 05_performance.R: Time-dep AUC with cluster bootstrap, calibration, KM
- 06_sensitivity.R: 6 sensitivity analyses (static ML, alt assoc, alt scores, etc.)
- 07_antibiotic_simulation.R: Antibiotic strategy simulation (Sankey + UpSet)
- 08_treatment_heterogeneity.R: Double ML ATE/CATE in TRACE cohort
- 09_github_upload.R: Git operations and push
- run_all.R: Master pipeline script
- reports/report.Rmd: R Markdown analysis report
- shiny_app/app.R: IPN risk calculator Shiny app"

system2("git", c("commit", "-m", commit_msg))

# ── Push to GitHub ────────────────────────────────────────────────────
message("\nPushing to GitHub ...")

# Check for PAT in environment
github_pat <- Sys.getenv("GITHUB_PAT")
if (github_pat == "") {
  github_pat <- Sys.getenv("GITHUB_TOKEN")
}

if (github_pat != "") {
  # Set up credential helper with PAT
  message("GITHUB_PAT found. Configuring authentication ...")
  # Construct URL with token for push
  push_url <- sub("https://", paste0("https://", github_pat, "@"), repo_url)

  push_result <- system2("git", c("push", "-u", push_url, "main"),
                         stdout = TRUE, stderr = TRUE)
  message(paste(push_result, collapse = "\n"))

  # Reset remote URL to clean version (without token)
  system2("git", c("remote", "set-url", "origin", repo_url))

} else {
  # Try SSH or existing credential helper
  message("GITHUB_PAT not set. Attempting push with existing credentials ...")
  push_result <- system2("git", c("push", "-u", "origin", "main"),
                         stdout = TRUE, stderr = TRUE)

  push_status <- attr(push_result, "status")
  if (!is.null(push_status) && push_status != 0) {
    message("\n===========================================================")
    message("  AUTOMATIC PUSH FAILED")
    message("===========================================================")
    message("  The code could not push to GitHub automatically.")
    message("")
    message("  Manual push instructions:")
    message("  1. Ensure you have a GitHub Personal Access Token (PAT)")
    message("     with 'repo' scope from: https://github.com/settings/tokens")
    message("  2. Set the environment variable:")
    message("     Sys.setenv(GITHUB_PAT = \"your_token_here\")")
    message("  3. Run: git push -u origin main")
    message("")
    message("  Repository URL: ", repo_url)
    message("===========================================================\n")
  } else {
    message("Push successful!")
  }
}

cat("=== GitHub upload complete ===\n")
