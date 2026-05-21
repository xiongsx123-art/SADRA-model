# 04_joint_model.R — Bayesian joint model (JMbayes2) for longitudinal SOFA + survival
cat("\n=== Step 4: Bayesian Joint Model (SADRA) ===\n")

train_long <- readRDS("data/train_long.rds")
selected_features <- readRDS("output/selected_features.rds")

# ── Configuration ─────────────────────────────────────────────────────
testing_mode <- Sys.getenv("SADRA_TESTING", unset = "FALSE") == "TRUE"
n_iter  <- if (testing_mode) 5000  else 30000
n_burn  <- if (testing_mode) 2500  else 15000
n_thin  <- if (testing_mode) 2     else 4
n_chains <- 4

if (testing_mode) {
  message("*** TESTING MODE: Reduced MCMC iterations (", n_iter, ") ***")
}

# ── Prepare data ──────────────────────────────────────────────────────
stopifnot(all(c("ID", "day", "SOFA", "event_time", "event_type") %in% colnames(train_long)))

# Survival data: one row per subject
surv_data <- train_long[!duplicated(train_long$ID), ]
surv_data$event_time <- pmax(surv_data$event_time, 0.1)

# Baseline covariates (from feature selection, intersect with data)
baseline_covars <- intersect(selected_features, colnames(surv_data))
if (length(baseline_covars) < 3) {
  baseline_covars <- intersect(c("Cr", "HR", "RR", "CTSI", "ALB"), colnames(surv_data))
  if (length(baseline_covars) == 0) baseline_covars <- c("Age", "Sex")
}

message(sprintf("Baseline covariates: %s", paste(baseline_covars, collapse = ", ")))
message(sprintf("Subjects: %d, Longitudinal observations: %d",
                length(unique(train_long$ID)), nrow(train_long)))

# ── Build formulas ────────────────────────────────────────────────────
# Longitudinal submodel: SOFA ~ ns(day, df=2) + baseline_covars
# with random intercept and random slope for day
lme_formula_str <- paste0(
  "SOFA ~ ns(day, df = 2) + ",
  paste(baseline_covars, collapse = " + ")
)
lme_formula <- as.formula(lme_formula_str)
message("Longitudinal formula: ", lme_formula_str)

# Survival submodel: Surv(event_time, event_type) ~ baseline_covars
cox_formula_str <- paste0(
  "Surv(event_time, event_type) ~ ",
  paste(baseline_covars, collapse = " + ")
)
cox_formula <- as.formula(cox_formula_str)
message("Survival formula: ", cox_formula_str)

# Random effects: random intercept + random slope for day
random_formula <- ~ day | ID

# ── Fit sub-models for JMbayes2 ───────────────────────────────────────
message("\nFitting longitudinal mixed model (nlme::lme) ...")
lme_fit <- tryCatch({
  nlme::lme(
    fixed = lme_formula,
    data = train_long,
    random = random_formula,
    control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200),
    na.action = na.omit
  )
}, error = function(e) {
  warning("nlme::lme failed: ", e$message, ". Trying simpler random effects.")
  tryCatch({
    nlme::lme(
      fixed = lme_formula,
      data = train_long,
      random = ~ 1 | ID,
      control = nlme::lmeControl(opt = "optim", maxIter = 200),
      na.action = na.omit
    )
  }, error = function(e2) {
    warning("Simplified lme also failed: ", e2$message)
    NULL
  })
})

message("\nFitting Cox survival model ...")
cox_fit <- tryCatch({
  survival::coxph(cox_formula, data = surv_data, model = TRUE, x = TRUE)
}, error = function(e) {
  warning("coxph failed: ", e$message)
  NULL
})

# ── Fit Bayesian joint model ──────────────────────────────────────────
message("\nFitting Bayesian joint model (", n_iter, " iter, ", n_chains, " chains) ...")
message("This may take 30-120 minutes for full MCMC ...")

jm_fit <- NULL

if (!is.null(lme_fit) && !is.null(cox_fit)) {
  jm_fit <- tryCatch({

    # Define priors explicitly
    # Fixed effects: N(0, 10^2) for each coefficient
    n_fixed <- length(fixef(lme_fit))
    prior_fixed <- lapply(seq_len(n_fixed), function(i) {
      list(mean = 0, Tau = 1 / (10^2))  # precision parameterization
    })

    # Random effects covariance: inverse Wishart with df=3, scale = diag(0.1, 0.01)
    n_random <- ncol(as.matrix(ranef(lme_fit)[[1]]))
    prior_random <- list(
      df = n_random + 1,
      scale = diag(c(0.1, 0.01)[seq_len(n_random)])
    )

    # Residual variance: inverse Gamma(0.01, 0.01)
    prior_residual <- list(shape = 0.01, rate = 0.01)

    # Construct prior list for JMbayes2
    jm_priors <- list(
      "fixed"  = prior_fixed,
      "random" = prior_random,
      "sigma"  = prior_residual
    )

    jm <- JMbayes2::jm(
      Surv_object  = cox_fit,
      Mixed_objects = lme_fit,
      time_var     = "day",
      functional_forms = ~ value(SOFA),
      n_iter     = n_iter,
      n_burnin   = n_burn,
      n_thin     = n_thin,
      n_chains   = n_chains,
      cores      = min(n_chains, max(1, parallel::detectCores() - 1)),
      priors     = jm_priors,
      seed       = 2025
    )
    message("Joint model fitted successfully.")
    jm

  }, error = function(e) {
    warning("Joint model fitting failed: ", e$message,
            "\nTrying with default priors ...")
    tryCatch({
      jm2 <- JMbayes2::jm(
        Surv_object  = cox_fit,
        Mixed_objects = lme_fit,
        time_var     = "day",
        functional_forms = ~ value(SOFA),
        n_iter     = min(n_iter, 10000),
        n_burnin   = min(n_burn, 5000),
        n_thin     = 2,
        n_chains   = 2,
        cores      = 2,
        seed       = 2025
      )
      message("Joint model fitted with default priors.")
      jm2
    }, error = function(e2) {
      warning("Joint model with default priors also failed: ", e2$message)
      NULL
    })
  })
} else {
  warning("Cannot fit joint model: sub-models unavailable.")
}

# ── Convergence diagnostics ───────────────────────────────────────────
if (!is.null(jm_fit)) {
  message("\nChecking convergence (R-hat) ...")
  rhats <- tryCatch({
    summary(jm_fit)$Outcome1$Rhat
  }, error = function(e) {
    warning("Could not extract R-hat: ", e$message)
    NULL
  })

  if (!is.null(rhats)) {
    rhat_ok <- all(rhats <= 1.05, na.rm = TRUE)
    message(sprintf("R-hat <= 1.05 for all parameters: %s",
                    ifelse(rhat_ok, "YES (converged)", "NO — some parameters may not have converged")))
    if (!rhat_ok) {
      high_rhat <- rhats[rhats > 1.05 & !is.na(rhats)]
      message(sprintf("  Parameters with R-hat > 1.05: %d", length(high_rhat)))
    }
  }
}

# ── Save model ────────────────────────────────────────────────────────
saveRDS(jm_fit, "output/joint_model.rds")
message("Joint model saved to output/joint_model.rds")

# ── Parameter tables (Table S2–S4) ────────────────────────────────────
if (!is.null(jm_fit)) {
  message("\nExtracting model parameters ...")

  jm_summary <- tryCatch({
    summary(jm_fit)
  }, error = function(e) {
    warning("Could not summarize model: ", e$message)
    NULL
  })

  if (!is.null(jm_summary)) {
    # Table S2: Longitudinal submodel parameters
    sink("output/Table_S2_longitudinal_params.txt")
    cat("Table S2. Longitudinal Submodel Parameters (SOFA trajectory)\n")
    cat("==========================================================\n\n")
    cat("Fixed effects:\n")
    print(jm_summary$Outcome1$statistics$Mean[grep("Y_", names(jm_summary$Outcome1$statistics$Mean))])
    cat("\nRandom effects covariance:\n")
    print(jm_summary$Outcome1$D)
    cat("\nResidual variance:\n")
    print(jm_summary$Outcome1$sigma)
    sink()
    message("Table S2 saved.")

    # Table S3: Survival submodel parameters
    sink("output/Table_S3_survival_params.txt")
    cat("Table S3. Survival Submodel Parameters (time-to-IPN)\n")
    cat("=====================================================\n\n")
    if (!is.null(jm_summary$Outcome2)) {
      print(jm_summary$Outcome2$statistics)
    } else {
      cat("Survival parameters embedded in full summary.\n")
      print(jm_summary)
    }
    sink()
    message("Table S3 saved.")

    # Table S4: Association parameter
    sink("output/Table_S4_association_params.txt")
    cat("Table S4. Association Parameter (Current-Value SOFA)\n")
    cat("====================================================\n\n")
    cat("Association structure: current value of SOFA\n")
    cat("Parameter interpretation: log hazard ratio per unit increase in expected SOFA\n\n")
    if (!is.null(jm_summary$Outcome2)) {
      assoc_idx <- grep("value", rownames(jm_summary$Outcome2$statistics), ignore.case = TRUE)
      if (length(assoc_idx) > 0) {
        print(jm_summary$Outcome2$statistics[assoc_idx, , drop = FALSE])
      }
    }
    sink()
    message("Table S4 saved.")

    # Save full summary as CSV for easy inspection
    write.csv(as.data.frame(jm_summary$Outcome1$statistics),
              "output/Table_S2_longitudinal_params.csv")
  }
}

cat("=== Joint model complete ===\n")
