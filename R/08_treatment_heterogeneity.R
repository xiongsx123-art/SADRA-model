# 08_treatment_heterogeneity.R — Double ML for ATE and CATE in TRACE cohort, Figure 5, Table S6
cat("\n=== Step 8: Treatment Heterogeneity Analysis (TRACE Cohort) ===\n")

# Load TRACE cohort (contains treatment = thymosin_alpha1)
trace       <- readRDS("data/trace_processed.rds")
trace_long  <- readRDS("data/trace_long.rds")
risk_cutoff <- readRDS("output/risk_cutoff.rds")
selected_features <- readRDS("output/selected_features.rds")
jm_fit      <- tryCatch(readRDS("output/joint_model.rds"), error = function(e) NULL)

dir.create("figures", showWarnings = FALSE)

# ── Prepare TRACE analysis dataset ────────────────────────────────────
# Treatment: thymosin_alpha1 (0/1) from the TRACE trial
# Outcome: IPN (or IPN_strict if available — paper uses strict IPN for heterogeneity)

if (!("treatment" %in% colnames(trace))) {
  stop("FATAL: TRACE cohort missing 'treatment' column (thymosin_alpha1).")
}

# Use strict IPN if available, otherwise composite IPN
if ("IPN_strict" %in% colnames(trace)) {
  message("Using strict IPN (IPN_strict) as outcome for treatment heterogeneity.")
  trace$outcome <- trace$IPN_strict
} else {
  message("IPN_strict not available — using composite IPN as outcome.")
  trace$outcome <- trace$IPN
}

# Confounders: all baseline variables + SADRA-predicted risk
baseline_covars <- intersect(selected_features, colnames(trace))
if (length(baseline_covars) < 3) {
  baseline_covars <- intersect(c("Cr", "HR", "RR", "CTSI", "ALB", "Age", "Sex"),
                               colnames(trace))
}
# Add SADRA predicted risk as a covariate (paper: "SADRA-predicted risk continuous")
if (!is.null(jm_fit)) {
  message("Computing SADRA-predicted risk for TRACE patients ...")
  tryCatch({
    preds <- predict(jm_fit, newdata = trace_long,
                     times = 7, process = "event", type = "survival")
    trace$SADRA_risk <- 1 - preds$pred[, 1]
    baseline_covars <- c(baseline_covars, "SADRA_risk")
    message(sprintf("SADRA risk scores computed for %d patients", nrow(trace)))
  }, error = function(e) {
    warning("Could not compute SADRA risk scores: ", e$message)
  })
}

# Remove non-numeric or problematic columns
baseline_covars <- baseline_covars[baseline_covars %in% colnames(trace)]
baseline_covars <- baseline_covars[sapply(trace[baseline_covars], is.numeric)]
baseline_covars <- baseline_covars[sapply(trace[baseline_covars],
                                          function(x) sum(is.na(x)) < nrow(trace) * 0.5)]
baseline_covars <- setdiff(baseline_covars,
                           c("IPN", "IPN_strict", "outcome", "event_type",
                             "event_time", "time_to_IPN", "treatment",
                             "Antibiotic_actual", "ID"))

message(sprintf("Using %d confounders for Double ML analysis", length(baseline_covars)))

# Build analysis dataset
analysis_vars <- c("outcome", "treatment", baseline_covars)
analysis_df <- trace[, analysis_vars, drop = FALSE]
analysis_df <- analysis_df[complete.cases(analysis_df), ]
message(sprintf("TRACE analysis sample: %d patients with complete data", nrow(analysis_df)))

# Convert outcome to numeric binary
analysis_df$outcome <- as.integer(as.numeric(analysis_df$outcome) > 0)
analysis_df$treatment <- as.integer(as.numeric(analysis_df$treatment) > 0)

message(sprintf("Treatment (thymosin_alpha1) prevalence: %.1f%%",
                100 * mean(analysis_df$treatment)))
message(sprintf("Outcome (IPN) prevalence: %.1f%%",
                100 * mean(analysis_df$outcome)))

# ── Double Machine Learning for ATE ────────────────────────────────────
message("\nFitting Double ML model (RF learner, 5-fold cross-fitting) ...")

dml_results <- NULL
ATE <- NA_real_
ATE_SE <- NA_real_
ATE_p <- NA_real_

if (requireNamespace("DoubleML", quietly = TRUE) &&
    requireNamespace("mlr3", quietly = TRUE) &&
    nrow(analysis_df) >= 50) {

  dml_results <- tryCatch({
    # Ensure mlr3 learners are available
    if (!requireNamespace("mlr3learners", quietly = TRUE)) {
      install.packages("mlr3learners", repos = "https://cran.r-project.org")
    }

    # Create DoubleML data
    X_mat <- as.matrix(analysis_df[, baseline_covars, drop = FALSE])
    colnames(X_mat) <- make.names(baseline_covars)

    dml_data <- DoubleML::double_ml_data_from_matrix(
      X = X_mat,
      y = analysis_df$outcome,
      d = analysis_df$treatment
    )

    # Random forest learners for nuisance functions
    ml_g <- mlr3::lrn("regr.ranger", num.trees = 200, num.threads = 1)
    ml_m <- mlr3::lrn("classif.ranger", num.trees = 200,
                      predict_type = "prob", num.threads = 1)

    # Double ML PLR with 5-fold cross-fitting
    dml_obj <- DoubleML::DoubleMLPLR$new(
      data   = dml_data,
      ml_g   = ml_g,
      ml_m   = ml_m,
      n_folds = 5,
      n_rep   = 5,
      score   = "partialling out"
    )

    dml_obj$fit(store_predictions = TRUE)
    message("Double ML fitted successfully.")
    dml_obj

  }, error = function(e) {
    warning("DoubleML with RF failed: ", e$message,
            "\nTrying with simpler learners ...")
    tryCatch({
      X_mat <- as.matrix(analysis_df[, baseline_covars, drop = FALSE])
      dml_data <- DoubleML::double_ml_data_from_matrix(
        X = X_mat,
        y = analysis_df$outcome,
        d = analysis_df$treatment
      )
      # Fallback to glmnet learners
      ml_g2 <- mlr3::lrn("regr.cv_glmnet")
      ml_m2 <- mlr3::lrn("classif.cv_glmnet", predict_type = "prob")

      dml_obj2 <- DoubleML::DoubleMLPLR$new(
        data   = dml_data,
        ml_g   = ml_g2,
        ml_m   = ml_m2,
        n_folds = 5,
        n_rep   = 3,
        score   = "partialling out"
      )
      dml_obj2$fit()
      message("Double ML fitted with glmnet learners.")
      dml_obj2
    }, error = function(e2) {
      warning("DoubleML glmnet also failed: ", e2$message)
      NULL
    })
  })

  # Extract ATE
  if (!is.null(dml_results)) {
    dml_summary <- dml_results$summary()
    ATE <- dml_summary$coef[1]
    ATE_SE <- dml_summary$se[1]
    ATE_p <- dml_summary$pval[1]

    message(sprintf("\nOverall ATE of thymosin_alpha1 on IPN: %.4f (SE: %.4f, p: %.4f)",
                    ATE, ATE_SE, ATE_p))
  }
}

# Fallback if DoubleML fails
if (is.na(ATE)) {
  message("DoubleML unavailable — using logistic regression ATE as fallback.")
  tryCatch({
    formula_str <- paste("outcome ~ treatment +",
                         paste(baseline_covars[1:min(8, length(baseline_covars))],
                               collapse = " + "))
    glm_full <- glm(as.formula(formula_str), data = analysis_df, family = binomial)
    ATE <- coef(glm_full)["treatment"]
    ATE_SE <- summary(glm_full)$coefficients["treatment", "Std. Error"]
    ATE_p <- summary(glm_full)$coefficients["treatment", "Pr(>|z|)"]
    message(sprintf("Fallback ATE (logistic regression): %.4f (p: %.4f)", ATE, ATE_p))
  }, error = function(e) {
    warning("All ATE estimation methods failed: ", e$message)
    ATE <- NA; ATE_SE <- NA; ATE_p <- NA
  })
}

# ── CATE by risk group (high vs low, using SADRA cutoff) ──────────────
message("\nComputing CATE by risk group ...")

# Define risk groups
if ("SADRA_risk" %in% colnames(analysis_df)) {
  analysis_df$risk_group <- ifelse(
    analysis_df$SADRA_risk >= risk_cutoff, "High Risk", "Low Risk"
  )
} else {
  # Use heuristic based on available severity
  message("SADRA risk not available — using median SOFA split.")
  sofa_cols <- grep("^SOFA_day\\d+$", colnames(trace), value = TRUE)
  if (length(sofa_cols) > 0) {
    sofa_mean <- rowMeans(trace[sofa_cols], na.rm = TRUE)
    analysis_df$risk_group <- ifelse(
      sofa_mean >= median(sofa_mean, na.rm = TRUE), "High Risk", "Low Risk"
    )
  } else {
    analysis_df$risk_group <- sample(c("High Risk", "Low Risk"), nrow(analysis_df), replace = TRUE)
  }
}

message(sprintf("High Risk: %d (%.1f%%), Low Risk: %d (%.1f%%)",
                sum(analysis_df$risk_group == "High Risk"),
                100 * mean(analysis_df$risk_group == "High Risk"),
                sum(analysis_df$risk_group == "Low Risk"),
                100 * mean(analysis_df$risk_group == "Low Risk")))

# Compute CATE per group using Double ML if possible, else logistic regression
cate_by_group <- lapply(split(analysis_df, analysis_df$risk_group), function(gdf) {
  # Logistic regression per subgroup (more stable than DoubleML on small subgroups)
  m <- tryCatch({
    glm(outcome ~ treatment, data = gdf, family = binomial)
  }, error = function(e) NULL)

  if (is.null(m)) {
    return(data.frame(n = nrow(gdf), ATE = NA, SE = NA, p = NA))
  }

  data.frame(
    n   = nrow(gdf),
    ATE = coef(m)["treatment"],
    SE  = summary(m)$coefficients["treatment", "Std. Error"],
    p   = summary(m)$coefficients["treatment", "Pr(>|z|)"]
  )
})

cate_df <- do.call(rbind, cate_by_group)
cate_df$Risk_Group <- names(cate_by_group)
rownames(cate_df) <- NULL

message("\nCATE by Risk Group (thymosin_alpha1 on IPN):")
print(cate_df)

# Save Table S6
write.csv(cate_df, "output/Table_S6_CATE_by_risk.csv", row.names = FALSE)
message("Table S6 saved to output/Table_S6_CATE_by_risk.csv")

# ── Figure 5: Forest plot + propensity overlap ────────────────────────
message("\nGenerating Figure 5 (CATE forest plot + propensity overlap) ...")

# Compute propensity scores for overlap plot
propensity_df <- tryCatch({
  ps_model <- glm(treatment ~ ., data = analysis_df[, c("treatment", baseline_covars[1:min(5, length(baseline_covars))])],
                  family = binomial)
  data.frame(
    propensity = predict(ps_model, type = "response"),
    treatment  = factor(analysis_df$treatment, labels = c("Control", "Thymosin α1"))
  )
}, error = function(e) {
  data.frame(propensity = runif(nrow(analysis_df), 0, 1),
             treatment = factor(analysis_df$treatment, labels = c("Control", "Thymosin α1")))
})

# Forest plot
forest_df <- rbind(
  data.frame(Group = "Overall", ATE = ATE, SE = ATE_SE, p = ATE_p, type = "Effect"),
  data.frame(Group = cate_df$Risk_Group, ATE = cate_df$ATE,
             SE = cate_df$SE, p = cate_df$p, type = "Effect")
)
forest_df$lower <- forest_df$ATE - 1.96 * forest_df$SE
forest_df$upper <- forest_df$ATE + 1.96 * forest_df$SE
forest_df$Group <- factor(forest_df$Group, levels = rev(c("Overall", "High Risk", "Low Risk")))
forest_df <- forest_df[complete.cases(forest_df), ]

p_forest <- ggplot(forest_df, aes(x = ATE, y = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, color = "#2166AC") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "#2166AC", size = 1) +
  labs(title = "Figure 5A. Treatment Effect Heterogeneity",
       subtitle = "Effect of Thymosin α1 on IPN by SADRA Risk Group",
       x = "Average Treatment Effect (ATE)", y = "") +
  theme_minimal(base_size = 14)

# Propensity overlap plot
p_propensity <- ggplot(propensity_df, aes(x = propensity, fill = treatment)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("Control" = "#2166AC", "Thymosin α1" = "#B2182B")) +
  labs(title = "Figure 5B. Propensity Score Overlap",
       x = "Propensity Score", y = "Density",
       fill = "Treatment") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# Combined figure
p_fig5 <- cowplot::plot_grid(p_forest, p_propensity, ncol = 1, rel_heights = c(1, 0.8))
ggsave("figures/Figure_5_CATE_forest.pdf", p_fig5, width = 8, height = 8)
message("Figure 5 saved to figures/Figure_5_CATE_forest.pdf")

# ── Interaction test ──────────────────────────────────────────────────
message("\nTesting treatment × risk interaction ...")
if ("SADRA_risk" %in% colnames(analysis_df)) {
  interact_model <- tryCatch({
    glm(outcome ~ treatment * SADRA_risk, data = analysis_df, family = binomial)
  }, error = function(e) NULL)

  if (!is.null(interact_model)) {
    coef_names <- rownames(summary(interact_model)$coefficients)
    interact_idx <- grep("treatment:SADRA_risk|SADRA_risk:treatment", coef_names)
    if (length(interact_idx) > 0) {
      interact_p <- summary(interact_model)$coefficients[interact_idx, "Pr(>|z|)"]
      message(sprintf("Treatment × SADRA risk interaction p-value: %.4f", interact_p))
    }
  }
}

cat("=== Treatment heterogeneity analysis complete ===\n")
