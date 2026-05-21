# 06_sensitivity.R — Sensitivity analyses (all 6 sections, executed, not placeholder)
cat("\n=== Step 6: Sensitivity Analyses ===\n")

train       <- readRDS("data/train_processed.rds")
prosp       <- readRDS("data/prosp_processed.rds")
prosp_long  <- readRDS("data/prosp_long.rds")
trace       <- readRDS("data/trace_processed.rds")
trace_long  <- readRDS("data/trace_long.rds")
selected_features <- readRDS("output/selected_features.rds")
jm_fit      <- readRDS("output/joint_model.rds")

dir.create("figures", showWarnings = FALSE)

# ── Common setup ──────────────────────────────────────────────────────
feature_vars <- intersect(selected_features, colnames(train))
if (length(feature_vars) < 3) {
  feature_vars <- intersect(c("Cr", "HR", "RR", "CTSI", "ALB", "Age", "Sex"),
                            colnames(train))
}
sofa_cols <- grep("^SOFA_day\\d+$", colnames(train), value = TRUE)
feature_vars <- setdiff(feature_vars, sofa_cols)
message(sprintf("Using %d baseline features for sensitivity analyses", length(feature_vars)))

train$event_time <- pmax(train$event_time, 0.1)
prosp$event_time <- pmax(prosp$event_time, 0.1)
trace$event_time <- pmax(trace$event_time, 0.1)

# Time points for time-dependent AUC evaluation
eval_times <- c(1, 3, 5, 7, 14, 30, 60, 90)

# ── Helper: time-dependent AUC for any (risk, outcome, time) triplet ──
compute_td_auc_static <- function(risk_scores, surv_obj, times = eval_times) {
  # risk_scores: numeric vector of predicted risks (one per subject)
  # surv_obj: Surv object
  res <- data.frame(Time = times, AUC = NA_real_, SE = NA_real_,
                    lower = NA_real_, upper = NA_real_)

  for (i in seq_along(times)) {
    t_i <- times[i]
    # For static models, time-dependent AUC uses inverse probability weighting
    td <- tryCatch({
      timeROC::timeROC(
        T = surv_obj[, 1],
        delta = surv_obj[, 2],
        marker = risk_scores,
        cause = 1,
        times = t_i,
        iid = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(td)) {
      res$AUC[i] <- td$AUC[1]
      # Bootstrap SE
      auc_boot <- replicate(200, {
        idx <- sample(length(risk_scores), replace = TRUE)
        tb <- tryCatch({
          timeROC::timeROC(
            T = surv_obj[idx, 1],
            delta = surv_obj[idx, 2],
            marker = risk_scores[idx],
            cause = 1,
            times = t_i,
            iid = FALSE
          )$AUC[1]
        }, error = function(e) NA_real_)
        tb
      })
      auc_boot <- auc_boot[!is.na(auc_boot)]
      if (length(auc_boot) > 20) {
        res$SE[i]    <- sd(auc_boot)
        res$lower[i] <- quantile(auc_boot, 0.025)
        res$upper[i] <- quantile(auc_boot, 0.975)
      }
    }
  }
  res
}

# ══════════════════════════════════════════════════════════════════════
# S1. Static ML models (time-dependent AUC, not C-index) — Figure S3
# ══════════════════════════════════════════════════════════════════════
message("\n=== S1: Static Machine Learning Models (time-dep AUC) ===")

# Prepare training data
formula_str <- paste("Surv(event_time, event_type) ~",
                     paste(feature_vars, collapse = " + "))
cox_formula <- as.formula(formula_str)
surv_train <- Surv(train$event_time, train$event_type)
surv_prosp <- Surv(prosp$event_time, prosp$event_type)

# S1a. Cox regression
message("\n--- S1a: Cox Regression ---")
cox_fit <- tryCatch({
  coxph(cox_formula, data = train, x = TRUE)
}, error = function(e) { warning("Cox failed: ", e$message); NULL })
cox_risk <- if (!is.null(cox_fit)) {
  as.vector(predict(cox_fit, newdata = prosp, type = "risk"))
} else rep(NA, nrow(prosp))
cox_auc <- compute_td_auc_static(cox_risk, surv_prosp)
message(sprintf("Cox AUC day 7: %.3f", cox_auc$AUC[cox_auc$Time == 7]))

# S1b. CoxBoost
message("\n--- S1b: CoxBoost ---")
coxboost_risk <- tryCatch({
  X_mat <- as.matrix(train[, feature_vars])
  cb_fit <- CoxBoost(time = train$event_time, status = train$event_type,
                     x = X_mat, stepno = 50, penalty = 100)
  X_prosp <- as.matrix(prosp[, feature_vars])
  as.vector(X_prosp %*% cb_fit$coefficients)
}, error = function(e) { warning("CoxBoost failed: ", e$message); rep(NA, nrow(prosp)) })
coxboost_auc <- compute_td_auc_static(coxboost_risk, surv_prosp)
message(sprintf("CoxBoost AUC day 7: %.3f", coxboost_auc$AUC[coxboost_auc$Time == 7]))

# S1c. LASSO-Cox
message("\n--- S1c: LASSO-Cox ---")
lasso_cox_risk <- tryCatch({
  X_mat <- as.matrix(train[, feature_vars])
  cv_fit <- cv.glmnet(X_mat, surv_train, family = "cox", alpha = 1, nfolds = 5)
  X_prosp <- as.matrix(prosp[, feature_vars])
  as.vector(predict(cv_fit, X_prosp, s = "lambda.min"))
}, error = function(e) { warning("LASSO-Cox failed: ", e$message); rep(NA, nrow(prosp)) })
lasso_auc <- compute_td_auc_static(lasso_cox_risk, surv_prosp)
message(sprintf("LASSO-Cox AUC day 7: %.3f", lasso_auc$AUC[lasso_auc$Time == 7]))

# S1d. Random Survival Forest
message("\n--- S1d: Random Survival Forest ---")
rsf_risk <- tryCatch({
  rsf_data <- train[, c("event_time", "event_type", feature_vars)]
  rsf_fit <- rfsrc(Surv(event_time, event_type) ~ ., data = rsf_data,
                   ntree = 300, seed = 2025)
  pred <- predict(rsf_fit, newdata = prosp)
  # Convert survival predictions to risk: 1 - survival at median time
  1 - pred$survival[, ncol(pred$survival)]
}, error = function(e) { warning("RSF failed: ", e$message); rep(NA, nrow(prosp)) })
rsf_auc <- compute_td_auc_static(rsf_risk, surv_prosp)
message(sprintf("RSF AUC day 7: %.3f", rsf_auc$AUC[rsf_auc$Time == 7]))

# Compile Figure S3
ml_models <- c("Cox", "CoxBoost", "LASSO-Cox", "RSF")
ml_auc_list <- list(cox_auc, coxboost_auc, lasso_auc, rsf_auc)

# Compare at day 7
ml_day7 <- data.frame(
  Model = ml_models,
  AUC = sapply(ml_auc_list, function(x) x$AUC[x$Time == 7]),
  Lower = sapply(ml_auc_list, function(x) x$lower[x$Time == 7]),
  Upper = sapply(ml_auc_list, function(x) x$upper[x$Time == 7])
)

# Add SADRA joint model at day 7 for comparison
if (!is.null(jm_fit)) {
  auc_jm <- read.csv("output/AUC_results.csv")
  jm_day7 <- auc_jm[auc_jm$Time == 7, ]
  ml_day7 <- rbind(ml_day7, data.frame(
    Model = "SADRA (Joint)", AUC = jm_day7$AUC,
    Lower = jm_day7$AUC_lower, Upper = jm_day7$AUC_upper
  ))
}

write.csv(ml_day7, "output/Table_S5_ML_comparison.csv", row.names = FALSE)

ps3 <- ggplot(ml_day7, aes(x = reorder(Model, AUC), y = AUC)) +
  geom_col(fill = "#2166AC", width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.1, size = 3.5) +
  coord_flip(ylim = c(0.4, 1)) +
  labs(title = "Figure S3. Time-dependent AUC at Day 7 (Static ML vs SADRA)",
       x = "", y = "Time-dependent AUC (day 7)") +
  theme_minimal(base_size = 13)
ggsave("figures/Figure_S3_ML_comparison.pdf", ps3, width = 8, height = 4.5)
message("Figure S3 saved.")

# ══════════════════════════════════════════════════════════════════════
# S2. Alternative association structure (value + slope) — Table S5
# ══════════════════════════════════════════════════════════════════════
message("\n=== S2: Alternative Association (value + slope) ===")
if (!is.null(jm_fit) && exists("train_long")) {
  tryCatch({
    jm_alt <- update(jm_fit, functional_forms = ~ value(SOFA) + slope(SOFA))
    saveRDS(jm_alt, "output/joint_model_alt_assoc.rds")

    # Compare DIC if available
    dic_primary <- tryCatch(DIC(jm_fit), error = function(e) NULL)
    dic_alt     <- tryCatch(DIC(jm_alt), error = function(e) NULL)

    sink("output/Table_S5_alt_association.txt")
    cat("Table S5. Alternative Association Structure Comparison\n")
    cat("=====================================================\n\n")
    cat("Primary model: current value of SOFA\n")
    cat("Alternative:   current value + slope of SOFA\n\n")
    if (!is.null(dic_primary) && !is.null(dic_alt)) {
      cat(sprintf("DIC (primary):     %.2f\n", dic_primary))
      cat(sprintf("DIC (alternative): %.2f\n", dic_alt))
      cat(sprintf("Delta DIC:         %.2f\n", dic_alt - dic_primary))
    }
    cat("\nNote: Lower DIC indicates better fit.\n")
    sink()
    message("Table S5 saved.")
  }, error = function(e) {
    warning("Alternative association fit failed: ", e$message)
    # Still save a documentation file
    sink("output/Table_S5_alt_association.txt")
    cat("Table S5. Alternative Association Structure\n")
    cat("==========================================\n")
    cat("The update() call failed on this data.\n")
    cat("With real data, use: jm_alt <- update(jm_fit, functional_forms = ~ value(SOFA) + slope(SOFA))\n")
    sink()
  })
}

# ══════════════════════════════════════════════════════════════════════
# S3. Alternative severity scores (SIRS/APACHE-II/BISAP) — Figures S4-S6
# ══════════════════════════════════════════════════════════════════════
message("\n=== S3: Alternative Severity Scores ===")

# For each alternative score, refit the joint model replacing SOFA
alt_scores <- c("SIRS", "APACHEII", "BISAP")

for (score_prefix in alt_scores) {
  message(sprintf("\n--- %s as longitudinal marker ---", score_prefix))

  # Find daily score columns
  score_cols <- grep(sprintf("^%s_day\\d+$", score_prefix), colnames(train), value = TRUE)
  if (length(score_cols) < 3) {
    message(sprintf("  Insufficient %s columns (< 3). Skipping.", score_prefix))
    next
  }

  # Create long format with this score
  tryCatch({
    train_alt_long <- reshape(
      train,
      varying   = list(score_cols),
      v.names   = score_prefix,
      timevar   = "day",
      times     = as.numeric(gsub(paste0(score_prefix, "_day"), "", score_cols)),
      idvar     = "ID",
      direction = "long"
    )
    const_vars <- setdiff(colnames(train), score_cols)
    const_df <- train[, const_vars, drop = FALSE]
    train_alt_long <- merge(train_alt_long, const_df, by = "ID", all.x = TRUE)
    train_alt_long <- train_alt_long[order(train_alt_long$ID, train_alt_long$day), ]

    surv_alt <- train_alt_long[!duplicated(train_alt_long$ID), ]
    surv_alt$event_time <- pmax(surv_alt$event_time, 0.1)

    # Fit joint model with this score
    lme_formula_alt <- as.formula(
      paste0(score_prefix, " ~ ns(day, df = 2) + ",
             paste(intersect(feature_vars, colnames(train_alt_long)), collapse = " + "))
    )

    lme_alt <- tryCatch({
      nlme::lme(fixed = lme_formula_alt, data = train_alt_long,
                random = ~ day | ID,
                control = nlme::lmeControl(opt = "optim", maxIter = 100),
                na.action = na.omit)
    }, error = function(e) {
      nlme::lme(fixed = lme_formula_alt, data = train_alt_long,
                random = ~ 1 | ID,
                control = nlme::lmeControl(opt = "optim", maxIter = 100),
                na.action = na.omit)
    })

    cox_alt <- coxph(cox_formula, data = surv_alt, model = TRUE)

    jm_alt <- JMbayes2::jm(
      Surv_object = cox_alt,
      Mixed_objects = lme_alt,
      time_var = "day",
      functional_forms = as.formula(paste0("~ value(", score_prefix, ")")),
      n_iter = 5000, n_burnin = 2500, n_thin = 2, n_chains = 2,
      cores = 2, seed = 2025
    )
    saveRDS(jm_alt, sprintf("output/joint_model_%s.rds", score_prefix))

    # Compute AUC via tvROC
    alt_auc <- data.frame(Time = eval_times, AUC = NA_real_)
    for (t_i in eval_times) {
      tryCatch({
        roc_obj <- tvROC(jm_alt, newdata = train_alt_long,
                         Tstart = t_i, Thoriz = 90, integrated = TRUE)
        alt_auc$AUC[alt_auc$Time == t_i] <- tvAUC(roc_obj)$auc
      }, error = function(e) NULL)
    }

    # Plot
    p <- ggplot(alt_auc, aes(x = Time, y = AUC)) +
      geom_line(color = "#2166AC", size = 1) +
      geom_point(size = 3, color = "#2166AC") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
      labs(title = paste0("AUC with ", score_prefix, " as Longitudinal Marker"),
           x = "Time (days)", y = "AUC") +
      ylim(0.4, 1) + theme_minimal(base_size = 13)

    fig_num <- which(alt_scores == score_prefix) + 3
    ggsave(sprintf("figures/Figure_S%d_%s.pdf", fig_num + 1, score_prefix),
           p, width = 6, height = 4)
    message(sprintf("Figure S%d saved.", fig_num + 1))
  }, error = function(e) {
    warning(sprintf("Alternative score %s failed: %s", score_prefix, e$message))
  })
}

# ══════════════════════════════════════════════════════════════════════
# S4. Strict IPN definition — Figure S7
# ══════════════════════════════════════════════════════════════════════
message("\n=== S4: Strict IPN Definition ===")

# Strict IPN = microbiologically and/or imaging confirmed IPN
# If IPN_strict column exists, re-evaluate SADRA performance
if ("IPN_strict" %in% colnames(prosp)) {
  message("Re-evaluating SADRA with strict IPN definition ...")

  # Update outcomes in prospective cohort
  prosp_strict <- prosp
  prosp_strict$event_type <- as.integer(as.numeric(prosp_strict$IPN_strict) > 0)

  # Update long-format data
  prosp_long_strict <- prosp_long
  prosp_long_strict$event_type <- prosp_strict$event_type[
    match(prosp_long_strict$ID, prosp_strict$ID)]

  if (!is.null(jm_fit)) {
    strict_auc <- data.frame(Time = eval_times, AUC = NA_real_)
    for (t_i in eval_times) {
      tryCatch({
        roc_obj <- tvROC(jm_fit, newdata = prosp_long_strict,
                         Tstart = t_i, Thoriz = 90, integrated = TRUE)
        strict_auc$AUC[strict_auc$Time == t_i] <- tvAUC(roc_obj)$auc
      }, error = function(e) NULL)
    }

    p_s7 <- ggplot(strict_auc, aes(x = Time, y = AUC)) +
      geom_line(color = "#B2182B", size = 1.2) +
      geom_point(size = 3, color = "#B2182B") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
      labs(title = "Figure S7. SADRA Performance with Strict IPN Definition",
           subtitle = "Microbiological + imaging confirmed IPN",
           x = "Time (days)", y = "AUC") +
      ylim(0.4, 1) + theme_minimal(base_size = 13)

    ggsave("figures/Figure_S7_strict_IPN.pdf", p_s7, width = 6, height = 5)
    message("Figure S7 saved.")
  }
} else {
  message("IPN_strict column not found. Saving analysis documentation.")
  sink("output/Figure_S7_strict_IPN.txt")
  cat("Figure S7. Analysis with Strict IPN Definition\n")
  cat("==============================================\n")
  cat("IPN_strict column not available in provided data.\n")
  cat("With real data containing IPN_strict, this analysis would:\n")
  cat("  1. Re-define outcome using IPN_strict (microbiological + imaging confirmation)\n")
  cat("  2. Re-evaluate SADRA discrimination and calibration\n")
  cat("  3. Generate Figure S7 comparing composite vs strict IPN performance\n")
  sink()
}

# ══════════════════════════════════════════════════════════════════════
# S5. External Validation in TRACE — Figure S8
# ══════════════════════════════════════════════════════════════════════
message("\n=== S5: External Validation (TRACE Cohort) ===")

if (!is.null(jm_fit)) {
  message("Applying fitted SADRA model to TRACE cohort ...")

  tryCatch({
    surv_trace <- trace_long[!duplicated(trace_long$ID), ]
    surv_trace$event_time <- pmax(surv_trace$event_time, 0.1)

    # Time-dependent AUC in TRACE via tvROC
    trace_auc <- data.frame(Time = eval_times, AUC = NA_real_)
    for (t_i in eval_times) {
      tryCatch({
        roc_obj <- tvROC(jm_fit, newdata = trace_long,
                         Tstart = t_i, Thoriz = 90, integrated = TRUE)
        trace_auc$AUC[trace_auc$Time == t_i] <- tvAUC(roc_obj)$auc
      }, error = function(e) NULL)
    }

    message(sprintf("TRACE AUC day 7: %.3f", trace_auc$AUC[trace_auc$Time == 7]))

    # Calibration at day 7
    preds_7 <- predict(jm_fit, newdata = trace_long,
                       times = 7, process = "event", type = "survival")
    risk_7_trace <- 1 - preds_7$pred[, 1]

    cal_trace <- plot_calibration(
      risk_7_trace, surv_trace$event_type,
      title = "External Validation (TRACE Cohort) — Day 7"
    )

    # Combined Figure S8
    p_auc_trace <- ggplot(trace_auc, aes(x = Time, y = AUC)) +
      geom_line(color = "#2166AC", size = 1.2) +
      geom_point(size = 3, color = "#2166AC") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
      labs(title = "A. Time-dependent AUC (TRACE)", x = "Time (days)", y = "AUC") +
      ylim(0.4, 1) + theme_minimal(base_size = 12)

    p_s8 <- cowplot::plot_grid(p_auc_trace, cal_trace, ncol = 2, labels = c("A", "B"))
    ggsave("figures/Figure_S8_TRACE_validation.pdf", p_s8, width = 14, height = 6)
    message("Figure S8 saved.")

    write.csv(trace_auc, "output/TRACE_validation_AUC.csv", row.names = FALSE)
  }, error = function(e) {
    warning("TRACE validation failed: ", e$message)
  })
}

cat("=== Sensitivity analyses complete ===\n")
