# 05_performance.R — Time-dependent AUC (cluster bootstrap), calibration, Youden cutoff, KM curves
cat("\n=== Step 5: Performance Evaluation ===\n")

jm_fit      <- readRDS("output/joint_model.rds")
selected_features <- readRDS("output/selected_features.rds")
train       <- readRDS("data/train_processed.rds")
train_long  <- readRDS("data/train_long.rds")
prosp       <- readRDS("data/prosp_processed.rds")
prosp_long  <- readRDS("data/prosp_long.rds")

dir.create("figures", showWarnings = FALSE)
eval_times <- c(1, 3, 5, 7, 14, 30, 60, 90)

# ── Helper: build model formulas from jm fit ─────────────────────────
extract_jm_formulas <- function(jm_obj, train_long, surv_data) {
  baseline_covars <- selected_features
  if (!exists("baseline_covars") || length(baseline_covars) < 3) {
    baseline_covars <- intersect(c("Cr", "HR", "RR", "CTSI", "ALB"), colnames(surv_data))
  }
  baseline_covars <- intersect(baseline_covars, colnames(surv_data))

  lme_formula_str <- paste0("SOFA ~ ns(day, df = 2) + ", paste(baseline_covars, collapse = " + "))
  cox_formula_str <- paste0("Surv(event_time, event_type) ~ ", paste(baseline_covars, collapse = " + "))

  list(
    lme_formula = as.formula(lme_formula_str),
    cox_formula = as.formula(cox_formula_str),
    random_formula = ~ day | ID,
    baseline_covars = baseline_covars
  )
}

# ── Refitting bootstrap helper ───────────────────────────────────────
refit_jm_bootstrap <- function(boot_ids, orig_long, orig_surv, formulas, seed_base) {
  # Build bootstrap long dataset with unique pseudo-IDs
  boot_rows <- do.call(rbind, lapply(seq_along(boot_ids), function(i) {
    rows <- orig_long[orig_long$ID == boot_ids[i], , drop = FALSE]
    rows$ID <- i
    rows
  }))
  boot_surv <- boot_rows[!duplicated(boot_rows$ID), ]
  boot_surv$event_time <- pmax(boot_surv$event_time, 0.1)

  lme_boot <- tryCatch({
    nlme::lme(fixed = formulas$lme_formula, data = boot_rows,
              random = formulas$random_formula,
              control = nlme::lmeControl(opt = "optim", maxIter = 100),
              na.action = na.omit)
  }, error = function(e) {
    nlme::lme(fixed = formulas$lme_formula, data = boot_rows,
              random = ~ 1 | ID,
              control = nlme::lmeControl(opt = "optim", maxIter = 100),
              na.action = na.omit)
  })

  cox_boot <- coxph(formulas$cox_formula, data = boot_surv, model = TRUE, x = TRUE)

  jm_boot <- JMbayes2::jm(
    Surv_object = cox_boot,
    Mixed_objects = lme_boot,
    time_var = "day",
    functional_forms = ~ value(SOFA),
    n_iter = 2000, n_burnin = 1000, n_thin = 1, n_chains = 1,
    cores = 1, seed = seed_base
  )

  list(jm = jm_boot, long = boot_rows, surv = boot_surv)
}

# ── Time-dependent AUC with tvROC and refitting bootstrap ────────────
message("\nComputing time-dependent AUC (tvROC + refitting bootstrap) ...")

compute_td_auc <- function(jm_obj, newdata_long, newdata_surv_raw,
                           times = eval_times, Thoriz = 90, n_boot = 200) {

  results <- data.frame(
    Time      = times,
    AUC       = NA_real_,
    AUC_SE    = NA_real_,
    AUC_lower = NA_real_,
    AUC_upper = NA_real_
  )

  # Point estimates via tvROC
  for (i in seq_along(times)) {
    t_i <- times[i]
    message(sprintf("  Point estimate at Tstart = %d days (Thoriz = %d) ...", t_i, Thoriz))

    roc_obj <- tryCatch({
      tvROC(jm_obj, newdata = newdata_long, Tstart = t_i,
            Thoriz = Thoriz, integrated = TRUE)
    }, error = function(e) {
      warning(sprintf("tvROC failed at Tstart %d: %s", t_i, e$message))
      NULL
    })

    if (!is.null(roc_obj)) {
      auc_val <- tryCatch(tvAUC(roc_obj)$auc, error = function(e) NA_real_)
      results$AUC[i] <- auc_val
      message(sprintf("    AUC = %.3f", auc_val))
    }
  }

  # Refitting bootstrap
  if (n_boot > 0 && !is.null(jm_obj)) {
    message(sprintf("\n  Refitting bootstrap (B = %d) ...", n_boot))
    formulas <- extract_jm_formulas(jm_obj, newdata_long, newdata_surv_raw)

    unique_ids <- unique(newdata_long$ID)
    n_patients <- length(unique_ids)

    # Parallel over bootstrap samples
    `%dofuture%` <- if (requireNamespace("doFuture", quietly = TRUE)) {
      doFuture::`%dofuture%`
    } else {
      foreach::`%do%`
    }

    boot_aucs <- foreach::foreach(
      b = seq_len(n_boot),
      .combine = rbind,
      .packages = c("JMbayes2", "nlme", "survival", "splines"),
      .errorhandling = "remove"
    ) %dofuture% {
      sampled_ids <- sample(unique_ids, n_patients, replace = TRUE)
      boot_fit <- refit_jm_bootstrap(sampled_ids, newdata_long,
                                      newdata_surv_raw, formulas, 2025 + b)

      sapply(times, function(t_i) {
        roc_b <- tryCatch({
          tvROC(boot_fit$jm, newdata = boot_fit$long,
                Tstart = t_i, Thoriz = Thoriz, integrated = TRUE)
        }, error = function(e) NULL)
        if (!is.null(roc_b)) {
          tryCatch(tvAUC(roc_b)$auc, error = function(e) NA_real_)
        } else NA_real_
      })
    }

    if (!is.null(boot_aucs) && nrow(boot_aucs) > 10) {
      for (j in seq_along(times)) {
        col_j <- boot_aucs[, j]
        col_j <- col_j[!is.na(col_j)]
        if (length(col_j) > 10) {
          results$AUC_SE[j]    <- sd(col_j) / sqrt(length(col_j))
          results$AUC_lower[j] <- quantile(col_j, 0.025, na.rm = TRUE)
          results$AUC_upper[j] <- quantile(col_j, 0.975, na.rm = TRUE)
        }
      }
      message(sprintf("  Bootstrap completed: %d valid samples", nrow(boot_aucs)))
    }
  }

  return(results)
}

# ── Helper: get one risk prediction per subject at a given horizon ────
predict_single_time <- function(jm_obj, newdata_long, t_horizon) {
  preds <- predict(jm_obj, newdata = newdata_long,
                   times = t_horizon, process = "event", type = "survival")
  risk <- 1 - preds$pred[, 1]
  unique_ids <- unique(newdata_long$ID)
  surv_ids <- newdata_long[!duplicated(newdata_long$ID), ]
  surv_ids$predicted_risk <- NA_real_
  if (length(unique_ids) == length(risk)) {
    surv_ids$predicted_risk <- risk
  }
  surv_ids$ID <- unique_ids
  return(surv_ids)
}

if (!is.null(jm_fit)) {
  auc_results <- compute_td_auc(jm_fit, prosp_long, prosp,
                                times = eval_times, n_boot = 200)

  write.csv(auc_results, "output/AUC_results.csv", row.names = FALSE)
  message("AUC results saved to output/AUC_results.csv")

  # Integrated AUC
  iauc <- compute_auc_jm(data.frame(time = auc_results$Time, AUC = auc_results$AUC))
  message(sprintf("Integrated AUC (days 1-90): %.3f", iauc))
} else {
  message("Skipping AUC: joint model not available.")
}

# ── Calibration plots (Figure 2) ──────────────────────────────────────
message("\nGenerating calibration plots ...")

if (!is.null(jm_fit)) {
  tryCatch({
    # Calibration at multiple time points
    cal_times <- c(3, 5, 7)
    cal_plots <- list()

    for (t_cal in cal_times) {
      surv_df <- predict_single_time(jm_fit, prosp_long, t_cal)
      surv_raw <- prosp[!duplicated(prosp$ID), ]
      surv_df$event_type <- surv_raw$event_type[match(surv_df$ID, surv_raw$ID)]

      valid <- complete.cases(surv_df$predicted_risk, surv_df$event_type)
      cal_plots[[as.character(t_cal)]] <- plot_calibration(
        surv_df$predicted_risk[valid],
        surv_df$event_type[valid],
        title = paste0("Calibration at Day ", t_cal)
      )
    }

    p_cal <- cowplot::plot_grid(
      cal_plots[["3"]], cal_plots[["5"]], cal_plots[["7"]],
      ncol = 3, labels = c("A", "B", "C")
    )

    ggsave("figures/Figure_2_calibration.pdf", p_cal,
           width = 15, height = 5, device = "pdf")
    message("Figure 2 (calibration) saved.")
  }, error = function(e) {
    warning("Calibration plot failed: ", e$message)
  })
}

# ── Youden index → risk cutoff ─────────────────────────────────────────
message("\nComputing optimal risk cutoff (Youden index, day 7) ...")

if (!is.null(jm_fit)) {
  tryCatch({
    roc_obj <- tvROC(jm_fit, newdata = prosp_long,
                     Tstart = 7, Thoriz = 90, integrated = FALSE)
    sens <- roc_obj$TPR
    spec <- 1 - roc_obj$FPR
    youden_idx <- which.max(sens + spec - 1)
    optimal_cutoff <- roc_obj$thresholds[youden_idx]
    optimal_sens  <- sens[youden_idx]
    optimal_spec  <- spec[youden_idx]

    message(sprintf("Optimal cutoff (Youden): %.4f (Sens: %.3f, Spec: %.3f)",
                    optimal_cutoff, optimal_sens, optimal_spec))
    saveRDS(optimal_cutoff, "output/risk_cutoff.rds")
  }, error = function(e) {
    warning("ROC analysis failed: ", e$message, ". Using default cutoff 0.32.")
    optimal_cutoff <- 0.32
    saveRDS(optimal_cutoff, "output/risk_cutoff.rds")
  })
} else {
  optimal_cutoff <- 0.32
  saveRDS(optimal_cutoff, "output/risk_cutoff.rds")
}

# ── KM curves stratified by risk group (Figure 3) ─────────────────────
message("\nGenerating Kaplan-Meier curves (Figure 3) ...")

if (!is.null(jm_fit)) {
  tryCatch({
    surv_df_7 <- predict_single_time(jm_fit, prosp_long, 7)
    surv_raw <- prosp[!duplicated(prosp$ID), ]
    surv_df_7$event_type <- surv_raw$event_type[match(surv_df_7$ID, surv_raw$ID)]
    surv_df_7$event_time <- surv_raw$event_time[match(surv_df_7$ID, surv_raw$ID)]

    valid <- complete.cases(surv_df_7$predicted_risk, surv_df_7$event_type)
    surv_plot <- surv_df_7[valid, ]
    surv_plot$risk_group <- ifelse(
      surv_plot$predicted_risk >= optimal_cutoff, "High Risk", "Low Risk"
    )
    surv_plot$event_time <- pmax(surv_plot$event_time, 0.1)

    km_fit <- survfit(Surv(event_time, event_type) ~ risk_group, data = surv_plot)

    # Log-rank test
    lr_test <- survdiff(Surv(event_time, event_type) ~ risk_group, data = surv_plot)
    lr_pval <- 1 - pchisq(lr_test$chisq, df = 1)
    message(sprintf("Log-rank test p-value: %.6f", lr_pval))

    km_plot <- ggsurvplot(
      km_fit,
      data = surv_plot,
      pval = TRUE,
      pval.method = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      palette = c("#2166AC", "#B2182B"),
      xlab = "Time (days)",
      ylab = "Event-free Survival",
      title = "Figure 3. KM Curves by SADRA Risk Group",
      legend.title = "Risk Group",
      legend.labs = c("High Risk", "Low Risk"),
      ggtheme = theme_minimal(base_size = 13)
    )

    pdf("figures/Figure_3_KM_curves.pdf", width = 8, height = 8)
    print(km_plot)
    dev.off()
    message("Figure 3 (KM curves) saved.")
  }, error = function(e) {
    warning("KM curve generation failed: ", e$message)
  })
}

# ── AUC curve plot (Figure 2 combined) ─────────────────────────────────
message("\nGenerating time-dependent AUC curve ...")
if (exists("auc_results") && !all(is.na(auc_results$AUC))) {
  p_auc <- ggplot(auc_results, aes(x = Time, y = AUC)) +
    geom_ribbon(aes(ymin = AUC_lower, ymax = AUC_upper), fill = "#2166AC", alpha = 0.2) +
    geom_line(color = "#2166AC", size = 1.2) +
    geom_point(size = 3, color = "#2166AC") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    labs(title = "Time-dependent AUC (95% CI, Cluster Bootstrap)",
         x = "Time (days)", y = "AUC") +
    ylim(0.4, 1) +
    theme_minimal(base_size = 13)
  ggsave("figures/Figure_2_AUC_curve.pdf", p_auc, width = 7, height = 5)
  message("Figure 2 AUC curve saved.")
}

cat("=== Performance evaluation complete ===\n")
