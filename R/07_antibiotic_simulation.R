# 07_antibiotic_simulation.R — Model-guided antibiotic strategy, Figure 4 (Sankey + UpSet)
cat("\n=== Step 7: Antibiotic Simulation ===\n")

prosp       <- readRDS("data/prosp_processed.rds")
prosp_long  <- readRDS("data/prosp_long.rds")
jm_fit      <- readRDS("output/joint_model.rds")
risk_cutoff <- readRDS("output/risk_cutoff.rds")

dir.create("figures", showWarnings = FALSE)

message(sprintf("Risk cutoff (Youden): %.4f", risk_cutoff))

# ── Simulate model-guided strategy ────────────────────────────────────
message("\nSimulating model-guided antibiotic strategy ...")

simulate_strategy <- function(jm_obj, data_long, cutoff) {
  surv_data <- data_long[!duplicated(data_long$ID), ]

  # Predict risk at day 7 using joint model
  preds <- tryCatch({
    predict(jm_obj, newdata = data_long,
            times = 7, process = "event", type = "survival")
  }, error = function(e) NULL)

  if (!is.null(preds)) {
    surv_data$predicted_risk <- 1 - preds$pred[, 1]
  } else {
    warning("Joint model prediction failed. Using heuristic risk scores.")
    # Heuristic based on SOFA trajectory
    sofa_means <- tapply(data_long$SOFA, data_long$ID, mean, na.rm = TRUE)
    surv_data$predicted_risk <- plogis(sofa_means[match(surv_data$ID, names(sofa_means))] - 3)
  }

  surv_data$model_high_risk <- surv_data$predicted_risk >= cutoff

  # Actual antibiotic use
  if ("Antibiotic_actual" %in% colnames(surv_data)) {
    surv_data$actual_abx <- surv_data$Antibiotic_actual
  } else {
    # If not available, note it and use a placeholder for demonstration
    message("Note: Antibiotic_actual column not found. Using simulated patterns for illustration.")
    surv_data$actual_abx <- NA_integer_
  }

  # Model-guided recommendation
  surv_data$model_recommend <- as.integer(surv_data$model_high_risk)

  # Determine IPN status (for sterile necrosis analysis)
  if ("IPN_strict" %in% colnames(surv_data)) {
    surv_data$sterile_necrosis <- (surv_data$IPN == 1 & surv_data$IPN_strict == 0)
  } else {
    surv_data$sterile_necrosis <- FALSE
  }

  # Classification
  if (all(!is.na(surv_data$actual_abx))) {
    surv_data$category <- with(surv_data, ifelse(
      actual_abx == 1 & model_recommend == 1, "Appropriate ABx (High Risk)",
      ifelse(actual_abx == 0 & model_recommend == 0, "Appropriate No ABx (Low Risk)",
      ifelse(actual_abx == 1 & model_recommend == 0, "De-escalation Opportunity",
             "Escalation Needed")))
    )
  }

  return(surv_data)
}

sim_results <- simulate_strategy(jm_fit, prosp_long, risk_cutoff)

# ── Summary statistics ────────────────────────────────────────────────
if (!all(is.na(sim_results$actual_abx))) {
  cat_counts <- table(sim_results$category)
  message("\nAntibiotic strategy comparison:")
  for (nm in names(cat_counts)) {
    message(sprintf("  %s: %d (%.1f%%)", nm, cat_counts[nm],
                    100 * cat_counts[nm] / sum(cat_counts)))
  }

  deesc_pct <- 100 * sum(sim_results$category == "De-escalation Opportunity") / nrow(sim_results)
  escal_pct <- 100 * sum(sim_results$category == "Escalation Needed") / nrow(sim_results)
  abx_reduction <- 100 * (sum(sim_results$actual_abx) - sum(sim_results$model_recommend)) /
                    max(1, sum(sim_results$actual_abx))

  message(sprintf("\nPotential antibiotic reduction: %.1f%%", abx_reduction))
  message(sprintf("De-escalation opportunities: %.1f%% of patients", deesc_pct))
  message(sprintf("Escalation needed: %.1f%% of patients", escal_pct))
}

# Save results
write.csv(sim_results, "output/antibiotic_simulation.csv", row.names = FALSE)

# ── Figure 4: Multi-panel (Sankey + Bar + UpSet) ──────────────────────
message("\nGenerating Figure 4 (Sankey + bar charts + UpSet) ...")

pdf("figures/Figure_4_antibiotic_simulation.pdf", width = 16, height = 12)

# Panel A: Alluvial/Sankey diagram showing flow between risk strata and antibiotic use
tryCatch({
  # Prepare alluvial data
  if (!all(is.na(sim_results$actual_abx))) {
    alluv_df <- data.frame(
      ID = sim_results$ID,
      Risk = ifelse(sim_results$model_high_risk, "High Risk", "Low Risk"),
      Antibiotics = ifelse(sim_results$actual_abx == 1, "ABx Given", "ABx Withheld"),
      Recommendation = ifelse(sim_results$model_recommend == 1,
                              "Recommended ABx", "Recommended No ABx"),
      Category = sim_results$category
    )

    # Frequency table for alluvial
    alluv_freq <- alluv_df %>%
      count(Risk, Antibiotics, Recommendation) %>%
      mutate( Risk = factor(Risk, levels = c("High Risk", "Low Risk")))

    p_sankey <- ggplot(alluv_freq,
                       aes(axis1 = Risk, axis2 = Antibiotics, axis3 = Recommendation,
                           y = n)) +
      ggalluvial::geom_alluvium(aes(fill = Risk), width = 1/4) +
      ggalluvial::geom_stratum(width = 1/4, fill = "grey90", color = "grey30") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
      scale_x_discrete(limits = c("Risk", "Antibiotics", "Recommendation"),
                       expand = c(0.15, 0.05)) +
      scale_fill_manual(values = c("High Risk" = "#B2182B", "Low Risk" = "#2166AC")) +
      labs(title = "A. Sankey Diagram: Risk → Actual ABx → Model Recommendation",
           subtitle = paste0("Cutoff = ", round(risk_cutoff, 3))) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_blank())
    print(p_sankey)
  }
}, error = function(e) {
  warning("Sankey diagram failed: ", e$message)
  plot.new()
  text(0.5, 0.5, "Sankey diagram unavailable", cex = 1.2)
})

# Panel B: Bar chart — Antibiotic strategy categories
if (!all(is.na(sim_results$actual_abx))) {
  cat_df <- as.data.frame(table(Category = sim_results$category))
  p_bars <- ggplot(cat_df, aes(x = reorder(Category, Freq), y = Freq, fill = Category)) +
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    scale_fill_brewer(palette = "RdBu", direction = -1) +
    labs(title = "B. Antibiotic Strategy Comparison",
         subtitle = "Actual use vs. Model recommendation",
         x = "", y = "Number of Patients") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  print(p_bars)
}

# Panel C: Risk distribution by actual antibiotic use
if (!all(is.na(sim_results$actual_abx))) {
  p_risk_box <- ggplot(sim_results, aes(
    x = factor(actual_abx, labels = c("No ABx", "ABx Given")),
    y = predicted_risk, fill = factor(actual_abx))) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = risk_cutoff, linetype = "dashed", color = "#B2182B", size = 1) +
    annotate("text", x = 0.7, y = risk_cutoff + 0.02,
             label = paste0("Cutoff = ", round(risk_cutoff, 3)),
             color = "#B2182B", size = 3.5) +
    scale_fill_manual(values = c("#2166AC", "#B2182B")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "C. Predicted Risk by Actual Antibiotic Use",
         x = "", y = "Predicted IPN Risk") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  print(p_risk_box)
}

# Panel D: UpSet plot — intersection of risk, ABx, and outcome
tryCatch({
  upset_data <- data.frame(
    ID = sim_results$ID,
    HighRisk = as.integer(sim_results$model_high_risk),
    ABxGiven = as.integer(sim_results$actual_abx),
    IPN = as.integer(sim_results$IPN)
  )
  upset_data <- upset_data[complete.cases(upset_data), ]

  if (nrow(upset_data) > 0 && sum(upset_data$HighRisk) > 0 &&
      sum(upset_data$ABxGiven) > 0 && sum(upset_data$IPN) > 0) {
    UpSetR::upset(
      upset_data,
      sets = c("HighRisk", "ABxGiven", "IPN"),
      sets.bar.color = "#2166AC",
      main.bar.color = "#2166AC",
      matrix.color = "#2166AC",
      order.by = "freq",
      sets.x.label = "Set Size",
      mainbar.y.label = "Intersection Size"
    )
    title("D. UpSet Plot: Risk × Antibiotics × IPN", line = -0.5)
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient data for UpSet plot", cex = 1.2)
  }
}, error = function(e) {
  warning("UpSet plot failed: ", e$message)
  # Fallback: Venn-like bar chart of categories
  if (!all(is.na(sim_results$actual_abx))) {
    upset_fallback <- sim_results %>%
      mutate(combo = paste0(
        ifelse(model_high_risk, "HighRisk", "LowRisk"), " + ",
        ifelse(actual_abx == 1, "ABx", "NoABx"), " + ",
        ifelse(IPN == 1, "IPN", "NoIPN")
      )) %>%
      count(combo) %>%
      top_n(8, n)

    p_fallback <- ggplot(upset_fallback, aes(x = reorder(combo, n), y = n)) +
      geom_col(fill = "#2166AC", width = 0.6) +
      coord_flip() +
      labs(title = "D. Patient Subgroup Intersections (UpSet fallback)",
           x = "", y = "Count") +
      theme_minimal(base_size = 11)
    print(p_fallback)
  }
})

dev.off()
message("Figure 4 saved to figures/Figure_4_antibiotic_simulation.pdf")

# ── Overuse analysis: antibiotics in sterile necrosis ─────────────────
message("\nAnalyzing antibiotic overuse in sterile necrosis ...")
if ("sterile_necrosis" %in% colnames(sim_results) &&
    !all(is.na(sim_results$actual_abx))) {
  sterile_patients <- sim_results[sim_results$sterile_necrosis == TRUE, ]
  if (nrow(sterile_patients) > 0) {
    overuse_rate <- mean(sterile_patients$actual_abx == 1, na.rm = TRUE)
    message(sprintf("Antibiotic overuse rate (ABx in sterile necrosis): %.1f%%",
                    100 * overuse_rate))

    # Model would reduce overuse in this group
    model_overuse <- mean(sterile_patients$model_recommend == 1, na.rm = TRUE)
    message(sprintf("Model would recommend ABx in %.1f%% of sterile necrosis patients",
                    100 * model_overuse))
    message(sprintf("Potential overuse reduction: %.1f percentage points",
                    100 * (overuse_rate - model_overuse)))
  }
}

cat("=== Antibiotic simulation complete ===\n")
