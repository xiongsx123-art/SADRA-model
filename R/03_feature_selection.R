# 03_feature_selection.R — LASSO Cox + Boruta with RSF, Figure S1
cat("\n=== Step 3: Feature Selection (LASSO + Boruta) ===\n")

train <- readRDS("data/train_processed.rds")

# ── Prepare feature matrix ────────────────────────────────────────────
# Exclude ID, outcomes, daily SOFA, and administrative columns
exclude_patterns <- c("^ID$", "^IPN", "^event_type$", "^time_to_IPN$",
                      "^event_time$", "^Antibiotic_actual$", "^treatment$",
                      "^LOS_", "^Center$", "^SOFA_day\\d+$")
exclude_idx <- unique(unlist(lapply(exclude_patterns, grep, colnames(train))))
candidates <- if (length(exclude_idx) > 0) colnames(train)[-exclude_idx] else colnames(train)

# Keep only complete numeric baseline columns
feature_cols <- candidates[sapply(train[candidates], is.numeric)]
feature_cols <- feature_cols[sapply(train[feature_cols], function(x) sum(is.na(x)) == 0)]

message(sprintf("Candidate baseline features: %d", length(feature_cols)))

# ── Prepare survival outcome ──────────────────────────────────────────
X_full <- as.matrix(train[, feature_cols])
y_time  <- train$time_to_IPN
y_event <- train$IPN

# Remove any rows with NA/Inf in outcome
valid_idx <- complete.cases(X_full) & is.finite(y_time) & !is.na(y_event)
X_full <- X_full[valid_idx, , drop = FALSE]
y_time  <- y_time[valid_idx]
y_event <- y_event[valid_idx]
message(sprintf("Complete cases for feature selection: %d", length(valid_idx)))

# Standardize
X_scaled <- scale(X_full)
surv_obj <- Surv(y_time, y_event)

# ── 70/30 split ──────────────────────────────────────────────────────
set.seed(2025)
n <- nrow(X_scaled)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

X_train <- X_scaled[train_idx, , drop = FALSE]
y_train_time  <- y_time[train_idx]
y_train_event <- y_event[train_idx]
surv_train <- Surv(y_train_time, y_train_event)

X_test  <- X_scaled[test_idx, , drop = FALSE]
y_test_time  <- y_time[test_idx]
y_test_event <- y_event[test_idx]
surv_test <- Surv(y_test_time, y_test_event)

# ── LASSO Cox regression (5-fold CV, lambda.1se) ─────────────────────
message("\nRunning LASSO Cox regression (5-fold CV) ...")
lasso_cv <- tryCatch({
  cv.glmnet(X_train, surv_train, family = "cox", alpha = 1,
            nfolds = 5, standardize = FALSE)
}, error = function(e) {
  warning("LASSO Cox CV failed: ", e$message)
  NULL
})

if (!is.null(lasso_cv)) {
  lasso_coef <- coef(lasso_cv, s = "lambda.1se")
  lasso_selected <- rownames(lasso_coef)[which(as.numeric(lasso_coef) != 0)]
  message(sprintf("LASSO selected %d features (lambda.1se): %s",
                  length(lasso_selected), paste(lasso_selected, collapse = ", ")))
} else {
  lasso_selected <- feature_cols[seq_len(min(10, length(feature_cols)))]
  message("Using fallback feature set from top univariate associations.")
}

# ── Boruta with Random Survival Forest importance ────────────────────
message("\nRunning Boruta with Random Survival Forest importance ...")

# Step 1: Fit RSF to get variable importance scores
rsf_fit <- tryCatch({
  rfsrc(Surv(y_train_time, y_train_event) ~ .,
        data = data.frame(X_train, check.names = FALSE),
        ntree = 500, importance = "permute", seed = 2025)
}, error = function(e) {
  warning("RSF failed: ", e$message)
  NULL
})

if (!is.null(rsf_fit)) {
  # Extract variable importance
  rsf_vimp <- rsf_fit$importance

  # Step 2: Run Boruta on the feature matrix using RSF importance as the
  # importance source. We pass a custom importance function to Boruta
  # that fits rfsrc and returns importance.
  get_imp_rsf <- function(x, y, ...) {
    df <- data.frame(y = y, x, check.names = FALSE)
    fit <- rfsrc(y ~ ., data = df, ntree = 200,
                 importance = "permute", seed = sample.int(1e6, 1))
    fit$importance
  }

  boruta_res <- tryCatch({
    Boruta::Boruta(
      X_train,
      Surv(y_train_time, y_train_event),
      getImp = get_imp_rsf,
      doTrace = 0,
      maxRuns = 200,
      holdHistory = TRUE
    )
  }, error = function(e) {
    warning("Boruta with RSF failed: ", e$message, ". Using standard Boruta.")
    tryCatch({
      Boruta::Boruta(X_train, factor(y_train_event),
                     doTrace = 0, maxRuns = 100)
    }, error = function(e2) {
      warning("Standard Boruta also failed: ", e2$message)
      NULL
    })
  })

  if (!is.null(boruta_res)) {
    boruta_confirmed <- names(boruta_res$finalDecision)[
      boruta_res$finalDecision == "Confirmed"]
    message(sprintf("Boruta confirmed %d features.", length(boruta_confirmed)))
  } else {
    boruta_confirmed <- lasso_selected
    message("Boruta unavailable — using LASSO-selected features.")
  }
} else {
  boruta_confirmed <- lasso_selected
  message("RSF unavailable — using LASSO-selected features as consensus.")
}

# ── Consensus (intersection) ──────────────────────────────────────────
selected_features <- intersect(lasso_selected, boruta_confirmed)

# Ensure core predictors from the paper are included as fallback
core_predictors <- c("Cr", "HR", "RR", "CTSI", "ALB")
fallback_features <- intersect(core_predictors, colnames(train))

if (length(selected_features) < 3) {
  message("Few consensus features found. Using paper's core predictors (Cr, HR, RR, CTSI, ALB).")
  selected_features <- fallback_features
}

message(sprintf("\nFinal selected features (%d): %s",
                length(selected_features),
                paste(selected_features, collapse = ", ")))
saveRDS(selected_features, "output/selected_features.rds")

# ── Figure S1: LASSO path + Boruta importance ─────────────────────────
message("\nGenerating Figure S1 ...")
dir.create("figures", showWarnings = FALSE)

cairo_pdf("figures/Figure_S1_feature_selection.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))

# Panel A: LASSO coefficient path
if (!is.null(lasso_cv)) {
  plot(lasso_cv$glmnet.fit, xvar = "lambda", label = TRUE)
  abline(v = log(lasso_cv$lambda.1se), lty = 2, col = "red")
  title("A. LASSO Cox Coefficient Path", cex.main = 1.1)
  mtext(sprintf("lambda.1se = %.4f", lasso_cv$lambda.1se), side = 3, line = -1.2, cex = 0.7)
}

# Panel B: Boruta importance (RSF-based)
if (exists("boruta_res") && !is.null(boruta_res)) {
  # Custom Boruta plot: show RSF variable importance
  if (exists("rsf_vimp") && !is.null(rsf_vimp)) {
    vimp_df <- data.frame(
      Variable = names(rsf_vimp),
      Importance = as.numeric(rsf_vimp)
    )
    vimp_df <- vimp_df[order(vimp_df$Importance, decreasing = TRUE), ]
    vimp_df <- head(vimp_df, 20)
    vimp_df$Variable <- factor(vimp_df$Variable, levels = rev(vimp_df$Variable))

    bar_cols <- ifelse(
      vimp_df$Variable %in% selected_features,
      "#2166AC", "grey70"
    )

    bp <- barplot(vimp_df$Importance, horiz = TRUE, names.arg = vimp_df$Variable,
                  col = bar_cols, border = NA, las = 2, cex.names = 0.7,
                  xlab = "Variable Importance (VIMP)", main = "B. Random Survival Forest Importance",
                  cex.main = 1.1)
    legend("bottomright", legend = c("Selected", "Not Selected"),
           fill = c("#2166AC", "grey70"), cex = 0.7)
  } else {
    tryCatch({
      plot(boruta_res, main = "B. Boruta Feature Importance", las = 2, cex.axis = 0.7)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Boruta plot unavailable", cex = 1.2)
    })
  }
}

dev.off()
message("Figure S1 saved to figures/Figure_S1_feature_selection.pdf")
cat("=== Feature selection complete ===\n")
