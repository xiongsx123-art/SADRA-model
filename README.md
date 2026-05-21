# SADRA Model: Sequential Assessment of Dynamic Risk in Acute Pancreatitis

Complete R analysis pipeline reproducing the SADRA model paper — a Bayesian joint model framework for predicting infected pancreatic necrosis (IPN) risk in acute pancreatitis patients.

## Overview

The SADRA model integrates longitudinal SOFA scores with baseline clinical characteristics (Cr, HR, RR, CTSI, ALB) using a Bayesian joint model (JMbayes2) to dynamically predict 90-day IPN risk. The pipeline includes:

- **Feature selection**: LASSO Cox regression + Boruta with Random Survival Forest
- **Joint model**: Bayesian joint model with current-value association, natural splines, explicit priors
- **Performance evaluation**: Time-dependent AUC with cluster bootstrap, calibration, Kaplan-Meier curves
- **Sensitivity analyses**: Static ML comparison, alternative association structures, alternative severity scores (SIRS/APACHE-II/BISAP), strict IPN definition, external validation (TRACE cohort)
- **Antibiotic simulation**: Model-guided antibiotic strategy with Sankey diagrams and UpSet plots
- **Treatment heterogeneity**: Double Machine Learning (DoubleML) for ATE/CATE estimation in the TRACE trial
- **Shiny app**: Interactive IPN risk calculator
- **R Markdown report**: Full analysis report

## Quick Start

### 1. Place Data Files
Copy your CSV files into the `data/` folder:
```
data/
  training_cohort.csv      (retrospective, n=812)
  prospective_cohort.csv   (prospective, n=465)
  trace_cohort.csv         (TRACE trial, n=437)
```

### 2. Set GitHub Token (Optional)
```r
Sys.setenv(GITHUB_PAT = "your_personal_access_token")
```

### 3. Run Pipeline
```bash
Rscript run_all.R
```

Or from RStudio:
```r
source("run_all.R")
```

### For faster testing:
```bash
# Windows
set SADRA_TESTING=TRUE && Rscript run_all.R

# Linux/Mac
SADRA_TESTING=TRUE Rscript run_all.R
```

## Required Data Format

### Required Columns (exact spelling)

**Baseline predictors:**
- `ID`, `Age`, `Sex`, `Cr`, `HR`, `RR`, `CTSI`, `ALB`

**Daily SOFA scores:**
- `SOFA_day1`, `SOFA_day2`, ..., `SOFA_day7` (minimum; extended range up to SOFA_day14 supported)

**Outcome:**
- `IPN` (composite, 0/1)
- `time_to_IPN` (days, max 90)

**TRACE cohort additional:**
- `treatment` (thymosin_alpha1, 0/1)

**Optional (for sensitivity analyses):**
- `IPN_strict`, `Antibiotic_actual`
- `SIRS_day1` through `SIRS_day7`
- `APACHEII_day1` through `APACHEII_day7`
- `BISAP_day1` through `BISAP_day7`

## Output Structure

```
output/
  selected_features.rds        Selected predictor set
  joint_model.rds              Fitted Bayesian joint model
  AUC_results.csv              Time-dependent AUC values
  risk_cutoff.rds              Optimal Youden cutoff
  antibiotic_simulation.csv    Simulation results
  TRACE_validation_AUC.csv     External validation results
  Table_S2_*.csv/.txt          Longitudinal submodel parameters
  Table_S3_*.txt               Survival submodel parameters
  Table_S4_*.txt               Association parameters
  Table_S5_*.csv/.txt          ML comparison / Alt association
  Table_S6_*.csv               CATE by risk group

figures/
  Figure_2_AUC_curve.pdf       Time-dependent AUC curve
  Figure_2_calibration.pdf     Calibration plots
  Figure_3_KM_curves.pdf       KM curves by risk group
  Figure_4_antibiotic_simulation.pdf  Sankey + bar + UpSet
  Figure_5_CATE_forest.pdf     Forest plot + propensity overlap
  Figure_S1_feature_selection.pdf     LASSO path + Boruta importance
  Figure_S3_ML_comparison.pdf         Static ML comparison
  Figure_S4-S6_*.pdf                  Alternative severity scores
  Figure_S7_strict_IPN.pdf            Strict IPN definition
  Figure_S8_TRACE_validation.pdf      External validation
```

## Reproducibility

- Set seed: `set.seed(2025)` for all analyses
- MCMC: 4 chains, 30,000 iterations, 15,000 burn-in, thin = 4
- MICE: m = 5 imputations
- Bootstrap: 1,000 cluster resamples (by patient ID)
- Cross-fitting: 5-fold, repeated 5 times (Double ML)

## Methods Summary

| Step | Method | Package |
|------|--------|---------|
| Missing data | MICE (m=5) for baseline; NA retained for SOFA | `mice` |
| Feature selection | LASSO Cox + Boruta with RSF importance | `glmnet`, `Boruta`, `randomForestSRC` |
| Joint model | Bayesian joint model, current-value association | `JMbayes2`, `nlme`, `survival` |
| Performance | Time-dependent AUC, cluster bootstrap | `pROC`, `timeROC` |
| Calibration | Decile-based calibration plots | `ggplot2` |
| Sensitivity | Static ML (Cox, CoxBoost, LASSO-Cox, RSF) | `glmnet`, `CoxBoost`, `randomForestSRC` |
| Antibiotic simulation | Sankey diagram, UpSet plot | `ggalluvial`, `UpSetR` |
| Treatment heterogeneity | Double ML with random forest | `DoubleML`, `mlr3`, `ranger` |

## Shiny App

```r
shiny::runApp("shiny_app")
```

The Shiny app provides:
- Individual IPN risk calculation with gauge plot
- Dynamic SOFA trajectory visualization
- Batch prediction with CSV upload
- Risk decomposition showing baseline + SOFA contributions

## Report

Knit `reports/report.Rmd` to generate an HTML report with all figures and tables.

```r
rmarkdown::render("reports/report.Rmd")
```

## Citation

This code accompanies the SADRA model paper. Please cite the corresponding publication when using this code or the SADRA model.
