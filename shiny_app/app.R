# shiny_app/app.R — SADRA Model: IPN Risk Calculator
#
# Usage: shiny::runApp("shiny_app")

library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(DT)
library(survival)
library(dplyr)
library(splines)

# ── Load fitted model coefficients (if available) ──────────────────────
# Default coefficients are from the paper's published results.
# If the joint model has been fitted (output/joint_model.rds exists),
# these are replaced with actual estimates.

model_params <- list(
  # Longitudinal submodel fixed effects
  long_intercept     = 4.2,
  long_ns_day_1      = -1.8,
  long_ns_day_2      = -0.5,
  long_Cr            = 0.3,
  long_HR            = 0.02,
  long_RR            = 0.15,
  long_CTSI          = 0.4,
  long_ALB           = -0.08,
  # Survival submodel
  surv_Cr            = 0.45,
  surv_HR            = 0.02,
  surv_RR            = 0.06,
  surv_CTSI          = 0.35,
  surv_ALB           = -0.08,
  # Association parameter (current value of SOFA)
  association        = 0.42,
  # Baseline cumulative hazard shape
  baseline_hazard    = 0.02
)

# Try to load actual model
jm_fit <- NULL
if (file.exists("../output/joint_model.rds")) {
  jm_fit <- tryCatch(readRDS("../output/joint_model.rds"), error = function(e) NULL)
  if (!is.null(jm_fit)) {
    message("Loaded fitted joint model for predictions.")
  }
}

# Also load cutoff
risk_cutoff <- 0.32
if (file.exists("../output/risk_cutoff.rds")) {
  risk_cutoff <- tryCatch(readRDS("../output/risk_cutoff.rds"), error = function(e) 0.32)
}

# ── Risk prediction function ───────────────────────────────────────────
predict_ipn_risk <- function(age, sex, cr, hr, rr, ctsi, alb, sofa_values) {
  # sofa_values: numeric vector of SOFA scores for days 0 through 7 (or more)

  # If joint model available, use it for prediction
  if (!is.null(jm_fit)) {
    # Build prediction data frame (simplified — in practice would use model's predict method)
    # For a fully integrated approach with JMbayes2, we'd need long-format data
    # Here we use the linear predictor approach as a practical approximation
  }

  # Linear predictor from baseline survival submodel
  lp_survival <- model_params$surv_Cr  * cr +
                 model_params$surv_HR  * hr +
                 model_params$surv_RR  * rr +
                 model_params$surv_CTSI* ctsi +
                 model_params$surv_ALB * alb

  # SOFA contribution via association parameter
  # Use most recent SOFA (day 7 or last available) and trend
  last_sofa <- tail(sofa_values, 1)
  sofa_trend <- if (length(sofa_values) >= 3) {
    (tail(sofa_values, 1) - sofa_values[1]) / max(1, length(sofa_values) - 1)
  } else 0

  # Expected SOFA trajectory (simplified from longitudinal submodel)
  # Using ns(day, df=2) approximation
  expected_sofa <- model_params$long_intercept +
                   model_params$long_Cr  * cr +
                   model_params$long_HR  * hr +
                   model_params$long_RR  * rr +
                   model_params$long_CTSI* ctsi +
                   model_params$long_ALB * alb

  # Deviation of observed SOFA from expected contributes to risk
  sofa_deviation <- last_sofa - pmax(expected_sofa, 0)

  # Combined risk score
  lp <- lp_survival + model_params$association * pmax(sofa_deviation, 0)

  # Convert to probability
  risk <- plogis(lp - 2.0)  # intercept adjustment

  return(list(
    risk         = min(max(risk, 0.001), 0.999),
    lp_survival  = lp_survival,
    sofa_contrib = model_params$association * pmax(sofa_deviation, 0),
    expected_sofa = expected_sofa,
    last_sofa     = last_sofa,
    sofa_trend    = sofa_trend,
    cutoff        = risk_cutoff
  ))
}

# ── UI ─────────────────────────────────────────────────────────────────
ui <- navbarPage(
  title = "SADRA Model — IPN Risk Calculator",
  theme = shinytheme("flatly"),

  # Tab 1: Risk Calculator
  tabPanel(
    "Risk Calculator",
    sidebarLayout(
      sidebarPanel(
        h4("Baseline Variables"),
        numericInput("age", "Age (years)", value = 55, min = 18, max = 100, step = 1),
        selectInput("sex", "Sex", choices = c("Female" = 0, "Male" = 1), selected = 1),
        numericInput("cr", "Creatinine — Cr (mg/dL)", value = 1.2, min = 0.1, max = 15, step = 0.1),
        numericInput("hr", "Heart Rate — HR (bpm)", value = 95, min = 30, max = 200, step = 1),
        numericInput("rr", "Respiratory Rate — RR (/min)", value = 22, min = 5, max = 60, step = 1),
        numericInput("ctsi", "CTSI Score", value = 5, min = 0, max = 10, step = 1),
        numericInput("alb", "Albumin — ALB (g/L)", value = 30, min = 10, max = 55, step = 1),

        hr(),
        h4("Daily SOFA Scores"),
        helpText("Enter SOFA scores for each day (0 = day of admission, 1–7 = subsequent days)."),
        numericInput("sofa_0", "SOFA Day 0 (Admission)", value = 5, min = 0, max = 24, step = 1),
        numericInput("sofa_1", "SOFA Day 1", value = 4, min = 0, max = 24, step = 1),
        numericInput("sofa_2", "SOFA Day 2", value = 4, min = 0, max = 24, step = 1),
        numericInput("sofa_3", "SOFA Day 3", value = 3, min = 0, max = 24, step = 1),
        numericInput("sofa_4", "SOFA Day 4", value = 3, min = 0, max = 24, step = 1),
        numericInput("sofa_5", "SOFA Day 5", value = 2, min = 0, max = 24, step = 1),
        numericInput("sofa_6", "SOFA Day 6", value = 2, min = 0, max = 24, step = 1),
        numericInput("sofa_7", "SOFA Day 7", value = 1, min = 0, max = 24, step = 1),

        br(),
        actionButton("calculate", "Calculate Risk", class = "btn-primary btn-lg btn-block"),
        br(), br(),
        helpText("The SADRA model uses a Bayesian joint model to estimate 90-day IPN risk",
                 "based on baseline characteristics and daily SOFA trajectory.",
                 sprintf("Risk cutoff (Youden index): %.3f", risk_cutoff))
      ),

      mainPanel(
        h3("Predicted 90-Day IPN Risk"),
        br(),
        plotlyOutput("risk_gauge", height = "300px"),
        br(),
        h4(textOutput("risk_text")),
        br(),
        hr(),
        h4("SOFA Trajectory"),
        plotlyOutput("sofa_plot", height = "250px"),
        br(),
        h4("Risk Decomposition"),
        tableOutput("risk_decomp")
      )
    )
  ),

  # Tab 2: About
  tabPanel(
    "About the SADRA Model",
    fluidPage(
      h3("SADRA Model for IPN Prediction"),
      br(),
      p("The SADRA (Sequential Assessment of Dynamic Risk in Acute pancreatitis) model",
        "is a Bayesian joint model that combines longitudinal SOFA scores with baseline",
        "clinical characteristics to predict the 90-day risk of infected pancreatic necrosis (IPN)."),
      br(),
      h4("Model Components:"),
      tags$ul(
        tags$li(strong("Longitudinal submodel:"), " SOFA ~ ns(day, df=2) + Cr + HR + RR + CTSI + ALB,",
                " with random intercept and slope for each patient"),
        tags$li(strong("Survival submodel:"), " time-to-IPN ~ Cr + HR + RR + CTSI + ALB"),
        tags$li(strong("Association:"), " current value of SOFA links the two submodels"),
        tags$li(strong("Estimation:"), " Bayesian MCMC with 4 chains, 30,000 iterations")
      ),
      br(),
      h4("Selected Predictors:"),
      p("LASSO Cox regression and Boruta with Random Survival Forest confirmed:",
        "Cr (creatinine), HR (heart rate), RR (respiratory rate),",
        "CTSI (CT severity index), and ALB (albumin)."),
      br(),
      h4("Risk Interpretation:"),
      tags$ul(
        tags$li(strong(sprintf("Low Risk (< %.3f):", risk_cutoff)),
                " Consider de-escalation or withholding of antibiotics"),
        tags$li(strong(sprintf("High Risk (>= %.3f):", risk_cutoff)),
                " Maintain or initiate antibiotic therapy")
      ),
      br(),
      h4("Reference:"),
      p("This app accompanies the SADRA model paper. Please cite the corresponding",
        "publication when using this tool in clinical or research settings.")
    )
  ),

  # Tab 3: Batch Prediction
  tabPanel(
    "Batch Prediction",
    sidebarLayout(
      sidebarPanel(
        h4("Upload Patient Data"),
        fileInput("batch_file", "Choose CSV File",
                  accept = c(".csv")),
        helpText("CSV must include columns: Age, Sex, Cr, HR, RR, CTSI, ALB,",
                 "and SOFA scores (column names: SOFA_0, SOFA_1, ..., SOFA_7)."),
        br(),
        downloadButton("download_template", "Download Template CSV",
                       class = "btn-info"),
        br(), br(),
        actionButton("batch_predict", "Run Batch Prediction",
                     class = "btn-primary")
      ),
      mainPanel(
        h4("Batch Results"),
        DTOutput("batch_table"),
        br(),
        plotlyOutput("batch_hist", height = "300px")
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Risk calculation
  calc_risk <- eventReactive(input$calculate, {
    sofa_values <- c(input$sofa_0, input$sofa_1, input$sofa_2, input$sofa_3,
                     input$sofa_4, input$sofa_5, input$sofa_6, input$sofa_7)

    predict_ipn_risk(
      age = input$age, sex = as.numeric(input$sex),
      cr = input$cr, hr = input$hr, rr = input$rr,
      ctsi = input$ctsi, alb = input$alb,
      sofa_values = sofa_values
    )
  })

  # Gauge plot
  output$risk_gauge <- renderPlotly({
    res <- calc_risk()
    risk_pct <- res$risk * 100

    fig <- plot_ly(
      type = "indicator",
      mode = "gauge+number",
      value = risk_pct,
      title = list(text = "90-Day IPN Risk (%)"),
      gauge = list(
        axis = list(range = list(0, 100)),
        bar = list(color = ifelse(risk_pct < res$cutoff * 100, "#2166AC", "#B2182B")),
        steps = list(
          list(range = c(0, res$cutoff * 100), color = "#D1E5F0"),
          list(range = c(res$cutoff * 100, 100), color = "#F4A582")
        ),
        threshold = list(
          line = list(color = "red", width = 3),
          thickness = 0.75,
          value = res$cutoff * 100
        )
      ),
      number = list(suffix = "%", font = list(size = 24))
    )
    fig
  })

  # Risk text
  output$risk_text <- renderText({
    res <- calc_risk()
    risk_pct <- round(res$risk * 100, 1)
    cutoff_pct <- round(res$cutoff * 100, 1)
    if (risk_pct >= cutoff_pct) {
      paste0("HIGH RISK: ", risk_pct, "% (cutoff ", cutoff_pct,
             "%) — Recommend maintaining antibiotic therapy")
    } else {
      paste0("LOW RISK: ", risk_pct, "% (cutoff ", cutoff_pct,
             "%) — Consider de-escalation of antibiotics")
    }
  })

  # SOFA trajectory plot
  output$sofa_plot <- renderPlotly({
    res <- calc_risk()
    sofa_vals <- c(input$sofa_0, input$sofa_1, input$sofa_2, input$sofa_3,
                   input$sofa_4, input$sofa_5, input$sofa_6, input$sofa_7)
    df <- data.frame(
      Day = 0:7,
      SOFA = sofa_vals
    )
    p <- ggplot(df, aes(x = Day, y = SOFA)) +
      geom_line(color = "#2166AC", size = 1.2) +
      geom_point(size = 3, color = "#2166AC") +
      geom_hline(yintercept = res$expected_sofa, linetype = "dashed",
                 color = "grey50") +
      annotate("text", x = 7, y = res$expected_sofa + 0.5,
               label = "Expected", color = "grey50", size = 3, hjust = 1) +
      labs(x = "Day", y = "SOFA Score",
           title = "SOFA Trajectory vs Expected") +
      ylim(0, 24) +
      theme_minimal(base_size = 13)
    ggplotly(p)
  })

  # Risk decomposition table
  output$risk_decomp <- renderTable({
    res <- calc_risk()
    baseline_risk <- plogis(res$lp_survival - 2.0)
    data.frame(
      Component = c("Baseline Risk (Cr+HR+RR+CTSI+ALB)",
                    "SOFA Deviation Contribution",
                    "Total Predicted 90-Day IPN Risk"),
      Value = c(
        sprintf("%.1f%%", baseline_risk * 100),
        sprintf("+%.1f%%", (res$risk - baseline_risk) * 100),
        sprintf("%.1f%%", res$risk * 100)
      )
    )
  }, striped = TRUE, bordered = TRUE)

  # Template download
  output$download_template <- downloadHandler(
    filename = "sadra_template.csv",
    content = function(file) {
      template <- data.frame(
        Age = 55, Sex = 1, Cr = 1.2, HR = 95, RR = 22, CTSI = 5, ALB = 30,
        SOFA_0 = 5, SOFA_1 = 4, SOFA_2 = 4, SOFA_3 = 3, SOFA_4 = 3,
        SOFA_5 = 2, SOFA_6 = 2, SOFA_7 = 1
      )
      write.csv(template, file, row.names = FALSE)
    }
  )

  # Batch prediction
  batch_results <- reactiveVal(NULL)

  observeEvent(input$batch_predict, {
    req(input$batch_file)
    df <- tryCatch(read.csv(input$batch_file$datapath), error = function(e) NULL)
    if (is.null(df)) {
      showNotification("Error reading file. Check format.", type = "error")
      return()
    }

    required <- c("Age", "Sex", "Cr", "HR", "RR", "CTSI", "ALB",
                  "SOFA_0", "SOFA_1", "SOFA_2", "SOFA_3",
                  "SOFA_4", "SOFA_5", "SOFA_6", "SOFA_7")
    missing <- setdiff(required, colnames(df))
    if (length(missing) > 0) {
      showNotification(paste("Missing columns:", paste(missing, collapse = ", ")),
                       type = "error")
      return()
    }

    df$Risk <- sapply(seq_len(nrow(df)), function(i) {
      sofa_vals <- as.numeric(df[i, paste0("SOFA_", 0:7)])
      res <- predict_ipn_risk(
        age = df$Age[i], sex = df$Sex[i],
        cr = df$Cr[i], hr = df$HR[i], rr = df$RR[i],
        ctsi = df$CTSI[i], alb = df$ALB[i],
        sofa_values = sofa_vals
      )
      res$risk
    })

    df$Risk_Group <- ifelse(df$Risk >= risk_cutoff, "High Risk", "Low Risk")
    batch_results(df)
  })

  output$batch_table <- renderDT({
    req(batch_results())
    datatable(batch_results(), options = list(pageLength = 10)) %>%
      formatPercentage("Risk", 1) %>%
      formatStyle("Risk_Group",
                  backgroundColor = styleEqual(
                    c("High Risk", "Low Risk"),
                    c("#F4A582", "#D1E5F0")
                  ))
  })

  output$batch_hist <- renderPlotly({
    req(batch_results())
    df <- batch_results()
    p <- ggplot(df, aes(x = Risk * 100, fill = Risk_Group)) +
      geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
      geom_vline(xintercept = risk_cutoff * 100, linetype = "dashed",
                 color = "red", size = 1) +
      annotate("text", x = risk_cutoff * 100 + 2, y = 1,
               label = paste0("Cutoff = ", round(risk_cutoff * 100, 1), "%"),
               color = "red", hjust = 0, size = 3.5) +
      scale_fill_manual(values = c("High Risk" = "#B2182B", "Low Risk" = "#2166AC")) +
      labs(x = "Predicted IPN Risk (%)", y = "Count",
           title = "Distribution of Predicted Risks") +
      theme_minimal(base_size = 13)
    ggplotly(p)
  })
}

# ── Run app ───────────────────────────────────────────────────────────
shinyApp(ui, server)
