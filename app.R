# ============================================================
#  ELISA Kd Calculator — R Shiny App
#  Model: One-site specific binding  Y = Bmax * X / (Kd + X)
#  Plate layout:
#    Rows A/B  = Novel binding protein replicate 1
#    Rows C/D  = Novel binding protein replicate 2
#    Rows E/F  = Positive control replicate 1
#    Rows G/H  = Positive control replicate 2
#  Odd rows = signal, Even rows = blank (subtracted)
# ============================================================

library(shiny)
library(rhandsontable)
library(ggplot2)
library(dplyr)
library(tidyr)
library(minpack.lm)   # robust NLS (Levenberg-Marquardt)

# ── Default absorbance data — rows = A-H, columns = 1-12 (matches machine output) ──
default_abs <- data.frame(
  Row = c("A","B","C","D","E","F","G","H"),
  `1`  = c(0.8096, 0.0798, 0.6853, 0.0642, 0.7711, 0.1247, 0.6712, 0.1017),
  `2`  = c(0.7588, 0.0612, 0.7489, 0.0602, 0.7351, 0.0702, 0.7392, 0.0773),
  `3`  = c(0.7612, 0.0620, 0.5425, 0.0602, 0.6753, 0.0686, 0.7538, 0.0664),
  `4`  = c(0.7627, 0.0584, 0.4365, 0.0611, 0.7100, 0.0659, 0.7906, 0.0660),
  `5`  = c(0.7485, 0.0587, 0.2819, 0.0590, 0.6063, 0.0640, 0.7094, 0.0625),
  `6`  = c(0.6667, 0.0586, 0.2403, 0.0612, 0.6451, 0.0617, 0.6275, 0.0618),
  `7`  = c(0.5548, 0.0591, 0.1830, 0.0615, 0.4808, 0.0634, 0.4594, 0.0687),
  `8`  = c(0.5002, 0.0602, 0.1156, 0.0643, 0.4526, 0.0639, 0.4128, 0.0656),
  `9`  = c(0.3802, 0.0611, 0.1291, 0.0627, 0.4384, 0.0678, 0.3629, 0.0669),
  `10` = c(0.3091, 0.0636, 0.0990, 0.0622, 0.2504, 0.0654, 0.2856, 0.0728),
  `11` = c(0.2629, 0.0663, 0.1011, 0.0645, 0.1579, 0.0648, 0.2054, 0.0673),
  `12` = c(0.1536, 0.0624, 0.0793, 0.0622, 0.1213, 0.0709, 0.1271, 0.0658),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ── One-site specific binding ─────────────────────────────────────────────────
one_site <- function(x, Bmax, Kd) (Bmax * x) / (Kd + x)

fit_kd <- function(conc, signal) {
  # Remove non-positive concentrations / NA
  keep <- is.finite(conc) & is.finite(signal) & conc > 0
  x <- conc[keep]; y <- signal[keep]
  if (length(x) < 4) return(NULL)
  
  Bmax0 <- max(y, na.rm = TRUE)
  Kd0   <- median(x)
  
  tryCatch(
    nlsLM(y ~ one_site(x, Bmax, Kd),
          start = list(Bmax = Bmax0, Kd = Kd0),
          lower = c(Bmax = 0, Kd = 0),
          control = nls.lm.control(maxiter = 200)),
    error = function(e) NULL
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', Arial, sans-serif; background: #f5f7fa; }
    .well { background: #ffffff; border: 1px solid #dde3ea; border-radius: 8px; box-shadow: none; }
    h2 { color: #1a3a5c; font-weight: 700; margin-bottom: 4px; }
    h4 { color: #2c5f8a; font-weight: 600; }
    .section-title { color: #1a3a5c; font-weight: 600; font-size: 15px;
                     border-bottom: 2px solid #2c8adf; padding-bottom: 4px; margin: 18px 0 10px; }
    .result-box { background: #eaf4ff; border-left: 4px solid #2c8adf; border-radius: 6px;
                  padding: 12px 16px; margin: 8px 0; }
    .result-box.positive { background: #eaffea; border-left-color: #2ca02c; }
    .result-label { font-size: 12px; color: #555; margin: 0; }
    .result-value { font-size: 22px; font-weight: 700; color: #1a3a5c; margin: 0; }
    .result-se    { font-size: 12px; color: #666; }
    .hot { margin-bottom: 6px; }
    .btn-primary { background-color: #2c8adf; border-color: #1a6ab0; }
    .btn-primary:hover { background-color: #1a6ab0; }
    .note { font-size: 12px; color: #777; }
  "))),
  
  titlePanel(
    div(h2("ELISA Dissociation Constant (K\u1d48) Calculator"),
        p("One-site specific binding model  |  Y = B\u2098\u2090\u02e3 \u00b7 [X] / (K\u1d48 + [X])",
          style = "color:#555; font-size:13px; margin:0;"))
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 wellPanel(
                   div(class = "section-title", "\u2699\ufe0f  Dilution Series"),
                   numericInput("start_conc", "Highest concentration (µM):", value = 10, min = 0, step = 0.1),
                   numericInput("fold",       "Fold dilution:",               value = 2,  min = 1.01, step = 0.5),
                   p(class = "note", "Column 1 = highest, Column 12 = lowest."),
                   hr(),
                   div(class = "section-title", "\ud83d\udcc8  Display"),
                   checkboxInput("show_rep",   "Show individual replicates", value = TRUE),
                   checkboxInput("show_resid", "Show residual plot",         value = FALSE),
                   hr(),
                   actionButton("calc", "Calculate K\u1d48", class = "btn-primary btn-block",
                                icon = icon("calculator"))
                 )
    ),
    
    mainPanel(width = 9,
              # ── Concentration preview ──────────────────────────────────────────────
              wellPanel(
                div(class = "section-title", "Concentrations (µM) — auto-generated from dilution settings"),
                verbatimTextOutput("conc_preview")
              ),
              
              # ── Data entry table ───────────────────────────────────────────────────
              wellPanel(
                div(class = "section-title", "Raw Absorbance Data"),
                p(icon("paste"),
                  tags$strong("You can copy and paste your data directly from the plate reader software into the table below."),
                  " Select the 8 x 12 block of values from your spreadsheet, click the top-left data cell (row A, column 1), and paste.",
                  style = "font-size:13px; color:#444; margin-bottom:8px;"),
                p(class = "note", style = "margin-bottom:8px;",
                  "Layout: rows A, C, E, G = signal | rows B, D, F, H = blank. The 'Row' column is read-only."),
                rHandsontableOutput("hot_table"),
                br(),
                p(class = "note",
                  icon("info-circle"),
                  " Blank subtraction is applied automatically: Net signal = signal row minus paired blank row.")
              ),
              
              # ── Results ───────────────────────────────────────────────────────────
              uiOutput("results_panel"),
              
              # ── Plots ─────────────────────────────────────────────────────────────
              uiOutput("plot_panel")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Reactive: concentration vector ────────────────────────────────────────
  conc_vec <- reactive({
    req(input$start_conc, input$fold)
    input$start_conc / (input$fold ^ (0:11))
  })
  
  output$conc_preview <- renderText({
    cv <- round(conc_vec(), 6)
    paste(sprintf("Col%2d: %s µM", 1:12, formatC(cv, format = "g", digits = 4)), collapse = "   ")
  })
  
  # ── Reactive: rhandsontable ────────────────────────────────────────────────
  rv <- reactiveVal(default_abs)
  
  output$hot_table <- renderRHandsontable({
    rhandsontable(rv(),
                  rowHeaders = NULL,
                  stretchH   = "all",
                  height     = 240) %>%
      hot_col("Row", readOnly = TRUE, width = 40) %>%
      hot_cols(format = "0.0000") %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })
  
  observe({
    if (!is.null(input$hot_table))
      rv(hot_to_r(input$hot_table))
  })
  
  # ── Reactive: processed data (blank-subtracted, net signals) ──────────────
  processed <- eventReactive(input$calc, {
    df   <- rv()
    conc <- conc_vec()
    
    # Helper: extract a named row as numeric vector across cols 1-12
    get_row <- function(row_label) {
      as.numeric(df[df$Row == row_label, as.character(1:12)])
    }
    
    # Net = signal row - blank row
    novel_net <- rbind(
      data.frame(rep = "Rep 1", conc = conc, net = get_row("A") - get_row("B")),
      data.frame(rep = "Rep 2", conc = conc, net = get_row("C") - get_row("D"))
    )
    ctrl_net <- rbind(
      data.frame(rep = "Rep 1", conc = conc, net = get_row("E") - get_row("F")),
      data.frame(rep = "Rep 2", conc = conc, net = get_row("G") - get_row("H"))
    )
    # Clamp negatives to 0 (below-blank artefacts)
    novel_net$net <- pmax(novel_net$net, 0)
    ctrl_net$net  <- pmax(ctrl_net$net,  0)
    
    # Mean per concentration
    novel_mean <- novel_net %>% group_by(conc) %>% summarise(mean_net = mean(net), .groups="drop")
    ctrl_mean  <- ctrl_net  %>% group_by(conc) %>% summarise(mean_net = mean(net), .groups="drop")
    
    list(novel_rep  = novel_net,
         ctrl_rep   = ctrl_net,
         novel_mean = novel_mean,
         ctrl_mean  = ctrl_mean)
  })
  
  # ── Reactive: fits ─────────────────────────────────────────────────────────
  fits <- eventReactive(input$calc, {
    p <- processed()
    
    # CHANGED: Fit Kd using individual replicates rather than the mean data
    list(
      novel = fit_kd(p$novel_rep$conc, p$novel_rep$net),
      ctrl  = fit_kd(p$ctrl_rep$conc,  p$ctrl_rep$net)
    )
  })
  
  # ── Results panel ─────────────────────────────────────────────────────────
  output$results_panel <- renderUI({
    req(fits())
    f <- fits()
    
    make_box <- function(fit, label, cls) {
      if (is.null(fit)) {
        div(class = paste("result-box", cls),
            p(class = "result-label", label),
            p(class = "result-value", "Fit failed"),
            p(class = "result-se", "Check data — ensure sufficient signal above blank"))
      } else {
        co   <- summary(fit)$coefficients
        Kd   <- co["Kd",   "Estimate"]
        KdSE <- co["Kd",   "Std. Error"]
        Bmax <- co["Bmax", "Estimate"]
        BmSE <- co["Bmax", "Std. Error"]
        r2   <- 1 - sum(residuals(fit)^2) / sum((fitted(fit) + residuals(fit) - mean(fitted(fit) + residuals(fit)))^2)
        div(class = paste("result-box", cls),
            p(class = "result-label", label),
            p(class = "result-value",
              sprintf("K\u1d48 = %.4f µM", Kd)),
            p(class = "result-se",
              sprintf("\u00b1 %.4f µM (SE)   |   B\u2098\u2090\u02e3 = %.4f \u00b1 %.4f   |   R\u00b2 = %.4f",
                      KdSE, Bmax, BmSE, r2))
        )
      }
    }
    
    wellPanel(
      div(class = "section-title", "\ud83d\udcca Results"),
      fluidRow(
        column(6, make_box(f$novel, "Novel Binding Protein", "")),
        column(6, make_box(f$ctrl,  "Positive Control",      "positive"))
      )
    )
  })
  
  # ── Plot panel ────────────────────────────────────────────────────────────
  output$plot_panel <- renderUI({
    req(fits())
    plots <- list(plotOutput("binding_plot", height = "380px"))
    if (input$show_resid)
      plots <- c(plots, list(plotOutput("resid_plot", height = "250px")))
    wellPanel(div(class = "section-title", "\ud83d\udcc9 Binding Curves"), do.call(tagList, plots))
  })
  
  output$binding_plot <- renderPlot({
    req(processed(), fits())
    p <- processed()
    f <- fits()
    
    x_seq <- seq(0, max(conc_vec()) * 1.05, length.out = 400)
    
    pred_df <- bind_rows(
      if (!is.null(f$novel)) data.frame(x = x_seq, y = predict(f$novel, newdata = list(x = x_seq)), Sample = "Novel Binding Protein") else NULL,
      if (!is.null(f$ctrl))  data.frame(x = x_seq, y = predict(f$ctrl,  newdata = list(x = x_seq)), Sample = "Positive Control") else NULL
    )
    
    mean_df <- bind_rows(
      data.frame(p$novel_mean, Sample = "Novel Binding Protein"),
      data.frame(p$ctrl_mean,  Sample = "Positive Control")
    )
    
    # CHANGED: Added facet_wrap and removed legend
    g <- ggplot(mean_df, aes(x = conc, y = mean_net, colour = Sample, fill = Sample)) +
      geom_point(size = 3.5, shape = 16) +
      labs(x = "Concentration (µM)", y = "Net Absorbance (blank-subtracted)",
           title = "One-Site Specific Binding Fit") +
      scale_colour_manual(values = c("Novel Binding Protein" = "#2c8adf",
                                     "Positive Control"      = "#2ca02c")) +
      scale_fill_manual(  values = c("Novel Binding Protein" = "#2c8adf44",
                                     "Positive Control"      = "#2ca02c44")) +
      facet_wrap(~ Sample) + 
      theme_bw(base_size = 13) +
      theme(legend.position = "none",
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", colour = "#1a3a5c"))
    
    if (input$show_rep) {
      rep_df <- bind_rows(
        data.frame(p$novel_rep, Sample = "Novel Binding Protein"),
        data.frame(p$ctrl_rep,  Sample = "Positive Control")
      )
      g <- g + geom_point(data = rep_df, aes(x = conc, y = net),
                          shape = 1, size = 2.5, alpha = 0.55, stroke = 1)
    }
    
    if (nrow(pred_df) > 0)
      g <- g + geom_line(data = pred_df, aes(x = x, y = y), linewidth = 1.1)
    
    # CHANGED: Use a dataframe to add Kd lines securely within respective facets
    kd_df <- bind_rows(
      if (!is.null(f$novel)) data.frame(Sample = "Novel Binding Protein", kd = coef(f$novel)["Kd"]) else NULL,
      if (!is.null(f$ctrl))  data.frame(Sample = "Positive Control",      kd = coef(f$ctrl)["Kd"]) else NULL
    )
    
    if (nrow(kd_df) > 0) {
      y_max <- max(mean_df$mean_net, na.rm = TRUE)
      y_pos <- y_max * 0.95 
      
      g <- g +
        geom_vline(data = kd_df, aes(xintercept = kd, colour = Sample), 
                   linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
        geom_label(data = kd_df, aes(x = kd, y = y_pos, label = sprintf("Kd = %.3f µM", kd), colour = Sample),
                   fill = "white", label.size = 0.3, hjust = -0.08, size = 4.5, fontface = "bold", show.legend = FALSE)
    }
    
    g
  })
  
  output$resid_plot <- renderPlot({
    req(processed(), fits())
    p <- processed(); f <- fits()
    
    resid_df <- bind_rows(
      if (!is.null(f$novel)) data.frame(
        fitted = fitted(f$novel), resid = residuals(f$novel), Sample = "Novel Binding Protein"
      ) else NULL,
      if (!is.null(f$ctrl)) data.frame(
        fitted = fitted(f$ctrl), resid = residuals(f$ctrl), Sample = "Positive Control"
      ) else NULL
    )
    
    # CHANGED: Added facet_wrap and removed legend
    ggplot(resid_df, aes(x = fitted, y = resid, colour = Sample)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_point(size = 3) +
      facet_wrap(~ Sample) +
      scale_colour_manual(values = c("Novel Binding Protein" = "#2c8adf",
                                     "Positive Control"      = "#2ca02c")) +
      labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", colour = "#1a3a5c"))
  })
}

shinyApp(ui, server)