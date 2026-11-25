# ============================================================
#  ELISA Kd Calculator - Guided Workflow (R Shiny)
#  Model: One-site specific binding  Y = Bmax * X / (Kd + X)
#
#  Workflow:
#    Step 1  Background   - what a Kd is + interactive explorer
#    Step 2  Enter data   - paste plate data, view raw QC plot
#    Step 3  Interpret    - prompted questions (gate: all must be answered)
#    Step 4  Results      - fit, Kd, comparison to prediction, reflection
#
#  Plate layout:
#    Rows A/B = Novel binding protein replicate 1 (A = signal, B = blank)
#    Rows C/D = Novel binding protein replicate 2
#    Rows E/F = Positive control replicate 1
#    Rows G/H = Positive control replicate 2
# ============================================================

library(shiny)
library(rhandsontable)
library(ggplot2)
library(dplyr)
library(tidyr)
library(minpack.lm)   # robust NLS (Levenberg-Marquardt)

# ── Settings ──────────────────────────────────────────────────────────────────
MIN_CHARS   <- 20     # minimum characters for each written response
STEP_LABELS <- c("Background", "Enter data", "Interpret", "Results")

# ── Default absorbance data: rows = A-H, columns = 1-12 ───────────────────────
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
  keep <- is.finite(conc) & is.finite(signal) & conc > 0
  x <- conc[keep]; y <- signal[keep]
  if (length(x) < 4) return(NULL)
  
  Bmax0 <- max(y, na.rm = TRUE)
  Kd0   <- median(x)
  
  tryCatch(
    nlsLM(y ~ one_site(x, Bmax, Kd),
          start   = list(Bmax = Bmax0, Kd = Kd0),
          lower   = c(Bmax = 0, Kd = 0),
          control = nls.lm.control(maxiter = 200)),
    error = function(e) NULL
  )
}

# ── Written questions (Step 3) ────────────────────────────────────────────────
questions <- list(
  list(id = "q1", title = "Signal vs. concentration",
       text = paste("Does the signal increase with concentration for each sample?",
                    "Does it look like it is levelling off (saturating) at the highest concentrations?")),
  list(id = "q2", title = "The blank rows",
       text = paste("Look at the dashed blank rows. Are they flat and low, or do they change with concentration?",
                    "What would a blank that rises with concentration suggest about your assay?")),
  list(id = "q3", title = "Novel protein vs. positive control",
       text = paste("Which sample do you predict has the lower Kd (tighter binding), and why?",
                    "What features of the plots support your prediction?")),
  list(id = "q4", title = "Replicates and outliers",
       text = paste("Do the two replicates for each sample agree with each other?",
                    "Do any individual points look like outliers, and what might have caused them?"))
)
gate_text_ids <- c("q1", "q2", "q3", "q4", "q5")
all_text_ids  <- c(gate_text_ids, "r1", "r2")

# ── CSS / JS ──────────────────────────────────────────────────────────────────
app_css <- "
  body { font-family: 'Segoe UI', Arial, sans-serif; background: #f5f7fa; }
  .well { background: #ffffff; border: 1px solid #dde3ea; border-radius: 8px; box-shadow: none; }
  h2 { color: #1a3a5c; font-weight: 700; margin-bottom: 4px; }
  h3, h4 { color: #2c5f8a; font-weight: 600; }
  .section-title { color: #1a3a5c; font-weight: 600; font-size: 15px;
                   border-bottom: 2px solid #2c8adf; padding-bottom: 4px; margin: 18px 0 10px; }
  .result-box { background: #eaf4ff; border-left: 4px solid #2c8adf; border-radius: 6px;
                padding: 12px 16px; margin: 8px 0; }
  .result-box.positive { background: #eaffea; border-left-color: #2ca02c; }
  .result-label { font-size: 12px; color: #555; margin: 0; }
  .result-value { font-size: 22px; font-weight: 700; color: #1a3a5c; margin: 0; }
  .result-se    { font-size: 12px; color: #666; }
  .btn-primary { background-color: #2c8adf; border-color: #1a6ab0; }
  .btn-primary:hover { background-color: #1a6ab0; }
  .note { font-size: 12px; color: #777; }

  /* progress indicator */
  .progress-steps { display: flex; gap: 8px; margin: 6px 0 18px; }
  .step-pill { flex: 1; display: flex; align-items: center; gap: 8px; padding: 8px 12px;
               border-radius: 22px; font-size: 13px; font-weight: 600;
               border: 1px solid #dde3ea; background: #fff; color: #8a97a6; }
  .step-pill .step-num { width: 24px; height: 24px; border-radius: 50%; background: #dde3ea;
                         color: #fff; display: flex; align-items: center; justify-content: center;
                         font-size: 12px; flex-shrink: 0; }
  .step-pill.active { border-color: #2c8adf; background: #eaf4ff; color: #1a3a5c; }
  .step-pill.active .step-num { background: #2c8adf; }
  .step-pill.done { border-color: #b7e1b7; background: #f3fff3; color: #2ca02c; }
  .step-pill.done .step-num { background: #2ca02c; }

  /* background page */
  .intro-wrap { max-width: 950px; margin: 0 auto; }
  .intro-wrap p, .intro-wrap li { font-size: 15px; line-height: 1.6; color: #333; }
  .eq-box { background: #f1f5f9; border-left: 4px solid #2c5f8a; border-radius: 6px;
            padding: 8px 16px; margin: 12px 0; text-align: center; }
  .key-idea { background: #fff8e5; border-left: 4px solid #f0a500; border-radius: 6px;
              padding: 10px 16px; margin: 12px 0; }

  /* questions */
  .question-card { background: #f9fbfd; border: 1px solid #dde3ea; border-radius: 8px;
                   padding: 12px 16px; margin-bottom: 14px; }
  .q-title { font-weight: 700; color: #1a3a5c; margin-bottom: 2px; }
  .question-card textarea { resize: vertical; }
  .count { font-size: 12px; color: #b36b00; }
  .count.ok { color: #2ca02c; }

  /* navigation */
  .nav-row { display: flex; gap: 10px; align-items: center; margin: 10px 0 30px; }
  .reading-guide { font-size: 13px; color: #444; background: #f1f5f9; border-radius: 6px;
                   padding: 8px 14px; margin-top: 8px; }
"

app_js <- "
  Shiny.addCustomMessageHandler('scrollTop', function(x) { window.scrollTo(0, 0); });
"

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML(app_css)), tags$script(HTML(app_js))),
  uiOutput("layout_css"),   # hides the sidebar on the Background step
  
  titlePanel(
    div(h2("ELISA Dissociation Constant (K\u1d48) Calculator"),
        p("A guided workflow  |  One-site specific binding model",
          style = "color:#555; font-size:13px; margin:0;")),
    windowTitle = "ELISA Kd Calculator"
  ),
  
  uiOutput("progress"),
  
  fluidRow(
    # Persistent sidebar (always in the DOM so inputs keep their values)
    div(class = "col-sm-3", id = "sidebar_col",
        wellPanel(
          div(class = "section-title", "\u2699\ufe0f  Dilution Series"),
          numericInput("start_conc", "Highest concentration (µM):", value = 10, min = 0, step = 0.1),
          numericInput("fold",       "Fold dilution:",               value = 2,  min = 1.01, step = 0.5),
          p(class = "note", "Column 1 = highest, Column 12 = lowest."),
          hr(),
          div(class = "section-title", "\ud83d\udcc8  Display"),
          radioButtons("x_scale", "X-Axis Scale:", choices = c("Linear", "Log10"), inline = TRUE),
          checkboxInput("show_rep",   "Show individual replicates", value = TRUE),
          checkboxInput("show_resid", "Show residual plot",         value = FALSE)
        )
    ),
    div(class = "col-sm-9", id = "main_col",
        uiOutput("step_ui")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ══════════════════════════════════════════════════════════════════════════
  #  Shared state
  # ══════════════════════════════════════════════════════════════════════════
  current_step <- reactiveVal(1)
  
  # Plate data (survives navigation because it lives here, not in the UI)
  rv <- reactiveVal(default_abs)
  
  # Student answers (survive navigation; UI is re-created on each step change)
  answers <- do.call(
    reactiveValues,
    c(setNames(as.list(rep("", length(all_text_ids))), all_text_ids),
      list(pred_novel = NA_real_, pred_ctrl = NA_real_))
  )
  
  conc_vec <- reactive({
    req(input$start_conc, input$fold)
    req(input$start_conc > 0, input$fold > 1)
    input$start_conc / (input$fold ^ (0:11))
  })
  
  # Is the plate data complete and numeric?
  data_ok <- reactive({
    df   <- rv()
    cols <- as.character(1:12)
    if (!all(cols %in% names(df))) return(FALSE)
    if (!setequal(df$Row, LETTERS[1:8])) return(FALSE)
    all(vapply(df[cols], function(x) all(is.finite(suppressWarnings(as.numeric(x)))), logical(1)))
  })
  
  # Are all gate questions answered?
  n_answered <- reactive({
    txt <- sum(vapply(gate_text_ids[1:4],
                      function(id) nchar(trimws(answers[[id]])) >= MIN_CHARS, logical(1)))
    pred_ok <- is.finite(answers$pred_novel) && answers$pred_novel > 0 &&
      is.finite(answers$pred_ctrl)  && answers$pred_ctrl  > 0 &&
      nchar(trimws(answers$q5)) >= MIN_CHARS
    txt + as.integer(pred_ok)
  })
  all_answered <- reactive(n_answered() == 5)
  
  # ══════════════════════════════════════════════════════════════════════════
  #  Navigation
  # ══════════════════════════════════════════════════════════════════════════
  observeEvent(input$start,      current_step(2))
  observeEvent(input$back_1,     current_step(1))
  observeEvent(input$to_q,       { req(data_ok());     current_step(3) })
  observeEvent(input$back_2,     current_step(2))
  observeEvent(input$to_results, { req(all_answered()); current_step(4) })
  observeEvent(input$back_3,     current_step(3))
  
  observeEvent(current_step(), session$sendCustomMessage("scrollTop", TRUE))
  
  output$layout_css <- renderUI({
    if (current_step() == 1)
      tags$style(HTML("#sidebar_col { display: none !important; }
                       #main_col { width: 100% !important; }"))
  })
  
  output$progress <- renderUI({
    s <- current_step()
    pills <- lapply(seq_along(STEP_LABELS), function(i) {
      cls <- if (i < s) "done" else if (i == s) "active" else "todo"
      div(class = paste("step-pill", cls),
          span(class = "step-num", if (i < s) "\u2713" else i),
          span(STEP_LABELS[i]))
    })
    div(class = "progress-steps", pills)
  })
  
  # ══════════════════════════════════════════════════════════════════════════
  #  Keep answers in sync with the text boxes
  # ══════════════════════════════════════════════════════════════════════════
  lapply(all_text_ids, function(id) {
    observeEvent(input[[id]], { answers[[id]] <- input[[id]] }, ignoreNULL = TRUE)
  })
  observeEvent(input$pred_novel, { answers$pred_novel <- input$pred_novel }, ignoreNULL = TRUE)
  observeEvent(input$pred_ctrl,  { answers$pred_ctrl  <- input$pred_ctrl  }, ignoreNULL = TRUE)
  
  # Live character counters
  lapply(gate_text_ids, function(id) {
    output[[paste0("count_", id)]] <- renderUI({
      n <- nchar(trimws(answers[[id]]))
      if (n >= MIN_CHARS) span(class = "count ok", sprintf("\u2713 %d characters", n))
      else                span(class = "count",    sprintf("%d / %d characters minimum", n, MIN_CHARS))
    })
  })
  
  # ══════════════════════════════════════════════════════════════════════════
  #  Step container
  #  NOTE: nothing reactive is read directly in here (answers are isolated),
  #  so this only re-renders when the step changes.
  # ══════════════════════════════════════════════════════════════════════════
  output$step_ui <- renderUI({
    switch(as.character(current_step()),
           "1" = step1_ui(),
           "2" = step2_ui(),
           "3" = step3_ui(),
           "4" = step4_ui())
  })
  
  # ──────────────────────────────────────────────────────────────────────────
  #  STEP 1: Background
  # ──────────────────────────────────────────────────────────────────────────
  step1_ui <- function() {
    withMathJax(
      div(class = "intro-wrap",
          wellPanel(
            h3("What is a dissociation constant (K\u1d48)?"),
            p("When a binding protein meets its target, the two stick together and come apart again.",
              "The ", strong("dissociation constant, K\u1d48"), ", is a single number that describes how tightly they stick."),
            div(class = "key-idea",
                strong("Key idea: "),
                "K\u1d48 is the concentration of binding protein at which ", strong("half"),
                " of the available binding sites are occupied. It is measured in concentration units (here, µM).",
                tags$br(),
                strong("A lower K\u1d48 means tighter binding"),
                ", because less protein is needed to fill half the sites."),
            
            h4("How is it calculated from an ELISA?"),
            p("In an ELISA, the target is stuck to the bottom of a plate and the binding protein is added in a ",
              strong("serial dilution"), ": the highest concentration goes in column 1, and each column to the right",
              "has the previous concentration divided by a fixed fold-dilution. Bound protein is detected with an enzyme-linked",
              "reagent that produces colour, which the plate reader measures as ", strong("absorbance"),
              ". More bound protein gives more signal."),
            p("Each signal row has a paired ", strong("blank"), " row that measures background (non-specific) signal.",
              "Subtracting the blank from the signal gives the ", em("net signal"), " that is due to specific binding."),
            p("Net signal plotted against concentration gives a curve that rises and then levels off. We describe it with the",
              strong("one-site specific binding"), " equation:"),
            
            div(class = "eq-box",
                "$$Y = \\frac{B_{max}\\,[X]}{K_d + [X]}$$"),
            
            tags$ul(
              tags$li(strong("Y"), ": the net signal (blank-subtracted absorbance)"),
              tags$li(strong("[X]"), ": the concentration of binding protein in the well"),
              tags$li(strong("B"), tags$sub("max"), ": the maximum signal, reached when every binding site is occupied"),
              tags$li(strong("K"), tags$sub("d"), ": the concentration at which Y = B", tags$sub("max"), "/2")
            ),
            p("When [X] is much lower than K\u1d48 the curve rises almost in a straight line. When [X] is much higher than K\u1d48 it flattens out near B",
              tags$sub("max"), ". The app finds the values of B", tags$sub("max"),
              " and K\u1d48 that make the curve sit as close as possible to your data (non-linear least squares)."),
            
            div(class = "section-title", "Try it: how do B\u2098\u2090\u2093 and K\u1d48 change the curve?"),
            fluidRow(
              column(4,
                     sliderInput("sim_kd",   "K\u1d48 (µM)", min = 0.1, max = 10, value = 1.5, step = 0.1),
                     sliderInput("sim_bmax", "B\u2098\u2090\u2093 (signal)", min = 0.5, max = 2, value = 1.2, step = 0.1),
                     checkboxInput("sim_log", "Log10 x-axis", value = FALSE),
                     p(class = "note", "The blue dashed lines mark K\u1d48 on the x-axis and half of B\u2098\u2090\u2093 on the y-axis. Notice where they meet the curve.")
              ),
              column(8, plotOutput("sim_plot", height = "320px"))
            ),
            
            hr(),
            actionButton("start", "Start: enter my data", class = "btn-primary btn-lg",
                         icon = icon("arrow-right"))
          )
      )
    )
  }
  
  output$sim_plot <- renderPlot({
    req(input$sim_kd, input$sim_bmax)
    kd <- input$sim_kd; bm <- input$sim_bmax
    log_x <- isTRUE(input$sim_log)
    x <- if (log_x) exp(seq(log(0.01), log(100), length.out = 300)) else seq(0, 20, length.out = 300)
    df <- data.frame(x = x, y = one_site(x, bm, kd))
    
    g <- ggplot(df, aes(x, y)) +
      geom_hline(yintercept = bm,     linetype = "dotted", colour = "grey50") +
      geom_hline(yintercept = bm / 2, linetype = "dashed", colour = "#2c8adf") +
      geom_vline(xintercept = kd,     linetype = "dashed", colour = "#2c8adf") +
      geom_line(linewidth = 1.2, colour = "#1a3a5c") +
      annotate("point", x = kd, y = bm / 2, size = 4, colour = "#2c8adf") +
      coord_cartesian(ylim = c(0, 2.1)) +
      labs(x = "Binding protein concentration [X] (µM)", y = "Net signal (Y)",
           title = sprintf("Kd = %.1f µM,  Bmax = %.1f", kd, bm)) +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", colour = "#1a3a5c"))
    if (log_x) g <- g + scale_x_log10(labels = function(x) format(x, scientific = FALSE))
    g
  })
  
  # ──────────────────────────────────────────────────────────────────────────
  #  STEP 2: Enter data + raw plots
  # ──────────────────────────────────────────────────────────────────────────
  step2_ui <- function() {
    tagList(
      wellPanel(
        div(class = "section-title", "Raw Absorbance Data"),
        p(icon("paste"),
          tags$strong("You can copy and paste your data directly from the plate reader software into the table below."),
          " Select the 8 x 12 block of values from your spreadsheet, click the top-left data cell (Row A, Col 1), and paste.",
          style = "font-size:13px; color:#444; margin-bottom:8px;"),
        p(class = "note", style = "margin-bottom:8px;",
          "Layout: rows A, C, E, G = signal | rows B, D, F, H = blank. The 'Row' column and the 'Conc' row are read-only.",
          " Set your dilution series in the panel on the left."),
        rHandsontableOutput("hot_table"),
        uiOutput("data_msg")
      ),
      wellPanel(
        div(class = "section-title", "\ud83d\udd0d Raw Data Quality Control"),
        plotOutput("raw_plot_s2", height = "380px"),
        raw_plot_guide()
      ),
      uiOutput("nav_2")
    )
  }
  
  raw_plot_guide <- function() {
    div(class = "reading-guide",
        strong("How to read this plot: "),
        "Each panel is one replicate pair. Solid lines are the signal row and dashed lines are the paired blank row.",
        " Look at the shape of the solid lines (do they rise and level off?) and at whether the dashed lines stay flat and low.",
        " Take a few minutes to look at your data before moving on; you'll be asked about it next.")
  }
  
  output$hot_table <- renderRHandsontable({
    df <- isolate(rv())   # isolate: only re-render when concentrations change, not on every edit
    
    conc_vals <- round(conc_vec(), 4)
    conc_row  <- data.frame(Row = "Conc", t(conc_vals), stringsAsFactors = FALSE)
    colnames(conc_row) <- colnames(df)
    display_df <- rbind(conc_row, df)
    
    rhandsontable(display_df,
                  rowHeaders = NULL,
                  stretchH   = "all",
                  height     = 275) %>%
      hot_cols(format = "0.0000", renderer = "
        function(instance, td, row, col, prop, value, cellProperties) {
          if (col === 0) {
            Handsontable.renderers.TextRenderer.apply(this, arguments);
          } else {
            Handsontable.renderers.NumericRenderer.apply(this, arguments);
          }
          if (row === 0) {
            cellProperties.readOnly = true;
            td.style.background = '#f1f5f9';
            td.style.color = '#2c5f8a';
            td.style.fontWeight = 'bold';
          }
          if (col === 0 && row > 0) {
            cellProperties.readOnly = true;
            td.style.background = '#f9f9f9';
          }
        }
      ")
  })
  
  observeEvent(input$hot_table, {
    new_df   <- hot_to_r(input$hot_table)
    clean_df <- new_df[new_df$Row != "Conc", ]
    rv(clean_df)
  })
  
  output$data_msg <- renderUI({
    if (!data_ok())
      div(style = "color:#b30000; font-size:13px; margin-top:8px;",
          icon("exclamation-triangle"),
          " Some cells are empty or non-numeric. Every cell in the 8 x 12 block needs a number before you can continue.")
  })
  
  output$nav_2 <- renderUI({
    div(class = "nav-row",
        actionButton("back_1", "Back", icon = icon("arrow-left")),
        if (data_ok())
          actionButton("to_q", "Continue to questions", class = "btn-primary", icon = icon("arrow-right"))
        else
          actionButton("to_q", "Continue to questions", class = "btn-primary", disabled = "disabled"))
  })
  
  # ── Raw data plot (pre-subtraction) ───────────────────────────────────────
  raw_plot_obj <- reactive({
    req(data_ok(), conc_vec())
    df   <- rv()
    conc <- conc_vec()
    
    raw_long <- df %>%
      pivot_longer(cols = -Row, names_to = "Column", values_to = "Absorbance") %>%
      mutate(Column = as.numeric(Column))
    
    conc_df  <- data.frame(Column = 1:12, Concentration = conc)
    raw_long <- raw_long %>% left_join(conc_df, by = "Column")
    
    raw_long <- raw_long %>%
      mutate(
        Type   = ifelse(Row %in% c("A", "C", "E", "G"), "Signal", "Blank"),
        Sample = ifelse(Row %in% c("A", "B", "C", "D"), "Novel Binding Protein", "Positive Control"),
        Pair   = case_when(
          Row %in% c("A", "B") ~ "Novel Protein - Rep 1 (Rows A & B)",
          Row %in% c("C", "D") ~ "Novel Protein - Rep 2 (Rows C & D)",
          Row %in% c("E", "F") ~ "Positive Control - Rep 1 (Rows E & F)",
          Row %in% c("G", "H") ~ "Positive Control - Rep 2 (Rows G & H)"
        )
      )
    
    g_raw <- ggplot(raw_long, aes(x = Concentration, y = Absorbance, color = Sample, linetype = Type)) +
      geom_point(size = 2.5) +
      geom_line(linewidth = 1) +
      facet_wrap(~ Pair, ncol = 2) +
      scale_linetype_manual(values = c("Signal" = "solid", "Blank" = "dashed"), name = "Row Type") +
      scale_color_manual(
        values = c("Novel Binding Protein" = "#2c8adf", "Positive Control" = "#2ca02c"),
        guide = "none"
      ) +
      labs(x = "Concentration (µM)", y = "Raw Absorbance (OD)",
           title = "Raw Binding vs. Negative Control Background") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", colour = "#1a3a5c"),
            legend.position = "right")
    
    if (input$x_scale == "Log10")
      g_raw <- g_raw + scale_x_log10(labels = function(x) format(x, scientific = FALSE))
    
    g_raw
  })
  
  output$raw_plot_s2 <- renderPlot(raw_plot_obj())
  output$raw_plot_s3 <- renderPlot(raw_plot_obj())
  
  # ──────────────────────────────────────────────────────────────────────────
  #  STEP 3: Interpret your data (the gate)
  # ──────────────────────────────────────────────────────────────────────────
  question_card <- function(number, title, text, id) {
    div(class = "question-card",
        div(class = "q-title", sprintf("Question %d: %s", number, title)),
        p(text, style = "margin-bottom:6px;"),
        textAreaInput(id, label = NULL, value = isolate(answers[[id]]), rows = 3, width = "100%",
                      placeholder = "Type your response here..."),
        uiOutput(paste0("count_", id)))
  }
  
  step3_ui <- function() {
    qcards <- lapply(seq_along(questions), function(i) {
      q <- questions[[i]]
      question_card(i, q$title, q$text, q$id)
    })
    
    q5 <- div(class = "question-card",
              div(class = "q-title", "Question 5: Estimate the half-maximal concentration"),
              p("Using the raw plots, make a rough estimate of the concentration at which the signal reaches",
                "half of its maximum for each sample. Then explain how you made your estimate.",
                style = "margin-bottom:6px;"),
              fluidRow(
                column(6, numericInput("pred_novel", "Novel binding protein (µM):",
                                       value = isolate(answers$pred_novel), min = 0, step = 0.1)),
                column(6, numericInput("pred_ctrl", "Positive control (µM):",
                                       value = isolate(answers$pred_ctrl), min = 0, step = 0.1))
              ),
              textAreaInput("q5", label = NULL, value = isolate(answers$q5), rows = 3, width = "100%",
                            placeholder = "How did you arrive at these estimates?"),
              uiOutput("count_q5"))
    
    tagList(
      wellPanel(
        div(class = "section-title", "Your raw data"),
        plotOutput("raw_plot_s3", height = "380px"),
        raw_plot_guide()
      ),
      wellPanel(
        div(class = "section-title", "\u270d\ufe0f Interpret your data"),
        p(class = "note",
          sprintf("Answer all five questions (at least %d characters each) to unlock the K\u1d48 calculation. ", MIN_CHARS),
          "There are no trick questions: the aim is to make your own predictions before you see the fit."),
        qcards,
        q5
      ),
      uiOutput("nav_3")
    )
  }
  
  output$nav_3 <- renderUI({
    n <- n_answered()
    div(class = "nav-row",
        actionButton("back_2", "Back", icon = icon("arrow-left")),
        if (n == 5)
          actionButton("to_results", "Calculate K\u1d48", class = "btn-primary", icon = icon("calculator"))
        else
          tagList(
            actionButton("to_results", "Calculate K\u1d48", class = "btn-primary", disabled = "disabled"),
            span(class = "note", sprintf(" %d of 5 questions complete.", n))
          ))
  })
  
  # ──────────────────────────────────────────────────────────────────────────
  #  STEP 4: Results
  # ──────────────────────────────────────────────────────────────────────────
  step4_ui <- function() {
    tagList(
      uiOutput("results_boxes"),
      uiOutput("comparison"),
      wellPanel(
        div(class = "section-title", "\ud83d\udcc9 Binding Curves"),
        plotOutput("binding_plot", height = "380px"),
        uiOutput("resid_ui")
      ),
      wellPanel(
        div(class = "section-title", "\ud83e\udde0 Reflect"),
        p(class = "note", "These reflection boxes are optional and are not saved or sent anywhere. Copy your text somewhere if you want to keep it."),
        div(class = "question-card",
            div(class = "q-title", "Your prediction vs. the fit"),
            p("How close were your estimates to the fitted K\u1d48 values? If they differ, what might explain it?",
              style = "margin-bottom:6px;"),
            textAreaInput("r1", label = NULL, value = isolate(answers$r1), rows = 3, width = "100%",
                          placeholder = "Type your reflection here...")),
        div(class = "question-card",
            div(class = "q-title", "Quality of the fit"),
            p("Does the curve describe your data well? Tick 'Show residual plot' in the sidebar: what do the residuals tell you?",
              style = "margin-bottom:6px;"),
            textAreaInput("r2", label = NULL, value = isolate(answers$r2), rows = 3, width = "100%",
                          placeholder = "Type your reflection here..."))
      ),
      div(class = "nav-row", actionButton("back_3", "Back to questions", icon = icon("arrow-left")))
    )
  }
  
  # ── Processed data (blank-subtracted, net signals) ────────────────────────
  processed <- reactive({
    req(current_step() == 4, all_answered(), data_ok())
    df   <- rv()
    conc <- conc_vec()
    
    get_row <- function(row_label) as.numeric(df[df$Row == row_label, as.character(1:12)])
    
    novel_net <- rbind(
      data.frame(rep = "Rep 1", conc = conc, net = get_row("A") - get_row("B")),
      data.frame(rep = "Rep 2", conc = conc, net = get_row("C") - get_row("D"))
    )
    ctrl_net <- rbind(
      data.frame(rep = "Rep 1", conc = conc, net = get_row("E") - get_row("F")),
      data.frame(rep = "Rep 2", conc = conc, net = get_row("G") - get_row("H"))
    )
    # NOTE: negatives are deliberately NOT clamped to 0, so non-specific binders
    # don't appear to have artificially low Kd values.
    
    novel_mean <- novel_net %>% group_by(conc) %>% summarise(mean_net = mean(net), .groups = "drop")
    ctrl_mean  <- ctrl_net  %>% group_by(conc) %>% summarise(mean_net = mean(net), .groups = "drop")
    
    list(novel_rep = novel_net, ctrl_rep = ctrl_net,
         novel_mean = novel_mean, ctrl_mean = ctrl_mean)
  })
  
  fits <- reactive({
    p <- processed()
    list(
      novel = fit_kd(p$novel_rep$conc, p$novel_rep$net),
      ctrl  = fit_kd(p$ctrl_rep$conc,  p$ctrl_rep$net)
    )
  })
  
  # ── Kd result boxes ───────────────────────────────────────────────────────
  output$results_boxes <- renderUI({
    f <- fits()
    
    make_box <- function(fit, label, cls) {
      if (is.null(fit)) {
        div(class = paste("result-box", cls),
            p(class = "result-label", label),
            p(class = "result-value", "Fit failed"),
            p(class = "result-se", "Check your data: make sure there is sufficient signal above the blank."))
      } else {
        co   <- summary(fit)$coefficients
        Kd   <- co["Kd",   "Estimate"]; KdSE <- co["Kd",   "Std. Error"]
        Bmax <- co["Bmax", "Estimate"]; BmSE <- co["Bmax", "Std. Error"]
        y    <- fitted(fit) + residuals(fit)
        r2   <- 1 - sum(residuals(fit)^2) / sum((y - mean(y))^2)
        div(class = paste("result-box", cls),
            p(class = "result-label", label),
            p(class = "result-value", sprintf("K\u1d48 = %.4f µM", Kd)),
            p(class = "result-se",
              sprintf("\u00b1 %.4f µM (SE)   |   B\u2098\u2090\u2093 = %.4f \u00b1 %.4f   |   R\u00b2 = %.4f",
                      KdSE, Bmax, BmSE, r2)))
      }
    }
    
    wellPanel(
      div(class = "section-title", "\ud83d\udcca Results"),
      fluidRow(
        column(6, make_box(f$novel, "Novel Binding Protein", "")),
        column(6, make_box(f$ctrl,  "Positive Control",      "positive"))
      ),
      p(class = "note",
        "The fit uses every individual replicate point (both replicates, blank-subtracted), not the replicate means. ",
        "The SE therefore reflects the scatter between replicates as well as the curve shape.")
    )
  })
  
  # ── Prediction vs. fitted Kd ──────────────────────────────────────────────
  output$comparison <- renderUI({
    f <- fits()
    
    row_for <- function(label, est, fit) {
      if (is.null(fit))
        return(tags$tr(tags$td(label), tags$td(sprintf("%.3f", est)), tags$td("Fit failed"), tags$td("-")))
      kd <- coef(fit)["Kd"]
      tags$tr(tags$td(label),
              tags$td(sprintf("%.3f", est)),
              tags$td(sprintf("%.3f", kd)),
              tags$td(sprintf("%+.0f%%", (est - kd) / kd * 100)))
    }
    
    wellPanel(
      div(class = "section-title", "\ud83d\udd2e Your prediction vs. the fitted K\u1d48"),
      tags$table(class = "table table-condensed", style = "max-width:600px;",
                 tags$thead(tags$tr(tags$th("Sample"), tags$th("Your estimate (µM)"),
                                    tags$th("Fitted K\u1d48 (µM)"), tags$th("Estimate vs. fit"))),
                 tags$tbody(
                   row_for("Novel Binding Protein", answers$pred_novel, f$novel),
                   row_for("Positive Control",      answers$pred_ctrl,  f$ctrl)
                 )),
      p(class = "note", "A positive percentage means your estimate was higher than the fitted K\u1d48.")
    )
  })
  
  # ── Binding plot ──────────────────────────────────────────────────────────
  output$binding_plot <- renderPlot({
    req(processed(), fits())
    p <- processed()
    f <- fits()
    
    if (input$x_scale == "Log10") {
      x_seq <- exp(seq(log(min(conc_vec()) * 0.1), log(max(conc_vec()) * 1.5), length.out = 400))
    } else {
      x_seq <- seq(0, max(conc_vec()) * 1.05, length.out = 400)
    }
    
    pred_df <- bind_rows(
      if (!is.null(f$novel)) data.frame(x = x_seq, y = predict(f$novel, newdata = list(x = x_seq)), Sample = "Novel Binding Protein") else NULL,
      if (!is.null(f$ctrl))  data.frame(x = x_seq, y = predict(f$ctrl,  newdata = list(x = x_seq)), Sample = "Positive Control") else NULL
    )
    
    mean_df <- bind_rows(
      data.frame(p$novel_mean, Sample = "Novel Binding Protein"),
      data.frame(p$ctrl_mean,  Sample = "Positive Control")
    )
    
    g <- ggplot(mean_df, aes(x = conc, y = mean_net, colour = Sample, fill = Sample)) +
      geom_hline(yintercept = 0, linetype = "solid", colour = "grey75", linewidth = 1) +
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
    
    kd_df <- bind_rows(
      if (!is.null(f$novel)) data.frame(Sample = "Novel Binding Protein", kd = coef(f$novel)["Kd"]) else NULL,
      if (!is.null(f$ctrl))  data.frame(Sample = "Positive Control",      kd = coef(f$ctrl)["Kd"]) else NULL
    )
    
    if (nrow(kd_df) > 0) {
      y_pos <- max(mean_df$mean_net, na.rm = TRUE) * 0.95
      g <- g +
        geom_vline(data = kd_df, aes(xintercept = kd, colour = Sample),
                   linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
        geom_label(data = kd_df, aes(x = kd, y = y_pos, label = sprintf("Kd = %.3f µM", kd), colour = Sample),
                   fill = "white", label.size = 0.3, hjust = -0.08, size = 4.5,
                   fontface = "bold", show.legend = FALSE)
    }
    
    if (input$x_scale == "Log10")
      g <- g + scale_x_log10(labels = function(x) format(x, scientific = FALSE))
    
    g
  })
  
  # ── Residual plot (toggle lives in its own uiOutput so the step doesn't re-render)
  output$resid_ui <- renderUI({
    if (isTRUE(input$show_resid)) plotOutput("resid_plot", height = "250px")
  })
  
  output$resid_plot <- renderPlot({
    req(processed(), fits())
    f <- fits()
    
    resid_df <- bind_rows(
      if (!is.null(f$novel)) data.frame(
        fitted = fitted(f$novel), resid = residuals(f$novel), Sample = "Novel Binding Protein") else NULL,
      if (!is.null(f$ctrl)) data.frame(
        fitted = fitted(f$ctrl),  resid = residuals(f$ctrl),  Sample = "Positive Control") else NULL
    )
    
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