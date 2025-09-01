# app.R
# shiny application for semiparametrically efficient stepped wedge designs

library(shiny)
library(ggplot2)

# source optimization functions
source("util_optimization.R")

# ui definition
ui <- fluidPage(
  withMathJax(),
  tags$head(tags$style(HTML("
    body { font-family: 'Times New Roman', serif; }
    .well { background-color: #f9f9f9; border: 1px solid #e3e3e3; }
    h3 { color: #333; margin-bottom: 15px; }
    .equation-display { text-align: center; margin: 20px 0; font-size: 18px; }
  "))),
  
  titlePanel("Semiparametrically Efficient Stepped Wedge Designs"),
  
  tabsetPanel(
    # individually randomized trial tab
    tabPanel("Individually Randomized Trial",
             br(),
             fluidRow(
               column(12,
                      wellPanel(
                        h3("Model Specification"),
                        div(class = "equation-display",
                            HTML("\\(Y^{\\text{obs}}_{ij} = \\beta_j + (\\tau_{j0} \\cdot Z_{ij} + \\delta_{j1} \\cdot Z_{i, j - 1} + \\cdots + \\delta_{j\\ell} \\cdot Z_{i, j - \\ell}) + \\upsilon_i \\cdot Z_{ij} + \\varepsilon_{ij}\\)")
                        ),
                        tags$ul(style = "list-style-type: none; margin-left: 20px;",
                                tags$li(HTML("• \\(Y_{ij}^{\\text{obs}}\\) is the observed outcome for unit \\(i\\) in period \\(j\\)")),
                                tags$li(HTML("• \\(Z_{ij}\\) is the binary treatment indicator for unit \\(i\\) in period \\(j\\)")),
                                tags$li(HTML("• \\(\\beta_j\\) is the time fixed effect")),
                                tags$li(HTML("• \\(\\delta_{js} = \\tau_{js} - \\tau_{j, s - 1}\\) for \\(s \\in \\{1, \\ldots, \\ell\\}\\) captures incremental carryover effects")),
                                tags$li(HTML("• \\(\\upsilon_i \\sim \\mathcal{N}(0, \\sigma_\\upsilon^2)\\) is a random unit-specific treatment effect")),
                                tags$li(HTML("• \\(\\boldsymbol{\\varepsilon}_i = (\\varepsilon_{i1}, \\cdots, \\varepsilon_{iJ})' \\sim \\mathcal{N}(0, \\sigma_\\epsilon^2\\mathbf{M})\\) with \\(\\mathbf{M}_{jj'} = r^{|j - j'|}\\)"))
                        )
                      )
               )
             ),
             
             sidebarLayout(
               sidebarPanel(
                 numericInput("J", 
                              label = HTML("Number of periods \\(J\\):"), 
                              value = 12, min = 3, max = 30, step = 1),
                 
                 numericInput("ell", 
                              label = HTML("Maximum duration of carryover effects \\(\\ell\\):"), 
                              value = 0, min = 0, max = 10, step = 1),
                 
                 radioButtons("estimand", 
                              label = "Target estimand:",
                              choiceNames = list(
                                HTML("\\(\\tau^{\\text{GATE}}\\)"),
                                HTML("\\(\\tau_s^{\\text{CTE}}\\) (s = 0, ..., \\(\\ell-1\\))"),
                                HTML("Both \\(\\tau^{\\text{GATE}}\\) and \\(\\tau_s^{\\text{CTE}}\\) (s = 0, ..., \\(\\ell-1\\))")
                              ),
                              choiceValues = list("gate", "cte", "both"),
                              selected = "gate",
                              inline = FALSE),
                 
                 helpText("Note: When ℓ = 0, all estimands are equivalent as treatment effect is constant."),
                 
                 sliderInput("r", 
                             label = HTML("Serial correlation \\(r\\):"), 
                             min = 0, max = 1, value = 0.9, step = 0.01),
                 
                 numericInput("quo", 
                              label = HTML("Variance ratio \\(\\sigma_{\\upsilon}^2/\\sigma_{\\epsilon}^2\\):"), 
                              value = 0, min = 0, max = 10, step = 0.1),
                 
                 actionButton("run", "Run optimization", class = "btn-primary"),
                 
                 br(), br(),
                 
                 downloadButton("downloadData", "Download optimal design")
               ),
               
               mainPanel(
                 plotOutput("piPlot", width = "600px", height = "600px"),
                 br(),
                 verbatimTextOutput("optSummary")
               )
             )
    ),
    
    # cluster randomized trial tab
    tabPanel("Cluster Randomized Trial",
             br(),
             tabsetPanel(
               # cross-section design tab
               tabPanel("Cross-section",
                        br(),
                        fluidRow(
                          column(12,
                                 wellPanel(
                                   h3("Model Specification"),
                                   p("Different individuals are observed in each cluster over time."),
                                   div(class = "equation-display",
                                       HTML("\\(Y^{\\text{obs}}_{ijk} = \\beta_j + (\\tau_{j0} \\cdot Z_{ij} + \\delta_{j1} \\cdot Z_{i, j - 1} + \\cdots + \\delta_{j\\ell} \\cdot Z_{i, j - \\ell}) + \\upsilon_i \\cdot Z_{ij} + \\gamma_{ij} + \\varepsilon_{ijk}\\)")
                                   ),
                                   tags$ul(style = "list-style-type: none; margin-left: 20px;",
                                           tags$li(HTML("• \\(Y_{ijk}^{\\text{obs}}\\) is the observed outcome for individual \\(k\\) in cluster \\(i\\) in period \\(j\\)")),
                                           tags$li(HTML("• \\(Z_{ij}\\) is the binary treatment indicator for cluster \\(i\\) in period \\(j\\)")),
                                           tags$li(HTML("• \\(\\beta_j\\) is the time fixed effect")),
                                           tags$li(HTML("• \\(\\delta_{js} = \\tau_{js} - \\tau_{j, s - 1}\\) for \\(s \\in \\{1, \\ldots, \\ell\\}\\) captures incremental carryover effects")),
                                           tags$li(HTML("• \\(\\upsilon_i \\sim \\mathcal{N}(0, \\sigma_\\upsilon^2)\\) is a random cluster-specific treatment effect")),
                                           tags$li(HTML("• \\(\\boldsymbol{\\gamma}_i = (\\gamma_{i1}, \\ldots, \\gamma_{iJ})' \\sim \\mathcal{N}(0, \\sigma_\\gamma^2\\mathbf{M})\\) with \\(\\mathbf{M}_{jj'} = r^{|j - j'|}\\)")),
                                           tags$li(HTML("• \\(\\varepsilon_{ijk} \\sim_{\\text{i.i.d.}} \\mathcal{N}(0, \\sigma_\\varepsilon^2)\\) is the individual-level residual error"))
                                   )
                                 )
                          )
                        ),
                        
                        sidebarLayout(
                          sidebarPanel(
                            numericInput("J_crt", 
                                         label = HTML("Number of periods \\(J\\):"), 
                                         value = 12, min = 3, max = 30, step = 1),
                            
                            numericInput("ell_crt", 
                                         label = HTML("Maximum duration of carryover effects \\(\\ell\\):"), 
                                         value = 0, min = 0, max = 10, step = 1),
                            
                            radioButtons("estimand_crt", 
                                         label = "Target estimand:",
                                         choiceNames = list(
                                           HTML("\\(\\tau^{\\text{GATE}}\\)")
                                         ),
                                         choiceValues = list("gate"),
                                         selected = "gate",
                                         inline = FALSE),
                            
                            sliderInput("r_crt", 
                                        label = HTML("Serial correlation \\(r\\):"), 
                                        min = 0, max = 1, value = 0.9, step = 0.01),
                            
                            numericInput("sigma_gamma_sq", 
                                         label = HTML("\\(\\sigma_{\\gamma}^2\\):"), 
                                         value = 1, min = 0, max = 10, step = 0.1),
                            
                            numericInput("sigma_epsilon_sq", 
                                         label = HTML("\\(\\sigma_{\\epsilon}^2\\):"), 
                                         value = 1, min = 0, max = 10, step = 0.1),
                            
                            numericInput("sigma_upsilon_sq", 
                                         label = HTML("\\(\\sigma_{\\upsilon}^2\\):"), 
                                         value = 0, min = 0, max = 10, step = 0.1),
                            
                            numericInput("K_ij", 
                                         label = HTML("Cluster-period size \\(K_{ij}\\):"), 
                                         value = 10, min = 1, max = 1000, step = 10),
                            
                            actionButton("run_crt", "Run optimization", class = "btn-primary"),
                            
                            br(), br(),
                            
                            downloadButton("downloadData_crt", "Download optimal design")
                          ),
                          
                          mainPanel(
                            plotOutput("piPlot_crt", width = "600px", height = "600px"),
                            br(),
                            verbatimTextOutput("optSummary_crt")
                          )
                        )
               ),
               
               # closed-cohort design tab
               tabPanel("Closed-cohort",
                        br(),
                        fluidRow(
                          column(12,
                                 wellPanel(
                                   h3("Model Specification"),
                                   p("Individuals are identified at the start of the trial and scheduled for repeated outcome assessment."),
                                   div(class = "equation-display",
                                       HTML("\\(Y^{\\text{obs}}_{ijk} = \\beta_j + (\\tau_{j0} \\cdot Z_{ij} + \\delta_{j1} \\cdot Z_{i, j - 1} + \\cdots + \\delta_{j\\ell} \\cdot Z_{i, j - \\ell}) + \\upsilon_i \\cdot Z_{ij} + \\gamma_{ij} + \\phi_{ik} + \\varepsilon_{ijk}\\)")
                                   ),
                                   tags$ul(style = "list-style-type: none; margin-left: 20px;",
                                           tags$li(HTML("• \\(Y_{ijk}^{\\text{obs}}\\) is the observed outcome for individual \\(k\\) in cluster \\(i\\) in period \\(j\\)")),
                                           tags$li(HTML("• \\(Z_{ij}\\) is the binary treatment indicator for cluster \\(i\\) in period \\(j\\)")),
                                           tags$li(HTML("• \\(\\beta_j\\) is the time fixed effect")),
                                           tags$li(HTML("• \\(\\delta_{js} = \\tau_{js} - \\tau_{j, s - 1}\\) for \\(s \\in \\{1, \\ldots, \\ell\\}\\) captures incremental carryover effects")),
                                           tags$li(HTML("• \\(\\upsilon_i \\sim \\mathcal{N}(0, \\sigma_\\upsilon^2)\\) is a random cluster-specific treatment effect")),
                                           tags$li(HTML("• \\(\\boldsymbol{\\gamma}_i = (\\gamma_{i1}, \\ldots, \\gamma_{iJ})' \\sim \\mathcal{N}(0, \\sigma_\\gamma^2\\mathbf{M})\\) with \\(\\mathbf{M}_{jj'} = r^{|j - j'|}\\)")),
                                           tags$li(HTML("• \\(\\phi_{ik} \\sim \\mathcal{N}(0, \\sigma_\\phi^2)\\) is the random effect for repeated measures from individual \\(k\\) in cluster \\(i\\)")),
                                           tags$li(HTML("• \\(\\varepsilon_{ijk} \\sim_{\\text{i.i.d.}} \\mathcal{N}(0, \\sigma_\\varepsilon^2)\\) is the individual-level residual error"))
                                   )
                                 )
                          )
                        ),
                        
                        sidebarLayout(
                          sidebarPanel(
                            numericInput("J_crt_cc", 
                                         label = HTML("Number of periods \\(J\\):"), 
                                         value = 12, min = 3, max = 30, step = 1),
                            
                            numericInput("ell_crt_cc", 
                                         label = HTML("Maximum duration of carryover effects \\(\\ell\\):"), 
                                         value = 0, min = 0, max = 10, step = 1),
                            
                            radioButtons("estimand_crt_cc", 
                                         label = "Target estimand:",
                                         choiceNames = list(
                                           HTML("\\(\\tau^{\\text{GATE}}\\)")
                                         ),
                                         choiceValues = list("gate"),
                                         selected = "gate",
                                         inline = FALSE),
                            
                            sliderInput("r_crt_cc", 
                                        label = HTML("Serial correlation \\(r\\):"), 
                                        min = 0, max = 1, value = 0.9, step = 0.01),
                            
                            numericInput("sigma_gamma_sq_cc", 
                                         label = HTML("\\(\\sigma_{\\gamma}^2\\):"), 
                                         value = 1, min = 0, max = 10, step = 0.1),
                            
                            numericInput("sigma_epsilon_sq_cc", 
                                         label = HTML("\\(\\sigma_{\\epsilon}^2\\):"), 
                                         value = 1, min = 0, max = 10, step = 0.1),
                            
                            numericInput("sigma_upsilon_sq_cc", 
                                         label = HTML("\\(\\sigma_{\\upsilon}^2\\):"), 
                                         value = 0, min = 0, max = 10, step = 0.1),
                            
                            numericInput("sigma_phi_sq_cc", 
                                         label = HTML("\\(\\sigma_{\\phi}^2\\):"), 
                                         value = 0, min = 0, max = 10, step = 0.1),
                            
                            numericInput("K_ij_cc", 
                                         label = HTML("Cluster-period size \\(K_{ij}\\):"), 
                                         value = 10, min = 1, max = 1000, step = 10),
                            
                            actionButton("run_crt_cc", "Run optimization", class = "btn-primary"),
                            
                            br(), br(),
                            
                            downloadButton("downloadData_crt_cc", "Download optimal design")
                          ),
                          
                          mainPanel(
                            plotOutput("piPlot_crt_cc", width = "600px", height = "600px"),
                            br(),
                            verbatimTextOutput("optSummary_crt_cc")
                          )
                        )
               )
             )
    )
  )
)

# server logic
server <- function(input, output, session) {
  
  # irt tab functionality
  observeEvent(input$J, {
    updateNumericInput(session, "ell", max = input$J - 2)
  })
  
  optimiseRes <- eventReactive(input$run, {
    optimize_irt(J = input$J, 
                 ell = input$ell, 
                 r = input$r, 
                 quo = input$quo, 
                 estimand = input$estimand)
  })
  
  output$piPlot <- renderPlot({
    res <- optimiseRes()
    req(res)
    df <- data.frame(Period = seq_len(length(res$pi)), Pi = res$pi)
    ggplot(df, aes(Period, Pi)) +
      geom_line(linewidth = 1) + 
      geom_point(size = 3) +
      scale_x_continuous(breaks = seq_len(max(df$Period))) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      labs(x = "Period (j)", y = expression(pi[j]), title = "Optimal Allocation Probabilities") +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
            panel.grid.major = element_line(colour = "grey85", linewidth = 0.25), 
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))
  })
  
  output$optSummary <- renderPrint({
    res <- optimiseRes()
    req(res)
    if(input$ell == 0) {
      cat("Target estimand: All equivalent (ℓ = 0, treatment effect is constant)\n")
    } else {
      cat("Target estimand: ", 
          switch(input$estimand,
                 "gate" = "τ^GATE",
                 "cte" = "τ_s^CTE (s = 0, ..., ℓ-1)",
                 "both" = "Both τ^GATE and τ_s^CTE (s = 0, ..., ℓ-1)"),
          "\n", sep = "")
    }
    cat("Function evaluations: ", res$iterations, "\n", sep = "")
    cat("Objective value: ", signif(res$objective_value, 4), "\n\n", sep = "")
    
    cat("Optimal allocation probabilities (pi):\n")
    for (i in seq_along(res$pi)) {
      cat(sprintf("  pi[%2d] = %.2f\n", i, res$pi[i]))
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("optimal_design_IRT_", input$estimand, "_J", input$J, "_ell", input$ell, "_r", input$r, "_quo", input$quo, ".csv")
    },
    content = function(file) {
      res <- optimiseRes()
      req(res)
      df <- data.frame(
        period = seq_len(length(res$pi)),
        pi = res$pi
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # crt cross-section functionality
  observeEvent(input$J_crt, {
    updateNumericInput(session, "ell_crt", max = input$J_crt - 2)
  })
  
  optimiseRes_crt <- eventReactive(input$run_crt, {
    optimize_crt(design_type = "cross_section",
                 J = input$J_crt,
                 ell = input$ell_crt,
                 r = input$r_crt,
                 sigma_gamma_sq = input$sigma_gamma_sq,
                 sigma_epsilon_sq = input$sigma_epsilon_sq,
                 sigma_upsilon_sq = input$sigma_upsilon_sq,
                 K = input$K_ij)
  })
  
  output$piPlot_crt <- renderPlot({
    res <- optimiseRes_crt()
    req(res)
    df <- data.frame(Period = seq_len(length(res$pi)), Pi = res$pi)
    ggplot(df, aes(Period, Pi)) +
      geom_line(linewidth = 1) + 
      geom_point(size = 3) +
      scale_x_continuous(breaks = seq_len(max(df$Period))) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      labs(x = "Period (j)", y = expression(pi[j]), title = "Optimal Allocation Probabilities") +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
            panel.grid.major = element_line(colour = "grey85", linewidth = 0.25), 
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))
  })
  
  output$optSummary_crt <- renderPrint({
    res <- optimiseRes_crt()
    req(res)
    cat("Target estimand: τ^GATE\n")
    cat("Function evaluations: ", res$iterations, "\n", sep = "")
    cat("Objective value: ", signif(res$objective_value, 4), "\n\n", sep = "")
    
    cat("Optimal allocation probabilities (pi):\n")
    for (i in seq_along(res$pi)) {
      cat(sprintf("  pi[%2d] = %.2f\n", i, res$pi[i]))
    }
  })
  
  output$downloadData_crt <- downloadHandler(
    filename = function() {
      paste0("optimal_design_CRT_gate_J", input$J_crt, "_ell", input$ell_crt, "_r", input$r_crt, 
             "_sg", input$sigma_gamma_sq, "_se", input$sigma_epsilon_sq, "_su", input$sigma_upsilon_sq, "_K", input$K_ij, ".csv")
    },
    content = function(file) {
      res <- optimiseRes_crt()
      req(res)
      df <- data.frame(
        period = seq_len(length(res$pi)),
        pi = res$pi
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # crt closed-cohort functionality
  observeEvent(input$J_crt_cc, {
    updateNumericInput(session, "ell_crt_cc", max = input$J_crt_cc - 2)
  })
  
  optimiseRes_crt_cc <- eventReactive(input$run_crt_cc, {
    optimize_crt(design_type = "closed_cohort",
                 J = input$J_crt_cc,
                 ell = input$ell_crt_cc,
                 r = input$r_crt_cc,
                 sigma_gamma_sq = input$sigma_gamma_sq_cc,
                 sigma_epsilon_sq = input$sigma_epsilon_sq_cc,
                 sigma_upsilon_sq = input$sigma_upsilon_sq_cc,
                 K = input$K_ij_cc,
                 sigma_phi_sq = input$sigma_phi_sq_cc)
  })
  
  output$piPlot_crt_cc <- renderPlot({
    res <- optimiseRes_crt_cc()
    req(res)
    df <- data.frame(Period = seq_len(length(res$pi)), Pi = res$pi)
    ggplot(df, aes(Period, Pi)) +
      geom_line(linewidth = 1) + 
      geom_point(size = 3) +
      scale_x_continuous(breaks = seq_len(max(df$Period))) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      labs(x = "Period (j)", y = expression(pi[j]), title = "Optimal Allocation Probabilities") +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
            panel.grid.major = element_line(colour = "grey85", linewidth = 0.25), 
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))
  })
  
  output$optSummary_crt_cc <- renderPrint({
    res <- optimiseRes_crt_cc()
    req(res)
    cat("Target estimand: τ^GATE\n")
    cat("Function evaluations: ", res$iterations, "\n", sep = "")
    cat("Objective value: ", signif(res$objective_value, 4), "\n\n", sep = "")
    
    cat("Optimal allocation probabilities (pi):\n")
    for (i in seq_along(res$pi)) {
      cat(sprintf("  pi[%2d] = %.2f\n", i, res$pi[i]))
    }
  })
  
  output$downloadData_crt_cc <- downloadHandler(
    filename = function() {
      paste0("optimal_design_CRT_CC_gate_J", input$J_crt_cc, "_ell", input$ell_crt_cc, "_r", input$r_crt_cc, 
             "_sg", input$sigma_gamma_sq_cc, "_se", input$sigma_epsilon_sq_cc, "_su", input$sigma_upsilon_sq_cc, 
             "_sp", input$sigma_phi_sq_cc, "_K", input$K_ij_cc, ".csv")
    },
    content = function(file) {
      res <- optimiseRes_crt_cc()
      req(res)
      df <- data.frame(
        period = seq_len(length(res$pi)),
        pi = res$pi
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
}

# run the application
shinyApp(ui, server)