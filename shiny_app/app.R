library(shiny)
library(ggplot2)
library(nloptr)

# utility functions for stepped wedge design optimization

# ============================================================================
# individually randomized trial (IRT) optimization functions
# ============================================================================

# minimize variance bound for gate
minimize_variance_bound_nloptr <- function(J, ell, r, quo) {
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (quo + 1) / pi[j - ell] + 1 / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]) - r ** abs(jp - j) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
          }
        }
      }
    }
    (term1 + term2)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.1, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints,
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective,
       status = result$status, iterations = result$iterations)
}

# minimize variance for cte
minimize_variance_cte <- function(J, ell, r, quo) {
  
  objective_fn <- function(pi) {
    for(i in 1:(J-1)) {
      if(pi[i+1] - pi[i] < 1e-8) {
        return(1e10)
      }
    }
    
    term_cte <- 0
    if (ell > 0) {
      for (s in 0:(ell - 1)) {
        term3 <- 0
        for (j in (s + 1):J) {
          if (j - s - 1 == 0) {
            pi_prev <- 0
          } else {
            pi_prev <- pi[j - s - 1]
          }
          term3 <- term3 + (quo + 1) / (pi[j - s] - pi_prev) + 1 / (1 - pi[j])
        }
        term4 <- 0
        if ((s + 1) <= (J - 1)) {
          for (j in (s + 1):(J - 1)) {
            for (jp in (j + 1):J) {
              if (j >= jp - s) {
                term4 <- term4 + 2 * (quo + r ** abs(jp - j)) / (1 - pi[j])
              }
            }
          }
        }
        term_cte <- term_cte + (term3 + term4) / ((J - s)^2)
      }
    }
    return(term_cte)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.01, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints, 
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective, 
       status = result$status, iterations = result$iterations)
}

# minimize variance for both gate and cte
minimize_variance_gate_cte <- function(J, ell, r, quo) {
  
  objective_fn <- function(pi) {
    for(i in 1:(J-1)) {
      if(pi[i+1] - pi[i] < 1e-8) {
        return(1e10)
      }
    }
    
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (quo + 1) / pi[j - ell] + 1 / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]) - r ** abs(jp - j) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
          }
        }
      }
    }
    
    term_cte <- 0
    if (ell > 0) {
      for (s in 0:(ell - 1)) {
        term3 <- 0
        for (j in (s + 1):J) {
          if (j - s - 1 == 0) {
            pi_prev <- 0
          } else {
            pi_prev <- pi[j - s - 1]
          }
          term3 <- term3 + (quo + 1) / (pi[j - s] - pi_prev) + 1 / (1 - pi[j])
        }
        term4 <- 0
        if ((s + 1) <= (J - 1)) {
          for (j in (s + 1):(J - 1)) {
            for (jp in (j + 1):J) {
              if (j >= jp - s) {
                term4 <- term4 + 2 * (quo + r ** abs(jp - j)) / (1 - pi[j])
              }
            }
          }
        }
        term_cte <- term_cte + (term3 + term4) / ((J - s)^2)
      }
    }
    
    return((term1 + term2) / ((J - ell)^2) + term_cte)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.01, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints, 
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective, 
       status = result$status, iterations = result$iterations)
}

# ============================================================================
# cluster randomized trial (CRT) optimization functions
# ============================================================================

# minimize variance bound for crt cross-section design
minimize_variance_crt_cross_section <- function(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, K) {
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (sigma_gamma_sq + sigma_epsilon_sq / K + sigma_upsilon_sq) / pi[j - ell] + (sigma_gamma_sq + sigma_epsilon_sq / K) / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j)) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j)) / (1 - pi[j]) - (sigma_gamma_sq * r ** abs(jp - j)) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
          }
        }
      }
    }
    (term1 + term2)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.1, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints,
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective,
       status = result$status, iterations = result$iterations)
}

# minimize variance bound for crt closed-cohort design
minimize_variance_crt_closed_cohort <- function(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, sigma_phi_sq, K) {
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (sigma_gamma_sq + sigma_epsilon_sq / K + sigma_phi_sq / K + sigma_upsilon_sq) / pi[j - ell] + (sigma_gamma_sq + sigma_epsilon_sq / K + sigma_phi_sq / K) / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq + sigma_phi_sq / K) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j) + sigma_phi_sq / K) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq + sigma_phi_sq / K) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j) + sigma_phi_sq / K) / (1 - pi[j]) - (sigma_gamma_sq * r ** abs(jp - j) + sigma_phi_sq / K) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
          }
        }
      }
    }
    (term1 + term2)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.1, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints,
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective,
       status = result$status, iterations = result$iterations)
}

# ============================================================================
# wrapper function for irt optimization
# ============================================================================

optimize_irt <- function(J, ell, r, quo, estimand = "gate") {
  # validate inputs
  if(ell >= J - 1) stop("ell must be less than J - 1")
  
  # when ell = 0, all estimands are equivalent
  if(ell == 0) {
    return(minimize_variance_bound_nloptr(J, ell, r, quo))
  }
  
  # choose optimization based on estimand
  if(estimand == "gate") {
    minimize_variance_bound_nloptr(J, ell, r, quo)
  } else if(estimand == "cte") {
    minimize_variance_cte(J, ell, r, quo)
  } else if(estimand == "both") {
    minimize_variance_gate_cte(J, ell, r, quo)
  } else {
    stop("Invalid estimand. Choose 'gate', 'cte', or 'both'")
  }
}

# ============================================================================
# wrapper function for crt optimization
# ============================================================================

optimize_crt <- function(design_type = "cross_section", J, ell, r, sigma_gamma_sq, 
                         sigma_epsilon_sq, sigma_upsilon_sq, K, sigma_phi_sq = NULL) {
  # validate inputs
  if(ell >= J - 1) stop("ell must be less than J - 1")
  
  if(design_type == "cross_section") {
    minimize_variance_crt_cross_section(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, K)
  } else if(design_type == "closed_cohort") {
    if(is.null(sigma_phi_sq)) stop("sigma_phi_sq is required for closed cohort design")
    minimize_variance_crt_closed_cohort(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, sigma_phi_sq, K)
  } else {
    stop("Invalid design_type. Choose 'cross_section' or 'closed_cohort'")
  }
}

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
                 
                 numericInput("J", 
                              label = HTML("Number of periods \\(J\\):"), 
                              value = 12, min = 3, max = 30, step = 1),
                 
                 numericInput("ell", 
                              label = HTML("Maximum duration of carryover effects \\(\\ell\\):"), 
                              value = 0, min = 0, max = 10, step = 1),
                 
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
                            radioButtons("estimand_crt", 
                                         label = "Target estimand:",
                                         choiceNames = list(
                                           HTML("\\(\\tau^{\\text{GATE}}\\)")
                                         ),
                                         choiceValues = list("gate"),
                                         selected = "gate",
                                         inline = FALSE),
                            
                            numericInput("J_crt", 
                                         label = HTML("Number of periods \\(J\\):"), 
                                         value = 12, min = 3, max = 30, step = 1),
                            
                            numericInput("ell_crt", 
                                         label = HTML("Maximum duration of carryover effects \\(\\ell\\):"), 
                                         value = 0, min = 0, max = 10, step = 1),
                            
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
                            radioButtons("estimand_crt_cc", 
                                         label = "Target estimand:",
                                         choiceNames = list(
                                           HTML("\\(\\tau^{\\text{GATE}}\\)")
                                         ),
                                         choiceValues = list("gate"),
                                         selected = "gate",
                                         inline = FALSE),
                            
                            numericInput("J_crt_cc", 
                                         label = HTML("Number of periods \\(J\\):"), 
                                         value = 12, min = 3, max = 30, step = 1),
                            
                            numericInput("ell_crt_cc", 
                                         label = HTML("Maximum duration of carryover effects \\(\\ell\\):"), 
                                         value = 0, min = 0, max = 10, step = 1),
                            
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