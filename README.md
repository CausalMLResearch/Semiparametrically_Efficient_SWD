# Semiparametrically Efficient Stepped Wedge Designs
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Shiny App](https://img.shields.io/badge/Shiny-Interactive%20App-blue)](https://f07k8s-hao-wang.shinyapps.io/Semiparametrically_Efficient_SWD/)

This repository contains the R code to reproduce the results presented in [Semiparametrically Efficient Stepped Wedge Designs](TBD).

### Overview

An illustration of various designs, including $\mathbf{Z}_{\text{ff}}$, $\mathbf{Z}_{\text{ba}}$, $\mathbf{Z}_{\text{ffba}}$, $\mathbf{Z}_{\text{linear}}$, and $\mathbf{Z}_{\text{efficient}}$.

<img width="2239" height="631" alt="designs" src="https://github.com/user-attachments/assets/bc78d9fd-f133-49ac-b165-f910b78e2362" />

### Quickstart
To reproduce the results, please download this repo on a machine with R, run each R script in the [`codes`](codes) directory without modification, and then the results are saved in [`figures`](figures), [`results`](results), and [`tables`](tables). All the R scripts can be run standalone. To run the R scripts, you do not need to set any pathnames; everything is relative. 

Required R packages: doRNG, doSNOW, dplyr, foreach, ggpattern, ggplot2, glmmTMB, lme4, MASS, nloptr, parallel, patchwork, RColorBrewer, RSSthemes, tidyverse, xtable.

### Solve Optimal Design
Run our [`Shiny App`](https://f07k8s-hao-wang.shinyapps.io/Semiparametrically_Efficient_SWD/) to find optimal stepped wedge designs:
- **Choose trial type:** Individually randomized trial or cluster randomized trial (cross-sectional or closed-cohort)
- **Select target estimand(s):**
  - $\tau^{\text{GATE}}$ - Global average treatment effect
  - $\tau_s^{\text{CTE}}$ for $s = 0, ..., \ell-1$ - Cumulative treatment effect  
  - Both $\tau^{\text{GATE}}$ and $\tau_s^{\text{CTE}}$ for $s = 0, ..., \ell-1$
- **Specify trial settings:** Number of clusters and periods, maximum duration of carryover effects ($\ell$), correlation parameters, etc.
- **Run optimization:** Click "Run optimization" and download your optimal design

### Generate Illustrative Figures

#### A Motivating Example

- Run [`figure_1_sim.R`](codes/figure_1/figure_1_sim.R) and then run [`figure_1.py`](codes/figure_1/figure_1.py)
  - run a simulation for a simple two-period stepped wedge trial and compares two estimators: difference-in-means and AIPW
  - generate [`Figure 1`](figures/figure_1/2_period_mse.pdf)

#### Impact of Serial Correlations in Outcomes

- Run [`figure_1_sim.R`](codes/figure_1/figure_1_sim.R) and then run [`figure_1.py`](codes/figure_1/figure_1.py)
  - run a simulation for a simple two-period stepped wedge trial and compares two estimators: difference-in-means and AIPW
  - generate [`Figure 1`](figures/figure_1/2_period_mse.pdf)

#### Impact of Maximum Duration of Carryover Effects

#### Impact of Treatment vs. Control Outcome Variance

#### Impact of Target Causal Estimand

### Application to the Mobile Health Study

### Helper Functions
The following scripts in [`codes`](codes) contain helper functions used throughout the analysis. These are automatically sourced by the main scripts - you don't need to run them separately:
