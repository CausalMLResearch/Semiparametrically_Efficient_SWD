# Semiparametrically Efficient Stepped Wedge Designs
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Shiny App](https://img.shields.io/badge/Shiny-Interactive%20App-blue)](https://f07k8s-hao-wang.shinyapps.io/Semiparametrically_Efficient_SWD/)

This repository contains the R code to reproduce the results presented in [Semiparametrically Efficient Stepped Wedge Designs](TBD).

### Overview

An illustration of various designs, including ![Z_ff](https://latex.codecogs.com/svg.image?\mathbf{Z}_{\text{ff}}), ![Z_ba](https://latex.codecogs.com/svg.image?\mathbf{Z}_{\text{ba}}), ![Z_ffba](https://latex.codecogs.com/svg.image?\mathbf{Z}_{\text{ffba}}), ![Z_linear](https://latex.codecogs.com/svg.image?\mathbf{Z}_{\text{linear}}), and ![Z_efficient](https://latex.codecogs.com/svg.image?\mathbf{Z}_{\text{efficient}}).

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

- Run [`figure_2a_J2.R`](codes/figure_2/figure_2a_J2.R), [`figure_2a_J4.R`](codes/figure_2/figure_2a_J4.R), and [`figure_2a_J10.R`](codes/figure_2/figure_2a_J10.R)
  - optimal designs for $\tau_0$ when $J = 2, 4, 10$, $r \in [0, 1]$, $\ell = 0$, and $\sigma_\upsilon^2 = 0$
  - generate [`J_2`](figures/figure_2/J_2.pdf), [`J_4`](figures/figure_2/J_4.pdf), and [`treated_fraction_r.pdf`](figures/figure_2/treated_fraction_r.pdf)

#### Impact of Maximum Duration of Carryover Effects

- Run [`figure_2b.R`](codes/figure_2/figure_2b.R)
  - optimal designs for $\tau_0$ for $\ell = 0$ and $\tau^{\text{GATE}}$ for $\ell = 1, 2, 3$ when $J = 10$, $r = 0.9$, and $\sigma_\upsilon^2 = 0$
  - generate [`treated_fraction_ell.pdf`](figures/figure_2/treated_fraction_ell.pdf)

#### Impact of Treatment vs. Control Outcome Variance

- Run [`figure_2c.R`](codes/figure_2/figure_2c.R)
  - optimal designs for $\tau^{\text{GATE}}$ when $J = 10$, $r = 0.9$, $\ell = 0$, and $\sigma_\upsilon^2/\sigma_\gamma^2 = 0, 1, 2$
  - generate [`treated_fraction_upsilon.pdf`](figures/figure_2/treated_fraction_upsilon.pdf)

#### Impact of Target Causal Estimand

- Run [`figure_2c_J2.R`](codes/figure_2/figure_2c_J2.R) and [`figure_2c_J10.R`](codes/figure_2/figure_2c_J10.R)
  - optimal designs for $\tau^{\text{GATE}}$, $\tau_0^{\text{CTE}}$, and both when $J = 2, 10$, $\ell = 1$, and $\sigma_\upsilon^2 = 0$
  - generate [`multiple_estimands.pdf`](figures/figure_2/multiple_estimands.pdf) and [`treated_fraction_vary_estimand.pdf`](figures/figure_2/treated_fraction_vary_estimand.pdf)

#### Different Designs

- Run [`figure_3.R`](codes/figure_3/figure_3.R)
  - benchmark designs and the efficient design
  - generate [`designs.pdf`](figures/figure_3/designs.pdf)
 
### Simulation Study

#### Base Case

- Run [`figure_4a_r_high_ell_1_sim`](codes/figure_4/figure_4a_r_high_ell_1_sim.R) and then run [`generate_tables.R`](codes/figure_4/generate_tables.R)
  - performance of AIPW, IPW, and OR under different designs when $I = 32$, $J = 10$, $\ell = 1$, $\tau_0^{\text{CTE}} = 1$ and $\tau^{\text{GATE}} = 2$
  - generate [`table_1.tex`](tables/table_1.tex)

#### Varying Maximum Duration of Carryover Effects

- Run [`figure_4a_r_high_ell_0_sim`](codes/figure_4/figure_4a_r_high_ell_0_sim.R) and then run [`figure_4a_r_high_ell_0.R`](codes/figure_4/figure_4a_r_high_ell_0.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $\ell = 0$, $\tau_0 = 2$
  - generate [`vary_ell_0.pdf`](figures/figure_4/vary_ell_0.pdf)
 
- Run [`figure_4a_r_high_ell_1_sim`](codes/figure_4/figure_4a_r_high_ell_1_sim.R) and then run [`figure_4a_r_high_ell_1.R`](codes/figure_4/figure_4a_r_high_ell_1.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $\ell = 1$, $\tau_0^{\text{CTE}} = 1$, $\tau^{\text{GATE}} = 2$
  - generate [`vary_ell_1.pdf`](figures/figure_4/vary_ell_1.pdf)
 
- Run [`figure_4a_r_high_ell_2_sim`](codes/figure_4/figure_4a_r_high_ell_2_sim.R) and then run [`figure_4a_r_high_ell_2.R`](codes/figure_4/figure_4a_r_high_ell_2.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $\ell = 2$, $\tau_0^{\text{CTE}} = 0.5$, $\tau_0^{\text{CTE}} = 1$, $\tau^{\text{GATE}} = 2$
  - generate [`vary_ell_2.pdf`](figures/figure_4/vary_ell_2.pdf)

#### Varying Serial Correlations

#### Varying Experimental Durations

#### Impact of Target Causal Estimands

### Application to the Mobile Health Study

#### Patient Cohort

#### Caregiver Cohort

#### A Dyad of Patient and Caregiver

### Helper Functions
The following scripts in [`codes`](codes) contain helper functions used throughout the analysis. These are automatically sourced by the main scripts - you don't need to run them separately:
