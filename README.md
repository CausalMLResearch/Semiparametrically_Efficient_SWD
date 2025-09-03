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
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $r = 0.9$, $\ell = 0$, $\tau_0 = 2$
  - generate [`vary_ell_0.pdf`](figures/figure_4/vary_ell_0.pdf)
 
- Run [`figure_4a_r_high_ell_1_sim`](codes/figure_4/figure_4a_r_high_ell_1_sim.R) and then run [`figure_4a_r_high_ell_1.R`](codes/figure_4/figure_4a_r_high_ell_1.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $r = 0.9$, $\ell = 1$, $\tau_0^{\text{CTE}} = 1$, $\tau^{\text{GATE}} = 2$
  - generate [`vary_ell_1.pdf`](figures/figure_4/vary_ell_1.pdf)
 
- Run [`figure_4a_r_high_ell_2_sim`](codes/figure_4/figure_4a_r_high_ell_2_sim.R) and then run [`figure_4a_r_high_ell_2.R`](codes/figure_4/figure_4a_r_high_ell_2.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $r = 0.9$, $\ell = 2$, $\tau_0^{\text{CTE}} = 0.5$, $\tau_0^{\text{CTE}} = 1$, $\tau^{\text{GATE}} = 2$
  - generate [`vary_ell_2.pdf`](figures/figure_4/vary_ell_2.pdf)

#### Varying Serial Correlations

- Run [`figure_4b_r_high_10_sim.R`](codes/figure_4/figure_4b_r_high_10_sim.R) and then run [`figure_4b_r_high_10.R`](codes/figure_4/figure_4b_r_high_10.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $r = 0.9$, $\ell = 0$, $\tau_0 = 2$
  - generate [`vary_r_0.9.pdf`](figures/figure_4/vary_r_0.9.pdf)
 
- Run [`figure_4b_r_med_10_sim.R`](codes/figure_4/figure_4b_r_med_10_sim.R) and then run [`figure_4b_r_med_10.R`](codes/figure_4/figure_4b_r_med_10.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $r = 0.5$, $\ell = 0$, $\tau_0 = 2$
  - generate [`vary_r_0.5.pdf`](figures/figure_4/vary_r_0.5.pdf)

#### Varying Experimental Durations

- Run [`figure_4b_r_high_8_sim.R`](codes/figure_4/figure_4b_r_high_8_sim.R), [`figure_4b_r_high_10_sim.R`](codes/figure_4/figure_4b_r_high_10_sim.R), [`figure_4b_r_high_12_sim.R`](codes/figure_4/figure_4b_r_high_12_sim.R), [`figure_4b_r_high_14_sim.R`](codes/figure_4/figure_4b_r_high_14_sim.R), [`figure_4b_r_high_16_sim.R`](codes/figure_4/figure_4b_r_high_16_sim.R), and then run [`figure_4c.R`](codes/figure_4/figure_4c.R)
  - performance of AIPW under different designs when $I = 32$, $J = 8, 10, 12, 14, 16$, $r = 0.9$, $\ell = 0$, $\tau_0 = 2$
  - generate [`vary_J.pdf`](figures/figure_4/vary_J.pdf)

#### Impact of Target Causal Estimands

- Run [`figure_4d_sim.R`](codes/figure_4/figure_4d_sim.R) and then run [`figure_4d.R`](codes/figure_4/figure_4d.R)
  - performance of AIPW under different designs when $I = 32$, $J = 10$, $r = 0.9$, $\ell = 1$, $\tau_0^{\text{CTE}} = 1$, $\tau^{\text{GATE}} = 2$
  - generate [`multiple_cte.pdf`](figures/figure_4/multiple_cte.pdf), [`multiple_gate.pdf`](figures/figure_4/multiple_gate.pdf), and [`multiple_both.pdf`](figures/figure_4/multiple_both.pdf)
 
#### Supplementary Materials

- Run scripts in [`supp_figure_1`](codes/supp_figure_1) to replicate all the simulations above by replacing $I = 32$ with $I = 1000$
 
### Extension to Cluster Randomized Trials

- Run [`supp_table_6_unequal_ell_0.R`](codes/supp_table_6/supp_table_6_unequal_ell_0.R), [`supp_table_6_unequal_ell_1.R`](codes/supp_table_6/supp_table_6_unequal_ell_1.R) and then run [`generate_supp_tables.R`](codes/supp_figure_1/generate_supp_tables.R)
  - performance of AIPW under different designs for SW-CRTs with unequal cluster sizes
  - generate [`supp_table_6_ell_0.tex`](tables/supp_table_6_ell_0.tex) and [`supp_table_6_ell_1.tex`](tables/supp_table_6_ell_1.tex)

### Application to the Mobile Health Study

#### Correlation Parameters Estimation

- Run [`patient_para_16.R`](codes/supp_figure_2/patient_para_16.R), [`caregiver_para_16.R`](codes/supp_figure_2/caregiver_para_16.R) and then run [`supp_figure_2.R`](codes/supp_figure_2/supp_figure_2.R)
  - estimated correlation parameter in the model health study when $I = 30$ and $H_{\text{hist}} = 16$
  - generate [`mobile_parameters.pdf`](figures/supp_figure_2/mobile_parameters.pdf)

##### Patient Cohort or Caregiver Cohort

- Run [`caregiver_ell_0_8_sim.R`](codes/figure_5/caregiver_ell_0_8_sim.R), [`caregiver_ell_0_10_sim.R`](codes/figure_5/caregiver_ell_0_10_sim.R), [`caregiver_ell_0_12_sim.R`](codes/figure_5/caregiver_ell_0_12_sim.R), [`caregiver_ell_0_14_sim.R`](codes/figure_5/caregiver_ell_0_14_sim.R), [`caregiver_ell_0_16_sim.R`](codes/figure_5/caregiver_ell_0_16_sim.R), [`caregiver_ell_1_8_sim.R`](codes/figure_5/caregiver_ell_1_8_sim.R), [`caregiver_ell_1_10_sim.R`](codes/figure_5/caregiver_ell_1_10_sim.R), [`caregiver_ell_1_12_sim.R`](codes/figure_5/caregiver_ell_1_12_sim.R), [`caregiver_ell_1_14_sim.R`](codes/figure_5/caregiver_ell_1_14_sim.R), [`caregiver_ell_1_16_sim.R`](codes/figure_5/caregiver_ell_1_16_sim.R), and then run [`figure_5.R`](codes/figure_5/figure_5.R)
  - performance of AIPW under different designs when $I = 40$, $J = 8, 10, 12, 14, 16$, $\ell = 0, 1$ for patient cohort and caregiver cohort in the mobile health study
  - generate [`mobile.pdf`](figures/figure_5/mobile.pdf)

#### A Dyad of Patient and Caregiver

- Run [`cluster_ell_0_8_sim.R`](codes/supp_figure_3/cluster_ell_0_8_sim.R), [`cluster_ell_0_10_sim.R`](codes/supp_figure_3/cluster_ell_0_10_sim.R), [`cluster_ell_0_12_sim.R`](codes/supp_figure_3/cluster_ell_0_12_sim.R), [`cluster_ell_0_14_sim.R`](codes/supp_figure_3/cluster_ell_0_14_sim.R), [`cluster_ell_0_16_sim.R`](codes/supp_figure_3/cluster_ell_0_16_sim.R) and then run [`supp_figure_3.R`](codes/supp_figure_3/supp_figure_3.R)
  - performance of AIPW under different designs when $I = 40$, $J = 8, 10, 12, 14, 16$, $\ell = 0, 1$ for dyads of patient cohort and caregiver cohort in the mobile health study
  - generate [`mobile_cluster.pdf`](figures/supp_figure_3/mobile_cluster.pdf)

### Helper Functions
The following scripts in [`codes`](codes) contain helper functions used throughout the analysis. These are automatically sourced by the main scripts - you don't need to run them separately:
