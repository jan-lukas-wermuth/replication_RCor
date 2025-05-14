
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication: “Statistical Inference for Rank Correlations”

<!-- badges: start -->
<!-- badges: end -->

This repository contains replication material for the paper “Statistical
Inference for Rank Correlations” by Marc-Oliver Pohle, Jan-Lukas Wermuth
and Christian H. Weiß. The corresponding `R` package is available in the
repository [RCor](https://github.com/jan-lukas-wermuth/RCor).

## Installation

In order to run the replication code, please install the `R` package
[RCor](https://github.com/jan-lukas-wermuth/RCor) first.

``` r
install.packages("devtools")
library(devtools)
install_github("jan-lukas-wermuth/RCor")
library(RCor)
```

It may be necessary to install further packages from
[CRAN](https://cran.r-project.org). All the necessary packages are
listed at the beginning of each `R` script.

## code

This folder contains the `R` scripts that are needed to calculate all
the results in the paper. I give a short overview over the several files
in the folder.

- True_Taus.R: This file simulates the values for $\tau$ for each
  continuous DGP that is described in the paper and a grid of dependence
  parameters $\alpha$. For the discrete DGPs, the attainable
  generalization $\gamma$ is computed. The $\alpha$-grid is chosen such
  that coverage plots (with alphas_CIs) and coverage tables (with
  alphas_CIs_short) can be computed. alphas_PowerGraph is needed for the
  power plots because their x-axis features $\gamma$ (which is
  equivalent to $\tau$ in the continuous case) and not $\alpha$, the
  parameter used in the simulations.

- True_Rhos.R: This script performs computations similar to True_Taus.R,
  where $\rho$ replaces $\tau$ and $\rho_b$ replaces $\gamma$. Since the
  power plots use $\gamma$ as their dependence measure, only alphas_CIs
  and alphas_CIs_short is needed.

- CIs_Coverage_Kendall.R: This script simulates empirical coverage rates
  of the confidence intervals for $\tau$ and $\gamma$ for the same DGPs
  as above.

- CIs_Coverage_Spearman.R: This script is the equivalent to
  CIs_Coverage_Kendall.R with $\rho$ instead of $\tau$ and $\rho_b$
  instead of $\gamma$.

- SizePower_IndepTest.R: This file simulates all the size and power
  values for all the independence tests in the paper. It uses the same
  DGPs as above.

- Applications.R: This file contains the data for both applications. It
  also includes the computation of all relevant correlation coefficients
  together with their confidence intervals and p-values for the
  uncorrelatedness and independence tests.

The files in the subfolder **functions** are called automatically by the
previous scripts and do not have to be run explicitly.

## results

### plots

This folder contains all plots.

### simulations

This folder contains all the raw simulation results in four subfolders:
**coverages**, **size_power**, **true_rhos** and **true_taus**. The
subfolder coverages is again splitted in two nested subfolders *kendall*
and *spearman*, where the former contains the coverage simulations for
$\tau$ and $\gamma$ and the latter the coverage simulations for $\rho$
and $\rho_b$.
