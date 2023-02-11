
# Mortality Benefits Associated with Decreasing PM<sub>2.5</sub>

Assessing the EPA's recent recommendation to reduce PM<sub>2.5</sub> from 12 &mu;g/m<sup>3</sup> to 10, 9, and 8 &mu;g/m<sup>3</sup>.

## Contents

### Analysis

The scripts described below are presented in the order with which they should be run.

1. [`Analysis/data_process.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/data_process.R): Creates stratified data sets from [`National Causal`](https://github.com/NSAPH/National-Causal-Analysis) data.
2. [`Analysis/fit_weights.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_weights.R): Applies functions from [`R/calibrate.R`](https://github.com/kevjosey/pm-risk/blob/main/R/calibrate.R) to fit nuisance inverse-probability/balancing weight models.
3. [`Analysis/fit_models.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_models.R): Applies functions from [`R/gam_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/gam_models.R) to fit nuisance outcome models.
4. [`Analysis/fit_erf.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_erf.R): Uses code from [`R/erf_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/erf_models.R) and the output from (2) and (3) to fit the exposure response functions.
5. [`Analysis/fit_boot.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_boot.R): Script for generating bootstrap samples used to estimate standard errors of the ERCs from (2)-(4).
6. [`Analysis/results`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/results/plot.R): Contains scripts for generating plots and tables from data generated in (1) and the models fit in (4)-(5).
7. [`Analysis/sensitivity`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/sensitivity): Code for fitting G-Computation and GPS as a regressor sensitivity analyses. Also included is code to visualize covariate balance for various weighting methods. 

### R

- [`R/calibrate.R`](https://github.com/kevjosey/pm-risk/blob/main/R/calibrate.R): Generic calibration function for fitting balancing weights. Implemented in [`Analysis/fit_weights.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_weights.R).</li>
- [`R/gam_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/gam_models.R): Wrapper functions for fitting spline outcome models (formerly generalized additive outcome models) assuming a quasi-Poisson likelihood. Functions output components used to construct the doubly-robust pseudo-outcome applied in [`R/erf_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/erf_models.R). Implemented in [`Analysis/fit_models.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_models.R).</li>
- [`R/erf_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/erf_models.R): Nonparametric doubly-robust estimator of the exposure response curve using the pseudo-outcomes fit with components found in [`R/calibrate.R`](https://github.com/kevjosey/pm-risk/blob/main/R/calibrate.R) and [`R/gam_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/gam_models.R). Implemented in [`Analysis/fit_erf.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_erf.R).