
# Mortality Benefits Associated with Decreasing PM<sub>2.5</sub>

Assessing the EPA's recent recommendation to reduce PM<sub>2.5</sub> from 12 &mu;g/m<sup>3</sup> to 10, 9, and 8 &mu;g/m<sup>3</sup>

###Order of Operations:

1. [`Analysis/data_process.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/data_process.R) – Creates stratified data sets from Xiao data.
2. [`Analysis/fit_weights.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_weights.R) – uses [`R/calibrate.R`](https://github.com/kevjosey/pm-risk/blob/main/R/calibrate.R) to fit nuisance weight models.
3. [`Analysis/fit_models.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_models.R) – uses [`R/gam_models.R`](https://github.com/kevjosey/pm-risk/blob/main/R/gam_models.R) to fit nuisance outcome models.
4. [`Analysis/fit_erf.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/fit_erf.R) – uses [`R/erf_models.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/erf_models.R) and the output from (2) and (3) to fit the exposure response functions.
5. [`Analysis/bootstrap.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/bootstrap.R) – yee ole bootstrap code to find standard errors of the ERFs from (2)-(4).
6. [`Analysis/results/plot.R`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/results/plot.R) – generates plots from models fit in (4)-(5).
7. [`Analysis/sensitivity`](https://github.com/kevjosey/pm-risk/blob/main/Analysis/sensitivity) – Code for fitting G-Computation and GPS as regressor sensitivity analyses. Also included is code to visualize covariate balance for various weighting methods. 

