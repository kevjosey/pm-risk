
# Risk Assesment of Increasing PM_{2.5}

### Code Order of Operations:

1. [`Analysis/data_process.R`] – Creates stratified data sets from Xiao data.
2. [`Analysis/fit_weights.R`] – uses [`R/calibrate.R`] to fit nuisance weight models.
3. [`Analysis/fit_models.R`] – uses [`R/gam_models.R`] to fit nuisance outcome models.
4. [`Analysis/fit_erf.R`] – uses [`R/erf_models.R`] and the output from (2) and (3) to fit the exposure response functions.
5. [`Analysis/bootstrap.R`] – yee ole bootstrap code to find standard errors of the ERFs from (2)-(4).
6. [`Analysis/results/plots.R`] – generates plots from models fit in (4).
7. [`Analysis/results/compare_plot.R`] – Code to compare Xiao's output with our methods.
8. [`Analysis/sensitivity.R`] – A place to talk about our feelings. Just kidding, there is no love in statistics. Only more analyses to prove to yourself that you are not an idiot!

Don’t worry about the other code yet. I still need to clean some of it up, this is for creating a standard Table 1.
