
# Risk Assesment of Increasing PM_{2.5}

### Order of Operations

1. Analysis/data_process.R – Creates stratified data sets from Xiao data
2. Analysis/fit_models.R – uses R/gam_models.R to fit nuisance models.
3. Analysis/fit_erf.R – uses R/erf_models.R and the output from (2) to fit the exposure response functions.
4. Analysis/plots.R – generates plots from models fit in (3)
5. Analysis/jasa_comparison_plot.R – Xiao’s code comparing our methods.

Don’t worry about descriptives.R, this is for creating a standard Table 1.
