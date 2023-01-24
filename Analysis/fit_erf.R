library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(KernSmooth)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), race = c("all", "white", "black"))
# scenarios <- expand.grid(dual = c("high", "low"), race = c("asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

# Load/Save models
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data_New/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All_New/'

filenames <- list.files(dir_mod, full.names = TRUE)
fnames <- list.files(dir_mod, full.names = FALSE)

## Run Models

for (i in c(6)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_mod, scenario$dual, "_", scenario$race, ".RData"))
  
  # 10-fold CV to find bandwidth
  if (i == 1) {
    
    # estimate effect curve with optimal bandwidth
    est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
    
    wts <- do.call(c, lapply(split(exp(model_data$log.pop),
                                   individual_data$id), sum))
    mat.list <- split(cbind(exp(log.pop), model_data$resid.lm,
                            model_data$muhat.mat), individual_data$id)
    
    # Aggregate by ZIP-code-year
    agg <- do.call(rbind, lapply(mat.list, function(vec) {
      mat <- matrix(vec, ncol = length(a.vals) + 2)
      colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
    } ))
    
    agg.new <- data.frame(wts = wts, id = names(mat.list), agg)
    mhat.vals <- colMeans(agg.new[,-c(1:3)], na.rm = TRUE)
    resid.dat <- inner_join(agg.new[,1:3], data.frame(a = zip_data$pm25, id = zip_data$id), by = "id")
    resid.dat <- resid.dat[sample(1:nrow(resid.dat), 10000, replace = FALSE),] # subset for speed
    
    resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
    
    # Pseudo-Outcome
    resid.dat$psi.lm <- with(resid.dat, X1 + mhat)
    
    # cross validation for bandwidth
    risk.est <- sapply(seq(0.2, 5, by = 0.2), risk.fn, a.vals = a.vals,
                       psi = resid.dat$psi.lm, a = resid.dat$a)
    bw <<- bw.seq[which.min(risk.est)]
    
    rm(wts, mat.list, agg, agg.new, resid.dat); gc()
    
  }
  
  # fit exposure response curves
  target <- count_erf(resid.lm = model_data$resid.lm, resid.cal = model_data$resid.cal, 
                      muhat.mat = model_data$muhat.mat, log.pop = model_data$log.pop,
                      w.id = individual_data$id, a = zip_data$pm25, x.id = zip_data$id, 
                      bw = bw, a.vals = a.vals, phat.vals = phat.vals, se.fit = TRUE)
  
  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal))
  
  extra <- list(lm = target$fit.lm, spl.lm = target$spl.lm,
                cal = target$fit.cal, spl.cal = target$spl.cal)
  
  print(paste0("Fit Complete: Scenario ", i))
  save(individual_data, zip_data, est_data, extra,
       file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  rm(individual_data, zip_data, model_data, est_data, target, extra); gc()
  
}
