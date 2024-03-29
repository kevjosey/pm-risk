library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(KernSmooth)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), 
                         race = c("all", "white", "black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)
bw.seq <- seq(0.1, 3, by = 0.1)
# bw <- c(2,2,2)

### Fit Exposure Responses from Pseudo Outcomes

# Load/Save models
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR/'

## Run Models

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_mod, scenario$dual, "_", scenario$race, ".RData"))
  
  # # leave-one-out cross validation
  if (i == 1) {
    
    # Separate Data into List
    mat.list <- with(model_data, split(cbind(exp(log.pop), resid.lm, resid.cal,
                                             resid.cal_trunc, muhat.mat), id))
    
    # Aggregate by ZIP-code-year
    mat <- do.call(rbind, lapply(mat.list, function(vec) {
      mat <- matrix(vec, ncol = length(a.vals) + 4)
      colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
    } ))
    
    mat.pool <- data.frame(id = names(mat.list), mat)
    mhat.vals <- colMeans(mat.pool[,-(1:4)], na.rm = TRUE)
    resid.dat <- inner_join(mat.pool[,(1:4)], data.frame(a = zip_data$pm25, id = zip_data$id), by = "id")
    resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
    
    # Pseudo-Outcomes
    resid.dat$psi.lm <- with(resid.dat, X1 + mhat)
    resid.dat$psi.cal <- with(resid.dat, X2 + mhat)
    resid.dat$psi.cal_trunc <- with(resid.dat, X3 + mhat)
    
    # grid search bandwidth
    risk.est.lm <- sapply(bw.seq, risk.fn, a.vals = a.vals,
                          psi = resid.dat$psi.lm, a = resid.dat$a)
    risk.est.cal <- sapply(bw.seq, risk.fn, a.vals = a.vals,
                           psi = resid.dat$psi.cal, a = resid.dat$a)
    risk.est.cal_trunc <- sapply(bw.seq, risk.fn, a.vals = a.vals,
                                 psi = resid.dat$psi.cal_trunc, a = resid.dat$a)
    bw <<- c(bw.seq[which.min(risk.est.lm)],
             bw.seq[which.min(risk.est.cal)],
             bw.seq[which.min(risk.est.cal_trunc)])
    
    rm(mat, mat.list, mat.pool, resid.dat); gc()
    
  }
  
  # fit exposure response curves
  target <- count_erf(resid.lm = model_data$resid.lm, resid.cal = model_data$resid.cal,
                      resid.cal_trunc = model_data$resid.cal_trunc, 
                      muhat.mat = model_data$muhat.mat, log.pop = model_data$log.pop, 
                      w.id = model_data$id, a = zip_data$pm25, x.id = zip_data$id,
                      bw = bw, a.vals = a.vals, phat.vals = phat.vals, se.fit = TRUE)
  
  # extract estimates
  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal),
                         estimate.cal_trunc = target$estimate.cal_trunc, se.cal_trunc = sqrt(target$variance.cal_trunc))
  
  print(paste0("Fit Complete: Scenario ", i))
  print(Sys.time())
  
  save(individual_data, zip_data, est_data,
       file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  rm(individual_data, zip_data, model_data, est_data, target); gc()
  
}