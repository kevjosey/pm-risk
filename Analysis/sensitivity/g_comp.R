library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), 
                         race = c("all", "white", "black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

### G-Computation Poisson Model

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GComp_All/'

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  w <- new_data$w
  x <- new_data$x
  
  # merge in ZIP-level covariates
  wx <- inner_join(w, x, by = c("zip", "year"))
  
  # remove collinear terms and identifiers
  if (scenario$dual == "both" & scenario$race == "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, id))
  } else if (scenario$dual == "both" & scenario$race != "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, race, id))
  } else if (scenario$dual != "both" & scenario$race == "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, id))
  } else if (scenario$dual != "both" & scenario$race != "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, race, id))
  }
  
  # outcome rates
  wx$ybar <- wx$dead/wx$time_count
  
  # estimate nuisance outcome model with glm + splines
  mumod <- glm(ybar ~ ns(a, 6) + . - a, weights = wx$time_count, model = FALSE,
               data = data.frame(ybar = wx$ybar, a = wx$pm25, w.tmp), family = quasipoisson())
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w.tmp)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # Separate Data into List
  mat.list <- split(cbind(wx$time_count, muhat.mat), wx$id)
  
  # Aggregate by ZIP-code-year
  mat <- do.call(rbind, lapply(mat.list, function(vec) {
    
    mat <- matrix(vec, ncol = length(a.vals) + 1)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
    
  } ))
  
  # Marginalize covariates
  mhat.vals <- colMeans(mat, na.rm = TRUE)
  
  print(paste0("Fit Complete: Scenario ", i))
  print(Sys.time())
  
  save(mumod, mhat.vals, a.vals, file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  rm(w, x, wx, mumod, muhat.mat, mat.list, mat); gc()
  
}