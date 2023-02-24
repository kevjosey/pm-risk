library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"),
                         race = c("all", "white", "black"))
# scenarios <- expand.grid(dual = c("high", "low"),
#                          race = c("white", "black", "asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

### Fit Outcome Models

# Load/Save models
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  w <- new_data$w
  x <- new_data$x
  phat.vals <- new_data$phat.vals
  
  # merge in ZIP-level covariates
  wx <- inner_join(w, x, by = c("zip", "year"))
  
  # remove collinear terms and identifiers
  if (scenario$dual == "both" & scenario$race == "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                    id, ipw, cal, cal_trunc))
  } else if (scenario$dual == "both" & scenario$race != "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, race, 
                                    id, ipw, cal, cal_trunc))
  } else if (scenario$dual != "both" & scenario$race == "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual,
                                    id, ipw, cal, cal_trunc))
  } else if (scenario$dual != "both" & scenario$race != "all") {
    w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, race,
                                    id, ipw, cal, cal_trunc))
  }
  
  # fit gam outcome model
  model_data <- gam_models(y = wx$dead, a = wx$pm25, w = w.tmp, log.pop = log(wx$time_count), id = wx$id, 
                           ipw = wx$ipw, cal = wx$cal, cal_trunc = wx$cal_trunc, a.vals = a.vals)
  individual_data <- data.frame(wx)
  zip_data <- data.frame(x)
  
  # check progress
  print(paste0("Fit Complete: Scenario ", i))
  print(Sys.time())
  
  # save data
  save(model_data, individual_data, zip_data, phat.vals,
       file = paste0(dir_mod, scenario$dual, "_", scenario$race, ".RData"))
  
}
