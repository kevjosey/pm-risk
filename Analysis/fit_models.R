library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gnm)
# library(gam)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(race = c("white","black", "asian", "hispanic", "other"), dual = c(0, 1))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

# Load/Save models
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data_New/'

for (i in 1:nrow(scenarios)) {
  
  print(i)
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  phat.vals <- new_data$phat.vals
  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(new_data$w)
  
  # merge in ZIP-level covariates
  wx.tmp <- inner_join(w.tmp, x.tmp, by = c("zip", "year"))
  
  # extract data components
  w.id <- wx.tmp$id
  y <- wx.tmp$dead
  a <- wx.tmp$pm25
  cal <- wx.tmp$cal
  ipw <- wx.tmp$ipw
  log.pop <- log(wx.tmp$time_count)
  
  wx.tmp$year <- factor(wx.tmp$year)
  wx.tmp$age_break <- factor(wx.tmp$age_break)
  wx.tmp$followup_year <- factor(wx.tmp$followup_year)
  
  # remove identifiers
  w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count,
                                  race, dual, id, ipw, cal))
  
  # fit gnm outcome model
  model_data <- gnm_models(y = y, a = a, w = w, log.pop = log.pop, ipw = ipw, cal = cal, a.vals = a.vals)
  individual_data <- data.frame(wx.tmp)
  zip_data <- data.frame(x.tmp)
  
  dual_name <- ifelse(scenario$dual == 0, "high", "low")
  
  # save data
  save(model_data, individual_data, zip_data, phat.vals,
       file = paste0(dir_mod, dual_name, "_", scenario$race, ".RData"))
  
}
