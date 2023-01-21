library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(gam)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(race = c("white","black", "asian", "hispanic", "other"),
                         dual = c(0, 1), sex = c("male", "female"),
                         age_break = c("[65,75)","[75,85)","[85,95)","[95,125)"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age_break <- as.character(scenarios$age_break)
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
  w.tmp <- setDF(subset(w.tmp, age_break == scenario$age_break))
  
  if (scenario$sex == "male") {
    w.tmp <- setDF(subset(w.tmp, female == 0))
  } else {
    w.tmp <- setDF(subset(w.tmp, female == 1))
  }
  
  wx.tmp <- inner_join(w.tmp, x.tmp, by = c("zip", "year"))
  
  y <- wx.tmp$dead
  a <- wx.tmp$pm25
  cal <- wx.tmp$cal
  ipw <- wx.tmp$ipw
  log.pop <- log(wx.tmp$time_count)
  
  w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count,
                                  race, dual, female, age_break, 
                                  id, ipw, cal))
  
  model_data <- gam_models(y = y, a = a, w = w, log.pop = log.pop, ipw = ipw, cal = cal, a.vals = a.vals)
  individual_data <- data.frame(wx.tmp)
  zip_data <- data.frame(x.tmp)
  
  dual_name <- ifelse(scenario$dual == 0, "high", "low")
  
  save(model_data, individual_data, zip_data, phat.vals,
       file = paste0(dir_mod, dual_name, "_", scenario$race, "_", 
                     scenario$sex, "_", scenario$age_break, ".RData"))
  
}
