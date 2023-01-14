
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(ggplot2)
library(cobalt)

source('/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/pm-risk/R/erf.R')
source('/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/pm-risk/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(race = c("white","black", "asian", "hispanic", "other"),
                         dual = c(0, 1), sex = c("male", "female"),
                         age_break = c("[65,75)","[75,85)","[85,95)","[95,125)"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(4, 16, length.out = 241)

# Load/Save models
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/Output/DR_mod/'

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(new_data$w)
  w.tmp <- setDF(subset(w.tmp, age_break == scenario$age_break))
  
  if (scenario$sex == "male")
    w.tmp <- setDF(subset(w.tmp, sex == 0))
  else
    w.tmp <- setDF(subset(w.tmp, sex == 1))
  
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(wx.tmp$zip, wx.tmp$year, sep = "-")
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  y <- wx.tmp$dead
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count, age_break))
  
  model_data <- gam_models(y = y, a_w = a_w, w = w, w.id = w.id, log.pop = log.pop, 
                           a_x = a_x, x = x, x.id = x.id, a.vals = a.vals)
  
  individual_data <- data.frame(wx.tmp,
                                resid.lm = model_data$resid.lm,
                                weights.lm = model_data$weights.lm_w,
                                resid.cal = model_data$resid.cal,
                                weights.cal = model_data$weights.cal_w)
  
  zip_data <- data.frame(x.tmp, weights.lm = model_data$weights.lm_x, weights.cal = model_data$weights.cal_x)
  
  save(model_data, individual_data, zip_data, 
       file = paste0(dir_mod, scenario$dual, "_", scenario$race, "_", 
                     scenario$sex, "_", scenario$age_break, ".RData"))
  
}
