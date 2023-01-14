
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/erf.R')
source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all", "white","black"), sex = c(0, 1, 2),
                         age_break = c("all","[65,75)","[75,85)","[85,95)","[95,125)"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(4, 16, length.out = 241)
n.boot <- 1000

# Load/Save models
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_mod = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_mod/'

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(new_data$w)

  if (scenario$age_break != "all")
    w.tmp <- setDF(subset(w.tmp, age_break == scenario$age_break))
  if (scenario$sex != 2)
    w.tmp <- setDF(subset(w.tmp, sex == scenario$sex))
  
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
