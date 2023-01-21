library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# scenarios
scenarios <- expand.grid(dual = c("both", "low", "high"), race = c("all","white", "black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All_New/'

tea <- data.frame()

# contrast indexes
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))

### Create Data

## Hazard Ratio
for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  # hazard ratio
  hr_tmp_8 <- c(as.numeric(est_data[,5])/as.numeric(est_data[idx8,5]))
  hr_tmp_9 <- c(as.numeric(est_data[,5])/as.numeric(est_data[idx9,5]))
  hr_tmp_10 <- c(as.numeric(est_data[,5])/as.numeric(est_data[idx10,5]))
  hr_tmp_11 <- c(as.numeric(est_data[,5])/as.numeric(est_data[idx11,5]))
  hr_tmp_12 <- c(as.numeric(est_data[,5])/as.numeric(est_data[idx12,5]))
  
  hr_var_8 <- (hr_tmp_8^2)*(c(est_data[,6]^2)/c(est_data[,5]^2) +
                              c(est_data[idx8,6]^2)/c(est_data[idx8,5]^2))
  hr_var_9 <- (hr_tmp_9^2)*(c(est_data[,6]^2)/c(est_data[,5]^2) +
                              c(est_data[idx9,6]^2)/c(est_data[idx9,5]^2))
  hr_var_10 <- (hr_tmp_10^2)*(c(est_data[,6]^2)/c(est_data[,5]^2) +
                                c(est_data[idx10,6]^2)/c(est_data[idx10,5]^2))
  hr_var_11 <- (hr_tmp_11^2)*(c(est_data[,6]^2)/c(est_data[,5]^2) +
                                c(est_data[idx11,6]^2)/c(est_data[idx11,5]^2))
  hr_var_12 <- (hr_tmp_12^2)*(c(est_data[,6]^2)/c(est_data[,5]^2) +
                                c(est_data[idx12,6]^2)/c(est_data[idx12,5]^2))
  
  m8 <- 1 - 1/hr_tmp_8
  m9 <- 1 - 1/hr_tmp_9
  m10 <- 1 - 1/hr_tmp_10
  m11 <- 1 - 1/hr_tmp_11
  m12 <- 1 - 1/hr_tmp_12
  
  v8 <- (hr_var_8)/(hr_tmp_8^4)
  v9 <- (hr_var_9)/(hr_tmp_9^4)
  v10 <- (hr_var_10)/(hr_tmp_10^4)
  v11 <- (hr_var_11)/(hr_tmp_11^4)
  v12 <- (hr_var_12)/(hr_tmp_12^4)
  
  tea_tmp <- data.frame(m8 = m8, m9 = m9, m10 = m10, m11 = m11, m12 = m12,
                        v8 = v8, v9 = v9, v10 = v10, v11 = v11, v12 = v12,
                        pm_closest = rep(est_data[,1], 5))
  
  idx <- sapply(zip_data$pm25, function(pm, a.vals) which.min(abs(a.vals - pm)), a.val = a.vals)
  x <- cbind(subset(zip_data, select = c(id, zip, year, pm25)), tea_tmp[idx,])
  
  w <- data.table(id = individual_data$id, dead = individual_data$dead, time_count = individual_data$time_count)[
    ,lapply(.SD, sum), by = c("id")]
  
  events <- inner_join(x, w, by = c("id"))
  
  tea_8 <- with(subset(events, pm25 > 8), sum(m8*dead))
  tea_9 <- with(subset(events, pm25 > 9), sum(m9*dead))
  tea_10 <- with(subset(events, pm25 > 10), sum(m10*dead))
  tea_11 <- with(subset(events, pm25 > 11), sum(m11*dead))
  tea_12 <- with(subset(events, pm25 > 12), sum(m12*dead))
  
  tea_var_8 <- with(subset(events, pm25 > 8), sum(v8*dead^2))
  tea_var_9 <- with(subset(events, pm25 > 9), sum(v9*dead^2))
  tea_var_10 <- with(subset(events, pm25 > 10), sum(v10*dead^2))
  tea_var_11 <- with(subset(events, pm25 > 11), sum(v11*dead^2))
  tea_var_12 <- with(subset(events, pm25 > 12), sum(v12*dead^2))
  
  tea_lower_8 <- tea_8 - 1.96*sqrt(tea_var_8)
  tea_lower_9 <- tea_9 - 1.96*sqrt(tea_var_9)
  tea_lower_10 <- tea_10 - 1.96*sqrt(tea_var_10)
  tea_lower_11 <- tea_11 - 1.96*sqrt(tea_var_11)
  tea_lower_12 <- tea_12 - 1.96*sqrt(tea_var_12)
  
  tea_upper_8 <- tea_8 + 1.96*sqrt(tea_var_8)
  tea_upper_9 <- tea_9 + 1.96*sqrt(tea_var_9)
  tea_upper_10 <- tea_10 + 1.96*sqrt(tea_var_10)
  tea_upper_11 <- tea_11 + 1.96*sqrt(tea_var_11)
  tea_upper_12 <- tea_12 + 1.96*sqrt(tea_var_12)
  
  n_8 <- with(subset(events, pm25 > 8), sum(dead))
  n_9 <- with(subset(events, pm25 > 9), sum(dead))
  n_10 <- with(subset(events, pm25 > 10), sum(dead))
  n_11 <- with(subset(events, pm25 > 11), sum(dead))
  n_12 <- with(subset(events, pm25 > 12), sum(dead))
  
  tea <- rbind(tea, data.frame(dual = scenario$dual, race = scenario$race,
                               tea_8 = tea_8, prop_8 = tea_8/n_8, tea_lower_8 = tea_lower_8, tea_upper_8 = tea_upper_8,
                               tea_9 = tea_9, prop_9 = tea_9/n_9, tea_lower_9 = tea_lower_9, tea_upper_9 = tea_upper_9,
                               tea_10 = tea_10, prop_10 = tea_10/n_10, tea_lower_10 = tea_lower_10, tea_upper_10 = tea_upper_10,
                               tea_11 = tea_11, prop_11 = tea_11/n_11, tea_lower_11 = tea_lower_11, tea_upper_11 = tea_upper_11,
                               tea_12 = tea_12, prop_12 = tea_12/n_12, tea_lower_12 = tea_lower_12, tea_upper_12 = tea_upper_12))
  
}

save(tea, file = '~/Data/tea.RData')