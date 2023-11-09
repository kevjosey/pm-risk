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
dir_out <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR/'

tea_8 <- tea_9 <- tea_10 <- tea_11 <- tea_12 <- data.frame()

# contrast indexes
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))

### Predict Total Events Avoided

## Hazard Ratio
for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  load(paste0(dir_out, scenario$dual, "_", scenario$race, "_boot.RData"))
  
  u.zip <- unique(individual_data$zip)
  m <- 1/log(length(u.zip)) # for m out of n bootstrap
  est_data$se.cal_trunc <- apply(boot_mat, 2, sd, na.rm = T)*sqrt(m)
  
  w <- data.table(id = individual_data$id, dead = individual_data$dead, time_count = individual_data$time_count)[
    ,lapply(.SD, sum), by = c("id")]
  
  events <- inner_join(subset(zip_data, select = c(id,pm25)), w, by = c("id"))
  
  tea_tmp_8 <- with(subset(events, pm25 > 8), sum(dead - time_count*as.numeric(est_data[idx8,6])))
  tea_tmp_9 <- with(subset(events, pm25 > 9), sum(dead - time_count*as.numeric(est_data[idx9,6])))
  tea_tmp_10 <- with(subset(events, pm25 > 10), sum(dead - time_count*as.numeric(est_data[idx10,6])))
  tea_tmp_11 <- with(subset(events, pm25 > 11), sum(dead - time_count*as.numeric(est_data[idx11,6])))
  tea_tmp_12 <- with(subset(events, pm25 > 12), sum(dead - time_count*as.numeric(est_data[idx12,6])))
  
  tea_var_8 <- with(subset(events, pm25 > 8), sum(m*var(boot_mat[,idx8])*time_count^2))
  tea_var_9 <- with(subset(events, pm25 > 9), sum(m*var(boot_mat[,idx9])*time_count^2))
  tea_var_10 <- with(subset(events, pm25 > 10), sum(m*var(boot_mat[,idx10])*time_count^2))
  tea_var_11 <- with(subset(events, pm25 > 11), sum(m*var(boot_mat[,idx11])*time_count^2))
  tea_var_12 <- with(subset(events, pm25 > 12), sum(m*var(boot_mat[,idx12])*time_count^2))
  
  tea_lower_8 <- tea_tmp_8 - 1.96*sqrt(tea_var_8)
  tea_lower_9 <- tea_tmp_9 - 1.96*sqrt(tea_var_9)
  tea_lower_10 <- tea_tmp_10 - 1.96*sqrt(tea_var_10)
  tea_lower_11 <- tea_tmp_11 - 1.96*sqrt(tea_var_11)
  tea_lower_12 <- tea_tmp_12 - 1.96*sqrt(tea_var_12)
  
  tea_upper_8 <- tea_tmp_8 + 1.96*sqrt(tea_var_8)
  tea_upper_9 <- tea_tmp_9 + 1.96*sqrt(tea_var_9)
  tea_upper_10 <- tea_tmp_10 + 1.96*sqrt(tea_var_10)
  tea_upper_11 <- tea_tmp_11 + 1.96*sqrt(tea_var_11)
  tea_upper_12 <- tea_tmp_12 + 1.96*sqrt(tea_var_12)
  
  n_8 <- with(subset(events, pm25 > 8), sum(dead))
  n_9 <- with(subset(events, pm25 > 9), sum(dead))
  n_10 <- with(subset(events, pm25 > 10), sum(dead))
  n_11 <- with(subset(events, pm25 > 11), sum(dead))
  n_12 <- with(subset(events, pm25 > 12), sum(dead))
  
  tea_8 <- rbind(tea_8, data.frame(dual = scenario$dual, race = scenario$race,
                                   tea_8 = tea_tmp_8, n_8 = n_8, prop_8 = tea_tmp_8/n_8,
                                   tea_lower_8 = tea_lower_8, tea_upper_8 = tea_upper_8))
  tea_9 <- rbind(tea_9, data.frame(dual = scenario$dual, race = scenario$race,
                                   tea_9 = tea_tmp_9, n_9 = n_9, prop_9 = tea_tmp_9/n_9, 
                                   tea_lower_9 = tea_lower_9, tea_upper_9 = tea_upper_9))
  tea_10 <- rbind(tea_10, data.frame(dual = scenario$dual, race = scenario$race,
                                     tea_10 = tea_tmp_10, n_10 = n_10, prop_10 = tea_tmp_10/n_10, 
                                     tea_lower_10 = tea_lower_10, tea_upper_10 = tea_upper_10))
  tea_11 <- rbind(tea_11, data.frame(dual = scenario$dual, race = scenario$race,
                                     tea_11 = tea_tmp_11, n_11 = n_11, prop_11 = tea_tmp_11/n_11, 
                                     tea_lower_11 = tea_lower_11, tea_upper_11 = tea_upper_11))
  tea_12 <- rbind(tea_12, data.frame(dual = scenario$dual, race = scenario$race,
                                     tea_12 = tea_tmp_12, n_12 = n_12, prop_12 = tea_tmp_12/n_12,
                                     tea_lower_12 = tea_lower_12, tea_upper_12 = tea_upper_12))
  
}

save(tea_8, tea_9, tea_10, tea_11, tea_12, file = '~/Data/tea.RData')