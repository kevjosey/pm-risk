library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

a.vals <- seq(2, 31, length.out = 146)

# scenarios
scenarios <- expand.grid(dual = c("high", "low", "both"),
                         race = c("white", "black", "all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

### G-Computation Poisson Model

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GComp_All/'
load(paste0(dir_data,"aggregate_data.RData"))

for (i in 1:nrow(scnearios)) {
  
  scenario <- scenarios[i,]
  
  if (scenario$dual == "high") {
    dual0 <- 0
  } else if (scenario$dual == "low") {
    dual0 <- 1
  } else if (scenario$dual == "both") {
    dual0 <- c(0,1)
  }
  
  if (scenario$race == "white") {
    race0 <- 1
  } else if (scenario$race == "black") {
    race0 <- 2
  } else if (scenario$race == "asian") {
    race0 <- 4
  } else if (scenario$race == "hispanic") {
    race0 <- 5
  } else if (scenario$race == "other") {
    race0 <- 3
  } else if (scenario$race == "all") {
    race0 <- c(1,2,3,4,5)
  }
  
  sub_data <- subset(aggregate_data, race %in% race0 & dual %in% dual0)
  
  # Individual Data
  w <- data.table(zip = sub_data$zip, year = sub_data$year, race = sub_data$race,
                  female = sub_data$female, dual = sub_data$dual, entry_age_break = sub_data$entry_age_break,
                  followup_year = sub_data$followup_year, dead = sub_data$dead, time_count = sub_data$time_count)[
                    ,lapply(.SD, sum), by = c("zip", "year", "race", "female", "dual", "entry_age_break", "followup_year")]
  
  w$zip <- factor(w$zip)
  w$year <- factor(w$year)
  w$female <- as.numeric(w$female)
  w$race <- factor(w$race)
  w$dual <- as.numeric(w$dual)
  w$entry_age_break <- factor(w$entry_age_break)
  w$followup_year <- factor(w$followup_year)
  
  rm(sub_data); gc()
  
  # Neighborhood Data
  zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
            "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")
  x <- data.table(zip = aggregate_data$zip, year = aggregate_data$year,
                  model.matrix(~ ., data = aggregate_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]
  
  x$zip <- factor(x$zip)
  x$year <- factor(x$year)
  x$id <- paste(x$zip, x$year, sep = "-")
  
  # merge in ZIP-level covariates
  wx <- inner_join(w, x, by = c("zip", "year"))
  
  # extract data components
  y <- wx$dead
  a <- wx$pm25
  w.id <- wx$id
  log.pop <- log(wx$time_count)
  
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
  
  rm(w, x, wx); gc()
  
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm + splines
  mumod <- glm(ybar ~ ns(a, 6) + . - a, weights = exp(log.pop), model = FALSE,
               data = data.frame(ybar = ybar, a = a, w.tmp), family = quasipoisson())
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w.tmp)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # Separate Data into List
  mat.list <- split(cbind(exp(log.pop), muhat.mat), w.id)
  
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
  
  rm(mumod, muhat.mat, mat.list, mat); gc()
  
}