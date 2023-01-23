library(parallel)
library(data.table)
library(tidyr)
library(dplyr)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

a.vals <- seq(2, 31, length.out = 146)

### Fit Weights

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
load(paste0(dir_data,"aggregate_data.RData"))

# initialize data

zip_cov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
             "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")
zip.tmp <- data.table(zip = aggregate_data$zip, year = aggregate_data$year,
                      model.matrix(~ ., data = aggregate_data[,zip_cov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]

a <- zip.tmp$pm25
x.tmp <- subset(zip.tmp, select = -c(zip, pm25))
x.tmp$year <- factor(x.tmp$year)

## LM GPS
pimod <- lm(a ~ ., data = data.frame(a = a, x.tmp))
pimod.vals <- c(pimod$fitted.values)
pimod.sd <- sigma(pimod)

# nonparametric density
a.std <- c(a - pimod.vals) / pimod.sd
dens <- density(a.std)
pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd

# get ipw numerator
pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  std <- c(a.tmp - pimod.vals) / pimod.sd
  approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
})

phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
phat <- predict(smooth.spline(a.vals, phat.vals), x = a)$y
phat[phat < 0] <- .Machine$double.eps
ipw <- phat/pihat

## Calibration Weights
x <- x.tmp %>% mutate_if(is.numeric, scale)
x.mat <- model.matrix(~ ., data = data.frame(x))
astar <- c(a - mean(a))/var(a)
astar2 <- c((a - mean(a))^2/var(a) - 1)
mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                 target = c(length(a), rep(0, ncol(x.mat) + 1)))
cal <- mod$weights

# truncation
trunc0 <- quantile(ipw, 0.005)
trunc1 <- quantile(ipw, 0.995)
ipw[ipw < trunc0] <- trunc0
ipw[ipw > trunc1] <- trunc1

trunc0 <- quantile(cal, 0.005)
trunc1 <- quantile(cal, 0.995)
cal[cal < trunc0] <- trunc0
cal[cal > trunc1] <- trunc1

x <- cbind(zip.tmp, ipw = ipw, cal = cal)

### Create Strata Data

create_strata <- function(data, x, phat.vals,
                          dual = c(0,1,2),
                          race = c("white", "black", "asian", 
                                          "hispanic", "other","all")) {
  
  if (dual == 0) {
    dual0 <- 0
  } else if (dual == 1) {
    dual0 <- 1
  } else if (dual == 2) {
    dual0 <- c(0,1)
  }
  
  if (race == "white") {
    race0 <- 1
  } else if (race == "black") {
    race0 <- 2
  } else if (race == "asian") {
    race0 <- 4
  } else if (race == "hispanic") {
    race0 <- 5
  } else if (race == "other") {
    race0 <- 3
  } else if (race == "all") {
    race0 <- c(1,2,3,4,5)
  }
  
  sub_data <- subset(data, race %in% race0 & dual %in% dual0)
  
  # Covariates and Outcomes
  w <- data.table(zip = sub_data$zip, year = sub_data$year, race = sub_data$race,
                  female = sub_data$female, dual = sub_data$dual, age_break = sub_data$age_break,
                  followup_year = sub_data$followup_year, dead = sub_data$dead, time_count = sub_data$time_count)[
                    ,lapply(.SD, sum), by = c("zip", "year", "race", "female", "dual", "age_break", "followup_year")]
  
  return(list(w = w, x = x, phat.vals = phat.vals))
  
}

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("white", "black", "asian", "hispanic", "other","all"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# format variables
aggregate_data$region <- factor(aggregate_data$region)
aggregate_data$zip <- factor(aggregate_data$zip)
aggregate_data$year <- factor(aggregate_data$year)
aggregate_data$female <- as.numeric(aggregate_data$female)
aggregate_data$race <- factor(aggregate_data$race)
aggregate_data$dual <- as.numeric(aggregate_data$dual)
aggregate_data$age_break <- factor(aggregate_data$age_break)
aggregate_data$followup_year <- factor(aggregate_data$followup_year)

x$zip <- factor(x$zip)
x$year <- factor(x$year)
x$id <- paste(x$zip, x$year, sep = "-")

# collate
lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data, x = x, phat.vals = phat.vals,
                            dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_data, "qd/", scenario$dual, "_", scenario$race, ".RData"))
  
})
