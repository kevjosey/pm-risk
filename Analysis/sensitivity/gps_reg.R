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

### Generalized Propensity Score as a Regressor

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GPSReg_All/'

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  w <- new_data$w
  x <- new_data$x
  a <- x$pm25
  x.id <- x$id
  x.tmp <- subset(x, select = -c(zip, pm25, id,
                                 ipw, cal, cal_trunc))
  x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
  
  # LM GPS
  pimod <- lm(a ~ ., data = data.frame(a = a, x.tmp))
  pimod.vals <- c(pimod$fitted.values)
  pimod.sd <- sigma(pimod)
  
  # nonparametric GPS
  a.std <- c(a - pimod.vals) / pimod.sd
  dens <- density(a.std)
  x$gps <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  # merge in ZIP-level covariates
  wx <- inner_join(w, x, by = c("zip", "year"))
  
  # extract data components
  w.id <- wx$id
  y <- wx$dead
  a <- wx$pm25
  gps <- wx$gps
  log.pop <- log(wx$time_count)
  
  # remove collinear terms and identifiers
  if (scenario$dual == "both" & scenario$race == "all") {
    w.tmp <- subset(wx, select = c(id, female, entry_age_break, followup_year, year, dual, race))
  } else if (scenario$dual == "both" & scenario$race != "all") {
    w.tmp <- subset(wx, select = c(id, female, entry_age_break, followup_year, year, dual))
  } else if (scenario$dual != "both" & scenario$race == "all") {
    w.tmp <- subset(wx, select = c(id, female, entry_age_break, followup_year, year, race))
  } else if (scenario$dual != "both" & scenario$race != "all") {
    w.tmp <- subset(wx, select = c(id, female, entry_age_break, followup_year, year))
  }
  
  rm(w, x, wx); gc()
  
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm + splines
  mumod <- glm(ybar ~ ns(a, 6)*gps + . - a, weights = exp(log.pop), model = FALSE,
               data = data.frame(ybar = ybar, a = a, gps = gps, w.tmp[,-1]), family = quasipoisson())
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    std <- c(a.tmp - pimod.vals) / pimod.sd
    gps.tmp <- approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
    wa.tmp <- inner_join(data.frame(a = a.tmp, gps = gps.tmp, id = x.id), w.tmp, by = "id")
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
  
  save(pimod, mumod, mhat.vals, a.vals, file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  rm(mumod, muhat.mat, mat.list, mat); gc()
  
}