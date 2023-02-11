library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)

library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(KernSmooth)

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
boot.iter <- 1000 # bootstrap iterations

dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GPSReg_All/'

# function for getting cluster bootstrap data
# need to tweak to be more general
bootstrap_data <- function(data, index, u.zip) {
  
  n.zip <- length(u.zip)
  boot <- data.frame()
  
  aa <- u.zip[index]
  aa <- aa[which(aa %in% data$zip)]
  bb <- table(aa)
  
  for (j in 1:max(bb)) {
    
    cc <- data[data$zip %in% names(bb[bb == j]),]
    
    for (k in 1:j) {
      cc$boot.id <- paste(cc$id, k, sep = "-")
      boot <- rbind(boot, cc)
    }
    
  }
  
  return(boot)
  
}

### Generalized Propensity Score as a Regressor Bootstrap

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GPSReg_All/'

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  zip_data <- new_data$x
  individual_data <- new_data$w
  individual_data$id <- paste(individual_data$zip, individual_data$year, sep = "-")
  u.zip <- unique(individual_data$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  
  boot_list <- mclapply(1:boot.iter, function(b, ...) {
    
    index <- sample(1:length(u.zip), m, replace = TRUE)  # initialize bootstrap  index
    
    # zip-code data
    x <- bootstrap_data(data = zip_data, index = index, u.zip = u.zip)
    x.tmp <- subset(x, select = -c(zip, pm25, id, boot.id, ipw, cal, cal_trunc))
    x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
    
    # LM GPS mean
    pimod <- lm(a ~ ., data = data.frame(a = x$pm25, x.tmp))
    pimod.vals <- c(pimod$fitted.values)
    pimod.sd <- sigma(pimod)
    
    # nonparametric GPS
    a.std <- c(x$pm25 - pimod.vals) / pimod.sd
    dens <- density(a.std)
    x$gps <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
    
    # individual-level data
    w <- bootstrap_data(data = individual_data, index = index, u.zip = u.zip)
    wx <- inner_join(w, subset(x, select = -c(id, zip, year), by = "boot.id"))
    
    # select relevant individual level factors
    if (scenario$dual == "both" & scenario$race == "all") {
      w.tmp <- subset(wx, select = c(boot.id, female, entry_age_break, followup_year, year, dual, race))
    } else if (scenario$dual == "both" & scenario$race != "all") {
      w.tmp <- subset(wx, select = c(boot.id, female, entry_age_break, followup_year, year, dual))
    } else if (scenario$dual != "both" & scenario$race == "all") {
      w.tmp <- subset(wx, select = c(boot.id, female, entry_age_break, followup_year, year, race))
    } else if (scenario$dual != "both" & scenario$race != "all") {
      w.tmp <- subset(wx, select = c(boot.id, female, entry_age_break, followup_year, year))
    }
    
    # outcome rates
    wx$ybar <- wx$dead/wx$time_count
    
    # estimate nuisance outcome model with glm + splines
    mumod <- glm(ybar ~ ns(a, 6)*gps + . - a, weights = wx$time_count, model = FALSE,
                 data = data.frame(ybar = wx$ybar, a = wx$pm25, gps = wx$gps, w.tmp[,-1]), 
                 family = quasipoisson())
    
    # predictions along a.vals
    muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
      
      std <- c(a.tmp - pimod.vals) / pimod.sd
      gps.tmp <- approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
      wa.tmp <- inner_join(data.frame(a = a.tmp, gps = gps.tmp, boot.id = x$boot.id), w.tmp, by = "boot.id")
      predict(mumod, newdata = wa.tmp, type = "response")
      
    })
    
    # Separate Data into List
    mat.list <- split(cbind(wx$time_count, muhat.mat), wx$boot.id)
    
    # Aggregate by ZIP-code-year
    mat <- do.call(rbind, lapply(mat.list, function(vec) {
      
      mat <- matrix(vec, ncol = length(a.vals) + 1)
      colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
      
    } ))
    
    # Marginalize covariates
    mhat.vals <- colMeans(mat, na.rm = TRUE)
    return(mhat.vals)
    
  }, mc.cores = 8)
  
  boot_mat <- do.call(rbind, boot_list)
  colnames(boot_mat) <- a.vals
  
  # save output
  save(boot_mat, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_boot.RData"))
  
}
