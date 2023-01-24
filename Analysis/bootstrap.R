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

### bootstrap

dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data_New'

a.vals <- seq(2, 31, length.out = 146)
boot.iter <- 500 # bootstrap iterations

# filenames
filenames <- list.files(dir_mod, full.names = TRUE)
fnames <- list.files(dir_mod, full.names = FALSE)

# scenarios
scenarios <- expand.grid(dual = c("both","high", "low"), race = c("all","white","black"))[-(2:3),]
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

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

# RUN IT!
for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_mod, scenario$dual, "_", scenario$race, ".RData"))
  
  boot_list <- mcapply(1:boot.iter, function(b, ...) {
    
    m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
    index <- sample(1:length(u.zip), m, replace = TRUE)  # initialize bootstrap  index
    
    ## GPS Model
    
    # zip-code data
    x <- bootstrap_data(data = zip_data, index = index, u.zip = u.zip)
    a <- x$pm25
    x.tmp <- subset(x, select = -c(zip, id, boot.id, pm25))
    x.tmp$year <- factor(x.tmp$year)
    
    # fit GPS model
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
    
    # truncation
    trunc0 <- quantile(ipw, 0.005)
    trunc1 <- quantile(ipw, 0.995)
    ipw[ipw < trunc0] <- trunc0
    ipw[ipw > trunc1] <- trunc1
    x$ipw <- ipw
    
    ## Outcome Model
    
    # individual level data
    w.tmp <- bootstrap_data(data = individual_data, index = index, u.zip = u.zip)
    wx <- inner_join(subset(w.tmp, select = -c(ipw, cal)),
                     data.frame(boot.id = x$boot.id, ipw = x$ipw), by = "boot.id")
    
    # factor strata variables
    wx$year <- factor(wx$year)
    wx$entry_age_break <- factor(wx$entry_age_break)
    wx$followup_year <- factor(wx$followup_year)
    
    # remove collinear terms and identifiers
    if (scenario$dual == "both" & scenario$race == "all") {
      wx$race <- factor(wx$race)
      w <- subset(wx, select = -c(zip, pm25, dead, time_count, id, boot.id, ipw))
    } else if (scenario$dual == "both" & scenario$race != "all") {
      w <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                  race, id, boot.id, ipw))
    } else if (scenario$dual != "both" & scenario$race == "all") {
      wx$race <- factor(wx$race)
      w <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                  dual, id, boot.id, ipw))
    } else if (scenario$dual != "both" & scenario$race != "all") {
      w <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                  dual, race, id, boot.id, ipw))
    }
    
    model_data <- gam_models_lm(y = wx$dead, a = wx$pm25, w = w,
                                log.pop = log(wx$time_count), 
                                ipw = wx$ipw, a.vals = a.vals)
    
    # set bandwidth from whole data
    target <- count_erf_lm(resid.lm = model_data$resid.lm, muhat.mat = model_data$muhat.mat,
                           log.pop = model_data$log.pop, w.id = wx$boot.id, a = x$pm25, x.id = x$boot.id,
                           bw = 1, a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)
    
    print(paste("Completed Scenario: ", i))
    return(target$estimate.lm)
    
  })
  
  boot_mat <- do.call(rbind, boot_list)
  colnames(boot_mat) <- a.vals
  
  # save output
  save(boot_mat, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_boot.RData"))
  rm(model_data, target, muhat.mat); gc()
  
}