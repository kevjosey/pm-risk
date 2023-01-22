library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(gam)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

## bootstrap

dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data_New'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All_New/'

load(paste0(dir_data,"aggregate_data.RData"))
a.vals <- seq(2, 31, length.out = 146)
boot.iter <- 100 # bootstrap iterations

# scenarios
scenarios <- expand.grid(dual = c("", "high", "low"), race = c("", "white","black"))
scen_names <- expand.grid(dual = c("both","high", "low"), race = c("all","white","black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# ZIP Code Data
zip_cov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
             "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")
z_data <- data.table(zip = aggregate_data$zip, year = aggregate_data$year,
                     model.matrix(~ ., data = aggregate_data[,zip_cov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]
z_data$id <- paste(z_data$zip, z_data$year, sep = "-")
u.zip <- unique(z_data$zip)

# filenames
filenames <- list.files(dir_mod, full.names = TRUE)
fnames <- list.files(dir_mod, full.names = FALSE)

# function for getting bootstrap data
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
  sname <- scen_names[i,]
  
  boot_mat <- sapply(1:boot.iter, function(b, ...) {
    
    # m <- 2*sqrt(length(u.zip)) # for m out of n bootstrap 
    
    # initialize bootstrap data
    index <- sample(1:length(u.zip), length(u.zip), replace = TRUE)
    x <- bootstrap_data(data = z_data, index = index, u.zip = u.zip)
    a <- x$pm25
    x.tmp <- subset(x, select = -c(zip, id, boot.id, pm25))
    x.tmp$year <- factor(x.tmp$year)
    
    ## GPS Model
    
    # fit models
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
    
    grep1 <- grep(scenario[1], fnames)
    grep2 <- grep(scenario[2], fnames)
    idx <- 1:length(filenames)
    fn <- filenames[(idx %in% grep1) & (idx %in% grep2)]
    
    log.pop <- muhat.mat <- resid.lm <- w.id <- NULL
    
    # Fit strata-specific outcome models
    for (j in 1:length(fn)) {
      
      load(paste0(fn[j]))
      print(j)
      
      w.tmp <- bootstrap_data(data = individual_data, index = index, u.zip = u.zip)
      wx.tmp <- inner_join(subset(w.tmp, select = -c(ipw, cal, resid.lm, resid.cal)),
                           data.frame(boot.id = x$boot.id, ipw = x$ipw), by = "boot.id")
      
      # extract data components
      y <- wx.tmp$dead
      a <- wx.tmp$pm25
      ipw <- wx.tmp$ipw
      id <- wx.tmp$boot.id
      log.pop <- log(wx.tmp$time_count)
      
      # remove identifiers
      w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count,
                                      race, dual, female, age_break, 
                                      boot.id, id, ipw))
      
      model_data <- gam_models_lm(y = y, a = a, w = w, log.pop = log.pop, ipw = ipw, a.vals = a.vals)
      
      # concatenate
      log.pop <- c(log.pop, model_data$log.pop)
      muhat.mat <- rbind(muhat.mat, model_data$muhat.mat)
      resid.lm <- c(resid.lm, model_data$resid.lm)
      w.id <- c(w.id, id)
      
    }
    
    a <- x$pm25
    x.id <- x$boot.id
    
    # set bandwidth form whole data
    target <- count_erf_lm(resid.lm = resid.lm, muhat.mat = muhat.mat, log.pop = log.pop, w.id = w.id, 
                           a = a, x.id = x.id, bw = 1, a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)
    
    print(paste("Completed Scenario: ", i))
    return(target$estimate.lm)
    
  })
  
  # save output
  save(boot_mat, file = paste0(dir_out, sname$dual, "_", sname$race, "_boot.RData"))
  
}
