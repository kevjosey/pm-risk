
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(ggplot2)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c("", "high", "low"), race = c("", "white","black"))
scen_names <- expand.grid(dual = c("both","high", "low"), race = c("all","white","black"))
# scenarios <- expand.grid(dual = c("high", "low"), race = c("asian","hispanic", "other"))
# scen_names <- expand.grid(dual = c("high", "low"), race = c("asian","hispanic","other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals.old <- seq(0.00783038, 30.92493, length.out = 201)
a.vals <- a.vals.old[c(T,F)]
rm_idx <- which(a.vals.old %in% a.vals)

# Load/Save models
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All/'

filenames <- list.files(dir_mod, full.names = TRUE)
fnames <- list.files(dir_mod, full.names = FALSE)

## Run Models

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  sname <- scen_names[i,]
  
  grep1 <- grep(scenario[1], fnames)
  grep2 <- grep(scenario[2], fnames)
  idx <- 1:length(filenames)
  
  fn <- filenames[(idx %in% grep1) & (idx %in% grep2)]
  
  w.id <- log.pop <- nval <- NULL
  muhat.mat <- phat.tmp <- NULL
  resid.lm <- resid.cal <- NULL
  i_data <- z_data <- NULL
  
  for (j in 1:length(fn)) {
    
    load(paste0(fn[j]))
    print(j)
    
    w.id <- c(w.id, model_data$w.id)
    log.pop <- c(log.pop, model_data$log.pop)
    nval <- c(nval, sum(exp(model_data$log.pop)))
    
    muhat.mat <- rbind(muhat.mat, model_data$muhat.mat[,rm_idx])
    phat.tmp <- rbind(phat.tmp, model_data$phat.vals[rm_idx])
    
    resid.lm <- c(resid.lm, model_data$resid.lm)
    resid.cal <- c(resid.cal, model_data$resid.cal)
    
    i_data <- rbind(i_data, individual_data)
    z_data <- rbind(z_data, zip_data)
    
  }
  
  # summary data
  phat.vals <- apply(phat.tmp, 2, weighted.mean, w = nval)
  zip_data <- z_data[!duplicated(paste(z_data$zip, z_data$year, sep = "-")),]
  individual_data <- i_data
  
  rm(z_data, i_data, phat.tmp, nval); gc()
  
  x.id <- paste(zip_data$zip, zip_data$year, sep = "-")
  a_x <- zip_data$pm25
  
  # 10-fold CV to find bandwidth
  # if (i == 1) {
  # 
  #   wts <- do.call(c, lapply(split(exp(log.pop), w.id), sum))
  #   list.cal <- split(data.frame(psi = psi.cal, wts = exp(log.pop)) , w.id)
  #   psi.cal.new <- data.frame(psi = do.call(c, lapply(list.cal, function(df) sum(df$psi*df$wts)/sum(df$wts))),
  #                             wts = wts, id = names(list.cal))
  #   cal.dat <- inner_join(psi.cal.new, data.frame(a = a_x, id = x.id), by = "id")
  #   cal.dat <- cal.dat[sample(1:nrow(cal.dat), 10000, replace = FALSE),]
  # 
  #   bw <<- cv_bw(a = cal.dat$a, psi = cal.dat$psi, weights = cal.dat$wts,
  #                bw.seq = seq(0.1, 4, by = 0.1), folds = 10)
  # 
  #   rm(wts, list.cal, psi.cal.new, cal.dat); gc()
  # 
  # }

  # fit exposure response curves
  target <- count_erf(resid.lm = resid.lm, resid.cal = resid.cal, muhat.mat = muhat.mat, log.pop = log.pop, w.id = w.id, 
                      a = a_x, x.id = x.id, bw = 1.8, a.vals = a.vals, phat.vals = phat.vals, se.fit = TRUE)
  
  print(paste0("Fit Complete: Scenario ", i))
  
  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm), n.lm = target$n.lm,
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal), n.cal = target$n.cal,
                         linear.lm = predict(target$fit.lm, newdata = data.frame(a = a.vals)),
                         linear.cal = predict(target$fit.cal, newdata = data.frame(a = a.vals)))
  
  extra <- list(lm.coef = target$fit.lm$coefficients,
                cal.coef = target$fit.cal$coefficients,
                lm.vcov = vcov(target$fit.lm),
                cal.vcov = vcov(target$fit.cal))
  
  save(individual_data, zip_data, est_data, extra,
       file = paste0(dir_out, sname$dual, "_", sname$race, ".RData"))
  
}
