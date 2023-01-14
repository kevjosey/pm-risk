
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(SuperLearner)
library(xgboost)
library(ggplot2)
library(cobalt)

source('/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/pm-risk/R/erf.R')
source('/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/pm-risk/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c("0", "1", ""), race = c("white","black","hispanic","asian",""))
scen_names <- expand.grid(dual = c("0", "1", "2"), race = c("white","black","hispanic","asian","all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(4, 16, length.out = 241)
n.boot <- 1000

# Load/Save models
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/Output/DR_mod/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey_pm25-mortality-np_erc_strata/Output/DR_all/'

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
  resid.lm <- resid.sl <- resid.cal <- NULL
  ind_data <- z_data <- NULL
  
  for (j in 1:length(fn)) {
    
    load(paste0(fn[j]))
    
    w.id <- c(w.id, model_data$w.id)
    log.pop <- c(log.pop, model_data$log.pop)
    nval <- c(nval, sum(exp(model_data$log.pop)))
    
    muhat.mat <- rbind(muhat.mat, model_data$muhat.mat)
    phat.tmp <- rbind(phat.tmp, model_data$phat.vals)
    
    resid.lm <- c(resid.lm, model_data$resid.lm)
    resid.cal <- c(resid.cal, model_data$resid.cal)
    
    ind_data <- rbind(ind_data, individual_data)
    z_data <- rbind(z_data, zip_data)
    
  }
  
  # summary data
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = exp(log.pop))
  phat.vals <- apply(phat.tmp, 2, weighted.mean, w = nval)
  mhat <- predict(smooth.spline(a.vals, mhat.vals), x = ind_data$pm25)$y
  zip_data <- z_data[!duplicated(paste(z_data$zip, z_data$year, sep = "-")),]
  individual_data <- ind_data
  rm(z_data,ind_data,phat.tmp)
  
  # integration matrix
  mhat.mat <- matrix(rep(mhat.vals, nrow(muhat.mat)), byrow = TRUE, nrow = nrow(muhat.mat))
  phat.mat <- matrix(rep(phat.vals, nrow(muhat.mat)), byrow = TRUE, nrow = nrow(muhat.mat))
  int.mat <- (muhat.mat - mhat.mat)*phat.mat
  
  psi.lm <- resid.lm + mhat
  psi.cal <- resid.cal + mhat
  
  x.id <- paste(zip_data$zip, zip_data$year, sep = "-")
  a_x <- zip_data$pm25
  
  target <- count_erf(psi.lm = psi.lm, psi.cal = psi.cal, w.id = w.id, x.id = x.id, a = a_x, 
                      log.pop = log.pop, int.mat = int.mat, bw = 1, a.vals = a.vals, se.fit = TRUE)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal),
                         linear.lm = predict(target$fit.lm, newdata = data.frame(a = a.vals)),
                         linear.cal = predict(target$fit.cal, newdata = data.frame(a = a.vals)))
  
  extra <- list(lm.coef = target$fit.lm$coefficients,
                cal.coef = target$fit.cal$coefficients,
                lm.vcov = vcov(target$fit.lm),
                sl.vcov = vcov(target$fit.sl),
                cal.vcov = vcov(target$fit.cal))
  
  print(paste0("Fit Complete: Scenario ", i))
  save(individual_data, zip_data, est_data, extra,
       file = paste0(dir_out, sname$dual, "_", sname$race,
                     sname$sex, "_", sname$age_break, ".RData"))
  
}