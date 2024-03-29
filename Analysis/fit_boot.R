library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(KernSmooth)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), 
                         race = c("all", "white", "black"))
# scenarios <- expand.grid(dual = c("both", "high", "low"), 
#                          race = c("asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)
boot.iter <- 1000 # bootstrap iterations
bw <- 1.8 # KWLS bandwidth

### M-out-of-N Bootstrap

dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR/'

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
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  zip_data <- new_data$x
  individual_data <- new_data$w
  individual_data$id <- paste(individual_data$zip, individual_data$year, sep = "-")
  u.zip <- unique(individual_data$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  
  boot_list <- mclapply(1:boot.iter, function(b, ...) {
    
    index <- sample(1:length(u.zip), m, replace = TRUE)  # initialize bootstrap  index
    
    ## GPS Model
    
    # zip-code data
    x <- bootstrap_data(data = zip_data, index = index, u.zip = u.zip)
    x.tmp <- subset(x, select = -c(zip, pm25, id, boot.id, ipw, cal, cal_trunc))
    x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
    
    # fit calibration weights
    x.mat <- model.matrix(~ ., data = data.frame(x.tmp))
    astar <- c(x$pm25 - mean(x$pm25))/var(x$pm25)
    astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
    mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                     target = c(nrow(x), rep(0, ncol(x.mat) + 1)))
    x$cal <- x$cal_trunc <- mod$weights
    
    # truncation
    trunc0 <- quantile(x$cal, 0.005)
    trunc1 <- quantile(x$cal, 0.995)
    x$cal_trunc[x$cal < trunc0] <- trunc0
    x$cal_trunc[x$cal > trunc1] <- trunc1
    
    ## Outcome Model
    
    # individual-level data
    w <- bootstrap_data(data = individual_data, index = index, u.zip = u.zip)
    wx <- inner_join(w, subset(x, select = -c(id, zip, year), by = "boot.id"))
    
    # remove collinear terms and identifiers
    if (scenario$dual == "both" & scenario$race == "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, 
                                      id, boot.id, ipw, cal, cal_trunc))
    } else if (scenario$dual == "both" & scenario$race != "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, race, 
                                      id, boot.id, ipw, cal, cal_trunc))
    } else if (scenario$dual != "both" & scenario$race == "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, 
                                      id, boot.id, ipw, cal, cal_trunc))
    } else if (scenario$dual != "both" & scenario$race != "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, race, 
                                      id, boot.id, ipw, cal, cal_trunc))
    }
    
    # fit gam model
    z <- gam_models_boot(y = wx$dead, a = wx$pm25, log.pop = log(wx$time_count), 
                         weights = wx$cal_trunc, id = wx$boot.id, w = w.tmp, a.vals = a.vals)
    
    # set bandwidth from original analysis
    target <- count_erf_boot(resid = z$resid, muhat.mat = z$muhat.mat, log.pop = z$log.pop,
                             w.id = z$id, x.id = x$boot.id, a = x$pm25, bw = bw,
                             a.vals = a.vals, phat.vals = phat.vals)
    
    print(paste("Completed Scenario: ", i))
    return(target$estimate)
    
  }, mc.cores = 8)
  
  boot_mat <- do.call(rbind, boot_list)
  colnames(boot_mat) <- a.vals
  
  # save output
  save(boot_mat, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_boot.RData"))
  
}