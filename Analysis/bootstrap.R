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
                         race = c("all", "white", "black"))[-(2:3),]
# scenarios <- expand.grid(dual = c("both", "high", "low"), 
#                          race = c("asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)
boot.iter <- 500 # bootstrap iterations

### M-out-of-N Bootstrap

dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All/'

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
  
  boot_list <- mclapply(1:boot.iter, function(b, ...) {
    
    m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
    index <- sample(1:length(u.zip), m, replace = TRUE)  # initialize bootstrap  index
    
    ## GPS Model
    
    # zip-code data
    x <- bootstrap_data(data = zip_data, index = index, u.zip = u.zip)
    a <- x$pm25
    x.tmp <- subset(x, select = -c(zip, id, boot.id, pm25,
                                   ipw, cal, cal_trunc))
    x.tmp$year <- factor(x.tmp$year)
    x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
    
    # fit calibration weights
    x.mat <- model.matrix(~ ., data = data.frame(x.tmp))
    astar <- c(a - mean(a))/var(a)
    astar2 <- c((a - mean(a))^2/var(a) - 1)
    mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                     target = c(length(a), rep(0, ncol(x.mat) + 1)))
    x$cal <- mod$weights
    
    # truncation
    trunc0 <- quantile(x$cal, 0.005)
    trunc1 <- quantile(x$cal, 0.995)
    x$cal[x$cal < trunc0] <- trunc0
    x$cal[x$cal > trunc1] <- trunc1
    
    ## Outcome Model
    
    # individual-level data
    w.tmp <- bootstrap_data(data = individual_data, index = index, u.zip = u.zip)
    wx <- inner_join(w.tmp, subset(x, select = -c(id, zip, year, ipw, cal_trunc), by = "boot.id"))
    
    # remove collinear terms and identifiers
    if (scenario$dual == "both" & scenario$race == "all") {
      w <- subset(wx, select = -c(zip, pm25, dead, time_count, id, boot.id, cal))
    } else if (scenario$dual == "both" & scenario$race != "all") {
      w <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                  race, id, boot.id, cal))
    } else if (scenario$dual != "both" & scenario$race == "all") {
      w <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                  dual, id, boot.id, cal))
    } else if (scenario$dual != "both" & scenario$race != "all") {
      w <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                  dual, race, id, boot.id, cal))
    }
    
    # fit gam model
    z <- gam_models_boot(y = wx$dead, a = wx$pm25, log.pop = log(wx$time_count), 
                         weights = wx$cal, id = wx$boot.id, w = w, a.vals = a.vals)
    
    # set bandwidth from whole data
    target <- count_erf_boot(resid = z$resid, muhat.mat = z$muhat.mat, log.pop = z$log.pop,
                             w.id = z$id, x.id = x$boot.id, a = a, bw = 2,
                             a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)
    
    print(paste("Completed Scenario: ", i))
    return(target$estimate)
    
  }, mc.cores = 8)
  
  boot_mat <- do.call(rbind, boot_list)
  colnames(boot_mat) <- a.vals
  
  # save output
  save(boot_mat, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_boot.RData"))
  
}