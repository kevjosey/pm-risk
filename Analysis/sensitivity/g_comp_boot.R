library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), 
                         race = c("all", "white", "black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)
boot.iter <- 1000 # bootstrap iterations

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

### G-Computation Bootstrap

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GComp_All/'

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
    
    # individual-level data
    w <- bootstrap_data(data = individual_data, index = index, u.zip = u.zip)
    wx <- inner_join(w, subset(x, select = -c(id, zip, year, ipw, cal, cal_trunc), by = "boot.id"))
    
    # remove collinear terms and identifiers
    if (scenario$dual == "both" & scenario$race == "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, id, boot.id))
    } else if (scenario$dual == "both" & scenario$race != "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, race, id, boot.id))
    } else if (scenario$dual != "both" & scenario$race == "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, id, boot.id))
    } else if (scenario$dual != "both" & scenario$race != "all") {
      w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count, dual, race, id, boot.id))
    }
    
    # outcome rates
    wx$ybar <- wx$dead/wx$time_count
    
    # estimate nuisance outcome model with glm + splines
    mumod <- glm(ybar ~ ns(a, 6) + . - a, weights = wx$time_count, model = FALSE,
                 data = data.frame(ybar = wx$ybar, a = wx$pm25, w.tmp), family = quasipoisson())
    
    # predictions along a.vals
    muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
      
      wa.tmp <- data.frame(a = a.tmp, w.tmp)
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
