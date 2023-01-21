dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
load(paste0(dir_data, "qd/", scenario$dual, "_", scenario$race, ".RData"))
     

for (j in 1:scenarios) {
  
  
  x <- new.data$x
  a_x <- x$pm25
  x.tmp <- zip.tmp[,-c(1,3)]
  x.tmp$year <- factor(x.tmp$year)
  
  x <- x.tmp %>% mutate_if(is.numeric, scale)
  x.mat <- model.matrix(~ ., data = data.frame(x))
  astar <- c(a_x - mean(a_x))/var(a_x)
  astar2 <- c((a - mean(a))^2/var(a_x) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                   target = c(length(a), rep(0, ncol(x.mat) + 1)))
  cal <- mod$weights
  
  
}



