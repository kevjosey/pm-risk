# estimate gam outcome model and lm/calibrated weights
gam_models <- function(y, a_w, w, w.id, a_x, x, x.id, 
                       a.vals, log.pop = NULL, trunc = 0.01, ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  m <- nrow(w)
  x <- data.frame(x)
  w <- data.frame(w)
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - 1e-6
  
  # estimate nuisance outcome model with glm
  # need better coding to generalize
  mumod <- gam(ybar ~ s(a, 5) + . - a + a:(regionWEST + regionNORTHEAST + regionSOUTH), weights = exp(log.pop),
               data = data.frame(ybar = ybar, a = a_w, w), family = quasipoisson())
  muhat <- mumod$fitted.values
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  ## LM
  pimod.lm <- lm(a ~ ., data = data.frame(a = a_x, x))
  pimod.vals.lm <- c(pimod.lm$fitted.values, predict(pimod.lm, newdata = data.frame(w)))
  pimod.sd.lm <- sigma(pimod.lm)
  
  # nonparametric denisty
  a.std.lm <- c(c(a_x, a_w) - pimod.vals.lm) / pimod.sd.lm
  dens.lm <- density(a.std.lm[1:n])
  pihat.lm <- approx(x = dens.lm$x, y = dens.lm$y, xout = a.std.lm)$y / pimod.sd.lm
  
  pihat.mat.lm <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals.lm) / pimod.sd.lm
    approx(x = dens.lm$x, y = dens.lm$y, xout = std)$y / pimod.sd.lm
  })
  
  phat.vals.lm <- colMeans(pihat.mat.lm[1:n,], na.rm = TRUE)
  phat.lm <- predict(smooth.spline(a.vals, phat.vals.lm), x = c(a_x, a_w))$y
  phat.lm[phat.lm < 0] <- .Machine$double.eps
  
  # truncation
  ipw.lm <- phat.lm/pihat.lm
  trunc0 <- quantile(ipw.lm[1:n], trunc)
  trunc1 <- quantile(ipw.lm[1:n], 1 - trunc)
  ipw.lm[ipw.lm < trunc0] <- trunc0
  ipw.lm[ipw.lm > trunc1] <- trunc1
  
  ## Calibration Weights
  x <- x %>% mutate_if(is.numeric, scale)
  x.mat <- model.matrix(~ ., data = data.frame(x))
  astar <- c(a_x - mean(a_x))/var(a_x)
  astar2 <- c((a_x - mean(a_x))^2/var(a_x) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                   target = c(n, rep(0, ncol(x.mat) + 1)))
  
  ipw.cal_x <- mod$weights
  ipw.cal_mat <- inner_join(x = data.frame(id = w.id), 
                            y = data.frame(id = x.id, wts = ipw.cal_x), 
                            by = "id")
  ipw.cal_w <- ipw.cal_mat$wts
  ipw.cal <- c(ipw.cal_x, ipw.cal_w)
  
  if (any(ipw.cal_mat$id != w.id))
    stop("id is getting scrambled!")
  
  # pseudo outcome
  resid.lm <- c(ybar - muhat)*ipw.lm[-(1:n)]
  resid.cal <- c(ybar - muhat)*ipw.cal[-(1:n)]
  
  out <- list(resid.lm = resid.lm, weights.lm_x = ipw.lm[1:n], weights.lm_w = ipw.lm[-(1:n)],
              resid.cal = resid.cal, weights.cal_x = ipw.cal_x, weights.cal_w = ipw.cal_w, 
              muhat.mat = muhat.mat, phat.vals = phat.vals.lm, w.id = w.id, log.pop = log.pop)
  
  return(out)
  
}
