# estimate gam outcome model and combine with BOTH lm and calibrated weights
gam_models <- function(y, a, w, ipw, cal, a.vals, log.pop = NULL, trunc = 0.01, ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  w <- data.frame(w)
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm
  # need better coding to generalize
  mumod <- glm(ybar ~ ns(a, 6) + . - a, weights = exp(log.pop),
               data = data.frame(ybar = ybar, a = a, w), family = quasipoisson())
  muhat <- mumod$fitted.values
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # pseudo outcome
  resid.lm <- c(ybar - muhat)*ipw
  resid.cal <- c(ybar - muhat)*cal
  
  out <- list(resid.lm = resid.lm,  resid.cal = resid.cal,
              muhat.mat = muhat.mat, log.pop = log.pop)
  
  return(out)
  
}

# estimate gam and combine ONLY with LM
gam_models_lm <- function(y, a, w, ipw, a.vals, log.pop = NULL, trunc = 0.01, ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  # get ybar = y/n
  w <- data.frame(w)
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm
  # need better coding to generalize
  mumod <- glm(ybar ~ ns(a, 6) + . - a, weights = exp(log.pop),
               data = data.frame(ybar = ybar, a = a, w), family = quasipoisson())
  muhat <- mumod$fitted.values
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # pseudo outcome
  resid.lm <- c(ybar - muhat)*ipw
  
  out <- list(resid.lm = resid.lm, muhat.mat = muhat.mat, log.pop = log.pop)
  
  return(out)
  
}
