# estimate gam outcome model and combine with BOTH lm and calibrated weights
gam_models <- function(y, a, w, ipw, cal, cal_trunc, id,
                       a.vals, log.pop = NULL, trunc = 0.01, ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  # need to model residuals as rates
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm + splines
  mumod <- gam(ybar ~ s(a, df = 6) + . + (a:.) - a, weights = exp(log.pop), model = FALSE,
               data = data.frame(ybar = ybar, a = a, w), family = quasipoisson())
  muhat <- mumod$fitted.values
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # pseudo outcome
  resid.lm <- c(ybar - muhat)*ipw
  resid.cal <- c(ybar - muhat)*cal
  resid.cal_trunc <- c(ybar - muhat)*cal_trunc
  
  out <- list(resid.lm = resid.lm, resid.cal = resid.cal, 
              resid.cal_trunc = resid.cal_trunc, muhat.mat = muhat.mat,
              id = id, log.pop = log.pop)
  
  return(out)
  
}

# estimate glm and combine with one weight
gam_models_boot <- function(y, a, w, weights, id, a.vals, log.pop = NULL, trunc = 0.01, ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm + splines
  mumod <- gam(ybar ~ s(a, df = 6) + . + (a:.) - a, weights = exp(log.pop), model = FALSE,
               data = data.frame(ybar = ybar, a = a, w), family = quasipoisson())
  muhat <- mumod$fitted.values
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # pseudo outcome
  resid <- c(ybar - muhat)*weights
  
  out <- list(resid = resid, muhat.mat = muhat.mat, 
              id = id, log.pop = log.pop)
  
  return(out)
  
}
