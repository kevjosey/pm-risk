# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
count_erf <- function(resid.lm, resid.cal, log.pop, muhat.mat, w.id, a, x.id, bw = 1,
                       a.vals = seq(min(a), max(a), length.out = 100), phat.vals = NULL, se.fit = TRUE) {	
  
  # marginalize psi within zip-year
  wts <- do.call(c, lapply(split(exp(log.pop), w.id), sum))
  
  # LM GPS
  list.lm <- split(data.frame(resid = resid.lm, wts = exp(log.pop)), w.id)
  resid.lm.new <- data.frame(resid = do.call(c, lapply(list.lm, function(df) sum(df$resid*df$wts)/sum(df$wts))), 
                           wts = wts, id = names(list.lm))
  lm.dat <- inner_join(resid.lm.new, data.frame(a = a, id = x.id), by = "id")
  
  # Calibration Weights
  list.cal <- split(data.frame(resid = resid.cal, wts = exp(log.pop)) , w.id)
  resid.cal.new <- data.frame(resid = do.call(c, lapply(list.cal, function(df) sum(df$resid*df$wts)/sum(df$wts))), 
                              wts = wts, id = names(list.cal))
  cal.dat <- inner_join(resid.cal.new, data.frame(a = a, id = x.id), by = "id")
  
  rm(list.lm, resid.lm.new, list.cal, resid.cal.new); gc()
  
  # Prediction Matrix
  mat.list <- split(cbind(exp(log.pop), muhat.mat), w.id)

  muhat.mat.new <- do.call(rbind, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = length(a.vals) + 1)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mhat.vals <- apply(muhat.mat.new, 2, weighted.mean, w = wts)
  mhat <- predict(smooth.spline(a.vals, mhat.vals), x = cal.dat$pm25)$y
  
  if (is.null(phat.vals))
    phat.vals <- rep(1, length(a.vals))
  
  # Integration Matrix
  mhat.mat <- matrix(rep(mhat.vals, nrow(muhat.mat.new)), byrow = TRUE, nrow = nrow(muhat.mat.new))
  phat.mat <- matrix(rep(phat.vals, nrow(muhat.mat.new)), byrow = TRUE, nrow = nrow(muhat.mat.new))
  int.mat <- (muhat.mat.new - mhat.mat)*phat.mat
  
  rm(phat.mat, phat.vals, mhat.mat, mhat.vals); gc()
  
  psi.lm <- lm.dat$resid + mhat
  psi.cal <- cal.dat$resid + mhat
  
  # KWLS regression
  out.lm <- sapply(a.vals, kern_est, psi = psi.lm, a = lm.dat$a, bw = bw,
                   a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)

  out.cal <- sapply(a.vals, kern_est, psi = psi.cal, a = cal.dat$a, bw = bw,
                    a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
  
  # linear model approx
  fit.lm <- lm(psi ~ a, weights = lm.dat$wts, data = lm.dat)
  fit.cal <- lm(psi ~ a, weights = cal.dat$wts, data = cal.dat)
  
  if (se.fit) {
    
    estimate.lm <- out.lm[1,]
    variance.lm <- out.lm[2,]
    n.lm <- out.lm[3,]
    estimate.cal <- out.cal[1,]
    variance.cal <- out.cal[2,]
    n.cal <- out.cal[3,]
    
    return(list(estimate.lm = estimate.lm, variance.lm = variance.lm, n.lm = n.lm, fit.lm = fit.lm,
                estimate.cal = estimate.cal, variance.cal = variance.cal, n.cal = n.cal, fit.cal = fit.cal))
    
  } else {
    
    estimate.lm <- out.lm[1,]
    n.lm <- out.lm[2,]
    estimate.cal <- out.cal[1,]
    n.cal <- out.cal[2,]
    
    return(list( estimate.lm = out.lm, n.lm = n.lm, fit.lm = fit.lm, 
                 estimate.cal = out.cal, n.cal = n.cal, fit.cal = fit.cal))
    
  }
  
}

# cross validated bandwidth
cv_bw <- function(a, psi, weights = NULL, folds = 5, bw.seq = seq(0.05, 1, by = 0.05)) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  fdx <- sample(x = folds, size = n, replace = TRUE)
  
  cv.mat <- sapply(bw.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a[fdx == k], kern_est, psi = psi[fdx != k], a = a[fdx != k], 
                      weights = weights[fdx != k], bw = h, se.fit = FALSE)[1,]
      cv.vec[k] <- weighted.mean(abs(psi[fdx == k] - preds), w = weights[fdx == k], na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  bw <- bw.seq[which.min(cv.err)]
  
  return(bw)
  
}
