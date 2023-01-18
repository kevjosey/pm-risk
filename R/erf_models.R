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
  
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = wts)
  mhat <- predict(smooth.spline(a.vals, mhat.vals), x = cal.dat$pm25)$y
  
  if (is.null(phat.vals))
    phat.vals <- rep(1, length(a.vals))
  
  # Integration Matrix
  mhat.mat <- matrix(rep(mhat.vals, nrow(muhat.mat)), byrow = TRUE, nrow = nrow(muhat.mat))
  phat.mat <- matrix(rep(phat.vals, nrow(muhat.mat)), byrow = TRUE, nrow = nrow(muhat.mat))
  int.mat <- (muhat.mat - mhat.mat)*phat.mat
  
  rm(phat.mat, phat.vals, mhat.mat, mhat.vals); gc()
  
  psi.lm <- lm.dat$resid + mhat
  psi.cal <- cal.dat$resid + mhat
  
  # KWLS regression
  out.lm <- sapply(a.vals, kern_est, psi = lm.dat$psi, a = lm.dat$a,
                   bw = bw, a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)

  out.cal <- sapply(a.vals, kern_est, psi = cal.dat$psi, a = cal.dat$a,
                    bw = bw, a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
  
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

# Kernel weighted least squares
kern_est <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  n.std <- sum(weights*k.std) 
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- c(g.std %*% b)
      
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.vals * kern.mat * int.mat
    int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1], n = n.std))
    
  } else
    return(c(mu = mu, n = n.std))
  
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
