# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
count_erf <- function(resid.lm, resid.cal, resid.cal_trunc, log.pop, 
                      muhat.mat, w.id, a, x.id, phat.vals = NULL,
                      a.vals = seq(min(a), max(a), length.out = 100), 
                      bw = c(1,2,3), se.fit = TRUE) {	
  
  # Separate Data into List
  mat.list <- split(cbind(exp(log.pop), resid.lm, resid.cal,
                          resid.cal_trunc, muhat.mat), w.id)
  
  # Aggregate by ZIP-code-year
  mat <- do.call(rbind, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = length(a.vals) + 4)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mat.pool <- data.frame(id = names(mat.list), mat)
  mhat.vals <- colMeans(mat.pool[,-c(1:4)], na.rm = TRUE)
  resid.dat <- inner_join(mat.pool[,1:4], data.frame(a = a, id = x.id), by = "id")
  resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
  
  # Pseudo-Outcomes
  resid.dat$psi.lm <- with(resid.dat, X1 + mhat)
  resid.dat$psi.cal <- with(resid.dat, X2 + mhat)
  resid.dat$psi.cal_trunc <- with(resid.dat, X3 + mhat)
  
  if (length(bw) != 3) {
    warning("length(bandwidth) != 3")
    bw <- rep(bw[1], 3)
  }
  
  if (se.fit) {
    
    if (is.null(phat.vals)) {
      warning("Setting phat.vals = rep(1, length(a.vals))")
      phat.vals <- rep(1, length(a.vals))
    }
    
    # Integration Matrix
    mhat.mat <- matrix(rep(mhat.vals, nrow(mat.pool[,-(1:4)])), byrow = TRUE, 
                       nrow = nrow(mat.pool[,-(1:4)]))
    phat.mat <- matrix(rep(phat.vals, nrow(mat.pool[,-(1:4)])), byrow = TRUE,
                       nrow = nrow(mat.pool[,-(1:4)]))
    int.mat <- (mat.pool[,-(1:4)] - mhat.mat)*phat.mat
    
    # KWLS Regression
    out.lm <- sapply(a.vals, kern_est, psi = resid.dat$psi.lm, a = resid.dat$a, 
                     bw = bw[1], a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    out.cal <- sapply(a.vals, kern_est, psi = resid.dat$psi.cal, a = resid.dat$a, 
                      bw = bw[2], a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    out.cal_trunc <- sapply(a.vals, kern_est, psi = resid.dat$psi.cal_trunc, a = resid.dat$a, 
                            bw = bw[3], a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    
    estimate.lm <- out.lm[1,]
    variance.lm <- out.lm[2,]
    estimate.cal <- out.cal[1,]
    variance.cal <- out.cal[2,]
    estimate.cal_trunc <- out.cal_trunc[1,]
    variance.cal_trunc <- out.cal_trunc[2,]
    
    return(list(estimate.lm = estimate.lm, variance.lm = variance.lm,
                estimate.cal = estimate.cal, variance.cal = variance.cal,
                estimate.cal_trunc = estimate.cal_trunc, variance.cal_trunc = variance.cal_trunc))
    
  } else {
    
    # KWLS Regression
    out.lm <- approx(locpoly(resid.dat$a, resid.dat$psi.lm, bandwidth = bw[1]), xout = a.vals)$y
    out.cal <- approx(locpoly(resid.dat$a, resid.dat$psi.cal, bandwidth = bw[2]), xout = a.vals)$y
    out.cal_trunc <- approx(locpoly(resid.dat$a, resid.dat$psi.cal_trunc, bandwidth = bw[3]), xout = a.vals)$y
    
    return(list(estimate.lm = out.lm, estimate.cal = out.cal,
                estimate.cal_trunc = out.cal_trunc))
    
  }
  
}

# count_erf where we consider only one weight implementation
count_erf_boot <- function(resid, log.pop, muhat.mat, w.id, a, x.id, phat.vals = NULL,
                           a.vals = seq(min(a), max(a), length.out = 100), 
                           bw = 1, se.fit = TRUE) {	
  
  # Separate Data into List
  mat.list <- split(cbind(exp(log.pop), resid, muhat.mat), w.id)
  
  # Aggregate by ZIP-code-year
  mat <- do.call(rbind, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = length(a.vals) + 2)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mat.pool <- data.frame(id = names(mat.list), mat)
  mhat.vals <- colMeans(mat.pool[,-(1:2)], na.rm = TRUE)
  resid.dat <- inner_join(mat.pool[,1:2], data.frame(a = a, id = x.id), by = "id")
  resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
  
  # Pseudo-Outcomes
  resid.dat$psi <- with(resid.dat, X1 + mhat)
  
  if (se.fit) {
    
    if (is.null(phat.vals)) {
      warning("Setting phat.vals = rep(1, length(a.vals))")
      phat.vals <- rep(1, length(a.vals))
    }
    
    # Integration Matrix
    mhat.mat <- matrix(rep(mhat.vals, nrow(mat.pool[,-(1:2)])), byrow = TRUE, 
                       nrow = nrow(mat.pool[,-(1:2)]))
    phat.mat <- matrix(rep(phat.vals, nrow(mat.pool[,-(1:2)])), byrow = TRUE,
                       nrow = nrow(mat.pool[,-(1:2)]))
    int.mat <- (mat.pool[,-(1:2)] - mhat.mat)*phat.mat
    
    # KWLS Regression
    out <- sapply(a.vals, kern_est, psi = resid.dat$psi, a = resid.dat$a, 
                  bw = bw[1], a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    
    estimate <- out[1,]
    variance <- out[2,]
    
    return(list(estimate = estimate, variance = variance))
    
  } else {
    
    out <- approx(locpoly(resid.dat$a, resid.dat$psi, bandwidth = bw[1]), xout = a.vals)$y
    
    return(list(estimate = out))
    
  }
  
}

## kernel estimation
kern_est <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ a.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  n.std <- sum(weights*k.std) 
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- c(g.std %*% b)
    
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    g.mat <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.mat * kern.mat * int.mat
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

## k-fold cross-validated bandwidth
cv_bw <- function(a, psi, weights = NULL, folds = 5, bw.seq = seq(0.1, 5, by = 0.1)) {
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  fdx <- sample(x = folds, size = n, replace = TRUE)
  
  cv.mat <- sapply(bw.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a[fdx == k], kern_est, psi = psi[fdx != k], a = a[fdx != k], 
                      weights = weights[fdx != k], bw = h, se.fit = FALSE)[1,]
      cv.vec[k] <- sqrt(mean((psi[fdx == k] - preds)^2, na.rm = TRUE))
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  bw <- bw.seq[which.min(cv.err)]
  
  return(bw)
  
}

## Leave-one-out cross-validated bandwidth
w.fn <- function(h, a, a.vals) {
  
  w.avals <- sapply(a.vals, function(a.tmp, ...) {
    a.std <- (a - a.tmp) / h
    k.std <- dnorm(a.std) / h
    return(mean(a.std^2 * k.std) * (dnorm(0) / h) /
             (mean(k.std) * mean(a.std^2 * k.std) - mean(a.std * k.std)^2))
  })
  
  return(w.avals / length(a))
  
}

hatvals <- function(h, a, a.vals) {
  approx(a.vals, w.fn(h = h, a = a, a.vals = a.vals), xout = a)$y
}

cts.eff.fn <- function(psi, a, h) {
  approx(locpoly(a, psi, bandwidth = h), xout = a)$y
}

risk.fn <- function(h, psi, a, a.vals) {
  hats <- hatvals(h = h, a = a, a.vals = a.vals)
  sqrt(mean(((psi - cts.eff.fn(psi = psi, a = a, h = h)) / (1 - hats))^2, na.rm = TRUE))
}
