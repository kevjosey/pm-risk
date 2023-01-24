# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
count_erf <- function(resid.lm, resid.cal, log.pop, muhat.mat, w.id, a, x.id,
                      a.vals = seq(min(a), max(a), length.out = 100), phat.vals = NULL, 
                      bw = 1, se.fit = TRUE) {	
  
  # Separate Data into List
  wts <- do.call(c, lapply(split(exp(log.pop), w.id), sum))
  mat.list <- split(cbind(exp(log.pop), resid.lm, resid.cal, muhat.mat), w.id)
  
  # Aggregate by ZIP-code-year
  mat <- do.call(rbind, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = length(a.vals) + 3)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mat.new <- data.frame(wts = wts, id = names(mat.list), mat)
  muhat.mat.new <- mat.new[,-c(1:4)]
  mhat.vals <- colMeans(muhat.mat.new, na.rm = TRUE)
  resid.dat <- inner_join(mat.new[,1:4], data.frame(a = a, id = x.id), by = "id")
  resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
  
  # Pseudo-Outcomes
  resid.dat$psi.lm <- with(resid.dat, X1 + mhat)
  resid.dat$psi.cal <- with(resid.dat, X2 + mhat)
  
  if (is.null(phat.vals))
    phat.vals <- rep(1, length(a.vals))
  
  # Integration Matrix
  mhat.mat <- matrix(rep(mhat.vals, nrow(muhat.mat.new)), byrow = TRUE, nrow = nrow(muhat.mat.new))
  phat.mat <- matrix(rep(phat.vals, nrow(muhat.mat.new)), byrow = TRUE, nrow = nrow(muhat.mat.new))
  int.mat <- (muhat.mat.new - mhat.mat)*phat.mat
  
  # Spline Regression Sensitivity
  spl.lm <- lm(psi.lm ~ ns(a, 6), data = resid.dat)
  spl.cal <- lm(psi.cal ~ ns(a, 6), data = resid.dat)
  
  # Least Squares Approximations
  fit.lm <- lm(psi.lm ~ a, data = resid.dat)
  fit.cal <- lm(psi.cal ~ a, data = resid.dat)
  
  if (se.fit) {
    
    # KWLS Regression
    out.lm <- sapply(a.vals, kern_est, psi = resid.dat$psi.lm, a = resid.dat$a, bw = bw,
                     a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    
    out.cal <- sapply(a.vals, kern_est, psi = resid.dat$psi.cal, a = resid.dat$a, bw = bw,
                      a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    
    estimate.lm <- out.lm[1,]
    variance.lm <- out.lm[2,]
    estimate.cal <- out.cal[1,]
    variance.cal <- out.cal[2,]
    
    return(list(estimate.lm = estimate.lm, variance.lm = variance.lm, 
                fit.lm = fit.lm, spl.lm = spl.lm,
                estimate.cal = estimate.cal, variance.cal = variance.cal,
                fit.cal = fit.cal, spl.cal = spl.cal))
    
  } else {
    
    out.lm <- approx(locpoly(resid.dat$a, resid$psi.lm, bandwidth = bw), xout = a.vals)$y
    out.cal <- approx(locpoly(resid.dat$a, resid$psi.cal, bandwidth = bw), xout = a.vals)$y
    
    return(list(estimate.lm = out.lm, fit.lm = fit.lm, spl.lm = spl.lm,
                estimate.cal = out.cal, fit.cal = fit.cal, spl.cal = spl.cal))
    
  }
  
}

# count_erf where we ONLY consider KWLS + LM GPS
count_erf_lm <- function(resid.lm, log.pop, muhat.mat, w.id, a, x.id,
                         a.vals = seq(min(a), max(a), length.out = 100), phat.vals = NULL, 
                         bw = 1, se.fit = TRUE) {	
  
  # Separate Data into List
  wts <- do.call(c, lapply(split(exp(log.pop), w.id), sum))
  mat.list <- split(cbind(exp(log.pop), resid.lm, muhat.mat), w.id)
  
  # Aggregate by ZIP-code-year
  mat <- do.call(rbind, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = length(a.vals) + 2)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mat.new <- data.frame(wts = wts, id = names(mat.list), mat)
  muhat.mat.new <- mat.new[,-c(1:3)]
  mhat.vals <- colMeans(muhat.mat.new, na.rm = TRUE)
  resid.dat <- inner_join(mat.new[,1:3], data.frame(a = a, id = x.id), by = "id")
  resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
  
  # Pseudo-Outcomes
  resid.dat$psi.lm <- with(resid.dat, X1 + mhat)
  
  if (is.null(phat.vals))
    phat.vals <- rep(1, length(a.vals))
  
  # Integration Matrix
  mhat.mat <- matrix(rep(mhat.vals, nrow(muhat.mat.new)), byrow = TRUE, nrow = nrow(muhat.mat.new))
  phat.mat <- matrix(rep(phat.vals, nrow(muhat.mat.new)), byrow = TRUE, nrow = nrow(muhat.mat.new))
  int.mat <- (muhat.mat.new - mhat.mat)*phat.mat
  
  if (se.fit) {
    
    # KWLS Regression
    out.lm <- sapply(a.vals, kern_est, psi = resid$psi.lm, a = resid.dat$a, bw = bw,
                     a.vals = a.vals, se.fit = se.fit, int.mat = int.mat)
    
    estimate.lm <- out.lm[1,]
    variance.lm <- out.lm[2,]
    
    return(list(estimate.lm = estimate.lm, variance.lm = variance.lm))
    
  } else {
    
    out.lm <- approx(locpoly(resid.dat$a, resid$psi.lm, bandwidth = bw), xout = a.vals)$y
    
    return(list(estimate.lm = out.lm))
    
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

## k-fold cross validated bandwidth
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
  mean(((psi - cts.eff.fn(psi = psi, a = a, h = h)) / (1 - hats))^2)
}
