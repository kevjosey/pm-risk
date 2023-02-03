# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
count_erf <- function(psi, w.id, log.pop, int.mat, x.id, a_x, loess = FALSE, 
                      a.vals = seq(min(a), max(a), length.out = 100), bw = NULL, 
                      bw.seq = seq(0.1, 2, by = 0.1), folds = 5, se.fit = TRUE) {	
  
  # select bw if null
  if (is.null(bw))
    bw <- cv_bw(a = a_w, psi = psi, weights = exp(log.pop), folds = folds, bw.seq = bw.seq)
  
  # marginalize psi within zip-year
  wts <- do.call(c, lapply(split(exp(log.pop), w.id), sum))
  
  list.cal <- split(data.frame(psi = psi, wts = exp(log.pop)) , w.id)
  psi.new <- data.frame(psi = do.call(c, lapply(list.cal, function(df) sum(df$psi*df$wts)/sum(df$wts))), 
                        wts = wts, id = names(list.cal))
  dat <- inner_join(psi.new, data.frame(a = a_x, id = x.id), by = "id")
  
  mat.list <- split(cbind(exp(log.pop), int.mat), w.id)
  int.mat.new <- do.call(rbind, lapply(split(cbind(exp(log.pop), int.mat), w.id), 
                                       function(vec) { mat <- matrix(vec, ncol = length(a.vals) + 1);
                                       colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1]) } ))
  
  # KWLS regression
  out <- sapply(a.vals, kern_est, psi = cal.dat$psi, a = cal.dat$a, weights = cal.dat$wts,
                loess = loess, bw = bw, se.fit = se.fit, int.mat = int.mat.new, a.vals = a.vals)
  
  # linear model approx
  fit <- lm(psi ~ a, weights = dat$wts, data = dat)
  
  if (se.fit) {
    
    estimate <- out[1,]
    variance <- out[2,]
    
    return(list(estimate = estimate, variance = variance, fit = fit))
    
  } else {
    
    return(list(estimate = out, fit = fit))
    
  }
  
}

# estimate glm outcome model (Keith - adapt code to fit a Poisson outcome as an excercise)
gam_est <- function(a_w, y, w, w.id, a_x, x, x.id, 
                    a.vals, log.pop = NULL, trunc = 0.01,
                    sl.lib = c("SL.mean","SL.glm"), ...) {
  
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
  mumod <- gam(ybar ~ s(a, 5) + . - a + a:(regionWEST + regionNORTHEAST + regionSOUTH), # need better coding to generalize
               weights = exp(log.pop),
               data = data.frame(ybar = ybar, a = a_w, w), 
               family = quasipoisson())
  muhat <- mumod$fitted.values
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = exp(log.pop))
  
  # LM
  pimod.lm <- lm(a ~ ., data = data.frame(a = a_x, x))
  pimod.vals.lm <- c(pimod.lm$fitted.values, predict(pimod.lm, newdata = data.frame(w)))
  pimod.sd.lm <- sigma(pimod.lm)
  
  # nonparametric denisty - LM
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
  
  # full calibration weights
  x <- x %>% mutate_if(is.numeric, scale)
  x.mat <- model.matrix(~ ., data = data.frame(x))
  astar <- c(a_x - mean(a_x))/var(a_x)
  astar2 <- c((a_x - mean(a_x))^2/var(a_x) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                   target = c(n, rep(0, ncol(x.mat) + 1)))
  
  ipw_x <- mod$weights
  ipw_mat <- inner_join(x = data.frame(id = w.id), 
                        y = data.frame(id = x.id, wts = ipw_x), by = "id")
  ipw_w <- ipw_mat$wts
  ipw <- c(ipw_x, ipw_w)
  
  if (any(ipw_mat$id != w.id))
    stop("id is getting scrambled!")
  
  # pseudo outcome
  resid <- c(ybar - muhat)*ipw[-(1:n)]
  
  out <- list(resid = resid.cal, weights_x = ipw_x, weights_w = ipw_w, 
              muhat.mat = muhat.mat, phat.vals = phat.vals, 
              w.id = w.id, log.pop = log.pop)
  
  return(out)
  
}

# Kernel weighted least squares
kern_est <- function(a.new, a, psi, bw, weights, se.fit = FALSE, 
                     int.mat = NULL, a.vals = NULL, loess = FALSE) {
  
  n <- length(a)
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Tricube Weights
  if (loess) {
    
    a.std <- (a - a.new) / bw
    idx <- which(abs(a.std) <= 1)
    k <- sum(idx)
    
    # subset
    a.std <- a.std[idx]
    psi <- psi[idx]
    weights <- weights[idx]
    int.mat <- int.mat[idx,]
    
    # construct kernel weight
    k.std <- c((1 - abs(a.std)^3)^3)
    g.std <- cbind(1, a.std)
    
  } else {   # Gaussian Kernel
    
    a.std <- (a - a.new) / bw
    k.std <- dnorm(a.std) / bw
    g.std <- cbind(1, a.std)
    
  }
  
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
    
    eta <- c(g.std %*% b)
    
    # LOESS
    if (loess) {
      
      kern.mat <- matrix(rep(c((1 - abs((a.vals - a.new)/bw)^3)^3), k), byrow = T, nrow = k)
      kern.mat[matrix(rep(abs(a.vals - a.new)/bw, k), byrow = T, nrow = k) > 1] <- 0
      g.vals <- matrix(rep(c(a.vals - a.new)/bw, k), byrow = T, nrow = k)
      
      intfn1.mat <- kern.mat * int.mat
      intfn2.mat <- g.vals * kern.mat * int.mat
      int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
      int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
      
    } else {  # Gaussian Kernel
      
      kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
      g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
      
      intfn1.mat <- kern.mat * int.mat
      intfn2.mat <- g.vals * kern.mat * int.mat
      int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
      int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
      
    }
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(mu)
  
}

# k-fold cross validation to select bw
cv_bw <- function(a, psi, weights = NULL, folds = 5, bw.seq = seq(0.1, 2, by = 0.1)) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  idx <- sample(x = folds, size = n, replace = TRUE)
  
  cv.mat <- sapply(bw.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a[idx == k], kern_est, psi = psi[idx != k], a = a[idx != k], 
                      weights = weights[idx != k], bw = h, se.fit = FALSE)
      cv.vec[k] <- mean((psi[idx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  bw <- bw.seq[which.min(cv.err)]
  
  return(bw)
  
}

opt_fun <- function(param, psi, g.std, k.std) {
  
  sum(k.std*(psi - exp(c(g.std %*% param)))^2)
  
}
