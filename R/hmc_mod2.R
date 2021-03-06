

betaLocHmcUpdate <- function(y.l, w, x.l, beta.loc, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta.loc)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta.loc)
  pStar <- p0 - 0.5 * delta * dU_beta_loc(y.l, w, x.l, beta.loc)
  
  # first full step for position
  bStar <- bcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_beta_loc(y.l, w, x.l, bStar)
    
    # position
    bStar <- bStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_beta_loc(y.l, w, x.l, bStar)
  
  # evaluate energies
  U0 <- U_beta_loc(y.l, w, x.l, bcurr)
  UStar <- U_beta_loc(y.l, w, x.l, bStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    alpha <- 0
  }
  
  if (runif(1, 0, 1) < alpha){
    bnext <- bStar
    accept <- 1
  } else {
    bnext <- bcurr
    accept <- 0
  }
  
  out <- list()
  out$b <- bnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


dU_beta_loc <- function(y.l, w, x.l, beta.loc){
  
  grad <- array(0, c(length(beta.loc), 1))
  
  # location contribution
  loc_pred <- x.l %*% beta.loc + w
  for (j in 1:length(beta.loc)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x.l[i,j] * (y.l[i] - expit(loc_pred[i]))
    }
  }
  
  # prior contribution
  prior_var_inv <- diag(rep(1/100, length(beta.loc)))
  grad <- grad + t(-t(beta.loc) %*% prior_var_inv)
  return(-grad)
  
}


U_beta_loc <- function(y.l, w, x.l, beta.loc){
  
  # likelihood
  logd <- 0
  lin_preds <- x.l %*% beta.loc + w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  prior_beta_mean <- rep(0, length(beta.loc))
  prior_beta_var <- diag(rep(100, length(beta.loc)))
  
  logd <- logd + dmvnorm(as.numeric(beta.loc), prior_beta_mean, prior_beta_var, log=T)
  
  return(-logd)
  
}


Uw_v2 <- function(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, w, sigma, loc.stats, offset){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- x.l %*% beta.loc + w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- offset + x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # likelihood: control counts
  count_pred <- offset + x.c %*% beta.co + alpha.co * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.co)){
    logd <- logd + dpois(y.co[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_v2 <- function(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, w, sigma.inv, loc.stats, offset){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  loc.pred <- x.l %*% beta.loc + w
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(loc.pred[i])
  }
  
  # case contribution
  lin.count <- offset + x.c %*% beta.ca
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # control contribution
  lin.count <- offset + x.c %*% beta.co
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.co * (y.co[i] - exp(lin.count[i] + alpha.co * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateV2 <- function(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, w, sigma, sigma.inv, loc.stats, delta, L, offset){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_v2(y.l, x.c, x.l,  y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, wcurr, sigma.inv, loc.stats, offset)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_v2(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats, offset)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_v2(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats, offset)
  
  # evaluate energies
  U0 <- Uw_v2(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, wcurr, sigma, loc.stats, offset)
  UStar <- Uw_v2(y.l, x.c, x.l, y.ca, alpha.ca, beta.ca, beta.loc, y.co, alpha.co, beta.co, wStar, sigma, loc.stats, offset)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    alpha <- 0
  }
  
  if (runif(1, 0, 1) < alpha){
    wnext <- wStar
    accept <- 1
  } else {
    wnext <- wcurr
    accept <- 0
  }
  
  out <- list()
  out$w <- wnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}
