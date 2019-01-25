

Uw_multi <- function(locs, case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, w, sigma){
  
  logd <- 0
  
  loc_pred <- w
  for (d in 1:length(locs)){
    
    # likelihood: locations
    y.d <- locs[[d]]$status
    for (i in 1:length(y.d)){
      logd <- logd + dbinom(y.d[i], size=1, prob=expit(loc_pred[i]), log=T)
    }
    
    # likelihood: case counts
    x.d <- case.data[[d]]$x.standardised
    y.d <- case.data[[d]]$y
    w.sub.d <- w[as.logical(locs[[d]]$status)]
    count_pred <- x.d %*% beta.ca[[d]] + alpha.ca[[d]] * w.sub.d
    rates <- exp(count_pred)
    for (i in 1:length(y.d)){
      logd <- logd + dpois(y.d[i], lambda=rates[i], log=T)
    }
    
    # likelihood: control counts
    y.d <- ctrl.data[[d]]$y
    count_pred <- x.d %*% beta.co[[d]] + alpha.co[[d]] * w.sub.d
    rates <- exp(count_pred)
    for (i in 1:length(y.d)){
      logd <- logd + dpois(y.d[i], lambda=rates[i], log=T)
    }
    
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_multi <- function(case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, w, sigma.inv, locs){
  
  grad <- array(0, c(length(w), 1))
  
  for (d in 1:length(locs)){
    
    # location contribution
    y.l <- locs[[d]]$status
    for (i in 1:length(w)){
      grad[i] <- grad[i] + y.l[i] - expit(w[i])
    }
    
    # case contribution
    x.c <- case.data[[d]]$x.standardised
    y.ca <- case.data[[d]]$y
    lin.count <- x.c %*% beta.ca[[d]]
    for (i in 1:length(locs[[d]]$ids)){
      id.i <- locs[[d]]$ids[i]
      grad[id.i] <- grad[id.i] + alpha.ca[[d]] * (y.ca[i] - exp(lin.count[i] + alpha.ca[[d]] * w[id.i]))
    }
    
    # control contribution
    y.co <- ctrl.data[[d]]$y
    lin.count <- x.c %*% beta.co[[d]]
    for (i in 1:length(locs[[d]]$ids)){
      id.i <- locs[[d]]$ids[i]
      grad[id.i] <- grad[id.i] + alpha.co[[d]] * (y.co[i] - exp(lin.count[i] + alpha.co[[d]] * w[id.i]))
    }
    
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateMulti <- function(case.data, ctrl.data, alpha.ca, beta.ca,
                            alpha.co, beta.co, w, sigma, sigma.inv, 
                            locs, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_multi(case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, wcurr, sigma.inv, locs)

  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_multi(case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, wStar, sigma.inv, locs)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_multi(case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, wStar, sigma.inv, locs)
  
  # evaluate energies
  U0 <- Uw_multi(locs, case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, wcurr, sigma)
  UStar <- Uw_multi(locs, case.data, alpha.ca, beta.ca, ctrl.data, alpha.co, beta.co, wStar, sigma)
  
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