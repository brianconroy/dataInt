

K <- function(p){
  return(t(p) %*% p/2)
}


betaHmcUpdateTemporal <- function(count.data_time, locs_time, w, u, beta, alpha, delta_c, L_c, offset){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_c * dU_beta_time(count.data_time, locs_time, w, u, bcurr, alpha, offset)
  
  # first full step for position
  bStar <- bcurr + delta_c*pStar
  
  # full steps
  for (jL in 1:c(L_c-1)){
    # momentum
    pStar <- pStar - delta_c * dU_beta_time(count.data_time, locs_time, w, u, bStar, alpha, offset)
    
    # position
    bStar <- bStar + delta_c*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_c * dU_beta_time(count.data_time, locs_time, w, u, bStar, alpha, offset)
  
  # evaluate energies
  U0 <- Ubeta_time(count.data_time, locs_time, w, u, bcurr, alpha, offset)
  UStar <- Ubeta_time(count.data_time, locs_time, w, u, bStar, alpha, offset)
  
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
  out$beta <- bnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


dU_beta_time <- function(count.data_time, locs_time, w, u, beta, alpha, offset){
  
  grad <- array(0, c(length(beta), 1))
  
  for (t in 1:length(count.data_time)){
    
    dat.t <- count.data_time[[t]]
    w.sub.t <- w[locs_time[[t]]$ids]
    u.t <- u[t]
    x.t <- dat.t$x
    y.t <- dat.t$y
    lin_preds <- offset + x.t %*% beta + alpha * (w.sub.t + u.t)
    
    for (j in 1:length(beta)){
      for (i in 1:length(w.sub.t)){
        grad[j] <- grad[j] + x.t[i, j] * (y.t[i] - exp(lin_preds[i,]))
      }
    }
    
  }
  
  if (length(beta) > 1){
    prior_var_inv <- diag(rep(1/100, length(beta)))
  } else{
    prior_var_inv <- 1/100
  }
  
  grad <- grad + t(-t(beta) %*% prior_var_inv)

  return(-grad)
  
}


dU_alpha_time <- function(count.data_time, locs_time, w, u, beta, alpha, prior_mean, prior_var, offset){
  
  grad <- 0
  
  for (t in 1:length(count.data_time)){
    
    dat.t <- count.data_time[[t]]
    w.sub.t <- w[locs_time[[t]]$ids]
    u.t <- u[t]
    x.t <- dat.t$x
    y.t <- dat.t$y
    lin_preds <- offset + x.t %*% beta + alpha * (w.sub.t + u.t)
    
    for (i in 1:length(w.sub.t)){
      grad <- grad + (w.sub.t[i] + u.t) * (y.t[i] - exp(lin_preds[i,]))
    }
    
  }
  
  prior_var_inv <- 1/prior_var
  grad <- grad - (alpha - prior_mean) * prior_var_inv
  
  return(-grad)
  
}


Ubeta_time <- function(count.data_time, locs_time, w, u, beta, alpha, offset){
  
  # likelihood
  logd <- 0
  for (t in 1:length(count.data_time)){
    
    dat.t <- count.data_time[[t]]
    w.sub.t <- w[locs_time[[t]]$ids]
    u.t <- u[t]
    x.t <- dat.t$x
    y.t <- dat.t$y
    lin_preds <- offset + x.t %*% beta + alpha * (w.sub.t + u.t)
    
    for (i in 1:length(y.t)){
      logd <- logd + dpois(y.t[i], lambda=exp(lin_preds[i]), log=T)
    }
    
  }
  
  # prior
  prior_beta_mean <- rep(0, length(beta))
  if (length(beta) > 1){
    prior_beta_var <- diag(rep(100, length(beta)))
  } else {
    prior_beta_var <- matrix(100)
  }
  
  logd <- logd + dmvnorm(as.numeric(beta), prior_beta_mean, prior_beta_var, log=T)
  
  return(-logd)
  
}


Ualpha_time <- function(count.data_time, locs_time, w, u, beta, alpha, prior_mean, prior_var, offset){
  
  # likelihood
  logd <- 0
  
  for (t in 1:length(count.data_time)){
    
    dat.t <- count.data_time[[t]]
    w.sub.t <- w[locs_time[[t]]$ids]
    u.t <- u[t]
    x.t <- dat.t$x
    y.t <- dat.t$y
    lin_preds <- offset + x.t %*% beta + alpha * (w.sub.t + u.t)
    
    for (i in 1:length(y.t)){
      logd <- logd + dpois(y.t[i], lambda=exp(lin_preds[i]), log=T)
    }
  }
  
  # prior
  logd <- logd + dnorm(alpha, prior_mean, prior_var, log=T)
  
  return(-logd)
  
}



alphaHmcUpdateTemporal <- function(count.data_time, locs_time, w, u, beta, alpha, delta_a, prior_mean, prior_var, L_a, offset){
  
  
  # sample random momentum
  p0 <- rnorm(1)
  
  # simulate Hamiltonian dynamics
  acurr <- alpha
  pStar <- p0 - 0.5 * delta_a * dU_alpha_time(count.data_time, locs_time, w, u, beta, acurr, prior_mean, prior_var, offset)
  
  # first full step for position
  aStar <- acurr + delta_a*pStar
  
  # full steps
  for (jL in 1:c(L_a-1)){
    # momentum
    pStar <- pStar - delta_a * dU_alpha_time(count.data_time, locs_time, w, u, beta, aStar, prior_mean, prior_var, offset)
    
    # position
    aStar <- aStar + delta_a*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_a * dU_alpha_time(count.data_time, locs_time, w, u, beta, aStar, prior_mean, prior_var, offset)
  
  # evaluate energies
  U0 <- Ualpha_time(count.data_time, locs_time, w, u, beta, acurr, prior_mean, prior_var, offset)
  UStar <- Ualpha_time(count.data_time, locs_time, w, u, beta, aStar, prior_mean, prior_var, offset)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  a <- min(1, exp((U0 + K0) - (UStar + KStar)))
  
  if (is.na(a)){
    a <- 0
  }
  
  if (runif(1, 0, 1) < a){
    anext <- aStar
    accept <- 1
  } else {
    anext <- acurr
    accept <- 0
  }
  
  out <- list()
  out$alpha <- anext
  out$accept <- accept
  out$a <- a
  return(out)
  
}


initialize_tuning <- function(m, target){
  return(list(
    M_adapt=m,
    delta_curr=0.05,
    delta_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, accepted){
  
  if (accepted){
    tuning$reject_streak <- 0
  } else{
    tuning$reject_streak <- tuning$reject_streak + 1
  }
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    delta_bar <- tuning$delta_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a)
    log_delta_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_delta_bar <- i^{-kappa} * log_delta_curr + (1 - i^{-kappa}) * log(delta_bar)
    delta_curr <- exp(log_delta_curr)
    delta_bar <- exp(log_delta_bar)
    
    tuning$delta_bar <- delta_bar
    tuning$HbarM <- HbarM
    tuning$delta_curr <- delta_curr
    
  } else{
    
    tuning$delta_curr <- tuning$delta_bar
    if (tuning$reject_streak > 1000) {
      tuning$delta_curr <- 0.90 * tuning$delta_curr
      tuning$delta_bar <- 0.90 * tuning$delta_bar
    }
    
  }
  
  return(tuning)
  
}


Uw_time <- function(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, w, u, sigma, offset){
  
  logd <- 0
  
  for (t in 1:length(case.data_time)){
    
    # likelihood: locations
    u.t <- u[t]
    loc_pred <- w + u.t
    y.l_t <- locs_time[[t]]$status
    for (i in 1:length(y.l_t)){
      logd <- logd + dbinom(y.l_t[i], size=1, prob=expit(loc_pred[i]), log=T)
    }
    
    # likelihood: case counts
    y.ca_t <- case.data_time[[t]]$y
    x.t <- case.data_time[[t]]$x
    w.sub.t <- w[as.logical(locs_time[[t]]$status)]
    count_pred <- offset + x.t %*% beta.ca + alpha.ca * (w.sub.t + u.t)
    rates <- exp(count_pred)
    for (i in 1:length(y.ca_t)){
      logd <- logd + dpois(y.ca_t[i], lambda=rates[i], log=T)
    }
    
    # likelihood: control counts
    y.co_t <- ctrl.data_time[[t]]$y
    count_pred <- offset + x.t %*% beta.co + alpha.co * (w.sub.t + u.t)
    rates <- exp(count_pred)
    for (i in 1:length(y.co_t)){
      logd <- logd + dpois(y.co_t[i], lambda=rates[i], log=T)
    }
    
  }

  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_time <- function(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, w, u, sigma.inv, offset){
  
  grad <- array(0, c(length(w), 1))
  
  for (t in 1:length(case.data_time)){
    
    y.l_t <- locs_time[[t]]$status
    u.t <- u[t]
    x.t <- case.data_time[[t]]$x
    
    # location contribution
    for (i in 1:length(w)){
      grad[i] <- grad[i] + y.l_t[i] - expit(w[i] + u.t)
    }
    
    # case contribution
    y.ca_t <- case.data_time[[t]]$y
    lin.count <- offset + x.t %*% beta.ca
    for (i in 1:length(locs_time[[t]]$ids)){
      id.i <- locs_time[[t]]$ids[i]
      grad[id.i] <- grad[id.i] + alpha.ca * (y.ca_t[i] - exp(lin.count[i] + alpha.ca * (w[id.i] + u.t)))
    }
    
    # control contribution
    y.co_t <- ctrl.data_time[[t]]$y
    lin.count <- offset + x.t %*% beta.co
    for (i in 1:length(locs_time[[t]]$ids)){
      id.i <- locs_time[[t]]$ids[i]
      grad[id.i] <- grad[id.i] + alpha.co * (y.co_t[i] - exp(lin.count[i] + alpha.co * (w[id.i] + u.t)))
    }
    
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateTemporal <- function(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, w, u, sigma, sigma.inv, delta, L, offset){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_time(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, wcurr, u, sigma.inv, offset)

  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_time(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, wStar, u, sigma.inv, offset)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_time(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, wStar, u, sigma.inv, offset)
  
  # evaluate energies
  U0 <- Uw_time(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, wcurr, u, sigma, offset)
  UStar <- Uw_time(case.data_time, ctrl.data_time, locs_time, alpha.ca, beta.ca, alpha.co, beta.co, wStar, u, sigma, offset)
  
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
