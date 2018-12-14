

Uw <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- x.c %*% beta.c + alpha * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.c)){
    logd <- logd + dpois(y.c[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


Ucase <- function(y, w, x, beta, alpha){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
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


Ualpha <- function(y, w, x, beta, alpha, prior_mean, prior_var){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dnorm(alpha, prior_mean, prior_var, log=T)
  
  return(-logd)
  
}


Ualpha_gamma <- function(y, w, x, beta, alpha, shape, scale, type){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  if (type == 'case'){
    logprior <- dgamma(alpha, shape=shape, scale=scale, log=T)
  } else {
    logprior <- dgamma(-alpha, shape=shape, scale=scale, log=T)
  }
  
  logd <- logd + logprior
  
  return(-logd)
  
}


Ualpha_flat <- function(y, w, x, beta, alpha, prior_lower_bound, prior_upper_bound){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  log_prior <- 0
  if (!is.na(alpha)){
    if (alpha < prior_lower_bound || alpha > prior_upper_bound){
      log_prior <- -Inf
    }
  } else {
    log_prior <- -Inf
  }
  logd <- logd + log_prior
  
  return(-logd)
  
}


dU_w <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma.inv, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(w[i])
  }
  
  # count contribution
  lin.count <- x.c %*% beta.c
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha * (y.c[i] - exp(lin.count[i] + alpha * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


dU_case <- function(y, w, x, beta, alpha){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + alpha * w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - exp(lin_preds[i,]))
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


dU_alpha <- function(y, w, x, beta, alpha, prior_mean, prior_var){
  
  grad <- 0
  lin_preds <- x %*% beta + alpha * w
  
  for (i in 1:length(w)){
    grad <- grad + w[i] * (y[i] - exp(lin_preds[i,]))
  }
  
  prior_var_inv <- 1/prior_var
  
  grad <- grad - (alpha - prior_mean) * prior_var_inv
  return(-grad)
  
}


dU_alpha_gamma <- function(y, w, x, beta, alpha, shape, scale, type){
  
  grad <- 0
  lin_preds <- x %*% beta + alpha * w
  
  for (i in 1:length(w)){
    grad <- grad + w[i] * (y[i] - exp(lin_preds[i,]))
  }
  
  if (type=='case'){
    prior_grad <- (shape - 1)/alpha - 1/scale
  } else{
    prior_grad <- (shape - 1)/(-alpha) - 1/scale
  }
  
  grad <- grad + prior_grad
  return(-grad)
  
}


dU_alpha_flat <- function(y, w, x, beta, alpha){
  
  grad <- 0
  lin_preds <- x %*% beta + alpha * w
  
  for (i in 1:length(w)){
    grad <- grad + w[i] * (y[i] - exp(lin_preds[i,]))
  }
  
  return(-grad)
  
}


K <- function(p){
  return(t(p) %*% p/2)
}


wHmcUpdate <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma, sigma.inv, loc.stats, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w(y.l, x.c, y.c, alpha, beta.c, wcurr, sigma.inv, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w(y.l, x.c, y.c, alpha, beta.c, wStar, sigma.inv, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w(y.l, x.c, y.c, alpha, beta.c, wStar, sigma.inv, loc.stats)
  
  # evaluate energies
  U0 <- Uw(y.l, x.c, y.c, alpha, beta.c, wcurr, sigma, loc.stats)
  UStar <- Uw(y.l, x.c, y.c, alpha, beta.c, wStar, sigma, loc.stats)
  
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



caseHmcUpdate <- function(y, w, x, beta, alpha, delta_c, L_c){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_c * dU_case(y, w, x, bcurr, alpha)
  
  # first full step for position
  bStar <- bcurr + delta_c*pStar
  
  # full steps
  for (jL in 1:c(L_c-1)){
    # momentum
    pStar <- pStar - delta_c * dU_case(y, w, x, bStar, alpha)
    
    # position
    bStar <- bStar + delta_c*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_c * dU_case(y, w, x, bStar, alpha)
  
  # evaluate energies
  U0 <- Ucase(y, w, x, bcurr, alpha)
  UStar <- Ucase(y, w, x, bStar, alpha)
  
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


alphaHmcUpdate <- function(y, w, x, beta, alpha, delta_a, prior_mean, prior_var, L_a){
  
  
  # sample random momentum
  p0 <- rnorm(1)
  
  # simulate Hamiltonian dynamics
  acurr <- alpha
  pStar <- p0 - 0.5 * delta_a * dU_alpha(y, w, x, beta, acurr, prior_mean, prior_var)
  
  # first full step for position
  aStar <- acurr + delta_a*pStar
  
  # full steps
  for (jL in 1:c(L_a-1)){
    # momentum
    pStar <- pStar - delta_a * dU_alpha(y, w, x, beta, aStar, prior_mean, prior_var)
    
    # position
    aStar <- aStar + delta_a*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_a * dU_alpha(y, w, x, beta, aStar, prior_mean, prior_var)
  
  # evaluate energies
  U0 <- Ualpha(y, w, x, beta, acurr, prior_mean, prior_var)
  UStar <- Ualpha(y, w, x, beta, aStar, prior_mean, prior_var)
  
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


alphaHmcUpdate_flat <- function(y, w, x, beta, alpha, delta_a, prior_lower_bound, prior_upper_bound, L_a){
  
  
  # sample random momentum
  p0 <- rnorm(1)
  
  # simulate Hamiltonian dynamics
  acurr <- alpha
  pStar <- p0 - 0.5 * delta_a * dU_alpha_flat(y, w, x, beta, acurr)
  
  # first full step for position
  aStar <- acurr + delta_a*pStar
  
  # full steps
  for (jL in 1:c(L_a-1)){
    # momentum
    pStar <- pStar - delta_a * dU_alpha_flat(y, w, x, beta, aStar)
    
    # position
    aStar <- aStar + delta_a*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_a * dU_alpha_flat(y, w, x, beta, aStar)
  
  # evaluate energies
  U0 <- Ualpha_flat(y, w, x, beta, acurr, prior_lower_bound, prior_upper_bound)
  UStar <- Ualpha_flat(y, w, x, beta, aStar, prior_lower_bound, prior_upper_bound)
  
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


alphaHmcUpdate_gamma <- function(y, w, x, beta, alpha, delta_a, shape, scale, L_a, type){
  
  
  # sample random momentum
  p0 <- rnorm(1)
  
  # simulate Hamiltonian dynamics
  acurr <- alpha
  pStar <- p0 - 0.5 * delta_a * dU_alpha_gamma(y, w, x, beta, acurr, shape, scale, type)
  
  # first full step for position
  aStar <- acurr + delta_a*pStar
  
  # full steps
  for (jL in 1:c(L_a-1)){
    # momentum
    pStar <- pStar - delta_a * dU_alpha_gamma(y, w, x, beta, aStar, shape, scale, type)
    
    # position
    aStar <- aStar + delta_a*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_a * dU_alpha_gamma(y, w, x, beta, aStar, shape, scale, type)
  
  # evaluate energies
  U0 <- Ualpha_gamma(y, w, x, beta, acurr, shape, scale, type)
  UStar <- Ualpha_gamma(y, w, x, beta, aStar, shape, scale, type)
  
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


Uw_cc <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # likelihood: control counts
  count_pred <- x.c %*% beta.co + alpha.co * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.co)){
    logd <- logd + dpois(y.co[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_cc <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma.inv, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(w[i])
  }
  
  # case contribution
  lin.count <- x.c %*% beta.ca
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # control contribution
  lin.count <- x.c %*% beta.co
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.co * (y.co[i] - exp(lin.count[i] + alpha.co * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateCC <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma, sigma.inv, loc.stats, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wcurr, sigma.inv, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats)
  
  # evaluate energies
  U0 <- Uw_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wcurr, sigma, loc.stats)
  UStar <- Uw_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma, loc.stats)
  
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
