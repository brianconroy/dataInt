


logisticGp <- function(y, d, n.sample, burnin, L, proposal.sd.theta=0.3,
                      w_initial=NULL, theta_initial=NULL, phi_initial=NULL,
                      prior_phi, prior_theta){
  
  
  # initial values
  N.w <- ncol(d)
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 10, 15)
  } else {
    phi.i <- phi_initial
  }
  
  # storage
  accept <- rep(0, 2)
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))

  # dual averaging quantities
  w_tuning <- initialize_tuning(m=700, target=0.75)
  deltas_w <- c()
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateLogit(y, w.i, sigma.i, sigma.inv.i, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      accept[2] <- accept[2] + theta.out$accept

    }
    
    w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
    deltas_w <- c(deltas_w, w_tuning$delta_curr)
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w

  return(output)
  
  
}


Uw_logit <- function(y, w, sigma){
  
  logd <- 0
  
  # likelihood
  probs <- expit(w)
  for (i in 1:length(y)){
    logd <- logd + dbinom(y[i], size=1, prob=probs[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_logit <- function(y, w, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  
  # likelihood contribution
  for (i in 1:length(y)){
    grad[i] <- grad[i] + (y[i] - expit(w[i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateLogit <- function(y, w, sigma, sigma.inv, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_logit(y, wcurr, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_logit(y, wStar, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_logit(y, wStar, sigma.inv)
  
  # evaluate energies
  U0 <- Uw_logit(y, wcurr, sigma)
  UStar <- Uw_logit(y, wStar, sigma)
  
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


view_logistic_output <- function(output){
  
  par(mfrow=c(1,3))
  plot(y=colMeans(output$samples.w), x=W, ylab='Estimated W', main='A)'); abline(0, 1, col=2)
  padded_plot(output$samples.theta, Theta, title='B)')
  padded_plot(output$samples.phi, Phi, title='C)')
  
}
