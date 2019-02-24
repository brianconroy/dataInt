

logisticMVGP <- function(y, d, n.sample, burnin, L, proposal.sd.theta=0.3,
                         w_initial=NULL, theta_initial=NULL, t_initial=NULL,
                         prior_t, prior_theta){
  
  
  N.d <- length(y)
  
  # initial values
  N.w <- ncol(d)
  if (is.null(w_initial)){
    w.i <- rnorm(N.w*N.d)
  } else {
    w.i <- w_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(t_initial)){
    t.i <- matrix(c(3, 0, 0, 3), nrow=2)
  } else {
    t.i <- t_initial
  }
  H.i <- Exponential(d, range=theta.i, phi=1)
  H.inv.i <- solve(H.i)
  
  Omega <- prior_t$scale
  r <- prior_t$df
  
  # storage
  accept <- rep(0, 2)
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w * N.d))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.t <- array(NA, c(n.keep, N.d * N.d))
  
  # dual averaging quantities
  w_tuning <- initialize_tuning(m=700, target=0.75)
  deltas_w <- c()
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- kronecker(H.i, t.i)
    sigma.inv.i <- kronecker(solve(H.i), solve(t.i))
    w.out.i <- wHmcUpdateMVGPLogit(w.i, sigma.i, sigma.inv.i, y, w_tuning$delta_curr, L)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMVGPupdate(H.i, t.i, w.i, d, theta.i, proposal.sd.theta, prior_theta)
    theta.i <- theta.out$theta
    H.i <- Exponential(d, range=theta.i, phi=1)
    H.inv.i <- solve(H.i)
    
    ## sample from T
    r_ <- r + N.w
    Omega_ <- Omega
    w1.i <- w.i[seq(1, length(w.i), by=N.d)]
    w2.i <- w.i[seq(2, length(w.i), by=N.d)]
    for (a in 1:N.w){
      for (b in 1:N.w){
        Omega_ <- Omega_ + H.inv.i[a, b] * matrix(c(w1.i[b], w2.i[b])) %*% t(matrix(c(w1.i[a], w2.i[a])))
      }
    }
    t.i <- riwish(r_, Omega_)
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.theta[j,] <- theta.i
      samples.t[j,] <- t.i
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
  output$samples.t <- samples.t
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$burnin <- burnin
  output$n.sample <- n.sample
  output$prior_t <- prior_t
  output$prior_theta <- prior_theta
  output$L <- L
  output$proposal.sd.theta <- proposal.sd.theta
  
  return(output)
  
  
}

#' logisticGp
#' 
#' Fits a spatial logistic regression model.
#' Updates random effects by Hamiltonian Monte Carlo,
#' spatial range by MHRW, and marginal variance by Gibbs. 
#'
#' @param y (vector) vector of outcomes
#' @param d (matrix) distance matrix
#' @param n.sample (numeric) number of samples
#' @param burnin (numeric) mcmc burnin
#' @param L (numeric) hamiltonian monte carlo length parameter
#' @param proposal.sd.theta (numeric) spatial range proposal standard deviation
#' @param w_initial (numeric) initial value of random effects
#' @param theta_initial (numeric) initial value of spatial range
#' @param phi_initial (numeric) initial value of spatial variance
#' @param prior_phi (numeric) shape and scale parameters of phi prior
#' @param prior_theta (numeric) shape and scale parameters of theta prior
#'
#' @return (list) mcmc samples
#' @export
#'
#' @examples
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
    w.out.i <- wHmcUpdateLogit(y, w.i, sigma.i, sigma.inv.i, w_tuning$delta_curr, L)
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
  output$burnin <- burnin
  output$n.sample <- n.sample
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$L <- L
  output$proposal.sd.theta <- proposal.sd.theta

  return(output)
  
  
}


burnin_logisticGp_mcmc <- function(output, n.burn){
  
  n.curr <- output$n.sample - output$burnin
  i.start <- n.burn + 1
  output$burnin <- output$burnin + n.burn
  output$samples.w <- output$samples.w[i.start:n.curr,]
  output$samples.phi <- output$samples.phi[i.start:n.curr]
  output$samples.theta <- output$samples.theta[i.start:n.curr]
  return(output)
  
}


continue_logisticGp_mcmc <- function(data, output, n.sample){
  
  # get initial values
  n.sample.old <- nrow(output$samples.w)
  w_initial <- output$samples.w[n.sample.old,]
  theta_initial <- output$samples.theta[n.sample.old]
  phi_initial <- output$samples.phi[n.sample.old]
  
  # get tuning parameters
  delta_w <- tail(output$deltas_w, 1)
  L <- output$L
  proposal.sd.theta <- output$proposal.sd.theta
  prior_phi <- output$prior_phi
  prior_theta <- output$prior_theta
  
  more_output <- logisticGp(data$y, d, n.sample, burnin=0, L, proposal.sd.theta=proposal.sd.theta,
                            w_initial=w_initial, theta_initial=theta_initial, phi_initial=phi_initial,
                            prior_phi=prior_phi, prior_theta=prior_theta)
  
  # merge samples
  new_output <- output
  new_output$samples.phi <- c(new_output$samples.phi, more_output$samples.phi)
  new_output$samples.theta <- c(new_output$samples.theta, more_output$samples.theta)
  new_output$samples.w <- rbind(new_output$samples.w, more_output$samples.w)
  new_output$deltas_w <- c(new_output$deltas_w, more_output$deltas_w)
  new_output$n.sample <- new_output$n.sample + n.sample
  new_output$accept <- (new_output$n.sample * new_output$accept + n.sample * more_output$accept)/(new_output$n.sample + n.sample)
  return(new_output)
  
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


Uw_mvgp_logit <- function(y, w, sigma){
  
  logd <- 0
  
  n_species <- length(y)
  for (s in 1:n_species){
    
    w_d <- w[seq(s, ncol(sigma), by=n_species)]
    # likelihood: locations
    y.d <- y[[s]]
    for (i in 1:length(y.d)){
      logd <- logd + dbinom(y.d[i], size=1, prob=expit(w_d[i]), log=T)
    }
    
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_mvgp_logit <- function(w, sigma.inv, y){
  
  grad <- array(0, c(length(w), 1))
  
  n_species <- length(y)
  for (s in 1:n_species){
    
    w_d <- w[seq(s, ncol(sigma.inv), by=n_species)]
    d_seq <- seq(s, ncol(sigma.inv), by=n_species)
    
    # location contribution
    y.l <- locs[[s]]$status
    grad[d_seq] <- y.l - expit(w_d)
    
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


wHmcUpdateMVGPLogit <- function(w, sigma, sigma.inv, y, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_mvgp_logit(wcurr, sigma.inv, y)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_mvgp_logit(wStar, sigma.inv, y)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_mvgp_logit(wStar, sigma.inv, y)
  
  # evaluate energies
  U0 <- Uw_mvgp_logit(y, wcurr, sigma)
  UStar <- Uw_mvgp_logit(y, wStar, sigma)
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
