

#####################################
# Fit the Besag-York-Mollie intrinsic
# autoregressive model
#####################################


carBYM <- function(data, proposal.sd.beta, proposal.sd.phi, proposal.sd.theta,
                    prior.var.beta=NULL, prior.tau2=NULL,
                    beta.initial=NULL, phi.initial=NULL, tau2.initial=NULL){
  
  
  #### setup
  K <- length(data$y)
  p <- data$p
  X <- data$x
  X.standardised <- data$x.standardised
  Y <- data$y
  
  
  #### priors
  prior.mean.beta <- rep(0, p)
  if (is.null(prior.var.beta)){ prior.var.beta <- rep(1000, p) }
  if (is.null(prior.tau2)){ prior.tau2 <- c(1, 0.01) }
  prior.sigma2 <- c(1, 0.01)
  
  
  #### Initial values
  mod.glm <- glm(Y~X.standardised-1, offset=NULL, family="quasipoisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  if (is.null(beta.initial)){
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  } else {
    beta <- beta.initial
  }
  
  log.Y <- log(Y)
  log.Y[Y==0] <- -0.1  
  res.temp <- log.Y - X.standardised %*% beta.mean
  res.sd <- sd(res.temp, na.rm=TRUE)
  if (is.null(phi.initial)){
    phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
  } else {
    phi <- phi.initial
  }
  if (is.null(tau2.initial)){
    tau2 <- var(phi)/8
  } else {
    tau2 <- tau2.initial
  }
  theta <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
  sigma2 <- var(theta)
  
  
  beta.initial <- beta
  phi.initial <- phi
  tau2.initial <- tau2
  
  
  #### MCMC quantities    
  n.sample <- 125000
  burnin <- 10000
  thin <- 1
  
  ## Storage
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.tau2 <- array(NA, c(n.keep, 1))
  samples.sigma2 <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, K))
  samples.theta <- array(NA, c(n.keep, K))
  
  ## Metropolis quantities
  accept <- rep(0, 3)
  tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
  sigma2.posterior.shape <- prior.sigma2[1] + 0.5 * K
  
  ## CAR quantities
  W.quants <- common.Wcheckformat(W, K)
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta
    offset.beta <- phi + theta
    beta.out <- betaPoissonUpdate(X.standardised, Y, beta, offset.beta,
                                  prior.mean.beta, prior.var.beta, proposal.sd.beta)
    beta <- beta.out$beta
    
    
    ## Sample from phi
    offset.phi <- X.standardised %*% beta + theta
    phi.out <- phiPoissonUpdate(Y, W, nNeighbors, phi, tau2, offset.phi, proposal.sd.phi, rho=1)
    phi <- phi.out$phi
    phi <- phi - mean(phi)
    
    
    ## Sample from theta
    offset.theta <- X.standardised %*% beta + phi
    theta.out <- thetaPoissonUpdate(Y, theta, offset.theta, proposal.sd.theta, sigma2)
    theta <- theta.out$theta
    theta <- theta - mean(theta)
    
    
    ## Sample from tau2
    temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, 1)
    tau2.posterior.scale <- temp2 + prior.tau2[2]
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
    
    
    ## Sample from sigma2
    sigma2.posterior.scale <- prior.sigma2[2] + 0.5*sum(theta^2)
    sigma2 <- 1/rgamma(1, sigma2.posterior.shape, scale=(1/sigma2.posterior.scale)) 
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta[ele,] <- beta
      accept[1] <- accept[1] + beta.out$accept
      samples.phi[ele,] <- phi
      accept[2] <- accept[2] + phi.out$accept
      samples.tau2[ele,] <- tau2
      samples.sigma2[ele,] <- sigma2
      samples.theta[ele,] <- theta
      accept[3] <- accept[3] + theta.out$accept
    }
    
    
  }
  
  
  accept <- accept/n.keep
  accept[2] <- accept[2]/K
  accept[3] <- accept[3]/K
  
  out <- list()
  out$samples.beta <- samples.beta
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$samples.sigma2 <- samples.sigma2
  out$accept <- accept
  out$beta.initial <- beta.initial
  out$phi.initial <- phi.initial
  out$tau2.initial <- tau2
  return(out)
  
  
}
