

dinvgamma <- function (x, shape, scale=1, log=FALSE){
  
  if (shape <= 0 | scale <= 0) {
    stop("Shape or scale parameter negative in dinvgamma().\n")
  }
  alpha <- shape
  beta <- scale
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  if (log){
    return(log.density)
  } else {
    return(log.density)
  }

}


wLogisticUpdate <- function(y, w, mu, sigma, proposal.sd){
  
  
  w.new <- w
  accept <- 0
  
  for(i in 1:length(w)){
    
    w.prop <- rnorm(1, w[i], proposal.sd)
    lp.curr <- w[i] + mu[i]
    lp.prop <- w.prop + mu[i]
    
    like.curr <- dbinom(y[i], 1, expit(lp.curr), log=T)
    like.prop <- dbinom(y[i], 1, expit(lp.prop), log=T)
    like.diff <- like.prop - like.curr
    
    # prior.sd <- 1#sqrt(sigma[1])
    # prior.curr <- dnorm(w[i], 0, prior.sd, log=T)
    # prior.prop <- dnorm(w.prop, 0, prior.sd, log=T)
    # prior.diff <- prior.prop - prior.curr
    
    # prior mean
    sigma_sub <- subdiv_sigmas[[i]]
    ids_sub <- subdiv_ids[[i]]
    w_sub <- w[subdiv_ids[[i]]]
    i_char <- as.character(i)
    i_not <- as.character(ids_sub[!ids_sub == i])
    
    sigma12 <- matrix(sigma_sub[i_char, i_not], nrow=1)
    sigma22 <- sigma[i_not, i_not]
    sigma22.inv <- solve(sigma22)
    w_minus <- w[i_not]
    prior_mean <- sigma12 %*% sigma22.inv %*% w_minus
    
    # prior variance
    sigma11 <- sigma_sub[i_char, i_char]
    sigma21 <- matrix(sigma_sub[i_not, i_char], ncol=1)
    prior_var <- sigma11 - sigma12 %*% sigma22.inv %*% sigma21
    
    prior.curr <- dnorm(w[i], prior_mean, sqrt(prior_var), log=T)
    prior.prop <- dnorm(w.prop, prior_mean, sqrt(prior_var), log=T)
    prior.diff <- prior.prop - prior.curr
    
    acceptance <- exp(like.diff + prior.diff)
    if(runif(1) <= acceptance) {
      w.new[i] <- w.prop
      accept <- accept + 1
    }
    
  }
  
  out <- list()
  out$w <- w.new
  out$accept <- accept
  return(out)
  
  
}


rangeMVGPupdate <- function(h.i, t.i, w.i, d, theta.i, proposal.sd, prior){
  
  w.i <- as.numeric(w.i)
  
  # proposal
  theta.new <- rlnorm(1, meanlog=log(theta.i), sdlog=proposal.sd)
  q.new <- dlnorm(theta.new, meanlog=log(theta.i), sdlog=proposal.sd, log=T)
  q.old <- dlnorm(theta.i, meanlog=log(theta.new), sdlog=proposal.sd, log=T)
  
  # covariance matrices
  h.new <- Exponential(d, range=theta.new, phi=1)
  
  # likelihoods
  loglik.curr <- dmvnorm(w.i, sigma=kronecker(h.i, t.i), log=T)
  loglik.new <- dmvnorm(w.i, sigma=kronecker(h.new, t.i), log=T)
  like.diff <- loglik.new - loglik.curr
  
  # priors
  shape <- prior[1]
  scale <- prior[2]
  prior.curr <- dgamma(theta.i, shape=shape, scale=scale, log=T)
  prior.new <- dgamma(theta.new, shape=shape, scale=scale, log=T)
  prior.diff <- prior.new - prior.curr
  
  out <- list()
  acceptance <- exp(like.diff + prior.diff + q.old - q.new)
  if(runif(1) <= acceptance) {
    out$theta <- theta.new
    out$accept <- 1
  } else { 
    out$theta <- theta.i
    out$accept <- 0
  }
  
  return(out)
  
}


wPoissonUpdate <- function(y, w, mu, sigma, proposal.sd){
  
  w.new <- w
  accept <- 0
  
  for(i in 1:length(w)){
    
    w.prop <- rnorm(1, w.new[i], proposal.sd)
    
    lp.curr <- w.new[i] + mu[i]
    lp.prop <- w.prop + mu[i]
    
    like.curr <- dpois(y[i], exp(lp.curr), log=T)
    like.prop <- dpois(y[i], exp(lp.prop), log=T)
    like.diff <- like.prop - like.curr
    
    prior.sd <- sqrt(sigma[1])
    prior.curr <- dnorm(w[i], 0, prior.sd, log=T)
    prior.prop <- dnorm(w.prop, 0, prior.sd, log=T)
    prior.diff <- prior.prop - prior.curr
    
    acceptance <- exp(like.diff + prior.diff)
    if(runif(1) <= acceptance) {
      w.new[i] <- w.prop
      accept <- accept + 1
    }
    
  }
  
  out <- list()
  out$w <- w.new
  out$accept <- accept
  return(out)
  
}


rangeMhUpdate <- function(theta, w, D, phi, proposal.sd, a, b){
  
  # proposal
  theta.new <- rlnorm(1, meanlog=log(theta), sdlog=proposal.sd)
  q.new <- dlnorm(theta.new, meanlog=log(theta), sdlog=proposal.sd, log=T)
  q.old <- dlnorm(theta, meanlog=log(theta.new), sdlog=proposal.sd, log=T)
  
  # covariance matrices
  sigma.curr <- Exponential(D, range=theta, phi=phi)
  sigma.new <- Exponential(D, range=theta.new, phi=phi)
  
  # likelihoods
  loglik.curr <- dmvnorm(w, sigma=sigma.curr, log=T)
  loglik.new <- dmvnorm(w, sigma=sigma.new, log=T)
  like.diff <- loglik.new - loglik.curr
  
  # priors
  prior.curr <- dgamma(theta, shape=a, scale=b, log=T)
  prior.new <- dgamma(theta.new, shape=a, scale=b, log=T)
  prior.diff <- prior.new - prior.curr
  
  out <- list()
  acceptance <- exp(like.diff + prior.diff + q.old - q.new)
  if(runif(1) <= acceptance) {
    out$theta <- theta.new
    out$accept <- 1
  } else { 
    out$theta <- theta
    out$accept <- 0
  }
  
  return(out)
  
}


muPoissonUpdate <- function(y, mu, w, proposal.sd){
  
  mu.prop <- rnorm(n=1, mu, proposal.sd)
  like.curr <- 0
  like.prop <- 0
  
  for (i in 1:length(y)){
    
    lp.curr <- w[i] + mu
    lp.prop <- w[i] + mu.prop
    like.curr <- like.curr + dpois(y[i], exp(lp.curr), log=T)
    like.prop <- like.prop + dpois(y[i], exp(lp.prop), log=T)
    
  }
  
  # prior? 
  like.diff <- like.prop - like.curr
  acceptance <- exp(like.diff)
  out <- list()
  if(runif(1) <= acceptance) {
    out$mu <- mu.prop
    out$accept <- 1
  } else {
    out$mu <- mu
    out$accept <- 0
  }
  
  return(out)
  
}


gpMCMC <- function(y, d, prior.theta, prior.phi, n_iter=50000, burnin=5000,
                   proposal.sd.w=0.1, proposal.sd.mu=0.1, proposal.sd.theta=0.5){
  
  nsamp <- length(y)
  n.keep <- n_iter - burnin
  samples.mu <- array(NA, n.keep)
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.w <- array(NA, c(n.keep, nsamp))
  
  theta.i <- runif(1, 1, 5)
  phi.i <- runif(1, 1, 5)
  mu.i <- runif(1, 2, 5)
  w.i <- mvrnorm(n=1, mu=rep(0, nsamp), diag(rep(1, nsamp)))
  
  accept <- c(0, 0, 0)
  cat("Generating", n.keep, "post burnin samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n_iter)
  
  for (i in 1:n_iter){
    
    # update mu
    mu.out <- muPoissonUpdate(y, mu.i, w.i, proposal.sd.mu)
    mu.i <- mu.out$mu
    accept[3] <- accept[3] + mu.out$accept
    
    # update w(s)
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    w.out <- wPoissonUpdate(y, w.i, mu.i, sigma.i, proposal.sd.w)
    w.i <- w.out$w
    w.i <- w.i - mean(w.i)
    accept[1] <- accept[1] + w.out$accept
    
    # update range
    theta.out <- rangeMhUpdate(theta.i, w.i, d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
    theta.i <- theta.out$theta
    accept[2] <- accept[2] + theta.out$accept
    
    # update tau2
    R.i <- sigma.i/phi.i
    phi.i <- 1 / rgamma(1, nsamp/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
    
    if (i > burnin){
      
      j <- i - burnin
      samples.mu[j] <- mu.i
      samples.w[j,] <- w.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      
      #k <- i/100
      #if(ceiling(k)==floor(k)){
      #  proposal.sd.w <- tuners.sd.1(accept[1], i*nsamp, proposal.sd.w, 75, 80)
      #  proposal.sd.mu <- tuners.sd.1(accept[3], i, proposal.sd.mu, 40, 50)
      #}
      
    }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n_iter)
    }
    
  }
  
  accept <- accept/n.keep
  accept[1] <- accept[1]/nsamp
  
  out <- list()
  out$accept <- accept
  out$samples.mu <- samples.mu
  out$samples.theta <- samples.theta
  out$samples.phi <- samples.phi
  out$samples.w <- samples.w
  
  return(out)
  
}


summarizeGp <- function(output){
  
  mu.hat <- mean(output$samples.mu)
  mu.pbias <- pbias(Mu, mu.hat)
  
  theta.hat <- mean(output$samples.theta)
  theta.pbias <- pbias(Theta, theta.hat)
  
  phi.hat <- mean(output$samples.phi)
  phi.pbias <- pbias(Phi, phi.hat)
  
  w.hat <- colMeans(output$samples.w)
  w.pbias <- mean(pbias(W, w.hat))
  
  return(list(
    mu.pbias=mu.pbias,
    theta.pbias=theta.pbias,
    phi.pbias=phi.pbias,
    w.mpbias=w.pbias
  ))
  
}


# for each posterior w sample in output
# generate realization at new points
# return those samples, and their means
krigeW <- function(output, d, ids){
  
  n.new <- nrow(d) - length(ids)
  nsamp.old <- nrow(output$samples.w)
  samples.new <- array(NA, c(nrow(output$samples.w), n.new))
  cat("generating", nsamp.old, "kriged samples\n")
  progressBar <- txtProgressBar(style=3)
  percentage.points<-round((1:100/100)*nsamp.old)
  
  for (i in 1:nsamp.old){
    
    theta.i <- output$samples.theta[i]
    phi.i <- output$samples.phi[i]
    covmat.i <- Exponential(d, range=theta.i, phi=phi.i)
    
    omega.11 <- covmat.i[-ids, -ids]
    omega.12 <- covmat.i[-ids, ids]
    omega.21 <- covmat.i[ids, -ids]
    omega.22 <- covmat.i[ids, ids]
    
    w.i <- output$samples.w[i,]
    e.cond <- omega.12 %*% solve(omega.22) %*% w.i
    var.cond <- omega.11 - omega.12 %*% solve(omega.22) %*% omega.21
    w.pred <- mvrnorm(n=1, mu=e.cond, var.cond)
    samples.new[i,] <- w.pred
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/nsamp.old)
    }
    
  }
  
  out <- list(
    mu.new=colMeans(samples.new),
    samples.new=samples.new
  )
  
  return(out)
  
}
