library(Matrix)


###########################################
## fit the Leroux et al. (2000) CAR model
## Metropolis Hastings random walk solution

## carLerouxIntegration
#   Fits the integrated Leroux CAR model
#   on multiple simulated datasets. Uses
#   Uses separate CAR models for each
#   dataset. No thinning.

## carLerouxIntegrationThinned
#   Fits the integrated Leroux CAR model 
#   on multiple thinned datasets. Uses 
#   separate CAR models for each dataset.

## carLerouxIntegrationUnif
#   Fits the integrated Leroux CAR model 
#   on multiple simulated datasets, also
#   integrates random survey counts. No
#   thinning.

## carLerouxThinnedUnif
#   Fits the integrated Leroux CAR model 
#   on multiple thinned simulated datasets,
#   also integrates cell counts from 
#   unthinned surveys.
###########################################


######################
# carLerouxIntegration
######################
#   Fits the integrated Leroux CAR model on multiple simulated datasets. 
#   Uses separate CAR models for each dataset. No thinning.
## arguments
#   data: list of data objects (outputs of simLeroux)
#   proposal.sd.beta: beta proposal distribution standard deviation 
#   proposal.sd.phi: phi proposal distribution standard deviation
#   n.sample: number of mcmc samples, not including burning
#   prior.var.beta: variance of the beta prior distribution
#   beta.initial: initial mcmc beta value
#   phi.initial: initial mcmc phi value
#   tau2.initial: initial mcmc tau2 value
carLerouxIntegration <- function(data, proposal.sd.beta, proposal.sd.phi, proposal.sd.rho=NULL, n.sample=125000,
                                 prior.var.beta=NULL, prior.tau2=NULL,
                                 beta.initial=NULL, phi.initial=NULL, tau2.initial=NULL, rho.initial=NULL,
                                 fix.rho=TRUE, rho.list=NULL, self.tune=FALSE){
  
  
  if (fix.rho & is.null(rho.list)) stop("you must supply rho.list if fix.rho=TRUE")
  
  
  #### setup
  D <- length(data)
  K <- length(data[[1]]$y)
  p <- data[[1]]$p
  X.list <- lapply(data, function(z){ z$x.standardised})
  Y.list <- lapply(data, function(z){ z$y})
  X.all <- X.list[[1]]
  Y.all <- Y.list[[1]]
  if (D > 1){
    for (d in 2:D){ 
      X.all <- rbind(X.all, X.list[[d]])
      Y.all <- c(Y.all, Y.list[[d]])
    }
  }

  
  #### priors
  prior.mean.beta <- rep(0, p)
  if (is.null(prior.var.beta)){ prior.var.beta <- rep(1000, p) }
  if (is.null(prior.tau2)){ prior.tau2 <- c(1, 0.01) }


  #### Initial values
  mod.glm <- glm(Y.all~X.all-1, family="poisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  if (is.null(beta.initial)){
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  } else {
    beta <- beta.initial
  }

  
  beta.initial <- beta
  phi.list <- list()
  tau2.list <- list()
  for (d in 1:D){
    Y.d <- Y.list[[d]]
    X.d <- X.list[[d]]
    log.Y <- log(Y.d)
    log.Y[Y.d==0] <- -0.1  
    res.temp <- log.Y - X.d %*% beta.mean
    res.sd <- sd(res.temp, na.rm=TRUE)
    if (is.null(phi.initial)){
      phi.d <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    } else {
      phi.d <- phi.initial[d]
    }
    if (is.null(tau2.initial)){
      tau2.d <- var(phi.d)/8
    } else {
      tau2.d <- tau2.initial[d]
    }
    phi.list[[d]] <- phi.d
    tau2.list[[d]] <- tau2.d
  }
  if (!fix.rho){
    rho.list <- lapply(1:D, function(x) {runif(1)})
  }
  
  
  beta.initial <- beta
  phi.initials <- phi.list
  tau2.initials <- tau2.list
  
  
  #### MCMC quantities    
  burnin <- 10000
  thin <- 1
  
  ## Storage
  n.keep <- floor((n.sample - burnin)/ thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.tau2 <- array(NA, c(n.keep, D))
  samples.phi <- array(NA, c(D, n.keep, K))
  if (!fix.rho){ samples.rho <- array(NA, c(D, n.keep, 1))}
  samples.deviance <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 2)
  if (!fix.rho){ accept <- c(accept, 0)}
  tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
  
  ## CAR quantities
  # ToDo: deal with rho
  W.quants <- common.Wcheckformat.leroux(W, K, fix.rho=FALSE, rho.1)
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  if (!fix.rho){
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.func <- function(x){
      0.5 * sum(log((x * Wstar.val + (1-x))))
    }
    det.Q.list <- lapply(rho.list, det.func)
  }
  
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta
    offsets.beta <- phi.list
    beta.out <- betaIntegratedUpdate(X.list, Y.list, beta, offsets.beta,
                                     prior.mean.beta, prior.var.beta, proposal.sd.beta)
    beta <- beta.out$beta
    
    
    ## Sample from phi
    phi.accept <- 0
    for (d in 1:D){
      y.d <- Y.list[[d]]
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      tau2.d <- tau2.list[[d]]
      offsets.phi.d <- X.list[[d]] %*% beta
      phi.out.d <- phiPoissonUpdate(y.d, W, nNeighbors, phi.d, tau2.d, offsets.phi.d, proposal.sd.phi, rho=rho.d)
      phi.list[[d]] <- phi.out.d$phi - mean(phi.out.d$phi)
      phi.accept <- phi.accept + phi.out.d$accept
    }
    
    
    ## Sample from tau2
    temp.list <- list()
    for (d in 1:D){
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, rho.d)
      temp.list[[d]] <- temp2
      tau2.posterior.scale <- temp2 + prior.tau2[2]
      tau2.d <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
      tau2.list[[d]] <- tau2.d
    }
    
    
    ## sample from rho
    if(!fix.rho){
      rho.accept <- 0
      for (d in 1:D){
        rho.d <- rho.list[[d]]
        phi.d <- phi.list[[d]]
        tau2.d <- tau2.list[[d]]
        temp2.d <- temp.list[[d]]
        det.Q.d <- det.Q.list[[d]]
        proposal.rho.d <- rtruncnorm(n=1, a=0, b=1, mean=rho.d, sd=proposal.sd.rho)
        temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, proposal.rho.d)
        det.Q.proposal <- det.func(proposal.rho.d)
        logprob.current <- det.Q.d - temp2.d / tau2.d
        logprob.proposal <- det.Q.proposal - temp3 / tau2.d
        prob <- exp(logprob.proposal - logprob.current)
      
        if (prob > runif(1)){
          rho.list[[d]] <- proposal.rho.d
          det.Q.list[[d]] <- det.Q.proposal
          rho.accept <- rho.accept + 1           
        }           
      }
    }
    
    
    ## Calculate the deviance
    lp <- as.numeric(X.all %*% beta) + unlist(phi.list)
    fitted <- exp(lp)
    deviance.all <- dpois(x=as.numeric(Y.all), lambda=fitted, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE) 
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta[ele,] <- beta
      accept[1] <- accept[1] + beta.out$accept
      accept[2] <- accept[2] + phi.accept
      samples.deviance[ele,] <- deviance
      for (d in 1:D){
        samples.phi[d,ele,] <- phi.list[[d]]
        samples.tau2[ele,][d] <- tau2.list[[d]]
        if (!fix.rho){
          samples.rho[d,ele,] <- rho.list[[d]]
          accept[3] <- accept[3] + rho.accept
        }
      }
      
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p>2){
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 40, 50)
          } else {
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 30, 40)
          }
          proposal.sd.phi <- tuners.sd.1(accept[2], ele*D*K, proposal.sd.phi, 55, 65)
          if(!fix.rho){
            proposal.sd.rho <- tuners.sd.2(accept[3], ele*D, proposal.sd.rho, 40, 50, 0.75)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  accept[2] <- accept[2]/(D*K)
  if (!fix.rho){ accept[3] <- accept[3]/D}
  
  
  #### Deviance information criterion (DIC)
  median.beta <- apply(samples.beta, 2, median)
  median.phi <- unlist(lapply(1:D, function(x) {apply(samples.phi[x,,], 2, median)}))
  fitted.median <- exp(X.all %*% median.beta + median.phi)
  deviance.fitted <- -2 * sum(dpois(x=Y.all, lambda=fitted.median, log=TRUE), na.rm=TRUE)
  DIC <- 2 * median(samples.deviance) - deviance.fitted
  
  
  #### Compute the % deviance explained
  fit.null <- glm(Y.all~1, family="poisson")$fitted.values
  deviance.null <- -2 * sum(dpois(x=Y.all, lambda=fit.null, log=TRUE), na.rm=TRUE)
  percent_dev_explained <- 100 * (deviance.null - deviance.fitted) / deviance.null

    
  out <- list()
  out$samples.beta <- samples.beta
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$accept <- accept
  out$beta.initial <- beta.initial
  out$phi.initial <- phi.initials
  out$tau2.initial <- tau2.initials
  out$samples.deviance <- samples.deviance
  out$fit <- list(DIC=DIC, percent_dev=percent_dev_explained)
  out$final.tune.beta <- proposal.sd.beta
  out$final.tune.phi <- proposal.sd.phi
  if (!fix.rho){
    out$final.tune.rho <- proposal.sd.rho
    out$samples.rho <- samples.rho
  }
  return(out)
  
  
}


##########################
# carLerouxIntegrationUnif
##########################
#   Fits the integrated Leroux CAR model on multiple simulated datasets,
#   and integrates a set of survey counts from random locations.
#   Uses separate CAR models for each dataset. No thinning.
## arguments
#   data: nested list with elements 'spatial' and 'survey'
#   proposal.sd.beta: beta proposal distribution standard deviation 
#   proposal.sd.phi: phi proposal distribution standard deviation
#   n.sample: number of mcmc samples, not including burning
#   prior.var.beta: variance of the beta prior distribution
#   beta.initial: initial mcmc beta value
#   phi.initial: initial mcmc phi value
#   tau2.initial: initial mcmc tau2 value
carLerouxIntegrationUnif <- function(data, proposal.sd.beta, proposal.sd.phi, proposal.sd.rho=NULL, n.sample=125000,
                                     prior.var.beta=NULL, prior.tau2=NULL,
                                     beta.initial=NULL, phi.initial=NULL, tau2.initial=NULL, 
                                     rho.initial=NULL, fix.rho=TRUE, rho.list=NULL, self.tune=FALSE){
  
  if (fix.rho & is.null(rho.list)) stop("you must supply rho.list if fix.rho=TRUE")
  
  #### setup
  spatial <- data$spatial
  survey <- data$survey
  S <- length(data$survey$y)
  D <- length(spatial)
  K <- length(spatial[[1]]$y)
  p <- spatial[[1]]$p
  X.list.sp <- lapply(spatial, function(z){ z$x.standardised})
  Y.list.sp <- lapply(spatial, function(z){ z$y})
  X.list <- list(spatial=X.list.sp, survey=survey$x.standardised)
  Y.list <- list(spatial=Y.list.sp, survey=survey$y)
  X.all <- X.list$survey
  Y.all <- Y.list$survey
  for (d in 1:D){ 
    X.all <- rbind(X.all, X.list$spatial[[d]])
    Y.all <- c(Y.all, Y.list$spatial[[d]])
  }
  
  
  #### priors
  prior.mean.beta <- rep(0, p)
  if (is.null(prior.var.beta)){ prior.var.beta <- rep(1000, p) }
  if (is.null(prior.tau2)){ prior.tau2 <- c(1, 0.01) }
  
  
  #### Initial values
  mod.glm <- glm(Y.all~X.all-1, family="poisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  if (is.null(beta.initial)){
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  } else {
    beta <- beta.initial
  }
  
  
  beta.initial <- beta
  phi.list <- list()
  tau2.list <- list()
  for (d in 1:D){
    Y.d <- Y.list$spatial[[d]]
    X.d <- X.list$spatial[[d]]
    log.Y <- log(Y.d)
    log.Y[Y.d==0] <- -0.1  
    res.temp <- log.Y - X.d %*% beta.mean
    res.sd <- sd(res.temp, na.rm=TRUE)
    if (is.null(phi.initial)){
      phi.d <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    } else {
      phi.d <- phi.initial[d]
    }
    if (is.null(tau2.initial)){
      tau2.d <- var(phi.d)/8
    } else {
      tau2.d <- tau2.initial[d]
    }
    phi.list[[d]] <- phi.d
    tau2.list[[d]] <- tau2.d
  }
  if (!fix.rho){
    rho.list <- lapply(1:D, function(x) {runif(1)})
  }
  
  
  beta.initial <- beta
  phi.initials <- phi.list
  tau2.initials <- tau2.list
  
  
  #### MCMC quantities    
  burnin <- 10000
  thin <- 1
  
  ## Storage
  n.keep <- floor((n.sample - burnin)/ thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.tau2 <- array(NA, c(n.keep, D))
  samples.phi <- array(NA, c(D, n.keep, K))
  if (!fix.rho){ samples.rho <- array(NA, c(D, n.keep, 1))}
  samples.deviance <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 2)
  if (!fix.rho){ accept <- c(accept, 0)}
  tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
  
  ## CAR quantities
  # ToDo: deal with rho
  W.quants <- common.Wcheckformat.leroux(W, K, fix.rho=FALSE, rho.1)
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  if (!fix.rho){
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.func <- function(x){
      0.5 * sum(log((x * Wstar.val + (1-x))))
    }
    det.Q.list <- lapply(rho.list, det.func)
  }
  
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta
    offsets.beta <- phi.list
    beta.out <- betaIntegratedUnifUpdate(X.list, Y.list, beta, offsets.beta,
                                         prior.mean.beta, prior.var.beta, proposal.sd.beta)
    beta <- beta.out$beta
    
    
    ## Sample from phi
    phi.accept <- 0
    for (d in 1:D){
      y.d <- Y.list$spatial[[d]]
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      tau2.d <- tau2.list[[d]]
      offsets.phi.d <- X.list$spatial[[d]] %*% beta
      phi.out.d <- phiPoissonUpdate(y.d, W, nNeighbors, phi.d, tau2.d, offsets.phi.d, proposal.sd.phi, rho=rho.d)
      phi.list[[d]] <- phi.out.d$phi - mean(phi.out.d$phi)
      phi.accept <- phi.accept + phi.out.d$accept
    }
    
    
    ## Sample from tau2
    temp.list <- list()
    for (d in 1:D){
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, rho.d)
      temp.list[[d]] <- temp2
      tau2.posterior.scale <- temp2 + prior.tau2[2]
      tau2.d <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
      tau2.list[[d]] <- tau2.d
    }
    
    
    ## sample from rho
    if(!fix.rho){
      rho.accept <- 0
      for (d in 1:D){
        rho.d <- rho.list[[d]]
        phi.d <- phi.list[[d]]
        tau2.d <- tau2.list[[d]]
        temp2.d <- temp.list[[d]]
        det.Q.d <- det.Q.list[[d]]
        proposal.rho.d <- rtruncnorm(n=1, a=0, b=1, mean=rho.d, sd=proposal.sd.rho)
        temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, proposal.rho.d)
        det.Q.proposal <- det.func(proposal.rho.d)
        logprob.current <- det.Q.d - temp2.d / tau2.d
        logprob.proposal <- det.Q.proposal - temp3 / tau2.d
        prob <- exp(logprob.proposal - logprob.current)
        
        if (prob > runif(1)){
          rho.list[[d]] <- proposal.rho.d
          det.Q.list[[d]] <- det.Q.proposal
          rho.accept <- rho.accept + 1
        }
      }
    }
    
    
    ## Calculate the deviance
    lp <- as.numeric(X.all %*% beta) + c(rep(0, S), unlist(phi.list))
    fitted <- exp(lp)
    deviance.all <- dpois(x=as.numeric(Y.all), lambda=fitted, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE) 
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta[ele,] <- beta
      accept[1] <- accept[1] + beta.out$accept
      accept[2] <- accept[2] + phi.accept
      samples.deviance[ele,] <- deviance
      for (d in 1:D){
        samples.phi[d,ele,] <- phi.list[[d]]
        samples.tau2[ele,][d] <- tau2.list[[d]]
        if (!fix.rho){
          samples.rho[d,ele,] <- rho.list[[d]]
          accept[3] <- accept[3] + rho.accept
        }
      }
      
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p>2){
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 40, 50)
          } else {
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 30, 40)
          }
          proposal.sd.phi <- tuners.sd.1(accept[2], ele*D*K, proposal.sd.phi, 55, 65)
          if(!fix.rho){
            proposal.sd.rho <- tuners.sd.2(accept[3], ele*D, proposal.sd.rho, 40, 50, 1)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  accept[2] <- accept[2]/(D*K)
  if (!fix.rho){ accept[3] <- accept[3]/D}
  
  
  #### Deviance information criterion (DIC)
  median.beta <- apply(samples.beta, 2, median)
  median.phi <- unlist(lapply(1:D, function(x) {apply(samples.phi[x,,], 2, median)}))
  median.lp <- as.numeric(X.all %*% median.beta) + c(rep(0, S), unlist(median.phi))
  fitted.median <- exp(median.lp)
  deviance.fitted <- -2 * sum(dpois(x=Y.all, lambda=fitted.median, log=TRUE), na.rm=TRUE)
  DIC <- 2 * median(samples.deviance) - deviance.fitted
  
  
  #### Compute the % deviance explained
  fit.null <- glm(Y.all~1, family="poisson")$fitted.values
  deviance.null <- -2 * sum(dpois(x=Y.all, lambda=fit.null, log=TRUE), na.rm=TRUE)
  percent_dev_explained <- 100 * (deviance.null - deviance.fitted) / deviance.null
  
  
  out <- list()
  out$samples.beta <- samples.beta
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$accept <- accept
  out$beta.initial <- beta.initial
  out$phi.initial <- phi.initials
  out$tau2.initial <- tau2.initials
  out$samples.deviance <- samples.deviance
  out$fit <- list(DIC=DIC, percent_dev=percent_dev_explained)
  out$final.tune.beta <- proposal.sd.beta
  out$final.tune.phi <- proposal.sd.phi
  if (!fix.rho){
    out$final.tune.rho <- proposal.sd.rho
    out$samples.rho <- samples.rho
  }
  return(out)
  
  
}


#############################
# carLerouxIntegrationThinned
#############################
#   Fits the integrated Leroux CAR model on multiple thinned simulated datasets. 
#   Uses separate CAR models for each dataset.
## arguments
#   data: list of outputs from simLerouxThinned
#   proposal.sd.beta: beta proposal distribution standard deviation 
#   proposal.sd.phi: phi proposal distribution standard deviation
#   n.sample: number of mcmc samples, not including burning
#   prior.var.beta: variance of the beta prior distribution
#   beta.initial: initial mcmc beta value
#   phi.initial: initial mcmc phi value
#   tau2.initial: initial mcmc tau2 value
carLerouxIntegrationThinned <- function(data, proposal.sd.beta, proposal.sd.phi, proposal.sd.rho=NULL,
                                        prior.var.beta=NULL, prior.tau2=NULL, beta.initial=NULL, phi.initial=NULL, 
                                        tau2.initial=NULL, rho.initial=NULL, fix.rho=TRUE, rho.list=NULL, n.sample=NULL, self.tune=FALSE){
  
  
  if (fix.rho & is.null(rho.list)) stop("you must supply rho.list if fix.rho=TRUE")
  
  
  #### setup
  D <- length(data)
  K <- length(data[[1]]$y)
  p <- data[[1]]$p
  X.list <- list()
  for (d in 1:D){
    p <- p + data[[d]]$p.z
    x.d <- data[[d]]$x.standardised
    z.d <- data[[d]]$z.standardised
    if (D > 1){
      if (d==1){
        x.d <- cbind(x.d, z.d, rep(0, nrow(x.d)))
      } else {
        x.d <- cbind(x.d, rep(0, nrow(x.d)), z.d)
      }
    }
    X.list[[d]] <- x.d
  }
  Y.list <- lapply(data, function(z){ z$y})
  X.all <- X.list[[1]]
  Y.all <- Y.list[[1]]
  if (D > 1){
    for (d in 2:D){ 
      X.all <- rbind(X.all, X.list[[d]])
      Y.all <- c(Y.all, Y.list[[d]])
    }
  }
  
  
  #### priors
  prior.mean.beta <- rep(0, p)
  if (is.null(prior.var.beta)){ prior.var.beta <- rep(1000, p) }
  if (is.null(prior.tau2)){ prior.tau2 <- c(1, 0.01) }
  
  
  #### Initial values
  mod.glm <- glm(Y.all~X.all-1, family="poisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  if (is.null(beta.initial)){
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  } else {
    beta <- beta.initial
  }
  
  
  beta.initial <- beta
  phi.list <- list()
  tau2.list <- list()
  for (d in 1:D){
    Y.d <- Y.list[[d]]
    X.d <- X.list[[d]]
    log.Y <- log(Y.d)
    log.Y[Y.d==0] <- -0.1  
    res.temp <- log.Y - X.d %*% beta.mean
    res.sd <- sd(res.temp, na.rm=TRUE)
    if (is.null(phi.initial)){
      phi.d <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    } else {
      phi.d <- phi.initial[d]
    }
    if (is.null(tau2.initial)){
      tau2.d <- var(phi.d)/8
    } else {
      tau2.d <- tau2.initial[d]
    }
    phi.list[[d]] <- phi.d
    tau2.list[[d]] <- tau2.d
  }
  if (!fix.rho){
    rho.list <- lapply(1:D, function(x) {runif(1)})
  }
  
  
  beta.initial <- beta
  phi.initials <- phi.list
  tau2.initials <- tau2.list
  
  
  #### MCMC quantities    
  if (is.null(n.sample)) {n.sample <- 125000}
  burnin <- 10000
  thin <- 1
  
  ## Storage
  n.keep <- floor((n.sample - burnin)/ thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.tau2 <- array(NA, c(n.keep, D))
  samples.phi <- array(NA, c(D, n.keep, K))
  if (!fix.rho){ samples.rho <- array(NA, c(D, n.keep, 1))}
  samples.deviance <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 2)
  if (!fix.rho){ accept <- c(accept, 0)}
  tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
  
  ## CAR quantities
  # ToDo: deal with rho
  W.quants <- common.Wcheckformat.leroux(W, K, fix.rho=FALSE, rho.1)
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  if (!fix.rho){
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.func <- function(x){
      0.5 * sum(log((x * Wstar.val + (1-x))))
    }
    det.Q.list <- lapply(rho.list, det.func)
  }
  
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta
    offsets.beta <- phi.list
    beta.out <- betaIntegratedUpdate(X.list, Y.list, beta, offsets.beta,
                                     prior.mean.beta, prior.var.beta, proposal.sd.beta)
    beta <- beta.out$beta
    
    
    ## Sample from phi
    phi.accept <- 0
    for (d in 1:D){
      y.d <- Y.list[[d]]
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      tau2.d <- tau2.list[[d]]
      offsets.phi.d <- X.list[[d]] %*% beta
      phi.out.d <- phiPoissonUpdate(y.d, W, nNeighbors, phi.d, tau2.d, offsets.phi.d, proposal.sd.phi, rho=rho.d)
      phi.list[[d]] <- phi.out.d$phi - mean(phi.out.d$phi)
      phi.accept <- phi.accept + phi.out.d$accept
    }
    
    
    ## Sample from tau2
    temp.list <- list()
    for (d in 1:D){
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, rho.d)
      temp.list[[d]] <- temp2
      tau2.posterior.scale <- temp2 + prior.tau2[2]
      tau2.d <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
      tau2.list[[d]] <- tau2.d
    }
    
    
    ## sample from rho
    if(!fix.rho){
      rho.accept <- 0
      for (d in 1:D){
        rho.d <- rho.list[[d]]
        phi.d <- phi.list[[d]]
        tau2.d <- tau2.list[[d]]
        temp2.d <- temp.list[[d]]
        det.Q.d <- det.Q.list[[d]]
        proposal.rho.d <- rtruncnorm(n=1, a=0, b=1, mean=rho.d, sd=proposal.sd.rho)
        temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, proposal.rho.d)
        det.Q.proposal <- det.func(proposal.rho.d)
        logprob.current <- det.Q.d - temp2.d / tau2.d
        logprob.proposal <- det.Q.proposal - temp3 / tau2.d
        prob <- exp(logprob.proposal - logprob.current)
        
        if (prob > runif(1)){
          rho.list[[d]] <- proposal.rho.d
          det.Q.list[[d]] <- det.Q.proposal
          rho.accept <- rho.accept + 1
        }
      }
    }
    
    
    ## Calculate the deviance
    lp <- as.numeric(X.all %*% beta) + unlist(phi.list)
    fitted <- exp(lp)
    deviance.all <- dpois(x=as.numeric(Y.all), lambda=fitted, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE)  
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta[ele,] <- beta
      accept[1] <- accept[1] + beta.out$accept
      accept[2] <- accept[2] + phi.accept
      samples.deviance[ele,] <- deviance
      for (d in 1:D){
        samples.phi[d,ele,] <- phi.list[[d]]
        samples.tau2[ele,][d] <- tau2.list[[d]]
        if (!fix.rho){
          samples.rho[d,ele,] <- rho.list[[d]]
          accept[3] <- accept[3] + rho.accept
        }
      }
      
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p>2){
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 40, 50)
          } else {
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 30, 40)
          }
          proposal.sd.phi <- tuners.sd.1(accept[2], ele*D*K, proposal.sd.phi, 55, 65)
          if(!fix.rho){
            proposal.sd.rho <- tuners.sd.2(accept[3], ele*D, proposal.sd.rho, 40, 50, 0.75)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  accept[2] <- accept[2]/(D*K)
  if (!fix.rho){ accept[3] <- accept[3]/D}
  
  
  #### Deviance information criterion (DIC)
  median.beta <- apply(samples.beta, 2, median)
  median.phi <- unlist(lapply(1:D, function(x) {apply(samples.phi[x,,], 2, median)}))
  fitted.median <- exp(X.all %*% median.beta + median.phi)
  deviance.fitted <- -2 * sum(dpois(x=Y.all, lambda=fitted.median, log=TRUE), na.rm=TRUE)
  DIC <- 2 * median(samples.deviance) - deviance.fitted 
  
  
  #### Compute the % deviance explained
  fit.null <- glm(Y.all~1, family="poisson")$fitted.values
  deviance.null <- -2 * sum(dpois(x=Y.all, lambda=fit.null, log=TRUE), na.rm=TRUE)
  percent_dev_explained <- 100 * (deviance.null - deviance.fitted) / deviance.null
  
  
  out <- list()
  out$samples.beta <- samples.beta
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$accept <- accept
  out$beta.initial <- beta.initial
  out$phi.initial <- phi.initials
  out$tau2.initial <- tau2.initials
  out$samples.deviance <- samples.deviance
  out$fit <- list(DIC=DIC, percent_dev=percent_dev_explained)
  out$final.tune.beta <- proposal.sd.beta
  out$final.tune.phi <- proposal.sd.phi
  if (!fix.rho){
    out$final.tune.rho <- proposal.sd.rho
    out$samples.rho <- samples.rho
  }
  return(out)
  
  
}


######################
# carLerouxThinnedUnif
######################
#   Fits the integrated Leroux CAR model on multiple thinned simulated datasets.
#   Also integrates cell counts from unthinned surveys.
#   Uses separate CAR models for each dataset.
## arguments
#   data: nested list with elements 'spatial' (output from simLerouxThinned)
#     and 'survey' (output of simUniform)
#   proposal.sd.beta: beta proposal distribution standard deviation 
#   proposal.sd.phi: phi proposal distribution standard deviation
#   n.sample: number of mcmc samples, not including burning
#   prior.var.beta: variance of the beta prior distribution
#   beta.initial: initial mcmc beta value
#   phi.initial: initial mcmc phi value
#   tau2.initial: initial mcmc tau2 value
carLerouxThinnedUnif <- function(data, proposal.sd.beta, proposal.sd.phi, proposal.sd.rho=NULL,
                                 prior.var.beta=NULL, prior.tau2=NULL,
                                 beta.initial=NULL, phi.initial=NULL, 
                                 tau2.initial=NULL, n.sample=NULL, rho.initial=NULL,
                                 fix.rho=TRUE, rho.list=NULL, self.tune=FALSE){

  if (fix.rho & is.null(rho.list)) stop("you must supply rho.list if fix.rho=TRUE")  

  #### setup
  spatial <- data$spatial
  survey <- data$survey
  S <- length(survey$y)
  D <- length(spatial)
  K <- length(spatial[[1]]$y)
  p <- spatial[[1]]$p
  p.z <- 0
  X.list.sp <- list()
  for (d in 1:D){
    p <- p + spatial[[d]]$p.z
    p.z <- p.z + spatial[[d]]$p.z
    x.d <- spatial[[d]]$x.standardised
    z.d <- spatial[[d]]$z.standardised
    if (D > 1){
      if (d==1){
        x.d <- cbind(x.d, z.d, rep(0, nrow(x.d)))
      } else {
        x.d <- cbind(x.d, rep(0, nrow(x.d)), z.d)
      }
    } else {
      x.d <- cbind(x.d, z.d)
    }
    X.list.sp[[d]] <- x.d
  }
  Y.list.sp <- lapply(spatial, function(z){ z$y})
  X.list <- list(spatial=X.list.sp, survey=cbind(survey$x.standardised, array(0, c(S, p.z))))
  Y.list <- list(spatial=Y.list.sp, survey=survey$y)
  X.all <- X.list$survey
  Y.all <- Y.list$survey
  for (d in 1:D){
    X.all <- rbind(X.all, X.list.sp[[d]])
    Y.all <- c(Y.all, Y.list.sp[[d]])
  }
  
  
  #### priors
  prior.mean.beta <- rep(0, p)
  if (is.null(prior.var.beta)){ prior.var.beta <- rep(1000, p) }
  if (is.null(prior.tau2)){ prior.tau2 <- c(1, 0.01) }
  
  
  #### Initial values
  mod.glm <- glm(Y.all~X.all-1, family="poisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  if (is.null(beta.initial)){
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  } else {
    beta <- beta.initial
  }
  
  
  beta.initial <- beta
  phi.list <- list()
  tau2.list <- list()
  for (d in 1:D){
    Y.d <- Y.list$spatial[[d]]
    X.d <- X.list$spatial[[d]]
    log.Y <- log(Y.d)
    log.Y[Y.d==0] <- -0.1  
    res.temp <- log.Y - X.d %*% beta.mean
    res.sd <- sd(res.temp, na.rm=TRUE)
    if (is.null(phi.initial)){
      phi.d <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    } else {
      phi.d <- phi.initial[d]
    }
    if (is.null(tau2.initial)){
      tau2.d <- var(phi.d)/8
    } else {
      tau2.d <- tau2.initial[d]
    }
    phi.list[[d]] <- phi.d
    tau2.list[[d]] <- tau2.d
  }
  if (!fix.rho){
    rho.list <- lapply(1:D, function(x) {runif(1)})
  }
  
  
  beta.initial <- beta
  phi.initials <- phi.list
  tau2.initials <- tau2.list
  
  
  #### MCMC quantities    
  if (is.null(n.sample)) {n.sample <- 125000}
  burnin <- 10000
  thin <- 1
  
  ## Storage
  n.keep <- floor((n.sample - burnin)/ thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.tau2 <- array(NA, c(n.keep, D))
  samples.phi <- array(NA, c(D, n.keep, K))
  if (!fix.rho){ samples.rho <- array(NA, c(D, n.keep, 1))}
  samples.deviance <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 2)
  if (!fix.rho){ accept <- c(accept, 0)}
  tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
  
  ## CAR quantities
  # ToDo: deal with rho
  W.quants <- common.Wcheckformat.leroux(W, K, fix.rho=FALSE, rho.1)
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  if (!fix.rho){
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.func <- function(x){
      0.5 * sum(log((x * Wstar.val + (1-x))))
    }
    det.Q.list <- lapply(rho.list, det.func)
  }
  
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta
    offsets.beta <- phi.list
    beta.out <- betaIntegratedUnifUpdate(X.list, Y.list, beta, offsets.beta,
                                         prior.mean.beta, prior.var.beta, proposal.sd.beta)
    beta <- beta.out$beta
    
    
    ## Sample from phi
    phi.accept <- 0
    for (d in 1:D){
      y.d <- Y.list$spatial[[d]]
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      tau2.d <- tau2.list[[d]]
      offsets.phi.d <- X.list$spatial[[d]] %*% beta
      phi.out.d <- phiPoissonUpdate(y.d, W, nNeighbors, phi.d, tau2.d, offsets.phi.d, proposal.sd.phi, rho=rho.d)
      phi.list[[d]] <- phi.out.d$phi - mean(phi.out.d$phi)
      phi.accept <- phi.accept + phi.out.d$accept
    }
    
    
    ## Sample from tau2
    temp.list <- list()
    for (d in 1:D){
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, rho.d)
      temp.list[[d]] <- temp2
      tau2.posterior.scale <- temp2 + prior.tau2[2]
      tau2.d <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
      tau2.list[[d]] <- tau2.d
    }
    
    
    ## sample from rho
    if(!fix.rho){
      rho.accept <- 0
      for (d in 1:D){
        rho.d <- rho.list[[d]]
        phi.d <- phi.list[[d]]
        tau2.d <- tau2.list[[d]]
        temp2.d <- temp.list[[d]]
        det.Q.d <- det.Q.list[[d]]
        proposal.rho.d <- rtruncnorm(n=1, a=0, b=1, mean=rho.d, sd=proposal.sd.rho)
        temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi.d, phi.d, proposal.rho.d)
        det.Q.proposal <- det.func(proposal.rho.d)
        logprob.current <- det.Q.d - temp2.d / tau2.d
        logprob.proposal <- det.Q.proposal - temp3 / tau2.d
        prob <- exp(logprob.proposal - logprob.current)
        
        if (prob > runif(1)){
          rho.list[[d]] <- proposal.rho.d
          det.Q.list[[d]] <- det.Q.proposal
          rho.accept <- rho.accept + 1
        }
      }
    }
    
    
    ## Calculate the deviance
    lp <- as.numeric(X.all %*% beta) + c(rep(0, S), unlist(phi.list))
    fitted <- exp(lp)
    deviance.all <- dpois(x=as.numeric(Y.all), lambda=fitted, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE) 
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta[ele,] <- beta
      accept[1] <- accept[1] + beta.out$accept
      accept[2] <- accept[2] + phi.accept
      samples.deviance[ele,] <- deviance
      for (d in 1:D){
        samples.phi[d,ele,] <- phi.list[[d]]
        samples.tau2[ele,][d] <- tau2.list[[d]]
        if (!fix.rho){
          samples.rho[d,ele,] <- rho.list[[d]]
          accept[3] <- accept[3] + rho.accept
        }
      }
      
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p>2){
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 40, 50)
          } else {
            proposal.sd.beta <- tuners.sd.1(accept[1], ele, proposal.sd.beta, 30, 40)
          }
          proposal.sd.phi <- tuners.sd.1(accept[2], ele*D*K, proposal.sd.phi, 55, 65)
          if(!fix.rho){
            proposal.sd.rho <- tuners.sd.2(accept[3], ele*D, proposal.sd.rho, 40, 50, 0.75)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  accept[2] <- accept[2]/(D*K)
  if (!fix.rho){ accept[3] <- accept[3]/D}
  
  
  #### Deviance information criterion (DIC)
  median.beta <- apply(samples.beta, 2, median)
  median.phi <- unlist(lapply(1:D, function(x) {apply(samples.phi[x,,], 2, median)}))
  median.lp <- as.numeric(X.all %*% median.beta) + c(rep(0, S), unlist(median.phi))
  fitted.median <- exp(median.lp)
  deviance.fitted <- -2 * sum(dpois(x=Y.all, lambda=fitted.median, log=TRUE), na.rm=TRUE)
  DIC <- 2 * median(samples.deviance) - deviance.fitted
  
  
  #### Compute the % deviance explained
  fit.null <- glm(Y.all~1, family="poisson")$fitted.values
  deviance.null <- -2 * sum(dpois(x=Y.all, lambda=fit.null, log=TRUE), na.rm=TRUE)
  percent_dev_explained <- 100 * (deviance.null - deviance.fitted) / deviance.null
  
  
  out <- list()
  out$samples.beta <- samples.beta
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$accept <- accept
  out$beta.initial <- beta.initial
  out$phi.initial <- phi.initials
  out$tau2.initial <- tau2.initials
  out$samples.deviance <- samples.deviance
  out$fit <- list(DIC=DIC, percent_dev=percent_dev_explained)
  out$final.tune.phi <- proposal.sd.phi
  out$final.tune.beta <- proposal.sd.beta
  if (!fix.rho){
    out$final.tune.rho <- proposal.sd.rho
    out$samples.rho <- samples.rho
  }
  return(out)
  
  
}


############################
# carLerouxIntegrationShared
############################
#   Fits the integrated Leroux CAR model on multiple simulated datasets. 
#   Shares CAR parameters among datasets. No thinning.
## arguments
#   data: list of outputs from simLerouxThinned
#   proposal.sd.beta: beta proposal distribution standard deviation 
#   proposal.sd.phi: phi proposal distribution standard deviation
#   n.sample: number of mcmc samples, not including burning
#   prior.var.beta: variance of the beta prior distribution
#   beta.initial: initial mcmc beta value
#   phi.initial: initial mcmc phi value
#   tau2.initial: initial mcmc tau2 value
carLerouxIntegrationShared <- function(data, proposal.sd.beta, proposal.sd.phi,
                                       prior.var.beta=NULL, prior.tau2=NULL,
                                       beta.initial=NULL, phi.initial=NULL, tau2.initial=NULL){
  
  
  #### setup
  D <- length(data)
  K <- length(data[[1]]$y)
  p <- data[[1]]$p
  X.list <- lapply(data, function(z){ z$x.standardised})
  Y.list <- lapply(data, function(z){ z$y})
  X.all <- X.list[[1]]
  Y.all <- Y.list[[1]]
  if (D > 1){
    for (d in 2:D){ 
      X.all <- rbind(X.all, X.list[[d]])
      Y.all <- c(Y.all, Y.list[[d]])
    }
  }
  
  
  #### priors
  prior.mean.beta <- rep(0, p)
  if (is.null(prior.var.beta)){ prior.var.beta <- rep(1000, p) }
  if (is.null(prior.tau2)){ prior.tau2 <- c(1, 0.01) }
  
  
  #### Initial values
  mod.glm <- glm(Y.all~X.all-1, family="poisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  if (is.null(beta.initial)){
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  } else {
    beta <- beta.initial
  }
  
  
  beta.initial <- beta
  phi.list <- list()
  tau2.vals <- c()
  for (d in 1:D){
    Y.d <- Y.list[[d]]
    X.d <- X.list[[d]]
    log.Y <- log(Y.d)
    log.Y[Y.d==0] <- -0.1  
    res.temp <- log.Y - X.d %*% beta.mean
    res.sd <- sd(res.temp, na.rm=TRUE)
    if (is.null(phi.initial)){
      phi.d <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    } else {
      phi.d <- phi.initial[d]
    }
    if (is.null(tau2.initial)){
      tau2.d <- var(phi.d)/8
    } else {
      tau2.d <- tau2.initial[d]
    }
    phi.list[[d]] <- phi.d
    tau2.vals <- c(tau2.vals, tau2.d)
  }
  tau2 <- mean(tau2.vals)
  
  
  beta.initial <- beta
  phi.initials <- phi.list
  tau2.initial <- tau2
  
  
  #### MCMC quantities    
  n.sample <- 125000
  burnin <- 10000
  thin <- 1
  
  
  ## Storage
  n.keep <- floor((n.sample - burnin)/ thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.tau2 <- array(NA, c(n.keep, D))
  samples.phi <- array(NA, c(D, n.keep, K))
  
  ## Metropolis quantities
  accept <- rep(0, 2)
  tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
  
  ## CAR quantities
  ## block diagonalize the W matrix
  K.all <- K*D
  W.all <- W
  if (D > 1){
    for (i in 2:D){ W.all <- as.matrix(bdiag(W.all, W))}
  }
  W.quants <- common.Wcheckformat.leroux(W.all, K.all, fix.rho=FALSE, rho)
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta
    offsets.beta <- phi.list
    beta.out <- betaIntegratedUpdate(X.list, Y.list, beta, offsets.beta,
                                     prior.mean.beta, prior.var.beta, proposal.sd.beta)
    beta <- beta.out$beta
    
    
    ## Sample from phi
    ## phi random effects updated separately
    for (d in 1:D){
      y.d <- Y.list[[d]]
      phi.d <- phi.list[[d]]
      offsets.phi.d <- X.list[[d]] %*% beta
      phi.out.d <- phiPoissonUpdate(y.d, W, nNeighbors, phi.d, tau2, offsets.phi.d, proposal.sd.phi, rho=1)
      phi.list[[d]] <- phi.out.d$phi - mean(phi.out.d$phi)
    }
    
    
    ## Sample from tau2
    phi.combined <- unlist(phi.list)
    temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K.all, phi.combined, phi.combined, rho)
    tau2.posterior.scale <- temp2 + prior.tau2[2]
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))

    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta[ele,] <- beta
      accept[1] <- accept[1] + beta.out$accept
      samples.tau2[ele,] <- tau2
      for (d in 1:D){
        samples.phi[d,ele,] <- phi.list[[d]]
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep

  
  out <- list()
  out$samples.beta <- samples.beta
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$accept <- accept
  out$beta.initial <- beta.initial
  out$phi.initial <- phi.initials
  out$tau2.initial <- tau2.initial
  return(out)
  
  
}
