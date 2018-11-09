
#############
## prefSample
#############
#   fits the preferential sampling
#   model by the MH-RW algorithm.
#   contains no gaussian process.
## arguments
#   data: output of simLocCond
#   n.sample: number of mcmc samples
#   burnin: burnin
prefSample <- function(data, n.sample=75000, burnin=10000, thin=1,
                       proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.01,
                       proposal.sd.alpha=0.05, self.tune=TRUE){
  
  
  #### Setup
  loc <- data$loc
  cond <- data$conditional
  K <- length(cond$y)
  p <- cond$p
  p.l <- loc$p
  X.c <- cond$x.standardised
  X.l <- loc$x.scaled
  if (ncol(X.l) > 1){
    X.l.sub <- X.l[as.logical(loc$status),]
  } else{
    X.l.sub <- matrix(X.l[as.logical(loc$status),])
  }
  Y.c <- cond$y
  Y.l <- loc$status
  
  
  #### Priors
  prior.mean.beta <- rep(0, p + p.l)
  prior.var.beta <- rep(1000, p + p.l)
  
  
  #### Initial values
  alpha <- runif(1, 0, 10)
  
  mod.c.glm <- glm(Y.c ~ X.c - 1, family="poisson")
  beta.c.mean <- mod.c.glm$coefficients
  beta.c.sd <- sqrt(diag(summary(mod.c.glm)$cov.scaled))
  beta.c <- rnorm(n=length(beta.c.mean), mean=beta.c.mean, sd=beta.c.sd)
  
  mod.l.glm <- glm(Y.l ~ X.l - 1, family='binomial')
  beta.l.mean <- mod.l.glm$coefficients
  beta.l.sd <- sqrt(diag(summary(mod.l.glm)$cov.scaled))
  beta.l <- rnorm(n=length(beta.l.mean), mean=beta.l.mean, sd=beta.l.sd)
  
  beta.c.initial <- beta.c
  beta.l.initial <- beta.l
  alpha.initial <- alpha
  
  
  #### Storage
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta.c <- array(NA, c(n.keep, p))
  samples.beta.l <- array(NA, c(n.keep, p.l))
  samples.alpha <- array(NA, c(n.keep, 1))
  samples.deviance <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 3)
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta (locational)
    beta.out.l <- betaLogisticUpdate(X.l, Y.l, beta.l, proposal.sd.beta.l)
    beta.l <- beta.out.l$beta
    
    
    ## Sample from beta (conditional)
    loc.pred <- expit(X.l.sub %*% beta.l)
    offset.beta.c <- alpha * loc.pred
    beta.out.c <- betaPoissonUpdate(X.c, Y.c, beta.c, offset.beta.c,
                                    prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.c <- beta.out.c$beta
    
    
    #  Sample from alpha
    alpha.out <- alphaPoissonUpdate(X.c, Y.c, alpha, loc.pred, beta.c, proposal.sd.alpha)
    alpha <- alpha.out$alpha
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta.c[ele,] <- beta.c
      samples.beta.l[ele,] <- beta.l
      samples.alpha[ele,] <- alpha
      accept[1] <- accept[1] + beta.out.l$accept
      accept[2] <- accept[2] + beta.out.c$accept
      accept[3] <- accept[3] + alpha.out$accept
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p.l>2){
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 40, 50)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 40, 50)
            proposal.sd.alpha <- tuners.sd.1(accept[3], ele, proposal.sd.alpha, 40, 50)
          } else {
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 30, 40)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 30, 40)
            proposal.sd.alpha <- tuners.sd.1(accept[3], ele, proposal.sd.alpha, 30, 40)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  
  
  out <- list()
  out$samples.beta <- samples.beta.c
  out$samples.beta.l <- samples.beta.l
  out$samples.alpha <- samples.alpha
  out$accept <- accept
  out$beta.c.initial <- beta.c.initial
  out$beta.l.initial <- beta.l.initial
  out$alpha.initial <- alpha.initial
  out$final.tune.beta.c <- proposal.sd.beta.c
  out$final.tune.beta.l <- proposal.sd.beta.l
  return(out)
  
  
}


# updates locational parameters from the 
# joint location-conditional likelihood,
# rather than the locational likelihood
# alone
prefSampleNew <- function(data, n.sample=75000, burnin=10000, thin=1,
                          proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.01,
                          proposal.sd.alpha=0.05, self.tune=TRUE){
  
  
  #### Setup
  loc <- data$loc
  cond <- data$conditional
  K <- length(cond$y)
  p <- cond$p
  p.l <- loc$p
  X.c <- cond$x.standardised
  X.l <- loc$x.scaled
  if (ncol(X.l) > 1){
    X.l.sub <- X.l[as.logical(loc$status),]
  } else{
    X.l.sub <- matrix(X.l[as.logical(loc$status),])
  }
  Y.c <- cond$y
  Y.l <- loc$status
  
  
  #### Priors
  prior.mean.beta <- rep(0, p + p.l)
  prior.var.beta <- rep(1000, p + p.l)
  
  
  #### Initial values
  alpha <- runif(1, 0, 10)
  
  mod.c.glm <- glm(Y.c ~ X.c - 1, family="poisson")
  beta.c.mean <- mod.c.glm$coefficients
  beta.c.sd <- sqrt(diag(summary(mod.c.glm)$cov.scaled))
  beta.c <- rnorm(n=length(beta.c.mean), mean=beta.c.mean, sd=beta.c.sd)
  
  mod.l.glm <- glm(Y.l ~ X.l - 1, family='binomial')
  beta.l.mean <- mod.l.glm$coefficients
  beta.l.sd <- sqrt(diag(summary(mod.l.glm)$cov.scaled))
  beta.l <- rnorm(n=length(beta.l.mean), mean=beta.l.mean, sd=beta.l.sd)
  
  beta.c.initial <- beta.c
  beta.l.initial <- beta.l
  alpha.initial <- alpha
  
  
  #### Storage
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta.c <- array(NA, c(n.keep, p))
  samples.beta.l <- array(NA, c(n.keep, p.l))
  samples.alpha <- array(NA, c(n.keep, 1))
  samples.deviance <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 3)
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta (locational)
    beta.out.l <- betaLocationalUpdate(X.l, Y.l, X.l.sub, X.c, Y.c, beta.l, beta.c, alpha, proposal.sd.beta.l)
    beta.l <- beta.out.l$beta
    
    
    ## Sample from beta (conditional)
    loc.pred <- expit(X.l.sub %*% beta.l)
    offset.beta.c <- alpha * loc.pred
    beta.out.c <- betaPoissonUpdate(X.c, Y.c, beta.c, offset.beta.c,
                                    prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.c <- beta.out.c$beta
    
    
    #  Sample from alpha
    alpha.out <- alphaPoissonUpdate(X.c, Y.c, alpha, loc.pred, beta.c, proposal.sd.alpha)
    alpha <- alpha.out$alpha
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta.c[ele,] <- beta.c
      samples.beta.l[ele,] <- beta.l
      samples.alpha[ele,] <- alpha
      accept[1] <- accept[1] + beta.out.l$accept
      accept[2] <- accept[2] + beta.out.c$accept
      accept[3] <- accept[3] + alpha.out$accept
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p.l>2){
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 40, 50)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 40, 50)
            proposal.sd.alpha <- tuners.sd.1(accept[3], ele, proposal.sd.alpha, 40, 50)
          } else {
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 30, 40)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 30, 40)
            proposal.sd.alpha <- tuners.sd.1(accept[3], ele, proposal.sd.alpha, 30, 40)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  
  
  out <- list()
  out$samples.beta <- samples.beta.c
  out$samples.beta.l <- samples.beta.l
  out$samples.alpha <- samples.alpha
  out$accept <- accept
  out$beta.c.initial <- beta.c.initial
  out$beta.l.initial <- beta.l.initial
  out$alpha.initial <- alpha.initial
  out$final.tune.beta.c <- proposal.sd.beta.c
  out$final.tune.beta.l <- proposal.sd.beta.l
  return(out)
  
  
}


# case control data
# cases and controls have the same covariates
# estimates alpha for cases and controls separately
prefSampleCC <- function(data, n.sample=75000, burnin=10000, thin=1,
                         proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.01,
                         proposal.sd.alpha=0.05, self.tune=TRUE){
  
  
  #### Setup
  loc <- data$loc
  case <- data$conditional.case
  control <- data$conditional.control
  K <- length(case$y)
  p <- case$p
  p.l <- loc$p
  X.c <- case$x.standardised
  X.l <- loc$x.scaled
  if (ncol(X.l) > 1){
    X.l.sub <- X.l[as.logical(loc$status),]
  } else{
    X.l.sub <- matrix(X.l[as.logical(loc$status),])
  }
  Y.case <- case$y
  Y.control <- control$y
  Y.c.list <- list(Y.case, Y.control)
  Y.l <- loc$status
  
  
  #### Priors
  prior.mean.beta <- rep(0, p + p.l)
  prior.var.beta <- rep(1000, p + p.l)
  
  
  #### Initial values
  alpha.case <- runif(1, 0, 10)
  alpha.control <- runif(1, -10, 0)
  
  mod.case.glm <- glm(Y.case ~ X.c - 1, family="poisson")
  beta.case.mean <- mod.case.glm$coefficients
  beta.case.sd <- sqrt(diag(summary(mod.case.glm)$cov.scaled))
  beta.case <- rnorm(n=length(beta.case.mean), mean=beta.case.mean, sd=beta.case.sd)
  
  mod.control.glm <- glm(Y.control ~ X.c - 1, family="poisson")
  beta.control.mean <- mod.control.glm$coefficients
  beta.control.sd <- sqrt(diag(summary(mod.control.glm)$cov.scaled))
  beta.control <- rnorm(n=length(beta.control.mean), mean=beta.control.mean, sd=beta.control.sd)
  
  mod.l.glm <- glm(Y.l ~ X.l - 1, family='binomial')
  beta.l.mean <- mod.l.glm$coefficients
  beta.l.sd <- sqrt(diag(summary(mod.l.glm)$cov.scaled))
  beta.l <- rnorm(n=length(beta.l.mean), mean=beta.l.mean, sd=beta.l.sd)
  
  beta.case.initial <- beta.case
  beta.control.initial <- beta.control
  beta.l.initial <- beta.l
  alpha.case.initial <- alpha.case
  alpha.control.initial <- alpha.control
  
  
  #### Storage
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta.case <- array(NA, c(n.keep, p))
  samples.beta.control <- array(NA, c(n.keep, p))
  samples.beta.l <- array(NA, c(n.keep, p.l))
  samples.alpha.case <- array(NA, c(n.keep, 1))
  samples.alpha.control <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 4)
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta (locational)
    beta.out.l <- betaLogisticUpdate(X.l, Y.l, beta.l, proposal.sd.beta.l)
    beta.l <- beta.out.l$beta
    
    
    ## Sample from beta case
    loc.pred <- expit(X.l.sub %*% beta.l)
    offset.beta.case <- alpha.case * loc.pred
    beta.out.case <- betaPoissonUpdate(X.c, Y.case, beta.case, offset.beta.case,
                                       prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.case <- beta.out.case$beta
    
    
    ## Sample from beta control
    offset.beta.control <- alpha.control * loc.pred
    beta.out.control <- betaPoissonUpdate(X.c, Y.control, beta.control, offset.beta.control,
                                          prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.control <- beta.out.control$beta
    
    
    #  Sample from alpha case
    alpha.case.out <- alphaPoissonUpdate(X.c, Y.case, alpha.case, loc.pred, beta.case, proposal.sd.alpha)
    alpha.case <- alpha.case.out$alpha
    
    
    #  Sample from alpha control
    alpha.control.out <- alphaPoissonUpdate(X.c, Y.control, alpha.control, loc.pred, beta.control, proposal.sd.alpha)
    alpha.control <- alpha.control.out$alpha
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta.case[ele,] <- beta.case
      samples.beta.control[ele,] <- beta.control
      samples.beta.l[ele,] <- beta.l
      samples.alpha.case[ele,] <- alpha.case
      samples.alpha.control[ele,] <- alpha.control
      accept[1] <- accept[1] + beta.out.l$accept
      accept[2] <- accept[2] + beta.out.case$accept
      accept[3] <- accept[3] + alpha.case.out$accept
      accept[4] <- accept[4] + alpha.control.out$accept
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p.l>2){
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 40, 50)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 40, 50)
            proposal.sd.alpha <- tuners.sd.1((accept[3] + accept[4])/2, ele, proposal.sd.alpha, 40, 50)
          } else {
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 30, 40)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 30, 40)
            proposal.sd.alpha <- tuners.sd.1((accept[3] + accept[4])/2, ele, proposal.sd.alpha, 30, 40)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  
  
  out <- list()
  out$samples.beta.case <- samples.beta.case
  out$samples.beta.control <- samples.beta.control
  out$samples.beta.l <- samples.beta.l
  out$samples.alpha.case <- samples.alpha.case
  out$samples.alpha.control <- samples.alpha.control
  out$accept <- accept
  out$beta.case.initial <- beta.case.initial
  out$beta.control.initial <- beta.control.initial
  out$beta.l.initial <- beta.l.initial
  out$alpha.case.initial <- alpha.case.initial
  out$alpha.control.initial <- alpha.control.initial
  out$final.tune.beta.c <- proposal.sd.beta.c
  out$final.tune.beta.l <- proposal.sd.beta.l
  return(out)
  
  
}


# case control data
# cases and controls have the same covariates
prefSampleConstrainedCC <- function(data, n.sample=75000, burnin=10000, thin=1,
                                  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.01,
                                  proposal.sd.alpha=0.05, self.tune=TRUE){
  
  
  #### Setup
  loc <- data$loc
  case <- data$conditional.case
  control <- data$conditional.control
  K <- length(case$y)
  p <- case$p
  p.l <- loc$p
  X.c <- case$x.standardised
  X.l <- loc$x.scaled
  X.l.sub <- matrix(X.l[as.logical(loc$status),])
  Y.case <- case$y
  Y.control <- control$y
  Y.c.list <- list(Y.case, Y.control)
  Y.l <- loc$status
  
  
  #### Priors
  prior.mean.beta <- rep(0, p + p.l)
  prior.var.beta <- rep(1000, p + p.l)
  
  
  #### Initial values
  alpha <- runif(1, 0, 10)
  
  mod.case.glm <- glm(Y.case ~ X.c - 1, family="poisson")
  beta.case.mean <- mod.case.glm$coefficients
  beta.case.sd <- sqrt(diag(summary(mod.case.glm)$cov.scaled))
  beta.case <- rnorm(n=length(beta.case.mean), mean=beta.case.mean, sd=beta.case.sd)
  
  mod.control.glm <- glm(Y.control ~ X.c - 1, family="poisson")
  beta.control.mean <- mod.control.glm$coefficients
  beta.control.sd <- sqrt(diag(summary(mod.control.glm)$cov.scaled))
  beta.control <- rnorm(n=length(beta.control.mean), mean=beta.control.mean, sd=beta.control.sd)
  
  mod.l.glm <- glm(Y.l ~ X.l - 1, family='binomial')
  beta.l.mean <- mod.l.glm$coefficients
  beta.l.sd <- sqrt(diag(summary(mod.l.glm)$cov.scaled))
  beta.l <- rnorm(n=length(beta.l.mean), mean=beta.l.mean, sd=beta.l.sd)
  
  beta.case.initial <- beta.case
  beta.control.initial <- beta.control
  beta.l.initial <- beta.l
  alpha.initial <- alpha
  
  
  #### Storage
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta.case <- array(NA, c(n.keep, p))
  samples.beta.control <- array(NA, c(n.keep, p))
  samples.beta.l <- array(NA, c(n.keep, p.l))
  samples.alpha <- array(NA, c(n.keep, 1))
  
  ## Metropolis quantities
  accept <- rep(0, 3)
  
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  
  
  #### Create the MCMC samples
  for(i in 1:n.sample){
    
    
    ## Sample from beta (locational)
    beta.out.l <- betaLogisticUpdate(X.l, Y.l, beta.l, proposal.sd.beta.l)
    beta.l <- beta.out.l$beta
    
    
    ## Sample from beta case
    loc.pred <- expit(X.l.sub %*% beta.l)
    offset.beta.case <- alpha * loc.pred
    beta.out.case <- betaPoissonUpdate(X.c, Y.case, beta.case, offset.beta.case,
                                       prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.case <- beta.out.case$beta
    
    
    ## Sample from beta control
    offset.beta.control <- alpha * loc.pred
    beta.out.control <- betaPoissonUpdate(X.c, Y.control, beta.control, offset.beta.control,
                                          prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.control <- beta.out.control$beta
    
    
    #  Sample from alpha
    alpha.out <- alphaPoissonUpdateCC(X.c, Y.case, Y.control, alpha, loc.pred, 
                                      beta.case, beta.control, proposal.sd.alpha)
    alpha <- alpha.out$alpha
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta.case[ele,] <- beta.case
      samples.beta.control[ele,] <- beta.control
      samples.beta.l[ele,] <- beta.l
      samples.alpha[ele,] <- alpha
      accept[1] <- accept[1] + beta.out.l$accept
      accept[2] <- accept[2] + beta.out.case$accept
      accept[3] <- accept[3] + alpha.out$accept
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p.l>2){
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 40, 50)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 40, 50)
            proposal.sd.alpha <- tuners.sd.1(accept[3], ele, proposal.sd.alpha, 40, 50)
          } else {
            proposal.sd.beta.l <- tuners.sd.1(accept[1], ele, proposal.sd.beta.l, 30, 40)
            proposal.sd.beta.c <- tuners.sd.1(accept[2], ele, proposal.sd.beta.c, 30, 40)
            proposal.sd.alpha <- tuners.sd.1(accept[3], ele, proposal.sd.alpha, 30, 40)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  
  
  out <- list()
  out$samples.beta.case <- samples.beta.case
  out$samples.beta.control <- samples.beta.control
  out$samples.beta.l <- samples.beta.l
  out$samples.alpha <- samples.alpha
  out$accept <- accept
  out$beta.case.initial <- beta.case.initial
  out$beta.control.initial <- beta.control.initial
  out$beta.l.initial <- beta.l.initial
  out$alpha.initial <- alpha.initial
  out$final.tune.beta.c <- proposal.sd.beta.c
  out$final.tune.beta.l <- proposal.sd.beta.l
  return(out)
  
  
}


# integrates case/control preferential sampling data with
# case only opportunistic counts.
prefSampleIntegration <- function(data, n.sample=75000, burnin=10000, thin=1,
                                  proposal.sd.beta.c=0.01, proposal.sd.beta.ctr=0.01,
                                  proposal.sd.beta.l=0.01,
                                  proposal.sd.alpha=0.05, proposal.sd.phi=0.03,
                                  proposal.sd.rho=0.05, self.tune=TRUE, fix.rho=FALSE, rho=NULL){
  
  
  #### setup
  ps <- data$ps
  spat <- data$spat
  loc <- ps$loc
  case <- ps$conditional.case
  control <- ps$conditional.control
  
  p <- case$p
  p.l <- loc$p
  
  Y.l <- loc$status
  Y.ps.case <- case$y
  Y.ps.control <- control$y
  Y.ps.list <- list(Y.ps.case, Y.ps.control)
  
  X.ps <- case$x.standardised
  X.l <- loc$x.scaled
  X.l.sub <- matrix(X.l[as.logical(loc$status),])
  
  D <- length(spat)
  K <- length(spat[[1]]$y)
  X.list <- lapply(spat, function(z){ z$x.standardised})
  Y.list <- lapply(spat, function(z){ z$y})
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
  prior.var.beta <- rep(1000, p)
  prior.tau2 <- c(1, 0.01)
  
  
  #### Initial values
  alpha.case <- runif(1, 0, 10)
  alpha.control <- runif(1, -10, 0)
  
  mod.case.glm <- glm(Y.ps.case ~ X.ps - 1, family="poisson")
  beta.case.mean <- mod.case.glm$coefficients
  beta.case.sd <- sqrt(diag(summary(mod.case.glm)$cov.scaled))
  beta.case <- rnorm(n=length(beta.case.mean), mean=beta.case.mean, sd=beta.case.sd)
  
  mod.control.glm <- glm(Y.ps.control ~ X.ps - 1, family="poisson")
  beta.control.mean <- mod.control.glm$coefficients
  beta.control.sd <- sqrt(diag(summary(mod.control.glm)$cov.scaled))
  beta.control <- rnorm(n=length(beta.control.mean), mean=beta.control.mean, sd=beta.control.sd)
  
  mod.l.glm <- glm(Y.l ~ X.l - 1, family='binomial')
  beta.l.mean <- mod.l.glm$coefficients
  beta.l.sd <- sqrt(diag(summary(mod.l.glm)$cov.scaled))
  beta.l <- rnorm(n=length(beta.l.mean), mean=beta.l.mean, sd=beta.l.sd)
  
  
  phi.list <- list()
  tau2.list <- list()
  for (d in 1:D){
    Y.d <- Y.list[[d]]
    X.d <- X.list[[d]]
    log.Y <- log(Y.d)
    log.Y[Y.d==0] <- -0.1  
    res.temp <- log.Y - X.d %*% beta.case
    res.sd <- sd(res.temp, na.rm=TRUE)
    phi.d <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    tau2.d <- var(phi.d)/4
    phi.list[[d]] <- phi.d
    tau2.list[[d]] <- tau2.d
  }
  if (!fix.rho){
    rho.list <- lapply(1:D, function(x) {runif(1)})
  }
  
  
  beta.case.initial <- beta.case
  beta.control.initial <- beta.control
  beta.l.initial <- beta.l
  alpha.case.initial <- alpha.case
  alpha.control.initial <- alpha.control
  phi.initials <- phi.list
  tau2.initials <- tau2.list
  
  #### Storage
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta.case <- array(NA, c(n.keep, p))
  samples.beta.control <- array(NA, c(n.keep, p))
  samples.beta.l <- array(NA, c(n.keep, p.l))
  samples.alpha.case <- array(NA, c(n.keep, 1))
  samples.alpha.control <- array(NA, c(n.keep, 1))
  samples.tau2 <- array(NA, c(n.keep, D))
  samples.phi <- array(NA, c(D, n.keep, K))
  if (!fix.rho){ samples.rho <- array(NA, c(D, n.keep, 1))}
  
  ## Metropolis quantities
  accept <- rep(0, 5)
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
    
    
    ## Sample from beta (locational)
    beta.out.l <- betaLogisticUpdate(X.l, Y.l, beta.l, proposal.sd.beta.l)
    beta.l <- beta.out.l$beta

    
    ## Sample from beta case
    loc.pred <- expit(X.l.sub %*% beta.l)
    offset.beta.case <- alpha.case * loc.pred
    offsets.beta.spat <- phi.list
    beta.out.case <- betaPsIntegratedUpdate(X.list, Y.list, X.ps, Y.ps.case, beta.case, offsets.beta.spat, 
                                            offset.beta.case, prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
    beta.case <- beta.out.case$beta

    
    ## Sample from beta control
    offset.beta.control <- alpha.control * loc.pred
    beta.out.control <- betaPoissonUpdate(X.ps, Y.ps.control, beta.control, offset.beta.control,
                                          prior.mean.beta, prior.var.beta, proposal.sd.beta.ctr)
    beta.control <- beta.out.control$beta

    
    #  Sample from alpha case
    alpha.case.out <- alphaPoissonUpdate(X.ps, Y.ps.case, alpha.case, loc.pred, beta.case, proposal.sd.alpha)
    alpha.case <- alpha.case.out$alpha
    
    
    #  Sample from alpha control
    alpha.control.out <- alphaPoissonUpdate(X.ps, Y.ps.control, alpha.control, loc.pred, beta.control, proposal.sd.alpha)
    alpha.control <- alpha.control.out$alpha
    
    
    ## Sample from phi
    phi.accept <- 0
    for (d in 1:D){
      y.d <- Y.list[[d]]
      phi.d <- phi.list[[d]]
      rho.d <- rho.list[[d]]
      tau2.d <- tau2.list[[d]]
      offsets.phi.d <- X.list[[d]] %*% beta.case
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
    
    
    ## Save the results
    if(i > burnin & (i-burnin)%%thin==0){
      ele <- (i - burnin)/thin
      samples.beta.l[ele,] <- beta.l
      samples.beta.case[ele,] <- beta.case
      samples.beta.control[ele,] <- beta.control
      samples.alpha.case[ele,] <- alpha.case
      samples.alpha.control[ele,] <- alpha.control
      accept[1] <- accept[1] + beta.out.case$accept
      accept[2] <- accept[2] + beta.out.control$accept
      accept[3] <- accept[3] + alpha.case.out$accept
      accept[4] <- accept[4] + alpha.control.out$accept
      accept[5] <- accept[5] + phi.accept
      for (d in 1:D){
        samples.phi[d,ele,] <- phi.list[[d]]
        samples.tau2[ele,][d] <- tau2.list[[d]]
        if (!fix.rho){
          samples.rho[d,ele,] <- rho.list[[d]]
          accept[6] <- accept[6] + rho.accept
        }
      }
      
      
      ## tune the proposal standard deviations
      if (self.tune){
        k <- i/100
        if(ceiling(k)==floor(k)){
          if (p>2){
            proposal.sd.beta.c <- tuners.sd.1(accept[1], ele, proposal.sd.beta.c, 40, 50)
            proposal.sd.beta.ctr <- tuners.sd.1(accept[2], ele, proposal.sd.beta.ctr, 40, 50)
          } else {
            proposal.sd.beta.c <- tuners.sd.1(accept[1], ele, proposal.sd.beta.c, 30, 40)
            proposal.sd.beta.ctr <- tuners.sd.1(accept[2], ele, proposal.sd.beta.ctr, 30, 40)
          }
          proposal.sd.alpha <- tuners.sd.1((accept[3] + accept[4])/2, ele, proposal.sd.alpha, 40, 50)
          proposal.sd.phi <- tuners.sd.1(accept[5], ele*D*K, proposal.sd.phi, 55, 65)
          if(!fix.rho){
            proposal.sd.rho <- tuners.sd.2(accept[6], ele*D, proposal.sd.rho, 40, 50, 2)
          }
        }
      }
    }
    
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
    
  }
  
  
  accept <- accept/n.keep
  accept[5] <- accept[5]/(D*K)
  if (!fix.rho){ accept[6] <- accept[6]/D}
  
  
  out <- list()
  out$samples.beta.l <- samples.beta.l
  out$samples.beta.case <- samples.beta.case
  out$samples.beta.control <- samples.beta.control
  out$samples.alpha.case <- samples.alpha.case
  out$samples.alpha.control <- samples.alpha.control
  out$samples.phi <- samples.phi
  out$samples.tau2 <- samples.tau2
  out$accept <- accept
  out$beta.initial <- beta.case.initial
  out$phi.initial <- phi.initials
  out$tau2.initial <- tau2.initials
  out$final.tune.beta.case <- proposal.sd.beta.c
  out$final.tune.beta.ctrl <- proposal.sd.beta.ctr
  out$final.tune.phi <- proposal.sd.phi
  if (!fix.rho){
    out$final.tune.rho <- proposal.sd.rho
    out$samples.rho <- samples.rho
  }
  return(out)
  
  
}
