library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


###############################################
# Tunes the thinned integrated Leroux CAR
# model on multiple simulated datasets, 
# including randomly sampled survey counts.

# To start, excludes systematic bias from thinning.

# integration, unthinned, uniform, separate CARs
###############################################


#######
# Setup
#######


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc)
W <- simRegion$W
caWc.disc <- simRegion$raster
caWc.disc <- caWc.disc[[c(1)]]
nNeighbors <- simRegion$nNeighbors


#### Sampling covariates
par(mfrow=c(1,3))
samp.disc.1 <- simRegion$raster[[c(3)]]
samp.disc.2 <- simRegion$raster[[c(4)]]
plot(caWc.disc[[1]])
plot(samp.disc.1)
plot(samp.disc.2)


#### True parameter values
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75 
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1)
gamma.1 <- -0.5
gamma.2 <- -0.25
gamma.true <- c(gamma.1, gamma.2)


#### Simulate data
data.1 <- simLerouxThinned(caWc.disc, samp.disc.1, beta.true, gamma.1, tau.sq.1, rho.1, nNeighbors, seed=401)
data.2 <- simLerouxThinned(caWc.disc, samp.disc.2, beta.true, gamma.2, tau.sq.2, rho.2, nNeighbors, seed=401)
data.unif <- simUniform(caWc.disc, beta.true, nsamp=25, seed=505)
data <- list(spatial=list(data.1, data.2), survey=data.unif)


par(mfrow=c(3,1))
hist(data.1$y)
hist(data.2$y)
hist(data.unif$y)


#### reference
glm(data.1$y ~ data.1$x.standardised + data.1$z.standardised -1, family='poisson')
glm(data.2$y ~ data.2$x.standardised + data.2$z.standardised -1, family='poisson')
glm(data.unif$y ~ data.unif$x.standardised-1, family='poisson')


#########################
## Tune Proposal & Priors
#########################


tune.beta <- c(0.01)
tune.phi <- c(0.5)
tune.rho <- c(0.05)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p + 2*data.1$p.z))


par(mfrow=c(1,1))
tuneResults <- tuneMCMC(
  data,
  carLerouxThinnedUnif,
  tune.beta,
  tune.phi,
  tune.rho,
  tune.tau.prior=tune.tau.prior,
  tune.beta.prior=tune.beta.prior,
  n.sample=150000,
  plot=TRUE,
  self.tune=TRUE,
  fix.rho=FALSE
)


par(mfrow=c(3,3))
params <- list(
  beta0=beta0,
  beta1=beta1,
  gamma.1=gamma.1,
  gamma.2=gamma.2,
  tau.sq.1=tau.sq.1,
  tau.sq.2=tau.sq.2
)
best <- bestTune(tuneResults, params)
overall.main <- paste(c("D:", length(data), "| beta:", 
                        beta.true, "| rho's:", rho.list), collapse=" ")
viewOutput(best$output, type='thinned', params, overall.main)
viewOutput(tuneResults[[1]], type='thinned', params, overall.main)


#############
## Traceplots
#############


overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", c(rho.1, rho.2)), collapse=" ")
par(mfrow=c(2, 2), oma=c(0,0,2,0))
viewTracesInt(tuneResults, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main))
viewTracesInt(tuneResults, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main))
viewTracesInt(tuneResults, 'beta2', line=gamma.1, overall.title=paste("gamma1 |",overall.main))
viewTracesInt(tuneResults, 'beta3', line=gamma.2, overall.title=paste("gamma2 |",overall.main))
viewTracesInt(tuneResults, 'tau2', line=tau.sq.1, d=1, overall.title=paste("tau2 (d=1) |",overall.main))
viewTracesInt(tuneResults, 'tau2', line=tau.sq.2, d=2, overall.title=paste("tau2 (d=2) |",overall.main))
viewTracesInt(tuneResults, 'rho.1', line=rho.1, d=1, overall.title=paste("rho (d=1) |",overall.main))
viewTracesInt(tuneResults, 'rho.2', line=rho.2, d=2, overall.title=paste("rho (d=2) |",overall.main))


viewHistsInt(tuneResults, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main))
viewHistsInt(tuneResults, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main))
viewHistsInt(tuneResults, 'beta2', line=gamma.1, overall.title=paste("beta2 |",overall.main))
viewHistsInt(tuneResults, 'tau2', line=tau.sq.1, d=1, overall.title=paste("tau2 (d=1) |",overall.main))
viewHistsInt(tuneResults, 'tau2', line=tau.sq.2, d=2, overall.title=paste("tau2 (d=2) |",overall.main))


################
## Multiple Runs
################


t.beta <- 0.01
t.phi <- 1.5
prior.tau <- c(1, 0.01)
prior.beta <- list(rep(1000, data.1$p))


## Random GLM Starting Values
n.runs <- 4
par(mfrow=c(2,2))
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLerouxIntegrationUnif(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau)
  results[[i]] <- output
  plot(output$samples.beta[,2], type='l', main='beta1')
  abline(h=beta1, col='2')
}


overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho.list), collapse=" ")
par(mfrow=c(2, 2), oma=c(0,0,2,0))
viewTracesInt(results, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main), title.type='title')
viewTracesInt(results, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main), title.type='title')
viewTracesInt(results, 'beta2', line=beta2, overall.title=paste("beta2 |",overall.main), title.type='title')
viewTracesInt(results, 'beta3', line=beta3, overall.title=paste("beta3 |",overall.main), title.type='title')
viewTracesInt(results, 'beta4', line=beta4, overall.title=paste("beta4 |",overall.main), title.type='title')
viewTracesInt(results, 'tau2', d=1, line=tau.sq.1, overall.title=paste("tau2 (d=1) |",overall.main), title.type='title')
viewTracesInt(results, 'tau2', d=2, line=tau.sq.2, overall.title=paste("tau2 (d=2) |",overall.main), title.type='title')


viewHistsInt(results, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main))
viewHistsInt(results, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main))
viewHistsInt(results, 'beta2', line=beta2, overall.title=paste("beta2 |",overall.main))
viewHistsInt(results, 'beta3', line=beta2, overall.title=paste("beta3 |",overall.main))
viewHistsInt(results, 'beta4', line=beta2, overall.title=paste("beta4 |",overall.main))
viewHistsInt(results, 'tau2', line=tau.sq.1, d=1, overall.title=paste("tau2 (d=1) |",overall.main))
viewHistsInt(results, 'tau2', line=tau.sq.2, d=2, overall.title=paste("tau2 (d=2) |",overall.main))


# pool samples
pool.results <- function(results){
  results.new <- list()
  new.beta <- results[[1]]$samples.beta
  new.tau2 <- results[[1]]$samples.tau2
  
  for (i in 2:length(results)){
    new.beta <- rbind(new.beta, results[[i]]$samples.beta)  
    new.tau2 <- rbind(new.tau2, results[[i]]$samples.tau2)
  }
  
  results.new$samples.beta <- new.beta
  results.new$samples.tau2 <- new.tau2
  
  return(results.new)
  
}


overall.main <- paste(c("D:", length(data), "| beta tune:", t.beta, "| phi tune:", t.phi), collapse=" ")
par(mfrow=c(3, 2), oma=c(0,0,2,0))
pooled <- pool.results(results)
hist(pooled$samples.beta[,1], main='beta0', xlab='value')
abline(v=beta0, col='2')
hist(pooled$samples.beta[,2], main='beta1', xlab='value')
abline(v=beta1, col='2')
hist(pooled$samples.beta[,3], main='beta2', xlab='value')
abline(v=beta2, col='2')
hist(pooled$samples.tau2[,1], main='tau2 (d=1)', xlab='value')
abline(v=tau.sq.1, col='2')
hist(pooled$samples.tau2[,2], main='tau2 (d=2)', xlab='value')
abline(v=tau.sq.2, col='2')
title(overall.main, outer=TRUE)
