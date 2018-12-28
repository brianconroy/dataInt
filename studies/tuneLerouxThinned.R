library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


##########################################
# Tunes the unthinned integrated Leroux CAR
# model on multiple simulated datasets. As
# a start, the CAR models of datasets are 
# analyzed separately.

# integration, unthinned, separate CARs
#########################################


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


#### True parameter values
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75 
# rho.list <- c(rho.1, rho.2)
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1)
gamma.1 <- -0.5
gamma.2 <- -0.25
gamma.true <- c(gamma.1, gamma.2)


#### Sampling covariates
par(mfrow=c(1,3))
samp.disc.1 <- simRegion$raster[[c(3)]]
samp.disc.2 <- simRegion$raster[[c(4)]]
plot(caWc.disc[[1]])
plot(samp.disc.1)
plot(samp.disc.2)


#### Simulate data
data.1 <- simLerouxThinned(caWc.disc, samp.disc.1, beta.true, gamma.1, tau.sq.1, rho.1, nNeighbors, seed=401)
data.2 <- simLerouxThinned(caWc.disc, samp.disc.2, beta.true, gamma.2, tau.sq.2, rho.2, nNeighbors, seed=401)
data <- list(data.1, data.2)
par(mfrow=c(1, 2))
hist(data.1$y)
hist(data.2$y)


#### reference
glm(data.1$y ~ data.1$x.standardised + data.1$z.standardised-1, family='poisson')
glm(data.2$y ~ data.2$x.standardised + data.2$z.standardised-1, family='poisson')


#########################
## Tune Proposal & Priors
#########################


tune.beta <- c(0.01, 0.075)
tune.phi <- c(1)
tune.rho <- c(0.05)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p + length(data)*data.1$p.z))


par(mfrow=c(2,1))
tuneResults <- tuneMCMC(
  data,
  carLerouxIntegrationThinned,
  tune.beta,
  tune.phi,
  tune.rho,
  tune.tau.prior=tune.tau.prior,
  tune.beta.prior=tune.beta.prior,
  n.sample=65000,
  plot=TRUE,
  self.tune=TRUE,
  fix.rho=FALSE
)


#############
## Traceplots
#############

overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", c(rho.1, rho.2)), collapse=" ")
par(mfrow=c(2,1), oma=c(0,0,2,0))
viewTracesInt(tuneResults, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main))
viewTracesInt(tuneResults, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main))
viewTracesInt(tuneResults, 'tau2', line=tau.sq.1, d=1, overall.title=paste("tau2 (d=1) |",overall.main))
viewTracesInt(tuneResults, 'tau2', line=tau.sq.2, d=2, overall.title=paste("tau2 (d=2) |",overall.main))
viewTracesInt(tuneResults, 'rho.1', line=rho.1, d=1, overall.title=paste("rho (d=1) |",overall.main))
viewTracesInt(tuneResults, 'rho.2', line=rho.2, d=2, overall.title=paste("rho (d=2) |",overall.main))

viewHistsInt(results, 'beta0', line=beta0)
viewHistsInt(results, 'beta1', line=beta1)
viewHistsInt(results, 'beta2', line=gamma.1)
viewHistsInt(results, 'beta2', line=gamma.2)
viewHistsInt(results, 'tau2', line=tau.sq.1, d=1)
viewHistsInt(results, 'tau2', line=tau.sq.1, d=2)


for (i in 1:length(tune.beta)){
  for(j in 1:length(tune.phi)){
    for (k in 1:length(tune.tau.prior)){
      for (l in 1:length(tune.beta.prior)){
        
        
        print("***************")
        t.beta <- tune.beta[i]
        t.phi <- tune.phi[j]
        print(t.beta)
        print(t.phi)
      }
    }
  }
}

################
## Multiple Runs
################


t.beta <- 0.05
t.phi <- 1.75
prior.tau <- NULL
prior.beta <- NULL


## Random GLM Starting Values
n.runs <- 4
par(mfrow=c(2,2))
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLerouxIntegration(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau)
  results[[i]] <- output
  
  plot(output$samples.beta[,2], type='l', main='beta1')
  abline(h=beta1, col='2')
}


overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho.list), collapse=" ")
par(mfrow=c(2, 2), oma=c(0,0,2,0))
viewTracesInt(results, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main), title.type='title')
viewTracesInt(results, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main), title.type='title')
viewTracesInt(results, 'beta2', line=beta2, overall.title=paste("beta2 |",overall.main), title.type='title')
viewTracesInt(results, 'tau2', d=1, line=tau.sq.1, overall.title=paste("tau2 (d=1) |",overall.main), title.type='title')
viewTracesInt(results, 'tau2', d=2, line=tau.sq.1, overall.title=paste("tau2 (d=2) |",overall.main), title.type='title')


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


viewHists(results, 'beta0', line=beta0)
viewHists(results, 'beta1', line=beta1)
viewHists(results, 'beta2', line=beta2)
viewHists(results, 'tau2', line=tau.sq)


## Write results to csv
df <- c()
for (r in results){
  
  beta0.bias <- round(mean(r$samples.beta[,1]) - beta0, 3)
  beta1.bias <- round(mean(r$samples.beta[,2]) - beta1, 3)
  tau2.bias <- round(mean(r$samples.tau2) - tau.sq, 3)
  
  beta0.bias.prc <- round(100*(mean(r$samples.beta[,1]) - beta0)/beta0, 3)
  beta1.bias.prc <- round(100*(mean(r$samples.beta[,2]) - beta1)/beta1, 3)
  tau2.bias.prc <- round(100*(mean(r$samples.tau2) - tau.sq)/tau.sq, 3)
  
  row <- c(r$tune.beta, r$tune.phi,
           r$tune.tau.prior[1], r$tune.beta.prior[1],
           round(r$accept[1], 3), round(r$accept[2], 3),
           beta0.bias, beta1.bias, tau2.bias,
           beta0.bias.prc, beta1.bias.prc, tau2.bias.prc)
  df <- rbind(df, row)
  
}
df <- data.frame(df)
names(df) <- c('tune.beta', 'tune.phi', 'tune.tau.prior',
               'tune.beta.prior', 'accept.beta', 'accept.phi', 
               'beta0.bias', 'beta1.bias', 'tau2.bias',
               'beta0.bias.prc', 'beta1.bias.prc', 'tau2.bias.prc')
write.table(df, file='/Users/brianconroy/Documents/research/dataInt/tune_results.csv', 
            sep=',', row.names=F)
