library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


##########################
# compare parameter bias & 
# model BIC of thinned vs
# unthinned models.

# iterate over different 
# strengths of thinning

# create
#   a table of thinning strengths & 
#   resulting param biases.
##########################


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
rho.list <- c(rho.1, rho.2)
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1)


#### Sampling covariates
par(mfrow=c(1,3))
samp.disc.1 <- simRegion$raster[[c(3)]]
samp.disc.2 <- simRegion$raster[[c(4)]]
plot(caWc.disc[[1]], main='environmental surface')
plot(samp.disc.1, main='sampling surface 1')
plot(samp.disc.2, main='sampling surface 2')


# how do you deal with each level of 
# thinning needing its own tuning?
# tune each setup separately & save results


gamma.list <- list(-2, -1, -0.25)
results.thinned.all <- list()
results.ignore.all <- list()
tune.params.thinned <- list()
tune.params.ignore <- list()
counter <- 1


i <- 1
gamma.1 <- gamma.list[[i]]
gamma.2 <- gamma.list[[i]]
cat("****************\n")
cat("gamma 1: ", gamma.1, "\n")
cat("gamma 2: ", gamma.2, "\n")


#### Simulate data
data.1 <- simLerouxThinned(
  caWc.disc, samp.disc.1, beta.true, gamma.1, 
  tau.sq.1, rho.1, nNeighbors, seed=401)
data.2 <- simLerouxThinned(
  caWc.disc, samp.disc.2, beta.true, gamma.2, 
  tau.sq.2, rho.2, nNeighbors, seed=401)
data <- list(data.1, data.2)


cat("running the thinned integrated analyses\n")
params <- c('beta0', 'beta1', 'tau2.1', 'tau2.2')
truevals <- c(beta0, beta1, tau.sq.1, tau.sq.2)
tune.beta <- c(0.025, 0.05)
tune.phi <- c(0.5, 1.5)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p + length(data)*data.1$p.z))
tuneResults <- tuneMCMC(
  data,
  carLerouxIntegrationThinned,
  tune.beta,
  tune.phi,
  tune.tau.prior=tune.tau.prior,
  tune.beta.prior=tune.beta.prior,
  n.sample=85000)
best <- bestTune(tuneResults, params, truevals)
best$output[['gamma.1']] <- gamma.1
best$output[['gamma.2']] <- gamma.2
viewOutput(best)
print(best$best.configs)
tune.params.thinned[[counter]] <- list(
  gamma.1=gamma.1,
  gamma.2=gamma.2,
  model='thinned integrated',
  tune.beta=best$best.configs$tune.beta,
  tune.phi=best$best.configs$tune.phi
)
results.thinned.all[[counter]] <- best$output


cat("running the integrated analyses ignoring thinning\n")
tuneResults <- tuneMCMC(
  data,
  carLerouxIntegration,
  tune.beta=c(0.025),
  tune.phi=c(0.5),
  tune.tau.prior=tune.tau.prior,
  tune.beta.prior=tune.beta.prior,
  n.sample=85000)
best <- bestTune(tuneResults, params, truevals)
best$output[['gamma.1']] <- gamma.1
best$output[['gamma.2']] <- gamma.2
print(best$best.configs)
viewBest(best, type='unthinned')
tune.params.ignore[[counter]] <- list(
  gamma.1=gamma.1,
  gamma.2=gamma.2,
  model='thinned integrated',
  tune.beta=best$best.configs$tune.beta,
  tune.phi=best$best.configs$tune.phi
)
results.ignore.all[[counter]] <- best$output


counter <- counter + 1


# save tuning params
list2csv(tune.params.thinned, 'tune_params_thinned.csv')
list2csv(tune.params.ignore, 'tune_params_ignore.csv')


# report some traceplots
j <- 2
truth <- list(
  beta0=beta0,
  beta1=beta1,
  gamma.1=gamma.list[[j]],
  gamma.2=gamma.list[[j]],
  tau.sq.1=tau.sq.1,
  tau.sq.2=tau.sq.2
)
viewOutput(results.thinned.all[[j]], truevals=truth)


## summarize each model
summaries <- list()
counter <- 1
for (r in results.thinned.all){
  r.summary <- summarize(r, params, truevals)
  r.summary[['gamma.1']] <- gamma.list[[counter]]
  r.summary[['gamma.2']] <- gamma.list[[counter]]
  r.summary[['model']] <- 'integrated thinned'
  summaries[[counter]] <- r.summary
  counter <- counter + 1
}

for (r in results.ignore.all){
  r.summary <- summarize(r, params, truevals)
  r.summary[['gamma.1']] <- gamma.list[[counter-3]]
  r.summary[['gamma.2']] <- gamma.list[[counter-3]]
  r.summary[['model']] <- 'integrated'
  summaries[[counter]] <- r.summary
  counter <- counter + 1
}

list2csv(summaries, "thin_comparison.csv")


## summarize absolute differences in bias
diffs <- list()
counter <- 0
for (i in 1:(length(summaries)/2)){
  summary.thin <- summaries[[i]]
  summary.unthin <- summaries[[i+3]]
  diffs.i <- list()
  
  diffs.i[['gamma.1']] <- summary.thin$gamma.1
  diffs.i[['gamma.2']] <- summary.thin$gamma.2
  diffs.i[['beta0']] <- abs(summary.unthin$bias.beta0) - abs(summary.thin$bias.beta0)
  diffs.i[['beta1']] <- abs(summary.unthin$bias.beta1) - abs(summary.thin$bias.beta1)
  diffs.i[['tau2.1']] <- abs(summary.unthin$bias.tau2.1) - abs(summary.thin$bias.tau2.1)
  diffs.i[['tau2.2']] <- abs(summary.unthin$bias.tau2.2) - abs(summary.thin$bias.tau2.2)
  
  diffs.i[['perc.beta0']] <- abs(summary.unthin$percbias.beta0) - abs(summary.thin$percbias.beta0)
  diffs.i[['perc.beta1']] <- abs(summary.unthin$percbias.beta1) - abs(summary.thin$percbias.beta1)
  diffs.i[['perc.tau2.1']] <- abs(summary.unthin$percbias.tau2.1) - abs(summary.thin$percbias.tau2.1)
  diffs.i[['perc.tau2.2']] <- abs(summary.unthin$percbias.tau2.2) - abs(summary.thin$percbias.tau2.2)
  diffs[[i]] <- diffs.i
  counter <- counter + 1
}

list2csv(diffs, "thin_comparison_diffs.csv")
