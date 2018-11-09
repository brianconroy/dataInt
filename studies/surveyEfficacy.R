library(plyr)
library(R.utils)
library(ggplot2)
sourceDirectory('Documents/research/dataInt/R/')


##########################################
# study 1:biased areal counts + unbiased surveys.
#   vary: type of bias (intercept only vs intercept 
#   + covariate driven)
#     vary: degree of bias (low, medium, high)
#       vary: nsamp
#########################################


#########################
# case 0: 
#   verify your alg converges
#   if given bias covariates
#########################


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc)
W <- simRegion$W
caWc.disc <- simRegion$raster
caWc.disc <- caWc.disc[[c(1, 4)]]
nNeighbors <- simRegion$nNeighbors


#### True parameter values
tau.sq.1 <- 0.35
rho.1 <- 0.5
beta0 <- 1
beta1 <- 2
beta2 <- -1
beta.true <- c(beta0, beta1, beta2)


#### Simulate data
data.1 <- simLeroux(caWc.disc, beta.true, tau.sq.1, rho.1, nNeighbors, seed=405)
data.unif <- simUniform(caWc.disc, beta.true, nsamp=25, seed=505)
data <- list(spatial=list(data.1), survey=data.unif)


par(mfrow=c(1,2))
hist(data.1$y)
hist(data.unif$y)
summary(data.1$y)
summary(data.unif$y)


glm(data.1$y ~ data.1$x.standardised-1, family='poisson')
glm(data.unif$y ~ data.unif$x.standardised-1, family='poisson')


#### Tune mcmc
tune.beta <- c(0.01, 0.05)
tune.phi <- c(1)
tune.rho <- c(0.05)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p))


par(mfrow=c(2,1))
tuneResults <- tuneMCMC(
  data,
  carLerouxIntegrationUnif,
  tune.beta,
  tune.phi,
  tune.rho,
  tune.tau.prior=tune.tau.prior,
  tune.beta.prior=tune.beta.prior,
  n.sample=100000,
  plot=TRUE,
  self.tune=TRUE,
  fix.rho=FALSE
)


overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho.1), collapse=" ")
params <- c('beta0', 'beta1', 'beta2', 'tau2.1', 'rho.1')
par(mfrow=c(2, 1), oma=c(0,0,2,0))
viewResults(tuneResults, params, overall.main)


#### Multiple runs
t.beta <- 0.05
t.phi <- 1
t.rho <- 0.05
prior.tau <- c(1, 0.01)
prior.beta <- rep(1000, data.1$p)

n.runs <- 4
par(mfrow=c(2,2))
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLerouxIntegrationUnif(
    data, t.beta, t.phi, t.rho, prior.var.beta=prior.beta, 
    prior.tau2=prior.tau, fix.rho=FALSE, n.sample=100000, self.tune=TRUE)
  results[[i]] <- output
  plot(output$samples.beta[,1], type='l', main='beta0')
  abline(h=beta0, col='2')
}

overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho.1), collapse=" ")
params <- c('beta0', 'beta1', 'beta2', 'tau2.1', 'rho.1')
par(mfrow=c(2, 2), oma=c(0,0,2,0))
viewResults(results, params, overall.main)


#### Pool samples
pooled <- pool.results(results)
overall.main <- "Pooled Results"
params <- c('beta0', 'beta1', 'beta2', 'tau2.1', 'rho.1')
par(mfrow=c(5, 2), oma=c(0,0,2,0))
viewResults(list(pooled), params, overall.main)

param.list <- list(
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  tau.sq.1=tau.sq.1,
  rho.1=rho.1
)
pooled.summary <- summarize(pooled, param.list)
list2csv(pooled.summary, 'sE_case0_summary.csv')


#########################
# case 1: 
#   covariate driven bias
#   single biased dataset
#########################


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
caWc.disc <- caWc.disc[[c(1, 4)]]
caWc.disc.unb <- caWc.disc[[1]]
nNeighbors <- simRegion$nNeighbors


#### True parameter values
tau.sq.1 <- 0.35
rho.1 <- 0.5
beta0 <- 1
beta1 <- 2
beta2 <- -1
beta.true <- c(beta0, beta1, beta2)
beta.unb <- c(beta0, beta1)


#### Simulate thinned data
data.1 <- simLeroux(caWc.disc, beta.true, tau.sq.1, rho.1, nNeighbors, seed=405)
data.1$x.standardised <- data.1$x.standardised[,1:2]
data.1$x <- data.1$x[,1:2]
data.1$p <- 2
data.unif <- simUniform(caWc.disc.unb, beta.unb, nsamp=25, seed=505)
data <- list(spatial=list(data.1), survey=data.unif)


par(mfrow=c(1,2))
hist(data.1$y)
hist(data.unif$y)
summary(data.1$y)
summary(data.unif$y)


glm(data.1$y ~ data.1$x.standardised-1, family='poisson')
glm(data.unif$y ~ data.unif$x.standardised-1, family='poisson')


#### Tune mcmc
tune.beta <- c(0.01, 0.05)
tune.phi <- c(1)
tune.rho <- c(0.05)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p))


par(mfrow=c(2,1))
tuneResults <- tuneMCMC(
  data,
  carLerouxIntegrationUnif,
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


overall.main <- paste(c("beta:", beta.true, "| rho's:", rho.1), collapse=" ")
params <- c('beta0', 'beta1', 'tau2.1', 'rho.1')
par(mfrow=c(2, 1), oma=c(0,0,2,0))
viewResults(tuneResults, params, overall.main)


#### Multiple runs
t.beta <- 0.05
t.phi <- 1
t.rho <- 0.05
prior.tau <- c(1, 0.01)
prior.beta <- rep(1000, data.1$p)

n.runs <- 4
par(mfrow=c(2,2))
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLerouxIntegrationUnif(
    data, t.beta, t.phi, t.rho, prior.var.beta=prior.beta, 
    prior.tau2=prior.tau, fix.rho=FALSE, n.sample=100000, self.tune=TRUE)
  results[[i]] <- output
  plot(output$samples.beta[,1], type='l', main='beta0')
  abline(h=beta0, col='2')
}

overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho.1), collapse=" ")
params <- c('beta0', 'beta1', 'tau2.1', 'rho.1')
par(mfrow=c(2, 2), oma=c(0,0,2,0))
viewResults(results, params, overall.main)


#### Pool samples
pooled <- pool.results(results)
overall.main <- "Pooled Results"
params <- c('beta0', 'beta1', 'tau2.1', 'rho.1')
par(mfrow=c(4, 2), oma=c(0,0,2,0))
viewResults(list(pooled), params, overall.main)

param.list <- list(
  beta0=beta0,
  beta1=beta1,
  tau.sq.1=tau.sq.1,
  rho.1=rho.1
)
pooled.summary <- summarize(pooled, param.list)
list2csv(pooled.summary, 'sE_case1_summary.csv')


#########################
# case 2: 
#   systematic bias
#########################


# 2.1
#   start with a single areal dataset
#   iterate over the thinning value 
#   plot effects of percent bias


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
tau.sq.1 <- 0.35
rho.1 <- 0.5
beta0 <- 1
beta1 <- 2
gamma <- -0.5
beta.true <- c(beta0, beta1)


#### Simulate data
data.1 <- simLerouxBias(caWc.disc, beta.true, gamma, tau.sq.1, rho.1, nNeighbors, seed=405)
data <- list(data.1)


glm(data.1$y ~ data.1$x.standardised-1, family='poisson')


#### Tune mcmc
tune.beta <- c(0.01, 0.05)
tune.phi <- c(1)
tune.rho <- c(0.05)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p))


par(mfrow=c(2,1))
tuneResults <- tuneMCMC(
  data,
  carLerouxIntegration,
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


overall.main <- paste(c("beta:", beta.true, "| rho's:", rho.1), collapse=" ")
params <- c('beta0', 'beta1', 'tau2.1', 'rho.1')
par(mfrow=c(2, 1), oma=c(0,0,2,0))
viewResults(tuneResults, params, overall.main)


pooled <- pool.results(tuneResults)
overall.main <- paste("Pooled Results (gamma = ", gamma, ")", sep="")
params <- c('beta0', 'beta1', 'tau2.1', 'rho.1')
par(mfrow=c(4, 2), oma=c(0,0,2,0))
viewResults(list(pooled), params, overall.main)

param.list <- list(
  beta0=beta0,
  beta1=beta1,
  tau.sq.1=tau.sq.1,
  rho.1=rho.1
)
pooled.summary <- summarize(pooled, param.list)
list2csv(pooled.summary, 'sE_case2_summary1.csv')


# 2.2
#   introduce uniform (unbiased) surveys
#   low / med / thinning
#   iterate over number of surveys


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
tau.sq.1 <- 0.35
rho.1 <- 0.5
beta0 <- 1
beta1 <- 2
beta.true <- c(beta0, beta1)


#### Monte carlo glm estimates
nsim <- 500
samps <- c(5, 10, 25, 50, 65)
est.glm <- lapply(samps, function(x){array(NA, c(nsim, 2))})
for (i in 1:nsim){
  
  counter.samp <- 1
  for (nsamp in samps){
    data.unif <- simUniform(caWc.disc, beta.true, nsamp=nsamp)
    mod <- glm(data.unif$y ~ data.unif$x.standardised - 1, family='poisson')
    est.glm[[counter.samp]][i,] <- mod$coefficients
    counter.samp <- counter.samp + 1
  }
  
}


results.glm <- list()
counter.glm <- 1
for (i in 1:length(est.glm)){
  
  est <- est.glm[[i]]
  mu.beta0 <- mean(est[,1])
  sd.beta0 <- sd(est[,1])
  mu.beta1 <- mean(est[,2])
  sd.beta1 <- sd(est[,2])
  
  summary.n <- list()
  summary.n$glm.beta0 <- round(mu.beta0, 3)
  summary.n$glm.sd.beta0 <- round(sd.beta0, 3)
  summary.n$glm.beta1 <- round(mu.beta1, 3)
  summary.n$glm.sd.beta1 <- round(sd.beta1, 3)
  summary.n$glm.pbias.beta0 <- round(100*(mu.beta0-beta0)/beta0, 3)
  summary.n$glm.pbias.beta1 <- round(100*(mu.beta1-beta1)/beta1, 3)
  summary.n$nsamp <- samps[i]
  results.glm[[counter.glm]] <- summary.n
  counter.glm <- counter.glm + 1
  
}
df.glm <- ldply(results.glm, data.frame)


#### Iterate over thinning
counter.gamma <- 1
results.gamma <- list()
gammas <- c(0, -0.25, -0.5, -1)
for (gamma in gammas){
  print(paste("gamma:", gamma))
  
  
  #### True parameter values
  tau.sq.1 <- 0.35
  rho.1 <- 0.5
  beta0 <- 1
  beta1 <- 2
  beta.true <- c(beta0, beta1)
  
  
  #### Iterate over nsamp
  counter <- 1
  samps <- c(5, 10, 25, 50, 65)
  results.iterated = list()
  for (nsamp in samps){
    print(paste("nsamp =", nsamp))
    
    
    #### Simulate data
    data.1 <- simLerouxBias(caWc.disc, beta.true, gamma, tau.sq.1, rho.1, nNeighbors, seed=405)
    data.unif <- simUniform(caWc.disc, beta.true, nsamp=nsamp)
    data <- list(spatial=list(data.1), survey=data.unif)
    
    
    print(glm(data.1$y ~ data.1$x.standardised-1, family='poisson'))
    print(glm(data.unif$y ~ data.unif$x.standardised-1, family='poisson'))
    
    
    #### Multiple runs
    t.beta <- 0.05
    t.phi <- 1
    t.rho <- 0.05
    prior.tau <- c(1, 0.01)
    prior.beta <- rep(1000, data.1$p)
    
    n.runs <- 3
    par(mfrow=c(1,3))
    results <- list()
    for (i in 1:n.runs){
      print(paste("iteration", i))
      output <- carLerouxIntegrationUnif(
        data, t.beta, t.phi, t.rho, prior.var.beta=prior.beta, 
        prior.tau2=prior.tau, fix.rho=FALSE, n.sample=70000, self.tune=TRUE)
      output$nsamp <- nsamp
      output$gamma <- gamma
      results[[i]] <- output
      plot(output$samples.beta[,1], type='l', main='beta0')
      abline(h=beta0, col='2')
    }
    
    pooled <- pool.results(results)
    results.iterated[[counter]] <- pooled
    counter <- counter + 1
    
  }
  
  
  results.gamma[[counter.gamma]] <- results.iterated
  counter.gamma <- counter.gamma + 1
  
  
}


param.list <- list(
  beta0=beta0,
  beta1=beta1,
  tau.sq.1=tau.sq.1,
  rho.1=rho.1
)
biases <- list()
j <- 1
for (h in 1:length(results.gamma)){
  for (i in 1:length(results.gamma[[h]])){
    pooled.summary <- summarize(results.gamma[[h]][[i]], param.list)
    for (r in pooled.summary){
      r$nsamp <- samps[[i]]
      r$gamma <- gammas[[h]]
      biases[[j]] <- r
      j <- j + 1
    }
  }
}
df <- ldply(biases, data.frame)
list2csv(biases, 'sE_case3_summary_all.csv')
df.b0 <- df[df$parameter=='beta0',]
df.b1 <- df[df$parameter=='beta1',]
df.b0.sub <- df.b0[,names(df.b0) %in% c('nsamp', 'gamma', 'posterior.mean', 'percbias', 'posterior.sd')]
df.b1.sub <- df.b1[,names(df.b1) %in% c('nsamp', 'gamma', 'posterior.mean', 'percbias', 'posterior.sd')]
comp.b0 <- merge(df.b0.sub, df.glm, by='nsamp')
comp.b1 <- merge(df.b1.sub, df.glm, by='nsamp')
comp.b0 <- comp.b0[with(comp.b0, order(-gamma)),]
comp.b1 <- comp.b1[with(comp.b1, order(-gamma)),]
comp.b0 <- comp.b0[, !names(comp.b0) %in% c('glm.beta1', 'glm.pbias.beta1')]
comp.b1 <- comp.b1[, !names(comp.b1) %in% c('glm.beta0', 'glm.pbias.beta0')]
comp.b0$abs.pb.diff <- abs(comp.b0$percbias) - abs(comp.b0$glm.pbias.beta0)
comp.b1$abs.pb.diff <- abs(comp.b1$percbias) - abs(comp.b1$glm.pbias.beta1)
comp.b0 <- comp.b0[c('gamma', 'nsamp', 'posterior.mean', 'posterior.sd',
                     'glm.beta0', 'glm.sd.beta0', 'abs.pb.diff')]
comp.b1 <- comp.b1[c('gamma', 'nsamp', 'posterior.mean', 'posterior.sd',
                     'glm.beta1', 'glm.sd.beta1', 'abs.pb.diff')]
write.table(comp.b0, file='Documents/research/dataInt/sE_case3_summaryb0.csv', sep=',', row.names=F)
write.table(comp.b1, file='Documents/research/dataInt/sE_case3_summaryb1.csv', sep=',', row.names=F)


par(mfrow=c(1,1))
df.b0$gamma <- as.factor(df.b0$gamma)
ggplot(data=df.b0, aes(x=nsamp, y=percbias, group=gamma))+
  geom_line(aes(color=gamma))+
  geom_point(aes(color=gamma))+ 
  labs(x ='number of surveys')+
  labs(y='percent bias')


comp.b0$gamma <- as.factor(comp.b0$gamma)
ggplot(data=comp.b0, aes(x=nsamp, y=abs.pb.diff, group=gamma))+
  geom_line(aes(color=gamma))+
  geom_point(aes(color=gamma))+ 
  labs(x ='number of surveys')+
  labs(y='APBD')


comp.b1$gamma <- as.factor(comp.b1$gamma)
ggplot(data=comp.b1, aes(x=nsamp, y=abs.pb.diff, group=gamma))+
  geom_line(aes(color=gamma))+
  geom_point(aes(color=gamma))+ 
  labs(x ='number of surveys')+
  labs(y='APBD')


comp.b0$gamma <- as.factor(comp.b0$gamma)
ggplot(data=comp.b0, aes(x=nsamp, y=comp.b0$glm.sd.beta0 - comp.b0$posterior.sd, group=gamma))+
  geom_line(aes(color=gamma))+
  geom_point(aes(color=gamma))+ 
  labs(x ='number of surveys')+
  labs(y='SD Difference')


comp.b1$gamma <- as.factor(comp.b1$gamma)
ggplot(data=comp.b1, aes(x=nsamp, y=comp.b1$glm.sd.beta1 - comp.b1$posterior.sd, group=gamma))+
  geom_line(aes(color=gamma))+
  geom_point(aes(color=gamma))+ 
  labs(x ='number of surveys')+
  labs(y='SD Difference')
