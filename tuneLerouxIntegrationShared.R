library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#########################################
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
caWc.disc <- caWc.disc[[c(1, 3)]]
nNeighbors <- simRegion$nNeighbors


#### True parameter values
tau.sq <- 0.5
rho <- 0.75
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1, beta2)


#### Simulate data
data.1 <- simLeroux(caWc.disc, beta.true, tau.sq, rho, nNeighbors, seed=405)
data.2 <- simLeroux(caWc.disc, beta.true, tau.sq, rho, nNeighbors, seed=400)
data <- list(data.1, data.2)
hist(data.1$y)
hist(data.2$y)


#### reference
glm(data.1$y ~ data.1$x.standardised-1, family='poisson')
glm(data.2$y ~ data.2$x.standardised-1, family='poisson')


#########################
## Tune Proposal & Priors
#########################


tune.beta <- c(0.01, 0.3)
tune.phi <- c(0.05, 0.005)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data.1$p))

par(mfrow=c(2, 2))
results <- list()
counter <- 1
for (i in 1:length(tune.beta)){
  for(j in 1:length(tune.phi)){
    for (k in 1:length(tune.tau.prior)){
      for (l in 1:length(tune.beta.prior)){
        
        
        print("***************")
        t.beta <- tune.beta[i]
        t.phi <- tune.phi[j]
        prior.tau <- tune.tau.prior[[k]]
        prior.beta <- tune.beta.prior[[l]]
        print(paste("tuning beta prop:", t.beta, "phi prop:", t.phi))
        print(paste("phi prior params:", prior.tau[1], prior.tau[2]))
        print(paste("beta prior params:", prior.beta[1], prior.beta[2]))
        
        output <- carLerouxIntegrationShared(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau)
        output$tune.beta <- t.beta
        output$tune.phi <- t.phi
        output$tune.tau.prior <- prior.tau
        output$tune.beta.prior <- prior.beta
        
        results[[counter]] <- output
        counter <- counter + 1
        
        print("acceptance rates (beta, phi)")
        print(output$accept)
        print(paste("beta0", mean(output$samples.beta[,1])))
        print(paste("beta1", mean(output$samples.beta[,2])))
        print(paste("tau2", mean(output$samples.tau2)))
        
        
        plot(output$samples.beta[,1], type='l', 
             main=paste('beta0 | ', 'phi tune:', t.phi, '| beta tune:', t.beta,
                        '| beta prior:', prior.beta[1], '| tau2 prior:', prior.tau[1]))
        abline(h=beta0, col='2')
      }
    }
  }
}


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


#############
## Traceplots
#############

overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho), collapse=" ")
par(mfrow=c(3, 3), oma=c(0,0,2,0))
viewTracesInt(results, 'beta0', line=beta0, overall.title=paste("beta0 |",overall.main))
viewTracesInt(results, 'beta1', line=beta1, overall.title=paste("beta1 |",overall.main))
viewTracesInt(results, 'beta2', line=beta2, overall.title=paste("beta2 |",overall.main))
viewTracesInt(results, 'tau2', line=tau.sq.1, d=1, overall.title=paste("tau2 (d=1) |",overall.main))
viewTracesInt(results, 'tau2', line=tau.sq.2, d=2, overall.title=paste("tau2 (d=2) |",overall.main))
viewTracesInt(results, 'tau2', line=tau.sq.2, d=3, overall.title=paste("tau2 (d=3) |",overall.main))

viewHists(results, 'beta0', line=beta0)
viewHists(results, 'beta1', line=beta1)
viewHists(results, 'beta2', line=beta2)
viewHists(results, 'tau2', line=tau.sq)


################
## Multiple Runs
################


# same CAR params
# t.beta <- 0.01
# t.phi <- 0.25
# prior.tau <- NULL
# prior.beta <- NULL


# different CAR params
# t.beta <- 0.025
# t.phi <- 0.1
# prior.tau <- NULL
# prior.beta <- NULL


# different CAR params, d=3
t.beta <- 0.025
t.phi <- 0.15
prior.tau <- NULL
prior.beta <- NULL


## Random GLM Starting Values
n.runs <- 6
par(mfrow=c(3,2))
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLerouxIndepIntegration(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau)
  results[[i]] <- output
  
  plot(output$samples.beta[,2], type='l', main='beta1')
  abline(h=beta1, col='2')
}


overall.main <- paste(c("D:", length(data), "| beta:", beta.true, "| rho's:", rho.list), collapse=" ")
par(mfrow=c(3, 2), oma=c(0,0,2,0))
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

