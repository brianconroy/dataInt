library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#######################################
# Tunes the BYM model on simulated data
#######################################


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
caWc.disc <- caWc.disc[[1]]
nNeighbors <- simRegion$nNeighbors


#### True parameter values
tau.sq <- 0.5
sigma.sq <- 0.5
beta0 <- 0.5
beta1 <- 3
beta.true <- c(beta0, beta1)


#### Simulate data
data <- simBYM(caWc.disc, beta.true, tau.sq, sigma.sq, nNeighbors)


#########################
## Tune Proposal & Priors
#########################


tune.theta <- c(0.015, 0.07)
tune.beta <- c(0.015, 0.025)
tune.phi <- c(1, 1.5)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data$p))

par(mfrow=c(4, 2))
results <- list()
counter <- 1
for (i in 1:length(tune.beta)){
  for(j in 1:length(tune.phi)){
    for (k in 1:length(tune.tau.prior)){
      for (l in 1:length(tune.beta.prior)){
        for (m in 1:length(tune.theta)){
        
        
          print("***************")
          t.beta <- tune.beta[i]
          t.phi <- tune.phi[j]
          t.theta <- tune.theta[m]
          prior.tau <- tune.tau.prior[[k]]
          prior.beta <- tune.beta.prior[[l]]
          print(paste("tuning beta prop:", t.beta, "phi prop:", t.phi, "theta prop:", t.theta))
          print(paste("phi prior params:", prior.tau[1], prior.tau[2]))
          print(paste("beta prior params:", prior.beta[1], prior.beta[2]))
          
          output <- carBYM(data, t.beta, t.phi, t.theta, prior.var.beta=prior.beta, prior.tau2=prior.tau)
          output$tune.beta <- t.beta
          output$tune.phi <- t.phi
          output$tune.tau.prior <- prior.tau
          output$tune.beta.prior <- prior.beta
          output$tune.theta <- t.theta
          
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
}


## Write results to csv
df <- c()
for (r in results){
  
  beta0.mean <- round(mean(r$samples.beta[,1]), 3)
  beta1.mean <- round(mean(r$samples.beta[,2]), 3)
  tau2.mean <- round(mean(r$samples.tau2), 3)
  sigma2.mean <- round(mean(r$samples.sigma2), 3)
  t.plus.s <- tau2.mean + sigma2.mean
  
  beta0.bias <- round(beta0.mean - beta0, 3)
  beta1.bias <- round(beta1.mean - beta1, 3)
  tau2.bias <- round(tau2.mean- tau.sq, 3)
  
  beta0.bias.prc <- round(100*beta0.bias/beta0, 3)
  beta1.bias.prc <- round(100*beta1.bias/beta1, 3)
  tau2.bias.prc <- round(100*tau2.bias/tau.sq, 3)
  
  row <- c(r$tune.beta, r$tune.phi, r$tune.theta,
           r$tune.tau.prior[1], r$tune.beta.prior[1],
           round(r$accept[1], 3), round(r$accept[2], 3),
           beta0.mean, beta1.mean, tau2.mean, sigma2.mean, t.plus.s,
           beta0.bias, beta1.bias, tau2.bias,
           beta0.bias.prc, beta1.bias.prc, tau2.bias.prc)
  df <- rbind(df, row)
  
}
df <- data.frame(df)
names(df) <- c('tune.beta', 'tune.phi', 'tune.theta', 'tune.tau.prior',
               'tune.beta.prior', 'accept.beta', 'accept.phi', 
               'beta0', 'beta1', 'tau2', 'sigma2', 'tau2 + sigma2',
               'beta0.bias', 'beta1.bias', 'tau2.bias',
               'beta0.bias.prc', 'beta1.bias.prc', 'tau2.bias.prc')
write.table(df, file='/Users/brianconroy/Documents/research/dataInt/tune_results.csv', 
            sep=',', row.names=F)


## Acceptance rates
for (r in results){
  print("***************")
  print("tuning params")
  print(paste("beta tune: ", r$tune.beta))
  print(paste("phi tune: ", r$tune.phi))
  print(paste("theta tune: ", r$tune.theta))
  print(paste("tau prior: ", r$tune.tau.prior))
  print(paste("beta prior: ", r$tune.beta.prior))
  print("acceptance rates (beta, phi, theta)")
  print(r$accept)
  print(paste("beta0", mean(r$samples.beta[,1])))
  print(paste("beta1", mean(r$samples.beta[,2])))
  print(paste("tau2", mean(r$samples.tau2)))
}


#############
## Traceplots
#############


par(mfrow=c(4, 2))
viewTraces(results, 'beta0', line=beta0)
viewTraces(results, 'beta1', line=beta1)
viewTraces(results, 'tau2', line=tau.sq)
viewTraces(results, 'sigma2', line=sigma.sq)

viewHists(results, 'beta0', line=beta0)
viewHists(results, 'beta1', line=beta1)
viewHists(results, 'beta2', line=beta2)
viewHists(results, 'tau2', line=tau.sq)

