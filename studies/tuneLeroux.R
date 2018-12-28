library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


##############################################
# Tunes the Leroux CAR model on simulated data
##############################################


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
rho <- 0.5
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1, beta2)


#### Simulate data
data <- simLeroux(caWc.disc, beta.true, tau.sq, rho, nNeighbors, seed=401)
glm(data$y ~ data$x.standardised-1, family='poisson')


#########################
## Tune Proposal & Priors
#########################


# tune.beta <- c(0.005, 0.05, 0.07)
# tune.phi <- c(0.5, 1, 1.5)


tune.beta <- c(0.020, 0.025)
tune.phi <- c(0.5)
tune.rho <- c(1)
tune.tau.prior <- list(c(1, 0.01))
tune.beta.prior <- list(rep(1000, data$p))
rho.ini <- rho


par(mfrow=c(3, 2))
results <- list()
counter <- 1
for (i in 1:length(tune.beta)){
  for(j in 1:length(tune.phi)){
    for (k in 1:length(tune.tau.prior)){
      for (l in 1:length(tune.beta.prior)){
        for (m in 1:length(tune.rho)){
        
        
          print("***************")
          t.beta <- tune.beta[i]
          t.phi <- tune.phi[j]
          t.rho <- tune.rho[m]
          prior.tau <- tune.tau.prior[[k]]
          prior.beta <- tune.beta.prior[[l]]
          print(paste("tuning beta prop:", t.beta, "phi prop:", t.phi, "rho prop:", t.rho))
          print(paste("phi prior params:", prior.tau[1], prior.tau[2]))
          print(paste("beta prior params:", prior.beta[1], prior.beta[2]))
          
          output <- carLeroux(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau,
                              fix.rho=TRUE, rho=rho, n.sample=70000)
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


## Acceptance rates
for (r in results){
  print("***************")
  print("tuning params")
  print(paste("beta tune: ", r$tune.beta))
  print(paste("phi tune: ", r$tune.phi))
  print(paste("tau prior: ", r$tune.tau.prior))
  print(paste("beta prior: ", r$tune.beta.prior))
  print("acceptance rates (beta, phi, rho)")
  print(r$accept)
  print(paste("beta0", mean(r$samples.beta[,1])))
  print(paste("beta1", mean(r$samples.beta[,2])))
  print(paste("tau2", mean(r$samples.tau2)))
}


#############
## Traceplots
#############


par(mfrow=c(2, 1))
viewTraces(results, 'beta0', line=beta0)
viewTraces(results, 'beta1', line=beta1)
viewTraces(results, 'beta2', line=beta2)
viewTraces(results, 'tau2', line=tau.sq)
viewTraces(results, 'rho', line=rho)

viewHists(results, 'beta0', line=beta0)
viewHists(results, 'beta1', line=beta1)
viewHists(results, 'beta2', line=beta2)
viewHists(results, 'tau2', line=tau.sq)
viewHists(results, 'rho', line=rho)


###############
## Best Configs
###############


# best configs (rho = 0.75)
# beta tune 0.05	
# tau tune 3
# tau prior 1	
# beta prior 1000

# best configs (rho = 1)
# beta tune 0.015	
# tau tune 2.25	
# tau prior 1	
# beta prior 1000

# best configs (multi params, rho=0.75)
# beta tune 0.015
# phi tune 1
# tau prior 1	
# beta prior 1000


################
## Multiple Runs
################
t.beta <- 0.015
t.phi <- 1.5
t.rho <- 1
prior.tau <- NULL
prior.beta <- NULL


## Random GLM Starting Values
n.runs <- 4
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLeroux(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau,
                      fix.rho=FALSE, proposal.sd.rho=t.rho, rho=rho.ini)
  results[[i]] <- output
}


par(mfrow=c(2,2))
viewTraces(results, 'beta0', line=beta0)
viewTraces(results, 'beta1', line=beta1)
viewTraces(results, 'beta2', line=beta2)
viewTraces(results, 'tau2', line=tau.sq)
viewTraces(results, 'rho', line=rho)

viewHists(results, 'beta0', line=beta0)
viewHists(results, 'beta1', line=beta1)
viewHists(results, 'beta2', line=beta2)
viewHists(results, 'tau2', line=tau.sq)
viewHists(results, 'rho', line=rho)


##################################
## Different Data Sim Realizations
##################################
t.beta <- 0.015
t.phi <- 2.25
prior.tau <- NULL
prior.beta <- NULL


n.runs <- 4
results <- list()
par(mfrow=c(2,2))
for (i in 1:n.runs){
  print(paste("iteration", i))
  data <- simLeroux(caWc.disc, beta.true, tau.sq, rho, nNeighbors)
  output <- carLeroux(data, t.beta, t.phi, prior.var.beta=prior.beta, prior.tau2=prior.tau)
  
  results[[i]] <- output
}


par(mfrow=c(2,2))
viewTraces(results, 'beta0', line=beta0)
viewTraces(results, 'beta1', line=beta1)
viewTraces(results, 'tau2', line=tau.sq)

viewHists(results, 'beta0', line=beta0)
viewHists(results, 'beta1', line=beta1)
viewHists(results, 'tau2', line=tau.sq)
