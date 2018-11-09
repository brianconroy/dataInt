library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
sourceDirectory('Documents/research/dataInt/R/')


##############################
## survey locations determined
## by Gaussian process
##############################


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 3
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
w <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)


hist(w)
w.disc <- caWc.disc[[1]]
w.disc[][!is.na(w.disc[])] <- w
plot(w.disc)


#### Simulate locations
beta.samp <- c(-1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=w, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


##############################
## First test - just update w
## beta fixed, cov params fixed
##############################


probLogisticUpdate <- function(logit_prob, y, proposal_sd, x, beta, w){
  
  logit_prob.new <- logit_prob
  accept <- 0
  prior_mean <- x %*% beta + w
  
  for(i in 1:length(prob)){
    
    logit_prob.prop <- rnorm(1, logit_prob[i], proposal_sd)
    prob.prop <- expit(logit_prob.prop)
    prob.curr <- expit(logit_prob[i])
    
    ll.prop <- dbinom(y[i], 1, prob.prop, log=TRUE)
    ll.curr <- dbinom(y[i], 1, prob.curr, log=TRUE)
    
    prior_mean_i <- prior_mean[i,]
    pr.prop <- dnorm(logit_prob.prop, prior_mean_i, 2, log=TRUE)
    pr.curr <- dnorm(logit_prob[i], prior_mean_i, 2, log=TRUE)
    
    acceptance <- exp(ll.prop + pr.prop - ll.curr - pr.curr)
    if(runif(1) <= acceptance) {
      logit_prob.new[i] <- logit_prob.prop
      accept <- accept + 1
    }
    
  }
  
  out <- list()
  out$logit_prob <- logit_prob.new
  out$accept <- accept
  return(out)
  
}


subdiv_matrix <- function(dmat, sigma, i, n){
  
  top <- sort(dmat[i,])[0:(n+1)]
  ids <- as.numeric(names(top))
  return(ids)
 
}


subdiv_ids <- list()
for (i in 1:length(w)){
  subdiv_ids[[i]] <- subdiv_matrix(d, sigma, i, 10)
}


subdiv_sigmas <- list()
for (i in 1:length(w)){
  ids_i <- subdiv_ids[[i]]
  subdiv_sigmas[[i]] <- sigma[ids_i, ids_i]
}


wNormalUpdate <- function(logit_prob, w, mu, sigma, proposal.sd){
  
  w.new <- w
  accept <- 0
  
  for(i in 1:length(w)){
    
    w.prop <- rnorm(1, w[i], proposal.sd)
    lp.curr <- w[i] + mu[i]
    lp.prop <- w.prop + mu[i]
    
    like.curr <- dnorm(logit_prob[i], lp.curr, sd=1, log=T)
    like.prop <- dnorm(logit_prob[i], lp.prop, sd=1, log=T)
    like.diff <- like.prop - like.curr
    
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


#### Fit
x.l <- locs$x.scaled
x.l.sub <- x.l[as.logical(locs$status),]
y.l <- locs$status

n.sample <- 15000
N <- nrow(d)

beta.l <- beta.samp
w.i <- w
prob <- runif(405, 0, 1)
logit_prob <- logit(prob)

proposal.sd.beta.l <- 0.2
proposal.sd.w <- 0.5
proposal_sd.prob <- 0.5

samples.beta.l <- array(NA, c(n.sample, length(beta.samp)))
samples.w <- array(NA, c(n.sample, length(w)))
samples.logit <- array(NA, c(n.sample, length(logit_prob)))

accept <- c(0, 0, 0)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  # sample from logit prob
  logit.out <- probLogisticUpdate(logit_prob, y.l, proposal_sd.prob, x.l, beta.l, w.i)
  logit_prob <- logit.out$logit_prob
  
  ## Sample from w
  sigma.i <- sigma
  mu.i <- x.l %*% beta.l
  w.out.i <- wNormalUpdate(logit_prob, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w.out.i$w
  w.i <- w.i - mean(w.i)
  
  accept[1] <- accept[1] + logit.out$accept
  accept[2] <- accept[2] + w.out.i$accept
  samples.logit[i,] <- logit_prob
  samples.w[i,] <- w.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
  k <- i/100
  if(ceiling(k)==floor(k) & k < 2000){
    proposal_sd.prob <- tuners.sd.1(accept[1], i * N, proposal_sd.prob, 45, 55)
    proposal.sd.w <- tuners.sd.1(accept[2], i * N, proposal.sd.w, 45, 55)
  }
  
}

accept <- accept/(n.sample * N)
print(accept)

logit.hat <- colMeans(samples.logit)
hist(logit.hat)

w.hat <- colMeans(samples.w)
hist(w.hat)
plot(w, w.hat)
abline(0,1, col='2')

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=w[j*i], col='2')
}

true.logit <- x.l %*% beta.l + w
par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.logit[,j*i], type='l'); abline(h=true.logit[j*i], col=2)
}

# try expanding n
# try adaptive samping (in burnin period)
# try hamiltonian mcmc


##############################
## Test 2) - just update w
## do you really have to 
## MHRW the logit quantity?

## gives about the same estimated
## ws, slightly worse
##############################


#### Fit
x.l <- locs$x.scaled
x.l.sub <- x.l[as.logical(locs$status),]
y.l <- locs$status

beta.l <- beta.samp
w.i <- w
prob <- runif(405, 0, 1)
logit_prob <- logit(prob)

proposal.sd.beta.l <- 0.2
proposal.sd.w <- 0.5

samples.beta.l <- array(NA, c(n.sample, length(beta.samp)))
samples.w <- array(NA, c(n.sample, length(w)))

n.sample <- 15000
accept <- c(0, 0, 0)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  sigma.i <- sigma
  mu.i <- x.l %*% beta.l
  w.out.i <- wLogisticUpdate(y.l, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w.out.i$w
  w.i <- w.i - mean(w.i)
  
  accept[2] <- accept[2] + w.out.i$accept
  samples.logit[i,] <- logit_prob
  samples.w[i,] <- w.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}


w.hat <- colMeans(samples.w)
hist(w.hat)
plot(w, w.hat)
abline(0,1, col='2')

par(mfrow=c(4, 4))
j <- 2
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=w[j*i], col='2')
}



##############################
## Test 3) - update covariance
## parameters
##############################


#### Fit
n.sample <- 15000
accept <- c(0, 0, 0)
N <- nrow(d)

x.l <- locs$x.scaled
x.l.sub <- x.l[as.logical(locs$status),]
y.l <- locs$status

beta.l <- beta.samp
w.i <- w
prob <- runif(405, 0, 1)
logit_prob <- logit(prob)
phi.i <- Phi
theta.i <- Theta

prior.phi <- c(3, 6)
prior.theta <- c(2.2, 2.2)

proposal.sd.beta.l <- 0.2
proposal.sd.w <- 0.5
proposal_sd.prob <- 0.5
proposal.sd.theta <- 0.1

samples.beta.l <- array(NA, c(n.sample, length(beta.samp)))
samples.w <- array(NA, c(n.sample, length(w)))
samples.logit <- array(NA, c(n.sample, length(logit_prob)))
samples.phi <- array(NA, c(n.sample, 1))
samples.theta <- array(NA, c(n.sample, 1))

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  # sample from logit prob
  logit.out <- probLogisticUpdate(logit_prob, y.l, proposal_sd.prob, x.l, beta.l, w.i)
  logit_prob <- logit.out$logit_prob
  
  ## sample from w
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  mu.i <- x.l %*% beta.l
  w.out.i <- wNormalUpdate(logit_prob, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w.out.i$w
  w.i <- w.i - mean(w.i)
  
  ## sample from Theta
  theta.out <- rangeMhUpdate(theta.i, w.i, d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  samples.phi[i,] <- phi.i
  
  accept[1] <- accept[1] + logit.out$accept
  accept[2] <- accept[2] + w.out.i$accept
  accept[3] <- accept[3] + theta.out$accept
  samples.logit[i,] <- logit_prob
  samples.w[i,] <- w.i
  samples.theta[i,] <- theta.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

accept <- accept/n.sample
accept[1:2] <- accept[1:2]/(N)
print(accept)

logit.hat <- colMeans(samples.logit)
hist(logit.hat)

w.hat <- colMeans(samples.w)
hist(w.hat)
plot(w, w.hat)
abline(0,1, col='2')

par(mfrow=c(4, 4))
j <- 2
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=count.data$w[j*i], col='2')
}

plot(samples.theta, type='l')
abline(h=Theta, col='2')


## Check Theta
theta.i <- Theta
prior.theta <- c(2.2, 2.2)
proposal.sd.theta <- 0.1
samples.theta <- array(NA, c(n.sample, 1))

w.i <- w
phi.i <- Phi

n.sample <- 15000
accept <- c(0, 0, 0)
N <- nrow(d)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  theta.out <- rangeMhUpdate(theta.i, w.i, d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  samples.theta[i,] <- theta.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}


