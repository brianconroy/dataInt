############################
# Fit the shared latent 
# process model

# no locational covariates
# or intercept, just w
############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 4
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, beta.samp, seed=42) # simBernoulliLocGp(loc.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


#### Simulate counts given locations
Alpha <- 2
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
count.data <- simConditionalGp2(cov.disc, locs, beta.case, beta.samp, Alpha, W)
glm(count.data$y ~ count.data$x.standardised-1, family='poisson')
prior_alpha_mean <- 1
prior_alpha_var <- 6
# beta.case <- c(1)
# w.sub <- W[locs$ids]
# x <- array(1, c(sum(locs$status), 1))
# rates <- exp(x %*% beta.case + Alpha * w.sub)
# counts <- sapply(rates, function(x){rpois(n=1, x)})
# count.data <- list(x.standardised=x, y=counts)
# glm(counts ~ x-1, family='poisson')


U <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- x.c %*% beta.c + alpha * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.c)){
    logd <- logd + dpois(y.c[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


Ucase <- function(y, w, x, beta, alpha){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  prior_beta_mean <- rep(0, length(beta))
  if (length(beta) > 1){
    prior_beta_var <- diag(rep(100, length(beta)))
  } else {
    prior_beta_var <- matrix(100)
  }
  
  logd <- logd + dmvnorm(as.numeric(beta), prior_beta_mean, prior_beta_var, log=T)
  
  return(-logd)
  
}


Ualpha <- function(y, w, x, beta, alpha){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  prior_mean <- 1
  prior_var <- prior_alpha_mean
  
  logd <- logd + dnorm(alpha, prior_mean, prior_var, log=T)
  
  return(-logd)
  
}


dU <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma.inv, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(w[i])
  }
  
  # count contribution
  lin.count <- x.c %*% beta.c
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha * (y.c[i] - exp(lin.count[i] + alpha * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


dU.case <- function(y, w, x, beta, alpha){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + alpha * w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - exp(lin_preds[i,]))
    }
  }
  
  if (length(beta) > 1){
    prior_var_inv <- diag(rep(1/5, length(beta)))
  } else{
    prior_var_inv <- 1/5
  }
  
  grad <- grad + t(-t(beta) %*% prior_var_inv)
  return(-grad)
  
}


dU.alpha <- function(y, w, x, beta, alpha){
  
  grad <- 0
  lin_preds <- x %*% beta + alpha * w
  
  for (i in 1:length(w)){
    grad <- grad + w[i] * (y[i] - exp(lin_preds[i,]))
  }
  
  prior_var_inv <- 1/prior_alpha_var
  
  grad <- grad + alpha * prior_var_inv
  return(-grad)
  
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(y.l, x.c, y.c, alpha, beta.c, w, sigma, sigma.inv, loc.stats){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(y.l, x.c, y.c, alpha, beta.c, wcurr, sigma.inv, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(y.l, x.c, y.c, alpha, beta.c, wStar, sigma.inv, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(y.l, x.c, y.c, alpha, beta.c, wStar, sigma.inv, loc.stats)
  
  # evaluate energies
  U0 <- U(y.l, x.c, y.c, alpha, beta.c, wcurr, sigma, loc.stats)
  UStar <- U(y.l, x.c, y.c, alpha, beta.c, wStar, sigma, loc.stats)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    wnext <- wcurr
    accept <- 0
  } else {
    if (runif(1, 0, 1) < alpha){
      wnext <- wStar
      accept <- 1
    } else {
      wnext <- wcurr
      accept <- 0
    }
  }
  
  out <- list()
  out$w <- wnext
  out$accept <- accept
  return(out)
  
}


caseHmcUpdate <- function(y, w, x, beta, alpha){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_c * dU.case(y, w, x, bcurr, alpha)
  
  # first full step for position
  bStar <- bcurr + delta_c*pStar
  
  # full steps
  for (jL in 1:c(L_c-1)){
    # momentum
    pStar <- pStar - delta_c * dU.case(y, w, x, bStar, alpha)
    
    # position
    bStar <- bStar + delta_c*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_c * dU.case(y, w, x, bStar, alpha)
  
  # evaluate energies
  U0 <- Ucase(y, w, x, bcurr, alpha)
  UStar <- Ucase(y, w, x, bStar, alpha)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    bnext <- bcurr
    accept <- 0
  } else {
    if (runif(1, 0, 1) < alpha){
      bnext <- bStar
      accept <- 1
    } else {
      bnext <- bcurr
      accept <- 0
    }
  }
  
  out <- list()
  out$beta <- bnext
  out$accept <- accept
  return(out)
  
}


alphaHmcUpdate <- function(y, w, x, beta, alpha){
  
  
  # sample random momentum
  p0 <- rnorm(1)
  
  # simulate Hamiltonian dynamics
  acurr <- alpha
  pStar <- p0 - 0.5 * delta_a * dU.alpha(y, w, x, beta, acurr)
  
  # first full step for position
  aStar <- acurr + delta_a*pStar
  
  # full steps
  for (jL in 1:c(L_a-1)){
    # momentum
    pStar <- pStar - delta_a * dU.alpha(y, w, x, beta, aStar)
    
    # position
    aStar <- aStar + delta_a*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_a * dU.alpha(y, w, x, beta, aStar)
  
  # evaluate energies
  U0 <- Ualpha(y, w, x, beta, acurr)
  UStar <- Ualpha(y, w, x, beta, aStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  a <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (runif(1, 0, 1) < a){
    anext <- aStar
    accept <- 1
  } else {
    anext <- acurr
    accept <- 0
  }
  
  out <- list()
  out$alpha <- anext
  out$accept <- accept
  return(out)
  
}


X.c <- count.data$x.standardised
Y.c <- count.data$y
Y.l <- locs$status
N.w <- length(locs$status)

## starting values
w.i <- rnorm(length(W))
mod <- glm(Y.c ~ X.c + w.i[locs$ids] - 1, family='poisson')
alpha.i <- runif(1, 1, 3)
beta.c <- as.numeric(coef(mod)[1:2])
theta.i <- runif(1, 0, 4)
phi.i <- runif(1, 0, 4)
p.c <- length(beta.c)

n.sample <- 8000
accept <- rep(0, 4)

samples.w <- array(NA, c(n.sample, length(Y.l)))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.alpha <- array(NA, c(n.sample, 1))
samples.beta.c <- array(NA, c(n.sample, length(beta.c)))

delta <- 0.05
L <- 20
delta_c <- 0.01
L_c <- 10
delta_a <- 0.01
L_a <- 10

proposal.sd.theta <- 0.3
proposal.sd.alpha <- 0.05

prior.phi <- c(4, 12)
prior.theta <- c(2.5, 2.5)
prior.mean.beta <- rep(0, p.c)
prior.var.beta <- rep(1000, p.c)

progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*n.sample)

for (i in 1:n.sample){
  
  ## sample from w
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  sigma.inv.i <- solve(sigma.i)
  w.out.i <- wHmcUpdate(Y.l, X.c, Y.c, alpha.i, beta.c, w.i, sigma.i, sigma.inv.i, locs)
  w.i <- w.out.i$w

  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1/rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  
  ## sample from beta.c
  w.i.sub <- w.i[locs$ids]
  beta.out.c <- caseHmcUpdate(Y.c, w.i[locs$ids], X.c, beta.c, alpha.i)
  beta.c <- beta.out.c$beta
  
  ## sample from alpha
  alpha.out <- alphaHmcUpdate(Y.c, w.i.sub, X.c, beta.c, alpha.i)
  alpha.i <- alpha.out$alpha

  samples.beta.c[i,] <- beta.c
  samples.alpha[i,] <- alpha.i
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.w[i,] <- t(w.i)
  
  accept[1] <- accept[1] + w.out.i$accept
  accept[2] <- accept[2] + theta.out$accept
  accept[3] <- accept[3] + beta.out.c$accept
  accept[4] <- accept[4] + alpha.out$accept
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}


accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
glm(Y.c ~ X.c + w.hat[locs$ids] -1, family='poisson')

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

par(mfrow=c(1,2))
plot(samples.theta, type='l', ylab='Theta'); abline(h=Theta, col='2')
hist(samples.theta, xlab='Theta', main=''); abline(v=Theta, col='2')
print(mean(samples.theta)); print(Theta)

plot(samples.phi, type='l', ylab='Phi'); abline(h=Phi, col='2')
hist(samples.phi, xlab='Phi', main=''); abline(v=Phi, col='2')
print(mean(samples.phi)); print(Phi)

plot(samples.beta.c[,1], type='l'); abline(h=beta.case[1], col='2')
lines(apply(samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
plot(apply(samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')

par(mfrow=c(1,2))
plot(samples.beta.c[,1], type='l', ylab='beta0'); abline(h=beta.case[1], col='2')
plot(samples.beta.c[,2], type='l', ylab='beta1'); abline(h=beta.case[2], col='2')

hist(samples.beta.c[,1]); abline(v=beta.case[1], col='2')
plot(samples.beta.c[,2], type='l'); abline(h=beta.case[2], col='2')
hist(samples.beta.c[,2]); abline(v=beta.case[2], col='2')
print(colMeans(samples.beta.c))

par(mfrow=c(1,2))
plot(samples.alpha, type='l', xlab='iteration', ylab='Alpha'); abline(h=Alpha, col=2)
hist(samples.alpha, xlab='Alpha'); abline(v=Alpha, col=2)
print(mean(samples.alpha))
