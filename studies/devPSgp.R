############################
# Hamiltonian MCMC to sample
# from case only preferential
# sampling model.

# gaussian process shared
# between model components.
# single alpha
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
# Theta <- 5
# Phi <- 1
# cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
# coords <- xyFromCell(caWc.disc, cell=cells.all)
# d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
# Sigma <- Exponential(d, range=Theta, phi=Phi)
# W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
# hist(W)

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
# beta.samp <- c(-1, 2)
# loc.disc <- caWc.disc[[c(1)]]
# locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W)
# sum(locs$status)
# hist(locs$w)
# plot(loc.disc)
# points(locs$coords)
beta.samp <- c(1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


#### Simulate counts given locations
Alpha <- 1
# beta.case <- c(1, 1)
beta.case <- c(1)
cov.disc <- caWc.disc[[c(1)]]
# count.data <- simConditionalGp2(cov.disc, locs, beta.case, beta.samp, Alpha, W)
w.sub <- W[locs$ids]
x <- array(1, c(sum(locs$status), 1))
rates <- exp(x %*% beta.case + Alpha * w.sub)
counts <- sapply(rates, function(x){rpois(n=1, x)})

count.data <- list(x.standardised=x, y=counts)
glm(counts ~ x-1, family='poisson')


U <- function(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, w, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- x.l %*% beta.l + w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- x.c %*% beta.c + alpha * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.c)){
    # logd <- logd + dpois(y.c[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU <- function(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, w, sigma.inv, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  lin.loc <- x.l %*% beta.l
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(lin.loc[i] + w[i])
  }
  
  # count contribution
  lin.count <- x.c %*% beta.c
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    # grad[id.i] <- grad[id.i] + alpha * (y.c[i] - exp(lin.count[i] + alpha * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdate <- function(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, w, sigma, sigma.inv, loc.stats){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, wcurr, sigma.inv, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, wStar, sigma.inv, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, wStar, sigma.inv, loc.stats)
  
  # evaluate energies
  U0 <- U(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, wcurr, sigma, loc.stats)
  UStar <- U(x.l, y.l, x.c, y.c, alpha, beta.l, beta.c, wStar, sigma, loc.stats)
  
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


X.c <- count.data$x.standardised
X.l <- locs$x.scaled
Y.c <- count.data$y
Y.l <- locs$status
N.w <- length(locs$status)
if (ncol(X.l) > 1){
  X.l.sub <- X.l[as.logical(locs$status),]
} else{
  X.l.sub <- matrix(X.l[as.logical(locs$status),])
}

alpha.i <- Alpha # runif(1, 0, 4)
theta.i <- runif(1, 0, 4) # Theta
phi.i <- runif(1, 0, 4) # Phi
w.i <- rnorm(length(cells.all))
beta.c <- beta.case # rnorm(ncol(X.c), 0, 2)
beta.l <- beta.samp # rnorm(ncol(X.l), 0, 2) # 
p.c <- length(beta.c)
p.l <- length(beta.l)

n.sample <- 5000
accept <- rep(0, 6)

samples.w <- array(NA, c(n.sample, length(Y.l)))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.alpha <- array(NA, c(n.sample, 1))
samples.beta.l <- array(NA, c(n.sample, length(beta.l)))
samples.beta.c <- array(NA, c(n.sample, length(beta.c)))

delta <- 0.2
L <- 20
proposal.sd.theta <- 0.3
proposal.sd.alpha <- 0.05
proposal.sd.beta.l <- 0.1
proposal.sd.beta.l.i <- 1
proposal.sd.beta.c <- 0.04
proposal.sd.beta.c.i <- 0.4

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
  w.out.i <- wHmcUpdate(X.l, Y.l, X.c, Y.c, alpha.i, beta.l, beta.c, w.i, sigma.i, sigma.inv.i, locs)
  w.i <- w.out.i$w

  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  
  ## sample from beta.l
  beta.out.l.i <- betaLogisticIntercept(X.l, Y.l, beta.l, proposal.sd.beta.l.i, w.i)
  beta.l <- beta.out.l.i$beta
  # beta.out.l.s <- betaLogisticSlope(X.l, Y.l, beta.l, proposal.sd.beta.l, w.i)
  # beta.l <- beta.out.l.s$beta
  # beta.l[2] <- beta.samp[2]
  
  ## sample from alpha
  # w.i.sub <- w.i[locs$ids]
  # alpha.out <- alphaPoissonUpdate(X.c, Y.c, alpha.i, w.i.sub, beta.c, proposal.sd.alpha)
  # alpha.i <- alpha.out$alpha
  
  ## sample from beta.c
  w.i.sub <- w.i[locs$ids]
  # beta.out.c <- betaPoissonUpdate(X.c, Y.c, beta.c, Alpha * w.i.sub, prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
  # beta.c <- beta.out.c$beta
  # beta.out.c.i <- betaPoissonIntercept(X.c, Y.c, beta.c, Alpha * w.i.sub, prior.mean.beta, prior.var.beta, proposal.sd.beta.c.i)
  # beta.c <- beta.out.c.i$beta
  #beta.out.c.s <- betaPoissonSlope(X.c, Y.c, beta.c, Alpha * w.i.sub, prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
  #beta.c <- beta.out.c.s$beta
  beta.c <- beta.case

  samples.beta.l[i,] <- beta.l
  samples.beta.c[i,] <- beta.c
  samples.alpha[i,] <- alpha.i
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.w[i,] <- t(w.i)
  
  accept[1] <- accept[1] + w.out.i$accept
  accept[2] <- accept[2] + theta.out$accept
  accept[3] <- accept[3] + beta.out.l.i$accept
  # accept[4] <- accept[4] + beta.out.l.s$accept
  # accept[5] <- accept[5] + beta.out.c.i$accept
  # accept[6] <- accept[6] + beta.out.c.s$accept
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}


glm(Y.c ~ X.c + w.hat[locs$ids] -1, family='poisson')

accept <- accept/n.sample
print(accept)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

plot(samples.theta, type='l'); abline(h=Theta, col='2')
hist(samples.theta); abline(v=Theta, col='2')
print(mean(samples.theta))

plot(samples.phi, type='l'); abline(h=Phi, col='2')
hist(samples.phi); abline(v=Phi, col='2')
print(mean(samples.phi))

ids.neg <- c(1:length(W))[W<0]
for(i in ids.neg[1:16]){
  plot(samples.w[,i], type='l'); abline(h=W[i], col='2')
}

plot(samples.alpha, type='l'); abline(h=Alpha, col='2')
hist(samples.alpha); abline(v=Alpha, col='2')
print(mean(samples.alpha))

par(mfrow=c(1,2))
plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col='2')
plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col='2')


plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col='2')
lines(apply(samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')



hist(samples.beta.l[,1])
abline(v=beta.samp[1], col='2')

hist(samples.beta.l[,2])
abline(v=beta.samp[2], col='2')

plot(samples.beta.c[,1], type='l'); abline(h=beta.case[1])
lines(apply(samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
hist(samples.beta.c[,1]); abline(v=beta.case[1], col=2)



plot(samples.beta.c[,2], type='l'); abline(h=beta.case[2], col='2')

print(colMeans(samples.beta.c))
