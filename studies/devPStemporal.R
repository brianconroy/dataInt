############################
# Fit the spatiotemporal
# shared latent process model
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
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, beta=0, seed=42) # simBernoulliLocGp(loc.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)


#### Simulate counts given locations
## Case Counts
Alpha.case <- 2
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W)
glm(case.data$y ~ case.data$x.standardised-1, family='poisson')
prior_alpha_ca_mean <- 2
prior_alpha_ca_var <- 6

Alpha.ctrl <- -2
beta.ctrl <- c(2, 0.5)
cov.disc <- caWc.disc[[c(1)]]
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W)
glm(ctrl.data$y ~ ctrl.data$x.standardised-1, family='poisson')
prior_alpha_co_mean <- -2
prior_alpha_co_var <- 6


Uw_cc <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma, loc.stats){
  
  logd <- 0
  
  # likelihood: locations
  loc_pred <- w
  for (i in 1:length(y.l)){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(loc_pred[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- w[as.logical(loc.stats$status)]
  count_pred <- x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # likelihood: control counts
  count_pred <- x.c %*% beta.co + alpha.co * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.co)){
    logd <- logd + dpois(y.co[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


dU_w_cc <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma.inv, loc.stats){
  
  grad <- array(0, c(length(w), 1))
  
  # location contribution
  for (i in 1:length(w)){
    grad[i] <- grad[i] + y.l[i] - expit(w[i])
  }
  
  # case contribution
  lin.count <- x.c %*% beta.ca
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # control contribution
  lin.count <- x.c %*% beta.co
  for (i in 1:length(loc.stats$ids)){
    id.i <- loc.stats$ids[i]
    grad[id.i] <- grad[id.i] + alpha.co * (y.co[i] - exp(lin.count[i] + alpha.co * w[id.i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


wHmcUpdateCC <- function(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, w, sigma, sigma.inv, loc.stats, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wcurr, sigma.inv, loc.stats)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma.inv, loc.stats)
  
  # evaluate energies
  U0 <- Uw_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wcurr, sigma, loc.stats)
  UStar <- Uw_cc(y.l, x.c, y.ca, alpha.ca, beta.ca, y.co, alpha.co, beta.co, wStar, sigma, loc.stats)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  if (is.na(alpha)){
    alpha <- 0
  }
  
  if (runif(1, 0, 1) < alpha){
    wnext <- wStar
    accept <- 1
  } else {
    wnext <- wcurr
    accept <- 0
  }
  
  out <- list()
  out$w <- wnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


X.c <- case.data$x.standardised
Y.ca <- case.data$y
Y.co <- ctrl.data$y
Y.l <- locs$status
N.w <- length(locs$status)

## starting values
# w.i <- W
# alpha.ca.i <- Alpha.case
# beta.ca <- beta.case
# alpha.co.i <- Alpha.ctrl
# beta.co <- beta.ctrl
# theta.i <- Theta
# phi.i <- Phi
# p.c <- length(beta.ca)

w.i <- W + rnorm(length(W))
beta.ca <- beta.case + rnorm(length(beta.case))
theta.i <- Theta + rnorm(1)
phi.i <- Phi + rnorm(1)
p.c <- length(beta.ca)

# w.i <- rnorm(length(W))
# alpha.ca.i <- runif(1, 1, 3)
# beta.ca <- rnorm(2, 1)
# alpha.co.i <- runif(1, -3, -1)
# beta.co <- rnorm(2, 1)
# theta.i <- runif(1, 0, 4)
# phi.i <- runif(1, 0, 4)
# p.c <- length(beta.ca)

n.sample <- 12000
accept <- rep(0, 6)

samples.w <- array(NA, c(n.sample, length(Y.l)))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.beta.ca <- array(NA, c(n.sample, length(beta.ca)))

L <- 20
L_ca <- 25


w_tuning <- initialize_tuning(m=700, target=0.75)
a_ca_tuning <- initialize_tuning(m=2000, target=0.75)
a_co_tuning <- initialize_tuning(m=2000, target=0.75)
ca_tuning <- initialize_tuning(m=700, target=0.75)
co_tuning <- initialize_tuning(m=700, target=0.75)


deltas_blyat <- c()
deltas_blyat_ca <- c()
deltas_blyat_co <- c()
deltas_blyat_aca <- c()
deltas_blyat_aco <- c()

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
  w.out.i <- wHmcUpdateCC(Y.l, X.c, Y.ca, alpha.ca.i, beta.ca, Y.co,
                          alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L)
  w.i <- w.out.i$w

  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1/rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  
  ## sample from beta.case
  w.i.sub <- w.i[locs$ids]
  beta.out.ca <- caseHmcUpdate(Y.ca, w.i[locs$ids], X.c, beta.ca, alpha.ca.i, ca_tuning$delta_curr, L_ca)
  beta.ca <- beta.out.ca$beta
  
  ## sample from alpha case
  alpha.out.ca <- alphaHmcUpdate(Y.ca, w.i.sub, X.c, beta.ca, alpha.ca.i, 
                                 a_ca_tuning$delta_curr, prior_alpha_ca_mean, prior_alpha_ca_var, L_a_ca)
  alpha.ca.i <- alpha.out.ca$alpha
  
  ## sample from beta.ctrl
  beta.out.co <- caseHmcUpdate(Y.co, w.i[locs$ids], X.c, beta.co, alpha.co.i, co_tuning$delta_curr, L_co)
  beta.co <- beta.out.co$beta
  
  ## sample from alpha control
  alpha.out.co <- alphaHmcUpdate(Y.co, w.i.sub, X.c, beta.co, alpha.co.i, 
                                 a_co_tuning$delta_curr, prior_alpha_co_mean, prior_alpha_co_var, L_a_co)
  alpha.co.i <- alpha.out.co$alpha
  
  samples.beta.ca[i,] <- beta.ca
  samples.beta.co[i,] <- beta.co
  samples.alpha.ca[i,] <- alpha.ca.i
  samples.alpha.co[i,] <- alpha.co.i
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.w[i,] <- t(w.i)
  
  accept[1] <- accept[1] + w.out.i$accept
  accept[2] <- accept[2] + theta.out$accept
  accept[3] <- accept[3] + beta.out.ca$accept
  accept[4] <- accept[4] + beta.out.co$accept
  accept[5] <- accept[5] + alpha.out.ca$accept
  accept[6] <- accept[6] + alpha.out.co$accept

  w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
  deltas_blyat <- c(deltas_blyat, w_tuning$delta_curr)

  a_ca_tuning <- update_tuning(a_ca_tuning, alpha.out.ca$a, i, alpha.out.ca$accept)
  deltas_blyat_aca <- c(deltas_blyat_aca, a_ca_tuning$delta_curr)
  
  a_co_tuning <- update_tuning(a_co_tuning, alpha.out.co$a, i, alpha.out.co$accept)
  deltas_blyat_aco <- c(deltas_blyat_aco, a_co_tuning$delta_curr)

  ca_tuning <- update_tuning(ca_tuning, beta.out.ca$a, i, beta.out.ca$accept)
  deltas_blyat_ca <- c(deltas_blyat_ca, ca_tuning$delta_curr)
  
  co_tuning <- update_tuning(co_tuning, beta.out.co$a, i, beta.out.co$accept)
  deltas_blyat_co <- c(deltas_blyat_co, co_tuning$delta_curr)
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}


plot(deltas_blyat)
plot(deltas_blyat_ca)
plot(deltas_blyat_co)
plot(deltas_blyat_aca)
plot(deltas_blyat_aco)


accept <- accept/n.sample
print(accept)


plot(apply(samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

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

plot(samples.beta.ca[,1], type='l', ylab='beta0'); abline(h=beta.case[1], col='2')
plot(samples.beta.ca[,2], type='l', ylab='beta1'); abline(h=beta.case[2], col='2')

plot(samples.beta.co[,1], type='l', ylab='beta0'); abline(h=beta.ctrl[1], col='2')
plot(samples.beta.co[,2], type='l', ylab='beta1'); abline(h=beta.ctrl[2], col='2')

par(mfrow=c(1,2))
plot(samples.alpha.ca, type='l', xlab='iteration', ylab='Alpha'); abline(h=Alpha.case, col=2)
hist(samples.alpha.ca, xlab='Alpha'); abline(v=Alpha.case, col=2)
print(mean(samples.alpha.ca))

plot(samples.alpha.co, type='l', xlab='iteration', ylab='Alpha'); abline(h=Alpha.ctrl, col=2)
hist(samples.alpha.co, xlab='Alpha'); abline(v=Alpha.ctrl, col=2)
print(mean(samples.alpha.co))

