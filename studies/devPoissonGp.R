############################
# Fit the case control 
# shared latent process model

# no locational covariates
# or intercept, just w.

# tunes step sizes by dual
# averaging (Hoffman and 
# Gelman, 2014)
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
N.samp <- 70
Theta <- 6
Phi <- 4
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
cells.samp <- sample(cells.all, size=N.samp)
coords <- xyFromCell(caWc.disc, cell=cells.samp)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, N.samp), Sigma)


#### Simulate counts 
Beta <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
X <- cov.disc[][cells.samp]
X.standardised <- matrix((X - mean(X))/sd(X))
X.standardised <- cbind(1, X.standardised)
rates <- exp(X.standardised %*% Beta + W)
Y <- sapply(rates, function(x){rpois(n=1, x)})


Uw_bench <- function(x, y, beta, w, sigma){
  
  logd <- 0
  
  # likelihood: counts
  count_pred <- x %*% beta + w
  rates <- exp(count_pred)
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), rep(0, length(w)), sigma, log=T)
  
  return(-logd)
  
}


Ubeta_bench <- function(y, w, x, beta){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + w
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


dU_w_bench <- function(x, y, beta, w, sigma.inv){
  
  grad <- array(0, c(length(w), 1))
  
  # count contribution
  lin.count <- x %*% beta
  for (i in 1:length(y)){
    grad[i] <- grad[i] + (y[i] - exp(lin.count[i] + w[i]))
  }
  
  # prior contribution
  grad <- grad + t(-t(w) %*% sigma.inv)
  
  return(-grad)
  
}


dU_beta_bench <- function(y, w, x, beta){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - exp(lin_preds[i,]))
    }
  }
  
  if (length(beta) > 1){
    prior_var_inv <- diag(rep(1/100, length(beta)))
  } else{
    prior_var_inv <- 1/100
  }
  
  grad <- grad + t(-t(beta) %*% prior_var_inv)
  return(-grad)
  
}


wHmcUpdateBench <- function(x, y, beta, w, sigma, sigma.inv, delta, L){
  
  # sample random momentum
  p0 <- matrix(rnorm(length(w)))
  
  # simulate Hamiltonian dynamics
  wcurr <- matrix(w)
  pStar <- p0 - 0.5 * delta * dU_w_bench(x, y, beta, wcurr, sigma.inv)
  
  # first full step for position
  wStar <- wcurr + delta*pStar
  
  # full steps
  for (jL in 1:c(L-1)){
    # momentum
    pStar <- pStar - delta * dU_w_bench(x, y, beta, wStar, sigma.inv)
    
    # position
    wStar <- wStar + delta*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta * dU_w_bench(x, y, beta, wStar, sigma.inv)
  
  # evaluate energies
  U0 <- Uw_bench(x, y, beta, wcurr, sigma)
  UStar <- Uw_bench(x, y, beta, wStar, sigma)
  
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


caseHmcUpdateBench <- function(y, w, x, beta, delta_c, L_c){
  
  
  # sample random momentum
  p0 <- matrix(rnorm(length(beta)))
  
  # simulate Hamiltonian dynamics
  bcurr <- matrix(beta)
  pStar <- p0 - 0.5 * delta_c * dU_beta_bench(y, w, x, bcurr)
  
  # first full step for position
  bStar <- bcurr + delta_c*pStar
  
  # full steps
  for (jL in 1:c(L_c-1)){
    # momentum
    pStar <- pStar - delta_c * dU_beta_bench(y, w, x, bStar)
    
    # position
    bStar <- bStar + delta_c*pStar
  }
  
  # last half step
  pStar <- pStar - 0.5 * delta_c * dU_beta_bench(y, w, x, bStar)
  
  # evaluate energies
  U0 <- Ubeta_bench(y, w, x, bcurr)
  UStar <- Ubeta_bench(y, w, x, bStar)
  
  K0 <- K(p0)
  KStar <- K(pStar)
  
  # accept/reject
  alpha <- min(1, exp((U0 + K0) - (UStar + KStar)))
  
  if (is.na(alpha)){
    alpha <- 0
  } 
  
  if (runif(1, 0, 1) < alpha){
    bnext <- bStar
    accept <- 1
  } else {
    bnext <- bcurr
    accept <- 0
  }
  
  
  out <- list()
  out$beta <- bnext
  out$accept <- accept
  out$a <- alpha
  return(out)
  
}


N.w <- ncol(d)
w.i <- c(W + rnorm(length(W)))
beta.i <- Beta + rnorm(length(Beta))
theta.i <- Theta + rnorm(1)
phi.i <- Phi + rnorm(1)
p.c <- length(beta.i)

n.sample <- 25000
accept <- rep(0, 3)

samples.w <- array(NA, c(n.sample, N.w))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))
samples.beta <- array(NA, c(n.sample, p.c))

L <- 20
L_b <- 22
w_tuning <- initialize_tuning(m=700, target=0.75)
b_tuning <- initialize_tuning(m=2000, target=0.75)

deltas_w <- c()
deltas_b <- c()

proposal.sd.theta <- 0.3
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
  w.out.i <- wHmcUpdateBench(X.standardised, Y, beta.i, w.i, sigma.i, sigma.inv.i, w_tuning$delta_curr, L)
  w.i <- w.out.i$w
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  
  ## sample from phi
  R.i <- sigma.i/phi.i
  phi.i <- 1/rgamma(1, N.w/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  
  ## sample from beta.case
  beta.out.i <- caseHmcUpdateBench(Y, w.i, X.standardised, beta.i, b_tuning$delta_curr, L_b)
  beta.i <- beta.out.i$beta
  
  samples.beta[i,] <- beta.i
  samples.theta[i,] <- theta.i
  samples.phi[i,] <- phi.i
  samples.w[i,] <- t(w.i)

  accept[1] <- accept[1] + w.out.i$accept
  accept[2] <- accept[2] + theta.out$accept
  accept[3] <- accept[3] + beta.out.i$accept
  
  w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
  deltas_w <- c(deltas_w, w_tuning$delta_curr)
  
  b_tuning <- update_tuning(b_tuning, beta.out.i$a, i, beta.out.i$accept)
  deltas_b <- c(deltas_b, b_tuning$delta_curr)
  
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

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

plot(samples.phi, type='l'); abline(h=Phi, col='2')
plot(samples.theta, type='l'); abline(h=Theta, col='2')

par(mfrow=c(1,2))
plot(samples.beta[,1], type='l'); abline(h=Beta[1], col='2')
plot(samples.beta[,2], type='l'); abline(h=Beta[2], col='2')


# blyat
# X.l <- cov.disc[][!is.na(cov.disc[])]
# X.l <- matrix((X.l - mean(X.l))/sd(X.l))
# X.l <- cbind(1, X.l)
# true_log_odds <- X.l %*% beta.case - X.l %*% beta.ctrl
# true_log_risk <- true_log_odds/(1 - true_log_odds)
# 
# model_ca <- glm(case.data$y ~ case.data$x.standardised-1, family='poisson')
# model_co <- glm(ctrl.data$y ~ ctrl.data$x.standardised-1, family='poisson')
# 
# est_case <- as.numeric(coefficients(model_ca))
# est_ctrl <- as.numeric(coefficients(model_co))
# est_log_odds <- X.l %*% est_case - X.l %*% est_ctrl
# est_log_risk <- est_log_odds/(1 - est_log_odds)
# 
# plot(x=true_log_odds, y=est_log_odds); abline(0, 1)
# plot(x=true_log_risk, y=est_log_risk); abline(0, 1)