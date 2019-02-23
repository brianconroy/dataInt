######################
# develop a separable 
# Multivariate Gaussian 
# Process Model
######################

library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=9)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate a 2d multivariate gaussian process 
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Tmat <- matrix(c(8, 2, 2, 6), nrow=2)
Theta <- 6
H <- Exponential(d, range=Theta, phi=1)
Sigma <- kronecker(H, Tmat)
set.seed(42)
W <-  mvrnorm(n=1, mu=rep(0, ncol(Sigma)), Sigma)
W1 <- W[seq(1, ncol(Sigma), by=2)]
W2 <- W[seq(2, ncol(Sigma), by=2)]

# ultimately...
# update W by HMC 
# update Theta by MHRW
# update Tmat by Gibbs sampling

# to start...
# fix W and Theta at true values
# update Tmat by gibbs sampling

## storage
n <- length(W1)
n.sample <- 500
burnin <- 0
n.keep <- n.sample - burnin
samples.t <- array(NA, c(n.keep, ncell(Tmat)))

## initial values 
T.i <- matrix(c(3, 0, 0, 3), nrow=2)
H.i <- H
H.inv.i <- solve(H.i)

## priors
Omega <- matrix(c(5, 0, 0, 5), nrow=2)
r <- 4
print(Omega/(r - ncol(Tmat) - 1))

cat("Generating", n.sample, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*n.sample)

for (i in 1:n.sample){
  
  ## sample from T
  r_ <- r + n
  Omega_ <- Omega
  for (a in 1:n){
    for (b in 1:n){
      Omega_ <- Omega_ + H.inv.i[a, b] * matrix(c(W1[b], W2[b])) %*% t(matrix(c(W1[a], W2[a])))
    }
  }
  T.i <- riwish(r_, Omega_)
  
  ## store samples
  samples.t[i,] <- T.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
}

matrix(colMeans(samples.t), nrow=2)
par(mfrow=c(2,2))
plot(samples.t[,1], type='l'); abline(h=Tmat[1,1], col=2)
plot(samples.t[,2], type='l'); abline(h=Tmat[2,1], col=2)
plot(samples.t[,3], type='l'); abline(h=Tmat[2,1], col=2)
plot(samples.t[,4], type='l'); abline(h=Tmat[2,2], col=2)


# next step...
# update range by MHRW
# fix T at its true value


rangeMVGPupdate <- function(h.i, t.i, w.i, d, theta.i, proposal.sd, prior){
  
  # proposal
  theta.new <- rlnorm(1, meanlog=log(theta.i), sdlog=proposal.sd)
  q.new <- dlnorm(theta.new, meanlog=log(theta.i), sdlog=proposal.sd, log=T)
  q.old <- dlnorm(theta.i, meanlog=log(theta.new), sdlog=proposal.sd, log=T)
  
  # covariance matrices
  h.new <- Exponential(d, range=theta.new, phi=1)

  # likelihoods
  loglik.curr <- dmvnorm(w.i, sigma=kronecker(h.i, t.i), log=T)
  loglik.new <- dmvnorm(w.i, sigma=kronecker(h.new, t.i), log=T)
  like.diff <- loglik.new - loglik.curr
  
  # priors
  shape <- prior[1]
  scale <- prior[2]
  prior.curr <- dgamma(theta.i, shape=shape, scale=scale, log=T)
  prior.new <- dgamma(theta.new, shape=shape, scale=scale, log=T)
  prior.diff <- prior.new - prior.curr
  
  out <- list()
  acceptance <- exp(like.diff + prior.diff + q.old - q.new)
  if(runif(1) <= acceptance) {
    out$theta <- theta.new
    out$accept <- 1
  } else { 
    out$theta <- theta.i
    out$accept <- 0
  }
  
  return(out)
  
}

## storage
n <- length(W1)
n.sample <- 200
burnin <- 0
n.keep <- n.sample - burnin
accept <- list(theta=0)
samples.theta <- array(NA, c(n.keep, 1))

## initial values 
T.i <- Tmat
theta.i <- 3
H.i <- Exponential(d, range=theta.i, phi=1)
H.inv.i <- solve(H.i)
w.i <- W

proposal.sd.theta <- 0.15

## priors
Omega <- matrix(c(5, 0, 0, 5), nrow=2)
r <- 4
print(Omega/(r - ncol(Tmat) - 1))
prior.theta <- c(3, 2)

cat("Generating", n.sample, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*n.sample)

for (i in 1:n.sample){
  
  ## sample from theta
  theta.out <- rangeMVGPupdate(H.i, T.i, w.i, d, theta.i, proposal.sd.theta, prior.theta)
  theta.i <- theta.out$theta
  # r_ <- r + n
  # Omega_ <- Omega
  # for (a in 1:n){
  #   for (b in 1:n){
  #     Omega_ <- Omega_ + H.inv.i[a, b] * matrix(c(W1[b], W2[b])) %*% t(matrix(c(W1[a], W2[a])))
  #   }
  # }
  # T.i <- riwish(r_, Omega_)
  
  ## store samples
  if (i > burnin){
    
    j <- i - burnin
    samples.theta[j,] <- theta.i
    accept$theta <- accept$theta + theta.out$accept
    # samples.t[i,] <- T.i
    
  }
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
}
accept$theta <- accept$theta/n.keep

print(accept)
plot(samples.theta, type='l'); abline(h=Theta, col=2)


# next step...
# update range and T matrix
## storage
n <- length(W1)
n.sample <- 200
burnin <- 0
n.keep <- n.sample - burnin
accept <- list(theta=0)
samples.theta <- array(NA, c(n.keep, 1))
samples.t <- array(NA, c(n.keep, ncell(Tmat)))

## initial values 
T.i <- matrix(c(3, 0, 0, 3), nrow=2)
theta.i <- 3
H.i <- Exponential(d, range=theta.i, phi=1)
H.inv.i <- solve(H.i)
w.i <- W

proposal.sd.theta <- 0.2

## priors
Omega <- matrix(c(5, 0, 0, 5), nrow=2)
r <- 4
print(Omega/(r - ncol(Tmat) - 1))
prior.theta <- c(3, 2)

cat("Generating", n.sample, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*n.sample)

for (i in 1:n.sample){
  
  ## sample from Theta
  theta.out <- rangeMVGPupdate(H.i, T.i, w.i, d, theta.i, proposal.sd.theta, prior.theta)
  theta.i <- theta.out$theta
  
  ## sample from T
  H.i <- Exponential(d, range=theta.i, phi=1)
  H.inv.i <- solve(H.i)
  r_ <- r + n
  Omega_ <- Omega
  for (a in 1:n){
    for (b in 1:n){
      Omega_ <- Omega_ + H.inv.i[a, b] * matrix(c(W1[b], W2[b])) %*% t(matrix(c(W1[a], W2[a])))
    }
  }
  T.i <- riwish(r_, Omega_)
  
  ## store samples
  if (i > burnin){
    
    j <- i - burnin
    samples.theta[j,] <- theta.i
    samples.t[j,] <- T.i
    accept$theta <- accept$theta + theta.out$accept
    
  }
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
}
accept$theta <- accept$theta/n.keep

print(accept)
plot(samples.theta, type='l'); abline(h=Theta, col=2)

print(matrix(colMeans(samples.t), nrow=2))
par(mfrow=c(2,2))
plot(samples.t[,1], type='l'); abline(h=Tmat[1,1], col=2)
plot(samples.t[,2], type='l'); abline(h=Tmat[2,1], col=2)
plot(samples.t[,3], type='l'); abline(h=Tmat[2,1], col=2)
plot(samples.t[,4], type='l'); abline(h=Tmat[2,2], col=2)

# next step...
# update W from poisson counts
# with theta and phi fixed
# ?
# merge Y1 and Y2 into Y
# merge W1 and W2 into W
# update W with your old sampler
# split it back to W1 and W2
# do all your other shit

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=11)
locs2 <- simLocW(W2, r, beta=0, seed=42)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r)
points(locs1$coords)
plot(r)
points(locs2$coords)
locs=list(locs1, locs2)

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
Alpha.case2 <- 1
Alpha.ctrl2 <- -1
beta.case1 <- c(0.25, 0.75, -0.50)
beta.ctrl1 <- c(2.75, 0.5, 0.5)
beta.case2 <- c(2, 0.8, -0.5)
beta.ctrl2 <- c(2, 0.5, 0.5)

case.data1 <- simConditionalGp2(caPr.disc, locs1, beta.case1, Alpha.case1, W1, seed=42)
ctrl.data1 <- simConditionalGp2(caPr.disc, locs1, beta.ctrl1, Alpha.ctrl1, W1, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))

case.data2 <- simConditionalGp2(caPr.disc, locs2, beta.case2, Alpha.case2, W2, seed=42)
ctrl.data2 <- simConditionalGp2(caPr.disc, locs2, beta.ctrl2, Alpha.ctrl2, W2, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))

case.data=list(case.data1, case.data2)
ctrl.data=list(ctrl.data1, ctrl.data2)

# alpha.ca <- list(Alpha.case1, Alpha.case2)
# beta.ca <- list(beta.case1, beta.case2)
# alpha.co <- list(Alpha.ctrl1, Alpha.ctrl2)
# beta.co <- list(beta.ctrl1, beta.ctrl2)
# w <- W
# sigma <- Sigma
# sigma.inv <- solve(sigma)
L_w <- 8
L_ca <- c(8, 8)
L_co <- c(8, 8)
L_a_ca <- c(8, 8)
L_a_co <- c(8, 8)
proposal.sd.theta <- 0.15

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

target_aca <- 0.65
target_aco <- 0.65
target_ca <- 0.65
target_co <- 0.65
target_w <- 0.65

n <- length(W1)
n.sample <- 2000
burnin <- 0
n.keep <- n.sample - burnin
accept <- list(w=0)
samples.w <- array(NA, c(n.keep, ncol(Sigma)))

## initial values 
T.i <- Tmat
theta.i <- Theta
H.i <- Exponential(d, range=theta.i, phi=1)
H.inv.i <- solve(H.i)
w.i <- rep(0, length(W))
alpha.ca.i <- list(Alpha.case1, Alpha.case2)
alpha.co.i <- list(Alpha.ctrl1, Alpha.ctrl2)
beta.ca <- list(beta.case1, beta.case2)
beta.co <- list(beta.ctrl1, beta.ctrl2)

self_tune_w <- T
deltas_w <- c()
w_tuning <- initialize_tuning(m=500, target=0.65)

cat("Generating", n.sample, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*n.sample)

sigma.i <- kronecker(H.i, T.i)
sigma.inv.i <- kronecker(solve(H.i), solve(T.i))

for (i in 1:n.sample){
  
  ## sample from w
  w.out.i <- wHmcUpdateMVGP(case.data, ctrl.data, alpha.ca.i, beta.ca,
                            alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w)
  w.i <- w.out.i$w
  
  samples.w[i,] <- t(w.i)
  
  accept$w <- accept$w + w.out.i$accept
  
  if (self_tune_w){
    w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
    deltas_w <- c(deltas_w, w_tuning$delta_curr)
  }
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
}

accept$w <- accept$w/n.keep
print(accept)
plot(deltas_w)
w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
view_tr_w(samples.w, w_true=W)
plot(apply(samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
