library(MCMCpack)
library(R.utils)
library(fields)
library(mvrnorm)
library(mvtnorm)
sourceDirectory('Documents/research/dataInt/R/')


####################################
# Test an MCMC algorithm for Poisson 
# counts with Gaussian Random Field
####################################


wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)
cov.disc <- caWc[[1]]

# observation sites
N <- 50
cells <- sample(c(1:length(cov.disc[]))[!is.na(cov.disc[])], N)
coords <- xyFromCell(cov.disc, cell=cells)
plot(cov.disc)
points(coords)

# distance matrix
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

# covariance matrix
Theta <- 5
Phi <- 3
Sigma <- Exponential(D, range=Theta, phi=Phi)
plot(as.im(Sigma))

# spatial process
W <- mvrnorm(n=1, mu=rep(0, N), Sigma)
hist(W)

# poisson counts
Mu <- 6
rates <- exp(Mu + W)
Y <- sapply(rates, function(x){rpois(n=1, x)})

curve(dgamma(x, 0.5, 0.005))



# fit model
output <- gpMCMC(Y, D, prior.theta=c(2.2, 2.2), prior.phi=c(3, 6),
                 proposal.sd.w=0.1)

print(mean(output$samples.mu))
plot(output$samples.mu, type='l')
abline(h=Mu, col=2)

print(mean(output$samples.theta))
plot(output$samples.theta, type='l')
abline(h=Theta, col='2')

what <- colMeans(output$samples.w)
plot(W, what)
abline(0, 1, col=2)
hist(what, xlim=c(-3, 3))
hist(W, col=2, add=T)

par(mfrow=c(4,4))
j <- 1
for (i in 1:16){
  plot(output$samples.w[,j*i], type='l')
  abline(h=W[j*i], col='2')
}

print(mean(output$samples.phi))
plot(output$samples.phi, type='l')
abline(h=Phi, col='2')

print(output$accept)


################
# test mu update
################


n_iter <- 50000
accept <- c(0, 0, 0)
samples.mu <- array(NA, n_iter)
proposal.sd.mu <- 1
mu.i <- runif(1, 2, 5)

for (i in 1:n_iter){
  
  mu.out <- muPoissonUpdate(Y, mu.i, W, proposal.sd.mu)
  mu.i <- mu.out$mu
  accept[3] <- accept[3] + mu.out$accept
  samples.mu[i] <- mu.i
  
}


print(mean(mu.i))
plot(samples.mu, type='l')
abline(h=Mu, col=2)



##################
# test w MH update
##################


n_iter <- 50000
accept <- c(0, 0)
samples.w <- array(NA, c(n_iter, N))
proposal.sd.w <- 0.5
w.i <- mvrnorm(n=1, mu=rep(0, N), diag(rep(1, N)))

for (i in 1:n_iter){
  
  # update w(s)
  w.out <- wPoissonUpdate(Y, w.i, Mu, Sigma, proposal.sd.w)
  w.i <- w.out$w
  w.i <- w.i - mean(w.i)
  accept[1] <- accept[1] + w.out$accept
  samples.w[i,] <- w.i
  
  #k <- i/1000
  #if(ceiling(k)==floor(k)){
  #  proposal.sd.w <- tuners.sd.1(accept[1], i*N, proposal.sd.w, 40, 50)
  #}
  
}

what <- colMeans(samples.w)
print(what)
hist(what, xlim=c(-5, 5))
hist(W, col=rgb(0,1,0.2), add=T)
hist(100*(W - what)/W)

plot(W, what)
abline(0,1, col='2')
abline(h=0)

par(mfrow=c(1,1))
plot(x=W, y=abs(100*(W - what)/W), ylim=c(0, 500))

par(mfrow=c(4,4))
j <- 2
for (i in 1:16){
  plot(samples.w[,j*i], type='l')
  abline(h=W[j*i], col='2')
}


###################
# test theta update
###################


n_iter <- 50000
accept <- c(0, 0)
samples.theta <- array(NA, c(n_iter, 1))
proposal.sd.theta <- 0.5
theta.i <- runif(1, 0, 4)
prior.theta <- c(4, 2)


for (i in 1:n_iter){
  
  # update range
  theta.out <- rangeMhUpdate(theta.i, W, D, Phi, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  accept[2] <- accept[2] + theta.out$accept
  samples.theta[i,] <- theta.i
  
}


print(colMeans(samples.theta))
plot(samples.theta, type='l')
abline(h=Theta, col='2')


#########################
# test theta and w update
#########################


n_iter <- 50000
accept <- c(0, 0)

samples.theta <- array(NA, c(n_iter, 1))
proposal.sd.theta <- 0.5
theta.i <- Theta
prior.theta <- c(4, 2)

samples.w <- array(NA, c(n_iter, N))
proposal.sd.w <- 0.1
w.i <- mvrnorm(n=1, mu=rep(0, N), diag(rep(1, N)))

for (i in 1:n_iter){
  
  # update w(s)
  w.out <- wPoissonUpdate(Y, w.i, Mu, Sigma, proposal.sd.w)
  w.i <- w.out$w
  w.i <- w.i - mean(w.i)
  accept[1] <- accept[1] + w.out$accept
  samples.w[i,] <- w.i
  
  # update range
  theta.out <- rangeMhUpdate(theta.i, w.i, D, Phi, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  accept[2] <- accept[2] + theta.out$accept
  samples.theta[i,] <- theta.i
  
  k <- i/100
  if(ceiling(k)==floor(k)){
    proposal.sd.w <- tuners.sd.1(accept[1], i*N, proposal.sd.w, 75, 80)
  }
  
}


print(mean(samples.theta))
plot(samples.theta, type='l')
abline(h=Theta, col='2')

what <- colMeans(samples.w)
plot(W, what)
abline(0, 1, col=2)
hist(what, xlim=c(-3, 3))
hist(W, col=2, add=T)

par(mfrow=c(4,4))
j <- 1
for (i in 1:16){
  plot(samples.w[,j*i], type='l')
  abline(h=W[j*i], col='2')
}


############
# update all
############


n_iter <- 50000
accept <- c(0, 0, 0)

samples.mu <- array(NA, n_iter)
proposal.sd.mu <- 0.1
mu.i <- runif(1, 2, 5)

samples.theta <- array(NA, c(n_iter, 1))
proposal.sd.theta <- 0.5
theta.i <- runif(1, 1, 5)
prior.theta <- c(4, 2)

samples.phi <- array(NA, c(n_iter, 1))
phi.i <- runif(1, 1, 5)
prior.phi <- c(2, 4)

samples.w <- array(NA, c(n_iter, N))
proposal.sd.w <- 0.1
w.i <- mvrnorm(n=1, mu=rep(0, N), diag(rep(1, N)))

for (i in 1:n_iter){
  
  # update mu
  mu.out <- muPoissonUpdate(Y, mu.i, w.i, proposal.sd.mu)
  mu.i <- mu.out$mu
  accept[3] <- accept[3] + mu.out$accept
  samples.mu[i] <- mu.i
  
  # update w(s)
  sigma.i <- Exponential(D, range=theta.i, phi=phi.i)
  w.out <- wPoissonUpdate(Y, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w.out$w
  w.i <- w.i - mean(w.i)
  accept[1] <- accept[1] + w.out$accept
  samples.w[i,] <- w.i
  
  # update range
  theta.out <- rangeMhUpdate(theta.i, w.i, D, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.i <- theta.out$theta
  accept[2] <- accept[2] + theta.out$accept
  samples.theta[i,] <- theta.i
  
  # update tau2
  R.i <- sigma.i/phi.i
  phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  samples.phi[i,] <- phi.i
  
  k <- i/100
  if(ceiling(k)==floor(k)){
    proposal.sd.w <- tuners.sd.1(accept[1], i*N, proposal.sd.w, 75, 80)
  }
  
}

print(mean(mu.i))
plot(samples.mu, type='l')
abline(h=Mu, col=2)

print(mean(samples.theta))
plot(samples.theta, type='l')
abline(h=Theta, col='2')

what <- colMeans(samples.w)
plot(W, what)
abline(0, 1, col=2)
hist(what, xlim=c(-3, 3))
hist(W, col=2, add=T)

par(mfrow=c(4,4))
j <- 1
for (i in 1:16){
  plot(samples.w[,j*i], type='l')
  abline(h=W[j*i], col='2')
}

print(mean(samples.phi))
plot(samples.phi, type='l')
abline(h=Phi, col='2')

accept <- accept/n_iter
accept[1] <- accept[1]/N
print(accept)


################
# Iterate over N
################


pbias <- function(true, est){
  return(round(100*(true - est)/true))
}


summaryGp <- data.frame()
for (n in c(10, 15, 20, 35, 50, 70, 100)){
  print(paste(n, "observation sites"))
  
  # observation sites
  cells <- sample(c(1:length(cov.disc[]))[!is.na(cov.disc[])], n)
  coords <- xyFromCell(cov.disc, cell=cells)
  plot(cov.disc)
  points(coords)
  
  # distance matrix
  D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  
  # covariance matrix
  Theta <- 4
  Phi <- 3
  Sigma <- Exponential(D, range=Theta, phi=Phi)

  # spatial process
  W <- mvrnorm(n=1, mu=rep(0, n), Sigma)
  hist(W)
  
  # poisson counts
  Mu <- 5
  rates <- exp(Mu + W)
  Y <- sapply(rates, function(x){rpois(n=1, x)})
  
  # fit model
  output <- gpMCMC(Y, D)
  summ.n <- summarizeGp(output)
  summ.n$n <- n
  summaryGp <- rbind(summaryGp, summ.n)
  
}


par(mfrow=c(2,2))
plot(x=summaryGp$n, y=summaryGp$mu.pbias, type='l', col='2')
abline(h=0)

plot(x=summaryGp$n, y=summaryGp$theta.pbias, type='l', col='2')
abline(h=0)

plot(x=summaryGp$n, y=summaryGp$phi.pbias, type='l', col='2')
abline(h=0)

plot(x=summaryGp$n, y=summaryGp$w.mpbias, type='l', col='2')
abline(h=0)


#########
# Kriging
#########


## test conditional normal formula


# all cells
N <- 100
cells <- sample(c(1:length(cov.disc[]))[!is.na(cov.disc[])], N)
coords <- xyFromCell(cov.disc, cell=cells)

# distance matrix
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

# covariance matrix
Theta <- 2
Phi <- 3
Sigma <- Exponential(D, range=Theta, phi=Phi)
plot(as.im(Sigma))

# spatial process
W <- mvrnorm(n=1, mu=rep(0, N), Sigma)
hist(W)

# sample from W
nsamp <- 25
ids <- sample(1:N, nsamp)
W.samp <- W[ids]

coords.samp <- coords[ids,]
plot(cov.disc)
points(coords)
points(coords.samp, col=2, pch=16)

# rearrange cov matrix to block format
omega.11 <- Sigma[-ids, -ids]
omega.12 <- Sigma[-ids, ids]
omega.21 <- Sigma[ids, -ids]
omega.22 <- Sigma[ids, ids]

E.cond <- omega.12 %*% solve(omega.22) %*% W.samp
Var.cond <- omega.11 - omega.12 %*% solve(omega.22) %*% omega.21

W.pred <- mvrnorm(n=1, mu=E.cond, Var.cond)
hist(W.pred - W[-ids])
hist((W[-ids] - W.pred)/W[-ids])
summary(W[-ids] - W.pred)
summary((W[-ids] - W.pred)/W[-ids])


## full run


wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)
cov.disc <- caWc[[1]]

# all sites
N <- 100
cells <- sample(c(1:length(cov.disc[]))[!is.na(cov.disc[])], N)
coords <- xyFromCell(cov.disc, cell=cells)
plot(cov.disc)
points(coords)

# observed sites
N.obs <- 50
ids <- sample(1:N, N.obs)

# distance matrix
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

# covariance matrix
Theta <- 2
Phi <- 3
Sigma <- Exponential(D, range=Theta, phi=Phi)
plot(as.im(Sigma))

# spatial process
W <- mvrnorm(n=1, mu=rep(0, N), Sigma)
hist(W)

# poisson counts
Mu <- 5
rates <- exp(Mu + W)
Y <- sapply(rates, function(x){rpois(n=1, x)})

# quantities for observed sites
W.obs <- W[ids]
Sigma.obs <- Sigma[ids, ids]
Y.obs <- Y[ids]
D.obs <- D[ids, ids]

# fit model
output <- gpMCMC(Y.obs, D.obs, n_iter=100000)


krg <- krigeW(output, D, ids)
hist((W[-ids] - krg$mu.new)/W[-ids])
summary((W[-ids] - krg$mu.new)/W[-ids])
plot(W[-ids], krg$mu.new)
abline(0, 1)
par(mfrow=c(4,4))
for (i in 1:16){
  plot(krg$samples.new[,i], type='l')
  abline(h=W[-ids][i], col=2)
}
