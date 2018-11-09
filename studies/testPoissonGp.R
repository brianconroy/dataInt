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
N.samp <- 100
Theta <- 6
Phi <- 10
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
cells.samp <- sample(cells.all, size=N.samp)
samp.inds <- cells.all %in% cells.samp
samp.ids <- c(1:length(cells.all))[samp.inds]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d_full <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
d <- d_full[samp.inds, samp.inds]
Sigma <- Exponential(d_full, range=Theta, phi=Phi)
W_full <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
W <- W_full[samp.inds]


#### Simulate counts 
Beta <- c(2, 1)
cov.disc <- caWc.disc[[c(1)]]
X <- cov.disc[][cells.samp]
X.standardised <- matrix((X - mean(X))/sd(X))
X.standardised <- cbind(1, X.standardised)
rates <- exp(X.standardised %*% Beta + W)
Y <- sapply(rates, function(x){rpois(n=1, x)})


#### Fit model
output <- poissonGp(X.standardised, Y, d, n.sample=10000, burnin=0, L_w=30, L_b=30, proposal.sd.theta=0.3,
                    beta_initial=Beta+rnorm(2), w_initial=W+rnorm(length(W)), 
                    phi_initial=Phi+rnorm(1), theta_initial=Theta+rnorm(1),
                    prior_phi=c(3, 40))

print(output$accept)

plot(apply(output$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(output$samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

par(mfrow=c(1,2))
plot(output$samples.phi, type='l'); abline(h=Phi, col='2')
plot(output$samples.theta, type='l'); abline(h=Theta, col='2')
print(mean(output$samples.phi))
print(mean(output$samples.theta))


plot(output$samples.beta[,1], type='l'); abline(h=Beta[1], col='2')
plot(output$samples.beta[,2], type='l'); abline(h=Beta[2], col='2')


#### test kriging
kriged_w <- krigeW(output, d_full, samp.ids)
true_W <- W_full[!samp.inds]
plot(x=true_W, y=kriged_w$mu.new, xlab='true W', ylab='kriged W', main='random sampling'); abline(0, 1, col='2')


##############################
# Preliminary study
# Preferential sampling effect
##############################


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
Phi <- 20
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d_full <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d_full, range=Theta, phi=Phi)
W_full <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)

N.samp <- 100
samp.inds <- W_full %in% sort(W_full, decreasing=T)[1:N.samp]
cells.samp <- cells.all[samp.inds]
samp.ids <- c(1:length(cells.all))[samp.inds]
d <- d_full[samp.inds, samp.inds]
W <- W_full[samp.inds]


#### Simulate counts 
Beta <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
X <- cov.disc[][cells.samp]
X.standardised <- matrix((X - mean(X))/sd(X))
X.standardised <- cbind(1, X.standardised)
rates <- exp(X.standardised %*% Beta + W)
Y <- sapply(rates, function(x){rpois(n=1, x)})


#### Fit model
output <- poissonGp(X.standardised, Y, d, n.sample=5000, burnin=0, L_w=20, L_b=22, proposal.sd.theta=0.3,
                    beta_initial=Beta+rnorm(2), w_initial=W+rnorm(length(W)), 
                    phi_initial=Phi+rnorm(1), theta_initial=Theta+rnorm(1),
                    prior_phi=c(3, 40))


print(output$accept)

plot(apply(output$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(output$samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}

par(mfrow=c(1,2))
plot(output$samples.phi, type='l'); abline(h=Phi, col='2')
plot(output$samples.theta, type='l'); abline(h=Theta, col='2')

plot(output$samples.beta[,1], type='l'); abline(h=Beta[1], col='2')
plot(output$samples.beta[,2], type='l'); abline(h=Beta[2], col='2')


#### test kriging
kriged_w <- krigeW(output, d_full, samp.ids)
true_W <- W_full[!samp.inds]
plot(x=true_W, y=kriged_w$mu.new, xlab='true W', ylab='kriged W', main='preferential sampling'); abline(0, 1, col='2')
