############################
# Fit the shared latent 
# process model

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
locs <- simLocW(W, beta=0, seed=42)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)


#### Simulate counts given locations
Alpha <- 2
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
count.data <- simConditionalGp2(cov.disc, locs, beta.case, beta.samp, Alpha, W)
glm(count.data$y ~ count.data$x.standardised-1, family='poisson')
prior_alpha_mean <- 2
prior_alpha_var <- 6


#### Combine data
data <- list(loc=locs, conditional=count.data)


#### Run model
output <- prefSampleGp(data, n.sample=11000, burnin=2000, L_w=20, L_c=20, L_a=14, proposal.sd.theta=0.3, target_a=0.95)
output <- prefSampleGp(data, n.sample=15000, burnin=4000, L_w=20, L_c=20, L_a=14, proposal.sd.theta=0.3, target_a=0.95,
                       target_w=0.9)

plot(output$deltas_w)
plot(output$deltas_a)
plot(output$deltas_c)
print(output$accept)


## Inspect w
plot(apply(output$samples.w, 1, mean), type='l', col='2', ylab='mean w'); abline(h=mean(W), col='2')
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
glm(Y.c ~ X.c + w.hat[locs$ids] -1, family='poisson')

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(output$samples.w[,j*i], type='l'); abline(h=W[j*i], col='2')
}


## Inspect other parameters
par(mfrow=c(1,2))
plot(output$samples.theta, type='l', ylab='Theta'); abline(h=Theta, col='2')
hist(output$samples.theta, xlab='Theta', main=''); abline(v=Theta, col='2')
print(mean(output$samples.theta)); print(Theta)

plot(output$samples.phi, type='l', ylab='Phi'); abline(h=Phi, col='2')
hist(output$samples.phi, xlab='Phi', main=''); abline(v=Phi, col='2')
print(mean(output$samples.phi)); print(Phi)

plot(output$samples.beta.c[,1], type='l', ylab='beta0'); abline(h=beta.case[1], col='2')
plot(output$samples.beta.c[,2], type='l', ylab='beta1'); abline(h=beta.case[2], col='2')

hist(output$samples.beta.c[,1]); abline(v=beta.case[1], col='2')
hist(output$samples.beta.c[,2]); abline(v=beta.case[2], col='2')
print(colMeans(output$samples.beta.c))

par(mfrow=c(1,2))
plot(output$samples.alpha, type='l', xlab='iteration', ylab='Alpha'); abline(h=Alpha, col=2)
hist(output$samples.alpha, xlab='Alpha'); abline(v=Alpha, col=2)
print(mean(output$samples.alpha))
