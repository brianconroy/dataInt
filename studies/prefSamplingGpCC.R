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
## Case Counts
Alpha.case <- 2
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W)
glm(case.data$y ~ case.data$x.standardised-1, family='poisson')
prior_alpha_ca_mean <- 2
prior_alpha_ca_var <- 6

## Control Counts
Alpha.ctrl <- -2
beta.ctrl <- c(2, 0.5)
cov.disc <- caWc.disc[[c(1)]]
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W)
glm(ctrl.data$y ~ ctrl.data$x.standardised-1, family='poisson')
prior_alpha_co_mean <- -2
prior_alpha_co_var <- 6

data <- list(
  loc=locs,
  case.data=case.data,
  ctrl.data=ctrl.data
)


## starting values
# w.i <- W
# alpha.ca.i <- Alpha.case
# beta.ca <- beta.case
# alpha.co.i <- Alpha.ctrl
# beta.co <- beta.ctrl
# theta.i <- Theta
# phi.i <- Phi
# p.c <- length(beta.ca)

# w.i <- W + rnorm(length(W))
# alpha.ca.i <- Alpha.case + rnorm(1)
# beta.ca <- beta.case + rnorm(length(beta.ca))
# alpha.co.i <- Alpha.ctrl + rnorm(1)
# beta.co <- beta.ctrl + rnorm(length(beta.co))
# theta.i <- Theta + rnorm(1)
# phi.i <- Phi + rnorm(1)
# p.c <- length(beta.ca)

# w.i <- rnorm(length(W))
# alpha.ca.i <- runif(1, 1, 3)
# beta.ca <- rnorm(2, 1)
# alpha.co.i <- runif(1, -3, -1)
# beta.co <- rnorm(2, 1)
# theta.i <- runif(1, 0, 4)
# phi.i <- runif(1, 0, 4)
# p.c <- length(beta.ca)


n.sample <- 4000
burnin <- 0
L <- 20
L_ca <- 25
L_co <- 25
L_a_ca <- 22
L_a_co <- 22


output <- prefSampleGpCC(data, n.sample, burnin, 
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=0.3,
                         m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, 
                         target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, 
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                         delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                         beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL)


plot(output$deltas_w)
plot(output$deltas_ca)
plot(output$deltas_co)
plot(output$deltas_aca)
plot(output$deltas_aco)

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
plot(output$samples.theta, type='l', ylab='Theta'); abline(h=Theta, col='2')
hist(output$samples.theta, xlab='Theta', main=''); abline(v=Theta, col='2')
print(mean(output$samples.theta)); print(Theta)

plot(output$samples.phi, type='l', ylab='Phi'); abline(h=Phi, col='2')
hist(output$samples.phi, xlab='Phi', main=''); abline(v=Phi, col='2')
print(mean(output$samples.phi)); print(Phi)

plot(output$samples.beta.ca[,1], type='l', ylab='beta0'); abline(h=beta.case[1], col='2')
plot(output$samples.beta.ca[,2], type='l', ylab='beta1'); abline(h=beta.case[2], col='2')

plot(output$samples.beta.co[,1], type='l', ylab='beta0'); abline(h=beta.ctrl[1], col='2')
plot(output$samples.beta.co[,2], type='l', ylab='beta1'); abline(h=beta.ctrl[2], col='2')

par(mfrow=c(1,2))
plot(output$samples.alpha.ca, type='l', xlab='iteration', ylab='Alpha'); abline(h=Alpha.case, col=2)
hist(output$samples.alpha.ca, xlab='Alpha'); abline(v=Alpha.case, col=2)
print(mean(output$samples.alpha.ca))

plot(output$samples.alpha.co, type='l', xlab='iteration', ylab='Alpha'); abline(h=Alpha.ctrl, col=2)
hist(output$samples.alpha.co, xlab='Alpha'); abline(v=Alpha.ctrl, col=2)
print(mean(output$samples.alpha.co))
