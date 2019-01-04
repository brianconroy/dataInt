############################
# Multispecies simulation:
# multiple preferentially sampled
# rodent species. 
#   2 to begin
############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Simulation parameters
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
Alpha.case2 <- 1
Alpha.ctrl2 <- -1
beta.case1 <- c(0.0, 0.75, -0.50)
beta.ctrl1 <- c(2.75, 1.00, 0.50)
beta.case2 <- c(-0.25, 0.8, -1)
beta.ctrl2 <- c(1, 1.25, 1)


Theta <- 6
Phi <- 12


prior_theta <- c(6, 1)
prior_phi <- c(18, 204)
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4
prior_alpha_ca_mean1 <- Alpha.case1
prior_alpha_co_mean1 <- Alpha.ctrl1
prior_alpha_ca_mean2 <- Alpha.case2
prior_alpha_co_mean2 <- Alpha.ctrl2


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate gaussian process
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)
hist(W)


#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W, r, beta=0, seed=11) # 42
locs2 <- simLocW(W, r, beta=0, seed=42)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r)
points(locs1$coords)
plot(r)
points(locs2$coords)


#### Simulate counts given locations
cov.disc <- caPr.disc
case.data1 <- simConditionalGp2(cov.disc, locs1, beta.case1, Alpha.case1, W, seed=42)
ctrl.data1 <- simConditionalGp2(cov.disc, locs1, beta.ctrl1, Alpha.ctrl1, W, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))

case.data2 <- simConditionalGp2(cov.disc, locs2, beta.case2, Alpha.case2, W, seed=42)
ctrl.data2 <- simConditionalGp2(cov.disc, locs2, beta.ctrl2, Alpha.ctrl2, W, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))


data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)




