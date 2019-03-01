################################
# Set simulation parameters
# Generate datasets
# zero correlation, medium, high
################################

library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


######################
#### Zero Correlation
######################

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

#### Simulate a 2d multivariate gaussian process 
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Tmat <- matrix(c(8, 0, 0, 9), nrow=2)
Theta <- 6
H <- Exponential(D, range=Theta, phi=1)
Sigma <- kronecker(H, Tmat)
set.seed(41)
W <-  mvrnorm(n=1, mu=rep(0, ncol(Sigma)), Sigma)
W1 <- W[seq(1, ncol(Sigma), by=2)]
W2 <- W[seq(2, ncol(Sigma), by=2)]

par(mfrow=c(1,2))
hist(W1); hist(W2)

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=11)
locs2 <- simLocW(W2, r, beta=0, seed=14)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r, main='A)')
points(locs1$coords, pch=16)
plot(r, main='B)')
points(locs2$coords, pch=16)
locs=list(locs1, locs2)

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
Alpha.case2 <- 1
Alpha.ctrl2 <- -1
beta.case1 <- c(-0.25, 0.75, -0.50)
beta.ctrl1 <- c(3.5, 0.5, 0.5)
beta.case2 <- c(2, 0.8, -0.5)
beta.ctrl2 <- c(2.5, 0.5, 0.5)

case.data1 <- simConditionalGp2(caPr.disc, locs1, beta.case1, Alpha.case1, W1, seed=42)
ctrl.data1 <- simConditionalGp2(caPr.disc, locs1, beta.ctrl1, Alpha.ctrl1, W1, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))
print(case.data1$y); print(ctrl.data1$y)

case.data2 <- simConditionalGp2(caPr.disc, locs2, beta.case2, Alpha.case2, W2, seed=42)
ctrl.data2 <- simConditionalGp2(caPr.disc, locs2, beta.ctrl2, Alpha.ctrl2, W2, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))
print(case.data2$y); print(ctrl.data2$y)

case.data=list(case.data1, case.data2)
ctrl.data=list(ctrl.data1, ctrl.data2)

data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)

save_output(data, "simMVGP_comparison_data_none.json")

params <- list(
  Alpha.cases=list(Alpha.case1, Alpha.case2),
  Alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
  beta.cases=list(beta.case1, beta.case2),
  beta.ctrls=list(beta.ctrl1, beta.ctrl2),
  Theta=Theta,
  Tmat=Tmat,
  W=W
)
save_output(params, "simMVGP_comparison_params_none.json")


######################
#### Medium Correlation
######################

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

#### Simulate a 2d multivariate gaussian process 
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Tmat <- matrix(c(8, 3, 3, 9), nrow=2)
Theta <- 6
H <- Exponential(D, range=Theta, phi=1)
Sigma <- kronecker(H, Tmat)
set.seed(40)
W <-  mvrnorm(n=1, mu=rep(0, ncol(Sigma)), Sigma)
W1 <- W[seq(1, ncol(Sigma), by=2)]
W2 <- W[seq(2, ncol(Sigma), by=2)]

par(mfrow=c(1,2))
hist(W1); hist(W2)

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=11)
locs2 <- simLocW(W2, r, beta=0, seed=14)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r, main='A)')
points(locs1$coords, pch=16)
plot(r, main='B)')
points(locs2$coords, pch=16)
locs=list(locs1, locs2)

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
Alpha.case2 <- 1
Alpha.ctrl2 <- -1
beta.case1 <- c(1, 0.75, -0.50)
beta.ctrl1 <- c(3.5, 0.5, 0.5)
beta.case2 <- c(2.25, 0.8, -0.5)
beta.ctrl2 <- c(2.5, 0.5, 0.5)

case.data1 <- simConditionalGp2(caPr.disc, locs1, beta.case1, Alpha.case1, W1, seed=42)
ctrl.data1 <- simConditionalGp2(caPr.disc, locs1, beta.ctrl1, Alpha.ctrl1, W1, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))
print(case.data1$y); print(ctrl.data1$y)

case.data2 <- simConditionalGp2(caPr.disc, locs2, beta.case2, Alpha.case2, W2, seed=42)
ctrl.data2 <- simConditionalGp2(caPr.disc, locs2, beta.ctrl2, Alpha.ctrl2, W2, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))
print(case.data2$y); print(ctrl.data2$y)

case.data=list(case.data1, case.data2)
ctrl.data=list(ctrl.data1, ctrl.data2)

data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)

save_output(data, "simMVGP_comparison_data_medium.json")

params <- list(
  Alpha.cases=list(Alpha.case1, Alpha.case2),
  Alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
  beta.cases=list(beta.case1, beta.case2),
  beta.ctrls=list(beta.ctrl1, beta.ctrl2),
  Theta=Theta,
  Tmat=Tmat,
  W=W
)
save_output(params, "simMVGP_comparison_params_medium.json")


######################
#### High Correlation
######################

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

#### Simulate a 2d multivariate gaussian process 
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Tmat <- matrix(c(8, 6, 6, 9), nrow=2)
Theta <- 6
H <- Exponential(D, range=Theta, phi=1)
Sigma <- kronecker(H, Tmat)
set.seed(40)
W <-  mvrnorm(n=1, mu=rep(0, ncol(Sigma)), Sigma)
W1 <- W[seq(1, ncol(Sigma), by=2)]
W2 <- W[seq(2, ncol(Sigma), by=2)]

par(mfrow=c(1,2))
hist(W1); hist(W2)

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=11)
locs2 <- simLocW(W2, r, beta=0, seed=14)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r, main='A)')
points(locs1$coords, pch=16)
plot(r, main='B)')
points(locs2$coords, pch=16)
locs=list(locs1, locs2)

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
Alpha.case2 <- 1
Alpha.ctrl2 <- -1
beta.case1 <- c(1.5, 0.75, -0.25)
beta.ctrl1 <- c(2.5, 0.5, 0.5)
beta.case2 <- c(0.15, 0.5, -0.5)
beta.ctrl2 <- c(2.75, 0.5, 0.5)

case.data1 <- simConditionalGp2(caPr.disc, locs1, beta.case1, Alpha.case1, W1, seed=42)
ctrl.data1 <- simConditionalGp2(caPr.disc, locs1, beta.ctrl1, Alpha.ctrl1, W1, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))
print(case.data1$y); print(ctrl.data1$y)

case.data2 <- simConditionalGp2(caPr.disc, locs2, beta.case2, Alpha.case2, W2, seed=42)
ctrl.data2 <- simConditionalGp2(caPr.disc, locs2, beta.ctrl2, Alpha.ctrl2, W2, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))
print(case.data2$y); print(ctrl.data2$y)

case.data=list(case.data1, case.data2)
ctrl.data=list(ctrl.data1, ctrl.data2)

data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)

save_output(data, "simMVGP_comparison_data_high.json")

params <- list(
  Alpha.cases=list(Alpha.case1, Alpha.case2),
  Alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
  beta.cases=list(beta.case1, beta.case2),
  beta.ctrls=list(beta.ctrl1, beta.ctrl2),
  Theta=Theta,
  Tmat=Tmat,
  W=W
)
save_output(params, "simMVGP_comparison_params_high.json")
