################################
# Simulate data for the study
# testing the separability 
# assumption. 

# Coregionalization model
#   exponential covariance functions
#   mean zero
#   marginal variance 1
# 3 datasets at varying differences
# of spatial range
# 2 responses at each site
################################

library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


#########################
#### low range difference
#########################

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

#### Simulate multivariate spatial process from
#### a linear coregionalization model
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Theta1 <- 2
Theta2 <- 4
Sigma1 <- Exponential(D, range=Theta1, phi=3)
Sigma2 <- Exponential(D, range=Theta2, phi=3)
set.seed(40)
w1 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma1)
w2 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma2)
A <- matrix(c(0.25, 0.25, -1, -0.5), nrow=2)

W <- c()
for (s in 1:length(w1)){
  w_s <- A %*% matrix(c(w1[s], w2[s]))
  W <- c(W, w_s)
}

W1 <- W[seq(1, length(W), by=2)]
W2 <- W[seq(2, length(W), by=2)]
hist(W1)
hist(W2)

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=33)
locs2 <- simLocW(W2, r, beta=0, seed=33)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r, main='A)')
points(locs1$coords, pch=16)
plot(r, main='B)')
points(locs2$coords, pch=16)
locs=list(locs1, locs2)

# for (nrgk in 10:40){
#   r <- caPr.disc[[1]]
#   locs1 <- simLocW(W1, r, beta=0, seed=nrgk)
#   locs2 <- simLocW(W2, r, beta=0, seed=nrgk)
#   print(paste("seed:", nrgk))
#   print(sum(locs1$status))
#   print(sum(locs2$status))
#   print("***")
# }

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.25
Alpha.case2 <- 0.45
Alpha.ctrl2 <- 0.15
beta.case1 <- c(2.2, 0.75, -0.50)
beta.ctrl1 <- c(3.5, 0.5, 0.5)
beta.case2 <- c(1.5, 0.8, -0.5)
beta.ctrl2 <- c(3.5, 0.5, 0.5)

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

save_output(data, "simMVGP_sensitivity_data_1.json")

params <- list(
  Alpha.cases=list(Alpha.case1, Alpha.case2),
  Alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
  beta.cases=list(beta.case1, beta.case2),
  beta.ctrls=list(beta.ctrl1, beta.ctrl2),
  Theta1=Theta1,
  Theta2=Theta2,
  W=W
)
save_output(params, "simMVGP_sensitivity_params_1.json")


#########################
#### medium range difference
#########################


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

#### Simulate multivariate spatial process from
#### a linear coregionalization model
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Theta1 <- 2
Theta2 <- 6
Sigma1 <- Exponential(D, range=Theta1, phi=3)
Sigma2 <- Exponential(D, range=Theta2, phi=3)
set.seed(40)
w1 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma1)
w2 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma2)
A <- matrix(c(0.25, 0.25, -1, -0.5), nrow=2)

W <- c()
for (s in 1:length(w1)){
  w_s <- A %*% matrix(c(w1[s], w2[s]))
  W <- c(W, w_s)
}

W1 <- W[seq(1, length(W), by=2)]
W2 <- W[seq(2, length(W), by=2)]
hist(W1)
hist(W2)

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=33)
locs2 <- simLocW(W2, r, beta=0, seed=33)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r, main='A)')
points(locs1$coords, pch=16)
plot(r, main='B)')
points(locs2$coords, pch=16)
locs=list(locs1, locs2)

# for (nrgk in 10:40){
#   r <- caPr.disc[[1]]
#   locs1 <- simLocW(W1, r, beta=0, seed=nrgk)
#   locs2 <- simLocW(W2, r, beta=0, seed=nrgk)
#   print(paste("seed:", nrgk))
#   print(sum(locs1$status))
#   print(sum(locs2$status))
#   print("***")
# }

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.25
Alpha.case2 <- 0.45
Alpha.ctrl2 <- 0.15
beta.case1 <- c(2.45, 0.75, -0.50)
beta.ctrl1 <- c(3.5, 0.5, 0.5)
beta.case2 <- c(1.75, 0.8, -0.5)
beta.ctrl2 <- c(3.5, 0.5, 0.5)

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

save_output(data, "simMVGP_sensitivity_data_2.json")

params <- list(
  Alpha.cases=list(Alpha.case1, Alpha.case2),
  Alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
  beta.cases=list(beta.case1, beta.case2),
  beta.ctrls=list(beta.ctrl1, beta.ctrl2),
  Theta1=Theta1,
  Theta2=Theta2,
  W=W
)
save_output(params, "simMVGP_sensitivity_params_2.json")


#########################
#### high range difference
#########################


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

#### Simulate multivariate spatial process from
#### a linear coregionalization model
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Theta1 <- 5
Theta2 <- 15
Sigma1 <- Exponential(D, range=Theta1, phi=3)
Sigma2 <- Exponential(D, range=Theta2, phi=3)
set.seed(40)
w1 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma1)
w2 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma2)
A <- matrix(c(0.25, 0.25, -1, -0.5), nrow=2)

W <- c()
for (s in 1:length(w1)){
  w_s <- A %*% matrix(c(w1[s], w2[s]))
  W <- c(W, w_s)
}

W1 <- W[seq(1, length(W), by=2)]
W2 <- W[seq(2, length(W), by=2)]
hist(W1)
hist(W2)

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=33)
locs2 <- simLocW(W2, r, beta=0, seed=33)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r, main='A)')
points(locs1$coords, pch=16)
plot(r, main='B)')
points(locs2$coords, pch=16)
locs=list(locs1, locs2)

# for (nrgk in 10:40){
#   r <- caPr.disc[[1]]
#   locs1 <- simLocW(W1, r, beta=0, seed=nrgk)
#   locs2 <- simLocW(W2, r, beta=0, seed=nrgk)
#   print(paste("seed:", nrgk))
#   print(sum(locs1$status))
#   print(sum(locs2$status))
#   print("***")
# }

#### Simulate counts given locations
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.25
Alpha.case2 <- 0.45
Alpha.ctrl2 <- 0.15
beta.case1 <- c(2.65, 0.75, -0.50)
beta.ctrl1 <- c(3.5, 0.5, 0.5)
beta.case2 <- c(1.75, 0.8, -0.5)
beta.ctrl2 <- c(3.5, 0.5, 0.5)

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

save_output(data, "simMVGP_sensitivity_data_3.json")

params <- list(
  Alpha.cases=list(Alpha.case1, Alpha.case2),
  Alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
  beta.cases=list(beta.case1, beta.case2),
  beta.ctrls=list(beta.ctrl1, beta.ctrl2),
  Theta1=Theta1,
  Theta2=Theta2,
  W=W
)
save_output(params, "simMVGP_sensitivity_params_3.json")