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


#### Simulation parameters for each level
n_sims <- 25
agg_factor <- 12
dst <- "/Users/brianconroy/Documents/research/dataInt/output/sim_mvgp/"


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


Alpha.case1 <- 0.5
Alpha.ctrl1 <- 0.1
Alpha.case2 <- 1
Alpha.ctrl2 <- 0.25
beta.case1 <- c(0.1, 0.75, -0.50)
beta.ctrl1 <- c(3.5, 0.5, 0.5)
beta.case2 <- c(1, 0.8, -0.5)
beta.ctrl2 <- c(4, 0.5, 0.5)
Theta <- 7
Phi <- 12
Tmat <- matrix(c(8, 0, 0, 9), nrow=2)
H <- Exponential(d, range=Theta, phi=1)

prevalences <- array(NA, c(n_sims, 2))
ps_contribs <- array(NA, c(n_sims, 2))
obs_cells <- array(NA, c(n_sims, 2))
n_specimen <- array(NA, c(n_sims, 2))
for (i in 1:n_sims){
  print(paste("dataset", i))
  
  #### Simulate gaussian process
  Sigma <- kronecker(H, Tmat)
  W <-  mvrnorm(n=1, mu=rep(0, ncol(Sigma)), Sigma)
  W1 <- W[seq(1, ncol(Sigma), by=2)]
  W2 <- W[seq(2, ncol(Sigma), by=2)]
  
  #### Simulate locations
  locs1 <- simLocW(W1, caPr.disc[[1]], beta=0)
  locs2 <- simLocW(W2, caPr.disc[[1]], beta=0)
  obs_cells[i,1] <- sum(locs1$status)
  obs_cells[i,2] <- sum(locs2$status)
  
  ps_contribs[i,1] <- calc_ps_contribution(caPr.disc, locs1, beta.case1, Alpha.case1, beta.ctrl1, Alpha.ctrl1, W1)
  ps_contribs[i,2] <- calc_ps_contribution(caPr.disc, locs2, beta.case2, Alpha.case2, beta.ctrl2, Alpha.ctrl2, W2)
  
  #### Simulate counts given locations
  case.data1 <- simConditionalGp2(caPr.disc, locs1, beta.case1, Alpha.case1, W1)
  ctrl.data1 <- simConditionalGp2(caPr.disc, locs1, beta.ctrl1, Alpha.ctrl1, W1)
  
  case.data2 <- simConditionalGp2(caPr.disc, locs2, beta.case2, Alpha.case2, W2)
  ctrl.data2 <- simConditionalGp2(caPr.disc, locs2, beta.ctrl2, Alpha.ctrl2, W2)
  
  prevalences[i,1] <- sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y)
  prevalences[i,2] <- sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y)
  n_specimen[i,1] <- sum(case.data1$y + ctrl.data1$y)
  n_specimen[i,2] <- sum(case.data2$y + ctrl.data2$y)
  
  params <- list(
    beta.cases=list(beta.case1, beta.case2),
    beta.ctrls=list(beta.ctrl1, beta.ctrl2),
    alpha.cases=list(Alpha.case1, Alpha.case2),
    alpha.ctrls=list(Alpha.ctrl1, Alpha.ctrl2),
    Tmat=Tmat,
    Theta=Theta,
    Phi=Phi,
    W=combine_ws(W1, W2)
  )
  
  save_output(params, paste("params_", i, ".json", sep=""), dst=dst)
  
  data <- list(
    locs=list(locs1, locs2),
    case.data=list(case.data1, case.data2),
    ctrl.data=list(ctrl.data1, ctrl.data2)
  )
  save_output(data, paste("data_", i, ".json", sep=""), dst=dst)
  
}

print(colMeans(prevalences))
print(colMeans(ps_contribs))
print(colMeans(obs_cells))
print(colMeans(n_specimen))
