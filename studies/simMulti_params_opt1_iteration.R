###########################
# Set simulation parameters
###########################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')

# iteration
#   goal: assess multispecies model performance under increasing
#   number of species.
#     does it work with more?
#     does it get better with more?
#   performance measures:
#     log odds accuracy
#     bias
#     w accuracy
#   computational efficiency:
#     number of samples drawn
#     burnin
#     runtime?
#   species numbers: c(2, 4, 6, 8)
# w: same


n_species <- 8


Theta <- 6
Phi <- 12

caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)

cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)

Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
beta.case1 <- c(0.0, 0.75, -0.50)
beta.ctrl1 <- c(2.75, 1.00, 0.50)

alpha.cases <- list(Alpha.case1)
alpha.ctrls <- list(Alpha.ctrl1)
beta.cases <- list(beta.case1)
beta.ctrls <- list(beta.ctrl1)

set.seed(42)
for (i in 1:c(n_species-1)){
  alpha.cases[[i+1]] <- round(alpha.cases[[1]] + abs(rnorm(1, sd=0.25)), 2)
  alpha.ctrls[[i+1]] <- round(alpha.ctrls[[1]] - abs(rnorm(1, sd=0.25)), 2)
  beta.cases[[i+1]] <- round(beta.cases[[1]] + rnorm(3, sd=0.25), 2)
  beta.ctrls[[i+1]] <- round(beta.ctrls[[1]] + rnorm(3, sd=0.25), 2)
}

#### Simulate locations and counts
seeds <- 1:n_species
r <- caPr.disc[[1]]
cov.disc <- caPr.disc
locs <- list()
case.datas <- list()
ctrl.datas <- list()
for (i in 1:n_species){
  locs[[i]] <- simLocW(W, r, beta=0, seed=seeds[i])
  print(sum(locs[[i]]$status))
  case.datas[[i]] <- simConditionalGp2(cov.disc, locs[[i]], beta.cases[[i]], alpha.cases[[i]], W, seed=42)
  ctrl.datas[[i]] <- simConditionalGp2(cov.disc, locs[[i]], beta.ctrls[[i]], alpha.ctrls[[i]], W, seed=40)
  print(sum(case.datas[[i]]$y)/sum(case.datas[[i]]$y + ctrl.datas[[i]]$y))
}

dat <- list(
  case.datas=case.datas,
  ctrl.datas=ctrl.datas,
  locs=locs,
  beta.cases=beta.cases,
  beta.ctrls=beta.ctrls,
  alpha.cases=alpha.cases,
  alpha.ctrls=alpha.ctrls,
  phi=Phi,
  theta=Theta,
  W=W
)
save_output(dat, 'simMulti_data_iteration.json')
