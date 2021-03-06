library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')




#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=3)
n_values(caPr.disc[[1]])
plot(caPr.disc)

cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=6, phi=12)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
r <- caPr.disc[[1]]
r[][!is.na(r[])] <- W
plot(r)











#### Level: Low preferential sampling
# Median preferential sampling contribution: 17.48%
prevalences <- c()
ps_contribs <- c()
obs_cells <- c()
n_specimen <- c()

Alpha.case <- 0.5
Alpha.ctrl <- 0.3
beta.case <- c(1, 0.75, 0.25)
beta.ctrl <- c(3, 1, 0.5)
Theta <- 7
Phi <- 12
for (i in 1:n_sims){
  
  #### Simulate gaussian process
  cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
  coords <- xyFromCell(caPr.disc, cell=cells.all)
  d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  Sigma <- Exponential(d, range=Theta, phi=Phi)
  W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
  N <- length(W)
  
  #### Simulate locations
  locs <-simLocW(W, caPr.disc[[1]], beta=0)
  if (sum(locs$status < 30)){
    locs <-simLocW(W, caPr.disc[[1]], beta=0)
  }
  obs_cells <- c(obs_cells, sum(locs$status))

  ps_contribs <- c(ps_contribs, 
                   calc_ps_contribution(caPr.disc, locs, beta.case, Alpha.case, beta.ctrl, Alpha.ctrl, W))
  
  
  #### Simulate counts given locations
  case.data <- simConditionalGp2(caPr.disc, locs, beta.case, Alpha.case, W)
  ctrl.data <- simConditionalGp2(caPr.disc, locs, beta.ctrl, Alpha.ctrl, W)
  prevalences <- c(prevalences, sum(case.data$y)/sum(case.data$y + ctrl.data$y))
  n_specimen <- c(n_specimen, sum(case.data$y + ctrl.data$y))
  
  params <- list(
    sampling="low",
    beta.case=beta.case,
    beta.ctrl=beta.ctrl,
    alpha.case=Alpha.case,
    alpha.ctrl=Alpha.ctrl,
    Theta=Theta,
    Phi=Phi,
    W=W
  )
  
  dst <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration_v1/"
  save_output(params, paste("params_low_", i, ".json", sep=""), dst=dst)

  data <- list(
    case.data=case.data,
    ctrl.data=ctrl.data,
    locs=locs
  )
  save_output(data, paste("data_low_", i, ".json", sep=""), dst=dst)

}

print(summary(prevalences))
print(summary(ps_contribs))
print(summary(obs_cells))
print(summary(n_specimen))


#### Level: High preferential sampling
# Median preferential sampling contribution: 17.48%
prevalences <- c()
ps_contribs <- c()
obs_cells <- c()
n_specimen <- c()

Alpha.case <- 1
Alpha.ctrl <- -0.5
beta.case <- c(-1.5, 0.25, 0.25)
beta.ctrl <- c(3.5, 1, 0.5)
Theta <- 7
Phi <- 12
for (i in 1:n_sims){
  
  #### Simulate gaussian process
  cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
  coords <- xyFromCell(caPr.disc, cell=cells.all)
  d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  Sigma <- Exponential(d, range=Theta, phi=Phi)
  W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
  N <- length(W)
  
  #### Simulate locations
  locs <-simLocW(W, caPr.disc[[1]], beta=0)
  if (sum(locs$status < 30)){
    locs <-simLocW(W, caPr.disc[[1]], beta=0)
  }
  obs_cells <- c(obs_cells, sum(locs$status))
  
  ps_contribs <- c(ps_contribs, 
                   calc_ps_contribution(caPr.disc, locs, beta.case, Alpha.case, beta.ctrl, Alpha.ctrl, W))
  
  
  #### Simulate counts given locations
  case.data <- simConditionalGp2(caPr.disc, locs, beta.case, Alpha.case, W)
  ctrl.data <- simConditionalGp2(caPr.disc, locs, beta.ctrl, Alpha.ctrl, W)
  prevalences <- c(prevalences, sum(case.data$y)/sum(case.data$y + ctrl.data$y))
  n_specimen <- c(n_specimen, sum(case.data$y + ctrl.data$y))
  
  params <- list(
    sampling="high",
    beta.case=beta.case,
    beta.ctrl=beta.ctrl,
    alpha.case=Alpha.case,
    alpha.ctrl=Alpha.ctrl,
    Theta=Theta,
    Phi=Phi,
    W=W
  )
  
  dst <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration_v1/"
  save_output(params, paste("params_high_", i, ".json", sep=""), dst=dst)

  data <- list(
    case.data=case.data,
    ctrl.data=ctrl.data,
    locs=locs
  )
  save_output(data, paste("data_high_", i, ".json", sep=""), dst=dst)
  
}

print(summary(prevalences))
print(summary(ps_contribs))
print(summary(obs_cells))
print(summary(n_specimen))
