############################
# Generate data for project 3
# simulation 1
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


agg_factor <- 12
n_sims <- 25


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=agg_factor) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])
print(n_values(caPr.disc_all[[1]][[1]]))
print(mean(area(caPr.disc_all[[1]][[1]])[]))


#### Alternating trend
level <- "alternating"
sim_name <- paste("simPsTemporal", level, sep="_")


Theta <- 7
Phi <- 12
U <- c(-1.5, 0, -1, 0.5, -1.25, 1, -0.6)
Alpha.case <- 0.75
Alpha.ctrl <- 0.20
beta.case <- c(0.75, 0.5, 0.5)
beta.ctrl <- c(3.75, 0.75, 0.75)
prevalences <- array(NA, c(n_sims, length(U)))
ps_contribs <- array(NA, c(n_sims, length(U)))
obs_cells <- array(NA, c(n_sims, length(U)))
n_specimen <- array(NA, c(n_sims, length(U)))
for (d in 1:n_sims){
  print(d)
  
  #### Simulate gaussian process 
  cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
  coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
  D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  Sigma <- Exponential(D, range=Theta, phi=Phi)
  W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
  N <- length(W)
  
  #### Simulate locational data over time
  r <- caPr.disc_all[[1]][[1]]
  locs_time <- list()
  counter <- 1
  for (h in U){
    locs <- simLocW(h + W, r, beta=0)
    while(sum(locs$status) < 10){
      locs <- simLocW(h + W, r, beta=0)
    }
    locs_time[[counter]] <- locs
    obs_cells[d,counter] <- sum(locs$status)
    ps_contribs[d,counter] <- calc_ps_contribution(caPr.disc_all[[counter]], locs, beta.case, Alpha.case, beta.ctrl, Alpha.ctrl, W)
    counter <- counter + 1
  }
  
  #### Simulate case and control data
  case.data_time <- list()
  ctrl.data_time <- list()
  counter <- 1
  for (h in 1:length(U)){
    locs_i <- locs_time[[h]]
    u_i <- U[h]
    cov.disc_i <- caPr.disc_all[[h]]
    case.data <- simConditionalGp2(cov.disc_i, locs_i, beta.case, alpha.case, u_i + W, center=FALSE)
    ctrl.data <- simConditionalGp2(cov.disc_i, locs_i, beta.ctrl, alpha.ctrl, u_i + W, center=FALSE)
    case.data_time[[h]] <- case.data
    ctrl.data_time[[h]] <- ctrl.data
    prevalences[d,counter] <- sum(case.data$y)/sum(case.data$y + ctrl.data$y)
    n_specimen[d,counter] <- sum(case.data$y + ctrl.data$y)
    counter <- counter + 1
  }
  
  data <- list(
    case.data=case.data_time,
    ctrl.data=ctrl.data_time,
    locs=locs_time
  )
  dst <- "/Users/brianconroy/Documents/research/dataInt/output/sim_temporal/"
  save_output(data, paste("data", level, "_", d, ".json", sep=""), dst=dst)
  
  params <- list(
    Alpha.case=alpha.case,
    Alpha.ctrl=alpha.ctrl,
    beta.case=beta.case,
    beta.ctrl=beta.ctrl,
    Theta=Theta,
    Phi=Phi,
    W=W,
    U=U
  )
  save_output(params, paste("params_", level, "_", d, ".json", sep=""), dst=dst)
  
}
print(colMeans(prevalences))
print(colMeans(ps_contribs))
print(colMeans(obs_cells))
print(apply(obs_cells, 2, median))
print(colMeans(n_specimen))


#### Decreasing trend
level <- "decreasing"
sim_name <- paste("simPsTemporal", level, sep="_")


Theta <- 7
Phi <- 12
U <- c(1.0, 0.0, -0.5, -1.0, -2.0, -2.5, -3.0)
Alpha.case <- 1
Alpha.ctrl <- 0.20
beta.case <- c(0.75, 0.5, 0.5)
beta.ctrl <- c(4.25, 0.75, 0.75)
prevalences <- array(NA, c(n_sims, length(U)))
ps_contribs <- array(NA, c(n_sims, length(U)))
obs_cells <- array(NA, c(n_sims, length(U)))
n_specimen <- array(NA, c(n_sims, length(U)))
for (d in 1:n_sims){
  print(d)
  
  #### Simulate gaussian process 
  cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
  coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
  D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  Sigma <- Exponential(D, range=Theta, phi=Phi)
  
  W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
  N <- length(W)
  
  #### Simulate locational data over time
  r <- caPr.disc_all[[1]][[1]]
  locs_time <- list()
  counter <- 1
  for (h in U){
    locs <- simLocW(h + W, r, beta=0)
    locs_time[[counter]] <- locs
    obs_cells[d,counter] <- sum(locs$status)
    ps_contribs[d,counter] <- calc_ps_contribution(caPr.disc_all[[counter]], locs, beta.case, Alpha.case, beta.ctrl, Alpha.ctrl, W)
    counter <- counter + 1
  }
  
  #### Simulate case and control data
  case.data_time <- list()
  ctrl.data_time <- list()
  counter <- 1
  for (h in 1:length(U)){
    locs_i <- locs_time[[h]]
    u_i <- U[h]
    cov.disc_i <- caPr.disc_all[[h]]
    case.data <- simConditionalGp2(cov.disc_i, locs_i, beta.case, Alpha.case, u_i + W, center=FALSE)
    ctrl.data <- simConditionalGp2(cov.disc_i, locs_i, beta.ctrl, Alpha.ctrl, u_i + W, center=FALSE)
    case.data_time[[h]] <- case.data
    ctrl.data_time[[h]] <- ctrl.data
    prevalences[d,counter] <- sum(case.data$y)/sum(case.data$y + ctrl.data$y)
    n_specimen[d,counter] <- sum(case.data$y + ctrl.data$y)
    counter <- counter + 1
  }
  
  data <- list(
    case.data=case.data_time,
    ctrl.data=ctrl.data_time,
    locs=locs_time
  )
  dst <- "/Users/brianconroy/Documents/research/dataInt/output/sim_temporal/"
  save_output(data, paste("data_", level, "_", d, ".json", sep=""), dst=dst)
  
  params <- list(
    Alpha.case=Alpha.case,
    Alpha.ctrl=Alpha.ctrl,
    beta.case=beta.case,
    beta.ctrl=beta.ctrl,
    Theta=Theta,
    Phi=Phi,
    W=W,
    U=U
  )
  save_output(params, paste("params_", level, "_", d, ".json", sep=""), dst=dst)
  
}
print(colMeans(prevalences))
print(colMeans(ps_contribs))
print(colMeans(obs_cells))
print(apply(obs_cells, 2, median))
print(colMeans(n_specimen))
