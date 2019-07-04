############################
# Generate data for project 3
# simulation 1
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


level <- "alternating"
sim_name <- paste("simPsTemporal", level, sep="_")


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=7) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])
print(n_values(caPr.disc_all[[1]][[1]]))
print(mean(area(caPr.disc_all[[1]][[1]])[]))


#### Simulate gaussian process 
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(D, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Specify temporal parameters
# U <- c(-3, -2.5, -2, -1, -0.5, 0, 1)
# U <- c(1.0, 0.0, -0.5, -1.0, -2.0, -2.5, -3.0)
# U <- rep(0, 7)
U <- c(-1.5, 0, -1, 0.5, -1.25, 1, -0.6)


#### Simulate locational data over time
r <- caPr.disc_all[[1]][[1]]
locs_time <- list()
counter <- 1
for (h in U){
  locs <- simLocW(h + W, r, beta=0, seed=11)
  locs_time[[counter]] <- locs
  counter <- counter + 1
}


for (l in locs_time){
  print(sum(l$status))
}


#### Simulate case and control data
alpha.case <- 0.75
alpha.ctrl <- 0.25
beta.case <- c(1, 0.5, 0.5)
beta.ctrl <- c(3, 0.75, 0.75)
case.data_time <- list()
ctrl.data_time <- list()
for (h in 1:length(U)){
  locs_i <- locs_time[[h]]
  u_i <- U[h]
  cov.disc_i <- caPr.disc_all[[h]]
  case.data <- simConditionalGp2(cov.disc_i, locs_i, beta.case, alpha.case, u_i + W, seed=42, center=FALSE)
  ctrl.data <- simConditionalGp2(cov.disc_i, locs_i, beta.ctrl, alpha.ctrl, u_i + W, seed=40, center=FALSE)
  case.data_time[[h]] <- case.data
  ctrl.data_time[[h]] <- ctrl.data
}


for (h in 1:length(U)){
  print(sum(case.data_time[[h]]$y)/sum(case.data_time[[h]]$y + ctrl.data_time[[h]]$y))
}


data <- list(
  case.data=case.data_time,
  ctrl.data=ctrl.data_time,
  locs=locs_time
)


save_output(data, paste(sim_name, "data.json", sep="_"))

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
save_output(params, paste(sim_name, "params.json", sep="_"))
