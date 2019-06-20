############################
# Fit the spatiotemporal
# shared latent process model
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (t in years){
  caPr_y <- load_prism_pcs_time(t)
  caPr.disc_y <- aggregate(caPr_y, fact=9) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])


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
U <- c(-3, -2.5, -2, -1, -0.5, 0, 1)


#### Simulate locational data over time
r <- caPr.disc_all[[1]][[1]]
locs_time <- list()
counter <- 1
for (t in U){
  locs <- simLocW(t + W, r, beta=0, seed=11)
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
for (i in 1:length(U)){
  locs_i <- locs_time[[i]]
  u_i <- U[i]
  cov.disc_i <- caPr.disc_all[[i]]
  case.data <- simConditionalGp2(cov.disc_i, locs_i, beta.case, alpha.case, u_i + W, seed=42, center=FALSE)
  ctrl.data <- simConditionalGp2(cov.disc_i, locs_i, beta.ctrl, alpha.ctrl, u_i + W, seed=40, center=FALSE)
  case.data_time[[i]] <- case.data
  ctrl.data_time[[i]] <- ctrl.data
}


for (i in 1:length(U)){
  print(sum(case.data_time[[i]]$y)/sum(case.data_time[[i]]$y + ctrl.data_time[[i]]$y))
}
for (i in 1:length(U)){
  print(sum(case.data_time[[i]]$y))
}
for (i in 1:length(U)){
  print(sum(ctrl.data_time[[i]]$y))
}
for (i in 1:length(U)){
  print(case.data_time[[i]]$y)
}
for (h in case.data_time){
  print(hist(h$x[,3]))
}


data <- list(
  case.data=case.data_time,
  ctrl.data=ctrl.data_time,
  locs=locs_time
)

w_initial=W + 2*rnorm(length(W))
beta_ca_initial=rep(0, 3)
beta_co_initial=rep(0, 3)
alpha_ca_initial=0
alpha_co_initial=0
theta_initial=Theta
phi_initial=Phi
u_initial=U

m_aca=2000
m_aco=2000
m_ca=700
m_co=700
m_w=700
m_u=700
target_aca=0.75
target_aco=0.75
target_ca=0.75
target_co=0.75
target_w=0.75
target_u=0.75
self_tune_w=TRUE
self_tune_aca=TRUE
self_tune_aco=TRUE
self_tune_ca=TRUE
self_tune_co=TRUE
self_tune_u=TRUE

L_w=8
L_ca=8
L_co=8
L_a_ca=8
L_a_co=8
offset=1

prior_alpha_ca_mean=0
prior_alpha_ca_var=4
prior_alpha_co_mean=0
prior_alpha_co_var=4

n.sample=2000
burnin=200

proposal.sd.theta=0.2
prior_theta=c(3,2)
d=D
prior_phi <- c(3, 40)

par(mfrow=c(2,3))
plot(samples.beta.ca[,1], type='l'); abline(h=beta.case[1], col=2)
plot(samples.beta.ca[,2], type='l'); abline(h=beta.case[2], col=2)
plot(samples.beta.ca[,3], type='l'); abline(h=beta.case[3], col=2)

plot(samples.beta.co[,1], type='l'); abline(h=beta.ctrl[1], col=2)
plot(samples.beta.co[,2], type='l'); abline(h=beta.ctrl[2], col=2)
plot(samples.beta.co[,3], type='l'); abline(h=beta.ctrl[3], col=2)

par(mfrow=c(1,2))
plot(samples.alpha.ca, type='l'); abline(h=alpha.case, col=2)
plot(samples.alpha.co, type='l'); abline(h=alpha.ctrl, col=2)

plot(samples.theta, type='l'); abline(h=Theta, col=2)
plot(samples.phi, type='l'); abline(h=Phi, col=2)

view_tr_w(samples.w, W)
w.hat <- colMeans(samples.w)
plot(x=w.hat, y=W); abline(0,1,col=2)
