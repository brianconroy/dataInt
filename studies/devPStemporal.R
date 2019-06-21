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
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
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
# U <- c(1.0, 0.0, -0.5, -1.0, -2.0, -2.5, -3.0)

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


#### Specify MCMC parameters
w_initial=rep(0, length(W))#rnorm(length(W)) # W + 2*rnorm(length(W))
beta_ca_initial=rep(0, 3)
beta_co_initial=rep(0, 3)
alpha_ca_initial=0
alpha_co_initial=0
theta_initial=Theta
phi_initial=Phi
u_initial=rep(0, length(U))

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
L_u=8

prior_alpha_ca_mean=0
prior_alpha_ca_var=4
prior_alpha_co_mean=0
prior_alpha_co_var=4
prior_u_mean=0
prior_u_var=4

n.sample=1500
burnin=200

proposal.sd.theta=0.2
prior_theta=c(3,2)
prior_phi <- c(3, 40)


#### Fit model
output <- prefSampleTemporal(data, D, n.sample, burnin,
                             L_w, L_ca, L_co, L_a_ca, L_a_co,
                             proposal.sd.theta=proposal.sd.theta,
                             m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, m_u=m_u,
                             target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_u=target_u,
                             self_tune_w=self_tune_w, self_tune_aca=self_tune_aca, self_tune_aco=self_tune_aco, self_tune_ca=self_tune_ca, self_tune_co=self_tune_co, self_tune_u=self_tune_u,
                             beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                             theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial, u_initial=u_initial,
                             prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var,
                             prior_u_mean=prior_u_mean, prior_u_var=prior_u_var,
                             offset=1)


#### View results
par(mfrow=c(2,3))
plot(output$samples.beta.ca[,1], type='l'); abline(h=beta.case[1], col=2)
plot(output$samples.beta.ca[,2], type='l'); abline(h=beta.case[2], col=2)
plot(output$samples.beta.ca[,3], type='l'); abline(h=beta.case[3], col=2)
plot(output$samples.beta.co[,1], type='l'); abline(h=beta.ctrl[1], col=2)
plot(output$samples.beta.co[,2], type='l'); abline(h=beta.ctrl[2], col=2)
plot(output$samples.beta.co[,3], type='l'); abline(h=beta.ctrl[3], col=2)

par(mfrow=c(1,2))
plot(output$samples.alpha.ca, type='l'); abline(h=alpha.case, col=2)
plot(output$samples.alpha.co, type='l'); abline(h=alpha.ctrl, col=2)

plot(output$samples.theta, type='l'); abline(h=Theta, col=2)
plot(output$samples.phi, type='l'); abline(h=Phi, col=2)

view_tr_w(output$samples.w, W)
w.hat <- colMeans(output$samples.w)
plot(x=w.hat, y=W); abline(0,1,col=2)

par(mfrow=c(2,4))
plot(output$samples.u[,1], type='l'); abline(h=U[1], col=2)
plot(output$samples.u[,2], type='l'); abline(h=U[2], col=2)
plot(output$samples.u[,3], type='l'); abline(h=U[3], col=2)
plot(output$samples.u[,4], type='l'); abline(h=U[4], col=2)
plot(output$samples.u[,5], type='l'); abline(h=U[5], col=2)
plot(output$samples.u[,6], type='l'); abline(h=U[6], col=2)
plot(output$samples.u[,7], type='l'); abline(h=U[7], col=2)

w.hat <- colMeans(output$samples.w)
u.hat <- colMeans(output$samples.u)
par(mfrow=c(2,4))
plot(x=W+U[1], y=w.hat+u.hat[1]); abline(0,1)
plot(x=W+U[2], y=w.hat+u.hat[2]); abline(0,1)
plot(x=W+U[3], y=w.hat+u.hat[3]); abline(0,1)
plot(x=W+U[4], y=w.hat+u.hat[4]); abline(0,1)
plot(x=W+U[5], y=w.hat+u.hat[5]); abline(0,1)
plot(x=W+U[6], y=w.hat+u.hat[6]); abline(0,1)
plot(x=W+U[7], y=w.hat+u.hat[7]); abline(0,1)

par(mfrow=c(2,4))
plot(x=W+U[1], y=w.hat); abline(0,1)
plot(x=W+U[2], y=w.hat); abline(0,1)
plot(x=W+U[3], y=w.hat); abline(0,1)
plot(x=W+U[4], y=w.hat); abline(0,1)
plot(x=W+U[5], y=w.hat); abline(0,1)
plot(x=W+U[6], y=w.hat); abline(0,1)
plot(x=W+U[7], y=w.hat); abline(0,1)
