############################
# Fit the spatiotemporal
# shared latent process model
# under differing prior
# variances
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


sim_name <- "simPsTemporal_sensitivity"


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=7) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
  plot(caPr_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])


#### Load simulation parameters
params <- load_output(paste(sim_name, '_params.json', sep=''))
alpha.case <- params$Alpha.case
alpha.ctrl <- params$Alpha.ctrl
beta.case <- params$beta.case
beta.ctrl <- params$beta.ctrl
W <- params$W
U <- params$U
Theta <- params$Theta
Phi <- params$Phi

cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


#### Load Data
data_ <- load_output(paste(sim_name, '_data.json', sep=''))
data <- reformat_saved_data(data_)

for (h in 1:7){
  print(sum(data$case.data[[h]]$y)/sum(data$case.data[[h]]$y + data$ctrl.data[[h]]$y))
}


################################
#### Fit spatiotemporal 
#### preferential sampling model
################################


#### Specify MCMC parameters
w_initial=rep(0, length(W))
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
prior_alpha_ca_var=10
prior_alpha_co_mean=0
prior_alpha_co_var=10
prior_u_mean=0
prior_u_var=10

n.sample=1000
burnin=500

proposal.sd.theta=0.2
prior_theta=get_gamma_prior(3, 10)
prior_phi=get_gamma_prior(15, 10)


#### Fit model
output <- prefSampleTemporal(data, D, n.sample, burnin,
                             L_w, L_ca, L_co, L_a_ca, L_a_co, L_u,
                             proposal.sd.theta=proposal.sd.theta,
                             m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, m_u=m_u,
                             target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_u=target_u,
                             self_tune_w=self_tune_w, self_tune_aca=self_tune_aca, self_tune_aco=self_tune_aco, self_tune_ca=self_tune_ca, self_tune_co=self_tune_co, self_tune_u=self_tune_u,
                             beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                             theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial, u_initial=u_initial,
                             prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var,
                             prior_u_mean=prior_u_mean, prior_u_var=prior_u_var)


#### Additional burnin
output <- burnin_after_temporal(output, n.burn=50)


#### Generate additional samples
output <- continue_mcmc_temporal(data, D, output, n.sample=2000)


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


output$description <- sim_name
save_output(output, paste("output_", sim_name, ".json", sep=""))


################################
#### Fit nontemporal 
#### preferential sampling model
################################


data_pooled <- pool_temporal_data(data)


## W initial value
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)
w_output <- logisticGp(y=data_pooled$locs$status, D, n.sample=1000, burnin=200, L=10,
                       prior_phi=prior_phi, prior_theta=prior_theta)
view_logistic_output(w_output)
save_output(w_output, paste("w_inival_output_", sim_name, ".json", sep=""))
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

# Beta & alpha initial values
ini_case <- glm(data_pooled$case.data$y ~ data_pooled$case.data$x + w_i[data_pooled$locs$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(data_pooled$ctrl.data$y ~ data_pooled$ctrl.data$x.standardised + w_i[data_pooled$locs$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

n.sample <- 5000
burnin <- 500
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

output_ps <- prefSampleGpCC(data_pooled, D, n.sample, burnin,
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=proposal.sd.theta,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                         target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                         delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                         beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                         theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                         prior_phi=prior_phi, prior_theta=prior_theta,
                         prior_alpha_ca_var, prior_alpha_co_var)


# optionally burnin the output more
output_ps <- burnin_after(output_ps, n.burn=2000)


# optionally continue running if necessary
output_ps <- continueMCMC(data_pooled, D, output_ps, n.sample=2000)


w.hat <- colMeans(output_ps$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)

par(mfrow=c(2, 3))
view_tr(output_ps$samples.beta.ca[,1], beta.case[1])
view_tr(output_ps$samples.beta.ca[,2], beta.case[2])
view_tr(output_ps$samples.beta.ca[,3], beta.case[3])
padded_plot(output_ps$samples.beta.co[,1], beta.ctrl[1])
padded_plot(output_ps$samples.beta.co[,2], beta.ctrl[2])
padded_plot(output_ps$samples.beta.co[,3], beta.ctrl[3])

par(mfrow=c(1,2))
view_tr(output_ps$samples.alpha.ca, alpha.case)
view_tr(output_ps$samples.alpha.co, alpha.ctrl)


## save results
tag <- paste(sim_name, "_aggregated_ps", sep="")
output_ps$description <- tag
save_output(output_ps, paste('output_', tag, ".json", sep=""))
