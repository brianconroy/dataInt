############################
# compare
  # normal alphas w/ uncalibrated initial vals
  # normal alphas w/ good initial vals
  # truncated normal alphas w/ uncalibrated initial vals
  # truncated normal alphas w/ calibrated initial vals

# run over fixed simulation parameters
############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


sim_name <- "priorcompare"


#### Load simulation parameters
params <- load_sim_params_priorcompare()
Alpha.case <- params$alpha.case
Alpha.ctrl <- params$alpha.ctrl
beta.case <- as.numeric(strsplit(params$beta.case, split=" ")[[1]])
beta.ctrl <- as.numeric(strsplit(params$beta.ctrl, split=" ")[[1]])


Theta <- 6
Phi <- 12


#### Priors
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)
prior_alpha_shape <- 0.05
prior_alpha_scale <- 10


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate gaussian process
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)
hist(W)


#### Simulate locations
r <- caPr.disc[[1]]
locs <- simLocW(W, r, beta=0, seed=11) # 42
sum(locs$status)
hist(W)
plot(r)
points(locs$coords)


#### Simulate counts given locations
cov.disc <- caPr.disc
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W, seed=42)
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40)
print(sum(case.data$y)/sum(case.data$y + ctrl.data$y))


save_true_params_general('priorcompare')


data <- list(
  loc=locs,
  case.data=case.data,
  ctrl.data=ctrl.data
)


################################
#### Normal priors
#### Uncalibrated initial values
################################


#### Tuning parameters
n.sample <- 2000
burnin <- 0
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

set.seed(241)
beta_ca_i <- abs(rnorm(3))
beta_co_i <- abs(rnorm(3))
alpha_ca_i <- runif(1, 2, 3)
alpha_co_i <- runif(1, -3, -2)
theta_i <- runif(1, 9, 10)
phi_i <- runif(1, 6, 8)
w_i <- rep(0, length(W))

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4

output_uncal <- prefSampleGpCC(data, n.sample, burnin,
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=proposal.sd.theta,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                         target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                         delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                         beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                         theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                         prior_phi=prior_phi, prior_theta=prior_theta,
                         prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)

# optionally burnin the output more
output_uncal <- burnin_after(output_uncal, n.burn=6000)


# optionally continue running if necessary
output_uncal <- continueMCMC(data, output_uncal, n.sample=2000)


plot(output_uncal$deltas_w)
plot(output_uncal$deltas_ca)
plot(output_uncal$deltas_co)
plot(output_uncal$deltas_aca)
plot(output_uncal$deltas_aco)
print(output_uncal$accept)

plot(apply(output_uncal$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_uncal$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output_uncal$samples.w, w_true=W)

view_tr(output_uncal$samples.theta, Theta)
print(mean(output_uncal$samples.theta)); print(Theta)

view_tr(output_uncal$samples.phi, Phi)
print(mean(output_uncal$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
view_tr(output_uncal$samples.beta.ca[,1], beta.case[1], ylim=c(0,3))
view_tr(output_uncal$samples.beta.ca[,2], beta.case[2], title='B)')
view_tr(output_uncal$samples.beta.ca[,3], beta.case[3], title='C)')
print(colMeans(output_uncal$samples.beta.ca))

view_tr(output_uncal$samples.beta.co[,1], beta.ctrl[1], title='D)')
view_tr(output_uncal$samples.beta.co[,2], beta.ctrl[2], title='E)')
view_tr(output_uncal$samples.beta.co[,3], beta.ctrl[3], title='F)')
print(colMeans(output_uncal$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output_uncal$samples.alpha.ca, Alpha.case, title='A)')
print(mean(output_uncal$samples.alpha.ca))

view_tr(output_uncal$samples.alpha.co, Alpha.ctrl, title='B)')
print(mean(output_uncal$samples.alpha.co))

tag <- "_normal_uncal"
output_uncal$description <- paste(sim_name, tag, sep="")
save_output(output_uncal, paste("output_uncal_", sim_name, tag, ".json", sep=""))
save_params_psgp(paste("params_", sim_name, tag, ".json", sep=""))

w.hat <- colMeans(output_uncal$samples.w)
beta_ca_h <- colMeans(output_uncal$samples.beta.ca)
beta_co_h <- colMeans(output_uncal$samples.beta.co)
alpha_ca_h <- mean(output_uncal$samples.alpha.ca)
alpha_co_h <- mean(output_uncal$samples.alpha.co)
phi_h <- mean(output_uncal$samples.phi)
theta_h <- mean(output_uncal$samples.theta)


################################
#### Normal priors
#### Calibrated initial values
################################


# W initial value
# w_output <- logisticGp(y=locs$status, d, n.sample=1000, burnin=200, L=10,
#                       prior_phi=prior_phi, prior_theta=prior_theta)
# view_logistic_output(w_output)
# save_output(w_output, "w_inival_output_priorcompare.json")
w_output <- load_output("w_inival_output_priorcompare.json")
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

ini_case <- glm(case.data$y ~ case.data$x.standardised + w_i[locs$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data$y ~ ctrl.data$x.standardised + w_i[locs$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4

n.sample <- 2000
burnin <- 500
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

output_cal <- prefSampleGpCC(data, n.sample, burnin,
                                    L_w, L_ca, L_co, L_a_ca, L_a_co,
                                    proposal.sd.theta=proposal.sd.theta,
                                    m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                                    target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                                    self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                                    delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                                    beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                                    theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                                    prior_phi=prior_phi, prior_theta=prior_theta,
                                    prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)


# optionally burnin the output_cal more
output_cal <- burnin_after(output_cal, n.burn=500)


# optionally continue running if necessary
output_cal <- continueMCMC(data, output_cal, n.sample=2000)


plot(output_cal$deltas_w)
plot(output_cal$deltas_ca)
plot(output_cal$deltas_co)
plot(output_cal$deltas_aca)
plot(output_cal$deltas_aco)
print(output_cal$accept)

plot(apply(output_cal$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_cal$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output_cal$samples.w, w_true=W)

view_tr(output_cal$samples.theta, Theta)
print(mean(output_cal$samples.theta)); print(Theta)

view_tr(output_cal$samples.phi, Phi)
print(mean(output_cal$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
view_tr(output_cal$samples.beta.ca[,1], beta.case[1], title='A)')
view_tr(output_cal$samples.beta.ca[,2], beta.case[2], title='B)')
view_tr(output_cal$samples.beta.ca[,3], beta.case[3], title='C)')
print(colMeans(output_cal$samples.beta.ca))

view_tr(output_cal$samples.beta.co[,1], beta.ctrl[1], title='D)')
view_tr(output_cal$samples.beta.co[,2], beta.ctrl[2], title='E)')
view_tr(output_cal$samples.beta.co[,3], beta.ctrl[3], title='F)')
print(colMeans(output_cal$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output_cal$samples.alpha.ca, Alpha.case, title='A)')
print(mean(output_cal$samples.alpha.ca))

view_tr(output_cal$samples.alpha.co, Alpha.ctrl, title='B)')
print(mean(output_cal$samples.alpha.co))

tag <- "_normal_cal"
output_cal$description <- paste(sim_name, tag, sep="")
save_output(output_cal, paste("output_", sim_name, tag, ".json", sep=""))
save_params_psgp(paste("params_", sim_name, tag, ".json", sep=""), prior='normal')


################################
#### Truncated Normal priors
#### Uncalibrated initial values
################################


n.sample <- 2000
burnin <- 500
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

set.seed(241)
beta_ca_i <- abs(rnorm(3))
beta_co_i <- abs(rnorm(3))
alpha_ca_i <- runif(1, 2, 3)
alpha_co_i <- runif(1, -3, -2)
theta_i <- runif(1, 9, 10)
phi_i <- runif(1, 6, 8)
w_i <- rep(0, length(W))

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4
prior_alpha_ca <- c(prior_alpha_ca_mean, prior_alpha_ca_var)
prior_alpha_co <- c(prior_alpha_co_mean, prior_alpha_co_var)

output_tnormal_uncal <- prefSampleGpCC_truncnorm(data, n.sample, burnin,
                                       L_w, L_ca, L_co, L_a_ca, L_a_co,
                                       proposal.sd.theta=proposal.sd.theta,
                                       m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                                       target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                                       self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                                       delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                                       beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                                       theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                                       prior_phi=prior_phi, prior_theta=prior_theta,
                                       prior_alpha_ca=prior_alpha_ca, 
                                       prior_alpha_co=prior_alpha_co)


# optionally burnin the output_normal more
output_tnormal_uncal <- burnin_after(output_tnormal_uncal, n.burn=500)


# optionally continue running if necessary
output_tnormal_uncal <- continueMCMC(data, output_tnormal_uncal, n.sample=2000)


plot(output_tnormal_uncal$deltas_w)
plot(output_tnormal_uncal$deltas_ca)
plot(output_tnormal_uncal$deltas_co)
plot(output_tnormal_uncal$deltas_aca)
plot(output_tnormal_uncal$deltas_aco)
print(output_tnormal_uncal$accept)

plot(apply(output_tnormal_uncal$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_tnormal_uncal$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output_tnormal_uncal$samples.w, w_true=W)

view_tr(output_tnormal_uncal$samples.theta, Theta)
print(mean(output_tnormal_uncal$samples.theta)); print(Theta)

view_tr(output_tnormal_uncal$samples.phi, Phi)
print(mean(output_tnormal_uncal$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
padded_plot(output_tnormal_uncal$samples.beta.ca[,1], beta.case[1])
padded_plot(output_tnormal_uncal$samples.beta.ca[,2], beta.case[2])
padded_plot(output_tnormal_uncal$samples.beta.ca[,3], beta.case[3])
print(colMeans(output_tnormal_uncal$samples.beta.ca))

padded_plot(output_tnormal_uncal$samples.beta.co[,1], beta.ctrl[1])
padded_plot(output_tnormal_uncal$samples.beta.co[,2], beta.ctrl[2])
padded_plot(output_tnormal_uncal$samples.beta.co[,3], beta.ctrl[3])
print(colMeans(output_tnormal_uncal$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output_tnormal_uncal$samples.alpha.ca, Alpha.case, title='A)')
print(mean(output_tnormal_uncal$samples.alpha.ca))

view_tr(output_tnormal_uncal$samples.alpha.co, Alpha.ctrl, title='B)')
print(mean(output_tnormal_uncal$samples.alpha.co))

tag <- "_tnormal_uncal"
output_tnormal_uncal$description <- paste(sim_name, tag, sep="")
save_output(output_tnormal_uncal, paste("output_", sim_name, tag, ".json", sep=""))
save_params_psgp(paste("params_", sim_name, tag, ".json", sep=""), prior='normal')


################################
#### Truncated Normal priors
#### Calibrated initial values
################################


n.sample <- 2000
burnin <- 500
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

w_output <- load_output("w_inival_output_priorcompare.json")
view_logistic_output(w_output)
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

ini_case <- glm(case.data$y ~ case.data$x.standardised + w_i[locs$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data$y ~ ctrl.data$x.standardised + w_i[locs$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4
prior_alpha_ca <- c(prior_alpha_ca_mean, prior_alpha_ca_var)
prior_alpha_co <- c(prior_alpha_co_mean, prior_alpha_co_var)

output_tnormal_uncal <- prefSampleGpCC_truncnorm(data, n.sample, burnin,
                                                 L_w, L_ca, L_co, L_a_ca, L_a_co,
                                                 proposal.sd.theta=proposal.sd.theta,
                                                 m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                                                 target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                                                 self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                                                 delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                                                 beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                                                 theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                                                 prior_phi=prior_phi, prior_theta=prior_theta,
                                                 prior_alpha_ca=prior_alpha_ca, 
                                                 prior_alpha_co=prior_alpha_co)


# optionally burnin the output_normal more
output_tnormal_uncal <- burnin_after(output_tnormal_uncal, n.burn=500)


# optionally continue running if necessary
output_tnormal_uncal <- continueMCMC(data, output_tnormal_uncal, n.sample=2000)


plot(output_tnormal_uncal$deltas_w)
plot(output_tnormal_uncal$deltas_ca)
plot(output_tnormal_uncal$deltas_co)
plot(output_tnormal_uncal$deltas_aca)
plot(output_tnormal_uncal$deltas_aco)
print(output_tnormal_uncal$accept)

plot(apply(output_tnormal_uncal$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_tnormal_uncal$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output_tnormal_uncal$samples.w, w_true=W)

view_tr(output_tnormal_uncal$samples.theta, Theta)
print(mean(output_tnormal_uncal$samples.theta)); print(Theta)

view_tr(output_tnormal_uncal$samples.phi, Phi)
print(mean(output_tnormal_uncal$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
padded_plot(output_tnormal_uncal$samples.beta.ca[,1], beta.case[1])
padded_plot(output_tnormal_uncal$samples.beta.ca[,2], beta.case[2])
padded_plot(output_tnormal_uncal$samples.beta.ca[,3], beta.case[3])
print(colMeans(output_tnormal_uncal$samples.beta.ca))

padded_plot(output_tnormal_uncal$samples.beta.co[,1], beta.ctrl[1])
padded_plot(output_tnormal_uncal$samples.beta.co[,2], beta.ctrl[2])
padded_plot(output_tnormal_uncal$samples.beta.co[,3], beta.ctrl[3])
print(colMeans(output_tnormal_uncal$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output_tnormal_uncal$samples.alpha.ca, Alpha.case, title='A)')
print(mean(output_tnormal_uncal$samples.alpha.ca))

view_tr(output_tnormal_uncal$samples.alpha.co, Alpha.ctrl, title='B)')
print(mean(output_tnormal_uncal$samples.alpha.co))

tag <- "_tnormal_cal"
output_tnormal_uncal$description <- paste(sim_name, tag, sep="")
save_output(output_tnormal_uncal, paste("output_", sim_name, tag, ".json", sep=""))
save_params_psgp(paste("params_", sim_name, tag, ".json", sep=""), prior='normal')

