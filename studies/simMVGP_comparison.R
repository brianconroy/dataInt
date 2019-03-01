#######################################
# Compare the multivariate gaussian
# process preferential sampling (PS) model
# against separate analyses of 
# the single species PS model. 
#######################################

# fix at 2 species
# iterate over? differing levels of alpha? done that
# iterate over? low, medium, high correlation? Do this, for now
# no iteration? just comparison over multiple realizations of data?
# both? low, high correlation, within each level generate 3 different datasets

library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


sim <- "simMVGP_comparison"
level <- "high"
sim_name <- paste(sim, level, sep="_")


#### Load simulation parameters
params <- load_output(paste('simMVGP_comparison_params_', level, '.json', sep=''))
Alpha.case1 <- params$Alpha.case1
Alpha.case2 <- params$Alpha.case2
Alpha.ctrl1 <- params$Alpha.ctrl1
Alpha.ctrl2 <- params$Alpha.ctrl2
beta.case1 <- params$beta.case1
beta.case2 <- params$beta.case2
beta.ctrl1 <- params$beta.ctrl1
beta.ctrl2 <- params$beta.ctrl2
Theta <- params$Theta
Tmat <- params$Tmat
W <- params$W


# #### Write simulation parameters to LaTeX
# sim_config <- list(
#   list(parameter="Alpha (case) 1", value=as.character(Alpha.case1)),
#   list(parameter="Alpha (case) 2", value=as.character(Alpha.case2)),
#   list(parameter="Alpha (control) 1",value=as.character(Alpha.ctrl1)),
#   list(parameter="Alpha (control) 2",value=as.character(Alpha.ctrl2)),
#   list(parameter="Beta (case) 1", value=paste(beta.case1, collapse=", ")),
#   list(parameter="Beta (case) 2", value=paste(beta.case2, collapse=", ")),
#   list(parameter="Beta (control) 1", value=paste(beta.ctrl1, collapse=", ")),
#   list(parameter="Beta (control) 2", value=paste(beta.ctrl2, collapse=", ")),
#   list(parameter="Range", value=as.character(Theta)),
#   list(parameter="Phi", value=as.character(Phi))
# )
# write_latex_table(ldply(sim_config, "data.frame"), fname=paste("sim_params_", level, ".txt", sep=""), path="/Users/brianconroy/Documents/research/project2/simulation_1_comparison")


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


#### Load Data
data <- load_output(paste('simMVGP_comparison_data_', level, '.json', sep=''))
locs1 <- list(
  status=data$locs$status[[1]],
  cells=data$locs$cells[[1]],
  coords=data$locs$coords[[1]],
  ids=data$locs$ids[[1]]
)
locs2 <- list(
  status=data$locs$status[[2]],
  cells=data$locs$cells[[2]],
  coords=data$locs$coords[[2]],
  ids=data$locs$ids[[2]]
)
case.data1 <- list(
  y=data$case.data$y[[1]],
  x.standardised=data$case.data$x.standardised[[1]],
  x=data$case.data$x[[1]],
  p=data$case.data$p[[1]]
)
case.data2 <- list(
  y=data$case.data$y[[2]],
  x.standardised=data$case.data$x.standardised[[2]],
  x=data$case.data$x[[2]],
  p=data$case.data$p[[2]]
)
ctrl.data1 <- list(
  y=data$ctrl.data$y[[1]],
  x.standardised=data$ctrl.data$x.standardised[[1]],
  x=data$ctrl.data$x[[1]],
  p=data$ctrl.data$p[[1]]
)
ctrl.data2 <- list(
  y=data$ctrl.data$y[[2]],
  x.standardised=data$ctrl.data$x.standardised[[2]],
  x=data$ctrl.data$x[[2]],
  p=data$ctrl.data$p[[2]]
)
data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)

print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))

r <- caPr.disc[[1]]
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r)
points(locs1$coords)
plot(r)
points(locs2$coords)
locs=list(locs1, locs2)

#############
#### Fit MVGP
#############

#### Set Priors
prior_theta <- c(3, 2)
prior_alpha_ca_mean <- c(Alpha.case1, Alpha.case2)
prior_alpha_co_mean <- c(Alpha.ctrl1, Alpha.ctrl2)
prior_alpha_ca_var <- c(4, 4)
prior_alpha_co_var <- c(4, 4)
Omega <- matrix(c(5, 0, 0, 5), nrow=2)
r <- 4
print(Omega/(r - ncol(Tmat) - 1))

#### Fit initial values
y <- list(locs1$status, locs2$status)
prior_t=list(scale=matrix(c(5, 0, 0, 5), nrow=2), df=4)
w_output <- logisticMVGP(y, D, n.sample=1000, burnin=200, L=10, 
                         prior_t=prior_t, prior_theta=prior_theta)
w.hat <- colMeans(w_output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
save_output(w_output, paste("w_inival_output_mvgp_", level, ".json", sep=""))

w_initial=colMeans(w_output$samples.w)
theta_initial <- mean(w_output$samples.theta)
t_initial <- matrix(colMeans(w_output$samples.t), nrow=2)

alpha_ca_initial <- list()
alpha_co_initial <- list()
beta_ca_initial <- list()
beta_co_initial <- list()
N.d <- length(data$case.data)
for (k in 1:N.d){
  k_seq <- seq(k, length(w_initial), by=N.d)
  w_k <- w_initial[k_seq]
  
  ini_case <- glm(data$case.data[[k]]$y ~ data$case.data[[k]]$x.standardised + w_k[data$locs[[k]]$ids] - 1, family='poisson')
  alpha_ca_initial[[k]] <- unname(coefficients(ini_case)[4])
  beta_ca_initial[[k]] <- unname(coefficients(ini_case)[1:3])
  
  ini_ctrl <- glm(data$ctrl.data[[k]]$y ~ data$ctrl.data[[k]]$x.standardised + w_k[data$locs[[k]]$ids] - 1, family='poisson')
  alpha_co_initial[[k]] <- unname(coefficients(ini_ctrl)[4])
  beta_co_initial[[k]] <- unname(coefficients(ini_ctrl)[1:3])
}

n.sample <- 2500
burnin <- 500
L_w <- 8
L_ca <- c(8, 8)
L_co <- c(8, 8)
L_a_ca <- c(8, 8)
L_a_co <- c(8, 8)
proposal.sd.theta <- 0.10

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

self_tune_w=TRUE
self_tune_aca=TRUE
self_tune_aco=TRUE
self_tune_ca=TRUE
self_tune_co=TRUE

target_aca <- 0.65
target_aco <- 0.65
target_ca <- 0.65
target_co <- 0.65
target_w <- 0.65

output <- prefSampleMVGP(data, D, n.sample, burnin,
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=0.3,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                         target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, 
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                         delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                         beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                         theta_initial=theta_initial, t_initial=t_initial, w_initial=w_initial,
                         prior_phi, prior_theta, prior_alpha_ca_mean, prior_alpha_co_mean, prior_alpha_ca_var, prior_alpha_co_var,
                         prior_t)

# optional: continue running the Markov chain
output <- continueMCMC_mvgp(data, D, output, n.sample=500)

# optional: additional burnin
output <- burnin_mvgp(output, n.burn=500)

print(colMeans(output$samples.t))

plot(output$samples.theta, type='l'); abline(h=Theta, col=2)
plot(output$samples.alpha.ca[1,,], type='l', ylim=c(0,1)); abline(h=Alpha.case1, col=2)
plot(output$samples.alpha.ca[2,,], type='l', ylim=c(0.5,2)); abline(h=Alpha.case2, col=2)
plot(output$samples.alpha.co[1,,], type='l', ylim=c(-1, -0.4)); abline(h=Alpha.ctrl1, col=2)
plot(output$samples.alpha.co[2,,], type='l', ylim=c(-1.5, -0.75)); abline(h=Alpha.ctrl2, col=2)
plot(output$samples.beta.ca[1,,1], type='l'); abline(h=beta.case1[1], col=2)
plot(output$samples.beta.ca[1,,2], type='l'); abline(h=beta.case1[2], col=2)
plot(output$samples.beta.ca[1,,3], type='l'); abline(h=beta.case1[3], col=2)
plot(output$samples.beta.ca[2,,1], type='l'); abline(h=beta.case2[1], col=2)
plot(output$samples.beta.ca[2,,2], type='l'); abline(h=beta.case2[2], col=2)
plot(output$samples.beta.ca[2,,3], type='l'); abline(h=beta.case2[3], col=2)

par(mfrow=c(2,3))
plot(output$samples.beta.co[1,,1], type='l'); abline(h=beta.ctrl1[1], col=2)
plot(output$samples.beta.co[1,,2], type='l'); abline(h=beta.ctrl1[2], col=2)
plot(output$samples.beta.co[1,,3], type='l'); abline(h=beta.ctrl1[3], col=2)
plot(output$samples.beta.co[2,,1], type='l'); abline(h=beta.ctrl2[1], col=2)
plot(output$samples.beta.co[2,,2], type='l'); abline(h=beta.ctrl2[2], col=2)
plot(output$samples.beta.co[2,,3], type='l'); abline(h=beta.ctrl2[3], col=2)

print(output$accept)
par(mfrow=c(1,1))
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
w.hat1 <- w.hat[seq(1,ncol(D),by=2)]
w.hat2 <- w.hat[seq(2,ncol(D),by=2)]
plot(x=W[seq(1,ncol(D),by=2)], y=w.hat1); abline(0, 1, col=2)
plot(x=W[seq(2,ncol(D),by=2)], y=w.hat2); abline(0, 1, col=2)
view_tr_w(output$samples.w, w_true=W)
plot(apply(output$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')

plot(output$samples.t[,1], type='l'); abline(h=Tmat[1], col=2)
plot(output$samples.t[,2], type='l'); abline(h=Tmat[2], col=2)
plot(output$samples.t[,4], type='l'); abline(h=Tmat[4], col=2)

tag <- paste(sim_name, 'mvgp', sep="_")
output$description <- tag
save_output(output, paste("output_", tag, ".json", sep=""))


#################
# Separate Models
#################


###############
# First Species
###############


data1 <- list(
  case.data=case.data1,
  ctrl.data=ctrl.data1,
  locs=locs1
)


# calibrated initial values
prior_phi <- c(18, 204)
w_output_s1 <- logisticGp(y=locs1$status, D, n.sample=1000, burnin=200, L=10,
                          prior_phi=prior_phi, prior_theta=prior_theta)
plot(x=colMeans(w_output_s1$samples.w), y=W[seq(1, length(W), by=2)]); abline(0, 1, col=2)
view_tr(w_output_s1$samples.theta, Theta)
save_output(w_output_s1, paste("w_inival_output_comparison_species1_", level, ".json", sep=""))

# w_output <- load_output("w_inival_output_multi.json")
w_i <- colMeans(w_output_s1$samples.w)
theta_i <- mean(w_output_s1$samples.theta)
phi_i <- mean(w_output_s1$samples.phi)

ini_case <- glm(case.data1$y ~ case.data1$x.standardised + w_i[locs1$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data1$y ~ ctrl.data1$x.standardised + w_i[locs1$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

prior_alpha_ca_mean <- Alpha.case1
prior_alpha_co_mean <- Alpha.ctrl1
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4

n.sample <- 2500
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

output_s1 <- prefSampleGpCC(data1, D, n.sample, burnin,
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
output_s1 <- burnin_after(output_s1, n.burn=500)


# optionally continue running if necessary
output_s1 <- continueMCMC(data1, output_s1, n.sample=1000)


plot(output_s1$deltas_w)
plot(output_s1$deltas_ca)
plot(output_s1$deltas_co)
plot(output_s1$deltas_aca)
plot(output_s1$deltas_aco)
print(output_s1$accept)

plot(apply(output_s1$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_s1$samples.w)
plot(x=W1, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output_s1$samples.w, w_true=W)

view_tr(output_s1$samples.theta, Theta)
print(mean(output_s1$samples.theta)); print(Theta)

view_tr(output_s1$samples.phi, Phi)
print(mean(output_s1$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
view_tr(output_s1$samples.beta.ca[,1], beta.case1[1])
view_tr(output_s1$samples.beta.ca[,2], beta.case1[2])
view_tr(output_s1$samples.beta.ca[,3], beta.case1[3])
print(colMeans(output_s1$samples.beta.ca))

view_tr(output_s1$samples.beta.co[,1], beta.ctrl1[1])
view_tr(output_s1$samples.beta.co[,2], beta.ctrl1[2])
view_tr(output_s1$samples.beta.co[,3], beta.ctrl1[3])
print(colMeans(output_s1$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output_s1$samples.alpha.ca, Alpha.case1)
print(mean(output_s1$samples.alpha.ca))

view_tr(output_s1$samples.alpha.co, Alpha.ctrl1)
print(mean(output_s1$samples.alpha.co))

tag <- paste(sim_name, 'species1', sep="_")
output_s1$description <- tag
save_output(output_s1, paste("output_", tag, ".json", sep=""))


#################
## Second Species
#################


data2 <- list(
  case.data=case.data2,
  ctrl.data=ctrl.data2,
  locs=locs2
)

prior_phi <- c(18, 204)
w_output_s2 <- logisticGp(y=locs2$status, D, n.sample=1000, burnin=200, L=10,
                          prior_phi=prior_phi, prior_theta=prior_theta)
plot(x=colMeans(w_output_s2$samples.w), y=W[seq(2, length(W), by=2)]); abline(0, 1, col=2)
view_tr(w_output_s2$samples.theta, Theta)
save_output(w_output_s2, paste("w_inival_output_comparison_species2_", level, ".json", sep=""))

ini_case <- glm(case.data2$y ~ case.data2$x.standardised + w_i[locs2$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data2$y ~ ctrl.data2$x.standardised + w_i[locs2$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

prior_alpha_ca_mean <- Alpha.case2
prior_alpha_co_mean <- Alpha.ctrl2
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

output_s2 <- prefSampleGpCC(data2, D, n.sample, burnin,
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
output_s2 <- burnin_after(output_s2, n.burn=500)


# optionally continue running if necessary
output_s2 <- continueMCMC(data2, output_s2, n.sample=1000)


plot(output_s2$deltas_w)
plot(output_s2$deltas_ca)
plot(output_s2$deltas_co)
plot(output_s2$deltas_aca)
plot(output_s2$deltas_aco)
print(output_s2$accept)

plot(apply(output_s2$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_s2$samples.w)
plot(x=W2, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output_s2$samples.w, w_true=W)

view_tr(output_s2$samples.theta, Theta)
print(mean(output_s2$samples.theta)); print(Theta)

view_tr(output_s2$samples.phi, Phi)
print(mean(output_s2$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
view_tr(output_s2$samples.beta.ca[,1], beta.case2[1])
view_tr(output_s2$samples.beta.ca[,2], beta.case2[2])
view_tr(output_s2$samples.beta.ca[,3], beta.case2[3])
print(colMeans(output_s2$samples.beta.ca))

view_tr(output_s2$samples.beta.co[,1], beta.ctrl2[1])
view_tr(output_s2$samples.beta.co[,2], beta.ctrl2[2])
view_tr(output_s2$samples.beta.co[,3], beta.ctrl2[3])
print(colMeans(output_s2$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output_s2$samples.alpha.ca, Alpha.case2, title='A)')
print(mean(output_s2$samples.alpha.ca))

view_tr(output_s2$samples.alpha.co, Alpha.ctrl2, title='B)')
print(mean(output_s2$samples.alpha.co))

tag <- paste(sim_name, 'species2', sep="_")
output_s2$description <- tag
save_output(output_s2, paste("output_", tag, ".json", sep=""))
