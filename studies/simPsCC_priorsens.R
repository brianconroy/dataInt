############################
# Simulation study
# monitoring performance of
# proposed model under 
# different priors 

# fix prevalence at 
# medium, medium
############################

# usage
# specify sim name
# specify sim iteration
# load sim params

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


sampling <- "medium"
prevalence <- "medium"
param <- "alpha" # theta, phi, alpha
index <- 5 # 1, 3, 5
sim_name <- paste("iterate_prior", param, index, sep="_")


#### Load simulation parameters
params <- load_sim_params(sampling, prevalence)
Alpha.case <- params$alpha.case
Alpha.ctrl <- params$alpha.ctrl
beta.case <- as.numeric(strsplit(params$beta.case, split=" ")[[1]])
beta.ctrl <- as.numeric(strsplit(params$beta.ctrl, split=" ")[[1]])


#### Load priors
priors <- load_priors(param, index)
prior_phi <- priors$prior_phi
prior_theta <- priors$prior_theta
prior_alpha_ca_var <- priors$prior_alpha
prior_alpha_co_var <- priors$prior_alpha
prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl

Theta <- 6
Phi <- 12


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


r.w <- caPr.disc[[1]]
r.w[][!is.na(r.w[])] <- W
par(mfrow=c(1,3))
plot(r.w)
plot(caPr.disc[[1]])
plot(caPr.disc[[2]])


#### Simulate locations
r <- caPr.disc[[1]]
locs <- simLocW(W, r, beta=0, seed=11)
sum(locs$status)
hist(W)
plot(r)
points(locs$coords)

rw <- r
rw[][!is.na(rw[])] <- W
plot(rw)
points(locs$coords, pch=16)

#### Simulate counts given locations
cov.disc <- caPr.disc
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W, seed=42)
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40)
print(sum(case.data$y)/sum(case.data$y + ctrl.data$y))


save_true_params(sampling, prevalence)


data <- list(
  loc=locs,
  case.data=case.data,
  ctrl.data=ctrl.data
)


#### Preferential sampling model


#### Load tuning parameters
# tune_params_psgp <- load_params(paste("params_", sim_name, ".json", sep=""))
# n.sample <- tune_params_psgp$n.sample
# burnin <- tune_params_psgp$burnin
# 
# L <- tune_params_psgp$L
# L_ca <- tune_params_psgp$L_ca
# L_co <- tune_params_psgp$L_co
# L_a_ca <- tune_params_psgp$L_a_ca
# L_a_co <- tune_params_psgp$L_a_co
# 
# beta_ca_i <- tune_params_psgp$beta_ca_i
# beta_co_i <- tune_params_psgp$beta_co_i
# alpha_ca_i <- tune_params_psgp$alpha_ca_i
# alpha_co_i <- tune_params_psgp$alpha_co_i
# theta_i <-tune_params_psgp$theta_i
# phi_i <- tune_params_psgp$phi_i
# w_i <- tune_params_psgp$w_i

#### Or manually define them
n.sample <- 2000
burnin <- 500
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

#### Initial values


# W initial value
# w_output <- logisticGp(y=locs$status, d, n.sample=1000, burnin=200, L=10,
#                        prior_phi=prior_phi, prior_theta=prior_theta)
# view_logistic_output(w_output)
# save_output(w_output, "w_inival_output_priorsens.json")
w_output <- load_output("w_inival_output_priorsens.json")
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

# Beta & alpha initial values
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

output <- prefSampleGpCC(data, n.sample, burnin,
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
output <- burnin_after(output, n.burn=1000)


# optionally continue running if necessary
output <- continueMCMC(data, output, n.sample=2000)


plot(output$deltas_w)
plot(output$deltas_ca)
plot(output$deltas_co)
plot(output$deltas_aca)
plot(output$deltas_aco)
print(output$accept)

padded_plot(apply(output$samples.w, 1, mean), mean(W))
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output$samples.w, w_true=W)

view_tr(output$samples.theta, Theta)
print(mean(output$samples.theta)); print(Theta)

view_tr(output$samples.phi, Phi)
print(mean(output$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
padded_plot(output$samples.beta.ca[,1], beta.case[1])
padded_plot(output$samples.beta.ca[,2], beta.case[2])
padded_plot(output$samples.beta.ca[,3], beta.case[3])
print(colMeans(output$samples.beta.ca))

padded_plot(output$samples.beta.co[,1], beta.ctrl[1])
padded_plot(output$samples.beta.co[,2], beta.ctrl[2])
padded_plot(output$samples.beta.co[,3], beta.ctrl[3])
print(colMeans(output$samples.beta.co))

par(mfrow=c(1,2))
padded_plot(output$samples.alpha.ca, Alpha.case)
print(mean(output$samples.alpha.ca))

padded_plot(output$samples.alpha.co, Alpha.ctrl)
print(mean(output$samples.alpha.co))

output$description <- sim_name
save_output(output, paste("output_", sim_name, ".json", sep=""))
save_params_psgp(paste("params_", sim_name, ".json", sep=""))

w.hat <- colMeans(output$samples.w)
beta_ca_h <- colMeans(output$samples.beta.ca)
beta_co_h <- colMeans(output$samples.beta.co)
alpha_ca_h <- mean(output$samples.alpha.ca)
alpha_co_h <- mean(output$samples.alpha.co)
phi_h <- mean(output$samples.phi)
theta_h <- mean(output$samples.theta)


####################################
# save and summarize param estimates
# calculate and compare risk surface
####################################


X.standard <- load_x_standard(as.logical(locs$status))
lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W

lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

