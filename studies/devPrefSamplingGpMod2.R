library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
sourceDirectory('Documents/research/dataInt/R/')


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=9)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate gaussian process
Phi <- 12
Theta <- 6
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)
hist(W)


#### Simulate locations
beta.samp <- c(-2, 1, -0.5)
locs <- simBernoulliLocCov(caPr.disc, beta.samp, w=W, seed=56)
sum(locs$status)
hist(locs$w)
plot(caPr.disc[[1]])
points(locs$coords)


#### Simulate counts given locations
Alpha.case <- 1.75
Alpha.ctrl <- 0.25
beta.case <- c(0.25, 1, -0.5)
beta.ctrl <- c(4, 1, -0.5)

case.data <- simConditionalGp2(caPr.disc, locs, beta.case, Alpha.case, W, seed=42)
ctrl.data <- simConditionalGp2(caPr.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40)
data <- list(
  loc=locs,
  case.data=case.data,
  ctrl.data=ctrl.data
)

print(sum(case.data$y)/sum(case.data$y + ctrl.data$y))
print(sum(case.data$y))
print(sum(ctrl.data$y))


#### Fit
m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_loc <- 1000
m_w <- 1000

L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
L_loc <- 8
proposal.sd.theta <- 0.15

target_aca=0.75
target_aco=0.75
target_ca=0.75
target_co=0.75
target_w=0.75
target_loc=0.75

self_tune_w=TRUE
self_tune_aca=TRUE
self_tune_aco=TRUE
self_tune_ca=TRUE
self_tune_co=TRUE
self_tune_loc=TRUE

beta_loc_initial=c(0,0,0)
beta_ca_initial=beta.case
beta_co_initial=beta.ctrl
alpha_ca_initial=Alpha.case
alpha_co_initial=Alpha.ctrl
theta_initial=Theta
phi_initial=Phi
w_initial=W

prior_alpha_ca_mean <- Alpha.case
prior_alpha_ca_var <- 4
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_co_var <- 4
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)

n.sample <- 200
burnin <- 0

# ToDo: estimate initial values

output <- prefSampleGpV2(data, d, n.sample, burnin, 
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=proposal.sd.theta,
                           m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, 
                           target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_loc=target_loc,
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE, self_tune_loc=TRUE,
                           beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial, beta_loc_initial=beta_loc_initial,
                           theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                           prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)

par(mfrow=c(1,3))
plot(output$samples.beta.loc[,1], type='l'); abline(h=beta.samp[1], col=2)
plot(output$samples.beta.loc[,2], type='l'); abline(h=beta.samp[2], col=2)
plot(output$samples.beta.loc[,3], type='l'); abline(h=beta.samp[3], col=2)

par(mfrow=c(1,1))
plot(colMeans(output$samples.w), W); abline(0, 1, col=2)

par(mfrow=c(1,2))
plot(output$samples.alpha.ca, type='l'); abline(h=Alpha.case, col=2)
plot(output$samples.alpha.co, type='l'); abline(h=Alpha.ctrl, col=2)
