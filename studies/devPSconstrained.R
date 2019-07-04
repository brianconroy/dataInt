############################
# Fit the spatiotemporal
# shared latent process model
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=9) 


#### Simulate gaussian process 
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)
hist(W)


#### Simulate locational data over time
r <- caPr.disc[[1]]
locs <- simLocW(W, r, beta=0, seed=11) # 42
sum(locs$status)
hist(W)


#### Simulate case and control data
beta.case <- c(1, 0.5, 0.5)
beta.ctrl <- c(2, 0.25, 0.5)
Alpha.case <- 0.75
Alpha.ctrl <- 0.5
cov.disc <- caPr.disc
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W, seed=42)
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40)
print(sum(case.data$y)/sum(case.data$y + ctrl.data$y))


data <- list(
  loc=locs,
  case.data=case.data,
  ctrl.data=ctrl.data
)


# how to perform the constrained update?

# update one
# update the other, propose values > or < than the other
# alternate which you propose first
# want different tuning parameters

#### Specify MCMC parameters
w_initial=W + 2*rnorm(length(W))
beta_ca_initial=rep(0, 3)
beta_co_initial=rep(0, 3)
alpha_ca_initial=0
alpha_co_initial=0
theta_initial=Theta
phi_initial=Phi

m_aca=2000
m_aco=2000
m_ca=700
m_co=700
m_w=700
target_aca=0.75
target_aco=0.75
target_ca=0.75
target_co=0.75
target_w=0.75
self_tune_w=TRUE
self_tune_aca=TRUE
self_tune_aco=TRUE
self_tune_ca=TRUE
self_tune_co=TRUE

L_w=8
L_ca=8
L_co=8
L_a_ca=8
L_a_co=8

prior_alpha_ca_mean=0
prior_alpha_ca_var=4
prior_alpha_co_mean=0
prior_alpha_co_var=4

n.sample=2500
burnin=20

proposal.sd.theta=0.2
prior_theta=c(3,2)
prior_phi <- c(3, 40)


#### Fit model
output <- prefSampleGpCC_constrained(data, d, n.sample, burnin, 
                                     L_w, L_ca, L_co, L_a_ca, L_a_co,
                                     proposal.sd.theta=0.3,
                                     m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, 
                                     target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, 
                                     self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                                     delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                                     beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                                     theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                                     prior_phi, prior_theta, prior_alpha_ca_var, prior_alpha_co_var)

#### Additional burnin
output <- burnin_after_temporal(output, n.burn=50)


#### Generate additional samples
output <- continue_mcmc_temporal(data, D, output, n.sample=1000)


#### View results
par(mfrow=c(2,3))
plot(output$samples.beta.ca[,1], type='l'); abline(h=beta.case[1], col=2)
plot(output$samples.beta.ca[,2], type='l'); abline(h=beta.case[2], col=2)
plot(output$samples.beta.ca[,3], type='l'); abline(h=beta.case[3], col=2)
plot(output$samples.beta.co[,1], type='l'); abline(h=beta.ctrl[1], col=2)
plot(output$samples.beta.co[,2], type='l'); abline(h=beta.ctrl[2], col=2)
plot(output$samples.beta.co[,3], type='l'); abline(h=beta.ctrl[3], col=2)

par(mfrow=c(1,2))
plot(output$samples.alpha.ca, type='l'); abline(h=Alpha.case, col=2)
plot(output$samples.alpha.co, type='l'); abline(h=Alpha.ctrl, col=2)

plot(output$samples.theta, type='l'); abline(h=Theta, col=2)
plot(output$samples.phi, type='l'); abline(h=Phi, col=2)

view_tr_w(output$samples.w, W)
w.hat <- colMeans(output$samples.w)
plot(x=w.hat, y=W); abline(0,1,col=2)

