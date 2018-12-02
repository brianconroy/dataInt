############################
# Simulation study
# monitoring performance of
# proposed model under 
# different priors 
############################

# plans
  # iterate over phi
  # iterate over theta
  # iterate over alpha
# fix prevalence at medium, medium

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
param <- "alpha"
index <- 1
sim_name <- paste("iterate_prior", param, index, sep="_")


#### Load simulation parameters
params <- load_sim_params(sampling, prevalence)
Alpha.case <- params$alpha.case
Alpha.ctrl <- params$alpha.ctrl
beta.case <- as.numeric(strsplit(params$beta.case, split=" ")[[1]])
beta.ctrl <- as.numeric(strsplit(params$beta.ctrl, split=" ")[[1]])


#### Load priors
priors <- load_priors(param, index)
prior_alpha_ca_var <- priors$prior_alpha
prior_alpha_co_var <- priors$prior_alpha
prior_phi <- priors$prior_phi
prior_theta <- priors$prior_theta


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
prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl


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
n.sample <- 13500
burnin <- 9500
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
w_i <- W + rnorm(length(W))

m_aca <- 4000
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
                         prior_alpha_ca_var, prior_alpha_co_var)


# optionally burnin the output more
output <- burnin_after(output, n.burn=1000)


# optionally continue running if necessary
# output <- continueMCMC(data, output, n.sample=1000)


plot(output$deltas_w)
plot(output$deltas_ca)
plot(output$deltas_co)
plot(output$deltas_aca)
plot(output$deltas_aco)
print(output$accept)

plot(apply(output$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output$samples.w, w_true=W)

view_tr(output$samples.theta, Theta)
print(mean(output$samples.theta)); print(Theta)

view_tr(output$samples.phi, Phi)
print(mean(output$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
view_tr(output$samples.beta.ca[,1], beta.case[1], title='A)')
view_tr(output$samples.beta.ca[,2], beta.case[2], title='B)')
view_tr(output$samples.beta.ca[,3], beta.case[3], title='C)')
print(colMeans(output$samples.beta.ca))

view_tr(output$samples.beta.co[,1], beta.ctrl[1], title='D)')
view_tr(output$samples.beta.co[,2], beta.ctrl[2], title='E)')
view_tr(output$samples.beta.co[,3], beta.ctrl[3], title='F)')
print(colMeans(output$samples.beta.co))

par(mfrow=c(1,2))
view_tr(output$samples.alpha.ca, Alpha.case, title='A)')
print(mean(output$samples.alpha.ca))

view_tr(output$samples.alpha.co, Alpha.ctrl, title='B)')
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
lrisk.true <- lodds.true/(1 - lodds.true)

lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
lrisk.ps <- lodds.ps/(1-lodds.ps)
plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

surf.true <- overlay(lodds.true, cov.disc[[1]])
surf.ps <- overlay(lodds.ps, cov.disc[[1]])
par(mfrow=c(1,2))
pal <- colorRampPalette(c("blue","red"))
brk <- seq(-30, 14, by=4)
plot(surf.true, col=pal(12), breaks=brk, main='A)')
plot(surf.ps, col=pal(12), breaks=brk, main='B)')

rsurf.true <- overlay(lrisk.true, cov.disc[[1]])
rsurf.ps <- overlay(lrisk.ps, cov.disc[[1]])
par(mfrow=c(1,2))
plot(rsurf.true)
plot(rsurf.ps)

par(mfrow=c(1,1))
pal <- colorRampPalette(c("green","red"))
diff.ps <- overlay(lodds.ps - lodds.true, cov.disc[[1]])
brks <- seq(-25, 15, by=5)
plot(diff.ps, col=pal(9))

rmse.ps <- round(sqrt(mean((lodds.true - lodds.ps)^2)), 2)
mae.ps <- round(abs(mean(lodds.true - lodds.ps)), 2)
