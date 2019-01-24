#######################
# test the preferential sampling
# model with an offset
######################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)
offset_ <- mean(area(caPr.disc[[1]])[])
offset_ <- offset_/1000 # Megameters
# offset_ <- 1 # no offset


#### Simulation parameters
Alpha.case <- 1
Alpha.ctrl <- -1
beta.case <- c(-0.25,  0.75, -0.50)
beta.ctrl <- c(3.0, 1.0, 0.5)
Theta <- 6
Phi <- 12


#### Priors
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)


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
rw <- r
rw[][!is.na(rw[])] <- W
pal <- colorRampPalette(c("blue","red"))
plot(rw, col=pal(10))
points(locs$coords, pch=16)


#### Simulate counts given locations
case.data <- simConditionalGp2(caPr.disc, locs, beta.case, Alpha.case, W, seed=42, offset=offset_)
ctrl.data <- simConditionalGp2(caPr.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40, offset=offset_)
print(sum(case.data$y))
print(sum(case.data$y)/sum(case.data$y + ctrl.data$y))


data <- list(
  loc=locs,
  case.data=case.data,
  ctrl.data=ctrl.data
)


##############
#### Fit Model
##############


w_output <- load_output("w_inival_output_priorcompare.json")
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

ini_case <- glm(case.data$y ~ case.data$x.standardised + w_i[locs$ids] + offset(rep(log(offset_), length(case.data$y))) - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data$y ~ ctrl.data$x.standardised + w_i[locs$ids] + offset(rep(log(offset_), length(case.data$y))) - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]


#### Tuning parameters
n.sample <- 1000
burnin <- 0
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

prior_alpha_ca_mean <- Alpha.case
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4

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
                         prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var, offset=offset_)

# optionally burnin the output more
output <- burnin_after(output, n.burn=6000)


# optionally continue running if necessary
output <- continueMCMC(data, output, n.sample=2000)


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
view_tr(output$samples.beta.ca[,1], beta.case[1])
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
