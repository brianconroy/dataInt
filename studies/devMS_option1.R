############################
# Multispecies simulation:
# multiple preferentially sampled
# rodent species. 
#   2 to begin
############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Simulation parameters
Alpha.case1 <- 0.5
Alpha.ctrl1 <- -0.5
Alpha.case2 <- 1
Alpha.ctrl2 <- -1
beta.case1 <- c(0.0, 0.75, -0.50)
beta.ctrl1 <- c(2.75, 1.00, 0.50)
beta.case2 <- c(-0.25, 0.8, -1)
beta.ctrl2 <- c(1, 1.25, 1)


Theta <- 6
Phi <- 12


prior_theta <- c(6, 1)
prior_phi <- c(18, 204)
prior_alpha_ca_var <- 4
prior_alpha_co_var <- 4
prior_alpha_ca_mean1 <- Alpha.case1
prior_alpha_co_mean1 <- Alpha.ctrl1
prior_alpha_ca_mean2 <- Alpha.case2
prior_alpha_co_mean2 <- Alpha.ctrl2


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
locs1 <- simLocW(W, r, beta=0, seed=11) # 42
locs2 <- simLocW(W, r, beta=0, seed=42)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r)
points(locs1$coords)
plot(r)
points(locs2$coords)


#### Simulate counts given locations
cov.disc <- caPr.disc
case.data1 <- simConditionalGp2(cov.disc, locs1, beta.case1, Alpha.case1, W, seed=42)
ctrl.data1 <- simConditionalGp2(cov.disc, locs1, beta.ctrl1, Alpha.ctrl1, W, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))

case.data2 <- simConditionalGp2(cov.disc, locs2, beta.case2, Alpha.case2, W, seed=42)
ctrl.data2 <- simConditionalGp2(cov.disc, locs2, beta.ctrl2, Alpha.ctrl2, W, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))


data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)


# calibrated initial values
w_output <- logisticGp(y=locs1$status, d, n.sample=1000, burnin=200, L=10,
                      prior_phi=prior_phi, prior_theta=prior_theta)
view_logistic_output(w_output)

save_output(w_output, "w_inival_output_priorcompare.json")

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
L_ca <- c(8, 8)
L_co <- c(8, 8)
L_a_ca <- c(8, 8)
L_a_co <- c(8, 8)
proposal.sd.theta <- 0.15

m_aca <- c(1000, 1000)
m_aco <- c(1000, 1000)
m_ca <- c(1000, 1000)
m_co <- c(1000, 1000)
m_w <- 1000

target_aca <- c(0.65, 0.65)
target_aco <- c(0.65, 0.65)
target_ca <- c(0.65, 0.65)
target_co <- c(0.65, 0.65)
target_w <- 0.65

prior_alpha_ca_var <- c(Alpha.case1, Alpha.case2)
prior_alpha_ca_var <- c(Alpha.case1, Alpha.case2)


output <- prefSampleMulti_1(data, n.sample, burnin, 
                            L_w, L_ca, L_co, L_a_ca, L_a_co,
                            proposal.sd.theta=0.3,
                            m_aca=NULL, m_aco=NULL, m_ca=NULL, m_co=NULL, m_w=NULL, 
                            target_aca=NULL, target_aco=NULL, target_ca=NULL, target_co=NULL, target_w=NULL, 
                            self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                            delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                            beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                            theta_initial=NULL, phi_initial=NULL, w_initial=NULL,
                            prior_phi, prior_theta, prior_alpha_ca_var, prior_alpha_co_var)
