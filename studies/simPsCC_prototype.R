############################
# Simulation study
# Compare the shared latent
# process model against 
# poisson and spatial poisson
# models in estimating risk
# surface
############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


sampling <- "none"
prevalence <- "low"
sim_name <- gen_sim_name(sampling, prevalence)


#### Load simulation parameters
params <- load_sim_params(sampling, prevalence)
Alpha.case <- params$alpha.case
Alpha.ctrl <- params$alpha.ctrl
beta.case <- as.numeric(strsplit(params$beta.case, split=" ")[[1]])
beta.ctrl <- as.numeric(strsplit(params$beta.ctrl, split=" ")[[1]])


#### Or manually define them
# Alpha.case <- 1
# beta.case <- c(-0.5, 0.25, -0.5)
# Alpha.ctrl <- -1
# beta.ctrl <- c(2, 1, 0.5)

Theta <- 6
Phi <- 12


prior_alpha_ca_mean <- Alpha.case
prior_alpha_ca_var <- 6
prior_alpha_co_mean <- Alpha.ctrl
prior_alpha_co_var <- 6
prior_theta <- c(2.5, 2.5)
prior_phi <- c(9.33, 100)


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


#### Simulate locations
r <- caPr.disc[[1]]
locs <- simLocW(W, r, beta=0, seed=42)
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
n.sample <- 3000
burnin <- 750
L <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

set.seed(241)
beta_ca_i <- abs(rnorm(3))
beta_co_i <- abs(rnorm(3))
alpha_ca_i <- runif(1, 1, 2)
alpha_co_i <- runif(1, -2, -1)
theta_i <- runif(1, 9, 10)
phi_i <- runif(1, 6, 8)
w_i <- rnorm(length(W))

m_aca <- 2000
m_aco <- 2000
m_ca <- 2000
m_co <- 2000
m_w <- 2000

output <- prefSampleGpCC(data, n.sample, burnin,
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=proposal.sd.theta,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                         target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                         delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                         beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                         theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                         prior_phi=prior_phi, prior_theta=prior_theta)


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
alpha_ca_h <- colMeans(output$samples.alpha.ca)
alpha_co_h <- colMeans(output$samples.alpha.co)
phi_h <- colMeans(output$samples.phi)
theta_h <- colMeans(output$samples.theta)

# parameter # true value # estimate # bias
smry <- list()
smry[[1]] <- summarize_param('beta +, 0', beta.case[1], beta_ca_h[1])
smry[[2]] <- summarize_param('beta +, 1', beta.case[2], beta_ca_h[2])
smry[[3]] <- summarize_param('beta -, 0', beta.ctrl[1], beta_co_h[1])
smry[[4]] <- summarize_param('beta -, 1', beta.ctrl[2], beta_co_h[2])
smry[[5]] <- summarize_param('alpha + ', Alpha.case, alpha_ca_h)
smry[[6]] <- summarize_param('alpha - ', Alpha.ctrl, alpha_co_h)
smry[[7]] <- summarize_param('phi', Phi, phi_h)
smry[[8]] <- summarize_param('theta', Theta, theta_h)
df_smry <- ldply(smry, data.frame)


#### Spatial poisson regression
## Cases
X.ca <- case.data$x.standardised
Y.ca <- case.data$y
d.sub <- d[as.logical(locs$status), as.logical(locs$status)]


#### load tuning parameters or set them manually
set.seed(314)
beta_ca_i_ <- beta.case + rnorm(length(beta.case))
w_i_ <- rnorm(nrow(d.sub))
phi_i_ <- Phi + rnorm(1)
theta_i_ <- Theta + rnorm(1)

n.sample_ <- 25000
burnin_ <- 3000
L_w_ <- 8
L_b_ <- 8

prior_phi_ <- c(3, 40)
prior_theta_ <- c(2.5, 2.5)

output.sp_ca <- poissonGp(X.ca, Y.ca, d.sub,
                          n.sample=n.sample_, burnin=burnin_, proposal.sd.theta=0.3,
                          L_w=L_w_, L_b=L_b_,
                          beta_initial=beta_ca_i_, w_initial=w_i_, 
                          phi_initial=phi_i_, theta_initial=theta_i_,
                          prior_phi=prior_phi_, prior_theta=prior_theta_)

print(output.sp_ca$accept)

plot(apply(output.sp_ca$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_ca$samples.w)

view_tr(output.sp_ca$samples.beta[,1])
view_tr(output.sp_ca$samples.beta[,2])
view_tr(output.sp_ca$samples.theta)
view_tr(output.sp_ca$samples.phi)
print(colMeans(output.sp_ca$samples.beta))
print(colMeans(output.sp_ca$samples.theta))
print(colMeans(output.sp_ca$samples.phi))

w.hat_spca <- colMeans(output.sp_ca$samples.w)
beta_ca_sp <- colMeans(output.sp_ca$samples.beta)
kriged_w_ca <- krigeW(output.sp_ca, d, locs$ids)
w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(locs$status))

output.sp_ca$description <- paste("spatial_poisson_case", sampling, prevalence, sep="_")
save_output(output.sp_ca, paste("output.sp_ca_", sim_name, ".json", sep=""))
save_params_psc(paste("params.sp_ca_", sim_name, ".json", sep=""))
save_output(kriged_w_ca, paste("output.krige_ca_", sim_name, ".json", sep=""))

## Controls
X.co <- ctrl.data$x.standardised
Y.co <- ctrl.data$y

n.sample__ <- 25000
burnin__ <- 3000
L_w__ <- 8
L_b__ <- 8

set.seed(314)
beta_co_i__ <- beta.ctrl + rnorm(length(beta.ctrl))
w_i__ <- rnorm(nrow(d.sub))
phi_i__ <- Phi + rnorm(1)
theta_i__ <- Theta + rnorm(1)

prior_phi__ <- c(3, 40)
prior_theta__ <- c(2.5, 2.5)

output.sp_co <- poissonGp(X.co, Y.co, d.sub,
                          n.sample=n.sample__, burnin=burnin__, 
                          L_w=L_w__, L_b=L_b__, proposal.sd.theta=0.3,
                          beta_initial=beta_co_i__, w_initial=w_i__, 
                          phi_initial=phi_i__, theta_initial=theta_i__,
                          prior_phi=prior_phi__, prior_theta=prior_theta__)

print(output.sp_co$accept)

plot(apply(output.sp_co$samples.w, 1, mean), type='l', col='2')
view_tr_w(output.sp_co$samples.w)

view_tr(output.sp_co$samples.beta[,1])
view_tr(output.sp_co$samples.beta[,2])
view_tr(output.sp_co$samples.theta)
view_tr(output.sp_co$samples.phi)
print(colMeans(output.sp_co$samples.beta))
print(colMeans(output.sp_co$samples.theta))
print(colMeans(output.sp_co$samples.phi))

w.hat_spco <- colMeans(output.sp_co$samples.w)
beta_co_sp <- colMeans(output.sp_co$samples.beta)
kriged_w_co <- krigeW(output.sp_co, d, locs$ids)
w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(locs$status))

output.sp_co$description <- paste("spatial_poisson_ctrl", sampling, prevalence, sep="_")
save_output(output.sp_co, paste("output.sp_co_", sim_name, ".json", sep=""))
save_params_psco(paste("params.sp_co_", sim_name, ".json", sep=""))
save_output(kriged_w_co, paste("output.krige_co_", sim_name, ".json", sep=""))

#### Poisson regression
## Cases
rmodel.ca <- glm(Y.ca ~ X.ca-1, family='poisson')
beta_ca_r <- coefficients(rmodel.ca)

## Controls
rmodel.co <- glm(Y.co ~ X.co-1, family='poisson')
beta_co_r <- coefficients(rmodel.co)


save_estimates_pr(beta_ca_r, beta_co_r, paste("estimates_poisson_", sim_name, ".json", sep=""))


####################################
# save and summarize param estimates
# calculate and compare risk surface
####################################


# X <- cov.disc[][!is.na(cov.disc[])]
# X.standard <- matrix((X - mean(X))/sd(X))
# X.standard <- cbind(1, X.standard)

X.standard <- load_x_standard(as.logical(locs$status))
lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
lrisk.true <- lodds.true/(1 - lodds.true)

par(mfrow=c(1, 3))
lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
lrisk.ps <- lodds.ps/(1-lodds.ps)
plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

lodds.sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
lrisk.sp <- lodds.sp/(1-lodds.sp)
plot(x=lodds.true, y=lodds.sp, main='B)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

lodds.r <- X.standard %*% beta_ca_r - X.standard %*% beta_co_r
lrisk.r <- lodds.r/(1-lodds.r)
plot(x=lodds.true, y=lodds.r, main='C)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')

surf.true <- overlay(lodds.true, cov.disc)
surf.ps <- overlay(lodds.ps, cov.disc)
surf.sp <- overlay(lodds.sp, cov.disc)
surf.r <- overlay(lodds.r, cov.disc)
par(mfrow=c(2,2))
pal <- colorRampPalette(c("blue","red"))
brk <- seq(-30, 14, by=4)
plot(surf.true, col=pal(12), breaks=brk, main='A)')
plot(surf.ps, col=pal(12), breaks=brk, main='B)')
plot(surf.sp, col=pal(12), breaks=brk, main='C)')
plot(surf.r, col=pal(12), breaks=brk, main='D)')

rsurf.true <- overlay(lrisk.true, cov.disc)
rsurf.ps <- overlay(lrisk.ps, cov.disc)
rsurf.sp <- overlay(lrisk.sp, cov.disc)
rsurf.r <- overlay(lrisk.r, cov.disc)
par(mfrow=c(2,2))
plot(rsurf.true)
plot(rsurf.ps)
plot(rsurf.sp)
plot(rsurf.r)

par(mfrow=c(1,3))
pal <- colorRampPalette(c("green","red"))
diff.ps <- overlay(lodds.ps - lodds.true, cov.disc)
diff.sp <- overlay(lodds.sp - lodds.true, cov.disc)
diff.r <- overlay(lodds.r - lodds.true, cov.disc)
brks <- seq(-25, 15, by=5)
plot(diff.ps, col=pal(9))
plot(diff.sp, col=pal(9))
plot(diff.r, col=pal(9))

rmse.ps <- round(sqrt(mean((lodds.true - lodds.ps)^2)), 2)
rmse.sp <- round(sqrt(mean((lodds.true - lodds.sp)^2)), 2)
rmse.r <- round(sqrt(mean((lodds.true - lodds.r)^2)), 2)
mae.ps <- round(abs(mean(lodds.true - lodds.ps)), 2)
mae.sp <- round(abs(mean(lodds.true - lodds.sp)), 2)
mae.r <- round(abs(mean(lodds.true - lodds.r)), 2)

metrics <- list(
  list(model='preferential sampling', 'rmse'=rmse.ps, 'mae'=mae.ps),
  list(model='spatial regression', 'rmse'=rmse.sp, 'mae'=mae.sp),
  list(model='poisson regression', 'rmse'=rmse.r, 'mae'=mae.r)
)
df <- ldply(metrics, data.frame)

dir <- "/Users/brianconroy/Documents/research/dataInt/output/"
write.table(lodds.true, paste(dir, 'lodds_true_med.txt', sep=''), row.names=F, col.names=F)
write.table(lodds.ps, paste(dir,'lodds_ps_med.txt', sep=''), row.names=F, col.names=F)
write.table(lodds.sp, paste(dir,'lodds_sp_med.txt', sep=''), row.names=F, col.names=F)
write.table(lodds.r, paste(dir,'lodds_r_med.txt', sep=''), row.names=F, col.names=F)
