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
sourceDirectory('Documents/research/dataInt/R/')


sim <- "simMVGP_comparison"
level <- "low"
sim_name <- paste(sim, '_', level, sep="_")


#### Load simulation parameters
all_params <- load_sim_params_general('simMulti_params_opt1_comparison.txt')
params <- all_params[all_params$perturb == level,]
Alpha.case1 <- params$alpha.case1
Alpha.case2 <- params$alpha.case2
Alpha.ctrl1 <- params$alpha.ctrl1
Alpha.ctrl2 <- params$alpha.ctrl2
beta.case1 <- as.numeric(strsplit(params$beta.case1, split=" ")[[1]])
beta.case2 <- as.numeric(strsplit(params$beta.case2, split=" ")[[1]])
beta.ctrl1 <- as.numeric(strsplit(params$beta.ctrl1, split=" ")[[1]])
beta.ctrl2 <- as.numeric(strsplit(params$beta.ctrl2, split=" ")[[1]])
Theta <- params$Theta
Phi <- params$Phi


#### Write simulation parameters to LaTeX
sim_config <- list(
  list(parameter="Alpha (case) 1", value=as.character(Alpha.case1)),
  list(parameter="Alpha (case) 2", value=as.character(Alpha.case2)),
  list(parameter="Alpha (control) 1",value=as.character(Alpha.ctrl1)),
  list(parameter="Alpha (control) 2",value=as.character(Alpha.ctrl2)),
  list(parameter="Beta (case) 1", value=paste(beta.case1, collapse=", ")),
  list(parameter="Beta (case) 2", value=paste(beta.case2, collapse=", ")),
  list(parameter="Beta (control) 1", value=paste(beta.ctrl1, collapse=", ")),
  list(parameter="Beta (control) 2", value=paste(beta.ctrl2, collapse=", ")),
  list(parameter="Range", value=as.character(Theta)),
  list(parameter="Phi", value=as.character(Phi))
)
write_latex_table(ldply(sim_config, "data.frame"), fname=paste("sim_params_", level, ".txt", sep=""), path="/Users/brianconroy/Documents/research/project2/simulation_1_comparison")


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
W1 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W1)
hist(W1)

if (w_difference == 'different'){
  set.seed(44)
  W2 <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
} else {
  W2 <- W1
}

#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W1, r, beta=0, seed=11)
locs2 <- simLocW(W2, r, beta=0, seed=42)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r)
points(locs1$coords)
plot(r)
points(locs2$coords)


#### Simulate counts given locations
cov.disc <- caPr.disc
case.data1 <- simConditionalGp2(cov.disc, locs1, beta.case1, Alpha.case1, W1, seed=42)
ctrl.data1 <- simConditionalGp2(cov.disc, locs1, beta.ctrl1, Alpha.ctrl1, W1, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))

case.data2 <- simConditionalGp2(cov.disc, locs2, beta.case2, Alpha.case2, W2, seed=42)
ctrl.data2 <- simConditionalGp2(cov.disc, locs2, beta.ctrl2, Alpha.ctrl2, W2, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))


data <- list(
  locs=list(locs1, locs2),
  case.data=list(case.data1, case.data2),
  ctrl.data=list(ctrl.data1, ctrl.data2)
)

save_true_params_multi(sim_name)

####################
# Multispecies Model
####################

# calibrated initial values
# w_output <- logisticGp(y=locs1$status, d, n.sample=1000, burnin=200, L=10,
#                       prior_phi=prior_phi, prior_theta=prior_theta)
# view_logistic_output(w_output)
# save_output(w_output, "w_inival_output_multi.json")

w_output <- load_output("w_inival_output_multi.json")
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

alpha_ca_i <- list()
alpha_co_i <- list()
beta_ca_i <- list()
beta_co_i <- list()
N.d <- length(data$case.data)
for (k in 1:N.d){
  ini_case <- glm(data$case.data[[k]]$y ~ data$case.data[[k]]$x.standardised + w_i[data$locs[[k]]$ids] - 1, family='poisson')
  alpha_ca_i[[k]] <- unname(coefficients(ini_case)[4])
  beta_ca_i[[k]] <- unname(coefficients(ini_case)[1:3])
  
  ini_ctrl <- glm(data$ctrl.data[[k]]$y ~ data$ctrl.data[[k]]$x.standardised + w_i[data$locs[[k]]$ids] - 1, family='poisson')
  alpha_co_i[[k]] <- unname(coefficients(ini_ctrl)[4])
  beta_co_i[[k]] <- unname(coefficients(ini_ctrl)[1:3])
}

prior_alpha_ca_mean <- c(Alpha.case1, Alpha.case2)
prior_alpha_co_mean <- c(Alpha.ctrl1, Alpha.ctrl2)
prior_alpha_ca_var <- c(4, 4)
prior_alpha_co_var <- c(4, 4)
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)

n.sample <- 2000
burnin <- 500
L_w <- 8
L_ca <- c(8, 8)
L_co <- c(8, 8)
L_a_ca <- c(8, 8)
L_a_co <- c(8, 8)
proposal.sd.theta <- 0.15

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

target_aca <- 0.65
target_aco <- 0.65
target_ca <- 0.65
target_co <- 0.65
target_w <- 0.65

output <- prefSampleMulti_1(data, n.sample, burnin, 
                            L_w, L_ca, L_co, L_a_ca, L_a_co,
                            proposal.sd.theta=0.3,
                            m_aca=m_aca, m_aco=m_aca, m_ca=m_aca, m_co=m_aca, m_w=m_aca, 
                            target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, 
                            self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                            delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                            beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                            theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                            prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_mean=prior_alpha_ca_mean,
                            prior_alpha_co_mean=prior_alpha_co_mean,
                            prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)


# output
output <- burnin_after(output, n.burn=500)


# continue mcmc
output <- continueMCMC_multi(data, output, n.sample=2000)


plot(apply(output$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W1), col='2')
w.hat <- colMeans(output$samples.w)
plot(x=W1, y=w.hat); abline(0, 1, col=2)
summary(100*(W1-w.hat)/W1)
view_tr_w(output$samples.w, w_true=W)

view_tr(output$samples.theta, Theta)
print(mean(output$samples.theta)); print(Theta)

view_tr(output$samples.phi, Phi)
print(mean(output$samples.phi)); print(Phi)

par(mfrow=c(2,2))
padded_plot(output$samples.alpha.ca[1,,], Alpha.case1)
padded_plot(output$samples.alpha.ca[2,,], Alpha.case2)

padded_plot(output$samples.alpha.co[1,,], Alpha.ctrl1)
padded_plot(output$samples.alpha.co[2,,], Alpha.ctrl2)

par(mfrow=c(2,3))
padded_plot(output$samples.beta.ca[1,,1], beta.case1[1])
padded_plot(output$samples.beta.ca[1,,2], beta.case1[2])
padded_plot(output$samples.beta.ca[1,,3], beta.case1[3])

padded_plot(output$samples.beta.co[1,,1], beta.ctrl1[1])
padded_plot(output$samples.beta.co[1,,2], beta.ctrl1[2])
padded_plot(output$samples.beta.co[1,,3], beta.ctrl1[3])

padded_plot(output$samples.beta.ca[2,,1], beta.case2[1])
padded_plot(output$samples.beta.ca[2,,2], beta.case2[2])
padded_plot(output$samples.beta.ca[2,,3], beta.case2[3])

padded_plot(output$samples.beta.co[2,,1], beta.ctrl2[1])
padded_plot(output$samples.beta.co[2,,2], beta.ctrl2[2])
padded_plot(output$samples.beta.co[2,,3], beta.ctrl2[3])

tag <- "_multi"
output$description <- paste(sim_name, tag, sep="")
save_output(output, paste("output_", sim_name, tag, ".json", sep=""))

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

output_s1 <- prefSampleGpCC(data1, n.sample, burnin,
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

tag <- "_separate1"
output_s1$description <- paste(sim_name, tag, sep="")
save_output(output_s1, paste("output_", sim_name, tag, ".json", sep=""))


#################
## Second Species
#################

data2 <- list(
  case.data=case.data2,
  ctrl.data=ctrl.data2,
  locs=locs2
)

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

output_s2 <- prefSampleGpCC(data2, n.sample, burnin,
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
output_s2 <- burnin_after(output_s2, n.burn=700)


# optionally continue running if necessary
output_s2 <- continueMCMC(data2, output_s2, n.sample=2000)


plot(output_s2$deltas_w)
plot(output_s2$deltas_ca)
plot(output_s2$deltas_co)
plot(output_s2$deltas_aca)
plot(output_s2$deltas_aco)
print(output_s2$accept)

plot(apply(output_s2$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_s2$samples.w)
plot(x=W1, y=w.hat); abline(0, 1, col=2)
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
view_tr(output_s2$samples.alpha.ca, Alpha.case1, title='A)')
print(mean(output_s2$samples.alpha.ca))

view_tr(output_s2$samples.alpha.co, Alpha.ctrl1, title='B)')
print(mean(output_s2$samples.alpha.co))

tag <- "_separate_species2"
output_s2$description <- paste(sim_name, tag, sep="")
save_output(output_s2, paste("output_", sim_name, tag, ".json", sep=""))


##################
# Unispecies Model
# With Pooled Data
##################

## Pool data
status <- as.numeric(locs1$status | locs2$status)
cells <- sort(unique(c(locs1$cells, locs2$cells)))
ids <- sort(unique(c(locs1$ids, locs2$ids)))
locs_pooled <- list(status=status, cells=cells, ids=ids)

df1 <- data.frame(cbind(locs1$ids, case.data1$y))
df2 <- data.frame(cbind(locs2$ids, case.data2$y))
df <- merge(df1, df2, by='X1', all=T)
df[is.na(df[,2]),2] <- 0
df[is.na(df[,3]),3] <- 0
y.ca <- df$X2.x + df$X2.y

df1 <- data.frame(cbind(locs1$ids, ctrl.data1$y))
df2 <- data.frame(cbind(locs2$ids, ctrl.data2$y))
df <- merge(df1, df2, by='X1', all=T)
df[is.na(df[,2]),2] <- 0
df[is.na(df[,3]),3] <- 0
y.co <- df$X2.x + df$X2.y

df1 <- data.frame(cbind(locs1$ids, case.data1$x.standardised))
df2 <- data.frame(cbind(locs2$ids, case.data2$x.standardised))
df <- merge(df1, df2, by='V1', all=T)
for (i in 2:7){
  df[is.na(df[,i]),i] <- -Inf
}
x.standardised <- array(NA, c(nrow(df), 3))
for (i in 1:nrow(x.standardised)){
  x.standardised[i,1] <- max(df[i,2], df[i,5])
  x.standardised[i,2] <- max(df[i,3], df[i,6])
  x.standardised[i,3] <- max(df[i,4], df[i,7])
}

case_pooled <- list(x.standardised=x.standardised, y=y.ca)
ctrl_pooled <- list(x.standardised=x.standardised, y=y.co)
  
data_pooled <- list(
  case.data=case_pooled,
  ctrl.data=ctrl_pooled,
  locs=locs_pooled
)

ini_case <- glm(case_pooled$y ~ case_pooled$x.standardised + w_i[locs_pooled$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl_pooled$y ~ ctrl_pooled$x.standardised + w_i[locs_pooled$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

prior_alpha_ca_mean <- mean(Alpha.case1, Alpha.case2)
prior_alpha_co_mean <- mean(Alpha.ctrl1, Alpha.ctrl2)
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

output_pooled <- prefSampleGpCC(data_pooled, n.sample, burnin,
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
output_pooled <- burnin_after(output_pooled, n.burn=500)


# optionally continue running if necessary
output_pooled <- continueMCMC(data2, output_pooled, n.sample=1000)


plot(output_pooled$deltas_w)
plot(output_pooled$deltas_ca)
plot(output_pooled$deltas_co)
plot(output_pooled$deltas_aca)
plot(output_pooled$deltas_aco)
print(output_pooled$accept)

plot(apply(output_pooled$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output_pooled$samples.w)
plot(x=W1, y=w.hat); abline(0, 1, col=2)
summary(100*(W1-w.hat)/W1)
view_tr_w(output_pooled$samples.w, w_true=W1)

view_tr(output_pooled$samples.theta, Theta)
print(mean(output_pooled$samples.theta)); print(Theta)

view_tr(output_pooled$samples.phi, Phi)
print(mean(output_pooled$samples.phi)); print(Phi)

par(mfrow=c(2, 3))
padded_plot(output_pooled$samples.beta.ca[,1], beta.case2[1])
padded_plot(output_pooled$samples.beta.ca[,2], beta.case2[2])
padded_plot(output_pooled$samples.beta.ca[,3], beta.case2[3])
print(colMeans(output_pooled$samples.beta.ca))

padded_plot(output_pooled$samples.beta.co[,1], beta.ctrl2[1])
padded_plot(output_pooled$samples.beta.co[,2], beta.ctrl2[2])
padded_plot(output_pooled$samples.beta.co[,3], beta.ctrl2[3])
print(colMeans(output_pooled$samples.beta.co))

par(mfrow=c(1,2))
padded_plot(output_pooled$samples.alpha.ca, Alpha.case1, title='A)')
print(mean(output_pooled$samples.alpha.ca))

padded_plot(output_pooled$samples.alpha.co, Alpha.ctrl1, title='B)')
print(mean(output_pooled$samples.alpha.co))

tag <- "_pooled"
output_pooled$description <- paste(sim_name, tag, sep="")
save_output(output_pooled, paste("output_", sim_name, tag, ".json", sep=""))

