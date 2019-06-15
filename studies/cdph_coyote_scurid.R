
library(mvtnorm)
library(R.utils)
library(plyr)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')

###############
# models to fit
#   MVGP: rodents + coyotes
#   Preferential sampling: coyotes
#   Preferential sampling: rodents (already done)
###############


analysis_name <- "cdph_coyote_scurid"
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=6)
N <- n_values(caPr.disc[[1]])
plot(caPr.disc)

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
coyotes <- read.csv(paste(src, "CDPH_coyote_recoded_full.csv", sep=""), header=T, sep=",")

loc.disc <- caPr.disc[[1]]
data_rodents <- assemble_data(rodents, loc.disc, caPr.disc)
data_coyotes <- assemble_data_coyotes(coyotes, caPr.disc)

cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

plot(caPr[[1]])
points(data_rodents$loc$coords, col='red', pch=16)
points(data_coyotes$loc$coords, col='blue')

xy_rodents <- cbind(rodents$Lon_Add_Fix, rodents$Lat_Add_Fix)
xy_coyotes <- cbind(coyotes$Long_QC, coyotes$Lat_QC)
plot(caPr[[1]])
points(xy_rodents, col=rgb(1, 0, 0, 0.25), pch=16)
points(xy_coyotes, col=rgb(0, 0, 1, 0.25))

# qa
par(mfrow=c(1,2))
plot(loc.disc)
points(data_rodents$loc$coords)
plot(loc.disc)
points(data_coyotes$loc$coords)

count_sum <- list()
count_sum[[1]] <- list(
  Species='Rodents',
  Cases=sum(data_rodents$case.data$y),
  Controls=sum(data_rodents$ctrl.data$y),
  Prevalence=round(sum(data_rodents$case.data$y)/sum(data_rodents$case.data$y + data_rodents$ctrl.data$y), 3)
)
count_sum[[2]] <- list(
  Species='Coyotes',
  Cases=sum(data_coyotes$case.data$y),
  Controls=sum(data_coyotes$ctrl.data$y),
  Prevalence=round(sum(data_coyotes$case.data$y)/sum(data_coyotes$case.data$y + data_coyotes$ctrl.data$y), 3)
)
count_df <- ldply(count_sum, 'data.frame')
print(count_df)
write_latex_table(count_df, 'latex_count_summary.txt', '/Users/brianconroy/Documents/research/Project2/cdph_coyotes_rodents/')

data <- list(
  locs=list(data_rodents$loc, data_coyotes$loc),
  case.data=list(data_rodents$case.data, data_coyotes$case.data),
  ctrl.data=list(data_rodents$ctrl.data, data_coyotes$ctrl.data)
)
save_output(data, paste(analysis_name, '_data.json', sep=''))

#############
#### Fit MVGP
#############

#### Set Priors
prior_theta <- c(3, 2)
prior_alpha_ca_mean <- c(0, 0)
prior_alpha_co_mean <- c(0, 0)
prior_alpha_ca_var <- c(4, 4)
prior_alpha_co_var <- c(4, 4)
Omega <- matrix(c(5, 0, 0, 5), nrow=2)
r <- 4
print(Omega/(r - length(data$case.data) - 1))

#### Fit initial values
y <- list(data_rodents$loc$status, data_coyotes$loc$status)
locs <- list(data_rodents$loc, data_coyotes$loc)
prior_t=list(scale=matrix(c(5, 0, 0, 5), nrow=2), df=4)
w_output <- logisticMVGP(y, D, n.sample=1000, burnin=200, L=10, 
                         prior_t=prior_t, prior_theta=prior_theta)

w_initial=colMeans(w_output$samples.w)
theta_initial <- mean(w_output$samples.theta)
t_initial <- matrix(colMeans(w_output$samples.t), nrow=2)

plot(w_output$samples.theta, type='l')
plot(w_output$samples.t[,1], type='l')
plot(w_output$samples.t[,2], type='l')
plot(w_output$samples.t[,4], type='l')

save_output(w_output, "cdph_w_inival_output_mvgp.json")

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

output <- continueMCMC_mvgp(data, D, output, n.sample=2500)

output <- burnin_mvgp(output, n.burn=400)

print(output$accept)

plot(output$samples.theta, type='l')
plot(output$samples.t[,1], type='l')
plot(output$samples.t[,2], type='l')
plot(output$samples.t[,4], type='l')

plot(output$samples.alpha.ca[1,,], type='l')
plot(output$samples.alpha.co[1,,], type='l')
plot(output$samples.alpha.ca[2,,], type='l')
plot(output$samples.alpha.co[2,,], type='l')

tag <- paste(analysis_name, "mvgp", sep="_")
output$description <- tag
save_output(output, paste("output_", tag, ".json", sep=""))

###################
#### Fit Single
#### Species Model
#### Coyotes
###################

# W initial value
prior_theta <- c(3, 2)
prior_phi <- c(18, 204)
w_output_c <- logisticGp(y=data_coyotes$loc$status, D, n.sample=1000, burnin=500, L=10,
                         prior_phi=prior_phi, prior_theta=prior_theta,
                         proposal.sd.theta=0.20)
w_output_c$accept
view_tr_w(w_output_c$samples.w)
view_tr(w_output_c$samples.theta)
view_tr(w_output_c$samples.phi)

hist(colMeans(w_output_c$samples.w))
save_output(w_output_c, "w_inival_output_c_cdph.json")

w_i <- colMeans(w_output_c$samples.w)
theta_i <- mean(w_output_c$samples.theta)
phi_i <- mean(w_output_c$samples.phi)

ini_case <- glm(data_coyotes$case.data$y ~ data_coyotes$case.data$x.standardised + w_i[data_coyotes$loc$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(data_coyotes$ctrl.data$y ~ data_coyotes$ctrl.data$x.standardised + w_i[data_coyotes$loc$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

# set prior theta and phi to the values estimated
# in the previous step
# iterate_g_variance(theta_i, 1, 15)
# iterate_ig_variance(phi_i, 400, 500)

prior_theta <- c(1.136, 3)
prior_phi <- c(28.5, 500)

g_var(1.136, 3)
ig_var(28.5, 500)

prior_alpha_ca_mean <- 0
prior_alpha_co_mean <- 0
prior_alpha_ca_var <- 6
prior_alpha_co_var <- 6

n.sample <- 4000
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

output_c <- prefSampleGpCC(data_coyotes, D, n.sample, burnin,
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

output_c <- burnin_after(output_c, n.burn=200)

output_c <- continueMCMC(data, output_c, n.sample=3000)

plot(apply(output_c$samples.w, 1, mean), type='l')
view_tr_w(output_c$samples.w)

view_tr(output_c$samples.alpha.ca)
view_tr(output_c$samples.alpha.co)

par(mfrow=c(2,3))
view_tr(output_c$samples.beta.ca[,1])
view_tr(output_c$samples.beta.ca[,2])
view_tr(output_c$samples.beta.ca[,3])

view_tr(output_c$samples.beta.co[,1])
view_tr(output_c$samples.beta.co[,2])
view_tr(output_c$samples.beta.co[,3])

par(mfrow=c(1,1))
view_tr(output_c$samples.theta)
view_tr(output_c$samples.phi)

tag <- paste(analysis_name, "coyote", sep="_")
output_c$description <- tag
save_output(output_c, paste("output_", tag, ".json", sep=""))


###################
#### Fit Single
#### Species Model
#### Rodents
###################

# W initial value
prior_theta <- c(3, 2)
prior_phi <- c(18, 204)
w_output_r <- logisticGp(y=data_rodents$loc$status, D, n.sample=1000, burnin=500, L=10,
                         prior_phi=prior_phi, prior_theta=prior_theta,
                         proposal.sd.theta=0.20)
w_output_r$accept
view_tr_w(w_output_r$samples.w)
view_tr(w_output_r$samples.theta)
view_tr(w_output_r$samples.phi)

hist(colMeans(w_output_r$samples.w))
save_output(w_output_r, "w_inival_output_r_cdph.json")

w_i <- colMeans(w_output_r$samples.w)
theta_i <- mean(w_output_r$samples.theta)
phi_i <- mean(w_output_r$samples.phi)

ini_case <- glm(data_rodents$case.data$y ~ data_rodents$case.data$x.standardised + w_i[data_rodents$loc$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(data_rodents$ctrl.data$y ~ data_rodents$ctrl.data$x.standardised + w_i[data_rodents$loc$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

# set prior theta and phi to the values estimated
# in the previous step
# iterate_g_variance(theta_i, 1, 15)
# iterate_ig_variance(phi_i, 400, 500)

prior_theta <- c(1.136, 3)
prior_phi <- c(28.5, 500)

g_var(1.136, 3)
ig_var(28.5, 500)

prior_alpha_ca_mean <- 0
prior_alpha_co_mean <- 0
prior_alpha_ca_var <- 6
prior_alpha_co_var <- 6

n.sample <- 4000
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

output_r <- prefSampleGpCC(data_rodents, D, n.sample, burnin,
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

output_r <- burnin_after(output_r, n.burn=700)

output_r <- continueMCMC(data_rodents, D, output_r, n.sample=500)

plot(apply(output_r$samples.w, 1, mean), type='l')
view_tr_w(output_r$samples.w)

view_tr(output_r$samples.alpha.ca)
view_tr(output_r$samples.alpha.co)

par(mfrow=c(2,3))
view_tr(output_r$samples.beta.ca[,1])
view_tr(output_r$samples.beta.ca[,2])
view_tr(output_r$samples.beta.ca[,3])

view_tr(output_r$samples.beta.co[,1])
view_tr(output_r$samples.beta.co[,2])
view_tr(output_r$samples.beta.co[,3])

par(mfrow=c(1,1))
view_tr(output_r$samples.theta)
view_tr(output_r$samples.phi)

tag <- paste(analysis_name, "rodent", sep="_")
output_r$description <- tag
save_output(output_r, paste("output_", tag, ".json", sep=""))
