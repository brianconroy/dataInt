#######################
# Fit the discrete time
# spatiotemporal model
# to Sciurid data
#######################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


analysis_name <- "cdph_temporal_analysis"
agg_factor <- 7
src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=agg_factor) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
  plot(caPr_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])
loc.disc_y <- caPr.disc_all[[1]]


#### Assemble data
locs <- list()
case.data <- list()
ctrl.data <- list()
bin_width <- 5

for (i in 1:length(years)){
  y <- years[i]
  y_ub <- y + (bin_width - 1)/2
  y_lb <- y - (bin_width - 1)/2
  
  rodents_y <- rodents[rodents$Year <= y_ub & rodents$Year >= y_lb,]
  loc.disc_y <- caPr.disc_all[[i]][[1]]
  caPr.disc_y <- caPr.disc_all[[i]]
  data_y <- assemble_data(rodents_y, loc.disc_y, caPr.disc_y)
  
  locs[[i]] <- data_y$loc
  case.data[[i]] <- data_y$case.data
  ctrl.data[[i]] <- data_y$ctrl.data
}


cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
data <- list(locs=locs, case.data=case.data, ctrl.data=ctrl.data)


###########
# fit model
###########


#### Specify MCMC parameters
m_aca=2000
m_aco=2000
m_ca=700
m_co=700
m_w=700
m_u=700
target_aca=0.75
target_aco=0.75
target_ca=0.75
target_co=0.75
target_w=0.75
target_u=0.75
self_tune_w=TRUE
self_tune_aca=TRUE
self_tune_aco=TRUE
self_tune_ca=TRUE
self_tune_co=TRUE
self_tune_u=TRUE

L_w=8
L_ca=8
L_co=8
L_a_ca=8
L_a_co=8
L_u=8


## Initial values
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=7)
data_pooled <- assemble_data(rodents, loc.disc_y, caPr.disc)

prior_theta <- c(6, 1)
prior_phi <- c(18, 204)
w_output <- logisticGp(y=data_pooled$loc$status, D, n.sample=1000, burnin=200, L=10,
                       prior_phi=prior_phi, prior_theta=prior_theta)
# w_output <- burnin_logisticGp_mcmc(w_output, n.burn=50)
plot(w_output$samples.theta, type='l')
plot(w_output$samples.phi, type='l')

save_output(w_output, paste("w_inival_output_", analysis_name, ".json", sep=""))

w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)


# Beta & alpha initial values
ini_case <- glm(data_pooled$case.data$y ~ data_pooled$case.data$x + w_i[data_pooled$loc$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(data_pooled$ctrl.data$y ~ data_pooled$ctrl.data$x.standardised + w_i[data_pooled$loc$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

w_initial=w_i
beta_ca_initial=beta_ca_i
beta_co_initial=beta_co_i
alpha_ca_initial=alpha_ca_i
alpha_co_initial=alpha_co_i
theta_initial=theta_i
phi_initial=phi_i
u_initial=rep(0, length(years))

prior_theta=get_gamma_prior(theta_i, 3)
prior_phi=get_igamma_prior(phi_initial, 3)
prior_alpha_ca_mean=alpha_ca_initial
prior_alpha_ca_var=2
prior_alpha_co_mean=alpha_ca_initial
prior_alpha_co_var=2
prior_u_mean=0
prior_u_var=2

n.sample=1000
burnin=500
proposal.sd.theta=0.2

output <- prefSampleTemporal(data, D, n.sample, burnin,
                             L_w, L_ca, L_co, L_a_ca, L_a_co, L_u,
                             proposal.sd.theta=proposal.sd.theta,
                             m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, m_u=m_u,
                             target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_u=target_u,
                             self_tune_w=self_tune_w, self_tune_aca=self_tune_aca, self_tune_aco=self_tune_aco, self_tune_ca=self_tune_ca, self_tune_co=self_tune_co, self_tune_u=self_tune_u,
                             beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                             theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial, u_initial=u_initial,
                             prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var,
                             prior_u_mean=prior_u_mean, prior_u_var=prior_u_var)


#### Additional burnin
output <- burnin_after_temporal(output, n.burn=500)


#### Generate additional samples
output <- continue_mcmc_temporal(data, D, output, n.sample=2000)


#### View results
par(mfrow=c(2,3))
plot(output$samples.beta.ca[,1], type='l')
plot(output$samples.beta.ca[,2], type='l')
plot(output$samples.beta.ca[,3], type='l')
plot(output$samples.beta.co[,1], type='l')
plot(output$samples.beta.co[,2], type='l')
plot(output$samples.beta.co[,3], type='l')

par(mfrow=c(1,2))
plot(output$samples.alpha.ca, type='l')
plot(output$samples.alpha.co, type='l')

plot(output$samples.theta, type='l')
plot(output$samples.phi, type='l')

w.hat <- colMeans(output$samples.w)
u.hat <- colMeans(output$samples.u)
par(mfrow=c(2,4))
hist(w.hat + u.hat[1])
hist(w.hat + u.hat[2])
hist(w.hat + u.hat[3])
hist(w.hat + u.hat[4])
hist(w.hat + u.hat[5])
hist(w.hat + u.hat[6])
hist(w.hat + u.hat[7])

par(mfrow=c(1,2))
plot(mean(w.hat) + u.hat, type='l')
prevalences <- c()
for (i in 1:length(years)){
  prevalences <- c(prevalences,
                   sum(data$case.data[[i]]$y)/(sum(data$case.data[[i]]$y) + sum(data$ctrl.data[[i]]$y)))
  
}
plot(prevalences, type='l')


#### Save output
output$description <- analysis_name
save_output(output, paste("output_", analysis_name, ".json", sep=""))


#### Calculate temporal risk surfaces
temporal_risks <- calc_temporal_risks(output, caPr.disc_all, caPr_all, agg_factor, years)

rescaled <- equalize_scales4(temporal_risks$r_risks_high)
par(mfrow=c(2,4))
for (i in 1:length(years)){
  plot(rescaled[[i]], main=years[i])
}


#####################################
# fit reference model
# (nontemporal preferential sampling)
#####################################


caPr <- load_prism_pcs2()
loc.disc <- caPr.disc[[1]]
caPr.disc <- aggregate(caPr, fact=agg_factor)
data_pooled <- assemble_data(rodents, loc.disc, caPr.disc)
print(sum(data_pooled$case.data$y + data_pooled$ctrl.data$y))
print(sum(data_pooled$case.data$y)/sum(data_pooled$case.data$y + data_pooled$ctrl.data$y))

## W initial value
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)
# w_output <- logisticGp(y=data_pooled$locs$status, D, n.sample=1000, burnin=200, L=10,
#                        prior_phi=prior_phi, prior_theta=prior_theta)
# view_logistic_output(w_output)
# save_output(w_output, paste("w_inival_output_", analysis_name, ".json", sep=""))
w_output <- load_output(paste("w_inival_output_", analysis_name, ".json", sep=""))
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

# Beta & alpha initial values
ini_case <- glm(data_pooled$case.data$y ~ data_pooled$case.data$x + w_i[data_pooled$loc$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(data_pooled$ctrl.data$y ~ data_pooled$ctrl.data$x.standardised + w_i[data_pooled$loc$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

n.sample <- 8000
burnin <- 1000
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

prior_theta=get_gamma_prior(theta_i, 3)
prior_phi=get_igamma_prior(phi_i, 3)
prior_alpha_ca_mean=alpha_ca_i
prior_alpha_ca_var=2
prior_alpha_co_mean=alpha_co_i
prior_alpha_co_var=2

output_ps <- prefSampleGpCC(data_pooled, D, n.sample, burnin,
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
output_ps <- burnin_after(output_ps, n.burn=1000)


# optionally continue running if necessary
output_ps <- continueMCMC(data_pooled, D, output_ps, n.sample=2000)


w.hat <- colMeans(output_ps$samples.w)
par(mfrow=c(1,2))
view_tr(output_ps$samples.theta)
view_tr(output_ps$samples.phi)

par(mfrow=c(2, 3))
view_tr(output_ps$samples.beta.ca[,1])
view_tr(output_ps$samples.beta.ca[,2])
view_tr(output_ps$samples.beta.ca[,3])
view_tr(output_ps$samples.beta.co[,1])
view_tr(output_ps$samples.beta.co[,2])
view_tr(output_ps$samples.beta.co[,3])

par(mfrow=c(1,2))
view_tr(output_ps$samples.alpha.ca)
view_tr(output_ps$samples.alpha.co)

## save results
tag <- paste(analysis_name, "_aggregated_ps", sep="")
output_ps$description <- tag
save_output(output_ps, paste('output_', tag, ".json", sep=""))

