############################
# Fit the spatiotemporal
# shared latent process model
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_temporal/"
agg_factor <- 12
n_sims <- 25


#### Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=agg_factor) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])
print(n_values(caPr.disc_all[[1]][[1]]))
print(mean(area(caPr.disc_all[[1]][[1]])[]))

cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


#### Specify trend
level <- "alternating"
sim_name <- paste("simPsTemporal", level, sep="_")
###################
#### Proposed model
###################
failures <- c(6, 8, 10, 24)
for (i in failures){
  
  print(paste("dataset", i))
  data <- reformat_saved_data(load_output(paste("data_", level, "_", i, ".json", sep=""), src=src))
  params <- load_output(paste("params_", level, "_", i, ".json", sep=""), src=src)
  
  caPr <- load_prism_pcs2()
  caPr.disc <- aggregate(caPr, fact=agg_factor)

  data_pooled <- pool_temporal_data(data)
  prior_theta <- get_gamma_prior(params$Theta, 5)
  prior_phi <- get_igamma_prior(params$Phi, 5)
  w_output <- logisticGp(y=data_pooled$loc$status, D, n.sample=1000, burnin=200, L=10,
                         prior_phi=prior_phi, prior_theta=prior_theta)
  plot(w_output$samples.theta, type='l'); abline(h=params$Theta)
  plot(w_output$samples.phi, type='l'); abline(h=params$Phi)
  plot(x=params$W, y=colMeans(w_output$samples.w));abline(0,1)
  
  w_initial <- colMeans(w_output$samples.w)
  theta_initial <- mean(w_output$samples.theta)
  phi_initial <- mean(w_output$samples.phi)
  
  
  # Beta & alpha initial values
  ini_case <- glm(data_pooled$case.data$y ~ data_pooled$case.data$x + w_initial[data_pooled$loc$ids] - 1, family='poisson')
  alpha_ca_initial <- coefficients(ini_case)[4]
  beta_ca_initial <- coefficients(ini_case)[1:3]
  
  ini_ctrl <- glm(data_pooled$ctrl.data$y ~ data_pooled$ctrl.data$x.standardised + w_initial[data_pooled$loc$ids] - 1, family='poisson')
  alpha_co_initial <- coefficients(ini_ctrl)[4]
  beta_co_initial <- coefficients(ini_ctrl)[1:3]
  
  u_initial=rep(0, length(years))
  prior_theta=get_gamma_prior(theta_initial, 2)
  prior_phi=get_igamma_prior(phi_initial, 2)
  prior_alpha_ca_mean=alpha_ca_initial
  prior_alpha_ca_var=2
  prior_alpha_co_mean=alpha_co_initial
  prior_alpha_co_var=2
  prior_u_mean=0
  prior_u_var=2
  
  # Tuning parameters
  m_aca=2000
  m_aco=2000
  m_ca=2000
  m_co=2000
  m_w=2000
  m_u=2000
  target_aca=0.65
  target_aco=0.65
  target_ca=0.65
  target_co=0.65
  target_w=0.65
  target_u=0.65
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
  
  n.sample=8000
  burnin=3000
  
  proposal.sd.theta=0.2
  
  #### Fit model
  output <- prefSampleTemporal(data, D, n.sample, burnin,
                               L_w, L_ca, L_co, L_a_ca, L_a_co, L_u,
                               proposal.sd.theta=proposal.sd.theta,
                               m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w, m_u=m_u,
                               target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, target_u=target_u,
                               self_tune_w=self_tune_w, self_tune_aca=self_tune_aca, self_tune_aco=self_tune_aco, self_tune_ca=self_tune_ca, self_tune_co=self_tune_co, self_tune_u=self_tune_u,
                               beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                               delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, delta_u=NULL,
                               theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial, u_initial=u_initial,
                               prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var,
                               prior_u_mean=prior_u_mean, prior_u_var=prior_u_var)
  
  #### Check estimated log odds
  # par(mfrow=c(2,4))
  # for (y_idx in 1:length(years)){
  #   year <- years[y_idx]
  #   params$alpha.case=params$Alpha.case
  #   params$alpha.ctrl=params$Alpha.ctrl
  #   lodds_true <- calc_true_lodds_time(params, data, year, y_idx, agg_factor=agg_factor)
  #   lodds_est_tmp <- calc_est_lodds_time(output, data, year, y_idx, agg_factor=agg_factor)
  #   n <- sum(data$locs[[y_idx]]$status)
  #   plot(x=lodds_true, y=lodds_est_tmp, main=paste("N: ", n)); abline(0,1,col=2)
  # }
  # par(mfrow=c(1,1))
  # plot(x=lodds_true, y=lodds_est_tmp, main=paste("N: ", n)); abline(0,1,col=2)
  w.hat <- colMeans(output$samples.w)
  u.hat <- colMeans(output$samples.u)
  par(mfrow=c(2,4))
  for (y_idx in 1:length(years)){
    n <- sum(data$locs[[y_idx]]$status)
    year <- years[y_idx]
    plot(x=params$U[y_idx] + params$W, y=u.hat[y_idx]+w.hat, main=paste("dataset: ", i, " N: ", n)); abline(0,1,col=2)
  }
  par(mfrow=c(1,1))
  
  #### Save output
  output$description <- paste(sim_name, "_", i, sep="")
  save_output(output, paste("output_", sim_name, "_", i, ".json", sep=""), dst=src)
  
}


# ####################
# #### Reference model
# ####################
level <- "decreasing"
sim_name <- paste("simPsTemporal", level, sep="_")

for (i in 1:n_sims){

  print(paste("dataset", i))
  data <- reformat_saved_data(load_output(paste("data_", level, "_", i, ".json", sep=""), src=src))
  params <- load_output(paste("params_", level, "_", i, ".json", sep=""), src=src)
  data_pooled <- pool_temporal_data(data)

  ## W initial value
  prior_theta <- get_gamma_prior(params$Theta, 5)
  prior_phi <- get_igamma_prior(params$Phi, 5)
  w_output <- logisticGp(y=data_pooled$locs$status, D, n.sample=1000, burnin=500, L=10,
                         prior_phi=prior_phi, prior_theta=prior_theta)
  w_i <- colMeans(w_output$samples.w)
  theta_i <- mean(w_output$samples.theta)
  phi_i <- mean(w_output$samples.phi)

  ## Beta & alpha initial values
  ini_case <- glm(data_pooled$case.data$y ~ data_pooled$case.data$x + w_i[data_pooled$locs$ids] - 1, family='poisson')
  alpha_ca_i <- coefficients(ini_case)[4]
  beta_ca_i <- coefficients(ini_case)[1:3]

  ini_ctrl <- glm(data_pooled$ctrl.data$y ~ data_pooled$ctrl.data$x.standardised + w_i[data_pooled$locs$ids] - 1, family='poisson')
  alpha_co_i <- coefficients(ini_ctrl)[4]
  beta_co_i <- coefficients(ini_ctrl)[1:3]

  ## Tuning parameters
  m_aca <- 2000
  m_aco <- 2000
  m_ca <- 2000
  m_co <- 2000
  m_w <- 2000
  proposal.sd.theta <- 0.15
  L_w <- 8
  L_ca <- 8
  L_co <- 8
  L_a_ca <- 8
  L_a_co <- 8

  n.sample <- 8000
  burnin <- 3000
  L_w <- 8
  L_ca <- 8
  L_co <- 8
  L_a_ca <- 8
  L_a_co <- 8
  
  prior_alpha_ca_var <- 3
  prior_alpha_co_var <- 3
  prior_alpha_co_mean <- alpha_co_i
  prior_alpha_ca_mean <- alpha_ca_i

  ## Run fit
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

  ## Check estimated log odds

  w.hat <- colMeans(output_ps$samples.w)
  beta_ca_h <- colMeans(output_ps$samples.beta.ca)
  beta_co_h <- colMeans(output_ps$samples.beta.co)
  alpha_ca_h <- mean(output_ps$samples.alpha.ca)
  alpha_co_h <- mean(output_ps$samples.alpha.co)
  phi_h <- mean(output_ps$samples.phi)
  theta_h <- mean(output_ps$samples.theta)
  
  Alpha.case <- params$Alpha.case
  Alpha.ctrl <- params$Alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  X.standard <- load_x_standard2(as.logical(data_pooled$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  
  par(mfrow=c(2,4))
  for (y_idx in 1:length(years)){
    year <- years[y_idx]
    params$alpha.case=Alpha.case
    params$alpha.ctrl=Alpha.ctrl
    u_idx <- params$U[y_idx]
    lodds.true <- X.standard %*% beta.case + Alpha.case * (u_idx + W) - X.standard %*% beta.ctrl - Alpha.ctrl * (u_idx + W)
    lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
    n <- sum(data$locs[[y_idx]]$status)
    plot(x=lodds.true, y=lodds.ps, main=paste("N: ", n)); abline(0,1,col=2)
  }
  par(mfrow=c(1,1))
  plot(x=lodds.true, y=lodds.ps, main=paste("N: ", n)); abline(0,1,col=2)
  
  #### Save output
  output_ps$description <- paste(sim_name, "_", i, "_pooled", sep="")
  save_output(output_ps, paste("output_pooled_", sim_name, "_", i, ".json", sep=""), dst=src)

}


### Check u + w estimation

level <- "alternating"
sim_name <- paste("simPsTemporal", level, sep="_")
for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- reformat_saved_data(load_output(paste("data_", level, "_", i, ".json", sep=""), src=src))
  params <- load_output(paste("params_", level, "_", i, ".json", sep=""), src=src)
  
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  w.hat <- colMeans(output$samples.w)
  u.hat <- colMeans(output$samples.u)
  par(mfrow=c(2,4))
  for (y_idx in 1:length(years)){
    n <- sum(data$locs[[y_idx]]$status)
    year <- years[y_idx]
    plot(x=params$U[y_idx] + params$W, y=u.hat[y_idx]+w.hat, main=paste("dataset: ", i, " N: ", n)); abline(0,1,col=2)
  }
  par(mfrow=c(1,1))
  
}
