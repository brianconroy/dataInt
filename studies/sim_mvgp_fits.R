############################
# Fit the spatiotemporal
# shared latent process model
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_mvgp/"
agg_factor <- 12
n_sims <- 25
sim_name <- "sim_mvgp"


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


###################
#### Proposed model
###################
for (i in 1:n_sims){
  
  print(paste("dataset", i))
  data <- reformat_saved_mvgp(load_output(paste("data_", i, ".json", sep=""), src=src))
  params <- load_output(paste("params_", i, ".json", sep=""), src=src)
  
  # Initial Values
  y <- list(data$locs[[1]]$status, data$locs[[2]]$status)
  locs <- list(data$locs[[1]], data$locs[[2]])
  prior_t=list(scale=matrix(c(5, 0, 0, 5), nrow=2), df=4)
  prior_theta=get_gamma_prior(prior_mean=params$Theta, prior_var=5)
  w_output <- logisticMVGP(y, locs, d, n.sample=700, burnin=250, L=10, prior_t=prior_t, prior_theta=prior_theta)
  
  w_initial <- colMeans(w_output$samples.w)
  theta_initial <- mean(w_output$samples.theta)
  t_initial <- matrix(colMeans(w_output$samples.t), nrow=2)
  
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
  
  # Tuning parameters
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
  
  # Priors
  prior_alpha_ca_mean=c(alpha_ca_initial[[1]], alpha_ca_initial[[2]])
  prior_alpha_ca_var=list(4, 4)
  prior_alpha_co_mean=c(alpha_co_initial[[1]], alpha_co_initial[[2]])
  prior_alpha_co_var=list(4, 4)
  prior_u_mean=0
  prior_theta=get_gamma_prior(theta_initial, 5)

  n.sample=5000
  burnin=2000
  
  # Fit model
  output <- prefSampleMVGP(data, d, n.sample, burnin,
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
  
  # Plot log odds
  par(mfrow=c(1,2))
  for (i in 1:2){
    lodds <- calc_lodds_mvgp(output, data, species=i, agg_factor=agg_factor)
    lodds_true <- calc_lodds_true_multi(params, data, species=i, agg_factor=agg_factor)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', main=paste("species", i)); abline(0, 1, col=2)
  }
  
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
