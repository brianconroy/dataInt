library(jsonlite)


get_gamma_prior <- function(prior_mean, prior_var){
  
  shape <- prior_mean^2/prior_var
  scale <- prior_var/prior_mean
  return(c(shape, scale))
  
}


get_igamma_prior <- function(prior_mean, prior_var){
  
  scale <- prior_mean^3 * (1/prior_var + 1/prior_mean^2)
  shape <- scale/prior_mean + 1
  return(c(shape, scale))
  
}


multispecies_plot <- function(samples, species, parameter, beta_id=NULL){
  
  if (grepl('Alpha.case', parameter)){
    true <- data_store$alpha.cases[species]
  } else if (grepl('Alpha.ctrl', parameter)){
    true <- data_store$alpha.ctrl[species]
  } else if (grepl('beta.case', parameter)){
    true <- data_store$beta.cases[species, beta_id]
  } else if (grepl('beta.ctrl', parameter)){
    true <- data_store$beta.ctrls[species, beta_id]
  } 
  
  if (!is.null(beta_id)){
    padded_plot(samples[species,,beta_id], true)
  } else {
    padded_plot(samples[species,,], true)
  }
  
}


summarize_mcmc_pscc <- function(output, model_desc){
  
  rows <- list()
  
  rows[[1]] <- list(
    model=model_desc,
    parameter='n.sample',
    value=output$n.sample
  )
  rows[[2]] <- list(
    model=model_desc,
    parameter='burnin',
    value=output$burnin
  )
  rows[[3]] <- list(
    model=model_desc,
    parameter='proposal.sd.theta',
    value=output$proposal.sd.theta
  )
  rows[[4]] <- list(
    model=model_desc,
    parameter='L (w)',
    value=output$L_w
  )
  rows[[5]] <- list(
    model=model_desc,
    parameter='L (beta case)',
    value=output$L_ca
  )
  rows[[6]] <- list(
    model=model_desc,
    parameter='L (beta control)',
    value=output$L_co
  )
  rows[[7]] <- list(
    model=model_desc,
    parameter='L (alpha case)',
    value=output$L_a_ca
  )
  rows[[8]] <- list(
    model=model_desc,
    parameter='L (alpha control)',
    value=output$L_a_co
  )
  rows[[9]] <- list(
    model=model_desc,
    parameter='m',
    value=1000
  )
  rows[[10]] <- list(
    model=model_desc,
    parameter='target acceptance',
    value=0.65
  )
  
  return(rows)
  
}


load_priors <- function(param, index){
  
  f1 <- 'simParams_alpha_prior.txt'
  f2 <- 'simParams_theta_prior.txt'
  f3 <- 'simParams_phi_prior.txt'
  priors_alpha <- load_output(f1)
  priors_theta <- load_output(f2)
  priors_phi <- load_output(f3)
  priors <- list()
  if (param == 'alpha'){
    priors$prior_alpha <- priors_alpha[index]
    priors$prior_phi <- priors_phi[2,]
    priors$prior_theta <- priors_theta[3,]
  } else if (param == 'phi'){
    priors$prior_alpha <- 4
    priors$prior_phi <- priors_phi[index,]
    priors$prior_theta <- c(2.5, 2.5)
  } else if (param == 'theta'){
    priors$prior_alpha <- 4
    priors$prior_phi <- c(14, 156)
    priors$prior_theta <- priors_theta[index,]
  }
  return(priors)
  
}


calc_log_odds_output <- function(output, true_params){
  
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  
  lodds.ps <- x_standard %*% beta_ca_h + alpha_ca_h * w.hat - x_standard %*% beta_co_h - alpha_co_h * w.hat
  return(lodds.ps)
  
}


calc_log_odds_multi <- function(output, true_params, species, output_type){
  
  location_indicators <- as.logical(true_params$data$locs$status[[species]])
  x_standard <- load_x_standard(location_indicators)
  
  w.hat <- colMeans(output$samples.w)
  if (output_type == '_multi'){
    beta_ca_h <- colMeans(output$samples.beta.ca[species,,])
    beta_co_h <- colMeans(output$samples.beta.co[species,,])
    alpha_ca_h <- mean(output$samples.alpha.ca[species,,])
    alpha_co_h <- mean(output$samples.alpha.co[species,,])
  } else {
    beta_ca_h <- colMeans(output$samples.beta.ca)
    beta_co_h <- colMeans(output$samples.beta.co)
    alpha_ca_h <- mean(output$samples.alpha.ca)
    alpha_co_h <- mean(output$samples.alpha.co)
  }
  
  lodds.ps <- x_standard %*% beta_ca_h + alpha_ca_h * w.hat - x_standard %*% beta_co_h - alpha_co_h * w.hat
  return(lodds.ps)
  
}


calc_log_odds_true_multi <- function(true_params, species){
  
  location_indicators <- as.logical(true_params$data$locs$status[[species]])
  x_standard <- load_x_standard(location_indicators)
  if (species == 1){
    beta.case <- true_params$beta.case1
    beta.ctrl <- true_params$beta.ctrl1
    Alpha.case <- true_params$Alpha.case1
    Alpha.ctrl <- true_params$Alpha.ctrl1
    W <- true_params$W1
  } else if (species == 2){
    beta.case <- true_params$beta.case2
    beta.ctrl <- true_params$beta.ctrl2
    Alpha.case <- true_params$Alpha.case2
    Alpha.ctrl <- true_params$Alpha.ctrl2
    W <- true_params$W2
  }
  
  lodds.true <- x_standard %*% beta.case + Alpha.case * W - x_standard %*% beta.ctrl - Alpha.ctrl * W
  return(lodds.true)
  
}


#' calc_log_odds_true_multi_general
#' 
#' general multispecies output where 
#' true_params store values as matrices for 
#' different species.
#'
#' @param true_params (list)
#' @param species (numeric) 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds_true_multi_general <- function(true_params, species){
  
  location_indicators <- as.logical(true_params$data$locs$status[[species]])
  x_standard <- load_x_standard(location_indicators)
  beta.case <- true_params$beta.cases[species,]
  beta.ctrl <- true_params$beta.ctrls[species,]
  Alpha.case <- true_params$alpha.cases[species,]
  Alpha.ctrl <- true_params$alpha.ctrls[species,]
  W <- true_params$W
  lodds.true <- x_standard %*% beta.case + Alpha.case * W - x_standard %*% beta.ctrl - Alpha.ctrl * W
  return(lodds.true)
  
}


get_output_priorsens <- function(outputs, param, number){
  
  output_target <- list()
  tag <- paste(param, number, sep="_")
  for (o in outputs){
    if (grepl(tag, o$description)){
      output_target <- o
      break
    }
  }
  return(output_target)
  
}


#' load_sim_outputs_priorsens
#' 
#' load simulation outputs for the prior sensitivity study
#'
#' @return
#' @export
#'
#' @examples
load_sim_outputs_priorsens <- function(){
  
  output_list <- list()
  
  counter <- 1
  for (f in list.files('/Users/brianconroy/Documents/research/dataInt/output/')){
    if (grepl('output_iterate_prior', f)) {
      output_list[[counter]] <- load_output(f)
      counter <- counter + 1
    }
  }
  
  return(output_list)
  
}


#' table_params
#'
#' @param outputs (list) mcmc outputs
#' @param sampling (character) strength of sampling
#' @param prevalence (character) disease prevalence ("low", "medium", "high")
#'
#' @return (list) list of lists containing parameter
#' estimates, true values and biases of different models
#' @export
#'
#' @examples
table_params <- function(outputs, sampling, prevalence){
  
  output_ps <- get_output(outputs, sampling, prevalence, 'prefSampleGpCC')
  output_sp_ca <- get_output(outputs, sampling, prevalence, 'spatial_poisson_case')
  output_sp_co <- get_output(outputs, sampling, prevalence, 'spatial_poisson_ctrl')
  betas <- load_params(paste('estimates_poisson_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
  betas$description <- "Poisson Regression"
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  
  rows <- list()
  shared <- c("Beta 0 (case)", "Beta 1 (case)", "Beta 2 (case)", 
              "Beta 0 (control)", "Beta 1 (control)", "Beta 2 (control)")
  counter <- 1
  for (p in shared){
    rows[[counter]] <- make_row(
      prevalence,
      sampling,
      "PS",
      p,
      true_params,
      output_ps
    )
    rows[[counter + 1]] <- make_row(
      prevalence,
      sampling,
      "SP",
      p,
      true_params,
      output_sp_ca
    )
    rows[[counter + 1]] <- make_row(
      prevalence,
      sampling,
      "PR",
      p,
      true_params,
      betas
    )
    counter <- counter + 3
  }
  
  ps <- c("Alpha (case)", "Alpha (control)", "Theta", "Phi")
  for (p in ps){
    rows[[counter]] <- make_row(
      prevalence,
      sampling,
      "PS",
      p,
      true_params,
      output_ps
    )
    counter <- counter + 1
  }
  
  return(rows)
}


make_row <- function(prevalence, sampling, model, parameter, true_params, output){
  
  est <- get_estimate(output, parameter)
  true <- get_true_val(true_params, parameter)
  row <- list(
    prevalence=prevalence,
    sampling=sampling,
    model=model,
    parameter=parameter,
    estimate=est,
    true=true,
    bias=est-true
  )
  return(row)
  
}


get_true_val <- function(true_params, parameter){
  
  if (parameter == "Beta 0 (case)"){
    return(true_params$beta.case[1])
  } else if (parameter == "Beta 1 (case)"){
    return(true_params$beta.case[2])
  } else if (parameter == "Beta 2 (case)"){
    return(true_params$beta.case[3])
  } else if (parameter == "Beta 0 (control)"){
    return(true_params$beta.ctrl[1])
  } else if (parameter == "Beta 1 (control)"){
    return(true_params$beta.ctrl[2])
  } else if (parameter == "Beta 2 (control)"){
    return(true_params$beta.ctrl[3])
  } else if (parameter == "Alpha (case)"){
    return(true_params$Alpha.case)
  } else if (parameter == "Alpha (control)"){
    return(true_params$Alpha.ctrl)
  } else if (parameter == "W"){
    return(true_params$W)
  } else if (parameter == "Theta"){
    return(true_params$Theta)
  } else if (parameter == "Phi"){
    return(true_params$Phi)
  }
  
}


get_estimate <- function(output, parameter){
  
  if (grepl("prefSampleGpCC", output$description)){
    type <- "PS"
  } else if (grepl("Poisson Regression", output$description)){
    type <- "PR"
  } else {
    type <- "SP"
  }
  
  if (parameter == "Beta 0 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,1]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,1]
    } else {
      target_samples <- output$case[1]
    }
  } else if (parameter == "Beta 1 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,2]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,2]
    } else {
      target_samples <- output$case[2]
    }
  } else if (parameter == "Beta 2 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,3]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,3]
    } else {
      target_samples <- output$case[3]
    }
  } else if (parameter == "Beta 0 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,1]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,1]
    } else {
      target_samples <- output$ctrl[1]
    }
  } else if (parameter == "Beta 1 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,2]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,2]
    } else {
      target_samples <- output$ctrl[2]
    }
  } else if (parameter == "Beta 2 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,3]
    } else if (type == "SP") {
      target_samples <- output$samples.beta[,3]
    } else {
      target_samples <- output$ctrl[3]
    }
  } else if (parameter == "Alpha (case)"){
    target_samples <- output$samples.alpha.ca
  } else if (parameter == "Alpha (control)"){
    target_samples <- output$samples.alpha.co
  } else if (parameter == "Theta"){
    target_samples <- output$samples.theta
  } else if (parameter == "Phi"){
    target_samples <- output$samples.phi
  }
  
  est <- mean(target_samples)
  return(round(est, 3))
  
}


get_pad <- function(val, typ){
  if (typ == 'lower'){
    if (val < 0){
      return(1.15 * val)
    } else if (val == 0) {
      return(-0.15)
    } else {
      return(0.75 * val)
    } 
  } else {
    if (val < 0){
      return(0.75 * val)
    } else if (val == 0) {
      return(0.15)
    } else {
      return(1.15 * val)
    } 
  }
}


padded_plot <- function(y, trueval, ylab='', title=''){
  
  ymax <- max(y)
  ymin <- min(y)
  if (trueval < ymin){
    lb <- get_pad(trueval, 'lower')
    ub <- get_pad(ymax, 'upper')
    plot(y, typ='l', ylab=ylab, ylim=c(lb, ub), main=title); abline(h=trueval, col='2')
  } else if (trueval > ymax){
    ub <- get_pad(trueval, 'upper')
    lb <- get_pad(ymin, 'lower')
    plot(y, typ='l', ylab=ylab, ylim=c(lb, ub), main=title); abline(h=trueval, col='2')
  } else{
    plot(y, typ='l', ylab=ylab, main=title); abline(h=trueval, col='2')
  }
  
}


padded_plot2 <- function(y, ylab=''){
  
  ymax <- max(y)
  ymin <- min(y)
  if (ymax > 0){
    ub <- 2 * ymax
    if (ymin > 0){
      lb <- 0
    } else{
      lb <- 2 * ymin
    }
  } else {
    lb <- 2 * ymin
    ub <- 2 * abs(ub)
  }
  plot(y, ylim=c(lb, ub), type='l')
  
}


#' plot_traces
#' 
#' plots traceplots for the preferential sampling model
#'
#' @param outputs (list) mcmc preferential sampling outputs
#' @param sampling (character) strength of preferential sampling
#' @param prevalence (character) disease prevalence
#'
#' @return
#' @export
#'
#' @examples
plot_traces <- function(outputs, sampling, prevalence){
  
  output <- get_output(outputs, sampling, prevalence, 'prefSampleGpCC')
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  
  par(mfrow=c(3,4))
  padded_plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', true_params$beta.case[1])
  padded_plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', true_params$beta.case[2])
  padded_plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', true_params$beta.case[3])
  padded_plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', true_params$beta.ctrl[1])
  padded_plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', true_params$beta.ctrl[2])
  padded_plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', true_params$beta.ctrl[3])
  padded_plot(output$samples.alpha.ca, ylab='Alpha (case)', true_params$Alpha.case[1])
  padded_plot(output$samples.alpha.co, ylab='Alpha (control)', true_params$Alpha.ctrl[1])
  padded_plot(output$samples.theta, ylab='Range', true_params$Theta[1])
  padded_plot(output$samples.phi, ylab='Marginal Variance', true_params$Phi)
  par(mfrow=c(1,1))
  
}


plot_traces_general <- function(output, true_params=NULL){
  
  par(mfrow=c(3,4))
  if (!is.null(true_params)){
    padded_plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', true_params$beta.case[1])
    padded_plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', true_params$beta.case[2])
    padded_plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', true_params$beta.case[3])
    padded_plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', true_params$beta.ctrl[1])
    padded_plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', true_params$beta.ctrl[2])
    padded_plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', true_params$beta.ctrl[3])
    padded_plot(output$samples.alpha.ca, ylab='Alpha (case)', true_params$Alpha.case[1])
    padded_plot(output$samples.alpha.co, ylab='Alpha (control)', true_params$Alpha.ctrl[1])
    padded_plot(output$samples.theta, ylab='Range', true_params$Theta[1])
    padded_plot(output$samples.phi, ylab='Marginal Variance', true_params$Phi)
  } else {
    plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', type='l')
    plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', type='l')
    plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', type='l')
    plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', type='l')
    plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', type='l')
    plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', type='l')
    plot(output$samples.alpha.ca, ylab='Alpha (case)', type='l')
    plot(output$samples.alpha.co, ylab='Alpha (control)', type='l')
    plot(output$samples.theta, ylab='Range', type='l')
    plot(output$samples.phi, ylab='Marginal Variance', type='l')
  }
  par(mfrow=c(1,1))
  
}


plot_traces_multi <- function(output, true_params, species, output_type){
  
  if (species == 1){
    beta.case <- true_params$beta.case1
    beta.ctrl <- true_params$beta.ctrl1
    alpha.case <- true_params$Alpha.case1
    alpha.ctrl <- true_params$Alpha.ctrl1
  } else if (species == 2){
    beta.case <- true_params$beta.case2
    beta.ctrl <- true_params$beta.ctrl2
    alpha.case <- true_params$Alpha.case2
    alpha.ctrl <- true_params$Alpha.ctrl2
  }
  
  par(mfrow=c(3,4))
  if (output_type == '_multi'){
    padded_plot(output$samples.beta.ca[species,,1], ylab='Beta 0 (case)', beta.case[1])
    padded_plot(output$samples.beta.ca[species,,2], ylab='Beta 1 (case)', beta.case[2])
    padded_plot(output$samples.beta.ca[species,,3], ylab='Beta 2 (case)', beta.case[3])
    padded_plot(output$samples.beta.co[species,,1], ylab='Beta 0 (control)', beta.ctrl[1])
    padded_plot(output$samples.beta.co[species,,2], ylab='Beta 1 (control)', beta.ctrl[2])
    padded_plot(output$samples.beta.co[species,,3], ylab='Beta 2 (control)', beta.ctrl[3])
    padded_plot(output$samples.alpha.ca[species,,1], ylab='Alpha (case)', Alpha.case[1])
    padded_plot(output$samples.alpha.co[species,,1], ylab='Alpha (control)', Alpha.ctrl[1])
    padded_plot(output$samples.theta, ylab='Range', true_params$Theta[1])
    padded_plot(output$samples.phi, ylab='Marginal Variance', true_params$Phi)
  } else {
    padded_plot(output$samples.beta.ca[,1], ylab='Beta 0 (case)', beta.case[1])
    padded_plot(output$samples.beta.ca[,2], ylab='Beta 1 (case)', beta.case[2])
    padded_plot(output$samples.beta.ca[,3], ylab='Beta 2 (case)', beta.case[3])
    padded_plot(output$samples.beta.co[,1], ylab='Beta 0 (control)', beta.ctrl[1])
    padded_plot(output$samples.beta.co[,2], ylab='Beta 1 (control)', beta.ctrl[2])
    padded_plot(output$samples.beta.co[,3], ylab='Beta 2 (control)', beta.ctrl[3])
    padded_plot(output$samples.alpha.ca, ylab='Alpha (case)', Alpha.case[1])
    padded_plot(output$samples.alpha.co, ylab='Alpha (control)', Alpha.ctrl[1])
    padded_plot(output$samples.theta, ylab='Range', true_params$Theta[1])
    padded_plot(output$samples.phi, ylab='Marginal Variance', true_params$Phi)
  }
  
  par(mfrow=c(1,1))
  
}


get_output <- function(outputs, sampling, prevalence, type){
  
  output_target <- list()
  tag <- paste(sampling, prevalence, sep="_")
  for (o in outputs){
    if (grepl(tag, o$description) & grepl(type, o$description)){
      output_target <- o
      break
    }
  }
  return(output_target)
  
}


get_output_general <- function(outputs, tag){
  
  output_target <- list()
  for (o in outputs){
    if (grepl(tag, o$description)){
      output_target <- o
      break
    }
  }
  return(output_target)
  
}


#' calc_log_odds
#' 
#' calculates the estimated log odds from the preferential sampling model
#'
#' @param outputs 
#' @param sampling 
#' @param prevalence 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds <- function(outputs, sampling, prevalence){
  
  output_target <- list()
  tag <- paste(sampling, prevalence, sep="_")
  for (o in outputs){
    if (grepl(tag, o$description) & grepl('prefSampleGpCC', o$description)){
      output_target <- o
      break
    }
  }
  
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  
  w.hat <- colMeans(output_target$samples.w)
  beta_ca_h <- colMeans(output_target$samples.beta.ca)
  beta_co_h <- colMeans(output_target$samples.beta.co)
  alpha_ca_h <- mean(output_target$samples.alpha.ca)
  alpha_co_h <- mean(output_target$samples.alpha.co)
  
  lodds.ps <- x_standard %*% beta_ca_h + alpha_ca_h * w.hat - x_standard %*% beta_co_h - alpha_co_h * w.hat
  return(lodds.ps)
  
}


#' calc_log_odds_true
#' 
#' calculates the true log odds
#'
#' @param sampling 
#' @param prevalence 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds_true <- function(sampling, prevalence){
  
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  beta.case <- true_params$beta.case
  beta.ctrl <- true_params$beta.ctrl
  Alpha.case <- true_params$Alpha.case
  Alpha.ctrl <- true_params$Alpha.ctrl
  W <- true_params$W
  
  lodds.true <- x_standard %*% beta.case + Alpha.case * W - x_standard %*% beta.ctrl - Alpha.ctrl * W
  return(lodds.true)
  
}


calc_log_odds_true_general <- function(true_params){
  
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  beta.case <- true_params$beta.case
  beta.ctrl <- true_params$beta.ctrl
  Alpha.case <- true_params$Alpha.case
  Alpha.ctrl <- true_params$Alpha.ctrl
  W <- true_params$W
  
  lodds.true <- x_standard %*% beta.case + Alpha.case * W - x_standard %*% beta.ctrl - Alpha.ctrl * W
  return(lodds.true)
  
}


#' calc_log_odds_sp
#' 
#' calculates the log odds from the spatial poisson models
#'
#' @param outputs 
#' @param sampling 
#' @param prevalence 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds_sp <- function(outputs, sampling, prevalence){
  
  tag <- paste(sampling, prevalence, sep='_')
  output_ca <- list()
  output_co <- list()
  for (o in outputs){
    if (grepl(tag, o$description) & grepl('spatial_poisson_case', o$description)){
      output_ca <- o
    } else if (grepl(tag, o$description) & grepl('spatial_poisson_ctrl', o$description)){
      output_co <- o
    }
  }
  
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  location_ids <- true_params$location_ids
  x_standard <- load_x_standard(location_indicators)
  
  w.hat_spca <- colMeans(output_ca$samples.w)
  beta_ca_sp <- colMeans(output_ca$samples.beta)
  kriged_w_ca <- load_output(paste('output.krige_ca_prefSampleGpCC_', sampling, '_',prevalence, '.json', sep=''))
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, location_indicators)
  
  w.hat_spco <- colMeans(output_co$samples.w)
  beta_co_sp <- colMeans(output_co$samples.beta)
  kriged_w_co <- load_output(paste('output.krige_co_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, location_indicators)
  
  lodds <- x_standard %*% beta_ca_sp + w_ca_est - x_standard %*% beta_co_sp - w_co_est
  return(lodds)
  
}


#' calc_log_odds_pr
#' 
#' calculates the log odds from the poisson regression models
#'
#' @param sampling 
#' @param prevalence 
#'
#' @return
#' @export
#'
#' @examples
calc_log_odds_pr <- function(sampling, prevalence){
  
  betas <- load_params(paste('estimates_poisson_prefSampleGpCC_', sampling, '_', prevalence, '.json', sep=''))
  beta_ca_r <- betas$case
  beta_co_r <- betas$ctrl
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  location_indicators <- true_params$location_indicators
  x_standard <- load_x_standard(location_indicators)
  lodds <- x_standard %*% beta_ca_r - x_standard %*% beta_co_r
  return(lodds)
  
}


#' save_output
#'
#' @param output (list) output of an MCMC function
#' @param fname (character) file name
#'
#' @return writes output as json file
#'
save_output <- function(output, fname){
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  write(toJSON(output), path)
  
}


#' load_output
#'
#' @param file (character) name of output JSON file to be loaded
#'
#' @return (list) loaded output
#' @export
#'
#' @examples load_output("prefSampleGpCC_med.json")
load_output <- function(fname){
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  return(fromJSON(path))
  
}


#' save_params
#' 
#' saves parameters of the preferential sampling case control model
#'
#' @param fname (character) name of JSON file to save parameters
#'
#' @return
#' @export
#'
#' @examples
save_params_psgp <- function(fname, prior="flat"){
  
  store=list()

  if (prior != "flat"){
    store$prior_alpha_co_mean=prior_alpha_co_mean
    store$prior_alpha_co_var=prior_alpha_co_var
  }
  store$prior_phi=prior_phi
  store$prior_theta=prior_theta
  
  store$n.sample=n.sample
  store$burnin=burnin
  store$L=L_w
  store$L_ca=L_ca
  store$L_co=L_co
  store$L_a_ca=L_a_ca
  store$L_a_co=L_a_co
  store$proposal.sd.theta=proposal.sd.theta
  
  store$beta_ca_i=beta_ca_i
  store$beta_co_i=beta_co_i
  store$alpha_ca_i=alpha_ca_i
  store$alpha_co_i=alpha_co_i
  store$theta_i=theta_i
  store$phi_i=phi_i
  store$w_i=w_i
  
  store$m_aca=m_aca
  store$m_aco=m_aco
  store$m_ca <- m_ca
  store$m_co <- m_co
  store$m_w <- m_w
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  write(toJSON(store), path)
  
}


#' save_params_psc
#' 
#' saves parameters of the spatial poisson regression case model
#'
#' @param fname (character) file name
#'
#' @return
#' @export
#'
#' @examples
save_params_psc <- function(fname){
  
  store=list()
  
  store$prior_phi_=prior_phi_
  store$prior_theta_=prior_theta_
  
  store$n.sample_=n.sample_
  store$burnin_=burnin_
  store$L_w_=L_w_
  store$L_b_=L_b_
  
  store$beta_ca_i_=beta_ca_i_
  store$theta_i_=theta_i_
  store$phi_i_=phi_i_
  store$w_i_=w_i_
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  write(toJSON(store), path)
  
}


#' save_params_psco
#' 
#' saves parameters of the spatial poisson regression control model
#'
#' @param fname fname (character) file name
#'
#' @return
#' @export
#'
#' @examples
save_params_psco <- function(fname){
  
  store=list()
  
  store$prior_phi__=prior_phi__
  store$prior_theta__=prior_theta__
  
  store$n.sample__=n.sample__
  store$burnin__=burnin__
  store$L_w__=L_w__
  store$L_b__=L_b__
  
  store$beta_co_i__=beta_co_i__
  store$theta_i__=theta_i__
  store$phi_i__=phi_i__
  store$w_i__=w_i__
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  write(toJSON(store), path)
  
}


#' save_estimates_pr
#' 
#' save parameter estimates of 2 (nonspatial) poisson models
#'
#' @param beta_ca (numeric) estimates
#' @param beta_ctrl (numeric) estimates
#' @param fname (character) file name
#'
#' @return
#' @export
#'
#' @examples
save_estimates_pr <- function(beta_ca, beta_ctrl, fname){
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  store <- list(
    case=beta_ca_r,
    ctrl=beta_co_r,
    description="Poisson Regression"
  )
  write(toJSON(store), path)
  
}


#' load_params
#'
#' @param fname (character) file name
#'
#' @return
#' @export
#'
#' @examples
load_params <- function(fname){
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  return(fromJSON(path))
  
}


save_true_params <- function(sampling, prevalence){
  
  fname <- paste('true_params_', sampling, '_', prevalence, '.json', sep='')
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  store <- list()
  store$Alpha.case <- Alpha.case
  store$Alpha.ctrl <- Alpha.ctrl
  store$beta.case <- beta.case
  store$beta.ctrl <- beta.ctrl
  store$W <- W
  store$Theta <- Theta
  store$Phi <- Phi
  store$location_indicators <- as.logical(locs$status)
  store$location_ids <- locs$ids
  write(toJSON(store), path)
  
  
}


save_true_params_multi <- function(tag){
  
  fname <- paste('true_params_', tag, '.json', sep='')
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  store <- list()
  store$Alpha.case1 <- Alpha.case1
  store$Alpha.ctrl1 <- Alpha.ctrl1
  store$Alpha.case2 <- Alpha.case2
  store$Alpha.ctrl2 <- Alpha.ctrl2
  store$beta.case1 <- beta.case1
  store$beta.ctrl1 <- beta.ctrl1
  store$beta.case2 <- beta.case2
  store$beta.ctrl2 <- beta.ctrl2
  store$W1 <- W1
  store$W2 <- W2
  store$Theta <- Theta
  store$Phi <- Phi
  store$data <- data
  write(toJSON(store), path)
  
}


save_true_params_general <- function(tag){
  
  fname <- paste('true_params_', tag, '.json', sep='')
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  store <- list()
  store$Alpha.case <- Alpha.case
  store$Alpha.ctrl <- Alpha.ctrl
  store$beta.case <- beta.case
  store$beta.ctrl <- beta.ctrl
  store$W <- W
  store$Theta <- Theta
  store$Phi <- Phi
  store$location_indicators <- as.logical(locs$status)
  store$location_ids <- locs$ids
  write(toJSON(store), path)
  
  
}



load_x_standard <- function(location_indicators){
  
  caPr <- load_prism_pcs()
  caPr.disc <- aggregate(caPr, fact=8)
  
  x_1 <- caPr.disc[[1]][]
  x_1 <- x_1[][!is.na(x_1[])]
  x_2 <- caPr.disc[[2]][]
  x_2 <- x_2[][!is.na(x_2[])]
  mu_1 <- mean(x_1[location_indicators])
  sd_1 <- sd(x_1[location_indicators])
  mu_2 <- mean(x_2[location_indicators])
  sd_2 <- sd(x_2[location_indicators])
  
  x_1_std <- (x_1 - mu_1)/sd_1
  x_2_std <- (x_2 - mu_2)/sd_2
  x_std <- array(1, c(length(x_1), 1))
  x_std <- cbind(x_std, x_1_std, x_2_std)
  return(x_std)

}


load_x_ca <- function(factor=NULL){
  
  caPr <- load_prism_pcs()
  if (!is.null(factor)){
    caPr <- aggregate(caPr, fact=factor)
    x_1 <- caPr[[1]][]
    x_1 <- x_1[][!is.na(x_1[])]
    x_2 <- caPr.disc[[2]][]
    x_2 <- x_2[][!is.na(x_2[])]
  } else {
    pc1 <- mask(caPr[[1]], caPr[[2]])
    pc2 <- caPr[[2]]
    x_1 <- pc1[]
    x_1 <- x_1[][!is.na(x_1[])]
    x_2 <- pc2[]
    x_2 <- x_2[][!is.na(x_2[])]
  }

  mu_1 <- mean(x_1)
  sd_1 <- sd(x_1)
  mu_2 <- mean(x_2)
  sd_2 <- sd(x_2)
  
  x_1_std <- (x_1 - mu_1)/sd_1
  x_2_std <- (x_2 - mu_2)/sd_2
  x_std <- array(1, c(length(x_1), 1))
  x_std <- cbind(x_std, x_1_std, x_2_std)
  return(x_std)
  
}



# load mcmc outputs
load_sim_outputs <- function(tag=''){
  
  output_list <- list()
  
  counter <- 1
  for (f in list.files('/Users/brianconroy/Documents/research/dataInt/output/')){
    if (grepl('output', f) & !grepl('krige', f) & grepl(tag, f)) {
      output_list[[counter]] <- load_output(f)
      counter <- counter + 1
    }
  }
  
  return(output_list)
  
}
