library(jsonlite)

# prevalence # param 	# model		# estimate	# true	# bias

table_params <- function(outputs, sampling, prevalence){
  
  output_ps <- get_output(outputs, sampling, prevalence, 'prefSampleGpCC')
  output_sp_ca <- get_output(outputs, sampling, prevalence, 'spatial_poisson_case')
  output_sp_co <- get_output(outputs, sampling, prevalence, 'spatial_poisson_ctrl')
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  
  rows <- list()
  
  list(
    prevalence=prevalence,
    sampling=sampling,
    model='PS',
    parameter='Beta 0 (case)',
    estimate=get_estimate(output, parameter),
    true=true_params$beta.case[1]
  )
  
}


make_row <- function(prevalence, sampling, model, parameter, output){
  
  est <- get_estimate(output, parameter)
  
  
}


get_estimate <- function(output, parameter){
  
  if (grepl("prefSampleGpCC", output$description)){
    type <- "PS"
  } else {
    type <- "SP"
  }
  
  if (parameter == "Beta 0 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,1]
    } else {
      target_samples <- output$samples.beta[,1]
    }
  } else if (parameter == "Beta 1 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,2]
    } else {
      target_samples <- output$samples.beta[,2]
    }
  } else if (parameter == "Beta 2 (case)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.ca[,3]
    } else {
      target_samples <- output$samples.beta[,3]
    }
  } else if (parameter == "Beta 0 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,1]
    } else {
      target_samples <- output$samples.beta[,1]
    }
  } else if (parameter == "Beta 1 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,2]
    } else {
      target_samples <- output$samples.beta[,2]
    }
  } else if (parameter == "Beta 2 (control)"){
    if (type == "PS"){
      target_samples <- output$samples.beta.co[,3]
    } else {
      target_samples <- output$samples.beta[,3]
    }
  }
  
  est <- mean(target_samples)
  return(round(est, 3))
  
}


plot_traces <- function(outputs, sampling, prevalence){
  
  output <- get_output(outputs, sampling, prevalence)
  true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
  
  par(mfrow=c(3,4))
  plot(output$samples.beta.ca[,1], typ='l', ylab='Beta 0 (case)'); abline(h=true_params$beta.case[1], col='2')
  plot(output$samples.beta.ca[,2], typ='l', ylab='Beta 1 (case)'); abline(h=true_params$beta.case[2], col='2')
  plot(output$samples.beta.ca[,3], typ='l', ylab='Beta 2 (case)'); abline(h=true_params$beta.case[3], col='2')
  
  plot(output$samples.beta.co[,1], typ='l', ylab='Beta 0 (control)'); abline(h=true_params$beta.ctrl[1], col='2')
  plot(output$samples.beta.co[,2], typ='l', ylab='Beta 1 (control)'); abline(h=true_params$beta.ctrl[2], col='2')
  plot(output$samples.beta.co[,3], typ='l', ylab='Beta 2 (control)'); abline(h=true_params$beta.ctrl[3], col='2')
  
  plot(output$samples.alpha.ca, typ='l', ylab='Alpha (case)'); abline(h=true_params$Alpha.case[1], col='2')
  plot(output$samples.alpha.co, typ='l', ylab='Alpha (control)'); abline(h=true_params$Alpha.ctrl[1], col='2')
  
  plot(output$samples.theta, typ='l', ylab='Range'); abline(h=true_params$Theta, col='2')
  plot(output$samples.phi, typ='l', ylab='Marginal Variance'); abline(h=true_params$Phi, col='2')
  par(mfrow=c(1,1))
  
}


get_output <- function(outputs, sampling, prevalence, type){
  
  output_target <- list()
  for (o in outputs){
    if (grepl(sampling, o$description) & grepl(prevalence, o$description)
        & grepl(type, o$description)){
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
  for (o in outputs){
    if (grepl(sampling, o$description) & grepl(prevalence, o$description) 
        & grepl('prefSampleGpCC', o$description)){
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
  alpha_ca_h <- colMeans(output_target$samples.alpha.ca)
  alpha_co_h <- colMeans(output_target$samples.alpha.co)
  
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
  
  output_ca <- list()
  output_co <- list()
  for (o in outputs){
    if (grepl(sampling, o$description) & grepl(prevalence, o$description) 
        & grepl('spatial_poisson_case', o$description)){
      output_ca <- o
    } else if (grepl(sampling, o$description) & grepl(prevalence, o$description) 
               & grepl('spatial_poisson_ctrl', o$description)){
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
save_params_psgp <- function(fname){
  
  store=list()

  store$prior_alpha_co_mean=prior_alpha_co_mean
  store$prior_alpha_co_var=prior_alpha_co_var
  store$prior_phi=prior_phi
  store$prior_theta=prior_theta
  
  store$n.sample=n.sample
  store$burnin=burnin
  store$L=L
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
    ctrl=beta_co_r
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


# load mcmc outputs
load_sim_outputs <- function(){
  
  output_list <- list()
  
  counter <- 1
  for (f in list.files('/Users/brianconroy/Documents/research/dataInt/output/')){
    if (grepl('output', f) & !grepl('krige', f)) {
      output_list[[counter]] <- load_output(f)
      counter <- counter + 1
    }
  }
  
  return(output_list)
  
}
