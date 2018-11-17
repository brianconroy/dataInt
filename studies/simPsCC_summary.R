############################
# Summarizes the preferential
# sampling x prevalence
# simulation results
############################

library(plyr)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)


outputs <- load_sim_outputs()


sampling <- "none"
prevalence <- "low"


calc_log_odds <- function(outputs, sampling, prevalence){
  
  output_target <- list()
  for (o in outputs){
    if (grepl(sampling, o$description) & grepl(prevalence, o$description)){
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


lodds <- calc_log_odds(outputs, sampling, prevalence)
lodds_true <- calc_log_odds_true(sampling, prevalence)
lodds_sp <- calc_log_odds_sp(outputs, sampling, prevalence)
lodds_pr <- calc_log_odds_pr(sampling, prevalence)

plot(x=lodds_true, y=lodds); abline(0, 1, col=2)
plot(x=lodds_sp, y=lodds); abline(0, 1, col=2)
plot(x=lodds_pr, y=lodds); abline(0, 1, col=2)


# log odds
  # scatterplots
  # box & whisker plots
  # maps

# param traceplots
# param tables
# plot of param biases
