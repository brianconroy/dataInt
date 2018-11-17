library(jsonlite)


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
