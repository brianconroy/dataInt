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
#' @param fname (character) name of JSON file to save parameters
#'
#' @return
#' @export
#'
#' @examples
save_params <- function(fname){
  
  store=list()
  store$Theta=Theta
  store$Phi=Phi
  store$Alpha.case=Alpha.case
  store$beta.case=beta.case
  store$prior_alpha_ca_mean=prior_alpha_ca_mean
  store$prior_alpha_ca_var=prior_alpha_ca_var
  store$Alpha.ctrl=Alpha.ctrl
  store$beta.ctrl=beta.ctrl
  
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
  
  store$beta_ca_i=beta_ca_i
  store$beta_co_i=beta_co_i
  store$alpha_ca_i=alpha_ca_i
  store$alpha_co_i=alpha_co_i
  store$theta_i=theta_i
  store$phi_i=phi_i
  store$w_i=w_i
  
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  write(toJSON(store), path)
  
}
