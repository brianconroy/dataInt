

summarize_ps_params <- function(output){
  
  params <- list()
  params[[1]] <- list(
    Name='Beta 0 (case)',
    Posterior_Mean=round(mean(output$samples.beta.ca[,1]), 3),
    Posterior_Median=round(median(output$samples.beta.ca[,1]), 3),
    Posterior_Variance=round(var(output$samples.beta.ca[,1]), 3)
  )
  params[[2]] <- list(
    Name='Beta 1 (case)',
    Posterior_Mean=round(mean(output$samples.beta.ca[,2]), 3),
    Posterior_Median=round(median(output$samples.beta.ca[,2]), 3),
    Posterior_Variance=round(var(output$samples.beta.ca[,2]), 3)
  )
  params[[3]] <- list(
    Name='Beta 2 (case)',
    Posterior_Mean=round(mean(output$samples.beta.ca[,3]), 3),
    Posterior_Median=round(median(output$samples.beta.ca[,3]), 3),
    Posterior_Variance=round(var(output$samples.beta.ca[,3]), 3)
  )
  params[[4]] <- list(
    Name='Beta 0 (control)',
    Posterior_Mean=round(mean(output$samples.beta.co[,1]), 3),
    Posterior_Median=round(median(output$samples.beta.co[,1]), 3),
    Posterior_Variance=round(var(output$samples.beta.co[,1]), 3)
  )
  params[[5]] <- list(
    Name='Beta 1 (control)',
    Posterior_Mean=round(mean(output$samples.beta.co[,2]), 3),
    Posterior_Median=round(median(output$samples.beta.co[,2]), 3),
    Posterior_Variance=round(var(output$samples.beta.co[,2]), 3)
  )
  params[[6]] <- list(
    Name='Beta 2 (control)',
    Posterior_Mean=round(mean(output$samples.beta.co[,3]), 3),
    Posterior_Median=round(median(output$samples.beta.co[,3]), 3),
    Posterior_Variance=round(var(output$samples.beta.co[,3]), 3)
  )
  params[[7]] <- list(
    Name='Alpha (case)',
    Posterior_Mean=round(mean(output$samples.alpha.ca), 3),
    Posterior_Median=round(median(output$samples.alpha.ca), 3),
    Posterior_Variance=round(var(output$samples.alpha.ca), 3)
  )
  params[[8]] <- list(
    Name='Alpha (control)',
    Posterior_Mean=round(mean(output$samples.alpha.co), 3),
    Posterior_Median=round(median(output$samples.alpha.co), 3),
    Posterior_Variance=round(var(output$samples.alpha.co), 3)
  )
  params[[9]] <- list(
    Name='Range',
    Posterior_Mean=round(mean(output$samples.theta), 3),
    Posterior_Median=round(median(output$samples.theta), 3),
    Posterior_Variance=round(var(output$samples.theta), 3)
  )
  params[[10]] <- list(
    Name='Marginal Variance',
    Posterior_Mean=round(mean(output$samples.phi), 3),
    Posterior_Median=round(median(output$samples.phi), 3),
    Posterior_Variance=round(var(output$samples.phi), 3)
  )
  return(ldply(params, 'data.frame'))
  
}


#' compare_params
#'
#' @param beta.ca.hat_p (numeric) vector of parameter estimates from case poisson model
#' @param beta.co.hat_p (numeric) vector of parameter estimates from control poisson model
#' @param beta.ca.hat (numeric) vector of parameter estimates from case preferential sampling model
#' @param beta.co.hat (numeric) vector of parameter estimates from control preferential sampling model
#'
#' @return
#' @export
#'
#' @examples
compare_params <- function(beta.ca.hat_p, beta.co.hat_p, beta.ca.hat, beta.co.hat){
  
  param_comp[[1]] <- list(
    Parameter='Beta 0 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[1], 3)
  )
  param_comp[[2]] <- list(
    Parameter='Beta 0 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[1], 3)
  )
  param_comp[[3]] <- list(
    Parameter='Beta 1 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[2], 3)
  )
  param_comp[[4]] <- list(
    Parameter='Beta 1 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[2], 3)
  )
  param_comp[[5]] <- list(
    Parameter='Beta 2 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[3], 3)
  )
  param_comp[[6]] <- list(
    Parameter='Beta 2 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[3], 3)
  )
  
  param_comp[[7]] <- list(
    Parameter='Beta 0 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[1], 3)
  )
  param_comp[[8]] <- list(
    Parameter='Beta 0 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[1], 3)
  )
  param_comp[[9]] <- list(
    Parameter='Beta 1 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[2], 3)
  )
  param_comp[[10]] <- list(
    Parameter='Beta 1 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[2], 3)
  )
  param_comp[[11]] <- list(
    Parameter='Beta 2 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[3], 3)
  )
  param_comp[[12]] <- list(
    Parameter='Beta 2 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[3], 3)
  )
  return(ldply(param_comp, 'data.frame'))
  
}
