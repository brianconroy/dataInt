

prefSampleTemporal <- function(data, d, n.sample, burnin, 
                               L_w, L_ca, L_co, L_a_ca, L_a_co,
                               proposal.sd.theta=0.3,
                               m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, m_u=700,
                               target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, target_u=0.75,
                               self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE, self_tune_u=TRUE,
                               delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, delta_u=NULL,
                               beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                               theta_initial=NULL, phi_initial=NULL, w_initial=NULL, u_initial=NULL,
                               prior_phi, prior_theta, prior_alpha_ca_var, prior_alpha_co_var, prior_u_mean, prior_u_var){

  
  ## setup
  case.data_time <- data$case.data
  ctrl.data_time <- data$ctrl.data
  locs_time <- data$locs
  p <- ncol(case.data_time[[1]]$x)
  N.w <- length(locs_time[[1]]$status)
  N.t <- length(locs_time)
  
  
  ## starting values
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(beta_ca_initial)){
    beta.ca <- rnorm(p)
  } else {
    beta.ca <- beta_ca_initial
  }
  if (is.null(beta_co_initial)){
    beta.co <- rnorm(p)
  } else {
    beta.co <- beta_co_initial
  }
  if (is.null(alpha_ca_initial)){
    alpha.ca.i <- runif(1, 1, 3)
  } else {
    alpha.ca.i <- alpha_ca_initial
  }
  if (is.null(alpha_co_initial)){
    alpha.co.i <- runif(1, -3, -1)
  } else {
    alpha.co.i <- alpha_co_initial
  }
  if (is.null(theta_initial)){
    theta.i <- runif(1, 5, 7)
  } else {
    theta.i <- theta_initial
  }
  if (is.null(phi_initial)){
    phi.i <- runif(1, 3.5, 4.5)
  } else {
    phi.i <- phi_initial
  }
  if (is.null(u_initial)){
    u.i <- rnorm(N.t)
  } else {
    u.i <- u_initial
  }
  

  # storage
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w))
  samples.u <- array(NA, c(n.keep, N.t))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.alpha.ca <- array(NA, c(n.keep, 1))
  samples.alpha.co <- array(NA, c(n.keep, 1))
  samples.beta.ca <- array(NA, c(n.keep, length(beta.ca)))
  samples.beta.co <- array(NA, c(n.keep, length(beta.co)))
  
  
  if (self_tune_w){
    w_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_aca){
    a_ca_tuning <- initialize_tuning(m=m_aca, target=target_aca)
  } else {
    a_ca_tuning <- list(delta_curr=delta_aca)
  }
  if (self_tune_aco){
    a_co_tuning <- initialize_tuning(m=m_aco, target=target_aco)
  } else {
    a_co_tuning <- list(delta_curr=delta_aco)
  }
  if (self_tune_ca){
    ca_tuning <- initialize_tuning(m=m_ca, target=target_ca)
  } else {
    ca_tuning <- list(delta_curr=delta_ca)
  }
  if (self_tune_co){
    co_tuning <- initialize_tuning(m=m_co, target=target_co)
  } else {
    co_tuning <- list(delta_curr=delta_co)
  }
  if (self_tune_u){
    u_tuning <- list()
    for (i in 1:N.t){
      u_tuning[[i]] <- initialize_tuning(m=m_u, target=target_u)
    }
  } else {
    u_tuning <- list(delta_curr=delta_u)
  }
  
  deltas_w <- c()
  deltas_ca <- c()
  deltas_co <- c()
  deltas_aca <- c()
  deltas_aco <- c()
  deltas_u <- array(NA, c(n.sample, N.t))
  
  accept <- list(
    w=0,
    beta_ca=0,
    beta_co=0,
    alpha_ca=0,
    alpha_co=0,
    u=rep(0, N.t),
    theta=0
  )
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateTemporal(case.data_time, ctrl.data_time, locs_time, alpha.ca.i, beta.ca,
                                  alpha.co.i, beta.co, w.i, u.i, sigma.i, sigma.inv.i, w_tuning$delta_curr, L_w, offset=0)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N.w/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])

    ## sample from beta case
    beta.out.ca <- betaHmcUpdateTemporal(case.data_time, locs_time, w.i, u.i, beta.ca, alpha.ca.i, ca_tuning$delta_curr, L_ca, offset=0)
    beta.ca <- beta.out.ca$beta

    ## sample from alpha case
    alpha.out.ca <- alphaHmcUpdateTemporal(case.data_time, locs_time, w.i, u.i, beta.ca, alpha.ca.i,
                                           a_ca_tuning$delta_curr, prior_alpha_ca_mean, prior_alpha_ca_var, L_a_ca, offset=0)
    alpha.ca.i <- alpha.out.ca$alpha
    
    ## sample from beta control
    beta.out.co <- betaHmcUpdateTemporal(ctrl.data_time, locs_time, w.i, u.i, beta.co, alpha.co.i, co_tuning$delta_curr, L_co, offset=0)
    beta.co <- beta.out.co$beta

    ## sample from alpha control
    alpha.out.co <- alphaHmcUpdateTemporal(ctrl.data_time, locs_time, w.i, u.i, beta.co, alpha.co.i, a_co_tuning$delta_curr, prior_alpha_co_mean, prior_alpha_co_var, L_a_co, offset=0)
    alpha.co.i <- alpha.out.co$alpha
    
    ## sample from u
    u.accept <- rep(0, N.t)
    u.alphas <- rep(0, N.t)
    for(t in 1:N.t){
      case.data_t <- case.data_time[[t]]
      ctrl.data_t <- ctrl.data_time[[t]]
      locs_t <- locs_time[[t]]
      u_tuning_t <- u_tuning[[t]]
      u_t <- u.i[t]
      u.out.t <- uHmcUpdate(case.data_t, ctrl.data_t, locs_t, alpha.ca.i, beta.ca, alpha.co.i, beta.co, w.i, u_t, u_tuning_t$delta_curr, L_u, prior_u_mean, prior_u_var)
      u.i[t] <- u.out.t$u
      u.accept[t] <- u.out.t$accept
      u.alphas[t] <- u.out.t$a
    }
    
    if (i > burnin){
      
      j <- i - burnin
      
      samples.beta.ca[j,] <- beta.ca
      samples.beta.co[j,] <- beta.co
      samples.alpha.ca[j,] <- alpha.ca.i
      samples.alpha.co[j,] <- alpha.co.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      samples.u[j,] <- u.i
      
      accept$w <- accept$w + w.out.i$accept
      accept$theta <- accept$theta + theta.out$accept
      accept$beta_ca <- accept$beta_ca + beta.out.ca$accept
      accept$beta_co <- accept$beta_co + beta.out.co$accept
      accept$alpha_ca <- accept$alpha_ca + alpha.out.ca$accept
      accept$alpha_co <- accept$alpha_co + alpha.out.co$accept
      accept$u <- accept$u + u.accept
      
    }
    
    if (self_tune_w){
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    if (self_tune_aca){
      a_ca_tuning <- update_tuning(a_ca_tuning, alpha.out.ca$a, i, alpha.out.ca$accept)
      deltas_aca <- c(deltas_aca, a_ca_tuning$delta_curr)
    }
    if (self_tune_aco){
      a_co_tuning <- update_tuning(a_co_tuning, alpha.out.co$a, i, alpha.out.co$accept)
      deltas_aco <- c(deltas_aco, a_co_tuning$delta_curr)
    }
    if (self_tune_ca){
      ca_tuning <- update_tuning(ca_tuning, beta.out.ca$a, i, beta.out.ca$accept)
      deltas_ca <- c(deltas_ca, ca_tuning$delta_curr)
    }
    if (self_tune_co){
      co_tuning <- update_tuning(co_tuning, beta.out.co$a, i, beta.out.co$accept)
      deltas_co <- c(deltas_co, co_tuning$delta_curr)
    }
    if (self_tune_u){
      for (t in 1:N.t){
        u_tuning[[t]] <- update_tuning(u_tuning[[t]], u.alphas[t], i, u.accept[t])
        deltas_u[i,t] <- u_tuning[[t]]$delta_curr
      }
    }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  for (n in names(accept)){
    accept[[n]] <- accept[[n]]/n.keep
  }
  
  output <- list()
  output$accept <- accept
  output$samples.beta.ca <- samples.beta.ca
  output$samples.beta.co <- samples.beta.co
  output$samples.alpha.ca <- samples.alpha.ca
  output$samples.alpha.co <- samples.alpha.co
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$samples.u <- samples.u
  output$deltas_w <- deltas_w
  output$deltas_aca <- deltas_aca
  output$deltas_aco <- deltas_aco
  output$deltas_ca <- deltas_ca
  output$deltas_co <- deltas_co
  output$deltas_u <- deltas_u
  output$L_w <- L_w
  output$L_ca <- L_ca
  output$L_co <- L_co
  output$L_a_ca <- L_a_ca 
  output$L_a_co <- L_a_co
  output$L_u <- L_u
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$prior_u_mean <- prior_u_mean
  output$prior_u_var <- prior_u_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}


burnin_after <- function(output, n.burn){
  
  n.curr <- output$n.sample - output$burnin
  i.start <- n.burn + 1
  output$burnin <- output$burnin + n.burn
  if (length(dim(output$samples.alpha.ca)) == 2 | class(output$samples.alpha.ca) == 'numeric'){
    output$samples.alpha.ca <- output$samples.alpha.ca[i.start:n.curr]
    output$samples.alpha.co <- output$samples.alpha.co[i.start:n.curr]
    output$samples.beta.ca <- output$samples.beta.ca[i.start:n.curr,]
    output$samples.beta.co <- output$samples.beta.co[i.start:n.curr,]
  } else {
    n_s <- dim(output$samples.alpha.ca)[1]
    new.aca <- array(NA, c(n_s,n.curr-i.start+1,1))
    new.aco <- array(NA, c(n_s,n.curr-i.start+1,1))
    for (j in 1:n_s){
      new.aca[j,,] <- matrix(output$samples.alpha.ca[j,i.start:n.curr,])
      new.aco[j,,] <- matrix(output$samples.alpha.co[j,i.start:n.curr,])
    }
    output$samples.alpha.ca <- new.aca
    output$samples.alpha.co <- new.aco
    output$samples.beta.ca <- output$samples.beta.ca[,i.start:n.curr,]
    output$samples.beta.co <- output$samples.beta.co[,i.start:n.curr,]
  }
  output$samples.w <- output$samples.w[i.start:n.curr,]
  output$samples.phi <- output$samples.phi[i.start:n.curr]
  output$samples.theta <- output$samples.theta[i.start:n.curr]
  return(output)
  
}


#' continueMCMC
#' 
#' continues running an MCMC chain from the output of prefSampleGpCC
#'
#' @param output (list) output of prefSampleGpCC
#' @param n.sample (numeric) number of additional samples to generate
#'
#' @return
#' @export
#'
#' @examples
continueMCMC <- function(data, D, output, n.sample){
  
  # get initial values
  n.sample.old <- nrow(output$samples.beta.ca)
  beta_ca_initial <- output$samples.beta.ca[n.sample.old,]
  beta_co_initial <- output$samples.beta.co[n.sample.old,]
  alpha_ca_initial <- output$samples.alpha.ca[n.sample.old]
  alpha_co_initial <- output$samples.alpha.co[n.sample.old]

  w_initial <- output$samples.w[n.sample.old,]
  theta_initial <- output$samples.theta[n.sample.old]
  phi_initial <- output$samples.phi[n.sample.old]
  
  # get tuning parameters
  delta_w <- tail(output$deltas_w, 1)
  delta_aca <- tail(output$deltas_aca, 1)
  delta_aco <- tail(output$deltas_aco, 1)
  delta_ca <- tail(output$deltas_ca, 1)
  delta_co <- tail(output$deltas_co, 1)
  L_w <- output$L_w
  L_ca <- output$L_ca
  L_co <- output$L_co
  L_a_ca <- output$L_a_ca 
  L_a_co <- output$L_a_co
  proposal.sd.theta <- output$proposal.sd.theta
  
  # get priors
  prior_phi <- output$prior_phi
  prior_theta <- output$prior_theta
  prior_alpha_ca_var <- output$prior_alpha_ca_var
  prior_alpha_co_var <- output$prior_alpha_co_var
  more_output <- prefSampleGpCC(data, D, n.sample, burnin=0, 
                                  L_w=L_w, L_ca=L_ca, L_co=L_co, L_a_ca=L_a_ca, L_a_co=L_a_co,
                                  proposal.sd.theta=proposal.sd.theta,
                                  self_tune_w=FALSE, self_tune_aca=FALSE, self_tune_aco=FALSE, self_tune_ca=FALSE, self_tune_co=FALSE,
                                  delta_w=delta_w, delta_aca=delta_aca, delta_aco=delta_aco, delta_ca=delta_ca, delta_co=delta_co, 
                                  beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, 
                                  alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                                  theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                                  prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_var=prior_alpha_ca_var, 
                                  prior_alpha_co_var=prior_alpha_ca_var)
  
  # combine outputs
  new_output <- output
  new_output$samples.alpha.ca <- c(new_output$samples.alpha.ca, more_output$samples.alpha.ca)
  new_output$samples.alpha.co <- c(new_output$samples.alpha.co, more_output$samples.alpha.co)
  new_output$samples.beta.ca <- rbind(new_output$samples.beta.ca, more_output$samples.beta.ca)
  new_output$samples.beta.co <- rbind(new_output$samples.beta.co, more_output$samples.beta.co)
  new_output$samples.phi <- c(new_output$samples.phi, more_output$samples.phi)
  new_output$samples.theta <- c(new_output$samples.theta, more_output$samples.theta)
  new_output$samples.w <- rbind(new_output$samples.w, more_output$samples.w)
  new_output$deltas_aca <- c(new_output$deltas_aca, more_output$deltas_aca)
  new_output$deltas_aco <- c(new_output$deltas_aco, more_output$deltas_aco)
  new_output$deltas_co <- c(new_output$deltas_co, more_output$deltas_co)
  new_output$deltas_ca <- c(new_output$deltas_ca, more_output$deltas_ca)
  new_output$deltas_w <- c(new_output$deltas_w, more_output$deltas_w)
  new_output$n.sample <- new_output$n.sample + n.sample
  new_output$accept <- (new_output$n.sample * new_output$accept + n.sample * more_output$accept)/(new_output$n.sample + n.sample)
  
  return(new_output)
  
}
