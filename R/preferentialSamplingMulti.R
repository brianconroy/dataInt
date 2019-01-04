

prefSampleMulti_1 <- function(data, n.sample, burnin, 
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=2000, m_aco=2000, m_ca=700, m_co=700, m_w=700, 
                           target_aca=0.75, target_aco=0.75, target_ca=0.75, target_co=0.75, target_w=0.75, 
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                           delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                           beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                           theta_initial=NULL, phi_initial=NULL, w_initial=NULL,
                           prior_phi, prior_theta, prior_alpha_ca_var, prior_alpha_co_var){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
  X.c <- case.data[[1]]$x.standardised
  locs <- data$locs
  N.w <- length(locs[[1]]$status)
  N.d <- length(locs)
  N.p <- ncol(case.data[[1]]$x.standardised)

  ## starting values
  if (is.null(w_initial)){
    w.i <- rnorm(N.w)
  } else {
    w.i <- w_initial
  }
  if (is.null(beta_ca_initial)){
    beta.ca <- list()
    for (i in 1:N.d){
      beta.ca[[i]] <- rnorm(N.p)
    }
  } else {
    beta.ca <- beta_ca_initial
  }
  if (is.null(beta_co_initial)){
    beta.co <- list()
    for (i in 1:N.d){
      beta.co[[i]] <- rnorm(N.p)
    }
  } else {
    beta.co <- beta_co_initial
  }
  if (is.null(alpha_ca_initial)){
    alpha.ca.i <- list()
    for (i in 1:N.d){
      alpha.ca.i[[i]] <- runif(1, 1, 3)
    }
  } else {
    alpha.ca.i <- alpha_ca_initial
  }
  if (is.null(alpha_co_initial)){
    alpha.co.i <- list()
    for (i in 1:N.d){
      alpha.co.i[[i]] <- runif(1, -3, -1)
    }
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

  
  # storage
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.phi <- array(NA, c(n.keep, 1))
  samples.alpha.ca <- array(NA, c(N.d, n.keep, 1))
  samples.alpha.co <- array(NA, c(N.d, n.keep, 1))
  samples.beta.ca <- array(NA, c(N.d, n.keep, N.p))
  samples.beta.co <- array(NA, c(N.d, n.keep, N.p))
  
  
  if (self_tune_w){
    w_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_aca){
    a_ca_tuning <- list()
    for (i in 1:N.d){
      a_ca_tuning[[i]] <- initialize_tuning(m=m_aca, target=target_aca)
    }
  } else {
    a_ca_tuning <- list()
    for (i in 1:N.d){
      a_ca_tuning[[i]] <- list(delta_curr=delta_aca)
    }
  }
  if (self_tune_aco){
    a_co_tuning <- list()
    for (i in 1:N.d){
      a_co_tuning[[i]] <- initialize_tuning(m=m_aco, target=target_aco)
    }
  } else {
    a_co_tuning <- list()
    for (i in 1:N.d){
      a_co_tuning[[i]] <- list(delta_curr=delta_aco)
    }
  }
  if (self_tune_ca){
    ca_tuning <- list()
    for (i in 1:N.d){
      ca_tuning[[i]] <- initialize_tuning(m=m_ca, target=target_ca)
    }
  } else {
    ca_tuning <- list()
    for (i in 1:N.d){
      ca_tuning[[i]] <- list(delta_curr=delta_ca)
    }
  }
  if (self_tune_co){
    co_tuning <- list()
    for (i in 1:N.d){
      co_tuning[[i]] <- initialize_tuning(m=m_co, target=target_co)
    }
  } else {
    co_tuning <- list()
    for (i in 1:N.d){
      co_tuning[[i]] <- list(delta_curr=delta_co)
    }
  }
  
  deltas_w <- c()
  deltas_ca <- array(NA, c(n.sample, N.p))
  deltas_co <- array(NA, c(n.sample, N.p))
  deltas_aca <- array(NA, c(n.sample, 1))
  deltas_aco <- array(NA, c(n.sample, 1))
  
  prior.mean.beta <- rep(0, p.c) # ToDo change
  prior.var.beta <- rep(1000, p.c)
  
  accept <- rep(0, 6) # ToDo change
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from w
    sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
    sigma.inv.i <- solve(sigma.i)
    w.out.i <- wHmcUpdateMulti(case.data, ctrl.data, alpha.ca.i, beta.ca,
                            alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w
    
    ## sample from theta
    theta.out <- rangeMhUpdate(theta.i, as.numeric(w.i), d, phi.i, proposal.sd.theta, a=prior_theta[1], b=prior_theta[2])
    theta.i <- theta.out$theta
    
    ## sample from phi
    R.i <- sigma.i/phi.i
    phi.i <- 1/rgamma(1, N/2 + prior_phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior_phi[2])
    
    ## sample from beta.case
    # w.i.sub <- w.i[locs$ids]
    # beta.out.ca <- caseHmcUpdate(Y.ca, w.i[locs$ids], X.c, beta.ca, alpha.ca.i, ca_tuning$delta_curr, L_ca)
    # beta.ca <- beta.out.ca$beta
    # 
    # ## sample from alpha case
    # alpha.out.ca <- alphaHmcUpdate(Y.ca, w.i.sub, X.c, beta.ca, alpha.ca.i, 
    #                                a_ca_tuning$delta_curr, prior_alpha_ca_mean, prior_alpha_ca_var, L_a_ca)
    # alpha.ca.i <- alpha.out.ca$alpha
    # 
    # ## sample from beta.ctrl
    # beta.out.co <- caseHmcUpdate(Y.co, w.i[locs$ids], X.c, beta.co, alpha.co.i, co_tuning$delta_curr, L_co)
    # beta.co <- beta.out.co$beta
    # 
    # ## sample from alpha control
    # alpha.out.co <- alphaHmcUpdate(Y.co, w.i.sub, X.c, beta.co, alpha.co.i, 
    #                                a_co_tuning$delta_curr, prior_alpha_co_mean, prior_alpha_co_var, L_a_co)
    # alpha.co.i <- alpha.out.co$alpha
    
    if (i > burnin){
      
      j <- i - burnin
      
      # samples.beta.ca[j,] <- beta.ca
      # samples.beta.co[j,] <- beta.co
      # samples.alpha.ca[j,] <- alpha.ca.i
      # samples.alpha.co[j,] <- alpha.co.i
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept[1] <- accept[1] + w.out.i$accept
      # accept[2] <- accept[2] + theta.out$accept
      # accept[3] <- accept[3] + beta.out.ca$accept
      # accept[4] <- accept[4] + beta.out.co$accept
      # accept[5] <- accept[5] + alpha.out.ca$accept
      # accept[6] <- accept[6] + alpha.out.co$accept
      
    }
    
    if (self_tune_w){
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    # if (self_tune_aca){
    #   a_ca_tuning <- update_tuning(a_ca_tuning, alpha.out.ca$a, i, alpha.out.ca$accept)
    #   deltas_aca <- c(deltas_aca, a_ca_tuning$delta_curr)
    # }
    # if (self_tune_aco){
    #   a_co_tuning <- update_tuning(a_co_tuning, alpha.out.co$a, i, alpha.out.co$accept)
    #   deltas_aco <- c(deltas_aco, a_co_tuning$delta_curr)
    # }
    # if (self_tune_ca){
    #   ca_tuning <- update_tuning(ca_tuning, beta.out.ca$a, i, beta.out.ca$accept)
    #   deltas_ca <- c(deltas_ca, ca_tuning$delta_curr)
    # }
    # if (self_tune_co){
    #   co_tuning <- update_tuning(co_tuning, beta.out.co$a, i, beta.out.co$accept)
    #   deltas_co <- c(deltas_co, co_tuning$delta_curr)
    # }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  
  output <- list()
  output$accept <- accept
  output$samples.beta.ca <- samples.beta.ca
  output$samples.beta.co <- samples.beta.co
  output$samples.alpha.ca <- samples.alpha.ca
  output$samples.alpha.co <- samples.alpha.co
  output$samples.theta <- samples.theta
  output$samples.phi <- samples.phi
  output$samples.w <- samples.w
  output$deltas_w <- deltas_w
  output$deltas_aca <- deltas_aca
  output$deltas_aco <- deltas_aco
  output$deltas_ca <- deltas_ca
  output$deltas_co <- deltas_co
  output$L_w <- L_w
  output$L_ca <- L_ca
  output$L_co <- L_co
  output$L_a_ca <- L_a_ca 
  output$L_a_co <- L_a_co
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}
