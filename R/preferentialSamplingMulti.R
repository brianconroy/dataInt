

prefSampleMulti_1 <- function(data, n.sample, burnin, 
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=NULL, m_aco=NULL, m_ca=NULL, m_co=NULL, m_w=NULL, 
                           target_aca=NULL, target_aco=NULL, target_ca=NULL, target_co=NULL, target_w=NULL, 
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
  deltas_ca <- array(NA, c(n.sample, N.d))
  deltas_co <- array(NA, c(n.sample, N.d))
  deltas_aca <- array(NA, c(n.sample, N.d))
  deltas_aco <- array(NA, c(n.sample, N.d))
  
  accept <- list(
    w=0,
    theta=0,
    beta.ca=rep(0, N.d),
    beta.co=rep(0, N.d),
    alpha.ca=rep(0, N.d),
    alpha.co=rep(0, N.d)
  )
  
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
    
    # sample from beta.cases, beta.ctrls, and alphas
    beta.ca.accepts <- c()
    beta.ca.as <- c()
    beta.co.accepts <- c()
    beta.co.as <- c()
    alpha.ca.accepts <- c()
    alpha.ca.as <- c()
    alpha.co.accepts <- c()
    alpha.co.as <- c()
    for (k in 1:N.d){
      w.i.sub <- w.i[locs[[k]]$ids]
      x.k <- case.data[[k]]$x.standardised
      beta.out.ca.k <- betaHmcUpdate(case.data[[k]]$y, w.i.sub, x.k, beta.ca[[k]], alpha.ca.i[[k]], ca_tuning[[k]]$delta_curr, L_ca[k])
      beta.ca.accepts[k] <- beta.out.ca.k$accept
      beta.ca.as[k] <- beta.out.ca.k$a
      beta.ca[[k]] <- beta.out.ca.k$beta
      
      beta.out.co.k <- betaHmcUpdate(ctrl.data[[k]]$y, w.i.sub, x.k, beta.co[[k]], alpha.co.i[[k]], co_tuning[[k]]$delta_curr, L_co[k])
      beta.co.accepts[k] <- beta.out.co.k$accept
      beta.co.as[k] <- beta.out.co.k$a
      beta.co[[k]] <- beta.out.co.k$beta
      
      alpha.out.ca <- alphaHmcUpdate(case.data[[k]]$y, w.i.sub, x.k, beta.ca[[k]], alpha.ca.i[[k]],
                                     a_ca_tuning[[k]]$delta_curr, prior_alpha_ca_mean[[k]], prior_alpha_ca_var[[k]], L_a_ca[k])
      alpha.ca.accepts[k] <- alpha.out.ca$accept
      alpha.ca.as[k] <- alpha.out.ca$a
      alpha.ca.i[[k]] <- alpha.out.ca$alpha
      
      alpha.out.co <- alphaHmcUpdate(ctrl.data[[k]]$y, w.i.sub, x.k, beta.co[[k]], alpha.co.i[[k]],
                                     a_co_tuning[[k]]$delta_curr, prior_alpha_co_mean[[k]], prior_alpha_co_var[[k]], L_a_co[k])
      alpha.co.accepts[k] <- alpha.out.co$accept
      alpha.co.as[k] <- alpha.out.co$a
      alpha.co.i[[k]] <- alpha.out.co$alpha
    }
    
    if (i > burnin){
      
      j <- i - burnin
      
      for (k in 1:N.d){
        samples.beta.ca[k,j,] <- beta.ca[[k]]
        samples.beta.co[k,j,] <- beta.co[[k]]
        samples.alpha.ca[k,j,] <- alpha.ca.i[[k]]
        samples.alpha.co[k,j,] <- alpha.co.i[[k]]
      }
      samples.theta[j,] <- theta.i
      samples.phi[j,] <- phi.i
      samples.w[j,] <- t(w.i)
      
      accept$w <- accept$w + w.out.i$accept
      accept$theta <- accept$theta + theta.out$accept
      for (k in 1:N.d){
        accept$beta.ca[k] <-  accept$beta.ca[k] + beta.ca.accepts[k]
        accept$beta.co[k] <-  accept$beta.co[k] + beta.co.accepts[k]
        accept$alpha.ca[k] <- accept$alpha.ca[k] + alpha.ca.accepts[k]
        accept$alpha.co[k] <- accept$alpha.co[k] + alpha.co.accepts[k]
      }
    }
    
    if (self_tune_w){
      w_tuning <- update_tuning(w_tuning, w.out.i$a, i, w.out.i$accept)
      deltas_w <- c(deltas_w, w_tuning$delta_curr)
    }
    if (self_tune_aca){
      for (k in 1:N.d){
        a_ca_tuning[[k]] <- update_tuning(a_ca_tuning[[k]], alpha.ca.as[k], i, alpha.ca.accepts[k])
        deltas_aca[i,k] <- a_ca_tuning[[k]]$delta_curr
      }
    }
    if (self_tune_aco){
      for (k in 1:N.d){
        a_co_tuning[[k]] <- update_tuning(a_co_tuning[[k]], alpha.co.as[k], i, alpha.co.accepts[k])
        deltas_aco[i,k] <- a_co_tuning[[k]]$delta_curr
      }
    }
    if (self_tune_ca){
      for (k in 1:N.d){
        ca_tuning[[k]] <- update_tuning(ca_tuning[[k]], beta.ca.as[k], i, beta.ca.accepts[k])
        deltas_ca[i,k] <- ca_tuning[[k]]$delta_curr
      }
    }
    if (self_tune_co){
      for (k in 1:N.d){
        co_tuning[[k]] <- update_tuning(co_tuning[[k]], beta.co.as[k], i, beta.co.accepts[k])
        deltas_co[i,k] <- co_tuning[[k]]$delta_curr
      }
    }
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/n.sample)
    }
    
  }
  
  accept <- accept/n.keep
  for (h in 1:length(accept)){
    accept[[h]] <- accept[[h]]/n.keep
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
