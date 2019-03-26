

#' prefSampleMVGP
#' 
#' Preferential sampling model for multiple species
#' with a multivariate Gaussian process.
#'
#' @param data 
#' @param d 
#' @param n.sample 
#' @param burnin 
#' @param L_w 
#' @param L_ca 
#' @param L_co 
#' @param L_a_ca 
#' @param L_a_co 
#' @param proposal.sd.theta 
#' @param m_aca 
#' @param m_aco 
#' @param m_ca 
#' @param m_co 
#' @param m_w 
#' @param target_aca 
#' @param target_aco 
#' @param target_ca 
#' @param target_co 
#' @param target_w 
#' @param self_tune_w 
#' @param self_tune_aca 
#' @param self_tune_aco 
#' @param self_tune_ca 
#' @param self_tune_co 
#' @param delta_w 
#' @param delta_aca 
#' @param delta_aco 
#' @param delta_ca 
#' @param delta_co 
#' @param beta_ca_initial 
#' @param beta_co_initial 
#' @param alpha_ca_initial 
#' @param alpha_co_initial 
#' @param theta_initial 
#' @param t_initial 
#' @param w_initial 
#' @param prior_phi 
#' @param prior_theta 
#' @param prior_alpha_ca_mean 
#' @param prior_alpha_co_mean 
#' @param prior_alpha_ca_var 
#' @param prior_alpha_co_var 
#' @param prior_t 
#'
#' @return
#' @export
#'
#' @examples
prefSampleMVGP <- function(data, d, n.sample, burnin,
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=NULL, m_aco=NULL, m_ca=NULL, m_co=NULL, m_w=NULL, 
                           target_aca=NULL, target_aco=NULL, target_ca=NULL, target_co=NULL, target_w=NULL, 
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                           delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                           beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                           theta_initial=NULL, t_initial=NULL, w_initial=NULL,
                           prior_phi, prior_theta, prior_alpha_ca_mean, prior_alpha_co_mean, prior_alpha_ca_var, prior_alpha_co_var,
                           prior_t){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
  locs <- data$locs
  N.w <- length(locs[[1]]$status)
  N.d <- length(locs)
  N.p <- ncol(case.data[[1]]$x.standardised)
  
  ## starting values
  if (is.null(w_initial)){
    w.i <- rnorm(N.w*N.d)
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
  if (is.null(t_initial)){
    t.i <- matrix(c(3, 0, 0, 3), nrow=2)
  } else {
    t.i <- t_initial
  }
  H.i <- Exponential(d, range=theta.i, phi=1)
  H.inv.i <- solve(H.i)
  
  Omega <- prior_t$scale
  r <- prior_t$df
  
  # storage
  n.keep <- n.sample - burnin
  samples.w <- array(NA, c(n.keep, N.w * N.d))
  samples.theta <- array(NA, c(n.keep, 1))
  samples.t <- array(NA, c(n.keep, N.d * N.d))
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
      a_ca_tuning[[i]] <- list(delta_curr=delta_aca[i])
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
      a_co_tuning[[i]] <- list(delta_curr=delta_aco[i])
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
      ca_tuning[[i]] <- list(delta_curr=delta_ca[i])
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
      co_tuning[[i]] <- list(delta_curr=delta_co[i])
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
    sigma.i <- kronecker(H.i, t.i)
    sigma.inv.i <- kronecker(solve(H.i), solve(t.i))
    w.out.i <- wHmcUpdateMVGP(case.data, ctrl.data, alpha.ca.i, beta.ca,
                              alpha.co.i, beta.co, w.i, sigma.i, sigma.inv.i, locs, w_tuning$delta_curr, L_w)
    w.i <- w.out.i$w

    ## sample from theta
    theta.out <- rangeMVGPupdate(H.i, t.i, w.i, d, theta.i, proposal.sd.theta, prior_theta)
    theta.i <- theta.out$theta
    H.i <- Exponential(d, range=theta.i, phi=1)
    H.inv.i <- solve(H.i)
    
    ## sample from T
    r_ <- r + N.w
    Omega_ <- Omega
    w1.i <- w.i[seq(1, length(w.i), by=N.d)]
    w2.i <- w.i[seq(2, length(w.i), by=N.d)]
    for (a in 1:N.w){
      for (b in 1:N.w){
        Omega_ <- Omega_ + H.inv.i[a, b] * matrix(c(w1.i[b], w2.i[b])) %*% t(matrix(c(w1.i[a], w2.i[a])))
      }
    }
    t.i <- riwish(r_, Omega_)
    
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
      w.i.k <- w.i[seq(k, length(w.i), by=N.d)]
      w.i.sub <- w.i.k[locs[[k]]$ids]
      x.k <- case.data[[k]]$x.standardised
      beta.out.ca.k <- betaHmcUpdate(case.data[[k]]$y, w.i.sub, x.k, beta.ca[[k]], alpha.ca.i[[k]], ca_tuning[[k]]$delta_curr, L_ca[k], offset=0)
      beta.ca.accepts[k] <- beta.out.ca.k$accept
      beta.ca.as[k] <- beta.out.ca.k$a
      beta.ca[[k]] <- beta.out.ca.k$beta
      
      beta.out.co.k <- betaHmcUpdate(ctrl.data[[k]]$y, w.i.sub, x.k, beta.co[[k]], alpha.co.i[[k]], co_tuning[[k]]$delta_curr, L_co[k], offset=0)
      beta.co.accepts[k] <- beta.out.co.k$accept
      beta.co.as[k] <- beta.out.co.k$a
      beta.co[[k]] <- beta.out.co.k$beta
      
      alpha.out.ca <- alphaHmcUpdate(case.data[[k]]$y, w.i.sub, x.k, beta.ca[[k]], alpha.ca.i[[k]],
                                     a_ca_tuning[[k]]$delta_curr, prior_alpha_ca_mean[[k]], prior_alpha_ca_var[[k]], L_a_ca[k], offset=0)
      alpha.ca.accepts[k] <- alpha.out.ca$accept
      alpha.ca.as[k] <- alpha.out.ca$a
      alpha.ca.i[[k]] <- alpha.out.ca$alpha
      
      alpha.out.co <- alphaHmcUpdate(ctrl.data[[k]]$y, w.i.sub, x.k, beta.co[[k]], alpha.co.i[[k]],
                                     a_co_tuning[[k]]$delta_curr, prior_alpha_co_mean[[k]], prior_alpha_co_var[[k]], L_a_co[k], offset=0)
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
      samples.t[j,] <- t.i
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
  output$samples.t <- samples.t
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
  output$m_w <- m_w
  output$m_aca <- m_aca
  output$m_aco <- m_aco
  output$m_ca <- m_ca
  output$m_co <- m_co
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_t <- prior_t
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}


prefSampleMVGP2 <- function(data, d, n.sample, burnin,
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=NULL, m_aco=NULL, m_ca=NULL, m_co=NULL, m_w=NULL, 
                           target_aca=NULL, target_aco=NULL, target_ca=NULL, target_co=NULL, target_w=NULL, 
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                           delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                           beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                           theta_initial=NULL, t_initial=NULL, w1_initial=NULL, w2_initial=NULL,
                           prior_phi, prior_theta, prior_alpha_ca_mean, prior_alpha_co_mean, prior_alpha_ca_var, prior_alpha_co_var,
                           prior_t){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
  locs <- data$locs
  N.w <- length(locs[[1]]$status)
  N.d <- length(locs)
  N.p <- ncol(case.data[[1]]$x.standardised)
  Y1.l <- locs[[1]]$status
  X1.c <- case.data[[1]]$x.standardised
  Y1.ca <- case.data[[1]]$y
  Y1.co <- ctrl.data[[1]]$y
  Y2.l <- locs[[2]]$status
  X2.c <- case.data[[2]]$x.standardised
  Y2.ca <- case.data[[2]]$y
  Y2.co <- ctrl.data[[2]]$y
  
  ## starting values
  if (is.null(g_initial)){
    g.i <- rnorm(N.w*N.d)
  } else {
    g.i <- g_initial
  }
  if (is.null(tau_initial)){
    tau.i <- runif(1, 5, 10)
  } else {
    tau.i <- tau_initial
  }
  if (is.null(w1_initial)){
    w1.i <- rnorm(N.w)
  } else {
    w1.i <- w1_initial
  }
  if (is.null(w2_initial)){
    w2.i <- rnorm(N.w)
  } else {
    w2.i <- w2_initial
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
  if (is.null(theta1_initial)){
    theta1.i <- runif(1, 5, 7)
  } else {
    theta1.i <- theta1_initial
  }
  if (is.null(theta2_initial)){
    theta2.i <- runif(1, 5, 7)
  } else {
    theta2.i <- theta2_initial
  }
  if (is.null(phi1_initial)){
    phi1.i <- runif(1, 5, 7)
  } else {
    phi1.i <- phi1_initial
  }
  if (is.null(phi2_initial)){
    phi2.i <- runif(1, 5, 7)
  } else {
    phi2.i <- phi2_initial
  }
  if (is.null(t_initial)){
    t.i <- 3 * diag(4)
  } else {
    t.i <- t_initial
  }
  H.i <- Exponential(d, range=tau.i, phi=1)
  H.inv.i <- solve(H.i)
  
  Omega <- prior_t$scale
  r <- prior_t$df
  
  # storage
  n.keep <- n.sample - burnin
  samples.g <- array(NA, c(n.keep, N.w * N.d))
  samples.t <- array(NA, c(n.keep, N.d * N.d))
  samples.tau <- array(NA, c(n.keep, 1))
  samples.w1 <- array(NA, c(n.keep, N.w))
  samples.w2 <- array(NA, c(n.keep, N.w))
  samples.theta1 <- array(NA, c(n.keep, 1))
  samples.theta2 <- array(NA, c(n.keep, 1))
  samples.phi1 <- array(NA, c(n.keep, 1))
  samples.phi2 <- array(NA, c(n.keep, 1))
  samples.alpha.ca <- array(NA, c(N.d, n.keep, 1))
  samples.alpha.co <- array(NA, c(N.d, n.keep, 1))
  samples.beta.ca <- array(NA, c(N.d, n.keep, N.p))
  samples.beta.co <- array(NA, c(N.d, n.keep, N.p))
  
  
  if (self_tune_g){
    g_tuning <- initialize_tuning(m=m_g, target=target_g)
  } else{
    g_tuning <- list(delta_curr=delta_g)
  }
  if (self_tune_w1){
    w1_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w1_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_w2){
    w2_tuning <- initialize_tuning(m=m_w, target=target_w)
  } else {
    w2_tuning <- list(delta_curr=delta_w)
  }
  if (self_tune_aca){
    a_ca_tuning <- list()
    for (i in 1:N.d){
      a_ca_tuning[[i]] <- initialize_tuning(m=m_aca, target=target_aca)
    }
  } else {
    a_ca_tuning <- list()
    for (i in 1:N.d){
      a_ca_tuning[[i]] <- list(delta_curr=delta_aca[i])
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
      a_co_tuning[[i]] <- list(delta_curr=delta_aco[i])
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
      ca_tuning[[i]] <- list(delta_curr=delta_ca[i])
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
      co_tuning[[i]] <- list(delta_curr=delta_co[i])
    }
  }
  
  deltas_g <- c()
  deltas_w1<- c()
  deltas_w2<- c()
  deltas_ca <- array(NA, c(n.sample, N.d))
  deltas_co <- array(NA, c(n.sample, N.d))
  deltas_aca <- array(NA, c(n.sample, N.d))
  deltas_aco <- array(NA, c(n.sample, N.d))
  
  accept <- list(
    g=0,
    tau=0,
    w1=0,
    w2=0,
    theta1=0,
    theta2=0,
    beta.ca=rep(0, N.d),
    beta.co=rep(0, N.d),
    alpha.ca=rep(0, N.d),
    alpha.co=rep(0, N.d)
  )
  
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*n.sample)
  
  for (i in 1:n.sample){
    
    ## sample from gamma
    sigma.g.i <- kronecker(H.i, t.i)
    sigma.g.inv.i <- kronecker(solve(H.i), solve(t.i))
    w.out.i <- wHmcUpdateMVGP2(case.data, ctrl.data, 
                               alpha.ca.i, beta.ca,
                               alpha.co.i, beta.co,
                               g.i,
                               w1.i, w2.i, 
                               sigma.g.i, sigma.g.inv.i, 
                               locs, g_tuning$delta_curr, L_g)
    w.i <- w.out.i$w
    
    tau.i <- Tau
    ## sample from tau 
    # theta.out <- rangeMVGPupdate(H.i, t.i, w.i, d, theta.i, proposal.sd.theta, prior_theta)
    # theta.i <- theta.out$theta
    # H.i <- Exponential(d, range=theta.i, phi=1)
    # H.inv.i <- solve(H.i)
    
    t.i <- Tmat
    ## sample from T
    # r_ <- r + N.w
    # Omega_ <- Omega
    # w1.i <- w.i[seq(1, length(w.i), by=N.d)]
    # w2.i <- w.i[seq(2, length(w.i), by=N.d)]
    # for (a in 1:N.w){
    #   for (b in 1:N.w){
    #     Omega_ <- Omega_ + H.inv.i[a, b] * matrix(c(w1.i[b], w2.i[b])) %*% t(matrix(c(w1.i[a], w2.i[a])))
    #   }
    # }
    # t.i <- riwish(r_, Omega_)
    
    ## sample from w1
    sigma1.i <- Exponential(d, range=theta1.i, phi=phi1.i)
    sigma1.inv.i <- solve(sigma1.i)
    w1.out.i <- wHmcUpdateCC(Y1.l, X1.c, Y1.ca, alpha.ca.i[[1]], beta.ca[[1]], Y1.co,
                            alpha.co.i[[1]], beta.co[[1]], w1.i, sigma1.i, sigma1.inv.i, locs[[1]], 
                            w1_tuning$delta_curr, L_w, offset=0)
    w1.i <- w1.out.i$w
    
    ## sample from w2
    sigma2.i <- Exponential(d, range=theta2.i, phi=phi2.i)
    sigma2.inv.i <- solve(sigma2.i)
    w2.out.i <- wHmcUpdateCC(Y2.l, X2.c, Y2.ca, alpha.ca.i[[2]], beta.ca[[2]], Y2.co,
                             alpha.co.i[[2]], beta.co[[2]], w2.i, sigma2.i, sigma2.inv.i, locs[[2]], 
                             w2_tuning$delta_curr, L_w, offset=0)
    w2.i <- w2.out.i$w
    
    ## sample from theta1
    theta1.out <- rangeMhUpdate(theta1.i, as.numeric(w1.i), d, phi1.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
    theta1.i <- theta1.out$theta
    
    ## sample from theta2
    theta2.out <- rangeMhUpdate(theta2.i, as.numeric(w2.i), d, phi2.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
    theta2.i <- theta2.out$theta
    
    ## sample from phi1
    R1.i <- sigma1.i/phi1.i
    phi1.i <- 1/rgamma(1, N.w/2 + prior.phi[1], t(w1.i) %*% solve(R1.i) %*% w1.i/2 + prior.phi[2])
    
    ## sample from phi2
    R2.i <- sigma2.i/phi2.i
    phi2.i <- 1/rgamma(1, N.w/2 + prior.phi[1], t(w2.i) %*% solve(R2.i) %*% w2.i/2 + prior.phi[2])
    
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
      #w.i.k <- w.i[seq(k, length(w.i), by=N.d)]
      #w.i.sub <- w.i.k[locs[[k]]$ids]
      if (k==1){
        w.i.k <- w1.i
      } else{
        w.i.k <- w2.i
      }
      w.i.sub <- w.i.k[locs[[k]]$ids]
      
      x.k <- case.data[[k]]$x.standardised
      beta.out.ca.k <- betaHmcUpdate(case.data[[k]]$y, w.i.sub, x.k, beta.ca[[k]], alpha.ca.i[[k]], ca_tuning[[k]]$delta_curr, L_ca[k], offset=0)
      beta.ca.accepts[k] <- beta.out.ca.k$accept
      beta.ca.as[k] <- beta.out.ca.k$a
      beta.ca[[k]] <- beta.out.ca.k$beta
      
      beta.out.co.k <- betaHmcUpdate(ctrl.data[[k]]$y, w.i.sub, x.k, beta.co[[k]], alpha.co.i[[k]], co_tuning[[k]]$delta_curr, L_co[k], offset=0)
      beta.co.accepts[k] <- beta.out.co.k$accept
      beta.co.as[k] <- beta.out.co.k$a
      beta.co[[k]] <- beta.out.co.k$beta
      
      alpha.out.ca <- alphaHmcUpdate(case.data[[k]]$y, w.i.sub, x.k, beta.ca[[k]], alpha.ca.i[[k]],
                                     a_ca_tuning[[k]]$delta_curr, prior_alpha_ca_mean[[k]], prior_alpha_ca_var[[k]], L_a_ca[k], offset=0)
      alpha.ca.accepts[k] <- alpha.out.ca$accept
      alpha.ca.as[k] <- alpha.out.ca$a
      alpha.ca.i[[k]] <- alpha.out.ca$alpha
      
      alpha.out.co <- alphaHmcUpdate(ctrl.data[[k]]$y, w.i.sub, x.k, beta.co[[k]], alpha.co.i[[k]],
                                     a_co_tuning[[k]]$delta_curr, prior_alpha_co_mean[[k]], prior_alpha_co_var[[k]], L_a_co[k], offset=0)
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
      samples.theta1[j,] <- theta1.i
      # samples.t[j,] <- t.i
      samples.w1[j,] <- t(w1.i)
      samples.w2[j,] <- t(w2.i)
      
      accept$w1 <- accept$w1 + w1.out.i$accept
      accept$w2 <- accept$w2 + w2.out.i$accept
      accept$theta1 <- accept$theta1 + theta1.out$accept
      accept$theta2 <- accept$theta2 + theta2.out$accept
      for (k in 1:N.d){
        accept$beta.ca[k] <-  accept$beta.ca[k] + beta.ca.accepts[k]
        accept$beta.co[k] <-  accept$beta.co[k] + beta.co.accepts[k]
        accept$alpha.ca[k] <- accept$alpha.ca[k] + alpha.ca.accepts[k]
        accept$alpha.co[k] <- accept$alpha.co[k] + alpha.co.accepts[k]
      }
    }
    
    if (self_tune_w1){
      w1_tuning <- update_tuning(w1_tuning, w1.out.i$a, i, w1.out.i$accept)
      deltas_w1 <- c(deltas_w1, w1_tuning$delta_curr)
    }
    if (self_tune_w2){
      w2_tuning <- update_tuning(w2_tuning, w2.out.i$a, i, w2.out.i$accept)
      deltas_w2 <- c(deltas_w2, w2_tuning$delta_curr)
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
  
  for (h in 1:length(accept)){
    accept[[h]] <- accept[[h]]/n.keep
  }
  
  output <- list()
  output$accept <- accept
  output$samples.beta.ca <- samples.beta.ca
  output$samples.beta.co <- samples.beta.co
  output$samples.alpha.ca <- samples.alpha.ca
  output$samples.alpha.co <- samples.alpha.co
  output$samples.theta1 <- samples.theta1
  output$samples.theta2 <- samples.theta2
  # output$samples.t <- samples.t
  output$samples.w1 <- samples.w1
  output$samples.w2 <- samples.w2
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
  output$m_w <- m_w
  output$m_aca <- m_aca
  output$m_aco <- m_aco
  output$m_ca <- m_ca
  output$m_co <- m_co
  output$proposal.sd.theta <- proposal.sd.theta
  # output$prior_t <- prior_t
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}


continueMCMC_mvgp <- function(data, d, output, n.sample){
  
  # get initial values
  n.sample.old <- nrow(output$samples.beta.ca[1,,])
  n_s <- dim(output$samples.alpha.ca)[1]
  beta_ca_initial <- list()
  beta_co_initial <- list()
  alpha_ca_initial <- list()
  alpha_co_initial <- list()
  for (j in 1:dim(output$samples.alpha.ca)[1]){
    beta_ca_initial[[j]] <- output$samples.beta.ca[j,n.sample.old,]
    beta_co_initial[[j]] <- output$samples.beta.co[j,n.sample.old,]
    alpha_ca_initial[[j]] <- output$samples.alpha.ca[j,n.sample.old,]
    alpha_co_initial[[j]] <- output$samples.alpha.co[j,n.sample.old,]
  }
  w_initial <- output$samples.w[n.sample.old,]
  theta_initial <- output$samples.theta[n.sample.old]
  t_initial <- matrix(output$samples.t[n.sample.old,], ncol=2)
  
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
  prior_t <- output$prior_t
  more_output <- prefSampleMVGP(data, d, n.sample, burnin=0, 
                                   L_w, L_ca, L_co, L_a_ca, L_a_co,
                                   proposal.sd.theta=proposal.sd.theta,
                                   m_aca=m_aca, m_aco=m_aca, m_ca=m_aca, m_co=m_aca, m_w=m_aca, 
                                   target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, 
                                   self_tune_w=FALSE, self_tune_aca=FALSE, self_tune_aco=FALSE, self_tune_ca=FALSE, self_tune_co=FALSE,
                                   delta_w=delta_w, delta_aca=delta_aca, delta_aco=delta_aco, delta_ca=delta_ca, delta_co=delta_co, 
                                   beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                                   theta_initial=theta_initial, t_initial=t_initial, w_initial=w_initial,
                                   prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_mean=prior_alpha_ca_mean,
                                   prior_alpha_co_mean=prior_alpha_co_mean,
                                   prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var,
                                   prior_t=prior_t)
  
  # combine outputs
  new_output <- output
  new_output$n.sample <- new_output$n.sample + n.sample
  samples.alpha.ca <- array(NA, c(n_s, n.sample.old + n.sample, 1))
  samples.alpha.co <- array(NA, c(n_s, n.sample.old + n.sample, 1))
  samples.beta.ca <- array(NA, c(n_s, n.sample.old + n.sample, 3))
  samples.beta.co <- array(NA, c(n_s, n.sample.old + n.sample, 3))
  for (j in 1:dim(output$samples.alpha.ca)[1]){
    samples.alpha.ca[j,,] <- matrix(c(new_output$samples.alpha.ca[j,,], more_output$samples.alpha.ca[j,,]))
    samples.alpha.co[j,,] <- matrix(c(new_output$samples.alpha.co[j,,], more_output$samples.alpha.co[j,,]))
    samples.beta.ca[j,,] <- rbind(new_output$samples.beta.ca[j,,], more_output$samples.beta.ca[j,,])
    samples.beta.co[j,,] <- rbind(new_output$samples.beta.co[j,,], more_output$samples.beta.co[j,,])
  }
  new_output$samples.alpha.ca <- samples.alpha.ca
  new_output$samples.alpha.co <- samples.alpha.co
  new_output$samples.beta.ca <- samples.beta.ca
  new_output$samples.beta.co <- samples.beta.co
  
  new_output$samples.t <- rbind(new_output$samples.t, more_output$samples.t)
  new_output$samples.theta <- matrix(c(new_output$samples.theta, more_output$samples.theta))
  new_output$samples.w <- rbind(new_output$samples.w, more_output$samples.w)
  # new_output$deltas_aca <- new_output$deltas_aca
  # new_output$deltas_aco <- new_output$deltas_aco
  # new_output$deltas_co <- new_output$deltas_co
  # new_output$deltas_ca <- new_output$deltas_ca
  # new_output$deltas_w <- new_output$deltas_w
  
  new_output$accept$w <- (new_output$n.sample * new_output$accept$w + n.sample * more_output$accept$w)/(new_output$n.sample + n.sample)
  new_output$accept$theta <- (new_output$n.sample * new_output$accept$theta + n.sample * more_output$accept$theta)/(new_output$n.sample + n.sample)
  for (k in 1:dim(output$samples.alpha.ca)[1]){
    new_output$accept$beta.ca[k] <- (new_output$n.sample * new_output$accept$beta.ca[k] + n.sample * more_output$accept$beta.ca[k])/(new_output$n.sample + n.sample)
    new_output$accept$beta.co[k] <- (new_output$n.sample * new_output$accept$beta.co[k] + n.sample * more_output$accept$beta.co[k])/(new_output$n.sample + n.sample)
    new_output$accept$alpha.ca[k] <- (new_output$n.sample * new_output$accept$alpha.ca[k] + n.sample * more_output$accept$alpha.ca[k])/(new_output$n.sample + n.sample)
    new_output$accept$alpha.co[k] <- (new_output$n.sample * new_output$accept$alpha.co[k] + n.sample * more_output$accept$alpha.co[k])/(new_output$n.sample + n.sample)
  }
  
  return(new_output)
  
}


burnin_mvgp <- function(output, n.burn){
  
  n.curr <- output$n.sample - output$burnin
  i.start <- n.burn + 1
  output$burnin <- output$burnin + n.burn
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
  output$samples.w <- output$samples.w[i.start:n.curr,]
  output$samples.t <- output$samples.t[i.start:n.curr,]
  output$samples.theta <- output$samples.theta[i.start:n.curr]
  return(output)
  
}


prefSampleMulti_1 <- function(data, n.sample, burnin, 
                           L_w, L_ca, L_co, L_a_ca, L_a_co,
                           proposal.sd.theta=0.3,
                           m_aca=NULL, m_aco=NULL, m_ca=NULL, m_co=NULL, m_w=NULL, 
                           target_aca=NULL, target_aco=NULL, target_ca=NULL, target_co=NULL, target_w=NULL, 
                           self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                           delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                           beta_ca_initial=NULL, beta_co_initial=NULL, alpha_ca_initial=NULL, alpha_co_initial=NULL,
                           theta_initial=NULL, phi_initial=NULL, w_initial=NULL,
                           prior_phi, prior_theta, prior_alpha_ca_mean, prior_alpha_co_mean, prior_alpha_ca_var, prior_alpha_co_var){
  
  
  ## setup
  case.data <- data$case.data
  ctrl.data <- data$ctrl.data
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
      a_ca_tuning[[i]] <- list(delta_curr=delta_aca[i])
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
      a_co_tuning[[i]] <- list(delta_curr=delta_aco[i])
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
      ca_tuning[[i]] <- list(delta_curr=delta_ca[i])
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
      co_tuning[[i]] <- list(delta_curr=delta_co[i])
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
  output$m_w <- m_w
  output$m_aca <- m_aca
  output$m_aco <- m_aco
  output$m_ca <- m_ca
  output$m_co <- m_co
  output$proposal.sd.theta <- proposal.sd.theta
  output$prior_phi <- prior_phi
  output$prior_theta <- prior_theta
  output$prior_alpha_ca_var <- prior_alpha_ca_var
  output$prior_alpha_co_var <- prior_alpha_co_var
  output$n.sample <- n.sample
  output$burnin <- burnin
  
  return(output)
  
}


continueMCMC_multi <- function(data, output, n.sample){
  
  # get initial values
  n.sample.old <- nrow(output$samples.beta.ca[1,,])
  n_s <- dim(output$samples.alpha.ca)[1]
  beta_ca_initial <- list()
  beta_co_initial <- list()
  alpha_ca_initial <- list()
  alpha_co_initial <- list()
  for (j in 1:dim(output$samples.alpha.ca)[1]){
    beta_ca_initial[[j]] <- output$samples.beta.ca[j,n.sample.old,]
    beta_co_initial[[j]] <- output$samples.beta.co[j,n.sample.old,]
    alpha_ca_initial[[j]] <- output$samples.alpha.ca[j,n.sample.old,]
    alpha_co_initial[[j]] <- output$samples.alpha.co[j,n.sample.old,]
  }
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
  more_output <- prefSampleMulti_1(data, n.sample, burnin=0, 
                                   L_w, L_ca, L_co, L_a_ca, L_a_co,
                                   proposal.sd.theta=0.3,
                                   m_aca=m_aca, m_aco=m_aca, m_ca=m_aca, m_co=m_aca, m_w=m_aca, 
                                   target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, 
                                   self_tune_w=FALSE, self_tune_aca=FALSE, self_tune_aco=FALSE, self_tune_ca=FALSE, self_tune_co=FALSE,
                                   delta_w=delta_w, delta_aca=delta_aca, delta_aco=delta_aco, delta_ca=delta_ca, delta_co=delta_co, 
                                   beta_ca_initial=beta_ca_initial, beta_co_initial=beta_co_initial, alpha_ca_initial=alpha_ca_initial, alpha_co_initial=alpha_co_initial,
                                   theta_initial=theta_initial, phi_initial=phi_initial, w_initial=w_initial,
                                   prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_mean=prior_alpha_ca_mean,
                                   prior_alpha_co_mean=prior_alpha_co_mean,
                                   prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)
  
  # combine outputs
  new_output <- output
  new_output$n.sample <- new_output$n.sample + n.sample
  samples.alpha.ca <- array(NA, c(n_s, n.sample.old + n.sample, 1))
  samples.alpha.co <- array(NA, c(n_s, n.sample.old + n.sample, 1))
  samples.beta.ca <- array(NA, c(n_s, n.sample.old + n.sample, 3))
  samples.beta.co <- array(NA, c(n_s, n.sample.old + n.sample, 3))
  for (j in 1:dim(output$samples.alpha.ca)[1]){
    samples.alpha.ca[j,,] <- matrix(c(new_output$samples.alpha.ca[j,,], more_output$samples.alpha.ca[j,,]))
    samples.alpha.co[j,,] <- matrix(c(new_output$samples.alpha.co[j,,], more_output$samples.alpha.co[j,,]))
    samples.beta.ca[j,,] <- rbind(new_output$samples.beta.ca[j,,], more_output$samples.beta.ca[j,,])
    samples.beta.co[j,,] <- rbind(new_output$samples.beta.co[j,,], more_output$samples.beta.co[j,,])
  }
  new_output$samples.alpha.ca <- samples.alpha.ca
  new_output$samples.alpha.co <- samples.alpha.co
  new_output$samples.beta.ca <- samples.beta.ca
  new_output$samples.beta.co <- samples.beta.co
  new_output$samples.phi <- matrix(c(new_output$samples.phi, more_output$samples.phi))
  new_output$samples.theta <- matrix(c(new_output$samples.theta, more_output$samples.theta))
  new_output$samples.w <- rbind(new_output$samples.w, more_output$samples.w)
  new_output$deltas_aca <- new_output$deltas_aca
  new_output$deltas_aco <- new_output$deltas_aco
  new_output$deltas_co <- new_output$deltas_co
  new_output$deltas_ca <- new_output$deltas_ca
  new_output$deltas_w <- new_output$deltas_w
  
  new_output$accept$w <- (new_output$n.sample * new_output$accept$w + n.sample * more_output$accept$w)/(new_output$n.sample + n.sample)
  new_output$accept$theta <- (new_output$n.sample * new_output$accept$theta + n.sample * more_output$accept$theta)/(new_output$n.sample + n.sample)
  for (k in 1:dim(output$samples.alpha.ca)[1]){
    new_output$accept$beta.ca[k] <- (new_output$n.sample * new_output$accept$beta.ca[k] + n.sample * more_output$accept$beta.ca[k])/(new_output$n.sample + n.sample)
    new_output$accept$beta.co[k] <- (new_output$n.sample * new_output$accept$beta.co[k] + n.sample * more_output$accept$beta.co[k])/(new_output$n.sample + n.sample)
    new_output$accept$alpha.ca[k] <- (new_output$n.sample * new_output$accept$alpha.ca[k] + n.sample * more_output$accept$alpha.ca[k])/(new_output$n.sample + n.sample)
    new_output$accept$alpha.co[k] <- (new_output$n.sample * new_output$accept$alpha.co[k] + n.sample * more_output$accept$alpha.co[k])/(new_output$n.sample + n.sample)
  }
  
  return(new_output)
  
}
