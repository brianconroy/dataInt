
############################
# Hamiltonian MCMC to sample
# from a bivariate normal

# No-U-Turn sampler. 
# Gooo Nutssss!

# not tuning delta here.
############################


mu <- c(0, 0)
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
sigma.inv <- solve(sigma)


L.norm <- function(x){
  return(- 0.5 * t(x) %*% sigma.inv %*% x)
}


dL.norm <- function(x){
  return(- t(x) %*% sigma.inv)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    n.prime <- as.numeric( u <= exp(L(theta.prime) - 0.5 * t(r.prime) %*% r.prime ) )
    s.prime <- as.numeric( L(theta.prime) - 0.5 * t(r.prime) %*% r.prime > log(u) - 1000 )
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
    
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime
      )
    )
  }
}


M <- 2000
epsilon <- 0.05
theta.m <- matrix(c(0, 6))
samples.theta <- array(NA, c(M, 2))

nrgks <- c()
for (m in 1:M){
  
  r0 <- matrix(rnorm(2))
  u <- runif(1, 0, exp(L.norm(theta.m) - 0.5 * t(r0) %*% r0))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r0
  r.p <- r0
  j <- 0
  n <- 1
  s <- 1
  
  nrgk <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, L=L.norm, dL=dL.norm)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, L=L.norm, dL=dL.norm)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
      }
    }
    
    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    
    nrgk <- nrgk + 1
    
  }
  nrgks <- c(nrgks, nrgk)
  samples.theta[m,] <- t(theta.m)
  
}

colMeans(samples.theta)
cov(samples.theta[,1], samples.theta[,2])
cov(samples.theta[,1], samples.theta[,1])
plot(samples.theta[,1], type='l')
plot(samples.theta[,2], type='l')
plot(samples.theta); abline(h=0, col=2); abline(v=0, col=2)


#####################
# NUTS sampler
# with dual averaging

# bivariate normal
#####################


mu <- c(0, 0)
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
sigma.inv <- solve(sigma)


L.norm <- function(x){
  return(- 0.5 * t(x) %*% sigma.inv %*% x)
}


dL.norm <- function(x){
  return(- t(x) %*% sigma.inv)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


K <- function(p){
  return(t(p) %*% p/2)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u < exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


initialize_tuning <- function(m, target){
  return(list(
    M_adapt=m,
    epsilon_curr=0.05,
    epsilon_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, n.a){
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    epsilon_bar <- tuning$epsilon_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a/n.a)
    log_epsilon_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_epsilon_bar <- i^{-kappa} * log_epsilon_curr + (1 - i^{-kappa}) * log(epsilon_bar)
    epsilon_curr <- exp(log_epsilon_curr)
    epsilon_bar <- exp(log_epsilon_bar)
    
    tuning$epsilon_bar <- epsilon_bar
    tuning$HbarM <- HbarM
    tuning$epsilon_curr <- epsilon_curr
    
  } else{
    
    tuning$epsilon_curr <- tuning$epsilon_bar
    if (tuning$reject_streak > 1000) {
      tuning$epsilon_curr <- 0.90 * tuning$epsilon_curr
      tuning$epsilon_bar <- 0.90 * tuning$epsilon_bar
    }
    
  }
  
  return(tuning)
  
}


M <- 5000
M.adapt <- 1000
tuning <- initialize_tuning(M.adapt, target=0.65)
theta.m <- matrix(c(0, 6))
samples.theta <- array(NA, c(M, 2))

Ls <- c()
epsilons <- c()

accept <- 0
n.accept <- 0
for (m in 1:M){
  
  epsilon <- tuning$epsilon_curr
  r.0 <- matrix(rnorm(2))
  u <- runif(1, 0, exp(L.norm(theta.m) - K(r.0)))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r.0
  r.p <- r.0
  j <- 0
  n <- 1
  s <- 1
  
  L.curr <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L.norm, dL=dL.norm)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L.norm, dL=dL.norm)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
        accept <- accept + 1
      }
    }

    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    
    L.curr <- L.curr + 1
    
  }
  
  n.accept <- n.accept + L.curr
  Ls <- c(Ls, L.curr)
  epsilons <- c(epsilons, epsilon)
  samples.theta[m,] <- t(theta.m)
  tuning <- update_tuning(tuning, alpha, m, n.alpha)

}

colMeans(samples.theta)
cov(samples.theta[,1], samples.theta[,2])
sd(samples.theta[,1])
sd(samples.theta[,2])
plot(samples.theta[,1], samples.theta[,2])
plot(samples.theta[,1], type='l'); abline(h=0, col='2')
plot(samples.theta[,2], type='l'); abline(h=0, col='2')

hist(Ls)
plot(epsilons)

######################
# Hamiltonian MCMC for 
# Poisson regression
######################


Beta <- c(2, 2)
N <- 100
X <- cbind(rep(1, N), rnorm(N))
rates <- exp(X %*% Beta)
Y <- sapply(rates, function(x){rpois(n=1, x)})

prior_mu <- rep(0, length(Beta))
prior_var <- rep(1000, length(Beta))
prior_sigma <- diag(prior_var)
prior_sigma_inv <- solve(prior_sigma)


L.beta_pois <- function(b){
  
  # likelihood
  logd <- 0
  lin_preds <- X %*% b
  for (i in 1:N){
    logd <- logd + dpois(Y[i], lambda=exp(lin_preds[i]), log=TRUE)
  }
  
  # prior
  for (i in 1:length(b)){
    logd <- logd + dnorm(b[i], mean=prior_mu[i], sd=sqrt(prior_var[i]), log=T)
  }
  
  return(logd)
}


dL.beta_pois <- function(b){
  
  grad <- array(0, c(length(Beta), 1))
  lin_preds <- X %*% b
  for (i in 1:length(Y)){
    for (j in 1:length(b)){
      grad[j] <- grad[j] + X[i,j]*Y[i] - X[i,j]*exp(lin_preds[i,])
    }
  }
  
  # prior contribution
  grad <- grad - t(t(b) %*% prior_sigma_inv)
  
  return(t(grad))
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


initialize_tuning <- function(m, target, eps){
  return(list(
    M_adapt=m,
    epsilon_curr=eps,
    epsilon_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, n.a){
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    epsilon_bar <- tuning$epsilon_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a/n.a)
    log_epsilon_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_epsilon_bar <- i^{-kappa} * log_epsilon_curr + (1 - i^{-kappa}) * log(epsilon_bar)
    epsilon_curr <- exp(log_epsilon_curr)
    epsilon_bar <- exp(log_epsilon_bar)
    
    tuning$epsilon_bar <- epsilon_bar
    tuning$HbarM <- HbarM
    tuning$epsilon_curr <- epsilon_curr
    
  } else{
    
    tuning$epsilon_curr <- tuning$epsilon_bar
    if (tuning$reject_streak > 1000) {
      tuning$epsilon_curr <- 0.90 * tuning$epsilon_curr
      tuning$epsilon_bar <- 0.90 * tuning$epsilon_bar
    }
    
  }
  
  return(tuning)
  
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(nrow(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


M <- 1000
M.adapt <- 500
theta.m <- matrix(c(-4, 6))
epsilon <- find_initial_epsilon(theta.m, L=L.beta_pois, dL=dL.beta_pois)


nuts_sampler <- function(M, L, dL, epsilon, theta.m){
  
  tuning <- initialize_tuning(M.adapt, target=0.65, eps=epsilon)
  samples.theta <- array(NA, c(M, 2))
  
  
  Ls <- c()
  epsilons <- c()
  
  accept <- 0
  n.accept <- 0
  cat("Generating", M, " samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*M)
  for (m in 1:M){
    
    r.0 <- matrix(rnorm(2))
    u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
    theta.n <- theta.m
    theta.p <- theta.m
    r.n <- r.0
    r.p <- r.0
    j <- 0
    n <- 1
    s <- 1
    
    L.curr <- 0
    while (s == 1) {
      v.j <- sample(c(-1, 1), 1)
      if (v.j == -1){
        output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
        theta.n <- output$theta.n
        r.n <- output$r.n
        theta.prime <- output$theta.prime
        n.prime <- output$n.prime
        s.prime <- output$s.prime
        alpha <- output$alpha
        n.alpha <- output$n.p.alpha
      } else {
        output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
        theta.p <- output$theta.p
        r.p <- output$r.p
        theta.prime <- output$theta.prime
        n.prime <- output$n.prime
        s.prime <- output$s.prime
        alpha <- output$alpha
        n.alpha <- output$n.p.alpha
      }
      if (s.prime == 1){
        if (runif(1, 0, 1) < min(1, n.prime/n)){
          theta.m <- theta.prime
          accept <- accept + 1
        }
      }
      
      n <- n + n.prime
      s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      j <- j + 1
      
      L.curr <- L.curr + 1
      
    }
    
    n.accept <- n.accept + L.curr
    Ls <- c(Ls, L.curr)
    epsilons <- c(epsilons, epsilon)
    samples.theta[m,] <- t(theta.m)
    tuning <- update_tuning(tuning, alpha, m, n.alpha)
    
    if(m %in% percentage.points){
      setTxtProgressBar(progressBar, m/M)
    }
    
  }
  
  accept/n.accept
  output <- list()
  output$accept <- accept
  output$samples.theta <- samples.theta
  output$Ls <- Ls
  return(output)
  
}

out <- nuts_sampler(M, L=L.beta_pois, dL=dL.beta_pois, epsilon, theta.m)

colMeans(out$samples.theta)
sd(out$samples.theta[,1])
sd(out$samples.theta[,2])
plot(out$samples.theta[,1], type='l'); abline(h=Beta[1], col='2')
plot(out$samples.theta[,2], type='l'); abline(h=Beta[2], col='2')

hist(out$Ls)
summary(out$Ls)


######################
# Hamiltonian MCMC for 
# Logistic regression
######################


logit <- function(x){
  return(log(x) - log(1-x))
}


expit <- function(x){
  return(1/(1 + exp(-x)))
}


Beta <- c(-2, 2)
N <- 100
X <- cbind(rep(1, N), rnorm(N))
probs <- expit(X %*% Beta)
Y <- sapply(probs, function(x){rbinom(n=1, size=1, x)})


prior_mu <- rep(0, length(Beta))
prior_var <- rep(1000, length(Beta))
prior_sigma <- diag(prior_var)
prior_sigma_inv <- solve(prior_sigma)


L.beta_logit <- function(b){
  
  # likelihood
  logd <- 0
  lin_preds <- X %*% b
  for (i in 1:N){
    logd <- logd + dbinom(Y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  for (i in 1:length(b)){
    logd <- logd + dnorm(b[i], mean=prior_mu[i], sd=sqrt(prior_var[i]), log=T)
  }
  
  return(logd)
}


dL.beta_logit <- function(b){
  
  grad <- array(0, c(length(Beta), 1))
  lin_preds <- X %*% b
  for (i in 1:length(Y)){
    for (j in 1:length(b)){
      grad[j] <- grad[j] + X[i,j]*Y[i] - X[i,j]*expit(lin_preds[i,])
    }
  }
  
  # prior contribution
  grad <- grad - t(t(b) %*% prior_sigma_inv)
  
  return(t(grad))
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


initialize_tuning <- function(m, target, eps){
  return(list(
    M_adapt=m,
    epsilon_curr=eps,
    epsilon_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, n.a){
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    epsilon_bar <- tuning$epsilon_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a/n.a)
    log_epsilon_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_epsilon_bar <- i^{-kappa} * log_epsilon_curr + (1 - i^{-kappa}) * log(epsilon_bar)
    epsilon_curr <- exp(log_epsilon_curr)
    epsilon_bar <- exp(log_epsilon_bar)
    
    tuning$epsilon_bar <- epsilon_bar
    tuning$HbarM <- HbarM
    tuning$epsilon_curr <- epsilon_curr
    
  } else{
    
    tuning$epsilon_curr <- tuning$epsilon_bar
    if (tuning$reject_streak > 1000) {
      tuning$epsilon_curr <- 0.90 * tuning$epsilon_curr
      tuning$epsilon_bar <- 0.90 * tuning$epsilon_bar
    }
    
  }
  
  return(tuning)
  
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(nrow(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


nuts_sampler <- function(M, L, dL, epsilon, theta.m){
  
  tuning <- initialize_tuning(M.adapt, target=0.65, eps=epsilon)
  samples.theta <- array(NA, c(M, 2))
  
  
  Ls <- c()
  epsilons <- c()
  
  accept <- 0
  n.accept <- 0
  cat("Generating", M, " samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*M)
  for (m in 1:M){
    
    r.0 <- matrix(rnorm(2))
    u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
    theta.n <- theta.m
    theta.p <- theta.m
    r.n <- r.0
    r.p <- r.0
    j <- 0
    n <- 1
    s <- 1
    
    L.curr <- 0
    while (s == 1) {
      v.j <- sample(c(-1, 1), 1)
      if (v.j == -1){
        output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
        theta.n <- output$theta.n
        r.n <- output$r.n
        theta.prime <- output$theta.prime
        n.prime <- output$n.prime
        s.prime <- output$s.prime
        alpha <- output$alpha
        n.alpha <- output$n.p.alpha
      } else {
        output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
        theta.p <- output$theta.p
        r.p <- output$r.p
        theta.prime <- output$theta.prime
        n.prime <- output$n.prime
        s.prime <- output$s.prime
        alpha <- output$alpha
        n.alpha <- output$n.p.alpha
      }
      if (s.prime == 1){
        if (runif(1, 0, 1) < min(1, n.prime/n)){
          theta.m <- theta.prime
          accept <- accept + 1
        }
      }
      
      n <- n + n.prime
      s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      j <- j + 1
      
      L.curr <- L.curr + 1
      
    }
    
    n.accept <- n.accept + L.curr
    Ls <- c(Ls, L.curr)
    epsilons <- c(epsilons, epsilon)
    samples.theta[m,] <- t(theta.m)
    tuning <- update_tuning(tuning, alpha, m, n.alpha)
    
    if(m %in% percentage.points){
      setTxtProgressBar(progressBar, m/M)
    }
    
  }
  
  accept/n.accept
  output <- list()
  output$accept <- accept
  output$samples.theta <- samples.theta
  output$Ls <- Ls
  return(output)
  
}


M <- 2000
M.adapt <- 500
theta.m <- matrix(c(0, 0))
epsilon <- find_initial_epsilon(theta.m, L=L.beta_logit, dL=dL.beta_logit)


out <- nuts_sampler(M, L=L.beta_logit, dL=dL.beta_logit, epsilon, theta.m)

colMeans(out$samples.theta)
sd(out$samples.theta[,1])
sd(out$samples.theta[,2])
plot(out$samples.theta[,1], type='l'); abline(h=Beta[1], col='2')
plot(out$samples.theta[,2], type='l'); abline(h=Beta[2], col='2')
hist(out$Ls)
summary(out$Ls)


##############################
# Gaussian process in logistic
# model. Known parameters 
##############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
Sigma_inv <- solve(Sigma)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, loc.disc, seed=56)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)
Y <- locs$status


L.w <- function(w){
  
  # likelihood
  logd <- 0
  lin_preds <- w
  for (i in 1:N){
    logd <- logd + dbinom(Y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), mean=rep(0, length(w)), sigma=Sigma, log=T)
  
  return(logd)
}


dL.w <- function(w){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- w
  for (i in 1:length(Y)){
    grad[i] <- grad[i] + Y[i] - expit(lin_preds[i])
  }
  
  # prior contribution
  grad <- grad - t(t(w) %*% Sigma_inv)
  
  return(t(grad))
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


initialize_tuning <- function(m, target, eps){
  return(list(
    M_adapt=m,
    epsilon_curr=eps,
    epsilon_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, n.a){
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    epsilon_bar <- tuning$epsilon_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a/n.a)
    log_epsilon_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_epsilon_bar <- i^{-kappa} * log_epsilon_curr + (1 - i^{-kappa}) * log(epsilon_bar)
    epsilon_curr <- exp(log_epsilon_curr)
    epsilon_bar <- exp(log_epsilon_bar)
    
    tuning$epsilon_bar <- epsilon_bar
    tuning$HbarM <- HbarM
    tuning$epsilon_curr <- epsilon_curr
    
  } else{
    
    tuning$epsilon_curr <- tuning$epsilon_bar
    if (tuning$reject_streak > 1000) {
      tuning$epsilon_curr <- 0.90 * tuning$epsilon_curr
      tuning$epsilon_bar <- 0.90 * tuning$epsilon_bar
    }
    
  }
  
  return(tuning)
  
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(nrow(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


nuts_update <- function(L, dL, epsilon, theta.m){
  
  r.0 <- matrix(rnorm(nrow(theta.m)))
  u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r.0
  r.p <- r.0
  j <- 0
  n <- 1
  s <- 1
  
  accept <- 0
  L.curr <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
        accept <- accept + 1
      }
    }
    
    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    
    L.curr <- L.curr + 1
    
  }
  
  out <- list(
    theta.m=theta.m,
    accept=accept,
    L.curr=L.curr,
    alpha=alpha,
    n.alpha=n.alpha
  )
  return(out)
  
  
}


nuts_sampler <- function(M, L, dL, epsilon, theta.m){
  
  tuning <- initialize_tuning(M.adapt, target=0.65, eps=epsilon)
  samples.theta <- array(NA, c(M, nrow(theta.m)))
  
  Ls <- c()
  epsilons <- c()
  
  accept <- 0
  n.accept <- 0
  cat("Generating", M, " samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*M)
  for (m in 1:M){
    
    r.0 <- matrix(rnorm(nrow(theta.m)))
    u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
    theta.n <- theta.m
    theta.p <- theta.m
    r.n <- r.0
    r.p <- r.0
    j <- 0
    n <- 1
    s <- 1
    
    L.curr <- 0
    while (s == 1) {
      v.j <- sample(c(-1, 1), 1)
      if (v.j == -1){
        output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
        theta.n <- output$theta.n
        r.n <- output$r.n
        theta.prime <- output$theta.prime
        n.prime <- output$n.prime
        s.prime <- output$s.prime
        alpha <- output$alpha
        n.alpha <- output$n.p.alpha
      } else {
        output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
        theta.p <- output$theta.p
        r.p <- output$r.p
        theta.prime <- output$theta.prime
        n.prime <- output$n.prime
        s.prime <- output$s.prime
        alpha <- output$alpha
        n.alpha <- output$n.p.alpha
      }
      if (s.prime == 1){
        if (runif(1, 0, 1) < min(1, n.prime/n)){
          theta.m <- theta.prime
          accept <- accept + 1
        }
      }
      
      n <- n + n.prime
      s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      j <- j + 1
      
      L.curr <- L.curr + 1
      
    }
    
    n.accept <- n.accept + L.curr
    Ls <- c(Ls, L.curr)
    epsilons <- c(epsilons, epsilon)
    samples.theta[m,] <- t(theta.m)
    tuning <- update_tuning(tuning, alpha, m, n.alpha)
    
    if(m %in% percentage.points){
      setTxtProgressBar(progressBar, m/M)
    }
    
  }
  
  output <- list()
  output$accept <- accept/n.accept
  output$samples.theta <- samples.theta
  output$Ls <- Ls
  return(output)
  
}


nuts_sampler_new <- function(M, L, dL, epsilon, theta.m){
  
  Ls <- c()
  epsilons <- c()
  samples.theta <- array(NA, c(M, nrow(theta.m)))
  
  accept <- 0
  n.accept <- 0
  cat("Generating", M, " samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100)*M)
  for (m in 1:M){
    
    theta.out <- nuts_update(L, dL, epsilon, theta.m)
    theta.m <- theta.out$theta.m
    accept <- accept + theta.out$accept
    n.accept <- n.accept + theta.out$L.curr
    Ls <- c(Ls, theta.out$L.curr)
    samples.theta[m,] <- t(theta.m)

    if(m %in% percentage.points){
      setTxtProgressBar(progressBar, m/M)
    }
    
  }
  
  output <- list()
  output$accept <- accept/n.accept
  output$samples.theta <- samples.theta
  output$Ls <- Ls
  return(output)
  
}

M <- 1500
M.adapt <- 500
theta.m <- matrix(rnorm(length(W)))
epsilon <- find_initial_epsilon(theta.m, L=L.w, dL=dL.w)
epsilon <- 0.5

out <- nuts_sampler_new(M, L=L.w, dL=dL.w, epsilon, theta.m)

out$accept
hist(out$Ls)
summary(out$Ls)
w.hat <- colMeans(out$samples.theta)
plot(x=W, y=w.hat); abline(1, 1, col='2')
view_tr_w(out$samples.theta, W)


##############################
# Gaussian process in logistic
# model. Estimate other GP
# parameters. 
##############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
Sigma_inv <- solve(Sigma)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, loc.disc, seed=56)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)
Y <- locs$status


L.w <- function(w){
  
  # likelihood
  logd <- 0
  lin_preds <- w
  for (i in 1:N){
    logd <- logd + dbinom(Y[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), mean=rep(0, length(w)), sigma=sigma.m, log=T)
  
  return(logd)
}


dL.w <- function(w){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- w
  for (i in 1:length(Y)){
    grad[i] <- grad[i] + Y[i] - expit(lin_preds[i])
  }
  
  # prior contribution
  grad <- grad - t(t(w) %*% sigma.inv.m)
  
  return(t(grad))
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


initialize_tuning <- function(m, target, eps){
  return(list(
    M_adapt=m,
    epsilon_curr=eps,
    epsilon_bar=1,
    delta_tar=target,
    mu=log(10*0.05),
    HbarM=0,
    gamma=0.05,
    t0=10,
    kappa=0.75,
    reject_streak=0
  ))
}


update_tuning <- function(tuning, a, i, n.a){
  
  if (i <= tuning$M_adapt){
    
    delta_tar <- tuning$delta_tar
    epsilon_bar <- tuning$epsilon_bar
    HbarM <- tuning$HbarM
    gamma <- tuning$gamma
    mu <- tuning$mu
    kappa <- tuning$kappa
    t0 <- tuning$t0
    
    HbarM <- (1 - 1/(i + t0)) * HbarM + (1/(i + t0)) * (delta_tar - a/n.a)
    log_epsilon_curr <- mu - (sqrt(i)/gamma) * HbarM
    log_epsilon_bar <- i^{-kappa} * log_epsilon_curr + (1 - i^{-kappa}) * log(epsilon_bar)
    epsilon_curr <- exp(log_epsilon_curr)
    epsilon_bar <- exp(log_epsilon_bar)
    
    tuning$epsilon_bar <- epsilon_bar
    tuning$HbarM <- HbarM
    tuning$epsilon_curr <- epsilon_curr
    
  } else{
    
    tuning$epsilon_curr <- tuning$epsilon_bar
    if (tuning$reject_streak > 1000) {
      tuning$epsilon_curr <- 0.90 * tuning$epsilon_curr
      tuning$epsilon_bar <- 0.90 * tuning$epsilon_bar
    }
    
  }
  
  return(tuning)
  
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(nrow(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


nuts_update <- function(L, dL, epsilon, theta.m){
  
  r.0 <- matrix(rnorm(nrow(theta.m)))
  u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r.0
  r.p <- r.0
  j <- 0
  n <- 1
  s <- 1
  
  accept <- 0
  L.curr <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
        accept <- accept + 1
      }
    }
    
    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    
    L.curr <- L.curr + 1
    
  }
  
  out <- list(
    theta.m=theta.m,
    accept=accept,
    L.curr=L.curr,
    alpha=alpha,
    n.alpha=n.alpha
  )
  return(out)
  
  
}


M <- 3000
N <- length(W)

set.seed(422)
w.m <- matrix(rnorm(length(W)))
theta.m <- Theta + 2 * rnorm(1)
phi.m <- Phi + 2 * rnorm(1)

epsilon <- find_initial_epsilon(w.m, L=L.w, dL=dL.w)
epsilon <- 0.5

proposal.sd.theta <- 0.25
prior.phi <- c(3, 24)
prior.theta <- c(3, 2.5)

Ls <- c()
samples.w <- array(NA, c(M, nrow(w.m)))
samples.theta <- array(NA, c(M, 1))
samples.phi <- array(NA, c(M, 1))

accept <- rep(0, 2)
n.accept <- 0
cat("Generating", M, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*M)
for (m in 1:M){
  
  ## sample from w
  sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
  sigma.inv.m <- solve(sigma.m)
  w.out <- nuts_update(L.w, dL.w, epsilon, w.m)
  w.m <- w.out$theta.m
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.m, as.numeric(w.m), d, phi.m, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.m <- theta.out$theta
  
  ## sample from phi
  R.m <- sigma.m/phi.m
  phi.m <- 1/rgamma(1, N/2 + prior.phi[1], t(w.m) %*% solve(R.m) %*% w.m/2 + prior.phi[2])
  
  accept[1] <- accept[1] + w.out$accept
  accept[2] <- accept[2] + theta.out$accept
  n.accept <- n.accept + w.out$L.curr
  Ls <- c(Ls, w.out$L.curr)
  
  samples.w[m,] <- t(w.m)
  samples.phi[m,] <- phi.m
  samples.theta[m,] <- theta.m
  
  if(m %in% percentage.points){
    setTxtProgressBar(progressBar, m/M)
  }
  
}

accept[1] <- accept[1]/n.accept
accept[2] <- accept[2]/M
print(accept)

hist(Ls)
summary(Ls)

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(1, 1, col='2')
view_tr_w(samples.w, W)

view_tr(samples.theta, Theta)
view_tr(samples.phi, Phi)


##############################
# Gaussian process in logistic
# model. Estimate other GP
# parameters. 
# Joint update w/ case data
##############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
Sigma_inv <- solve(Sigma)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, loc.disc, seed=56)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)
Y.l <- locs$status


## Case Counts
Alpha.case <- 1
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W)
glm(case.data$y ~ case.data$x.standardised-1, family='poisson')
prior_alpha_ca_mean <- 2
prior_alpha_ca_var <- 6
y.ca <- case.data$y
x.c <- case.data$x.standardised


L.w <- function(w, y.l, x.c, beta.ca, y.ca, alpha.ca, sigma.m){
  
  alpha.ca <- as.numeric(alpha.ca)
  
  # likelihood
  logd <- 0
  lin_preds <- w
  for (i in 1:N){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- as.numeric(w)[as.logical(locs$status)]
  count_pred <- x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), mean=rep(0, length(w)), sigma=sigma.m, log=T)
  
  return(logd)
}


dL.w <- function(w, y.l, x.c, beta.ca, y.ca, alpha.ca, sigma.inv.m){
  
  alpha.ca <- as.numeric(alpha.ca)
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- w
  for (i in 1:length(Y.l)){
    grad[i] <- grad[i] + y.l[i] - expit(lin_preds[i])
  }
  
  # case contribution
  lin.count <- x.c %*% beta.ca
  for (i in 1:length(locs$ids)){
    id.i <- locs$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # prior contribution
  grad <- grad - t(t(w) %*% sigma.inv.m)
  
  return(t(grad))
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(nrow(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


nuts_update <- function(L, dL, epsilon, theta.m){
  
  r.0 <- matrix(rnorm(nrow(theta.m)))
  u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r.0
  r.p <- r.0
  j <- 0
  n <- 1
  s <- 1
  
  accept <- 0
  L.curr <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
        accept <- accept + 1
      }
    }
    
    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    
    L.curr <- L.curr + 1
    
  }
  
  out <- list(
    theta.m=theta.m,
    accept=accept,
    L.curr=L.curr,
    alpha=alpha,
    n.alpha=n.alpha
  )
  return(out)
  
  
}


M <- 2000
N <- length(W)

set.seed(422)
w.m <- matrix(rnorm(length(W)))
theta.m <- Theta + 2 * rnorm(1)
phi.m <- Phi + 2 * rnorm(1)
alpha.m <- Alpha.case
beta.m <- beta.case

sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
sigma.inv.m <- solve(sigma.m)
# epsilon <- find_initial_epsilon(w.m, L=L.w, dL=dL.w)
epsilon <- 0.08

proposal.sd.theta <- 0.25
prior.phi <- c(3, 24)
prior.theta <- c(3, 2.5)

Ls <- c()
samples.w <- array(NA, c(M, nrow(w.m)))
samples.theta <- array(NA, c(M, 1))
samples.phi <- array(NA, c(M, 1))

accept <- rep(0, 2)
n.accept <- 0
cat("Generating", M, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*M)
for (m in 1:M){
  
  ## sample from w
  sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
  sigma.inv.m <- solve(sigma.m)
  L.w_wrap <- function(w){
    return(L.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.m))
  }
  dL.w_wrap <- function(w){
    return(dL.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.inv.m))
  }
  w.out <- nuts_update(L.w_wrap, dL.w_wrap, epsilon, w.m)
  w.m <- w.out$theta.m
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.m, as.numeric(w.m), d, phi.m, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.m <- theta.out$theta
  
  ## sample from phi
  R.m <- sigma.m/phi.m
  phi.m <- 1/rgamma(1, N/2 + prior.phi[1], t(w.m) %*% solve(R.m) %*% w.m/2 + prior.phi[2])
  
  accept[1] <- accept[1] + w.out$accept
  accept[2] <- accept[2] + theta.out$accept
  n.accept <- n.accept + w.out$L.curr
  Ls <- c(Ls, w.out$L.curr)
  
  samples.w[m,] <- t(w.m)
  samples.phi[m,] <- phi.m
  samples.theta[m,] <- theta.m
  
  if(m %in% percentage.points){
    setTxtProgressBar(progressBar, m/M)
  }
  
}

accept[1] <- accept[1]/n.accept
accept[2] <- accept[2]/M
print(accept)

hist(Ls)
summary(Ls)
plot(Ls, type='l')

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(1, 1, col='2')
view_tr_w(samples.w, W)

view_tr(samples.theta, Theta)
view_tr(samples.phi, Phi)


##############################
# Gaussian process in logistic
# model. Estimate other GP
# parameters. 
# Joint update w/ case data.
# Estimate betas. 
##############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
Sigma_inv <- solve(Sigma)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
beta.samp <- c(-1, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, loc.disc, seed=56)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)
Y.l <- locs$status


## Case Counts
Alpha.case <- 1
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W)
glm(case.data$y ~ case.data$x.standardised-1, family='poisson')
prior_alpha_ca_mean <- 1
prior_alpha_ca_var <- 6
y.ca <- case.data$y
x.c <- case.data$x.standardised


L.w <- function(w, y.l, x.c, beta.ca, y.ca, alpha.ca, sigma.m){
  
  # likelihood
  logd <- 0
  lin_preds <- w
  for (i in 1:N){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- as.numeric(w)[as.logical(locs$status)]
  count_pred <- x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), mean=rep(0, length(w)), sigma=sigma.m, log=T)
  
  return(logd)
}


dL.w <- function(w, y.l, x.c, beta.ca, y.ca, alpha.ca, sigma.inv.m){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- w
  for (i in 1:length(Y.l)){
    grad[i] <- grad[i] + y.l[i] - expit(lin_preds[i])
  }
  
  # case contribution
  lin.count <- x.c %*% beta.ca
  for (i in 1:length(locs$ids)){
    id.i <- locs$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # prior contribution
  grad <- grad - t(t(w) %*% sigma.inv.m)
  
  return(t(grad))
}


L.beta <- function(w, x.c, beta.ca, y.ca, alpha.ca){
  
  # likelihood
  logd <- 0
  lin_preds <- x.c %*% beta.ca + alpha.ca * w
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior 
  logd <- logd + dmvnorm(as.numeric(beta.ca), mean=rep(0, nrow(beta.ca)), sigma=diag(rep(10, nrow(beta.ca))), log=T)
  
  return(logd)
}


dL.beta <- function(w, x, beta, y, alpha){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + alpha * w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - exp(lin_preds[i,]))
    }
  }
  
  if (length(beta) > 1){
    prior_var_inv <- diag(rep(1/10, length(beta)))
  } else{
    prior_var_inv <- 1/10
  }
  
  grad <- grad + t(-t(beta) %*% prior_var_inv)
  return(t(grad))
  
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(nrow(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


nuts_update <- function(L, dL, epsilon, theta.m){
  
  r.0 <- matrix(rnorm(nrow(theta.m)))
  u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r.0
  r.p <- r.0
  j <- 0
  n <- 1
  s <- 1
  
  accept <- 0
  L.curr <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
        accept <- accept + 1
      }
    }
    
    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    
    L.curr <- L.curr + 1
    
  }
  
  out <- list(
    theta.m=theta.m,
    accept=accept,
    L.curr=L.curr,
    alpha=alpha,
    n.alpha=n.alpha
  )
  return(out)
  
  
}


M <- 1000
N <- length(W)

set.seed(422)
w.m <- matrix(rnorm(length(W)))
theta.m <- Theta + 2 * rnorm(1)
phi.m <- Phi + 2 * rnorm(1)
alpha.m <- Alpha.case
beta.m <- matrix(rnorm(2))

sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
sigma.inv.m <- solve(sigma.m)
L.w_wrap <- function(w){
  return(L.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.m))
}
dL.w_wrap <- function(w){
  return(dL.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.inv.m))
}
find_initial_epsilon(w.m, L=L.w_wrap, dL=dL.w_wrap)
epsilon <- 0.08

w.sub <- w.m[as.logical(locs$status)]
L.beta_wrap <- function(beta){
  return(L.beta(w.sub, x.c, beta, y.ca, alpha.m))
}
dL.beta_wrap <- function(beta){
  return(dL.beta(w.sub, x.c, beta, y.ca, alpha.m))
}
find_initial_epsilon(beta.m, L=L.beta_wrap, dL=dL.beta_wrap)
epsilon.b <- 0.03

proposal.sd.theta <- 0.25
prior.phi <- c(3, 24)
prior.theta <- c(3, 2.5)

Ls <- c()
Ls.beta <- c()
samples.w <- array(NA, c(M, nrow(w.m)))
samples.theta <- array(NA, c(M, 1))
samples.phi <- array(NA, c(M, 1))
samples.beta <- array(NA, c(M, length(beta.case)))

accept <- rep(0, 3)
n.accept <- rep(0, 3)
cat("Generating", M, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*M)
for (m in 1:M){
  
  ## sample from w
  sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
  sigma.inv.m <- solve(sigma.m)
  L.w_wrap <- function(w){
    return(L.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.m))
  }
  dL.w_wrap <- function(w){
    return(dL.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.inv.m))
  }
  w.out <- nuts_update(L.w_wrap, dL.w_wrap, epsilon, w.m)
  w.m <- w.out$theta.m
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.m, as.numeric(w.m), d, phi.m, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.m <- theta.out$theta
  
  ## sample from phi
  R.m <- sigma.m/phi.m
  phi.m <- 1/rgamma(1, N/2 + prior.phi[1], t(w.m) %*% solve(R.m) %*% w.m/2 + prior.phi[2])
  
  ## sample from beta
  w.sub <- w.m[as.logical(locs$status)]
  L.beta_wrap <- function(beta){
    return(L.beta(w.sub, x.c, beta, y.ca, Alpha.case))
  }
  dL.beta_wrap <- function(beta){
    return(dL.beta(w.sub, x.c, beta, y.ca, Alpha.case))
  }
  beta.out.ca <- nuts_update(L.beta_wrap, dL.beta_wrap, epsilon.b, beta.m)
  beta.m <- beta.out.ca$theta.m
  
  accept[1] <- accept[1] + w.out$accept
  accept[2] <- accept[2] + theta.out$accept
  accept[3] <- accept[3] + beta.out.ca$accept
  n.accept[1] <- n.accept[1] + w.out$L.curr
  n.accept[2] <- n.accept[2] + 1
  n.accept[3] <- n.accept[3] + beta.out.ca$L.curr
  Ls <- c(Ls, w.out$L.curr)
  Ls.beta <- c(Ls.beta, beta.out.ca$L.curr)
  
  samples.w[m,] <- t(w.m)
  samples.phi[m,] <- phi.m
  samples.theta[m,] <- theta.m
  samples.beta[m,] <- t(beta.m)
  
  if(m %in% percentage.points){
    setTxtProgressBar(progressBar, m/M)
  }
  
}

accept <- accept/n.accept
print(accept)

hist(Ls)
summary(Ls)
plot(Ls, type='l')

hist(Ls.beta)
summary(Ls.beta)
plot(Ls.beta, type='l')

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(1, 1, col='2')
view_tr_w(samples.w, W)

view_tr(samples.theta, Theta)
view_tr(samples.phi, Phi)

view_tr(samples.beta[,1], beta.case[1])
view_tr(samples.beta[,2], beta.case[2])


##############################
# Gaussian process in logistic
# model. Estimate other GP
# parameters. 
# Joint update w/ case data.
# Estimate betas. Estimate 
# alpha. 
##############################


library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
Sigma_inv <- solve(Sigma)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
loc.disc <- caWc.disc[[c(6)]]
locs <- simLocW(W, loc.disc, seed=56)
sum(locs$status)
hist(W)
plot(loc.disc)
points(locs$coords)
Y.l <- locs$status


## Case Counts
Alpha.case <- 1
beta.case <- c(1, 1)
cov.disc <- caWc.disc[[c(1)]]
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W)
glm(case.data$y ~ case.data$x.standardised-1, family='poisson')
prior_alpha_ca_mean <- 1
prior_alpha_ca_var <- 10
y.ca <- case.data$y
x.c <- case.data$x.standardised


L.w <- function(w, y.l, x.c, beta.ca, y.ca, alpha.ca, sigma.m){
  
  # likelihood
  logd <- 0
  lin_preds <- w
  for (i in 1:N){
    logd <- logd + dbinom(y.l[i], size=1, prob=expit(lin_preds[i]), log=T)
  }
  
  # likelihood: case counts
  w.sub <- as.numeric(w)[as.logical(locs$status)]
  count_pred <- x.c %*% beta.ca + alpha.ca * w.sub
  rates <- exp(count_pred)
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=rates[i], log=T)
  }
  
  # prior
  logd <- logd + dmvnorm(as.numeric(w), mean=rep(0, length(w)), sigma=sigma.m, log=T)
  
  return(logd)
}


dL.w <- function(w, y.l, x.c, beta.ca, y.ca, alpha.ca, sigma.inv.m){
  
  grad <- array(0, c(length(w), 1))
  lin_preds <- w
  for (i in 1:length(Y.l)){
    grad[i] <- grad[i] + y.l[i] - expit(lin_preds[i])
  }
  
  # case contribution
  lin.count <- x.c %*% beta.ca
  for (i in 1:length(locs$ids)){
    id.i <- locs$ids[i]
    grad[id.i] <- grad[id.i] + alpha.ca * (y.ca[i] - exp(lin.count[i] + alpha.ca * w[id.i]))
  }
  
  # prior contribution
  grad <- grad - t(t(w) %*% sigma.inv.m)
  
  return(t(grad))
}


L.beta <- function(w, x.c, beta.ca, y.ca, alpha.ca){
  
  # likelihood
  logd <- 0
  lin_preds <- x.c %*% beta.ca + alpha.ca * w
  for (i in 1:length(y.ca)){
    logd <- logd + dpois(y.ca[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior 
  logd <- logd + dmvnorm(as.numeric(beta.ca), mean=rep(0, nrow(beta.ca)), sigma=diag(rep(10, nrow(beta.ca))), log=T)
  
  return(logd)
}


dL.beta <- function(w, x, beta, y, alpha){
  
  grad <- array(0, c(length(beta), 1))
  lin_preds <- x %*% beta + alpha * w
  
  for (j in 1:length(beta)){
    for (i in 1:length(w)){
      grad[j] <- grad[j] + x[i, j] * (y[i] - exp(lin_preds[i,]))
    }
  }
  
  if (length(beta) > 1){
    prior_var_inv <- diag(rep(1/10, length(beta)))
  } else{
    prior_var_inv <- 1/10
  }
  
  grad <- grad + t(-t(beta) %*% prior_var_inv)
  return(t(grad))
  
}


L.alpha <- function(w, x, beta, y, alpha, prior_mean, prior_var){
  
  # likelihood
  logd <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(y)){
    logd <- logd + dpois(y[i], lambda=exp(lin_preds[i]), log=T)
  }
  
  # prior
  logd <- logd + dnorm(alpha, prior_mean, prior_var, log=T)
  
  return(logd)
  
}


dL.alpha <- function(w, x, beta, y, alpha, prior_mean, prior_var){
  
  grad <- 0
  lin_preds <- x %*% beta + alpha * w
  for (i in 1:length(w)){
    grad <- grad + w[i] * (y[i] - exp(lin_preds[i,]))
  }
  
  grad <- grad - (alpha) / prior_var
  
  return(grad)
  
}


K <- function(p){
  return(sum(t(p) %*% p)/2)
}


leapfrog <- function(theta, r, epsilon, dL){
  r <- r + (epsilon/2) * t(dL(theta))
  theta <- theta + epsilon * r
  r <- r + (epsilon/2) * t(dL(theta))
  output <- list(
    r=r,
    theta=theta
  )
  return(output)
}


build_tree <- function(theta, r, u, v, j, epsilon, theta.0, r.0, L, dL){
  
  if (j == 0) {
    lf <- leapfrog(theta, r, v * epsilon, dL=dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    L.prime <- L(theta.prime)
    K.prime <- K(r.prime)
    L.0 <- L(theta.0)
    K.0 <- K(r.0)
    n.prime <- as.numeric(u <= exp(L.prime - K.prime))
    s.prime <- as.numeric(u <= exp(1000 + L.prime - K.prime))
    alpha <- min(1, exp(L.prime - K.prime - L.0 + K.0))
    if (is.null(alpha)){
      alpha <- 0
    }
    return(
      list(
        theta.n=theta.prime,
        r.n=r.prime,
        theta.p=theta.prime,
        r.p=r.prime,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha,
        n.p.alpha=1
      )
    )
  } else{
    bt <- build_tree(theta, r, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
    theta.n <- bt$theta.n
    r.n <- bt$r.n
    theta.p <- bt$theta.p
    r.p <- bt$r.p
    theta.prime <- bt$theta.prime
    n.prime <- bt$n.prime
    s.prime <- bt$s.prime
    alpha.prime <- bt$alpha
    n.p.alpha <- bt$n.p.alpha
    if (s.prime == 1){
      
      if (v == -1){
        bt <- build_tree(theta.n, r.n, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.n <- bt$theta.n
        r.n <- bt$r.n
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      } else {
        bt <- build_tree(theta.p, r.p, u, v, j-1, epsilon, theta.0, r.0, L=L, dL=dL)
        theta.p <- bt$theta.p
        r.p <- bt$r.p
        theta.pp <- bt$theta.prime
        n.pp <- bt$n.prime
        s.pp <- bt$s.prime
        alpha.pp <- bt$alpha
        n.pp.alpha <- bt$n.p.alpha
      }
      if (n.pp + n.prime > 0){
        if (runif(1, 0, 1) < n.pp/(n.pp + n.prime)){
          theta.prime <- theta.pp
        }
      }
      alpha.prime <- alpha.prime + alpha.pp
      n.p.alpha <- n.p.alpha + n.pp.alpha
      s.prime <- s.pp * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
      n.prime <- n.prime + n.pp
      
    }
    return(
      list(
        theta.n=theta.n,
        r.n=r.n,
        theta.p=theta.p,
        r.p=r.p,
        theta.prime=theta.prime,
        n.prime=n.prime,
        s.prime=s.prime,
        alpha=alpha.prime,
        n.p.alpha=n.p.alpha
      )
    )
  }
}


p_ratio <- function(L, r1, theta1, r2, theta2){
  
  return(exp(L(theta1) - K(r1) - L(theta2) + K(r2)))
  
}


find_initial_epsilon <- function(theta, L, dL){
  
  epsilon <- 1
  r <- rnorm(length(theta))
  lf <- leapfrog(theta, r, epsilon, dL)
  theta.prime <- lf$theta
  r.prime <- lf$r
  pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  a <- 2 * as.numeric(pr > 0.5) - 1
  while ((pr^a) > 2^(-a)){
    epsilon <- epsilon * 2^a
    lf <- leapfrog(theta, r, epsilon, dL)
    theta.prime <- lf$theta
    r.prime <- lf$r
    pr <- p_ratio(L, r.prime, theta.prime, r, theta)
  }
  return(epsilon)
  
}


nuts_update <- function(L, dL, epsilon, theta.m){
  
  r.0 <- matrix(rnorm(length(theta.m)))
  u <- runif(1, 0, exp(L(theta.m) - K(r.0)))
  theta.n <- theta.m
  theta.p <- theta.m
  r.n <- r.0
  r.p <- r.0
  j <- 0
  n <- 1
  s <- 1
  
  accept <- 0
  L.curr <- 0
  while (s == 1) {
    v.j <- sample(c(-1, 1), 1)
    if (v.j == -1){
      output <- build_tree(theta.n, r.n, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.n <- output$theta.n
      r.n <- output$r.n
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    } else {
      output <- build_tree(theta.p, r.p, u, v.j, j, epsilon, theta.m, r.0, L=L, dL=dL)
      theta.p <- output$theta.p
      r.p <- output$r.p
      theta.prime <- output$theta.prime
      n.prime <- output$n.prime
      s.prime <- output$s.prime
      alpha <- output$alpha
      n.alpha <- output$n.p.alpha
    }
    if (is.na(s.prime)){
      s.prime <- 0
    }
    if (s.prime == 1){
      if (runif(1, 0, 1) < min(1, n.prime/n)){
        theta.m <- theta.prime
        accept <- accept + 1
      }
    }
    
    n <- n + n.prime
    s <- s.prime * as.numeric( t(theta.p - theta.n) %*% r.n >= 0 ) * as.numeric( t(theta.p - theta.n) %*% r.p >= 0 )
    j <- j + 1
    if (is.na(s)){
      s <- 0
    }
    
    L.curr <- L.curr + 1
    
  }
  
  out <- list(
    theta.m=theta.m,
    accept=accept,
    L.curr=L.curr,
    alpha=alpha,
    n.alpha=n.alpha
  )
  return(out)
  
  
}


M <- 1000
N <- length(W)

set.seed(422)
w.m <- matrix(rnorm(length(W)))
theta.m <- Theta + 2 * rnorm(1)
phi.m <- Phi + 2 * rnorm(1)
alpha.m <- Alpha.case
beta.m <- matrix(rnorm(2))

sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
sigma.inv.m <- solve(sigma.m)
L.w_wrap <- function(w){
  return(L.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.m))
}
dL.w_wrap <- function(w){
  return(dL.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.inv.m))
}
find_initial_epsilon(w.m, L=L.w_wrap, dL=dL.w_wrap)
epsilon <- 0.1

w.sub <- w.m[as.logical(locs$status)]
L.beta_wrap <- function(beta){
  return(L.beta(w.sub, x.c, beta, y.ca, alpha.m))
}
dL.beta_wrap <- function(beta){
  return(dL.beta(w.sub, x.c, beta, y.ca, alpha.m))
}
find_initial_epsilon(beta.m, L=L.beta_wrap, dL=dL.beta_wrap)
epsilon.b <- 0.03

L.alpha_wrap <- function(alpha){
  return(L.alpha(w.sub, x.c, beta.m, y.ca, alpha, prior_alpha_ca_mean, prior_alpha_ca_var))
}
dL.alpha_wrap <- function(alpha){
  return(L.alpha(w.sub, x.c, beta.m, y.ca, alpha, prior_alpha_ca_mean, prior_alpha_ca_var))
}
find_initial_epsilon(alpha.m, L=L.alpha_wrap, dL=dL.alpha_wrap)
epsilon.a <- 0.0005


proposal.sd.theta <- 0.25
prior.phi <- c(3, 24)
prior.theta <- c(3, 2.5)

Ls <- c()
Ls.beta <- c()
Ls.alpha <- c()
samples.w <- array(NA, c(M, nrow(w.m)))
samples.theta <- array(NA, c(M, 1))
samples.phi <- array(NA, c(M, 1))
samples.beta <- array(NA, c(M, length(beta.case)))
samples.alpha <- array(NA, c(M, 1))

accept <- rep(0, 4)
n.accept <- rep(0, 4)
cat("Generating", M, " samples\n")
progressBar <- txtProgressBar(style = 3)
percentage.points <- round((1:100/100)*M)
for (m in 1:M){
  
  ## sample from w
  sigma.m <- Exponential(d, range=theta.m, phi=phi.m)
  sigma.inv.m <- solve(sigma.m)
  L.w_wrap <- function(w){
    return(L.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.m))
  }
  dL.w_wrap <- function(w){
    return(dL.w(w, Y.l, x.c, beta.m, y.ca, alpha.m, sigma.inv.m))
  }
  w.out <- nuts_update(L.w_wrap, dL.w_wrap, epsilon, w.m)
  w.m <- w.out$theta.m
  
  ## sample from theta
  theta.out <- rangeMhUpdate(theta.m, as.numeric(w.m), d, phi.m, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  theta.m <- theta.out$theta
  
  ## sample from phi
  R.m <- sigma.m/phi.m
  phi.m <- 1/rgamma(1, N/2 + prior.phi[1], t(w.m) %*% solve(R.m) %*% w.m/2 + prior.phi[2])
  
  ## sample from beta
  w.sub <- w.m[as.logical(locs$status)]
  L.beta_wrap <- function(beta){
    return(L.beta(w.sub, x.c, beta, y.ca, Alpha.case))
  }
  dL.beta_wrap <- function(beta){
    return(dL.beta(w.sub, x.c, beta, y.ca, Alpha.case))
  }
  beta.out.ca <- nuts_update(L.beta_wrap, dL.beta_wrap, epsilon.b, beta.m)
  beta.m <- beta.out.ca$theta.m
  
  ## sample from alpha
  # L.alpha_wrap <- function(alpha){
  #   return(L.alpha(w.sub, x.c, beta.m, y.ca, alpha, prior_alpha_ca_mean, prior_alpha_ca_var))
  # }
  # dL.alpha_wrap <- function(alpha){
  #   return(L.alpha(w.sub, x.c, beta.m, y.ca, alpha, prior_alpha_ca_mean, prior_alpha_ca_var))
  # }
  # alpha.out.ca <- nuts_update(L.alpha_wrap, dL.alpha_wrap, epsilon.a, alpha.m)
  # alpha.m <- alpha.out.ca$theta.m
  alpha.out.ca <- alphaHmcUpdate(y.ca, w.sub, x.c, beta.m, alpha.m,
                                 delta_a=5e-4, prior_alpha_ca_mean, prior_alpha_ca_var, L_a=4)
  alpha.m <- alpha.out.ca$alpha
  # alpha.m <- Alpha.case
  
  accept[1] <- accept[1] + w.out$accept
  accept[2] <- accept[2] + theta.out$accept
  accept[3] <- accept[3] + beta.out.ca$accept
  accept[4] <- accept[4] + alpha.out.ca$accept
  n.accept[1] <- n.accept[1] + w.out$L.curr
  n.accept[2] <- n.accept[2] + 1
  n.accept[3] <- n.accept[3] + beta.out.ca$L.curr
  n.accept[4] <- n.accept[4] + 1 # alpha.out.ca$L.curr
  Ls <- c(Ls, w.out$L.curr)
  Ls.beta <- c(Ls.beta, beta.out.ca$L.curr)
  # Ls.alpha <- c(Ls.alpha, alpha.out.ca$L.curr)
  
  samples.w[m,] <- t(w.m)
  samples.phi[m,] <- phi.m
  samples.theta[m,] <- theta.m
  samples.beta[m,] <- t(beta.m)
  samples.alpha[m,] <- alpha.m
  
  if(m %in% percentage.points){
    setTxtProgressBar(progressBar, m/M)
  }
  
}

accept <- accept/n.accept
print(accept)

hist(Ls)
summary(Ls)
plot(Ls, type='l')

hist(Ls.beta)
summary(Ls.beta)
plot(Ls.beta, type='l')

w.hat <- colMeans(samples.w)
plot(x=W, y=w.hat); abline(1, 1, col='2')
view_tr_w(samples.w, W)

view_tr(samples.theta, Theta)
view_tr(samples.phi, Phi)

view_tr(samples.beta[,1], beta.case[1])
view_tr(samples.beta[,2], beta.case[2])

view_tr(samples.alpha)
 