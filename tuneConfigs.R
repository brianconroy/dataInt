

###################################
### D = 1, unthinned, estimate rho
###################################
t.beta <- 0.015
t.phi <- 1.5
t.rho <- 1
tau.sq <- 0.5
rho <- 0.75
beta0 <- 0.5
beta1 <- 3
beta2 <- -1


###################################
### D = 1, unthinned, fixed rho
###################################
t.beta <- 0.015
t.phi <- 1.5
t.rho <- 1
tau.sq <- 0.5
rho <- 0.75
beta0 <- 0.5
beta1 <- 3
beta2 <- -1


###################################
### D = 1, separate CARs, unthinned
###################################
tau.sq <- 0.5
rho.1 <- 0.75
beta0 <- 0.5
beta1 <- 3
tune.beta <- 0.025
tune.phi <- 1
tune.tau.prior <- list(c(1, 0.01))


###################################
### D = 2, separate CARs, unthinned
###################################
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
tune.beta <- 0.01 
tune.phi <- 0.5
tune.tau.prior <- list(c(1, 0.01))


###################################
### D = 2, separate CARs, unthinned,
### p = 4
###################################
caWc.disc <- caWc.disc[[c(1, 3, 4)]]
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75 
rho.list <- c(rho.1, rho.2)
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta3 <- -0.5
beta.true <- c(beta0, beta1, beta2, beta3)
tune.beta <- 0.025
tune.phi <- 1


###################################
### D = 3, separate CARs, unthinned
###################################
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
tau.sq.3 <- 0.75
rho.1 <- 0.25
rho.2 <- 0.5
rho.3 <- 0.75
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
t.beta <- 0.025
t.phi <- 0.15


###################################
### D = 1 (spatial) + 1 (survey),
#   unthinned, works well
###################################
caWc.disc <- caWc.disc[[c(1)]]
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75 
rho.list <- c(rho.1, rho.2)
beta0 <- 1
beta1 <- 2
beta.true <- c(beta0, beta1)
data.1 <- simLeroux(caWc.disc, beta.true, tau.sq.1, rho.1, nNeighbors, seed=405)
data.unif <- simUniform(caWc.disc, beta.true, nsamp=25, seed=505)
data <- list(spatial=list(data.1), survey=data.unif)
t.beta <- 0.01
t.phi <- 1.5


###################################
### D = 2 (spatial) + 1 (survey),
#   unthinned, works well
###################################
caWc.disc <- caWc.disc[[c(1)]]
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.5
rho.list <- c(rho.1, rho.2)
beta0 <- 1
beta1 <- 2
beta.true <- c(beta0, beta1)
data.1 <- simLeroux(caWc.disc, beta.true, tau.sq.1, rho.1, nNeighbors, seed=405)
data.2 <- simLeroux(caWc.disc, beta.true, tau.sq.2, rho.2, nNeighbors, seed=404)
data.unif <- simUniform(caWc.disc, beta.true, nsamp=50, seed=505)
data <- list(spatial=list(data.1, data.2), survey=data.unif)
t.beta <- 0.1
t.phi <- 0.5


###################################
### D = 1, thinned, sort of works
###################################
tau.sq.1 <- 0.5
rho.1 <- 0.75
beta0 <- 0.5
beta1 <- 2
gamma.1 <- -0.5
tune.beta <- c(0.025, 0.03, 0.04)
tune.phi <- c(0.025, 0.5)
tune.tau.prior <- list(c(1, 0.01))


###################################
### thinned, D=2, p=5, fixed rho,
#   pretty much works
###################################
caWc.disc <- caWc.disc[[c(1)]]
nNeighbors <- simRegion$nNeighbors
samp.disc.1 <- simRegion$raster[[c(3)]]
samp.disc.2 <- simRegion$raster[[c(4)]]
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75 
rho.list <- c(rho.1, rho.2)
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1)
gamma.1 <- -0.5
gamma.2 <- -0.25
gamma.true <- c(gamma.1, gamma.2)
data.1 <- simLerouxThinned(caWc.disc, samp.disc.1, beta.true, gamma.1, tau.sq.1, rho.1, nNeighbors, seed=401)
data.2 <- simLerouxThinned(caWc.disc, samp.disc.2, beta.true, gamma.2, tau.sq.2, rho.2, nNeighbors, seed=401)
tune.beta <- 0.05
tune.phi <- 1.5


########################################
### thinned, D=3 (2 spatial, 1 uniform), 
#   p=5, fixed rho, pretty much works
########################################
caWc.disc <- caWc.disc[[c(1)]]
samp.disc.1 <- simRegion$raster[[c(3)]]
samp.disc.2 <- simRegion$raster[[c(4)]]
tau.sq.1 <- 0.25
tau.sq.2 <- 0.5
rho.1 <- 0.25
rho.2 <- 0.75 
rho.list <- c(rho.1, rho.2)
beta0 <- 0.5
beta1 <- 3
beta2 <- -1
beta.true <- c(beta0, beta1)
gamma.1 <- -0.5
gamma.2 <- -0.25
gamma.true <- c(gamma.1, gamma.2)
data.1 <- simLerouxThinned(caWc.disc, samp.disc.1, beta.true, gamma.1, tau.sq.1, rho.1, nNeighbors, seed=401)
data.2 <- simLerouxThinned(caWc.disc, samp.disc.2, beta.true, gamma.2, tau.sq.2, rho.2, nNeighbors, seed=401)
data.unif <- simUniform(caWc.disc, beta.true, nsamp=25, seed=505)

