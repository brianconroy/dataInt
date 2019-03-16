library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
sourceDirectory('Documents/research/dataInt/R/')


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=2)
caWc.disc <- simRegion$raster


#### Simulate gaussian process
Theta <- 5
Phi <- 3
cells.all <- c(1:ncell(caWc.disc))[!is.na(values(caWc.disc[[1]]))]
coords <- xyFromCell(caWc.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
w <- mvrnorm(n=1, mu=rep(0, length(cells.all)), sigma)

ids_high <- as.numeric(names(sort(w, decreasing=T)[1:50]))
w_obs <- w[ids_high]
w_not <- w[-ids_high]
set.seed(42)
ids_random <- sample(1:length(w), size=50, replace=F)
w_obs2 <- w[ids_random]
w_not2 <- w[-ids_random]

conditional_mean <- function(ids, covmat, w_observed){
  
  omega.11 <- covmat[-ids, -ids]
  omega.12 <- covmat[-ids, ids]
  omega.21 <- covmat[ids, -ids]
  omega.22 <- covmat[ids, ids]
  
  e.cond <- omega.12 %*% solve(omega.22) %*% w_observed
  var.cond <- omega.11 - omega.12 %*% solve(omega.22) %*% omega.21
  
  return(e.cond)
  
}

w_est1 <- conditional_mean(ids_high, sigma, w_obs)
w_est2 <- conditional_mean(ids_random, sigma, w_obs2)

par(mfrow=c(1,2))
plot(x=w_not, y=w_est1, xlim=c(-3,3), main='A)', xlab='True W', ylab='Predicted W'); abline(0,1,col=2)
plot(x=w_not, y=w_est2, xlim=c(-3.5,3.5), main='B)', xlab='True W', ylab='Predicted W'); abline(0,1,col=2)
