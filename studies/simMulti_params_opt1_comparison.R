###########################
# Set simulation parameters
###########################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


# Scenario 1: 
#   identical gaussian process
#   when the disease covariates/sampling params
#   are moderately different

#### Load simulation parameters
perturb_factor <- 1.75 #1.05  1.75
perturb <- "high"
ws <- "different"
  
Alpha.case1 <- 0.5
Alpha.case2 <- Alpha.case1 * perturb_factor
Alpha.ctrl1 <- -0.5
Alpha.ctrl2 <- Alpha.ctrl1 * perturb_factor
beta.case1 <- c(0.0, 0.75, -0.50)
beta.case2 <- beta.case1 * perturb_factor
beta.ctrl1 <- c(2.75, 1.00, 0.50)
beta.ctrl2 <- c(beta.ctrl1[1], beta.ctrl1[2:3] * perturb_factor)


Theta <- 6
Phi <- 12


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)


#### Simulate gaussian process
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)
hist(W)


#### Simulate locations
r <- caPr.disc[[1]]
locs1 <- simLocW(W, r, beta=0, seed=11) # 42
locs2 <- simLocW(W, r, beta=0, seed=42)
sum(locs1$status)
sum(locs2$status)
par(mfrow=c(1,2))
plot(r)
points(locs1$coords)
plot(r)
points(locs2$coords)


#### Simulate counts given locations
cov.disc <- caPr.disc
case.data1 <- simConditionalGp2(cov.disc, locs1, beta.case1, Alpha.case1, W, seed=42)
ctrl.data1 <- simConditionalGp2(cov.disc, locs1, beta.ctrl1, Alpha.ctrl1, W, seed=40)
print(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y))

case.data2 <- simConditionalGp2(cov.disc, locs2, beta.case2, Alpha.case2, W, seed=42)
ctrl.data2 <- simConditionalGp2(cov.disc, locs2, beta.ctrl2, Alpha.ctrl2, W, seed=40)
print(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y))


new_row <- list(
  beta.case1=paste(as.character(beta.case1), collapse=' '),
  beta.case2=paste(as.character(beta.case2), collapse=' '),
  beta.ctrl1=paste(as.character(beta.ctrl1), collapse=' '),
  beta.ctrl2=paste(as.character(beta.ctrl2), collapse=' '),
  alpha.case1=Alpha.case1,
  alpha.case2=Alpha.case2,
  alpha.ctrl1=Alpha.ctrl1,
  alpha.ctrl2=Alpha.ctrl2,
  Theta=Theta,
  Phi=Phi,
  perturb=perturb,
  w1=W1,
  w2=W2,
  ws=ws
)

write.table(data.frame(new_row), 
            file='/Users/brianconroy/Documents/research/dataInt/output/simMulti_params_opt1_comparison.txt', 
            row.names=F, 
            append=T)
