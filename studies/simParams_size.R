#############################
# Choose simulation parameter
# values

# fix prevalence ~ 15-20%
# fix alphas = 0.5, -0.5
# record total case + control
#############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate gaussian process
Theta <- 6
Phi <- 12
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)


#### Simulate locations
r <- caPr.disc[[1]]
locs <- simLocW(W, r, beta=0, seed=11) # 42
sum(locs$status)
hist(W)
plot(r)
points(locs$coords)


Alpha.case <- 0.5
beta.case <- c(0, 0.75, -0.5)
Alpha.ctrl <- -0.5
beta.ctrl <- c(2.75, 1, 0.5)


#### Simulate counts given locations
cov.disc <- caPr.disc
case.data <- simConditionalGp2(cov.disc, locs, beta.case, Alpha.case, W, seed=42)
ctrl.data <- simConditionalGp2(cov.disc, locs, beta.ctrl, Alpha.ctrl, W, seed=40)


sum(case.data$y)/sum(case.data$y + ctrl.data$y)
sum(case.data$y)
sum(ctrl.data$y)
sum(case.data$y + ctrl.data$y)
case.data$y
ctrl.data$y

size <- '5495'

new_row <- list(
  size=size,
  beta.case=paste(as.character(beta.case), collapse=' '),
  beta.ctrl=paste(as.character(beta.ctrl), collapse=' '),
  alpha.case=Alpha.case,
  alpha.ctrl=Alpha.ctrl,
  total.y.ca=sum(case.data$y),
  total.y.co=sum(ctrl.data$y),
  prev=round(sum(case.data$y)/sum(case.data$y + ctrl.data$y), 2)
)


write.table(data.frame(new_row), 
            file='/Users/brianconroy/Documents/research/dataInt/output/simParams_size.txt', 
            row.names=F, 
            append=T)
