library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
sourceDirectory('Documents/research/dataInt/R/')


##############################
## survey locations determined
## by Gaussian process
##############################


##########################
## test the joint w update
##########################


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


hist(w)
w.disc <- caWc.disc[[1]]
w.disc[][!is.na(w.disc[])] <- w
plot(w.disc)


#### Simulate locations
beta.samp <- c(-2, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, w=w, seed=56)
sum(locs$status)
hist(locs$w)
plot(loc.disc)
points(locs$coords)


#### Chosen alpha
alpha <- 3


#### Disease covariate surface
cov.disc <- caWc.disc[[c(1)]]
beta.case <- c(1, 2)


#### Simulate counts given locations
count.data <- simConditionalGp(cov.disc, locs, beta.case, beta.samp, alpha, w=w, seed=43)
count.data$y
hist(count.data$w)
summary(count.data$w)
count.data$w - w[as.logical(locs$status)]


#### Fit
x.l <- locs$x.scaled
x.l.sub <- x.l[as.logical(locs$status),]
y.l <- locs$status

x.c <- count.data$x.standardised
y.c <- count.data$y

n.sample <- 25000

beta.c <- beta.case
beta.l <- beta.samp
w.i <- w

proposal.sd.beta.l <- 0.2
proposal.sd.beta.c <- 0.01
proposal.sd.w <- 0.5

p <- ncol(x.c)
prior.mean.beta <- rep(0, p)
prior.var.beta <- rep(1000, p)

samples.beta.l <- array(NA, c(n.sample, length(beta.samp)))
samples.beta.c <- array(NA, c(n.sample, length(beta.case)))
samples.w <- array(NA, c(n.sample, length(w)))

accept <- c(0, 0, 0)

for(i in 1:n.sample){
  
  ## Sample from beta (locational)
  beta.out.l <- betaLogisticUpdate(x.l, y.l, beta.l, proposal.sd.beta.l, offset=w.i)
  beta.l <- beta.out.l$beta
  
  
  ## Sample from beta (case)
  #loc.pred <- expit(x.l.sub %*% beta.samp)
  #offset.beta.c <- alpha * loc.pred + count.data$w
  #beta.out.c <- betaPoissonUpdate(x.c, y.c, beta.c, offset.beta.c,
  #                                prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
  beta.c <- beta.case #beta.out.c$beta
  
  
  ## Sample from w
  #mu.i <- alpha * loc.pred + x.c %*% beta.c
  #sigma.i <- count.data$sigma
  #w.out.i <- wPoissonUpdate(y.c, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w #w.out.i$w
  
  
  samples.beta.l[i,] <- beta.l
  samples.beta.c[i,] <- beta.c
  samples.w[i,] <- w.i
  
  accept[1] <- accept[1] + beta.out.l$accept
  #accept[2] <- accept[2] + beta.out.c$accept
  #accept[3] <- accept[3] + w.out.i$accept
  
}

accept <- accept/n.sample
accept[3] <- accept[3]/length(count.data$w)
print(accept)

plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col='2')
plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col='2')

plot(samples.beta.c[,1], type='l'); abline(h=beta.case[1], col='2')
plot(samples.beta.c[,2], type='l'); abline(h=beta.case[2], col='2')

what <- colMeans(samples.w)
plot(count.data$w, what)
abline(0, 1, col=2)

par(mfrow=c(4, 4))
j <- 2
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=count.data$w[j*i], col='2')
}


#########################
## test the case count
## model param estimation
#########################


#######
# Setup
#######


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=2)
W <- simRegion$W
caWc.disc <- simRegion$raster


#### Simulate survey locations
Theta <- 5
Phi <- 3
beta.samp <- c(-1.5, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=56)
sum(locs$status)


#### Disease covariate surface
cov.disc <- caWc.disc[[c(1)]]
beta.case <- c(1, 2)


#### Visualize surfaces
pal1 <- colorRampPalette(c("blue","red"))
pal2 <- colorRampPalette(c("blue","green"))


# raw covariates
par(mfrow=c(1,3))
plot(cov.disc, main='A)', col=pal1(8))
plot(loc.disc, main='B)', col=pal2(8))
plot(w.disc, main='C)', col=pal2(8))


# sampling intensity
par(mfrow=c(1,1))
prob.rast <- loc.disc
prob.rast[!is.na(values(loc.disc))] <- locs$probs
plot(prob.rast, col=pal2(8))
points(locs$coords, pch=16)


# chosen alpha
alpha <- 3



d <- as.matrix(dist(locs$coords, diag=TRUE, upper=TRUE))
sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
w <- mvrnorm(n=1, mu=rep(0, ncol(d)), sigma)


# Simulate counts given locations
count.data <- simConditionalGp(cov.disc, locs, beta.case, beta.samp, alpha, Theta, Phi, seed=43)
# count.data <- simConditionalGp(cov.disc, locs, beta.case, beta.samp, alpha, w=w, seed=43)
count.data$y
plot(as.im(count.data$sigma))
hist(count.data$w)
summary(count.data$w)


x.l <- locs$x.scaled
x.l.sub <- x.l[as.logical(locs$status),]
y.l <- locs$status

x.c <- count.data$x.standardised
y.c <- count.data$y

n.sample <- 500000

beta.c <- beta.case
beta.l <- beta.samp
w.i <- count.data$w

proposal.sd.beta.l <- 0.2
proposal.sd.beta.c <- 0.05
proposal.sd.w <- 0.5

p <- ncol(x.c)
prior.mean.beta <- rep(0, p)
prior.var.beta <- rep(1000, p)

samples.beta.l <- array(NA, c(n.sample, length(beta.samp)))
samples.beta.c <- array(NA, c(n.sample, length(beta.case)))
samples.w <- array(NA, c(n.sample, length(count.data$w)))

accept <- c(0, 0, 0)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  
  ## Sample from beta (locational)
  beta.out.l <- betaLogisticUpdate(x.l, y.l, beta.l, proposal.sd.beta.l)
  beta.l <- beta.out.l$beta
  
  
  ## Sample from beta (case)
  loc.pred <- expit(x.l.sub %*% beta.l)
  offset.beta.c <- alpha * loc.pred + w.i
  beta.out.c <- betaPoissonUpdate(x.c, y.c, beta.c, offset.beta.c,
                                  prior.mean.beta, prior.var.beta, proposal.sd.beta.c)
  beta.c <- beta.out.c$beta
  
  
  ## Sample from w
  mu.i <- alpha * loc.pred + x.c %*% beta.c
  sigma.i <- count.data$sigma
  w.out.i <- wPoissonUpdate(y.c, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w.out.i$w
  
  
  samples.beta.l[i,] <- beta.l
  samples.beta.c[i,] <- beta.c
  samples.w[i,] <- w.i
  
  accept[1] <- accept[1] + beta.out.l$accept
  accept[2] <- accept[2] + beta.out.c$accept
  accept[3] <- accept[3] + w.out.i$accept
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

accept <- accept/n.sample
accept[3] <- accept[3]/length(count.data$w)
print(accept)

plot(samples.beta.l[,1], type='l'); abline(h=beta.samp[1], col='2')
plot(samples.beta.l[,2], type='l'); abline(h=beta.samp[2], col='2')

plot(samples.beta.c[,1], type='l'); abline(h=beta.case[1], col='2')
plot(samples.beta.c[,2], type='l'); abline(h=beta.case[2], col='2')

what <- colMeans(samples.w)
plot(count.data$w, what)
abline(0, 1, col=2)

par(mfrow=c(4, 4))
j <- 2
for(i in 1:16){
  plot(samples.w[,j*i], type='l'); abline(h=count.data$w[j*i], col='2')
}


#########################
## test the locational
## model param estimation
## wLogisticUpdate
#########################


#### Worldclim data
wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=2)
W <- simRegion$W
caWc.disc <- simRegion$raster


#### Simulate survey locations
Theta <- 5
Phi <- 3
beta.samp <- c(-0.75, 2)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLocGp(loc.disc, beta.samp, Theta, Phi, seed=56)
sum(locs$status)

hist(locs$w)
w.disc <- loc.disc
w.disc[][!is.na(w.disc[])] <- locs$w


#### Fit


x.l <- locs$x.scaled
y.l <- locs$status
n.sample <- 25000

samples.beta <- array(NA, c(n.sample, length(beta.samp)))
samples.w <- array(NA, c(n.sample, length(locs$w)))
samples.theta <- array(NA, c(n.sample, 1))
samples.phi <- array(NA, c(n.sample, 1))

beta.l <- beta.samp
w.i <- locs$w
theta.i <- Theta
phi.i <- Phi

proposal.sd.beta.l <- 0.1
proposal.sd.theta <- 0.1
proposal.sd.w <- 0.5
d <- locs$d
N <- nrow(d)

prior.theta <- c(2.2, 2.2)
prior.phi <- c(3, 6)

accept <- c(0, 0, 0)

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

for(i in 1:n.sample){
  
  ## Sample from beta (locational)
  beta.out.l <- betaLogisticUpdate(x.l, y.l, beta.l, proposal.sd.beta.l, offset=w.i)
  beta.l <- beta.out.l$beta
  accept[1] <- accept[1] + beta.out.l$accept
  
  ## sample from w
  mu.i <- x.l %*% beta.l
  sigma.i <- Exponential(d, range=theta.i, phi=phi.i)
  w.out.i <- wLogisticUpdate(y.l, w.i, mu.i, sigma.i, proposal.sd.w)
  w.i <- w.out.i$w
  accept[2] <- accept[2] + w.out.i$accept
  
  ## sample from theta
  # theta.out <- rangeMhUpdate(theta.i, w.i, d, phi.i, proposal.sd.theta, a=prior.theta[1], b=prior.theta[2])
  # theta.i <- theta.out$theta
  # accept[3] <- accept[3] + theta.out$accept
  # 
  # ## sample from phi
  # R.i <- sigma.i/phi.i
  # phi.i <- 1 / rgamma(1, N/2 + prior.phi[1], t(w.i) %*% solve(R.i) %*% w.i/2 + prior.phi[2])
  # samples.phi[i,] <- phi.i
  
  samples.beta[i,] <- beta.l
  samples.w[i,] <- w.i
  #samples.theta[i,] <- theta.i
  
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.sample)
  }
  
}

print(accept[2]/(n.sample * length(locs$status)))

plot(samples.beta[,1], type='l')
abline(h=beta.samp[1], col=2)

plot(samples.beta[,2], type='l')
abline(h=beta.samp[2], col=2)

plot(samples.theta, type='l')
abline(h=Theta, col=2)

what <- colMeans(samples.w)
plot(x=locs$w, y=what, ylab='estimated w', xlab='true w')
abline(0, 1, col=2)

plot(x=locs$w, y=(locs$w - what)/locs$w)

par(mfrow=c(4,4))
j = 1
for (i in 237:253){
  plot(samples.w[,j*i], type='l')
  abline(h=locs$w[j*i], col=2)
}

beta.hat <- colMeans(samples.beta)
probs <- expit(x.l %*% beta.hat + what)
preds <- as.numeric(probs > 0.5)

sum(abs(preds - locs$status))




#################################
# compare occupancy predictions
# against a simple logistic model
# w/o gaussian process
#################################


mod <- glm(y.l ~ x.l, family='binomial')
predict(mod)






#### chosen alpha
alpha <- 3


#### Simulate counts given locations
count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha, seed=42)
data <- list(conditional=count.data, loc=locs)


par(mfrow=c(1,1))
plot(loc.disc)
points(locs$coords, cex=count.data$y/40)


plot(cov.disc)
points(locs$coords, cex=count.data$y/35)


hist(count.data$y)
mean(count.data$y)


##############################
## Preferential Sampling Model
##############################


output <- prefSampleNew(
  data, n.sample=500000, burnin=20000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
)


truevals <- list(
  beta0.loc=beta.samp[1],
  beta1.loc=beta.samp[2],
  beta0=beta.case[1],
  beta1=beta.case[2],
  alpha=alpha
)
print(output$accept)
par(mfrow=c(2,3))
viewOutput(output, type='preferential.alpha', truevals, view.hist=FALSE)


prefSummary <- summarize(output, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PS'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


###################################
# Compare against model which
# updates locational parameters
# only from the location likelihood
###################################


output.loc <- prefSample(
  data, n.sample=500000, burnin=20000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
)


print(output.loc$accept)
par(mfrow=c(2,2))
viewOutput(output.loc, type='preferential.alpha', truevals, view.hist=FALSE)


prefSummaryLoc <- summarize(output.loc, truevals, dic=FALSE)
df.ps.loc <- ldply(prefSummaryLoc, data.frame)
df.ps.loc$model <- 'PS loc'
df.ps.loc <- df.ps.loc[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps.loc) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


df.comp <- rbind(df.ps, df.ps.loc)
df.comp <- df.comp[with(df.comp, order(parameter)),]
write.table(df.comp, file='Documents/research/dataInt/pS_summary_comp1.csv', sep=',', row.names=F)


######################
## Iterate models over 
## various alpha
######################


summaries <- data.frame()
for (a in c(1, 2, 4, 6, 8)){
  
  print(a)
  count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, a)
  data <- list(conditional=count.data, loc=locs)
  
  output.a <- prefSampleNew(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  output.a.loc <- prefSample(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  truevals <- list(
    beta0.loc=beta.samp[1],
    beta1.loc=beta.samp[2],
    beta0=beta.case[1],
    beta1=beta.case[2],
    alpha=a
  )
  
  viewOutput(output.a, type='preferential.alpha', truevals)
  
  prefSummary <- summarize(output.a, truevals, dic=FALSE)
  df.ps <- ldply(prefSummary, data.frame)
  df.ps$alpha <- a
  df.ps$model <- 'PS'
  summaries <- rbind(summaries, df.ps)
  
  prefSummaryLoc <- summarize(output.a.loc, truevals, dic=FALSE)
  df.ps.loc <- ldply(prefSummaryLoc, data.frame)
  df.ps.loc$alpha <- a
  df.ps.loc$model <- 'PS loc'
  summaries <- rbind(summaries, df.ps.loc)
  
}


par(mfrow=c(2,3))
for (p in unique(summaries$parameter)){
  
  summaries.p <- summaries[summaries$parameter == p,]
  summaries.p.joint <- summaries.p[summaries.p$model == 'PS',]
  summaries.p.loc <- summaries.p[summaries.p$model == 'PS loc',]
  
  uy <- max(1.1 * max(summaries.p.joint$percbias),
            1.1 * max(summaries.p.loc$percbias),
            5)
  ly <- min(0,
            1.1 * min(summaries.p.joint$percbias), 
            1.1 * min(summaries.p.loc$percbias))
  
  plot(x=summaries.p.joint$alpha, y=summaries.p.joint$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.joint$alpha, y=summaries.p.joint$percbias, type='l')
  points(x=summaries.p.loc$alpha, y=summaries.p.loc$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.loc$alpha, y=summaries.p.loc$percbias, type='l', col='2')
  abline(h=0, lty=2)
  
}


#################################
# Now consider a shared covariate
#################################


#### Simulate survey locations
beta.samp <- c(-3, 1, 1)
loc.disc <- caWc.disc[[c(6, 15)]]
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
sum(locs$status)


#### Disease covariate surface
beta.case <- c(1, 1, 2)
cov.disc <- caWc.disc[[c(12, 15)]]


#### chosen alpha
alpha <- 2


#### Simulate counts given locations
count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha, seed=42)
data <- list(conditional=count.data, loc=locs)


output <- prefSampleNew(
  data, n.sample=500000, burnin=20000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
)


truevals <- list(
  beta0.loc=beta.samp[1],
  beta1.loc=beta.samp[2],
  beta2.loc=beta.samp[3],
  beta0=beta.case[1],
  beta1=beta.case[2],
  beta2=beta.case[3],
  alpha=alpha
)
print(output$accept)
par(mfrow=c(2,4))
viewOutput(output, type='preferential.alpha', truevals, view.hist=FALSE)


prefSummary <- summarize(output, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PS'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


#### Compare against original
output.loc <- prefSample(
  data, n.sample=500000, burnin=20000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
)


prefSummaryLoc <- summarize(output.loc, truevals, dic=FALSE)
df.ps.loc <- ldply(prefSummaryLoc, data.frame)
df.ps.loc$model <- 'PS loc'
df.ps.loc <- df.ps.loc[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps.loc) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


df.comp <- rbind(df.ps, df.ps.loc)
df.comp <- df.comp[with(df.comp, order(parameter)),]


######################
## Iterate models over 
## various alpha
######################


summaries <- data.frame()
for (a in c(1, 2, 4, 6, 8)){
  
  print(a)
  count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, a)
  data <- list(conditional=count.data, loc=locs)
  
  output.a <- prefSampleNew(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  output.a.loc <- prefSample(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  truevals <- list(
    beta0.loc=beta.samp[1],
    beta1.loc=beta.samp[2],
    beta0=beta.case[1],
    beta1=beta.case[2],
    alpha=a
  )
  
  viewOutput(output.a, type='preferential.alpha', truevals)
  
  prefSummary <- summarize(output.a, truevals, dic=FALSE)
  df.ps <- ldply(prefSummary, data.frame)
  df.ps$alpha <- a
  df.ps$model <- 'PS'
  summaries <- rbind(summaries, df.ps)
  
  prefSummaryLoc <- summarize(output.a.loc, truevals, dic=FALSE)
  df.ps.loc <- ldply(prefSummaryLoc, data.frame)
  df.ps.loc$alpha <- a
  df.ps.loc$model <- 'PS loc'
  summaries <- rbind(summaries, df.ps.loc)
  
}


par(mfrow=c(2,3))
for (p in unique(summaries$parameter)){
  
  summaries.p <- summaries[summaries$parameter == p,]
  summaries.p.joint <- summaries.p[summaries.p$model == 'PS',]
  summaries.p.loc <- summaries.p[summaries.p$model == 'PS loc',]
  
  uy <- max(1.1 * max(summaries.p.joint$percbias),
            1.1 * max(summaries.p.loc$percbias),
            5)
  ly <- min(0,
            1.1 * min(summaries.p.joint$percbias), 
            1.1 * min(summaries.p.loc$percbias))
  
  plot(x=summaries.p.joint$alpha, y=summaries.p.joint$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.joint$alpha, y=summaries.p.joint$percbias, type='l')
  points(x=summaries.p.loc$alpha, y=summaries.p.loc$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.loc$alpha, y=summaries.p.loc$percbias, type='l', col='2')
  abline(h=0, lty=2)
  
}


########################
## Iterate models over 
## locational intercepts
## (i.e. number of surveys)
########################


a <- 2
summaries <- data.frame()
for (beta.loc0 in c(-4, -3, -2, -1, 1)){
  
  print(beta.loc0)
  beta.samp <- c(beta.loc0, 1)
  locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
  nsamp <- sum(locs$status)
  
  count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, a)
  data <- list(conditional=count.data, loc=locs)
  
  output.a <- prefSampleNew(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  output.a.loc <- prefSample(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  truevals <- list(
    beta0.loc=beta.samp[1],
    beta1.loc=beta.samp[2],
    beta0=beta.case[1],
    beta1=beta.case[2],
    alpha=a
  )
  
  prefSummary <- summarize(output.a, truevals, dic=FALSE)
  df.ps <- ldply(prefSummary, data.frame)
  df.ps$n <- nsamp
  df.ps$model <- 'PS'
  summaries <- rbind(summaries, df.ps)
  
  prefSummaryLoc <- summarize(output.a.loc, truevals, dic=FALSE)
  df.ps.loc <- ldply(prefSummaryLoc, data.frame)
  df.ps.loc$n <- nsamp
  df.ps.loc$model <- 'PS loc'
  summaries <- rbind(summaries, df.ps.loc)
  
}


par(mfrow=c(2,3))
for (p in unique(summaries$parameter)){
  
  summaries.p <- summaries[summaries$parameter == p,]
  summaries.p.joint <- summaries.p[summaries.p$model == 'PS',]
  summaries.p.loc <- summaries.p[summaries.p$model == 'PS loc',]
  
  uy <- max(1.1 * max(summaries.p.joint$percbias),
            1.1 * max(summaries.p.loc$percbias),
            5)
  ly <- min(0,
            1.1 * min(summaries.p.joint$percbias), 
            1.1 * min(summaries.p.loc$percbias))
  
  plot(x=summaries.p.joint$n, y=summaries.p.joint$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.joint$n, y=summaries.p.joint$percbias, type='l')
  points(x=summaries.p.loc$n, y=summaries.p.loc$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.loc$n, y=summaries.p.loc$percbias, type='l', col='2')
  abline(h=0, lty=2)
  
}


########################
## and the same for the
## shared covariate 
## scenario
########################


#### Simulate survey locations
beta.samp <- c(-3, 1, 1)
loc.disc <- caWc.disc[[c(6, 15)]]


#### Disease covariate surface
beta.case <- c(1, 1, 2)
cov.disc <- caWc.disc[[c(12, 15)]]


a <- 2
summaries <- data.frame()
for (beta.loc0 in c(-4, -3, -2, -1, 1)){
  
  print(beta.loc0)
  beta.samp <- c(beta.loc0, 1, 1)
  locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
  nsamp <- sum(locs$status)
  
  count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, a)
  data <- list(conditional=count.data, loc=locs)
  
  output.a <- prefSampleNew(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  output.a.loc <- prefSample(
    data, n.sample=300000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  truevals <- list(
    beta0.loc=beta.samp[1],
    beta1.loc=beta.samp[2],
    beta0=beta.case[1],
    beta1=beta.case[2],
    beta2=beta.case[3],
    alpha=a
  )
  
  prefSummary <- summarize(output.a, truevals, dic=FALSE)
  df.ps <- ldply(prefSummary, data.frame)
  df.ps$n <- nsamp
  df.ps$model <- 'PS'
  summaries <- rbind(summaries, df.ps)
  
  prefSummaryLoc <- summarize(output.a.loc, truevals, dic=FALSE)
  df.ps.loc <- ldply(prefSummaryLoc, data.frame)
  df.ps.loc$n <- nsamp
  df.ps.loc$model <- 'PS loc'
  summaries <- rbind(summaries, df.ps.loc)
  
}


par(mfrow=c(2,3))
for (p in unique(summaries$parameter)){
  
  summaries.p <- summaries[summaries$parameter == p,]
  summaries.p.joint <- summaries.p[summaries.p$model == 'PS',]
  summaries.p.loc <- summaries.p[summaries.p$model == 'PS loc',]
  
  uy <- max(1.1 * max(summaries.p.joint$percbias),
            1.1 * max(summaries.p.loc$percbias),
            5)
  ly <- min(0,
            1.1 * min(summaries.p.joint$percbias), 
            1.1 * min(summaries.p.loc$percbias))
  
  plot(x=summaries.p.joint$n, y=summaries.p.joint$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.joint$n, y=summaries.p.joint$percbias, type='l')
  points(x=summaries.p.loc$n, y=summaries.p.loc$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p.loc$n, y=summaries.p.loc$percbias, type='l', col='2')
  abline(h=0, lty=2)
  
}
