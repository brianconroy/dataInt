library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


#########################
# Integrates areal counts
# w/ preferential samples
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
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
W <- simRegion$W
caWc.disc <- simRegion$raster
nNeighbors <- simRegion$nNeighbors


#### Set parameters
beta.case <- c(1, 2)
beta.control <- c(3, -1)
beta.samp <- -1
alpha.case <- 4
alpha.control <- -6
rho.1 <- 0.75
tau.sq.1 <- 0.5


#### Spatial covariate surfaces
cov.disc <- caWc.disc[[c(12)]]
loc.disc <- caWc.disc[[c(6)]]
values(loc.disc) <- values(loc.disc) + abs(min(values(loc.disc), na.rm=T))


#### Simulate survey locations
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
sum(locs$status)


prob.rast <- loc.disc
prob.rast[!is.na(values(loc.disc))] <- locs$probs
par(mfrow=c(1,3))
plot(cov.disc, main='A)')
plot(loc.disc, main='B)')
plot(prob.rast, main='C)')
points(locs$coords, pch=16)


par(mfrow=c(2,2))
cuts=c(-5, -2, 0, 2, 4, 6, 10)
pal1 <- colorRampPalette(c("blue","red"))
pal2 <- colorRampPalette(c("purple","yellow"))
intens.case.true <- calcIntens(cov.disc, beta.case)
plot(cov.disc, main='A)')
plot(loc.disc, main='B')
plot(intens.case.true, breaks=cuts, col=pal1(8), main='C)')
points(locs$coords, pch=16)
plot(prob.rast, col=pal2(8), main='D)')
points(locs$coords, pch=16)



#### Simulate counts given locations
count.case <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha.case, global.center=T, seed=505)
count.control <- simConditional(cov.disc, locs, beta.control, beta.samp, alpha.control, global.center=T, seed=505)
data.ps <- list(conditional.case=count.case, conditional.control=count.control, loc=locs)


par(mfrow=c(2,2))
plot(prob.rast, main='A)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.case$y/30)

plot(cov.disc, main='B)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.case$y/30)

plot(prob.rast, main='C)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.control$y/10)

plot(cov.disc, main='D)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.control$y/10)

par(mfrow=c(2,1))
hist(count.case$y, mar=c(4, 4, 3, 2) + 0.1)
hist(count.control$y, mar=c(4, 4, 3, 2) + 0.1)


#### Simulate areal data
data.1 <- simLeroux(cov.disc, beta.case, tau.sq.1, rho.1, nNeighbors, seed=105)
data.spat <- list(data.1)


hist(data.1$y)
mean(data.1$y)


data <- list(ps=data.ps, spat=data.spat)


###################
## Integrated Model
###################


output <- prefSampleIntegration(
  data, n.sample=85000, burnin=20000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.ctr=0.01,
  proposal.sd.beta.l=0.05, proposal.sd.alpha=0.05, 
  proposal.sd.phi=0.03, proposal.sd.rho=0.05, 
  self.tune=TRUE, fix.rho=FALSE
)
print(output$accept)


truevals <- list(
  beta1.loc=beta.samp,
  beta0.case=beta.case[1],
  beta1.case=beta.case[2],
  beta0.control=beta.control[1],
  beta1.control=beta.control[2],
  alpha.case=alpha.case,
  alpha.control=alpha.control
)
par(mfrow=c(2,4))
viewOutput(output, type='preferential.cc.diff', truevals, view.hist=F)
viewOutput(output, type='preferential.cc.diff', truevals, view.trace=F)


################
## Multiple Runs
################


n.runs <- 4
par(mfrow=c(2,4))
results <- list()
for (i in 1:n.runs){
  
  print(paste("iteration", i))
  output <- prefSampleIntegration(
    data, n.sample=85000, burnin=20000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.ctr=0.01,
    proposal.sd.beta.l=0.05, proposal.sd.alpha=0.05, 
    proposal.sd.phi=0.03, proposal.sd.rho=0.05, 
    self.tune=TRUE, fix.rho=FALSE
  )
  results[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
  
}


pool.results <- function(results){
  
  results.new <- list()
  new.beta.l <- results[[1]]$samples.beta.l
  new.alpha.case <- results[[1]]$samples.alpha.case
  new.alpha.control <- results[[1]]$samples.alpha.control
  new.beta.case <- results[[1]]$samples.beta.case
  new.beta.control <- results[[1]]$samples.beta.control
  new.accept <- results[[1]]$accept
  new.tau2 <- results[[1]]$samples.tau2

  D <- dim(results[[1]]$samples.rho)[1]
  n.keep <- dim(results[[1]]$samples.rho)[2]
  k <- dim(results[[1]]$samples.phi)[3]
  new.rho <- array(NA, c(D, length(results)*n.keep, 1))
  new.phi <- array(NA, c(D, length(results)*n.keep, k))
  for (d in 1:D){
    new.rho[d,1:n.keep,] <- results[[1]]$samples.rho[d,,]
    new.phi[d,1:n.keep,] <- results[[1]]$samples.phi[d,,]
  }


  for (i in 2:length(results)){
    new.beta.l <- rbind(new.beta.l, results[[i]]$samples.beta.l)  
    new.alpha.case <- rbind(new.alpha.case, results[[i]]$samples.alpha.case)
    new.alpha.control <- rbind(new.alpha.control, results[[i]]$samples.alpha.control)
    new.beta.case <- rbind(new.beta.case, results[[i]]$samples.beta.case)
    new.beta.control <- rbind(new.beta.control, results[[i]]$samples.beta.control)
    new.accept <- new.accept + results[[i]]$accept
    new.tau2 <- rbind(new.tau2, results[[i]]$samples.tau2)
    for (d in 1:D){
      new.rho[d,(n.keep*(i-1) + 1):(n.keep*i),] <- results[[i]]$samples.rho[d,,]
      new.phi[d,(n.keep*(i-1) + 1):(n.keep*i),] <- results[[i]]$samples.phi[d,,]
    }
  }
  
  new.accept <- new.accept/length(results)
  results.new$samples.beta.l <- new.beta.l
  results.new$samples.beta.case <- new.beta.case
  results.new$samples.beta.control <- new.beta.control
  results.new$samples.alpha.case <- new.alpha.case
  results.new$samples.alpha.control <- new.alpha.control
  results.new$accept <- new.accept
  results.new$samples.tau2 <- new.tau2
  results.new$samples.rho <- new.rho
  results.new$samples.phi <- new.phi
  
  return(results.new)
  
}


results.trimmed <- lapply(results, function(x){ trimSamps(x, 10000)})


pooled <- pool.results(results)
par(mfrow=c(2,4))
viewOutput(pooled, type='preferential.cc.diff', truevals, view.hist=F)
viewOutput(pooled, type='preferential.cc.diff', truevals, view.trace=F)
par(mfrow=c(1,2))
plot(pooled$samples.rho, type='l', main='rho')
abline(h=rho.1, col='2')
plot(pooled$samples.tau2, type='l', main='tau2')
abline(h=tau.sq.1, col='2')


mean.phi <- colMeans(pooled$samples.phi[1,,])
hist(mean.phi)
df.phi <- data.frame(mean.phi)
ggplot(data=df.phi, aes(df.phi$mean.phi)) + 
  geom_histogram(breaks=seq(-1.5, 1.5, by =0.125), 
                 col="blue", 
                 fill="green", 
                 alpha=0.2) + 
  labs(x="phi")
round(summary(mean.phi), 3)
phi.disc <- cov.disc
phi.disc[!is.na(values(phi.disc))] <- mean.phi
plot(phi.disc)


cuts=seq(-0.8, 1.25, by=0.5)
pal1 <- colorRampPalette(c("blue","green"))
plot(phi.disc, col=pal1(16))


prefSummary <- summarize(pooled, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PS CAR'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


### comparison
# integrated model summary
# integrated only
# pS only
# usual comps
# + risk surfaces - these should be best for the integrated model


#########
# pS only
#########


n.runs <- 8
par(mfrow=c(2,4))
results.ps <- list()
for (i in 1:n.runs){
  
  print(paste("iteration", i))
  output <- prefSampleDiffCC(
    data$ps, n.sample=150000, burnin=15000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  results.ps[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
  
}


pooled.ps <- pool.results.psccd(results.ps)
par(mfrow=c(2,4))
viewOutput(pooled.ps, type='preferential.cc.diff', truevals, view.hist=F)


prefSummary.only <- summarize(pooled.ps, truevals, dic=FALSE)
df.ps.only <- ldply(prefSummary.only, data.frame)
df.ps.only$model <- 'PS only'
df.ps.only <- df.ps.only[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps.only) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


##########
# CAR only
##########


n.runs <- 4
par(mfrow=c(2,2))
results.car <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- carLerouxIntegration(
    data$spat, 
    proposal.sd.beta=0.025, 
    proposal.sd.phi=1,
    proposal.sd.rho=0.05,
    prior.var.beta=NULL, 
    prior.tau2=NULL,
    self.tune=TRUE,
    fix.rho=FALSE,
    n.sample=90000)
  results.car[[i]] <- output
  plot(output$samples.beta[,2], type='l', main='beta1')
  abline(h=beta.case[2], col='2')
}


pooled.car <- pool.results.car(results.car)
params.car <- list(beta0=beta.case[1], beta1=beta.case[2])
carSummary <- summarize(pooled.car, params.car, dic=FALSE)
df.car <- ldply(carSummary, data.frame)
df.car$model <- 'CAR only'
df.car$parameter <- c('beta0.case', 'beta1.case')
df.car <- df.car[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.car) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


#############################
# Compare Parameter Estimates
#############################


df.comp <- rbind(df.ps, df.ps.only, df.car)
df.comp <- df.comp[with(df.comp, order(parameter)),]
df.comp <- df.comp[c('parameter', 'model', 'estimate', 'sd', 'pbias')]
write.table(df.comp, file='Documents/research/dataInt/output/pS_summary_int.csv', sep=',', row.names=F)


#########################
# Estimated Risk Surfaces
#########################


cuts=c(-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16)
pal <- colorRampPalette(c("blue","red"))


## True surface
intens.case.true <- calcIntens(cov.disc, beta.case)
intens.control.true <- calcIntens(cov.disc, beta.control)
diff.true <- rastDiff(intens.case.true, intens.control.true)
plot(diff.true, breaks=cuts, col=pal(8))


## PS Integrated estimate
beta.case.psi <- colMeans(pooled$samples.beta.case)
beta.control.psi <- colMeans(pooled$samples.beta.control)
intens.case.psi <- calcIntens(cov.disc, beta.case.psi)
intens.control.psi <- calcIntens(cov.disc, beta.control.psi)
diff.psi <- rastDiff(intens.case.psi, intens.control.psi)
plot(diff.psi, breaks=cuts, col=pal(8))


## PS Only estimate
beta.case.ps <- colMeans(pooled.ps$samples.beta.case)
beta.control.ps <- colMeans(pooled.ps$samples.beta.control)
intens.case.ps <- calcIntens(cov.disc, beta.case.ps)
intens.control.ps <- calcIntens(cov.disc, beta.control.ps)
diff.ps <- rastDiff(intens.case.ps, intens.control.ps)
plot(diff.ps, breaks=cuts, col=pal(8))


## CAR Only estimate
beta.case.car <- colMeans(pooled.car$samples.beta)
intens.case.car <- calcIntens(cov.disc, beta.case.car)
plot(intens.case.car, breaks=cuts, col=pal(8))


## Combined figure
par(mfrow=c(2,2))
plot(diff.true, breaks=cuts, col=pal(11), main='A)')
plot(diff.psi, breaks=cuts, col=pal(11), main='B)')
plot(diff.ps, breaks=cuts, col=pal(11), main='C)')
plot(intens.case.car, breaks=cuts, col=pal(11), main='D)')


summs <- list()
summs[[1]] <- summarizeSurf(diff.true, diff.psi, 'PS')
summs[[2]] <- summarizeSurf(diff.true, diff.ps, 'PS only')
summs[[3]] <- summarizeSurf(diff.true, intens.case.car, 'CAR only')
df.risks <- ldply(summs, data.frame)
write.table(df.risks, file='Documents/research/dataInt/output/pS_int_mse.csv', sep=',', row.names=F)
