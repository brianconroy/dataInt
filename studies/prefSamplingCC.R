library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


####################################
# Fits the (unconstrained) preferential 
# sampling Bayesian model, with case 
# and control counts
####################################


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
beta.samp <- -1
loc.disc <- caWc.disc[[c(6)]]
values(loc.disc) <- values(loc.disc) + abs(min(values(loc.disc), na.rm=T))
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=505)
sum(locs$status)
plot(loc.disc)
points(locs$coords, pch=16)


#### Disease covariate surface
cov.disc <- caWc.disc[[c(12)]]


#### Count parameter values
beta.case <- c(1, 2)
beta.control <- c(2, -1)


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


#### chosen alphas
alpha.case <- 4
alpha.control <- -6


#### Simulate counts given locations
count.case <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha.case)
count.control <- simConditional(cov.disc, locs, beta.control, beta.samp, alpha.control)
data <- list(conditional.case=count.case, conditional.control=count.control, loc=locs)


par(mfrow=c(2,1))
hist(count.case$y)
hist(count.control$y)


par(mfrow=c(2,2))
plot(prob.rast, main='A)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.case$y/15)

plot(cov.disc, main='B)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.case$y/15)

plot(prob.rast, main='C)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.control$y/5)

plot(cov.disc, main='D)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.control$y/5)


par(mfrow=c(2,1))
hist(count.case$y, mar=c(4, 4, 3, 2) + 0.1)
hist(count.control$y, mar=c(4, 4, 3, 2) + 0.1)


###########################################
## Differential Preferential Sampling Model
###########################################


output <- prefSampleDiffCC(
  data, n.sample=75000, burnin=10000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
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


n.runs <- 8
par(mfrow=c(2,4))
results <- list()
for (i in 1:n.runs){
  
  print(paste("iteration", i))
  output <- prefSampleDiffCC(
    data, n.sample=150000, burnin=15000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  results[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
  
}


pooled <- pool.results.psccd(results)
par(mfrow=c(2,4))
viewOutput(pooled, type='preferential.cc.diff', truevals, view.hist=F)
viewOutput(pooled, type='preferential.cc.diff', truevals, view.trace=F)


beta.case.psccd <- colMeans(pooled$samples.beta.case)
beta.control.psccd <- colMeans(pooled$samples.beta.control)


prefSummary <- summarize(pooled, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PSD'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


###################
# Comparison Models
###################


summarizeGlm <- function(mod, mod.name, tag, truevals){
  
  
  if (tag == 'case'){
    b0 <- truevals$beta0.case
    b1 <- truevals$beta1.case
  } else {
    b0 <- truevals$beta0.control
    b1 <- truevals$beta1.control
  }
  
  est <- round(unname(coefficients(mod)), 3)
  est.b0 <- est[1]
  est.b1 <- est[2]
  sds <- round(unname(sqrt(diag(summary(mod)$cov.scaled))), 3)
  sd.b0 <- sds[1]
  sd.b1 <- sds[2]
  
  summ <- list()
  summ.b0 <- list()
  summ.b0$model <- mod.name
  summ.b0$parameter <- paste('beta0', tag, sep='.')
  summ.b0$estimate <- est.b0
  summ.b0$sd <- sd.b0
  summ.b0$pbias <- round(100 * (est.b0 - b0)/b0, 3)
  summ[[1]] <- summ.b0
  
  summ.b1 <- list()
  summ.b1$model <- mod.name
  summ.b1$parameter <- paste('beta1', tag, sep='.')
  summ.b1$estimate <- est.b1
  summ.b1$sd <- sd.b1
  summ.b1$pbias <- round(100 * (est.b1 - b1)/b1, 3)
  summ[[2]] <- summ.b1
  
  
  return(ldply(summ, data.frame))
  
  
}


## glm fits
locs.x.sub <- data$loc$x.scaled[as.logical(data$loc$status)]
mod.1.ca <- glm(data$conditional.case$y ~ data$conditional.case$x.standardised - 1, family='poisson')
mod.2.ca <- glm(data$conditional.case$y ~ data$conditional.case$x.standardised + locs.x.sub - 1, family='poisson')
beta.case.glm1 <- unname(coefficients(mod.1.ca))
beta.case.glm2 <- unname(coefficients(mod.2.ca))[1:2]

mod.1.co <- glm(data$conditional.control$y ~ data$conditional.control$x.standardised - 1, family='poisson')
mod.2.co <- glm(data$conditional.control$y ~ data$conditional.control$x.standardised + locs.x.sub - 1, family='poisson')
beta.control.glm1 <- unname(coefficients(mod.1.co))
beta.control.glm2 <- unname(coefficients(mod.2.co))[1:2]

df.glm <- data.frame()
df.glm <- rbind(df.glm, summarizeGlm(mod.1.ca, 'GLM', 'case', truevals))
df.glm <- rbind(df.glm, summarizeGlm(mod.2.ca, 'GLM location', 'case', truevals))
df.glm <- rbind(df.glm, summarizeGlm(mod.1.co, 'GLM', 'control', truevals))
df.glm <- rbind(df.glm, summarizeGlm(mod.2.co, 'GLM location', 'control', truevals))


# legacy case control preferential sample
n.runs <- 8
par(mfrow=c(2,4))
results.legacy <- list()
for (i in 1:n.runs){
  
  print(paste("iteration", i))
  output <- prefSampleCC(
    data, n.sample=150000, burnin=15000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  results.legacy[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
  
}


pooled.legacy <- pool.results.pscc(results.legacy)
par(mfrow=c(2,3))
viewOutput(pooled.legacy, type='preferential.cc', truevals, view.hist=F)
viewOutput(pooled.legacy, type='preferential.cc', truevals, view.trace=F)


beta.case.pscc <- colMeans(pooled.legacy$samples.beta.case)
beta.control.pscc <- colMeans(pooled.legacy$samples.beta.control)


truevals.legacy <- list(
  beta1.loc=beta.samp,
  beta0.case=beta.case[1],
  beta1.case=beta.case[2],
  beta0.control=beta.control[1],
  beta1.control=beta.control[2]
)
prefSummary.legacy <- summarize(pooled.legacy, truevals.legacy, dic=FALSE)
df.ps.legacy <- ldply(prefSummary.legacy, data.frame)
df.ps.legacy$model <- 'PS'
df.ps.legacy <- df.ps.legacy[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps.legacy) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


df.comp <- rbind(df.ps, df.ps.legacy, df.glm)
df.comp <- df.comp[with(df.comp, order(parameter)),]
df.comp <- df.comp[c('parameter', 'model', 'estimate', 'sd', 'pbias')]
write.table(df.comp, file='Documents/research/dataInt/output/pS_diff_comp.csv', sep=',', row.names=F)


#######################
# Compare Risk surfaces
#######################


cuts=c(-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16)
pal <- colorRampPalette(c("blue","red"))


## True surface
intens.case.true <- calcIntens(cov.disc, beta.case)
intens.control.true <- calcIntens(cov.disc, beta.control)
diff.true <- rastDiff(intens.case.true, intens.control.true)
plot(diff.true, breaks=cuts, col=pal(11))


## PSD CC surface
intens.case.psccd <- calcIntens(cov.disc, beta.case.psccd)
intens.control.psccd <- calcIntens(cov.disc, beta.control.psccd)
diff.psccd <- rastDiff(intens.case.psccd, intens.control.psccd)
plot(diff.psccd, breaks=cuts, col=pal(11))


## PS CC surface
intens.case.pscc <- calcIntens(cov.disc, beta.case.pscc)
intens.control.pscc <- calcIntens(cov.disc, beta.control.pscc)
diff.pscc <- rastDiff(intens.case.pscc, intens.control.pscc)
plot(diff.pscc, breaks=cuts, col=pal(11))


## glm case control surface
intens.case.glm1 <- calcIntens(cov.disc, beta.case.glm1)
intens.control.glm1 <- calcIntens(cov.disc, beta.control.glm1)
diff.glmcc1 <- rastDiff(intens.case.glm1, intens.control.glm1)
plot(diff.glmcc1, breaks=cuts, col=pal(11))


## glm location case control surface
intens.case.glm2 <- calcIntens(cov.disc, beta.case.glm2)
intens.control.glm2 <- calcIntens(cov.disc, beta.control.glm2)
diff.glmcc2 <- rastDiff(intens.case.glm2, intens.control.glm2)
plot(diff.glmcc2, breaks=cuts, col=pal(11))


par(mfrow=c(2,3))
plot(diff.true, breaks=cuts, col=pal(11), main='A)')
plot(diff.psccd, breaks=cuts, col=pal(11), main='B)')
plot(diff.pscc, breaks=cuts, col=pal(11), main='C)')
plot(diff.glmcc1, breaks=cuts, col=pal(11), main='D)')
plot(diff.glmcc2, breaks=cuts, col=pal(11), main='E)')


summs <- list()
summs[[1]] <- summarizeSurf(diff.true, diff.psccd, 'PSD case control')
summs[[2]] <- summarizeSurf(diff.true, diff.pscc, 'PS case control')
summs[[3]] <- summarizeSurf(diff.true, diff.glmcc1, 'GLM case control')
summs[[4]] <- summarizeSurf(diff.true, diff.glmcc2, 'GLM location case control')
df.risks <- ldply(summs, data.frame)
write.table(df.risks, file='Documents/research/dataInt/output/pS_risk_comps.csv', sep=',', row.names=F)


##############
# Scatterplots
##############


riskScatter <- function(rast.true, rast.est, title=NULL){
  
  inds <- !is.na(values(rast.true))
  x <- values(rast.true)[inds]
  y <- values(rast.est)[inds]
  df <- data.frame(cbind(x, y))
  p <- ggplot(df, aes(x=x, y=y)) + 
    geom_point(shape=1) +
    labs(x ='true log risk') +
    labs(y='estimated log risk') + 
    ylim(-5, 10) + 
    geom_abline(slope=1, intercept=0)
  if (!is.null(title)){
    p <- p + ggtitle(title)
  }
  return(p)
}

s1 <- riskScatter(diff.true, diff.psccd, 'A)')
s2 <- riskScatter(diff.true, diff.pscc, 'B)')
s3 <- riskScatter(diff.true, diff.glmcc1, 'C)') 
s4 <- riskScatter(diff.true, diff.glmcc1, 'D)')
grid.arrange(s1, s2, s3, s4, ncol=2)
