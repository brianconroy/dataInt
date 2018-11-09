library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


####################################
# Fits the differential preferential 
# sampling Bayesian model, with case 
# and control counts

# Considers effect of preferential 
# sampling when location and
# conditional intensities share a 
# covariate
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
beta.samp <- c(-3, 1, 1)
loc.disc <- caWc.disc[[c(6, 15)]]
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
sum(locs$status)


#### Disease covariate surface
cov.disc <- caWc.disc[[c(12, 15)]]


#### Count parameter values
beta.case <- c(1.5, 1.5, 2)
beta.control <- c(1.5, -1.5, 2)


#### Visualize surfaces
pal1 <- colorRampPalette(c("blue","red"))
pal2 <- colorRampPalette(c("blue","green"))


# raw covariates
par(mfrow=c(2,2))
plot(cov.disc[[1]], main='A)', col=pal1(8))
plot(cov.disc[[2]], main='B)')
plot(loc.disc[[1]], main='C)', col=pal2(8))
plot(loc.disc[[2]], main='D)')


# sampling intensity
par(mfrow=c(1,1))
prob.rast <- loc.disc[[1]]
prob.rast[!is.na(values(loc.disc))] <- locs$probs
plot(prob.rast, col=pal2(8))
points(locs$coords, pch=16)


#### chosen alphas
alpha.case <- 3
alpha.control <- -1


#### Simulate counts given locations
count.case <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha.case)
count.control <- simConditional(cov.disc, locs, beta.control, beta.samp, alpha.control)
data <- list(conditional.case=count.case, conditional.control=count.control, loc=locs)


par(mfrow=c(2,1))
hist(count.case$y)
hist(count.control$y)


par(mfrow=c(1,2))
plot(prob.rast, main='A)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.case$y/15)


plot(prob.rast, main='B)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.control$y/15)


par(mfrow=c(2,1))
hist(count.case$y, mar=c(4, 4, 3, 2) + 0.1)
hist(count.control$y, mar=c(4, 4, 3, 2) + 0.1)


##############################
## Preferential Sampling Model
##############################


output <- prefSampleCC(
  data, n.sample=500000, burnin=10000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
)
print(output$accept)


truevals <- list(
  beta1.loc=beta.samp,
  beta0.case=beta.case[1],
  beta1.case=beta.case[2],
  beta2.case=beta.case[3],
  beta0.control=beta.control[1],
  beta1.control=beta.control[2],
  beta2.control=beta.control[3],
  alpha.case=alpha.case,
  alpha.control=alpha.control
)
par(mfrow=c(2,4))
viewOutput(output, type='preferential.cc.diff', truevals, view.hist=F)
viewOutput(output, type='preferential.cc.diff', truevals, view.trace=F)


beta.case.psccd <- colMeans(output$samples.beta.case)
beta.control.psccd <- colMeans(output$samples.beta.control)


prefSummary <- summarize(output, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PSD'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')


###################
# Comparison Models
###################


summarizeGlm <- function(mod, mod.name, tag, truevals, max.pos=3){
  
  est <- round(unname(coefficients(mod)), 3)[1:max.pos]
  sds <- round(unname(sqrt(diag(summary(mod)$cov.scaled))), 3)
  
  summ <- list()
  for (i in 1:length(est)){
    param <- paste('beta', i-1, '.', tag, sep='')
    true.i <- unlist(truevals[param])
    summ.i <- list()
    summ.i$model <- mod.name
    summ.i$parameter <- param
    summ.i$estimate <- est[i]
    summ.i$sd <- sds[i]
    summ.i$pbias <- round(100 * (est[i] - true.i)/true.i, 3)
    summ[[i]] <- summ.i
  }
  
  return(ldply(summ, data.frame))
  
}


## glm fits
locs.x.sub <- data$loc$x.scaled[as.logical(data$loc$status),]
locs.x.sub <- locs.x.sub[,2]
mod.1.ca <- glm(data$conditional.case$y ~ data$conditional.case$x.standardised - 1, family='poisson')
mod.2.ca <- glm(data$conditional.case$y ~ data$conditional.case$x.standardised + locs.x.sub - 1, family='poisson')
beta.case.glm1 <- unname(coefficients(mod.1.ca))
beta.case.glm2 <- unname(coefficients(mod.2.ca))[1:3]

mod.1.co <- glm(data$conditional.control$y ~ data$conditional.control$x.standardised - 1, family='poisson')
mod.2.co <- glm(data$conditional.control$y ~ data$conditional.control$x.standardised + locs.x.sub - 1, family='poisson')
beta.control.glm1 <- unname(coefficients(mod.1.co))
beta.control.glm2 <- unname(coefficients(mod.2.co))[1:3]


df.glm <- data.frame()
df.glm <- rbind(df.glm, summarizeGlm(mod.1.ca, 'GLM', 'case', truevals))
df.glm <- rbind(df.glm, summarizeGlm(mod.2.ca, 'GLM location', 'case', truevals))
df.glm <- rbind(df.glm, summarizeGlm(mod.1.co, 'GLM', 'control', truevals))
df.glm <- rbind(df.glm, summarizeGlm(mod.2.co, 'GLM location', 'control', truevals))

df.comp <- rbind(df.ps, df.glm)
df.comp <- df.comp[with(df.comp, order(parameter)),]
df.comp <- df.comp[c('parameter', 'model', 'estimate', 'sd', 'pbias')]
write.table(df.comp, file='Documents/research/dataInt/output/pS_diff_comp.csv', sep=',', row.names=F)


#######################
# Compare Risk surfaces
#######################


pal <- colorRampPalette(c("blue","red"))


## True surface
intens.case.true <- calcIntens(cov.disc, beta.case)
intens.control.true <- calcIntens(cov.disc, beta.control)
diff.true <- rastDiff(intens.case.true, intens.control.true)
plot(diff.true, col=pal(11))


## PS CC surface
intens.case.psccd <- calcIntens(cov.disc, beta.case.psccd)
intens.control.psccd <- calcIntens(cov.disc, beta.control.psccd)
diff.psccd <- rastDiff(intens.case.psccd, intens.control.psccd)
plot(diff.psccd, col=pal(11))


## glm case control surface
intens.case.glm1 <- calcIntens(cov.disc, beta.case.glm1)
intens.control.glm1 <- calcIntens(cov.disc, beta.control.glm1)
diff.glmcc1 <- rastDiff(intens.case.glm1, intens.control.glm1)
plot(diff.glmcc1, col=pal(11))


## glm location case control surface
intens.case.glm2 <- calcIntens(cov.disc, beta.case.glm2)
intens.control.glm2 <- calcIntens(cov.disc, beta.control.glm2)
diff.glmcc2 <- rastDiff(intens.case.glm2, intens.control.glm2)
plot(diff.glmcc2, col=pal(11))


par(mfrow=c(2,2))
plot(diff.true, col=pal(11), main='A)')
plot(diff.psccd, col=pal(11), main='B)')
plot(diff.glmcc1, col=pal(11), main='D)')
plot(diff.glmcc2, col=pal(11), main='E)')


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
s2 <- riskScatter(diff.true, diff.glmcc1, 'B)') 
s3 <- riskScatter(diff.true, diff.glmcc2, 'C)')
grid.arrange(s1, s2, s3, ncol=3)
