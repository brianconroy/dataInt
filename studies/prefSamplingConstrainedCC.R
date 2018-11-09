library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


#########################
# Tunes the preferential 
# sampling Bayesian model
# with case and control 
# counts
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


summarizeSamps <- function(samples, beta){
  
  resultsList <- list()
  
  mu.b0 <- round(mean(samples[,1]), 3)
  resultsList$mu.b0 <- mu.b0
  resultsList$sd.b0 <- round(sd(samples[,1]), 3)
  resultsList$pbias.b0 <- round(100*(mu.b0 - beta[1])/beta[1], 3)
  
  mu.b1 <- round(mean(samples[,2]), 3)
  resultsList$mu.b1 <- mu.b1
  resultsList$sd.b1 <- round(sd(samples[,2]), 3)
  resultsList$pbias.b1 <- round(100*(mu.b1 - beta[2])/beta[2], 3)
  
  resultsList
}


#### Iterate over alpha
results.glm <- list()
counter <- 1
for (a in c(0, 0.15, 0.5, 2, 10, 20)){
  
  M <- 500
  samples.1 <- array(NA, c(M, 2))
  samples.2 <- array(NA, c(M, 2 + length(names(loc.disc))))
  samples.loc <- array(NA, c(M, length(names(loc.disc))))
  for (i in 1:M){
    
    data.a <- simLocCond(cov.disc, loc.disc, beta.case, beta.samp, a)
    locs.x.sub <- data.a$loc$x.scaled[as.logical(data.a$loc$status)]
    mod.1 <- glm(data.a$conditional$y ~ data.a$conditional$x.standardised - 1, family='poisson')
    mod.2 <- glm(data.a$conditional$y ~ data.a$conditional$x.standardised + locs.x.sub - 1, family='poisson')
    mod.loc <- glm(data.a$loc$status ~ data.a$loc$x.scaled-1, family='binomial')
    samples.1[i,] <- mod.1$coefficients
    samples.2[i,] <- mod.2$coefficients
    samples.loc[i,] <- mod.loc$coefficients
    
  }
  results.1 <- summarizeSamps(samples.1, beta.case)
  results.1$model <- 'glm'
  results.1$a <- a
  results.glm[[counter]] <- results.1
  counter <- counter + 1
  
  results.2 <- summarizeSamps(samples.2, beta.case)
  results.2$model <- 'glm location'
  results.2$a <- a
  results.glm[[counter]] <- results.2
  counter <- counter + 1
  
}
df.glm <- ldply(results.glm, data.frame)
write.table(df.glm[c('a', 'model', 'mu.b0', 'sd.b0', 'pbias.b0',
                     'mu.b1', 'sd.b1', 'pbias.b1')],
            file='Documents/research/dataInt/pS_iterate_a.csv', sep=',', row.names=F)
df.glm.1 <- df.glm[df.glm$model == 'glm',]
df.glm.2 <- df.glm[df.glm$model == 'glm location', ]


reformat <- function(df){
  
  df.glm.b0 <- df[,c('a', 'mu.b0', 'sd.b0', 'pbias.b0')]
  df.glm.b1 <- df[,c('a', 'mu.b1', 'sd.b1', 'pbias.b1')]
  df.glm.b0$parameter <- 'beta0.cond'
  df.glm.b1$parameter <- 'beta1.cond'
  names(df.glm.b0) <- c('a', 'mu', 'sd', 'pbias', 'parameter')
  names(df.glm.b1) <- names(df.glm.b0)
  df.glm.new <- rbind(df.glm.b0, df.glm.b1)
  return(df.glm.new)
  
}

df.glm.new1 <- reformat(df.glm.1)
df.glm.new2 <- reformat(df.glm.2)

p1 <- ggplot(data=df.glm.new1, aes(x=a, y=df.glm.new1$pbias, group=parameter))+
  geom_line(aes(color=parameter))+
  geom_point(aes(color=parameter))+ 
  labs(x ='alpha')+
  labs(y='percent bias')+
  ggtitle('A)')

p2 <- ggplot(data=df.glm.new2, aes(x=a, y=df.glm.new2$pbias, group=parameter))+
  geom_line(aes(color=parameter))+
  geom_point(aes(color=parameter))+ 
  labs(x ='alpha')+
  labs(y='percent bias')+
  ggtitle('B)')
grid.arrange(p1, p2, ncol=2)


#### chosen alpha
alpha <- 2


#### Simulate counts given locations
count.case <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha)
count.control <- simConditional(cov.disc, locs, beta.control, beta.samp, alpha)
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
points(locs$coords, cex=count.control$y/15)

plot(cov.disc, main='D)', mar=c(4, 3, 3, 2) + 0.1)
points(locs$coords, pch=16)
points(locs$coords, cex=count.control$y/15)


par(mfrow=c(2,1))
hist(count.case$y, mar=c(4, 4, 3, 2) + 0.1)
hist(count.control$y, mar=c(4, 4, 3, 2) + 0.1)
hist(count.control$y)


##############################
## Preferential Sampling Model
##############################


output <- prefSampleCC(
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
  alpha=alpha
)
par(mfrow=c(2,3))
viewOutput(output, type='preferential.cc', truevals, view.hist=F)
viewOutput(output, type='preferential.cc', truevals, view.trace=F)


################
## Multiple Runs
################


n.runs <- 8
par(mfrow=c(2,4))
results <- list()
for (i in 1:n.runs){
  
  print(paste("iteration", i))
  output <- prefSampleCC(
    data, n.sample=150000, burnin=15000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  results[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
  
}


pooled <- pool.results.pscc(results)
par(mfrow=c(2,3))
viewOutput(pooled, type='preferential.cc', truevals, view.hist=F)
viewOutput(pooled, type='preferential.cc', truevals, view.trace=F)


beta.case.pscc <- colMeans(pooled$samples.beta.case)
beta.control.pscc <- colMeans(pooled$samples.beta.control)


prefSummary <- summarize(pooled, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PS'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')

df.glm.a1 <- df.glm.new1[df.glm.new1$a == alpha,]
df.glm.a2 <- df.glm.new2[df.glm.new2$a == alpha,]
df.glm.a1$parameter <- c('beta0.case', 'beta1.case')
df.glm.a2$parameter <- c('beta0.case', 'beta1.case')
df.glm.a1$model <- 'glm'
df.glm.a2$model <- 'glm location'
names(df.glm.a1) <- c('a', 'estimate', 'sd', 'pbias', 'parameter', 'model')
names(df.glm.a2) <- c('a', 'estimate', 'sd', 'pbias', 'parameter', 'model')
df.glm.a1 <- df.glm.a1[c('model', 'parameter', 'estimate', 'sd', 'pbias')]
df.glm.a2 <- df.glm.a2[c('model', 'parameter', 'estimate', 'sd', 'pbias')]

df.comp <- rbind(df.ps, df.glm.a1, df.glm.a2)
df.comp <- df.comp[with(df.comp, order(parameter)),]
df.comp <- df.comp[c('parameter', 'model', 'estimate', 'sd', 'pbias')]
write.table(df.comp, file='Documents/research/dataInt/pS_summary_comp3.csv', sep=',', row.names=F)


### Compare Risk surfaces
  ## The case control PS model
  ## The case only PS model
  ## glm case control (non locational)


data.co <- data
data.co$conditional <- data.co$conditional.case
n.runs <- 8
par(mfrow=c(2,4))
results.co <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- prefSample(
    data.co, n.sample=250000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  results.co[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
}

pooled.co <- pool.results.ps(results.co)
names(pooled.co) <- c("samples.beta.l", "samples.beta.case", "samples.alpha",  "accept")
beta.case.psco <- colMeans(pooled.co$samples.beta.case)


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


cuts=c(-5, -2, 0, 2, 4, 6, 10)
pal <- colorRampPalette(c("blue","red"))


## True surface
intens.case.true <- calcIntens(cov.disc, beta.case)
intens.control.true <- calcIntens(cov.disc, beta.control)
diff.true <- rastDiff(intens.case.true, intens.control.true)
plot(diff.true, breaks=cuts, col=pal(8))


## PS CC surface
intens.case.pscc <- calcIntens(cov.disc, beta.case.pscc)
intens.control.pscc <- calcIntens(cov.disc, beta.control.pscc)
diff.pscc <- rastDiff(intens.case.pscc, intens.control.pscc)
plot(diff.pscc, breaks=cuts, col=pal(8))


## PS case only surface
intens.case.psco <- calcIntens(cov.disc, beta.case.psco)
plot(intens.case.psco, breaks=cuts, col=pal(8))


## glm case control surface
intens.case.glm <- calcIntens(cov.disc, beta.case.glm1)
intens.control.glm <- calcIntens(cov.disc, beta.control.glm1)
diff.glmcc <- rastDiff(intens.case.glm, intens.control.glm)
plot(diff.glmcc, breaks=cuts, col=pal(8))


## glm location case control surface
intens.case.glm2 <- calcIntens(cov.disc, beta.case.glm2)
intens.control.glm2 <- calcIntens(cov.disc, beta.control.glm2)
diff.glmcc2 <- rastDiff(intens.case.glm2, intens.control.glm2)
plot(diff.glmcc2, breaks=cuts, col=pal(8))


par(mfrow=c(2,3))
plot(diff.true, breaks=cuts, col=pal(8), main='A)')
plot(diff.pscc, breaks=cuts, col=pal(8), main='B)')
plot(intens.case.psco, breaks=cuts, col=pal(8), main='C)')
plot(intens.case.glm, breaks=cuts, col=pal(8), main='D)')
plot(diff.glmcc, breaks=cuts, col=pal(8), main='E)')


summs <- list()
summs[[1]] <- summarizeSurf(diff.true, diff.pscc, 'PS case control')
summs[[2]] <- summarizeSurf(diff.true, intens.case.psco, 'PS case only')
summs[[3]] <- summarizeSurf(diff.true, intens.case.glm, 'GLM case only')
summs[[4]] <- summarizeSurf(diff.true, diff.glmcc, 'GLM case control')
summs[[5]] <- summarizeSurf(diff.true, diff.glmcc2, 'GLM location case control')
df.risks <- ldply(summs, data.frame)
write.table(df.risks, file='Documents/research/dataInt/pS_summary_comp6.csv', sep=',', row.names=F)


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

s1 <- riskScatter(diff.true, diff.pscc, 'A)')
s2 <- riskScatter(diff.true, intens.case.psco, 'B)') 
s3 <- riskScatter(diff.true, intens.case.glm, 'C)')
s4 <- riskScatter(diff.true, diff.glmcc, 'C)')
grid.arrange(s1, s2, s3, s4, ncol=2)


##############################
# Iterate over different 
# alpha_case and alpha_control
##############################


intens.case.true <- calcIntens(cov.disc, beta.case)
intens.control.true <- calcIntens(cov.disc, beta.control)
diff.true <- rastDiff(intens.case.true, intens.control.true)
plot(diff.true, breaks=cuts, col=pal(8))


alpha.iter <- list()
counter <- 1
for (a.ca in c(1, 2, 4, 6, 8)){
  for (a.co in c(8, 6, 4, 2, 1, -1, -2, -4, -6, -8)){

    
    count.case <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha=a.ca)
    count.control <- simConditional(cov.disc, locs, beta.control, beta.samp, alpha=a.co)
    data <- list(conditional.case=count.case, conditional.control=count.control, loc=locs)
    
    
    locs.x.sub <- data$loc$x.scaled[as.logical(data$loc$status)]
    mod.1.ca <- glm(data$conditional.case$y ~ data$conditional.case$x.standardised - 1, family='poisson')
    mod.2.ca <- glm(data$conditional.case$y ~ data$conditional.case$x.standardised + locs.x.sub - 1, family='poisson')
    beta.case.glm1 <- unname(coefficients(mod.1.ca))
    beta.case.glm2 <- unname(coefficients(mod.2.ca))[1:2]
    
    
    mod.1.co <- glm(data$conditional.control$y ~ data$conditional.control$x.standardised - 1, family='poisson')
    mod.2.co <- glm(data$conditional.control$y ~ data$conditional.control$x.standardised + locs.x.sub - 1, family='poisson')
    beta.control.glm1 <- unname(coefficients(mod.1.co))
    beta.control.glm2 <- unname(coefficients(mod.2.co))[1:2]
    
    
    intens.case.glm <- calcIntens(cov.disc, beta.case.glm1)
    intens.control.glm <- calcIntens(cov.disc, beta.control.glm1)
    diff.glmcc <- rastDiff(intens.case.glm, intens.control.glm)

    
    intens.case.glm2 <- calcIntens(cov.disc, beta.case.glm2)
    intens.control.glm2 <- calcIntens(cov.disc, beta.control.glm2)
    diff.glmcc2 <- rastDiff(intens.case.glm2, intens.control.glm2)
    
    
    summs1 <- summarizeSurf(diff.true, diff.glmcc, 'GLM case control')
    summs2 <- summarizeSurf(diff.true, diff.glmcc2, 'GLM location case control')
    summs1$alphacase <- a.ca
    summs1$alphacontrol <- a.co
    summs2$alphacase <- a.ca
    summs2$alphacontrol <- a.co
    alpha.iter[[counter]] <- summs1
    alpha.iter[[counter + 1]] <- summs2
    counter <- counter + 2
    
    
  }
}

df.alpha <- ldply(alpha.iter, data.frame)
df.1 <- df.alpha[df.alpha$model == 'GLM case control',]
df.2 <- df.alpha[df.alpha$model == 'GLM location case control',]


p1 <- ggplot(df.alpha, aes(x=alphacase, y=alphacontrol, size=mse, group=model)) +
  geom_point(shape=21, aes(color=model)) + 
  scale_size(range = c(1,25)) + 
  ggtitle('A)')
p2 <- ggplot(df.alpha, aes(x=alphacase, y=alphacontrol, size=mae, group=model)) +
  geom_point(shape=21, aes(color=model)) + 
  scale_size(range = c(1,10)) + 
  ggtitle('B)')
grid.arrange(p1, p2, ncol=2)
