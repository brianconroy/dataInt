library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


#########################
# Tunes the preferential 
# sampling Bayesian model

# Updates locational parameters
# based on the joint likelihood,
# not merely the locational 
# component
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
beta.samp <- c(-1.5, 1)
loc.disc <- caWc.disc[[c(6)]]
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
sum(locs$status)


#### Disease covariate surface
cov.disc <- caWc.disc[[c(12)]]


#### Visualize surfaces
pal1 <- colorRampPalette(c("blue","red"))
pal2 <- colorRampPalette(c("blue","green"))


# raw covariates
par(mfrow=c(1,2))
plot(cov.disc, main='A)', col=pal1(8))
plot(loc.disc, main='C)', col=pal2(8))


# sampling intensity
par(mfrow=c(1,1))
prob.rast <- loc.disc
prob.rast[!is.na(values(loc.disc))] <- locs$probs
plot(prob.rast, col=pal2(8))
points(locs$coords, pch=16)


#### True disease parameter values
beta.case <- c(2, 2)


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
            file='Documents/research/dataInt/pS_iterate_a_pos.csv', sep=',', row.names=F)
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
  ggtitle('C)')

p2 <- ggplot(data=df.glm.new2, aes(x=a, y=df.glm.new2$pbias, group=parameter))+
  geom_line(aes(color=parameter))+
  geom_point(aes(color=parameter))+ 
  labs(x ='alpha')+
  labs(y='percent bias')+
  ggtitle('D)')
grid.arrange(p1, p2, ncol=2)


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
