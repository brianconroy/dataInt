library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


#########################
# Tunes the preferential 
# sampling Bayesian model
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


#### 2 preference scenaries
#### caWc.disc[[6]] and disease is caWc.disc[[12]], alpha > 0
#### caWc.disc[[6]] and disease is caWc.disc[[1]], alpha < 0


#### Simulate survey locations
beta.samp <- -1
loc.disc <- caWc.disc[[c(6)]]
values(loc.disc) <- values(loc.disc) + abs(min(values(loc.disc), na.rm=T))
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
sum(locs$status)


#### Disease covariate surface
cov.disc <- caWc.disc[[c(1)]] # for negative preferential sampling
# cov.disc <- caWc.disc[[c(12)]] # for positive preferential sampling


prob.rast <- loc.disc
prob.rast[!is.na(values(loc.disc))] <- locs$probs
par(mfrow=c(1,3))
plot(cov.disc, main='A)')
plot(loc.disc, main='B)')
plot(prob.rast, main='C)')
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
alpha <- 2


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


output <- prefSample(
  data, n.sample=75000, burnin=10000, thin=1,
  proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
  proposal.sd.alpha=0.05, self.tune=TRUE
)


truevals <- list(
  beta1.loc=beta.samp,
  beta0.cond=beta.case[1],
  beta1.cond=beta.case[2],
  alpha=alpha
)
print(output$accept)
par(mfrow=c(2,4))
viewOutput(output, type='preferential.alpha', truevals)


################
## Multiple Runs
################


n.runs <- 4
par(mfrow=c(2,2))
results <- list()
for (i in 1:n.runs){
  print(paste("iteration", i))
  output <- prefSample(
    data, n.sample=250000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  results[[i]] <- output
  plot(output$samples.beta.l[,1], type='l', main='beta0')
  abline(h=beta.samp, col='2')
}


pool.results <- function(results){
  
  results.new <- list()
  new.beta.l <- results[[1]]$samples.beta.l
  new.beta.c <- results[[1]]$samples.beta.c
  new.alpha <- results[[1]]$samples.alpha
  new.accept <- results[[1]]$accept
  
  for (i in 2:length(results)){
    new.beta.l <- rbind(new.beta.l, results[[i]]$samples.beta.l)  
    new.beta.c <- rbind(new.beta.c, results[[i]]$samples.beta.c)
    new.alpha <- rbind(new.alpha, results[[i]]$samples.alpha)
    new.accept <- new.accept + results[[i]]$accept
  }
  
  new.accept <- new.accept/length(results)
  results.new$samples.beta.l <- new.beta.l
  results.new$samples.beta.c <- new.beta.c
  results.new$samples.alpha <- new.alpha
  results.new$accept <- new.accept
  
  return(results.new)
  
}


pooled <- pool.results(results)
par(mfrow=c(2,4), oma=c(1,1,1,1))
viewOutput(pooled, type='preferential.alpha', truevals)


prefSummary <- summarize(pooled, truevals, dic=FALSE)
df.ps <- ldply(prefSummary, data.frame)
df.ps$model <- 'PS'
df.ps <- df.ps[c('model', 'parameter', 'posterior.mean', 'posterior.sd', 'percbias')]
names(df.ps) <- c('model', 'parameter', 'estimate', 'sd', 'pbias')

df.glm.a1 <- df.glm.new1[df.glm.new1$a == alpha,]
df.glm.a2 <- df.glm.new2[df.glm.new2$a == alpha,]
df.glm.a1$model <- 'glm'
df.glm.a2$model <- 'glm location'
names(df.glm.a1) <- c('a', 'estimate', 'sd', 'pbias', 'parameter', 'model')
names(df.glm.a2) <- c('a', 'estimate', 'sd', 'pbias', 'parameter', 'model')
df.glm.a1 <- df.glm.a1[c('model', 'parameter', 'estimate', 'sd', 'pbias')]
df.glm.a2 <- df.glm.a2[c('model', 'parameter', 'estimate', 'sd', 'pbias')]

df.comp <- rbind(df.ps, df.glm.a1, df.glm.a2)
write.table(df.comp, file='Documents/research/dataInt/pS_summary_comp1.csv', sep=',', row.names=F)


#######################
## Iterate Preferential
## Sampling Model over 
## various alpha
#######################


summaries <- data.frame()
for (a in c(1, 2, 4, 6, 8)){
  
  print(a)
  count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, a)
  data <- list(conditional=count.data, loc=locs)
  
  output.a <- prefSample(
    data, n.sample=75000, burnin=10000, thin=1,
    proposal.sd.beta.c=0.01, proposal.sd.beta.l=0.05,
    proposal.sd.alpha=0.05, self.tune=TRUE
  )
  
  truevals <- list(
    beta0.loc=beta.samp,
    beta0=beta.case[1],
    beta1=beta.case[2],
    alpha=a
  )
  
  viewOutput(output.a, type='preferential.alpha', truevals)
  prefSummary <- summarize(output.a, truevals, dic=FALSE)
  df.ps <- ldply(prefSummary, data.frame)
  df.ps$alpha <- a
  summaries <- rbind(summaries, df.ps)
  
}


par(mfrow=c(2,2))
for (p in unique(summaries$parameter)){
  
  summaries.p <- summaries[summaries$parameter == p,]
  uy <- max(1.1 * max(summaries.p$percbias), 5)
  ly <- min(0, 1.1 * min(summaries.p$percbias))
  plot(x=summaries.p$alpha, y=summaries.p$percbias, ylim=c(ly, uy), main=p)
  lines(x=summaries.p$alpha, y=summaries.p$percbias, type='l')
  abline(h=0, col=2)
  
}
