library(plyr)
library(R.utils)
library(ggplot2)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


#########################
# Tunes the preferential 
# sampling Bayesian model

# Considers effect of 
# preferential sampling
# when location and
# conditional intensities
# share a covariate
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
beta.samp <- c(-3, 1, 1)
loc.disc <- caWc.disc[[c(6, 15)]]
locs <- simBernoulliLoc(loc.disc, beta.samp, seed=42)
sum(locs$status)


#### Disease covariate surface
beta.case <- c(1, 1, 2)
cov.disc <- caWc.disc[[c(12, 15)]]


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


#### Iterate over alpha
results.glm <- list()
counter <- 1
for (a in c(0, 0.15, 0.5, 2, 5, 10)){
  
  M <- 500
  samples.1 <- array(NA, c(M, length(beta.case)))
  samples.2 <- array(NA, c(M, length(beta.case) + 1))
  
  for (i in 1:M){
    
    data.a <- simLocCond(cov.disc, loc.disc, beta.case, beta.samp, a, global.center=TRUE)
    locs.x.sub <- data.a$loc$x.scaled[as.logical(data.a$loc$status),]
    locs.x.sub <- locs.x.sub[,2]
    mod.1 <- glm(data.a$conditional$y ~ data.a$conditional$x.standardised - 1, family='poisson')
    mod.2 <- glm(data.a$conditional$y ~ data.a$conditional$x.standardised + locs.x.sub - 1, family='poisson')
    samples.1[i,] <- mod.1$coefficients
    samples.2[i,] <- mod.2$coefficients

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
write.table(df.glm[c('a', 'model', 
                     'mu.b0', 'sd.b0', 'pbias.b0',
                     'mu.b1', 'sd.b1', 'pbias.b1',
                     'mu.b2', 'sd.b2', 'pbias.b2')],
            file='Documents/research/dataInt/output/pS_iterate_a_shared.csv', sep=',', row.names=F)
df.glm.1 <- df.glm[df.glm$model == 'glm',]
df.glm.2 <- df.glm[df.glm$model == 'glm location', ]


reformat <- function(df){
  
  df.glm.b0 <- df[,c('a', 'mu.b0', 'sd.b0', 'pbias.b0')]
  df.glm.b1 <- df[,c('a', 'mu.b1', 'sd.b1', 'pbias.b1')]
  df.glm.b2 <- df[,c('a', 'mu.b2', 'sd.b2', 'pbias.b2')]
  df.glm.b0$parameter <- 'beta0'
  df.glm.b1$parameter <- 'beta1'
  df.glm.b2$parameter <- 'beta2'
  names(df.glm.b0) <- c('a', 'mu', 'sd', 'pbias', 'parameter')
  names(df.glm.b1) <- names(df.glm.b0)
  names(df.glm.b2) <- names(df.glm.b0)
  df.glm.new <- rbind(df.glm.b0, df.glm.b1, df.glm.b2)
  return(df.glm.new)
  
}


df.glm.new1 <- reformat(df.glm.1)
df.glm.new2 <- reformat(df.glm.2)


p1 <- ggplot(data=df.glm.new1, aes(x=a, y=df.glm.new1$pbias, group=parameter))+
  geom_line(aes(color=parameter))+
  geom_point(aes(color=parameter))+ 
  labs(x ='alpha')+
  labs(y='percent bias')+
  ggtitle('A)')+ 
  ylim(-250, 250)

p2 <- ggplot(data=df.glm.new2, aes(x=a, y=df.glm.new2$pbias, group=parameter))+
  geom_line(aes(color=parameter))+
  geom_point(aes(color=parameter))+ 
  labs(x ='alpha')+
  labs(y='percent bias')+
  ggtitle('B)')+
  ylim(-250, 250)

grid.arrange(p1, p2, ncol=2)


#### chosen alpha
alpha <- 2


#### Simulate counts given locations
count.data <- simConditional(cov.disc, locs, beta.case, beta.samp, alpha, seed=42)
data <- list(conditional=count.data, loc=locs)


par(mfrow=c(1,1))
plot(prob.rast, col=pal2(8))
points(locs$coords, cex=count.data$y/10)


hist(count.data$y)
mean(count.data$y)


##############################
## Preferential Sampling Model
##############################


output <- prefSample(
  data, n.sample=300000, burnin=10000, thin=1,
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


df.glm.a1 <- df.glm.new1[df.glm.new1$a == alpha,]
df.glm.a2 <- df.glm.new2[df.glm.new2$a == alpha,]
df.glm.a1$model <- 'glm'
df.glm.a2$model <- 'glm location'
names(df.glm.a1) <- c('a', 'estimate', 'sd', 'pbias', 'parameter', 'model')
names(df.glm.a2) <- c('a', 'estimate', 'sd', 'pbias', 'parameter', 'model')
df.glm.a1 <- df.glm.a1[c('model', 'parameter', 'estimate', 'sd', 'pbias')]
df.glm.a2 <- df.glm.a2[c('model', 'parameter', 'estimate', 'sd', 'pbias')]


df.comp <- rbind(df.ps, df.glm.a1, df.glm.a2)
df.comp <- df.comp[with(df.comp, order(parameter)), ]
write.table(df.comp, file='Documents/research/dataInt/pS_summary_comp_shared.csv', sep=',', row.names=F)
