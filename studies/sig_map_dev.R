##############################
# New method to calculate high
# resolution significance maps
##############################

library(gstat)
library(mgcv)
library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


dst <- "/Users/brianconroy/Documents/research/project2/cdph_by_species/"
src <- "/Users/brianconroy/Documents/research/cdph/data/"
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=5)
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
output <- load_output("output_cdph_baseline.json")


#########################
#### Faster interpolaters 
#########################
all_ids <- 1:ncell(caPr.disc[[1]])
obs_ids <- all_ids[!is.na(caPr.disc[[1]][])]
coords <- xyFromCell(caPr.disc[[1]], obs_ids)
v <- data.frame(colMeans(output$samples.w))
names(v) <- "v"

dsp <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp, v)

# nearest neighbors
gs <- gstat(formula=v~1, locations=dsp, nmax=3, set=list(idp = 0))
nn <- interpolate(caPr[[1]], gs)
nnr <- mask(nn, caPr[[1]])
plot(nnr)

# idw
gs <- gstat(formula=v~1, locations=dsp)
idw <- interpolate(caPr[[1]], gs)
idr <- mask(idw, caPr[[1]])
plot(idr)

# gam
df <- data.frame(cbind(v, coords))
mod <- gam(formula=v~s(x,y), data=df, gamma=0.005)
df_new <- data.frame(xyFromCell(caPr[[1]], (1:ncell(caPr[[1]]))[!is.na(caPr[[1]][])]))
plot(overlay(predict(mod, newdata=df_new, type='response'), caPr[[1]]))
plot(overlay(predict(mod), caPr.disc[[1]]))
plot(overlay(colMeans(output$samples.w), caPr.disc[[1]]))

# fast tps
r <- overlay(colMeans(output$samples.w), caPr.disc[[1]])
xy <- data.frame(xyFromCell(r, 1:ncell(r)))
z <- getValues(r)
p <- raster(caPr[[2]])
tps <- fastTps(xy, z, theta=2, na.rm=T, lambda=0)
p <- interpolate(p, tps)
p <- mask(p, caPr[[1]])
plot(p)

# original tps
r <- overlay(colMeans(output$samples.w), caPr.disc[[1]])
xy <- data.frame(xyFromCell(r, 1:ncell(r)))
z <- getValues(r)
p <- raster(caPr[[2]])
tps <- Tps(xy, z)
p <- interpolate(p, tps)
p <- mask(p, caPr[[1]])
plot(p)

# kernel smoothing
library(np)
library(tictoc)

template <- caPr.disc[[1]]
txdat <- data.frame(xyFromCell(template, (1:ncell(template))[!is.na(template[])]))
z = colMeans(output$samples.w)
x = txdat[,1]
y = txdat[,2]
df_new <- data.frame(xyFromCell(caPr[[1]], (1:ncell(caPr[[1]]))[!is.na(caPr[[1]][])]))

# calculate bandwidth
bw <- npregbw(formula=z~x+y,
              regtype="lc",
              bwmethod="cv.ls",
              ckertype="gaussian")

tic("smoothing")
model.np <- npreg(bws=c(0.08, 0.09),
                  formula=z~x+y,
                  regtype="lc",
                  ckertype="gaussian")
pred <- predict(model.np, newdata=df_new)
toc()

plot(overlay(pred, caPr[[1]]))
plot(overlay(z, caPr.disc[[1]]))
plot(overlay(predict(hen), caPr.disc[[1]]))

# which regression type to use?
#   lc specifies a local-constant estimator (Nadaraya-Watson) 
#   ll specifies a local-linear estimator. 

# bw selection: which cv metric?
#   cv.aic specifies expected KullbackLeibler cross-validation (Hurvich, Simonoff, and Tsai (1998))
#   cv.ls specifies least-squares cross-validation.

# which kernel to use?
#   second order gaussian 

# test run
ns <- 100
samples.pred <- array(NA, c(ns, sum(!is.na(caPr[[1]][]))))
progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*ns)
tic("smoothing")
for (i in 1:ns){
  z = output$samples.w[i,]
  model.np = npreg(bws=c(0.08, 0.09),
                   formula=z~x+y,
                   regtype="lc",
                   ckertype="gaussian")
  pred <- predict(model.np,
                  newdata = df_new)
  samples.pred[i,] <- pred
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/ns)
  }
}
toc()

plot(overlay(colMeans(samples.pred), caPr[[1]]))
plot(overlay(apply(samples.pred, 2, var), caPr[[1]]))

# in parallel
library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

tic("smoothing")
finalMatrix <- foreach(i=1:ns, .combine=rbind) %dopar% {
  library(np)
  z = output$samples.w[i,]
  model.np = npreg(bws=c(0.08, 0.09),
                   formula=z~x+y,
                   regtype="lc",
                   ckertype="gaussian")
  pred <- predict(model.np,
                  newdata = df_new)
  pred
}
toc()

plot(overlay(colMeans(samples.pred), caPr[[1]]), main='A)')
plot(overlay(colMeans(finalMatrix), caPr[[1]]), main='B')

# implement 3 new functions
# 1 function to generate interpolated random field samples
# args: samples.w, bandwidths, r_train, r_pred
# have it save samples
# function to calc high res posterior samples
# function to calc sig maps
