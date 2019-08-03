###########################
# Summarizes results of the 
# by species (Scurid)
# cdph analysis
###########################

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
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")


######################
#### Significance maps
######################


analysis_name <- "cdph_temporal_analysis"
output_agg <- load_output(paste('output_', paste(analysis_name, "_aggregated_ps", sep=""), ".json", sep=""))
caPr.disc <- aggregate(caPr, fact=7)

samples.w <- output_agg$samples.w
samples.theta  <- output_agg$samples.theta
samples.phi <- output_agg$samples.phi
n.samples <- nrow(samples.w)

# knots
loc.disc <- caPr.disc[[1]]
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
all_coords <- xyFromCell(loc.disc, cell=all_ids)
m <- 100
km.knots <- kmeans(all_coords, m)$centers
knot_ids <- cellFromXY(loc.disc, km.knots)
w_ids <- (1:length(all_ids))[all_ids %in% knot_ids]
knot_coords <- xyFromCell(loc.disc, cell=knot_ids)
plot(overlay(colMeans(samples.w), loc.disc))
points(knot_coords, pch=16, cex=0.5)
d <- as.matrix(dist(knot_coords, diag=TRUE, upper=TRUE))

# prediction points
loc.disc_ <- caPr[[2]]
pred_ids <- c(1:length(loc.disc_[]))[!is.na(loc.disc_[])]
pred_coords <- xyFromCell(loc.disc_, cell=pred_ids)

# precompute distances from pred points to knot points
d_lk <- array(NA, c(length(pred_ids), length(knot_ids)))
observed <- c()
for (l in 1:length(pred_ids)){
  for (k in 1:length(knot_ids)){
    dis <- dist(rbind(pred_coords[l,], knot_coords[k,]))
    d_lk[l, k] <- dis
    if (round(dis, 4) == 0){
      observed <- rbind(observed, c(l, k))
    }
  }
}


n.s <- 400
samples.pred <- array(NA, c(n.s, length(pred_ids)))

progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.s)
for (i in 1:n.s){
  
  w.i <- matrix(samples.w[i,])
  w.i.knot <- w.i[w_ids]
  theta.i <- samples.theta[i]
  phi.i <- samples.phi[i]
  
  C <- Exponential(d, range=theta.i, phi=phi.i)
  C.inv <- solve(C)
  
  # iterate through locations
  for (l in 1:length(pred_ids)){
    
    if (!(l  %in% observed[,1])){
      
      # calculate covariance of point l with each knot point
      c <- matrix(Exponential(d_lk[l,], range=theta.i, phi=phi.i))
      a <- t(c) %*% C.inv
      mu <- a %*% w.i.knot
      var <- phi.i - a %*% c
      samples.pred[i, l] <- rnorm(1, mu, sd=sqrt(var))
      
    } else {
      
      idx <- observed[observed[,1]==l][2]
      samples.pred[i, l] <- w.i.knot[k]
      
    }
    
  }
    
  if(i %in% percentage.points){
    setTxtProgressBar(progressBar, i/n.s)
  }
  
}

pred_mu <- colMeans(samples.pred)
pred_var <- apply(samples.pred, 2, var)
plot(overlay(pred_mu, caPr[[2]]))
points(knot_coords, pch=16, cex=0.65)

w.hat <- colMeans(samples.w[1:n.s,])
plot(overlay(w.hat, caPr.disc[[2]]))
points(knot_coords, pch=16, cex=0.65)

