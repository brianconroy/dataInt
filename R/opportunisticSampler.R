library(USAboundaries)
library(spatstat)
library(maptools)
library(raster)
library(sp)


# st: SpatialPolygonsDataFrame for a state
# vals: numeric of raster values
makeSampRast <- function(numrow, numcol, st, vals){
  r <- raster(nrow=numrow,
              ncol=numcol,
              xmn=extent(ca)@xmin,
              xmx=extent(ca)@xmax,
              ymn=extent(ca)@ymin,
              ymx=extent(ca)@ymax)
  r <- crop(r, extent(st))
  r[] <- vals
  r <- mask(r, st)
  return(r)
}


# values from a North-South sampling gradient
# x: rasterLayer
# demarclat: latitude to demarcate the north south bound
sampGradNS <- function(x, demarclat){
  ext <- extent(x)
  fraclat <- (ext@ymax - demarclat)/(ext@ymax - ext@ymin)
  reps1 <- round(fraclat * nrow(x))
  reps2 <- nrow(x) - reps1
  vals <- c(rep(1, ncol(x)*reps1), rep(-1, ncol(x)*reps2))
  return(vals)
}


# r.mold: RasterBrick to serve as a "mold" for a new raster
# minval: minimum value of the sampling covariate
# maxval: maximum value of the sampling covariate
nsGradSmooth <- function(r.mold, minval, maxval, cov.name, seed=NULL){
  
  # setup the new raster
  r.sampling <- r.mold
  vals.samp <- values(r.sampling)
  num.vals <- length(vals.samp[!is.na(vals.samp)])
  increment <- (maxval - minval)/num.vals
  running.val <- minval
  
  # assign new values
  for (i in 1:length(vals.samp)){
    if (!is.na(vals.samp)[i]){
      vals.samp[i] <- running.val #+ noise[i]
      running.val <- running.val + increment
    }
  }
  
  if (!is.null(seed)){ set.seed(seed) }
  noise <- runif(num.vals, 0, 0.1*maxval)
  vals.samp[!is.na(vals.samp)] <- vals.samp[!is.na(vals.samp)] + noise
  
  r.sampling[] <- vals.samp
  names(r.sampling) <- cov.name
  return(r.sampling)
  
}


# produces values from a North-South sampling gradient
# ext: extent
# demarclat: latitude to demarcate the north south bound
sampGradSmooth <- function(numrow, numcol, minval, maxval){
  increments <- seq(minval, maxval, by=(maxval - minval)/numrow)
  vals <- c()
  for (i in 1:numrow){
    vals <- c(vals, rep(increments[i], numcol))
  }
  return(vals)
}


# x: rasterLayer or rasterBrick of sampling covariates
# betas: list of corresponding beta values
prObserve <- function(x, betas){
  if (class(covRast1)[1] == "RasterLayer"){
    biasProb = exp(covRast1 * betas[[1]])/(1 + exp(covRast1 * betas[[1]]))
    return(biasProb)
  }
}


# instead of an intensity + density approach just 
# define a probability surface
# biasDensity <- function(x){
#   x: rasterLayer of the bias intensity surface
#   cellAreas <- area(x)
#   normalizingConst <- sum(na.omit(values(cellAreas * values(x))))
#   densRast <- x/normalizingConst
#   return(densRast)
# }


# biasSurface: RasterLayer defining the bias density surface
# pp: planar point pattern
# returns: list of observation indicators for each point
thinPoints <- function(biasSurface, pp){
  densVals <- raster::extract(biasSurface, coords(pp))
  inds <- sapply(densVals, function(x){ rbinom(1, 1, x)})
  inds[is.na(inds)] <- 0
  thinned <- pp[as.logical(inds), ]
  return(thinned)
}
