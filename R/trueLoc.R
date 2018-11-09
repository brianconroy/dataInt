library(USAboundaries)
library(spatstat)
library(maptools)
library(raster)
library(sp)


# Generate true case locations according to an inhomogeneous Poisson point process
# with intensity function informed by bioclimatic (worldclim) covariates. 
# Point patterns sampled from the rpoispp function {spatstat}


readWC <- function(wcdir='Documents/research/dataInt/data/wc10_2/'){
  rasStack <- raster(paste(wcdir, 'wc2.0_bio_10m_01.tif', sep=''))
  for (f in list.files(wcdir)[2:19]){
    rasF <- raster(paste(wcdir, f, sep=''))
    rasStack <- stack(rasStack, rasF)
  }
  names(rasStack) <- paste('bio', seq(1,19), sep='')
  return(rasStack)
}


getState <- function(state){
  us <- readRDS('Documents/research/dataInt/data/GADM_2.8_USA_adm1.rds')
  st <- us[match(toupper(state), toupper(us$NAME_1)), ]
  return(st)
}


getStateWC <- function(state, wc){
  us <- readRDS('Documents/research/dataInt/data/GADM_2.8_USA_adm1.rds')
  datSub <- crop(wc, extent(state))
  datSt <- mask(datSub, ca)
  return(datSt)
}


intensFunc <- function(beta, rasterCovariates, intercept=NULL){
  intens <- beta[[1]]*rasterCovariates[[1]]
  if (length(beta) > 1){
    for (i in 2:length(beta)){
      intens <- intens + beta[[i]]*rasterCovariates[[i]] 
    }
  }
  if (!is.null(intercept)){
    intens <- intens + intercept
  }
  intensIm <- as.im.RasterLayer(exp(intens))
  return(intensIm)
}
