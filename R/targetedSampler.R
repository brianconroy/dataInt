
# how to define targeted covariates
# targeted locations are point processes with intensity
# functions informed by: covariates, and, the disease process

# place the covariate target grids
# sample survey locations
# collect points around each survey in some box

# better to sample from an IPP
# within (possibly multiple) prespecified
# extents of the state

injectTarget <- function(x, ext, spdf=NULL){
  # x: rasterLayer
  # ext: extent defining the injection rectangle
  # spdf: study area boundary
  x_tar <- crop(x, ext)
  vals <- rnorm(ncell(x_tar), 3, 1)
  x_tar[] <- vals
  x_merge <- merge(x_tar, x)
  if (!is.null(spdf)){
    x_merge <- mask(x_merge, spdf)
  }
  names(x_merge) <- names(x)
  return(x_merge)
}

# options
# what is the actual use case?
# state health department surveys areas thought to have plague
# - evidenced by previous cases of plague, or being a campsite







