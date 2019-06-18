library(prism)
library(sp)
library(raster)
library(RStoolbox)
library(maps)
library(rgeos)
library(prism)
library(stringr)


set.seed(42) # the answer to life, the universe and everything

years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
options(prism.path="/Users/brianconroy/Documents/research/dataInt/data_time")

# Download files
get_prism_annual_new(type= "tmean",
                     years=years,
                     keepZip=FALSE)
get_prism_annual_new(type= "tmax",
                     years=years,
                     keepZip=FALSE)
get_prism_annual_new(type= "tmin",
                     years=years,
                     keepZip=FALSE)
get_prism_annual_new(type= "ppt",
                     years=years,
                     keepZip=FALSE)
get_prism_annual_new(type= "vpdmin",
                     years=years
                     ,keepZip=FALSE)
get_prism_annual_new(type= "vpdmax",
                     years=years,
                     keepZip=FALSE)
get_prism_annual_new(type= "tdmean",
                     years=years,
                     keepZip=FALSE)


calc_prism_pcs <- function(year){
  
  # Convert to Rasters
  files <- ls_prism_data(absPath=T)[,2]
  ppt_rast <- get_rast(files, 'ppt', year)
  tdmean_rast <- get_rast(files, 'tdmean', year)
  tmax_rast <- get_rast(files, 'tmax', year)
  tmean_rast <- get_rast(files, 'tmean', year)
  tmin_rast <- get_rast(files, 'tmin', year)
  vpdmax_rast <- get_rast(files, 'vpdmax', year)
  vpdmin_rast <- get_rast(files, 'vpdmin', year)
  
  # Reproject PRISM Rasters
  crs_us <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  reproj_ppt_rast <- raster::projectRaster(ppt_rast, crs = crs(crs_us))
  reproj_tdmean_rast <- raster::projectRaster(tdmean_rast, crs = crs(crs_us))
  reproj_tmax_rast <- raster::projectRaster(tmax_rast, crs = crs(crs_us))
  reproj_tmean_rast <- raster::projectRaster(tmean_rast, crs = crs(crs_us))
  reproj_tmin_rast <- raster::projectRaster(tmin_rast, crs = crs(crs_us))
  reproj_vpdmax_rast <- raster::projectRaster(vpdmax_rast, crs = crs(crs_us))
  reproj_vpdmin_rast <- raster::projectRaster(vpdmin_rast, crs = crs(crs_us))
  
  # Scale Rasters by Range Transformation (can do other scaling)
  stand_reproj_ppt_rast <- (reproj_ppt_rast-min(na.omit(reproj_ppt_rast@data@values)))/(max(na.omit(reproj_ppt_rast@data@values))-min(na.omit(reproj_ppt_rast@data@values)))
  stand_reproj_tdmean_rast <- (reproj_tdmean_rast-min(na.omit(reproj_tdmean_rast@data@values)))/(max(na.omit(reproj_tdmean_rast@data@values))-min(na.omit(reproj_tdmean_rast@data@values)))
  stand_reproj_tmax_rast <- (reproj_tmax_rast-min(na.omit(reproj_tmax_rast@data@values)))/(max(na.omit(reproj_tmax_rast@data@values))-min(na.omit(reproj_tmax_rast@data@values)))
  stand_reproj_tmean_rast <- (reproj_tmean_rast-min(na.omit(reproj_tmean_rast@data@values)))/(max(na.omit(reproj_tmean_rast@data@values))-min(na.omit(reproj_tmean_rast@data@values)))
  stand_reproj_tmin_rast <- (reproj_tmin_rast-min(na.omit(reproj_tmin_rast@data@values)))/(max(na.omit(reproj_tmin_rast@data@values))-min(na.omit(reproj_tmin_rast@data@values)))
  stand_reproj_vpdmax_rast <- (reproj_vpdmax_rast-min(na.omit(reproj_vpdmax_rast@data@values)))/(max(na.omit(reproj_vpdmax_rast@data@values))-min(na.omit(reproj_vpdmax_rast@data@values)))
  stand_reproj_vpdmin_rast <-(reproj_vpdmin_rast-min(na.omit(reproj_vpdmin_rast@data@values)))/(max(na.omit(reproj_vpdmin_rast@data@values))-min(na.omit(reproj_vpdmin_rast@data@values)))
  
  ### Principal Component Analysis
  rasters_stand <- raster::stack(stand_reproj_ppt_rast, stand_reproj_tdmean_rast,stand_reproj_tmax_rast,stand_reproj_tmean_rast,stand_reproj_tmin_rast,stand_reproj_vpdmax_rast,stand_reproj_vpdmin_rast)
  
  # Spatial PCA
  pca1 <- RStoolbox::rasterPCA(rasters_stand)
  summary(pca1$model) # PCA components
  pca1$model$loadings # PCA loadings
  
  # Extract Bands from PCA
  pc1 <- pca1$map
  pc1_b1 <- pc1[[1]] # PC1
  pc1_b2 <- pc1[[2]] # PC2
  
  # Mask scaled rasters by study area (window)
  us <- raster::getData("GADM", country="USA", level=1)
  ca <- us[match(toupper("California"),toupper(us$NAME_1)),]
  # region union kills the data frame so don't overwrite 'wus'
  regs <- rgeos::gUnaryUnion(ca)
  # takes way too long to plot without simplifying the polygons
  regs <- rgeos::gSimplify(regs, 0.05, topologyPreserve = TRUE)
  # Add 0.1 degree buffer to capture all of PRISM
  ca_buffer <- rgeos::gBuffer(regs, width=0.1, byid=TRUE) # same projection as crs_us
  
  # Mask for California
  mask_pc1 <- mask(pc1_b1, ca_buffer)
  mask_pc2 <- mask(pc1_b2, ca_buffer)
  
  ca_pc1 <- crop(pc1_b1, ca_buffer)
  ca_pc1 <- mask(ca_pc1, ca_buffer)
  
  ca_pc2 <- crop(pc1_b2, ca_buffer)
  ca_pc2 <- mask(ca_pc2, ca_buffer)
  
  ca_pcs <- stack(ca_pc1, ca_pc2)
  return(ca_pcs)
    
}


get_rast <- function(files, type, year){
  
  files_type <- files[str_detect(files, type)]
  file_year <- files_type[str_detect(files_type, as.character(year))]
  return(raster(file_year))
  
}


for (y in years){
  print(y)
  pcs_y <- calc_prism_pcs(y)
  plot(pcs_y, main=y)
  writeRaster(pcs_y,
              file=paste("Documents/research/dataInt/data_time/prism_pcas_ca_", y, ".grd", sep=""),
              bandorder='BIL',
              overwrite=TRUE)
}
