library(prism)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')
options(prism.path="/Users/brianconroy/Documents/research/dataInt/data")

# 30 year tmean normals for July at a 4km grid resolution
get_prism_normals(type = 'tmean', resolution = '4km', mon = 12, keepZip = TRUE)
RS <- prism_stack(ls_prism_data()[18,1])
plot(RS)

# precipitation for the month of July from 1985 to 2016
get_prism_monthlys(type = 'ppt', years=1990:2016, mon = 7, keepZip = TRUE)
ls_prism_data(name=TRUE)
RS <- prism_stack(ls_prism_data()[1:10,1])
plot(RS[[1]])
par(mfrow=c(2,2))
for (i in 1:4){
  plot(RS[[i]])
}

# filter data by state
ca <- getState('california')
caWin <- as.owin(ca)
sr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
RS <- projectRaster(RS, crs=sr)
caPr <- getStateWC(ca, RS)
plot(caPr[[1]])
