############################
# develop area to point kriging
###############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


Theta <- 6
Phi <- 12


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=6)
n_values(caPr.disc[[1]])
plot(caPr.disc)


#### Simulate gaussian process at fine scale
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
N <- length(W)
hist(W)
W_r_fine <- caPr.disc[[1]]
W_r_fine[][!is.na(W_r_fine[])] <- W
plot(W_r_fine)


#### Aggregate the Gaussian process to a coarse scale
#   make a coarse scale raster
#   extract points of the fine scale one at the coarse scale centroids
caPr.disc_coarse <- aggregate(caPr, fact=8)
cells_coarse <- c(1:ncell(caPr.disc_coarse))[!is.na(values(caPr.disc_coarse[[1]]))]
coords_coarse <- xyFromCell(caPr.disc_coarse[[1]], cell=cells_coarse)
plot(caPr.disc_coarse[[1]])
points(coords_coarse)

plot(W_r_fine)
points(coords_coarse)
W_coarse <- extract(W_r_fine, coords_coarse)
W_r_coarse <- caPr.disc_coarse[[1]]
W_r_coarse[][!is.na(W_r_coarse[])] <- W_coarse

par(mfrow=c(1,2))
plot(W_r_fine)
plot(W_r_coarse)


#### Now Krige
# block to block covariance

# block to point covariance

# given a coarse cell id
#   get cell ids of fine grained points inside it

# given a coords of a coarse cell centroid
#   get coords of fine grained cell centroids which fall inside it

coarse_poly <- rasterToPolygons(W_r_coarse)
plot(W_r_fine)
plot(coarse_poly, add=T)

plot(W_r_fine)
plot(rasterToPolygons(disaggregate(W_r_coarse, fact=2)), add=T)

plot(coarse_poly)
plot(rasterToPolygons(disaggregate(W_r_coarse, fact=2)), add=T)

plot(coarse_poly)
points(coordinates(rasterToPolygons(disaggregate(W_r_coarse, fact=2))))

finer_poly <- rasterToPolygons(disaggregate(W_r_coarse, fact=2))
nrgk = extract(W_r_fine, coordinates(finer_poly))
blyat = W_r_fine



# disaggregate
template_finer <- disaggregate(W_r_coarse, fact=2)
plot(template_finer)

# extract fine scale raster values at spatial polygon cell centroids
cells.all <- c(1:ncell(template_finer))[!is.na(values(template_finer[[1]]))]
coords_f <- xyFromCell(template_finer, cell=cells.all)
vals <- extract(W_r_fine, coords_f)

plot(W_r_coarse)
points(coords_f)

# map those onto the disaggregated raster
template_finer[][!is.na(template_finer[])] <- vals
plot(template_finer)

par(mfrow=c(1,2))
plot(W_r_coarse)
plot(template_finer)


# why not just start fine then aggregate?


##########
# tapering?
##########

cells.all <- c(1:ncell(caPr[[1]]))[!is.na(values(caPr[[1]]))]
coords <- xyFromCell(caPr[[1]], cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))

Sigma <- Exponential(d, range=Theta, phi=Phi)
Sigma_tap <- Sigma * Wendland(d, theta=6, k=1, dimension=2)
plot(as.im(Sigma_tap))
nnzero(Sigma_tap)
nnzero(Sigma)
nrgk=solve(Sigma_tap)
nrgk=solve(Sigma)


#########
# splines
#########
library(fields)

r <- caPr[[1]]
ra <- aggregate(r, 10)

xy <- data.frame(xyFromCell(ra, 1:ncell(ra)))
v <- getValues(ra)

# Thin plate spline model
tps <- Tps(xy, v)
p <- raster(r)

# use model to predict values at all locations
p <- interpolate(p, tps)
p <- mask(p, r)
par(mfrow=c(1,2))
plot(r)
plot(p)
par(mfrow=c(1,1))
plot(x=getValues(r), y=getValues(p)); abline(0, 1, col=2)

# add your covariate values to the interpolation step
elevation <- (init(r, 'x') * init(r, 'y')) / 100000000
names(elevation) <- 'elev'
# elevation <- mask(elevation, r)
z <- extract(elevation, xy)
z[is.na(z)] <- mean(z, na.rm=T)

# add as another independent variable
xyz <- cbind(xy, z)
tps2 <- Tps(xyz, v)
p2 <- interpolate(elevation, tps2, xyOnly=FALSE)
p2 <- mask(p2, r)
par(mfrow=c(1,3))
plot(r)
plot(p)
plot(p2)

par(mfrow=c(1,2))
plot(x=getValues(r), y=getValues(p)); abline(0, 1, col=2)
plot(x=getValues(r), y=getValues(p2)); abline(0, 1, col=2)


###################
# simulate gaussian 
# field w
###################


#### Prism Principal Components
caPr <- load_prism_pcs()
caPr <- aggregate(caPr[[1]], 3)

#### Simulate gaussian process at fine scale
cells.all <- c(1:ncell(caPr))[!is.na(values(caPr[[1]]))]
coords <- xyFromCell(caPr, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
Sigma <- Exponential(d, range=Theta, phi=Phi)
set.seed(40)
W <- mvrnorm(n=1, mu=rep(0, length(cells.all)), Sigma)
hist(W)

rw <- caPr
rw[][!is.na(rw[])] <- W
ra <- aggregate(rw, 3)

par(mfrow=c(1,2))
plot(rw)
plot(ra)

xy <- data.frame(xyFromCell(ra, 1:ncell(ra)))
v <- getValues(ra)

# Thin plate spline model
tps <- Tps(xy, v)
p <- raster(rw)

p <- interpolate(p, tps)
p <- mask(p, rw)
par(mfrow=c(1,2))
plot(rw)
plot(p)
par(mfrow=c(1,1))
plot(x=getValues(rw), y=getValues(p)); abline(0, 1, col=2)
