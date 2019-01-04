library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')
ca <- getState('california')
caWin <- as.owin(ca)

#########
# coyotes
#########

src <- "/Users/brianconroy/Documents/research/cdph/data/"
coyotes <- read.csv(paste(src, "CDPH_coyote_recoded_full.csv", sep=""), header=T, sep=",")

names(coyotes)
xy <- matrix(coyotes$Long_QC)
xy <- cbind(xy, coyotes$Lat_QC)
plot(caWin)
points(xy, col='2')


par(mfrow=c(1,2))
plot(caWin, main='coyotes')
points(xy, col='2')

plot_rodents(rodents, 'rodents')


src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")

names(rodents)

# distribution of sampling sites by species
table(rodents$Short_Name)

rodents_cagsq <- rodents[rodents$Short_Name == 'CA G Sq',]
rodents_chips <- rodents[rodents$Short_Name == 'Chipmunk, S',]
rodents_chipy <- rodents[rodents$Short_Name == 'Chipmunk, YP',]
rodents_chipl <- rodents[rodents$Short_Name == 'Chipmunk, LP',]
rodents_gmgsq <- rodents[rodents$Short_Name == 'GM G Sq',]


plot_rodents <- function(df, title=''){
  xy <- matrix(df$Lon_Add_Fix)
  xy <- cbind(xy, df$Lat_Add_Fix)
  plot(caWin, main=title)
  points(xy, col='2')
}


par(mfrow=c(3,2))
plot_rodents(rodents_cagsq, 'CA G Sq')
plot_rodents(rodents_chips, 'Chipmunk, S')
plot_rodents(rodents_chipy, 'Chipmunk, YP')
plot_rodents(rodents_chipl, 'Chipmunk, LP')
plot_rodents(rodents_gmgsq, 'GM G Sq')



coords_neg <- matrix(rodents_neg$Lon_Add_Fix)
coords_neg <- cbind(coords_neg, rodents_neg$Lat_Add_Fix)
points(coords_neg)
coords_pos <- matrix(rodents_pos$Lon_Add_Fix)
coords_pos <- cbind(coords_pos, rodents_pos$Lat_Add_Fix)
points(coords_pos, col='2')




tr <- summary(rodents$Res)
print(tr)
print(tr[2]/(tr[1] + tr[2]))

ty <- table(as.character(rodents$Year))
plot(x=names(ty), y=as.numeric(ty), ylab='obs', xlab='year', type='l')
points(x=names(ty), y=as.numeric(ty), ylab='obs', xlab='year')

tc <- table(rodents$less_confident)
print(tc)
print(tc[2]/(tc[1] + tc[2]))

summary(rodents$Elevation)
hist(rodents$Elevation, main='', xlab='Elevation')

table(rodents$Agency_ID)
barplot(table(rodents$Agency_ID), names.arg='', ylab='# observations', xlab='agency')

table(rodents$Mamm_Group)

table(rodents$Sex_ID)

head(rodents$Latitude, n=100)
head(rodents$Lat_Add_Fix, n=100)
head(rodents$Longitude, n=100)

sum(is.na(rodents$Longitude))
sum(na.omit(abs(rodents$Longitude) > 10000))
sum(is.na(rodents$Lon_Add_Fix))
sum(rodents$Lon_Add_Fix == 0)

sum(is.na(rodents$Latitude))
sum(is.na(rodents$Lat_Add_Fix))
sum(rodents$Lat_Add_Fix == 0)

table(rodents$Specimen_Type_ID) # Collection method

table(rodents$gcode_confid)

library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')

ca <- getState('california')
caWin <- as.owin(ca)
plot(caWin)

rodents_pos <- rodents[rodents$Res == 'POS',]
rodents_neg <- rodents[rodents$Res == 'NEG',]

coords_neg <- matrix(rodents_neg$Lon_Add_Fix)
coords_neg <- cbind(coords_neg, rodents_neg$Lat_Add_Fix)
points(coords_neg)
coords_pos <- matrix(rodents_pos$Lon_Add_Fix)
coords_pos <- cbind(coords_pos, rodents_pos$Lat_Add_Fix)
points(coords_pos, col='2')


# how many collection events (distinct days) are there? 
# distinct date & location
collection_events <- unique(rodents[,c('Date', 'Location', 'Year')])
ty <- table(as.character(collection_events$Year))
plot(x=names(ty), y=as.numeric(ty), ylab='events', xlab='year', type='l', main='collection events')
points(x=names(ty), y=as.numeric(ty))

tl <- table(collection_events$Location)
barplot(tl, names.arg='', xlab='site', ylab='# surveys')
summary(as.numeric(tl))
sum(tl < 10)/nrow(collection_events)
sum(table(rodents$Location) < 2)/nrow(collection_events)


# rates over time
dfp <- data.frame(table(rodents_pos$Year))
names(dfp)[2] <- 'pos'
dfn <- data.frame(table(rodents_neg$Year))
names(dfn)[2] <- 'neg'
df <- merge(dfn, dfp, all=T)
df$rate <- df$pos/(df$neg + df$pos)
df$Var1 <- as.numeric(levels(df$Var1))[df$Var1]
plot(x=df$Var1, y=df$rate, xlab='year', ylab='disease rate', type='l')
points(x=df$Var1, y=df$rate, xlab='year', ylab='disease rate')
(max(df$rate) - min(df$rate))/min(df$rate)


# change in obs sites over time
par(mfrow=c(3,3))
for (i in 0:7){
  
  lb <- 1983 + i*4
  ub <- 1983 + (i+1)*4
  rodents_yr <- rodents[rodents$Year >= lb & rodents$Year < ub, ]

  plot(caWin, main = paste(lb, "to", ub))
  
  rodents_pos <- rodents_yr[rodents_yr$Res == 'POS',]
  rodents_neg <- rodents_yr[rodents_yr$Res == 'NEG',]
  
  coords_neg <- matrix(rodents_neg$Lon_Add_Fix)
  coords_neg <- cbind(coords_neg, rodents_neg$Lat_Add_Fix)
  points(coords_neg)
  coords_pos <- matrix(rodents_pos$Lon_Add_Fix)
  coords_pos <- cbind(coords_pos, rodents_pos$Lat_Add_Fix)
  points(coords_pos, col='2')
  
}


# summary of recent data
rodents_00 <- rodents[rodents$Year >= 2000,]
tr <- table(rodents_00$Res)
tr[2]/(tr[1] + tr[2])

rodents_pos <- rodents_00[rodents_00$Res == 'POS',]
rodents_neg <- rodents_00[rodents_00$Res == 'NEG',]

plot(caWin)
coords_neg <- matrix(rodents_neg$Lon_Add_Fix)
coords_neg <- cbind(coords_neg, rodents_neg$Lat_Add_Fix)
points(coords_neg)
coords_pos <- matrix(rodents_pos$Lon_Add_Fix)
coords_pos <- cbind(coords_pos, rodents_pos$Lat_Add_Fix)
points(coords_pos, col='2')

barplot(table(rodents_00$Agency_ID), names.arg='', ylab='# observations', xlab='agency')


####################
# Case only PS model
####################
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


wc <- readWC()
ca <- getState('california')
caWin <- as.owin(ca)
caWc <- getStateWC(ca, wc)


#### Discretize the study region
simRegion <- discretizeSimRegion(caWin, caWc, factor=3)
caWc.disc <- simRegion$raster

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")

rodents_pos <- rodents[rodents$Res == 'POS',]
coords_pos <- matrix(rodents_pos$Lon_Add_Fix)
coords_pos <- cbind(coords_pos, rodents_pos$Lat_Add_Fix)

loc.disc <- caWc.disc[[1]]
cells_pos <- cellFromXY(loc.disc, coords_pos)
cells_pos <- cells_pos[!is.na(cells_pos[])]
cells_obs <- unique(cells_pos)
cells_obs <- sort(cells_obs)

loc.disc[][!is.na(loc.disc[])] <- 0
loc.disc[][cells_obs] <- 5

plot(loc.disc)
plot(rasterToPolygons(loc.disc), add=T, border='black', lwd=1) 
points(coords_pos, col='2')


# counts at observed cells
counts_pos <- data.frame(table(cells_pos))
names(counts_pos) <- c('cell', 'count_pos')
counts_df <- counts_pos
counts_df <- counts_df[with(counts_df, order(cell)),]
counts_df$cell <- as.numeric(as.character(counts_df$cell))
points(coords_pos, col='2')


# location data
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]

locs <- list(
  cells=cells_obs,
  status=1 * c(all_ids %in% cells_obs),  
  coords=xyFromCell(loc.disc, cells_obs)
)
locs$ids <- c(1:length(all_ids))[as.logical(locs$status)]

plot(loc.disc)
points(locs$coords)

# disease data
cov.disc <- caWc.disc[[c(1)]] # Annual mean temperature

x <- cov.disc[][locs$cells]
x.standardised <- (x - mean(x))/sd(x)
x <- cbind(1, x)
x.standardised <- cbind(1, x.standardised)
y <- counts_df$count_pos

count.data <- list(
  y=y,
  x.standardised=x.standardised,
  x=x,
  p=2
)

# fit model
coords <- xyFromCell(caWc.disc, cell=all_ids)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
data <- list(loc=locs, conditional=count.data)
prior_alpha_mean <- 2
prior_alpha_var <- 6


output <- prefSampleGp(data, n.sample=20000, burnin=2000, L_w=25, L_c=23, L_a=20, proposal.sd.theta=0.3, #target_a=0.95,
                       target_w=0.9, target_c=0.75, m_c=700, beta_initial=c(1.5, 1), self_tune_a=F, delta_a=0.01,
                       alpha_initial=c(1))
output <- prefSampleGp(data, n.sample=30000, burnin=2000, L_w=25, L_c=23, L_a=20, proposal.sd.theta=0.3, target_a=0.95,
                       target_w=0.9, target_c=0.75, m_c=700, beta_initial=c(1.5, 1),
                       alpha_initial=c(1))
plot(output$deltas_w)
plot(output$deltas_a)
plot(output$deltas_c)
print(output$accept)

plot(apply(output$samples.w, 1, mean), type='l', col='2', ylab='mean w')
w.hat <- colMeans(output$samples.w)
hist(w.hat)

par(mfrow=c(4, 4))
j <- 1
for(i in 1:16){
  plot(output$samples.w[,j*i], type='l')
}

par(mfrow=c(1,2))
plot(output$samples.theta, type='l', ylab='Theta')
hist(output$samples.theta, xlab='Theta', main='')

plot(output$samples.phi, type='l', ylab='Phi')
hist(output$samples.phi, xlab='Phi', main='')

plot(output$samples.beta.c[,1], type='l', ylab='beta0')
plot(output$samples.beta.c[,2], type='l', ylab='beta1')

plot(output$samples.alpha, type='l', xlab='iteration', ylab='Alpha')
hist(output$samples.alpha, xlab='Alpha')
