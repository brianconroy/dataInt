
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=6)
N <- n_values(caPr.disc[[1]])
plot(caPr.disc)

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")

# - 'CA G Sq': CA Ground Squirrel
# - 'GM G Sq': Golden Mantle GS
# - 'Pine Squirrel': Douglas Squirrel
# - 'Chipmunk, YP': T. amoeunus (Yellow Pine)
# - 'Chipmunk, LP': T. speciousus (Lodgepole)??
# - 'Chipmunk, S': T. senex??
# - 'Chipmunk, M': T. merriami??
# - all but ds
# - CGS and 'Chipmunk, YP' combined

# all_species <- unique(rodents$Short_Name)
# species <- as.character(all_species[all_species != 'Pine Squirrel'])
# analysis_name <- 'analysis_all_but_ds'

species <- c('CA G Sq', 'Chipmunk, YP')
analysis_name <- gsub(',', '', gsub(' ', '_', paste('analysis', paste(species, collapse="_"), sep='_'), fixed=T))
rodents <- rodents[rodents$Short_Name %in% species,]

coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
loc.disc <- caPr.disc[[1]]
cells_all <- cellFromXY(loc.disc, coords_all)
cells_all <- cells_all[!is.na(cells_all[])]
cells_obs <- sort(unique(cells_all))

loc.disc[][!is.na(loc.disc[])] <- 0
loc.disc[][cells_obs] <- 5

plot(loc.disc)
plot(rasterToPolygons(loc.disc), add=T, border='black', lwd=1) 
points(coords_all, col='2')

# positive counts at each cell
rodents_pos <- rodents[rodents$Res == 'POS',]
coords_pos <- cbind(matrix(rodents_pos$Lon_Add_Fix), rodents_pos$Lat_Add_Fix)
cells_pos <- cellFromXY(loc.disc, coords_pos)
counts_pos <- data.frame(table(cells_pos))
names(counts_pos) <- c('cell', 'count_pos')

# negative counts at each cell
rodents_neg <- rodents[rodents$Res == 'NEG',]
coords_neg <- cbind(matrix(rodents_neg$Lon_Add_Fix), rodents_neg$Lat_Add_Fix)
cells_neg <- cellFromXY(loc.disc, coords_neg)
counts_neg <- data.frame(table(cells_neg))
names(counts_neg) <- c('cell', 'count_neg')

# combine counts
counts_all <- merge(counts_pos, counts_neg, by='cell', all=T)
counts_all$cell <- as.numeric(as.character(counts_all$cell))

counts_all[is.na(counts_all$count_pos),]$count_pos <- 0
if (sum(is.na(counts_all$count_neg)) > 0){
  counts_all[is.na(counts_all$count_neg),]$count_neg <- 0
}
counts_all <- counts_all[with(counts_all, order(cell)),]

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

# case data
cov.disc <- caPr.disc
x1 <- cov.disc[[1]][][locs$cells]
x2 <- cov.disc[[2]][][locs$cells]
x1.standardised <- (x1 - mean(x1))/sd(x1)
x2.standardised <- (x2 - mean(x2))/sd(x2)
x <- cbind(1, x1, x2)
x.standardised <- cbind(1, x1.standardised, x2.standardised)

case.data <- list(
  y=counts_all$count_pos,
  x.standardised=x.standardised,
  x=x,
  p=3
)

# control data
ctrl.data <- list(
  y=counts_all$count_neg,
  x.standardised=x.standardised,
  x=x,
  p=3
)

coords <- xyFromCell(caPr.disc, cell=all_ids)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
data <- list(loc=locs, case.data=case.data, ctrl.data=ctrl.data)

###########
# fit model
###########

# W initial value
prior_theta <- get_gamma_prior(prior_mean=7, prior_var=10)
prior_phi <- get_igamma_prior(prior_mean=15, prior_var=10)
w_output <- logisticGp(y=locs$status, d, n.sample=3000, burnin=1000, L=10,
                       theta_initial=prior_theta[1]*prior_theta[2],
                       phi_initial=prior_phi[2]/(prior_phi[1] - 1),
                       prior_phi=prior_phi, prior_theta=prior_theta,
                       proposal.sd.theta=0.15)

# optional extra burnin
w_output <- burnin_logisticGp_mcmc(w_output, n.burn=100)

w_output$accept
view_tr_w(w_output$samples.w)
view_tr(w_output$samples.theta)
view_tr(w_output$samples.phi)
hist(colMeans(w_output$samples.w))
save_output(w_output, paste("w_inival_output_cdph_", analysis_name, ".json", sep=""))

# w_output <- load_output(paste("w_inival_output_cdph_", analysis_name, ".json", sep=""))
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

prior_theta <- get_gamma_prior(prior_mean=theta_i, prior_var=10)
prior_phi <- get_igamma_prior(prior_mean=phi_i, prior_var=10)

ini_case <- glm(case.data$y ~ case.data$x.standardised + w_i[locs$ids] - 1, family='poisson')
alpha_ca_i <- coefficients(ini_case)[4]
beta_ca_i <- coefficients(ini_case)[1:3]

ini_ctrl <- glm(ctrl.data$y ~ ctrl.data$x.standardised + w_i[locs$ids] - 1, family='poisson')
alpha_co_i <- coefficients(ini_ctrl)[4]
beta_co_i <- coefficients(ini_ctrl)[1:3]

prior_alpha_ca_mean <- alpha_ca_i
prior_alpha_co_mean <- alpha_co_i
prior_alpha_ca_var <- 1
prior_alpha_co_var <- 1

n.sample <- 9000
burnin <- 1000
L_w <- 8
L_ca <- 8
L_co <- 8
L_a_ca <- 8
L_a_co <- 8
proposal.sd.theta <- 0.15

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

output <- prefSampleGpCC(data, d, n.sample, burnin,
                         L_w, L_ca, L_co, L_a_ca, L_a_co,
                         proposal.sd.theta=proposal.sd.theta,
                         m_aca=m_aca, m_aco=m_aco, m_ca=m_ca, m_co=m_co, m_w=m_w,
                         target_aca=0.65, target_aco=0.65, target_ca=0.65, target_co=0.65, target_w=0.65,
                         self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                         delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL,
                         beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                         theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                         prior_phi=prior_phi, prior_theta=prior_theta,
                         prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)

output <- burnin_after(output, n.burn=1000)

output <- continueMCMC(data, d, output, n.sample=3000)

plot(apply(output$samples.w, 1, mean), type='l')
view_tr_w(output$samples.w)

view_tr(output$samples.alpha.ca)
view_tr(output$samples.alpha.co)

par(mfrow=c(2,3))
view_tr(output$samples.beta.ca[,1])
view_tr(output$samples.beta.ca[,2])
view_tr(output$samples.beta.ca[,3])

view_tr(output$samples.beta.co[,1])
view_tr(output$samples.beta.co[,2])
view_tr(output$samples.beta.co[,3])

par(mfrow=c(1,1))
view_tr(output$samples.theta)
view_tr(output$samples.phi)

output$description <- paste("cdph", analysis_name, sep="_")
save_output(output, paste("cdph_", analysis_name,".json", sep=""))

###########
# downscale
###########

# interpolate w
w.hat <- colMeans(output$samples.w)
rw <- caPr.disc[[1]]
rw[][!is.na(rw[])] <- w.hat

xy <- data.frame(xyFromCell(rw, 1:ncell(rw)))
v <- getValues(rw)

tps <- Tps(xy, v)
p <- raster(caPr[[2]])
p <- interpolate(p, tps)
p <- mask(p, caPr[[1]])
w.hat_ds <- p[][!is.na(p[])]

par(mfrow=c(1,2))
plot(rw)
plot(p)

par(mfrow=c(1,2))
plot(rw)
plot(rw)
points(locs$coords, col=1, pch=16)

################
# calculate risk 
# surfaces
################

alpha.ca.hat <- mean(output$samples.alpha.ca)
alpha.co.hat <- mean(output$samples.alpha.co)
beta.ca.hat <- colMeans(output$samples.beta.ca)
beta.co.hat <- colMeans(output$samples.beta.co)

# low resolution
X_low <- load_x_ca2(factor=6)
lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * w.hat - X_low %*% beta.co.hat - alpha.co.hat * w.hat
risk_low <- calc_risk(lodds_low)

r_lodds_low <- caPr.disc[[1]]
r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
plot(r_lodds_low)

r_risk_low <- caPr.disc[[1]]
r_risk_low[][!is.na(r_risk_low[])] <- risk_low
plot(r_risk_low)

# high resolution
X_high <- load_x_ca2()
lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * w.hat_ds - X_high %*% beta.co.hat - alpha.co.hat * w.hat_ds
risk_high <- calc_risk(lodds_high)

r_lodds_high <- caPr[[2]]
r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
plot(r_lodds_high)

r_risk_high <- caPr[[2]]
r_risk_high[][!is.na(r_risk_high[])] <- risk_high

plot(r_risk_high)

r_cases <- rodents[rodents$Res == 'POS',]
r_coords <- cbind(r_cases$Lon_Add_Fix, r_cases$Lat_Add_Fix)
r_ctrls <- rodents[rodents$Res == 'NEG',]
r_coords_ctrl <- cbind(r_ctrls$Lon_Add_Fix, r_ctrls$Lat_Add_Fix)
points(r_coords_ctrl, col=rgb(0,0,1,0.05), pch=16, cex=0.65)
points(r_coords, col=rgb(1,0,0,0.2), pch=16, cex=0.65)

par(mfrow=c(1,2))
plot(r_lodds_low)
plot(r_lodds_high)
