#######################################
# Compare the multispecies preferential 
# sampling model against separate 
# analyses of the single species model. 

# iteration
#   goal: assess multispecies model performance under increasing
#   number of species.
#     does it work with more?
#     does it get better with more?
#   performance measures:
#     log odds accuracy
#     bias
#     w accuracy
#   computational efficiency:
#     number of samples drawn
#     burnin
#     runtime?
#   species numbers: c(2, 4, 6, 8)
# w: same
#######################################


library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')

n_species <- 8
sim <- "simMulti_opt1_iteration"
sim_name <- paste(sim, n_species, sep='_')

get_ca_data <- function(idx){
  y_idx <- data_store$case.data$y[[idx]]
  x_idx<- data_store$case.data$x.standardised[[idx]]
  return(list(
    y=y_idx,
    x.standardised=x_idx
  ))
}

get_co_data <- function(idx){
  y_idx <- data_store$ctrl.data$y[[idx]]
  x_idx<- data_store$ctrl.data$x.standardised[[idx]]
  return(list(
    y=y_idx,
    x.standardised=x_idx
  ))
}

get_locs <- function(idx){
  status_idx <- data_store$locs$status[[idx]]
  cells_idx <- data_store$locs$cells[[idx]]
  coords_idx <- data_store$locs$coords[[idx]]
  ids_idx <- data_store$locs$ids[[idx]]
  return(list(
    status=status_idx,
    cells=cells_idx,
    coords=coords_idx,
    ids=ids_idx
  ))
}

#### Prism Principal Components
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
n_values(caPr.disc[[1]])
plot(caPr.disc)

data_store <- load_output('simMulti_data_iteration.json')
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
N <- length(data_store$W)
W <- data_store$W
Theta <- data_store$theta
Phi <- data_store$phi
locs_species <- list()
case.data_species <- list()
ctrl.data_species <- list()
for (i in 1:n_species){
  locs_species[[i]] <- get_locs(i)
  case.data_species[[i]] <- get_ca_data(i)
  ctrl.data_species[[i]] <- get_co_data(i)
}

data <- list(
  locs=locs_species,
  case.data=case.data_species,
  ctrl.data=ctrl.data_species
)

#### Write simulation parameters to LaTeX
Alpha.case1 <- data$alpha.cases[1]
Alpha.ctrl1 <- data$alpha.ctrl[1]
beta.case1 <- data$beta.cases[1]
beta.ctrl1 <- data$beta.ctrls[1]
sim_config <- list(
  list(parameter="Alpha (case) 1", value=as.character(Alpha.case1)),
  list(parameter="Alpha (control) 1",value=as.character(Alpha.ctrl1)),
  list(parameter="Beta (case) 1", value=paste(beta.case1, collapse=", ")),
  list(parameter="Beta (control) 1", value=paste(beta.ctrl1, collapse=", ")),
  list(parameter="Range", value=as.character(Theta)),
  list(parameter="Phi", value=as.character(Phi))
)
write_latex_table(ldply(sim_config, "data.frame"), fname=paste("sim_params", ".txt", sep=""), path="/Users/brianconroy/Documents/research/project2/simulation_1_iteration")


####################
# Multispecies Model
####################

# calibrated initial values
# w_output <- logisticGp(y=locs1$status, d, n.sample=1000, burnin=200, L=10,
#                       prior_phi=prior_phi, prior_theta=prior_theta)
# view_logistic_output(w_output)
# save_output(w_output, "w_inival_output_multi.json")

w_output <- load_output("w_inival_output_multi.json")
w_i <- colMeans(w_output$samples.w)
theta_i <- mean(w_output$samples.theta)
phi_i <- mean(w_output$samples.phi)

alpha_ca_i <- list()
alpha_co_i <- list()
beta_ca_i <- list()
beta_co_i <- list()
N.d <- length(data$case.data)
for (k in 1:N.d){
  ini_case <- glm(data$case.data[[k]]$y ~ data$case.data[[k]]$x.standardised + w_i[data$locs[[k]]$ids] - 1, family='poisson')
  alpha_ca_i[[k]] <- unname(coefficients(ini_case)[4])
  beta_ca_i[[k]] <- unname(coefficients(ini_case)[1:3])
  
  ini_ctrl <- glm(data$ctrl.data[[k]]$y ~ data$ctrl.data[[k]]$x.standardised + w_i[data$locs[[k]]$ids] - 1, family='poisson')
  alpha_co_i[[k]] <- unname(coefficients(ini_ctrl)[4])
  beta_co_i[[k]] <- unname(coefficients(ini_ctrl)[1:3])
}

prior_alpha_ca_mean <- data_store$alpha.cases[1:n_species]
prior_alpha_co_mean <- data_store$alpha.ctrls[1:n_species]
prior_alpha_ca_var <- rep(4, n_species)
prior_alpha_co_var <- rep(4, n_species)
prior_theta <- c(6, 1)
prior_phi <- c(18, 204)

n.sample <- 4000
burnin <- 1000
L_w <- 8
L_ca <- rep(8, n_species)
L_co <- rep(8, n_species)
L_a_ca <- rep(8, n_species)
L_a_co <- rep(8, n_species)
proposal.sd.theta <- 0.15

m_aca <- 1000
m_aco <- 1000
m_ca <- 1000
m_co <- 1000
m_w <- 1000

target_aca <- 0.65
target_aco <- 0.65
target_ca <- 0.65
target_co <- 0.65
target_w <- 0.65

output <- prefSampleMulti_1(data, n.sample, burnin, 
                            L_w, L_ca, L_co, L_a_ca, L_a_co,
                            proposal.sd.theta=0.3,
                            m_aca=m_aca, m_aco=m_aca, m_ca=m_aca, m_co=m_aca, m_w=m_aca, 
                            target_aca=target_aca, target_aco=target_aco, target_ca=target_ca, target_co=target_co, target_w=target_w, 
                            self_tune_w=TRUE, self_tune_aca=TRUE, self_tune_aco=TRUE, self_tune_ca=TRUE, self_tune_co=TRUE,
                            delta_w=NULL, delta_aca=NULL, delta_aco=NULL, delta_ca=NULL, delta_co=NULL, 
                            beta_ca_initial=beta_ca_i, beta_co_initial=beta_co_i, alpha_ca_initial=alpha_ca_i, alpha_co_initial=alpha_co_i,
                            theta_initial=theta_i, phi_initial=phi_i, w_initial=w_i,
                            prior_phi=prior_phi, prior_theta=prior_theta, prior_alpha_ca_mean=prior_alpha_ca_mean,
                            prior_alpha_co_mean=prior_alpha_co_mean,
                            prior_alpha_ca_var=prior_alpha_ca_var, prior_alpha_co_var=prior_alpha_co_var)


# output
output <- burnin_after(output, n.burn=500)


# continue mcmc
output <- continueMCMC_multi(data, output, n.sample=500)


plot(apply(output$samples.w, 1, mean), type='l', col='2'); abline(h=mean(W), col='2')
w.hat <- colMeans(output$samples.w)
plot(x=W, y=w.hat); abline(0, 1, col=2)
summary(100*(W-w.hat)/W)
view_tr_w(output$samples.w, w_true=W)

view_tr(output$samples.theta, Theta)
print(mean(output$samples.theta)); print(Theta)

view_tr(output$samples.phi, Phi)
print(mean(output$samples.phi)); print(Phi)

par(mfrow=c(4,2))
for (i in 1:n_species){
  multispecies_plot(output$samples.alpha.ca, i, 'Alpha.case')
}
for (i in 1:n_species){
  multispecies_plot(output$samples.alpha.co, i, 'Alpha.ctrl')
}

par(mfrow=c(2,3))
multispecies_plot(output$samples.beta.ca, 1, 'beta.case', beta_id=1)
multispecies_plot(output$samples.beta.ca, 1, 'beta.case', beta_id=2)
multispecies_plot(output$samples.beta.ca, 1, 'beta.case', beta_id=3)

multispecies_plot(output$samples.beta.ca, 2, 'beta.case', beta_id=1)
multispecies_plot(output$samples.beta.ca, 2, 'beta.case', beta_id=2)
multispecies_plot(output$samples.beta.ca, 2, 'beta.case', beta_id=3)

multispecies_plot(output$samples.beta.co, 1, 'beta.ctrl', beta_id=1)
multispecies_plot(output$samples.beta.co, 1, 'beta.ctrl', beta_id=2)
multispecies_plot(output$samples.beta.co, 1, 'beta.ctrl', beta_id=3)

multispecies_plot(output$samples.beta.co, 2, 'beta.ctrl', beta_id=1)
multispecies_plot(output$samples.beta.co, 2, 'beta.ctrl', beta_id=2)
multispecies_plot(output$samples.beta.co, 2, 'beta.ctrl', beta_id=3)

output$description <- sim_name
save_output(output, paste("output_", sim_name, ".json", sep=""))
