#############################
# Choose simulation parameter
# values
#############################

library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


# Alpha variances
prior_alpha_var <- c(6, 12, 18, 24, 30)
write(toJSON(prior_alpha_var), 
      file='/Users/brianconroy/Documents/research/dataInt/output/simParams_alpha_prior.txt')


# target variances
# 4, 9, 12, 20, 40
iterate_ig_variance(Phi)
priors_phi <- list(
  c(32, 372),
  c(18, 204),
  c(14, 156),
  c(9.167, 98),
  c(5.58, 55)
)
write(toJSON(priors_phi), 
      file='/Users/brianconroy/Documents/research/dataInt/output/simParams_phi_prior.txt')

# target theta variances
# 3, 6, 12, 24, 48
iterate_g_variance(Theta, 0.25, 50)
priors_theta <- list(
  c(12, 0.5),
  c(6, 1),
  c(3, 2),
  c(1.5, 4),
  c(0.75, 8)
)
write(toJSON(priors_theta), 
      file='/Users/brianconroy/Documents/research/dataInt/output/simParams_theta_prior.txt')
