library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


#######################
# save priors and other 
# values in output
#######################


outputs <- load_sim_outputs()


for (s in c("none", "medium", "high")){
  for (p in c("low", "medium", "high")){
    output <- get_output(outputs, s, p, "prefSampleGpCC")
    params <- load_output(paste('params_prefSampleGpCC_', s, '_', p, '.json', sep=""))

    output$L_w <- params$L
    output$L_ca <- params$L_ca
    output$L_co <- params$L_co
    output$L_a_ca <- params$L_a_ca 
    output$L_a_co <- params$L_a_co
    output$proposal.sd.theta <- params$proposal.sd.theta
    output$prior_phi <- params$prior_phi
    output$prior_theta <- params$prior_theta
    output$prior_alpha_ca_var <- 6
    output$prior_alpha_co_var <- 6
    output$n.sample <- params$n.sample
    output$burnin <- params$burnin    
    
    sim_name <- gen_sim_name(s, p)
    save_output(output, paste("output_", sim_name, ".json", sep=""))
    
  }
}


#######################
# save priors and other 
# values in output
#######################


outputs <- load_sim_outputs()


# save priors and other values in output
for (s in c("none", "medium", "high")){
  for (p in c("low", "medium", "high")){
    output <- get_output(outputs, s, p, "prefSampleGpCC")
    print(output$n.sample)
    print(output$description)
  }
}


#######################
# save priors and other 
# values in output for
# priors sensitivity 
# studies 
#######################

outputs <- load_sim_outputs_priorsens()
p <- "theta"
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, p, i)
  fname <- paste("params_iterate_prior_", p, "_", i, ".json", sep="")
  tuning_params <- load_output(fname)
  o$L_w <- tuning_params$L
  o$L_ca <- tuning_params$L_ca
  o$L_co <- tuning_params$L_co
  o$L_a_ca <- tuning_params$L_a_ca
  o$L_a_co <- tuning_params$L_a_co
  o$proposal.sd.theta <- tuning_params$proposal.sd.theta
  o$n.sample <- tuning_params$n.sample
  o$burnin <- tuning_params$burnin
  output_name <- paste("output_iterate_prior_", p, "_", i, ".json", sep="")
  print(output_name)
  save_output(o, output_name)

}
