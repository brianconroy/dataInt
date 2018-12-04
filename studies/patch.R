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

