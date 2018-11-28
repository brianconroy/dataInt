############################
# Summarizes the preferential
# sampling x prevalence
# simulation results
############################

library(plyr)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs()


# log odds
  # scatterplots
    # for each sampling level export a 3x3 figure
    # iterated over the different prevalences
par(mfrow=c(3,3))
s <- "none"
prevs <- c("low", "medium", "high")
mains <- c("A)", "B)", "C)")
xl <- c(-8, 5)
yl <- c(-8, 5)
for (i in 1:3){
  p <- prevs[i]
  lodds <- calc_log_odds(outputs, s, p)
  lodds_true <- calc_log_odds_true(s, p)
  lodds_sp <- calc_log_odds_sp(outputs, s, p)
  lodds_pr <- calc_log_odds_pr(s, p)
  
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', main=mains[i], xlim=xl, ylim=yl); abline(0, 1, col=2)
  plot(x=lodds_true, y=lodds_sp, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)
  plot(x=lodds_true, y=lodds_pr, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)
}


s <- "none"
for (p in c("low", "medium", "high")){
  fname <- paste("traces_", s, "_", p, ".png", sep="")
  png(paste("/Users/brianconroy/Documents/research/project1/simulations/", fname, sep=""),
      width=900, height=700, res=100)
  plot_traces(outputs, s, p)
  dev.off()
}


tables <- list()
counter <- 1
s <- "none"
for (p in c("low", "medium", "high")){
  tables[[counter]] <- table_params(outputs, s, p)
  counter <- counter + 1
}
for (i in 1:length(tables)){
  df <- ldply(tables[[i]], 'data.frame')
  s <- df$sampling[1]
  p <- df$prevalence[1]
  fname <- paste("latex_param_estimates_", s, "_", p, ".txt", sep="")
  write_latex_table(df, fname)
}


# w scatterplots
  # wait till you fit other models

# log odds
  # scatterplots
  # box & whisker plots
  # maps

# param traceplots
# param tables
# plot of param biases
