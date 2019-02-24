############################
# Summarizes the preferential
# sampling x prevalence
# simulation results
############################

library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='prefSampleGpCC')
s <- "high"
dst <- "/Users/brianconroy/Documents/research/project1/simulations_iteration/"


# summary of simulation parameters
params_all <- list()
counter <- 1
for (s in c('none', 'medium', 'high')){
  for (p in c('low', 'medium', 'high')){
    params <- load_sim_params(s, p)
    params <- params[, !names(params) %in% c('total.y.ca', 'total.y.co')]
    params$beta.case <- paste('(', gsub(' ', ', ', params$beta.case), ')', sep='')
    params$beta.ctrl <- paste('(', gsub(' ', ', ', params$beta.ctrl), ')', sep='')
    params_all[[counter]] <- params
    counter <- counter + 1
  }
}
write_latex_table(ldply(params_all, 'data.frame'), 'iteration_sim_params.txt', dst)


# log odds
  # scatterplots
    # for each sampling level export a 3x3 figure
    # iterated over the different prevalences
par(mfrow=c(3,3))
prevs <- c("low", "medium", "high")
mains <- c("A)", "B)", "C)")
xl <- c(-44, 17)
yl <- xl
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


for (s in c("none", "medium", "high")){
  for (p in c("low", "medium", "high")){
    fname <- paste("traces_", s, "_", p, ".png", sep="")
    png(paste(dst, fname, sep=""),
        width=900, height=700, res=100)
    plot_traces(outputs, s, p)
    dev.off()
  }
}


for (s in c("none", "medium", "high")){
  tables <- list()
  counter <- 1
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
}


# w scatterplots
par(mfrow=c(3,3))
labels <- c("A)", "", "", "B)", "", "", "C)", "", "")
counter <- 1
for (s in c("none", "medium", "high")){
  for (p in c("low", "medium", "high")){
    lab <- labels[counter]
    output <- get_output(outputs, s, p, "prefSampleGpCC")
    truevals <- load_params(paste('true_params_', s, '_', p, '.json', sep=''))
    w.hat <- colMeans(output$samples.w)
    plot(x=truevals$W, y=w.hat, main=lab, xlab="True W", ylab="Estimated W"); abline(0, 1, col='2')
    counter <- counter + 1
  }
}
  

# log odds mse table
rows <- list()
counter <- 1
for (s in c("none", "medium", "high")){
  print(s)
  for (p in c("low", "medium", "high")){
    print(p)
    lodds <- calc_log_odds(outputs, s, p)
    lodds_true <- calc_log_odds_true(s, p)
    lodds_sp <- calc_log_odds_sp(outputs, s, p)
    lodds_pr <- calc_log_odds_pr(s, p)
    new_row <- list(
      sampling=s,
      prevalence=p,
      PS=round(sqrt(mean((lodds-lodds_true)^2)), 2),
      SP=round(sqrt(mean((lodds_sp-lodds_true)^2)), 2),
      PR=round(sqrt(mean((lodds_pr-lodds_true)^2)), 2)
    )
    rows[[counter]] <- new_row
    counter <- counter + 1
  }
}
df_rmse <- ldply(rows, 'data.frame')
write_latex_table(df_rmse, "latex_rmse.txt")


# log odds barplot
reformat <- function(sampling){
  df_s <- df_rmse[df_rmse$sampling == sampling,]
  rows_new <- list()
  i <- 1
  for (p in df_s$prevalence){
    r_p <- df_s[df_s$prevalence==p,]
    rows_new[[i]] <- list(
      prevalence=p,
      model='PS',
      rmse=r_p$PS
    )
    rows_new[[i+1]] <- list(
      prevalence=p,
      model='SP',
      rmse=r_p$SP
    )
    rows_new[[i+2]] <- list(
      prevalence=p,
      model='PR',
      rmse=r_p$PR
    )
    i <- i + 3
  }
  return(ldply(rows_new, 'data.frame'))
}


rmse_none <- reformat('none')
rmse_medium <- reformat('medium')
rmse_high <- reformat('high')

p1 <- ggplot(rmse_none,aes(x=prevalence,y=rmse,fill=factor(model)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(name="Model", values = c("PS" = "blue", "SP" = "orange", "PR" = "red"))+
  xlab("Prevalence")+ylab("RMSE")+
  ylim(0, 13)+
  ggtitle('A)')

p2 <- ggplot(rmse_medium,aes(x=prevalence,y=rmse,fill=factor(model)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(name="Model", values = c("PS" = "blue", "SP" = "orange", "PR" = "red"))+
  xlab("Prevalence")+ylab("RMSE")+
  ylim(0, 13)+
  ggtitle('B)')

p3 <- ggplot(rmse_high,aes(x=prevalence,y=rmse,fill=factor(model)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(name="Model", values = c("PS" = "blue", "SP" = "orange", "PR" = "red"))+
  xlab("Prevalence")+ylab("RMSE")+
  ylim(0, 13)+
  ggtitle('C)')

grid.arrange(p1, p2, p3)


# export tuning parameters to LaTeX
rows_all <- c()
for (s in c("none", "medium", "high")){
  for (p in c("low", "medium", "high")){
    o <- get_output(outputs, s, p, "prefSampleGpCC")
    rows_o <- summarize_mcmc_pscc(o, s)
    rows_all <- c(rows_o, rows_all)
  }
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_iteration_mcmc.txt", path=dst)


# parameter biases
for (s in c("none", "medium", "high")){
  biases <- table_params_wide(outputs, s)
  write_latex_table(biases, paste("latex_iteration_bias_", s, ".txt", sep=""), path=dst)
}


# model 1 specific parameter biases
for (s in c("none", "medium", "high")){
  rows <- list()
  for (p in c("low", "medium", "high")){
    output_ps <- get_output(outputs, s, p, 'prefSampleGpCC')
    true_params <- load_params(paste('true_params_', s, '_', p, '.json', sep=''))
    params <- c("Alpha (case)", "Alpha (control)", "Theta", "Phi")
    for (param in params){
      rows[[counter]] <- make_row(
        p,
        s,
        "PS",
        param,
        true_params,
        output_ps
      )
      counter <- counter + 1
    }
  }
  rows_df <- ldply(rows, 'data.frame')
  rows_df <- rows_df[,!colnames(rows_df) %in% c('model', 'true')]
  write_latex_table(rows_df, paste("latex_iteration_bias_model1_", s, ".txt", sep=""), path=dst)
}

