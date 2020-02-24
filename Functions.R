
######## Making coefficient plots from rstanarm objects 

names_rstanarm <- function(stan_object){
  draws_df <- as.data.frame(stan_object)
  names <- colnames(draws_df)
  return(names)
}


PPlot_COEF <- function(stan_object, par_select,new_par_names,title_string) {
  require(ggplot2)
  draws_df <- as.data.frame(stan_object)
  draws_list <- lapply(draws_df, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))
  lows95 <- rep(0,length(draws_list))
  lows50 <- rep(0,length(draws_list))
  his95 <- rep(0,length(draws_list))
  his50 <- rep(0,length(draws_list))
  meds <- rep(0, length(draws_list))
  params <- rep(NA,length(draws_list))
  for(i in 1:length(draws_list)){
    params[i] <- paste(colnames(draws_df[i]))
    lows95[i] <- draws_list[[i]][1]
    lows50[i] <- draws_list[[i]][2]
    meds[i] <- draws_list[[i]][3]
    his50[i] <- draws_list[[i]][4]
    his95[i] <- draws_list[[i]][5]
  }
  plot_df <- data.frame(par = as.factor(params), lo95 = lows95, lo50 = lows50, med = meds, hi50 = his50, hi95 = his95)
  plot_df_real <- plot_df[which(plot_df$par%in%par_select),]
  plot_df_real$pars <- factor(plot_df_real$par, levels = par_select)
  levels(plot_df_real$pars) <- new_par_names
  coef_plot <- ggplot(plot_df_real, aes(x = pars)) + geom_point(aes(y=med),size = 1.5) +
    geom_segment(aes(xend = pars, y = lo95, yend = hi95)) +
    geom_segment(aes(xend=pars,y=lo50,yend=hi50),size = 1.1) +
    geom_hline(aes(yintercept = 0),linetype = "dashed") +
    coord_flip() + ylab("Estimates") + xlab("") +
    scale_x_discrete(limits = rev(levels(plot_df_real$pars))) +
    ggtitle(paste(title_string)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"));
  
  print(coef_plot)
  return(coef_plot)
}


# automated pairwise comparisons 
# This takes advantage of the fact that MCMC returns joint density over all 
# parameters, and thus we can estimate differences seamlessly. 
# Use case: coefficients grouped into batches (say N responses by cultivar), and
# we want to compare all the cultivars to each other...etc. 

pairwise_function <- function(stan_model, par_select, new_par_names){
  
  draw_df2 <- as.data.frame(stan_model)
  
  diff_array <- array(0, dim = c(M,N,N))
  
  for(i in 1:N){
    for(j in 1:N){
      diff_array[,j,i] <- draw_df2[,i] - draw_df2[,j]
    }
  }
  
  nums2select <- which(colnames(draw_df2) %in% par_select)
  
  out_array <- diff_array[,nums2select,nums2select]
  dimnames(out_array)[[2]] <- new_par_names
  dimnames(out_array)[[3]] <- new_par_names
  
  mean_diff <- apply(out_array, c(2,3),mean)
  se_diff <- apply(out_array, c(2,3),sd)
  low95_diff <- apply(out_array, c(2,3),quantile, probs = c(0.025))
  hi95_diff <- apply(out_array,c(2,3),quantile, probs = c(0.975))
  
  diff_names <- matrix(0,nrow(mean_diff),ncol(mean_diff))
  
  for(i in 1:nrow(mean_diff)){
    for(j in 1:ncol(mean_diff)){
      diff_names[i,j] <- paste0(paste(rownames(mean_diff)[j]), paste("-"), paste(colnames(mean_diff)[i]))
    }
  }
  
  diff_names2 <- diff_names[upper.tri(mean_diff)]
  mean_diffs <- mean_diff[upper.tri(mean_diff)]
  low95_diffs <- low95_diff[upper.tri(mean_diff)]
  hi95_diffs <- hi95_diff[upper.tri(mean_diff)]
  
  diff_df <- data.frame(par_diff = as.factor(diff_names2),mean_diffs = mean_diffs, lo95 = low95_diffs, hi95 = hi95_diffs)
  
  coef_plot <- ggplot(diff_df, aes(x = par_diff)) + geom_point(aes(y=mean_diffs),size = 1.5) + geom_segment(aes(xend = par_diff, y = lo95, yend = hi95)) +
    geom_hline(aes(yintercept = 0),linetype = "dashed") + coord_flip() + ylab("Estimates") + xlab("")
  
  print(coef_plot)
  return(coef_plot)
  
}


