
#################################################################################
# functions for fitting the extended Elo rating model ###########################
#################################################################################

# function to run the extended Elo rating method for a single group and either return the log-likelihood or Elo scores and k values for each individual in an interaction
elo_one_group <- function(par, elo_data, burn_in = 100, select_par = c(1,0,0,0,0,0,0,0,0,0,0), return_elo_scores=F)
{    
  # initialize parameters; selected parameters are initialized with the set value, all other parameters are set to zero 
  k_par <- rep(0,11)
  k_par[select_par==1] <- par

  # get all individuals
  all_ids <- unique(c(as.character(elo_data$Winner),as.character(elo_data$Loser)))

  # current Elo scores
  currentELO <- rep(0,length(all_ids)); names(currentELO) <- all_ids
  
  # the date of the most recent interaction, needed for recentering Elo scores
  prev_date <- rep(0,length(all_ids)); names(prev_date) <- all_ids
  
  # initialize columns to save Elo scores and k values
  if (return_elo_scores==T) 
  {
    elo_data$elo_w <- NA  #Elo score of the winner
    elo_data$elo_l <- NA  #Elo score of the loser
    elo_data$k_w <- NA  #k value of the winner
    elo_data$k_l <- NA  #k value of the loser
  }
  
  L <- 0  #log-likelihood
  
  for(i in 1:nrow(elo_data)) 
  {       
    ind1 <- which(names(currentELO)==elo_data$Winner[i])  #Winner in this interaction
    ind2 <- which(names(currentELO)==elo_data$Loser[i])   #Loser in this interaction
    
    p_win <- 1/(1+exp(-.01*(currentELO[ind1] - currentELO[ind2]))) # winning chance of the winner
    
    if (i <= burn_in)   #during burn-in period all k values are fixed to 100
    {
      currentELO[ind1] <- currentELO[ind1] + 100 * (1 - p_win)  #new Elo score of the Winner
      currentELO[ind2] <- currentELO[ind2] - 100 * (1 - p_win)  #new Elo score of the Loser
    }
    else  #after the burn-in period k values are a function of selected predictor variables
    {
      k1 <- exp(sum(k_par[1:2]) + sum(k_par[3:4]) * elo_data$hybrid_w[i] + sum(k_par[5:6]) * elo_data$Age_w[i] + sum(k_par[7:8]) * elo_data$Lag_day_w[i] + sum(k_par[9:10]) * elo_data$aggr_index[i] + k_par[11] * elo_data$hybrid_w[i] * elo_data$aggr_index[i])
      k2 <- exp(k_par[1]        +        k_par[3] * elo_data$hybrid_l[i] +        k_par[5] * elo_data$Age_l[i] +        k_par[7] * elo_data$Lag_day_l[i] +         k_par[9] * elo_data$aggr_index[i] + k_par[11] * elo_data$hybrid_l[i] * elo_data$aggr_index[i])

      currentELO[ind1] <- currentELO[ind1] + k1 * (1 - p_win) #new Elo score of the Winner
      currentELO[ind2] <- currentELO[ind2] - k2 * (1 - p_win) #new Elo score of the Loser

      # save Elo scores and k values for output
      if (return_elo_scores==T) 
      {
        elo_data$elo_w[i] <- currentELO[ind1]
        elo_data$elo_l[i] <- currentELO[ind2]
        elo_data$k_w[i] <- k1
        elo_data$k_l[i] <- k2  
      }
      else L <- L + log(p_win)  #update log-likelihood 
    }
    
    #recenter Elo scores around mean 0, consider only individuals that interacted within the past 90 days
    currentELO[elo_data$date[i] - prev_date <= 90] <- currentELO[elo_data$date[i] - prev_date <= 90] - mean(currentELO[elo_data$date[i] - prev_date <= 90])

    #update date of last interaction
    prev_date[ind1] <- prev_date[ind2] <- elo_data$Date[i]
    
  }  
  if (return_elo_scores==T) return(elo_data)    #return data frame with added Elo scores and k values
  else return(-1*L)   #return negative log-likelihood
}




#################################################################################

# function to combine results of multiple group 
elo_multiple_groups <- function(par, ago_data, select_par, return_elo_scores)
{
  groups <- unique(ago_data$Group)  # get all group IDs
  
  if (return_elo_scores) #return Elo scores and k values
  {
    new_data <- data.frame()
    for (i in groups) new_data <- rbind(new_data, elo_one_group(par, subset(ago_data, Group == i), burn_in=100, select_par=select_par, return_elo_scores=return_elo_scores))
    return(new_data)
  }
  else  # return negative log-likelihood
  {
    log_lik <- 0
    for (i in groups) log_lik <- log_lik + elo_one_group(par, subset(ago_data, Group == i), burn_in=100, select_par=select_par, return_elo_scores=F)
    return(log_lik)
  }
}


elo_multiple_groups(c(4.6), ago_data = ago_data, select_par = c(1,0,0,0,0,0,0,0,0,0,0), return_elo_scores=F)
res <- elo_multiple_groups(c(4.6,0,0,0,0,0,0,0,0,0,0), ago_data = ago_data, select_par = c(1,1,1,1,1,1,1,1,1,1,1), return_elo_scores=T)




#################################################################################
# read data file  ###############################################################
#################################################################################


# please insert the correct path where the file 'ago_data.csv' was saved
ago_data <- read.csv('C:/../ago_data.csv', sep=';')



#################################################################################
# model fitting  ################################################################
#################################################################################


# fitting the full model
res_fit_full <- optim(c(4.6,0,0,0,0,0,0,0,0,0,0), elo_multiple_groups, ago_data = ago_data, select_par = c(1,1,1,1,1,1,1,1,1,1,1), return_elo_scores=F, method='BFGS')

# fitting the reduced model
res_fit_red <- optim(c(4.6,0,0,0,0,0,0), elo_multiple_groups, ago_data = ago_data, select_par = c(1,1,1,0,1,1,1,0,1,0,0), return_elo_scores=F, method='BFGS')




#################################################################################
#plot k values ##################################################################
#################################################################################


#get k values and Elo scores
res_data <- elo_multiple_groups(res_fit_red$par, ago_data = ago_data, select_par = c(1,1,1,0,1,1,1,0,1,0,0), return_elo_scores=T)
res_data <- subset(res_data, !is.na(elo_w)) #remove burn-in data points


#reformat data for plotting
res_data_all <- data.frame('winner' = c(rep('Winner', nrow(res_data)), rep('Loser', nrow(res_data))),  'k' = c(res_data$k_w, res_data$k_l) , 'hybrid' = c(res_data$hybrid_w, res_data$hybrid_l), 'age' = c(res_data$Age_w, res_data$Age_l), 'lag' = c(res_data$Lag_day_w, res_data$Lag_day_l), 'agg' = c(res_data$aggr_index, res_data$aggr_index), 'elo' = c(res_data$elo_w, res_data$elo_l), 'group' =  c(res_data$Group, res_data$Group), 'date' =  c(res_data$Date, res_data$Date))

# create categories for plotting
res_data_all$a <- ''
res_data_all$a[res_data_all$agg< mean(res_data_all$agg)] <- 'Low aggresssion intensity'
res_data_all$a[res_data_all$agg>=mean(res_data_all$agg)] <- 'High aggresssion intensity'
res_data_all$a <- factor(res_data_all$a, levels=c('Low aggresssion intensity','High aggresssion intensity'))
res_data_all$winner <- factor(res_data_all$winner, levels=c('Winner','Loser'))

#randomize rows (to avoid order effects in plots)
rand_data <- res_data_all[sample(nrow(res_data_all)),]

library(ggplot2)
p <- ggplot(rand_data) + geom_point(aes(age, log(k), col=hybrid), size=1.5) + xlab('Age') + ylab('log(k)') + scale_colour_continuous(name="Hybrid score", low=rgb(.95, .95, 0), high=rgb(0,.3,0))   + theme(panel.background = element_rect(fill = grey(.6)))
p + facet_grid(winner ~ a)




