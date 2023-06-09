library(rethinking)

###################################################
# import data
de <- read.csv("survival_dataset_all.csv")
df <- read.csv("survival_dataset_ff.csv")
dm <- read.csv("survival_dataset_fm.csv")

# make sure data has only complete cases
df <- df[complete.cases(df),]
dm <- dm[complete.cases(dm),]

# load fitted models and assign to new objects, because model names are the same across a three types of datasets
load("aft_models_all.RDATA")
# objects with "a" refer to all group members dataset
agg <- mgg
agr <- mgr
asg <- msg
asr <- msr
af <- mf

load("aft_models_ff.RDATA")
# objects with "f" refer to adult female dataset
fgg <- mgg
fgr <- mgr
fsg <- msg
fsr <- msr
ff <- mf

load("aft_models_fm.RDATA")
# objects with "m" refer to adult male dataset (no need to assign new objects)

###################################################

# model names that you need to loop through
models <- c("agg","agr","asg","asr","af",
            "fgg","fgr","fsg","fsr","ff",
            "mgg", "mgr", "msg", "msr", "mf")

# Fix age to mean age (standardized)
age <- 0
#  Fix rank: the mean rank (standardized)
rank <- 0
# fix mean group size (standardized)
group_size <- 0

# create an empty list to store sociality measures for each behavior
soc_stand <- list()
soc_stand[[1]] <- (de$gg_mean - mean(de$gg_mean) ) / sd(de$gg_mean)
soc_stand[[2]] <- (de$gr_mean - mean(de$gr_mean) ) / sd(de$gr_mean)
soc_stand[[3]] <- (de$sg_mean - mean(de$sg_mean) ) / sd(de$sg_mean)
soc_stand[[4]] <- (de$sr_mean - mean(de$sr_mean) ) / sd(de$sr_mean)
soc_stand[[5]] <- (de$f_mean - mean(de$f_mean) ) / sd(de$f_mean)

soc_stand[[6]] <- (df$gg_mean - mean(df$gg_mean) ) / sd(df$gg_mean)
soc_stand[[7]] <- (df$gr_mean - mean(df$gr_mean) ) / sd(df$gr_mean)
soc_stand[[8]] <- (df$sg_mean - mean(df$sg_mean) ) / sd(df$sg_mean)
soc_stand[[9]] <- (df$sr_mean - mean(df$sr_mean) ) / sd(df$sr_mean)
soc_stand[[10]] <- (df$f_mean - mean(df$f_mean) ) / sd(df$f_mean)

soc_stand[[11]] <- (dm$gg_mean - mean(dm$gg_mean) ) / sd(dm$gg_mean)
soc_stand[[12]] <- (dm$gr_mean - mean(dm$gr_mean) ) / sd(dm$gr_mean)
soc_stand[[13]] <- (dm$sg_mean - mean(dm$sg_mean) ) / sd(dm$sg_mean)
soc_stand[[14]] <- (dm$sr_mean - mean(dm$sr_mean) ) / sd(dm$sr_mean)
soc_stand[[15]] <- (dm$f_mean - mean(dm$f_mean) ) / sd(dm$f_mean)

# get min and max and make a sequence between them. This will be used for calculations.
soc_seq <- list()
for( i in 1:15){
  soc_seq[[i]] <- seq( min(soc_stand[[i]]), max(soc_stand[[i]]), length.out = 20 )
}

# convert sequence to the natural scale. This willl be used for the plot.
soc_nat <- list()
soc_nat[[1]] <- soc_seq[[1]] * sd(de$gg_mean) + mean(de$gg_mean)
soc_nat[[2]] <- soc_seq[[2]] * sd(de$gr_mean) + mean(de$gr_mean)
soc_nat[[3]] <- soc_seq[[3]] * sd(de$sg_mean) + mean(de$sg_mean)
soc_nat[[4]] <- soc_seq[[4]] * sd(de$sr_mean) + mean(de$sr_mean)
soc_nat[[5]] <- soc_seq[[5]] * sd(de$f_mean) + mean(de$f_mean)


soc_nat[[6]] <- soc_seq[[6]] * sd(df$gg_mean) + mean(df$gg_mean)
soc_nat[[7]] <- soc_seq[[7]] * sd(df$gr_mean) + mean(df$gr_mean)
soc_nat[[8]] <- soc_seq[[8]] * sd(df$sg_mean) + mean(df$sg_mean)
soc_nat[[9]] <- soc_seq[[9]] * sd(df$sr_mean) + mean(df$sr_mean)
soc_nat[[10]] <- soc_seq[[10]] * sd(df$f_mean) + mean(d$f_mean)

soc_nat[[11]] <- soc_seq[[11]] * sd(dm$gg_mean) + mean(dm$gg_mean)
soc_nat[[12]] <- soc_seq[[12]] * sd(dm$gg_mean) + mean(dm$gg_mean)
soc_nat[[13]] <- soc_seq[[13]] * sd(dm$gg_mean) + mean(dm$gg_mean)
soc_nat[[14]] <- soc_seq[[14]] * sd(dm$gg_mean) + mean(dm$gg_mean)
soc_nat[[15]] <- soc_seq[[15]] * sd(dm$gg_mean) + mean(dm$gg_mean)



# compute predictions:
for( i in models){
  
  # extract posterior samples
  post <- extract.samples( eval( parse( text=paste0( i ) ) ) )
  
  # social integration measure sequence 
  soc <- soc_seq[[ which(i==models) ]]
  
  # loop through social integration measures and compute predictions
    mu_link <- function(x) exp( post$base_rate + 
                                  post$b_sociality * x + 
                                  post$b_rank * rank + 
                                  post$b_age * age + 
                                  post$b_group_size * group_size )
    
    
    mu <- sapply( soc, mu_link ) 
    
    # number of days
    mu_mean <- apply( mu , 2 , mean )
    mu_HPDI <- apply( mu , 2 , HPDI , prob=0.89 )
    
    # daily death rate
    lambda_day <- 1/mu
    # approx. annual probability of dying
    lambda <- 1 - (1 - lambda_day)^365 
    
    # approx annual probability of dying
    lambda_mean <- apply( lambda , 2 , mean )
    lambda_HPDI <- apply( lambda , 2 , HPDI , prob=0.89 )
    
    # save results
    results_lambda_mean <- lambda_mean
    results_lambda_HPDI<- lambda_HPDI
  
  
  # store the results outside of the loop
  assign( paste0(i, "_mean"), results_lambda_mean )
  assign( paste0(i, "_HPDI"), results_lambda_HPDI )
  
}

# the data is ready for plotting



