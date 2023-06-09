
behavior <- c("gg", "gr","sg","sr","f")

d <- read.csv(paste0("survival_dataset_all.csv"), stringsAsFactors = FALSE )


# standardize predictors
d$age_s <- (d$age - mean(d$age)) / sd(d$age)
d$rank_s <- (d$rank - mean(d$rank)) / sd(d$rank)
d$group_size_s <- (d$group_size - mean(d$group_size)) / sd(d$group_size)


## standardize sociality measures
d$gg_mean_s <- (d$gg_mean - mean(d$gg_mean)) / sd(d$gg_mean)
d$gg_sd_s <- d$gg_sd  / sd(d$gg_mean)

d$gr_mean_s <- (d$gr_mean - mean(d$gr_mean)) / sd(d$gr_mean)
d$gr_sd_s <- d$gr_sd  / sd(d$gr_mean)

d$sg_mean_s <- (d$sg_mean - mean(d$sg_mean)) / sd(d$sg_mean)
d$sg_sd_s <- d$sg_sd  / sd(d$sg_mean)

d$sr_mean_s <- (d$sr_mean - mean(d$sr_mean)) / sd(d$sr_mean)
d$sr_sd_s <- d$sr_sd  / sd(d$sr_mean)

d$f_mean_s <- (d$f_mean - mean(d$f_mean)) / sd(d$f_mean)
d$f_sd_s <- d$f_sd  / sd(d$f_mean)


chains = 4
refresh = 100
iter <- 2400 #1200
warmup <- 1200 #600

sm <- stan_model(file="aft_survival_capuchins.stan")

dat <- list(
  
  I = length(unique(d$female)),
  N = length(d$female),
  G = length(unique(d$group)),
  
  index = 1:nrow( d ),
  female = d$female,
  group = d$group,
  days_survived = as.integer( d$days_survived_this_year_this_group ),
  died = as.integer( d$died ),
  
  # Prepare for individual mean and sd for sociality measure
  ind_rate_mean = NA,
  ind_rate_sd = NA,
  
  # other
  age_s = as.numeric( d$age_s ),
  rank_s = as.numeric( d$rank_s ),
  group_size_s = as.numeric( d$group_size_s )
  
)

for(i in 1:length(behavior) ){

# update sociality measure
dat$ind_rate_mean <- as.numeric( eval( parse( text=paste0("d$",behavior[i],"_mean_s") ) ) )
dat$ind_rate_sd <- as.numeric( eval( parse( text=paste0("d$",behavior[i],"_sd_s") ) ) )

#run model
res <- sampling(object = sm, data = dat, chains = chains, refresh=refresh, iter = iter, warmup = warmup)

assign( paste0("m",behavior[i]), res )
}

save.image("aft_models_all.RDATA")
