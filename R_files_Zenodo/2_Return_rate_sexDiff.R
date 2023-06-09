#### -------------------------------------------------------------------------------------- ####
#### Packages & functions
#### -------------------------------------------------------------------------------------- ####
{
  options(stringsAsFactors = FALSE)
  sapply(c('magrittr','data.table','brms','arm','lme4','car'),
         require, character.only = TRUE)
  
  ### VIF check on lm model
  vif.mer <- function (fit) {
    ## adapted from rms::vif
    
    v <- vcov(fit)
    nam <- names(fixef(fit))
    
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
      v <- v[-(1:ns), -(1:ns), drop = FALSE]
      nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
  }
  
}

#### -------------------------------------------------------------------------------------- ####
#### Analysis 2: Sex difference in return rate                                                                
#### -------------------------------------------------------------------------------------- ####

#### Load in data
ret.dm <- readRDS("/ds/grpkempenaers/EKwon/FIDELITY/R output/return_sexDiff_input_osf.rds") %>% data.table


#### Variables included in the dataset:
{ # [1] "popID"         "sourceID"      "family"        "mating.system" "scinam"        "subsp"         "species"      
  # [8] "male"          "female"        "dm.return"     "wingM"         "wingF"         "dm.wing"       "N"
  
#### Units
  # "male"           : the proportion of banded male birds returned, ranges from 0 to 1
  # "female"         : the proportion of banded female birds returned, ranges from 0 to 1
  # "dm.return"      : male return rate divided by the sum of male and female return rates
  # "wingM"          : species mean wing length for males in mm     
  # "wingF"          : species mean wing length for females in mm     
  # "dm.wing"        : log(wingM) - log(wingF)
  # "N"              : total number of birds banded for the population      
  
}

#### Data transformation
{
  ret.dm$dm.wing       <- scale(ret.dm$dm.wing)
  ret.dm$mating.system <- as.factor(ret.dm$mating.system)
  ret.dm$subsp         <- as.factor(ret.dm$subsp)
}   

#### Remove rows with missing values 
ret.dm[rowSums(is.na(ret.dm)) > 0,]
ret.dm <- na.omit(ret.dm)

#### Remove data from 'mixed' mating system   
ret.dm <- ret.dm[!mating.system=="mixed", ]
ret.dm$mating.system <- droplevels(ret.dm$mating.system)

#### This leaves us with 64 estimates of return rates from 33 species or 64 populations. 

#### Create a dummy variable with a binary value
ret.dm[male-female==0, bi.return := 1] # 1 if males return more or as same as female
ret.dm[male-female>0,  bi.return := 1]
ret.dm[male-female<0,  bi.return := 0] # 0 if females return more

#### In order to use the sample size (total number of banded birds) as the weight
#### on model parameterization, we first adjusted the outlier value (changed from 10000 to 1500),
#### and divided by the standard deviation
ret.dm[N==10000, N := 1500]
ret.dm$N <- ret.dm$N/sd(ret.dm$N)
#### scale the measure of sexual size dimorphism
ret.dm$dm.wing <- scale(ret.dm$dm.wing)


#### --------------------------------------------------------------------------- ####
#### MODEL A. Beta regression                                             
#### --------------------------------------------------------------------------- ####

#### MODEL A-1. Additive ------------------------------------------------------- ####
{
 ### PRIOR
  priors <- get_prior(dm.return | weights(N) ~ mating.system + dm.wing + (1|subsp), 
                      data=ret.dm, family="Beta")
  
  ### MODEL 1
  fm_GA_wt = brm(form    = dm.return | weights(N) ~ mating.system + dm.wing + (1|subsp), 
                 data    = ret.dm,   
                 family  = Beta("logit"),
                 cores   = 5,
                 chains  = 5,
                 control = list(adapt_delta = 0.99),
                 iter    = 50000,
                 thin    = 5,
                 save_pars = save_pars(all = TRUE),
                 sample_prior="yes",
                 prior   = priors 
                 )
  
  print(summary(fm_GA_wt))
  
}

#### MODEL A-2. Mating system only  -------------------------------------------- ####
{
  ### PRIOR
  priors <- get_prior(dm.return | weights(N) ~ mating.system + (1|subsp), 
                      data=ret.dm, family="Beta")
  
  ### MODEL 2
  fm_GM_wt = brm(form    = dm.return | weights(N) ~ mating.system + (1|subsp), 
                 data    = ret.dm,   
                 family  = Beta("logit"),
                 cores   = 5,
                 chains  = 5,
                 control = list(adapt_delta = 0.99),
                 iter    = 50000,
                 thin    = 5,
                 save_pars = save_pars(all = TRUE),
                 sample_prior="yes",
                 prior   = priors 
                 )
  
  print(summary(fm_GM_wt))
}
  
#### MODEL A-3. SSD only  ------------------------------------------------------ ####
{  
  ### PRIOR
  priors <- get_prior(dm.return | weights(N) ~ dm.wing + (1|subsp), 
                      data=ret.dm, family="Beta")
  
  ### MODEL 3
  fm_GS_wt = brm(form    = dm.return.prop | weights(N) ~ dm.wing + (1|subsp), 
                 data    = ret.dm,   
                 family  = Beta("logit"),
                 cores   = 5,
                 chains  = 5,
                 control = list(adapt_delta = 0.99),
                 iter    = 50000,
                 thin    = 5,
                 save_pars = save_pars(all = TRUE),
                 sample_prior="yes",
                 prior   = priors 
                 )
  
  print(summary(fm_GS_wt))
}
  
  
#### --------------------------------------------------------------------------- ####
#### MODEL B. Logistic regression                                             
#### --------------------------------------------------------------------------- ####

#### MODEL B-1. Additive ------------------------------------------------------- ####
{
  ### PRIOR
  priors <- get_prior(bi.return | weights(N) ~ mating.system + dm.wing + (1|subsp), 
                      data=ret.dm, 
                      family=bernoulli(link = "logit"))
  
  ### MODEL 4
  fm_BA_wt = brm(form    = bi.return | weights(N) ~ mating.system + dm.wing + (1|subsp), 
                 data    = ret.dm,   
                 family  = bernoulli(link = "logit"),
                 cores   = 5,
                 chains  = 5,
                 control = list(adapt_delta = 0.99),
                 iter    = 50000,
                 thin    = 5,
                 save_pars = save_pars(all = TRUE),
                 sample_prior="yes",
                 prior   = priors 
                 )
  
  print(summary(fm_BA_wt))
}

#### MODEL B-2. Mating system only  -------------------------------------------- ####
{  
  ### PRIOR
  priors <- get_prior(bi.return | weights(N) ~ mating.system + (1|subsp), 
                      data=ret.dm, 
                      family=bernoulli(link = "logit"))
  
  ### MODEL 5
  fm_BM_wt = brm(form    = bi.return | weights(N) ~ mating.system + (1|subsp), 
                 data    = ret.dm,   
                 family  = bernoulli(link = "logit"),
                 cores   = 5,
                 chains  = 5,
                 control = list(adapt_delta = 0.99),
                 iter    = 50000,
                 thin    = 5,
                 save_pars = save_pars(all = TRUE),
                 sample_prior="yes",
                 prior   = priors 
                 )
  
  print(summary(fm_BM_wt))
  
}

#### MODEL B-3. SSD only  ------------------------------------------------------ ####
{  
  ### PRIOR
  priors <- get_prior(bi.return | weights(N) ~ dm.wing + (1|subsp), 
                      data=ret.dm, 
                      family=bernoulli(link = "logit"))
  
  ### MODEL 6
  fm_BS_wt = brm(form    = bi.return | weights(N) ~ dm.wing + (1|subsp), 
                 data    = ret.dm,   
                 family  = bernoulli(link = "logit"),
                 cores   = 5,
                 chains  = 5,
                 control = list(adapt_delta = 0.99),
                 iter    = 50000,
                 thin    = 5,
                 save_pars = save_pars(all = TRUE),
                 sample_prior="yes",
                 prior   = priors 
                 )
  
  print(summary(fm_BS_wt))
}




#--------------------------------------------------------------------------------------------
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.2 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
# LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
# 
# locale:
# [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
# [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggspatial_1.1.4     RColorBrewer_1.1-2  rnaturalearth_0.1.0 glue_1.4.2          sf_0.9-7            GGally_2.1.2        geiger_2.0.7       
# [8] phytools_0.7-70     maps_3.3.0          stringr_1.4.0       rstudioapi_0.13     pals_1.6            leaflet_2.0.3       sp_1.4-5           
# [15] scales_1.1.1        tidyr_1.1.2         readxl_1.3.1        reshape2_1.4.4      doFuture_0.10.0     future_1.20.1       foreach_1.5.0      
# [22] showtext_0.9-2      showtextdb_3.0      sysfonts_0.8.3      paletteer_1.3.0     wesanderson_0.3.6   sjPlot_2.8.6        modelr_0.1.8       
# [29] tibble_3.1.2        tidybayes_2.3.1     car_3.0-10          carData_3.0-4       plyr_1.8.6          arm_1.11-2          lme4_1.1-25        
# [36] Matrix_1.3-3        MASS_7.3-54         ggplot2_3.3.4       ape_5.4-1           brms_2.14.4         Rcpp_1.0.6          data.table_1.13.2  
# [43] magrittr_2.0.1