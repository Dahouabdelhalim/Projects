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
#### Analysis 1: Return rate                                                                ####
#### -------------------------------------------------------------------------------------- ####

#### Load in data
ret <- readRDS(".../return_input_osf.rds") %>% data.table


#### Variables included in the dataset:
{# [1]  "popID"             "sourceID"          "family"            "mating.system"     "scinam"           
 # [6]  "species"           "sex"               "return.rate"       "N"                 "wing.length"      
 # [11] "longitude"         "latitude"          "study.area.size"   "distance.boundary" "distance.qt"      
 # [16] "range.span"        "relative.latitude" "start.year"        "end.year"          "study.duration"   
 # [21] "location"

#### Units
# "return.rate" : the proportion of banded birds returned, ranges from 0 to 1
# "N"           : number of birds      
# "wing.length" : mm     
# "longitude"   : in decimal degrees      
# "latitude"    : in decimal degrees      
# "study.area.size"   : sqkm
# "distance.boundary" : m
# "distance.qt"       : quantile
# "range.span"        : 1000km
# "relative.latitude" : 1000km
# "study.duration"    : years
}

#### Data transformation
{
  ### Log(e)-transformation: study area size, N2
  ### Log(10)-transformation: wing length
  ret$log.study.area <- log(ret$study.area.size)
  ret$log.N          <- log(ret$N)
  ret$log.wing       <- log10(ret$wing.length)
  
  ### SCALE (subtract mean and divided by the 1SD)
  ret$log.N       <- scale(ret$log.N)
  ret$log.area    <- scale(ret$log.study.area)
  ret$study.dur   <- scale(ret$study.duration)
  ret$start.yr    <- scale(ret$start.year)
  ret$rel.lat     <- scale(ret$relative.latitude)
  ret$dist.qt     <- scale(ret$distance.qt)
  ret$range.span  <- scale(ret$range.span)
  ret$log.wing    <- scale(ret$log.wing)
  
  ret$mating.system <- as.factor(ret$mating.system)
  ret$sex <- as.factor(ret$sex)
  ret$subsp <- as.factor(ret$subsp)
}   

#### Remove rows with missing values 
  ret[rowSums(is.na(ret)) > 0,]
  ret <- na.omit(ret)

#### Remove data from 'mixed' mating system   
  ret.s <- ret[!mating.system=="mixed", ]
  ret.s$mating.system <- droplevels(ret.s$mating.system)
  
#### This leaves us with 175 estimates of return rates from 49 species or 111 populations.  

#### Check for multicollinearity  
 
  lm1 <- lmer(return.rate ~ sex + mating.system + log.N + study.dur + log.area +  
                            rel.lat + dist.qt + start.yr + log.wing + 
                            (1|subsp), data=ret.s)
  vif.mer(lm1) # all < 1.5
  

#### Setting priors
#### Generate priors using the ‘get_prior’ function in ‘brms’, which sets non-informative priors 
#### for all slope coefficients and uses a Student’s t distribution for the intercept, 
#### standard deviation, and uses a gamma distribution for phi for the beta logistic model
  priors <- get_prior(return.rate ~ sex + mating.system + log.N + study.dur + log.area +  
                                    rel.lat + dist.qt + start.yr + log.wing + 
                                    (1|subsp), 
                      data=ret.s, family="Beta")
  
  
#### Fitting the model: beta logistic regression model with a logit-link function 
 
  fm1 = brm(
    form    = return.rate ~ sex + mating.system + log.N + study.dur + log.area +  
                            rel.lat + dist.qt + start.yr + log.wing + (1|subsp), 
    data    = ret.s,   
    family  = Beta("logit"),
    cores   = 5,
    chains  = 5,
    control = list(adapt_delta = 0.99),
    iter    = 50000,
    thin    = 5,
    sample_prior="yes",
    prior   = priors 
  )
  
  print(summary(fm1))
  
  mcmc_plot(fm1, type = "hist")
  mcmc_plot(fm1, type = "trace")
  mcmc_plot(fm1, type = "neff")
  mcmc_plot(fm1, type = "rhat")
  
  
  #--------------------------------------------------------------------------------------------
  # print(sessionInfo())
  
  # R version 4.1.0 (2021-05-18)
  # Platform: x86_64-pc-linux-gnu (64-bit)
  # Running under: Ubuntu 20.04.2 LTS
  # 
  # Matrix products: default
  # BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
  # LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
  # 
  # locale:
  # [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8   
  # [6] LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C        
  # [11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
  # 
  # attached base packages:
  # [1] stats     graphics  grDevices utils     datasets  methods   base     
  # 
  # other attached packages:
  # [1] ggspatial_1.1.4     RColorBrewer_1.1-2  rnaturalearth_0.1.0 glue_1.4.2          sf_0.9-7            GGally_2.1.2       
  # [7] geiger_2.0.7        phytools_0.7-70     maps_3.3.0          stringr_1.4.0       rstudioapi_0.13     pals_1.6           
  # [13] leaflet_2.0.3       sp_1.4-5            scales_1.1.1        tidyr_1.1.2         readxl_1.3.1        reshape2_1.4.4     
  # [19] doFuture_0.10.0     future_1.20.1       foreach_1.5.0       showtext_0.9-2      showtextdb_3.0      sysfonts_0.8.3     
  # [25] paletteer_1.3.0     wesanderson_0.3.6   sjPlot_2.8.6        modelr_0.1.8        tibble_3.1.2        tidybayes_2.3.1    
  # [31] car_3.0-10          carData_3.0-4       plyr_1.8.6          arm_1.11-2          lme4_1.1-25         Matrix_1.3-3       
  # [37] MASS_7.3-54         ggplot2_3.3.4       ape_5.4-1           brms_2.14.4         Rcpp_1.0.6          data.table_1.13.2  
  # [43] magrittr_2.0.1  