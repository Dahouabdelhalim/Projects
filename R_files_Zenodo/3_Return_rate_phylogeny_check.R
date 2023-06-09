#### -------------------------------------------------------------------------------------- ####
#### Packages & functions
#### -------------------------------------------------------------------------------------- ####
{
  options(stringsAsFactors = FALSE)
  sapply(c('magrittr','data.table','brms','arm','lme4','car',
           'ape','phytools','geiger','stringr','plyr'),
         require, character.only = TRUE)
  
  physet <- function(x, p = phy_subset) {
      ape::drop.tip(p, setdiff(p$tip.label, x$scinam) )
    }
  # function to trim the phylo tree fitting to the dataset
    
  phy2A <- function(phys) {
      INphylo <- MCMCglmm::inverseA(phys, scale = TRUE)
      A <- solve(INphylo$Ainv)
      rownames(A) <- rownames(INphylo$Ainv)
      A
    }
  # function to generate the matrix of phylogenetic relatedness
    
  firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
  # function to capitalize the first letter
    
  fitLambda<-function(tree,x,model="ER"){
      lik<-function(lambda,tree,x,model)
        logLik(ace(x,rescale(tree,model="lambda",lambda),
                   type="discrete",model=model))
      obj<-optimize(lik,c(0,1),tree=tree,x=x,model=model,maximum=TRUE)
      fit<-ace(x,rescale(tree,model="lambda",lambda=obj$maximum),
               type="discrete",model=model)
      I<-fit$index.matrix
      fitted.Q=matrix(fit$rates[I],dim(I)[1],dim(I)[2],
                      dimnames=list(dimnames(fit$lik.anc)[[2]],
                                    dimnames(fit$lik.anc)[[2]]))
      diag(fitted.Q)<--rowSums(fitted.Q,na.rm=TRUE)
      list(Q=fitted.Q,lambda=obj$maximum,logLik=logLik(fit))
    }
   
} 
    

#### -------------------------------------------------------------------------------------- ####
#### Check 1: Return rate                                                                   ####
#### -------------------------------------------------------------------------------------- ####
{
#### Load in data
  ## Phylogenetic tree
  phy = read.nexus("DATA/phylo_maxCladeCred.tre")
  ## Bayesian model 'fm1' from r_script "1_Return_rate" or the following
  fm1 <- readRDS("DATA/return.RDS")
  ## Raw data input for 'fm1'
  ret <- readRDS("DATA/return_input.rds") %>% data.table
  
  
#### Variables included in the dataset:
  { # [1]  "popID"             "sourceID"          "family"            "mating.system"     "scinam"           
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
  
  
#### Residuals from fm1
  fm1_res <- residuals(fm1, type='pearson') %>% data.table
  
  ret.sp <- cbind(ret.s, fm1_res) 
  
  ret.sp[species=="dunlin", scinam := "calidris alpina"]
  ret.sp[species=="stilt sandpiper", scinam := "micropalma himantopus"]
  
  ret.sp[, scinam := firstup(as.character(as.factor(scinam)))]
  ret.sp[, scinam := str_replace_all(string=scinam, pattern=" ", repl="_")]
  ret.sp$scinam <- as.character(as.factor(ret.sp$scinam))
  
  
  x <- ret.sp
  x <- x[!scinam %in% c("Charadrius_morinellus","Charadrius_novaeseelandiae","Charadrius_dubius"),]
  
  
  d <- c(phy$tip.label)
  d[!d %in% c(unique(x$scinam))]
  phy_s <- drop.tip(phy, tip=c(d[!d %in% c(unique(x$scinam))]))
  
  A = physet(x, phy_s)  %>% phy2A
  
  INphylo <- MCMCglmm::inverseA(phy_s, scale = TRUE)
  A <- solve(INphylo$Ainv)
  rownames(A) <- rownames(INphylo$Ainv)
  A
  
  x2 <- x[,c("Estimate","scinam")]
  x.s <- ddply(x2, .(scinam), summarise,
               spmean = mean(Estimate, na.rm=TRUE))
  
  x.ss <- setNames(x.s$spmean, x.s$scinam)
  
  phy_s$tip.label[!phy_s$tip.label %in% x.s$scinam]
  x.s$scinam[!x.s$scinam %in% phy_s$tip.label]
  
  
  fitLambda(phy_s, x.ss) 
  
  
#### Fitting intercept-only model to the residuals
#### WITH phylogeny
  priors <- get_prior(Estimate ~ 0 + Intercept + (1|species) + (1|gr(scinam, cov = A)), data=x)
  
  fm1.res.wP = brm(form    = Estimate ~ 0 + Intercept + (1|species) + (1|gr(scinam, cov = A)), 
                   data    = x,   
                   data2   = list(A = A),  #phylo_relatedness matrix
                   cores   = 5,
                   chains  = 5,
                   control = list(adapt_delta = 0.99),
                   iter    = 5000,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE), 
                   prior   = priors 
                   )
  
  summary(fm1.res.wP)
  
 
#### Fitting intercept-only model to the residuals
#### WITHOUT phylogeny
  priors <- get_prior(Estimate ~ 0 + Intercept + (1|species), data=x)
  
  fm1.res.woP = brm(form    = Estimate ~ 0 + Intercept + (1|species), 
                    data    = x,   
                    cores   = 5,
                    chains  = 5,
                    control = list(adapt_delta = 0.99),
                    iter    = 5000,
                    sample_prior="yes",
                    save_pars = save_pars(all = TRUE),
                    prior   = priors 
                    )
  
  summary(fm1.res.woP)
  
 
#### Model comparison
#### WITH vs. WITHOUT phylogeny
  
# 1. LOOic
  loo(fm1.res.wP)    
  loo(fm1.res.woP)   
  
  
# 2. Bayes factor
  bayes_factor(fm1.res.woP, fm1.res.wP)
  
  
# 3. Poterior probability  
  post_prob(fm1.res.wP, fm1.res.woP)
 
  
# 4. Bayesian hypothesis testing  
  hypothesis(fm1.res.wP, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sigma^2) = 0",  class = NULL)
  
}  
  
  
#### -------------------------------------------------------------------------------------- ####
#### Check 2: Sex difference in return rate                                                                   ####
#### -------------------------------------------------------------------------------------- ####
{  
#### Load in data
  ## Phylogenetic tree
  phy = read.nexus("DATA/phylo_maxCladeCred.tre")
  ## Bayesian model 'fm_GA_wt' from r_script "2_Return_rate_sexDiff" or the following
  fm_GA_wt <- readRDS("DATA/return_sexDiff.RDS")
  ## Raw data input for 'fm1'
  ret.dm <- readRDS("DATA/return_sexDiff_input_osf.rds") %>% data.table
  
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
  
  
  
#### Residuals from fm_GA_wt
  fm2_res <- residuals(fm_GA_wt, type='pearson') %>% data.table
  
  ret.dmp <- cbind(ret.dm, fm2_res) 
  
  ret.dmp[species=="dunlin", scinam := "calidris alpina"]
  ret.dmp[species=="stilt sandpiper", scinam := "micropalma himantopus"]
  
  ret.dmp[, scinam := firstup(as.character(as.factor(scinam)))]
  ret.dmp[, scinam := str_replace_all(string=scinam, pattern=" ", repl="_")]
  ret.dmp$scinam <- as.character(as.factor(ret.dmp$scinam))
  
  x <- ret.dmp
  x <- x[!scinam %in% c("Charadrius_morinellus","Charadrius_novaeseelandiae","Charadrius_dubius"),]
  
  d <- c(phy$tip.label)
  d[!d %in% c(unique(x$scinam))]
  phy_s <- drop.tip(phy, tip=c(d[!d %in% c(unique(x$scinam))]))
  
  A = physet(x, phy_s)  %>% phy2A
  
  INphylo <- MCMCglmm::inverseA(phy_s, scale = TRUE)
  A <- solve(INphylo$Ainv)
  rownames(A) <- rownames(INphylo$Ainv)
  A
  
  x2 <- x[,c("Estimate","scinam")]
  x.s <- ddply(x2, .(scinam), summarise,
               spmean = mean(Estimate, na.rm=TRUE))
  
  x.ss <- setNames(x.s$spmean, x.s$scinam)
  
  phy_s$tip.label[!phy_s$tip.label %in% x.s$scinam]
  x.s$scinam[!x.s$scinam %in% phy_s$tip.label]
  
  
  fitLambda(phy_s, x.ss)
  
  
  #### Fitting intercept-only model to the residuals
  #### WITH phylogeny
  priors <- get_prior(Estimate ~ 0 + Intercept + (1|species) + (1|gr(scinam, cov = A)), data=x)
  
  fm2.res.wP = brm(form    = Estimate ~ 0 + Intercept + (1|species) + (1|gr(scinam, cov = A)), 
                   data    = x,   
                   data2   = list(A = A),  #phylo_relatedness matrix
                   cores   = 5,
                   chains  = 5,
                   control = list(adapt_delta = 0.99),
                   iter    = 5000,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE), 
                   prior   = priors 
  )
  
  summary(fm2.res.wP)
  
  
  #### Fitting intercept-only model to the residuals
  #### WITHOUT phylogeny
  priors <- get_prior(Estimate ~ 0 + Intercept + (1|species), data=x)
  
  fm2.res.woP = brm(form    = Estimate ~ 0 + Intercept + (1|species), 
                    data    = x,   
                    cores   = 5,
                    chains  = 5,
                    control = list(adapt_delta = 0.99),
                    iter    = 5000,
                    sample_prior="yes",
                    save_pars = save_pars(all = TRUE),
                    prior   = priors 
  )
  
  summary(fm2.res.woP)
  
  
  #### Model comparison
  #### WITH vs. WITHOUT phylogeny
  
  # 1. LOOic
  loo(fm2.res.wP)    
  loo(fm2.res.woP)   
  
  
  # 2. Bayes factor
  bayes_factor(fm2.res.woP, fm2.res.wP)
  
  
  # 3. Poterior probability  
  post_prob(fm2.res.wP, fm2.res.woP)
  
  
  # 4. Bayesian hypothesis testing  
  hypothesis(fm2.res.wP, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sigma^2) = 0",  class = NULL)
  
}