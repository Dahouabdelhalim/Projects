
#=============================================================================================================#
# Script created by Antica Culina, a.culina@yahoo.com
# Script created in version R 1.4.1106
# This script is for running meta-analytic & meta-regression models on within season divorce in birds
# from Culina & Brouwer 2022:No evidence of immediate fitness benefits of within-season divorce in monogamous birds 
#=============================================================================================================#

# Set working directory to wherever you have downloaded the data and code
#session info  can be found at the end of this code

# inspiration: from Shinichi Nakagawa and Malgorzata Lagisz's code: https://github.com/mlagisz and Sanchez-Tojar https://osf.io/cwkxb/ 

##############################################################
# Description of script and Instructions
##############################################################

# This code is used to conduct a multilevel meta-analyses and meta-regressions
# on the breeding success and divorce in monogamous birds
# three separate meta-analysis were conducted:
#### A) on the breeding success before divorce and probability of divorce 
#### B) on the breeding success after divorce and occurrence of divorce
#### c) on the adaptivenes of divorce for females vs males


##############################################################
# Packages needed
##############################################################
# Clear memory
rm(list=ls())


# load libraries and packages

library(MCMCglmm)
library(ape)
library(phangorn) # needed to convert phil tree into an ulatrametic one

sessionInfo()


##############################################################
# Load  Data Sets
##############################################################

before <- read.table(file.choose(), header = T, sep = ";")  # breeding success before divorce and occurrence of divorce MA

after<- read.table(file.choose(), header = T, sep = ";") # occurrence of divorce and breeding success after divorce and MA

FvsM <- read.table(file.choose(), header = T, sep = ";") # between sexes benefits and costs of partner change MA

#loading  and adjusting phylogenetic tree for before and after MA
# change the path name to find your file

tree <- read.nexus("F:\\\\NIOO3\\\\My Projects\\\\within season divorce\\\\wtihin_season divorce\\\\paper drafting\\\\biology letters\\\\Revision\\\\data_code\\\\all.nex")

# tree is non ultrametric -> needs to be changed to ultrametirc
# force an utrametric tree

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\\n\\n")
  tree
}

ult.nnls<-force.ultrametric(tree[[1]]) ## default method
is.ultrametric(ult.nnls)

plot(tree[[1]], cex = 1, main="phylogenetic tree of species included in the overall analysis of divorce and breeding success")
Ainv<-inverseA(ult.nnls)$Ainv

prior1 <- list(G=list(G1=list(V=1, nu=0.002)), R=list(V=1, nu=0.002))  # defining prior

##################################################################
# Before MA: breeding success before divorce and occurence of divorce
################################################################## 

 ### scale countinuous predictors
before$year <- as.numeric(scale(before$year, scale = FALSE))

basic <- MCMCglmm(Z~1, mev=before$SE^2, burnin=200000, nitt=2000000,thin=1000, prior = prior1, data=before)

spphyl.1 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=before)
sp.1 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=before$SE^2, burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=before)
study.1 <- MCMCglmm(Z~1, random=~title, pr=TRUE, mev=before$SE^2, burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=before)

summary(basic) 
summary(spphyl.1) #for random effect it reports variance and credible intervals
summary(sp.1)
summary(study.1)

plot(spphyl.1$Sol)

# Run random model two more times
spphyl.2 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=before$SE^2, ginverse=list(latin_name=Ainv_b), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=before)
spphyl.3 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=before$SE^2, ginverse=list(latin_name=Ainv_b), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=before)


### Gelman-Rubin diagnostic for model convergence (shirnk factor needs to be <1.05)
gelman.diag(list(spphyl.1$Sol[,1],spphyl.2$Sol[,1],spphyl.3$Sol[,1]))
gelman.diag(list(spphyl.1$VCV[,c(1)],spphyl.2$VCV[,c(1)],spphyl.3$VCV[,c(1)]))
gelman.diag(list(spphyl.1$Deviance,spphyl.2$Deviance,spphyl.3$Deviance))

# All of these are either 1 or very close, indicating very good convergence between runs

# Check autocorrelation between estimates
# This should be less than 0.1 at the end of the run
autocorr(spphyl.1$Sol[,1])
autocorr(spphyl.1$VCV[,c(1,3)])
autocorr(spphyl.2$Sol[,1])
autocorr(spphyl.2$VCV[,c(1,3)])
autocorr(spphyl.3$Sol[,1])
autocorr(spphyl.3$VCV[,c(1,3)])



##########################        
###   Model heterogeneity ###
##########################
        
# Typical measurement error variance
Weight<-1/before$SE^2
MV<-sum(Weight*(length(Weight)-1))/(sum(Weight)^2-sum(Weight^2))

#phylogeny model

total_var <- (spphyl.1$VCV[,"latin_name"] + spphyl.1$VCV[,"units"]+MV)

        
# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 
      
# I2 for phylogeny
I2p<-spphyl.1$VCV[,"latin_name"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 
        
        
# I2 for units (residuals)
I2u<-spphyl.1$VCV[,"units"]/total_var
summary(I2u)
        
### fixed effects (overall intercept)
mean(spphyl.1$Sol[,1])
HPDinterval(spphyl.1$Sol[,1])

#species model

total_var <- (sp.1$VCV[,"latin_name"] + sp.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for species
I2p<-sp.1$VCV[,"latin_name"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 


# I2 for units (residuals)
I2u<-sp.1$VCV[,"units"]/total_var
summary(I2u)

#study model

total_var <- (study.1$VCV[,"title"] + study.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for phylogeny
I2p<-study.1$VCV[,"title"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 


# I2 for units (residuals)
I2u<-sp.1$VCV[,"units"]/total_var
summary(I2u)

##########################        
###   METHODOLGOCAL MODERATORS
##########################


## is the measure success dichotmized?
dich <- MCMCglmm (Z~dichotomisation-1, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=before)

## whether study only takes mating outcomes after failure or all
after_failure <- MCMCglmm (Z~is_renesting_after_failure-1, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=before)


summary(dich)
summary(after_failure)


########################### biological moderator: BS [breeding success] 

BS <- MCMCglmm (Z~BS2-1, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=before)
BS2 <- MCMCglmm (Z~dichotomisation + BS2 -1, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=before)

summary(BS)
summary(BS2)


########## frequency of extra pair offspring

model_EPP <- MCMCglmm (Z~EPP, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=before)

model_EPP2 <- MCMCglmm (Z~dichotomisation + EPP, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=before)

summary(model_EPP)
summary(model_EPP2)

######################################################  
######## after MA
########################### ########################### 

### scale continuous predictors
after$year <- as.numeric(scale(after$year, scale = FALSE))

basic <- MCMCglmm(Z~1, mev=after$SE^2, burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=after)
             
spphyl.1 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=after)
sp.1 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=after$SE^2, burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=after)
study.1 <- MCMCglmm(Z~1, random=~title, pr=TRUE, mev=after$SE^2, burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=after)

summary(basic) 
summary(spphyl.1) 
summary(sp.1)
summary(study.1) 


# Run random model two more times
spphyl.2 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=after)
spphyl.3 <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=after)


### Gelman-Rubin diagnostic for model convergence (shirnk factor needs to be <1.05)
gelman.diag(list(spphyl.1$Sol[,1],spphyl.2$Sol[,1],spphyl.3$Sol[,1]))
gelman.diag(list(spphyl.1$VCV[,c(1)],spphyl.2$VCV[,c(1)],spphyl.3$VCV[,c(1)]))
gelman.diag(list(spphyl.1$Deviance,spphyl.2$Deviance,spphyl.3$Deviance))

# All of these are either 1 or very close, indicating very good convergence between runs

# Check autocorrelation between estimates
# This should be less than 0.1 at the end of the run
autocorr(spphyl.1$Sol[,1])
autocorr(spphyl.1$VCV[,c(1,3)])
autocorr(spphyl.2$Sol[,1])
autocorr(spphyl.2$VCV[,c(1,3)])
autocorr(spphyl.3$Sol[,1])
autocorr(spphyl.3$VCV[,c(1,3)])



##########################        
###   Model heterogeneity ###
##########################

# Typical measurement error variance
Weight<-1/after$SE^2
MV<-sum(Weight*(length(Weight)-1))/(sum(Weight)^2-sum(Weight^2))


### phylogeny

total_var <- (spphyl.1$VCV[,"latin_name"] + spphyl.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for phylogeny
I2p<-spphyl.1$VCV[,"latin_name"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

### fixed effects (overall intercept)
mean(spphyl.1$Sol[,1])
HPDinterval(spphyl.1$Sol[,1])

# species

total_var <- (sp.1$VCV[,"latin_name"] + sp.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for species
I2p<-sp.1$VCV[,"latin_name"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 


# study

total_var <- (study.1$VCV[,"title"] + study.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for phylogeny
I2p<-study.1$VCV[,"title"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

##########################        
###   METHODOLGOCAL MODERATORS
##########################

## experimental or observational study
m.study <- MCMCglmm (Z~type_of_study - 1, random=~latin_name, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=after)

summary(m.study)



# is renesting after failure?

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


afterf <- completeFun(after, "is_renesting_after_failure")
ext <- MCMCglmm(Z~iraf2 - 1, random=~latin_name, pr=TRUE, mev=afterf$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=afterf)
summary(ext)

contr <- MCMCglmm(Z~ 1, random=~latin_name, pr=TRUE, mev=afterf$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=afterf)
summary(contr)


########################### biological moderator: BS [breeding success], sex

m.BS <- MCMCglmm (Z~BS2-1, random=~latin_name, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=after)

summary(m.BS)



########## frequency of extra pair offspring
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


afterepp <- completeFun(after, "EPP")
afterepp$EPP <- as.numeric(scale(afterepp$EPP, scale = FALSE))

model_EPP <- MCMCglmm (Z~EPP, random=~latin_name, mev=afterepp$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=afterepp)
basic <- MCMCglmm (Z~1, random=~latin_name, mev=afterepp$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=afterepp)

summary(model_EPP)
summary(basic)


######################################################################
#############  Between-sexes MA: Females vs M ########################
######################################################################

### scale continuous predictors
FvsM$year <- as.numeric(scale(FvsM$year, scale = FALSE))

basic <- MCMCglmm (Z~1, mev=FvsM$SE^2, thin=1000, burnin=200000, nitt=2000000, prior=prior1, data=FvsM)

spphil.1 <- MCMCglmm (Z~1, random=~latin_name, mev=FvsM$SE^2, ginverse=list(latin_name=Ainv), thin=1000, burnin=200000, nitt=2000000, prior=prior1, data=FvsM)
sp.1 <- MCMCglmm (Z~1, random=~latin_name, mev=FvsM$SE^2, thin=1000, burnin=200000, nitt=2000000, prior=prior1, data=FvsM)
study.1 <- MCMCglmm (Z~1, random=~title, mev=FvsM$SE^2, thin=1000, burnin=200000, nitt=2000000, prior=prior1, data=FvsM)

summary(basic)
summary(spphil.1)
summary(sp.1)
summary(study.1)


# study is the best

###   Model heterogeneity ###
##########################

# Run random model two more times
study.2 <- MCMCglmm (Z~1, random=~title, mev=FvsM$SE^2, thin=1000, burnin=200000, nitt=2000000, prior=prior1, data=FvsM)
study.3 <- MCMCglmm (Z~1, random=~title, mev=FvsM$SE^2, thin=1000, burnin=200000, nitt=2000000, prior=prior1, data=FvsM)


### Gelman-Rubin diagnostic for model convergence (shirnk factor needs to be <1.05)
gelman.diag(list(study.1$Sol[,1],study.2$Sol[,1],study.3$Sol[,1]))
gelman.diag(list(study.1$VCV[,c(1)],study.2$VCV[,c(1)],study.3$VCV[,c(1)]))
gelman.diag(list(study.1$Deviance,study.2$Deviance,study.3$Deviance))

# All of these are either 1 or very close, indicating very good convergence between runs

# Check autocorrelation between estimates
# This should be less than 0.1 at the end of the run
autocorr(study.1$Sol[,1])
autocorr(study.1$VCV[,c(1,3)])
autocorr(study.2$Sol[,1])
autocorr(study.2$VCV[,c(1,3)])
autocorr(study.3$Sol[,1])
autocorr(study.3$VCV[,c(1,3)])



##########################        
###   Model heterogeneity ###
##########################

# Typical measurement error variance
Weight<-1/FvsM$SE^2
MV<-sum(Weight*(length(Weight)-1))/(sum(Weight)^2-sum(Weight^2))


 # study

total_var <- (study.1$VCV[,"title"] + study.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for study
I2p<-study.1$VCV[,"title"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

# I2 for units
I2p<-study.1$VCV[,"units"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

### fixed effects (overall intercept)
mean(study.1$Sol[,1])
HPDinterval(study.1$Sol[,1])

# phylogeny

total_var <- (spphil.1$VCV[,"latin_name"] + spphil.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for phyl
I2p<-spphil.1$VCV[,"latin_name"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

# I2 for units
I2p<-spphil.1$VCV[,"units"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

# species

total_var <- (sp.1$VCV[,"latin_name"] + sp.1$VCV[,"units"]+MV)

# total I2
I2t<-(total_var - MV)/total_var
summary(I2t)
mean(I2t)
HPDinterval(I2t) 

# I2 for study
I2p<-sp.1$VCV[,"latin_name"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 

# I2 for units
I2p<-sp.1$VCV[,"units"]/total_var
summary(I2p)
mean(I2p)
HPDinterval(I2p) 



########################### methdological moderators, keepng study

m.study <- MCMCglmm (Z~type_of_study-1, random=~title, mev=FvsM$SE^2, burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=FvsM)

m.comparisson <- MCMCglmm (Z~comparisson-1, random=~title, mev=FvsM$SE^2, burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=FvsM)

summary(m.study)
summary(m.comparisson)


# MA with effect sizes where renesting is after failure

FvsM3 <- subset(FvsM, is_renesting_after_failure =='yes')
fail <- MCMCglmm(Z~1, random=~latin_name, pr=TRUE, mev=FvsM3$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin= 1000, prior=prior1, data=FvsM3)
summary(fail)

#### to look into BS, only consider BS mesured in the t+1 season 


BS <- MCMCglmm (Z~BS2-1, random=~title, mev=FvsM$SE^2,  burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=FvsM)
summary(BS)


### re-level to get the estimates for widowed
# levels(FvsM$comparisson)
# FvsM$comparisson <- as.factor(FvsM$comparisson)
# FvsM$comparisson <- factor(FvsM$comparisson, levels=c("within_widowed", "within_divorced"))


########################### EPP

EP <- FvsM[ which(FvsM$EPP !='NA'), ]

basic <- MCMCglmm (Z~ 1, random=~title, mev=EP$SE^2, burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=EP)
m.EPP<- MCMCglmm (Z~EPP, random=~title, mev=EP$SE^2, burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=EP)

summary(basic)
summary(m.EPP)


########################### backtransforming Zr into r


library(psych)
fisherz2r()

#===============================================================
# session info
#===============================================================

#R version 4.0.4 (2021-02-15)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 14393)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_Europe.1252  LC_CTYPE=English_Europe.1252    LC_MONETARY=English_Europe.1252
#[4] LC_NUMERIC=C                    LC_TIME=English_Europe.1252    

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] phangorn_2.7.0 MCMCglmm_2.32  ape_5.5        coda_0.19-4    Matrix_1.3-2  

#loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.6       codetools_0.2-18 quadprog_1.5-8   lattice_0.20-41  corpcor_1.6.9    grid_4.0.4      
#[7] nlme_3.1-152     magrittr_2.0.1   cubature_2.0.4.2 fastmatch_1.1-0  tools_4.0.4      igraph_1.2.6    
#[13] parallel_4.0.4   compiler_4.0.4   pkgconfig_2.0.3  tensorA_0.36.2  
#> 


   