
##############################################################
# Author: 
# Antica Culina
# NIOO-KNAW & University of Oxford 
# Email: a.culina@yahoo.com


##############################################################
# Description of script and Instructions
##############################################################

# This code is used to conduct publication bias tests, and produce Supplementary
# Figure S2

# Set working directory to wherever you have downloaded the data and code
#session info  can be found at the end of this code

rm(list=ls())


# load libraries

library(metafor)
library(rmeta)
library(ape)
library(phangorn)
library(MCMCglmm)

sessionInfo()

####### variables ###########
# Z = Fishers Zr, effect size used in meta-anlaysis
#SE = standard error of the effect size
# title = represents the study, random effect
# latin_name = represents species, random effect



#####################################
####trim and fill in the methaphor package, plus constructing funnel plots
####################################3


######BEFORE DATASET

before<- read.table(file.choose(), header = T, sep = ";")
before$year <- as.numeric(scale(before$year, scale = FALSE))

##### timelag bias

year <- MCMCglmm (Z~year, random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=before)

summary(year)

# first, do funnel plot and run the basic model on all the effect sizes
# this funnel plot corresponds to Figure S2, panel B)

res<- rma(Z, SE^2, method="REML", data=before)
summary(res)


funnel(res, yaxis="seinv", xlab="Effect size (Zr)", ylim=c(0.01,33), xlim=c(-1.1,1.1), cex=1, digits=1,  back="white", shade="white", hlines="white")

before$latin_name <- as.factor(before$latin_name)

points(x= before$Z, y=1/before$SE, pch = 16, cex = 1.5, col = c("blue", "blueviolet", "darkgoldenrod1", "black", "aquamarine4", "red", "brown", "aquamarine1", "gray", "chartreuse")[before$latin_name])
trimfill(res, maxit=50)


#do this but with randomly chosen effect size from each study

beforemiss.studies <-beforemiss.side <- beforemiss.b <- array(NA, dim=2000)

beforeall.output <- list()

for(i in 1:2000){

  uni <- sapply(unique(before$title), function(x) sample (which(before$title == x), 1))
  before2 <-before[ uni , ]
ressubset<- rma(Z, SE^2, method="REML", data=before2)
#summary(ressubset)
#funnel(ressubset)
beforetrim <-trimfill(ressubset, maxit=50)
#str(beforetrim)

beforemiss.studies[i] <-beforetrim$k0

beforemiss.side[i] <- beforetrim$side
beforemiss.b[i]<- beforetrim$b
beforeall.output[[i]] <-beforetrim }

beforetrim

# this histogram corresponds to Figure S2, panel A)

hist(ifelse(beforemiss.side=="right",beforemiss.studies,-beforemiss.studies), xlab='estimated number of missing studies', main='')
title(main='Histogram of estimated number of missing studies for r_before')
mean(beforemiss.studies)
sd(beforemiss.studies)
mean(beforemiss.b)   #### mean of the all estimates of the Effect size
sd(beforemiss.b)    #### sd of estimates of the Effect sizes


##### AFTER DATASET

after<- read.table(file.choose(), header = T, sep = ";")

# time-lag bias

after$year <- as.numeric(scale(after$year, scale = FALSE))

m.year <- MCMCglmm (Z~year, random=~latin_name, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, thin = 1000, prior=prior1, data=after)

summary(m.year)

# first, do funnel plot and run the basic model on all the effect sizes

afres<- rma(Z, SE^2, method="REML", data=after)

summary(afres)
after$latin_name <- as.factor(after$latin_name)

unique(after$latin_name)

# this funnel plot corresponds to Figure S2, panel D)

funnel(afres, yaxis="seinv", xlab="Effect size (Zr)", ylim=c(0.01,33), xlim=c(-1.1,1.1), digits=1,  back="white", shade="white", hlines="white",
       pch = 16, cex = 1.5, col = c("darkolivegreen4", "darkolivegreen2", "darkorange", "darkorange4", "black", "darkorchid", "darkred", "darksalmon", "deepskyblue2", "gold", "deeppink", "darkgreen", "blueviolet", "blue4")[after$latin_name])
trimfill(afres, maxit=50)



#afres2 <- rma.mv(Z, SE, random = ~ 1 | title, data=after)
#funnel(afres2, yaxis="seinv", xlab="Effect size (Zr)", ylim=c(0.01,20), xlim=c(-1.1,1.1), cex=0.75, digits=1,  back="white", shade="white", hlines="white")

#do this but with randomly choosen effect size from each study

miss.studies <- miss.side <- miss.b <- array(NA, dim=2000)

all.output <- list()

for(i in 1:1000){

  uni <- sapply(unique(after$title), function(x) sample (which(after$title == x), 1))
  after2 <-after[ uni , ]


afressubset<- rma(Z, SE^2, method="REML", data=after2, control=list(maxiter=2000))

trim <- trimfill(afressubset, maxit=50)
miss.studies[i] <- trim$k0
miss.b[i] <- trim$b
miss.side[i] <- trim$side
all.output[[i]] <- trim
}

# did not converge

# this histogram corresponds to Figure S2, panel C)

hist(ifelse(miss.side=="right",miss.studies,-miss.studies), xlab='estimated number of missing studies', main='Histogram of estimated number of missing studies for r')
mean(miss.studies)
median(miss.studies)
sd(miss.studies)
mean(miss.b)
sd(miss.b)


####### sexes (F vs M) dataset


FvsM<- read.table(file.choose(), header = T, sep = ";")

# time-lag bias

FvsM$year <- as.numeric(scale(FvsM$year, scale = FALSE))

m.year <- MCMCglmm (Z~year, random=~title, mev=FvsM$SE^2, burnin=200000, nitt=2000000, thin=1000, prior=prior1, data=FvsM)
summary(m.year)

# first, do funnel plot and run the basic model on all the effect sizes

FvsMres<- rma(Z, SE^2, method="REML", data=FvsM)
summary(FvsMres)

unique(FvsM$latin_name)
FvsM$latin_name <- as.factor(FvsM$latin_name)

# this funnel plot corresponds to Figure S2, panel F)

funnel(FvsMres, yaxis="seinv", xlab="Effect size (Zr)", ylim=c(0.01,33), xlim=c(-1.1,1.1), 
       pch = 16, cex = 1.5, col = c("blue", "green", "yellow", "black", "orange", "red", "brown", "gray", "darkorchid")[FvsM$latin_name], digits=1,  back="white", shade="white", hlines="white")
trimfill(FvsMres, maxit=50)

#do this but with randomly choosen effect size from each study

FvsMmiss.studies <-FvsMmiss.side <- FvsMmiss.b <- array(NA, dim=2000)

FvsMall.output <- list()


for(i in 1:2000){
  
  uni <- sapply(unique(FvsM$title), function(x) sample (which(FvsM$title == x), 1))
  FvsM2 <-FvsM[ uni , ]
  
  FvsMsubset<- rma(Z, SE^2, method="REML", data=FvsM2)
  FvsMtrim <- trimfill(FvsMsubset)
  
  FvsMmiss.studies[i] <-FvsMtrim$k0
  FvsMmiss.b[i] <- FvsMtrim$b
  FvsMmiss.side[i] <- FvsMtrim$side
  FvsMall.output[[i]] <-FvsMtrim 
}

# did not converge

# this histogram corresponds to Figure S2, panel E)

hist(ifelse(FvsMmiss.side=="right",FvsMmiss.studies,-FvsMmiss.studies), xlab='estimated number of missing studies', main='Histogram of estimated number of missing studies for r')
mean(FvsMmiss.studies)
sd(FvsMmiss.studies)
mean(FvsMmiss.b)
sd(FvsMmiss.b)


#############################################
###########FUNEL PLOT ASSYMETRY 
#############################################


#prior for the models

prior1 <- list(G=list(G1=list(V=1, nu=0.002)), R=list(V=1, nu=0.002))

#loading  and adjusting phylogenetic tree

tree <- read.nexus("F:\\\\NIOO3\\\\My Projects\\\\within season divorce\\\\wtihin_season divorce\\\\Anlysis\\\\within codes\\\\Trees_for_all_species\\\\all.nex")

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

Ainv<-inverseA(ult.nnls)$Ainv

# regressions for all three MA, keeping the best random effect structure

beforeegger <- MCMCglmm(Z~1+(SE^2), random=~latin_name, mev=before$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, prior=prior1, data=before)
summary(beforeegger)


aftereger <- MCMCglmm(Z~1+SE^2, random=~latin_name, pr=TRUE, mev=after$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, prior=prior1, data=after)
summary(aftereger)



FvsMeger <- MCMCglmm (Z~1+SE^2, random=~title, mev=FvsM$SE^2, ginverse=list(latin_name=Ainv), burnin=200000, nitt=2000000, prior=prior1, data=FvsM)
summary(FvsMeger)

#========================================================
# session info
#========================================================

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
#  [1] rmeta_3.0      metafor_2.4-0  phangorn_2.7.0 MCMCglmm_2.32  ape_5.5        coda_0.19-4    Matrix_1.3-2  

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.6       codetools_0.2-18 quadprog_1.5-8   lattice_0.20-41  corpcor_1.6.9    grid_4.0.4      
#[7] nlme_3.1-152     magrittr_2.0.1   cubature_2.0.4.2 fastmatch_1.1-0  tools_4.0.4      igraph_1.2.6    
#[13] parallel_4.0.4   compiler_4.0.4   pkgconfig_2.0.3  tensorA_0.36.2  
