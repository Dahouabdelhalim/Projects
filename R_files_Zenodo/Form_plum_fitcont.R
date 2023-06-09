###Fit continuous

#load packages
library(ape)
library(geiger)
library(nlme)
library(phytools)

####Males
plum.m<-read.csv("Form_plum_mean.2.0_m.csv",row.names = 1)
head(plum.m)


#remove columns not used (sex, RGB mean and sumPower)
plum.m<-plum.m[,-c(1,2:4,7,9:11,14,16:18,21)]
head(plum.m)

#log transformation
plum.m.ln<-log(plum.m)
head(plum.m.ln)


#Open tree
tree.m <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(tree.m)
tree.m

#Check tree and data
objm<-name.check(tree.m,plum.m.ln)
objm

#Remove tips that we don't have data
tree.m <- drop.tip(tree.m, objm$tree_not_data)
plot(tree.m)

##Remove data of taxa that we don't have in the tree 
plum.m.ln<-plum.m.ln[-c(25,54,56),]
plum.m.ln

#Check tree and data again
name.check(tree.m, plum.m.ln)



#2.4.3.2 Testing models 



#----------------------------------------------------

#   BROWNIAN MOTION

#---------------------------------------------------

brownFit <-  fitContinuous(tree.m,plum.m.ln)

aic.brown<-numeric(9) #9 plumage traits

for (i in 1:9) aic.brown[i]<-brownFit[[i]]$opt$aicc

lik.brown<-numeric(9)

for (i in 1:9) lik.brown[i]<-brownFit[[i]]$opt$lnL

par.brown<-numeric(9)

for (i in 1:9) par.brown[i]<-brownFit[[i]]$opt[1]


#----------------------------------------------------

#   OU MODEL: ALPHA

#---------------------------------------------------

ouFit<-fitContinuous(tree.m,plum.m.ln, model="OU")

aic.ou<-numeric(9)

for (i in 1:9) aic.ou[i]<-ouFit[[i]]$opt$aicc

lik.ou<-numeric(9)

for (i in 1:9) lik.ou[i]<-ouFit[[i]]$opt$lnL

par.ou<-numeric(9)

for (i in 1:9) par.ou[i]<-ouFit[[i]]$opt[1]


#----------------------------------------------------

#   Pagel lambda MODEL: LAMBDA

#---------------------------------------------------

lambdaFit<-fitContinuous(tree.m,plum.m.ln, model="lambda")

aic.lambda<-numeric(9)

for (i in 1:9) aic.lambda[i]<-lambdaFit[[i]]$opt$aicc

lik.lambda<-numeric(9)

for (i in 1:9) lik.lambda[i]<-lambdaFit[[i]]$opt$lnL

par.lambda<-numeric(9)

for (i in 1:9) par.lambda[i]<-lambdaFit[[i]]$opt[1]

#----------------------------------------------------

#   COMPARE ALL MODELS MULTIVARIATE

#---------------------------------------------------

for (i in 1:9)
  
{
  
  ## REMEMBER TO MODIFY FILE NAMES
  
  files<-c("Lum_dorsal_model_comparison_m2.txt", "MP_dorsal_model_comparison_m2.txt", "Cont_dorsal_model_comparison_m2.txt","Lum_ventral_model_comparison_m2.txt", "MP_ventral_model_comparison_m2.txt", "Cont_ventral_model_comparison_m2.txt","Lum_wing_model_comparison_m2.txt", "MP_wing_model_comparison_m2.txt", "Cont_wing_model_comparison_m2.txt")
  
  vars<-c("Lum_dorsal", "MP_dorsal", "Cont_dorsal","Lum_ventral", "MP_ventral", "Cont_ventral","Lum_wing", "MP_wing", "Cont_wing")
  
  liksVARS<-rbind(lik.brown[i], lik.ou[i], lik.lambda[i])
  
  paramsVARS<-rbind(par.brown[i], par.ou[i], par.lambda[i])
  
  aicVARS<-rbind(aic.brown[i], aic.ou[i], aic.lambda[i])
  
  foo<-function(x) x-x[which(x==min(x))]
  
  daic<-t(apply(aicVARS, 2, foo))
  
  tab<-cbind(paramsVARS,liksVARS,aicVARS,daic[1,])
  
  colnames(tab)<-c("Parameters","Likelihood","AICc","DeltaAICc")
  
  rownames(tab)<-c("Brownian", "OU", "Lambda")
  
  cat(vars[i], "\\n")
  
  print(tab, digits=9)
  
  cat("\\n")
  
  write.table(tab, file=files[i], sep=",")
  
}




######################################################################################
#-----------------------------------------------------------------------------------
######################################################################################
#-----------------------------------------------------------------------------------
######################################################################################
#-----------------------------------------------------------------------------------


#Females
plum.f<-read.csv("Form_plum_mean.2.0_f.csv",row.names = 1)
head(plum.f)

#remove columns not used (sex, RGB mean and sumPower)
plum.f<-plum.f[,-c(1,2:4,7,9:11,14,16:18,21)]
head(plum.f)

#log transformation
plum.f.ln<-log(plum.f)
head(plum.f.ln)


#Open tree
tree.f <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(tree.f)
tree.f

#Check tree and data
obj_f<-name.check(tree.f,plum.f.ln)
obj_f

#Remove tips that we don't have data
tree.f <- drop.tip(tree.f, obj_f$tree_not_data)
plot(tree.f)

##Remove data of taxa that we don't have in the tree   
plum.f.ln<-plum.f.ln[-c(25,53,55),]

#Check tree and data again
name.check(tree.f, plum.f.ln)


#########################################################
#2.4.3.2 Testing model of shape evolution



#----------------------------------------------------

#   BROWNIAN MOTION

#---------------------------------------------------

brownFit <-  fitContinuous(tree.f,plum.f.ln)

aic.brown<-numeric(9) #9 plumage traits

for (i in 1:9) aic.brown[i]<-brownFit[[i]]$opt$aicc

lik.brown<-numeric(9)

for (i in 1:9) lik.brown[i]<-brownFit[[i]]$opt$lnL

par.brown<-numeric(9)

for (i in 1:9) par.brown[i]<-brownFit[[i]]$opt[1]


#----------------------------------------------------

#   OU MODEL: ALPHA

#---------------------------------------------------

ouFit<-fitContinuous(tree.f,plum.f.ln, model="OU")

aic.ou<-numeric(9)

for (i in 1:9) aic.ou[i]<-ouFit[[i]]$opt$aicc

lik.ou<-numeric(9)

for (i in 1:9) lik.ou[i]<-ouFit[[i]]$opt$lnL

par.ou<-numeric(9)

for (i in 1:9) par.ou[i]<-ouFit[[i]]$opt[1]


#----------------------------------------------------

#   Pagel lambda MODEL: LAMBDA

#---------------------------------------------------

lambdaFit<-fitContinuous(tree.f,plum.f.ln, model="lambda")

aic.lambda<-numeric(9)

for (i in 1:9) aic.lambda[i]<-lambdaFit[[i]]$opt$aicc

lik.lambda<-numeric(9)

for (i in 1:9) lik.lambda[i]<-lambdaFit[[i]]$opt$lnL

par.lambda<-numeric(9)

for (i in 1:9) par.lambda[i]<-lambdaFit[[i]]$opt[1]


#----------------------------------------------------

#   COMPARE ALL MODELS MULTIVARIATE

#---------------------------------------------------

for (i in 1:9)
  
{
  
  ## REMEMBER TO MODIFY FILE NAMES
  
  files<-c("Lum_dorsal_model_comparison_f2.txt", "MP_dorsal_model_comparison_f2.txt", "Cont_dorsal_model_comparison_f2.txt","Lum_ventral_model_comparison_f2.txt", "MP_ventral_model_comparison_f2.txt", "Cont_ventral_model_comparison_f2.txt","Lum_wing_model_comparison_f2.txt", "MP_wing_model_comparison_f2.txt", "Cont_wing_model_comparison_f2.txt")
  
  vars<-c("Lum_dorsal", "MP_dorsal", "Cont_dorsal","Lum_ventral", "MP_ventral", "Cont_ventral","Lum_wing", "MP_wing", "Cont_wing")
  
  liksVARS<-rbind(lik.brown[i], lik.ou[i], lik.lambda[i])
  
  paramsVARS<-rbind(par.brown[i], par.ou[i], par.lambda[i])
  
  aicVARS<-rbind(aic.brown[i], aic.ou[i], aic.lambda[i])
  
  foo<-function(x) x-x[which(x==min(x))]
  
  daic<-t(apply(aicVARS, 2, foo))
  
  tab<-cbind(paramsVARS,liksVARS,aicVARS,daic[1,])
  
  colnames(tab)<-c("Parameters","Likelihood","AICc","DeltaAICc")
  
  rownames(tab)<-c("Brownian", "OU", "Lambda")
  
  cat(vars[i], "\\n")
  
  print(tab, digits=9)
  
  cat("\\n")
  
  write.table(tab, file=files[i], sep=",")
  
}

