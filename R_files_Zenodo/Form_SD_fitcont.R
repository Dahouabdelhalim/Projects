#load packages
library(ape)
library(geiger)
library(nlme)
library(phytools)

#Sexual dimorphism data
Form.SDeco<-read.csv("Form_SD_eco_9traits.csv",row.names = 1)
head(Form.SDeco)

#Remove columns not used (ecological and behavioral traits)
Form.SD<-Form.SDeco[,c(1,2)]
head(Form.SD)

library(tidyr)
#remove NAs
Form.SD<-drop_na(Form.SD)
head(Form.SD)


#Open tree
tree <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(tree)
tree

#Check tree and data
obj<-name.check(tree,Form.SD)
obj

##Remove tips that we don't have data
tree <- drop.tip(tree, obj$tree_not_data)
plot(tree)


#Check tree and data again
name.check(tree, Form.SD)


####################################################
#2.4.3.2 Testing model of shape evolution



#----------------------------------------------------

#   BROWNIAN MOTION

#---------------------------------------------------

brownFit <-  fitContinuous(tree,Form.SD)

aic.brown<-numeric(2) #2 traits

for (i in 1:2) aic.brown[i]<-brownFit[[i]]$opt$aicc

lik.brown<-numeric(2)

for (i in 1:2) lik.brown[i]<-brownFit[[i]]$opt$lnL

par.brown<-numeric(2)

for (i in 1:2) par.brown[i]<-brownFit[[i]]$opt[1]


#----------------------------------------------------

#   OU MODEL: ALPHA

#---------------------------------------------------

ouFit<-fitContinuous(tree,Form.SD, model="OU")

aic.ou<-numeric(2)

for (i in 1:2) aic.ou[i]<-ouFit[[i]]$opt$aicc

lik.ou<-numeric(2)

for (i in 1:2) lik.ou[i]<-ouFit[[i]]$opt$lnL

par.ou<-numeric(2)

for (i in 1:2) par.ou[i]<-ouFit[[i]]$opt[1]


#----------------------------------------------------

#   Pagel lambda MODEL: LAMBDA

#---------------------------------------------------

lambdaFit<-fitContinuous(tree,Form.SD, model="lambda")

aic.lambda<-numeric(2)

for (i in 1:2) aic.lambda[i]<-lambdaFit[[i]]$opt$aicc

lik.lambda<-numeric(2)

for (i in 1:2) lik.lambda[i]<-lambdaFit[[i]]$opt$lnL

par.lambda<-numeric(2)

for (i in 1:2) par.lambda[i]<-lambdaFit[[i]]$opt[1]


#----------------------------------------------------

#   COMPARE ALL MODELS MULTIVARIATE

#---------------------------------------------------

for (i in 1:2)
  
{
  
  ## REMEMBER TO MODIFY FILE NAMES
  
  files<-c("ED_plum2.txt", "ED_vocal2.txt")
  
  vars<-c("ED_plum",  "ED_vocal")
  
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
