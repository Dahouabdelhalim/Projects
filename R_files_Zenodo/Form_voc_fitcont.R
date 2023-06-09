###Fit continuous


#load packages
library(ape)
library(geiger)
library(nlme)
library(phytools)


####Males
voc.m<-read.csv("Form_vocal_mean_m_ok.csv",row.names = 1)
head(voc.m)
str(voc.m)

#log transformation and remove columns not used (sex and entropy)
voc.m.ln<-log(voc.m[-c(1,4,5,6,7)])
head(voc.m.ln)


#Open tree
tree.m <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(tree.m)
tree.m

#Check tree and data
objm<-name.check(tree.m,voc.m.ln)
objm

#Remove tips that we don't have data
tree.m <- drop.tip(tree.m, objm$tree_not_data)
plot(tree.m)

##Remove data of taxa that we don't have in the tree 
voc.m.ln<-voc.m.ln[-c(23,51,53),]
head(voc.m.ln)

#Check tree and data again
name.check(tree.m, voc.m.ln)



#2.4.3.2 Testing model of shape evolution



#----------------------------------------------------

#   BROWNIAN MOTION

#---------------------------------------------------

brownFit <-  fitContinuous(tree.m,voc.m.ln)

aic.brown<-numeric(8) #8 vocal traits

for (i in 1:8) aic.brown[i]<-brownFit[[i]]$opt$aicc

lik.brown<-numeric(8)

for (i in 1:8) lik.brown[i]<-brownFit[[i]]$opt$lnL

par.brown<-numeric(8)

for (i in 1:8) par.brown[i]<-brownFit[[i]]$opt[1]


#----------------------------------------------------

#   OU MODEL: ALPHA

#---------------------------------------------------

ouFit<-fitContinuous(tree.m,voc.m.ln, model="OU")

aic.ou<-numeric(8)

for (i in 1:8) aic.ou[i]<-ouFit[[i]]$opt$aicc

lik.ou<-numeric(8)

for (i in 1:8) lik.ou[i]<-ouFit[[i]]$opt$lnL

par.ou<-numeric(8)

for (i in 1:8) par.ou[i]<-ouFit[[i]]$opt[1]


#----------------------------------------------------

#   Pagel lambda MODEL: LAMBDA

#---------------------------------------------------

lambdaFit<-fitContinuous(tree.m,voc.m.ln, model="lambda")

aic.lambda<-numeric(8)

for (i in 1:8) aic.lambda[i]<-lambdaFit[[i]]$opt$aicc

lik.lambda<-numeric(8)

for (i in 1:8) lik.lambda[i]<-lambdaFit[[i]]$opt$lnL

par.lambda<-numeric(8)

for (i in 1:8) par.lambda[i]<-lambdaFit[[i]]$opt[1]


#----------------------------------------------------

#   COMPARE ALL MODELS MULTIVARIATE

#---------------------------------------------------

for (i in 1:8)
  
{
  
  ## REMEMBER TO MODIFY FILE NAMES
  
  files<-c("Notecount_model_comparison_m2.txt", "Notetype_model_comparison_m2.txt", "Duration_model_comparison_m2.txt","Peakfreq_model_comparison_m2.txt", "Notediversity_model_comparison_m2.txt", "Noterate_model_comparison_m2.txt","Songband_model_comparison_m2.txt", "Loudsongmodulationrate_model_comparison_m2.txt")
  
  vars<-c("Notecount", "Notetype", "Duration","Peakfreq", "Notediversity", "Noterate","Songband", "Loudsongmodulationrate")
  
  liksVARS<-rbind(lik.brown[i], lik.ou[i], lik.lambda[i])
  
  paramsVARS<-rbind(par.brown[i], par.ou[i], par.lambda[i])
  
  aicVARS<-rbind(aic.brown[i], aic.ou[i], aic.lambda[i])
  
  foo<-function(x) x-x[which(x==min(x))]
  
  daic<-t(apply(aicVARS, 2, foo))
  
  tab<-cbind(paramsVARS,liksVARS,aicVARS,daic[1,])
  
  colnames(tab)<-c("Parameters","Likelihood","AICc","DeltaAICc")
  
  rownames(tab)<-c("Brownian", "OU", "Lambda")
  
  cat(vars[i], "\\n")
  
  print(tab, digits=8)
  
  cat("\\n")
  
  write.table(tab, file=files[i], sep=",")
  
}


#####################################################################################################################################################################################################################################################################


####Females
voc.f<-read.csv("Form_vocal_mean_f_ok.csv",row.names = 1)
head(voc.f)
str(voc.f)


#log transformation and remove columns not used (sex and entropy)
voc.f.ln<-log(voc.f[-c(1,4:6,7)])
head(voc.f.ln)

#Open tree
tree.f <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(tree.f)
tree.f

#Check tree and data
objf<-name.check(tree.f,voc.f.ln)
objf


#Remove tips that we don't have data
tree.f <- drop.tip(tree.f, objf$tree_not_data)
plot(tree.f)

#Check tree and data again
name.check(tree.f,voc.f.ln)



#2.4.3.2 Testing model of shape evolution



#----------------------------------------------------

#   BROWNIAN MOTION

#---------------------------------------------------

brownFit <-  fitContinuous(tree.f,voc.f.ln)

aic.brown<-numeric(8) #8 vocal traits

for (i in 1:8) aic.brown[i]<-brownFit[[i]]$opt$aicc

lik.brown<-numeric(8)

for (i in 1:8) lik.brown[i]<-brownFit[[i]]$opt$lnL

par.brown<-numeric(8)

for (i in 1:8) par.brown[i]<-brownFit[[i]]$opt[1]


#----------------------------------------------------

#   OU MODEL: ALPHA

#---------------------------------------------------

ouFit<-fitContinuous(tree.f,voc.f.ln, model="OU")

aic.ou<-numeric(8)

for (i in 1:8) aic.ou[i]<-ouFit[[i]]$opt$aicc

lik.ou<-numeric(8)

for (i in 1:8) lik.ou[i]<-ouFit[[i]]$opt$lnL

par.ou<-numeric(8)

for (i in 1:8) par.ou[i]<-ouFit[[i]]$opt[1]


#----------------------------------------------------

#   Pagel lambda MODEL: LAMBDA

#---------------------------------------------------

lambdaFit<-fitContinuous(tree.f,voc.f.ln, model="lambda")

aic.lambda<-numeric(8)

for (i in 1:8) aic.lambda[i]<-lambdaFit[[i]]$opt$aicc

lik.lambda<-numeric(8)

for (i in 1:8) lik.lambda[i]<-lambdaFit[[i]]$opt$lnL

par.lambda<-numeric(8)

for (i in 1:8) par.lambda[i]<-lambdaFit[[i]]$opt[1]

#----------------------------------------------------

#   COMPARE ALL MODELS MULTIVARIATE

#---------------------------------------------------

for (i in 1:8)
  
{
  
  ## REMEMBER TO MODIFY FILE NAMES
  
  files<-c("Notecount_model_comparison_f2.txt", "Notetype_model_comparison_f2.txt", "Duration_model_comparison_f2.txt","Peakfreq_model_comparison_f2.txt", "Notediversity_model_comparison_f2.txt", "Noterate_model_comparison_f2.txt","Songband_model_comparison_f2.txt", "Loudsongmodulationrate_model_comparison_f2.txt")
  
  vars<-c("Notecount", "Notetype", "Duration","Peakfreq", "Notediversity", "Noterate","Songband", "Loudsongmodulationrate")
  
  liksVARS<-rbind(lik.brown[i], lik.ou[i], lik.lambda[i])
  
  paramsVARS<-rbind(par.brown[i], par.ou[i], par.lambda[i])
  
  aicVARS<-rbind(aic.brown[i], aic.ou[i], aic.lambda[i])
  
  foo<-function(x) x-x[which(x==min(x))]
  
  daic<-t(apply(aicVARS, 2, foo))
  
  tab<-cbind(paramsVARS,liksVARS,aicVARS,daic[1,])
  
  colnames(tab)<-c("Parameters","Likelihood","AICc","DeltaAICc")
  
  rownames(tab)<-c("Brownian", "OU", "Lambda")
  
  cat(vars[i], "\\n")
  
  print(tab, digits=8)
  
  cat("\\n")
  
  write.table(tab, file=files[i], sep=",")
  
}
