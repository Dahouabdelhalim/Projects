#final hypothesis-testing models performed over a tree block of 1,000 parrot trees (see Tables
#2 and S6)

#Miquel Vall-Llosera kindly provided us with example code using a loop to run the analyses over a tree block 
#which we've since adapted for our own use

library(Rmisc)
library(ape)
library(caper)
library(rcompanion)
library(DescTools)
library(tidyverse)

#read trees and data in 
parrottrees <- read.nexus("parrot_trees.nex")
parrot_data<-read.csv("parrot_comp_data.csv", header=TRUE)
attach(parrot_data)

#whole body SB and maximum feeding group size model

#set up objects to store the values for each term 
t<- c(rep(0, length(parrottrees)))#t stat
p<- c(rep(0, length(parrottrees)))#p value
r2<- c(rep(0, length(parrottrees))) # r2
lambda <- c(rep(0, length(parrottrees))) #lambda
#run loop. For each of the 1,000 parrot trees, run a PGLS model based on the comparative object and store 
#in the object 'mods'
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(BSB+1) ~ log(Max_feed_size), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]#extract the values for each parameter in turn and place within the previously
    p[i]<-summary(mods)$coef[2,4]#created objects for each parameter
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

#to get overall n and dfs
summary(mods)

#getting median values and CIs
#t
MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)
#r2
MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)
#lambda
MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)
#p
MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#whole body SB and communal roosting

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(BSB+1) ~ Communal_roost+Current_EE+Early_EE, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#whole body SB and food search

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(BSB+1) ~ Food_search+Stand_cage+Cap_diet_div, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#whole body SB and food handling

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(BSB~ sqrt(Food_handling)+sqrt(Prop_adult)+ Prop_female+Stand_cage, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#whole body SB and diet breadth

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(BSB~ Diet_breadth, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#whole body SB and habitat breadth

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(BSB~ Habitat_breadth, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now innovation rate

#have to filter dataset first so that only species within innovation rate regions are included
parrot_data_inn<-filter(parrot_data, Inn_region==1)
detach(parrot_data)
attach(parrot_data_inn)

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data_inn <- comparative.data(parrottrees[[i]], parrot_data_inn, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(BSB)~ Innovation+log(Research_effort)+Human_reared, parrot_comp_data_inn, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now whole body SB and brain volume model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(BSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)


MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now whole body SB and IUCN model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(BSB ~ IUCN_code, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and maximum feeding group size model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB ~ log(Max_feed_size), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and communal roosting model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB ~ Communal_roost+Current_EE+Early_EE, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and food search model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB ~ Food_search+Stand_cage+Cap_diet_div, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and food handling model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(OSB) ~ sqrt(Food_handling)+Prop_adult+ Prop_female+Stand_cage, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and diet breadth model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB~ Diet_breadth, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and habitat breadth model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB~ Habitat_breadth, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and innovation rate 

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data_inn <- comparative.data(parrottrees[[i]], parrot_data_inn, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(OSB) ~ Innovation+log(Research_effort)+Human_reared, parrot_comp_data_inn, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now oral SB and brain volume model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#oral SB and IUCN model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(OSB~ IUCN_code, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#FDB and maximum feeding group size model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB) ~ log(Max_feed_size), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#FDB and communal roosting model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB) ~ Communal_roost+Early_EE+Current_EE, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#FDB and food search model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB) ~ Food_search+Stand_cage+Cap_diet_div, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now FDB and food handling model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(FDB~ sqrt(Food_handling)+Prop_adult+ Prop_female+Stand_cage, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now FDB and diet breadth model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB)~ Diet_breadth, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now FDB and habitat breadth model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB)~ Habitat_breadth, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now FDB and habitat breadth model (with maximum foraging group size). Note we want to extract t and p values 
#for a correlated predictor (max feeding group size) so we set up extra empty objects to store those parameters in

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
t2<- c(rep(0, length(parrottrees)))
p2<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(FDB+1) ~ Habitat_breadth+log(Max_feed_size), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    t2[i]<-summary(mods)$coef[3,3]
    p2[i]<-summary(mods)$coef[3,4]   
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    t2[i]<-NA
    p2[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(t2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#FDB and innovation rate model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data_inn <- comparative.data(parrottrees[[i]], parrot_data_inn, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB)~ Innovation+log(Research_effort)+Human_reared, parrot_comp_data_inn, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now FDB and brain volume model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(FDB~ log(Brain_vol)+log(Body_mass)+log(Stand_cage+1), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now FDB and IUCN model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(sqrt(FDB)~ IUCN_code, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)


#hatch rate and maximum feeding group size model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~ log(Max_feed_size)+Nat_fecund, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now for hatch rate and communal roosting model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~ Communal_roost+Nat_fecund, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now for hatch rate and food search model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1)~sqrt(Food_search)+log(Nat_fecund), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#hatch rate and food search model (with brain volume). Note we want to extract t and p values for a 
#correlated predictor (brain vol) so we set up extra empty objects to store those parameters in

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
t2<- c(rep(0, length(parrottrees)))
p2<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~ sqrt(Food_search)+log(Brain_vol)+log(Nat_fecund)+log(Body_mass), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    t2[i]<-summary(mods)$coef[3,3]
    p2[i]<-summary(mods)$coef[3,4]   
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    t2[i]<-NA
    p2[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(t2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now for hatch rate and food handling model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1)~sqrt(Food_handling)+log(Nat_fecund), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now for hatch rate and diet breadth model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1)~ Diet_breadth+log(Nat_fecund), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now for hatch rate and habitat breadth model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1)~ Habitat_breadth+Nat_fecund, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now hatch rate and innovation rate

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data_inn <- comparative.data(parrottrees[[i]], parrot_data_inn, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~log(Innovation+1) +Nat_fecund+ log(Research_effort), parrot_comp_data_inn, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now hatch rate and brain volume model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~ log(Brain_vol)+log(Nat_fecund)+ log(Body_mass), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)


#now hatch rate and IUCN model

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~IUCN_code+Nat_fecund, parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now hatch rate, IUCN and habitat breadth model. Note we want to extract t and p values for a correlated predictor (max 
#habitat breadth) so we set up extra empty objects to store those parameters in

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
t2<- c(rep(0, length(parrottrees)))
p2<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1)~ IUCN_code+Habitat_breadth+log(Nat_fecund), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    t2[i]<-summary(mods)$coef[3,3]
    p2[i]<-summary(mods)$coef[3,4]   
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    t2[i]<-NA
    p2[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(t2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

#now hatch rate, IUCN and maximum foraging group size model. Note we want to extract t and p values for a correlated predictor (max 
#maximum foraging group size) so we set up extra empty objects to store those parameters in

t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
t2<- c(rep(0, length(parrottrees)))
p2<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees))) 
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~ IUCN_code+log(Max_feed_size)+log(Nat_fecund), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    t2[i]<-summary(mods)$coef[3,3]
    p2[i]<-summary(mods)$coef[3,4]   
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    t2[i]<-NA
    p2[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(t2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)


#now for hatch rates and n breeding pairs model
t<- c(rep(0, length(parrottrees)))
p<- c(rep(0, length(parrottrees)))
r2<- c(rep(0, length(parrottrees)))
lambda <- c(rep(0, length(parrottrees)))
for (i in 1:length(parrottrees)) {
  parrot_comp_data <- comparative.data(parrottrees[[i]], parrot_data, Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
  mods<- try(pgls(log(Hatch_rate+1) ~ log(Hatch_n_breed_pairs)+log(Nat_fecund), parrot_comp_data, lambda = "ML"), silent = TRUE)
  if(class(mods) != "try-error") {
    t[i]<-summary(mods)$coef[2,3]
    p[i]<-summary(mods)$coef[2,4]
    r2[i]<-summary(mods)$r.squared
    lambda[i]<-summary(mods)$param[2]
  }
  else {
    mods <- NA
    t[i]<-NA
    p[i]<-NA
    r2[i]<-NA
    lambda[i]<-NA
  }}

summary(mods)

MedianCI(t,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(r2,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(lambda,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

MedianCI(p,
         conf.level = 0.95,
         na.rm = TRUE,
         method = "exact",
         R = 10000)

