### Code for Babar Manuscript ###

noldbabs<-read.csv("~/Babar Experiment Data - RSOS.csv", header = T)
dim(noldbabs)
summary(noldbabs)
names(noldbabs)

group.membs<-table(noldbabs$model,noldbabs$pb.group)
group.membs       

library(ggplot2) # for graphs
library(pastecs) # for descriptives
library(lme4) 
library(lmerTest)
library(lmtest) # for bptest for heteroscedasity in models 
library(MuMIn)

###### HOUSE KEEPING ##########
noldbabs$focal.id<-as.factor(noldbabs$focal.id)    
noldbabs$rings<-as.factor(noldbabs$rings)        
noldbabs$experimenter<-as.factor(noldbabs$experimenter)  
noldbabs$pres.num<-as.factor(noldbabs$pres.num)  
noldbabs$model<-as.factor(noldbabs$model)   
noldbabs$pb.group<-as.factor(noldbabs$pb.group)  
noldbabs$pb.id<-as.factor(noldbabs$pb.id)   
noldbabs$out.time<-as.numeric(noldbabs$out.time)     
noldbabs$duration<-as.numeric(noldbabs$duration)      
noldbabs$ground.perch<-as.factor(noldbabs$ground.perch)  
noldbabs$focal.scold<-as.factor(noldbabs$focal.scold)   
noldbabs$focal.flight<-as.factor(noldbabs$focal.flight)  
noldbabs$flight.latency<-as.numeric(noldbabs$flight.latency)   
noldbabs$pr20.fl.num<-as.numeric(noldbabs$pr20.fl.num) 
noldbabs$pr40.fl.num<-as.numeric(noldbabs$pr40.fl.num) 
noldbabs$pr60.fl.num<-as.numeric(noldbabs$pr60.fl.num) 
noldbabs$po20.fl.num<-as.numeric(noldbabs$po20.fl.num) 
noldbabs$po40.fl.num<-as.numeric(noldbabs$po40.fl.num) 
noldbabs$po60.fl.num<-as.numeric(noldbabs$po60.fl.num) 
noldbabs$pr20.fl.dist<-as.numeric(noldbabs$pr20.fl.dist) 
noldbabs$pr40.fl.dist<-as.numeric(noldbabs$pr40.fl.dist) 
noldbabs$pr60.fl.dist<-as.numeric(noldbabs$pr60.fl.dist) 
noldbabs$po20.fl.dist<-as.numeric(noldbabs$po20.fl.dist) 
noldbabs$po40.fl.dist<-as.numeric(noldbabs$po40.fl.dist) 
noldbabs$po60.fl.dist<-as.numeric(noldbabs$po60.fl.dist) 
noldbabs$box.chick<-as.factor(noldbabs$box.chick)  
noldbabs$days.fledged<-as.numeric(noldbabs$days.fledged)  

# DISTANCE #
noldbabs$prdist60<-noldbabs$pr20.fl.dist+noldbabs$pr40.fl.dist+noldbabs$pr60.fl.dist
noldbabs$prdist40<-noldbabs$pr40.fl.dist+noldbabs$pr20.fl.dist
noldbabs$prdist20<-noldbabs$pr20.fl.dist
noldbabs$podist20<-noldbabs$po20.fl.dist
noldbabs$podist40<-noldbabs$po20.fl.dist+noldbabs$po40.fl.dist
noldbabs$podist60<-noldbabs$po20.fl.dist+noldbabs$po40.fl.dist+noldbabs$po60.fl.dist
#declare them
noldbabs$prdist60<-as.numeric(noldbabs$prdist60)
noldbabs$prdist40<-as.numeric(noldbabs$prdist40)
noldbabs$prdist20<-as.numeric(noldbabs$prdist20)
noldbabs$podist20<-as.numeric(noldbabs$podist20)
noldbabs$podist40<-as.numeric(noldbabs$podist40)
noldbabs$podist60<-as.numeric(noldbabs$podist60)

### Flight Number
noldbabs$n.pre60<-noldbabs$pr20.fl.num+noldbabs$pr40.fl.num+noldbabs$pr60.fl.num
noldbabs$n.pre40<-noldbabs$pr40.fl.num+noldbabs$pr20.fl.num
noldbabs$n.pre20<-noldbabs$pr20.fl.num
noldbabs$n.post20<-noldbabs$po20.fl.num
noldbabs$n.post40<-noldbabs$po20.fl.num+noldbabs$po40.fl.num
noldbabs$n.post60<-noldbabs$po20.fl.num+noldbabs$po40.fl.num+noldbabs$po60.fl.num
#declare them
noldbabs$n.pre60<-as.numeric(noldbabs$n.pre60)
noldbabs$n.pre40<-as.numeric(noldbabs$n.pre40)
noldbabs$n.pre20<-as.numeric(noldbabs$n.pre20)
noldbabs$n.post20<-as.numeric(noldbabs$n.post20)
noldbabs$n.post40<-as.numeric(noldbabs$n.post40)
noldbabs$n.post60<-as.numeric(noldbabs$n.post60)

# TRANSFORMATIONS #
noldbabs$sqrt.post60<-sqrt(noldbabs$n.post60)
noldbabs$sqrt.post60<-as.numeric(noldbabs$sqrt.post60)
hist(noldbabs$n.pre60,100)
noldbabs$sqrt.pre60<-sqrt(noldbabs$n.pre60)
noldbabs$sqrt.pre60<-as.numeric(noldbabs$sqrt.pre60)
hist(noldbabs$sqrt.pre60,100) 

# Change model name from Babar to Elephant for use in manuscript graphs (one must look mature & responsible)
noldbabs$model2 [noldbabs$model=="BABAR"]<- "Elephant"  
noldbabs$model2 [noldbabs$model=="FOX"]<- "Fox" 
noldbabs$model2<-as.factor(noldbabs$model2)

# Make a big 4 level variable for ease of plotting
noldbabs$pb.model [noldbabs$model2=="Elephant" & noldbabs$pb.group=="CONTACT"]<- "Eleph.Con"  
noldbabs$pb.model [noldbabs$model2=="Fox" & noldbabs$pb.group=="CONTACT"]<- "Fox.Con"
noldbabs$pb.model [noldbabs$model2=="Fox" & noldbabs$pb.group=="SCOLD"]<- "Fox.Scold"  
noldbabs$pb.model [noldbabs$model2=="Elephant" & noldbabs$pb.group=="SCOLD"]<- "Eleph.Scold"
noldbabs$pb.model<-as.factor(noldbabs$pb.model)
summary(noldbabs$pb.model)
# Create subsets for later use
babsco<-subset(noldbabs, model=="BABAR" & pb.group=="SCOLD")
dim(babsco)
babcon<-subset(noldbabs, model=="BABAR" & pb.group=="CONTACT")
dim(babcon)
foxsco<-subset(noldbabs, model=="FOX" & pb.group=="SCOLD")
dim(foxsco)
foxcon<-subset(noldbabs, model=="FOX" & pb.group=="CONTACT")
dim(foxcon)

babs<-subset(noldbabs,model=="BABAR")
foxy<-subset(noldbabs,model=="FOX")

scold<-subset(noldbabs,pb.group=="SCOLD")
contact<-subset(noldbabs,pb.group=="CONTACT")

# Create subsets of the the different presentation periods
pres1<-subset(noldbabs,pres.num==1)
summary(pres1)
pres2<-subset(noldbabs,pres.num==2)
summary(pres2)
pres3<-subset(noldbabs,pres.num==3)
summary(pres3)

babs1<-subset(babs,pres.num==1)
summary(babs1)
babs2<-subset(babs,pres.num==2)
summary(babs2)
babs3<-subset(babs,pres.num==3)
summary(babs3)
foxy1<-subset(foxy,pres.num==1)
summary(foxy1)
foxy2<-subset(foxy,pres.num==2)
summary(foxy2)
foxy3<-subset(foxy,pres.num==3)
summary(foxy3)

pres1con<-subset(pres1,pb.group=="CONTACT")
summary(pres1con)
pres2con<-subset(pres2,pb.group=="CONTACT")
summary(pres1con)
pres3con<-subset(pres3,pb.group=="CONTACT")
summary(pres3con)
pres1sco<-subset(pres1,pb.group=="SCOLD")
summary(pres1sco)
pres2sco<-subset(pres2,pb.group=="SCOLD")
summary(pres2sco)
pres3sco<-subset(pres3,pb.group=="SCOLD")
summary(pres3sco)

P1BCo<-subset(babs1,pb.group=="CONTACT")
P1BSc<-subset(babs1,pb.group=="SCOLD")
P3BCo<-subset(babs3,pb.group=="CONTACT")
P3BSc<-subset(babs3,pb.group=="SCOLD")
P1FCo<-subset(foxy1,pb.group=="CONTACT")
P1FSc<-subset(foxy1,pb.group=="SCOLD")
P3FCo<-subset(foxy3,pb.group=="CONTACT")
P3FSc<-subset(foxy3,pb.group=="SCOLD")

# How correlated where the 2 potential response variables measured
cor(noldbabs$n.post60, noldbabs$podist60)
plot(noldbabs$n.post60, noldbabs$podist60) # very highly correlated. Use n.flights, as this is less subjective

# Descriptive Stats
### Presentation 1 flight number ### 
#Subsets & descriptives
pres1.babs<-subset(pres1,model=="BABAR")
p1ScB<-subset(pres1.babs,pb.group=="SCOLD")
p1CoB<-subset(pres1.babs,pb.group=="CONTACT")
stat.desc(p1ScB$n.pre60)  #mean Babs Scold   n.pre60 (SE) = 3.33 (1.36)
stat.desc(p1CoB$n.pre60)  #mean babs Contact n.pre60 (SE) = 5.33 (1.64)
stat.desc(p1ScB$n.post60)  #mean Babs Scold   n.pre60 (SE) = 9.83 (1.92)
stat.desc(p1CoB$n.post60)  #mean babs Contact n.pre60 (SE) = 9.67 (2.16)


pres1.fox<-subset(pres1,model=="FOX")
p1ScF<-subset(pres1.fox,pb.group=="SCOLD")
p1CoF<-subset(pres1.fox,pb.group=="CONTACT")
stat.desc(p1ScF$n.pre60) # mean Fox Scold   n.pre60 (SE) = 6.17 (1.59)
stat.desc(p1CoF$n.pre60) # mean Fox Contact n.pre60 (SE) = 4.33 (1.41)
stat.desc(p1ScF$n.post60) # mean Fox Scold   n.post60 (SE) = 11.08 (2.80)
stat.desc(p1CoF$n.post60) # mean Fox Contact n.post60 (SE) = 8.25 (1.56)

### Presentation 2 flight number ### 
#Subsets & descriptives
pres2.babs<-subset(pres2,model=="BABAR")
p2ScB<-subset(pres2.babs,pb.group=="SCOLD")
p2CoB<-subset(pres2.babs,pb.group=="CONTACT")
stat.desc(p2ScB$n.pre60) # mean Babar Scold    n.pre60 (SE) = 3.83 (1.54)
stat.desc(p2CoB$n.pre60) # mean Babar Contact  n.pre60 (SE) = 4.58 (1.70)
stat.desc(p2ScB$n.post60) # mean Babar Scold   n.post60 (SE) = 10.75 (2.10)
stat.desc(p2CoB$n.post60) # mean Babar Contact n.post60 (SE) = 9.67 (2.35)

pres2.fox<-subset(pres2,model=="FOX")
p2ScF<-subset(pres2.fox,pb.group=="SCOLD")
p2CoF<-subset(pres2.fox,pb.group=="CONTACT")
stat.desc(p2ScF$n.pre60)  # mean Fox Scold   n.pre60 (SE) = 8.33 (1.93)
stat.desc(p2CoF$n.pre60)  # mean Fox Contact n.pre60 (SE) = 4.92 (2.15)
stat.desc(p2ScF$n.post60) # mean Fox Scold   n.post60 (SE) = 12.67 (2.42)
stat.desc(p2CoF$n.post60) # mean Fox Contact n.post60 (SE) = 5.50 (1.03)

### Presentation 3 flight number plots ### 
#Subsets & descriptives
pres3.babs<-subset(pres3,model=="BABAR")
p3ScB<-subset(pres3.babs,pb.group=="SCOLD")
p3CoB<-subset(pres3.babs,pb.group=="CONTACT")
stat.desc(p3ScB$n.pre60)  # mean Babar Scold    n.pre60 (SE) = 6.67 (1.61)
stat.desc(p3CoB$n.pre60)  # mean Babar Contact  n.pre60 (SE) = 5.25 (2.05)
stat.desc(p3ScB$n.post60) # mean Babar Scold   n.post60 (SE) = 7.42 (1.64)
stat.desc(p3CoB$n.post60) # mean Babar Contact n.post60 (SE) = 4.08 (1.90)

pres3.fox<-subset(pres3,model=="FOX")
p3ScF<-subset(pres3.fox,pb.group=="SCOLD")
p3CoF<-subset(pres3.fox,pb.group=="CONTACT")
stat.desc(p3ScF$n.pre60)  # mean Fox Scold   n.pre60 (SE) = 9.00 (2.24)
stat.desc(p3CoF$n.pre60)  # mean Fox Contact n.pre60 (SE) = 4.17 (1.19)
stat.desc(p3ScF$n.post60) # mean Fox Scold   n.post60 (SE) = 11.58 (2.18)
stat.desc(p3CoF$n.post60) # mean Fox Contact n.post60 (SE) = 7.33 (1.65)

# more general numbers
stat.desc(pres1.babs$n.post60) # Babs Pres1 mean n.post 60 (SE) = 9.75 (1.41)
stat.desc(pres1.fox$n.post60)  # Fox Pres1 mean n.post 60 (SE) = 9.66 (1.59)
stat.desc(pres3.babs$n.post60) # Babs Pres3 mean n.post 60 (SE) = 5.75 (1.28)
stat.desc(pres3.fox$n.post60)  # Fox Pres3 mean n.post 60 (SE) = 9.46 (1.41)

stat.desc(pres3sco$n.post60) # pres 3 scold mean (SE) = 9.5 (1.4)
stat.desc(pres3con$n.post60) # pres 3 contact mean (SE) = 5.7 (1.3)

stat.desc(scold$n.post60)
stat.desc(contact$n.post60)



### Analyses ###
# Is it possible to analyse the data with a poisson error structure
poiss1<-glmer(n.post60~(model*pb.group*pres.num)+(1|focal.id),family=poisson,control=glmerControl(optimizer="bobyqa"),data=noldbabs)
summary(poiss1)
# Run as a GLM to evaluate levels of overdispersion
test1<-glm(n.post60~(pres.num+model+pb.group+pres.num*model*pb.group),data=noldbabs,family=poisson)
summary(test1) # Poisson GLM converges ok, but huge overdispersion. Cannot use poisson
coeftest(test1, vcov = sandwich)
test2<-glm(n.post60~(pres.num+model+pb.group+pres.num*model*pb.group),data=noldbabs,family=quasipoisson)
summary(test2) # Dispersion parameter is 5.523

# Bolker style test for overdispersion (from: http://avesbiodiv.mncn.csic.es/estadistica/curso2011/regm26.pdf)
grpMeans<- with(noldbabs, tapply(n.post60, list(pb.group,model,pres.num),mean))
summary(grpMeans)
grpVars <- with(noldbabs, tapply(n.post60, list(pb.group,model,pres.num),var))
summary(grpVars)
plot(grpVars ~ grpMeans)
abline(c(0, 1), lty = 2) # points fall above the line, so data is overdispersed

### What does SQRT transformation do to n.post60 (will insert transformation code above) so that subsets don't have to be remade
hist(noldbabs$n.post60)
hist(noldbabs$sqrt.post60)
plot(n.post60~sqrt.post60, data=noldbabs)


# I will generate all of the models required for comparison and we can 
# go from there. The main influences of how we will proceed from now are:
# Richards et al http://link.springer.com/article/10.1007%2Fs00265-010-1035-8
# but no model averaging: http://onlinelibrary.wiley.com/doi/10.1890/14-1639.1/full
# using this as an example of application: http://rspb.royalsocietypublishing.org/content/282/1811/20151086

sqrt.sel1<-lmer(sqrt.post60~sqrt.pre60+model*pres.num*pb.group+(1|focal.id)+(1|pb.id),REML=F,control=lmerControl(optimizer="bobyqa"),data=noldbabs)
summary(sqrt.sel1) # So this runs fine, and again shows only two 2-way interactions and no 3 way, just like the others.
plot(sqrt.sel1) # fairly heretogeneous spread, but bounded by zero which makes it ugly, but it is the best avenue available 
AICc(sqrt.sel1)
r.squaredGLMM(sqrt.sel1)

sqrt.sel2<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+model*pres.num+model*pb.group+pb.group*pres.num)+(1|focal.id)+(1|pb.id),REML=F,lmerControl(optimizer="bobyqa"),data=noldbabs)
summary(sqrt.sel2) # Converges ok (optimizer required)
AICc(sqrt.sel2)
r.squaredGLMM(sqrt.sel2)

sqrt.sel3<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+model*pres.num+model*pb.group)+(1|focal.id)+(1|pb.id),REML=F,lmerControl(optimizer="bobyqa"),data=noldbabs)
summary(sqrt.sel3)
AICc(sqrt.sel3)
r.squaredGLMM(sqrt.sel3)

sqrt.sel4<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+model*pb.group+pb.group*pres.num)+(1|focal.id)+(1|pb.id),REML=F,lmerControl(optimizer="bobyqa"),data=noldbabs)
summary(sqrt.sel4) 
AICc(sqrt.sel4)
r.squaredGLMM(sqrt.sel4)

sqrt.sel5<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+model*pres.num+pb.group*pres.num)+(1|focal.id)+(1|pb.id),REML=F,lmerControl(optimizer="bobyqa"),data=noldbabs)
summary(sqrt.sel5)
AICc(sqrt.sel5)
r.squaredGLMM(sqrt.sel5)

sqrt.sel6<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+model*pres.num)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel6)
AICc(sqrt.sel6)
r.squaredGLMM(sqrt.sel6)

sqrt.sel7<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+model*pb.group)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel7)
AICc(sqrt.sel7)
r.squaredGLMM(sqrt.sel7)

sqrt.sel8<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num+pres.num*pb.group)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel8) 
AICc(sqrt.sel8)
r.squaredGLMM(sqrt.sel8)

sqrt.sel9<-lmer(sqrt.post60~(sqrt.pre60+model+pres.num+model*pres.num)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel9)
AICc(sqrt.sel9)
r.squaredGLMM(sqrt.sel9)

sqrt.sel10<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+model*pb.group)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel10)
AICc(sqrt.sel10)
r.squaredGLMM(sqrt.sel10)

sqrt.sel11<-lmer(sqrt.post60~(sqrt.pre60+pb.group+pres.num+pres.num*pb.group)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel11) 
AICc(sqrt.sel11)
r.squaredGLMM(sqrt.sel11)

sqrt.sel12<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group+pres.num)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel12)
AICc(sqrt.sel12)
r.squaredGLMM(sqrt.sel12)

sqrt.sel13<-lmer(sqrt.post60~(sqrt.pre60+model+pb.group)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel13)
AICc(sqrt.sel13)
r.squaredGLMM(sqrt.sel13)

sqrt.sel14<-lmer(sqrt.post60~(sqrt.pre60+model+pres.num)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel14)
AICc(sqrt.sel14)
r.squaredGLMM(sqrt.sel14)

sqrt.sel15<-lmer(sqrt.post60~(sqrt.pre60+pb.group+pres.num)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel15)
AICc(sqrt.sel15)
r.squaredGLMM(sqrt.sel15)

sqrt.sel16<-lmer(sqrt.post60~(sqrt.pre60+model)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel16)
AICc(sqrt.sel16)
r.squaredGLMM(sqrt.sel16)

sqrt.sel17<-lmer(sqrt.post60~(sqrt.pre60+pb.group)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel17) 
AICc(sqrt.sel17)
r.squaredGLMM(sqrt.sel17)

sqrt.sel18<-lmer(sqrt.post60~(sqrt.pre60+pres.num)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel18)
AICc(sqrt.sel18)
r.squaredGLMM(sqrt.sel18)

sqrt.sel19<-lmer(sqrt.post60~(sqrt.pre60)+(1|focal.id)+(1|pb.id),lmerControl(optimizer="bobyqa"),REML=F,data=noldbabs)
summary(sqrt.sel19)
AICc(sqrt.sel19)
r.squaredGLMM(sqrt.sel19)


# Use MuMIn package to get the R-Sq values of our GLMMs, and determine the top set for application of the 
# nesting rule as Richards et al (2011)suggest on p79
model.compare<-model.sel(sqrt.sel1,sqrt.sel2,sqrt.sel3,sqrt.sel4,sqrt.sel5,sqrt.sel6,sqrt.sel7,sqrt.sel8,
                         sqrt.sel9,sqrt.sel10,sqrt.sel11,sqrt.sel12,sqrt.sel13,sqrt.sel14,sqrt.sel15,
                         sqrt.sel16,sqrt.sel17,sqrt.sel18,sqrt.sel19,rank="AICc")
model.compare
write.table(model.compare, "C:\\\\RStudio files\\\\Babar Main Analysis Models Table.csv", sep=",")

# Lastly, can I obtain the Akaike Weights for the two models retained once the nesting rule
# was applied? 
(Weights(AICc(sqrt.sel6, sqrt.sel9)))








### Consider Presentation 3 only
p3.model1<-lm(sqrt.post60~sqrt.pre60+pb.group*model, data=pres3)
summary(p3.model1)
plot(p3.model1) # plots look ok
p3.model2<-lm(sqrt.post60~(sqrt.pre60+pb.group+model), data=pres3)
summary(p3.model2)
p3.model3<-lm(sqrt.post60~(sqrt.pre60+pb.group), data=pres3)
summary(p3.model3)
p3.model4<-lm(sqrt.post60~(sqrt.pre60+model), data=pres3)
summary(p3.model4)
p3.model5<-lm(sqrt.post60~(pb.group*model), data=pres3)
summary(p3.model5)
p3.model6<-lm(sqrt.post60~(pb.group+model), data=pres3)
summary(p3.model6)
p3.model7<-lm(sqrt.post60~(sqrt.pre60), data=pres3)
summary(p3.model7)
p3.model8<-lm(sqrt.post60~(model), data=pres3)
summary(p3.model8)
p3.model9<-lm(sqrt.post60~(pb.group), data=pres3)
summary(p3.model9)

model.compare<-model.sel(p3.model1,p3.model,p3.model3,p3.model4,
                         p3.model5,p3.model6,p3.model7,p3.model8,
                         p3.model9,rank="AICc")
model.compare

write.table(model.compare, "C:\\\\RStudio files\\\\omega.csv", sep=",")



### Graphs ###

# Figure 1a
palette<-c("grey97","grey45")
ggplot(data=babs, aes(x=pres.num,y=n.post60,fill=pb.group))+geom_boxplot()+ylim(0,30)+
  scale_fill_manual(values=palette,name="Playback Group")+theme(legend.position="right")+
  labs(x="Presentation Number",y="Number of flights")+ggtitle("a) Elephant Group")+theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(colour="black"))+
  theme(text = element_text(size=14,face= "bold",colour="black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.text=element_text(size=14,colour="black"))

# Figure 1b
palette<-c("darkorange1","darkorange4")
ggplot(data=foxy, aes(x=pres.num,y=n.post60,fill=pb.group))+geom_boxplot()+ylim(0,30)+
  scale_fill_manual(values=palette,name="Playback Group")+theme(legend.position="right")+
  labs(x="Presentation Number",y="Number of Flights")+ggtitle("b) Fox Group")+theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(colour="black"))+
  theme(text = element_text(size=14,face= "bold",colour="black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.text=element_text(size=14,colour="black"))  

stat.desc(pres3$n.pre60)
stat.desc(pres3$n.post60)

# Figure 2
palette<-c("grey75","grey45","darkorange1","darkorange4")
ggplot(pres3, aes(x=n.pre60,y=n.post60,colour=pb.model))+geom_point(size=3)+ylim(0,30)+xlim(0,22)+
  scale_fill_manual(values=palette)+scale_colour_manual(values=palette,name="Model & \\nPlayback Group")+
  labs(x="Number of flights in 60s prior to Presentation",y="Number of Flights in 60s post Presentation")+
  stat_smooth(method=lm,se=FALSE,fullrange=TRUE,size=1.6)+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(colour="black"))+
  theme(text = element_text(size=14,face= "bold",colour="black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.text=element_text(size=14,colour="black")) 


# SuppMat1
ggplot(babs,aes(x=pres.num,y=n.post60,group=focal.id,colour=pb.group))+geom_line(size=1.5)+
  labs(x="Presentation Number",y="Number of Flights in 60s after Presentation")+ggtitle("Change in individuals' response to the Elephant presentation")+
  scale_colour_manual(values=c("grey70","grey30"))+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(colour="black"))+
  theme(text = element_text(size=14,face= "bold",colour="black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.text=element_text(size=14,colour="black"))  
# SuppMat2
ggplot(foxy,aes(x=pres.num,y=n.post60,group=focal.id,colour=pb.group))+geom_line(size=1.5)+
  labs(x="Presentation Number",y="Number of Flights in 60s after Presentation")+ggtitle("Change in individuals' response to the Fox presentation")+
  scale_colour_manual(values=c("darkorange1","darkorange4"))+
  theme(panel.background = element_rect(fill='white'), axis.line = element_line(colour="black"))+
  theme(text = element_text(size=14,face= "bold",colour="black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.text=element_text(size=14,colour="black"))  
