############## Does the solving behaviour spread socially throughout the populations and do rates differ between social network structures? ########


library(NBDA) #downloaded from WH github

#Create NBDA object for all pens combined

load("cumulativeLearningPotentialNetworksChicksModularityCombined.RData")
#Here I made the 15 'big' matrices showing which birds went into the testing together each round for all 6 pens together with zero's between the different pens
#then I used the learning indicator array to get learning potential per round and then cumulative learning potential
#Id's are recoded so that pens: 4A=1-27, 4B=28-54, 5A=55-78, 5B=79-102, 7A=103-126, 7B=127-150)

dim(cumulativeLearningPotential)

#We do want 6 networks- with the connections from each group, and with zeroes for the other group, allowing us to test for differences in social
#transmission as before. This is easily done:

cumulativeLearningPotential2<-array(0,dim=c(150,150,6,15))

#WH remember in a stratified OADA the network for group 1 includes the individuals from other groups, but just has zero connections for any dyads where both individuals are not in group 1
#and the same for group 2 etc.

cumulativeLearningPotential2[1:27,1:27,1,]<-cumulativeLearningPotential[1:27,1:27,]
cumulativeLearningPotential2[28:54,28:54,2,]<-cumulativeLearningPotential[28:54,28:54,]
cumulativeLearningPotential2[55:78,55:78,3,]<-cumulativeLearningPotential[55:78,55:78,]
cumulativeLearningPotential2[79:102,79:102,4,]<-cumulativeLearningPotential[79:102,79:102,]
cumulativeLearningPotential2[103:126,103:126,6,]<-cumulativeLearningPotential[103:126,103:126,]
cumulativeLearningPotential2[127:150,127:150,6,]<-cumulativeLearningPotential[127:150,127:150,]

#here is the aquisition order for all new coded birds:
Aq_order_all<-c(65,84,61,14,95,81,35,12,21,29,34,47,39,62,77,85,96,98,106,7,8,30,51,43,44,49,31,68,79,86,100,115,1,10,2,15,24,28,38,72,75,78,88,83,91,99,133,4,18,26,32,54,53,67,93,97,127,131,11,13,19,41,45,56,71,3,22,69,76,87,102,17,27,36,57,89,116,25,40,55,70,90)

#here is the round those birds learnt in:
assMatrixIndexAll<-c(2,2,5,6,6,6,7,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,13,13,13,13,13,13,14,14,14,14,14,14,15,15,15,15,15)

#and the true ties list:
trueTiesAll<-list(c(8,9),c(11,12),c(22,23),c(33,34),c(36,37),c(44,45,46),c(51,52))

demons = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Individual level variables:
#WH in the NBDA package we specify each ILV as a matrix- this is to allow for the possibility we have time-dependent ILVs with multiple columns, one for each event
sex4A<-cbind(c(0,0,0,0,0,1,0,1,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0)) 
mass24A<-cbind(c(133,106,148,121,143,121,121,103,114,126,119,120,122,143,95,104,108,120,122,106,120,88,94,102,111,136,115))

sex4B<-cbind(c(0,0,0,0,1,1,1,0,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1,1,1,0,1))
mass24B<-cbind(c(124,94,121,120,90,107,105,121,121,104,123,120,102,103,101,79,102,96,143,107,92,94,84,104,107,128,114))

sex5A<-cbind(c(0,0,1,0,0,0,1,1,1,0,0,0,0,1,0,1,1,1,0,0,1,1,1,1))
mass25A<-cbind(c(124,124,116,120,130,137,111,111,134,133,124,156,123,122,129,131,134,142,142,133,136,128,110,139))

sex5B<-cbind(c(0,0,1,1,0,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,1,1,0,1))
mass25B<-cbind(c(118,132,137,124,120,133,138,114,142,132,125,122,132,136,142,135,129,114,133,144,115,122,117,143))

sex7A<-cbind(c(0,0,0,1,1,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,1,1,0,1))
mass27A<-cbind(c(105,110,134,128,108,116,121,129,120,118,133,128,130,119,119,116,115,112,112,128,142,122,119,117))

sex7B<-cbind(c(0,1,1,0,1,0,1,0,0,0,0,1,1,1,1,1,0,0,0,0,1,0,1,1))
mass27B<-cbind(c(127,137,146,120,121,121,125,118,114,127,116,131,128,148,125,129,98,132,116,118,120,127,124,119))

sex<-rbind(sex4A,sex4B,sex5A,sex5B,sex7A,sex7B)
mass2<-rbind(mass24A, mass24B,mass25A, mass25B, mass27A, mass27B)

#Usually in a multi-diffusion using TADA or stratified OADA we would want to have a set of binary ILVs (amounting to a factor) controlling for the possibility 
#that each group might have a different rate of asocial learning- which might look like social transmission.
#However, here is a little different. Individuals were assigned to groups at random, and there is no reason to assume that different groups would
#have different rates of asocial learning.

asoc_ilv<-c("sex", "mass2") 
int_ilv<-c("sex","mass2")


AllPens<-nbdaData(label="AllPens", assMatrix=cumulativeLearningPotential2,assMatrixIndex=assMatrixIndexAll,orderAcq=Aq_order_all, demons=demons,asoc_ilv=asoc_ilv,int_ilv=int_ilv,trueTies = trueTiesAll)

model_t<-oadaFit(AllPens)

#We are using the unconstrained model to allow for the possibility that sex/mass might have an effect on social learning that is different to its effect on asocial learning.

#In order to run the AIC tables we specify what I call a constraintsVectMatrix.
#Each line of the matrix specifies a model.
#Each entry in that line represents a parameter in the model
#There are 6 networks in the object, so the first 6 parameters are the s parameters corresponding to
#each network.
#The next 1 corresponds to the effect of sex on asocial learning
#The next 1 corresponds to the effect of mass on asocial learning
#The next 1 corresponds to the effect of sex on social learning
#The next 1 corresponds to the effect of mass on social learning
#When a 0 appears- it means that parameter is constrained to have a value of zero- i.e. the corresponding variable is dropped from the model
#When two parameters have the same non-zero number assigned to them they are constrained to have the same value.


# For no effect of the ILV we set both to 0.
#e.g. To fit a model in which sex affects only asocial learning we set the 7th entry to 1 and the next 3 to 0. 

#setting S parameter to look at whether the social transmission rate is equal for all 4 pens, differs in alll 4 pens or is dependent on grouping condition (Pen A and D are the same and C and B are the same)

constraintsVectMatrix<-rbind(
  #s1,s2,s3,s4,s5,s6,7sex-asocial,8mass-asocial, 9sex-social, 10mass-social
  #Asocial model
  #No ILVs
  c(0,0,0,0,0,0,0,0,0,0),
  #Asocial effect of sex
  c(0,0,0,0,0,0,1,0,0,0),
  #Asocial effect of mass
  c(0,0,0,0,0,0,0,1,0,0),
  #ILVs cannot have an effect when there is no social learning so we do not
  #need the two models with an effect of sex/mass on social learning
  
  #Same s Params each diffusion, none for the control group
  #No ILVs
  c(1,1,1,1,0,0,0,0,0,0),
  #Asocial effect of sex
  c(1,1,1,1,0,0,2,0,0,0),
  #Asocial effect of mass
  c(1,1,1,1,0,0,0,2,0,0),
  #Asocial effect of sex and mass
  c(1,1,1,1,0,0,2,3,0,0),
  #Social effect of sex
  c(1,1,1,1,0,0,0,0,2,0),
  #Social effect of mass
  c(1,1,1,1,0,0,0,0,0,2),
  #Social effect of sex and mass
  c(1,1,1,1,0,0,0,0,2,3),
  #Asocial and social effect of sex
  c(1,1,1,1,0,0,2,0,3,0),
  #Asocial and social effects of mass
  c(1,1,1,1,0,0,0,2,0,3),
  #Asocial and social effects of sex and mass
  c(1,1,1,1,0,0,2,3,4,5),
  
  #One s for high mod, another for low mod
  #No ILVs
  c(1,2,2,1,0,0,0,0,0,0),
  #Asocial effect of sex
  c(1,2,2,1,0,0,3,0,0,0),
  #Asocial effect of mass
  c(1,2,2,1,0,0,0,3,0,0),
  #Asocial effect of sex and mass
  c(1,2,2,1,0,0,3,4,0,0),
  #Social effect of sex
  c(1,2,2,1,0,0,0,0,3,0),
  #Social effect of mass
  c(1,2,2,1,0,0,0,0,0,3),
  #Social effect of sex and mass
  c(1,2,2,1,0,0,0,0,3,4),
  #Asocial and social effect of sex
  c(1,2,2,1,0,0,3,0,4,0),
  #Asocial and social effects of mass
  c(1,2,2,1,0,0,0,3,0,4),
  #Asocial and social effects of sex and mass
  c(1,2,2,1,0,0,3,4,5,6),
  
  
  #Different s Params each diffusion
  #No ILVs
  c(1,2,3,4,0,0,0,0,0,0),
  #Asocial effect of sex
  c(1,2,3,4,0,0,5,0,0,0),
  #Asocial effect of mass
  c(1,2,3,4,0,0,0,5,0,0),
  #Asocial effect of sex and mass
  c(1,2,3,4,0,0,5,6,0,0),
  #Social effect of sex
  c(1,2,3,4,0,0,0,0,5,0),
  #Social effect of mass
  c(1,2,3,4,0,0,0,0,0,5),
  #Social effect of sex and mass
  c(1,2,3,4,0,0,0,0,5,6),
  #Asocial and social effect of sex
  c(1,2,3,4,0,0,5,0,6,0),
  #Asocial and social effect of mass
  c(1,2,3,4,0,0,0,5,0,6),
  #Asocial and social effect of sex and mass
  c(1,2,3,4,0,0,5,6,7,8)
)


#This fits all the models defined above and orders them by AICc
aicTable1<-oadaAICtable(nbdadata=AllPens, constraintsVectMatrix= constraintsVectMatrix)

print(aicTable1)
save(aicTable1,file="StratOADAaicTable1.rdata")

#We can get the support for different social learning hypotheses here:
networksSupport(aicTable1)

#            support       numberOfModels 
#0:0:0:0:0:0 0.00003283236              3
#1:1:1:1:0:0 0.51365778533             10
#1:2:2:1:0:0 0.16697439070             10
#1:2:3:4:0:0 0.31933499160             10

#We can see what the parameter estimates are by requesting model averaged estimates, 
#but filtering by network combination by giving the argument netFilter

modelAverageEstimates(aicTable1, includeAsocial = F, netFilter = "1:1:1:1:0:0")
#s1          s2          s3          s4          s5          s6                         
#0.640486312  0.640486312  0.640486312  0.640486312  0.000000000  0.000000000 -0.066188529 -0.001900636  0.028949326 -0.001857035 

#out of interest
modelAverageEstimates(aicTable1, includeAsocial = F, netFilter = "1:2:2:1:0:0")
#s1       s2         s3       s4        s5     s6       ASOCIALsex  ASOCIALmass2 SOCIALsex  SOCIALmass2   #For the OADA based on total connections to informed individuals:
#0.818    0.67      0.67      0.818    0.00   0.00      -0.051       -0.002     -0.0019     -0.00215 


#We can also get the support for each variable in the model
variableSupport(aicTable1)
#The support for the s parameters does not mean a lot because there are so many more models with the
#s parameters in than there are without them, so here we are looking at the support for the ILVs
#we can see there is little support for an effect of sex or mass on asocial or social learning rate
#The support is <50% meaning these effects are more likely not to be in the best predictive modelthan they are to be in it

#           s1        s2        s3        s4     s5  s6   ASOC:sex ASOC:mass2 SOCIAL:sex SOCIAL:mass2
#support 0.9999672 0.9999672 0.9999672 0.9999672  0   0   0.2251027  0.2065243  0.3224017    0.2935539    


# The best model

print(aicTable1)[1,]   

#model   type    netCombo baseline CONS.s1 CONS.s2 CONS.s3 CONS.s4 CONS.s5 CONS.s6 CONS.ASOC.sex CONS.ASOC.mass2 CONS.SOCIAL.sex CONS.SOCIAL.mass2 OFFs1 OFFs2
#4     noILVs 1:1:1:1:0:0       NA       1       1       1       1       0       0             0               0               0                 0     0     0
#OFFs3 OFFs4 OFFs5 OFFs6 OFFASOC.sex OFFASOC.mass2 OFFSOCIAL.sex OFFSOCIAL.mass2 convergence  loglik        s1        s2        s3        s4  s5 s6 ASOCIALsex
# 0     0     0     0           0             0             0               0        TRUE    363.873  0.5281524 0.5281524 0.5281524 0.5281524  0  0          0
#ASOCIALmass2 SOCIALsex SOCIALmass2      SEs1      SEs2      SEs3      SEs4 SEs5 SEs6 SEasocialsex SEasocialmass2 SEsocialsex SEsocialmass2     aic    aicc
#     0         0           0       0.2262787 0.2262787 0.2262787 0.2262787    0    0            0              0           0             0 729.746 729.796
#deltaAICc RelSupport AkaikeWeight
#     0          1    0.1503701


#To get 95% CIs for s parameters using a profile likelihood method: We do this by creating new nbdadata objects with a constraintsVect (works as the constraintsVectMatrix did above)
bestModelData1<-constrainedNBDAdata(AllPens,constraintsVect=c(1,1,1,1,0,0,0,0,0,0))
bestModel<-oadaFit(bestModelData1)
#Check the AICc is the same as the model at the top of aicTable1
min(aicTable1@aicc)
bestModel@aicc

#We can then look at the estimates from the model as follows
data.frame(Parameter=bestModel@varNames,Estimate=bestModel@outputPar,SE=bestModel@se)

#Parameter              Estimate      SE             
#Social transmission 1 0.5281524 0.2262787   


#We can get CIs using the profile likelihood technique:

#The usual function plotProfLik does not work for trueTies so I have an experimental alternative function for getting CIs for trueTies:
plotProfLikTrueTies(which=1, model=bestModel,range=c(0,3)) 
#All the values of s which have a profile likelihood below the dotted line would not be rejected at the 5% significance level in
#a likelihood ratio test (LRT) and so are in the 95% CI
profLikCITrueTies(which=1, model=bestModel, interval=c(0,1))     #lower limit 0.202    
profLikCITrueTies(which=1, model=bestModel, interval=c(1,3))     #upper limit 1.158    

                                                       
#We can get the estimated proportion of events that occured by social transmission for our estimate of s and for the end points of the 95% CI
nbdaPropSolveByST(model=bestModel)
#ties strat by all:
#P(Network 1)  P(S offset) 
#0.3984       0.0000        #For the OADA based on total connections to informed individuals:
# s=0 corresponds to 0% events by social transmission

nbdaPropSolveByST(par=0.2020056,nbdadata =bestModelData1)
#0.259
nbdaPropSolveByST(par=1.157808,nbdadata =bestModelData1)
#0.50209


################# Does social structure influence the adoption of different solving techniques by the birds? ############


#Reading in the number of red and blue doors birds opened during the experiment
RBdata<-read.delim(file='ChicksModularity Number RB doors.txt', header = TRUE , sep = "\\t")

####  COMBINED HIGH MODULARITY NETWORKS (excluding those who did not solve so this does not affect assortment bias)
data<-read.delim(file='HighCombinedNetsNAs.txt', header=TRUE, row.names=1, sep="\\t")
AM_combined_high_NAs<-as.matrix(data, nrow=52, ncol=52, dimnames=c(1,1))

library(assortnet)

diag(AM_combined_high_NAs) <- 0
prior<-c(1,1) #could be 0.5,0.5

#subsetting RBdata to HIGH MOD
RB_High<-RBdata[which(RBdata$Treatment=='High'),]

shape1<-RB_High$R + prior[1]
shape2<-RB_High$B + prior[2]

true.pref<-matrix(nrow=length(RB_High$R), ncol=1000)
for(i in 1:1000){
  true.pref[,i]<-rbeta(n=length(RB_High$R), shape1=shape1, shape2=shape2)
}

assort.sim <-NA
for (i in 1:1000){
  assort.sim[i]<-assortment.continuous(AM_combined_high_NAs, true.pref[,i], weighted=TRUE)$r
}
hist(assort.sim)

mean_r <- mean(assort.sim)                  #-0.10
ci_r <- quantile(assort.sim,c(0.025,0.975)) #-0.25 - 0.07


#### COMBINED LOW MODULAIRTY NETWORKS (excluding those who did not solve)

data<-read.delim(file='LowCombinedNetsNAs.txt', header=TRUE, row.names=1, sep="\\t")
AM_combined_low_NAs<-as.matrix(data, nrow=52, ncol=52, dimnames=c(1,1))

library(assortnet)

diag(AM_combined_low_NAs) <- 0
prior<-c(1,1) #could be 0.5,0.5

#subsetting RBData to LOW MOD
RB_Low<-RBdata[which(RBdata$Treatment=='Low'),]

shape1<-RB_Low$R + prior[1]
shape2<-RB_Low$B + prior[2]

true.pref<-matrix(nrow=length(RB_Low$R), ncol=1000)
for(i in 1:1000){
  true.pref[,i]<-rbeta(n=length(RB_Low$R), shape1=shape1, shape2=shape2)
}

assort.sim <-NA
for (i in 1:1000){
  assort.sim[i]<-assortment.continuous(AM_combined_low_NAs, true.pref[,i], weighted=TRUE)$r
}
hist(assort.sim)

mean_r <- mean(assort.sim)                  #-0.003
ci_r <- quantile(assort.sim,c(0.025,0.975)) #-0.14 - 0.19



########### Does social learning and network modularity affect the likelihood of acquiring the solving behaviour? #######

library(survival)

dir()
#data <- read.delim(file="Chicks 4&5&7 survival.txt", header=TRUE, sep="\\t")
data <- read.delim(file="Chicks 4&5&7 survival D rem.txt", header=TRUE, sep="\\t")

head(data)

data$SurvObj<-with(data, Surv(Rounds.to.learn,Censored))

#survival curve split into scrounging treatments
plot(survfit(Surv(Rounds.to.learn,Censored)~Social,data=data),col=c('red','blue'))    #social condition pens (high and low mod = blue) vs asocial ( = red)
plot(survfit(Surv(Rounds.to.learn,Censored)~Modularity,data=data),col=c('red','blue', 'green'))


# Survival of birds in high modularity social, low modularity social and asocial pens
model2<-coxph(Surv(Rounds.to.learn,Censored)~
                Modularity,
              data=data)
summary(model2)

library(coxme) #with random effect of population(pen)
model3<-coxme(Surv(Rounds.to.learn,Censored)~
                Modularity + (1|Pen),
              data=data)
summary(model3)

#HIGH MOD
#95% Wald CI on log scale:
log(12.94)+c(-1,1)*1.96*0.44
# 1.697923 3.422723
#Back transformed:
exp(log(12.94)+c(-1,1)*1.96*.44)
#5.462591 30.652778

#LOW MOD
#95% Wald CI on log scale:
log(10.18)+c(-1,1)*1.96*0.44
# 1.458025 3.182825
#Back transformed:
exp(log(10.18)+c(-1,1)*1.96*.44)
# 4.297464 24.114782

#ASOCIAL
#95% Wald CI on log scale:
log(0.79)+c(-1,1)*1.96*0.29
# -0.8041223  0.3326777
#Back transformed:
exp(log(0.79)+c(-1,1)*1.96*.29)
# 0.4474805 1.3946977

anova(model3)

#Here we are making a new dataframe that models an individuals 'survival' when in social high, social low or asocial
new.data=data.frame(Modularity=c('High','Low','Asocial'))

new.data
a<-survfit(model2, newdata = new.data)
b<-summary(a)
b$surv 
b$std.err
b$upper
b$lower

ModularityConditionn_Survival<-data.frame(b$surv,b$std.err,b$upper,b$lower) #use this dataframe to then plot model predictions 


# Graph

dir()
data2<-read.delim(file= "Chicks4&5&7SurvivalModelOutputs D Rem.txt", header=TRUE, sep="\\t")
as.factor(data2$Modularity)

head(data2)

library(ggplot2)
df <- data.frame(data2$Round,data2$X1.Surv)

ggplot(data=data2, aes(x=Round, y=X1.Surv, group=Modularity))

Graph<-ggplot(data=data2,aes(x=Round, y=X1.Surv, group=Modularity, color=Modularity)) + 
  labs(x = "Round", y = "Cumulative probability of acquiring behaviour") +
  coord_cartesian(ylim=c(0, 1), xlim=c(0,13)) + 
  scale_y_continuous(breaks=seq(0, 1, 0.25)) +
  scale_x_continuous(breaks=seq(0,12,2)) +
  scale_color_manual(values = c("gold1", "cyan4","cyan2")) +
  scale_fill_manual(values = c("gold1", "cyan4","cyan2")) +
  theme(axis.line = element_line(colour="grey"), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent'), panel.border = element_rect(colour='grey', fill=NA)) +
  geom_line(colour="black") + geom_point(aes(shape=Modularity), colour="black", size=2) +
  geom_ribbon(data=data2,aes(ymin=X1.Lower,ymax=X1.Upper, fill=Modularity), alpha=0.5)

Graph
#ggsave(Graph,filename='Chicks457SurvivalGraph.png',bg='transparent')
