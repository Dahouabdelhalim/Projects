


library(MASS)
library(nlme)
library(lme4) 
library(car)
library(RVAideMemoire)

read.table(file.choose(), header=TRUE, sep=",", dec=",") ##### import the files .csv "Costantini_TOTeggs" and "dataset_Costantini"



#####################################################################
##### FEMALE INVESTMENT as number of COCOONS PER INDIVIDUAL
#####################################################################


str(Costantini_TOTeggs)


a1=lm(log10(N_cocoons_per_individual+1)~
        Level.of.sexual.maturity+
        factor(Density)+
        factor(Exp_replicate)+
        mating_opportunities+
        Level.of.sexual.maturity*mating_opportunities+
        factor(Density)*mating_opportunities, 
      subset(Costantini_TOTeggs, PrimaryLast2==1)) #### PrimaryLast2==1 is there because the mean egg production per individual was calculated per pair and is repeated for every line of the pair-->"primaryLast=1" works so that each bowl contributes only one value to the model. 
####Model reduction follows rules in Burnham, K.P., Anderson, D.R. & Huyvaert, K.P. AIC model selection and multimodel inference in behavioral ecology: some background, observations, and comparisons. Behav Ecol Sociobiol 65, 23???35 (2011). https://doi.org/10.1007/s00265-010-1029-6

plot(a1)

vif(a1) #test for multicollinearity --- to calculate on a1 after removing interactions

Anova (a1, type="3")
summary(a1)
AIC(a1)
emmeans(a1, pairwise~  mating_opportunities:Level.of.sexual.maturity)



#####################################################################
##### FEMALE INVESTMENT as number of EGGS PER COCOON
#####################################################################
Costantini_TOTeggs$line=c(1:length(Costantini_TOTeggs[,1]))


b1=glmer(Number_of_eggs~
           Level.of.sexual.maturity+
           factor(Density)+
           factor(Exp_replicate)+
           mating_opportunities+
           Level.of.sexual.maturity*mating_opportunities+
          factor(Density)*mating_opportunities+
           (1|newBowlID)+(1|line),
         family=poisson(link="log"),Costantini_TOTeggs)

overdisp.glmer(b1)

summary(b1)
Anova(b1, type="3")
print(VarCorr(b1),comp="Variance")


#####################################################################
##### OXYDATIVE STATUS MARKERS
#####################################################################

dataset_Costantini$batch<-as.factor(dataset_Costantini$batch)
dataset_Costantini$matingopp<-as.factor(dataset_Costantini$matingopp)
dataset_Costantini$volbowl<-as.factor(dataset_Costantini$volbowl)
dataset_Costantini$agel<-as.factor(dataset_Costantini$age)


##########OXY    OXY      OXY     OXY      OXY      OXY

c1=lm(log10(oxy)~matingopp + volbowl + age + batch + matingopp:volbowl + matingopp:age, data=dataset_Costantini)

Anova (c1, type="3")
summary(c1)
AIC(c1)
emmeans(c1, pairwise~  matingopp)
#########Thiols    Thiols      Thiols       Thiols      Thiols

d1=lm(log10(thiols)~matingopp + volbowl + age + batch + matingopp:volbowl+ matingopp:age, data=dataset_Costantini)

Anova (d1, type="3")
summary(d1)
AIC(d1)
emmeans(d1, pairwise~  matingopp)

#########SOD    SOD      SOD       SOD      SOD     SOD

e1=lm(log10(sod)~matingopp + volbowl + age + batch + matingopp:volbowl+ matingopp:age,
      data=dataset_Costantini)

Anova (e1, type="3")
summary(e1)
AIC(e1)
emmeans(e1, pairwise~  matingopp)

#########GPx    GPx      GPx       GPx      GPx      GPx

f1=lm(log10(gpx)~matingopp + volbowl + age + batch + matingopp:volbowl+ matingopp:age,
      data=dataset_Costantini)

Anova (f1, type="3")
summary(f1)
AIC(f1)
emmeans(e1, pairwise~  matingopp)



