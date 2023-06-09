#script for replicating statistical analyses

fem<-read.csv("female.choice.data.for.archiving.csv",stringsAsFactors = FALSE,na.strings="NA")

#Chi-squared test of outcome (mating or not) and cross type (interspecific or intraspecific)
c<-table(fem$cross.type,fem$outcome)
chisq.test(c)

#analysis of male and female behaviours
fem$courts<-fem$male.court.or.chase+fem$male.mate.attempt #creates a column of the total number of male courtship behaviours conducted in each trial

wilcox.test(fem$courts~fem$male.species,paired=FALSE)#Mann-Whitney U test to determine whether total number of courtships of H. pachinus and H. cydno males are drawn from the same distribution

#logistic regressions to determine effect of male species on courtship rate
crts<-fem[complete.cases(fem[ ,21]),] #creates a subset of data frame containing only observations that have no missing data for total number of courtships
z<-glm(outcome~male.species*courts,data=crts,family=binomial(link="logit")) #logistic regression with both male species and number of courtships as predictors
anova(z,test="Chisq") #likelihood ratio tests of coefficients in the full model
z1<-glm(outcome~male.species,data=crts,family=binomial(link="logit")) #reduced model of logistic regression with only male species as predictor
anova(z,z1,test="Chisq") #likelihood ratio test comparing full to reduced model
AIC(z) #Akaike's Information Criterion of full model
AIC(z1) #Akaike's Information Criterion of reduced model

#create columns for numbers of each female behaviour per male behaviour (i.e. divided by total number of male courtship behaviours) to correct for collinearity of number of male behaviours with number of female responses
fem$flutter<-fem$female.flutter/fem$courts
fem$open<-fem$female.wings.open/fem$courts
fem$closed<-fem$female.wings.closed/fem$courts
fem$fly<-fem$female.fly/fem$courts

#do female behaviours predict trial outcome?
intra<-subset(fem,fem$male.species=="cydno") #creates subset of data frame containing only intraspecific trials
#logistic regressions of trial outcome versus each female behaviour rate
z1<-glm(outcome~closed,data=intra,family=binomial(link="logit"))
z2<-glm(outcome~open,data=intra,family=binomial(link="logit"))
z3<-glm(outcome~flutter,data=intra,family=binomial(link="logit"))
z4<-glm(outcome~fly,data=intra,family=binomial(link="logit"))
#output of these logistic regression are accessed using summary(z)

#Mann-Whitney U tests to determine whether female have different behaviour rates towards conspecific and heterospecific males
wilcox.test(fem$closed~fem$male.species,paired=FALSE)
wilcox.test(fem$open~fem$male.species,paired=FALSE)
wilcox.test(fem$fly~fem$male.species,paired=FALSE)
wilcox.test(fem$flutter~fem$male.species,paired=FALSE)
