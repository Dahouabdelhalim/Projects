#Journal of applied eccology - script submission
#Script contains details of the final models

require(lme4) 
require(plyr) 
require(data.table)
require(Rmisc)
require(plyr)
require(MASS)


#load data
census<- read.csv("Arce.et.al_census.csv")
obs1  <- read.csv("Arce.et.al_observer1.csv")
obs2  <- read.csv("Arce.et.al_observer2.csv")
all   <- read.csv("Arce.et.al_forage_freq.csv")

census$queen<-as.factor(census$queen)
census$pair<-as.factor(census$pair)

# When wind excluded there are no missing variables giving 600 observations.
# removing the NA from the dataset allows model with and without wind to be compared using "anova" command
obs1b<- na.omit(obs1)
obs2b<- na.omit(obs2)
all.2<-na.omit(all) 

#t-test for sucrose consumed, ml/ day between treatment groups
control	  <-c(	639,	497.5,	380.5,	410,	529,	648,	381.5,	811.5,	770,	407.5)
treatment	<-c(	804.5,	718.5,	375,	519.5,	449.5,	445.5,	516.5,	688.5,	565.5,	552.5)
t.test(control,treatment)

#Comparing the original colony size (worker and pupae No.)
obs_W  <-glm(workers.st~obs,       family = poisson, data=census)
obs_P  <-glm(pupae.st  ~obs,       family = quasipoisson, data=census)
treat_W<-glm(workers.st~treatment, family = poisson, data=census)
treat_P<-glm(pupae.st  ~treatment, family = poisson, data=census)
treat_P<-glm(pupae.st  ~treatment, family = quasipoisson, data=census)

summary(obs_W)
summary(obs_P)
summary(treat_W)
summary(treat_P)

########################################################################################################################################################
#Table S5a and S8a
#Number of foragers per hour
S5a <-glmer(Freq~treatment+obs.hour
            +scale(wind)+scale(temperature)
            +treatment:obs.hour
            +(day|observer/pair/colony), #full model
            family=poisson,
            data=all.2,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

S8a<-glmer(Freq~treatment+obs.hour
            +scale(wind)+scale(temperature)
            +treatment:obs.hour+observer
            +(day|pair/colony), #full model
            family=poisson,
            data=all.2,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

source(system.file("utils", "allFit.R", package="lme4"))

S5a.all <- allFit(S5a)
S8a.all <- allFit(S8a)
ss_S5a <- summary(S5a.all)
ss_S8a <- summary(S8a.all)
ss_S5a
ss_S8a

summary(S5a)
summary(S8a)

########################################################################################################################################################
#Table S5b and S8b
#Proportion of pollen carriers per hour

y<-cbind(all.2$success,all.2$fail)# gives the response variable for foragers returning with pollen

S5b <-glmer(y~treatment+obs.hour+wind+temperature+treatment:obs.hour+(day|observer/pair/colony), #full model
           family=binomial(link='logit'),
           data=all.2,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

S8b<-glmer(y~treatment+obs.hour+wind+temperature+treatment:obs.hour+observer+(day|pair/colony), #full model
            family=binomial(link='logit'),
            data=all.2,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

S5b.all <-allFit(S5b)
S8b.all <-allFit(S8b)
ss_S5b  <-summary(S5b.all)
ss_S8b  <-summary(S8b.all)

summary(S5b)
summary(S8b)

#standard model fails to converge, however all available optimisers are practically equivalent so the error message is likely a 
#false positive see: help("convergence")

########################################################################################################################################################
#c) Table S5c:
#Average weight of pollen per hour

S5c<-lmer(p.ball~treatment+obs.hour+temperature+treatment:obs.hour+(day|pair/colony), #full model
          data=obs1, subset=(pollen=="1"))

summary(S5c)

########################################################################################################################################################
#d) Table S5d:
#Average area of pollen per hour

S5d<-lmer(p.ball~treatment+obs.hour+treatment:obs.hour+(day|pair/colony), #full model
           data=obs2,subset=(pollen=="1"))

S5d.all <- allFit(S5d)
ss_S5d <- summary(S5d.all)
ss_S5d
summary(S5d)

#standard model fails to converge, however all available optimisers are practically equivalent so the error message is likely a 
#false positive see: help("convergence")

########################################################################################################################################################
#Table S5e:
#total weight of pollen per hour
dt1 <- data.table(obs1)
sum.em<-dt1[,list(sum=sum(p.ball, na.rm=TRUE), temperature=mean(temperature), obs.hour=mean(obs.hour), wind=mean(wind)),by=c("pair","colony","day","treatment")]

S5e <-lmer(sum~treatment+obs.hour+treatment:obs.hour+(day|pair/colony),    
           data=sum.em)

summary(S5e)

########################################################################################################################################################
#Table S5f:
#total area of pollen per hour
dt2 <- data.table(obs2b)
sum.tom<-dt2[,list(sum=sum(p.ball, na.rm=TRUE), temperature=mean(temperature), obs.hour=mean(obs.hour), wind=mean(wind)),by=c("pair","colony","day","treatment")]

S5f<-lmer(sum~treatment+obs.hour+treatment:obs.hour+(day|pair/colony), 
          data=sum.tom)

summary(S5f)

########################################################################################################################################################
#Table S6a and Table 8c
#number of foragers per day
day2<-all$day^2

S6a<-glmer(Freq~treatment+day+day2+
           temperature+treatment:day+treatment:day2+
           treatment+
           (day|observer/pair/colony), #wind removed
           family=poisson,
           data=all,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))) # use the whole dataset if not including wind

S8c<-glmer(Freq~treatment+day+day2+
            temperature+treatment:day+treatment:day2+
            treatment+observer+
            (day|pair/colony), #wind removed
            family=poisson,
            data=all,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))) # use the whole dataset if not including wind

summary(S6a)
summary(S8c)


ss.CS <- transform(all, day2=scale(day2)) #<- this scales the variable to remove the warning messages
S6a <- update(S6a, data=ss.CS)
S8c <- update(S8c, data=ss.CS)

S6a.all <- allFit(S6a)
S8c.all <- allFit(S8c)

ss_S6a <- summary(S6a.all)
ss_S8c <- summary(S8c.all)

summary(S6a)
summary(S8c)

#standard model fails to converge, however fixef from all available optimisers are practically equivalent so the error message is likely a 
#false positive see: help("convergence")

########################################################################################################################################################
#Table S6b and Table S8d
#proportion of pollen carriers per day

all$day2<-all$day^2
y<-cbind(all$success,all$fail)# gives the response variable for foragers returning with pollen

S6b<-glmer(y~treatment+day+day2+wind+temperature+treatment:day+treatment:day2+(day|observer/pair/colony), #full model
           family=binomial(link='logit'),
           data=all,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
S8d<-glmer(y~treatment+day+day2+wind+temperature+treatment:day+treatment:day2+observer+(day|pair/colony), #full model
            family=binomial(link='logit'),
            data=all,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))


ss.CS <- transform(all, day2=scale(day2))

S6b  <- update(S6b, data=ss.CS)
S8d <- update(S8d, data=ss.CS)

S6b.all  <- allFit(S6b)
S8d.all  <- allFit(S8d)
ss_S6b <- summary(S6b.all)
ss_S8d <- summary(S8d.all)
ss_S6b
ss_S8d

summary(S6b)
summary(S8d)

#standard model fails to converge, however fixef from all available optimisers are practically equivalent so the error message is likely a 
#false positive see: help("convergence")


########################################################################################################################################################
#Table S6c
#Average weight of pollen per day - LMER, normal distribution

obs1$day2<-obs1$day^2
             
S6c<-lmer(p.ball~treatment+day+day2+treatment:day+treatment:day2+(day|pair/colony), #full model <-dont include temp
            data=obs1, subset=(pollen=="1"),control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

#model converges when temperature removed. 

summary(S6c)

#model fails to converge, with temperature, fixef from all available optimisers are not similar error message is likely a real issue with 
#convergence e.g. overspecification of model, removing temperature from the model. 
#false positive see: help("convergence")
#######################################################################################################################################################
#Table S6d
#d) average area of pollen ball area per day
obs2$day2<-as.integer(obs2$day^2)

S6d<-lmer(p.ball~treatment+day+temperature+treatment:day+(day|pair/colony), #full model
                        data=obs2, subset=(pollen=="1"))

summary(S6d)

########################################################################################################################################################
#Table S6e
#e) sum weight of pollen collected per day
dt1 <- data.table(obs1b)
sum.em <-dt1[,list(sum=sum(p.ball, na.rm=TRUE), temperature=mean(temperature), wind=mean(wind), 
        hour=mean(obs.hour)),by=c("pair","colony","day","treatment")]
sum.em$day2<-sum.em$day^2


S6e<-lmer(sum~treatment+day+day2+wind+treatment:day+treatment:day2+(day|pair/colony), #full model
           data=sum.em)

summary(S6e)

########################################################################################################################################################
#Table S6f
#Sum area of pollen collected per day
obs2$day2<-obs2$day^2
dt2 <- data.table(obs2)
sum.tom<-dt2[,list(sum=sum(p.ball, na.rm=TRUE), temperature=mean(temperature), wind=mean(wind), hour=mean(obs.hour)),by=c("pair","colony","day", "day2","treatment")]
             
S6f<-lmer(sum~treatment+day+day2+wind+treatment:day+treatment:day2+(day|pair/colony), #full model
          data=sum.tom, na.action = na.omit)
           
########################################################################################################################################################
#Table S7a

#S7 a)Change in colony weight
census$pair<-as.factor(census$pair)
S7a <-lmer(weight.increase~treatment+b.wt.st+(1|obs/pair), data=census) #full
summary(S7a)

#S7 b)Number of eggs
S7b <-glmer(eggs~treatment+b.wt.st+(1|obs)+(1|pair),family=poisson, data=census) #full
summary(S7b)
             
#S7 c)Number of larvae
S7c <-glmer(sum.larval.no~treatment+b.wt.st+(1|obs)+(1|pair), family=poisson, data=census)
summary(S7c)
             
#S7 d)Number of pupae
S7d <-glmer(pupae~+treatment+b.wt.st+(1|obs)+(1|pair), family=poisson, data=census)
summary(S7d)
             
#S7 e)Number of workers
S7e <-glmer(workers~+treatment+b.wt.st+(1|obs)+(1|pair), family=poisson, data=census)
summary(S7e)
             
#S7 f)Number of drones
S7f<-glmer(males~+treatment+b.wt.st+(1|obs)+(1|pair), family=poisson, data=census)
summary(S7f)           

#S7 g)Number of gynes
S7g<-glmer(gyne~+treatment+b.wt.st+(1|obs)+(1|pair), family=poisson, data=census)
summary(S7g)