#Bee foraging behaviour

#Load Packages
#require(sjPlot)
require(ggplot2)
require(lme4)
require(plyr)
require(lmerTest)

###########################
#Table 4 SI & Figure 1e: Proportions of visits to each concentration (colony level)

#set data folder (add path between ""):
data = ""

sum_all_t           <- read.csv(paste0(data, "colony_level_visits_&_feed_time.csv"))
sum_all_t$treatment <- factor(sum_all_t$treatment, levels=(c("0 ppb","2 ppb","11 ppb"))) # relevels the factor to C, L, H 

#fit the mixed effect model and view summary:
x                   <-cbind(sum_all_t$treat_visits,(sum_all_t$day_visits-sum_all_t$treat_visits))
TS4                 <-glmer(x~day*treatment+period+(day|colony), binomial, data=sum_all_t)  
summary(TS4)

#the random effects that we include in the model explain essentialy 0 variance
#however i'm keeping them in due to the design of the experiment (repeated measures on 
#the colonies)

#Produce the model predictions for amount of sucrose consumed.
#create an empty DF to use with "predict"

a<-data.frame(expand.grid(
  treatment=c("0 ppb", "2 ppb", "11 ppb"),
  period=c("P1"),
  day=c(1,2,4)
))

b<-data.frame(expand.grid(
  treatment=c("0 ppb", "2 ppb", "11 ppb"),
  period=c("P2"),
  day=c(6,8,10)
))

b<-rbind(a,b)

b$prop<- predict(TS4, b,type="response",re.form=NA)

ggplot(sum_all_t, aes(x=as.factor(day), y=((treat_visits/day_visits))))+ 
  theme_bw()+
  geom_boxplot(alpha=0)+  facet_wrap( ~ treatment)+
  geom_point(data=b, aes(x=as.factor(day), y=prop),shape=21, fill="red")+
  scale_x_discrete(name="Day") +
  scale_y_continuous(name="Proportion of visits", limits = c(0,0.6))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(size=24),
        legend.position="none"
  )+ 
  geom_hline(yintercept=0.33)+
  theme(aspect.ratio=6/3)

###########################

#Table 5 SI: Counts of workers visiting each concentration

TS5 <-glmer(treat_visits~day*treatment+period+(day|colony),  
            family = "poisson", data=sum_all_t
            )  

summary(TS5)

#Produce the model predictions for amount of sucrose consumed.
#create an empty DF to use with "predict"

a<-data.frame(expand.grid(
  treatment=c("0 ppb", "2 ppb", "11 ppb"),
  period=c("P1"),
  day=c(1,2,4)
))

b<-data.frame(expand.grid(
  treatment=c("0 ppb", "2 ppb", "11 ppb"),
  period=c("P2"),
  day=c(6,8,10)
))

b<-rbind(a,b)

b$prop<- predict(TS5, b,type="response",re.form=NA)




###########################

#Table S6 SI: Time spent at each concentration
#Mixed effects model (note - times were not recorded for Day 1 observations so these have been excluded from model)

TS6 <-lmer(av_feed~day*treatment+period+(day|colony),  data=sum_all_t[ which(sum_all_t$day !=1), ])  #full model, with RE data correctly specified. Is this all necessary?
summary(TS6)

#Produce the model predictions for amount of sucrose consumed.
#create an empty DF to use with "predict"

a<-data.frame(expand.grid(
  treatment=c("0 ppb", "2 ppb", "11 ppb"),
  period=c("P1"),
  day=c(1,2,4)
))

b<-data.frame(expand.grid(
  treatment=c("0 ppb", "2 ppb", "11 ppb"),
  period=c("P2"),
  day=c(6,8,10)
))

b<-rbind(a,b)

b$time<- predict(TS6, b,type="response",re.form=NA)

###########################

#Table S7 SI: Proportion of visits to each feeder (Only individuals, observed on >3 days
#and contribute >1% of observed foraging trips)

commited.forager = read.csv(paste0(data, "commited_foragers.csv"))

y                 = cbind(commited.forager$treat_visits, (commited.forager$day_visits-commited.forager$treat_visits) )

TS7  =glmer(y~day*treatment+period+(day|colony/UID), family = "binomial", data=commited.forager)
summary(TS7)

#plot the data(SI):

new_dat<-expand.grid (      #creates a data frame with the model predictors 
  day       = c(1,2,4,6,8,10),
  treatment = levels(commited.forager$treatment)
  
) 
new_dat$period  = rep( c(  rep("P1",3), rep("P2",3)  ), 3)

#to plot the graphs, first order the df by treatment and day (Fig S4):

new_dat=with(new_dat, new_dat[order(treatment,day ),])
new_dat$resp.pred<-predict(TS7,newdata=new_dat,re.form=~0, type = "response") #Prediction

#now produce the boxplots:

ggplot(commited.forager, aes(x=as.factor(day), y=((treat_visits/day_visits))))+ 
  theme_bw()+
  geom_boxplot(alpha=0)+  facet_wrap( ~ treatment)+
  geom_point(data=new_dat, aes(x=as.factor(day),  y=resp.pred),shape=21, fill="red")+
  scale_x_discrete(name="Day") +
  scale_y_continuous(name="Proportion of visits")+
  expand_limits(y=c(0,1))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #text = element_text(size=24),
        legend.position="none")+ 
  theme(aspect.ratio=6/3)









