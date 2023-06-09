#Analysing the volume of sucrose consumed by colonies across 10 days:

#load packages:
require(ggplot2)
require(lme4)
require(lmerTest)
require(plyr)

#Analysing the volume of sucrose consumed:

#set data folder (add path between ""):
data = ""

volume<-read.csv(paste0(data, "volume_data.csv"))
volume$treatment<- factor(volume$treatment, levels = c("0 ppb", "2 ppb", "11 ppb"))


# Mixed model for Table S3, average volume of sucrose consumed from each pesticide concentration
TS3<-lmer(tr_vol ~day*treatment+period+(day|colony), data=volume)
summary(TS3)

#Produce the model predictions for amount of sucrose consumed.
#create an empty DF to use with "predict"

a<-data.frame(expand.grid(
  treatment=c("0 ppb","2 ppb","11 ppb"),
  period=c("P1"),
  day=c(1,2,3,4,5)
))

b<-data.frame(expand.grid(
  treatment=c("0 ppb","2 ppb","11 ppb"),
  period=c("P2"),
  day=c(6,7,8,9,10)
))

b<-rbind(a,b)

#predict the average values predicted by the model:
b$volume   <- predict(TS3, b,type="response",re.form=NA)

#calculate bootstrapped CI for model paramaters
#sort the dataframe by treatment and day
b = with(b, b[order(treatment, day),]) 

mm     <-model.matrix(~day*treatment+period,data=b)      # model matrix
predFun<-function(.) mm%*%fixef(.)                       # predicttions for the fixed effects?  
bb     <-bootMer(TS3, FUN=predFun,nsim=200)              #bootstrapping for CI, do this 200 times
bb_se  <-apply(bb$t,2,function(x) x[order(x)][c(5,195)])
b$LC   <-bb_se[1,]
b$UC   <-bb_se[2,] 
b$pred <-predict(TS3,newdata=b,re.form=~0)


#Plot Fig 1d the raw data, the mean volume consumed from model TS3, and the CI

ggplot(data=volume, aes(y = tr_vol, x = day))+
  theme_bw()+
  geom_line(  data=volume, aes(y = tr_vol,
                                 x = day,
                                 group = colony,
                                 colour = colony
  ), alpha = 0.4)+
  
  geom_point(data=volume, aes(y = tr_vol,
                                x = day,
                                group = colony,
                                colour = colony
  ), alpha = 0.4)+
  geom_line  (data = b, aes(x = day, y = pred), size = 1)+
  geom_ribbon(data = b, aes(x = day,ymin =(b$LC), 
                                    ymax =(b$UC)
   ),inherit.aes=FALSE, alpha = 0.4)+
  facet_wrap( ~ treatment)+
  scale_x_continuous(breaks=seq(0, 10, 1))+   
  ylab("Sucrose consumption (ml)")+ 
  xlab("Experimental day")+ 
  theme(aspect.ratio=6/3)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=24),
        axis.line = element_line(colour = "black"),
        legend.position="none")
