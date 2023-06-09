##### 2015 Lab amoebocyte analyses #####

# 3-23-20 

data<-read.csv("Amoeb_data_2015_3-23-20.csv")  
head(data)
length(data$sample) # 411 
levels(data$sample) # 80: 10 colonies * 4 treatments * 2 time points

# Omitting colony 9 due to pre-existing infection
head(data)
data_no9<-subset(data,colony!="9")
head(data_no9)
summary(data_no9$colony)
length(data_no9$sample) # 374

levels(data_no9$sample) 

# Factors
class(data_no9$perc_amoeb)
data_no9$colony<-as.factor(data_no9$colony)
class(data_no9$colony)
data_no9$treatment<-as.factor(data_no9$treatment)
class(data_no9$treatment)
data_no9$time<-as.factor(data_no9$time)
class(data_no9$time) 
class(data_no9$onepath) # onepath = "parasite_single" in text

###### LINEAR MODELS 

# Normality: 
par(mfrow=c(1,1))
qqnorm(data_no9$perc_amoeb)
qqline(data_no9$perc_amoeb) # 
shapiro.test(data_no9$perc_amoeb) # not quite normal

# sqrt transformation works
qqnorm(sqrt(data_no9$perc_amoeb)) 
qqline(sqrt(data_no9$perc_amoeb)) # Yes 
shapiro.test(sqrt(data_no9$perc_amoeb)) # normal

sqrtperc9<-sqrt(data_no9$perc_amoeb)
qqnorm(sqrtperc9)
qqline(sqrtperc9)

qqnorm(data_no9$sqrtperc)
qqline(data_no9$sqrtperc)

# CONTRASTS
contrasts(data_no9$treatment) # treatment all vs. N
contrasts(data_no9$treatment) <- contr.sum # control "N" vs each treatment 

contrasts(data_no9$path_YN) 
contrasts(data_no9$colony) <- contr.sum 
contrasts(data_no9$time) <- contr.sum 
contrasts(data_no9$trtmt_M) <- cbind(c(-1,1)) 
contrasts(data_no9$trtmt_M) 
contrasts(data_no9$trtmt_A) <-  cbind(c(-1,1))
contrasts(data_no9$trtmt_A)

######### MODELS #########

library(lme4) # glmer in lme4

amoeb_null<-lmer(sqrtperc9~1+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_1<-lmer(sqrtperc9~trtmt_M*trtmt_A*time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_2<-lmer(sqrtperc9~trtmt_M*trtmt_A+time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_3<-lmer(sqrtperc9~trtmt_M+trtmt_A*time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_4<-lmer(sqrtperc9~trtmt_M*time+trtmt_A+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_5<-lmer(sqrtperc9~trtmt_M+trtmt_A+time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_6<-lmer(sqrtperc9~trtmt_M+trtmt_A+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_7<-lmer(sqrtperc9~trtmt_M+time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_8<-lmer(sqrtperc9~trtmt_A+time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_9<-lmer(sqrtperc9~trtmt_M+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_10<-lmer(sqrtperc9~trtmt_A+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_11<-lmer(sqrtperc9~path_YN*time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_12<-lmer(sqrtperc9~path_YN+time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_13<-lmer(sqrtperc9~path_YN+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_14<-lmer(sqrtperc9~onepath+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_15<-lmer(sqrtperc9~onepath+time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

amoeb_16<-lmer(sqrtperc9~onepath*time+(1|colony)+(1|sample),data=data_no9,REML=FALSE) 

library(AICcmodavg)
amoeb_list<-list(amoeb_null,amoeb_1,amoeb_2,amoeb_3,amoeb_4,amoeb_5,amoeb_6,amoeb_7,amoeb_8,amoeb_9,amoeb_10,amoeb_11,amoeb_12,amoeb_13,amoeb_14,amoeb_15,amoeb_16)
aictab(amoeb_list,modnames=c("null","1M*A*time","2M*A+time","3M+A*time","4M*time+A","5M+A+time","6M+A","7M+time","8A+time","9M","10A","11pathYN*time","12pathYN+time","13pathYN","14onepath","15onepath+time","16onepath*time"),second.ord=FALSE)

# Best model: 9M (just M again); actually pretty consistent with perc_amoeb_norm results 

library(lmtest)
lrtest(amoeb_null,amoeb_9) # 9 is better than the null & by AIC


# REML 
library(lmerTest)

amoeb_9_REML<-lmer(sqrtperc9~trtmt_M+(1|colony)+(1|sample),data=data_no9,REML=TRUE) 
summary(amoeb_9_REML)

library(MuMIn)
confint(amoeb_9_REML,level=0.95, method="Wald") # 

hist(resid(amoeb_9_REML)) # sufficiently normal
plot(predict(amoeb_9_REML),resid(amoeb_9_REML)) # good

tapply(data_no9$perc_amoeb,data_no9$trtmt_M,FUN=mean) # with and without copepod; 20.22 to 21.877
(21.877-20.226)/20.226 # 8.16%

# residuals vs. the independent variable also looks good 
library(ggplot2)
ggplot(data.frame(x1=data_no9$treatment,pearson=residuals(amoeb_9_REML,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw()




###### BAR PLOT 

# summarySE - http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/ .....5-23-17
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


## Lines: treatment 1SE ###
data_no9$treatment <- factor(data_no9$treatment, levels=c("N","A","M","MA"))

data_se <- summarySE(data_no9, measurevar="perc_amoeb", groupvars=c("treatment"))
data_se

ggplot(data_se, aes(x=treatment, y=perc_amoeb)) + 
  geom_errorbar(aes(ymin=perc_amoeb-1*se, ymax=perc_amoeb+1*se), width=.1) +
  geom_line() +
  geom_point() +
  xlab("Treatment")+
  ylab("Amoebocyte density (%)")+
  # labs(title = "% Change in amoebocytes by treatment")+ 
  ylim(17,25) # 



##### Plot rand. effects ###

data_no9$trtmt_M <- as.factor(data_no9$trtmt_M)


library(nlme)
amoeb_9l <- lme(sqrtperc9 ~ trtmt_M, data = data_no9, random = ~ 1|colony)

newdat <- expand.grid(trtmt_M=unique(data_no9$trtmt_M))

library(ggplot2)
p <- ggplot(data_no9, aes(x=trtmt_M, y=sqrtperc9, group=1)) +
  geom_point(size=1) +
  geom_line(aes(y=predict(amoeb_9l), group=colony,color=colony)) +
  geom_line(data=newdat, aes(y=predict(amoeb_9l, level=0, newdata=newdat))) +
  scale_size_manual(name="Predictions", values=c("Colonies"=0.5)) +
  theme_bw(base_size=15) +
  ylab("Amoebocyte density (%area)")+
  xlab("Copepod presence")
print(p) 


