##### 2015 Lab copepod severity #####

# Updated: 3-20-20

data<-read.csv("copepod_sev_3-20-20.csv")  
head(data)
length(data$time) # 40 - half of lab samples are in copepod or co-infection treatment

# Factors
data$colony<-as.factor(data$colony)
class(data$colony)
class(data$copYN)
class(data$time)
class(data$sex)
class(data$coinfectYN)

# response variable:
class(data$cop_sev)

###### CHECKING PLOTS 
plot(time~colony,data=data)
plot(coinfectYN~colony,data=data)
plot(coinfectYN~slide,data=data)

####### EXPLORATORY PLOTS 
plot(cop_sev~time,data=data)
plot(cop_sev~sex,data=data)
plot(cop_sev~colony,data=data) # lots of variation
plot(cop_sev~coinfectYN,data=data) # 

tapply(data$cop_sev,data$coinfectYN,FUN=mean) # 6.35 vs. 4.9 - not very different

######################
###### LINEAR MODELS 
####################

# Distribution: Use poisson- it's count data; but overdispersed 

######### MODELS #########

## Poisson ##
library(lme4) # glmer in lme4; REML = True at the very end with only best model

copsev_null<-glmer(cop_sev~1+(1|colony),family="poisson",data=data) 
copsev_1<-glmer(cop_sev~coinfectYN+(1|colony),family="poisson",data=data) 
copsev_2<-glmer(cop_sev~coinfectYN*time+(1|colony),family="poisson",data=data) 

library(AICcmodavg)
amoeb_list<-list(copsev_null,copsev_1,copsev_2)
aictab(amoeb_list,modnames=c("null","1Coinfect","2Coinfect*Time"),second.ord=FALSE)

# Null is the best
library(lmtest)
lrtest(copsev_1,copsev_null)

library("blmeco") 
dispersion_glmer(copsev_1) # should be between 0.75 and 1.4

# Other diagnostics
qqnorm(resid(copsev_1), main="normal qq-plot, residuals")
qqline(resid(copsev_1)) # good

plot(fitted(copsev_1), resid(copsev_1)) #residuals vs fitted
abline(h=0)

overdisp_fun_poisson <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun_poisson(copsev_1) # overdispersed

copsev_1g<-glm(cop_sev~coinfectYN,family="poisson",data=data) 
summary(copsev_1g) # also overdispersed

## Negative binomial##

library(lme4) # glmer in lme4; REML = True at the very end with only best model

copsev_null<-glmer.nb(cop_sev~1+(1|colony),data=data) 

copsev_1<-glmer.nb(cop_sev~coinfectYN+(1|colony),data=data) 

copsev_2<-glmer.nb(cop_sev~coinfectYN*time+(1|colony),data=data) 

library(AICcmodavg)
amoeb_list<-list(copsev_null,copsev_1,copsev_2)
aictab(amoeb_list,modnames=c("null","1Coinfect","2Coinfect*Time"),second.ord=FALSE)



#### PLOTS #####

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


#### Plot of 40 colonies ######
data_summarize<-summarySE(data,measurevar="cop_sev", groupvars=c("coinfectYN"))

library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize, aes(x=coinfectYN, y=cop_sev,fill="red")) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=cop_sev-2*se, ymax=cop_sev+2*se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("Copepod Severity")+xlab("Co-infection")

