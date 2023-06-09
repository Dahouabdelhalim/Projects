##Logistic regressions for urchin survival in lobster functional response experiments

library(ggplot2)
############################################
#Logistic regression of urchin survival by test diameter, purple urchins only

#load dataset: Survival by TD Logstic Reg Data.txt
#First column is urchin TD, second column is survival,  1=alive   0=dead

UrchTDsurv<-Survival_by_TD_Logistic_Reg_Data
colnames(UrchTDsurv) <- c("TD","Survival") #rename columns
UrchTDsurv$Treatment<-c("PurpleOnly")  #for later analysis

#Model for purple urchin only trials
glm.LogTD<-glm(Survival~TD, family=binomial(logit),data=UrchTDsurv)
summary(glm.LogTD) ##OUTPUT for purple only size-dependent mortality
par(mar=c(5,5,.5,.5))
plot(Survival~TD, data=UrchTDsurv, xlim=c(15,110), xlab="Test diameter (mm)",
     ylab="Survival probability", las=1, bty="L")
curve(predict(glm.LogTD, data.frame(TD=x), type="resp"), 
      add=TRUE, lwd=2)

#create 95% CI for this model fit
predict.data<-with(UrchTDsurv, data.frame(TD=seq(15,86, length=100)))#seq(min(TD), max(TD), length=100))) #set up data on which to predict
predictions<-predict(glm.LogTD, newdata=predict.data, type="link", se.fit=T)  #now a list w/components 'fit' and 'se.fit'
critval<-1.96   #from Gaussian 95% CI 
upr<-predictions$fit + (critval*predictions$se.fit)
lwr<-predictions$fit - (critval*predictions$se.fit)
fit<-predictions$fit

#now apply the inverse of (logit) link function
fit2<-glm.LogTD$family$linkinv(fit)
upr2<-glm.LogTD$family$linkinv(upr)
lwr2<-glm.LogTD$family$linkinv(lwr)

#compare with ggplot2
predict.data$lwr<- lwr2 
predict.data$upr <- upr2 
ggplot(data=UrchTDsurv, mapping=aes(x=TD,y=Survival)) + 
  geom_point() + 
  stat_smooth(method="glm", method.args=list(family=binomial)) + 
  geom_line(data=predict.data, mapping=aes(x=TD, y=upr), col="red") + 
  geom_line(data=predict.data, mapping=aes(x=TD, y=lwr), col="red")

#ggplot2 logistic regression: purple only
purple<-ggplot(UrchTDsurv, aes(x=TD, y=Survival))+
  geom_point(position=position_jitter(width=0.15, height=0.02), alpha=0.5, shape=21, size=2.5, col="darkorchid4")+
  stat_smooth(method="glm", method.args=list(family=binomial), color="black") +
  ylab("Survival probability") + xlab("Test diameter (mm)")
purple

##Figure S1A ##Final version
purple2<- purple + theme_bw() + scale_x_continuous(expand=c(0,0), limits=c(10,110)) +
  theme(panel.border = element_blank(), axis.line.y=element_line(color="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(axis.line=element_line(color="black")) +
  theme(axis.text.x=element_text(color="black", size=14)) +
  theme(axis.text.y=element_text(color="black", size=14)) +
  theme(axis.title.x=element_text(size=16)) +
  theme(axis.title.y=element_text(size=16))

purple2


################################################################################
#Logistic regression of urchin survival by test diameter, purple + red urchins
################################################################################
#load dataset: Survival by TD Logstic Reg DataBoth2.txt
#First column is urchin TD, second column is survival,   1=alive   0=dead

UrchTDsurvBoth<-Survival_by_TD_Logistic_Reg_Data_Both2
colnames(UrchTDsurvBoth) <- c("TD","Survival") #rename columns
UrchTDsurvBoth$Treatment<-c("Both") #for later analysis

#Model for both urchin species together: logistic regression
glm.LogTDboth<-glm(Survival~TD, family=binomial(logit),data=UrchTDsurvBoth)
summary(glm.LogTDboth)

plot(Survival~TD, data=UrchTDsurvBoth, xlim=c(15,110), xlab="Test diameter (mm)",
     ylab="Survival probability", col="blue")
curve(predict(glm.LogTDboth, data.frame(TD=x), type="resp"), 
      add=TRUE, col="black", lwd=2)

#create 95% CI for this model fit
predict.data<-with(UrchTDsurvBoth, data.frame(TD=seq(min(TD), max(TD), length=100))) #set up data on which to predict
predictions<-predict(glm.LogTDboth, newdata=predict.data, type="link", se.fit=T)  #now a list w/components 'fit' and 'se.fit'
uprb<-predictions$fit + (critval*predictions$se.fit)
lwrb<-predictions$fit - (critval*predictions$se.fit)
fitb<-predictions$fit

#now apply the inverse of (logit) link function
fitb2<-glm.LogTD$family$linkinv(fitb)
uprb2<-glm.LogTD$family$linkinv(uprb)
lwrb2<-glm.LogTD$family$linkinv(lwrb)

#ggplot logistic regression
purpred<-ggplot(UrchTDsurvBoth, aes(x=TD, y=Survival))+
  geom_point(position=position_jitter(width=0.15, height=0.02), alpha=0.5, shape=21, size=2.5, col="firebrick")+
  stat_smooth(method="glm", method.args=list(family=binomial), color="black") +
  ylab("Survival probability") + xlab("Test diameter (mm)")
purpred

##Plot for Fig. S1B  ###Final version
purpred2<- purpred + theme_bw() + scale_x_continuous(expand=c(0,0), limits=c(10,110)) +
  theme(panel.border = element_blank(), axis.line.y=element_line(color="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(axis.line=element_line(color="black")) +
  theme(axis.text.x=element_text(color="black", size=14)) +
  theme(axis.text.y=element_text(color="black", size=14)) +
  theme(axis.title.x=element_text(size=16)) +
  theme(axis.title.y=element_text(size=16))
purpred2  


#########################################################
#######Overall test of "treatment"######################
###ie, are purple only and both urchin trials different#
########################################################
#stack UrchTDsurv and UrchTDsurvBoth
bigD<-rbind(UrchTDsurv,UrchTDsurvBoth)
str(bigD); head(bigD, 8)
bigGLM<-glm(Survival~TD + Treatment + TD*Treatment, family=binomial(logit),data=bigD)
summary(bigGLM) #significant interaction btw. TD and treatment
