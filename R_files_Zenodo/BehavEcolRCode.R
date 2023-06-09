#Read in GRTI Data
GRTI <- read.csv("C:/Users/Jesee/Dropbox/DEJU-GRTI.Wakeup/CSV.Files/GreatTitRData.csv")

#GRTI Repeatability with Confidence Intervals
rep<-lmer(wakeup ~ inc.day2 + (1|nestid), data=GRTI)
summary(rep, correlation=T)
confint(rep)

#Linear Mixed Effects Model-GRTI
GRTIModel <-lmer(wakeup ~ egg1 + inc.day2 + apr.date + (1|nestid), data=GRTI)
summary(GRTIModel, correlation=T)
anova(GRTIModel)

#Read in WWJU Data
WWJU <- read.csv("C:/Users/Jesee/Dropbox/DEJU-GRTI.Wakeup/CSV.Files/JuncoRData.csv")

#WWJU Repeatability with Confidence Intervals
rep2 <-lmer(wakeup ~ incdate2 + (1|nestid), data=WWJU)
summary(rep2, correlation=T)
confint(rep2)

#Linear Mixed Effects Model With Year as RE
WWJUModel <-lmer(wakeup ~ egg1 + incdate2 + date + (1|nestid) + (1|Year), data=WWJU)
summary(WWJUModel, correlation=T)
anova(WWJUModel)

#Light Intensity at sunrise
WWJUIntensity<-lmer(wakeup ~ intensity + egg1 + incdate2 + date + (1|nestid), data=WWJU)
summary(WWJUIntensity, correlation=T)
anova(WWJUIntensity)


#Linear Mixed Effects Model w/Amb Temp & Year as RE
WWJUModelTemp <-lmer(wakeup ~ egg1 + incdate2 + date + ambtemp + (1|nestid) + (1|Year), data=WWJU)
summary(WWJUModelTemp, correlation=T)







#Read in WWJU Data
WWJU2 <- read.csv("C:/Users/Jesee/Dropbox/DEJU-GRTI.Wakeup/CSV.Files/JuncoBivariateData.csv")

#Linear Mixed Effects Model With Year as RE
WWJUBi <-MCMCglmm(cbind(wakeup, egg1) ~ trait:date + trait:incdate - 1, random = ~units:nestid, rcov = ~us(trait):units, prior = prior, family = c("gaussian", "gaussian"), data=WWJU2)
summary(WWJUBi, correlation=T)

#Linear Mixed Effects Model With Year as RE
WWJUBi <-lm(cbind(wakeup, egg1) ~ date + incdate + (1|nestid) - 1, data=WWJU2)
summary(WWJUBi, corr=TRUE)
