##########################################################################################
# Script to analyse data and produce the results presetned in                            #
# " Environmental drivers of Sphagnum growth in peatlands across the Holarctic region"   #
#                                                                                        #
# Contact: Gustaf.Granath@gmail.com                                                      #
##########################################################################################
library(ggplot2)
library(gridExtra)
# Load extra functions (R2, plot etc) ####
source("gsp_extra_func.R")
setwd("~/Dropbox/gsp_r/rev_analyses")
# Load data
dat <- read.csv("Bengtsson_etal_gsp_prod.csv")

# Remove two site in China (inconsistent sampling) and one site that do not have N concentration
rem <- which(dat$Site == "Hani"| dat$Site == "Mangui"| dat$Site == "Verh-Tarka")
dat <- dat[-rem, ]

#Make df in long format
# and fix variables
library(tidyr)
library(dplyr)
data_long <- gather(dat, year, LI, c(LI13, LI14), factor_key=TRUE)
levels(data_long$year)<-c("yr2013", "yr2014")

data_long$grow_days <- c(dat$days2013, dat$days2014)
data_long$LI_d <- c(dat$LI_per_day2013, dat$LI_per_day2014)
data_long$prod <- c(dat$Prod13, dat$Prod14)
data_long$prod_d <- c(dat$prod_per_day2013, dat$prod_per_day2014)
data_long$temp <- c(dat$temp13, dat$temp14)
data_long$evap <- c(dat$evap13, dat$evap14)
data_long$evap_d <- c(dat$ev_d13, dat$ev_d14)
data_long$prec <- c(dat$precip13, dat$precip14)
data_long$prec_d <- c(dat$prec_d13, dat$prec_d14)
data_long$norain <- c(dat$norain13.mean, dat$norain14.mean)
data_long$par <- c(dat$par13, dat$par14)
data_long$par_d <- c(dat$par_d13, dat$par_d14)
data_long$cover <- c(dat$cover13, dat$cover14)
data_long$hwt <- c(dat$HWT13, dat$HWT14)
data_long$BD <- c(dat$BD13, dat$BD14)
data_long$numden <- c(dat$num.den13, dat$num.den14)
data_long$ev_pre_d <- data_long$prec_d - data_long$evap_d # precipitation minus evaporation per day
data_long$ev_pre <- data_long$prec - data_long$evap # precipitation minus evaporation 

# Add site mean nitrogen concentration and NP ratio
data_long$id <- 1:nrow(data_long)
data_long <- data_long %>% 
  group_by(Site, Species) %>% 
  mutate(Nmean = mean(N_per, na.rm=TRUE),
         NPmean = mean(NP, na.rm=TRUE)) %>% arrange(id) %>% as.data.frame()

# Loop to: scale for each response variable after removing NA rows,
# and save each data frame in a list named 'resp. for later use
responses <- c("LI_d", "LI", "prod_d", "prod", "BD")
resp <-list()
for (i in 1:length(responses)) {
temp <- data_long[complete.cases(data_long[ , c(responses[i], "prec_d", "evap_d", "par_d",
                                                       "prec", "temp", "evap", "par", "hwt", "cover", "norain", 
                                                       "ndep_Lam13", "Nmean", "NPmean", "numden"),]),]

temp.sc <- scale(temp[ , c("prec_d", "evap_d", "par_d", "ev_pre", "ev_pre_d",
                                  "prec", "temp", "evap", "par", "hwt", "cover", "norain", 
                                  "ndep_Lam13", "Nmean", "NPmean", "numden"),])
colnames(temp.sc) <- c("pr_d", "ev_d", "pa_d","ev_pre","ev_pre_d",
                      "pr", "tem", "ev", "pa", "wt", "cov", "nor", "ndeL", "Nm", "NPm","numd")
resp[[i]] <- as.data.frame(cbind(temp, temp.sc))
}
names(resp) <- responses # name list

# Run models ####

# Mixed models. Sites as random, 
# Using lmer and Powertransform to transform data so residuals 
# are normaly distributed, predictors are scaled
library(lme4)
library(car) #needed for powerTransform()
library(nlme)

#___LI per day####
# first get lambda for Box-cox transformation, response is too small, so *100
# make all values positive for easier transformation
min.val <- min(resp[["LI_d"]]$LI_d, na.rm=T)*-1+0.0001 # get value to add for positive response
mod1 <- lmer(((LI_d+min.val)*100) ~ Species+year+ Species*pr_d  + Species*tem + Species*ev_d + 
               Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm+ Species*numd+
             (year-1|Site),
            data = resp[["LI_d"]])
plot(mod1) #increasing variance
summary(l1<-powerTransform(mod1, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.30

#transform response
resp[["LI_d"]]$LId_trans<-((resp[["LI_d"]]$LI_d+min.val)*100)^(l1$lambda)
mod1.2<-lmer(LId_trans ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + Species*pa_d + 
               Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+Species*numd+
               (year-1|Site),
         data = resp[["LI_d"]])

        #data_long.na2<- mod1.2@frame # save data w/o NAs
#Anova to check if interactions are interesting library(car)
Anova(mod1.2, test.statistic = "F")
summary(mod1.2)
plot(mod1.2)

#to get coefficients per species
mod1.2.1<-lmer(LId_trans ~ Species+year+Species/pr_d  + Species/tem + Species/ev_d + Species/pa_d + 
                 Species/wt + Species/cov + 
               Species/nor + Species/ndeL + Species/Nm+ Species/numd+
               (year-1|Site),
             data = resp[["LI_d"]])
summary(mod1.2.1)

#LI_d null model
mod1n<-lmer((LId_trans)~1 + (1|Site),
            data = resp[["LI_d"]])
#plot(mod1n)
summary(mod1n)

#______Excluding outliers from model with 2-way interactions ####
outliersLId<-resp[["LI_d"]][which(residuals(mod1.2)<(-1.5)),]#which are the outliers?
data_long1<-resp[["LI_d"]][which(residuals(mod1.2)>(-1.5)),]
modxx1<-lmer((LId_trans) ~ Species/year+ Species/pr_d  + Species/tem + Species/ev_d + Species/pa_d + 
               Species/wt + Species/cov + 
              Species/nor + Species/nde + Species/Nm + (1|Site), 
             data = data_long1, na.action = na.omit)

modx1<-lmer((LId_trans) ~ Species*year+ Species*pr_d  + Species*tem + Species*ev_d + Species*pa_d + 
              Species*wt + Species*cov + 
              Species*nor + Species*nde + Species*Nm + (1|Site),
            data = data_long1, na.action = na.omit)
plot(modx1)

#Anova to check if interactions are interesting library(car)
summary(modxx1)
Anova(modx1, test.statistic = "F")
summary(modx1)

#Model with only main effects
mod1.5<-lmer(LId_trans ~ Species + year + pr_d + tem +ev_d + pa_d + wt + cov + nor+ ndeL + Nm + numd + 
              (year-1|Site), 
               data = resp[["LI_d"]])
plot(mod1.5)
summary(mod1.5)
Anova(mod1.5, test.statistic = "F")

# Check residuals
res.test <- mod4.2@frame
res.test$resid <- residuals(mod4.2,type="pearson")
# numeric predictors
res.test %>%
  select_if(is.numeric) %>%
  gather(preds, value, 2:(ncol(.)-1)) %>% 
  ggplot(aes(x=value, y=resid)) +
  geom_point() +
  facet_wrap(~ preds, scales = "free")
# categorical predictors also
res.test %>%
  select(-Site) %>%
  gather(preds, value, 2:(ncol(.)-1)) %>% 
  ggplot(aes(x=value, y=resid)) +
  geom_point() +
  facet_wrap(~ preds, scales = "free")

# check influential points
ggplot(data.frame(lev=hatvalues(mod4.2),pearson=residuals(mod4.2,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()

#identify points with high leverage, set a cutoff
res.test[which(hatvalues(mod2.2) >= .25),]


#______Excluding outliers from model with no interactions####
outliersLId2 <- resp[["LI_d"]][which(residuals(mod1.5)<(-1.5)),]

data_long2<-resp[["LI_d"]][which(residuals(mod1.5)>(-1.5)),]

modx2<-lmer((LId_trans) ~ Species + year + pr_d  + tem + ev_d + pa_d + wt + cov + nor+ nde + Nm+
              (1|Site), 
            data = data_long2, na.action = na.omit)
plot(modx2)
summary(modx2)
Anova(modx2, test.statistic = "F")##end outlier test


#______Prec-evap mixed model####
# Test if a index works better than evap and precip as separate variables
min.val <- min(resp[["LI_d"]]$LI_d, na.rm=T)*-1+0.0001 # get value to add for positive response
mod1_evpre <- lmer(((LI_d+min.val)*100) ~ Species+year+ Species*ev_pre_d  + Species*tem + 
               Species*pa_d + Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+ Species*numd+
               (year-1|Site),
             data = resp[["LI_d"]])
plot(mod1_evpre) #increasing variance
summary(l1<-powerTransform(mod1_evpre, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.32

#transform response
resp[["LI_d"]]$LId_trans<-((resp[["LI_d"]]$LI_d+min.val)*100)^(l1$lambda)
mod2_evpre<-lmer(LId_trans ~ Species+year+Species*ev_pre_d  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+ Species*numd+
               (year-1|Site),
             data = resp[["LI_d"]])
#Anova to check if interactions are interesting library(car)
Anova(mod2_evpre, test.statistic = "F")

#to get coefficients per species
mod2.1_evpre<-lmer(LId_trans ~ Species+year+Species/ev_pre_d  + Species/tem +  Species/pa_d + 
                 Species/wt + Species/cov + 
                 Species/nor + Species/ndeL + Species/Nm+ Species/numd+
                 (year-1|Site),
               data = resp[["LI_d"]])
summary(mod2.1_evpre)

#Model with only main effects
mod3_evpre<-lmer(LId_trans ~ Species + year + ev_pre_d + tem + pa_d + wt + cov + nor+ ndeL + Nm + numd +
               (year-1|Site), 
             data = resp[["LI_d"]])
plot(mod3_evpre)
summary(mod3_evpre)
Anova(mod3_evpre, test.statistic = "F")


#___NPP per day####
# first get lambda for Box-cox transformation, response is too small, so *100
# make all values positive for easier transformation
min.val <- min(resp[["prod_d"]]$prod_d)*-1+0.001 # get value to add for positive response
mod2<-lmer(( (prod_d+min.val)*100) ~ Species+year+ Species*pr_d  + Species*tem + Species*ev_d + 
             Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm+Species*numd+
             (year-1|Site), data = resp[["prod_d"]])

plot(mod2)
summary(l2<-powerTransform(mod2, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.32

#transform response
resp[["prod_d"]]$prod.d_trans <- ((resp[["prod_d"]]$prod_d+min.val)*100)^(l2$lambda)
mod2.2<-lmer(prod.d_trans ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + 
               Species*pa_d + Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+Species*numd+
               (year-1|Site), 
             data= resp[["prod_d"]])

#data_long.prod.d.na<- mod2.2@frame # save data w/o NAs
plot(mod2.2)
#to get coefficients
mod2.2.1<-lmer(prod.d_trans ~ Species+year+Species/pr_d  + Species/tem + Species/ev_d + Species/pa_d + 
                 Species/wt + Species/cov + 
                 Species/nor + Species/ndeL + Species/Nm+Species/numd+
                 (year-1|Site), 
               data = resp[["prod_d"]])
summary(mod2.2.1)
Anova(mod2.2, test.statistic = "F")
summary(mod2.2)

#prod_d null model
mod2n<-lmer((prod.d_trans)~1 + (1|Site),
            data = resp[["prod_d"]])
#plot(mod2n)
#summary(mod2n)

#Run without interactions to get coefficients
mod2.5<-lmer(prod.d_trans ~ Species + year + pr_d  + tem + ev_d + pa_d + wt + cov + nor+ ndeL + Nm+numd+
               (year-1|Site), 
             data = resp[["prod_d"]])
plot(mod2.5)
summary(mod2.5)
Anova(mod2.5, test.statistic = "F")

#_______Prec-evap Mixed model####
min.val <- min(resp[["prod_d"]]$prod_d)*-1+0.001 # get value to add for positive response
modE_P<-lmer(( (prod_d+min.val)*100) ~ Species+year+ Species*ev_pre_d  + Species*tem + 
             Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm+Species*numd+
             (year-1|Site), data = resp[["prod_d"]])

plot(modE_P)
summary(l2<-powerTransform(modE_P, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.32

#transform response
resp[["prod_d"]]$prod.d_trans <- ((resp[["prod_d"]]$prod_d+min.val)*100)^(l2$lambda)

modE_P2<-lmer(prod.d_trans ~ Species+year+Species*ev_pre_d  + Species*tem + 
               Species*pa_d + Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+Species*numd+
               (year-1|Site), 
             data= resp[["prod_d"]])

Anova(modE_P2, test.statistic = "F")


modE_P4<-lmer(prod.d_trans ~ Species+year + Species/ev_pre_d  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm+Species/numd+
                (year-1|Site), data = resp[["prod_d"]])
summary(modE_P4)

#Run without interactions to get coefficients
modE_P5<-lmer(prod.d_trans ~ Species + year + ev_pre_d  + tem + pa_d + wt + cov + nor+ ndeL + Nm+ numd+
               (year-1|Site),
              data = resp[["prod_d"]])
plot(modE_P5)
summary(modE_P5)
Anova(modE_P5, test.statistic = "F")

#___LI per season ####

min.val <- min(resp[["LI"]]$LI)*-1+0.001 # get value to add for positive response
mod3<-lmer((LI+min.val) ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
            Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd+
           (year-1|Site), data = resp[["LI"]])
plot(mod3)

summary(l3<-powerTransform(mod3, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.38

resp[["LI"]]$LI_trans <- (resp[["LI"]]$LI+min.val)^(l3$lambda)

mod3.2<-lmer(LI_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd +
               (year-1|Site), data = resp[["LI"]])

# test if rain 2013 can predict 2014 LI
mod3.2.prec13<-lmer(LI_trans ~ Species +Species*scale(precip13)+ Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (1|Site), data = resp[["LI"]][resp[["LI"]]$year=="yr2014",])
summary(mod3.2.prec13)
#

#LI only PAR variable
mod3.2par<-lmer(LI_trans ~ Species+year + Species*pa +
               (year-1|Site), data = resp[["LI"]])
summary(mod3.2par)

#LI season null model
mod3n<-lmer((LI_trans)~1 + (year-1|Site),
            data = resp[["LI"]])
summary(mod3.2)
summary(mod3n)

#Calculate within and between site variance 
mod3n<-lmer((LI_trans)~1 + (1|Site),
            data = resp[["LI"]],subset = Species=="S.magellanicum")

# code to test if random slope is better
mod3.2b<-lmer(LI_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (1|Site), data = resp[["LI"]])
chistat <- anova(mod3.2b, mod3.2)[2,"Chisq"] 
0.5 * pchisq(chistat, 1, lower.tail=FALSE) + 0.5 * pchisq(chistat, 2, lower.tail=FALSE)

plot(mod3.2)

# get coefs
mod3.2.1<-lmer(LI_trans ~ Species+year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +Species/numd +
               (year-1|Site), data = resp[["LI"]])

summary(mod3.2.1)
Anova(mod3.2, test.statistic = "F")
summary(mod3.2)


#no interactions
mod3.4<-lmer(LI_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + ndeL + Nm + numd
             + (year-1|Site), 
             data = resp[["LI"]])

plot(mod3.4)
summary(mod3.4)
Anova(mod3.4, test.statistic = "F")

#______Excluding outliers from model with no interactions####
outliersLIseas<-resp[["LI"]][which(residuals(mod3.2)<(-1.5)),]

data_long2<-resp[["LI"]][which(residuals(mod3.2)>(-1.5)),]

modx3.1<-lmer(LI_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                (1|Site), data = data_long2, na.action = na.omit)
plot(modx3.1)
summary(modx3.1)
Anova(modx3.1, test.statistic = "F")

#no interactions
modx3<-lmer(LI_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + nde + Nm +
              (1|Site), data = data_long2, na.action = na.omit)
plot(modx3)
summary(modx3)
Anova(modx3, test.statistic = "F")


#______SEM models####

# save specific data frame
li.seas <- resp[["LI"]][, c("LI_trans", "Species","year","pr","tem","ev",
                            "pa","wt","nor","ndeL","Nm","cov", "numd", 
                            "Site")]
# If we want standardised coefs we need LI standardised as well
li.seas$LI_trans_st <- scale(li.seas$LI_trans)[,1]

                          
####
# use piecewiseSEM package etc
library(piecewiseSEM)

#factors dont work so need to have them as numerical
li.seas$Species <- as.numeric(li.seas$Species)-1
li.seas$year <- as.numeric(li.seas$year)-1

sems <- psem(
  lme(wt ~ Species + year + pr +tem + ev +  
        nor, 
      random = ~1|Site, data = li.seas, 
      weights=varExp(), control = list(maxIter=1500)),
  lme(cov ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
         Species*wt + Species*ndeL, 
       random = ~1|Site, data = li.seas, 
       weights=varExp()),
  lme(Nm ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
         Species*wt + Species*nor + Species*ndeL + Species*cov, 
       random = ~year-1|Site, data = li.seas, 
       weights=varExp(), na.action = na.omit),
  lme(LI_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
         Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm,
        random = ~year-1|Site, data = li.seas)
)
dd <- summary(sems)
# for plotting
#dd <- summary(sems.per.sp)
#cc <- dd$coefficients
#colnames(cc)[8] <- "estimate"

####

# SEM customised
wt.liseas <- lme(wt ~ year + pr + tem + ev + nor, 
                  random = ~1|Site, data = li.seas, 
                  weights=varExp())
cov.liseas <-  lme(cov ~ year + Species/pr + Species/tem + Species/ev +  
                      Species/wt + Species/ndeL, 
                    random = ~1|Site, data = li.seas, 
                    weights=varExp(), na.action = na.omit, control = list(maxIter = 500))
nm.liseas <-  lme(Nm ~ #Species/pa + Species/wt + 
            Species/ndeL + Species/cov, 
                   random = ~1|Site, data = li.seas, 
                   weights=varExp())
nd.liseas <-lme(numd ~ year + Species/pr + Species/tem + Species/ev +  
                   Species/wt + Species/cov,
                 random = ~ 1|Site, weights=varExp(), li.seas[!(is.na(li.seas$numd)),])
li.liseas <-  lme(LI_trans_st ~ year + Species/pr + Species/tem + Species/ev + Species/pa + 
                      Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm + 
                    Species/numd,
                    random = ~year-1|Site, data = li.seas[!(is.na(li.seas$numd)),])
#MuMIn::r.squaredGLMM(nd.nppseas)

# save coefs from models
wt.tab <- data.frame(summary(wt.liseas)$tTable[-1,])
wt.tab$Predictor <- rownames(wt.tab)
cov.tab <- data.frame(summary(cov.liseas)$tTable[-1,])
cov.tab$Predictor <- rownames(cov.tab)
nm.tab <- data.frame(summary(nm.liseas)$tTable[-1,])
nm.tab$Predictor <- rownames(nm.tab)
nd.tab <- data.frame(summary(nd.liseas)$tTable[-1,])
nd.tab$Predictor <- rownames(nd.tab)
li.tab <- data.frame(summary(li.liseas)$tTable[-1,])
li.tab$Predictor <- rownames(li.tab)
# put together data frame with coefs and change column names
sem.li <- data.frame(Response = c(rep("wt",nrow(wt.tab)),rep("cov",nrow(cov.tab)),
                                   rep("Nm",nrow(nm.tab)),rep("numd",nrow(nd.tab)),
                                   rep("LI",nrow(li.tab))), 
                      rbind(wt.tab, cov.tab, nm.tab, nd.tab, li.tab))
colnames(sem.li)[2] <- "estimate"
colnames(sem.li)[6] <- "P.Value"
sem.li <- sem.li[,c(1,7,2,3,4,5,6)]

# remove species and year effect and split up coefs on species
sem.li.mag <- sem.li[unique(c(grep("S.mag", sem.li$Predictor), grep("S.fus", sem.li$Predictor, invert = TRUE))),]
sem.li.mag <- sem.li.mag[!(sem.li.mag$Response == "wt"),]
sem.li.mag[2] <- lapply(sem.li.mag[2], gsub, pattern = "SpeciesS.magellanicum:", 
                         replacement = "", fixed = TRUE)
sem.li.mag <- sem.li.mag[!(sem.li.mag$Predictor == "SpeciesS.magellanicum") &
                             !(sem.li.mag$Predictor == "yearyr2014"),]

sem.li.fus <- sem.li[unique(c(grep("S.fus", sem.li$Predictor), grep("S.mag", sem.li$Predictor, invert = TRUE))),]
sem.li.fus <- sem.li.fus[!(sem.li.fus$Response == "wt"),]
sem.li.fus[2] <- lapply(sem.li.fus[2], gsub, pattern = "SpeciesS.fuscum:", 
                         replacement = "", fixed = TRUE)
sem.li.fus <- sem.li.fus[!(sem.li.fus$Predictor == "yearyr2014"),]

library(openxlsx) 
write.xlsx(sem.li, "sem.li_april2.xlsx")

# plot sem
#def.par <- par(no.readonly = TRUE) 
x11(width=6, height=12)
#pdf("LIseas.pdf",width=6, height=12)
#setEPS()
#postscript("LIplot.eps", width=6, height=12)
par(mfrow=c(2,1), xpd = NA)
sem.plot(coef.table=sem.li.mag,alpha = 0.1)
legend(-1.8,1.8,lty=c(1,2,1,1), col=c(1,1,1,2),
       legend=c("P<0.1","P>0.1","pos","neg"), xpd = TRUE, bty="n")
mtext("LI season\\n", line=2)
mtext("S.mag\\n")
sem.plot(coef.table=sem.li.fus,alpha = 0.1)
mtext("S.fuscum\\n")
dev.off()

# function to plot SEM
# modified fromm piecewiseSEM v1.2.1
sem.plot = function( modelList = NULL, data = NULL, coef.table = NULL, corr.errors = NULL, 
                      show.nonsig = TRUE, scaling = 10, alpha = 0.05, ...) {
    # Get coefficients
    if(!is.null(modelList) & !is.null(data) & is.null(coef.table))
    coef.table = sem.coefs(modelList, data, corr.errors = corr.errors, ...) else
      if(is.null(modelList) & is.null(data) & is.null(coef.table))
        stop("Please provide model list and data, or coefficient table!")

  # Prepare coef.table
  coef.table$corr.errors = grepl("~~", coef.table[, 1])
  coef.table[, 1] = gsub("~~ ", "", coef.table[, 1])
  coef.table[, 2] = gsub("~~ ", "", coef.table[, 2])

  # Strip transformations
  for(i in 1:2) coef.table[, i] = gsub(".*\\\\((.*)\\\\).*", "\\\\1", coef.table[, i])

  # Get vector of labels
  lbls = unlist(coef.table[, 1:2])
  lbls = as.character(unname(lbls[!duplicated(lbls)]))

  # # Shorten label names if necessary
  # if(any(sapply(lbls, function(x) nchar(x) > 10))) {
  # 
  #  new.lbls = gsub("a|e|i|o|u", "", lbls) } else new.lbls = lbls
  #  if(any(sapply(new.lbls, function(x) nchar(x) > 10))) 
  #    new.lbls = sapply(lbls, function(x) ifelse(nchar(x) > 10, substr(x, 1, 10), x))
  
  #names(new.lbls) = lbls
  new.lbls = lbls
  names(new.lbls) = lbls
  
    # Set graphical parameters
    par(mar = rep(5, 4), xpd = NA)

  # Initialize plot
  plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", ann = FALSE, axes = FALSE)

    # Prepare circle coef.table
    theta = seq(0, 2 * pi, length = 200)
  
    # Add to coef.table.frame
    circle = data.frame(
    x = cos(theta),
    y = sin(theta)
    )
  # Get labels
  row.n = seq(1, 200, by = 200/length(lbls))
  names(row.n) = lbls

  # Get linewidth scale
  if(is.na(scaling)) scl.fctr = rep(1, nrow(coef.table)) else {
    scaling = scaling / diff(range(coef.table[, "estimate"]))
    scl.fctr = abs(coef.table[, "estimate"]) * scaling
  }

  # Add arrows
  for(i in 1:nrow(coef.table)) {
    resp = row.n[names(row.n) == coef.table[i, 1]]
    pred = row.n[names(row.n) == coef.table[i, 2]]
    if(show.nonsig == FALSE & coef.table[i, "P.Value"] >= alpha) next else {
      arrows(
        x0 = circle[pred, "x"],
        y0 = circle[pred, "y"],
        x1 = circle[resp, "x"],
        y1 = circle[resp, "y"],
        code = ifelse(coef.table[i, "corr.errors"] == TRUE, 3, 2),
        #col = ifelse(coef.table[i, "P.Value"] < alpha & coef.table[i, "estimate"] > 0, "black", 
        #             ifelse(coef.table[i, "P.Value"] < alpha & coef.table[i, "estimate"] < 0, "red", "grey50")),
        col = ifelse(coef.table[i, "estimate"] > 0, "black", "red"),
        
        lwd = scl.fctr[i],
        lty = ifelse(coef.table[i, "P.Value"] < alpha, 1, 2)
      )
    }
  }

  # Plot labels
  for(i in row.n) {
    x = circle[i, "x"]
    y = circle[i, "y"]
    xjust = ifelse(x > 0, 0.5 - (x / 2), 0.5 + (-x / 2))
    yjust = ifelse(y > 0, 0, 1)
    legend(
      x, 
      y, 
      legend = new.lbls[names(new.lbls) == names(row.n[row.n == i])], 
      cex = 0.8,
      x.intersp = 0,
      xjust = xjust,
      yjust = yjust
    )
  }

  # Reset par
  par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = FALSE)
}

#______Prec-evap mixed model####
min.val <- min(resp[["LI"]]$LI)*-1+0.001 # get value to add for positive response
LI_P_E1<-lmer((LI+min.val) ~ Species+year + Species*ev_pre + Species*tem + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd +
             (year-1|Site), data = resp[["LI"]])
plot(LI_P_E1)

summary(l3<-powerTransform(LI_P_E1, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.34

resp[["LI"]]$LI_trans <- (resp[["LI"]]$LI+min.val)^(l3$lambda)

LI_P_E<-lmer(LI_trans ~ Species+year + Species*ev_pre  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm+Species*numd +
               (year-1|Site), data = resp[["LI"]])

plot(LI_P_E)

LI_P_E2<-lmer(LI_trans ~ Species+year + Species/ev_pre  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm+Species/numd +
                (year-1|Site), data = resp[["LI"]])
summary(LI_P_E2)
Anova(LI_P_E, test.statistic = "F")

#Run without interactions to get coefficients
LI_P_E3<-lmer(LI_trans ~ Species + year + ev_pre  + tem + pa_d + wt + cov + nor+ ndeL + Nm+ numd +
               (year-1|Site),
             data = resp[["LI"]])
plot(LI_P_E3)
summary(LI_P_E3)
Anova(LI_P_E3, test.statistic = "F")

#___NPP per season ####
min.val <- min(resp[["prod"]]$prod)*-1+0.001 # get value to add for positive response
mod4<-lmer((prod+min.val) ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd+
             (year-1|Site), data = resp[["prod"]])
plot(mod4)
summary(l4<-powerTransform(mod4, family = "bcPower"))#bcnPower neccessary with neg. values 
resp[["prod"]]$prod_trans<-(resp[["prod"]]$prod+min.val)^(l4$lambda)

mod4.2<-lmer(prod_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm+ Species*numd+
               (year-1|Site), data = resp[["prod"]])
plot(mod4.2)

# coefs
mod4.2.1<-lmer(prod_trans ~ Species+year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +Species/numd+
                 (year-1|Site), data = resp[["prod"]])

summary(mod4.2.1)
Anova(mod4.2, test.statistic = "F")
summary(mod4.2)

# test if rain 2013 can predict 2014 NPP
mod4.2.prec13<-lmer(prod_trans ~ Species + Species*scale(precip13)+ Species*tem + Species*ev + Species*pa + 
                      Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
                      (1|Site), data = resp[["prod"]][resp[["LI"]]$year=="yr2014",])
summary(mod4.2)
#



# no interactions
mod4.4<-lmer(prod_trans ~ Species + year + pr + tem + ev + pa + 
               wt + cov + nor + ndeL + Nm + numd+ (year-1|Site), 
             data = resp[["prod"]])

plot(mod4.4)
summary(mod4.4)
Anova(mod4.4, test.statistic = "F")

#prod season null model
mod4n<-lmer((prod_trans)~1 + (year-1|Site),
            data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
#plot(mod4n)
summary(mod4n)

#Calculate site variance
#Overall
mod4n<-lmer((prod_trans)~1 + (1|Site),
            data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
summary(mod4n)
mod4n<-lmer((prod_trans)~1 + (1|Site),
            data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),],subset = Species=="S.fuscum")

#______Excluding outliers from model with no interactions####
outliersprodseas<-resp[["prod"]][which(residuals(mod4.2)<(-1.5)),]

data_long4<-resp[["prod"]][which(residuals(mod4.2)>(-1.5)),]

modx4.1<-lmer(prod_trans ~ Species/year + Species/pr + Species/tem + Species/ev + 
                Species/pa + Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                (year-1|Site), 
              data = data_long4)
plot(modx4.1)
summary(modx4.1)
Anova(modx4.1, test.statistic = "F")

#no interactions
modx4<-lmer(prod_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + nde + Nm +
              (year-1|Site), 
            data = data_long4)
plot(modx4)
summary(modx4)
Anova(modx4, test.statistic = "F")#

#______SEM models####

# save specific data frame
#npp.seas <- resp[["prod"]] #[-(which(resp[["prod"]]$orgOrder ==177)),]
npp.seas <- resp[["prod"]][, c("prod_trans", "Species","year","pr","tem","ev",
                            "pa","wt","nor","ndeL","Nm","cov", "numd", 
                            "Site")]
# If we want standardised coefs we need LI standardised as well
npp.seas$prod_trans_st <- scale(npp.seas$prod_trans)[,1]

# use piecewiseSEM package
library(piecewiseSEM)
# factors have to be numeric in piecewiseSEM
npp.seas$Species <- as.numeric(npp.seas$Species)-1
npp.seas$year <- as.numeric(npp.seas$year)-1

sems.npp.seas <- psem(
  lme(wt ~ Species+year + pr +tem + ev +  
        nor, 
      random = ~year-1|Site, data = npp.seas, 
      weights=varExp(), na.action = na.omit),
  lme(cov ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
        Species*wt + Species*ndeL, 
      random = ~year-1|Site, data = npp.seas, 
      weights=varExp(), na.action = na.omit),
  lme(Nm ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
        Species*wt + Species*nor + Species*ndeL + Species*cov, 
      random = ~year-1|Site, data = npp.seas, 
      weights=varExp(), na.action = na.omit),
  lme(prod_trans_st ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
        Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm,
      random = ~year-1|Site, data = npp.seas)
)
summary(sems.npp.seas)

# SEM customised
wt.nppseas <- lme(wt ~ year + pr + tem + ev + nor, 
                 random = ~1|Site, data = npp.seas, 
                 weights=varExp())
cov.nppseas <-  lme(cov ~ year + Species/pr + Species/tem + Species/ev +  
                     Species/wt + Species/ndeL, 
                   random = ~1|Site, data = npp.seas, 
                   weights=varExp(), na.action = na.omit, control = list(maxIter = 500))
nm.nppseas <-  lme(Nm ~ #Species/pa + Species/wt + 
                    Species/ndeL + Species/cov, 
                  random = ~1|Site, data = npp.seas, 
                  weights=varExp())
nd.nppseas <-lme(numd ~ year + Species/pr + Species/tem + Species/ev +  
                 Species/wt + Species/cov,
                   random = ~ 1|Site, weights=varExp(), npp.seas[!(is.na(npp.seas$numd)),])
npp.nppseas <-  lme(prod_trans_st ~ year + Species/pr + Species/tem + Species/ev + Species/pa + 
                    Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm + 
                      Species/numd,
                  random = ~year-1|Site, data = npp.seas[!(is.na(npp.seas$numd)),])
#MuMIn::r.squaredGLMM(nd.nppseas)

# save coefs from models
wt.tab <- data.frame(summary(wt.nppseas)$tTable[-1,])
wt.tab$Predictor <- rownames(wt.tab)
cov.tab <- data.frame(summary(cov.nppseas)$tTable[-1,])
cov.tab$Predictor <- rownames(cov.tab)
nm.tab <- data.frame(summary(nm.nppseas)$tTable[-1,])
nm.tab$Predictor <- rownames(nm.tab)
nd.tab <- data.frame(summary(nd.nppseas)$tTable[-1,])
nd.tab$Predictor <- rownames(nd.tab)
npp.tab <- data.frame(summary(npp.nppseas)$tTable[-1,])
npp.tab$Predictor <- rownames(npp.tab)
# put together data frame with coefs and change column names
sem.npp <- data.frame(Response = c(rep("wt",nrow(wt.tab)),rep("cov",nrow(cov.tab)),
                              rep("Nm",nrow(nm.tab)),rep("numd",nrow(nd.tab)),
                              rep("NPP",nrow(npp.tab))), 
                 rbind(wt.tab, cov.tab, nm.tab, nd.tab, npp.tab))
colnames(sem.npp)[2] <- "estimate"
colnames(sem.npp)[6] <- "P.Value"
sem.npp <- sem.npp[,c(1,7,2,3,4,5,6)]


# remove species and year effect and split up coefs on species
sem.npp.mag <- sem.npp[unique(c(grep("S.mag", sem.npp$Predictor), grep("S.fus", sem.npp$Predictor, invert = TRUE))),]
sem.npp.mag <- sem.npp.mag[!(sem.npp.mag$Response == "wt"),]

sem.npp.mag[2] <- lapply(sem.npp.mag[2], gsub, pattern = "SpeciesS.magellanicum:", 
                    replacement = "", fixed = TRUE)
sem.npp.mag <- sem.npp.mag[!(sem.npp.mag$Predictor == "SpeciesS.magellanicum") &
                   !(sem.npp.mag$Predictor == "yearyr2014"),]


sem.npp.fus <- sem.npp[unique(c(grep("S.fus", sem.npp$Predictor), grep("S.mag", sem.npp$Predictor, invert = TRUE))),]
sem.npp.fus <- sem.npp.fus[!(sem.npp.fus$Response == "wt"),]
sem.npp.fus[2] <- lapply(sem.npp.fus[2], gsub, pattern = "SpeciesS.fuscum:", 
                    replacement = "", fixed = TRUE)
sem.npp.fus <- sem.npp.fus[!(sem.npp.fus$Predictor == "yearyr2014"),]
write.xlsx(sem.npp, "sem.npp_april2.xlsx")

 # plot sem
#def.par <- par(no.readonly = TRUE) 
x11(width=6, height=12)
pdf("NPPseastest.pdf",width=6, height=12)
#Save as .eps:
#setEPS()
#postscript("NPPplot.eps",width=6, height=12)
par(mfrow=c(2,1), xpd = NA)
sem.plot(coef.table=sem.npp.mag,alpha = 0.1)
legend(-1.8,1.8,lty=c(1,2,1,1), col=c(1,1,1,2),
       legend=c("P<0.1","P>0.1","pos","neg"), xpd = TRUE, bty="n")
mtext("NPP season\\n", line=2)
mtext("S.mag\\n")
sem.plot(coef.table=sem.npp.fus,alpha = 0.1)
mtext("S.fuscum\\n")
dev.off()



#______Prec-Evap mixed model####
min.val <- min(resp[["prod"]]$prod)*-1+0.001 # get value to add for positive response
NPP_E_P<-lmer((prod+min.val) ~ Species+year + Species*ev_pre + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd +
             (year-1|Site), data = resp[["prod"]])
plot(NPP_E_P)
summary(l4<-powerTransform(NPP_E_P, family = "bcPower"))#bcnPower neccessary with neg. values 
resp[["prod"]]$prod_trans<-(resp[["prod"]]$prod+min.val)^(l4$lambda)

NPP_E_P2<-lmer(prod_trans ~ Species+year + Species*ev_pre + Species*tem + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd +
               (year-1|Site), data = resp[["prod"]])
plot(NPP_E_P2)

NPP_E_P3<-lmer(prod_trans ~ Species+year + Species/ev_pre + Species/tem + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +Species/numd +
                 (year-1|Site), data = resp[["prod"]])

summary(NPP_E_P3)
Anova(NPP_E_P2, test.statistic = "F")

#no interactions
NPP_E_P4<-lmer(prod_trans ~ Species + year + ev_pre + tem + pa + 
               wt + cov + nor + ndeL + Nm + numd+ (year-1|Site), 
             data = resp[["prod"]])

plot(NPP_E_P4)
summary(NPP_E_P4)
Anova(NPP_E_P4, test.statistic = "F")


#___Bulk density####
mod5<-lmer(BD ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
             (1|Site), data = resp[["BD"]])
plot(mod5)
summary(l5<-powerTransform(mod5, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.27
resp[["BD"]]$BD_trans<-resp[["BD"]]$BD^(l5$lambda)

mod5.2<-lmer(BD_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (1|Site), data = resp[["BD"]])
plot(mod5.2)

mod5.2.1<-lmer(BD_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
                 (1|Site), data = resp[["BD"]])

summary(mod5.2.1)
Anova(mod5.2, test.statistic = "F")
summary(mod5.2)


#no interactions
mod5.4<-lmer(BD_trans ~ Species + year + pr + tem + ev + pa + 
               wt + cov + nor + ndeL + Nm + (1|Site), 
             data = resp[["BD"]])

plot(mod5.4)
summary(mod5.4)
Anova(mod5.4, test.statistic = "F")

#BD season null model
mod5n<-lmer((BD_trans)~1 + (1|Site),
            data = resp[["BD"]])
plot(mod5n)
summary(mod5n)


#______Excluding outliers from model with no interactions####
outliersBDseas<-resp[["BD"]][which(residuals(mod5.2)<(-1.0)),]

data_long5<-resp[["BD"]][which(residuals(mod5.2)>(-1.0)),]

modx5.1<-lmer(BD_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
                (1|Site), 
              data = data_long5)

mod5.1<-lmer(BD_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
                (1|Site), 
              data = data_long5)

plot(mod5.1)
summary(modx5.1)
Anova(mod5.1, test.statistic = "F")

#no inrteractions
modx5<-lmer(BD_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + ndeL + Nm +
              (1|Site), 
            data = data_long5, na.action = na.omit)
plot(modx5)
summary(modx5)
Anova(modx5, test.statistic = "F")


#MIXED MODEL with prec.-evap.
mod5.5<-lmer(BD_trans ~ Species*year + Species*ev_pre  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm+
               (1|Site),
             data = resp[["BD"]])
plot(mod5.5)

mod5.5x<-lmer(BD_trans ~ Species/year + Species/ev_pre  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm+
                (1|Site),
              data = resp[["BD"]])

summary(mod5.5)
Anova(mod5.5, test.statistic = "F")
summary(mod5.5x)

#Run without interactions to get coefficients
mod5.6<-lmer(BD_trans ~ Species + year + ev_pre  + tem + pa_d + wt + cov + nor+ ndeL + Nm+ 
               (1|Site),
             data = resp[["BD"]])
plot(mod5.6)
summary(mod5.6)
Anova(mod5.6, test.statistic = "F")
#data_long.seas2.na<- mod3.6@frame # save data w/o NAs

#___Numerical density####
mod6<-lmer(numden ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
             (1|Site), data = resp[["numden"]])
plot(mod6)
summary(l6<-powerTransform(mod6, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.1366
resp[["numden"]]$den_trans<-(resp[["numden"]]$numden*100)^(l6$lambda)

mod6.2<-lmer(den_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (1|Site), data = resp[["numden"]])
plot(mod6.2)

mod6.2.1<-lmer(den_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
                 (1|Site), resp[["numden"]])

summary(mod6.2.1)
Anova(mod6.2, test.statistic = "F")
summary(mod6.2)


#no interactions
mod6.4<-lmer(den_trans ~ Species + year + pr + tem + ev + pa + 
               wt + cov + nor + ndeL + Nm + (1|Site), 
             data = resp[["numden"]])

plot(mod6.4)
summary(mod6.4)
Anova(mod6.4, test.statistic = "F")

#den season null model
mod6n<-lmer((den_trans)~1 + (1|Site),
            data = resp[["numden"]])
plot(mod6n)
summary(mod6n)


#______Excluding outliers from model with no interactions####
outliersdenseas<-resp[["numden"]][which(residuals(mod6.2)<(-0.55)),]

data_long6<-data_long.na6[which(residuals(mod6.2)>(-0.55)),]

modx6.1<-lmer(den_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
                (1|Site), 
              data = data_long6)

mod6.1<-lmer(den_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (1|Site), 
             data = data_long6)

plot(mod6.1)
summary(modx6.1)
Anova(mod6.1, test.statistic = "F")

#no inrteractions
modx6<-lmer(den_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + ndeL + Nm +
              (1|Site), 
            data = data_long6)
plot(modx6)
summary(modx6)
Anova(modx6, test.statistic = "F")##end outlier test


#MIXED MODEL with prec.-evap.
mod6.5<-lmer(den_trans ~ Species*year + Species*ev_pre  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm+
               (1|Site),
             data = resp[["numden"]])
plot(mod6.5)

mod6.5x<-lmer(den_trans ~ Species/year + Species/ev_pre  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm+
                (1|Site),
              data = resp[["numden"]])

summary(mod6.5)
Anova(mod6.5, test.statistic = "F")
summary(mod6.5x)

#Run without interactions to get coefficients
mod6.6<-lmer(den_trans ~ Species + year + ev_pre  + tem + pa_d + wt + cov + nor+ ndeL + Nm+ 
               (1|Site),
             data = resp[["numden"]])
plot(mod6.6)
summary(mod6.6)
Anova(mod6.6, test.statistic = "F")

# Variances, R2 tables ####

# get R2 values for (1) between, and within site, (2) marginal (total variation explained by fixed effects)
# and conditional (total variation fixed + random effects)

#___LI season####
m1 <- lmer(LI_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm + Species*numd+ 
             (year-1|Site), data = resp[["LI"]])
m2 <- lmer(LI_trans ~ 1 +
             (year-1|Site), data = resp[["LI"]])
r2.fnc(mod3.2par, m2)
mod3.2par
r2.fnc(LI_P_E, m2)
#__Main effects only
m1<-lmer(LI_trans ~ Species+year+pr  + tem +ev +pa + 
           wt + cov + nor + ndeL + Nm+numd+
           (year-1|Site), data = resp[["LI"]])
m2 <- lmer(LI_trans ~ 1 +
             (year-1|Site), data = resp[["LI"]])
r2.fnc(m1, m2)
r2.fnc(LI_P_E3, m2)
# only PAR model
m1 <-lmer(LI_trans ~ Species+year + Species*pa +
                  (year|Site), data = resp[["LI"]])
m2 <-lmer(LI_trans ~ Species+year + 
            (year|Site), data = resp[["LI"]])
r2.fnc(m1, m2)

#___LI day####
m1<-lmer(LId_trans ~ Species+year+Species*pr_d + Species*tem + Species*ev_d + Species*pa_d + 
               Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm + Species*numd +
               (year-1|Site), data = resp[["LI_d"]])
m2 <- lmer(LId_trans ~ 1 +
             (year-1|Site), data = resp[["LI_d"]])
r2.fnc(m1, m2)
r2.fnc(mod2_evpre, m2)
#__Main effects only
m1<-lmer(LId_trans ~ Species+year+pr_d  + tem +ev_d +pa_d + 
          wt + cov + nor + ndeL + Nm + numd +
           (year-1|Site), data = resp[["LI_d"]])
m2 <- lmer(LId_trans ~ 1 +
             (year-1|Site), data = resp[["LI_d"]])
r2.fnc(m1, m2)
r2.fnc(mod3_evpre, m2)
#___NPP season####
m1 <-lmer(prod_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
            Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +Species*numd +
            (year-1|Site), data = resp[["prod"]])
m2 <-lmer(prod_trans ~ 1 +
            (year-1|Site), data = resp[["prod"]])
r2.fnc(m1, m2)
r2.fnc(NPP_E_P2, m2)
#__Main effects only
m1<-lmer(prod_trans ~ Species+year + pr + tem + ev + pa + 
           wt + cov + nor + ndeL + Nm +numd+
           (year-1|Site), data = resp[["prod"]])
m2 <- lmer(prod_trans ~ 1 +
             (year-1|Site), data = resp[["prod"]])
r2.fnc(m1, m2)
r2.fnc(NPP_E_P4, m2)
#___NPP day####
m1 <-lmer(prod.d_trans ~ Species + year + Species*pr_d + Species*tem + Species*ev_d + 
            Species*pa_d + Species*wt + Species*cov + 
            Species*nor + Species*ndeL + Species*Nm + Species*numd +
            (year-1|Site), 
          data= resp[["prod_d"]])
m2 <-lmer(prod.d_trans ~ 1 +
            (year-1|Site), data = resp[["prod_d"]])
r2.fnc(m1, m2)
r2.fnc(modE_P2, m2)

#__Main effects only
m1 <-lmer(prod.d_trans ~ Species +year + pr_d + tem + ev_d + 
            pa_d + wt + cov + 
            nor + ndeL + Nm + numd +
            (year-1|Site), 
          data= resp[["prod_d"]])
m2 <-lmer(prod.d_trans ~ 1 +
            (year-1|Site), data = resp[["prod_d"]])
r2.fnc(m1, m2)
r2.fnc(modE_P5, m2)
#__R2 CIs####
# run resampling R2 function to get CIs of R2 
# set sims=1000, this takes a while though!

ee <- r2.resamp(mod2.2, # m1 is the model of interest
                sims = 1000,  # number of simulations
                only.whole=FALSE) # TRUE=only whole model CIs, if FALSE also each variable(SLOW!!!)


# get Table with CI and median R2 values
ee.vec <- sapply(ee, function (x) apply(x,1:2, median))
#ee.vec[is.na(ee.vec)] <- 0
ee.ci <- apply(ee.vec,1, function (x) paste(round(quantile(x, c(0.05,0.50,0.95), na.rm =T), digits=3), 
                                            collapse = ",", sep = "") )
out = data.frame(row.names = as.character(rownames(ee[[1]])), 
                 b.site = ee.ci[1:nrow(ee[[1]])], 
                 w.site = ee.ci[(nrow(ee[[1]])+1):(2*nrow(ee[[1]]))], 
                 marginal = ee.ci[(2*nrow(ee[[1]])+1):(3*nrow(ee[[1]]))],
                 conditional = ee.ci[(3*nrow(ee[[1]])+1):(4*nrow(ee[[1]]))])
out



# Plot response ####
library(interactions)
library(jtools)

#__NPP####
toplot <- plot.model.fnc(mod=mod4.2,           # model name          
           back = TRUE,      # if you want to backtransform the response (FALSE if not) 
           x.scale.orig = TRUE,   # should x-axis (predictor) be on the original scale?
           orig.data = resp[["prod"]],  # orignal data frame used in fitting the model resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),] or resp[["LI"]]
           lam = l4$lambda,       # give transformation lambda value (lam=) for backtransformation, 
           # does not account for the *100 for daily values
           shift = min.val,  #min.val value given to get positive responses. If not set as 0.
           mult = 1,       # if response was multiplied for larger values. If not set as 1.
           ylab = element_blank(), #expression(NPP ~ (g ~m^{-2} ~yr^{-1})), # set y-axis label 
           vars.plot = c("Species", "year",  # you always need these # IF all variables in model, set vars=NULL 
                        "pr","tem", "ev", "wt", "ndeL", "pa","nor","cov", "Nm","numd"))

png("npp.png", width=22, height=26, units="cm", res=300)
grid.arrange(grobs=toplot, ncol=2) # plot the plots
dev.off()

#__LI####
toplot <- plot.model.fnc(mod=mod3.2,           # model name          
                         back = TRUE,      # if you want to backtransform the response (FALSE if not) 
                         x.scale.orig = TRUE,   # should x-axis (predictor) be on the original scale?
                         orig.data = resp[["LI"]],  # orignal data frame used in fitting the model resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),] or resp[["LI"]]
                         lam = l3$lambda,       # give transformation lambda value (lam=) for backtransformation, 
                         # does not account for the *100 for daily values
                         shift = min.val,  #min.val value given to get positive responses. If not set as 0.
                         mult = 1,       # if response was multiplied for larger values. If not set as 1.
                         ylab = element_blank(), #ylab = expression(LI ~ (mm ~yr^{-1})), # set y-axis label LI,  LI ~ (mm ~yr^{-1})
                         vars.plot = c("Species", "year",  # you always need these # IF all variables in model, set vars=NULL 
                                       "pr", "tem", "ev", "wt", "ndeL", "pa","nor","cov", "Nm","numd")) # LI

png("li.png", width=22, height=26, units="cm", res=300)
grid.arrange(grobs=toplot, ncol=2) # plot the plots
dev.off()

# Plot Biome and sites####
# Need to add worldclim variables to data file
#The biome plot using worldclim data for our sites
#names(dat3)
#tem per month 167:178
#prec per month 179:190
#calculate rowMeans, so av. per year
temp_year<-rowMeans(subset(dat, select=167:178),na.rm = TRUE)
prec_year<-(rowSums(subset(dat, select=179:190),na.rm = TRUE))/10
temp_year_site <- aggregate(temp_year ~ site.verified, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)
prec_year_site<- aggregate(prec_year ~ site.verified, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)
year_site<-merge(prec_year_site, temp_year_site)
#plot(temp_year_site$temp_year,prec_year_site$prec_year/10, ylim=c(16,440), xlim=c(-17, 33))


#install.packages("devtools")
# you need the 'plotbiomes' package from github
#devtools::install_github("valentinitnelav/plotbiomes")
library(ggplot2)
library(plotbiomes)
whittaker_base_plot() +
  geom_point(data = year_site,
             aes(x = temp_year,
                 y = prec_year),
             size   = 3,
             shape  = 21,
             stroke = 1,
             alpha  = 0.5) +
  theme_bw()

#Site level data and tables ####
# aggregate data frame by returning means for numeric variables
# means based on site and species

# only include bulk density data from patches where there are NPP data as well.
dat$BD13 <- ifelse(is.na(dat$Prod13), NA,  dat$BD13)
dat$BD14 <- ifelse(is.na(dat$Prod14), NA,  dat$BD14)


Avs <- aggregate(cbind(Beging_season_vasc_2013, End_season_vasc_2013, HWT_begin_season_2013, HWT_end_season_2013, Beging_season_vasc_2014,End_season_vasc_2014,HWT_begin_season_2014,HWT_end_season_2014,
                       Patch_coord_lat, Patch_coord_lon, days2013, days2014, LI13, LI14, Prod13, Prod14, LI_per_day2013, LI_per_day2014, prod_per_day2013, prod_per_day2014, BD13, BD14,
                       cover13, cover14, HWT13, HWT14, N_per, C_per, CN, P_per, ndep_Lam13,
                       temp13, temp14, evap13, evap14, ev_d13, ev_d14, precip13, precip14, prec_d13, prec_d14, norain13.mean, norain14.mean, norain13.max, norain14.max, par13, par14, par_d13, par_d14,
                       num.den13, num.den14,Patch_coord_lat, googlemap_elevation_masl) ~ Site + Species, dat,na.rm=TRUE, na.action="na.pass", FUN=mean)


#Average per site
Avs_site <- aggregate(cbind(Beging_season_vasc_2013, End_season_vasc_2013, HWT_begin_season_2013, HWT_end_season_2013, 
                            Beging_season_vasc_2014,End_season_vasc_2014,
                            HWT_begin_season_2014,HWT_end_season_2014,
                            cover13, cover14, HWT13, HWT14,Patch_coord_lat, Patch_coord_lon, 
                            ndep_Lam13, temp13, temp14, num.den13, num.den14,
                            evap13, evap14, precip13, precip14, norain13.mean, norain14.mean, par13, par14,
                            googlemap_elevation_masl, LI13, LI14, Prod13, Prod14, ) ~ 
                      Site, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)


#Average per site, no matter which species (for table S1,S2)
Avs_site2 <- aggregate(cbind(ndep_Lam13, googlemap_elevation_masl) ~ site.verified, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)
#library(openxlsx)
#write.xlsx(Avs_site2,"Avs_site.xlsx") 

# Save site averages based on data long
av.long<-aggregate(cbind(LI_d, prod_d, LI, prod, prec_d, temp, evap_d, par_d,
                         hwt, cover, norain, ndep_Lamarque2013, Nmean) ~ 
                     Site, data_long, na.rm=TRUE,na.action="na.pass", FUN=mean)

# Correlations####
#_ndep and n tissue and ndep and vasc. plant cover####
cor.test(av.long$ndep_Lamarque2013,av.long$Nmean)#0.35, df 97
cor.test(av.long$ndep_Lamarque2013,av.long$cover)#0.25, df 95

Avs_site$mean_cover <- apply(Avs_site[, c("cover13","cover14")],1, mean, na.rm=T)#average cover between the years
cor.test(Avs_site$ndep_Lam13, Avs_site$mean_cover)#0.23, p=0.045, df=74
cor.test(Avs_site$ndep_Lam13, Avs_site$Nmean)

#_LI and NPP correlated?####
with(data_long[data_long$Species=="S.fuscum",], cor.test(prod, LI))#0.59, df 574
with(data_long[data_long$Species=="S.magellanicum",], cor.test(prod, LI))#0.57, df 598

#_Beginning and end season####
cor.test(data_long$HWT_begin_season_2013,data_long$HWT_end_season_2013)#cor=0.75, df=1106 
cor.test(data_long$HWT_begin_season_2014,data_long$HWT_end_season_2014)#=0.70, df 1132
cor.test(data_long$Beging_season_vasc_2013,data_long$End_season_vasc_2013)#0.92, df 1036
cor.test(data_long$Beging_season_vasc_2014,data_long$End_season_vasc_2014)#0.86, df 896

#On site level.
cor.test(Avs_site$HWT_begin_season_2013, Avs_site$HWT_end_season_2013)#r=0.79, df 83
cor.test(Avs_site$HWT_begin_season_2014, Avs_site$HWT_end_season_2014)#r=0.79, df 85
cor.test(Avs_site$Beging_season_vasc_2013, Avs_site$End_season_vasc_2013)#r=0.92, df 76
cor.test(Avs_site$Beging_season_vasc_2014, Avs_site$End_season_vasc_2014)#r=0.88, df 64

#_Within species correlations####
Av.fu<-Avs[which(Avs$Species== "S.fuscum"),]#S. fuscum df
levels(as.factor(Av.fu$Site))#85 Sites
cor.test(Av.fu$HWT13, Av.fu$HWT14)#r=0,74, df 67
cor.test(Av.fu$cover13, Av.fu$cover14)#r=0,89, df 61
cor.test(Av.fu$LI13, Av.fu$LI14)#r=0,68, df=68, p<0.0001
cor.test(Av.fu$Prod13, Av.fu$Prod14)#r=0,58, df=65, p<0.0001
cor.test(Av.fu$LI13, Av.fu$Prod13)#r=0,51, df=74, p<0.0001
cor.test(Av.fu$LI14, Av.fu$Prod14)#r=0,55, df=73, p<0.0001


Av.mg<-Avs[which(Avs$Species=='S.magellanicum'),]#S. magellanicum df
levels(sitmg<-as.factor(Av.mg$Site))#91
cor.test(Av.mg$HWT13, Av.mg$HWT14)#r=0,78, df 78
cor.test(Av.mg$cover13, Av.mg$cover14)#r=0,85, df 68
cor.test(Av.mg$LI13, Av.mg$LI14)#r=0,69, df=77, p<0.0001
cor.test(Av.mg$Prod13, Av.mg$Prod14)#r=0,48, df=73, p<0.0001
cor.test(Av.mg$LI13, Av.mg$Prod13)#r=0,56, df=79, p<0.0001
cor.test(Av.mg$LI14, Av.mg$Prod14)#r=0,57, df=82, p<0.0001


#_corr plots####
library(corrplot)
M <- cor(data.frame(Avs[c("days2013", "days2014", "LI13", "LI14", "Prod13", "Prod14", "LI_per_day2013", "LI_per_day2014",   
                                "prod_per_day2013", "prod_per_day2014", "BD13", "BD14", "cover13",  "cover14", "HWT13", "HWT14",  "N_per",  "C_per", "CN",
                                "P_per",  "ndep_Lam13", "temp13", "temp14", "evap13", "evap14", "ev_d13", "ev_d14", "precip13", "precip14", "prec_d13", "prec_d14", "norain13.mean",   
                                 "norain14.mean", "norain13.max", "norain14.max", "par13", "par14", "par_d13", "par_d14", "Patch_coord_lat")] ), use = "complete.obs")
par(mfrow=c(1,1))
corrplot(M, method = "number", number.cex = 0.7)

#Split per year and day/season to see better
M_per_d13 <- cor(data.frame(Avs[c("LI_per_day2013", "prod_per_day2013", "BD13", "cover13",  "HWT13", "N_per", "C_per", "CN",
                          "P_per",  "ndep_Lam13", "temp13", "ev_d13", "prec_d13", "norain13.mean",   
                           "norain13.max", "par_d13", "Patch_coord_lat")] ), use = "complete.obs")
M_per_d14 <- cor(data.frame(Avs[c( "LI_per_day2014", "prod_per_day2014", "BD14",   "cover14",  "HWT14",  "N_per",  "C_per", "CN",
                                  "P_per",  "ndep_Lam13",  "temp14",  "ev_d14",  "prec_d14", "norain14.mean",  "norain14.max", 
                                  "par_d14", "Patch_coord_lat")] ), use = "complete.obs")
M_season13 <- cor(data.frame(Avs[c("LI13", "Prod13",  "BD13", "cover13", "HWT13", "N_per",  "C_per", "CN",
                          "P_per",  "ndep_Lam13", "temp13", "evap13", "precip13", "norain13.mean",   
                           "norain13.max",  "par13", "Patch_coord_lat")] ), use = "complete.obs")
M_season14 <- cor(data.frame(Avs[c("days2014", "LI14", "Prod14", "BD14",  "cover14", "HWT14",  "N_per",  "C_per", "CN",
                                   "P_per",  "ndep_Lam13",  "temp14", "evap14", "precip14", "norain14.mean", "norain14.max", "par14", "Patch_coord_lat")] ), use = "complete.obs")

par(mfrow=c(1,1))
corrplot(M, method = "number", number.cex = 0.7)
corrplot(M_per_d13, method = "number", number.cex = 0.7)
corrplot(M_per_d14, method = "number", number.cex = 0.7)
corrplot(M_season13, method = "number", number.cex = 0.7)
corrplot(M_season14, method = "number", number.cex = 0.7)

#Pairs for each year
#2013
pairs(data.frame(Avs2[c("LI_per_day2013", "prod_per_day2013", "BD13",
                       "BD14", "temp13", "evap13",  "precip13",
                       "norain13.mean", "norain13.max", "par13", "ndep_Lam13", "Beging_season_vasc_2013","End_season_vasc_2013",
                       "HWT_begin_season_2013",  "HWT_end_season_2013", "CN", "P_per", "Patch_coord_lat")] ), main= "2013 excluding China & Japan")
#2014
pairs(data.frame(Avs[c("LI_per_day2014", "prod_per_day2014", 
                       "BD14",  "temp14",  "evap14",  "precip14",
                        "norain14.mean", "norain14.max", "par14", "ndep_Lam13", 
                       "HWT_begin_season_2014", "HWT_end_season_2014", "CN", "P_per", "Patch_coord_lat")] ),main= "2014 excluding China & Japan")


#Figure 4#######
# Save site averages based on data long
av_long<-aggregate(cbind(LI_d, prod_d, LI, prod, prec_d, temp, evap_d, par_d,
                         hwt, cover, norain, ndep_Lamarque2013, Nmean, numden, BD) ~ 
                     Species + Site +year, data_long, na.rm=TRUE,na.action="na.pass", FUN=mean)

library(ggplot2)

av_long$YEAR<-factor(av_long$year,levels = c("yr2013","yr2014"))
av_long$SPECIES<-factor(av_long$Species,levels = c("S.fuscum","S.magellanicum"))
av_long$sp.year<-paste(av_long$SPECIES,av_long$YEAR,sep="_")
av_long$sp.year<-factor(av_long$sp.year,levels = c("S.fuscum_yr2013","S.fuscum_yr2014","S.magellanicum_yr2013","S.magellanicum_yr2014"))

p1 <- ggplot(data=av_long,aes(x=sp.year,y=LI, fill=YEAR))+
  scale_fill_grey(start=0.8, end=0.4)+
  ylab("LI"~ (mm ~yr^{-1}))+
  xlab(element_blank())+
  scale_shape_identity()+
  geom_boxplot()+
  stat_summary(colour="white")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     axis.title=element_text(size=20), axis.text=element_text(size=16),axis.text.x = element_blank(), axis.ticks = element_blank(),
                     legend.title = element_text(size = rel(1.5)),legend.text = element_text(size=rel(1.3)),legend.key.size =  unit(0.3, "in"))
p2 <- ggplot(data=av_long,aes(x=sp.year,y=BD,fill=YEAR))+
  scale_fill_grey(start=0.8, end=0.4)+
  ylab("BD"~ (kg ~m^{-3}))+
  xlab(element_blank())+
  scale_shape_identity()+
  geom_boxplot()+
  stat_summary(colour="white")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     axis.title=element_text(size=20), axis.text=element_text(size=16),axis.text.x = element_blank(), axis.ticks = element_blank(),
                     legend.title = element_text(size = rel(1.5)),legend.text = element_text(size=rel(1.3)),legend.key.size =  unit(0.3, "in"))

p3 <- ggplot(data=av_long,aes(x=sp.year,y=prod,fill=YEAR))+
  scale_fill_grey(start=0.8, end=0.4)+
  ylab("NPP"~ (g ~m^{-2} ~yr^{-1}))+
  xlab(element_blank())+
  scale_shape_identity()+
  geom_boxplot()+
  stat_summary(colour="white")+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     axis.title=element_text(size=20), axis.text=element_text(size=16),axis.text.x = element_blank(), axis.ticks = element_blank(),
                     legend.title = element_text(size = rel(1.5)),legend.text = element_text(size=rel(1.3)),legend.key.size =  unit(0.3, "in"))


png("box.png", width=15, height=27, units="cm", res=600)
grid.arrange(grobs=list(p1, p2, p3), ncol=1) # plot the plots
dev.off()


#Tables with descriptive stats for Supplementary material, from averaged per site df####
std<-function(x) sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
names(Avs)

numden13.av<-tapply(Avs$num.den13,list(Avs$Species),mean,na.rm=T)
numden13.SE<-tapply(Avs$num.den13,list(Avs$Species),std)
numden13.min<-tapply(Avs$num.den13,list(Avs$Species),min, na.rm=T)
numden13.max<-tapply(Avs$num.den13,list(Avs$Species),max, na.rm=T)

LI13average<-tapply(Avs$LI13,list(Avs$Species),mean,na.rm=T)
LI13SE<-tapply(Avs$LI13,list(Avs$Species),std)
LI13min<-tapply(Avs$LI13,list(Avs$Species),min, na.rm=T)
LI13max<-tapply(Avs$LI13,list(Avs$Species),max, na.rm=T)
NPP13average<-tapply(Avs$Prod13,list(Avs$Species),mean,na.rm=T)
NPP13SE<-tapply(Avs$Prod13,list(Avs$Species),std)
NPP13min<-tapply(Avs$Prod13,list(Avs$Species),min, na.rm=T)
NPP13max<-tapply(Avs$Prod13,list(Avs$Species),max, na.rm=T)
Prec13average<-tapply(Avs$precip13,list(Avs$Species),mean,na.rm=T)
Prec13SE<-tapply(Avs$precip13,list(Avs$Species),std)
Prec13min<-tapply(Avs$precip13,list(Avs$Species),min, na.rm=T)
Prec13max<-tapply(Avs$precip13,list(Avs$Species),max, na.rm=T)
Temp13average<-tapply(Avs$temp13,list(Avs$Species),mean,na.rm=T)
Temp13SE<-tapply(Avs$temp13,list(Avs$Species),std)
Temp13min<-tapply(Avs$temp13,list(Avs$Species),min, na.rm=T)
Temp13max<-tapply(Avs$temp13,list(Avs$Species),max, na.rm=T)
Evap13average<-tapply(Avs$evap13,list(Avs$Species),mean,na.rm=T)
Evap13SE<-tapply(Avs$evap13,list(Avs$Species),std)
Evap13min<-tapply(Avs$evap13,list(Avs$Species),min, na.rm=T)
Evap13max<-tapply(Avs$evap13,list(Avs$Species),max, na.rm=T)
PAR13average<-tapply(Avs$par13,list(Avs$Species),mean,na.rm=T)
PAR13SE<-tapply(Avs$par13,list(Avs$Species),std)
PAR13min<-tapply(Avs$par13,list(Avs$Species),min, na.rm=T)
PAR13max<-tapply(Avs$par13,list(Avs$Species),max, na.rm=T)
HWT13average<-tapply(Avs$HWT13,list(Avs$Species),mean,na.rm=T)
HWT13SE<-tapply(Avs$HWT13,list(Avs$Species),std)
HWT13min<-tapply(Avs$HWT13,list(Avs$Species),min, na.rm=T)
HWT13max<-tapply(Avs$HWT13,list(Avs$Species),max, na.rm=T)
Cover13average<-tapply(Avs$cover13,list(Avs$Species),mean,na.rm=T)
Cover13SE<-tapply(Avs$cover13,list(Avs$Species),std)
Cover13min<-tapply(Avs$cover13,list(Avs$Species),min, na.rm=T)
Cover13max<-tapply(Avs$cover13,list(Avs$Species),max, na.rm=T)
Nor13average<-tapply(Avs$norain13.mean,list(Avs$Species),mean,na.rm=T)
Nor13SE<-tapply(Avs$norain13.mean,list(Avs$Species),std)
Nor13min<-tapply(Avs$norain13.mean,list(Avs$Species),min, na.rm=T)
Nor13max<-tapply(Avs$norain13.mean,list(Avs$Species),max, na.rm=T)

numden14.av<-tapply(Avs$num.den14,list(Avs$Species),mean,na.rm=T)
numden14.SE<-tapply(Avs$num.den14,list(Avs$Species),std)
numden14.min<-tapply(Avs$num.den14,list(Avs$Species),min, na.rm=T)
numden14.max<-tapply(Avs$num.den14,list(Avs$Species),max, na.rm=T)

LI14average<-tapply(Avs$LI14,list(Avs$Species),mean,na.rm=T)
LI14SE<-tapply(Avs$LI14,list(Avs$Species),std)
LI14min<-tapply(Avs$LI14,list(Avs$Species),min, na.rm=T)
LI14max<-tapply(Avs$LI14,list(Avs$Species),max, na.rm=T)
NPP14average<-tapply(Avs$Prod14,list(Avs$Species),mean,na.rm=T)
NPP14SE<-tapply(Avs$Prod14,list(Avs$Species),std)
NPP14min<-tapply(Avs$Prod14,list(Avs$Species),min, na.rm=T)
NPP14max<-tapply(Avs$Prod14,list(Avs$Species),max, na.rm=T)
Prec14average<-tapply(Avs$precip14,list(Avs$Species),mean,na.rm=T)
Prec14SE<-tapply(Avs$precip14,list(Avs$Species),std)
Prec14min<-tapply(Avs$precip14,list(Avs$Species),min, na.rm=T)
Prec14max<-tapply(Avs$precip14,list(Avs$Species),max, na.rm=T)
Temp14average<-tapply(Avs$temp14,list(Avs$Species),mean,na.rm=T)
Temp14SE<-tapply(Avs$temp14,list(Avs$Species),std)
Temp14min<-tapply(Avs$temp14,list(Avs$Species),min, na.rm=T)
Temp14max<-tapply(Avs$temp14,list(Avs$Species),max, na.rm=T)
Evap14average<-tapply(Avs$evap14,list(Avs$Species),mean,na.rm=T)
Evap14SE<-tapply(Avs$evap14,list(Avs$Species),std)
Evap14min<-tapply(Avs$evap14,list(Avs$Species),min, na.rm=T)
Evap14max<-tapply(Avs$evap14,list(Avs$Species),max, na.rm=T)
PAR14average<-tapply(Avs$par14,list(Avs$Species),mean,na.rm=T)
PAR14SE<-tapply(Avs$par14,list(Avs$Species),std)
PAR14min<-tapply(Avs$par14,list(Avs$Species),min, na.rm=T)
PAR14max<-tapply(Avs$par14,list(Avs$Species),max, na.rm=T)
HWT14average<-tapply(Avs$HWT14,list(Avs$Species),mean,na.rm=T)
HWT14SE<-tapply(Avs$HWT14,list(Avs$Species),std)
HWT14min<-tapply(Avs$HWT14,list(Avs$Species),min, na.rm=T)
HWT14max<-tapply(Avs$HWT14,list(Avs$Species),max, na.rm=T)
Cover14average<-tapply(Avs$cover14,list(Avs$Species),mean,na.rm=T)
Cover14SE<-tapply(Avs$cover14,list(Avs$Species),std)
Cover14min<-tapply(Avs$cover14,list(Avs$Species),min, na.rm=T)
Cover14max<-tapply(Avs$cover14,list(Avs$Species),max, na.rm=T)
Nor14average<-tapply(Avs$norain14.mean,list(Avs$Species),mean,na.rm=T)
Nor14SE<-tapply(Avs$norain14.mean,list(Avs$Species),std)
Nor14min<-tapply(Avs$norain14.mean,list(Avs$Species),min, na.rm=T)
Nor14max<-tapply(Avs$norain14.mean,list(Avs$Species),max, na.rm=T)

Ndepaverage<-tapply(Avs$ndep_Lam13,list(Avs$Species),mean,na.rm=T)
NdepSE<-tapply(Avs$ndep_Lam13,list(Avs$Species),std)
Ndepmin<-tapply(Avs$ndep_Lam13,list(Avs$Species),min, na.rm=T)
Ndepmax<-tapply(Avs$ndep_Lam13,list(Avs$Species),max, na.rm=T)
Ntissaverage<-tapply(Avs$N_per,list(Avs$Species),mean,na.rm=T)
NtissSE<-tapply(Avs$N_per,list(Avs$Species),std)
Ntissmin<-tapply(Avs$N_per,list(Avs$Species),min, na.rm=T)
Ntissmax<-tapply(Avs$N_per,list(Avs$Species),max, na.rm=T)
Eleaverage<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),mean,na.rm=T)
EleSE<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),std)
Elemin<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),min, na.rm=T)
Elemax<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),max, na.rm=T)

tLI<-cbind(LI13average, LI13SE, LI13min, LI13max, LI14average, LI14SE, LI14min, LI14max)
tNPP<-cbind(NPP13average, NPP13SE, NPP13min, NPP13max, NPP14average, NPP14SE, NPP14min, NPP14max)
tPrec<-cbind(Prec13average, Prec13SE, Prec13min, Prec13max, Prec14average, Prec14SE, Prec14min, Prec14max)
tTemp<-cbind(Temp13average, Temp13SE, Temp13min, Temp13max, Temp14average, Temp14SE, Temp14min, Temp14max)
tEvap<-cbind(Evap13average, Evap13SE, Evap13min, Evap13max, Evap14average, Evap14SE, Evap14min, Evap14max)
tPAR<-cbind(PAR13average, PAR13SE, PAR13min, PAR13max, PAR14average, PAR14SE, PAR14min, PAR14max)
tHWT<-cbind(HWT13average, HWT13SE, HWT13min, HWT13max, HWT14average, HWT14SE, HWT14min, HWT14max)
tCover<-cbind(Cover13average, Cover13SE, Cover13min, Cover13max, Cover14average, Cover14SE, Cover14min, Cover14max)
tNOR<-cbind(Nor13average, Nor13SE, Nor13min, Nor13max, Nor14average, Nor14SE, Nor14min, Nor14max)
tNdep<-cbind(Ndepaverage, NdepSE, Ndepmin, Ndepmax)
tNtiss<-cbind(Ntissaverage, NtissSE, Ntissmin, Ntissmax)
tElevation<-cbind(Eleaverage, EleSE, Elemin, Elemax)

desc.table<-cbind(tLI, tNPP, tPrec, tTemp, tEvap, tPAR, tHWT, tCover, tNOR, tNdep, tNtiss, tElevation)
#library(openxlsx)
#write.xlsx(desc.table,"Descriptive.xlsx") 
#desc.table<-cbind(tNdep, tNtiss, tElevation)

#For how many sites do we have data?####
#Cover
length(which(Avs_site$cover13>0))# & Avs_site$cover13>0)) 90 2013,  2014, 83 
length(which(Avs_site$cover13>0 | Avs_site$cover14 >0)) #total 76 with data for both years, either 97

#HWT
length(which(Avs_site$HWT13>0))# 94 2013, 90 2014
length(which(Avs_site$HWT13>0 & Avs_site$HWT14>0)) #total 85 with data for both years, 97 either

#Data for both spring and fall?
#Cover 
length(which(Avs_site$Beging_season_vasc_2013>0 & Avs_site$End_season_vasc_2013>0)) #total 78 with data for both spring and fall in 2013
91-78#=13 sites with data that only has data for one period 
length(which(Avs_site$End_season_vasc_2013>0))#87 begining, 82 end season

length(which(Avs_site$Beging_season_vasc_2014>0 & Avs_site$End_season_vasc_2014>0)) #total 66 with data for both spring and fall in 2013
83-66#= 17 sites with data for only one period 
length(which(Avs_site$Beging_season_vasc_2014>0))#76 begining, 73 end season

#HWT
length(which(Avs_site$HWT_begin_season_2013>0 & Avs_site$HWT_end_season_2013>0)) #total 85 with data for both spring and fall in 2013
95-85#=10 sites that only has data for one period 
length(which(Avs_site$HWT_begin_season_2013>0))#91 begining, 89 end season

length(which(Avs_site$HWT_begin_season_2014>0 & Avs_site$HWT_end_season_2014>0)) #total 87 with data for both spring and fall in 2013
90-87#= 3 sites with data for only one period 
length(which(Avs_site$HWT_end_season_2014>0))#89 begining, 88 end season

#Overall Avs_site : averages between species per site calculated 
mean(Avs_site$LI13, na.rm=T)#16.95
mean(Avs_site$LI14, na.rm=T)#18.12
mean(Avs_site$Prod13, na.rm=T)#194.70
mean(Avs_site$Prod14, na.rm=T)#182.10
std<-function(x) sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
std(Avs_site$LI13)#1.09
std(Avs_site$LI14)#1.19
std(Avs_site$Prod13)#12.13
std(Avs_site$Prod14)#11.12

# GGs BRMS stuff etc ####
mm0 <- (lme(LI_d ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
              Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm, random = ~ 1|Site,
            data = resp[["LI_d"]], control = list(maxIter=50)))
mm <- (lme(scale(LI_d*10) ~ Species+year+ Species*pr_d  + Species*tem + Species*ev_d + 
             Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm, random = ~ year-1|Site,
           data = resp[["LI_d"]], weights = varExp(form=~LI_d), 
           control = list(maxIter=50))) #
summary(mm)
summary(mm0)
anova(mm0, mm)
plot(mm)
plot(resid(mm, type="p") ~ resp[["prod"]]$Species)

library(brms)
bmgam <- brm(bf((LI_d+0.02826901)*100 ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + Species*pa_d + 
                  Species*wt + Species*cov + 
                  Species*nor + Species*ndeL + Species*Nm+
                  (year-1|Site)),family=lognormal(link="identity", link_sigma = "log"),,#brmsfamily("Gamma", link="identity", link_sigma = "log"),
             data = resp[["LI_d"]],
             iter = 3000, chains = 4, cores = 2)
summary(bmgam)
bb = bf( (LI_d+0.02826901)*10 ~ Species+year +Nm + (1|Site))
mod = brm(bb, #prior=bprior,
          data = resp[["LI_d"]], family=lognormal(link="identity", link_sigma = "log"),#family="gaussian", #link_sigma = "log", #family = 'gaussian',
          iter = 3000, chains = 2, cores = 2)


bb = bf(LI_d*10 ~ Species+year +Nm + (1|Site), 
        sigma ~  LI_d)
bb = bf(LI_d ~ Species+year +Nm + (1|Site)) + 
  nlf(sigma ~  exp(2*t*LI_d), 
      t ~ 1)
bprior <- prior(normal(0, 10), nlpar = t, lb = 0)
mod = brm(bb, #prior=bprior,
          data = resp[["LI_d"]], family=brmsfamily("gaussian", link_sigma = "log"),#family="gaussian", #link_sigma = "log", #family = 'gaussian',
          iter = 3000, chains = 2, cores = 2)
summary(mod)
plot(fitted(mod)[,1], residuals(mod, type="pearson",method = "fitted")[,1])

mo <- (lme(LI_d*10 ~ Species+year +Nm, random = ~ 1|Site,
           data = resp[["LI_d"]], control = list(maxIter=50), weights = varExp(form=~LI_d)))
summary(mo)
plot(mo)

mo.ga <- (glmer((LI_d+min.val) ~ Species+year +Nm + (1|Site),
                data = resp[["LI_d"]], family=Gamma(link = "log")))
plot(mo.ga)
summary(mo.ga)

mo.va <- (lme( (LI_d+0.02826901)*10 ~ Species+year +Nm, random = ~ 1|Site,
               data = resp[["LI_d"]], control = list(maxIter=50), weights = varExp()))
summary(mo.va)
mo.e <- (lme(log((LI_d+0.04826901)*10) ~ Species+year +Nm, random = ~ 1|Site,
             data = resp[["LI_d"]], control = list(maxIter=50)))
summary(mo.e)
plot(mo.e)
fitted_values <- fitted(mod)
## plot fitted means against actual response
dat <- as.data.frame(cbind(Y = standata(mod)$Y, fitted_values))
ggplot(dat) + geom_point(aes(x = Estimate, y = Y))


