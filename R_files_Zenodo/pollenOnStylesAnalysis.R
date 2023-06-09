# The goal of this script is to run 2 analyses that examine the relationship between pollen
# collected off of styles and the floral neighborhood predictors

# 1. How does the floral neighborhood influence conspecific pollen count on styles? 
# # numSpro ~ floral neighborhood predictors + observation day , data = ssffCon1
# # (could be overdispersed), 

# 2. How does the floral neighborhood influence heterospecific pollen count on styles?
# # hetAbun ~ floral neighborhood predictors + observation day , data = ssffCon1
# load packages                                                      ####
library(lme4)
library(tidyverse)
library(vegan)
library(DHARMa)
library(bbmle) 
library(MuMIn)
library(effects)
library(car)
library(glmmTMB)
set.seed(66)
# get datasets                                                       ####


# Starting Data Frames
ssHet <- read.csv('ssHet.csv', stringsAsFactors = F)
ssCon <- read.csv('ssCon.csv', stringsAsFactors = F)
ff <- read.csv('ff.csv', stringsAsFactors = F)
rr <- read.csv('rr.csv', stringsAsFactors = F)

# fixing ff
ff[ff$plantObsID == "aa-2-z",]
ff[ff$plantObsID == "aa-2-z","plantObsID"]
ff[ff$plantObsID == "aa-2-z","plantObsID"] <- "aa-2-e"
ff[ff$plantObsID == "aa-2-e",]

# Creating ssffCon
str(ssCon) # 337 obs. of  6 variables
str(ff) # 224 obs. of  53 variables

ssffCon <- merge(ssCon, ff, by = 'plantObsID', all = T)
str(ssffCon) # 347 obs. of  60 variables
ssffCon$plantObsID <- as.factor(ssffCon$plantObsID)
ssffCon$styleID <- as.factor(ssffCon$styleID)
ssffCon$fncDay <- as.factor(ssffCon$fncDay)
ssffCon$site <- as.factor(ssffCon$site)
ssffCon$tag <- as.factor(ssffCon$tag)

str(ssffCon)
head(ssffCon)
names(ssffCon) # preStyle = 1 is pre-obs style
ssffCon0 <- ssffCon[ssffCon$preStyle == 0,] # Here we are only interested in post-obs styles

str(ssffCon0) # 285 obs. of  60 variables
ssffCon1 <- na.omit(ssffCon0)
ssffCon1 <- ssffCon1 %>%
  mutate(fncDay = as.factor(ssffCon1$fncDay)) %>%
  mutate(infl.ecan = as.integer(infl.ecan + 1)) %>% # add one ech infl for the focal
  mutate(inflCt = rowSums(select(.,starts_with('infl.')))) %>%
  mutate(rich = specnumber(select(.,starts_with('infl.')))) %>%
  mutate(logAbun = log(inflCt))

# ssffCon0 = Data frame with merged ssCon (conspecific pollen) & ff (floral neighborhood) including only styles collected after pollinator observations were completed.

#ssffCon1a <- ssffCon0[ , 1:19]
#str(ssffCon1a)
#
# Option to create ssCon with only the day3 styles in rr.csv
# We did this, but currently (April 16, 2019) that is not what we are using to do the analyses below.

str(rr) # 136 obs
head(rr)

#want<-rr$styleID
#look<-ssffCon0[ssffCon0$styleID %in% want,]
#str(look) # 120 obs.
#tapply(look$styleID, look$fncDay, length)

#look2<-ssffCon0[ssffCon0$fncDay == 1,]
#str(look2) # 62 obs.

#head(look)
#head(look2)
#?rbind
#look3<-rbind(look,look2)
#str(look3)
#head(look3)
#tapply(look3$styleID, look3$fncDay, length)

#ssffCon<-look3
#str(ssffCon)
###ssffCon = Conspecific pollen on styles combined with floral neighborhood data for styles collected after pollinator observations & including only styles from focal plants with style persistance data.

# Creating ssffHet
str(ssHet) # 391 obs. of  4 variables
str(ff) # 224 obs. of  53 variables

ssffHet <- merge(ssHet, ff, by = 'plantObsID', all = T)
str(ssffHet) # 401 obs. of  56 variables
ssffHet$plantObsID <- as.factor(ssffHet$plantObsID)
ssffHet$styleID <- as.factor(ssffHet$styleID)
ssffHet$fncDay <- as.factor(ssffHet$fncDay)
ssffHet$site <- as.factor(ssffHet$site)
ssffHet$tag <- as.factor(ssffHet$tag)

# Create ssffHet with only the day3 styles in rr.csv

str(rr) # 136 obs
head(rr)

# want<-ssffCon2$styleID
# look4<-ssffHet[ssffHet$styleID %in% want,]
# str(look4) # 120 obs.
# tapply(look4$styleID, look4$fncDay, length)

# ssffHet<-look4
# str(ssffHet)
###ssffHet = Heterospecific pollen on styles combined with floral neighborhood data for styles collected after pollinator observations & including only styles from focal plants with style persistance data.

# 
# 
# # step 1: visualize                                                  ####
#
# #    look at pairwise relationships between fixed predictors 
# plot(jitter(ssffCon1$rich), ssffCon1$nn2)
# plot(jitter(ssffCon1$rich), jitter(as.integer(ssffCon1$fncDay)))
# plot(jitter(ssffCon1$rich), ssffCon1$logAbun)
# plot(ssffCon1$nn2, jitter(as.integer(ssffCon1$fncDay)))
# plot(ssffCon1$nn2, ssffCon1$logAbun)
# plot(ssffCon1$logAbun, jitter(as.integer(ssffCon1$fncDay)))
# 
# # #    look at site, the random effect 
# aggregate(plantObsID ~ site         , ssffCon1, length)
# aggregate(plantObsID ~ site + fncDay, ssffCon1, length)  # spp only on FNC day 3
# 
# # #    look at histograms of each fixed predictors per site 
# 
# ssffCon1 %>%
#   ggplot(aes(site, rich))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(site, nn2))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(site, logAbun))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# 
# #    look at pairwise relationships between preds w/in sites 
# 
# ssffCon1 %>%
#   ggplot(aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(rich, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(nn2, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(fncDay, logAbun))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(fncDay, nn2))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ssffCon1 %>%
#   ggplot(aes(fncDay, rich))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# # there don't seem to be any problems, so let's move forward....
# 
# 1. How does the floral neighborhood influence conspecific pollen count on styles? 
# # numSpro ~ floral neighborhood predictors + observation day , data = ssffCon1
# # (could be overdispersed), 

# How does the floral neighborhood influence conspecific pollen count on styles? ####
# assign models                                                      ####
ssffCon2 <- ssffCon1
ssffCon1 <- ssffCon2 %>%
  group_by(plantObsID)%>%
  mutate(hetAbun = round(mean(hetAbun))) %>%
  mutate(numSpro = round(mean(numSpro)))%>%
  filter(!duplicated(plantObsID))%>%
  mutate(rich = rich - mean(ssffCon2$rich)) %>%
  mutate(logAbun = logAbun - mean(ssffCon2$logAbun)) %>%
  mutate(nn2 = nn2 - mean(ssffCon2$nn2))
contrasts(ssffCon1$fncDay) <- contr.sum(3)


sproNOEA <- c('infl.acmi',
              'infl.ciar',
              'infl.cifl',
              'infl.copa',
              'infl.erst',
              'infl.phpi',
              'infl.somi')
stylesWNoea<- ssffHet[ssffHet$polCode %in% 'NOEA', 'styleID']
stylesToRem1 <- ssffCon[(rowSums(ssffCon1[,sproNOEA]) >0), 'styleID']
stylesToRemFinal <- stylesToRem1[!stylesToRem1 %in% stylesWNoea]
ssffCon3 <- ssffCon1[!ssffCon1$styleID %in% stylesToRemFinal, ]
#ssffCon1 <- ssffCon3 # you can toggle this on and off to look at results with and without
# sproNOEA pollen- which shows that the results are biologically equivalent even
# when we remove floral neighborhoods with lookalike pollen

x.1 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.2 <- glmmTMB(numSpro ~ nn2 + fncDay + logAbun +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.3 <- glmmTMB(numSpro ~ rich +  fncDay + logAbun +
               (1 | site), data = ssffCon1,
             family = nbinom2)

# because x.3 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.3.1 <- glmmTMB(numSpro ~    fncDay + logAbun +
                 (1 | site), data = ssffCon1,
               family = nbinom2)

x.3.2 <- glmmTMB(numSpro ~ rich +    logAbun +
                 (1 | site), data = ssffCon1,
               family = nbinom2)

x.3.3 <- glmmTMB(numSpro ~ rich +  fncDay +  
                 (1 | site), data = ssffCon1,
               family = nbinom2)

x.4 <- glmmTMB(numSpro ~ rich + nn2 + logAbun +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.5 <- glmmTMB(numSpro ~ rich + nn2 + fncDay +
               (1 | site), data = ssffCon1,
             family = nbinom2)

# because x.5 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.5.1 <- glmmTMB(numSpro ~      nn2 + fncDay +
                 (1 | site), data = ssffCon1,
               family = nbinom2)

x.5.2 <- glmmTMB(numSpro ~ rich +       fncDay +
                 (1 | site), data = ssffCon1,
               family = nbinom2)

x.5.3 <- glmmTMB(numSpro ~ rich + nn2 + 
                 (1 | site), data = ssffCon1,
               family = nbinom2)

x.6 <- glmmTMB(numSpro ~ rich + nn2 +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.7 <- glmmTMB(numSpro ~  fncDay + logAbun +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.8 <- glmmTMB(numSpro ~ rich + logAbun +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.9 <- glmmTMB(numSpro ~  nn2 + fncDay +
               (1 | site), data = ssffCon1,
             family = nbinom2)

x.10 <- glmmTMB(numSpro ~  rich + fncDay + 
                (1 | site), data = ssffCon1,
              family = nbinom2)

# because x.10 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.10.1 <- glmmTMB(numSpro ~       fncDay + 
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.10.2 <- glmmTMB(numSpro ~  rich +     
                  (1 | site), data = ssffCon1,
                family = nbinom2)


x.11 <- glmmTMB(numSpro ~  nn2 + logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.12 <- glmmTMB(numSpro ~  logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.13 <- glmmTMB(numSpro ~ rich +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.14 <- glmmTMB(numSpro ~  nn2 + 
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.15 <- glmmTMB(numSpro ~  fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.16 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                rich*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.17 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                rich*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.18 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.19 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

# because x.19 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.19.1 <- glmmTMB(numSpro ~        nn2 + fncDay + logAbun +
                  nn2*logAbun +
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.19.2 <- glmmTMB(numSpro ~ rich +       fncDay + logAbun +
                  nn2*logAbun +
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.19.3 <- glmmTMB(numSpro ~ rich + nn2 +          logAbun +
                  nn2*logAbun +
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.20 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.21 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun+
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.22 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.23 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.24 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + nn2*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.25 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 + rich*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.26 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 + rich*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.27 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                rich*fncDay + rich*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.28 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*nn2 +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.29 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.30 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.31 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay +  rich*nn2 + 
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.32 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + rich*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.33 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + rich*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.34 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*nn2 +
                (1 | site), data = ssffCon1,
              family = nbinom2)

# because x.34 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.34.1 <- glmmTMB(numSpro ~       nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*nn2 +
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.34.2 <- glmmTMB(numSpro ~ rich +      fncDay + logAbun +
                  nn2*logAbun + rich*nn2 +
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.34.3 <- glmmTMB(numSpro ~ rich + nn2 +         logAbun +
                  nn2*logAbun + rich*nn2 +
                  (1 | site), data = ssffCon1,
                family = nbinom2)

x.35 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*fncDay +
                (1 | site), data = ssffCon1,
              family = nbinom2)

x.36 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*logAbun +
                (1 | site), data = ssffCon1,
              family = nbinom2)
x.37 <- glmmTMB(numSpro ~ 
                  (1 | site), data = ssffCon1,
                family = nbinom2)

# compare models                                                     ####
i <- model.sel(x.1,x.2,x.3,x.4,x.5,x.6,x.7, x.8, x.9, x.10, x.11, x.12, 
            x.13, x.14,x.15,x.16, x.17, x.18, x.19, x.20, x.21, x.22, 
            x.23, x.24, x.25, x.26, x.27, x.28, x.29, x.30, x.31,
            x.32, x.33, x.34, x.35, x.36, x.37, x.19.1
            ) 

ms <- i

i <- 1:length(i$fncDay) # indices of columns with model terms
response <- "a"

res <- as.data.frame(ms)
v <- names(ms)[i]
v[v == "(Intercept)"] <- 1


mnames <- sapply(attr(ms, "modelList"), function(x) deparse(formula(x)))


res$model <- mnames
res1 <- as.data.frame(cbind(model = sapply(res$model, "[[", 1), 
                            df = print(res$df),
                            delta = print(res$delta), 
                            weight = print(res$weight)))
res1 <- as.data.frame(res1)
rownames(res1) <- c()
res1$delta1 <- unname(res1$delta)
res1$delta1 <- round(as.numeric(res1$delta1), digits = 2)
# make table of model selection

# sw(subset(i, delta <= 2))
# avg1 <- (model.avg(i, subset = delta < 2, fit = T))
# avgmod.95p <- model.avg(i, cumsum(weight) <= .95) 
# confint(avg1)


# model diagnostics                                                  ####
# some diagnostics, I followed the procedure here:
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

simulationOutput <- simulateResiduals(fittedModel = x.19.1, n =1000)
plot(simulationOutput)
testResiduals(simulationOutput)
plotResiduals(simulationOutput, ssffCon1$rich)
plotResiduals(simulationOutput, ssffCon1$fncDay)
testDispersion(simulationOutput)
testQuantiles(simulationOutput)

simulationOutput$scaledResiduals



Anova(x.19.1)
r.squaredGLMM(x.19.1)
plot(allEffects(x.19.1))
summary(allEffects(x.19.1))

# The final model is x.10.4 because that deals wthe issue of residuals vs. predicted variables (see dispformula component).

#


rich <- seq(0,9,1)
nn2 <- seq(0, 11, .2)
logAbun <- seq(0.2,8.5, 0.2)
fncDay <- as.factor(unique(ssffCon1$fncDay))
makeLook <- expand.grid(nn2, logAbun, unique(ssffCon1$site), fncDay, rich)
names(makeLook) <- c('nn2' , 'logAbun', 'site', 'fncDay','rich')

makeLook <- makeLook %>%
  mutate(rich = makeLook$rich - mean(ssffCon2$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun - mean(ssffCon2$logAbun)) %>%
  mutate(nn2 = nn2 - mean(ssffCon2$nn2))
contrasts(makeLook$fncDay) <- contr.sum(3)
fit <- predict(avg1, newdata = makeLook, type="response")
makeLook$numSpro <- fit
fit1 <- predict(x.19.1, newdata = makeLook, type="response")
makeLook$numSpro1 <- fit1

# Make Figure S4 ####
makeLook %>%
  filter(fncDay %in% '3') %>%
  mutate(rich = rich + mean(ssffCon2$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun + mean(ssffCon2$logAbun)) %>%
  mutate(nn2 = nn2 + mean(ssffCon2$nn2)) %>%
  ggplot(aes(nn2, logAbun)) +
  #facet_wrap(~fncDay) +
  geom_raster(aes(fill = numSpro1), interpolate = T) +
  theme_classic() +
  geom_jitter(data = ssffCon2, pch = 21, aes(fill = numSpro)) +
  scale_fill_gradient(low = 'white', high = 'black')+
  labs(x = 'Distance (m) to second nearest flowering conspecific',
       y = 'Count of inflorescences',
       fill = 'Conspecific\\npollen\\ncount')+
  scale_y_continuous(labels = c(round(exp(0)), round(exp(2)), 
                                round(exp(4)), round(exp(6)), 
                                round(exp(8))))




  



# 2. How does the floral neighborhood influence heterospecific pollen count on styles? ####
# # hetAbun ~ floral neighborhood predictors + observation day , data = ssffCon1

# Analysis of heterospecific pollen abundance (hetAbun)
names(ssffHet)
# hetAbun: number of heterospecific pollen grains on styles

hist(ssffCon1$hetAbun) # might do a zeroinflation model
par(mfcol=c(1,1), las=2)
mean(ssffCon1$hetAbun)
median(ssffCon1$hetAbun)



# assign models                                                      ####

y.1 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                 (1 | site), zi=~rich + nn2 + fncDay + logAbun, data =
                 ssffCon1, family = nbinom2)

y.2 <- glmmTMB(hetAbun ~ nn2 + fncDay + logAbun +
                 (1 | site), zi=~nn2 + fncDay + logAbun, data = ssffCon1,
               family = nbinom2)

y.3 <- glmmTMB(hetAbun ~ rich +  fncDay + logAbun +
                 (1 | site), zi=~rich + fncDay + logAbun, data = ssffCon1,
               family = nbinom2)

y.4 <- glmmTMB(hetAbun ~ rich + nn2 + logAbun +
                 (1 | site), zi=~rich + nn2 + logAbun, data = ssffCon1,
               family = nbinom2)

y.5 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay +
                 (1 | site), zi=~rich + nn2 + fncDay, data = ssffCon1,
               family = nbinom2)

y.6 <- glmmTMB(hetAbun ~ rich + nn2 +
                 (1 | site), zi=~rich + nn2, data = ssffCon1,
               family = nbinom2)

y.7 <- glmmTMB(hetAbun ~  fncDay + logAbun +
                 (1 | site), zi=~fncDay + logAbun, data = ssffCon1,
               family = nbinom2)

y.8 <- glmmTMB(hetAbun ~ rich + logAbun +
                 (1 | site), zi=~rich + logAbun, data = ssffCon1,
               family = nbinom2)

y.9 <- glmmTMB(hetAbun ~  nn2 + fncDay +
                 (1 | site), zi=~ nn2 + fncDay , data = ssffCon1,
               family = nbinom2)

y.10 <- glmmTMB(hetAbun ~  rich + fncDay + 
                  (1 | site), zi=~rich  + fncDay, data = ssffCon1,
                family = nbinom2)

y.11 <- glmmTMB(hetAbun ~  nn2 + logAbun +
                  (1 | site), zi=~nn2 + logAbun, data = ssffCon1,
                family = nbinom2)

y.12 <- glmmTMB(hetAbun ~  logAbun +
                  (1 | site), zi=~logAbun, data = ssffCon1,
                family = nbinom2)

y.13 <- glmmTMB(hetAbun ~ rich +
                  (1 | site), zi=~rich, data = ssffCon1,
                family = nbinom2)

y.14 <- glmmTMB(hetAbun ~  nn2 + 
                  (1 | site), zi=~nn2, data = ssffCon1,
                family = nbinom2)

y.15 <- glmmTMB(hetAbun ~  fncDay +
                  (1 | site), zi=~fncDay, data = ssffCon1,
                family = nbinom2)

y.16 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  rich*logAbun, data = ssffCon1,
                family = nbinom2)

y.17 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  rich*fncDay, data = ssffCon1,
                family = nbinom2)

y.18 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  rich*nn2, data = ssffCon1,
                family = nbinom2)


y.19 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*logAbun, data = ssffCon1,
                family = nbinom2)

y.20 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*fncDay, data = ssffCon1,
                family = nbinom2)

y.21 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun, data = ssffCon1,
                family = nbinom2)

y.22 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay, data = ssffCon1,
                family = nbinom2)

y.23 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*logAbun, data = ssffCon1,
                family = nbinom2)

y.24 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + nn2*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + nn2*logAbun, data = ssffCon1,
                family = nbinom2)

y.25 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*fncDay, data = ssffCon1,
                family = nbinom2)

y.26 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*logAbun, data = ssffCon1,
                family = nbinom2)

y.27 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*fncDay + rich*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  rich*fncDay + rich*logAbun, data = ssffCon1,
                family = nbinom2)

y.28 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*nn2 +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*nn2, data = ssffCon1,
                family = nbinom2)

y.29 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*fncDay, data = ssffCon1,
                family = nbinom2)

y.30 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*logAbun, data = ssffCon1,
                family = nbinom2)

y.31 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +  rich*nn2 + 
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +  rich*nn2, data = ssffCon1,
                family = nbinom2)

y.32 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*fncDay, data = ssffCon1,
                family = nbinom2)

y.33 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*logAbun, data = ssffCon1,
                family = nbinom2)

y.34 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*nn2 +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*nn2, data = ssffCon1,
                family = nbinom2)

y.35 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*fncDay +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*fncDay, data = ssffCon1,
                family = nbinom2)

y.36 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*logAbun +
                  (1 | site), zi=~rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*logAbun, data = ssffCon1,
                family = nbinom2)






# all of those look fine

y.n <- glmmTMB(hetAbun ~  1 + 
                 (1 | site), data = ssffCon1,
               family = nbinom2)

y.n.1 <- glmmTMB(hetAbun ~  1 + 
                   (1 | site), zi= ~1, data = ssffCon1,
                 family = nbinom2)
# compare models ####
g <- model.sel(y.1,y.2,y.3,y.4,y.5,y.6,y.7, y.8, y.9, y.10, y.11, y.12, 
            y.13, y.14,y.15,y.16, y.17, y.18, y.19, y.20, y.21, y.22, 
            y.23, y.24, y.25, y.26, y.27, y.28, y.29, y.30, y.31,
            y.32, y.33, y.34, y.35, y.36, y.n) 
g$prAIC <- round(exp((min(g$AICc, na.rm = T)-g$AICc)/2),4)

ms <- g

i <- 1:length(g$fncDay) # indices of columns with model terms
response <- "a"

res <- as.data.frame(ms)
v <- names(ms)[i]
v[v == "(Intercept)"] <- 1


mnames <- sapply(attr(ms, "modelList"), function(x) deparse(formula(x)))


res$model <- mnames
res1 <- as.data.frame(cbind(model = sapply(res$model, "[[", 1), 
                            df = print(res$df),
                            delta = print(res$delta), 
                            weight = print(res$weight)))
res1 <- as.data.frame(res1)
rownames(res1) <- c()
res1$delta1 <- unname(res1$delta)
res1$delta1 <- round(as.numeric(res1$delta1), digits = 2)
# Write model selection results table
# as.data.frame(res1) %>%
#   select(model, df, delta, weight) %>%
#   write.csv(file = '/Users/learichardson/Desktop/sffHet.csv')

Anova(y.13)
r.squaredGLMM(y.13)
plot(allEffects(y.13))
summary(allEffects(y.13))
avg1 <- (model.avg(g, subset = delta < 2, fit = T))



# I guess it doesn't need a figure then.


