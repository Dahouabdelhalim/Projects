# The goal of this script is to run 2 analyses that examine the relationship between pollen
# collected off of pollinators bodies and the floral neighborhood predictors

# 1. How does the floral neighborhood influence heterospecific pollen count on pollinators? 
# # hetAbun ~ floral neighborhood predictors + observation day 
# # (could be overdispersed), 

# 2. How does the floral neighborhood influence conspecific pollen count on pollinators?
# # numSpro ~ floral neighborhood predictors + observation day 

# load packages                                                      ####
library(lme4)
library(tidyverse)
library(vegan)
library(DHARMa)
library(bbmle) 
library(glmmTMB)
library(MuMIn)
library(effects)
library(ggpubr)
library(car)
set.seed(66)

# get datasets                                                       ####
ff <- read.csv('ff.csv')
bb <- read.csv('bb.csv')
abunRich <-read.csv('ppCon.csv')
levels <- levels(abunRich$sp)
levels[length(levels) + 1] <- "None"
abunRich$sp <- factor(abunRich$sp, levels = levels)
abunRich$sp <- as.character(abunRich$sp)
abunRich[is.na(abunRich$sp), "sp"] <- "unidentified"
abunRich <- na.omit(abunRich)
pp1 <-abunRich[abunRich['totAbun']!=0,]
pp2 <-pp1[pp1['totAbun']!=1,]
pp3 <-pp2[pp2['sp']!='Coleoptera sp. 2',]
pp4 <-pp3[pp3['sp']!="Coelioxys rufitarsis",]
pp5 <-pp4[pp4['sp']!="Megachile sp. 1",]
pp6 <-pp5[pp5['sp']!="Megachile sp. 2",]
pp7 <-pp6[pp6['sp']!="unidentified",]
pp7$sp<-factor(pp7$sp)
pp <- pp7
remove(pp1,pp2,pp3,pp4,pp5,pp6,pp7,abunRich)

# match the bees to the floral neighborhoods they came from
bff <- merge(bb, ff, by = 'plantObsID')

# make count and richness columns for floral neighborhood
bff <- bff %>%
  mutate(fncDay = as.factor(bff$fncDay)) %>% #make day a factor
  mutate(infl.ecan = as.integer(infl.ecan + 1)) %>% # add one ech infl for the focal
  mutate(inflCt = rowSums(select(.,starts_with('infl.')))) %>%
  mutate(rich = specnumber(select(.,starts_with('infl.')))) %>%
  mutate(logAbun = log(inflCt)) %>%
  select(!site.x) %>%
  rename(site = 'site.y')
  
ppff <-  merge(bff, pp, by = 'specimenID')


# step 1: visualize                                                  ####


# 
# 
# ppff %>%
#   ggplot(aes(site, rich))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(site, nn2))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(site, logAbun))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# 
# #    look at pairwise relationships between preds w/in sites 
# 
# 
# 
# ppff %>%
#   ggplot(aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(rich, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(nn2, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(fncDay, logAbun))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(fncDay, nn2))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# ppff %>%
#   ggplot(aes(fncDay, rich))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()

# there don't seem to be any problems, so let's move forward....

# assign models             ####

##################################################
####### MOD 1, heterospecific pollen count #######
##################################################

# ppff <- ppff[!ppff$specimenID %in% specsToRemoveFinal$specimenID, ] # You can use
# this line of code to remove specimens collected from sites with NOEA species in the 
# floral neighborhood. But to run this like you need to source the script:
# pollinatorCommunity.R

ppff <- ppff %>%
  group_by(plantObsID)%>%
  mutate(hetAbun = round(mean(hetAbun)))%>%
  mutate(numSpro = round(mean(numSpro)))%>%
  filter(!duplicated(plantObsID))

ppff1 <- ppff

ppff <- ppff %>%
  mutate(rich = rich - mean(bff$rich)) %>%
  mutate(logAbun = logAbun - mean(bff$logAbun)) %>%
  mutate(nn2 = nn2 - mean(bff$nn2))
contrasts(ppff$fncDay) <- contr.sum(3)

x.1 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)

x.2 <- glmmTMB(hetAbun ~ nn2 + fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)
x.3 <- glmmTMB(hetAbun ~ rich +  fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)
x.4 <- glmmTMB(hetAbun ~ rich + nn2 + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)
x.5 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay +
                 (1 | site), data = ppff,
               family = nbinom2)

x.6 <- glmmTMB(hetAbun ~ rich + nn2 +
                 (1 | site), data = ppff,
               family = nbinom2)

x.7 <- glmmTMB(hetAbun ~  fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)

x.8 <- glmmTMB(hetAbun ~ rich + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)

x.9 <- glmmTMB(hetAbun ~  nn2 + fncDay +
                 (1 | site), data = ppff,
               family = nbinom2)

x.10 <- glmmTMB(hetAbun ~  rich + fncDay + 
                  (1 | site), data = ppff,
                family = nbinom2)

x.11 <- glmmTMB(hetAbun ~  nn2 + logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.12 <- glmmTMB(hetAbun ~  logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.13 <- glmmTMB(hetAbun ~ rich +
                  (1 | site), data = ppff,
                family = nbinom2)

x.14 <- glmmTMB(hetAbun ~  nn2 + 
                  (1 | site), data = ppff,
                family = nbinom2)

x.15 <- glmmTMB(hetAbun ~  fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.16 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.17 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.18 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 +
                  (1 | site), data = ppff,
                family = nbinom2)

x.19 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.20 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.21 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), data = ppff,
                family = nbinom2)


x.22 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.23 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.24 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + nn2*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.25 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.26 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.27 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  rich*fncDay + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.28 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*nn2 +
                  (1 | site), data = ppff,
                family = nbinom2)

x.29 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.30 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.31 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +  rich*nn2 + 
                  (1 | site), data = ppff,
                family = nbinom2)

x.32 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

x.33 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.34 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*nn2 +
                  (1 | site), data = ppff,
                family = nbinom2)

x.35 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

#take out nn2
x.35.1 <- glmmTMB(hetAbun ~ rich + fncDay + logAbun +
                    nn2*logAbun + rich*fncDay +
                    (1 | site), data = ppff,
                  family = nbinom2)

#take out nn2 and nn2*logabun (leave logAbun)
x.35.2 <- glmmTMB(hetAbun ~ rich + fncDay + logAbun +
                    rich*fncDay +
                    (1 | site), data = ppff,
                  family = nbinom2)

# ONLY richness*day interaction
x.35.3 <- glmmTMB(hetAbun ~ rich + fncDay + 
                    rich*fncDay +
                    (1 | site), data = ppff,
                  family = nbinom2)

x.36 <- glmmTMB(hetAbun ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

x.37 <- glmmTMB(hetAbun ~ (1 | site), data = ppff,
                family = nbinom2)
# compare models                                                     ####

h <- model.sel(x.1,x.2,x.3,x.4,x.5,x.6,x.7, x.8, x.9, x.10, x.11, x.12, 
               x.13, x.14,x.15,x.16, x.17, x.18, x.19, x.20, x.21, x.22, 
               x.23, x.24, x.25, x.26, x.27, x.28, x.29, x.30, x.31,
               x.32, x.33, x.34, x.35, x.36, x.35.1, x.35.2, x.35.3, x.37)

sw(subset(h, delta <= 2))
avg <- model.avg(h, subset = delta <= 2, fit = TRUE)
avgmod.95p <- model.avg(h, cumsum(weight) <= .95) 
confint(avgmod.95p)

ms <- h

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
# Make model selection table 
# as.data.frame(res1) %>%
#   select(model, df, delta, weight) %>%
#   write.csv(file = '/Users/learichardson/Desktop/pffHet.csv')



summary(x.35.2)
Anova(x.35.2)
r.squaredGLMM(x.35.2)
plot(allEffects(x.35.2))
summary(allEffects(x.35.2))


simulationOutput <- simulateResiduals(fittedModel = x.35.2, n =1000)
plot(simulationOutput)
testResiduals(simulationOutput)
# great!



makeLook1 %>%
  group_by(fncDay) %>%
  summarize_at(vars(rich, logAbun, nn2), funs(min, max))
r.squaredGLMM(x.35.2, x.37)
plot(allEffects(x.35.2))
summary(allEffects(x.35.2))


nn2 <- seq(0.2, 11, 1)
rich <- seq(1, 10 , .5)
logAbun <- seq(0.2,8.5, 1)
fncDay <- as.factor(c(1,2,3))
makeLook <- expand.grid(rich, unique(ppff$site), fncDay, nn2, logAbun)
names(makeLook) <- c('rich' , 'site', 'fncDay', 'nn2', 'logAbun')
makeLook <- makeLook %>%
  mutate(rich = rich - mean(ppff1$rich)) %>%
  mutate(logAbun = logAbun - mean(ppff1$logAbun)) %>%
  mutate(nn2 = nn2 - mean(ppff1$nn2))
contrasts(makeLook$fncDay) <- contr.sum(3)
fit <- predict(x.35.2, newdata = makeLook, type="response")
makeLook$hetAbun <- fit

makeLook1 <- makeLook %>%
  filter(fncDay %in% '1' & between(rich, -1.68, 5.32) |
           fncDay %in% '2' & between(rich,-2.68, 2.32) | 
           fncDay %in% '3' & between(rich, -1.68, 2.32))
# Make Figure 3b ####
a <- makeLook1 %>%
  mutate(rich = rich + mean(ppff1$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun + mean(ppff1$logAbun)) %>%
  mutate(nn2 = nn2 + mean(ppff1$nn2)) %>%
  ggplot(aes(rich, hetAbun, linetype=fncDay)) +
  geom_jitter(data = ppff1, aes(shape=fncDay, alpha = 0.8))+
  geom_smooth(color = 'black', method = 'loess') + 
  scale_linetype_manual(values=c(1,4,5))+
  scale_shape_manual(values=c(1,3,4))+
  scale_x_continuous(breaks=c(2,4,6,8))+
  theme_classic()+
  labs(x = 'Floral species richness',
       y = 'Heterospecific pollen count')






###############################################
####### MOD 2, conspecific pollen count #######
###############################################

z.1 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)

z.2 <- glmmTMB(numSpro ~ nn2 + fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)
z.3 <- glmmTMB(numSpro ~ rich +  fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)
z.4 <- glmmTMB(numSpro ~ rich + nn2 + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)
z.5 <- glmmTMB(numSpro ~ rich + nn2 + fncDay +
                 (1 | site), data = ppff,
               family = nbinom2)

z.6 <- glmmTMB(numSpro ~ rich + nn2 +
                 (1 | site), data = ppff,
               family = nbinom2)

z.7 <- glmmTMB(numSpro ~  fncDay + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)

z.8 <- glmmTMB(numSpro ~ rich + logAbun +
                 (1 | site), data = ppff,
               family = nbinom2)

z.9 <- glmmTMB(numSpro ~  nn2 + fncDay +
                 (1 | site), data = ppff,
               family = nbinom2)

z.10 <- glmmTMB(numSpro ~  rich + fncDay + 
                  (1 | site), data = ppff,
                family = nbinom2)

z.11 <- glmmTMB(numSpro ~  nn2 + logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.12 <- glmmTMB(numSpro ~  logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.13 <- glmmTMB(numSpro ~ rich +
                  (1 | site), data = ppff,
                family = nbinom2)

z.14 <- glmmTMB(numSpro ~  nn2 + 
                  (1 | site), data = ppff,
                family = nbinom2)

z.15 <- glmmTMB(numSpro ~  fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.16 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.17 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.18 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 +
                  (1 | site), data = ppff,
                family = nbinom2)

z.19 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.20 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.21 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), data = ppff,
                family = nbinom2)


z.22 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.23 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.24 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + nn2*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.25 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.26 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  rich*nn2 + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.27 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  rich*fncDay + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.28 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*nn2 +
                  (1 | site), data = ppff,
                family = nbinom2)

z.29 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.30 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.31 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay +  rich*nn2 + 
                  (1 | site), data = ppff,
                family = nbinom2)

z.32 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

z.33 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*fncDay + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)

z.34 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*nn2 +
                  (1 | site), data = ppff,
                family = nbinom2)

z.35 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*fncDay +
                  (1 | site), data = ppff,
                family = nbinom2)

#take out nn2
z.35.1 <- glmmTMB(numSpro ~ rich + fncDay + logAbun +
                    nn2*logAbun + rich*fncDay +
                    (1 | site), data = ppff,
                  family = nbinom2)

#take out nn2 and nn2*logabun (leave logAbun)
z.35.2 <- glmmTMB(numSpro ~ rich + fncDay + logAbun +
                    rich*fncDay +
                    (1 | site), data = ppff,
                  family = nbinom2)

# ONLY richness*day interaction
z.35.3 <- glmmTMB(numSpro ~ rich + fncDay + 
                    rich*fncDay +
                    (1 | site), data = ppff,
                  family = nbinom2)

z.36 <- glmmTMB(numSpro ~ rich + nn2 + fncDay + logAbun +
                  nn2*logAbun + rich*logAbun +
                  (1 | site), data = ppff,
                family = nbinom2)
z.37 <- glmmTMB(numSpro ~ 
                  (1 | site), data = ppff,
                family = nbinom2)

# compare models                                                     ####

i <- model.sel(z.1,z.2,z.3,z.4,z.5,z.6,z.7, z.8, z.9, z.10, z.11, z.12, 
            z.13, z.14,z.15,z.16, z.17, z.18, z.19, z.20, z.21, z.22, 
            z.23, z.24, z.25, z.26, z.27, z.28, z.29, z.30, z.31,
            z.32, z.33, z.34, z.35, z.36, z.35.1, z.35.2, z.35.3, z.37) 

ms <- i

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
# Make model selection table 
# as.data.frame(res1) %>%
#   select(model, df, delta, weight) %>%
#   write.csv(file = '/Users/learichardson/Desktop/pffCon.csv')

z.10
summary(z.10)

sw(subset(i, delta <= 2))
avgi <- model.avg(i, subset = delta < 2, fit = TRUE)
avgmod.95p <- model.avg(h, cumsum(weight) <= .95) 
confint(avgi)


r.squaredGLMM(z.10)
plot(allEffects(z.10))
summary(allEffects(z.10))



nn2 <- seq(0.2, 11, 1)
rich <- seq(1, 10 , .5)
logAbun <- seq(0.2,8.5, 1)
fncDay <- as.factor(c(1,2,3))
makeLook <- expand.grid(rich, unique(ppff$site), fncDay, nn2, logAbun)
names(makeLook) <- c('rich' , 'site', 'fncDay', 'nn2', 'logAbun')
makeLook <- makeLook %>%
  mutate(rich = rich - mean(ppff1$rich)) %>%
  mutate(logAbun = logAbun - mean(ppff1$logAbun)) %>%
  mutate(nn2 = nn2 - mean(ppff1$nn2))
contrasts(makeLook$fncDay) <- contr.sum(3)
fit <- predict(z.10, newdata = makeLook, type="response")
makeLook$numSpro <- fit

makeLook1 <- makeLook %>%
  filter(fncDay %in% '1' & between(rich, -1.68, 5.32) |
           fncDay %in% '2' & between(rich,-2.68, 2.32) | 
           fncDay %in% '3' & between(rich, -1.68, 2.32))
# Make Figure 3a ####
b <- makeLook1 %>%
  mutate(rich = rich + mean(ppff1$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun + mean(ppff1$logAbun)) %>%
  mutate(nn2 = nn2 + mean(ppff1$nn2)) %>%
  ggplot(aes(rich, numSpro, linetype=fncDay)) +
  geom_jitter(data = ppff1, aes(shape=fncDay, alpha = 0.8))+
  geom_smooth(color = 'black', method = 'loess') + 
  scale_linetype_manual(values=c(1,4,5))+
  scale_shape_manual(values=c(1,3,4))+
  theme_classic()+
  scale_x_continuous(breaks=c(2,4,6,8))+
  labs(x = 'Floral species richness',
       y = 'Conspecific pollen count')

# Make overall Figure 3 ####
ggarrange(b, a, label.x = 'A.', label.y = 'B.',
          ncol = 2, common.legend = T)
