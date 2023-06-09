# The goal of this script is to run 2 analyses that examine the relationship between 
# style persistence and the floral neighborhood predictors

# 1. Style persistence ~ floral neighborhood predictors + observation day

# load packages                                                      ####
library(lme4)
library(tidyverse)
library(vegan)
library(DHARMa)
library(bbmle) 
library(MuMIn)
library(effects)
library(car)
set.seed(66)

# get datasets                                                       ####

rr <- read.csv("rr.csv", stringsAsFactors = F)
ff <- read.csv('ff.csv', stringsAsFactors = F)
rrff <- merge(rr[, c('plantObsID', 'spLevel')], ff)
rff<- rrff
rrff <- rrff %>%
  mutate(fncDay = as.factor(rrff$fncDay)) %>%
  mutate(infl.ecan = as.integer(infl.ecan + 1)) %>% # add one ech infl for the focal
  mutate(inflCt = rowSums(select(.,starts_with('infl.')))) %>%
  mutate(rich = specnumber(select(.,starts_with('infl.')))) %>%
  mutate(logAbun = log(inflCt))

rff <- rrff
rrff <- rrff %>%
  mutate(rich = rich - mean(rrff$rich)) %>%
  mutate(logAbun = logAbun - mean(rrff$logAbun)) %>%
  mutate(nn2 = nn2 - mean(rrff$nn2))
contrasts(rrff$fncDay) <- contr.sum(2)

# step 1: visualize                                                  ####



# #    how does logAbun look? 
# hist(rrff$logAbun, 20)
# 
# #    look at pairwise relationships between fixed predictors 
# plot(jitter(rrff$rich), rrff$nn2)
# plot(jitter(rrff$rich), jitter(as.integer(rrff$fncDay)))
# plot(jitter(rrff$rich), rrff$logAbun)
# plot(rrff$nn2, jitter(as.integer(rrff$fncDay)))
# plot(rrff$nn2, rrff$logAbun)
# plot(rrff$logAbun, jitter(as.integer(rrff$fncDay)))
# 
# # #    look at site, the random effect 
# aggregate(plantObsID ~ site         , rrff, length)
# aggregate(plantObsID ~ site + fncDay, rrff, length)
# 
# # #    look at histograms of each fixed predictors per site 
# 
# rrff %>%
#   ggplot(aes(site, rich))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(site, nn2))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(site, logAbun))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# 
# #    look at pairwise relationships between preds w/in sites 
# 
# rrff %>%
#   ggplot(aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(rich, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(nn2, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(fncDay, logAbun))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(fncDay, nn2))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# rrff %>%
#   ggplot(aes(fncDay, rich))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()

# there don't seem to be any problems, so let's move forward....

# assign models                                                      ####

x.1 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.2 <- glmer(spLevel ~ nn2 + fncDay + logAbun +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.3 <- glmer(spLevel ~ rich +  fncDay + logAbun +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)
x.4 <- glmer(spLevel ~ rich + nn2 + logAbun +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)
x.5 <- glmer(spLevel ~ rich + nn2 + fncDay +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.6 <- glmer(spLevel ~ rich + nn2 +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.7 <- glmer(spLevel ~  fncDay + logAbun +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.8 <- glmer(spLevel ~ rich + logAbun +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.9 <- glmer(spLevel ~  nn2 + fncDay +
               (1 | site), data = rrff,
             family = binomial, nAGQ = 2)

x.10 <- glmer(spLevel ~  rich + fncDay + 
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.11 <- glmer(spLevel ~  nn2 + logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.12 <- glmer(spLevel ~  logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.13 <- glmer(spLevel ~ rich +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.14 <- glmer(spLevel ~  nn2 + 
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.15 <- glmer(spLevel ~  fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.16 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                rich*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.17 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 5)

x.18 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.19 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.20 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 5)

x.21 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun+
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

# because x.21 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.21.1 <- glmer(spLevel ~        nn2 + fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 2)

x.21.2 <- glmer(spLevel ~ rich       + fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 2)

x.21.3 <- glmer(spLevel ~             fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 2)

x.21.4 <- glmer(spLevel ~             fncDay + logAbun +
                  fncDay*logAbun+
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 5)

x.22 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 3)

# because x.22 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important

x.22.1 <- glmer(spLevel ~        nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 2)

x.22.2 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)
x.22.3 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 1)

x.22.4 <- glmer(spLevel ~        nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 1)

x.22.5 <- glmer(spLevel ~        nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 3)

x.22.6 <- glmer(spLevel ~        nn2 + fncDay + logAbun +
                  fncDay*logAbun + nn2*fncDay +
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 5)

x.23 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 3)

x.24 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + nn2*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.25 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 + rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)



x.27 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                rich*fncDay + rich*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)

x.28 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*nn2 +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.29 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 1)
x.29.1 <- glmer(spLevel ~ rich +     fncDay + logAbun +
                fncDay*logAbun + rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 1)

x.29.2 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + 
                (1 | site), data = rrff,
              family = binomial, nAGQ = 1)

x.29.3 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                                 rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 1)

x.29.4 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + 
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 2)

x.29.5 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                  fncDay*logAbun + 
                  (1 | site), data = rrff,
                family = binomial, nAGQ = 4)

x.30 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)

x.31 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay +  rich*nn2 + 
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.32 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 1)

x.33 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + rich*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)

x.34 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*nn2 +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)

x.35 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*fncDay +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 2)

x.36 <- glmer(spLevel ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*logAbun +
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)

x.37 <- glmer(spLevel ~ 
                (1 | site), data = rrff,
              family = binomial, nAGQ = 0)
# compare models                                                     ####
# initially I compared models 1-36, but after seeing that the top models
# were x.22, x.21, x.29, x.23, and x.28, I went back to add a few extra models 
# that would help us figure out which predictors are really the most important

g <- model.sel(x.1,x.2,x.3,x.4,x.5,x.6,x.7, x.8, x.9, x.10, x.11, x.12, 
            x.13, x.14,x.15,x.16, x.17, x.18, x.19, x.20, x.21, x.22, 
            x.23, x.24, x.25, x.27, x.28, x.29, x.30, x.31,
            x.32, x.33, x.34, x.35, x.36,
             x.21.1, x.21.2, x.21.3, x.37) #x.21.4,
            # x.22.1, x.22.2, x.22.3, x.22.4, x.22.5, x.22.6,
            # x.29.1, x.29.2, x.29.3, x.29.4, x.29.5
            # ) 

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
# Make model selection table 
# as.data.frame(res1) %>%
#   select(model, df, delta, weight) %>%
#   write.csv(file = '/Users/learichardson/Desktop/rff.csv')


# sw(subset(i, delta <= 2))
# avgi <- model.avg(g, subset = delta < 2, fit = TRUE)
# avgmod.95p <- model.avg(h, cumsum(weight) <= .95)
# confint(avgi)

Anova(x.22)
r.squaredGLMM(x.22)
plot(allEffects(x.22))
summary(allEffects(x.22))

g$prAIC <- round(exp((min(g$dAIC)-g$dAIC)/2),4)

# model diagnostics                                                  ####
# some diagnostics, I followed the procedure here:
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# simulationOutput <- simulateResiduals(fittedModel = x.22, n =1000)
# plot(simulationOutput)
# testResiduals(simulationOutput)

# make Figure 5 ####
nn2 <- seq(0, 11.1, .2)
rich <- seq(-0.2, 10, .2)
logAbun <- seq(0,8.3, .2)
fncDay <- as.factor(c(2,3))
makeLook <- expand.grid(nn2, rich, unique(ff$site), logAbun, fncDay)
names(makeLook) <- c('nn2' , 'rich', 'site', 'logAbun', 'fncDay')
makeLook <- makeLook %>%
  mutate(rich = rich - mean(rff$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun - mean(rff$logAbun)) %>%
  mutate(nn2 = nn2 - mean(rff$nn2, na.rm = T))
contrasts(makeLook$fncDay) <- contr.sum(2)
fit <- predict(x.22, newdata = makeLook, type="response")
makeLook$styleP <- fit


makeLook %>%
  mutate(rich = rich + mean(rff$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun + mean(rff$logAbun)) %>%
  mutate(nn2 = nn2 + mean(rff$nn2)) %>%
  ggplot(aes(nn2, rich)) +
  facet_wrap(~fncDay) +
  geom_raster(aes(fill = styleP), interpolate = T) +
  theme_classic() +
  geom_jitter(data = rff, pch = 21, 
              col = ifelse(rrff$fncDay == '2', 'white', 'black'), 
              aes(fill = spLevel)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  labs(x = 'Distance (m) to second nearest flowering conspecific',
       y = 'Floral species richness',
       fill = 'Style\\npersistence')+
  scale_y_continuous(breaks=c(0,3,6,9))


                                                                   



