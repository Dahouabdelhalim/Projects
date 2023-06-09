# to what extent do elements of fl neighborhood (richess, nn2, fncDay, abundance, and site)
# predict style persistence?

# load packages                                                      ####
library(lme4) # version 1.1.23
library(tidyverse) # version 1.3.0
library(vegan) # version 2.5.6
library(DHARMa) # version 0.3.1
library(bbmle) # version 1.0.23.1
library(effects)
library(MuMIn)
set.seed(66)
# created using R version 3.6.3


# get datasets                                                       ####
ff <- read.csv('ff.csv', stringsAsFactors = F)
bb <- read.csv('bb.csv', stringsAsFactors = F)
bff <- merge(bb, ff, by = c('plantObsID', 'site'))
bff <- bff %>%
  mutate(fncDay = as.factor(bff$fncDay)) %>%
  mutate(infl.ecan = as.integer(infl.ecan + 1)) %>% # add one ech infl for the focal
  mutate(inflCt = rowSums(select(.,starts_with('infl.')))) %>%
  mutate(rich = specnumber(select(.,starts_with('infl.')))) %>%
  mutate(logAbun = log(inflCt))

# step 1: visualize                                                  ####
# 
# #    look at histograms of each fixed predictors per site 
# 
# bff %>%
#   ggplot(aes(site, rich))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(site, nn2))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(site, logAbun))+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()


# this calls for some tidyverse magic!

#    look at pairwise relationships between preds w/in sites 

# bff %>%
#   ggplot(aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(rich, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(nn2, logAbun))+
#   facet_wrap(~site)+
#   geom_point() +
#   geom_smooth(method = lm)+
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(fncDay, logAbun))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(fncDay, nn2))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()
# 
# bff %>%
#   ggplot(aes(fncDay, rich))+
#   facet_wrap(~site)+
#   geom_violin()+
#   geom_jitter(color = 'black', alpha = 1/5) +
#   theme_classic()

# there don't seem to be any problems, so let's move forward....

# assign models                                                      ####

bff1 <- bff %>% 
  group_by(plantObsID) %>% 
  mutate(totalPols = sum(polObserved)) %>%
  add_tally() %>%
  mutate(totalObs = n) %>%
  mutate(totalFails = n - totalPols) %>%
  select(totalFails, totalPols, rich, nn2, fncDay, logAbun, site, plantObsID, polObserved) %>%
  distinct() %>%
  mutate(rich = rich - mean(bff$rich)) %>%
  mutate(logAbun = logAbun - mean(bff$logAbun)) %>%
  mutate(nn2 = nn2 - mean(bff$nn2))
contrasts(bff1$fncDay) <- contr.sum(3)

y <- cbind(bff1$totalPols, bff1$totalFails)

x.1 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)

x.2 <- glmer(y ~ nn2 + fncDay + logAbun +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)
x.3 <- glmer(y ~ rich +  fncDay + logAbun +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)
x.4 <- glmer(y ~ rich + nn2 + logAbun +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)
x.5 <- glmer(y ~ rich + nn2 + fncDay +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)

x.6 <- glmer(y ~ rich + nn2 +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)

x.7 <- glmer(y ~  fncDay + logAbun +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)

x.8 <- glmer(y ~ rich + logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.9 <- glmer(y ~  nn2 + fncDay +
                 (1 | site), data = bff1,
               family = binomial, nAGQ = 2)

x.10 <- glmer(y ~  rich + fncDay + 
                 (1 | site), data = bff1,
               family = binomial, nAGQ = 2)

x.11 <- glmer(y ~  nn2 + logAbun +
                 (1 | site), data = bff1,
               family = binomial, nAGQ = 2)
               
x.12 <- glmer(y ~  logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.13 <- glmer(y ~ rich +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.14 <- glmer(y ~  nn2 + 
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.15 <- glmer(y ~  fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.16 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
               rich*logAbun +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 2)

x.17 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
               rich*fncDay +
               (1 | site), data = bff1,
             family = binomial, nAGQ = 5)

x.18 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

# because x.18 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.18.1 <- glmer(y ~ rich + nn2 +        logAbun +
                  rich*nn2 +
                  (1 | site), data = bff1,
                family = binomial, nAGQ = 2)

x.18.2 <- glmer(y ~ rich + nn2 + fncDay  +
                  rich*nn2 +
                  (1 | site), data = bff1,
                family = binomial, nAGQ = 5)

x.18.3 <- glmer(y ~ rich + nn2 +  
                  rich*nn2 +
                  (1 | site), data = bff1,
                family = binomial, nAGQ = 2)

x.19 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.20 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 5)

x.21 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun+
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

# because x.21 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.21.1 <- glmer(y ~        nn2 + fncDay + logAbun +
                fncDay*logAbun+
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.21.2 <- glmer(y ~ rich       + fncDay + logAbun +
                fncDay*logAbun+
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.21.3 <- glmer(y ~             fncDay + logAbun +
                fncDay*logAbun+
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.22 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 3)

x.23 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + nn2*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 3)

x.24 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + nn2*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.25 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 + rich*fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 0)


x.26 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                rich*nn2 + rich*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

# because x.26 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.26.1 <- glmer(y ~ rich + nn2  + logAbun +
                rich*nn2 + rich*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.27 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                rich*fncDay + rich*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 0)

x.28 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*nn2 +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)
 
x.29 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 1)

x.30 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                fncDay*logAbun + rich*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 0)

x.31 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay +  rich*nn2 + 
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

# because x.31 is ranked relatively highly by AIC, lets add a few models
# to try and determine which predictors are most important
x.31.1 <- glmer(y ~ rich + nn2 + fncDay +
                  nn2*fncDay +  rich*nn2 +
                  (1 | site), data = bff1,
                family = binomial, nAGQ = 2)
x.31.2 <- glmer(y ~ rich + nn2 + fncDay +
                  rich*nn2 +
                  (1 | site), data = bff1,
                family = binomial, nAGQ = 2)

x.32 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
               nn2*fncDay + rich*fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 1)

x.33 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*fncDay + rich*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 0)

x.34 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*nn2 +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 0)

x.35 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*fncDay +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 2)

x.36 <- glmer(y ~ rich + nn2 + fncDay + logAbun +
                nn2*logAbun + rich*logAbun +
                (1 | site), data = bff1,
              family = binomial, nAGQ = 0)

x.37 <- glmer(y ~ (1 | site), data = bff1,
              family = binomial, nAGQ = 0)
# compare models                                                     ####

g <- model.sel(x.1,x.2,x.3,x.4,x.5,x.6,x.7, x.8, x.9, x.10, x.11, x.12,
            x.13, x.14,x.15,x.16, x.17, x.18, x.19, x.20, x.21, x.22,
            x.23, x.24, x.25, x.26, x.27, x.28, x.29, x.30, x.31,
            x.32, x.33, x.34, x.35, x.36,
            x.18.1, x.18.2, x.18.3,
            x.21.1, x.21.2, x.21.3,
            x.26.1,
            x.31.1, x.31.2, x.37)


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
# Make model selection table #
# as.data.frame(res1) %>%
#   select(model, df, delta, weight) %>%
#   write.csv(file = '/Users/learichardson/Desktop/bff.csv')
# 
# 
# 
sw(subset(g, delta <= 2))
avg <- model.avg(g, subset = delta < 2)
avgmod.95p <- model.avg(g, cumsum(weight) <= .95)
confint(avg)



# model diagnostics                                                  ####
# some diagnostics, I followed the procedure here:
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

simulationOutput <- simulateResiduals(fittedModel = x.28, n =1000)
plot(simulationOutput)
testResiduals(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = x.31.1, n =1000)
plot(simulationOutput)
testResiduals(simulationOutput)

# things look good 


# visualize predicted values                                         ####



## visualize x.31.1: the relationship between nn2, richness, and pollinator visitation ####

# look at predicted pollinator visitation across the gradient of nn2 and rich
# nn2 <- seq(0, 11.1, 0.07)
# rich <- seq(-0.2, 10, 0.07)
# logAbun <- 
# makeLook <- expand.grid(nn2, rich, unique(bff$site))
# names(makeLook) <- c('nn2' , 'rich', 'site')
# makeLook$fncDay <- as.factor(3)
# 
# makeLook %>%
#   mutate(rich = rich - mean(bff$rich)) %>%
#   mutate(logAbun = logAbun - mean(bff$logAbun)) %>%
#   mutate(nn2 = nn2 - mean(bff$nn2))
# contrasts(bff1$fncDay) <- contr.sum(3)
# 
# fit <- predict(x.31.1, newdata = makeLook, type="response")
# makeLook$fit <- fit
# 
# 
# makeLook$colors <- gray(makeLook$fit)
# 
# ggplot(makeLook, aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point(color = makeLook$colors)
# 
# makeLook$fncDay <- as.factor(2)
# fit <- predict(x.28.1, newdata = makeLook, type="response")
# makeLook$fit <- fit
# 
# ggplot(makeLook, aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point(color = makeLook$colors)
# 
# makeLook$fncDay <- as.factor(3)
# fit <- predict(x.28.1, newdata = makeLook, type="response")
# makeLook$fit <- fit
# 
# ggplot(makeLook, aes(nn2, rich))+
#   facet_wrap(~site)+
#   geom_point(color = makeLook$colors)


# all of thes plots show that despite subtle differences by day or by site, 
# the most obvious pattern is the one we already identified in the 
# previous analysis, that isolated Echinacea surrounded by high richness 
# and Echinacea with nearby mates surrounded by low richness have the 
# highest probability of pollinator visitation

# make a plot showing relationship between nn2 and rich faceted by day
# Make Figure 2 ####
nn2 <- seq(0, 11.1, .2)
rich <- seq(-0.2, 10, .2)
logAbun <- seq(0,8.3, .2)
fncDay <- as.factor(c(1, 2,3))
makeLook <- expand.grid(nn2, rich, unique(ff$site), logAbun, fncDay)
names(makeLook) <- c('nn2' , 'rich', 'site', 'logAbun', 'fncDay')
makeLook <- makeLook %>%
mutate(rich = rich - mean(bff$rich)) %>%
  mutate(logAbun = logAbun - mean(bff$logAbun)) %>%
  mutate(nn2 = nn2 - mean(bff$nn2))
contrasts(makeLook$fncDay) <- contr.sum(3)
fit <- predict(x.31.1, newdata = makeLook, type="response")
makeLook$fit <- fit


makeLook %>%
  mutate(rich = rich + mean(bff$rich, na.rm = T)) %>%
  mutate(logAbun = logAbun + mean(bff$logAbun)) %>%
  mutate(nn2 = nn2 + mean(bff$nn2)) %>%
  ggplot(aes(nn2, rich)) +
  geom_raster(aes(fill = fit), interpolate = T) +
  theme_classic() +
  geom_jitter(data = bff, pch = 21, 
              col = ifelse(bff$polObserved == '1', 'white', 'black'), 
              aes(fill = polObserved)) +
  scale_fill_gradient(low = 'white', high = 'black') +
  labs(x = 'Distance (m) to second nearest flowering conspecific',
       y = 'Floral species richness',
       fill = 'Probability\\nof\\npollinator\\nvisitation')+
  scale_y_continuous(breaks=c(0,3,6,9))

