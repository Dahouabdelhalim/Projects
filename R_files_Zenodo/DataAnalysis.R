rm(list=ls())

library(tidyverse)
library(readxl)
library(lubridate)
library(lme4)
library(unmarked)
library(pbkrtest)
library(cowplot)
library(AICcmodavg)


# Importing data ----------------------------------------------------------


#Set the working directory
setwd("...")

#Read in the spreadsheets
stand = read_excel("Data.xlsx", sheet='stand')
site = read_excel("Data.xlsx", sheet='site')
dawnSurvey = read_excel("Data.xlsx", sheet='dawnSurvey') %>% 
  mutate(dateSurveyed = as_date(dateSurveyed))
prdRecording = read_excel("Data.xlsx", sheet='prdRecording') %>% 
  mutate(mamuDetectCount = as.integer(mamuDetectCount))


# Analyzing recording data 2016 -----------------------------------------------------


modelData = stand %>% 
  select(standName, treatment, historicStatus) %>% 
  full_join(site, by='standName') %>% 
  select(standName, siteNumber, treatment, historicStatus) %>% 
  full_join(prdRecording, by='siteNumber') %>% 
  filter(date < as.Date('2017-01-01'))%>% #A few haphazard recordings were made in 2017
  group_by(standName, treatment, date, historicStatus) %>% 
  summarise(calls = sum(mamuDetectCount, na.rm=T),
            minutes = sum(lengthReviewedNearestMinute, na.rm=T)) %>% 
  ungroup() %>% 
  arrange(minutes) %>% 
  filter(minutes >= 100) %>% 
  mutate(week = cut(date, breaks = 'week', start.on.monday=F),
         detected = ifelse(calls > 0, 1, 0),
         doy = scale(yday(as.Date(date)), center=T, scale=T),
         doy2 = doy^2)


#What proportion of days did we detect MAMU at treatment and control stands?
modelData %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(detected))

#Fitting 3 models
full = glmer(detected ~ treatment + doy + doy2 +
               treatment*doy + treatment*doy2 + (1|standName),
             data=modelData, family=binomial(link='logit'))

noInteraction = glmer(detected ~ treatment + doy + doy2 + (1|standName),
                      data=modelData, family=binomial(link='logit'))

reduced = glmer(detected ~ doy + doy2 + (1|standName),
                data=modelData, family=binomial(link='logit'))

#Comparing the models using AICc
aictab(list('full'=full, 'noInteraction'=noInteraction, 'null'=reduced))


#Does this full model fit?  To assess this we are simulating datasets
#from the fitted model and then fitting a model to those.  We record
#the chisquare at each and calculate the proportion of the time
#this value is > the chisquare value for the actual data.

nSims = 10  #We simulated 500 for the manuscript, but that takes a while
simResults = rep(NA, nSims)

tmp1 = simulate(full, nsim=nSims)
for(i in 1:nSims){
  tmp2 = modelData %>%
    mutate(detected = tmp1[,i])
  tmp3 = glmer(detected ~ treatment + doy + doy2 +
          treatment*doy + treatment*doy2 + (1|standName),
        data=tmp2, family=binomial(link='logit'))
  o = tmp2$detected
  e = predict(tmp3, type='response')
  simResults[i] = sum((o-e)^2/e, na.rm=T)
}

simResults = data.frame(simulation = 1:nSims, chiSquare = simResults)


modChisq = sum((modelData$detected - predict(full, type='response'))^2/predict(full, type='response'))

#Plotting the distribution of chiSquare values
ggplot(simResults, aes(x=chiSquare))+
  geom_histogram(fill='white', color='black')+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_vline(xintercept=modChisq, color='red')

sum(simResults[,2] > modChisq)/length(simResults[,2])

#A p-value in the upper or lower 0.025 quantile would indicate a lack of fit.



#Bootstrapping confidence intervals around expected odds and probabilities of recording
#murrelets at treatment and control stands
tmp = data.frame('doyUnscale' = seq(min(yday(modelData$date)), max(yday(modelData$date)), 5)) %>% 
  mutate(doy = (doyUnscale-mean(yday(modelData$date)))/sd(yday(modelData$date)))
newData = tmp %>% 
  mutate(treatment='Playback') %>% 
  rbind(tmp %>% 
          mutate(treatment='Control')) %>% 
  mutate(doy2 = doy^2,
         date = as.Date(doyUnscale, origin=as.Date('2016-01-01')))
mm = model.matrix(~treatment+doy+doy2+treatment*doy+treatment*doy2, newData)
newData$odds = exp(mm%*%fixef(full))
newData$detected = plogis(mm%*%fixef(full))

predFun = function(.) mm %*% fixef(.)
bb = bootMer(full, FUN=predFun, nsim=10) #Used 500 for the manuscript
newData$lclOdds = exp(apply(bb$t, 2, function(x) quantile(x, probs=c(0.025, 0.975))))[1,]
newData$uclOdds = exp(apply(bb$t, 2, function(x) quantile(x, probs=c(0.025, 0.975))))[2,]
newData$lclDet = plogis(apply(bb$t, 2, function(x) quantile(x, probs=c(0.025, 0.975))))[1,]
newData$uclDet = plogis(apply(bb$t, 2, function(x) quantile(x, probs=c(0.025, 0.975))))[2,]
newData = newData %>% 
  mutate(treatment = ifelse(treatment=='Playback', 'Treatment', 'Control'))

#newData will be used to look at predicted values for treatment and control stands.
#newData2 will be used to look at the difference between treatment and control stands.
newData2 = newData %>% 
  filter(treatment=='Treatment') %>% 
  select(doyUnscale, doy, doy2, date, odds, detected) %>% 
  rename('oddsTreat'=odds, 'detectedTreat'=detected) %>% 
  left_join(newData %>% 
              filter(treatment=='Control') %>% 
              select(doyUnscale, doy, doy2, date, odds, detected) %>% 
              rename('oddsControl'=odds, 'detectedControl'=detected), by=c('doyUnscale', 'doy', 'doy2', 'date')) %>% 
  mutate(oddsRatio = oddsTreat/oddsControl)

newData2$oddsRatioLcl = apply(exp(bb$t[,which(newData$treatment=='Treatment')])/exp(bb$t[,which(newData$treatment=='Control')]),2, function(x) quantile(x, probs=c(0.025, 0.975)))[1,]
newData2$oddsRatioUcl = apply(exp(bb$t[,which(newData$treatment=='Treatment')])/exp(bb$t[,which(newData$treatment=='Control')]),2, function(x) quantile(x, probs=c(0.025, 0.975)))[2,]
newData2$diff = apply(plogis(bb$t[,which(newData$treatment=='Treatment')]) - plogis(bb$t[,which(newData$treatment=='Control')]), 2, function(x) mean(x))
newData2$diffLcl = apply(plogis(bb$t[,which(newData$treatment=='Treatment')]) - plogis(bb$t[,which(newData$treatment=='Control')]),2, function(x) quantile(x, probs=c(0.025, 0.975)))[1,]
newData2$diffUcl = apply(plogis(bb$t[,which(newData$treatment=='Treatment')]) - plogis(bb$t[,which(newData$treatment=='Control')]),2, function(x) quantile(x, probs=c(0.025, 0.975)))[2,]

#Odds ratio (treatment/control) over time
newData2 %>% select(date, oddsRatio, oddsRatioLcl, oddsRatioUcl)

#Creating Figure 2 Building a figure to show expected values with confidence intervals

a = ggplot()+
  geom_ribbon(data=newData, aes(x=date, y=detected, color=treatment, fill=treatment,ymin=lclDet, ymax=uclDet), alpha=0.3, linetype='blank')+
  geom_line(data=newData, aes(x=date, y=detected, color=treatment), linetype='solid', size=2)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  ylab('Probability of recording\\na murrelet')+
  theme(legend.title=element_blank())+
  theme(legend.position=c(0.2, 0.8))+
  guides(color=F)+
  scale_color_brewer(palette='Set2')+
  scale_fill_brewer(palette='Set2')+
  xlab('Date')+
  theme(axis.title.y=element_text(size=10))

b = ggplot()+
  geom_ribbon(data=newData2, aes(x=date, y=diff, ymin=diffLcl, ymax=diffUcl), alpha=0.3, fill='#9999CC')+
  geom_line(data=newData2, aes(x=date, y=diff), size=2, color='#9999CC')+
  geom_abline(intercept=0, slope=0, color='black', linetype='dashed')+
  theme_bw()+
  theme(panel.grid=element_blank())+
  xlab('Date')+
  ylab('Difference in probability\\n(treatment - control)')+
  theme(axis.title.y=element_text(size=10))

plot_grid(a, b, ncol=1, labels=c('A', 'B'))



# Dawn survey data 2017 ---------------------------------------------------


#Looking at the raw data
dawnSurvey %>% 
  mutate(occupied = ifelse(outcome=='Occupied', 1, 0),
         presOrOccupied = ifelse(outcome=='No detections', 0, 1)) %>% 
  group_by(standName) %>% 
  summarise(occupied = max(occupied, na.rm=T),
            presOrOccupied = max(presOrOccupied, na.rm=T)) %>% 
  ungroup() %>% 
  left_join(stand %>% 
              select(standName, treatment), by='standName') %>% 
  group_by(treatment) %>% 
  summarise(totalstands = n(),
            numOccupied = sum(occupied, na.rm=T),
            propOccupied = mean(occupied, na.rm=T),
            numPresOrOccupied = sum(presOrOccupied, na.rm=T),
            propPresOrOccupied = mean(presOrOccupied, na.rm=T))


#### Presence dataset
y = dawnSurvey %>% 
  mutate(outcome = ifelse(outcome=='No detections', 0, 1)) %>% 
  select(standName, outcome, dateSurveyed) %>% 
  group_by(standName) %>% 
  mutate(surveyNumber = rank(dateSurveyed)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=standName, names_from=surveyNumber, values_from=outcome) %>% 
  arrange(standName) %>% 
  select(-standName)

standCovars = dawnSurvey %>% 
  select(standName) %>% 
  unique() %>% 
  left_join(stand %>%
              select(standName, treatment, historicStatus), by='standName') %>% 
  arrange(standName) %>% 
  mutate(treatment = ifelse(treatment=='Playback', 'Treatment', 'Control'))

obsCovars = dawnSurvey %>% 
  select(standName, dateSurveyed, observer) %>% 
  group_by(standName) %>% 
  mutate(surveyNumber=rank(dateSurveyed)) %>% 
  ungroup() %>% 
  full_join(expand.grid('standName' = unique(dawnSurvey$standName),
                        'surveyNumber' = 1:ncol(y)),
            by=c('standName', 'surveyNumber')) %>% 
  mutate(doy = scale(yday(dateSurveyed)), center=T, scale=T) %>% 
  mutate(doy2 = doy^2) %>% 
  arrange(standName, surveyNumber)
presenceData = unmarkedFrameOccu(y=y, siteCovs=standCovars, obsCovs=obsCovars)

#Fitting full and reduced occupancy models
presFull = occu(~doy+doy2~treatment, data=presenceData, method='BFGS')
presReduced = occu(~doy+doy2~1, data=presenceData, method='BFGS')

#Comparing models using AICc
aictab(list('reduced'=presReduced,
            'full'=presFull))


#Checking model fit
chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  #sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE) #This was from the Unmarked vignette, but it doesn't seem right to me. (chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://cran.r-project.org/web/packages/unmarked/vignettes/unmarked.pdf)
  sum((y-fv)^2/fv, na.rm=T)
}

#Model is a good fit if t0 does not fall in the upper or lower 0.025 quantile
(pb <- parboot(presReduced, statistic=chisq, nsim=500))

#Estimated odds ratio (treatment/control) from full model
exp(coef(presFull)[2]) #Odds ratio
exp(coef(presFull)[2] - 1.96*SE(presFull)[2]) #Odds ratio LCL
exp(coef(presFull)[2] + 1.96*SE(presFull)[2]) #Odds ratio UCL



#### Occupancy dataset
y = dawnSurvey %>% 
  mutate(outcome = ifelse(outcome=='Occupied', 1, 0)) %>% 
  select(standName, outcome, dateSurveyed) %>% 
  group_by(standName) %>% 
  mutate(surveyNumber = rank(dateSurveyed)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=standName, names_from=surveyNumber, values_from=outcome) %>% 
  arrange(standName) %>% 
  select(-standName)

dawnData = unmarkedFrameOccu(y=y, siteCovs=standCovars, obsCovs=obsCovars)

occFull = occu(~doy+doy2~treatment, data=dawnData, method='BFGS')
occReduced = occu(~doy+doy2~1, data=dawnData, method='BFGS')

#Comparing models using AICc
aictab(list('reduced'=occReduced,
            'full'=occFull))

#Model is a good fit if t0 does not fall in the upper or lower 0.025 quantile
(pb <- parboot(occFull, statistic=chisq, nsim=500)) #Model fits

#Estimated odds ratio (treatment/control) from full model
exp(coef(occFull)[2]) #Odds ratio
exp(coef(occFull)[2] - 1.96*SE(occFull)[2]) #Odds ratio LCL
exp(coef(occFull)[2] + 1.96*SE(occFull)[2]) #Odds ratio UCL


#Creating plots from these models

tmp = data.frame(treatment = c('Control', 'Treatment'))

results = predict(presFull, newdata = tmp, type = 'state') %>% 
  mutate(model = 'Presence', type = 'Presence', treatment=c('Control', 'Treatment')) %>% 
  rbind(predict(occFull, newdata = tmp, type = 'state') %>% 
          mutate(model = 'Breeding', type = 'Breeding', treatment=c('Control', 'Treatment'))) %>% 
  mutate(model = ifelse(model=='Presence', 'Presence', 'Occupancy')) %>% 
  mutate(model = factor(as.character(model), levels=c('Presence', 'Occupancy')))

a = ggplot(results, aes(x=model, y=Predicted))+
  geom_errorbar(aes(ymin=lower, ymax=upper, color=treatment), position=position_dodge(width=0.5), size=1)+
  geom_point(aes(color=treatment), position=position_dodge(width=0.5), size=4)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  ylim(0,1.1)+
  ylab('Estimated probability')+
  xlab('State variable')+
  theme(legend.title = element_blank())+
  scale_color_brewer(palette='Set2')+
  theme(axis.title = element_text(size=10))+
  theme(legend.position='bottom')+
  theme(plot.margin=unit(c(0,0,0,0), 'cm'))+
  theme(axis.text = element_text(size=8))

tmp = data.frame(date = seq(min(dawnSurvey$dateSurveyed), max(dawnSurvey$dateSurveyed),1)) %>% 
  mutate(doy = (yday(date) - mean(yday(dawnSurvey$dateSurveyed)))/sd(yday(dawnSurvey$dateSurveyed))) %>% 
  mutate(doy2 = doy^2)

results = predict(presFull, newdata = tmp, type='det', append=T) %>% 
  mutate(model = 'Presence', type = 'Detection') %>% 
  rbind(predict(occFull, newdata = tmp, type='det', append=T) %>% 
          mutate(model = 'Occupancy', type = 'Detection')) %>% 
  mutate(model = factor(as.character(model), levels=c('Presence', 'Occupancy')))


b = ggplot(results, aes(x=date, y=Predicted, linetype=model))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill='#9999CC', show.legend=F)+
  geom_line(size=2, color='#9999CC')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab('Date')+
  ylab('Estimated detection probability')+
  theme(legend.title=element_blank())+
  theme(axis.title = element_text(size=10))+
  theme(legend.position='bottom')+
  scale_linetype_manual(values = c('solid', '11'))+
  theme(plot.margin=unit(c(0,0,0,0), 'cm'))+
  theme(axis.text = element_text(size=8))+
  scale_x_date(limits = as.Date(c('2017-06-20','2017-08-02')))+
  theme(legend.key.width = unit(2, "cm"))


plot_grid(get_legend(a),
               a+theme(legend.position='none'),
               NULL,
               get_legend(b),
               b+theme(legend.position='none'),
               ncol=1, labels=c('A', '', '', 'B', ''), rel_heights=c(0.1,1,0.05,0.1,1))

