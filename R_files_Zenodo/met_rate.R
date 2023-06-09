# ---- GOAL: METABOLIC RATE DIFFERENCES 
# ---- install and load packages ----
library("lme4") #LM, LMM, GLM, GLMM
library("dplyr")
library("car") #anova for LM and LMM
library("multcomp")

library("Rmisc") #stats (not correlation)
library("Hmisc") #correlation significance

library("ggplot2") #for all graphs
library("ggpubr") #oublication ready graphs (also stat_cor for ggplot)
library("ggfortify") #diagnostic plots of models

library("corrplot") #ggcorr plots
library("RColorBrewer")
library("reshape2")

# loading Java everytime.
Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jdk-14")
system("java -version", intern = TRUE)
library("rJava") #glmulti needs rJava

# model selection
library("metafor") #rma function
library("glmulti") 
library("MuMIn") #multi-model inference


# ---- working directories shortcut ----
raw.dir <- ("G:/My Drive/R/met-rate/data-input") #raw datasheet
out.sheet <- ("G:/My Drive/R/met-rate/data-output") # datasheets after editing
gr.col <- ("G:/My Drive/R/met-rate/graphs/whole-colony") # 208 and 2019 dataset
gr.brd <- ("G:/My Drive/R/met-rate/graphs/brood-only") # 2019 brood
gr.adlt <- ("G:/My Drive/R/met-rate/graphs/adults-only") #2019 adults 
gr.data <- ("G:/My Drive/R/met-rate/graphs/data-aspect") #raw data
gr.qn <- ("G:/My Drive/R/met-rate/graphs/queens-only") #2019 groups of queens


# ---- loading raw per file data and assigning features----
setwd(raw.dir)
list.files() #visualize files in the folder

# loading compiled metabolic rate CSV file and adding variables
# whole colony measurement
df.col <- read.csv("mr-colony.csv")
df.col <- df.col %>%
  mutate(net.co2 = co2.mean - trash.co2.mean, #ppm. Identifying CO2 from only colony
         mass.diff= mass.before - mass.after, #grams. Change in weight during the experiment
         co2 = (net.co2*air.flow*60)/1000000, #1ppm = 10^-6 g/ml, airflow = ml/min, 60 minute = 1 hour. CO2 flow in grams per hour from whole colony
         co2.emission.rate = co2/mass.before, # ml/g*hr. CO2 emission rate per gram of colony
         o2.consmp = co2.emission.rate/rq, # ml/g*hr. #Oxygen consumption rate per gm of colony. Species specific respiratory quotient to convert CO2 emission rate to O2 consumption rate.
         oxy.joule = (16 + (5.164*rq)), #Energy rate. mathematical formula to calculate joule heat energy produced per unit of O2 consumed
         met.rate.j = (o2.consmp*oxy.joule)/3600, #Joules/sec at measured temperature per gram of colony
         met.rate.w = met.rate.j*1000000,#microwatts at measured temperature from whole colony
         met.rate.25 = met.rate.w/(2^((temp.mean-25)/10)), #micro watts per gram of colony at 25'C
         met.rate.25.col = met.rate.25*mass.before, #micro watts for whole colony at 25'C not weighted to colony's mass
         young.brood = eggs + young.larvae,
         old.brood = old.larvae + worker.pupae,
         total.brood = young.brood + old.brood,
         adults = workers + queens,
         colony.size = adults + total.brood,
         brood.prp = (total.brood/colony.size),#proportion of brood per colony
         yb.prp = (young.brood/colony.size), #proportion of eggs and young larvae per colony
         old.prp = (old.brood/colony.size), #proportion of late-instar larvae and worker pupae per colony
         adlt.prp = (adults/colony.size), #proportion of queens and worker per colony
         humid.mean.col = humid.mean - trash.humid.mean) #humidity due to background

# brood measurement
df.brd <- read.csv("mr-brood.csv")
df.brd <- df.brd %>%
  mutate(net.co2 = co2.mean - trash.co2.mean, #ppm. Identifying CO2 from only colony
         mass.diff= mass.before - mass.after, #grams. Change in weight during the experiment
         co2 = (net.co2*air.flow*60)/1000000, #1ppm = 10^-6 g/ml, airflow = ml/min, 60 minute = 1 hour. CO2 flow in grams per hour from whole colony
         co2.emission.rate = co2/mass.before, # ml/g*hr. CO2 emission rate per gram of colony
         o2.consmp = co2.emission.rate/rq, # ml/g*hr. #Oxygen consumption rate per gm of colony. Species specific respiratory quotient to convert CO2 emission rate to O2 consumption rate.
         oxy.joule = (16 + (5.164*rq)), #Energy rate. mathematical formula to calculate joule heat energy produced per unit of O2 consumed
         met.rate.j = (o2.consmp*oxy.joule)/3600, #Joules/sec at measured temperature per gram of colony
         met.rate.w = met.rate.j*1000000,#microwatts at measured temperature from whole colony
         met.rate.25 = met.rate.w/(2^((temp.mean-25)/10)), #micro watts per gram of colony at 25'C
         met.rate.25.brd = met.rate.25*mass.before, #micro watts for whole brood (not weighted by mass)
         young.brood = eggs + young.larvae,
         old.brood = old.larvae + worker.pupae,
         total.brood = young.brood + old.brood,
         yb.prp = (young.brood/total.brood),
         old.prp = (old.brood/total.brood)) #proportion of brood per colony

# brood measurement
df.qn <- read.csv("mr-queens.csv")
df.qn <- df.qn %>%
  mutate(co2 = (co2.mean*air.flow*60)/1000000, #1ppm = 10^-6 g/ml, airflow = ml/min, 60 minute = 1 hour. CO2 flow in grams per hour from whole colony
         co2.emission.rate = co2/mass, # ml/g*hr. CO2 emission rate per gram of colony
         o2.consmp = co2.emission.rate/rq, # ml/g*hr. #Oxygen consumption rate per gm of colony. Species specific respiratory quotient to convert CO2 emission rate to O2 consumption rate.
         oxy.joule = (16 + (5.164*rq)), #Energy rate. mathematical formula to calculate joule heat energy produced per unit of O2 consumed
         met.rate.j = (o2.consmp*oxy.joule)/3600, #Joules/sec at measured temperature per gram of colony
         met.rate.w = met.rate.j*1000000,#microwatts at measured temperature from whole colony
         met.rate.25 = met.rate.w/(2^((temp-25)/10)), #micro watts per gram of colony at 25'C
         met.rate.25.qn = met.rate.25*mass) #micro watts for whole brood (not weighted by mass)

# ---- loading processed files ----
setwd(out.sheet)

df.col <- read.csv ("metrate-colony.csv") #metabolic rates of whole colony 2018, 2019
df.brd <- read.csv("mr2019-brood.csv") #met rate of only brood 2019
df.adlt <- read.csv("mr2019-adult.csv") #adults2019 = colony2019-brood2019
df.qn <- read.csv("mr2019-queens.csv") #queen only
df.hmd <- read.csv("mr2019-humidity-colony.csv") #queen only

df.col2019 <- data.frame(subset(df.col, year == "2019"))
df.col2018 <- data.frame(subset(df.col, year == "2018"))

# ---- Importing video analysis data ----
setwd(raw.dir)

df.vid <- read.csv("met-rate-vid-activity.csv")

# merging average activity per colony witl metabolic rate data
#finding common columns
common_col_names <- intersect(names(df.col), names(df.vid))

#combining two dataframes with common column names
df.col.vid <- merge(df.col, df.vid, by=common_col_names, all.x=TRUE)

# ---- Creating dataset as per requirement----
# removing colony # 7 (uninfected, 2019) since it had high humidity during measurement due to leaky respirometer chamber for whole colony recording
df.col.vid <- data.frame(subset(df.col.vid, humid.mean < "9" ))
df.brd <- data.frame(subset(df.brd, colony != "7"))

# Humidity: same colony has been measured with and without humidity
df.hmd <- data.frame(subset(df.col.vid, year == "2019")) #testing effect of humidity on measurement
df.hmd <- df.hmd %>% 
  group_by(colony) %>% 
  filter(n() == 2 )

# updating datasheet to exclude colonies in 2019 that have been measured twice
df.col.vid <- data.frame(subset(df.col.vid,humidity == "yes" | year == "2018"))

df.col2019 <- data.frame(subset(df.col.vid, year == "2019"))
df.col2018 <- data.frame(subset(df.col.vid, year == "2018"))

# Adults(= whole colony - brood)
df.adlt <- df.col2019 %>%
  filter(colony %in% df.brd$colony) #same colonies analyzed adults and brood

# merge two dataframes.some values duplicated. x = whole colony values, y = only brood values
df.adlt <- merge(df.col2019,df.brd,
                  by=c("colony", "wolbachia", "queen.age"))

df.adlt$net.co2 = df.adlt$net.co2.x - df.adlt$net.co2.y
df.adlt$mass = df.adlt$mass.before.x - df.adlt$mass.before.y

df.adlt <- df.adlt %>%
  mutate(co2 = (net.co2*air.flow.x*60)/1000000, #1ppm = 10^-6 g/ml, airflow = ml/min, 60 minute = 1 hour. CO2 flow in grams per hour from whole colony
         co2.emission.rate = co2/mass, # ml/g*hr. CO2 emission rate per gram of colony
         o2.consmp = co2.emission.rate/rq.x,# ml/g*hr. #Oxygen consumption rate per gm of colony. Species specific respiratory quotient to convert CO2 emission rate to O2 consumption rate.
         oxy.joule = (16 + (5.164*rq.x)),
         met.rate.j = (o2.consmp*oxy.joule)/3600, #Joules/sec at measured temperature per gram of colony
         met.rate.w = met.rate.j*1000000,#microwatts at measured temperature from whole colony
         mr.25 = met.rate.w/(2^((temp.mean.x-25)/10)),#metabolic rate of adults per gram of adults
         mr.25.adlt = mr.25*mass) #metabolic rate from all adults, not weighted by mass

#removing columns to make adult only dataset.
df.adlt <- data.frame(subset(df.adlt, select = c(colony:queen.age,
                                                 queens:rq.x,
                                                 humid.mean.x,
                                                 temp.mean.x,
                                                 adults,
                                                 net.co2:mr.25.adlt)))
df.adlt <- df.adlt %>%
  rename(air.flow = air.flow.x,
         rq =rq.x,
         humid.mean = humid.mean.x,
         temp.mean = temp.mean.x)

# ---- saving this file in output directory -----
setwd(out.sheet)

write.csv(df.col.vid, "metrate-colony.csv") #all colonies 2019 and 2018 with average change in pixel from video
write.csv(df.brd, "mr2019-brood.csv") #2019 brood only data
write.csv(df.adlt, "mr2019-adult.csv") #2019 adult only data
write.csv(df.hmd, "mr2019-humidity-colony.csv") #2019 humidity comparisons of the same colony measurement
write.csv(df.qn, "mr2019-queens.csv") #2019 humidity comparisons of the same colony measurement

# ---- Sample size -----
aggregate(colony ~ wolbachia, df.col2019, 
          function(x) length(unique(x))) 

#sample size for whole dataset
nrow(df.brd)

# ---- Normality of raw data ####################
#Shapiro Wilks Test for sampling from normal distribution, p < 0.05: not drawn from normal distribution (use GLM
shapiro.test(df.col2019$met.rate.25)
shapiro.test(df.col2019$met.rate.25.col)

shapiro.test(df.brd$met.rate.25)
shapiro.test(df.brd$met.rate.25.brd) 

shapiro.test(df.col2018$met.rate.25)
shapiro.test(df.col2018$met.rate.25.col) 

#Batlett's test for unequal variance. p-value < 0.05 : variance can not be assumed equal
bartlett.test(met.rate.25 ~ wolbachia, data=df.col) #p < 0.05. unequal variance
bartlett.test(adlt.mr.25 ~ wolbachia, data = df.adlt) #p < 0.05. unequal variance
bartlett.test(met.rate.25 ~ wolbachia, data = df.hmd) #p > 0.05. equal variance
bartlett.test(met.rate.25 ~ wolbachia, data = df.brd) #p > 0.05. equal variance

#F-test for equal variance, p < 0.05: variance not equal
var.test(met.rate.25 ~ wolbachia, data=df.col, alternative = "two.sided") #p < 0.05. unequal variance
var.test(adlt.mr.25 ~ wolbachia, data = df.adlt, alternative = "two.sided") #p < 0.05. unequal variance
var.test(met.rate.25 ~ wolbachia, data = df.hmd, alternative = "two.sided") #p > 0.05. equal variance
var.test(met.rate.25 ~ wolbachia, data = df.brd, alternative = "two.sided") #p > 0.05. equal variance


# ---- Model Selection ----
# goal: 

# designing a function which takes model and dataset as input
# fits a random/mixed-effects meta-regression model to the given data using maximum likelihood estimation
# source : http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti_and_mumin
rma.glmulti <- function(formula, data, ...)
  rma (formula, vi, data=data, method="ML", ...)

# creating a model with all the variables that we think can affect
#level = 1 (no interaction, fixed factors)
#level = 2 (interactions) - although a lot of model combinations since many possible combinations
col2019 <- glmulti(met.rate.25.col~mass.before+humid.mean+queen.age+eggs+young.larvae+old.larvae+worker.pupae+workers+avg.chg.pxl+wolbachia,
                   data=df.col2019,level=1,fitfunc=lm, crit="aicc", confsetsize=2000)

col2019 <- glmulti(met.rate.25~humid.mean+queen.age+brood.prp+avg.chg.pxl+wolbachia,
                   data=df.col2019,level=1,fitfunc=lm, crit="aicc", confsetsize=2000)

brd2019 <- glmulti(met.rate.25.brd~mass.before +humid.mean+eggs+young.larvae+queen.age+old.larvae+worker.pupae+wolbachia,
                   data=df.brd,level=1,fitfunc=lm, crit="aicc", confsetsize=2000)

col2018 <- glmulti(met.rate.25.col~mass.before+queen.age+eggs + young.larvae+old.larvae+worker.pupae+workers+wolbachia,
                   data=df.col2018,level=1,fitfunc=lm, crit="aicc", confsetsize=2000)

model <- col2018
print(model)
top <- weightable(model)
top <- top[top$aicc <= min(top$aicc) +2,] #top models based on lowest AIC
top
plot(model, type="s") #model-averaged importance

# ---------------- Functions for evaluating GLM and producing test statistics --------------------
# Goodness of fit for GLM ----
# Chi-square test on the residual deviance and degree of freedon
# ran GLM with Poisson, quasipoisson and negative binomial model
# null hypothesis = family (Poisson, quasipoisson or negative binomial) is an adquate fit for the data
# if P > 0.05, family specified in the model is a good fit, i.e., accept null hypothesis
p.gof <- function(model){
  dev <- deviance(model)
  rdf <- df.residual(model)
  pchisq(dev,
         df=rdf,
         lower.tail=FALSE)
}

# Overdispersion of GLM -----
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# if overdispersed then use quasipoisson or negative binomial
# selected between the 2 by looking at 3 criteria in order of preference -
#(a) linear fit of residuals (q-q plot), (b) goodness of fit of modelm & (c) lack of over dispersion

# Model diagnostics----
# Q-Q plot of model 
qq <- function(model){
  res <- residuals(model)
  qqnorm(res)
  qqline(res, probs=c(0.25,0.75))
}

# Histogam of residuals ---
res <- function(model) {
  res1 <- residuals(model)
  hist(res1)
}

# Normality of the residuals ---
shp.tst <- function(model) {
  res <- residuals(model)
  shapiro.test(res)
}
 
# Test statistics for model ----
# for quasi distribution and 
stats.glm <- function(model){
  drop1(model,.~.,test="F")
}

stats.lm <- function(model) {
  car::Anova(model, ddf="Satterthwaite")
}
# ---- LMER to identify effects #########
# understanding of model terms : https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R
#http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/#:~:text=Regression%20diagnostics%20plots%20can%20be,creates%20a%20ggplot2%2Dbased%20graphics.&text=The%20diagnostic%20plots%20show%20residuals,check%20the%20linear%20relationship%20assumptions.
# ---- 2018+2019: HYPOMETRIC SCALING AND OUTLIER TEST #########
# hypometric scaling of microwatts with mass of colony
lm.scale <- lm(log(met.rate.25.col) ~ log(mass.before),
               data = df.col)

model <- lm.scale
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
outlierTest(model) #row # of data that is an outlier
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# removing outlier row #24 and row #34
df.col <- subset(df.col, colony !="4-a")
df.col <- subset(df.col, colony !="15")

# removing these colonies from brood since these are from 2019
df.brd <- subset(df.brd, colony !="4-a")
df.brd <- subset(df.brd, colony !="15")

# removing these colonies from adult since these are from 2019
df.adlt <- subset(df.adlt, colony !="4-a")
df.adlt <- subset(df.adlt, colony !="15")

# saving csv file
write.csv(df.col, #data frame
          file.path(out.sheet,# file path
                    "metrate-colony.csv")) #2018 + 2019: metabolic rates, census and video activity for whole colony

write.csv(df.brd, #data frame
          file.path(out.sheet,# file path
                    "mr2019-brood.csv")) #2019: metabolic rates and census for brood

write.csv(df.adlt, #data frame
          file.path(out.sheet,# file path
                    "mr2019-adult.csv")) #2019: metabolic rates and census for adults


# ---- 2018 + 2019: analyzed together ----
#wolbachia-by-year interaction
lm <- lm (met.rate.25.col ~ wolbachia*mass.before,
               data = df.col)

lm.mass <- lm (mass.before ~ wolbachia,
               data = df.col)

model <- lm.mass
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# mass
lm.mass <- lm (met.rate.25.col ~ mass.before,
               data = df.col)

lm.mass <- lm (mass.before ~ queen.age,
               data = df.col)

model <- lm.mass
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# queen age
lm.qn.age <- lm (met.rate.25 ~ queen.age,
               data = df.col)

model <- lm.qn.age
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# humidity
lm.hmd <- lm (met.rate.25.col ~ humidity,
                 data = df.col)

model <- lm.hmd
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---- 2018: whole colony ----
df.col2018 <- subset(df.col, year == 2018)

# multivariate
lm.2018 <- lm (met.rate.25.col~wolbachia + mass.before + colony.size + avg.chg.pxl,
               data = df.col2018)

model <- lm.2018
vif(model) # multicollinearity of factors
outlierTest(model) #row # which is an outlier
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# univariate
lm.2018 <- lm (mass.before~ wolbachia,
               data = df.col2018)

lm.2018 <- lm (met.rate.25~ workers,
               data = df.col2018)

lm.2018 <- lm (met.rate.25.col~ workers,
               data = df.col2018)

model <- lm.2018
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# proportion of brood analysis
glm.brd <- glm (cbind(adults,(old.brood+young.brood))~ wolbachia, 
                data = df.col2018,
                family=quasibinomial(link=logit))

model <- glm.brd
summary(model) #summary
stats.glm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# differences in counts of individuals - used p.gof and overdisp_fun to assess goodness of fit of negative binomial and quasipoisson models
glm.brd <- glm.nb (colony.size ~ wolbachia, 
                  data = df.col2018)


model <- glm.brd
summary(model) #summary
stats.glm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---- 2019 : whole colony ---- 
df.col2019 <- subset(df.col, year == 2019)

# multivairate
lm.19 <- lm (met.rate.25.col ~ wolbachia + mass.before + queen.age +
               colony.size + avg.chg.pxl + humid.mean, data = df.col2019)
model <- lm.19
vif(model) # multicollinearity of factors
outlierTest(model) #row # which is an outlier
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# univariate
lm.2019 <- lm (met.rate.25 ~ humid.mean,
               data = df.col2019)

lm.2019 <- lm (met.rate.25.col~ humid.mean,
               data = df.col2019)

model <- lm.2019
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---- 2019: brood ----
# multivariate
lm.19b <- lm (met.rate.25.brd~wolbachia+ queen.age + mass.before + total.brood + humid.mean,
              data = df.brd)

model <- lm.19b
vif(model) # multicollinearity of factors
outlierTest(model) #row # which is an outlier
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# univariate
lm.19b <- lm (met.rate.25~ humid.mean,
               data = df.brd)

lm.19b <- lm (met.rate.25.brd~worker.pupae,
               data = df.brd)

model <- lm.19b
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# proportion of brood analysis
glm.brd <- glm (cbind(young.brood, old.brood)~ queen.age, 
                data = df.brd,
                family=quasibinomial(link=logit))

model <- glm.brd
summary(model) #summary
stats.glm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---- 2019: queens only ----
# Interaction included since multiple sampleID per queen age
lm.qn <- lm(met.rate.25.qn~queen.age*wolbachia + mass,
              data = df.qn)
model <- lm.qn
vif(model)
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# univariate
lm.qn <- lm(met.rate.25~wolbachia*queen.age,
            data = df.qn)

lm.qn <- lm (met.rate.25.qn~mass,data = df.qn)

model <- lm.qn
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# spearman rank correlation between met.rate and queen.age
cor.test(~ met.rate.25.qn + queen.age,
          data=df.qn,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95)

# ---- 2019: Adults -------------
# multivariate
lm.adlt <- lm(mr.25 ~ queen.age+wolbachia+ adlt.prp,
              data = df.adlt)

lm.adlt <- lm (mr.25.adlt ~ mass+ queen.age + wolbachia + adults ,data = df.adlt)

model <- lm.adlt
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# univariate
lm.adlt <- lm(mr.25~ queen.age,
              data = df.adlt)

lm.adlt <- lm (mr.25.adlt~mass,data = df.adlt)

model <- lm.adlt
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---- for colonies measured twice in 2019 to compare effects of humidity -----
#check disrtibution of data
hist(df.hmd$humid.mean,
     xlab = "mean humidity during measurement",
     ylab = "frequency",
     main = "Distribution of mean humidity values")

lm.hmd <- lm (met.rate.25.col ~ humidity + wolbachia,
             data = df.hmd) #what about paired comparisons?

model <- lm.hmd
summary(model) #summary
stats.lm(model) #ANOVA test statistic
qq(model) #qq plot of residuals
shp.tst(model) # test for normality of residuals
res(model) # histograph of residuals
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---- Correlation Matrix ----
# loading processed files from the analysis below
df.corr2018 <- read.csv("correlation-metrate2018.csv")
df.corr2019 <- read.csv("correlation-metrate2019.csv")

# step 1: shortlisting columns for calculating correlation
# removing avg.pxl.activity with non-NA values
df.col18 <- data.frame(subset(df.col2018, avg.chg.pxl != "is.NA",
                              select = c(queen.age:workers, brood.prp,
                                         humid.mean,mass.before,
                                         met.rate.25:adlt.prp,
                                         avg.chg.pxl)))

df.col19 <- data.frame(subset(df.col2019,select = c(queen.age:workers, brood.prp,
                                                    humid.mean,mass.before,
                                                    met.rate.25:adlt.prp,
                                                    avg.chg.pxl)))

df.brd19 <- data.frame(subset(df.brd,select = c(queen.age:worker.pupae,mass.before,
                                                    humid.mean, met.rate.25:old.prp)))

# Step 2: renaming column names to generic names for graphs
colnames(df.col18) <- c("Queen age","Eggs",
                        "Young instar larvae","Late-instar larvae","Worker pupae",
                        "Queens","Workers","Brood proportion","Humidity",
                        "Mass","Mass-specific metabolic rate",
                        "Metabolic rate","Total younger brood","Total older brood",
                        "Total brood","Adults","Colony size","Young brood proportion",
                        "Older brood proportion","Adult proportion","Activity")

colnames(df.brd19) <- c("Queen age","Eggs",
                        "Young instar larvae","Late-instar larvae","Worker pupae",
                        "Mass","Humidity","Mass-specific metabolic rate",
                        "Metabolic rate","Total younger brood","Total older brood",
                        "Total brood","Young brood proportion",
                        "Older brood proportion")


# Step 3: Spearman Correlation by year (normally distributed data)
# creating a list where correlation values are rounded to 3 decimal points
corr18 <- round(cor(df.col18,
                       method = "spearman"), 3)

corr19 <- round(cor(df.col19,
                       method = "spearman"), 3)

corr19.brd <- round(cor(df.brd19,
                        method = "spearman"),3)

# Step 3: P-values of Spearman Correlation stored in a separate lists upto 3 decimal point
corr18m <- rcorr(as.matrix(corr18),
                       type = "spearman") #storing all results from Spearman correlation
corr18.p <- round(corr18m[["P"]], 3) #extracting only P-values

corr19m <- rcorr(as.matrix(corr19),
                       type = "spearman") #storing all results from Spearman correlation
corr19.p <- round(corr19m[["P"]], 3) #extracting only P-values

corr19b <- rcorr(as.matrix(corr19.brd),
                 type = "spearman") #storing all results from Spearman correlation
corr19b.p <- round(corr19b[["P"]], 3) #extracting only P-values

# Step 3: Re-shaping the data into long format along with renaming columns
melt_corr18 <- melt(corr18)
melt_corr18 <- melt_corr18 %>%
  dplyr::rename(trait1 =Var1,
                trait2 =Var2,
                spearman.corr=value)

melt_corr18.p <- melt(corr18.p)
melt_corr18.p <- melt_corr18.p %>%
  dplyr::rename(trait1 =Var1,
                trait2 =Var2,
                p.value=value)

melt_corr19 <- melt(corr19)
melt_corr19 <- melt_corr19 %>%
  dplyr::rename(trait1 =Var1,
                trait2 =Var2,
                spearman.corr=value)

melt_corr19.p <-melt(corr19.p)
melt_corr19.p <- melt_corr19.p %>%
  dplyr::rename(trait1 =Var1,
                trait2 =Var2,
                p.value=value)

melt_corr19b <- melt(corr19)
melt_corr19b <- melt_corr19b %>%
  dplyr::rename(trait1 =Var1,
                trait2 =Var2,
                spearman.corr=value)

melt_corr19b.p <-melt(corr19b.p)
melt_corr19b.p <- melt_corr19b.p %>%
  dplyr::rename(trait1 =Var1,
                trait2 =Var2,
                p.value=value)

# Step 4: Combining Spearman correlation and P-values per year
#finding common columns
common_col_names18 <- intersect(names(melt_corr18), names(melt_corr18.p))
common_col_names19 <- intersect(names(melt_corr19), names(melt_corr19.p))
common_col_names_brd <- intersect(names(melt_corr19), names(melt_corr19.p))

#combining two dataframes with common column names
df.corr2018 <- merge(melt_corr18, melt_corr18.p,
                     by=common_col_names18, all.x=TRUE)

df.corr2019 <- merge(melt_corr19, melt_corr19.p,
                     by=common_col_names19, all.x=TRUE)

df.corr2019.brd <- merge(melt_corr19, melt_corr19.p,
                     by=common_col_names19, all.x=TRUE)

# dropping dataframes with P-value = NA
df.corr2018 <- subset(df.corr2018, p.value != "is.NA")
df.corr2019 <- subset(df.corr2019, p.value != "is.NA")
df.corr2019.brd <- subset(df.corr2019.brd, p.value != "is.NA")

# Step 5: save these files
write.csv(df.corr2018, file.path(out.sheet,"corr-mr2018.csv"))
write.csv(df.corr2019, file.path(out.sheet,"corr-mr2019-colony.csv"))
write.csv(df.corr2019.brd, file.path(out.sheet,"corr-mr2019-brood.csv"))

# ---- PLOT : Box plot with raw values ---------
yr.color <-c ("#52969A","#004549")
my.color <-c ("darkorchid4","chocolate2" )

#Step 1: create universal variable and 
data <- df.col2019
data$y <- data$met.rate.25.col
data$x <- data$wolbachia
y.lab <- "Metabolic Rates (microwatts)"
x.lab <- "Wolbachia Infection"
title <- "Metabolic Rates of \\nWhole Colonies"
#subtitle <- "(1- -month-old Queens)"

# Step 2: re-order
data$x <- factor(data$x, 
                 levels = c("infected", "uninfected"))

# Step 3: calculate mean
data.mean <- summarySE(data, measurevar="y",
                       groupvars=c("x"), na.rm = T)

data.summary <- data %>%
  dplyr::group_by(x) %>%
  dplyr::summarise(mean = mean(y),
                   median = median(y),
                   sd = sd(y),
                   sample.size = n())

# Step 4: plot
ggplot(data , aes(x=x,
                  y = y,
                  fill = x)) + 
  geom_boxplot(alpha = 0.3,
               width = 0.4)+
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_dotplot(binaxis = "y", #which axis to bin along
               binwidth = 28,
               stackdir = "center", #centered
               shape = 21,
               dotsize = 0.8,
               position = position_dodge(0.2),
               alpha = 0.6)+ 
  stat_summary(fun = mean, 
               geom = "point",
               shape = 24, 
               size = 3, 
               color = "black",
               fill = "black") +
  scale_fill_manual(values = my.color)+
  scale_color_manual(values = my.color)+
  #coord_cartesian(ylim = c(330, 970))+
  #scale_y_continuous(breaks = seq(330, 650, 970)) +
  labs(y = y.lab,
       title = title)+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial",
                                  colour = "black",
                                  face = "italic",
                                  size = "12"),
        panel.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "14"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "14"),
        axis.title.x=element_blank(), #removing X-axis title
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "14", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))


#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "mr2019-colony-nomass.png",
       path = gr.col,
       dpi = 300, 
       width = 3,
       height = 3.5, 
       units = "in")


# ---- PLOT : Correlation (Spearman Rank) plot ---------
my.color <-c ("chocolate2","darkorchid4")
yr.color <-c ("#52969A","#004549")

#Step 1: make universal variables
data <- df.col
data$y <- data$met.rate.25.col
data$x <- data$mass.before
data$yr.factor <- as.factor(data$year)
y.lab <- "metabolic rates \\n(microwatts)"
x.lab <- "weight of the colony (gram)"
title <- "Effect of Colony Mass on \\nMetabolic Rates of \\nWhole Colonies (2018+2019)"

#Step 2: plot
ggplot(data, aes(x=x,
                 y=y,
                 fill = yr.factor)) +
  scale_fill_manual(values  = yr.color) +
  geom_point(shape=21,
             size = 3,
             alpha = 0.8)+
  geom_smooth(method=lm,
              formula = y~x,
              se = T,
              fullrange = F,
              alpha = 0.2,
              fill = "#52969A",
              color ="#52969A" )+
  labs(title=title,
       x = x.lab,
       y = y.lab) +
  stat_cor(#label.y = 1700,
           method = "spearman")+ #position of R-square
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial",
                                  colour = "black",
                                  face = "italic",
                                  size = "12"),
        panel.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "16"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "14"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "16", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "mr20182019-mW-humid.png",
       path = gr.col,
       dpi = 300, 
       width = 5,
       height =4.5, 
       units = "in")

# ---- PLOT : Log-Log Scatter plot ---------
# Resource: https://statisticsbyjim.com/regression/log-log-plots/
# Resource: https://stackoverflow.com/questions/29586858/what-is-the-difference-between-scale-transformation-and-coordinate-system-trans
yr.color <-c ("#52969A","#004549")

#Step 1: make universal variables
data <- df.brd
data$y <- log(data$met.rate.25.brd)
data$x <- log(data$mass.before)
data$yr.factor <- as.factor(data$year)
x.lab <- "Log transformed mass (gram)"
y.lab <- "Log transformed metabolic \\nrates (microwatts)"
title <- "Metabolic Rate Scaling with \\nMass of Brood"

#Step 2: plot
ggplot(data, aes(x=x,
                 y=y)) +
  #scale_fill_manual(values  = yr.color) +
  geom_point(shape=21,
             size = 3,
             fill = "#52969A",
             alpha = 0.8,)+
  geom_smooth(method=lm,
              formula = y~x,
              se = T,
              fullrange = F,
              alpha = 0.2,
              fill = "#004549",
              color ="#004549" )+
  labs(title=title,
       x = x.lab,
       y = y.lab) +
  stat_cor(label.y = 6,
           method = "spearman")+ #position of R-square
  stat_regline_equation(mapping = NULL, #inheriting above info
                        data = NULL, #inheriting above info
                        formula = y ~x,
                        geom ="text",
                        label.y = 5.8)+ #position of equation using Y-axis as reference
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial",
                                  colour = "black",
                                  face = "italic",
                                  size = "12"),
        panel.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "16"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "14"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "16", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "brd2019-mW-mass-loglog.png",
       path = gr.brd,
       dpi = 300, 
       width = 4.5,
       height =4, 
       units = "in")

# ---- PLOT : Dot plot with raw values ---------
my.color <-c ("chocolate2","darkorchid4")

df.hmd$m <- df.hmd$brood.prp
df.hmd$wolbachia<- factor(df.hmd$wolbachia, +
                            levels = c("uninfected", "infected"))

ggplot(df.hmd , aes(x=humidity,
                    y = m,
                    fill = wolbachia)) + 
  geom_point(aes(fill = wolbachia),
             shape = 21, 
             size = 2)+
  facet_wrap(~colony)+
  #scale_y_continuous(breaks = seq(from = 300, to = 5600 , by = 1500)) +
  scale_fill_manual(values = my.color)+
  labs(title ="Brood Proportion in Colonies \\nMeasured to test Effect of Humidty",
       y = "brood versus adults")+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial",
                                  colour = "black",
                                  face = "italic",
                                  size = "12"),
        panel.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #legend.position = c(0.08,0.80),
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "10"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "11"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "12"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "12", 
                                  margin = margin(t='1', r='0.5', b='0.5', l='1', unit = 'cm')))

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "humidity-colony2019-brdprp.png",
       path = gr.col,
       dpi = 300, 
       width = 4.5,
       height = 4.5, 
       units = "in")


# ---- PLOT : Correlation Matrix ----
# reference:https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html#visualization-of-a-correlation-matrix
# used the correlation matrix for this 

data <- corr18m

png(height= 8.5,
    width= 8,
    units = "in",
    res = 300,
    file="corr2018-sprm-colony.png")

corrplot(data$r, method = "color",
         type="lower", order="alphabet", #half of matrix and alphabetical ordering
         p.mat = data$P, sig.level = 0.05, insig = "blank", # p>0.05 are blank
         col = brewer.pal(n = 6, name = "BrBG"), #color scheme
         #outline = TRUE, #outline for each cell
         tl.col = "black", tl.srt = 45, #change label color and angle
         diag = FALSE, #hide correlation on the diaganol
         title = "Correlation Plot Colony with 2018 Dataset")
dev.off()

# good reference: https://www.guru99.com/r-pearson-spearman-correlation.html
#pick colors - https://htmlcolorcodes.com/

# plotting only correlations with Metabolic Rate
df.corr2018.mr <- data.frame(subset(df.corr2018, trait1 == "met.rate.25"))

# common variables
data = df.corr2018
x = data$trait1
y = data$trait2
value = data$spearman.corr
ggtitle = "Spearman Correlations \\n2019"

ggplot(data,aes(x,y))+
  geom_tile(aes(fill = value),
            color = "gray70")+
  ggtitle(ggtitle)+
  geom_text(aes(x, y, label = value),
            color = "black", size = 4) +
  scale_fill_gradient2(low = "#FF6347",# color scheme
                       mid = "white",
                       high = "cadetblue",
                       midpoint = 0,
                       limit = c(-1,1),
                       name="Correlations \\n(Pearson)",
                       guide = "colourbar") +
  theme_minimal()+ # minimal theme
  theme(plot.title = element_text(family = "Arial", colour = "black",
                                   face = "bold", size = "12", 
                                   margin = margin(t='0.5', r='0.5',
                                                   b='0.5', l='0.5', 
                                                   unit = 'cm')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "10"),
        legend.background = element_rect(fill = NA,linetype = NULL),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_blank(),
        axis.text = element_text(family = "Arial",
                                 colour = "black",
                                 size = 11),
        axis.text.x = element_text(angle = 45,vjust = 1,size = 11,hjust = 1))        
#axis.ticks = element_line(colour = "black", size = 1)

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(filename = "corr2019-sprm-mr.png",
       path = gr.col,
       plot = last_plot(),
       dpi = 300, 
       width = 3.5,
       height = 4, 
       units = "in")
