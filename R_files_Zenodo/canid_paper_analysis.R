#Code to fit models from:
#Responses of sympatric canids to human development revealed through citizen science
#Kellner KF, Hil JE, Gantchoff MG, Kramer DW, Bailey AM, Belant JL. 
#Ecology and Evolution

#Required libraries
library(tidyverse)
library(lme4)
library(car)
library(MuMIn)
library(psych)

#Overdispersion test function from:
#https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#Read in raw data
canid_data <- read_csv("data/canid_data.csv")

#Check correlation among covariates
cor(as.matrix(canid_data[6:8]))
pairs.panels(canid_data[6:8])

#Coyotes-----------------------------------------------------------------------

#Fit model
coy_fit <- glmer(coyote ~ ag + urbanL + para_mn + (1|year) + (1|id) +
             + (1|wmu) + offset(log(hours)), data=canid_data,
             family='poisson',na.action='na.fail')
#Save fit
saveRDS(coy_fit, "data/coyote_model.Rds")

#Look at output (Table 1)
summary(coy_fit)

#Check VIF
vif(coy_fit)

#Check fit
r.squaredGLMM(coy_fit)
overdisp_fun(coy_fit)
pchisq(deviance(coy_fit), df=df.residual(coy_fit), lower.tail=FALSE)

#Red fox-----------------------------------------------------------------------

#Get coyote abundance covariate
coy_sc <- c(scale(canid_data$coyote/canid_data$hours))

#Fit model
rf_fit <- glmer(redfox ~ ag + I(ag^2) + urbanL + para_mn +coy_sc + (1|year) + (1|id)
             + (1|wmu) + offset(log(hours)), data=canid_data,
             family='poisson',na.action='na.fail')
#Save fit
saveRDS(rf_fit, "data/redfox_model.Rds")

#Look at output (Table 2)
summary(rf_fit)

#Check VIF
vif(rf_fit)

#Check fit
r.squaredGLMM(rf_fit)
overdisp_fun(rf_fit)
pchisq(deviance(rf_fit), df=df.residual(rf_fit), lower.tail=FALSE)
