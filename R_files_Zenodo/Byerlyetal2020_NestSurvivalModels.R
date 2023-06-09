#Roseate Tern daily hatch survival, chick survival and predation probability generalized linear mized effects models (GLMER)
#from Byerly et al. 2020: Colony characteristics influence nest survival of Caribbean Roseate Terns

#Packages ----
library(ggplot2)
library(lme4)
library(dplyr)
library(tidyr)
library(car)
library(AICcmodavg)
library(gridExtra)
library(grid)
library(ggeffects)
library(scales)

pubtheme= theme_bw() +theme(panel.grid = element_blank(), axis.text = element_text(size = 16), text=element_text(size=16, family="serif"))

#read in data ---- 
rost=read.csv("Byerlyetal2020_dailysurvival.csv")

#standardize large covariates by centering (mean of 0) and standard deviation of 1 (scaling)
size <- scale(rost$n, center = TRUE, scale = TRUE)
hatchdate  <- scale(rost$HD, center = TRUE, scale = TRUE)
#rename rest of covariates
hatch_age  <- rost$nest_age
chick_age <- rost$chick_age
age <- rost$age
density  <- rost$dens
nearest <- rost$nn
cover <- rost$per_cover
str <- rost$max_str
lagu <- rost$nest_lagu
tern <- rost$nest_tern


#full model to check VIF
full=glmer(survive~size+hatchdate+age+density+nearest+cover+str+lagu+tern+veg+(1|island), family=binomial, data=rost)
car::vif(full)

#logistic exposure model
logexp <- function(exposure = 1) {
  linkfun <- function(mu) qlogis(mu^(1/exposure))
  linkinv <- function(eta) plogis(eta)^exposure
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta)>30,.Machine$double.eps, 
           exp(eta)/(1+exp(eta))^2)
  }
  mu.eta <- function(eta) {
    exposure * plogis(eta)^(exposure-1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),
            class = "link-glm")
}

#code to convert intercept to DSR
logit <- function(x){ log(x/(1-x))}
expit <- function(x){1/(1+exp(-x))}

####################################################################################################################################
#MODEL SET 1: HATCH SUCCESS

#MODEL SETS: random effect (island) + scaled covariates
#hatch_age
mod=list()
mod[[1]] <- glmer(Hatch ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(Hatch ~ hatch_age+hatchdate+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(Hatch ~ hatch_age+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(Hatch ~ hatchdate+(1|island), family=binomial(link=logit), data=rost)

#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 2: colony/site covariates (size, nesting tern, nesting LAGU, density)
mod=list()
mod[[1]] <- glmer(Hatch ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(Hatch ~ size+density+lagu+tern+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(Hatch ~ size+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(Hatch ~ density+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(Hatch ~ size+density+(1|island), family=binomial(link=logit), data=rost)
mod[[6]] <- glmer(Hatch ~ size+tern+(1|island), family=binomial(link=logit), data=rost)
mod[[7]] <- glmer(Hatch ~ size+lagu+(1|island), family=binomial(link=logit), data=rost)
mod[[8]] <- glmer(Hatch ~ density+lagu+(1|island), family=binomial(link=logit), data=rost)
mod[[9]] <- glmer(Hatch ~ density+tern+(1|island), family=binomial(link=logit), data=rost)
#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 3: nest covariates
mod=list()
mod[[1]] <- glmer(Hatch ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(Hatch ~ nearest+place+cover+str+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(Hatch ~ nearest+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(Hatch ~ place+cover+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(Hatch ~ str+(1|island), family=binomial(link=logit), data=rost)
mod[[6]] <- glmer(Hatch ~ nearest+cover+(1|island), family=binomial(link=logit), data=rost)
mod[[7]] <- glmer(Hatch ~ place+cover+str+(1|island), family=binomial(link=logit), data=rost)
mod[[8]] <- glmer(Hatch ~ cover+str+(1|island), family=binomial(link=logit), data=rost)
mod[[9]] <- glmer(Hatch ~ place+(1|island), family=binomial(link=logit), data=rost)
mod[[10]] <- glmer(Hatch ~ cover+(1|island), family=binomial(link=logit), data=rost)
##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

##Stage 4: top covariates from stages 1-3
mod=list()
mod[[1]] <- glmer(Hatch ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(Hatch ~ size+tern+hatch_age+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(Hatch ~ size+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(Hatch ~ hatch_age+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(Hatch ~ hatch_age+size+(1|island), family=binomial(link=logit), data=rost)
mod[[6]] <- glmer(Hatch ~ tern+(1|island), family=binomial(link=logit), data=rost)
##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#model averaging for top models
averaging<-as.data.frame(matrix(rep(NA,5),ncol=5, byrow=FALSE))
names(averaging)<-c("Param","Estimate","SE","Lower CI","Upper CI")
tot<-list(mod[[5]],mod[[2]])
Modnames<-c(5,2)
#model averaging: intercept (DSR), then expit() values reported to exponentiate
a<-modavg(cand.set = tot, parm="(Intercept)", Modnames, conf.level = 0.95, second.ord = T)
#model averaging: age
b<-modavg(cand.set = tot, parm="size", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)
#model averaging: size
c<-modavg(cand.set = tot, parm="hatch_age", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)
#model averaging: size
#back transform coefficients and CI with expit to get correct odds
d<-modavg(cand.set = tot, parm="tern", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)

#PLOT RESULTS
mod <- glmer(Hatch ~ n+hatch_age+(1|island), family=binomial(link=logit), data=rost)
#effects: https://strengejacke.github.io/ggeffects/
colony_n = ggpredict(mod, terms="n")
hatch=ggplot(colony_n, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0.87, 1) +
  labs(title ="a.", 
       x = "Colony size", y = "Daily survival rate\\n") + 
  pubtheme
hatch


hatchage = ggpredict(mod, terms="hatch_age[all]")
age=ggplot(hatchage, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0.85, 1) + xlim(0, 25) +
  labs(title ="b.", 
       x = "Nest age in days", y = "") + 
  pubtheme

grid.arrange(hatch, age, ncol=2) 
####################################################################################################################################

#MODEL SET 2: CHICK SUCCESS

#MODEL SETS: random effect (island) + scaled covariates
#Stage 1: hatch date
mod=list()
mod[[1]] <- glmer(chick ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(chick ~ hatchdate+chick_age+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(chick ~ chick_age+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(chick ~ hatchdate+(1|island), family=binomial(link=logit), data=rost)
#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 2: colony/site covariates (size, nesting tern, nesting LAGU, density)
mod=list()
mod[[1]] <- glmer(chick ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(chick ~ size+density+lagu+tern+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(chick ~ size+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(chick ~ density+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(chick ~ size+density+(1|island), family=binomial(link=logit), data=rost)
mod[[6]] <- glmer(chick ~ size+tern+(1|island), family=binomial(link=logit), data=rost)
mod[[7]] <- glmer(chick ~ size+lagu+(1|island), family=binomial(link=logit), data=rost)
mod[[8]] <- glmer(chick ~ density+lagu+(1|island), family=binomial(link=logit), data=rost)
mod[[9]] <- glmer(chick ~ density+tern+(1|island), family=binomial(link=logit), data=rost)
#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 3: nest covariates
mod=list()
mod[[1]] <- glmer(chick ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(chick ~ nearest+place+cover+str+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(chick ~ nearest+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(chick ~ place+cover+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(chick ~ str+(1|island), family=binomial(link=logit), data=rost)
mod[[6]] <- glmer(chick ~ nearest+cover+(1|island), family=binomial(link=logit), data=rost)
mod[[7]] <- glmer(chick ~ place+cover+str+(1|island), family=binomial(link=logit), data=rost)
mod[[8]] <- glmer(chick ~ cover+str+(1|island), family=binomial(link=logit), data=rost)
mod[[9]] <- glmer(chick ~ place+(1|island), family=binomial(link=logit), data=rost)
mod[[10]] <- glmer(chick ~ cover+(1|island), family=binomial(link=logit), data=rost)
##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 4: top covariates from stages 1-3
mod=list()
mod[[1]] <- glmer(chick ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(chick ~ size+nearest+chick_age+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(chick ~ size+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(chick ~ nearest+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(chick ~ size+chick_age+(1|island), family=binomial(link=logit), data=rost)
mod[[6]] <- glmer(chick ~ chick_age+(1|island), family=binomial(link=logit), data=rost)
##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#model averaging for top models
averaging<-as.data.frame(matrix(rep(NA,5),ncol=5, byrow=FALSE))
names(averaging)<-c("Param","Estimate","SE","Lower CI","Upper CI")
tot<-list(mod[[6]],mod[[2]])
Modnames<-c(6, 2)
#model averaging: intercept (DSR), then expit() values reported to exponentiate
a<-modavg(cand.set = tot, parm="(Intercept)", Modnames, conf.level = 0.95, second.ord = T)
#model averaging: age
b<-modavg(cand.set = tot, parm="chick_age", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)
#model averaging: size
c<-modavg(cand.set = tot, parm="nearest", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)
#model averaging: size
#back transform coefficients and CI with expit to get correct odds
d<-modavg(cand.set = tot, parm="size", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)

#PLOT RESULTS
mod2 <- glmer(chick ~ n+chick_age+(1|island), family=binomial(link=logit), data=rost)
#effects: https://strengejacke.github.io/ggeffects/
colony_n = ggpredict(mod2, terms="n")
chick=ggplot(colony_n, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0.6, 1) +
  labs(title ="c.", 
       x = "Colony size", y = "") + 
  pubtheme
chick
chick=chick + scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01))

grid.arrange(hatch, age, chick, ncol=3) 

####################################################################################################################################

#MODEL SET 3: DAILY PREDATION PROBABILITY

#MODEL SETS: random effect (island) + scaled covariates
#Stage 1: time
mod=list()
mod[[1]] <- glmer(pred ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(pred ~ hatchdate+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(pred ~ hatchdate+age+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(pred ~ age+(1|island), family=binomial(link=logit), data=rost)

#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 2: colony/site covariates (size, nesting tern, nesting LAGU, density)
mod=list()
mod[[1]] <- glmer(pred ~ 1+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[2]] <- glmer(pred ~ size+density+lagu+tern+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[3]] <- glmer(pred ~ size+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[4]] <- glmer(pred ~ density+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[5]] <- glmer(pred ~ size+density+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[6]] <- glmer(pred ~ size+tern+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[7]] <- glmer(pred ~ size+lagu+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[8]] <- glmer(pred ~ size+density+lagu+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[9]] <- glmer(pred ~ size+density+tern+(1|island), family=binomial(link=logexp(rost$t)), data=rost)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#Stage 3: nest covariates
#mod <- glmer(survive ~ col_nn+place+cover+str+(1|nest), family=binomial(link=logexp(rost$obvs)), data=rost)
mod=list()
mod[[1]] <- glmer(pred ~ 1+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[2]] <- glmer(pred ~ nearest+place+cover+str+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[3]] <- glmer(pred ~ nearest+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[4]] <- glmer(pred ~ place+cover+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[5]] <- glmer(pred ~ str+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[6]] <- glmer(pred ~ nearest+cover+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[7]] <- glmer(pred ~ nearest+str+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[8]] <- glmer(pred ~ cover+str+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[9]] <- glmer(pred ~ place+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
mod[[10]] <-glmer(pred ~ cover+(1|island), family=binomial(link=logexp(rost$t)), data=rost)
##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

##Stage 4: top covariates from stages 1-3
mod=list()
mod[[1]] <- glmer(pred ~ 1+(1|island), family=binomial(link=logit), data=rost)
mod[[2]] <- glmer(pred ~ size+age+(1|island), family=binomial(link=logit), data=rost)
mod[[3]] <- glmer(pred ~ size+density+(1|island), family=binomial(link=logit), data=rost)
mod[[4]] <- glmer(pred ~ size+(1|island), family=binomial(link=logit), data=rost)
mod[[5]] <- glmer(pred ~ age+(1|island), family=binomial(link=logit), data=rost)
##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#model averaging for top models
averaging<-as.data.frame(matrix(rep(NA,5),ncol=5, byrow=FALSE))
names(averaging)<-c("Param","Estimate","SE","Lower CI","Upper CI")
tot<-list(mod[[2]])
Modnames<-c(2)
#model averaging: intercept (DSR), then expit() values reported to exponentiate
a<-modavg(cand.set = tot, parm="(Intercept)", Modnames, conf.level = 0.95, second.ord = T)
#model averaging: age
b<-modavg(cand.set = tot, parm="age", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)
#model averaging: size
c<-modavg(cand.set = tot, parm="size", Modnames, c.hat = 1, conf.level = 0.95, second.ord = T)

mod <- glmer(pred ~ n+age+(1|island), family=binomial(link=logit), data=rost)

#PLOT RESULTS
pred1 = ggpredict(mod, terms="age[all]")
pred_age=ggplot(pred1, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 0.2) + xlim(0, 28) +
  labs(title ="a.", 
       x = "Nest age", y = "Daily predation probability\\n") +
  pubtheme
pred_age

pred2 = ggpredict(mod, terms="n")
pred_n=ggplot(pred2, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 0.09) + 
  labs(title ="b.", 
       x = "Colony size", y = "") +
  pubtheme
pred_n
#change number of decimal places on y axes to 2
pred=pred_n + scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01))
  
grid.arrange(pred_age,pred, ncol=2) 