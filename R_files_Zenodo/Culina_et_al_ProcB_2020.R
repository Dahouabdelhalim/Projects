##############################################################
# Author: 
# Antica Culina
# NIOO-KNAW & University of Oxford 
# Email: a.culina@yahoo.com
# pardon me for misspelling....

##############################################################
# Description of script and Instructions
##############################################################

# version of R used in the analysis: 3.6.3

# This code was used to conduct analyses on the influence of a period when a pair has met in year
# on their breeding success in the following season (section 1), and their probability to breed together again (section 2)
# Section 1: Two datasets are used, as year data collection protocol was different between these two sets of years
# one with the data from 2007/08 to 2009/10 years, 2008 to 2010 breeding seasons (section 1A)
# one with the data from 2011/12 to 2013/14 years, 2012 to 2014 breeding seasons (section 1B)

# description of the variables can be found in Data description document


rm(list=ls())

install.packages("lme4", dep = T)
library(MASS)
library(lattice)
require(ggplot2)
library(jtools)

######################## SECTION 1 #############################

####################### 1A: 2007-2010


timing <- read.table(file.choose(), header = T, sep = "") # load data, pairs_data_2007_10
str(timing)
timing$year <- as.factor(timing$year)

### supplementary analysis (includes figures S1, S3)

# Fig S1
#
# subseting for years
#

timing_08 <- subset(timing, timing$year == "2007")
timing_09 <- subset(timing, timing$year == "2008")
timing_10 <- subset(timing, timing$year == "2009")


par(mfrow = c(3,1))

hist(timing_08$first_seen_together, breaks=seq(0.5,8.5, by=1), 
     xlab="", main="2007/08 year", col="gray")
hist(timing_09$first_seen_together, breaks=seq(0.5,8.5, by=1), 
     xlab="", main="2008/09 year", col="gray")
hist(timing_10$first_seen_together, breaks=seq(0.5,8.5, by=1), 
     xlab="", main="2009/10 year", col="gray")


cor.test(timing$f_first_seen, timing$m_first_seen)
cor.test(timing$f_first_seen, timing$first_seen_together)
cor.test(timing$m_first_seen, timing$first_seen_together)


# define vairables to measure difference between SRI of a pair bond and a few other SRI measures, in the
# month when a pair has met


timing$paired1 <- timing$SRI_partner_1met - timing$X50q  # 50 quantile of the distribution for all birds

timing$paired2 <- timing$SRI_partner_1met - timing$X75q   # 75 quantile of the distribution for all birds

timing$paired3 <- timing$SRI_partner_1met - timing$avSRI_F_1met # avareage SRI of a F

timing$paired4 <- timing$SRI_partner_1met - timing$maxSRI_F_1met # max SRI of a F

timing$paired5 <- timing$SRI_partner_1met - timing$avSRI_M_f1met # average SRI of a M

timing$paired6 <- timing$SRI_partner_1met - timing$maxSRI_M_1met # max SRI of a M


# freq. of the values of the above defined var

a1<-as.data.frame(table(cut(timing$paired4,breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1,-0.002, 0),
                            labels=c("-0.5/-0.4","-0.4/-0.3","-0.3/0.2","-0.2/-0.1","-0.1/-0.002", "0"))) )
colnames(a1)<-c("numbers","Freq")
a1   # 2/3 of valuesare between -0.1 and 0,  and 37 values 0

summary(timing$paired4)  # median is - 0.06885; 


a1<-as.data.frame(table(cut(timing$paired6,breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1,-0.002, 0),
                            labels=c("-1/-0.4","-0.4/-0.3","-0.3/0.2","-0.2/-0.1","-0.1/-0.002", "0"))) )
colnames(a1)<-c("numbers","Freq")
a1  


# correlation between these variables

cor.test(timing$paired1, timing$paired2) # 0.991 [0.988-0.993]
cor.test(timing$paired3, timing$paired5) # 0.916 [0.884-0.939] correlation
cor.test(timing$paired4, timing$paired6) # 0.521 [0.388-0.632]


# models on what variables influence the bond strength

#Transformation function
library(car)

f.logit<-function(a)round(logit(round(a,5),adjust=0.025),5)


## does year SRI coorelate with the meeting month, f and male arival time?
# n of associates and relative numebr of months a bird was seen


cor.test(timing$m_first_seen, timing$f_first_seen)
cor.test(timing$X1met_degree_F, timing$X1met_degree_M) 

mod1 <- lm(f.logit(SRI_partner_winter) ~ year  + degree_F_winter + 
             degree_M_winter + f_first_seen *  first_seen_together * m_first_seen, data=timing)

summary(mod1) 


## does relative bond year SRI coorelate with the time of meeting, 
#  and contolleld for n of associates, 

timing$F.relSRI <- timing$SRI_partner_winter/timing$SRI_F_winter
timing$M.relSRI <- timing$SRI_partner_winter/timing$SRI_M_winter

cor.test(timing$F.relSRI, timing$M.relSRI) # 0.789 [0.716-0.844] correlation

mod1 <- lm(f.logit(F.relSRI) ~ year + first_seen_together*m_first_seen*f_first_seen + 
             degree_F_winter+ degree_M_winter, data=timing)

summary(mod1)  


# does bond SRI in the month they met depend on the n of associates 
# and rel n off associates in that month?

timing$relmetdg.F <- timing$X1met_degree_F/timing$degree_F_winter
timing$relmetdg.M <- timing$X1met_degree_M/timing$degree_M_winter

cor.test(timing$relmetdg.F, timing$relmetdg.M) # 0.596 [0.477-0.693]

mod1 <- lm(f.logit(SRI_partner_1met) ~ year + X1met_degree_F + X1met_degree_M + 
             first_seen_together * m_first_seen * f_first_seen, data=timing)
mod2 <- lm(f.logit(SRI_partner_1met) ~ year + relmetdg.F + relmetdg.M + 
             first_seen_together*m_first_seen*f_first_seen, data=timing)

summary(mod1) 
summary(mod2) 

# does meeting month correlate with the n of associates in that month, 
#year degree, and rel n of associatrs in that month?


cor.test(timing$X1met_degree_F, timing$X1met_degree_M) # 0.65 [0.538-0.733] corr

summary(timing$X1met_degree_F)
summary(timing$X1met_degree_M)


### histogram of degree of a F & a M in a month when a pair was first detected together in the same flock

par(mfrow= c(2,2))

hist(timing$X1met_degree_F, breaks = seq(0.5,115.5, 5), 
     main = "Females", xlab = "meeting month degree")
hist(timing$X1met_degree_M, breaks = seq(0.5, 125.5,5), 
     main = "Males", xlab = "meeting month degree", ylab="")


mod1 <- glm(first_seen_together ~ year + X1met_degree_F + X1met_degree_M + 
               f_first_seen*m_first_seen, data=timing, family= poisson)

summary(mod1)

mod1 <- glm(first_seen_together ~ year + degree_F_winter +degree_M_winter + 
               f_first_seen*m_first_seen, data=timing, family= poisson)

summary(mod1) 

mod1 <- glm(first_seen_together ~ year + relmetdg.F +relmetdg.M + 
               f_first_seen*m_first_seen, data=timing, family= poisson)

summary(mod1) 

#for those where they met after they both have arrived

dat <- subset(timing, timing$met_vs_arrival != 0)


mod1 <- lm(first_seen_together ~ f_first_seen*m_first_seen + year, data=dat)
mod2 <- lm(first_seen_together ~ f_first_seen*m_first_seen * year, data=dat)

AIC(mod1, mod2)
summary(mod1)


#### from above, we can conclude that meeting time does not depend on the N of associates in that month; or overal
# it only depends on f and m arival time
## pairs that meet earler have higher year SRI
### but to cover all possibilities control for (sensitivity analysis)
# 1) Meeting month SRI of a pair bond
# 2) year SRI of a pair bond
# 3) relative year SRI (SRI/OveralSRI)
# 4) SRI vs max strngth (of F and M) (paired4 & 6)
# 5) SRI vs average SRI (of F) (paried3)
# 6) SRI vs 50 perentile of the distribution for a month (paired1)


################## MAIN ANALYSIS, and sensitivity analysis ##########################

# models should not have more that around 7-14 parameters 
# for the sample size of 140 (e.g. Harrell's Regression Modelling strategies)


attach(timing)

######## laydate

basic <- glm (laydate ~  1 , family=gaussian)
m1 <- glm (laydate ~  pair_type, family=gaussian)
m3 <- glm (laydate ~  year+first_seen_together*pair_type, family=gaussian)
m4 <- glm (laydate ~  year+first_seen_together+pair_type, family=gaussian)
m5 <- glm (laydate ~  year*first_seen_together+pair_type, family=gaussian)
m6 <- glm (laydate ~  first_seen_together+year*pair_type, family=gaussian)
m7 <- glm (laydate ~  first_seen_together*year, family=gaussian)
m8 <- glm (laydate ~  first_seen_together+year, family=gaussian)
m9 <- glm (laydate ~  first_seen_together, family=gaussian)
m10 <- glm (laydate ~  year, family=gaussian)
m11 <- glm (laydate ~  first_seen_together*pair_type, family=gaussian)
m12 <- glm (laydate ~  first_seen_together+pair_type, family=gaussian)
m13 <- glm (laydate ~  year*pair_type, family=gaussian)
m14 <- glm (laydate ~  year+pair_type, family=gaussian)

AIC(basic, m1,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) # m3 and m8

# diagnostics

plot(m8)
plot(m3)

#### with polinomials

m8_q2 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + year, family=gaussian)
m8_q3 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + 
                 I(first_seen_together^3) + year, family=gaussian)
m3_q3 <- glm (laydate ~  first_seen_together*pair_type + year + 
                 I(first_seen_together^2):pair_type +  pair_type:I(first_seen_together^3), family=gaussian)
m3_q2 <- glm (laydate ~  year+first_seen_together*pair_type + 
                 pair_type:I(first_seen_together^2), family=gaussian)

AIC(m8, m8_q2, m8_q3, m3, m3_q2, m3_q3) # m8_q3 is the best (quadratci is 2 AIC lower)

summary(m8_q2)
summary(m8_q3)

# diagnostics

plot(m8_q2)
plot(m8_q3)


##########################################################
# lets see if we include female/male arival time ##############
##########################################################


F1 <- glm (laydate ~  f_first_seen + I(f_first_seen^2) + I(f_first_seen^3) + year, family=gaussian)
F2 <- glm (laydate ~  f_first_seen +first_seen_together+ I(f_first_seen^2) + 
              I(f_first_seen^3) + year, family=gaussian)
F3 <- glm (laydate ~  first_seen_together + f_first_seen + 
              I(first_seen_together^2) + I(first_seen_together^3) + year, family=gaussian)


M1 <- glm (laydate ~  m_first_seen + I(m_first_seen^2) + I(m_first_seen^3) + year, family=gaussian)
M2 <- glm (laydate ~  m_first_seen +first_seen_together+ I(m_first_seen^2) + 
              I(m_first_seen^3) + year, family=gaussian)
M3 <- glm (laydate ~  first_seen_together + m_first_seen + 
              I(first_seen_together^2) + I(first_seen_together^3) + year, family=gaussian)

AIC(m8_q3, F1, F2, F3, M1, M2, M3) 


### Sensitivity analysis

# 1) Meeting month SRI of a pair bond
# 2) year SRI of a pair bond
# 3) relative year SRI (SRI/OveralSRI)
# 4) SRI vs max strngth (of F and M) (paired4 & 6)
# 5) SRI vs average SRI (of F) (paried3)
# 6) SRI vs 75 perentile of the distribution for a month (paired2)
# 7) meeting degree of F and M

cor.test(X1met_degree_F, X1met_degree_M)


m_q3 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                year -1, family=gaussian)
m_q3_1 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                 year + SRI_partner_winter, family=gaussian)
m_q3_2 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                 year + SRI_partner_1met, family=gaussian)
m_q3_3 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                 year + F.relSRI, family=gaussian)
m_q3_4 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                 year+paired4+paired6, family=gaussian)
m_q3_5 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                 year + paired3, family=gaussian)
m_q3_6 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                 year + paired2, family=gaussian)
m_q3_7a <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                  year + X1met_degree_F, family=gaussian)
m_q3_7b <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3) + 
                  year + X1met_degree_M, family=gaussian)

AIC(m_q3, m_q3_1, m_q3_2, m_q3_3, m_q3_4, m_q3_5, m_q3_6, m_q3_7a, m_q3_7b)  # 

# model with female meeting degree is better supported than the model without
#model with partner SRI has equivalent support to the model without it

summary(m_q3)
summary(m_q3_1) # partner SRI does not reach singificance, while meeting time keeps the same effect
summary(m_q3_7a) # F degree has marginal significance + it is very low effect size, effect of meeting time stays equivalent to the model without  degree



####################
#################### PLOTTING
par(mfrow = c(1,1))

##### what are the values for different months

min(laydate)
coef(m_q3)
newdata2 <- with(timing, data.frame(first_seen_together = rep(seq(from = 1, to = 8, length.out = 8),2), 
                                    year = factor(rep(c("2007", "2008", "2009"), each=16))))

newdata3 <- cbind(newdata2, predict(m_q3, newdata = newdata2, type = "link",se = TRUE))
newdata3 <- within(newdata3, { PredictedProb <- (fit)
LL <- (fit - (1.96 * se.fit))
UL <- (fit + (1.96 * se.fit))})
newdata3


### plot Fig 2A

ggplot(newdata3, aes(x = first_seen_together, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
                                                                                    
                                                                                    ymax = UL, fill = year), alpha = 0.2) + geom_line(aes(colour = year),size = 1) + geom_point(data = timing,
                                                                                                                                                                                    mapping = aes(x = first_seen_together, y = laydate, color = year)) + geom_jitter(aes(colour = year), width = 0.2)

###############################################################################
################################ Clutch size ##################################
###############################################################################


basic <- glm (clutchsz ~  laydate, family=poisson)
m2 <- glm (clutchsz ~  laydate+ first_seen_together*year+pair_type, family=poisson)
m3 <- glm (clutchsz ~  laydate+ first_seen_together+year*pair_type, family=poisson)
m4 <- glm (clutchsz ~  laydate+ first_seen_together*pair_type+year, family=poisson)
m5 <- glm (clutchsz ~  laydate+ first_seen_together+year+pair_type, family=poisson)
m6 <- glm (clutchsz ~  laydate+ first_seen_together*year, family=poisson)
m7 <- glm (clutchsz ~  laydate+ first_seen_together+year, family=poisson)
m8 <- glm (clutchsz ~  laydate+year*pair_type, family=poisson)
m9 <- glm (clutchsz ~  laydate+ year+pair_type, family=poisson)
m10 <- glm (clutchsz ~  laydate+ first_seen_together*pair_type, family=poisson)
m11 <- glm (clutchsz ~  laydate+ first_seen_together+pair_type, family=poisson)
m12 <- glm (clutchsz ~  laydate+ year, family=poisson)
m13 <- glm (clutchsz ~  laydate+ first_seen_together, family=poisson)
m14 <- glm (clutchsz ~  laydate+ pair_type, family=poisson)

AIC(basic, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) #
anova(basic, m13, test="Chi")

plot(basic) #diagnostics more less ok

basic <- glm (clutchsz ~  laydate, family=quasipoisson)

summary(basic)

summary(m13) ################## the term first_seen_togheter does not reach stat sign

m_q2 <- glm (clutchsz ~  laydate+ first_seen_together + I(first_seen_together^2), family=poisson)
m_q3 <- glm (clutchsz ~  laydate+ first_seen_together+ I(first_seen_together^2)+ 
                I(first_seen_together^3), family=poisson)

AIC(m13, m_q2, m_q3) # quadratic and cubic prefrom worse thant m13

###### female / male arrival time controlled for

F1 <- glm (clutchsz ~  laydate+ first_seen_together + f_first_seen, family=poisson)
F2 <- glm (clutchsz ~  laydate+ f_first_seen, family=poisson)
M1 <- glm (clutchsz ~  laydate+ first_seen_together + m_first_seen, family=poisson)
M2 <- glm (clutchsz ~  laydate+ m_first_seen, family=poisson)

AIC(m13, F1, F2, M1, M2)  # all are close to m13, but none is better


###############################################################################
################################ N of chicks ##################################
###############################################################################


basic <- glm (n_chicks ~  clutchsz, family=poisson)
m2 <- glm (n_chicks ~  clutchsz+ first_seen_together*year+pair_type, family=poisson)
m3 <- glm (n_chicks ~  clutchsz+ first_seen_together+year*pair_type, family=poisson)
m4 <- glm (n_chicks ~  clutchsz+ first_seen_together*pair_type+year, family=poisson)
m5 <- glm (n_chicks ~  clutchsz+ first_seen_together+year+pair_type, family=poisson)
m6 <- glm (n_chicks ~  clutchsz+ first_seen_together*year, family=poisson)
m7 <- glm (n_chicks ~  clutchsz+ first_seen_together+year, family=poisson)
m8 <- glm (n_chicks ~  clutchsz+year*pair_type, family=poisson)
m9 <- glm (n_chicks ~  clutchsz+ year+pair_type, family=poisson)
m10 <- glm (n_chicks ~  clutchsz+ first_seen_together*pair_type, family=poisson)
m11 <- glm (n_chicks ~  clutchsz+ first_seen_together+pair_type, family=poisson)
m12 <- glm (n_chicks ~  clutchsz+ year, family=poisson)
m13 <- glm (n_chicks ~  clutchsz+ first_seen_together, family=poisson)
m14 <- glm (n_chicks ~  clutchsz+ pair_type, family=poisson)

AIC(basic, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) 

basic <- glm (n_chicks ~  clutchsz, family=quasipoisson)
summary(basic)



###############################################################################
################################ n fledged for all pairs ######################
###############################################################################


basic <- glm (n_fledged ~  n_chicks, family=poisson)
m2 <- glm (n_fledged ~  n_chicks+ first_seen_together*year+pair_type, family=poisson)
m3 <- glm (n_fledged ~  n_chicks+ first_seen_together+year*pair_type, family=poisson)
m4 <- glm (n_fledged ~  n_chicks+ first_seen_together*pair_type+year, family=poisson)
m5 <- glm (n_fledged ~  n_chicks+ first_seen_together+year+pair_type, family=poisson)
m6 <- glm (n_fledged ~  n_chicks+ first_seen_together*year, family=poisson)
m7 <- glm (n_fledged ~  n_chicks+ first_seen_together+year, family=poisson)
m8 <- glm (n_fledged ~  n_chicks+year*pair_type, family=poisson)
m9 <- glm (n_fledged ~  n_chicks+ year+pair_type, family=poisson)
m10 <- glm (n_fledged ~  n_chicks+ first_seen_together*pair_type, family=poisson)
m11 <- glm (n_fledged ~  n_chicks+ first_seen_together+pair_type, family=poisson)
m12 <- glm (n_fledged ~  n_chicks+ year, family=poisson)
m13 <- glm (n_fledged ~  n_chicks+ first_seen_together, family=poisson)
m14 <- glm (n_fledged ~  n_chicks+ pair_type, family=poisson)

AIC(basic, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) # several models have same support. 
#the simplest one (the lowest n of parameters) is nasic and meeting time one

summary(basic)

summary(m13)
# polinomial with m13

m13_q2 <- glm (n_fledged ~  n_chicks+ first_seen_together + I(first_seen_together^2), family=poisson)
m13_q3 <- glm (n_fledged ~  n_chicks+ first_seen_together + I(first_seen_together^2) + 
                  I(first_seen_together^3), family=poisson)

AIC(m13, m13_q2, m13_q3)  # simmialr support to all


### f/m arival time (add on m13) as a basis

F1 <- glm (n_fledged ~  n_chicks + f_first_seen, family=poisson)
M1 <- glm (n_fledged ~  n_chicks + m_first_seen, family=poisson)
F2 <- glm (n_fledged ~  n_chicks+  f_first_seen + first_seen_together, family=poisson)
M2 <- glm (n_fledged ~  n_chicks+ m_first_seen + first_seen_together, family=poisson)

AIC(m13, F1, F2, M1, M2) # m1 is best supported

M1 <- glm (n_fledged ~  n_chicks + m_first_seen, family=quasipoisson)


summary(M1)


############# n fledged, fledged more than 1

t2 <- subset(timing, timing$n_fledged > "0") # as of sample size, some of the models included in the
#previous fledgling model, were now excluded

basic<- glm (n_fledged ~  n_chicks, family=poisson, data=t2)
m2 <- glm (n_fledged ~  n_chicks+ first_seen_together*year+pair_type, family=poisson, data=t2)
m4 <- glm (n_fledged ~  n_chicks+ first_seen_together*pair_type+year, family=poisson, data=t2)
m5 <- glm (n_fledged ~  n_chicks+ first_seen_together+year+pair_type, family=poisson, data=t2)
m6 <- glm (n_fledged ~  n_chicks+ first_seen_together*year, family=poisson, data=t2)
m7 <- glm (n_fledged ~  n_chicks+ first_seen_together+year, family=poisson, data=t2)
m9 <- glm (n_fledged ~  n_chicks+ year+pair_type, family=poisson, data=t2)
m10 <- glm (n_fledged ~  n_chicks+ first_seen_together*pair_type, family=poisson, data=t2)
m11 <- glm (n_fledged ~  n_chicks+ first_seen_together+pair_type, family=poisson, data=t2)
m12 <- glm (n_fledged ~  n_chicks+ year, family=poisson, data=t2)
m13 <- glm (n_fledged ~  n_chicks+ first_seen_together, family=poisson, data=t2)
m14 <- glm (n_fledged ~  n_chicks+ pair_type, family=poisson, data=t2)

AIC(basic, m2, m4,m5,m6, m7, m9, m10, m11, m12, m13, m14) # intercept is the best supported


###############################################################################
################################ Hatching success #############################
###############################################################################

#################### laydate vs laydate + clutchsize

basic <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  laydate, family=binomial)
basic2 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  laydate+ clutchsz, family=binomial)

AIC(basic, basic2) # laydate only is better, so continue with that

m2 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  first_seen_together*year+pair_type+ laydate , family=binomial)
m3 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+year*pair_type+ laydate , family=binomial)
m4 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together*pair_type+year+ laydate, family=binomial)
m5 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+year+pair_type+ laydate, family=binomial)
m6 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together*year+ laydate, family=binomial)
m7 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+year+ laydate, family=binomial)
m8 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ year*pair_type+ laydate, family=binomial)
m9 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~ year+pair_type+ laydate, family=binomial)
m10 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ laydate, family=binomial)
m11 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+pair_type+ laydate, family=binomial)
m12 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~   year+ laydate, family=binomial)
m13 <-glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+ laydate, family=binomial)
m14 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  pair_type+ laydate, family=binomial)

AIC(basic, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14)

# m4  (AIC 375.57) and m10 are best supported

m4 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together*pair_type+year+ 
              laydate, family=quasibinomial)

m10 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
               laydate, family=quasibinomial)

summary(m4)
summary (m10) # the only significant term in both is laydate + both models are overdispersed


##### polinomials, only quadratic were considered for m4 and m10 
#####because cubic would intorduce to many paremeters to estimate vs sample size


m4_q2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ year + first_seen_together*pair_type+ 
                 laydate + I(first_seen_together^2):pair_type, family=binomial)


m7_q2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ year +
                  first_seen_together+ 
                  laydate + I(first_seen_together^2), family=binomial)
m7_q3 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ year +
                 first_seen_together+ 
                 laydate + I(first_seen_together^2)+ I(first_seen_together^3), family=binomial)
m10_q2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ laydate + 
                  I(first_seen_together^2):pair_type, family=binomial)

AIC(m4, m4_q2, m7, m7_q2, m7_q3, m10, m10_q2) # polinomial models do not preform any better

###### female / male arrival time controlled for

F1 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ 
             f_first_seen*pair_type + laydate, family=binomial)
F2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ first_seen_together + f_first_seen*pair_type+ 
             laydate, family=binomial)
F3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~
             first_seen_together*pair_type+ f_first_seen + laydate , family=binomial)


M1 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ 
              m_first_seen*pair_type + laydate, family=binomial)
M2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ first_seen_together + m_first_seen*pair_type+ 
              laydate, family=binomial)
M3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~
              first_seen_together*pair_type+ m_first_seen + laydate , family=binomial)

AIC(m10, F3, M3, M1, F1, M2, F2) 


### sensitivity analysis


m10_1 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                 laydate + SRI_partner_winter, family=binomial)
m10_2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                 laydate + SRI_partner_1met, family=binomial)
m10_3 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                 laydate + F.relSRI, family=binomial)
m10_4 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                 laydate + paired4 + paired6, family=binomial)
m10_5 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                 laydate + paired3, family=binomial)
m10_6 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                 laydate + paired2, family=binomial)
m10_7a <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                  laydate +  X1met_degree_F, family=binomial)
m10_7b <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
                  laydate +  X1met_degree_M, family=binomial)

AIC(m10 ,m10_1, m10_2, m10_3, m10_4, m10_5, m10_6, m10_7a, m10_7b) # m10_ 4 is far best supported


summary(m10_4) # paired 6 but not paired 4 has significant influence


#### plot the influence of max male partner SRI vs max sri (paired6), Fig S4

mplot1 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~laydate + pair_type*first_seen_together + 
                  paired6, family=quasibinomial)

summary(mplot1)

# for this model the estimate of the paired6 stayes the same as for m10_4


newdata1 <- with(timing, data.frame(paired6 = rep(seq(from = -1, to = 0, length.out = 140),), 
                                    pair_type = factor(rep(c("newnew", "newold", "oldnew", "oldold"), each=140)), laydate= mean(laydate),
                                    first_seen_together=mean(first_seen_together)))

newdata2 <- cbind(newdata1, predict(mplot1, newdata = newdata1, type = "link",se = TRUE))
newdata2 <- within(newdata2, { PredictedProb <- plogis(fit)
LL <- plogis(fit - (1.96 * se.fit))
UL <- plogis(fit + (1.96 * se.fit))

})
newdata2


ggplot(newdata2, aes(x = paired6, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
                                                                        
                                                                        ymax = UL, fill = pair_type), alpha = 0.2) + geom_line(aes(colour = pair_type),size = 1)


# no CI (this is used for Fig S4)

ggplot(newdata2, aes(x = paired6, y = PredictedProb)) + geom_line(aes(colour = pair_type),size = 1)



###############################################################################
################################ fledgling success ######################
###############################################################################


##
####### proportion of young fledged out of chicks
##

basic <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  laydate, family=binomial)
basic2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  laydate+ clutchsz, family=binomial)
basic3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  clutchsz, family=binomial)

AIC(basic,basic2, basic3)

summary(basic2)
summary(basic3)

anova(basic2, basic3, test = "Chisq")

# so keep only clutch size here


basic <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  clutchsz, family=binomial)
m2 <- glm (cbind(n_fledged, (n_chicks - n_fledged)) ~  first_seen_together*year+pair_type+ clutchsz, family=binomial)
m3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year*pair_type + clutchsz-1, family=binomial)
m4 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*pair_type+year + clutchsz - 1, family=binomial)
m5 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+pair_type+ clutchsz, family=binomial)
m6 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*year+ clutchsz, family=binomial)
m7 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year + clutchsz, family=binomial)
m8 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~ year*pair_type + clutchsz, family=binomial)
m9 <- glm (cbind(n_fledged,(n_chicks - n_fledged))  ~ year+pair_type + clutchsz, family=binomial)
m10 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~first_seen_together*pair_type + clutchsz, family=binomial)
m11 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+pair_type + clutchsz, family=binomial)
m12 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~   year+ clutchsz, family=binomial)
m13 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+ clutchsz, family=binomial)
m14 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  pair_type+ clutchsz, family=binomial)

AIC(basic, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) # m4 best supported

summary(m4) # overdispersed

# polinomial models, based on m5 


m5_q2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together + pair_type+
                 year + I(first_seen_together^2) + clutchsz, family=binomial)
m5_q3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together + pair_type+
                 year + I(first_seen_together^2)+ I(first_seen_together^3) + clutchsz, family=binomial)

AIC(m5_q2, m5_q3, m5) 


# female and male arrival time

F1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*pair_type+year + 
              f_first_seen + clutchsz - 1, family=binomial)
F2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  f_first_seen*pair_type+year + 
            + clutchsz - 1, family=binomial)
F3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  f_first_seen*pair_type+year + 
              first_seen_together + clutchsz - 1, family=binomial)
M1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*pair_type+year + 
              m_first_seen + clutchsz - 1, family=binomial)
M2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  m_first_seen*pair_type+year + 
              + clutchsz - 1, family=binomial)
M3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  m_first_seen*pair_type+year + 
              first_seen_together + clutchsz, family=binomial)


AIC(M3, F1, M2,  M1, F3, F2)

M3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  m_first_seen*pair_type+year + 
              first_seen_together + clutchsz, family=quasibinomial)
F1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*pair_type+year + 
              f_first_seen + clutchsz - 1, family=quasibinomial)

summary(M3)

summary(F1)



### sensitivity analysis

m4_1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  SRI_partner_1met + first_seen_together+ m_first_seen* pair_type+
                year + clutchsz, family=binomial)
m4_2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  SRI_partner_winter + first_seen_together + m_first_seen * pair_type+
                year + clutchsz, family=binomial)
m4_3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  F.relSRI + first_seen_together +m_first_seen* pair_type+
                year + clutchsz, family=binomial)
m4_4 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~ paired4 + paired6 + first_seen_together + m_first_seen* pair_type+
                year +  clutchsz, family=binomial)
m4_5 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  paired3 + first_seen_together + m_first_seen* pair_type+
                year + clutchsz, family=binomial)
m4_6 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  paired2 + first_seen_together + m_first_seen* pair_type+
                year + clutchsz, family=binomial)
m4_7a <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  X1met_degree_F + first_seen_together + m_first_seen* pair_type+
                 year +  clutchsz, family=binomial)
m4_7b <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  X1met_degree_M + first_seen_together + m_first_seen* pair_type+
                 year + clutchsz, family=binomial)

AIC(M3, m4_1, m4_2, m4_3, m4_4, m4_5,m4_6, m4_7a, m4_7b) # m4_7b with 876.73 AIC

# overdispersion accounted for

m4_7b <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  X1met_degree_M + first_seen_together + m_first_seen* pair_type+
                 year + clutchsz, family=quasibinomial)

summary(m4_7b) # no term significant

## the outputs of these models show that, after accounting for overdisp., meeting time is in no way sign

########## fl_succ, only those pairs that have fledged 

t3 <- subset(timing, timing$n_fledged>0)

basic <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  clutchsz, family=binomial, data = t3)
basic2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  laydate+ clutchsz, family=binomial, data = t3)
basic3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  laydate, family=binomial, data = t3)


AIC(basic, basic2, basic3)

basic <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  clutchsz, family=binomial, data = t3)
m2 <- glm (cbind(n_fledged, (n_chicks - n_fledged)) ~  first_seen_together*year+pair_type + clutchsz, family=binomial, data=t3)
m4 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*pair_type+year+ clutchsz - 1, family=binomial, data=t3)
m5 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+pair_type+ clutchsz, family=binomial, data=t3)
m6 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*year+ clutchsz, family=binomial, data=t3)
m7 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+ clutchsz, family=binomial, data=t3)
m9 <- glm (cbind(n_fledged,(n_chicks - n_fledged))  ~ year+pair_type+ clutchsz, family=binomial, data=t3)
m10 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~first_seen_together*pair_type+ clutchsz, family=binomial, data=t3)
m11 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+pair_type+ clutchsz, family=binomial, data=t3)
m12 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~   year+ clutchsz, family=binomial, data=t3)
m13 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+ clutchsz, family=binomial, data=t3)
m14 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  pair_type+ clutchsz, family=binomial, data=t3)

AIC(basic, m2, m4,m5,m6, m7, m9, m10, m11, m12, m13, m14) # m2 and m4 are the best

m2 <- glm (cbind(n_fledged, (n_chicks - n_fledged)) ~  first_seen_together*year+pair_type + 
              clutchsz, family=quasibinomial, data=t3)
m4 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together*pair_type+year+ 
              clutchsz, family=quasibinomial, data=t3)

summary(m2)
summary(m4) # when overdispersion avvounted for, meeting time looses sign and no interaction is sign

# polinomials, based on m5 as of number of paramenters vs sample size. Plus no interaction is
# sign if m2 and m4

m_q2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+ 
                pair_type+ clutchsz + I(first_seen_together^2), family=binomial, data=t3)
m_q3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+ 
                pair_type+ clutchsz + I(first_seen_together^2) + I(first_seen_together^3), family=binomial, data=t3)

AIC(m5, m_q2, m_q3) # no support

# with arivla times

F1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+pair_type+ 
              clutchsz + f_first_seen, family=binomial, data=t3)
F2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  f_first_seen+year+pair_type+ 
              clutchsz, family=binomial, data=t3)
M1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  first_seen_together+year+pair_type+ 
              clutchsz + m_first_seen, family=binomial, data=t3)
M2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  m_first_seen+year+pair_type+ 
              clutchsz, family=binomial, data=t3)


AIC(m5, M1, M2, F1, F2) # m5

### sensitivity analysis

m5_1 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  SRI_partner_1met + first_seen_together + pair_type+
                year + clutchsz, family=binomial, data=t3)
m5_2 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  SRI_partner_winter + first_seen_together + pair_type+
                year + clutchsz, family=binomial, data=t3)
m5_3 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  F.relSRI + first_seen_together + pair_type+
                year + clutchsz, family=binomial, data=t3)
m5_4 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~ paired4 + paired6 + first_seen_together + pair_type+
                year + clutchsz, family=binomial, data=t3)
m5_5 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  paired3 + first_seen_together + pair_type+
                year + clutchsz, family=binomial, data=t3)
m5_6 <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  paired2 + first_seen_together + pair_type+
                year + clutchsz, family=binomial,data=t3)
m5_7a <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  X1met_degree_F + first_seen_together + pair_type+
                 year + clutchsz, family=binomial, data=t3)
m5_7b <- glm (cbind(n_fledged, (n_chicks - n_fledged))  ~  X1met_degree_M + first_seen_together + pair_type+
                 year + clutchsz, family=binomial, data=t3)


AIC(m5, m5_1, m5_2, m5_3, m5_4, m5_5,m5_6, m5_7a, m5_7b) #m5


#################### binary - flegded vs not fledged ######

basic <- glm (succ_failure ~  first_seen_together*year+pair_type+ laydate, family=binomial)
basic2 <- glm (succ_failure ~  first_seen_together*year+pair_type+ clutchsz, family=binomial)
basic3 <- glm (succ_failure ~  first_seen_together*year+pair_type+ laydate + clutchsz, family=binomial)

AIC(basic, basic2, basic3) # keep laydat

m2 <- glm (succ_failure ~  first_seen_together*year+pair_type+ laydate, family=binomial)
m3 <- glm (succ_failure  ~  first_seen_together+year*pair_type+ laydate -1, family=binomial)
m4 <- glm (succ_failure  ~  first_seen_together*pair_type+year+ laydate - 1, family=binomial)
m5 <- glm (succ_failure ~  first_seen_together+year+pair_type+ laydate, family=binomial)
m6 <- glm (succ_failure  ~  first_seen_together*year+ laydate, family=binomial)
m7 <- glm (succ_failure  ~  first_seen_together+year+ laydate, family=binomial)
m8 <- glm (succ_failure ~ year*pair_type+ laydate, family=binomial)
m9 <- glm (succ_failure ~ year+pair_type+ laydate, family=binomial)
m10 <- glm (succ_failure ~first_seen_together*pair_type+ laydate, family=binomial)
m11 <- glm (succ_failure  ~  first_seen_together+pair_type+ laydate, family=binomial)
m12 <- glm (succ_failure ~   year+ laydate, family=binomial)
m13 <- glm (succ_failure ~  first_seen_together+ laydate, family=binomial)
m14 <- glm (succ_failure~  pair_type+ laydate, family=binomial)

AIC(basic, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) # m12 and m13 are the best

summary(m13) # terms for meeting time not significant and small

m13_q2 <- glm (succ_failure ~  first_seen_together+ laydate + I(first_seen_together^2), family=binomial)
m13_q3 <- glm (succ_failure ~  first_seen_together+ laydate + I(first_seen_together^2)+ I(first_seen_together^3), family=binomial)

AIC(m13, m13_q2, m13_q3) # no support for polinomial effects

# with male and female arival times

F1 <- glm (succ_failure ~  f_first_seen+ laydate, family=binomial)
M1 <- glm (succ_failure ~  m_first_seen+ laydate, family=binomial)
F2<- glm (succ_failure ~  first_seen_together+ f_first_seen + laydate, family=binomial)
M2 <- glm (succ_failure ~  first_seen_together + m_first_seen + laydate, family=binomial)

AIC(m13, F1, F2, M1, M2) # no best support for any




#######################################################################################################
#############################  SECTION 1B: 2011-2014  #################################################
#######################################################################################################
#######################################################################################################

bs <- read.table(file.choose(), header = T, sep = "") # load the dataset (pairs_data_2011_14)

str(bs)
bs$year <- as.factor(bs$year)


# Figure S2


f2011 <- subset(bs, bs$year =="2011")
f2012 <- subset(bs, bs$year =="2012")
f2013 <- subset(bs, bs$year =="2013")

par(mfrow=c(3,1))
hist(f2011$first_seen_together, breaks=seq(0.5,26.5, by=1), xlab="", main = "", col="gray")
title("2011/12 year", adj = 0.8, line = -0.2)
hist(f2012$first_seen_together, breaks=seq(0.5,26.5, by=1), xlab="", main="", col="gray")
title("2012/13 year", adj = 0.8, line = -0.2)
hist(f2013$first_seen_together, breaks=seq(0.5,26.5, by=1), xlab="weekend", main="", col="gray", ylim= c(0,15))
title("2013/14 year", adj = 0.8, line = -0.2)




### define some variables to descibe bond strenght: difference between parter SRI when firts met and a few other SRI measures

bs$paired1 <- bs$SRI_partner_1met - bs$X50q  # 50 quanile of the distribution for all birds

bs$paired2 <- bs$SRI_partner_1met - bs$X75q   # 75 quantile of the distribution for all birds

bs$paired3 <- bs$SRI_partner_1met - bs$avSRI_F_1met # average SRI of a F

bs$paired4 <- bs$SRI_partner_1met - bs$maxSRI_F_1met # max SRI of a F

bs$paired5 <- bs$SRI_partner_1met - bs$avSRI_M_f1met # average SRI of a M

bs$paired6 <- bs$SRI_partner_1met - bs$maxSRI_M_1met # max SRI of a M

# freq. of the values of the variables above

a1<-as.data.frame(table(cut(bs$paired4,breaks=c(-0.5,-0.4,-0.3,-0.2,-0.16,-0.002, 0),
                            labels=c("-0.5/-0.4","-0.4/-0.3","-0.3/0.2","-0.2/-0.1","-0.1/-0.002", "0"))) )
colnames(a1)<-c("numbers","Freq")
a1   # 1/2 of valuesare between -0.1 and 0,  and 64 values 0

summary(bs$paired4)  # median is - 0.093; 


a1<-as.data.frame(table(cut(bs$paired6,breaks=c(-0.5,-0.4,-0.3,-0.2,-0.16,-0.002, 0),
                            labels=c("-1/-0.4","-0.4/-0.3","-0.3/0.2","-0.2/-0.1","-0.1/-0.002", "0"))) )
colnames(a1)<-c("numbers","Freq")
a1   # 1/2 of valuesare between -0.1 and 0,  and 71 values 0

summary(bs$paired6)  # median is - 0.098; 

cor.test(bs$paired1, bs$paired2) #  0.99 cc
cor.test(bs$paired3, bs$paired5) # 0.89  cc
cor.test(bs$paired4, bs$paired6) # 0.37


### supplementary analysis (includes figur S3)


## does year SRI corelate with the time of meeting, f and male arival time? and n of associates


cor.test(bs$m_first_seen, bs$f_first_seen) # 0.508 [0.409-0.596]

mod1 <- lm(f.logit(SRI_partner_winter) ~ year +  + degree_F_winter + 
             degree_M_winter + f_first_seen *  first_seen_together * m_first_seen, data=bs)

summary(mod1) 

## does relative bond SRI in year coorelate on the time of meeting, and contolleld for n of associates

bs$F.relSRI <- bs$SRI_partner_winter/bs$SRI_F_winter
bs$M.relSRI <- bs$SRI_partner_winter/bs$SRI_M_winter

cor.test(bs$F.relSRI, bs$M.relSRI) # 0.92 [0.89-0.93]

mod1 <- lm(f.logit(M.relSRI) ~ year + first_seen_together*m_first_seen*f_first_seen + 
             degree_F_winter+ degree_M_winter, data=bs)

summary(mod1)  


# does bond SRI in the month they met depend on the n of associates & rel n of associates in that month?

bs$relmetdg.F <- bs$X1met_degree_F/bs$degree_F_winter
bs$relmetdg.M <- bs$X1met_degree_M/bs$degree_M_winter

cor.test(bs$relmetdg.F, bs$relmetdg.M) # 0.545 [0.451-0.628]

mod1 <- lm(f.logit(SRI_partner_1met) ~ year + X1met_degree_F + X1met_degree_M + 
             first_seen_together * m_first_seen * f_first_seen, data=bs)
mod2 <- lm(f.logit(SRI_partner_1met) ~ year + relmetdg.F + relmetdg.M + 
             first_seen_together*m_first_seen*f_first_seen, data=bs)

summary(mod1) 
summary(mod2) 

# does time of meeting depend on degree in that motnh, year degree, and rel degree in that month?

cor.test(bs$X1met_degree_F, bs$X1met_degree_M) # 0.68 [0.612-0.746] 

summary(bs$X1met_degree_F) 
summary(bs$X1met_degree_M)

# ranges between 1 and 69 - so quite a large spread, mean of 23

# Fig S3

par(mfrow= c(1,2))

hist(bs$X1met_degree_F, breaks = seq(1,80, 1), 
     main = "", xlab = "meeting weekend degree")
hist(bs$X1met_degree_M, breaks = seq(1,80, 1), 
     main = "", xlab = "meeting weekend degree", ylab="")

mod1 <- glm(first_seen_together ~ year + X1met_degree_F + X1met_degree_M + 
              f_first_seen*m_first_seen, data=bs, family= poisson)

summary(mod1) #

mod1 <- glm(first_seen_together ~ year + degree_F_winter +degree_M_winter + 
              f_first_seen*m_first_seen, data=bs, family= poisson)

summary(mod1) 

mod1 <- glm(first_seen_together ~ year + relmetdg.F +relmetdg.M + f_first_seen*m_first_seen, data=bs, family= poisson)

summary(mod1) 

mod1 <- lm(first_seen_together ~ m_first_seen*f_first_seen* year, data=bs)
mod2 <- lm(first_seen_together ~ m_first_seen*f_first_seen + year, data=bs)

AIC(mod1, mod2) 

summary(mod1)

#for those where they met after they both have arrived

dat <- subset(bs, bs$met_vs_arrival != 0)  # in 99 pair both partners arrived before they were first seen associated

mod1 <- lm(first_seen_together ~ m_first_seen*f_first_seen* year, data=dat)
mod2 <- lm(first_seen_together ~ m_first_seen*f_first_seen + year, data=dat)

AIC(mod1, mod2) 

summary(mod2)



#### Main model selection & sensitivity analysis

attach(bs)

############# laydate

basic <- glm (laydate ~  1, family=gaussian, data=bs)
m1 <- glm (laydate ~  pair_type, family=gaussian, data=bs)
m2 <- glm (laydate ~  first_seen_together*year + first_seen_together*pair_type, family=gaussian, data=bs)
m3 <- glm (laydate ~  year+first_seen_together*pair_type, family=gaussian, data=bs)
m4 <- glm (laydate ~  year+first_seen_together+pair_type, family=gaussian, data=bs)
m5 <- glm (laydate ~  year*first_seen_together+pair_type, family=gaussian, data=bs)
m6 <- glm (laydate ~  first_seen_together+year*pair_type, family=gaussian, data=bs)
m7 <- glm (laydate ~  first_seen_together*year, family=gaussian, data=bs)
m8 <- glm (laydate ~  first_seen_together+year, family=gaussian, data=bs)
m9 <- glm (laydate ~  first_seen_together, family=gaussian, data=bs)
m10 <- glm (laydate ~  year, family=gaussian, data=bs)
m11 <- glm (laydate ~  first_seen_together*pair_type, family=gaussian, data=bs)
m12 <- glm (laydate ~  first_seen_together+pair_type, family=gaussian, data=bs)
m13 <- glm (laydate ~  year*pair_type, family=gaussian, data=bs)
m14 <- glm (laydate ~  year+pair_type, family=gaussian, data=bs)
m15 <- glm (laydate ~  first_seen_together*year + year*pair_type, family=gaussian, data=bs)
m16 <- glm (laydate ~  pair_type*year + first_seen_together*pair_type, family=gaussian, data=bs)

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16) # m8 & m9 the best

summary(m8)
summary(m9) #  both give the similar estimate for the slope of the effect of meeting time

plot(m8)

# polinomial models

m8_q2 <- glm (laydate ~  first_seen_together+year + I(first_seen_together^2), family=gaussian, data=bs)
m8_q3 <- glm (laydate ~  first_seen_together+year + I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
m9_q2 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) , family=gaussian, data=bs)
m9_q3 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)

AIC(m8, m8_q2, m8_q3, m9, m9_q2, m9_q3) # m9_q3, followed by m8_q3 and m9_q2

summary(m9_q3)


# female and male arrival times


F1 <- glm (laydate ~  first_seen_together +f_first_seen + 
              I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
F2 <- glm (laydate ~  first_seen_together +f_first_seen + I(f_first_seen^2) + 
              I(f_first_seen^3), family=gaussian, data=bs)
F3 <- glm (laydate ~  f_first_seen + I(f_first_seen^2) + I(f_first_seen^3), family=gaussian, data=bs)

M1 <- glm (laydate ~  first_seen_together +m_first_seen + I(first_seen_together^2) + 
              I(first_seen_together^3), family=gaussian, data=bs)
M2 <- glm (laydate ~  first_seen_together +m_first_seen + I(m_first_seen^2) + 
              I(m_first_seen^3), family=gaussian, data=bs)
M3 <- glm (laydate ~  m_first_seen + I(m_first_seen^2) + I(m_first_seen^3), family=gaussian, data=bs)

AIC(m9_q3, M1, F1, M2, F2, M3, F3) # lower support for M added models


### sensitivity analysis


s1 <- glm (laydate ~  first_seen_together + I(first_seen_together^2) + 
             I(first_seen_together^3) + SRI_partner_1met, family=gaussian, data=bs)
s2 <- glm (laydate ~  first_seen_together + SRI_partner_winter + 
             I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s3 <- glm (laydate ~  first_seen_together + F.relSRI + 
             I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s4 <- glm (laydate ~  first_seen_together + paired4 + paired6 + 
             I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s5 <- glm (laydate ~  first_seen_together  + paired3 +
             I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s6 <- glm (laydate ~  first_seen_together + paired2 + 
             I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s7a <- glm (laydate ~  first_seen_together + X1met_degree_F + 
              I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s7b <- glm (laydate ~  first_seen_together + X1met_degree_M + 
              I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)
s7c <- glm (laydate ~  first_seen_together + relmetdg.F + 
              I(first_seen_together^2) + I(first_seen_together^3), family=gaussian, data=bs)

AIC(m9_q3, s1, s2,s3,s4,s5,s6,s7a,s7b,s7c) # best is s2

summary(m9_q3)
summary(s2)

## effect of meeting time stays the same thouth these models

# PLOT results

#ggplot(bs, aes(x = SRI_partner_winter, y = laydate)) + 
#  geom_point() +
#  stat_smooth(method = "lm", col = "black")

### Fig 2B

newdata2 <- with(bs, data.frame(first_seen_together = rep(seq(from = 1, to = 26, length.out = 26),2)))

newdata3 <- cbind(newdata2, predict(m9_q3, newdata = newdata2, type = "link",se = TRUE))
newdata3 <- within(newdata3, { PredictedProb <- (fit)
LL <- (fit - (1.96 * se.fit))
UL <- (fit + (1.96 * se.fit))})
newdata3

ggplot(newdata3, aes(x = first_seen_together, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
                                                                                    
                                                                                    ymax = UL), alpha = 0.2) + geom_line() + geom_point(data = bs,
                                                                                                  mapping = aes(x = first_seen_together, y = laydate)) + geom_jitter(width = 0.2)

################################################################################################
################################ CLUTCH size
################################################################################################

basic <- glm (clutchsz ~  laydate, family=poisson, data=bs)
m1 <- glm (clutchsz ~  laydate + pair_type, family=poisson, data=bs)
m2 <- glm (clutchsz ~  laydate + first_seen_together*year + first_seen_together*pair_type, family=poisson, data=bs)
m3 <- glm (clutchsz ~  laydate + year+first_seen_together*pair_type, family=poisson, data=bs)
m4 <- glm (clutchsz ~  laydate + year+first_seen_together+pair_type, family=poisson, data=bs)
m5 <- glm (clutchsz ~  laydate + year*first_seen_together+pair_type, family=poisson, data=bs)
m6 <- glm (clutchsz ~ laydate +  first_seen_together+year*pair_type, family=poisson, data=bs)
m7 <- glm (clutchsz ~  laydate + first_seen_together*year, family=poisson, data=bs)
m8 <- glm (clutchsz ~  laydate + first_seen_together+year, family=poisson, data=bs)
m9 <- glm (clutchsz ~  laydate + first_seen_together, family=poisson, data=bs)
m10 <- glm (clutchsz ~  laydate + year, family=poisson, data=bs)
m11 <- glm (clutchsz ~  laydate + first_seen_together*pair_type, family=poisson, data=bs)
m12 <- glm (clutchsz ~  laydate + first_seen_together+pair_type, family=poisson, data=bs)
m13 <- glm (clutchsz ~  laydate + year*pair_type, family=poisson, data=bs)
m14 <- glm (clutchsz ~  laydate + year+pair_type, family=poisson, data=bs)
m15 <- glm (clutchsz ~  laydate + first_seen_together*year + year*pair_type, family=poisson, data=bs)
m16 <- glm (clutchsz ~  laydate + pair_type*year + first_seen_together*pair_type, family=poisson, data=bs)

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16) 

basic <- glm (clutchsz ~  laydate, family=quasipoisson, data=bs)
summary(basic)

# polinomial

m8_q2 <- glm (clutchsz ~  laydate + first_seen_together + I(first_seen_together^2) + year, family=poisson, data=bs)
m8_q3 <- glm (clutchsz ~  laydate + first_seen_together + I(first_seen_together^2) + 
                 I(first_seen_together^3) +year, family=poisson, data=bs)

AIC(m8, m8_q2, m8_q3)

### m.f arival time

F1 <- glm (clutchsz ~  laydate + f_first_seen + year , family=poisson, data=bs)
M1 <- glm (clutchsz ~  laydate  + m_first_seen + year, family=poisson, data=bs)
F2 <- glm (clutchsz ~  laydate + first_seen_together + f_first_seen, family=poisson, data=bs)
M2 <- glm (clutchsz ~  laydate + m_first_seen + first_seen_together, family=poisson, data=bs)

AIC(m8, F1, F2, M1, M2) # m18 and m19 have simmilar support


################################################################################################
################################ Number of chicks
################################################################################################


basic <- glm (n_chicks ~  clutchsz, family=poisson, data=bs)
m1 <- glm (n_chicks ~  clutchsz + pair_type, family=poisson, data=bs)
m2 <- glm (n_chicks ~  clutchsz + first_seen_together*year + first_seen_together*pair_type, family=poisson, data=bs)
m3 <- glm (n_chicks ~  clutchsz + year+first_seen_together*pair_type, family=poisson, data=bs)
m4 <- glm (n_chicks ~  clutchsz + year+first_seen_together+pair_type, family=poisson, data=bs)
m5 <- glm (n_chicks ~  clutchsz + year*first_seen_together+pair_type, family=poisson, data=bs)
m6 <- glm (n_chicks ~ clutchsz +  first_seen_together+year*pair_type, family=poisson, data=bs)
m7 <- glm (n_chicks ~  clutchsz + first_seen_together*year, family=poisson, data=bs)
m8 <- glm (n_chicks ~  clutchsz + first_seen_together+year, family=poisson, data=bs)
m9 <- glm (n_chicks ~  clutchsz + first_seen_together, family=poisson, data=bs)
m10 <- glm (n_chicks ~  clutchsz + year, family=poisson, data=bs)
m11 <- glm (n_chicks ~  clutchsz + first_seen_together*pair_type, family=poisson, data=bs)
m12 <- glm (n_chicks ~  clutchsz + first_seen_together+pair_type, family=poisson, data=bs)
m13 <- glm (n_chicks ~  clutchsz + year*pair_type, family=poisson, data=bs)
m14 <- glm (n_chicks ~  clutchsz + year+pair_type, family=poisson, data=bs)
m15 <- glm (n_chicks ~  clutchsz + first_seen_together*year + year*pair_type, family=poisson, data=bs)
m16 <- glm (n_chicks ~  clutchsz + pair_type*year + first_seen_together*pair_type, family=poisson, data=bs)

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16) #


basic <- glm (n_chicks ~  clutchsz, family=quasipoisson, data=bs)

summary(basic)


#######################################
############# n fledged
#######################################

basic <- glm (n_fledged ~  n_chicks, family=poisson, data=bs)
m1 <- glm (n_fledged ~  n_chicks + pair_type, family=poisson, data=bs)
m2 <- glm (n_fledged ~  n_chicks + first_seen_together*year + 
              first_seen_together*pair_type, family=poisson, data=bs)
m3 <- glm (n_fledged ~  n_chicks + year+first_seen_together*pair_type, family=poisson, data=bs)
m4 <- glm (n_fledged ~  n_chicks + year+first_seen_together+pair_type, family=poisson, data=bs)
m5 <- glm (n_fledged ~  n_chicks + year*first_seen_together+pair_type, family=poisson, data=bs)
m6 <- glm (n_fledged ~ n_chicks +  first_seen_together+year*pair_type, family=poisson, data=bs)
m7 <- glm (n_fledged ~  n_chicks + first_seen_together*year, family=poisson, data=bs)
m8 <- glm (n_fledged ~  n_chicks + first_seen_together+year, family=poisson, data=bs)
m9 <- glm (n_fledged ~  n_chicks + first_seen_together, family=poisson, data=bs)
m10 <- glm (n_fledged ~  n_chicks + year, family=poisson, data=bs)
m11 <- glm (n_fledged ~  n_chicks + first_seen_together*pair_type, family=poisson, data=bs)
m12 <- glm (n_fledged ~  n_chicks + first_seen_together+pair_type, family=poisson, data=bs)
m13 <- glm (n_fledged ~  n_chicks + year*pair_type, family=poisson, data=bs)
m14 <- glm (n_fledged ~  n_chicks + year+pair_type, family=poisson, data=bs)
m15 <- glm (n_fledged ~  n_chicks + first_seen_together*year + 
               year*pair_type, family=poisson, data=bs)
m16 <- glm (n_fledged ~  n_chicks + pair_type*year + 
               first_seen_together*pair_type, family=poisson, data=bs)

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15,m16) 

basic <- glm (n_fledged ~  n_chicks, family=quasipoisson, data=bs)
summary(basic)

# polinomail

m9_q2 <- glm (n_fledged ~  n_chicks+first_seen_together + 
                 I(first_seen_together^2) + year, family=poisson, data=bs)
m9_q3 <- glm (n_fledged ~  n_chicks +first_seen_together + 
                 I(first_seen_together^2) + I(first_seen_together^3) + year, family=poisson, data=bs)

AIC(m9, m9_q2, m9_q3) # no support for polinomials

#
# F/M arival time
#

F1 <- glm (n_fledged ~  n_chicks +f_first_seen + year, family=poisson, data=bs)
M1 <- glm (n_fledged ~  n_chicks +m_first_seen + year, family=poisson, data=bs)
F2 <- glm (n_fledged ~  n_chicks +first_seen_together+ f_first_seen + year, family=poisson, data=bs)
M2 <- glm (n_fledged ~  n_chicks +first_seen_together + year + m_first_seen, family=poisson, data=bs)

AIC(m9, M1, M2, F1, F2) # 

summary(M1)

############# n fledged, fledged more than 1

t2 <- subset(bs, bs$n_fledged > "0")

basic <- glm (n_fledged ~  n_chicks, family=poisson, data=t2)
m1 <- glm (n_fledged ~  n_chicks + pair_type, family=poisson, data=t2)
m2 <- glm (n_fledged ~  n_chicks + first_seen_together*year + first_seen_together*pair_type, family=poisson, data=t2)
m3 <- glm (n_fledged ~  n_chicks + year+first_seen_together*pair_type, family=poisson, data=t2)
m4 <- glm (n_fledged ~  n_chicks + year+first_seen_together+pair_type, family=poisson, data=t2)
m5 <- glm (n_fledged ~  n_chicks + year*first_seen_together+pair_type, family=poisson, data=t2)
m6 <- glm (n_fledged ~ n_chicks +  first_seen_together+year*pair_type, family=poisson, data=t2)
m7 <- glm (n_fledged ~  n_chicks + first_seen_together*year, family=poisson, data=t2)
m8 <- glm (n_fledged ~  n_chicks + first_seen_together+year, family=poisson, data=t2)
m9 <- glm (n_fledged ~  n_chicks + first_seen_together, family=poisson, data=t2)
m10 <- glm (n_fledged ~  n_chicks + year, family=poisson, data=t2)
m11 <- glm (n_fledged ~  n_chicks + first_seen_together*pair_type, family=poisson, data=t2)
m12 <- glm (n_fledged ~  n_chicks + first_seen_together+pair_type, family=poisson, data=t2)
m13 <- glm (n_fledged ~  n_chicks + year*pair_type, family=poisson, data=t2)
m14 <- glm (n_fledged ~  n_chicks + year+pair_type, family=poisson, data=t2)

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) # intercept it the best supported

summary(basic)


#########################################################
################### hatching success
#########################################################

basic <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  1, family=binomial, data=bs)
basic1 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  laydate, family=binomial, data=bs)
basic2 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  laydate + clutchsz, family=binomial, data=bs)
basic3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  clutchsz, family=binomial, data=bs)

AIC(basic, basic1, basic2, basic3) # equivalent support to all

# control for laydate

m1 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  pair_type+ laydate, family=binomial, data=bs)
m2 <- glm (cbind(n_chicks, (clutchsz - n_chicks) )~  first_seen_together*year+ 
             first_seen_together*pair_type + laydate, family=binomial, data=bs)
m3 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together*pair_type+
             year+ laydate, family=binomial, data=bs)
m4 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+year+
             pair_type+ laydate, family=binomial, data=bs)
m5 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  first_seen_together*year+
             pair_type+ laydate , family=binomial, data=bs)
m6 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+
             year*pair_type+ laydate , family=binomial, data=bs)
m7 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together*year+ 
             laydate, family=binomial, data=bs)
m8 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+
             year+ laydate, family=binomial, data=bs)
m9 <-glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+ laydate, family=binomial, data=bs)
m10 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~   year+ laydate, family=binomial, data=bs)
m11 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~first_seen_together*pair_type+ 
              laydate, family=binomial, data=bs)
m12 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+pair_type+ 
              laydate, family=binomial, data=bs)
m13 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~ year*pair_type+ laydate, family=binomial, data=bs)
m14 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~ year+pair_type+ laydate, family=binomial, data=bs)
m15 <- glm (cbind(n_chicks, (clutchsz - n_chicks) )~  first_seen_together*year+ 
              year*pair_type + laydate, family=binomial, data=bs)
m16 <- glm (cbind(n_chicks, (clutchsz - n_chicks) )~  pair_type*year+ 
              first_seen_together*pair_type + laydate, family=binomial, data=bs)

AIC(basic1, m1, m2, m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16) # m4 & m5 best supported


# polinomial

m4_q2 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+year+ laydate + pair_type +
                I(first_seen_together^2), family=binomial, data=bs)
m4_q3 <- glm (cbind(n_chicks, (clutchsz - n_chicks))  ~  first_seen_together+year+ pair_type +
                I(first_seen_together^2) + I(first_seen_together^3)+laydate, family=binomial, data=bs)
m5_q2 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  first_seen_together*year+I(first_seen_together^2):year
              + pair_type+ laydate , family=binomial, data=bs)
m5_q3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  first_seen_together*year+
                pair_type+ laydate +I(first_seen_together^2):year + I(first_seen_together^3):year, family=binomial, data=bs)

AIC(m5, m5_q2, m5_q3, m4, m4_q2, m4_q3) # m_5_q2 amd q3 have best and equvialent support



# f / m arrival time

F1 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  f_first_seen + first_seen_together*year+I(first_seen_together^2):year
           + pair_type+ laydate , family=binomial, data=bs)
F2 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  first_seen_together+ f_first_seen*year+I(f_first_seen^2):year
           + pair_type+ laydate , family=binomial, data=bs)
F3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  f_first_seen*year+I(f_first_seen^2):year
           + pair_type+ laydate , family=binomial, data=bs)
M1 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  m_first_seen + first_seen_together*year+I(first_seen_together^2):year
           + pair_type+ laydate , family=binomial, data=bs)
M2 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  first_seen_together+ m_first_seen*year+I(m_first_seen^2):year
           + pair_type+ laydate , family=binomial, data=bs)
M3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  m_first_seen*year+I(m_first_seen^2):year
           + pair_type+ laydate , family=binomial, data=bs)

AIC(M3, M2, M1, F3, F2, F1, m5_q2) # 

M3 <- glm (cbind(n_chicks, (clutchsz - n_chicks)) ~  m_first_seen*year+I(m_first_seen^2):year
           + pair_type+ laydate , family=quasibinomial, data=bs)

summary(M3)

#######################################################################
########### Fledging success
######################################################################


#####################
##################### fl succ proportion
#####################

basic <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate, family=binomial, data=bs)
basic2 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate + clutchsz, family=binomial, data=bs)
basic3 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  clutchsz, family=binomial, data=bs)

AIC(basic, basic2, basic3) # contol for both cutch and laydate


##################### fl succ controled for laydate and clutchsize

m1 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ pair_type, family=binomial, data=bs)
m2 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together*year + 
              first_seen_together*pair_type, family=binomial, data=bs)
m3 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + 
              clutchsz+ year+first_seen_together*pair_type, family=binomial, data=bs)
m4 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + 
              clutchsz+ year+first_seen_together+pair_type, family=binomial, data=bs)
m5 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + 
              clutchsz+ year*first_seen_together+pair_type, family=binomial, data=bs)
m6 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~ laydate  + clutchsz+  
              first_seen_together+year*pair_type, family=binomial, data=bs)
m7 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ 
              first_seen_together*year, family=binomial, data=bs)
m8 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together+year, family=binomial, data=bs)
m9 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together, family=binomial, data=bs)
m10 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate + clutchsz + year, family=binomial, data=bs)
m11 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ 
               first_seen_together*pair_type, family=binomial, data=bs)
m12 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ 
               first_seen_together+pair_type, family=binomial, data=bs)
m13 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ 
               year*pair_type, family=binomial, data=bs)
m14 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year+
               pair_type, family=binomial, data=bs)
m15 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ 
               first_seen_together*year + year*pair_type, family=binomial, data=bs)
m16 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ pair_type*year + first_seen_together*pair_type, family=binomial, data=bs)

AIC(basic2, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16) # best m5

#polinomial using  m5

m5_q2 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*first_seen_together+
                 I(first_seen_together^2):year+pair_type, family=binomial, data=bs)
m5_q3 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*first_seen_together+
                 I(first_seen_together^2):year+I(first_seen_together^3):year + 
                 pair_type, family=binomial, data=bs)


AIC(m5, m5_q2, m5_q3) # m5_q2 is convincingly the best


# female and male arival times


F1 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*f_first_seen+
              I(f_first_seen^2):year+pair_type + first_seen_together, family=binomial, data=bs)
F2 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*f_first_seen+
              I(f_first_seen^2):year+pair_type, family=binomial, data=bs)
F3 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~ f_first_seen +  laydate  + clutchsz+ year*first_seen_together+
              I(first_seen_together^2):year+pair_type, family=binomial, data=bs)
M1 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*m_first_seen+
              I(m_first_seen^2):year+pair_type + first_seen_together, family=binomial, data=bs)
M2 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*m_first_seen+
              I(m_first_seen^2):year+pair_type, family=binomial, data=bs)
M3 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  m_first_seen + laydate  + clutchsz+ year*first_seen_together+
              I(first_seen_together^2):year+pair_type, family=binomial, data=bs)

AIC(m5_q2, F1, F2, F3, M3, M2, M1) 

F1 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*f_first_seen+
              I(f_first_seen^2):year+pair_type + first_seen_together, family=quasibinomial, data=bs)

summary(F1) # meeting time is NS when overdispersion accounted for


##################### fl succ, only those fledged mor than 0, controled for laydate and clutchsize


basic <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate + clutchsz, family=binomial, data=t2)
m1 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ pair_type, family=binomial, data=t2)
m2 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together*year + first_seen_together*pair_type, family=binomial, data=t2)
m3 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year+first_seen_together*pair_type, family=binomial, data=t2)
m4 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ 
              year+first_seen_together+pair_type, family=binomial, data=t2)
m5 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*first_seen_together+pair_type, family=binomial, data=t2)
m6 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~ laydate  + clutchsz+  first_seen_together+year*pair_type, family=binomial, data=t2)
m7 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together*year, family=binomial, data=t2)
m8 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together+year, family=binomial, data=t2)
m9 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together, family=binomial, data=t2)
m10 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate + clutchsz + year, family=binomial, data=t2)
m11 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together*pair_type, family=binomial, data=t2)
m12 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ first_seen_together+pair_type, family=binomial, data=t2)
m13 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year*pair_type, family=binomial, data=t2)
m14 <- glm (cbind(n_fledged, (n_chicks - n_fledged) ) ~  laydate  + clutchsz+ year+pair_type, family=binomial, data=t2)
m15 <- glm (succ_failure ~  laydate + first_seen_together*year + 
               year*pair_type, family=binomial, data=t2) # does not converge
m16 <- glm (succ_failure ~  laydate + pair_type*year + 
               first_seen_together*pair_type, family=binomial, data=t2) # does not converge

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14) # several models have the same support
#m13, m14, m3
summary(m4)

##### success vs failure

basic <- glm (succ_failure ~  laydate, family=binomial, data=bs)
m1 <- glm (succ_failure ~  laydate + pair_type, family=binomial, data=bs)
m2 <- glm (succ_failure ~  laydate + first_seen_together*year + first_seen_together*pair_type, family=binomial, data=bs)
m3 <- glm (succ_failure ~  laydate + year+first_seen_together*pair_type, family=binomial, data=bs)
m4 <- glm (succ_failure ~  laydate + year+first_seen_together+pair_type, family=binomial, data=bs)
m5 <- glm (succ_failure ~  laydate + year*first_seen_together+pair_type, family=binomial, data=bs)
m6 <- glm (succ_failure ~ laydate +  first_seen_together+year*pair_type, family=binomial, data=bs)
m7 <- glm (succ_failure ~  laydate + first_seen_together*year, family=binomial, data=bs)
m8 <- glm (succ_failure ~  laydate + first_seen_together+year, family=binomial, data=bs)
m9 <- glm (succ_failure ~  laydate + first_seen_together, family=binomial, data=bs)
m10 <- glm (succ_failure ~  laydate + year, family=binomial, data=bs)
m11 <- glm (succ_failure ~  laydate + first_seen_together*pair_type, family=binomial, data=bs)
m12 <- glm (succ_failure ~  laydate + first_seen_together+pair_type, family=binomial, data=bs)
m13 <- glm (succ_failure ~  laydate + year*pair_type, family=binomial, data=bs)
m14 <- glm (succ_failure ~  laydate + year+pair_type, family=binomial, data=bs)
m15 <- glm (succ_failure ~  laydate + first_seen_together*year + year*pair_type, family=binomial, data=bs)
m16 <- glm (succ_failure ~  laydate + pair_type*year + first_seen_together*pair_type, family=binomial, data=bs)

AIC(basic, m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12, m13, m14, m15 ,m16) 


# polinomial

m8_q2 <- glm (succ_failure ~  laydate + first_seen_together+year +I(first_seen_together^2), family=binomial, data=bs)
m8_q3 <- glm (succ_failure ~  laydate + first_seen_together+year +I(first_seen_together^2)  +I(first_seen_together^3), family=binomial, data=bs)

AIC(m8, m8_q2, m8_q3) # m8 has the best support

# F/M arrival time

F1 <- glm (succ_failure ~  laydate + first_seen_together+year + f_first_seen, family=binomial, data=bs)
F2 <- glm (succ_failure ~  laydate +year + f_first_seen, family=binomial, data=bs)
M1 <- glm (succ_failure ~  laydate + first_seen_together+year + m_first_seen, family=binomial, data=bs)
M2 <- glm (succ_failure ~  laydate +year + m_first_seen, family=binomial, data=bs)

AIC(M2,F2,M1, F1)



########################################################################################
########################################################################################
########################################################################################

############# SECTION 2

# Meeting time (in year prior to the breeding seson t) and DIVORCE (between breeding seasons t and t+1)

########################################################################################
########################################################################################
########################################################################################

###################### laoding the dataset - only new pairs that have been seen associating. also, there was on pair for 2011/12 for which laydate
###################### was unknown and that one is excluded from the dataset

divorce <- read.table(file.choose(), header = T, sep = "")

divorce$year <- as.factor(divorce$year)

x <- table(divorce$status_t.1, divorce$year) #Table 2
x

str(divorce)
######################2007to2009 years

div_1 <- subset(divorce, divorce$period =="2007to2010")
str(div_1)


m1 <- glm(status_bin ~ 1, data=div_1, family = "binomial")
m2 <- glm(status_bin ~ first_seen_together, data=div_1, family = "binomial")
m3 <- glm(status_bin ~ clutchsz, data=div_1, family = "binomial")
m4 <- glm(status_bin ~ first_seen_together+clutchsz, data=div_1, family = "binomial")
m5 <- glm(status_bin ~ first_seen_together+laydate, data=div_1, family = "binomial")
m6 <- glm(status_bin ~ laydate, data=div_1, family = "binomial")

AIC(m1,m2,m3,m4,m5,m6)

summary(m5)

############################### 2011

div_2 <- subset(divorce, divorce$period =="2011to2014")
str(div_2)

m1 <- glm(status_bin ~ 1, data=div_2, family = "binomial")
m2 <- glm(status_bin ~ first_seen_together, data=div_2, family = "binomial")
m3 <- glm(status_bin ~ clutchsz, data=div_2, family = "binomial")
m4 <- glm(status_bin ~ first_seen_together+clutchsz, data=div_2, family = "binomial")
m5 <- glm(status_bin ~ first_seen_together+laydate, data=div_2, family = "binomial")
m6 <- glm(status_bin ~ laydate, data=div_2, family = "binomial")

m7 <- glm(status_bin ~ year, data=div_2, family = "binomial")
m8 <- glm(status_bin ~ first_seen_together + year, data=div_2, family = "binomial")
m9 <- glm(status_bin ~ clutchsz + year, data=div_2, family = "binomial")
m10 <- glm(status_bin ~ first_seen_together+clutchsz + year, data=div_2, family = "binomial")
m11 <- glm(status_bin ~ first_seen_together+ laydate + year, data=div_2, family = "binomial")
m12 <- glm(status_bin ~ laydate + year, data=div_2, family = "binomial")

AIC(m1,m2,m3,m4,m5,m6, m7, m8, m9, m10, m11, m12)

