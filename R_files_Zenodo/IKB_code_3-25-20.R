##### IKB Gene expression from laboratory experiment 2015 ####

# 3-25-20

data<-read.csv("IKB_data_3-25-20.csv")
head(data)
length(data$sample) # 80 - 10 colonies*4 treatments* 2 timepoints

## Normality 

# IKB referred to as nfkb for shorthand in this file (really an inhibitor of NF-ÎºB - D. melanogaster) 

# nkfb (IKB)
hist(data$nfkb_ab) # not normal
qqnorm(data$nfkb_ab)
qqline(data$nfkb_ab)  

# Transform
qqnorm(sqrt(data$nfkb_ab))
qqline(sqrt(data$nfkb_ab))  

shapiro.test(sqrt(data$nfkb_ab)) # Normal when square root transformed

## Data classes
class(data$sample) # factor
class(data$treatment) # factor
data$colony<-as.factor(data$colony)
class(data$colony) # factor
class(data$time) #

class(data$nfkb_ab)

class(data$pathYN)
class(data$trtmt_M)
class(data$trtmt_A)
class(data$trtmt_MA)
class(data$onepath)

## Contrasts
contrasts(data$treatment) 
contrasts(data$treatment) <- contr.sum # contrasts the control to each treatment
contrasts(data$treatment) 

contrasts(data$colony) <- contr.sum
contrasts(data$time) <- cbind(c(-1,1))

contrasts(data$pathYN) <- cbind(c(-1,1))
contrasts(data$pathYN)
contrasts(data$trtmt_M) <- cbind(c(-1,1))
contrasts(data$trtmt_M) 
contrasts(data$trtmt_A) <- cbind(c(-1,1))
contrasts(data$trtmt_A) 
contrasts(data$trtmt_MA) <- cbind(c(1,0,-1),c(-1,1,0))
contrasts(data$trtmt_MA) 

## Exploratory plots
plot(colony~sample, data=data) # 1 each - good
plot(time~sample, data=data) 
plot(treatment~colony,data=data)
plot(treatment~time,data=data)
plot(colony~time,data=data)

# Response variable: nfkb_AB
sqrt_nfkb_AB<-sqrt(data$nfkb_ab)

# Homoscedasticity:
plot(sqrt_nfkb_AB~treatment, data=data) 
plot(sqrt_nfkb_AB~colony, data=data) # differences 
plot(sqrt_nfkb_AB~time, data=data) 
plot(sqrt_nfkb_AB~pathYN, data=data) 
plot(sqrt_nfkb_AB~trtmt_M, data=data) 
plot(sqrt_nfkb_AB~trtmt_A, data=data) 
plot(sqrt_nfkb_AB~trtmt_MA, data=data) 


########### Without colony 9  #######

data_no9<-subset(data,colony!="9")
summary(data_no9$colony)
length(data_no9$time) # 72 

hist(data_no9$nfkb_ab) # not normal

qqnorm(data_no9$nfkb_ab)
qqline(data_no9$nfkb_ab)  

# Transform
qqnorm(sqrt(data_no9$nfkb_ab))
qqline(sqrt(data_no9$nfkb_ab)) 

shapiro.test(sqrt(data_no9$nfkb_ab)) # Normal 

## Data classes
class(data_no9$sample) # factor
class(data_no9$treatment) # factor
data_no9$colony<-as.factor(data_no9$colony)
class(data_no9$colony) # factor
class(data_no9$time) #

class(data_no9$nfkb_ab)

class(data_no9$onepath)

## Contrasts
contrasts(data_no9$treatment) <- contr.sum
contrasts(data_no9$treatment)<-cbind(c(-1,0,0,1),c(-1,1,0,0),c(-1,0,1,0)) # relative to N
contrasts(data_no9$colony) <- contr.sum 
contrasts(data_no9$time) <- contr.sum 
contrasts(data_no9$onepath) <- cbind(c(1,0,-1),c(0,1,-1))

## Exploratory plots
plot(colony~sample, data=data_no9)
plot(time~sample, data=data_no9) 
plot(treatment~colony,data=data_no9)
plot(treatment~time,data=data_no9)
plot(colony~time,data=data_no9)
plot(onepath~treatment, data=data_no9)

# Response variable: nfkb_AB
sqrt_nfkb_AB<-sqrt(data_no9$nfkb_ab)
plot(sqrt_nfkb_AB~treatment, data=data_no9) 
plot(sqrt_nfkb_AB~colony, data=data_no9) #  
plot(sqrt_nfkb_AB~time, data=data_no9) # 

###################
##### Linear models ## 
###################

# Fixed effects: timing, treatments
# Random effects: colony and test for random slope 

library(lme4)

### Deciding on random slopes 
#1) Fit a global model with all predictors as both fixed effects and as random slopes. 
#2) See which random slopes are worth including, then use for testing for fixed effects

###### EF1 #########

ef_whole<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 
# test onepath separately because otherwise does not converge

# no pathYN as random slope
ef_b<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_M
ef_c<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_A
ef_d<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+time |colony),REML=FALSE,data=data_no9) 

# no time
ef_e<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A |colony),REML=FALSE,data=data_no9) 

# none
ef_f<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(1 |colony),REML=FALSE,data=data_no9) 

# onepath
ef_g<-lmer(ef1~ onepath+time +(1 |colony),REML=FALSE,data=data_no9) 
ef_g2<-lmer(ef1~onepath+time +(onepath |colony),REML=FALSE,data=data_no9) 
library(lmtest)
lrtest(ef_g,ef_g2) # not better with onepath as random slope

# LRT 
library(lmtest)
lrtest(ef_whole,ef_b) # pathYN not worth including
lrtest(ef_whole,ef_c) # trtmt_M not worth it
lrtest(ef_whole,ef_d) # trtmt_A not worth it
lrtest(ef_whole,ef_e) # time not worth it 

# All NS vs. one with no random slopes: none are needed then b/c all NS
lrtest(ef_f,ef_whole)
lrtest(ef_f,ef_b)
lrtest(ef_f,ef_c)
lrtest(ef_f,ef_d)
lrtest(ef_f,ef_e)

ef_f2<-lm(ef1~trtmt_M+trtmt_A + pathYN + time,data=data_no9) 

lrtest(ef_f,ef_f2) # colony as random intercept is worthwhile 

# None worth including as random slopes for EF1, reference gene 
# Use 1|colony for random intercept in EF1 models

### EF1 models #######

efnull<-lmer(ef1~1+(1|colony),REML=FALSE,data=data_no9) 

# M&A 

ef0<-lmer(ef1~trtmt_A*time+(1|colony),REML=FALSE,data=data_no9) 

ef1<-lmer(ef1~trtmt_M*time+(1|colony),REML=FALSE,data=data_no9)

ef2<-lmer(ef1~trtmt_M*trtmt_A+(1|colony),REML=FALSE,data=data_no9)

ef3<-lmer(ef1~trtmt_M*trtmt_A+time+(1|colony),REML=FALSE,data=data_no9)

ef4<-lmer(ef1~trtmt_M*trtmt_A*time+(1|colony),REML=FALSE,data=data_no9) 

ef5<-lmer(ef1~trtmt_M+trtmt_A*time+(1|colony),REML=FALSE,data=data_no9)

ef6<-lmer(ef1~trtmt_M*time+trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

ef7<-lmer(ef1~trtmt_M+time+trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

ef8<-lmer(ef1~trtmt_M+trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

ef9<-lmer(ef1~trtmt_A+time+(1|colony),REML=FALSE,data=data_no9) 

ef10<-lmer(ef1~trtmt_M+time+(1|colony),REML=FALSE,data=data_no9) 

ef11<-lmer(ef1~trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

ef12<-lmer(ef1~trtmt_M+(1|colony),REML=FALSE,data=data_no9) 

ef13<-lmer(ef1~pathYN*time+(1|colony),REML=FALSE,data=data_no9) 

ef14<-lmer(ef1~pathYN+time+(1|colony),REML=FALSE,data=data_no9) 

ef15<-lmer(ef1~pathYN+(1|colony),REML=FALSE,data=data_no9) 

ef16<-lmer(ef1~onepath+(1|colony),REML=FALSE,data=data_no9) 

ef17<-lmer(ef1~onepath+time+(1|colony),REML=FALSE,data=data_no9) 

ef18<-lmer(ef1~onepath*time+(1|colony),REML=FALSE,data=data_no9) 


# AIC table
AIC(efnull,ef0,ef1,ef2,ef3,ef4,ef5,ef6,ef7,ef8,ef9,ef10,ef11,ef12,ef13,ef14,ef15,ef16,ef17,ef18)

# Test whether models are better:13,18,0,4,5,1,14,10,9,17,6
lrtest(efnull,ef13) # 13 is better
lrtest(efnull,ef18) # 18 is better
lrtest(efnull,ef4) # 4 is better
lrtest(efnull,ef0) # 0 is better
lrtest(efnull,ef5) # 5 is better
lrtest(efnull,ef1) # 1 is better
lrtest(efnull,ef14) # null is better
lrtest(efnull,ef10) # null is better
lrtest(efnull,ef9) # null is better

# check summary for ones where null is NOT better
library(lmerTest) # adds p value

# 4
ef4_reml<-lmer(ef1~trtmt_M*trtmt_A*time+(1|colony),REML=TRUE,data=data_no9) 
summary(ef4_reml) # all NS

# 5
ef5_reml<-lmer(ef1~trtmt_M+trtmt_A*time+(1|colony),REML=TRUE,data=data_no9)
summary(ef5_reml) # all NS

# 0
ef0_reml<-lmer(ef1~trtmt_A*time+(1|colony),REML=TRUE,data=data_no9) 
summary(ef0_reml) # all NS

# 13
ef13_reml<-lmer(ef1~pathYN*time+(1|colony),REML=TRUE,data=data_no9) 
summary(ef13_reml) # all NS

# 18
ef18_reml<-lmer(ef1~onepath*time+(1|colony),REML=TRUE,data=data_no9) 
summary(ef18_reml) # all NS

# 1
ef1_reml<-lmer(ef1~trtmt_M*time+(1|colony),REML=TRUE,data=data_no9)
summary(ef1_reml) # all NS


###### NFKB #########

# Erika:
#1) Random slopes decision 
#2) Fixed effects models

# Random slopes
nfkb_whole<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) # fixed effect model matrix is rank deficient if include M*A in this model
# test onepath models separately for convergence, too

# no pathYN as random slope
nfkb_b<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_M
nfkb_c<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_A
nfkb_d<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+time |colony),REML=FALSE,data=data_no9) 

# no time
nfkb_e<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A |colony),REML=FALSE,data=data_no9) 

# none
nfkb_f<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(1 |colony),REML=FALSE,data=data_no9) 

# LRT 
library(lmtest)
lrtest(nfkb_whole,nfkb_b) # whole not better 
lrtest(nfkb_whole,nfkb_c) # whole not better
lrtest(nfkb_whole,nfkb_d) # whole not better
lrtest(nfkb_whole,nfkb_e) # whole not better

# All actually better than only colony?
lrtest(nfkb_f,nfkb_b) # no (omits pathYN)
lrtest(nfkb_f,nfkb_c) # yes (omits trtmt_M)
lrtest(nfkb_f,nfkb_d) # yes (omits trtmt_A)
lrtest(nfkb_f,nfkb_e) # almost (omits time)

# pathYN and time combos supported most
nfkb_g<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(pathYN |colony),REML=FALSE,data=data_no9) 

nfkb_h<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(time |colony),REML=FALSE,data=data_no9) 

nfkb_i<-lmer(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time +(pathYN+time|colony),REML=FALSE,data=data_no9) 

# LRT - better than colony alone?
lrtest(nfkb_f,nfkb_g) # yes: pathYN 
lrtest(nfkb_f,nfkb_h) # no: time
lrtest(nfkb_f,nfkb_i) # yes: pathYN+time
lrtest(nfkb_g,nfkb_i) # no - just pathYN is BEST

# c and d are not better
lrtest(nfkb_g,nfkb_c) 
lrtest(nfkb_g,nfkb_d)

# AIC
AIC(nfkb_g,nfkb_c,nfkb_d,nfkb_i,nfkb_f) # g is best with only pathYN as a random slope

nfkb_g2<-lm(sqrt_nfkb_AB~trtmt_M+trtmt_A + pathYN + time, data=data_no9) 

lrtest(nfkb_g,nfkb_g2) 
lrtest(nfkb_f,nfkb_g2) # colony worth including as intercept 

# onepath
nfkb_j<-lmer(sqrt_nfkb_AB~ onepath+ time+(1 |colony),REML=FALSE,data=data_no9) 
nfkb_j2<-lmer(sqrt_nfkb_AB~ onepath+ time+(onepath |colony),REML=FALSE,data=data_no9) 
lrtest(nfkb_j,nfkb_j2) # onepath random slope is worthwhile

# PathYN and onepath both supported as random slopes


### Test fixed effects: pathYN random slope ####
# models must have pathYN because it's a random slope 

nfkb_null<-lmer(sqrt_nfkb_AB~1+(pathYN|colony),REML=FALSE,data=data_no9) 

nfkb_1<-lmer(sqrt_nfkb_AB~pathYN+(pathYN|colony),REML=FALSE,data=data_no9)

nfkb_2<-lmer(sqrt_nfkb_AB~pathYN+time+(pathYN|colony),REML=FALSE,data=data_no9)

nfkb_3<-lmer(sqrt_nfkb_AB~pathYN*time+(pathYN|colony),REML=FALSE,data=data_no9)

nfkb_4<-lmer(sqrt_nfkb_AB~pathYN*time+trtmt_M+trtmt_A+(pathYN|colony),REML=FALSE,data=data_no9)
drop1(nfkb_4) # drop trtmt_M and trtmt_A 

nfkb_5<-lmer(sqrt_nfkb_AB~pathYN+time+trtmt_M+trtmt_A+(pathYN|colony),REML=FALSE,data=data_no9) 
drop1(nfkb_5) # drop trtmt_A
nfkb_5b<-lmer(sqrt_nfkb_AB~pathYN+time+trtmt_M+(pathYN|colony),REML=FALSE,data=data_no9)
drop1(nfkb_5b) # drop trtmt-M

nfkb_6<-lmer(sqrt_nfkb_AB~pathYN+trtmt_M*trtmt_A+time+(pathYN|colony),REML=FALSE,data=data_no9) # 
drop1(nfkb_6) # drop time
nfkb_6b<-lmer(sqrt_nfkb_AB~pathYN+trtmt_M*trtmt_A+(pathYN|colony),REML=FALSE,data=data_no9) # 
drop1(nfkb_6b)

nfkb_7<-lmer(sqrt_nfkb_AB~pathYN+trtmt_M*time+(pathYN|colony),REML=FALSE,data=data_no9) # 

nfkb_8<-lmer(sqrt_nfkb_AB~pathYN+trtmt_A*time+(pathYN|colony),REML=FALSE,data=data_no9) # 

# AIC
AIC(nfkb_null,nfkb_1,nfkb_2,nfkb_3,nfkb_4,nfkb_5,nfkb_5b,nfkb_6b,nfkb_7,nfkb_8)

# Only pathYN models
AIC(nfkb_null,nfkb_1,nfkb_2,nfkb_3)


# Best model = 3 by AIC
library(lmtest)
lrtest(nfkb_null,nfkb_3) # null is best by LRT

# REML with 3
nfkb_3_reml<-lmer(sqrt_nfkb_AB~pathYN*time+(pathYN|colony),REML=TRUE,data=data_no9)
summary(nfkb_3_reml) # only time is sig.

# % change pathYN
(0.0837-0.0660)/0.0660 # 26.8%


### Test fixed effects: onepath random slope###

nfkb_null<-lmer(sqrt_nfkb_AB~1+(onepath|colony),REML=FALSE,data=data_no9) 

nfkb_1<-lmer(sqrt_nfkb_AB~onepath+(onepath|colony),REML=FALSE,data=data_no9)

nfkb_2<-lmer(sqrt_nfkb_AB~onepath+time+(onepath|colony),REML=FALSE,data=data_no9)

nfkb_3<-lmer(sqrt_nfkb_AB~onepath*time+(onepath|colony),REML=FALSE,data=data_no9)


# AIC with all
AIC(nfkb_null,nfkb_1,nfkb_2,nfkb_3)

# Best model is null by AIC & LRT
library(lmtest)
lrtest(nfkb_null,nfkb_3) # null is best by LRT


####### PLOTS #####

# REACTION NORMS  
data$treatment <- as.character(data$treatment)
data$treatment <- factor(data$treatment, levels=c("N", "M", "MA","A")) # re-orders

# subset
one_col<-subset(data, data$colony=="1")
two_col<-subset(data, data$colony=="2")
three_col<-subset(data, data$colony=="3")
four_col<-subset(data, data$colony=="4")
five_col<-subset(data, data$colony=="5")
six_col<-subset(data, data$colony=="6")
seven_col<-subset(data, data$colony=="7")
eight_col<-subset(data, data$colony=="8")
nine_col<-subset(data, data$colony=="9")
ten_col<-subset(data, data$colony=="10")

# X axis
x<-c(1:4)

# IKB reaction norm
one_norm<-tapply(one_col$nfkb_ab,one_col$treatment,mean)
two_norm<-tapply(two_col$nfkb_ab,two_col$treatment,mean)
three_norm<-tapply(three_col$nfkb_ab,three_col$treatment,mean)
four_norm<-tapply(four_col$nfkb_ab,four_col$treatment,mean)
five_norm<-tapply(five_col$nfkb_ab,five_col$treatment,mean)
six_norm<-tapply(six_col$nfkb_ab,six_col$treatment,mean)
seven_norm<-tapply(seven_col$nfkb_ab,seven_col$treatment,mean)
eight_norm<-tapply(eight_col$nfkb_ab,eight_col$treatment,mean)
nine_norm<-tapply(nine_col$nfkb_ab,nine_col$treatment,mean)
ten_norm<-tapply(ten_col$nfkb_ab,ten_col$treatment,mean)


##### All 10 colonies

colony_norm<-data.frame(one_norm,two_norm,three_norm,four_norm,five_norm,six_norm,seven_norm,eight_norm,nine_norm,ten_norm)

colony_norm
col_mat<-as.matrix(colony_norm)
col_mat
### Plot
colnames<-c(1:4)
colors<-c("darkred","darkorange","gold","green","blue","purple","lightsalmon4","black","aquamarine3","deeppink")
plot(x, y= rep(-10000,4), ylim= c(0,0.2), # initiaties plot
     ylab= "nfkb_AB", xlab= "Treatment", main= "NFKB by Colony")    
for(i in c(1:10)){
  points(x, col_mat[,i], col=colors[i], type= "b") 
}


############## BAR PLOTS ############
# summarySE - http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/ .....5-23-17
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##### NO COLONY 9########
data_no9$treatment <- as.character(data_no9$treatment)
data_no9$treatment <- factor(data_no9$treatment, levels=c("N", "A", "M","MA")) # re-orders

# Treatment
data_summarize1<-summarySE(data_no9,measurevar="nfkb_ab", groupvars=c("treatment"))

library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize1, aes(x=treatment, y=nfkb_ab,fill="red")) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nfkb_ab-1*se, ymax=nfkb_ab+1*se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) 

# PathYN * time (best by AIC but not LRT)
data_summarize<-summarySE(data,measurevar="nfkb_ab", groupvars=c("time","pathYN"))

library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize, aes(x=time, y=nfkb_ab, fill=pathYN)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nfkb_ab-1*se, ymax=nfkb_ab+1*se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("NFKB expression (corrected to reference gene)")+
  xlab("pathYN")

# PathYN * time with all treatments shown
data_no9$treatment <- factor(data_no9$treatment, levels=c("N", "A", "M","MA")) # re-orders
data_summarize<-summarySE(data_no9,measurevar="nfkb_ab", groupvars=c("time","treatment"))

library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize, aes(x=time, y=nfkb_ab, fill=treatment)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nfkb_ab-1*se, ymax=nfkb_ab+1*se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("NFKB expression (corrected to reference gene)")+
  xlab("Treatment")

