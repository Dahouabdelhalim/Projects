##### T5A Gene expression from laboratory experiment 2015 ####

# Date: 3-23-20

# With efficiency calculations from PCR Miner (see manuscript)

data<-read.csv("ef1_t5A_data_3-23-20.csv")
head(data)
length(data$sample) # 80 - 10 colonies* 4 treatments* 2 times

### EF1 #####
hist(data$ef1)
shapiro.test(data$ef1) # normal

qqnorm(data$ef1)
qqline(data$ef1)  


##### Tachylectin 5A #####

## Normality 
hist(data$t5a_ab)
shapiro.test(data$t5a_ab) # not normal

qqnorm(data$t5a_ab)
qqline(data$t5a_ab)  

# Sqrt transformation
qqnorm(sqrt(data$t5a_ab))
qqline(sqrt(data$t5a_ab))
hist(sqrt(data$t5a_ab))

shapiro.test(sqrt(data$t5a_ab)) # Normal when square root transformed


## Data classes
class(data$sample) # factor
class(data$treatment) # factor
data$colony<-as.factor(data$colony)
class(data$colony) # factor
class(data$time) #

class(data$ef1)
class(data$t5a_ab) #

class(data$trtmt_M)
class(data$trtmt_A)
class(data$trtmt_MA)
class(data$onepath) # "parasite_single" in text
class(data$pathYN)

## Contrasts

contrasts(data$treatment) 
contrasts(data$treatment) <- contr.sum # control vs. each treatment - YES.

contrasts(data$colony) <- contr.sum 
contrasts(data$time) <- cbind(c(-1,1)) # early is lower

contrasts(data$trtmt_M) <- cbind(c(-1,1))
contrasts(data$trtmt_A) <- cbind(c(-1,1))
contrasts(data$onepath) <- cbind(c(1,0,-1),c(0,1,-1))
contrasts(data$onepath)
contrasts(data$pathYN) <- cbind(c(-1,1))

## Check data with plots
plot(colony~sample, data=data) # 1 each 
plot(treatment~sample, data=data) # 4 equal groups
plot(time~sample, data=data) # yes
plot(treatment~colony,data=data)
plot(treatment~time,data=data)
plot(colony~time,data=data)

plot(trtmt_M~treatment,data=data) # 2 groups
plot(trtmt_A~treatment,data=data) # 2 groups
plot(trtmt_MA~treatment,data=data) # 3 groups 

# homoscedasticity: ef1
plot(ef1~treatment, data=data) 
plot(ef1~colony, data=data) # differences by colony 
plot(ef1~time, data=data) 
plot(ef1~trtmt_M,data=data)
plot(ef1~trtmt_A,data=data)
plot(ef1~trtmt_MA,data=data)

# Response variable & homoscedasticity: t5a_AB
sqrt_t5a_ab<-sqrt(data$t5a_ab)
plot(sqrt_t5a_ab~treatment, data=data) 
plot(sqrt_t5a_ab~colony, data=data) # differences
plot(sqrt_t5a_ab~time, data=data) 

plot(t5a_ab~trtmt_M,data=data)
plot(t5a_ab~trtmt_A,data=data)
plot(t5a_ab~trtmt_MA,data=data)

##### Linear models #### 

# Fixed effects: colony, timing, treatment (infection)
# Random: colony as random intercept, random slopes TBD

library(lme4)

##################
# Without colony 9  
#############
data_no9<-subset(data, colony!="9")
summary(data_no9$colony)
length(data_no9$colony) # 72

hist(data_no9$ef1)
shapiro.test(data_no9$ef1) # normal

qqnorm(data_no9$ef1)
qqline(data_no9$ef1)  

class(data_no9$ef1)
class(data_no9$t5a_ab)
class(data_no9$treatment)
class(data_no9$pathYN)
class(data_no9$time)
class(data_no9$trtmt_M)
class(data_no9$trtmt_A)

## Contrasts
contrasts(data_no9$treatment) # default is treatment contrasts
contrasts(data_no9$treatment) <- contr.sum # control vs. each treatment
contrasts(data_no9$treatment) 

contrasts(data_no9$colony) <- contr.sum 

contrasts(data_no9$time) <- cbind(c(-1,1)) # early is lower

contrasts(data_no9$trtmt_M) <- cbind(c(-1,1))
contrasts(data_no9$trtmt_A) <- cbind(c(-1,1))
contrasts(data_no9$trtmt_MA) <- cbind(c(-1,0,1),c(0,1,-1))

contrasts(data_no9$pathYN) <- cbind(c(-1,1))

### Homoscedasticity
# ef1
plot(ef1~treatment, data=data_no9) 
plot(ef1~colony, data=data_no9) # differences by colony 
plot(ef1~time, data=data_no9) 
plot(ef1~trtmt_M,data=data_no9)
plot(ef1~trtmt_A,data=data_no9)
plot(ef1~trtmt_MA,data=data_no9)

# Response variable & homoscedasticity: t5a_AB
sqrt_t5a_ab<-sqrt(data_no9$t5a_ab)
plot(sqrt_t5a_ab~sample, data=data_no9) 

plot(sqrt_t5a_ab~treatment, data=data_no9) 
plot(sqrt_t5a_ab~colony, data=data_no9) # differences
plot(sqrt_t5a_ab~time, data=data_no9) 

plot(t5a_ab~trtmt_M,data=data)
plot(t5a_ab~trtmt_A,data=data)
plot(t5a_ab~trtmt_MA,data=data)


###### EF1: Determining random slopes #####

#1) Fit a  model with all predictors as both fixed effects and as random slopes
#2) See which random slopes are worth including & use for testing fixed effects 

# Whole model split into several because otherwise does not converge

### Whole 1
ef_whole<-lmer(ef1~trtmt_M+trtmt_A + time+(trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmtM as random slope
ef_wholeb<-lmer(ef1~trtmt_M+trtmt_A + time+(trtmt_A+time |colony),REML=FALSE,data=data_no9)

# no trtmtA
ef_wholec<-lmer(ef1~trtmt_M+trtmt_A + time+(trtmt_M+time |colony),REML=FALSE,data=data_no9) 

# no time
ef_wholed<-lmer(ef1~trtmt_M+trtmt_A + time+(trtmt_M+trtmt_A |colony),REML=FALSE,data=data_no9) 

# none
ef_wholee<-lmer(ef1~trtmt_M+trtmt_A + time+(1|colony),REML=FALSE,data=data_no9) 

# LRT 
library(lmtest)
lrtest(ef_whole,ef_wholeb) # trtmt_M not worth it
lrtest(ef_whole,ef_wholec) # trtmt_A not worth it
lrtest(ef_whole,ef_wholed) # time not worth it
lrtest(ef_whole,ef_wholee) # None are worth including as random slopes 

### Whole #2
ef_whole_path<-lmer(ef1~pathYN + time+(pathYN+time |colony),REML=FALSE,data=data_no9) #

# no pathYN as random slope
ef_whole_pathb<-lmer(ef1~pathYN + time+(time |colony),REML=FALSE,data=data_no9)

# no time
ef_whole_pathc<-lmer(ef1~pathYN + time+(pathYN |colony),REML=FALSE,data=data_no9) #

# none
ef_whole_pathd<-lmer(ef1~pathYN + time+(1|colony),REML=FALSE,data=data_no9)

# LRT 
library(lmtest)
lrtest(ef_whole_path,ef_whole_pathb) # pathYN not worth it
lrtest(ef_whole_path,ef_whole_pathc) #time not worth it
lrtest(ef_whole_path,ef_whole_pathd) # none are worth it again


#### Whole #3
ef_whole_onepath<-lmer(ef1~onepath + time+(onepath+time |colony),REML=FALSE,data=data_no9)

# no onepath
ef_whole_onepathb<-lmer(ef1~onepath + time+(time |colony),REML=FALSE,data=data_no9)

# no time
ef_whole_onepathc<-lmer(ef1~onepath + time+(onepath |colony),REML=FALSE,data=data_no9)

# none
ef_whole_onepathd<-lmer(ef1~onepath + time+(1 |colony),REML=FALSE,data=data_no9)
  
  
# LRT 
library(lmtest)
lrtest(ef_whole_onepath,ef_whole_onepathb) # onepath not worth it
lrtest(ef_whole_onepath,ef_whole_onepathc) #time not worth it
lrtest(ef_whole_onepath,ef_whole_onepathd) # none are worth it again


### No factors are worth including as random slopes for EF1, reference gene 
### Use 1|colony for random intercept in EF1 models

### EF1 models: fixed effects ####

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

# model ef13 lowest by AIC
lrtest(efnull,ef13) # Better than the null also by LRT

library(lmerTest) # adds p value
summary(ef13) # p values are not significant: no predictors of ef1, reference gene.


#########################
#### T5a no 9: test allmodels 5  ###
#########################

hist(data_no9$t5a_ab)
qqnorm(sqrt(data_no9$t5a_ab))
qqline(sqrt(data_no9$t5a_ab))

sqrt_t5a9<-sqrt(data_no9$t5a_ab)

shapiro.test(sqrt_t5a9) # good

###### MODELS - determine random effects structure

t5a_whole<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 
# All rank deficient if fit onepath as well- separate model
## onepath
t5a_g<-lmer(sqrt_t5a9~onepath+ time +(1 |colony),REML=FALSE,data=data_no9) 
t5a_g2<-lmer(sqrt_t5a9~onepath+ time +(onepath |colony),REML=FALSE,data=data_no9) 

lrtest(t5a_g,t5a_g2) #  not worth including onepath as random slope

# no pathYN as random slope
t5a_b<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_M
t5a_c<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_A
t5a_d<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+time |colony),REML=FALSE,data=data_no9) 

# no time
t5a_e<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A |colony),REML=FALSE,data=data_no9) 

# none
t5a_f<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(1 |colony),REML=FALSE,data=data_no9) 

# LRT 
library(lmtest)
lrtest(t5a_whole,t5a_b) # pathYN not worth it
lrtest(t5a_whole,t5a_c) # trtmtM not worth it
lrtest(t5a_whole,t5a_d) # trtmtA not worth it
lrtest(t5a_whole,t5a_e) # whole better WITH time

# All better than only colony
lrtest(t5a_f,t5a_b) #
lrtest(t5a_f,t5a_c) # 
lrtest(t5a_f,t5a_d) # 
lrtest(t5a_f,t5a_e) # unless time is absent. 

# compare pairs to time alone
t5a_e1<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(time+trtmt_A |colony),REML=FALSE,data=data_no9) 

t5a_e2<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(pathYN+time |colony),REML=FALSE,data=data_no9) 

t5a_e3<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(time+trtmt_M |colony),REML=FALSE,data=data_no9) 

t5a_e4<-lmer(sqrt_t5a9~trtmt_M+trtmt_A + pathYN + time +(time |colony),REML=FALSE,data=data_no9) 

AIC(t5a_e,t5a_e1,t5a_e2,t5a_e3,t5a_e4) # best is t5a_e4 with only time.

# none with 3 or 4 are better than time alone: all NS
lrtest(t5a_e4,t5a_whole) 
lrtest(t5a_e4,t5a_b) 
lrtest(t5a_e4,t5a_c) 
lrtest(t5a_e4,t5a_d) 

# none with 2 random slopes are better than time alone: all NS
lrtest(t5a_e4,t5a_e1)
lrtest(t5a_e4,t5a_e2)
lrtest(t5a_e4,t5a_e3)

#### T5A MODELS: random slope is time  ######
# Must have time as a fixed effect 

t5null<-lmer(sqrt_t5a9~1+(time|colony),REML = FALSE,data=data_no9) 

# M&A 
t5_0<-lmer(sqrt_t5a9~trtmt_A*time+(time|colony),REML=FALSE,data=data_no9) 

t5_1<-lmer(sqrt_t5a9~trtmt_M*time+(time|colony),REML=FALSE,data=data_no9)

#t5_2<-lmer(sqrt_t5a9~trtmt_M*trtmt_A+(time|colony),REML=FALSE,data=data_no9) #### time not a fixed effect

t5_3<-lmer(sqrt_t5a9~trtmt_M*trtmt_A+time+(time|colony),REML=FALSE,data=data_no9)

t5_4<-lmer(sqrt_t5a9~trtmt_M*trtmt_A*time+(time|colony),REML=FALSE,data=data_no9) 

t5_5<-lmer(sqrt_t5a9~trtmt_M+trtmt_A*time+(time|colony),REML=FALSE,data=data_no9) 

t5_6<-lmer(sqrt_t5a9~trtmt_M*time+trtmt_A+(time|colony),REML=FALSE,data=data_no9) 

t5_7<-lmer(sqrt_t5a9~trtmt_M+time+trtmt_A+(time|colony),REML=FALSE,data=data_no9) 

#t5_8<-lmer(sqrt_t5a9~trtmt_M+trtmt_A+(time|colony),REML=FALSE,data=data_no9) ## time not a fixed effect

t5_9<-lmer(sqrt_t5a9~trtmt_A+time+(time|colony),REML=FALSE,data=data_no9) 

t5_10<-lmer(sqrt_t5a9~trtmt_M+time+(time|colony),REML=FALSE,data=data_no9) 

#t5_11<-lmer(sqrt_t5a9~trtmt_A+(time|colony),REML=FALSE,data=data_no9)  ## time not a fixed effect

#t5_12<-lmer(sqrt_t5a9~trtmt_M+(time|colony),REML=FALSE,data=data_no9) ## time not a fixed effect

t5_13<-lmer(sqrt_t5a9~pathYN*time+(time|colony),REML=FALSE,data=data_no9) # can test b/c EF1 NS

t5_14<-lmer(sqrt_t5a9~pathYN+time+(time|colony),REML=FALSE,data=data_no9) 

#t5_15<-lmer(sqrt_t5a9~pathYN+(time|colony),REML=FALSE,data=data_no9) ## time not a fixed effect

#t5_16<-lmer(sqrt_t5a9~onepath+(time|colony),REML=FALSE,data=data_no9) ## time not a fixed effect

t5_17<-lmer(sqrt_t5a9~onepath+time+(time|colony),REML=FALSE,data=data_no9) 

t5_18<-lmer(sqrt_t5a9~onepath*time+(time|colony),REML=FALSE,data=data_no9) 


# AIC

AIC(t5null,t5_0,t5_1,t5_3,t5_4,t5_5,t5_6,t5_7,t5_9,t5_10,t5_13,t5_14,t5_17,t5_18)

# Model 13, then 14 are best by > 4 AIC 
library(lmtest)
lrtest(t5null,t5_14) # t5_14 better than null
lrtest(t5null,t5_13) # t5_13 better than null

lrtest(t5_14,t5_13) # t5_14 is best: intxn of 13 not justified

## Best model t5_14: pathYN+time
summary(t5_14) #
plot(t5_14) 

library(lmerTest)

t5_14_REML<-lmer(sqrt_t5a9~pathYN+time+(time|colony),REML = TRUE,data=data_no9) 
summary(t5_14_REML) # pathYN is significant, time is not

qqnorm(resid(t5_14_REML))
qqline(resid(t5_14_REML)) 

tapply(data$t5a_ab,data$pathYN,mean) # N: 0.162, Y: 0.225
((0.225-0.162)/0.162)*100 # 38.9% increase

library(MuMIn)
confint(t5_14_REML,level=0.95, method="Wald") # 


###################
####### PLOTS #####
#################


#### REACTION NORMS ###
data$treatment <- factor(data$treatment, levels=c("N","MA","M","A"))

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

#### 4 treatments 
# X axis
x<-c(1:4)

# t5a_AB reaction norm: 4 treatments
one_norm<-tapply(one_col$t5a_ab,one_col$treatment,mean)
two_norm<-tapply(two_col$t5a_ab,two_col$treatment,mean)
three_norm<-tapply(three_col$t5a_ab,three_col$treatment,mean)
four_norm<-tapply(four_col$t5a_ab,four_col$treatment,mean)
five_norm<-tapply(five_col$t5a_ab,five_col$treatment,mean)
six_norm<-tapply(six_col$t5a_ab,six_col$treatment,mean)
seven_norm<-tapply(seven_col$t5a_ab,seven_col$treatment,mean)
eight_norm<-tapply(eight_col$t5a_ab,eight_col$treatment,mean)
nine_norm<-tapply(nine_col$t5a_ab,nine_col$treatment,mean)
ten_norm<-tapply(ten_col$t5a_ab,ten_col$treatment,mean)


### Reaction norm: without colony 9 ###

# 4 treatments 
# X axis
x<-c(1:4)

colony_norm9<-data.frame(one_norm,two_norm,three_norm,four_norm,five_norm,six_norm,seven_norm,eight_norm,ten_norm)

colony_norm9
col_mat<-as.matrix(colony_norm9)
col_mat
### Plot
colnames<-c(1:4)
colors<-c("darkred","darkorange","gold","green","blue","purple","lightsalmon4","black","deeppink")
plot(x, y= rep(-10000,4), ylim= c(0,.5), # initiaties plot
     ylab= "t5a_AB", xlab= "Treatment", main= "T5A by Colony")    
for(i in c(1:9)){
  points(x, col_mat[,i], col=colors[i], type= "b") 
}



################################
################ BAR PLOTS with error bars ##############
##############################

data$treatment <- factor(data$treatment, levels=c("N","A","M","MA"))

## t5a_1 model
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


## WITHOUT COLONY 9 ####

# PathYN
data_no9$treatment <- factor(data_no9$treatment, levels=c("N","A","M","MA"))

data_summarize1_no9<-summarySE(data_no9,measurevar="t5a_ab", groupvars=c("treatment"))

library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize1_no9, aes(x=treatment, y=t5a_ab,fill="red")) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=t5a_ab-1*se, ymax=t5a_ab+1*se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("T5A expression (corrected to reference gene)") +
  xlab("Treatment")

