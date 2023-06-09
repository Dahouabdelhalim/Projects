##### Co-infection Gene expression experiment 2015 ####

# Date: 3-30-20

data<-read.csv("MMP_data_3-30-20.csv")
head(data)
length(data$sample) # 80 - 10 colonies*4 treatments* 2 timepoints

## Normality 
hist(data$mmp_ab) # not normal

qqnorm(data$mmp_ab)
qqline(data$mmp_ab)  

# Transform
qqnorm(sqrt(data$mmp_ab))
qqline(sqrt(data$mmp_ab))  

shapiro.test(sqrt(data$mmp_ab)) # Normal when square root transformed

## Data classes
class(data$sample) # factor
class(data$treatment) # factor
data$colony<-as.factor(data$colony)
class(data$colony) # factor
class(data$time) #
class(data$onepath)

class(data$mmp_ab)

## Contrasts
contrasts(data$treatment) <- contr.sum # contrasts control to each treatment 
contrasts(data$treatment)

contrasts(data$colony) <- contr.sum 
contrasts(data$time) <- contr.sum 

## Exploratory plots
plot(colony~sample, data=data) # 1 each 
plot(time~sample, data=data) 
plot(treatment~colony,data=data)
plot(treatment~time,data=data)
plot(colony~time,data=data)

# Response variable: mmp_AB
sqrt_mmpAB<-sqrt(data$mmp_ab)

# heteroscedastiscity
plot(sqrt_mmpAB~treatment, data=data) 
plot(sqrt_mmpAB~colony, data=data) # differences 
plot(sqrt_mmpAB~time, data=data) 


#### Without colony 9 ######
data_no9<-subset(data,colony!="9")
summary(data_no9$colony)
length(data_no9$colony) # 72 

hist(data_no9$mmp_ab)
shapiro.test(data_no9$mmp_ab) # not normal
qqnorm(data_no9$mmp_ab)
qqline(data_no9$mmp_ab)  

# Transform
qqnorm(sqrt(data_no9$mmp_ab))
qqline(sqrt(data_no9$mmp_ab))  

shapiro.test(sqrt(data_no9$mmp_ab)) # Normal when square root transformed

## Data classes
class(data_no9$sample) # factor
class(data_no9$treatment) # factor
data_no9$colony<-as.factor(data_no9$colony)
class(data_no9$colony) # factor
class(data_no9$time) #
class(data$trtmt_M) 
class(data$trtmt_A) 
class(data$onepath) 
class(data_no9$mmp_ab)

## Contrasts
contrasts(data_no9$treatment) <- contr.sum 
contrasts(data_no9$colony) <- contr.sum 
contrasts(data_no9$time) <- cbind(c(-1,1)) # early is lower
contrasts(data_no9$trtmt_M) <- cbind(c(-1,1))
contrasts(data_no9$trtmt_A) <- cbind(c(-1,1))

contrasts(data_no9$pathYN) <- cbind(c(-1,1))

contrasts(data_no9$onepath) <- cbind(c(1,0,-1),c(0,1,-1))
contrasts(data_no9$onepath)

## Exploratory plots
plot(colony~sample, data=data_no9) 
plot(treatment~sample, data=data_no9) 
plot(time~sample, data=data_no9) 
plot(treatment~colony,data=data_no9)
plot(treatment~time,data=data_no9)
plot(colony~time,data=data_no9)
plot(onepath~treatment,data=data_no9)

# Response variable: mmp_AB
sqrt_mmpAB9<-sqrt(data_no9$mmp_ab)
plot(sqrt_mmpAB9~sample, data=data_no9) 
plot(sqrt_mmpAB9~treatment, data=data_no9)  
plot(sqrt_mmpAB9~colony, data=data_no9) 
plot(sqrt_mmpAB9~time, data=data_no9) 

## Linear models ## 

# Fixed effects: timing, treatment
# Random effects: colony and test for random slopes

library(lme4)

###### EF1 #####

# Normality
hist(data_no9$ef1)
shapiro.test(data_no9$ef1) #  normal 

qqnorm(data_no9$ef1)
qqline(data_no9$ef1)  

## Steps of model selection 
#1) Test models for random slopes. 
#2) Use random slopes structure to determine best model with fixed effects

# Random slopes
ef_whole<-lmer(ef1~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) # fixed effect model matrix is rank deficient if include M*A in this model

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

# LRT 
library(lmtest)
lrtest(ef_whole,ef_b) # pathYN not worth including
lrtest(ef_whole,ef_c) # trtmt_M not worth it
lrtest(ef_whole,ef_d) # trtmt_A not worth it
lrtest(ef_whole,ef_e) # time not worth it 

# All NS vs. one with no random slopes: none are needed then
lrtest(ef_f,ef_whole)
lrtest(ef_f,ef_b)
lrtest(ef_f,ef_c)
lrtest(ef_f,ef_d)
lrtest(ef_f,ef_e)

# Onepath tested separately due to convergence issues if all in same model
ef_g<-lmer(ef1~onepath+ time +(1 |colony),REML=FALSE,data=data_no9) 
ef_g2<-lmer(ef1~onepath+ time +(onepath |colony),REML=FALSE,data=data_no9) 

lrtest(ef_g,ef_g2) # also not worth including

# None worth including as random slopes for EF1, reference gene 
# Use 1|colony for random intercept in EF1 models

### Determining fixed effects ###

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

# step 3 - AIC table

AIC(efnull,ef0,ef1,ef2,ef3,ef4,ef5,ef6,ef7,ef8,ef9,ef10,ef11,ef12,ef13,ef14,ef15,ef16,ef17,ef18)

# Best model 5 by AIC
library(lmtest)
lrtest(efnull,ef5) # not better than null by LRT 
lrtest(efnull,ef0) # not better than null

# Null is best model

###### MMP #########

## Models for random effects structure

mmp_whole<-lmer(sqrt_mmpAB9~trtmt_M*trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) # fixed effect model matrix is rank deficient if include M*A in this model

# no pathYN as random slope
mmp_b<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(trtmt_M+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_M
mmp_c<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_A+time |colony),REML=FALSE,data=data_no9) 

# no trtmt_A
mmp_d<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+time |colony),REML=FALSE,data=data_no9) 

# no time
mmp_e<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(pathYN+trtmt_M+trtmt_A |colony),REML=FALSE,data=data_no9) 

# none
mmp_f<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(1 |colony),REML=FALSE,data=data_no9) 

# LRT 
library(lmtest)
lrtest(mmp_whole,mmp_b) # whole not better
lrtest(mmp_whole,mmp_c) # whole not better
lrtest(mmp_whole,mmp_d) # whole not better
lrtest(mmp_whole,mmp_e) # whole not better

# None better than only colony
lrtest(mmp_f,mmp_b)
lrtest(mmp_f,mmp_c) 
lrtest(mmp_f,mmp_d) 
lrtest(mmp_f,mmp_e) 

mmp_g<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(time|colony),REML=FALSE,data=data_no9) 

mmp_h<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(pathYN |colony),REML=FALSE,data=data_no9) 

mmp_i<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(trtmt_M |colony),REML=FALSE,data=data_no9) 

mmp_j<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time +(trtmt_A |colony),REML=FALSE,data=data_no9) 

# Singles vs. just colony? Still all NS
lrtest(mmp_f,mmp_g)
lrtest(mmp_f,mmp_h) 
lrtest(mmp_f,mmp_i) 
lrtest(mmp_f,mmp_j) 

# onepath
mmp_g<-lmer(sqrt_mmpAB9~onepath + time +(1 |colony),REML=FALSE,data=data_no9) 
mmp_g2<-lmer(sqrt_mmpAB9~onepath + time +(onepath |colony),REML=FALSE,data=data_no9) 
lrtest(mmp_g,mmp_g2) # NS too 

# Colony random intercept
mmp_f2<-lm(sqrt_mmpAB9~trtmt_M+trtmt_A + pathYN + time,data=data_no9) 

lrtest(mmp_f2,mmp_f) # colony is worth including 

# Oly use colony as random intercept for MMP models

#### Fixed effects models ### 

mmp_null<-lmer(sqrt_mmpAB9~1+(1|colony),REML=FALSE,data=data_no9) 

mmp_1<-lmer(sqrt_mmpAB9~trtmt_M*trtmt_A*time+(1|colony),REML=FALSE,data=data_no9) 

mmp_2<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A*time+(1|colony),REML=FALSE,data=data_no9) 
mmp_2b<-lmer(sqrt_mmpAB9~trtmt_A*time+(1|colony),REML=FALSE,data=data_no9) 

mmp_3<-lmer(sqrt_mmpAB9~trtmt_M*time+trtmt_A+(1|colony),REML=FALSE,data=data_no9) 
mmp_3b<-lmer(sqrt_mmpAB9~trtmt_M*time+(1|colony),REML=FALSE,data=data_no9) 

mmp_4<-lmer(sqrt_mmpAB9~trtmt_M+time+trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

mmp_5<-lmer(sqrt_mmpAB9~trtmt_M+trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

mmp_6<-lmer(sqrt_mmpAB9~trtmt_A+time+(1|colony),REML=FALSE,data=data_no9) 

mmp_7<-lmer(sqrt_mmpAB9~trtmt_M+time+(1|colony),REML=FALSE,data=data_no9) 

mmp_8<-lmer(sqrt_mmpAB9~trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

mmp_9<-lmer(sqrt_mmpAB9~trtmt_M+(1|colony),REML=FALSE,data=data_no9) 

mmp_10<-lmer(sqrt_mmpAB9~pathYN*time+(1|colony),REML=FALSE,data=data_no9) 

mmp_11<-lmer(sqrt_mmpAB9~pathYN+time+(1|colony),REML=FALSE,data=data_no9) 

mmp_12<-lmer(sqrt_mmpAB9~pathYN+(1|colony),REML=FALSE,data=data_no9) 

mmp_13<-lmer(sqrt_mmpAB9~treatment*time+(1|colony),REML=FALSE,data=data_no9)

mmp_14<-lmer(sqrt_mmpAB9~trtmt_M*trtmt_A+(1|colony),REML=FALSE,data=data_no9) 

mmp_15<-lmer(sqrt_mmpAB9~trtmt_M*trtmt_A+time+(1|colony),REML=FALSE,data=data_no9) 

mmp_16<-lmer(sqrt_mmpAB9~onepath+(1|colony),REML=FALSE,data=data_no9) 

mmp_17<-lmer(sqrt_mmpAB9~onepath+time+(1|colony),REML=FALSE,data=data_no9) 

mmp_18<-lmer(sqrt_mmpAB9~onepath*time+(1|colony),REML=FALSE,data=data_no9) 


# AIC
AIC(mmp_null,mmp_1,mmp_2,mmp_2b,mmp_3,mmp_3b,mmp_4,mmp_5,mmp_6,mmp_7,mmp_8,mmp_9,mmp_10,mmp_11,mmp_12,mmp_13,mmp_14,mmp_15,mmp_16,mmp_17,mmp_18)

# Null is best model by AIC and LRT
library(lmtest)
lrtest(mmp_null,mmp_9)
lrtest(mmp_null,mmp_3b)
lrtest(mmp_null,mmp_10)


####### PLOTS #####

# REACTION NORMS 
data$treatment <- as.character(data$treatment)
data$treatment <- factor(data$treatment, levels=c("N", "A","M","MA")) # re-orders

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

# mmp_ab reaction norm
one_norm<-tapply(one_col$mmp_ab,one_col$treatment,mean)
two_norm<-tapply(two_col$mmp_ab,two_col$treatment,mean)
three_norm<-tapply(three_col$mmp_ab,three_col$treatment,mean)
four_norm<-tapply(four_col$mmp_ab,four_col$treatment,mean)
five_norm<-tapply(five_col$mmp_ab,five_col$treatment,mean)
six_norm<-tapply(six_col$mmp_ab,six_col$treatment,mean)
seven_norm<-tapply(seven_col$mmp_ab,seven_col$treatment,mean)
eight_norm<-tapply(eight_col$mmp_ab,eight_col$treatment,mean)
nine_norm<-tapply(nine_col$mmp_ab,nine_col$treatment,mean)
ten_norm<-tapply(ten_col$mmp_ab,ten_col$treatment,mean)


##### All 10 colonies

colony_norm<-data.frame(one_norm,two_norm,three_norm,four_norm,five_norm,six_norm,seven_norm,eight_norm,nine_norm,ten_norm)

colony_norm
col_mat<-as.matrix(colony_norm)
col_mat
### Plot
colnames<-c(1:4)
colors<-c("darkred","darkorange","gold","green","blue","purple","lightsalmon4","black","aquamarine3","deeppink")
plot(x, y= rep(-10000,4), ylim= c(0,0.08), # initiaties plot
     ylab= "mmp_ab", xlab= "Treatment", main= "MMP by Colony")    
for(i in c(1:10)){
  points(x, col_mat[,i], col=colors[i], type= "b") # redo colors b/c there are two black lines
}

# no 9
colony_norm_no9<-data.frame(one_norm,two_norm,three_norm,four_norm,five_norm,six_norm,seven_norm,eight_norm,ten_norm)

colony_norm_no9
col_mat<-as.matrix(colony_norm_no9)
col_mat
### Plot
colnames<-c(1:4)
colors<-c("darkred","darkorange","gold","green","blue","purple","lightsalmon4","black","deeppink")
plot(x, y= rep(-10000,4), ylim= c(0,0.08), # initiates plot
     ylab= "mmp_ab", xlab= "Treatment", main= "No9 MMP by Colony")    
for(i in c(1:10)){
  points(x, col_mat[,i], col=colors[i], type= "b") # redo colors b/c there are two black lines
}

# Variability among colonies
