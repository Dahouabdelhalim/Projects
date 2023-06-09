# Energy balance and feeding competition in female bonobos. 
# Model response is log of urinary C-peptide ng/mg crea level per female.
# Updated variables March 2018: "proportion fruit scans" replaced "sine/cosine" for seasonality index,
# and "proportion feed scans" replaced "proportion time in focal patches". Added control 
# for behavioural shift, because feeding may differ between AM and PM shifts.




# set up the working directory:

setwd("Y:/primint/Niina/2018_ucp_model_seasonality_update")
getwd()


# check the list to make sure its empty
ls()
rm(list=ls())

#######################

save.image("Y:/primint/Niina/2018_ucp_model_seasonality_update/main_ucp_model_workspace_05032018.Rdata")
load("Y:/primint/Niina/2018_ucp_model_seasonality_update/main_ucp_model_workspace_05032018.Rdata")


##############################


# read in C-peptide data from Primint as "xdata":
xdata=read.table(file="main_ucp_model_data_r_2018_seasonality_update_with_feed_scans.txt", header=T, sep="\\t")

str(xdata)


# --------------INSPECTION------------------------------------------------------------------------


# first inspect data as tables:
table(xdata$bono_id)
# Djulie   Gwen   Iris   Luna Martha   Nina   Olga  Paula  Polly    Rio   Susi    Uma  Wilma    Zoe 
#   5     10     28      7     15     19     27     18      9     10     18     14     16     22 

# 14 females with more or less balanced representation.


# check year/month/day entries:
table(xdata$year_coll)
# 2012 2013 2014 
# 64   88   66

table(xdata$month_coll)
# 1  2  3  4  5  6  7  8  9 10 11 12 
# 23 19 37 26 15 10 14 11 11  7 13 32


table(xdata$day_coll)
# all entries look fine


# check entries for lactation:
table(xdata$lactation)
# no yes 
# 28 190 


# check entries for rank:
table(xdata$rank)
# 0 0.143 0.286 0.429 0.571 0.714 0.857     1 
# 5     9    19    25    16    10    66    68


min(xdata$tot_party)
# 3
max(xdata$tot_party)
# 20.69



# checking assigments of bonobo IDs to ranks
xx=table(xdata$bono_id, xdata$rank)

range(apply(X=xx>0, MARGIN=1, FUN=sum))
# 1 1

# only one rank assigned to each female, looks good.


# TIME as R-TIME:

# time currently coded as hour and minute of sample collection
# make a new column "time elapsed since midnight"
# incase we need to control for urine sample collection time:

xdata$urine_coll_time=xdata$hour_coll*3600 + xdata$min_coll*60
hist(xdata$urine_coll_time)



# DISTRIBUTIONS:

# check the distribution of the response, urinary C-peptide value:
hist(xdata$cpep)

# distribution is quite left skewed, better to log-transform:
min(xdata$cpep)
# 0.387

hist(log(xdata$cpep))

# log-transformed looks okay, store them in a new column:
xdata$log_cpep=log(xdata$cpep)

# ----------------------

# check distribution of the main test predictor variables:
hist(xdata$tot_party)
# total party size looks well balanced. 

# -------------------

# check number of patches visited/hour:
hist(xdata$patches_per_h)

# number of patches visited/hour looks left skewed, better to log-transform: 
# check min value:
min(xdata$patches_per_h)
# 0.16

hist(log(xdata$patches_per_h))

# log-transformed look more balanced, store them in new column:
xdata$log_patches_per_h=log(xdata$patches_per_h)

# ------------------------

# check proportion of THV scans distribution:
hist(xdata$prop_thv_scans)

hist(sqrt(xdata$prop_thv_scans))

xdata$sq_root_thv=sqrt(xdata$prop_thv_scans)

hist(xdata$sq_root_thv)

# -----------------------


# check meters travelled per hour:
hist(xdata$meters_per_hour)

# looks pretty okay, leave as is.


# --------------------------

# check proportion of fruit scans:
hist(xdata$aver_month_fruit_scan)

# it looks okay

# ---------------------

# then the new proportion of feed scans over all 30min scans:
hist(xdata$prop_feed_scans)

# leave as is, looks okay

# -------------------------------




# DATE in to R-DATE:

# make day ID as random effect, currently d/m/y 
# as separate columns, create new variable:

xdate=paste(xdata$day_coll, xdata$month_coll, xdata$year_coll, sep="/")

xdate=as.Date(xdate, format="%d/%m/%Y")
xdate

xdata$new_date=xdate



str(xdata)

    
aggregate(cpep~bono_id, data=xdata, mean)
aggregate(aver_month_fruit_scan~, data=xdata, mean)
aggregate(patches_per_h~, data=xdata, mean)

mean(xdata$cpep~bono_id)


####################################

mean(xdata$aver_month_fruit_scan)
#0.6775076

min(xdata$aver_month_fruit_scan)
#0.3828289

max(xdata$aver_month_fruit_scan)
# 0.8841694

######################################

mean(xdata$patches_per_h)
# 0.6958716

min(xdata$patches_per_h)
# 0.16

max(xdata$patches_per_h)
# 2.89

#####################################

mean(xdata$prop_thv_scans) 
# 0.134403

min(xdata$prop_thv_scans) 
# 0

max(xdata$prop_thv_scans)
# 0.43
  
###################################


mean(xdata$prop_feed_scans) 
# 0.391789

min(xdata$prop_feed_scans) 
# 0

max(xdata$prop_feed_scans)
# 0.88

##########################



  
# RANDOM SLOPES

source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\diagnostic_fcns.R")

# check which ones log-transformed:
# log "UCP" (response), and log "number patches per hour", and square root THV.
# "prop feed scans" is okay, 'tot party" is okay and "meters per hour" is okay, "aver_month_fruit_scan" is okay.


rand_slope_data=fe.re.tab(fe.model="log_cpep~tot_party*rank + prop_feed_scans*rank + 
                          sq_root_thv + meters_per_hour + aver_month_fruit_scan + lactation + beh_shift", 
                          re="(1|bono_id)+(1|new_date)", data=xdata)


rand_slope_data$summary

#$`tot_party_within_bono_id (covariate)`     YES
#5  7  9 10 11 13 15 16 18 21 24 
#1  2  1  1  1  1  1  1  2  1  2            

#$`tot_party_within_new_date (covariate)`    NO
#1   2 
#110   5 

#$`rank_within_bono_id (covariate)          NO
#1 
#14 

#$`rank_within_new_date (covariate)       YES    
#1  2  3  4  5 
#57 40 12  5  1 

#$`prop_feed_scans_within_bono_id (covariate)`      NO
#4  6  9 10 12 13 16 18 20 24 
#1  1  2  1  2  2  1  1  1  2 

#$`prop_feed_scans_within_new_date (covariate)`     NO
#1   2 
#109   6 

#$`sq_root_thv_within_bono_id (covariate)`        YES
#5  6  8  9 10 12 13 16 
#2  3  1  2  2  1  2  1 

#$`sq_root_thv_within_new_date (covariate)`     NO
#1   2 
#110   5 

#$`meters_per_hour_within_bono_id (covariate)`    YES
#5  7  9 10 13 16 18 22 27 
#1  1  1  2  2  1  3  1  2 

#$`meters_per_hour_within_new_date (covariate)`     NO
#1   2 
#109   6 

#$`aver_month_fruit_scan_within_bono_id (covariate)`  YES
#2  5  7  8  9 10 11 12 14 
#1  2  1  3  2  1  1  1  2 

#$`aver_month_fruit_scan_within_new_date (covariate)`   NO
#1 
#115 

#$`lactation_within_bono_id (factor)`                NO
#1 
#14 

#$`lactation_within_new_date (factor)`          NO
#0  1  2 
#58 54  3 

#$`beh_shift_within_bono_id (factor)`           NO
#1  2 
#1 13 

#$`beh_shift_within_new_date (factor)`          NO
#0  1 
#54 61 

#$tot_party_rank_within_bono_id            NO
#tot_party rank n.bono_id
#1          5    1         1
#2          7    1         2
#3          9    1         1
#4         10    1         1
#5         11    1         1
#6         13    1         1
#7         15    1         1
#8         16    1         1
#9         18    1         2
#10        21    1         1
#11        24    1         2

#$tot_party_rank_within_new_date       ?
#tot_party rank n.new_date
#1         1    1         56
#2         2    1          1
#3         1    2         36
#4         2    2          4
#5         1    3         12
#6         1    4          5
#7         1    5          1

#$rank_prop_feed_scans_within_bono_id         ?
#rank prop_feed_scans n.bono_id
#1     1               4         1
#2     1               6         1
#3     1               9         2
#4     1              10         1
#5     1              12         2
#6     1              13         2
#7     1              16         1
#8     1              18         1
#9     1              20         1
#10    1              24         2

#$rank_prop_feed_scans_within_new_date          NO
#rank prop_feed_scans n.new_date
#1    1               1         56
#2    2               1         35
#3    3               1         12
#4    4               1          5
#5    5               1          1
#6    1               2          1
#7    2               2          5


# include following random slopes:
# within "bono id" :  "tot party", "prop_feed_scans, "log patches per h", "meter_per_h", "sq_root_thv", "aver_month_fruit_scan"
# within "new date" : "rank"




# Z_TRANSFORMATIONS: 

# z-transform predictor covariates:

xdata$z.tot_party=as.vector(scale(xdata$tot_party))

xdata$z.prop_feed_scans=as.vector(scale(xdata$prop_feed_scans))

xdata$z.meters_per_hour=as.vector(scale(xdata$meters_per_hour))

xdata$z.rank=as.vector(scale(xdata$rank))

xdata$z.sqrt_thv=as.vector(scale(xdata$sq_root_thv))

xdata$z.log_patches_per_hour=as.vector(scale(xdata$log_patches_per_h))

xdata$z.aver_month_fruit_scan=as.vector(scale(xdata$aver_month_fruit_scan))


# CORRELATIONS 

# prepare data-frame tdata:
tdata=data.frame(xdata)

vars=c("z.sqrt_thv", "z.tot_party", "z.meters_per_hour", "z.prop_feed_scans", "z.log_patches_per_hour", "z.aver_month_fruit_scan")

cor(tdata[, vars])


#                        z.sqrt_thv z.tot_party z.meters_per_hour z.prop_feed_scans z.log_patches_per_hour z.aver_month_fruit_scan
# z.sqrt_thv               1.00000000  0.06867110        0.41873778        0.18357868            -0.02262039             -0.16650330
# z.tot_party              0.06867110  1.00000000        0.01040873        0.24724569            -0.16478837             -0.01581122
# z.meters_per_hour        0.41873778  0.01040873        1.00000000       -0.13517193             0.19398051             -0.05555038
# z.prop_feed_scans        0.18357868  0.24724569       -0.13517193        1.00000000            -0.11084909              0.05498503
# z.log_patches_per_hour  -0.02262039 -0.16478837        0.19398051       -0.11084909             1.00000000             -0.18948033
# z.aver_month_fruit_scan -0.16650330 -0.01581122       -0.05555038        0.05498503            -0.18948033              1.00000000

# highest cor 0.42 between meters per hour and propoportion of THV scans over all scans


# RUN FULL MODEL:

library(lme4)


# prepare control object before running the full model:
control_model=lmerControl(optCtrl = list(maxfun=100000))

date.fac=as.factor(as.character(xdata$new_date))

res=lmer(log_cpep~z.tot_party * z.rank + z.prop_feed_scans * z.rank + z.log_patches_per_hour * z.rank + z.sqrt_thv + z.meters_per_hour + 
           z.aver_month_fruit_scan + lactation + beh_shift +
           (1 + z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.sqrt_thv + z.meters_per_hour + z.aver_month_fruit_scan||bono_id) + 
           (1 + z.rank||date.fac), 
         data=xdata, REML=F, control=control_model)



# MODEL ASSUMPTIONS/VIF:

source("diagnostic_fcns.r")

# check plots:
ranef.diagn.plot(res)  # ok
diagnostics.plot(res)  # ok 


# check VIFs (exclude random effects):
library(car)

vif_model=lm(log_cpep~z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.rank + z.sqrt_thv + 
               z.meters_per_hour + z.aver_month_fruit_scan + lactation + beh_shift , 
             data=xdata) 

vif(vif_model)

#if bigger than 10 then reason to worry.

# z.tot_party       z.prop_feed_scans  z.log_patches_per_hour                  z.rank              z.sqrt_thv       z.meters_per_hour 
# 1.161813                1.328669                1.139359                   1.174189                1.719049                1.482247 

# z.aver_month_fruit_scan               lactation               beh_shift 
# 1.180313                               1.216369                1.712017 

# Highest is 1.74, okay.



#ADD date.fac to xdata:

xdata$date.fac=date.fac



# MODEL STABILITY:

# source Rogers function:
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\glmm_stability.r")
ls()

head(glmm.model.stab)

model_stability=glmm.model.stab(model.res=res, contr=control_model, n.cores=6, para=T, data=xdata)
table(model_stability$detailed$warnings)

round(model_stability$summary[, -1], 3)

#                                     orig    min    max
# (Intercept)                         0.549  0.361  0.847
# z.tot_party                         0.060  0.016  0.086
# z.rank                             -0.048 -0.089  0.020
# z.prop_feed_scans                  -0.023 -0.056  0.015
# z.log_patches_per_hour             -0.063 -0.080 -0.049
# z.sqrt_thv                          0.054  0.017  0.085
# z.meters_per_hour                   0.050  0.006  0.122
# lactationyes                        0.062 -0.259  0.217
# z.aver_month_fruit_scan             0.176  0.149  0.200
# beh_shiftpm                         0.120  0.048  0.173
# z.tot_party:z.rank                 -0.024 -0.087  0.042
# z.rank:z.prop_feed_scans            0.039 -0.015  0.058
# z.rank:z.log_patches_per_hour       0.016  0.003  0.028
# date.fac@z.rank@NA                  0.000  0.000  0.000
# date.fac@(Intercept)@NA             0.228  0.179  0.277
# bono_id@z.aver_month_fruit_scan@NA  0.000  0.000  0.028
# bono_id@z.meters_per_hour@NA        0.049  0.000  0.078
# bono_id@z.sqrt_thv@NA               0.000  0.000  0.000
# bono_id@z.log_patches_per_hour@NA   0.000  0.000  0.000
# bono_id@z.prop_feed_scans@NA        0.073  0.000  0.109
# bono_id@z.tot_party@NA              0.144  0.040  0.186
# bono_id@(Intercept)@NA              0.190  0.125  0.214
# Residual                            0.452  0.419  0.476

m.stab.plot(model_stability$summary[, -1])
# looks ok, expect for lactation was a bit off



# NULL MODEL:
null_model=lmer(log_cpep~lactation + beh_shift +
                  (1 + z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.sqrt_thv + z.meters_per_hour + z.aver_month_fruit_scan||bono_id) + 
                  (1 + z.rank||date.fac), 
                data=xdata, REML=F, control=control_model)




# FULL - NULL COMPARISON:
anova(null_model, res, test="Chisq")

#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# null_model 13 392.94 436.94 -183.47   366.94                           
# res        23 392.91 470.75 -173.45   346.91 20.029     10    0.02898 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# the full null is significant.




#P-VALUES FULL MODEL
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\drop1_para.r")
ls()
head(drop1p)

drop1p(model.res=res, data=xdata, contr=control_model)

#                               logLik      AIC      Chisq Chi.Df   Pr..Chisq.   n n.opt.warnings n.fun.warnings
# none                          -173.4542 392.9084         NA     NA           NA 218             NA             NA
# z.sqrt_thv                    -173.9769 391.9538  1.0453391      1 0.3065829861 218              0              0
# z.meters_per_hour             -173.9201 391.8402  0.9317929      1 0.3343974696 218              0              0
# lactation                     -173.5101 391.0201  0.1117147      1 0.7382003326 218              0              0
# z.aver_month_fruit_scan       -179.1521 402.3042 11.3957410      1 0.0007361268 218              0              0
# beh_shift                     -174.0417 392.0835  1.1750641      1 0.2783631443 218              0              0
# z.tot_party:z.rank            -173.5407 391.0814  0.1729626      1 0.6774921815 218              0              0  NOT SIGNIF
# z.rank:z.prop_feed_scans      -173.8582 391.7163  0.8078975      1 0.3687425867 218              0              0  NOT SIGNIF
# z.rank:z.log_patches_per_hour -173.5470 391.0939  0.1855034      1 0.6666857859 218              0              0  NOT SIGNIF


summary(res)

# Number of obs: 218, groups:  date.fac, 115; bono_id, 14
# 
# Fixed effects:
#                                Estimate Std. Error t value
# (Intercept)                    0.54858    0.18897   2.903
# z.tot_party                    0.05957    0.05995   0.994
# z.rank                        -0.04814    0.06584  -0.731
# z.prop_feed_scans             -0.02320    0.04969  -0.467
# z.log_patches_per_hour        -0.06347    0.04208  -1.508
# z.sqrt_thv                     0.05423    0.05145   1.054
# z.meters_per_hour              0.05049    0.05064   0.997
# lactationyes                   0.06178    0.18259   0.338
# z.aver_month_fruit_scan        0.17639    0.04266   4.135
# beh_shiftpm                    0.11956    0.10825   1.104
# z.tot_party:z.rank            -0.02370    0.05659  -0.419
# z.rank:z.prop_feed_scans       0.03928    0.04329   0.907
# z.rank:z.log_patches_per_hour  0.01576    0.03597   0.438




#REDUCED MODEL:
res_nointer=lmer(log_cpep~z.rank + z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.sqrt_thv + z.meters_per_hour + 
                  z.aver_month_fruit_scan + lactation + beh_shift +
                  (1 + z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.sqrt_thv + z.meters_per_hour + z.aver_month_fruit_scan||bono_id) + 
                  (1 + z.rank||date.fac), 
                         data=xdata, REML=F, control=control_model)


#P-VALUES REDUCED MODEL:
drop1p(model.res=res_nointer, data=xdata, contr=control_model)


#                           logLik      AIC       Chisq Chi.Df   Pr..Chisq.   n n.opt.warnings n.fun.warnings
# none                    -174.0233 388.0465          NA     NA           NA 218             NA             NA
# z.rank                  -174.3333 386.6665  0.62000049      1 0.4310471050 218              0              0
# z.tot_party             -174.5360 387.0719  1.02543882      1 0.3112323708 218              0              0
# z.prop_feed_scans       -174.1118 386.2237  0.17719653      1 0.6737936477 218              0              0
# z.log_patches_per_hour  -175.0611 388.1223  2.07574996      1 0.1496563638 218              0              0
# z.sqrt_thv              -174.6426 387.2851  1.23859917      1 0.2657411986 218              0              0
# z.meters_per_hour       -174.4681 386.9362  0.88969202      1 0.3455607824 218              0              0
# z.aver_month_fruit_scan -179.6664 397.3327 11.28620865      1 0.0007808498 218              0              0    SIGNIF
# lactation               -174.0713 386.1427  0.09616349      1 0.7564830969 218              0              0
# beh_shift               -174.4908 386.9816  0.93513432      1 0.3335323206 218              0              0



#DIRECTION OF EFFECT:
summary(res_nointer)

# Fixed effects:
#                         Estimate Std. Error t value
# (Intercept)              0.56092    0.18623   3.012
# z.rank                  -0.05167    0.06524  -0.792
# z.tot_party              0.06253    0.06119   1.022
# z.prop_feed_scans       -0.02216    0.04943  -0.448
# z.log_patches_per_hour  -0.06131    0.04208  -1.457
# z.sqrt_thv               0.05872    0.05127   1.145
# z.meters_per_hour        0.04879    0.04988   0.978
# z.aver_month_fruit_scan  0.17402    0.04270   4.075   +tive effect
# lactationyes             0.05657    0.18043   0.314
# beh_shiftpm              0.10560    0.10715   0.986






#CONF. INTERVALS:
cbind(coefficients(res_nointer), confint(object=res_nointer))

#                         2.5 %               97.5 %    
#   .sig01                List,10 0           0.1499861 
# .sig02                  List,10 0.06399514  0.338225  
# .sig03                  List,10 0           0.1354985 
# .sig04                  List,10 0           0.1598104 
# .sig05                  List,10 0           0.1364289 
# .sig06                  List,10 0           0.1113634 
# .sig07                  List,10 0           0.2099749 
# .sig08                  List,10 0           0.3034142 
# .sig09                  List,10 0.08746256  0.3405814 
# .sigma                  List,10 0.3902123   0.5287163 
# (Intercept)             List,10 0.1691378   0.9373134 
# z.rank                  List,10 -0.1933772  0.08307101
# z.tot_party             List,10 -0.06223659 0.2005284 
# z.prop_feed_scans       List,10 -0.1202396  0.09173599
# z.log_patches_per_hour  List,10 -0.1457974  0.0222832 
# z.sqrt_thv              List,10 -0.04488449 0.1638396 
# z.meters_per_hour       List,10 -0.05314449 0.1603789 
# z.aver_month_fruit_scan List,10 0.08539206  0.2590736 
# lactationyes            List,10 -0.3080111  0.4520507 
# beh_shiftpm             List,10 -0.1103959  0.3182453



save.image("Y:/primint/Niina/2018_ucp_model_seasonality_update/main_ucp_model_workspace_05032018.Rdata")
load("Y:/primint/Niina/2018_ucp_model_seasonality_update/main_ucp_model_workspace_05032018.Rdata")

############################################################################################################
############################################################################################################


#UCP vs aver month PROPORTION FRUIT SCANS
tiff("D:/Users/niina_nurmi/Documents/research_2017_pan_manuscripts/manuscript_i_energy_balance/ms_i_behavioural_ecology/figure_5_tiff_ucp_fruit.tiff", width = 10, height = 7.4, units = 'in', res = 350)
plot.new()

xdata$lactation.code=as.numeric(xdata$lactation==levels(xdata$lactation)[2])
xdata$lactation.code=xdata$lactation.code-mean(xdata$lactation.code)
plot.res_nointer = lmer(log_cpep~z.rank + z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.meters_per_hour + z.sqrt_thv + 
                          z.aver_month_fruit_scan + lactation.code + beh_shift +
          (1 + z.tot_party + z.prop_feed_scans + z.log_patches_per_hour + z.meters_per_hour + z.sqrt_thv + z.aver_month_fruit_scan||bono_id) +
          (1 + z.rank||date.fac), 
                 data=xdata, REML=F, control=control_model)


par(mar=c(5.9, 6.3, 1.0, 1.5), mgp=c(4.6, 0.9, 0), tcl=-0.3)
#par(mar=c(6.4, 6.4, 0.5, 0.8), mgp=c(3.6, 0.8, 0.0))
plot(x=xdata$z.aver_month_fruit_scan, y=xdata$log_cpep, ylim=c(-1.0, 2.5), pch=19, cex=1.7, 
     col=adjustcolor("Black", alpha=0.5), xlab="Proportion of fruit feeding scans (monthly average)", 
     ylab="Urinary C-peptide levels (ng/mg Crea)", 
     xaxt="n", yaxt="n", las=1, cex.lab=2.2, cex.axis=1.9)
range(xdata$aver_month_fruit_scan)

x.labs=pretty(xdata$aver_month_fruit_scan)
axis(side=1, at=((x.labs)-mean(xdata$aver_month_fruit_scan))/sd(xdata$aver_month_fruit_scan), labels=x.labs, 
     cex.lab=2.3 , cex.axis=1.9)

y.labs=0.5*2^(0:4)
axis(side=2, at=log(y.labs), labels=y.labs, las=1, cex.lab=2.2 , cex.axis=1.9)

plot.coeffs<-fixef(plot.res_nointer)

plot.xvals<-seq(min(xdata$z.aver_month_fruit_scan), max(xdata$z.aver_month_fruit_scan), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.aver_month_fruit_scan"]*plot.xvals

lines(plot.xvals, plot.yvals, lwd=3, lty=2)

dev.off()








###########################################################################################################

#---OLDER GRAPH-SCRIPT HERE--------------------------------------------------------------------------
#GRAPH RESULTS:

#UCP vs PATCH OCCUPANCY
xdata$lactation.code=as.numeric(xdata$lactation==levels(xdata$lactation)[2])
xdata$lactation.code=xdata$lactation.code-mean(xdata$lactation.code)
plot.res_nointer = lmer(log_cpep ~ z.rank + z.tot_party + z.dur_feed_over_obs + z.log_patches_per_hour + z.meters_per_hour + z.sq_root_thv + 
                          lactation.code + z.prop_fruit_scans + 
                          (1 + z.tot_party + z.dur_feed_over_obs + z.log_patches_per_hour + z.meters_per_hour + z.sq_root_thv + z.prop_fruit_scans||bono_id) +
                          (1 + z.rank||date.fac), 
                        data=xdata, REML=F, control=control_model)

par(mar=c(6.4, 6.4, 0.5, 0.5), mgp=c(4.0, 0.9, 0))
plot(x=xdata$z.dur_feed_over_obs, y=xdata$log_cpep, ylim=c(-1.0, 2.5), pch=19, cex=1.4, 
     col=adjustcolor("Black", alpha=0.5), xlab="Proportion of time in focal food patches", ylab="Urinary C-peptide levels (ng/mg crea)", 
     xaxt="n", yaxt="n", las=1, cex.lab=1.6, cex.axis=1.4)
range(xdata$dur_feed_over_obs)

axis(side=1, at=(pretty(xdata$dur_feed_over_obs)-mean(xdata$dur_feed_over_obs))/sd(xdata$dur_feed_over_obs), labels=pretty(xdata$dur_feed_over_obs), 
     cex.lab=1.5 , cex.axis=1.4)

y.labs=0.5*2^(0:4)
axis(side=2, at=log(y.labs), labels=y.labs, las=1, cex.lab=1.6 , cex.axis=1.4)

plot.coeffs<-fixef(plot.res_nointer)

plot.xvals<-seq(min(xdata$z.dur_feed_over_obs), max(xdata$z.dur_feed_over_obs), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.dur_feed_over_obs"]*plot.xvals

lines(plot.xvals, plot.yvals, lwd=3, lty=2)




#--------skip confint-----------------------------------------
dur.data=data.frame(z.dur_feed_over_obs=seq(from=min(xdata$z.dur_feed_over_obs), to=max(xdata$z.dur_feed_over_obs),
                                            length.out=100), z.log_patches_per_hour=mean(xdata$z.log_patches_per_hour),
                    z.exp_prop_fruit_scans=mean(xdata$z.exp_prop_fruit_scans), z.sq_root_thv=mean(xdata$z.sq_root_thv),
                    z.rank=mean(xdata$z.rank), z.meters_per_hour=mean(xdata$z.meters_per_hour), z.tot_party=mean(xdata$z.tot_party))

ci.plot=predict.lm(object=res_nointer, newdata = dur.data, interval="confidence")
#-------------------------------------------------------------




#ucp vs thv
par(mar=c(6.4, 6.4, 0.5, 0.5), mgp=c(4.5, 0.9, 0))
plot(x=xdata$z.sq_root_thv, y=xdata$log_cpep, ylim=c(-1.0, 2.5), pch=19, cex=1.4, 
     col=adjustcolor("Black", alpha=0.4), xlab="Proportion of THV scans", ylab="Urinary C-peptide levels [ng/mg crea]", 
     xaxt="n", yaxt="n", las=1, cex.lab=1.5, cex.axis=1.4)
x.labs=pretty(xdata$prop_thv_scans)
axis(side=1, at=(sqrt(x.labs)-mean(xdata$sq_root_thv))/sd(xdata$sq_root_thv), labels=x.labs, 
     cex.lab=1.5 , cex.axis=1.4)
y.labs=0.5*2^(0:4)
axis(side=2, at=log(y.labs), labels=y.labs, las=1,
     cex.lab=1.5 , cex.axis=1.4)

plot.coeffs<-fixef(plot.res_nointer)

plot.xvals<-seq(min(xdata$z.sq_root_thv), max(xdata$z.sq_root_thv), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.sq_root_thv"]*plot.xvals
#plot.yvals<-exp(plot.yvals)/(1+exp(plot.yvals))

lines(plot.xvals, plot.yvals, lwd=3, lty=2)
#-----------------------------------------------------------------------------------------------------
