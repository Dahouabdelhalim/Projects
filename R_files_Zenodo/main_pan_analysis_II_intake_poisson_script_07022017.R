#-----MODEL 2. Feeding efficiency and feeding competition in bonobo females 14.11.2016-----
#The model predicts the feeding efficiency of females via their intake rate
#values as the response (Poisson). The main test predictors are rank and its interactions
#with number of individuals in patch, patch size and number of fruits in patch.
#To account for depletion effect, we have time elapsed counting from end of bout as an
#additional test predictor. As of 22.09, we rerun the model with a squared time elap term added.
#We have time of day (elapsed since midnight) as a fixed control predictor and bonobo ID,
#tree species as a random control predictor.

#4.11. rank squared and time until end of bout standardised as proportion into bout 


############################################################################
save.image("Y:/primint/Niina/Files_for_R_analysis_II_intake/MAIN_pan_analysis_II_intake_poisson_workspace_18112016.Rdata")
load("Y:/primint/Niina/Files_for_R_analysis_II_intake/MAIN_pan_analysis_II_intake_poisson_workspace_18112016.Rdata")

############################################################################################################################


#set up the working directory for MODEL 2:

setwd("Y:/primint/Niina/Files_for_R_analysis_II_intake")
getwd()


#check the list to make sure its empty
ls()


#the list needs to be cleared
rm(list=ls())

#then read in the intake rate data from Primint as "xdata":
xdata=read.table(file="MASTER_analyses_II_R_file_intake_f.txt", header=T, sep="\\t")

str(xdata)


#-----------INSPECT DATA---------------------------

#first inspect data for typos etc.
table(xdata$bonobo_id)


#14 females; more or less balanced (Luna underrepresented, Iris & Susi overrepresented)
#next check the year/month/day entries:

table(xdata$year)
# 2012 2013 2014 
# 1028 2994  528 

table(xdata$month)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31 
# 135  98  85 175  90 172 103 170 279 105 111  91 103 172 249 289  99 158 132 143 224 175  80 195  76 136 178 153  89 194  91 


table(xdata$day)


#all dates look fine; next we check ranks:
table(xdata$rank)
# 0  0.142857 0.2857143 0.4285714 0.5714286 0.7142857 0.8571429         1 
# 223       236       343      1024       267       120       951      1386 



#check food spp and dbh as tables:
table(xdata$foc_patch_spp)

range(xdata$dbh_cm)

#check the assignments per rank:
xx=table(xdata$bonobo_id, xdata$rank)

range(apply(X=xx>0, MARGIN=1, FUN=sum))

table(xdata$fruit_crop_categ)

table(xdata$ind_in_patch)


#--------------DISTRIBUTION/TRANSFORM--------------


#after inspecting for typos and inconsistencies, check the distribution of 
#our response, intake - should follow a Poisson distribution:

hist(xdata$intake_rate)


range(xdata$intake_rate)
#0 38

mean(xdata$intake_rate)
#7.821319




#looks well balanced, so no need to transform.
#next check the predictor distribution:
hist(xdata$ind_in_patch)

hist(xdata$dbh_cm)


hist(xdata$fruit_crop_categ)

aggregate(fruit_crop_categ~foc_patch_id, data=xdata, mean)

median(xdata$fruit_crop_categ)
#4

mean(xdata$fruit_crop_categ)
#3.787473

range(xdata$fruit_crop_categ)
#2 to 5


#log transform predictors:

xdata$log_dbh=log(xdata$dbh_cm)
xdata$log_ind_in=log(xdata$ind_in_patch)



#---------DATES, TIMES------------------------------------------



#next we need to transform the time, as entered in separate columns
#will make it a new column "time elapsed since midnight":

xdata$new_time=xdata$hour*3600 + xdata$min*60 + xdata$sec

hist(xdata$new_time)

#new time is now time elapsed since midnight in seconds.
#we also have the "end of bout time" in three colums,
#which we need to compare to our focal minute;

xdata$new_end_time=xdata$end_hour*3600 + xdata$end_min*60 + xdata$end_sec
str(xdata)



#----------Z-TRANSFORM------------------------

#next we z-transform the predictor covariates:
xdata$z.log_ind_in=as.vector(scale(xdata$log_ind_in))       
xdata$z.rank=as.vector(scale(xdata$rank))
xdata$z.log_dbh=as.vector(scale(xdata$log_dbh))
xdata$z.fruit_crop_categ=as.vector(scale(xdata$fruit_crop_categ))
xdata$z.time_elap_end=as.vector(scale(xdata$time_elap_end)) 
xdata$z.new_time=as.vector(scale(xdata$new_time)) 


#---------------------------------------------


#here create from z transf time ela since midnight a new squared term to account for satiation etc
#but first have to z transform before squaring the time!
xdata$z.sq_time=xdata$z.new_time*xdata$z.new_time

#create bout id by combining date to consequtive running nth focal tree of day:

xdata$bout_id=paste(xdata$year, xdata$month, xdata$day, xdata$foc_patch_no, sep="/")

table(xdata$bout_id)
sum(xdata$bout_id)
unique(xdata$bout)


#create proportion time into bout variable

#get start times per bout id
bout.start.times<-aggregate(xdata$new_time, list(xdata$bout_id), min)
#match into xdata based on bout id
xdata$bout.start.time<-bout.start.times$x[match(xdata$bout_id, bout.start.times$Group.1)]

sum(xdata$bout.start.time>xdata$new_end_time)
sum(xdata$new_time<xdata$bout.start.time)

xdata$prop.into.bout<-(xdata$new_time - xdata$bout.start.time)/(xdata$new_end_time - xdata$bout.start.time)

plot(xdata$prop.into.bout, xdata$new_end_time - xdata$bout.start.time, pch=19, col=adjustcolor("black", alpha=0.2))


xdata$z.prop.into.bout<-as.vector(scale(xdata$prop.into.bout))
xdata$z.prop.into.bout.squared<-xdata$z.prop.into.bout*xdata$z.prop.into.bout

range(xdata$prop.into.bout)
mean(xdata$prop.into.bout)
# 0.3837576


xdata$z.rank.squared<-xdata$z.rank*xdata$z.rank

str(xdata)
table(xdata$rank)


aggregate(ind_in_patch~foc_patch_id, data=xdata, mean)
aggregate(dbh_cm~foc_patch_id, data=xdata, mean)
aggregate(fruit_crop_categ~foc_patch_id, data=xdata, mean)


mean(xdata$prop.into.bout)
#0.3837576

min(xdata$prop.into.bout)
#0

max(xdata$prop.into.bout)
#1

#----------RANDOM SLOPES----------------------


#check for random slopes:
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\diagnostic_fcns.R")

rand_slope_data=fe.re.tab(fe.model="intake_rate~z.log_ind_in*(z.rank+z.rank.squared) + z.log_dbh*(z.rank+z.rank.squared) + z.fruit_crop_categ*(z.rank+z.rank.squared) + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time", re="(1|bonobo_id) + (1|foc_patch_spp) + (1|foc_patch_id) + (1|bout_id)", data=xdata)

rand_slope_data$summary



#check which bout id  has two values for DBH :            <---------two bout ids changed in the original excel-template
aggregate(xdata$dbh_cm, list(xdata$bout_id), function(x){length(unique(x))})





#according to the summary, random slopes need to be added for: see printed list


#------------------------RUN MODEL---------------------------


library(lme4)


#random slope within foc patch species with ind-in-patch and rank interaction provided error msg,
#need to create new variables of their interaction:
xdata$i_rank_dbh=xdata$z.rank*xdata$z.log_dbh
xdata$i_rank_sq_dbh=xdata$z.rank.squared*xdata$z.log_dbh

xdata$i_rank_fruit=xdata$z.rank*xdata$z.fruit_crop_categ
xdata$i_rank_sq_fruit=xdata$z.rank.squared*xdata$z.fruit_crop_categ

xdata$i_rank_ind=xdata$z.rank*xdata$z.log_ind_in
xdata$i_rank_sq_ind=xdata$z.rank.squared*xdata$z.log_ind_in



#Model II: res poisson with interactions manually coded for random slopes
#need a control object before running the full model:
control_model=glmerControl(optCtrl = list(maxfun=100000), optimizer="bobyqa", calc.derivs=F)

start.time<-Sys.time()
res=glmer(intake_rate~z.log_ind_in*(z.rank+z.rank.squared) + z.log_dbh*(z.rank+z.rank.squared) + z.fruit_crop_categ*(z.rank+z.rank.squared) + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time +

(1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bonobo_id) + 

(1 + i_rank_dbh + i_rank_sq_dbh + i_rank_fruit + i_rank_sq_fruit + i_rank_ind + i_rank_sq_ind + z.log_ind_in + z.rank + z.rank.squared + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_spp) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_id) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bout_id),

data=xdata, family=poisson, control=control_model)

end.time<-Sys.time()
end.time - start.time
#difference of 30.68011 mins
summary(res)





#-----------MODEL ASSUMPTIONS/STABILITY-----------------------

#check that the model assumptions are fulfilled:
ranef.diagn.plot(res)
overdisp.test(res)

#      chisq     df     P       dispersion.parameter
#1  2392.907    4534    1            0.5277695
#> 




#then we check VIFs:
library(car)

vif_model=lm(intake_rate~z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.rank + z.prop.into.bout + z.new_time, data=xdata) 

vif(vif_model)

#   z.log_ind_in          z.log_dbh z.fruit_crop_categ             z.rank   z.prop.into.bout         z##.new_time 
#          1.233332           1.134660           1.158586           1.020197           1.083380         #  1.057398 




#next we check the model stability by sourcing Rogers function:
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\glmm_stability.r")
ls()

#head(glmm.model.stab)


model_stability=glmm.model.stab(model.res=res, contr=control_model, n.cores=6, para=T, data=xdata)
table(model_stability$detailed$warnings)
round(model_stability$summary[, -1], 3)
m.stab.plot(model_stability$summary[, -1])


#based on the stability, check the model 



#---------------FULL NULL COMPARISON---------------------



#the null model:
   null_model=glmer(intake_rate~z.new_time + z.sq_time + 

   (1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bonobo_id) + 

(1 + i_rank_dbh + i_rank_sq_dbh + i_rank_fruit + i_rank_sq_fruit + i_rank_ind + i_rank_sq_ind + z.log_ind_in + z.rank + z.rank.squared + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_spp) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_id) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bout_id), 

     data=xdata, family=poisson, control=control_model)





#then the full null comparison:

anova(null_model, res, test="Chisq")

#          Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#null_model 39 19126 19377 -9524.0    19048                             
#res        52 19109 19443 -9502.4    19005 43.257     13  4.073e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





#-----------------P-VALUES--------------------------------


#the full null comparison is significant, next, p-values for the interactions:

source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\drop1_para.r")
ls()
head(drop1p)
drop1.results.int.sq<-drop1p(model.res=res, data=xdata, contr=control_model, to.del=c("z.rank.squared:z.fruit_crop_categ","z.rank.squared:z.log_dbh","z.log_ind_in:z.rank.squared"))

#results for sq-rank interactions not significant:
#drop1.res
#$drop1.res
#                                     logLik      AIC       Chisq Chi.Df Pr..Chisq.    n n.opt.warnings
#none                              -9502.406 19108.81          NA     NA         NA 4550             NA
#z.rank.squared:z.fruit_crop_categ -9503.319 19108.64 1.825785990      1  0.1766261 4550              0
#z.rank.squared:z.log_dbh          -9503.180 19108.36 1.547262310      1  0.2135401 4550              0
#z.log_ind_in:z.rank.squared       -9502.409 19106.82 0.006600014      1  0.9352507 4550              0
#                                  n.fun.warnings
#none                                          NA
#z.rank.squared:z.fruit_crop_categ              0
#z.rank.squared:z.log_dbh                       0
#z.log_ind_in:z.rank.squared                    0




summary(res)
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                        1.5781886  0.2096250   7.529 5.13e-14 ***
#  z.log_ind_in                      -0.0126450  0.0158591  -0.797  0.42526    
#z.rank                             0.0336738  0.0123851   2.719  0.00655 ** 
#  z.rank.squared                    -0.0025975  0.0131854  -0.197  0.84383    
#z.log_dbh                          0.0052507  0.0292708   0.179  0.85764    
#z.fruit_crop_categ                 0.0492313  0.0327758   1.502  0.13308    
#z.prop.into.bout                  -0.0710533  0.0082682  -8.594  < 2e-16 ***
#  z.prop.into.bout.squared           0.0277021  0.0099498   2.784  0.00537 ** 
#  z.new_time                         0.0739131  0.0177376   4.167 3.09e-05 ***
#  z.sq_time                         -0.0267381  0.0196017  -1.364  0.17254    
#z.log_ind_in:z.rank               -0.0093655  0.0101714  -0.921  0.35717    
#z.log_ind_in:z.rank.squared       -0.0009439  0.0112091  -0.084  0.93289    
#z.rank:z.log_dbh                  -0.0131087  0.0128653  -1.019  0.30824    
#z.rank.squared:z.log_dbh          -0.0160019  0.0127444  -1.256  0.20926    
#z.rank:z.fruit_crop_categ          0.0261008  0.0140266   1.861  0.06277 .  
#z.rank.squared:z.fruit_crop_categ  0.0186199  0.0131764   1.413  0.15762    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1








#rerun drop 1 on RES to check p-value of squared proportion into bout:
drop1.results.bout.sq<-drop1p(model.res=res, data=xdata, contr=control_model, to.del=c("z.prop.into.bout.squared"))


#the p-value for propoprtion into bout squared was highly significant so we keep in the model for now:

#$drop1.res
#                            logLik      AIC    Chisq Chi.Df  Pr..Chisq.    n n.opt.warnings n.fun#.warnings
#none                     -9502.406 19108.81       NA     NA          NA 4550             NA            # NA
#z.prop.into.bout.squared -9505.886 19113.77 6.959219      1 0.008338841 4550              0              0

#$model.results
#NULL



#rerun model by dropping the squared rank interactions but keeping into bout squared:
start.time<-Sys.time()

res_no_sq_int=glmer(intake_rate~z.log_ind_in*z.rank + z.log_dbh*z.rank + z.fruit_crop_categ*z.rank + z.rank.squared + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time +

(1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bonobo_id) + 

(1 + i_rank_dbh + i_rank_fruit + i_rank_ind + z.log_ind_in + z.rank + z.rank.squared + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_spp) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_id) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bout_id),

data=xdata, family=poisson, control=control_model)



#next we need p-values for the rank interactions and rank squared:

drop1.results.int.ranksq<-drop1p(model.res=res_no_sq_int, data=xdata, contr=control_model, to.del=c("z.rank:z.fruit_crop_categ","z.rank:z.log_dbh","z.log_ind_in:z.rank", "z.rank.squared"))



#$drop1.res
#                             logLik      AIC     Chisq Chi.Df Pr..Chisq.    n n.opt.warnings n.fun#.warnings
#none                      -9503.647 19099.29        NA     NA         NA 4550             NA           #  NA
#z.rank:z.fruit_crop_categ -9504.036 19098.07 0.7775934      1  0.3778782 4550              0           #   0
#z.rank:z.log_dbh          -9503.708 19097.42 0.1215628      1  0.7273457 4550              0           #   0
#z.log_ind_in:z.rank       -9503.856 19097.71 0.4170211      1  0.5184272 4550              0           #   0
#z.rank.squared            -9503.785 19097.57 0.2768260      1  0.5987891 4550              0           #   0
#
#$model.results
#NULL





#none of the p-values for the rank interactions are significant and neither is the rank squared, thus we run he 
#model again without the interactions


summary(res_no_sq_int)

#Random effects:
#Number of obs: 4550, groups:  bout_id, 404; foc_patch_id, 328; foc_patch_spp, 18; bonobo_id, 14

#Fixed effects:
#                         Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                1.569037   0.208791   7.515 5.70e-14 ***
#  z.log_ind_in              -0.014214   0.012013  -1.183 0.236722    
#z.rank                     0.039156   0.011832   3.309 0.000936 ***
#  z.log_dbh                 -0.012995   0.025966  -0.500 0.616735    
#z.fruit_crop_categ         0.067549   0.030104   2.244 0.024843 *  
#  z.rank.squared             0.006208   0.011003   0.564 0.572597    
#z.prop.into.bout          -0.071621   0.008284  -8.645  < 2e-16 ***
#  z.prop.into.bout.squared   0.027719   0.010119   2.739 0.006159 ** 
#  z.new_time                 0.073311   0.017992   4.075 4.61e-05 ***
#  z.sq_time                 -0.025646   0.019146  -1.340 0.180394    
#z.log_ind_in:z.rank       -0.008328   0.009594  -0.868 0.385387    
#z.rank:z.log_dbh          -0.004418   0.011069  -0.399 0.689818    
#z.rank:z.fruit_crop_categ  0.014816   0.012071   1.227 0.219685    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





#-----------------New reduced model for 17.11.-------------------------------------------------------------


res_drop_inter=glmer(intake_rate~z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.rank + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time +

(1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bonobo_id) + 

(1 + z.log_ind_in + z.rank + z.log_dbh + z.fruit_crop_categ + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_spp) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_id) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bout_id),

data=xdata, family=poisson, control=control_model)


 

#then the p-values for the new model:
drop1.results.final<-drop1p(model.res=res_drop_inter, data=xdata, contr=control_model)
drop1.results.final


#----------NEW P-values for the reduced model WITHOUT interactions 17.11.-----------------------------
round(drop1.results.final$drop1.res, 4)

#                            logLik      AIC   Chisq Chi.Df Pr..Chisq.    n n.opt.warnings n.fun.warnings
#none                     -9504.520 19085.04      NA     NA         NA 4550             NA            # NA
#z.log_ind_in             -9505.224 19084.45  1.4081      1     0.2354 4550              0             # 0
#z.log_dbh                -9504.635 19083.27  0.2300      1     0.6315 4550              0              0
#z.fruit_crop_categ       -9506.520 19087.04  4.0003      1     0.0455 4550              0              0
#z.rank                   -9507.124 19088.25  5.2089      1     0.0225 4550              0              0
#z.prop.into.bout         -9516.819 19107.64 24.5985      1     0.0000 4550              0              0
#z.prop.into.bout.squared -9507.882 19089.76  6.7252      1     0.0095 4550              0              0
#z.new_time               -9510.100 19094.20 11.1609      1     0.0008 4550              0              0
#z.sq_time                -9505.359 19084.72  1.6783      1     0.1951 4550              0              0
#>


#rank is significant (p=0.0225) with a positive correlation with intake
#Proportion in to bout is also significant as squared (non linear)term (p=0.0095)
#
#direction of effect:
summary(res_drop_inter)


#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )

#Fixed effects:
#                          Estimate Std. Error z value Pr(>|z|)    
#(Intercept)               1.574672   0.208683   7.546 4.50e-14 ***
#z.log_ind_in             -0.014294   0.011815  -1.210  0.22633    
#z.log_dbh                -0.012613   0.025428  -0.496  0.61987    
#z.fruit_crop_categ        0.067684   0.030143   2.245  0.02474 *  
#z.rank                    0.034236   0.010693   3.202  0.00137 ** 
#z.prop.into.bout         -0.071284   0.008328  -8.559  < 2e-16 ***
#z.prop.into.bout.squared  0.028434   0.010469   2.716  0.00660 ** 
#z.new_time                0.073023   0.018312   3.988 6.67e-05 ***
#z.sq_time                -0.025097   0.018244  -1.376  0.16893    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



intake_resid=resid(res_drop_inter)
plot(xdata$rank, intake_resid, ylab = "Residuals intake", xlab="Rank")
plot.coeffs<-fixef(res_drop_inter)

plot.xvals<-seq(min(xdata$z.rank), max(xdata$z.rank), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.rank"]*plot.xvals
plot.yvals<-exp(plot.yvals)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)


############################################################################################
############################################################################################
#PLOT INTAKE VS PROP INTO BOUT 25.01. :
#D:/Users/niina_nurmi/Documents/research_2017_pan_manuscripts/manuscript_i_energy_balance/ms_i_behavioural_ecology
tiff("D:/Users/niina_nurmi/Documents/research_2017_pan_manuscripts/manuscript_i_energy_balance/ms_i_behavioural_ecology/figure_intake_prop.tiff", width = 6, height = 6, units = 'in', res = 300)
plot.new()


dev.off()



tiff("D:/Users/niina_nurmi/Documents/research_2017_pan_manuscripts/manuscript_i_energy_balance/ms_i_behavioural_ecology/figure_intake_prop.tiff", width = 6, height = 6, units = 'in', res = 300)
plot.new()

par(mar=c(6.5, 6.5, 0.5, 0.5), mgp=c(4.1, 0.8, 0))
plot(x=xdata$z.prop.into.bout, y=xdata$intake_rate, pch=19, col=adjustcolor("blue", alpha=0.1), xlab="Proportion of time until end of bout", ylab="Intake rate per min", xaxt="n", las=1, cex.lab=2.1, cex.axis=1.5)
axis(side=1, at=(c(0, 0.20, 0.40, 0.60, 0.80, 1)-mean(xdata$prop.into.bout))/sd(xdata$prop.into.bout), labels=c(0, 0.20, 0.40, 0.60, 0.80, 1), cex.lab=2.1, cex.axis=1.5)

plot.coeffs<-fixef(res_drop_inter)

plot.xvals<-seq(min(xdata$z.prop.into.bout), max(xdata$z.prop.into.bout), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.prop.into.bout"]*plot.xvals+plot.coeffs["z.prop.into.bout.squared"]*(plot.xvals^2)
plot.yvals<-exp(plot.yvals)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)

dev.off()









par(mar=c(6.5, 6.5, 0.5, 0.5), mgp=c(4.1, 0.8, 0))
plot(x=xdata$z.prop.into.bout, y=xdata$intake_rate, pch=19, col=adjustcolor("blue", alpha=0.1), xlab="Proportion of time until end of bout", ylab="Intake rate per min", xaxt="n", las=1, cex.lab=2.1, cex.axis=1.5)
axis(side=1, at=(c(0, 0.20, 0.40, 0.60, 0.80, 1)-mean(xdata$prop.into.bout))/sd(xdata$prop.into.bout), labels=c(0, 0.20, 0.40, 0.60, 0.80, 1), cex.lab=2.1, cex.axis=1.5)

plot.coeffs<-fixef(res_drop_inter)

plot.xvals<-seq(min(xdata$z.prop.into.bout), max(xdata$z.prop.into.bout), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.prop.into.bout"]*plot.xvals+plot.coeffs["z.prop.into.bout.squared"]*(plot.xvals^2)
plot.yvals<-exp(plot.yvals)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)


#check range of proportion into bout, 
#too many unique, we need to bin values
range(xdata$z.prop.into.bout)

binwidth=0.08 


#note this was originally : binwidth=0.1 

xdata$time_bins=xdata$prop.into.bout %/% binwidth

xdata$time_bins=xdata$time_bins*binwidth

xdata$time_bins=xdata$time_bins+binwidth/2


unique(xdata$time_bins)


#[1] 0.05 0.15 0.25 0.55 0.65 0.35 0.45 0.75 0.85 0.95 1.05

#looks pretty good
#now need to find proportion of what is the response per bin


proportion_bin=aggregate(x=xdata$intake_rate, by=list(xdata$time_bins, xdata$foc_patch_spp), FUN=mean)  

head(proportion_bin)

#now turn the bin values to z space equalents:
proportion_bin$z_bin=(proportion_bin$Group.1-mean(xdata$prop.into.bout))/sd(xdata$prop.into.bout)


#for the size of data points create the object with number of obs per bin:
data_point_size=aggregate(xdata$intake_rate, by=list(xdata$time_bins, xdata$foc_patch_spp), FUN=length)

min(proportion_bin$z_bin)
max(proportion_bin$z_bin)


#plotting the results:
tiff("D:/Users/niina_nurmi/Documents/research_2017_pan_manuscripts/manuscript_i_energy_balance/ms_i_behavioural_ecology/figure_2b_tiff_intake_prop.tiff", width = 7, height = 7.2, units = 'in', res = 350)
plot.new()


par(mar=c(5.9, 6.3, 0.5, 1.5), mgp=c(4.6, 0.9, 0), tcl=-0.3)
plot(x=proportion_bin$z_bin, y=proportion_bin$x, pch=19, cex=data_point_size$x^(1/4), col=adjustcolor("black", alpha=0.5),
     xlab="Proportion of time until end of bout", ylab="Intake rate per minute", xaxt="n", las=1, cex.lab=2.3, cex.axis=1.9)
axis(side=1, at=(c(0, 0.20, 0.40, 0.60, 0.80, 1)-mean(xdata$prop.into.bout))/sd(xdata$prop.into.bout), 
     labels=c(0, 0.20, 0.40, 0.60, 0.80, 1), cex.lab=2.3, cex.axis=1.9)


plot.coeffs<-fixef(res_drop_inter)

plot.xvals<-seq(min(xdata$z.prop.into.bout), max(xdata$z.prop.into.bout), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.prop.into.bout"]*plot.xvals+plot.coeffs["z.prop.into.bout.squared"]*(plot.xvals^2)
plot.yvals<-exp(plot.yvals)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)

dev.off()


#savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\xxxxx.png", type="png")




#----------------------------------------------------------------------------------------------------------------------
#add species categories
xdata$sp.category<-factor(c("Pods","Fruits","Fruits","Pods","Terrestrial","Fruits","Leaves","Fruits","Fruits","Fruits","Fruits","Fruits","Fruits","Leaves","Dialium","Fruits","Pods","Terrestrial")[as.numeric(xdata$foc_patch_spp)])

proportion_bin$sp.category<-xdata$sp.category[match(proportion_bin$Group.2, xdata$foc_patch_spp)]
proportion_bin$color<-c("black","black","black","blue","black")[as.numeric(proportion_bin$sp.category)]

par(mar=c(6.5, 6.5, 0.5, 0.5), mgp=c(4.6, 0.9, 0))
plot(x=proportion_bin$z_bin, y=proportion_bin$x, pch=19, cex=data_point_size$x^(1/2.8), col=adjustcolor(proportion_bin$color, alpha=0.7), xlab="Proportion of time until end of FT bout", ylab="Intake rate per minute", xaxt="n", las=1, cex.lab=2.1, cex.axis=1.5)
axis(side=1, at=(c(0, 0.20, 0.40, 0.60, 0.80, 1)-mean(xdata$prop.into.bout))/sd(xdata$prop.into.bout), labels=c(0, 0.20, 0.40, 0.60, 0.80, 1), cex.lab=2.1, cex.axis=1.5)
lines(plot.xvals, plot.yvals, lwd=4, lty=2)


#legend
legend("topright", legend=c("Dialium","Fruits","Leaves" ,"Pods","Terrestrial"), col=adjustcolor(c("black","orange","green","brown","grey"), alpha=0.6), bty="n", pch=19, pt.cex=2, text.font=c(3,1,1,1,1))

#savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\xxxxxxxxx.png", type="png")




#------------------------------------------------------------------------------
###with species lines included :


ranef(res_drop_inter)
sp.colors<-proportion_bin[!duplicated(proportion_bin[,c("Group.2","sp.category","color")]),]

for(i in 1:nrow(ranef(res_drop_inter)$foc_patch_spp)){
    sp.yvals<-plot.coeffs["(Intercept)"] + ranef(res_drop_inter)$foc_patch_spp[i,"(Intercept)"] + plot.coeffs["z.prop.into.bout"]*plot.xvals + (plot.coeffs["z.prop.into.bout.squared"]+ranef(res_drop_inter)$foc_patch_spp[i,"z.prop.into.bout.squared"])*(plot.xvals^2)
    sp.yvals<-exp(sp.yvals)
    lines(plot.xvals, sp.yvals, lwd=2, lty=2, col=sp.colors$color[i])
}

savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_vs_prop_into_bout_spp_colours_lines.png", type="png")




#################################################################################
#################################################################################
#
#PLOT INTAKE VS RANK 25.01:

#plotting the results:

#create data point sizes:
data_point_size_rank=aggregate(xdata$intake_rate, by=list(xdata$bonobo_id, xdata$rank), FUN=length)
#lets look at it:
data_point_size_rank

str(data_point_size_rank)

data_point_size_rankmean=aggregate(xdata$intake_rate, by=list(xdata$bonobo_id, xdata$rank), FUN=mean)

#now turn the bin values to z space equalents:
data_point_size_rankmean$z_rank=(data_point_size_rankmean$Group.2-mean(xdata$rank))/sd(xdata$rank)



#then plot
tiff("D:/Users/niina_nurmi/Documents/research_2017_pan_manuscripts/manuscript_i_energy_balance/ms_i_behavioural_ecology/figure_2a_tiff_intake_rank.tiff", width = 7, height = 7.2, units = 'in', res = 350)
plot.new()

par(mar=c(5.9, 6.3, 0.5, 1.5), mgp=c(4.6, 0.9, 0), tcl=-0.3)
plot(x=data_point_size_rankmean$z_rank, y=data_point_size_rankmean$x, pch=19, cex=data_point_size_rank$x^(1/4), 
     col=adjustcolor("black", alpha=0.5), xlab="Increasing female dominance rank", ylab="Intake rate per minute", 
     xaxt="n", las=1, ylim=c(0, 14), cex.lab=2.3, cex.axis=1.9)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1), cex.lab=2.3, cex.axis=1.9)

plot.coeffs<-fixef(res_drop_inter)

plot.xvals<-seq(min(xdata$z.rank), max(xdata$z.rank), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.rank"]*plot.xvals
plot.yvals<-exp(plot.yvals)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)

dev.off()


#---------------------------------------------------------------
#by species and individual

par(mar=c(5.1, 5.8, 0.5, 0.5), mgp=c(3.7, 1.0, 0))
intake.per.sp.ind<-aggregate(xdata$intake_rate, list(xdata$bonobo_id, xdata$foc_patch_spp, xdata$rank), mean)
intake.per.sp.ind$z.rank<-(intake.per.sp.ind$Group.3-mean(xdata$rank))/sd(xdata$rank)

sample.size.intake.per.sp.ind<-aggregate(xdata$intake_rate, list(xdata$bonobo_id, xdata$foc_patch_spp, xdata$rank), length)

plot(x=intake.per.sp.ind$z.rank, y=intake.per.sp.ind$x, pch=19, cex=sample.size.intake.per.sp.ind$x^(1/3.5), col=adjustcolor("red", alpha=0.8), xlab="Female rank", ylab="Intake rate per minute", xaxt="n", las=1, cex.lab=1.9, cex.axis=1.5)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1), cex.lab=1.9, cex.axis=1.4)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)


#---------------------------------------------------------------------
#with colors by species category

intake.per.sp.ind$sp.category<-xdata$sp.category[match(intake.per.sp.ind$Group.2, xdata$foc_patch_spp)]
intake.per.sp.ind$color<-c("grey40","firebrick4","green4","goldenrod2","deepskyblue3")[as.numeric(intake.per.sp.ind$sp.category)]

plot(x=intake.per.sp.ind$z.rank, y=intake.per.sp.ind$x, pch=19, cex=sample.size.intake.per.sp.ind$x^(1/3.5), col=adjustcolor(intake.per.sp.ind$color, alpha=0.8), xlab="Female rank", ylab="Intake rate per minute", xaxt="n", las=1, cex.lab=1.9, cex.axis=1.4)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1), cex.lab=1.9, cex.axis=1.4)
lines(plot.xvals, plot.yvals, lwd=3, lty=2)

#legend
legend("topleft", legend=c("Dialium","Fruits","Leaves" ,"Pods","Terrestrial"), col=adjustcolor(c("grey40","firebrick4","green4","goldenrod2","deepskyblue3"), alpha=0.8), bty="n", pch=19, pt.cex=2, text.font=c(3,1,1,1,1))



for(i in 1:nrow(ranef(res_drop_inter)$foc_patch_spp)){
    current.species<-rownames(ranef(res_drop_inter)$foc_patch_spp)[i]
    sp.xvals<-seq(min(xdata$z.rank[xdata$foc_patch_spp==current.species]), max(xdata$z.rank[xdata$foc_patch_spp==current.species]), length.out=100)

    sp.yvals<-plot.coeffs["(Intercept)"] + ranef(res_drop_inter)$foc_patch_spp[i,"(Intercept)"] + (plot.coeffs["z.rank"] + ranef(res_drop_inter)$foc_patch_spp[i,"z.rank"])*plot.xvals
    sp.yvals<-exp(sp.yvals)
#    lines(plot.xvals, sp.yvals, lwd=2, lty=2, col=sp.colors$color[i])
  #  lines(sp.xvals, sp.yvals, lwd=2, lty=2, col=sp.colors$color[i])
    lines(sp.xvals, sp.yvals, lwd=2, lty=2, col=intake.per.sp.ind$color[i])
}




######################################################### oma vari punanen
 #######################################################
#with colors by species category

intake.per.sp.ind$sp.category<-xdata$sp.category[match(intake.per.sp.ind$Group.2, xdata$foc_patch_spp)]
intake.per.sp.ind$color<-c("black","black","black","red","black")[as.numeric(intake.per.sp.ind$sp.category)]

plot(x=intake.per.sp.ind$z.rank, y=intake.per.sp.ind$x, pch=19, cex=sample.size.intake.per.sp.ind$x^(1/3.5), col=adjustcolor(intake.per.sp.ind$color, alpha=0.7), xlab="Female rank", ylab="Intake rate per minute", xaxt="n", las=1, cex.lab=1.9, cex.axis=1.4)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1), cex.lab=1.9, cex.axis=1.4)
lines(plot.xvals, plot.yvals, lwd=3, lty=2)

#legend
legend("topleft", legend=c("Dialium","Fruits","Leaves" ,"Pods","Terrestrial"), col=adjustcolor(c("grey40","firebrick4","green4","goldenrod2","deepskyblue3"), alpha=0.8), bty="n", pch=19, pt.cex=2, text.font=c(3,1,1,1,1))



savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_rank_plot_color_new_species.png", type="png")

savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_rank_plot_color_new_species_lines_added.png", type="png")

savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_rank_plot_color_new_spp_lines_added_short.png", type="png")




##########################################################################################################################
#####################################################################################################################
#####################################################################################################################
#PLOT INTAKE vs FRUIT CROP SIZE 25.01


#create data point sizes:
data_point_size_crop=aggregate(xdata$intake_rate, by=list(xdata$bonobo_id, xdata$fruit_crop_categ), FUN=length)
#lets look at it:
data_point_size_crop

str(data_point_size_crop)

data_point_size_cropmean=aggregate(xdata$intake_rate, by=list(xdata$bonobo_id, xdata$fruit_crop_categ), FUN=mean)

#now turn the bin values to z space equalents:
data_point_size_cropmean$z_crop=(data_point_size_cropmean$Group.2-mean(xdata$fruit_crop_categ))/sd(xdata$fruit_crop_categ)

range(data_point_size_cropmean$x)

#then plot
par(mar=c(5.1, 5.8, 0.5, 2.5), mgp=c(3.8, 0.8, 0))

plot(x=data_point_size_cropmean$z_crop, y=data_point_size_cropmean$x, pch=19, cex=data_point_size_crop$x^(1/4), col=adjustcolor("darkorange", alpha=0.7), xlab="FT food crop size", ylab="Intake rate per minute", xaxt="n", las=1, ylim=c(0, 18), cex.lab=1.9, cex.axis=1.5)
axis(side=1, at=(c(2, 3, 4, 5)-mean(xdata$fruit_crop_categ))/sd(xdata$fruit_crop_categ), labels=c("10 - 99", "100 - 999", "1000 - 9999", "10000 - 99999"), cex.axis=1.1, cex.lab=1.8)

plot.coeffs<-fixef(res_drop_inter)

plot.xvals<-seq(min(xdata$z.fruit_crop_categ), max(xdata$z.fruit_crop_categ), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.fruit_crop_categ"]*plot.xvals
plot.yvals<-exp(plot.yvals)

lines(plot.xvals, plot.yvals, lwd=3, lty=2)





#---------------------------------------------------------------------------------------------------------------
#by FRUIT CROP CATEG by species and individual

intake.crop.per.sp.ind<-aggregate(xdata$intake_rate, list(xdata$bonobo_id, xdata$foc_patch_spp, xdata$fruit_crop_categ), mean)

intake.crop.per.sp.ind$z.crop<-(intake.crop.per.sp.ind$Group.3-mean(xdata$fruit_crop_categ))/sd(xdata$fruit_crop_categ)

sample.size.intake.crop.per.sp.ind<-aggregate(xdata$intake_rate, list(xdata$bonobo_id, xdata$foc_patch_spp, xdata$fruit_crop_categ), length)


intake.crop.per.sp.ind$sp.category<-xdata$sp.category[match(intake.crop.per.sp.ind$Group.2, xdata$foc_patch_spp)]

intake.crop.per.sp.ind$color<-c("black","orange","green","brown","grey")[as.numeric(intake.crop.per.sp.ind$sp.category)]


plot(x=intake.crop.per.sp.ind$z.crop, y=intake.crop.per.sp.ind$x, pch=19, cex=sample.size.intake.crop.per.sp.ind$x^(1/3), col=adjustcolor(intake.crop.per.sp.ind$color, alpha=0.6), xlab="Fruit crop category", ylab="Intake rate per minute", xaxt="n", las=1)
axis(side=1, at=(c(2, 3, 4, 5)-mean(xdata$fruit_crop_categ))/sd(xdata$fruit_crop_categ), labels=c(2, 3, 4, 5))
lines(plot.xvals, plot.yvals, lwd=3, lty=2)


#legend
legend("topleft", legend=c("Dialium","Fruits","Leaves" ,"Pods","Terrestrial"), col=adjustcolor(c("black","orange","green","brown","grey"), alpha=0.6), bty="n", pch=19, pt.cex=2, text.font=c(3,1,1,1,1))


#for the species lines, check fixed effects:
fixef(res_drop_inter)

ranef(res_drop_inter)
#no random slope so all same line above or below the main line...


for(i in 1:nrow(ranef(res_drop_inter)$foc_patch_spp)){
    current.species<-rownames(ranef(res_drop_inter)$foc_patch_spp)[i]
    sp.xvals<-seq(min(xdata$z.fruit_crop_categ[xdata$foc_patch_spp==current.species]), max(xdata$z.fruit_crop_categ[xdata$foc_patch_spp==current.species]), length.out=100)
    sp.yvals<-plot.coeffs["(Intercept)"] + ranef(res_drop_inter)$foc_patch_spp[i,"(Intercept)"] + (plot.coeffs["z.fruit_crop_categ"] + ranef(res_drop_inter)$foc_patch_spp[i,"z.fruit_crop_categ"])*sp.xvals
    sp.yvals<-exp(sp.yvals)
    lines(sp.xvals, sp.yvals, lwd=2, lty=2, col=sp.colors$color[i])
}



savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_crop_plot_color_new_spp_lines_added.png", type="png")

savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_crop_plot_color_new_spp_lines_added_short.png", type="png")







##################plot fruit crop category against residuals of model without it instead
no.crop<-glmer(intake_rate~z.log_ind_in + z.log_dbh + z.rank + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time +

(1 + z.log_ind_in + z.log_dbh + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bonobo_id) + 

(1 + z.log_ind_in + z.rank + z.log_dbh +z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_spp) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||foc_patch_id) + 

(1 + z.log_ind_in + z.prop.into.bout + z.prop.into.bout.squared + z.new_time + z.sq_time||bout_id),

data=xdata, family=poisson, control=control_model)


plot(x=xdata$fruit_crop_categ, y=residuals(no.crop), pch=19, col=adjustcolor("black", alpha=0.05))





simple.res<-glmer(intake_rate~z.fruit_crop_categ + 

(1 + z.fruit_crop_categ||bonobo_id) + 

(1 + z.fruit_crop_categ||foc_patch_spp) + 

(1|foc_patch_id) + 

(1|bout_id),

data=xdata, family=poisson, control=control_model)


fixef(simple.res)
fixef(res_drop_inter)

simple.res.no.rs<-glmer(intake_rate~z.fruit_crop_categ + 

(1 |bonobo_id) + 

(1|foc_patch_spp) + 

(1|foc_patch_id) + 

(1|bout_id),

data=xdata, family=poisson, control=control_model)

fixef(simple.res.no.rs)
fixef(res_drop_inter)





#confidence intervals for plot
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\boot_glmm.r")
rank.plot.ci<-boot.glmm.2(res_no_inter, use=c("z.rank"), para=T, n.boots=1000, contr=control_model)

#what to do when the term is squared?
time.until.end.plot.ci<-boot.glmm.2(res_no_inter, use=c("z.time_elap_end"), para=T, n.boots=1000, contr=control_model)








##########################
#----NEW plot intake vs rank different colours per species 10.11.16-----------

#new by species and individual:

intake.per.sp.ind<-aggregate(xdata$intake_rate, list(xdata$bonobo_id, xdata$foc_patch_spp[], xdata$rank), mean)
intake.per.sp.ind$z.rank<-(intake.per.sp.ind$Group.3-mean(xdata$rank))/sd(xdata$rank)

sample.size.intake.per.sp.ind<-aggregate(xdata$intake_rate, list(xdata$bonobo_id, xdata$foc_patch_spp, xdata$rank), length)

plot(x=intake.per.sp.ind$z.rank, y=intake.per.sp.ind$x, pch=19, cex=sample.size.intake.per.sp.ind$x^(1/3), col=adjustcolor("black", alpha=0.7), xlab="Female rank", ylab="Intake rate per minute", xaxt="n", las=1)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1))

lines(plot.xvals, plot.yvals, lwd=3, lty=2)

#legend
legend("topleft", legend=c("Dialium","Fruits","Leaves" ,"Pods","Terrestrial"), col=adjustcolor(c("gray1","tomato2","springgreen4","tan4","yellow3"), alpha=0.7), bty="n", pch=19, pt.cex=2, text.font=c(3,1,1,1,1))

#with colors by species category

intake.per.sp.ind$sp.category<-xdata$sp.category[match(intake.per.sp.ind$Group.2, xdata$foc_patch_spp)]
intake.per.sp.ind$color<-c("gray1","tomato2","springgreen4","tan4","yellow3")[as.numeric(intake.per.sp.ind$sp.category)]

plot(x=intake.per.sp.ind$z.rank, y=intake.per.sp.ind$x, pch=19, cex=sample.size.intake.per.sp.ind$x^(1/3), col=adjustcolor(intake.per.sp.ind$color, alpha=0.7), xlab="Female rank", ylab="Intake rate per minute", xaxt="n", las=1)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1))
lines(plot.xvals, plot.yvals, lwd=3, lty=2)

#legend
legend("topleft", legend=c("Dialium","Fruits","Leaves" ,"Pods","Terrestrial"), col=adjustcolor(c("gray1","tomato2","springgreen4","tan4","yellow3"), alpha=0.7), bty="n", pch=19, pt.cex=2, text.font=c(3,1,1,1,1))


source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\boot_glmm.r")

intake.per.sp.ind$sp.category<-xdata$sp.category[match(intake.per.sp.ind$Group.2, xdata$foc_patch_spp)]
intake.per.sp.ind$color<-c("gray1","tomato2","springgreen4","tan4","yellow3")[as.numeric(intake.per.sp.ind$sp.category)]
















savePlot("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_II_intake\\\\intake_rank_plot_species_diff_colours.png", type="png")






############################################################################
save.image("Y:/primint/Niina/Files_for_R_analysis_II_intake/MAIN_pan_analysis_II_intake_poisson_workspace_18112016.Rdata")
load("Y:/primint/Niina/Files_for_R_analysis_II_intake/MAIN_pan_analysis_II_intake_poisson_workspace_18112016.Rdata")