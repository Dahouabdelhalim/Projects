#-----MODEL 3. Feeding efficiency and feeding competition in bonobo females 26.09.2016-----
#This model predicts feeding efficiency of females via their movement (meters/min) in feeding patches
#as the response. The main test predictors are rank and its interactions
#with number of individuals in patch, patch size and number of fruits in patch.
#To account for depletion effect, we have time until end of bout as an
#additional test predictor, including its square term to account for non-linear depletion effect.
#We also have time of day (elapsed since midnight) as a fixed control predictor and bonobo ID,
#tree species and bout ID as a random control predictors.


#set up the working directory for MODEL 3:

setwd("Y:/primint/Niina/Files_for_R_analysis_III_movement")
getwd()


#check the list to make sure its empty
ls()
rm(list=ls())

#then read in the movement data from Primint as "xdata":
xdata=read.table(file="MASTER_analyses_III_R_file_movement_f.txt", header=T, sep="\\t")

str(xdata)


#-------------INSPECT DATA-----------------------------------------------------


table(xdata$bonobo_id)

#check the assignments per rank:
xx=table(xdata$bonobo_id, xdata$rank)

range(apply(X=xx>0, MARGIN=1, FUN=sum))

table(xdata$fruit_crop_categ)

table(xdata$ind_in_patch)

table(xdata$rank)



#--------------DISTRIBUTION/TRANSFORM--------------


#check the distribution of the response:

hist(xdata$change_poss_m, xlab="Meters moved during feeding focal minute", ylab="Frequency", col="grey")

#zero inflated, highly right skewed distribution

#as a table:

table(xdata$change_poss_m)
#   0    0.5   1   1.5    2    3    4    5    6    7    8   10   15   20   25 
# 3703   59  303    4   178   96   68   36   16    7    3   15    4    6    2 




hist(log(xdata$change_poss_m + 1))


#must make new response where zero is zero and 
#everything else becomes 1:

xdata$move_y_n=as.numeric(xdata$change_poss_m>0) 

head(xdata$move_y_n)
table(xdata$move_y_n)

#   0    1 
#3703  797

#looks good so the new response variable was created.


#these are the new rank categories corrected on the 1.11. this is why model was rerun for 3rd time on 1.11.:
# 0       0.1428571   0.2857143   0.4285714   0.5714286   0.7142857   0.8571429         1 
# 223       231       343         1010        266         118         934              1375 




#next check the predictor distribution:
hist(xdata$ind_in_patch)
hist(xdata$dbh_cm)
hist(xdata$fruit_crop_categ)


hist(log(xdata$ind_in_patch))
hist(log(xdata$dbh_cm))

#log transform predictors:

xdata$log_dbh=log(xdata$dbh_cm)
xdata$log_ind_in=log(xdata$ind_in_patch)




str(xdata)


#the time needs to be transformed, as entered in separate columns
#will make it a new column "time elapsed since midnight":

xdata$new_time=xdata$hour*3600 + xdata$min*60 + xdata$sec

hist(xdata$new_time)

#new time is now time elapsed since midnight in seconds.
#we also have the "end of bout time" in three colums,
#which we need to compare to our focal minute;

xdata$new_end_time=xdata$end_hour*3600 + xdata$end_min*60 + xdata$end_sec
str(xdata)


#then make the time elapsed variable:
xdata$time_elap_end=xdata$new_end_time-xdata$new_time

plot(xdata$time_elap_end)


#----------Z-TRANSFORM------------------------


#next we z-transform the predictor covariates:
xdata$z.log_ind_in=as.vector(scale(xdata$log_ind_in))       
xdata$z.rank=as.vector(scale(xdata$rank))
xdata$z.log_dbh=as.vector(scale(xdata$log_dbh))
xdata$z.fruit_crop_categ=as.vector(scale(xdata$fruit_crop_categ))
xdata$z.time_elap_end=as.vector(scale(xdata$time_elap_end)) 
xdata$z.new_time=as.vector(scale(xdata$new_time)) 



#then create from the z-transformed time elapsed since midnight a new 
#squared term to account for satiation.
#but first have to z-transform before squaring the time!
#create the same for the time until end of bout.
xdata$z.sq_time=xdata$z.new_time*xdata$z.new_time
xdata$z.sq_time_elap_end=xdata$z.time_elap_end*xdata$z.time_elap_end

str(xdata)


#create bout id by combining date to consecutive running nth focal tree of day:

xdata$bout_id=paste(xdata$year, xdata$month, xdata$day, xdata$foc_patch_no, sep="/")



#----------RANDOM SLOPES----------------------


#check for random slopes:
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\diagnostic_fcns.R")

rand_slope_data=fe.re.tab(fe.model="move_y_n~z.log_ind_in*z.rank + z.log_dbh*z.rank + z.fruit_crop_categ*z.rank + z.time_elap_end + z.sq_time_elap_end + z.new_time + z.sq_time", re="(1|bonobo_id) + (1|foc_patch_spp) + (1|foc_patch_id) + (1|bout_id)", data=xdata)

rand_slope_data$summary





#------------------------RUN MODEL---------------------------


library(lme4)

#need a control object before running the full model:
control_model=glmerControl(optCtrl = list(maxfun=100000), optimizer="bobyqa", calc.derivs=F)





#random slope within foc patch species with ind-in-patch and rank interaction provided error msg,
#need to create new variables of their interaction:
xdata$i_rank_dbh=xdata$z.rank*xdata$z.log_dbh
xdata$i_rank_fruit=xdata$z.rank*xdata$z.fruit_crop_categ
xdata$i_rank_ind=xdata$z.rank*xdata$z.log_ind_in



#Model III: Res WITH interactions manually coded for random slopes


start.time<-Sys.time()
res=glmer(move_y_n~z.log_ind_in*z.rank + z.log_dbh*z.rank + z.fruit_crop_categ*z.rank + z.time_elap_end + z.sq_time_elap_end + z.new_time + z.sq_time +
           
           (1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.time_elap_end + z.new_time + z.sq_time||bonobo_id) + 
           
           (1 + i_rank_dbh + i_rank_fruit + i_rank_ind + z.log_ind_in + z.rank + z.log_dbh + z.fruit_crop_categ + z.time_elap_end + z.new_time + z.sq_time||foc_patch_spp) + 
           
           (1 + z.log_ind_in + z.time_elap_end + z.new_time + z.sq_time||foc_patch_id) + 
           
           (1 + z.log_ind_in + z.time_elap_end + z.new_time + z.sq_time||bout_id),
         
         data=xdata, family=binomial, control=control_model)

end.time<-Sys.time()
end.time - start.time


summary(res)

#---------------------------IGNORE---------------------------------------------

#we include summary of original full model because we need to 
#compare these results to a rerun where we have added the CALC DERIVS term:

#Fixed effects:
#                          Estimate Std. Error z value Pr(>|z|)    
#(Intercept)               -1.39504    0.19757  -7.061 1.65e-12 ***
#z.log_id_in              -0.11649    0.09850  -1.183  0.23697    
#z.rank                    -0.17546    0.07777  -2.256  0.02407 *  
#z.log_dbh                  0.09970    0.14280   0.698  0.48507    
#z.fruit_crop_categ         0.04545    0.13541   0.336  0.73712    
#z.time_elap_end           -0.31411    0.10091  -3.113  0.00185 ** 
#z.sq_time_elap_end         0.07827    0.05649   1.385  0.16590    
#z.new_time                -0.12369    0.06574  -1.882  0.05990 .  
#z.sq_time                 -0.10243    0.09347  -1.096  0.27314    
#z.log_ind_in:z.rank        0.10191    0.05957   1.711  0.08713 .  
#z.rank:z.log_dbh           0.05614    0.06657   0.843  0.39903    
#z.rank:z.fruit_crop_categ -0.08628    0.06102  -1.414  0.15739    

#these are the fixed effects WITH CALC DERIVS to compare:
#Fixed effects:
#                          Estimate Std. Error z value Pr(>|z|)    
#(Intercept)               -1.39504    0.19510  -7.150 8.65e-13 ***
#z.log_ind_in              -0.11649    0.09485  -1.228 0.219420    
#z.rank                    -0.17546    0.07456  -2.353 0.018604 *  
#z.log_dbh                  0.09970    0.13328   0.748 0.454436    
#z.fruit_crop_categ         0.04545    0.12855   0.354 0.723653    
#z.time_elap_end           -0.31411    0.09471  -3.317 0.000911 ***
#z.sq_time_elap_end         0.07827    0.04618   1.695 0.090071 .  
#z.new_time                -0.12369    0.05962  -2.075 0.038013 *  
#z.sq_time                 -0.10243    0.09185  -1.115 0.264763    
#z.log_ind_in:z.rank        0.10191    0.05779   1.764 0.077815 .  
#z.rank:z.log_dbh           0.05614    0.06442   0.871 0.383508    
#z.rank:z.fruit_crop_categ -0.08628    0.05498  -1.569 0.116567   


#resuls are similar enough to continue with the CALC DERIVS model:


#------------------IGNORE results from rerun 28.10. here--------------------------
summary(res)
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                  -1.39675    0.19502  -7.162 7.94e-13 ***
#  z.log_ind_in               -0.11588    0.09483  -1.222 0.221713    
#z.rank                       -0.17570    0.07434  -2.364 0.018101 *  
#  z.log_dbh                  0.10082    0.13358   0.755 0.450421    
#z.fruit_crop_categ           0.04582    0.12873   0.356 0.721891    
#z.time_elap_end              -0.31997    0.09722  -3.291 0.000998 ***
#  z.sq_time_elap_end         0.08236    0.04631   1.778 0.075338 .  
#z.new_time                   -0.12265    0.05965  -2.056 0.039765 *  
#  z.sq_time                 -0.10283    0.09220  -1.115 0.264757    
#z.log_ind_in:z.rank          0.10250    0.05780   1.773 0.076171 .  
#z.rank:z.log_dbh             0.05654    0.06449   0.877 0.380646    
#z.rank:z.fruit_crop_categ  -0.08588    0.05499  -1.562 0.118373    
#-----------------------------------------------------------------------------------------

#-------------------MAIN fixed effects from 1.11. here:--------------------------------------------

summary(res)
#Fixed effects:
#                           Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                -1.39675    0.19502  -7.162 7.94e-13 ***
#  z.log_ind_in             -0.11589    0.09483  -1.222 0.221711    
#z.rank                     -0.17570    0.07434  -2.364 0.018101 *  
#  z.log_dbh                0.10081    0.13358   0.755 0.450426    
#z.fruit_crop_categ         0.04582    0.12873   0.356 0.721885    
#z.time_elap_end            -0.31997    0.09722  -3.291 0.000998 ***
#  z.sq_time_elap_end       0.08236    0.04631   1.778 0.075336 .  
#z.new_time                 -0.12265    0.05965  -2.056 0.039765 *  
#  z.sq_time                -0.10283    0.09220  -1.115 0.264759    
#z.log_ind_in:z.rank        0.10250    0.05780   1.773 0.076171 .  
#z.rank:z.log_dbh           0.05654    0.06449   0.877 0.380649    
#z.rank:z.fruit_crop_categ  -0.08588    0.05499  -1.562 0.118373    

#---------------------------------------------------------------------------------
#-----------MODEL ASSUMPTIONS/STABILITY-----------------------

#check that the model assumptions are fulfilled:
ranef.diagn.plot(res)



library(car)

vif_model=lm(move_y_n~z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.rank + z.time_elap_end + z.new_time, data=xdata) 


vif(vif_model)


#------IGNORE----------------------------------------------------
#  z.log_ind_in          z.log_dbh z.fruit_crop_categ             z.rank    z.time_elap_end 
  #        1.183293           1.129366           1.219009           1.020212           1.128141 
#       z.new_time        1.055642 
#------------------IGNORE from model rerun 28.10.--------------------------------------

#z.log_ind_in          z.log_dbh      z.fruit_crop_categ          z.rank        z.time_elap_end         z.new_time 
#1.183522               1.129396      1.218745                  1.020205           1.128287           1.055664 
#--------------------------------------------------------------------------------------------
#----------------MAIN VIFs 1.11. here:---------------------------------------------------------------------


#z.log_ind_in          z.log_dbh      z.fruit_crop_categ      z.rank    z.time_elap_end         z.new_time 
#1.183522             1.129396           1.218745           1.020205           1.128287           1.055664 

#---------------------------------------------------------------------------------------------------





#next we check the model stability by sourcing Rogers function:
source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\glmm_stability.r")
ls()

head(glmm.model.stab)


model_stability=glmm.model.stab(model.res=res, contr=control_model, n.cores=6, para=T, data=xdata)
table(model_stability$detailed$warnings)
round(model_stability$summary[, -1], 3)
m.stab.plot(model_stability$summary[, -1])



#---------------FULL NULL COMPARISON---------------------



#the null model:
   null_model=glmer(move_y_n~z.new_time + z.sq_time + 

   (1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.time_elap_end + z.sq_time_elap_end + z.new_time + z.sq_time||bonobo_id) + 

    (1 + i_rank_dbh + i_rank_fruit + i_rank_ind + z.log_ind_in + z.rank + z.log_dbh + z.fruit_crop_categ + z.time_elap_end + z.sq_time_elap_end +
    z.new_time + z.sq_time||foc_patch_spp) + 

    (1 + z.log_ind_in + z.time_elap_end + z.sq_time_elap_end + z.new_time + z.sq_time||foc_patch_id) + 

    (1 + z.log_ind_in + z.time_elap_end + z.sq_time_elap_end + z.new_time + z.sq_time||bout_id), 

     data=xdata, family=binomial, control=control_model)



#then the full null comparison significance:
anova(null_model, res, test="Chisq")


#------------IGNORE---------------------------------------------------

#          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#null_model 35 3863.7 4088.1 -1896.9   3793.7                             
#res        40 3852.5 4108.9 -1886.2   3772.5 21.263      5  0.0007225 ***

#-----------------IGNORE--------------------------------

# new p-values for the CALC DERIVS model, comply with the older model:

#          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#null_model 35 3863.7 4088.1 -1896.9   3793.7                             
#res        40 3852.5 4108.9 -1886.2   3772.5 21.263      5  0.0007225 ***

#---------------------------------------------------------------------------
#-----------------IGNORE older comparison 28.10.------------------------

#               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#null_model     35 3863.5 4088.0 -1896.8   3793.5                             
#res            40 3852.1 4108.6 -1886.1   3772.1 21.425      5  0.0006732 ***
#-----------------------------------------------------------------------------------------------------------

#-------------MAIN full null results from rerun 3.11.------------------------------------------------------

#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#null_model         35 3863.5 4088.0 -1896.8   3793.5                             
#res                40 3852.1 4108.6 -1886.1   3772.1 21.425      5  0.0006732 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#--------------------------------------------------------------------------------





#the full null comparison is significant, 
#next get p-values for the interactions:

source("X:\\\\Statistics\\\\R scripts\\\\Roger\\\\drop1_para.r")
ls()
head(drop1p)
drop1.results<-drop1p(model.res=res, data=xdata, contr=control_model)

#-------------IGNORE----------------------------------------------------------------------------------------
#                             AIC     Chisq Chi.Df Pr..Chisq.    n n.warnings warnings
#none                      3852.469        NA     NA         NA 4500          0         
#z.time_elap_end           3856.762 6.2935469      1 0.01211783 4500          0         
#z.sq_time_elap_end        3852.724 2.2553465      1 0.13315365 4500          0         
#z.new_time                3853.882 3.4128746      1 0.06468966 4500          0         
#z.sq_time                 3851.642 1.1726517      1 0.27885706 4500          0         
#z.log_ind_in:z.rank       3853.270 2.8008792      1 0.09421263 4500          0         
#z.rank:z.log_dbh          3851.148 0.6785935      1 0.41007148 4500          0         
#z.rank:z.fruit_crop_categ 3852.537 2.0680371      1 0.15041501 4500          0 
#--------ignore these older CALC DERIVS results---------------------------------------------------------------------------------

#results from new full model with CALC DERIVS, looks similar, so again we drop the 
#interaction terms:
#                               AIC     Chisq Chi.Df Pr..Chisq.    n n.warnings warnings
#none                      3852.469        NA     NA         NA 4500          0         
#z.time_elap_end           3856.762 6.2935469      1 0.01211783 4500          0         
#z.sq_time_elap_end        3852.724 2.2553465      1 0.13315365 4500          0         
#z.new_time                3853.882 3.4128746      1 0.06468966 4500          0         
#z.sq_time                 3851.642 1.1726517      1 0.27885706 4500          0         
#z.log_ind_in:z.rank       3853.270 2.8008792      1 0.09421263 4500          0         
#z.rank:z.log_dbh          3851.148 0.6785935      1 0.41007148 4500          0         
#z.rank:z.fruit_crop_categ 3852.537 2.0680371      1 0.15041501 4500          0  
#----------------------------------------------------------------
#---------------IGNORE these 28.10. rerun------------------------------------------
#                           logLik      AIC     Chisq Chi.Df Pr..Chisq.    n n.opt.warnings n.fun.warnings
#none                      -1886.057 3852.115        NA     NA         NA 4500             NA             NA
#z.time_elap_end           -1889.181 3856.362 6.2475418      1 0.01243658 4500              0              0
#z.sq_time_elap_end        -1887.231 3852.461 2.3462648      1 0.12558359 4500              0              0
#z.new_time                -1887.723 3853.447 3.3315826      1 0.06796145 4500              0              0
#z.sq_time                 -1886.646 3851.291 1.1763548      1 0.27809934 4500              0              0
#z.log_ind_in:z.rank       -1887.470 3852.940 2.8247658      1 0.09282043 4500              0              0
#z.rank:z.log_dbh          -1886.401 3850.803 0.6876424      1 0.40696744 4500              0              0
#z.rank:z.fruit_crop_categ -1887.087 3852.174 2.0592692      1 0.15128271 4500              0              0
#-----------------------------------------------------------------------------------------------------------------------

#---------------------------------MAIN results from the 3.11. rerun-----------------------------------------


#                           logLik      AIC     Chisq Chi.Df Pr..Chisq.    n n.opt.warnings n.fun.warnings
#none                      -1886.057 3852.115        NA     NA         NA 4500             NA             NA
#z.time_elap_end           -1889.181 3856.362 6.2475420      1 0.01243658 4500              0              0
#z.sq_time_elap_end        -1887.231 3852.461 2.3462650      1 0.12558358 4500              0              0
#z.new_time                -1887.723 3853.447 3.3316388      1 0.06795913 4500              0              0
#z.sq_time                 -1886.646 3851.291 1.1763556      1 0.27809918 4500              0              0
#z.log_ind_in:z.rank       -1887.470 3852.940 2.8247640      1 0.09282053 4500              0              0
#z.rank:z.log_dbh          -1886.401 3850.803 0.6876429      1 0.40696728 4500              0              0
#z.rank:z.fruit_crop_categ -1887.087 3852.174 2.0592688      1 0.15128275 4500              0              0
#
#$model.results
#NULL

#-----------------------------------------------------------------------------------------------------


#based on the rerun 3.11., interactions still not significant, hence we run a reduced model WITHOUT interactions:


#new reduced Model III: res WITHOUT interactions and without the SQUARE TERM time elapsed, rerun 4.11.:

start.time<-Sys.time()
res_without_inter=glmer(move_y_n~z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.rank + z.time_elap_end + z.new_time + z.sq_time +
           
           (1 + z.log_ind_in + z.log_dbh + z.fruit_crop_categ + z.time_elap_end + z.new_time + z.sq_time||bonobo_id) + 
           
           (1 + z.log_ind_in + z.rank + z.log_dbh + z.fruit_crop_categ + z.time_elap_end + z.new_time + z.sq_time||foc_patch_spp) + 
           
           (1 + z.log_ind_in + z.time_elap_end + z.new_time + z.sq_time||foc_patch_id) + 
           
           (1 + z.log_ind_in + z.time_elap_end + z.new_time + z.sq_time||bout_id),
         
         data=xdata, family=binomial, control=control_model)

end.time<-Sys.time()
end.time - start.time




#then we need the p-values for the new model:
drop1.results.final<-drop1p(model.res=res_without_inter, data=xdata, contr=control_model)
drop1.results.final


#--------------------IGNORE-----------------------------------------------------------------------
#                       AIC     Chisq Chi.Df Pr..Chisq.    n n.warnings warnings
#none               3848.304        NA     NA         NA 4500          0         
#z.log_ind_in       3847.073 0.7690227      1 0.38051954 4500          0         
#z.log_dbh          3846.710 0.4055190      1 0.52425271 4500          0         
#z.fruit_crop_categ 3846.423 0.1188181      1 0.73031992 4500          0         
#z.rank             3850.192 3.8876968      1 0.04864109 4500          0         
#z.time_elap_end    3852.906 6.6016936      1 0.01018818 4500          0         
#z.new_time         3850.062 3.7573578      1 0.05257560 4500          0         
#z.sq_time          3847.448 1.1436746      1 0.28487720 4500          0                                                                                                           #                                

#---------------IGNORE reduced model 31.10.--------------------------------------------------------

#                   logLik      AIC       Chisq   Chi.Df  Pr..Chisq.  n   n.opt.warnings n.fun.warnings
#none               -1891.034 3848.068        NA     NA          NA 4500             NA             NA
#z.log_ind_in       -1891.415 3846.830 0.7617287      1 0.382788067 4500              0              0
#z.log_dbh          -1891.239 3846.477 0.4091172      1 0.522417904 4500              0              0
#z.fruit_crop_categ -1891.094 3846.188 0.1202233      1 0.728792468 4500              0              0
#z.rank             -1892.980 3849.959 3.8911119      1 0.048542275 4500              0              0
#z.time_elap_end    -1894.357 3852.715 6.6463801      1 0.009935748 4500              0              0
#z.new_time         -1892.921 3849.843 3.7746689      1 0.052034197 4500              0              0
#z.sq_time          -1891.603 3847.206 1.1381822      1 0.286036764 4500              0              0
#-----------------------------------------------------------------------------------------------------

#--------------MAIN new P-values for 4.11. rerun-------------------------------------------------------
#                   logLik      AIC     Chisq Chi.Df  Pr..Chisq.    n n.opt.warnings n.fun.warnings
#none               -1891.034 3848.068        NA     NA          NA 4500             NA             NA
#z.log_ind_in       -1891.415 3846.830 0.7617287      1 0.382788060 4500              0              0
#z.log_dbh          -1891.239 3846.477 0.4091171      1 0.522417955 4500              0              0
#z.fruit_crop_categ -1891.094 3846.188 0.1202230      1 0.728792742 4500              0              0
#z.rank             -1892.980 3849.959 3.8911115      1 0.048542289 4500              0              0
#z.time_elap_end    -1894.357 3852.715 6.6463799      1 0.009935749 4500              0              0
#z.new_time         -1892.921 3849.843 3.7746687      1 0.052034204 4500              0              0
#z.sq_time          -1891.603 3847.206 1.1381821      1 0.286036787 4500              0              0


#--------------------------------------------------------------------------------------------------------







#we see that rank and time until end of bout are significant, next, get the direction
#of effect:

summary(res_without_inter)


#-------------------IGNORE----------------------------------------------------
#Fixed effects:
#                   Estimate Std. Error z value Pr(>|z|)    
#(Intercept)        -1.34778    0.19421  -6.940 3.93e-12 ***
#z.log_ind_in       -0.10213    0.10401  -0.982  0.32610    
#z.log_dbh           0.08788    0.13254   0.663  0.50730    
#z.fruit_crop_categ  0.04674    0.12835   0.364  0.71575    
#z.rank             -0.15043    0.06869  -2.190  0.02851 *  
#z.time_elap_end    -0.25759    0.08720  -2.954  0.00314 ** 
#z.new_time         -0.12778    0.05952  -2.147  0.03179 *  
#z.sq_time          -0.09185    0.08238  -1.115  0.26485    
#----------------------------------------------------------------------------------------
#----------------IGNORE rerun 31.10.-----------------------------------------
#Fixed effects:
#                   Estimate  Std. Error z value Pr(>|z|)    
#(Intercept)        -1.34889    0.19413  -6.948  3.7e-12 ***
#z.log_ind_in       -0.10165    0.10406  -0.977  0.32866    
#z.log_dbh           0.08834    0.13256   0.666  0.50513    
#z.fruit_crop_categ  0.04701    0.12834   0.366  0.71416    
#z.rank             -0.15046    0.06863  -2.192  0.02836 *  
#z.time_elap_end    -0.26059    0.08755  -2.976  0.00292 ** 
#z.new_time         -0.12813    0.05953  -2.152  0.03137 *  
#z.sq_time          -0.09180    0.08256  -1.112  0.26616    


#------------------MAIN direction of the effect 4.11.---------------------------------
#Fixed effects:
#                   Estimate Std. Error z value Pr(>|z|)    
#(Intercept)        -1.34889    0.19413  -6.948  3.7e-12 ***
#z.log_ind_in       -0.10165    0.10406  -0.977  0.32866    
#z.log_dbh           0.08835    0.13256   0.666  0.50513    
#z.fruit_crop_categ  0.04701    0.12834   0.366  0.71417    
#z.rank             -0.15046    0.06863  -2.192  0.02836 *  
#z.time_elap_end    -0.26059    0.08755  -2.976  0.00292 ** 
#z.new_time         -0.12813    0.05953  -2.152  0.03138 *  
#z.sq_time          -0.09180    0.08256  -1.112  0.26617    
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#-----------------------------------------------------------------------------------------


#------------RESULT SUMMARY-------------------------------------------------------------------------


#rank is significant (p=0.0485), the direction of the relationship is negative (-0.15)
#with movement meaning that the higher the rank, the lower the movement, tends towards zero 
#time until end of bout is signif with also a negative relationship (-0.26)
#hence the closer to end of bout we are the more there is movement


#-------------PLOT RESULTS---------------------------------------------------------------------------------

#plotting the results:
plot(x=xdata$z.time_elap_end, y=xdata$move_y_n, pch=19, col=adjustcolor("grey19", alpha=0.1), xlab="Time until end of bout", ylab="Movement", xaxt="n", las=1)
axis(side=1, at=(c()-mean(xdata$time_elap_end))/sd(xdata$time_elap_end), labels=c())


xdata$time_elap_end
unique(xdata$time_elap_end)

#too many unique, we need to bin values
range(xdata$z.time_elap_end)

binwidth=1000 

xdata$time_bins=xdata$time_elap_end %/% binwidth

xdata$time_bins=xdata$time_bins*binwidth

xdata$time_bins=xdata$time_bins+binwidth/2


unique(xdata$time_bins)

#looks pretty good
#now need to find proportion of what is the response per bin


proportion_bin=aggregate(x=xdata$move_y_n, by=list(xdata$time_bins), FUN=mean)  

head(proportion_bin)

#now turn the bin values to z space equalents:
proportion_bin$z_bin=(proportion_bin$Group.1-mean(xdata$time_elap_end))/sd(xdata$time_elap_end)


#for the size of data points create the object with number of obs per bin:
data_point_size=aggregate(xdata$move_y_n, by=list(xdata$time_bins), FUN=length)



plot(x=proportion_bin$z_bin, y=proportion_bin$x, pch=19, cex=data_point_size$x^(1/6), col="grey19", xlab="Minutes until end of bout", ylab="Probability of movement", xaxt="n", las=1)
axis(side=1, at=(c(0, 3600, 7200, 10800, 14400)-mean(xdata$time_elap_end))/sd(xdata$time_elap_end), labels=c(0, 3600/60, 7200/60, 10800/60, 14400/60))

min(xdata$time_elap_end)
max(xdata$time_elap_end)


#next put model line in:
library(lme4)

plot.coeffs<-fixef(res_without_inter)

plot.xvals<-seq(min(xdata$z.time_elap_end), max(xdata$z.time_elap_end), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.time_elap_end"]*plot.xvals
plot.yvals<-exp(plot.yvals)/(1+exp(plot.yvals))

lines(plot.xvals, plot.yvals, lwd=3, lty=2)


#lets instead bin original time elap values then transform in Z space so was redone above




##################################
#-----------PLOT Move vs Rank next:----------------------------------------


#----create data points with means per rank----------------------------------------------------
datapoint_size_rank=aggregate(xdata$move_y_n, by=list(xdata$bonobo_id, xdata$rank), FUN=length)
#lets look at it:
datapoint_size_rank

str(datapoint_size_rank)

datapoint_size_rankmean=aggregate(xdata$move_y_n, by=list(xdata$bonobo_id, xdata$rank), FUN=mean)
datapoint_size_rankmean


#now turn the bin values to z space equalents:
datapoint_size_rankmean$z_rank=(datapoint_size_rankmean$Group.2-mean(xdata$rank))/sd(xdata$rank)

#--------------------------------------------------------------------------------------------------



#plotting the results:
plot(x=datapoint_size_rankmean$z_rank, y=datapoint_size_rankmean$x, pch=19, cex=datapoint_size_rank$x^(1/4), col=adjustcolor("black", alpha=0.6), xlab="Female rank", ylab="Probability of movement", xaxt="n", las=1)
axis(side=1, at=(c(0, 0.5, 1)-mean(xdata$rank))/sd(xdata$rank), labels=c(0, 0.5, 1))





#then the MODEL LINE for Move vs Rank---------------------------------------------------------------

library(lme4)

plot.coeffs<-fixef(res_without_inter)

plot.xvals<-seq(min(xdata$z.rank), max(xdata$z.rank), length.out=100)

plot.yvals<-plot.coeffs["(Intercept)"]+plot.coeffs["z.rank"]*plot.xvals
plot.yvals<-exp(plot.yvals)/(1+exp(plot.yvals))

lines(plot.xvals, plot.yvals, lwd=3, lty=2)






q()








save.image("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_III_movement\\\\MAIN_pan_analysis_III_movement_workspace_07112016.Rdata")
load("Y:\\\\primint\\\\Niina\\\\Files_for_R_analysis_III_movement\\\\MAIN_pan_analysis_III_movement_workspace_07112016.Rdata")









