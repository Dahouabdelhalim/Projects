library(nlme)
library(lme4)
library(car)
library(effects)
library(lsmeans)
library(pbkrtest)
library(effects)
#install.packages("ggplot2")
library("r2glmm")
library(ggplot2)
library(RRPP)
library(plyr)
spectra=read.csv('7_27_20_new_spectra_smooth_combined_PCA_ready_LR_averaged.csv')
colspace<-read.csv('Childers_Jalmenus_perceptual_tetrahedral_colspace_and_stimulation_values_RARC_curated.csv')
spectraSt=subset(spectra, Illumination=="Standard")

spectraFWSt=subset(spectraSt, Wing=="FW")
spectraHWSt=subset(spectraSt, Wing=="HW")

colFWSt=subset(colspace, Wing=="FW")
colHWSt=subset(colspace, Wing=="HW")




RRPPaov=function(input_spectra,inputcolspace){
  y=data.matrix(input_spectra[1:26,7:457])
  y1=scale(y)
  colscores=data.matrix(inputcolspace[1:26,12:14])
  sex=input_spectra[1:26,4]
  rrppframe=rrpp.data.frame(y=y1,
                            colscores=colscores,
                            'Sex'=sex)
  rrppframe
  fit <- lm.rrpp(y ~ Sex,
                 SS.type = "I", data = rrppframe,
                 print.progress = FALSE,iter=10000)
  a=anova(fit, effect.type = "F")
  rawlist=c(a$table$Df[1],a$table$SS[1],a$table$Rsq[1],a$table$F[1],a$table$Z[1],a$table$`Pr(>F)`[1])
  
  fit2 <- lm.rrpp(colscores ~ Sex,
                  SS.type = "I", data = rrppframe,
                  print.progress = FALSE,iter=10000)
  
  c=anova(fit2, effect.type = "F")
  collist=c(c$table$Df[1],c$table$SS[1],c$table$Rsq[1],c$table$F[1],c$table$Z[1],c$table$`Pr(>F)`[1])
  
  returnframe=data.frame('Statistic'=c('Df','SS','Rsq','F','Z (effect size)','p-value'),
                         'Raw_Results'=rawlist,'Tetracolor_results'=collist)
  returnframe
  return(returnframe)
  
}


FWSt=RRPPaov(spectraFWSt,colFWSt)
HWSt=RRPPaov(spectraHWSt,colHWSt)


resultsframe=data.frame('FW_standard_Raw'=FWSt[2],'FW_standard_tetracolor'=FWSt[3],
                        'HW_standard_Raw'=HWSt[2],'HW_standard_tetracolor'=HWSt[3])


names(resultsframe)=c('FW_standard_Raw','FW_standard_tetracolor',
                      'HW_standard_Raw','HW_standard_tetracolor')

resultsframe1=t(resultsframe)
resultsframe1=data.frame(Treatment=c('FW_standard_Raw','FW_standard_tetracolor',
                                     'HW_standard_Raw','HW_standard_tetracolor')
                         ,'Variable'=c('Sex','Sex','Sex','Sex')
                         ,resultsframe1)
names(resultsframe1)=c('Treatment','Variable','Df','SS','Rsq','F','Z (Effect size)','p_value')
resultsframe1


resultsframe1['P_adjusted_BH']=p.adjust(resultsframe1$p_value,method='BH')
resultsframe1
# write.csv(resultsframe1,'RRPP_anova_results_for_rawdata_noPC_scores_tetrahedral_colspace_lm_standard_only.csv')
