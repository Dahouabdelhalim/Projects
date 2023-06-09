
#SEM derives from the CE and SE on final stability of grain yields.
library(sem)
library(dplyr)
library(DiagrammeR)
options (digits=1,scipen=300)
library(xlsx)
data1 <-"***/Wu et al., Data for R.xlsx"
dat<- read.xlsx(data1,1,header = TRUE)
dat
table(TSPB)
scale(dat,center=T,scale=F)
scale(dat,center=T,scale=T)
cor_num <- cor(dat) 
opt <- options(fit.indices = c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"))
model.kerch <- specifyModel(
  text = '
  
  CE -> TSPB , CE _ TSPB, NA
  SE -> TSPB , SE _ TSPB, NA
  CE -> TSPG , CE _ TSPG, NA
  SE -> TSPG , SE _ TSPG, NA
  TSPB -> TSPG , TSPB _ TSPG, NA
  CE -> SUBTSPB , CE _ SUBTSPB, NA
  SE -> SUBTSPB, SE _ SUBTSPB, NA
  CE -> DOMTSPB , CE _ DOMTSPB, NA
  SE -> DOMTSPB, SE _ DOMTSPB, NA
  DOMTSPB ->SUBTSPB,   DOMTSPB _  SUBTSPB, NA

  SUBTSPB -> TSPB , SUBTSPB _TSPB, NA
  DOMTSPB -> TSPB, DOMTSPB _TSPB, NA
  
  SE <-> SE  , SE _SE, NA
  CE <-> CE  , CE _CE, NA
  DOMTSPB <->DOMTSPB,   DOMTSPB _  DOMTSPB, NA
  SUBTSPB <-> SUBTSPB , SUBTSPB _ SUBTSPB, NA
  TSPB <->TSPB, TSPB _ TSPB, NA
  TSPG <->TSPG, TSPG _ TSPG, NA
  
  ')

dat[!is.na(dat)]
out_sem <- sem(model.kerch, cor_num, nrow(dat))
coef <- out_sem$coeff
coeff_name <- out_sem$semmod[,1]
standardizedCoefficients(out_sem)
summary(out_sem)
semPaths(out_sem, "Standardized", "Estimates")
pathDiagram(out_sem, edge.labels="values",edge.weight ="proportional", standardized =T,rank.direction="TB") 



#Final stable SEM derives from the CE and SE on final stability of grain yields.
library(sem)
library(dplyr)
library(DiagrammeR)
options (digits=1,scipen=300)
library(xlsx)
data1 <-"***/Wu et al., Data for R.xlsx"
dat<- read.xlsx(data1,1,header = TRUE)
attach(dat)
table(TSPB)
scale(dat,center=T,scale=F)
scale(dat,center=T,scale=T)
cor_num <- cor(dat) 
opt <- options(fit.indices = c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"))
model.kerch <- specifyModel(
  text = '
  
  #CE -> TSPB , CE _ TSPB, NA
  #SE -> TSPB , SE _ TSPB, NA
  #CE -> TSPG , CE _ TSPG, NA
  #SE -> TSPG , SE _ TSPG, NA
  TSPB -> TSPG , TSPB _ TSPG, NA
  CE -> SUBTSPB , CE _ SUBTSPB, NA
  #SE -> SUBTSPB, SE _ SUBTSPB, NA
  #CE -> DOMTSPB , CE _ DOMTSPB, NA
  #SE -> DOMTSPB, SE _ DOMTSPB, NA
  
  SUBTSPB -> TSPB , SUBTSPB _TSPB, NA
  #DOMTSPB -> TSPB, DOMTSPB _TSPB, NA
  
  #SE <-> SE  , SE _SE, NA
  CE <-> CE  , CE _CE, NA
  #DOMTSPB <->DOMTSPB,   DOMTSPB _  DOMTSPB, NA
  SUBTSPB <-> SUBTSPB , SUBTSPB _ SUBTSPB, NA
  TSPB <->TSPB, TSPB _ TSPB, NA
  TSPG <->TSPG, TSPG _ TSPG, NA
  
  ')

dat[!is.na(dat)]
out_sem <- sem(model.kerch, cor_num, nrow(dat))
coef <- out_sem$coeff
coeff_name <- out_sem$semmod[,1]
standardizedCoefficients(out_sem)
summary(out_sem)
semPaths(out_sem, "Standardized", "Estimates")
pathDiagram(out_sem, edge.labels="values",edge.weight ="proportional", standardized =T,rank.direction="TB")




#SEM derives from the RII of competitively dominant/subordination species on final stability of grain yields.
library(dplyr)
library(DiagrammeR)
options (digits=1,scipen=300)
library(xlsx)
data1 <-"***/Wu et al., Data for R.xlsx"
dat<- read.xlsx(data1,2,header = TRUE)
dat
cor_num <- cor(dat) 
opt <- options(fit.indices = c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"))
model.kerch <- specifyModel(
  text = '

  #SRII->TSPB, CE_TSPG, NA
  #SRII->TSPG, SRII_TSPG,NA
  #SRII->DOMTSPB,SRII_DOMTSPB,NA
  #SRII->SUBTSPB,SRII_SUBTSPB,NA
  
  #DOMRII->SRII, DOMRII_SRII, NA
  DOMRII->TSPB, DOMRII_TSPG, NA
  DOMRII->TSPG, DOMRII_TSPG,NA
  DOMRII->DOMTSPB,DOMRII_DOMTSPB,NA
  DOMRII->SUBTSPB,DOMRII_SUBTSPB,NA
  
  #SUBRII->SRII, SUBRII_SRII, NA
  SUBRII->TSPB, SUBRII_TSPG, NA
  SUBRII->TSPG, SUBRII_TSPG,NA
  SUBRII->DOMTSPB,SUBRII_DOMTSPB,NA
  SUBRII->SUBTSPB,SUBRII_SUBTSPB,NA
  
  #CE -> DOMRII , CE _ DOMRII, NA
  #SE -> DOMRII , SE _ DOMRII, NA
  #CE -> SUBRII , CE _ SUBRII, NA
  #SE -> SUBRII , SE _ SUBRII, NA
 
  #CE -> TSPB , CE _ TSPB, NA
  #SE -> TSPB , SE _ TSPB, NA
  #CE -> TSPG , CE _ TSPG, NA
  #SE -> TSPG , SE _ TSPG, NA
  TSPB -> TSPG , TSPB _ TSPG, NA
  #CE -> SUBTSPB , CE _ SUBTSPB, NA
  #SE -> SUBTSPB, SE _ SUBTSPB, NA
  #CE -> DOMTSPB , CE _ DOMTSPB, NA
  #SE -> DOMTSPB, SE _ DOMTSPB, NA
  DOMTSPB ->SUBTSPB,   DOMTSPB _  SUBTSPB, NA

  SUBTSPB -> TSPB , SUBTSPB _TSPB, NA
  DOMTSPB -> TSPB, DOMTSPB _TSPB, NA
  #SUBTSPB -> TSPG , SUBTSPB _TSPG, NA
  #DOMTSPB -> TSPG, DOMTSPB _TSPG, NA
  
  SUBRII<-> SUBRII, SUBRII_ SUBRII, NA
  DOMRII<-> DOMRII, DOMRII_ DOMRII, NA
  #SRII<->SRII,SRII_SRII,NA
  #SE <-> SE  , SE _SE, NA
  #CE <-> CE  , CE _CE, NA
  DOMTSPB <->DOMTSPB,   DOMTSPB _  DOMTSPB, NA
  SUBTSPB <-> SUBTSPB , SUBTSPB _ SUBTSPB, NA
  #TSPB <->TSPB, TSPB _ TSPB, NA
  #TSPG <->TSPG, TSPG _ TSPG, NA
  
  ')

out_sem <- sem(model.kerch, cor_num, nrow(dat))
coef <- out_sem$coeff
coeff_name <- out_sem$semmod[,1]
standardizedCoefficients(out_sem)
summary(out_sem)
semPaths(out_sem, "Standardized", "Estimates")
pathDiagram(out_sem, edge.labels="values",edge.weight ="proportional", standardized =T,rank.direction="TB") 


#Final stable SEM derives from the RII of competitively dominant/subordination species on final stability of grain yields.
library(sem)
library(dplyr)
library(DiagrammeR)
options (digits=1,scipen=300)
library(xlsx)
data1 <-"***/Wu et al., Data for R.xlsx"
dat<- read.xlsx(data1,2,header = TRUE)
dat
cor_num <- cor(dat) 
opt <- options(fit.indices = c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"))
model.kerch <- specifyModel(
  text = '
  
  #SRII->TSPB, CE_TSPG, NA
  #SRII->TSPG, SRII_TSPG,NA
  #SRII->DOMTSPB,SRII_DOMTSPB,NA
  #SRII->SUBTSPB,SRII_SUBTSPB,NA
  
  #DOMRII->SRII, DOMRII_SRII, NA
  #DOMRII->TSPB, DOMRII_TSPG, NA
  #DOMRII->TSPG, DOMRII_TSPG,NA
  #DOMRII->DOMTSPB,DOMRII_DOMTSPB,NA
  #DOMRII->SUBTSPB,DOMRII_SUBTSPB,NA
  
  
  #SUBRII->SRII, SUBRII_SRII, NA
  #SUBRII->TSPB, SUBRII_TSPG, NA
  #SUBRII->TSPG, SUBRII_TSPG,NA
  #SUBRII->DOMTSPB,SUBRII_DOMTSPB,NA
  SUBRII->SUBTSPB,SUBRII_SUBTSPB,NA
  
  #CE -> DOMRII , CE _ DOMRII, NA
  #SE -> DOMRII , SE _ DOMRII, NA
  #CE -> SUBRII , CE _ SUBRII, NA
  #SE -> SUBRII , SE _ SUBRII, NA
 
  
  #CE -> TSPB , CE _ TSPB, NA
  #SE -> TSPB , SE _ TSPB, NA
  #CE -> TSPG , CE _ TSPG, NA
  #SE -> TSPG , SE _ TSPG, NA
  TSPB -> TSPG , TSPB _ TSPG, NA
  #CE -> SUBTSPB , CE _ SUBTSPB, NA
  #SE -> SUBTSPB, SE _ SUBTSPB, NA
  #CE -> DOMTSPB , CE _ DOMTSPB, NA
  #SE -> DOMTSPB, SE _ DOMTSPB, NA
  #DOMTSPB ->SUBTSPB,   DOMTSPB _  SUBTSPB, NA


  SUBTSPB -> TSPB , SUBTSPB _TSPB, NA
  #DOMTSPB -> TSPB, DOMTSPB _TSPB, NA
  
  #SUBTSPB -> TSPG , SUBTSPB _TSPG, NA
  #DOMTSPB -> TSPG, DOMTSPB _TSPG, NA
  
  
  SUBRII<-> SUBRII, SUBRII_ SUBRII, NA
  #DOMRII<-> DOMRII, DOMRII_ DOMRII, NA
  #SRII<->SRII,SRII_SRII,NA
  #SE <-> SE  , SE _SE, NA
  #CE <-> CE  , CE _CE, NA
  #DOMTSPB <->DOMTSPB,   DOMTSPB _  DOMTSPB, NA
  SUBTSPB <-> SUBTSPB , SUBTSPB _ SUBTSPB, NA
  #TSPB <->TSPB, TSPB _ TSPB, NA
  #TSPG <->TSPG, TSPG _ TSPG, NA
  
  ')

out_sem <- sem(model.kerch, cor_num, nrow(dat))
coef <- out_sem$coeff
coeff_name <- out_sem$semmod[,1]
standardizedCoefficients(out_sem)
summary(out_sem)
semPaths(out_sem, "Standardized", "Estimates")
pathDiagram(out_sem, edge.labels="values",edge.weight ="proportional", standardized =T,rank.direction="TB") 




# The SEM (resource using as mediating variables) derives from the RII of competitively dominant/subordination species on final stability of grain yields.
library(sem)
library(dplyr)
library(DiagrammeR)
options (digits=1,scipen=300)
library(xlsx)
data1 <-"***/Wu et al., Data for R.xlsx"
dat<- read.xlsx(data1,3,header = TRUE)
table(TSPB)
table(W)
table(N)
table(DOMTSPB)
scale(dat,center=T,scale=F)
scale(dat,center=T,scale=T)
cor_num <- cor(dat) 
opt <- options(fit.indices = c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"))
model.kerch <- specifyModel(
  text = '
  
  DOMRII -> TSPB , DOMRII _ TSPB, NA
  W -> TSPB , W _ TSPB, NA 
  N -> TSPB, N _ TSPB, NA  
  P -> TSPB  , P _ TSPB, NA        					
  K -> TSPB  , K _ TSPB, NA 
  SUBTSPB -> TSPB , SUBTSPB _TSPB, NA
  DOMTSPB -> TSPB, DOMTSPB _TSPB, NA
  DOMTSPB -> SUBTSPB, DOMTSPB _SUBTSPB, NA
   
  W -> SUBTSPB , W _ SUBTSPB, NA 
  N -> SUBTSPB, N _ SUBTSPB, NA  
  P -> SUBTSPB  , P _ SUBTSPB, NA        					
  K -> SUBTSPB  , K _ SUBTSPB, NA 

  W -> DOMTSPB , W _ DOMTSPB, NA 
  N -> DOMTSPB, N _ DOMTSPB, NA  
  P -> DOMTSPB , P _ DOMTSPB, NA 
  K -> DOMTSPB , K _ DOMTSPB, NA 
  
  SUBRII ->W , SUBRII_W  , NA 
  SUBRII->N , SUBRII_N  , NA  
  SUBRII->P , SUBRII_P  , NA        					
  SUBRII->K , SUBRII_K  , NA
  #SUBRII ->SUBTSPB , SUBRII_SUBTSPB  , NA 
  #SUBRII ->SUBTSPB , SUBRII_SUBTSPB  , NA 
  SUBRII ->TSPB , SUBRII_TSPB  , NA 
  
  DOMRII ->W , DOMRII_W  , NA 
  DOMRII->N , DOMRII_N  , NA  
  DOMRII->P , DOMRII_P  , NA        					
  DOMRII->K , DOMRII_K  , NA
  #DOMRII ->DOMTSPB , DOMRII_DOMTSPB  , NA 
  #DOMRII ->SUBTSPB , DOMRII_SUBTSPB  , NA 
  DOMRII ->TSPB , DOMRII_TSPB  , NA 

  W ->  N,  W _  N , NA	
  W -> P, W _ P, NA  
  W -> K,  W _ K , NA
  
  W <-> W,  W _  W , NA	
  N <-> N, N _ N, NA  
  P <-> P,  P _ P , NA
  K <-> K,  K _ K , NA
  
  SUBRII <-> SUBRII  , SUBRII_ SUBRII, NA
  DOMRII <-> DOMRII  , DOMRII_ DOMRII, NA
  DOMTSPB <->DOMTSPB,   DOMTSPB _  DOMTSPB, NA
  SUBTSPB <-> SUBTSPB , SUBTSPB _ SUBTSPB, NA
  TSPB <->TSPB, TSPB _ TSPB, NA
  
  ')

dat[!is.na(dat)]
out_sem <- sem(model.kerch, cor_num, nrow(dat))
coef <- out_sem$coeff
coeff_name <- out_sem$semmod[,1]
standardizedCoefficients(out_sem)
summary(out_sem)
semPaths(out_sem, "Standardized", "Estimates")
pathDiagram(out_sem, edge.labels="values",edge.weight ="proportional", standardized =T,rank.direction="TB")


