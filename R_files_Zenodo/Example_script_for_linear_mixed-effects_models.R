library(usdm)
library(MuMIn) 
library(nlme)
library(lme4)
library(lmerTest)
library(arm)
library(reshape)
library(MASS)
library(vegan)
library(moments)

ind <- read.csv("All_region_ind.csv",sep=",", fileEncoding="UTF-8-BOM", row.names="Site",header=TRUE) 
#"ind" lists all the biological response metrics
env <- read.csv("All_region_ENV.csv",sep=",",row.names="Unique_sample_ID",header=TRUE)
# "env" lists all the environmental (i.e. intermittence, impacts etc.) variables

envind<-merge(env,ind, by ="row.names")
write.table(envind,"envind.xls",sep="\\t",dec=".")

# Are there missing values?
colSums(is.na(envind))
#no missing values

names(envind)
#find the col numbers to use below...
envind[,c(11:15)] = apply(envind[,c(11:15)], 2, function(x) as.numeric(as.character(x)));

allvar<-c("Latitude","Longitude","Mean_Temp","Mean_Precip","Mean_Impact","MeanNoFlow","PropNoFlow","TimeSince", "Aridity")
#here we select the variables we are interested in - everything (above - just to see) and a subset (below - what we use)
allvar2<-c("Mean_Impact","PropNoFlow","Aridity")

envindcor <- envind[,allvar] 
envindcor2 <- envind[,allvar2] 

######remove NAs
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

envindS<-envind
#create a copy of envind
envindS<-completeFun(envindS, allvarz2)
#modify by removing all rows which have NA in any of the three col
envindS[,allvar2]<-scale(envindS[,allvar2])
#the above scales the variables (centred around zero and divided by SD - standardizing the to vary within the same range)

vif(envindcor2)

FamSkew<-c(envind$fam_rich)
skewness(FamSkew)
FRedSkew<-c(envind$ffred)
skewness(FRedSkew)
FRichSkew<-c(envind$ffrich)
skewness(FRichSkew)
UK_BMWPSkew<-c(envind$BMWP)
skewness(UK_BMWPSkew)
UK_ASPTSkew<-c(envind$ASPT)
skewness(UK_ASPTSkew)
IBMWPSkew<-c(envind$IBMWP)
skewness(IBMWPSkew)
IASPTSkew<-c(envind$IASPT)
skewness(IASPTSkew)
WHPTSkew<-c(envind$WHPT_Total)
skewness(WHPTSkew)
WHPT_ASPTSkew<-c(envind$WHPT_ASPT)
skewness(WHPT_ASPTSkew)

options(na.action = "na.omit") # necessary to run dredge()

#################################
# multi-model inference    # 
#################################

envindDF<-as.data.frame(envindS)

names(envindDF)

#transform to improve the normality e.g.
sqrtIBMWP<-sqrt(envindDF$IBMWP)
envindDF["sqrtIBMWP"]<-sqrtIBMWP
sqrtIBMWPSkew<-c(envindDF$sqrtIBMWP)
skewness(sqrtIBMWPSkew)

###open "0_model_functions" script from Soria et al. https://doi.org/10.1111/1365-2664.13538 and run it in a different window

#below we are building a model...

mult_res<-multi.lm(envindDF, sel.resp=envindDF[,c(24:length(envindDF))], 
                   sel.fixed=c("Aridity" , "Mean_Impact" , "PropNoFlow"  , "Aridity:Mean_Impact", "Mean_Impact:PropNoFlow", "Aridity:PropNoFlow"),
                   sel.random="(1|Dataprov/Site.code)",
                   delta.set="d2",
                   #above builds modesl and averages best models - delta is to select within delta (i.e. AIC) 2 - <2 from the AIC.
                   av.method="subset",var.pois=c(F,rep(F,ncol(envindDF[,c(24:length(envindDF))]))),
                   plot_assum=T)

write.table(mult_res,"All_region_WC.xls",sep="\\t",dec=".")

#variance partitioning
library(variancePartition)

var<-c(24:length(envindDF))
sel.fixed=c("Aridity","Mean_Impact" , "PropNoFlow"  , "PropNoFlow:Mean_Impact","Mean_Impact:Aridity","PropNoFlow:Aridity")
matvarpart<-matrix(rep(NA,(length(var)*(length(sel.fixed)+3))),nrow=length(var),ncol=(length(sel.fixed)+3),dimnames=(list(colnames(envindDF[,c(24:length(envindDF))]),c("Dataprov","Site.code:Dataprov",sel.fixed,"Residuals"))))

for(v in 1:length(c(24:length(envindDF)))){
  try({V<-var[v]
  resp<-envindDF[,V]
  mod<-lmer(resp~Aridity + Mean_Impact + PropNoFlow + PropNoFlow:Mean_Impact + Mean_Impact:Aridity+ PropNoFlow:Aridity + (1|Dataprov/Site.code), REML=F, data=envindDF)
  varpart<-calcVarPart(mod,showWarnings=FALSE)
  matvarpart[v,]<-varpart
  }, silent=TRUE)
}

matvarpart
write.table(matvarpart,"results.xls",sep="\\t",dec=".")

