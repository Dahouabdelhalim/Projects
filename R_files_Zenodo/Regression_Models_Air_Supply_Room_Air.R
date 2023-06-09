#module load gcc/6.3
#module load R/4.0.0
library(data.table)



setwd("/lustre/ahome3/m/mr251/EMBER/Reintegrated_Data_Nov2020/RA_Air_Supply_Feb2021")
dd_n277<-read.table(file="AS_peak table based on 277 template for Matt with 9 reintegrated features.csv",header=TRUE,sep=",")
dd_n277<-read.table(file="RA_peak table based on 277 template for Matt with 9 reintegrated features.csv",header=TRUE,sep=",")



dd_gcgc_breath_feature<-dd_n277[,-c(1:4)]
dd_gcgc_breath_feature<-data.frame(dd_gcgc_breath_feature)
non_zero<-function(x){
(sum(as.numeric(x !=0),na.rm=TRUE)/length(x))*100
}


for(i in 1:nrow(dd_n277)){
if(as.numeric(substr(as.character(dd_n277$Date_Tray[i]),1,6)) >=170830 & as.numeric(substr(as.character(dd_n277$Date_Tray[i]),1,6)) <=171019) {dd_n277$Batch_ID[i]<-1}
if(as.numeric(substr(as.character(dd_n277$Date_Tray[i]),1,6)) >=171109 & as.numeric(substr(as.character(dd_n277$Date_Tray[i]),1,6)) <=180328) {dd_n277$Batch_ID[i]<-2}
if(as.numeric(substr(as.character(dd_n277$Date_Tray[i]),1,6)) >=180412 & as.numeric(substr(as.character(dd_n277$Date_Tray[i]),1,6)) <=181203) {dd_n277$Batch_ID[i]<-3}
}
dd_batch_effects<-read.table(file="~/EMBER/Final_Clinical_Data_and_Peak_Table_April2020/Peak_Tables_April_2020/Peak_tables_from_templates/Batch_Effects.csv",header=TRUE,sep=",")
dd_clinical<-read.table(file="~/EMBER/Reintegrated_Data_Nov2020/FINAL_DATASET.csv",header=TRUE,sep=",")
names(dd_clinical)[1]<-"Study_number"
tbl_gcgc<-setDT(dd_n277)
tbl_clinical<-setDT(dd_clinical)
setkey(tbl_gcgc,Patient.ID)
setkey(tbl_clinical,Study_number)
tbl_gcgc_clinical <-tbl_gcgc[tbl_clinical, nomatch=0]
tbl_batch<-setDT(dd_batch_effects)
setkey(tbl_batch,Patient.ID)

tbl_gcgc_clinical_batch <-tbl_gcgc_clinical[tbl_batch, nomatch=0]
dd_gcgc_clinical_batch<-data.frame(tbl_gcgc_clinical_batch)

for(i in 1:nrow(dd_gcgc_clinical_batch)){
if(dd_gcgc_clinical_batch$Diagnosis[i]==1) {dd_gcgc_clinical_batch$Grp[i]<-"Healthy"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==2) {dd_gcgc_clinical_batch$Grp[i]<-"Acute Asthma"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==3) {dd_gcgc_clinical_batch$Grp[i]<-"Acute COPD"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==4) {dd_gcgc_clinical_batch$Grp[i]<-"HF"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==5) {dd_gcgc_clinical_batch$Grp[i]<-"Pneumonia"}
}

gcgc_feature_counts<-apply(dd_gcgc_breath_feature,2,non_zero)
gcgc_feature_counts<-data.frame(gcgc_feature_counts)
names(gcgc_feature_counts)<-"Percent_Non_Zero"
sum(as.numeric(gcgc_feature_counts$Percent_Non_Zero) !=0,na.rm=TRUE)

prop_non_zero<-function(p_r){
sum(as.numeric(gcgc_feature_counts$Percent_Non_Zero) >= p_r,na.rm=TRUE)
}


dd_feature_counts<-data.frame(c(100,98,95,90,85,80,75,70,65,60,55,50),unlist(lapply(c(100,98,95,90,85,80,75,70,65,60,55,50),prop_non_zero)))
names(dd_feature_counts)<-c("Present in Percentage of Samples","N")
dd_feature_counts$`Present in Percentage of Samples`<-c("=100",">=98",">=95",">=90",">=85",">=80",">=75",">=70",">=65",">=60",">=55",">=50")
dd_feature_counts$Percent<-(dd_feature_counts$N/ncol(dd_gcgc_breath_feature))*100
names(dd_feature_counts)[2:3]<-c("Number of Features","Percent of Total Features")

names(dd_n277)[3]<-"Patient_ID"
dd_n277<-data.frame(dd_n277)
dd_gcgc<-dd_n277[,as.character(names(dd_n277)) %in% c("Patient_ID",rownames(gcgc_feature_counts)[which(gcgc_feature_counts$Percent_Non_Zero >80)])]
dd_gcgc<-data.frame(dd_gcgc[,-1])
dd_gcgc<-apply(dd_gcgc,2,log1p)
dd_gcgc<-data.frame(dd_gcgc)



library(limma)
library(mixOmics)
library(PMA)
library(sva)
library(pamr)
#Un-processed Data
edata<-t(dd_gcgc)
#pheno<-data.frame(as.character(dd_n277$Patient_ID),as.factor(dd_n277$Batch_ID))
pheno<-data.frame(as.character(dd_gcgc_clinical_batch$Patient.ID),as.factor(dd_gcgc_clinical_batch$Batch_ID),as.factor(dd_gcgc_clinical_batch$Operator),as.factor(dd_gcgc_clinical_batch$Time.Collected),as.factor(dd_gcgc_clinical_batch$Time.Stored.Wet),as.factor(dd_gcgc_clinical_batch$Time.Stored.Dry),as.factor(dd_gcgc_clinical_batch$Collection.volume),as.factor(dd_gcgc_clinical_batch$Diagnosis))
names(pheno)<-c("Sample_ID","Batch_ID","Operator","Time.Collected","Time.Stored.Wet","Time.Stored.Dry","Collection.volume","Diagnosis")
mod = model.matrix(~as.factor(Diagnosis), data=pheno)
mod0 = model.matrix(~1,data=pheno)
n.sv = num.sv(edata,mod,method="be")
svobj = sva(edata,mod,mod0,n.sv=n.sv)
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")
sum(as.numeric(qValues < 0.05))
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
sum(as.numeric(qValuesSv < 0.05))
dd_Sig_Features_latent<-data.frame(qValuesSv[qValuesSv<0.05])
dd_Sig_Features_no_latent<-data.frame(qValuesSv[qValues<0.05])

#Adjust for Batch Effect using ComBat
#Parametric empirical Bayesian adjustment

batch = pheno$Batch_ID
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
X_test<-t(combat_edata)
X_test<-data.frame(X_test)
dd_tmp<-cbind(pheno,X_test)

tbl_tmp<-setDT(dd_tmp)
tbl_clinical<-setDT(dd_clinical)
setkey(tbl_tmp,Sample_ID)
setkey(tbl_clinical,Study_number)
tbl_tmp_clinical <-tbl_tmp[tbl_clinical, nomatch=0]
dd_wadah_gcgc_clinical<-data.frame(tbl_tmp_clinical)


tbl_gcgc_clinical_batch <-tbl_gcgc_clinical[tbl_batch, nomatch=0]
dd_gcgc_clinical_batch<-data.frame(tbl_gcgc_clinical_batch)

for(i in 1:nrow(dd_gcgc_clinical_batch)){
if(dd_gcgc_clinical_batch$Diagnosis[i]==1) {dd_gcgc_clinical_batch$Grp[i]<-"Healthy"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==2) {dd_gcgc_clinical_batch$Grp[i]<-"Acute Asthma"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==3) {dd_gcgc_clinical_batch$Grp[i]<-"Acute COPD"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==4) {dd_gcgc_clinical_batch$Grp[i]<-"HF"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==5) {dd_gcgc_clinical_batch$Grp[i]<-"Pneumonia"}
}

for(i in 1:nrow(dd_wadah_gcgc_clinical)){
if(dd_wadah_gcgc_clinical$Diagnosis[i]==1) {dd_wadah_gcgc_clinical$Grp[i]<-"Healthy"}
else if(dd_wadah_gcgc_clinical$Diagnosis[i]==2) {dd_wadah_gcgc_clinical$Grp[i]<-"Acute Asthma"}
else if(dd_wadah_gcgc_clinical$Diagnosis[i]==3) {dd_wadah_gcgc_clinical$Grp[i]<-"Acute COPD"}
else if(dd_wadah_gcgc_clinical$Diagnosis[i]==4) {dd_wadah_gcgc_clinical$Grp[i]<-"HF"}
else if(dd_wadah_gcgc_clinical$Diagnosis[i]==5) {dd_wadah_gcgc_clinical$Grp[i]<-"Pneumonia"}
}

dd_gcgc_clinical_batch$Gender..1.1..0.0..<-as.factor(dd_gcgc_clinical_batch$Gender..1.1..0.0..)
dd_gcgc_clinical_batch$CXR.coding..0..normal..1..consolidation..pneumonia...2..septal.lines..heart.failure...3.hyperinflation..COPD.and.asthma.<-as.factor(dd_gcgc_clinical_batch$CXR.coding..0..normal..1..consolidation..pneumonia...2..septal.lines..heart.failure...3.hyperinflation..COPD.and.asthma.)
dd_gcgc_clinical_batch$Eosinophilic.asthma.and.COPD..Eos.0.3.<-as.factor(dd_gcgc_clinical_batch$Eosinophilic.asthma.and.COPD..Eos.0.3.)

#Train, Discovery
X<-X_test
X$Data_Set<-dd_gcgc_clinical_batch$Discovery.1
X$Grp<-as.character(dd_gcgc_clinical_batch$Grp)
#X$Grp<-as.character(dd_gcgc_clinical_batch$BreathlessVsnonbreathless..breathless..1..nonbreathless.0.)
X$Patient.ID<-dd_gcgc_clinical_batch$Patient.ID
Y<-as.factor(X$Grp)
X<-X[,-which(names(X) %in% c("Data_Set","Grp","Batch_ID","Patient.ID"))]




#EN Model
library(glmnetUtils)
library(glmnet)
library(plotmo)
library(caret)
library(doParallel)
library(foreach)
library(c060)
registerDoParallel(8)
setwd("/lustre/ahome3/m/mr251/EMBER/Reintegrated_Data_Nov2020/RA_Air_Supply_Feb2021/RA_EN")
library(stabs)
path_res<-stabpath(Y,as.matrix(X_red),weakness=0.7,mc.cores=2,family = "multinomial")
plot(path_res)
stabsel(path_res,error=0.05,type="pfer")

alpha_est<-0.5
#Repeated cross-validation, 5 clinical groups
for(i in 1:100){
set.seed(.Random.seed[i])
cvfit=cv.glmnet(as.matrix(X) ,Y, family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
m_coefs<-coef(cvfit, cvfit$lambda.1se)
dd_coefs_Healthy<-data.frame(m_coefs$Healthy[which(abs(m_coefs$Healthy[,1])>0),],i)
write.table(dd_coefs_Healthy,file="dd_coefs_Healthy_EN_grouped.csv",append=T,col.names=F,sep=",")
dd_coefs_Asthma<-data.frame(m_coefs$`Acute Asthma`[which(abs(m_coefs$`Acute Asthma`[,1])>0),],i)
write.table(dd_coefs_Asthma,file="dd_coefs_Asthma_EN_grouped.csv",append=T,col.names=F,sep=",")
dd_coefs_COPD<-data.frame(m_coefs$`Acute COPD`[which(abs(m_coefs$`Acute COPD`[,1])>0),],i)
write.table(dd_coefs_COPD,file="dd_coefs_COPD_EN_grouped.csv",append=T,col.names=F,sep=",")
dd_coefs_HF<-data.frame(m_coefs$HF[which(abs(m_coefs$HF[,1])>0),],i)
write.table(dd_coefs_HF,file="dd_coefs_HF_EN_grouped.csv",append=T,col.names=F,sep=",")
dd_coefs_Pneumonia<-data.frame(m_coefs$Pneumonia[which(abs(m_coefs$Pneumonia[,1])>0),],i)
write.table(dd_coefs_Pneumonia,file="dd_coefs_Pneumonia_EN_grouped.csv",append=T,col.names=F,sep=",")
}

library(sqldf)
dd_Feature_Healthy<-read.table(file="dd_coefs_Healthy_EN_grouped.csv",col.names=c("Feature","Coef","Sim_Num"),sep=",")
dd_Feature_set_Healthy<-sqldf("select Feature, count(Feature) as Feature_Count, sum(coef)/count(Feature) as av_coef from dd_Feature_Healthy group by Feature having Feature_Count >=80")
dd_Feature_set_Healthy$Feature<-gsub("-","_",dd_Feature_set_Healthy$Feature)
cat("dd_gcgc_clinical_batch$Healthy_Score[i]<-sum(",gsub("\\\\(Intercept)\\\\[i]\\\\*","",paste("(dd_gcgc_clinical_batch$",dd_Feature_set_Healthy$Feature,"[i]*",dd_Feature_set_Healthy$av_coef,")",",",sep=""),")"),")")

dd_Feature_Asthma<-read.table(file="dd_coefs_Asthma_EN_grouped.csv",col.names=c("Feature","Coef","Sim_Num"),sep=",")
dd_Feature_set_Asthma<-sqldf("select Feature, count(Feature) as Feature_Count, sum(coef)/count(Feature) as av_coef from dd_Feature_Asthma group by Feature having Feature_Count >=80")
dd_Feature_set_Asthma$Feature<-gsub("-","_",dd_Feature_set_Asthma$Feature)
cat("dd_gcgc_clinical_batch$Asthma_Score[i]<-sum(",gsub("\\\\(Intercept)\\\\[i]\\\\*","",paste("(dd_gcgc_clinical_batch$",dd_Feature_set_Asthma$Feature,"[i]*",dd_Feature_set_Asthma$av_coef,")",",",sep=""),")"),")")

dd_Feature_COPD<-read.table(file="dd_coefs_COPD_EN_grouped.csv",col.names=c("Feature","Coef","Sim_Num"),sep=",")
dd_Feature_set_COPD<-sqldf("select Feature, count(Feature) as Feature_Count, sum(coef)/count(Feature) as av_coef from dd_Feature_COPD group by Feature having Feature_Count >=80")
dd_Feature_set_COPD$Feature<-gsub("-","_",dd_Feature_set_COPD$Feature)
cat("dd_gcgc_clinical_batch$COPD_Score[i]<-sum(",gsub("\\\\(Intercept)\\\\[i]\\\\*","",paste("(dd_gcgc_clinical_batch$",dd_Feature_set_COPD$Feature,"[i]*",dd_Feature_set_COPD$av_coef,")",",",sep=""),")"),")")

dd_Feature_HF<-read.table(file="dd_coefs_HF_EN_grouped.csv",col.names=c("Feature","Coef","Sim_Num"),sep=",")
dd_Feature_set_HF<-sqldf("select Feature, count(Feature) as Feature_Count, sum(coef)/count(Feature) as av_coef from dd_Feature_HF group by Feature having Feature_Count >=80")
dd_Feature_set_HF$Feature<-gsub("-","_",dd_Feature_set_HF$Feature)
cat("dd_gcgc_clinical_batch$HF_Score[i]<-sum(",gsub("\\\\(Intercept)\\\\[i]\\\\*","",paste("(dd_gcgc_clinical_batch$",dd_Feature_set_HF$Feature,"[i]*",dd_Feature_set_HF$av_coef,")",",",sep=""),")"),")")

dd_Feature_Pneumonia<-read.table(file="dd_coefs_Pneumonia_EN_grouped.csv",col.names=c("Feature","Coef","Sim_Num"),sep=",")
dd_Feature_set_Pneumonia<-sqldf("select Feature, count(Feature) as Feature_Count, sum(coef)/count(Feature) as av_coef from dd_Feature_Pneumonia group by Feature having Feature_Count >=80")
dd_Feature_set_Pneumonia$Feature<-gsub("-","_",dd_Feature_set_Pneumonia$Feature)
cat("dd_gcgc_clinical_batch$Pneumonia_Score[i]<-sum(",gsub("\\\\(Intercept)\\\\[i]\\\\*","",paste("(dd_gcgc_clinical_batch$",dd_Feature_set_Pneumonia$Feature,"[i]*",dd_Feature_set_Pneumonia$av_coef,")",",",sep=""),")"),")")

dd_feature_space<-read.table(file="/home/m/mr251/EMBER/Reintegrated_Data_Nov2020/FeatureList139_green_amber_red_labels_ForMatt301120.csv",header=TRUE,sep=",") 
dd_feature_space<-dd_feature_space[dd_feature_space$feature_group!="Red",]

intersect(dd_Feature_set_Healthy$Feature,dd_feature_space$feature_list)
intersect(dd_Feature_set_Asthma$Feature,dd_feature_space$feature_list)
intersect(dd_Feature_set_COPD$Feature,dd_feature_space$feature_list)
intersect(dd_Feature_set_HF$Feature,dd_feature_space$feature_list)
intersect(dd_Feature_set_Pneumonia$Feature,dd_feature_space$feature_list)

