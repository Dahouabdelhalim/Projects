#module load gcc/6.3
#module load R/4.0.0
library(data.table)


#Main data set REINTEGRATED
setwd("/lustre/ahome3/m/mr251/EMBER/Reintegrated_Data_Nov2020")
dd_n277<-read.table(file="templaten277_peaktable277_vol1p4_SNR30_FILTERED_TO_278_FEATURES_with_9_FEATURES_REINTEGRATED.csv",header=TRUE,sep=",")


#Intitial Data Set
setwd("/lustre/ahome3/m/mr251/EMBER/Final_Clinical_Data_and_Peak_Table_April2020/Peak_Tables_April_2020/Peak_tables_from_templates")
dd_n277<-read.table(file="templaten277_peaktable277_vol1p4_SNR30.csv",header=TRUE,sep=",")

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
dd_gcgc<-dd_n277[,as.character(names(dd_n277)) %in% c("Patient_ID",rownames(gcgc_feature_counts)[which(gcgc_feature_counts$Percent_Non_Zero >=2)])]
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



tbl_gcgc_clinical_batch <-tbl_gcgc_clinical[tbl_batch, nomatch=0]
dd_gcgc_clinical_batch<-data.frame(tbl_gcgc_clinical_batch)

for(i in 1:nrow(dd_gcgc_clinical_batch)){
if(dd_gcgc_clinical_batch$Diagnosis[i]==1) {dd_gcgc_clinical_batch$Grp[i]<-"Healthy"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==2) {dd_gcgc_clinical_batch$Grp[i]<-"Acute Asthma"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==3) {dd_gcgc_clinical_batch$Grp[i]<-"Acute COPD"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==4) {dd_gcgc_clinical_batch$Grp[i]<-"HF"}
else if(dd_gcgc_clinical_batch$Diagnosis[i]==5) {dd_gcgc_clinical_batch$Grp[i]<-"Pneumonia"}
}


#Aditional Clinical Variables
dd_gcgc_clinical_batch$Gender..1.1..0.0..<-as.factor(dd_gcgc_clinical_batch$Gender..1.1..0.0..)
dd_gcgc_clinical_batch$CXR.coding..0..normal..1..consolidation..pneumonia...2..septal.lines..heart.failure...3.hyperinflation..COPD.and.asthma.<-as.factor(dd_gcgc_clinical_batch$CXR.coding..0..normal..1..consolidation..pneumonia...2..septal.lines..heart.failure...3.hyperinflation..COPD.and.asthma.)
dd_gcgc_clinical_batch$Eosinophilic.asthma.and.COPD..Eos.0.3.<-as.factor(dd_gcgc_clinical_batch$Eosinophilic.asthma.and.COPD..Eos.0.3.)

#Train, Discovery
X<-X_test
X$Data_Set<-dd_gcgc_clinical_batch$Discovery.1
X$Grp<-as.character(dd_gcgc_clinical_batch$Grp)
#X$Grp<-as.character(dd_gcgc_clinical_batch$BreathlessVsnonbreathless..breathless..1..nonbreathless.0.)
X$Patient.ID<-dd_gcgc_clinical_batch$Patient.ID
X$Age<-dd_gcgc_clinical_batch$Age
X$Gender<-dd_gcgc_clinical_batch$Gender..1.1..0.0..
X$BNP<-log1p(dd_gcgc_clinical_batch$V1.BNP)
X$CRP<-log1p(dd_gcgc_clinical_batch$V1.CRP)
X$Troponin<-log1p(dd_gcgc_clinical_batch$V1.Troponin)
X$EOS<-dd_gcgc_clinical_batch$Eosinophilic.asthma.and.COPD..Eos.0.3.
X$X_Ray<-dd_gcgc_clinical_batch$CXR.coding..0..normal..1..consolidation..pneumonia...2..septal.lines..heart.failure...3.hyperinflation..COPD.and.asthma.
X<-subset(X,X$Data_Set=="Discovery")
X<-X[X$Grp!="Healthy",]
Y<-as.factor(X$Grp)
X<-X[,-which(names(X) %in% c("Data_Set","Grp","Batch_ID","Patient.ID"))]


#Test, Replication
X_t<-X_test
X_t$Data_Set<-dd_gcgc_clinical_batch$Discovery.1
X_t$Grp<-as.character(dd_gcgc_clinical_batch$Grp)
#X_t$Grp<-as.character(dd_gcgc_clinical_batch$BreathlessVsnonbreathless..breathless..1..nonbreathless.0.)
X_t$Patient.ID<-dd_gcgc_clinical_batch$Patient.ID
X_t$Age<-dd_gcgc_clinical_batch$Age
X_t$Gender<-dd_gcgc_clinical_batch$Gender..1.1..0.0..
X_t$BNP<-log1p(dd_gcgc_clinical_batch$V1.BNP)
X_t$CRP<-log1p(dd_gcgc_clinical_batch$V1.CRP)
X_t$Troponin<-log1p(dd_gcgc_clinical_batch$V1.Troponin)
X_t$EOS<-dd_gcgc_clinical_batch$Eosinophilic.asthma.and.COPD..Eos.0.3.
X_t$X_Ray<-dd_gcgc_clinical_batch$CXR.coding..0..normal..1..consolidation..pneumonia...2..septal.lines..heart.failure...3.hyperinflation..COPD.and.asthma.
X_t<-subset(X_t,X_t$Data_Set=="Replication")
X_t<-X_t[X_t$Grp!="Healthy",]
Y_t<-as.factor(X_t$Grp)
X_t<-X_t[,-which(names(X_t) %in% c("Data_Set","Grp","Batch_ID","Patient.ID"))]


dd_feature_space<-read.table(file="FeatureList139_green_amber_red_labels_ForMatt301120.csv",header=TRUE,sep=",") 
dd_adjusted_Mike<-dd_wadah_gcgc_clinical[,names(dd_wadah_gcgc_clinical) %in% c(dd_feature_space$feature_list,"Grp","Discovery.1")]
write.table(dd_adjusted_Mike,file="Mike_Red_Green_Amber_Batch_Adjusted.csv",row.names=FALSE,sep=",")
dd_feature_space<-dd_feature_space[dd_feature_space$feature_group!="Red",]
dd_feature_space<-dd_feature_space[dd_feature_space$feature_group=="Green",]
dd_feature_space<-dd_feature_space[dd_feature_space$feature_group=="Amber",]
X<-as.data.frame(X)
X_red<-X[,names(X) %in% c(dd_feature_space$feature_list,"Grp","Patient.ID")]
X_red<-X[,names(X) %in% c(dd_feature_space$feature_list,"Age","Gender","BNP","CRP","Troponin","EOS","X_Ray")]
X_red<-X[,names(X) %in% c(dd_feature_space$feature_list)]
X_red<-X[,names(X) %in% c("Age","Gender","BNP","CRP","Troponin","EOS","X_Ray")]
X_red<-X[,names(X) %in% c(dd_feature_space$feature_list,"BNP","CRP","Troponin")]
X_red<-X[,names(X) %in% c("BNP","CRP","Troponin")]
X_t<-as.data.frame(X_t)
X_red_t<-X_t[,names(X_t) %in% c(dd_feature_space$feature_list,"Grp","Patient.ID")]
X_red_t<-X_t[,names(X_t) %in% c(dd_feature_space$feature_list,"Age","Gender","BNP","CRP","Troponin","EOS","X_Ray")]
X_red_t<-X_t[,names(X_t) %in% c(dd_feature_space$feature_list)]
X_red_t<-X_t[,names(X_t) %in% c("Age","Gender","BNP","CRP","Troponin","EOS","X_Ray")]
X_red_t<-X_t[,names(X_t) %in% c(dd_feature_space$feature_list,"BNP","CRP","Troponin")]
X_red_t<-X_t[,names(X_t) %in% c("BNP","CRP","Troponin")]

#Model using 773 features, overfitted
library(glmnetUtils)
library(glmnet)
library(plotmo)
library(caret)
library(doParallel)
library(foreach)
library(c060)
registerDoParallel(8)

alpha_est<-0
Y_bin<-as.numeric(Y=="Pneumonia")
Y_t_bin<-as.numeric(Y_t=="Pneumonia")
r_s<-sample(seq(1,222,by=1),111,replace=FALSE)
X_r<-X[r_s,]
X_rt<-X[-r_s,]
Y_r<-Y[r_s]
Y_rt<-Y[-r_s]
cfit_assess=cv.glmnet(as.matrix(X),Y_bin,family="binomial",type.measure="auc",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
assess_disc<-assess.glmnet(cfit_assess,newx=as.matrix(X),newy=Y_bin)
m_coefs<-coef(cfit_assess, cfit_assess$lambda.1se)
assess_rep<-assess.glmnet(cfit_assess,newx=as.matrix(X_t),newy=Y_t_bin)
dd_pred_disc<-data.frame(Y_bin,predict(cfit_assess,newx=as.matrix(X), s = c("lambda.min"),type=c("link")))
dd_pred_rep<-data.frame(Y_t_bin,predict(cfit_assess,newx=as.matrix(X_t), s = c("lambda.min"),type=c("link")))
library(pROC)
sc1<-roc(dd_pred_disc$Y_bin,dd_pred_disc$X1)
sc2<-roc(dd_pred_rep$Y_t_bin,dd_pred_rep$X1)
plot(sc1,print.auc=TRUE,print.auc.x=0.7, print.auc.y=0.95)
plot(sc2, add=TRUE,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.5,col=c("blue"))
legend("bottomright", legend=c("Discovery", "Replication"),col=c("black", "blue"),lwd=1.5)

cfit <- cv.glmnet(as.matrix(X), Y_bin, family = "binomial", type.measure = "auc",keep = TRUE,alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
rocs <- roc.glmnet(cfit$fit.preval, newy = Y_bin)
best <- cfit$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")
#

X_red$Sub_ID<-paste(X_red$Patient.ID,X_red$Grp)
X_red_t$Sub_ID<-paste(X_red_t$Patient.ID,X_red_t$Grp)
dd_tda_disc<-data.frame(t(X_red))
dd_tda_rep<-data.frame(t(X_red_t))
names(dd_tda_disc)<-dd_tda_disc[104,]
dd_tda_disc<-dd_tda_disc[c(-102,-103,-104),]
names(dd_tda_rep)<-dd_tda_rep[104,]
dd_tda_rep<-dd_tda_rep[c(-102,-103,-104),]
dd_tda_disc$Feature<-rownames(dd_tda_disc)
dd_tda_rep$Feature<-rownames(dd_tda_rep)
write.table(dd_tda_disc,file="Discovery_TDA_Features_Rows.csv",row.names=FALSE,sep=",")
write.table(dd_tda_rep,file="Replication_TDA_Features_Rows.csv",row.names=FALSE,sep=",")

X_red_train<-data.frame(X_red,Y)
X_red_train<-X_red_train[complete.cases(X_red_train),]
X_red_test<-data.frame(X_red_t,Y_t)
X_red_test<-X_red_test[complete.cases(X_red_test),]
names(X_red_test)[ncol(X_red_test)]<-"Y"

#EN Model
library(glmnetUtils)
library(glmnet)
library(plotmo)
library(caret)
library(doParallel)
library(foreach)
library(c060)
registerDoParallel(8)
setwd("/lustre/ahome3/m/mr251/EMBER/Reintegrated_Data_Nov2020/EN_Model_Green_Amber")
library(stabs)
path_res<-stabpath(Y,as.matrix(X_red),weakness=0.7,mc.cores=2,family = "multinomial")
plot(path_res)
stabsel(path_res,error=0.05,type="pfer")

#Single cross-validation, Breathlessness
breath_features<-read.table(file="Breathlessness_Features.csv",header=TRUE,sep=",")
X<-as.data.frame(X)
X_red<-X[,names(X) %in% c(rownames(breath_features))]
X_t<-as.data.frame(X_t)
X_red_t<-X_t[,names(X_t) %in% c(rownames(breath_features))]

alpha_fit<-cva.glmnet(as.matrix(X_red),Y, family="binomial")
plot(alpha_fit)

alpha_est<-0
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="binomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
m_coefs<-coef(cfit_assess, cfit_assess$lambda.1se)
write.table(data.frame(as.matrix(m_coefs)),file="Breathlessness_Features.csv",sep=",")
assess.glmnet(cfit_assess,newx=as.matrix(X_red),newy=Y)
assess.glmnet(cfit_assess,newx=as.matrix(X_red_t),newy=Y_t)
dd_disc<-data.frame(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min")),Y)
dd_rep<-data.frame(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min")),Y_t)
names(dd_disc)<-c("Breathlessness_Score","Grp")
names(dd_rep)<-c("Breathlessness_Score","Grp")
dd_disc_prob<-data.frame(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min"),type=c("response")),Y)
dd_rep_prob<-data.frame(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min"),type=c("response")),Y_t)
names(dd_disc_prob)<-c("Breathlessness_Prob","Grp")
names(dd_rep_prob)<-c("Breathlessness_Prob","Grp")

confusionMatrix(as.factor(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min"),type=c("class"))),Y)
confusionMatrix(as.factor(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min"),type=c("class"))),Y_t)

library(ggpubr)

p1<-ggviolin(dd_disc_prob, x = "Grp",
          y = c("Breathlessness_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Breathlessness Probability", 
          add = "boxplot")

p2<-ggviolin(dd_rep_prob, x = "Grp",
          y = c("Breathlessness_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Breathlessness Probability", 
          add = "boxplot")
tiff(filename = "Breathlessnsess_Probs.tiff" ,units="in", width=18, height=9, res=300,compression = "lzw")
ggarrange(p1,p2)
dev.off()

library(pROC)
tiff(filename = "Green_Amber_Features_Breathless2.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
breath_score_disc<-roc(as.numeric(dd_disc$Grp=="1"),dd_disc$Breathlessness_Score,ci = TRUE)
breath_score_rep<-roc(as.numeric(dd_rep$Grp=="1"),dd_rep$Breathlessness_Score,ci = TRUE)
plot(breath_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=1,print.auc.y=0.9)
plot(breath_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.8)
legend("bottomright", legend=c("Breathlessness Score Discovery", "Breathlessness Score Replication"),col=c("black", "blue"),lwd=0.7)
dev.off()
dd_roc_disc<-data.frame(reportROC(as.numeric(dd_disc$Grp=="1"),dd_disc$Breathlessness_Score))
dd_roc_rep<-data.frame(reportROC(as.numeric(dd_rep$Grp=="1"),dd_rep$Breathlessness_Score))
write.table(dd_roc_disc,file="Breathlessness_ROC_Discovery.csv",sep=",")
write.table(dd_roc_rep,file="Breathlessness_ROC_Replication.csv",sep=",")

X_un<-X_test
X_un$Data_Set<-dd_gcgc_clinical_batch$Discovery.1
X_un$Grp<-as.character(dd_gcgc_clinical_batch$BreathlessVsnonbreathless..breathless..1..nonbreathless.0.)
X_un$Grp2<-as.character(dd_gcgc_clinical_batch$Grp)
X_un$Patient.ID<-dd_gcgc_clinical_batch$Patient.ID
X_un$Diagnostic_uncertainty<-dd_gcgc_clinical_batch$Diagnostic.uncertainty.higher.than.upper.quartile..20.
X_un1<-X_un[X_un$Diagnostic_uncertainty==1 & !is.na(X_un$Diagnostic_uncertainty),]
X_un2<-X_un[X_un$Grp2=="Healthy",]
X_un<-rbind(X_un1,X_un2)
X_un<-data.frame(X_un$Grp,predict(cfit_assess,newx=as.matrix(X_un[,names(X_un) %in% c(dd_feature_space$feature_list)]), s = c("lambda.min")))
names(X_un)<-c("Grp","Breathlessness_Score")

library(pROC)
tiff(filename = "Green_Amber_Features_Breathless_UN.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
breath_score_un<-roc(as.numeric(X_un$Grp=="1"),X_un$Breathlessness_Score,ci = TRUE)
plot(breath_score_un, col=c("black"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.8, print.auc.y=0.8,main="Breathlessness Score Diagnostic Uncertainty High verus Healthy")
dev.off()
library(reportROC)
dd_roc_un<-data.frame(reportROC(as.numeric(X_un$Grp=="1"),X_un$Breathlessness_Score))
write.table(dd_roc_un,file="Breathlessness_ROC_UN.csv",sep=",")


predictions<-predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.1se"),type="class")
library(DescTools)
library(jmuOutlier)
my_acc<-c()
my_v2<-c()
my_jac<-c()
for(i in 1:1000){
#set.seed(i)
Y_r<-as.factor(sample(Y,length(Y),replace = FALSE))
my_v2[i]<-CramerV(as.matrix(data.frame(as.numeric(Y),as.numeric(Y_r))),method = c("fisheradj"),correct = FALSE)
my_acc[i]<-mean(predictions==Y_r)
}

tiff(filename = "Model_Accuracy_Breath.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(my_v2,my_acc,ylim=c(0:1),xlab="Cramer's V, Y, Y rand",ylab="Model Accuracy")
abline(h=mean(predictions==Y),col="red")
dev.off()

perm.test(mean(predictions==Y),my_acc)

predictions<-predict(cfit_assess,newx=as.matrix(X_red_t),s = c("lambda.1se"),type="class")
library(DescTools)
library(jmuOutlier)
my_acc<-c()
my_v2<-c()
my_jac<-c()
for(i in 1:1000){
#set.seed(i)
Y_r<-as.factor(sample(Y_t,length(Y_t),replace = FALSE))
my_v2[i]<-CramerV(as.matrix(data.frame(as.numeric(Y_t),as.numeric(Y_r))),method = c("fisheradj"),correct = FALSE)
my_acc[i]<-mean(predictions==Y_r)
}

tiff(filename = "Model_Accuracy_Breath_Rep.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(my_v2,my_acc,ylim=c(0:1),xlab="Cramer's V, Y, Y rand",ylab="Model Accuracy")
abline(h=mean(predictions==Y_t),col="red")
dev.off()

perm.test(mean(predictions==Y_t),my_acc)

requireNamespace("mlr3measures")
mlr3measures::bacc(truth, response)

m_err<-c()
for(i in 1:1000){
alpha_est<-0
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="binomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red),newy=Y)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.1se"),type="class"))
m_err[i]<-mlr3measures::bacc(Y, predictions)
}

m_err_r<-c()
for(i in 1:1000){
alpha_est<-0
Y_r<-as.factor(sample(Y,length(Y),replace = FALSE))
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="binomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err_r[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red),newy=Y_r)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.1se"),type="class"))
m_err_r[i]<-mlr3measures::bacc(Y_r, predictions)
}
tiff(filename = "Model_Accuracy_Breath_EN.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(seq(1:1000),m_err,ylim=c(0:1),xlab="Index",ylab="Balanced Accuracy",type="l",col="green")
lines(seq(1:1000),m_err_r,col="red")
dev.off()

m_err<-c()
for(i in 1:1000){
alpha_est<-0
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="binomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red_t),newy=Y_t)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red_t),s = c("lambda.1se"),type="class"))
m_err[i]<-mlr3measures::bacc(Y_t, predictions)
}

m_err_r<-c()
for(i in 1:1000){
alpha_est<-0
Y_r<-as.factor(sample(Y_t,length(Y_t),replace = FALSE))
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="binomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err_r[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red_t),newy=Y_r)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red_t),s = c("lambda.1se"),type="class"))
m_err_r[i]<-mlr3measures::bacc(Y_r, predictions)
}
tiff(filename = "Model_Accuracy_Breath_EN_Rep.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(seq(1:1000),m_err,ylim=c(0:1),xlab="Index",ylab="Balanced Accuracy",type="l",col="green")
lines(seq(1:1000),m_err_r,col="red")
dev.off()

#Adaptive Elastic Net
library(msaenet)
ad_fit<-msaenet(as.matrix(X_red),Y, family = c("binomial"),alphas = seq(0.1, 0.9, 0.1),nsteps=20L,tune.nsteps="ebic",seed=1968)
ad_fit$best.step
msaenet.nzv(ad_fit)
coef(ad_fit)
msaenet.nzv.all(ad_fit)
plot(ad_fit, label = TRUE)
plot(ad_fit, type = "dotplot", label = TRUE, label.cex = 1)

#Single cross-validation, 5 clinical groups
alpha_fit<-cva.glmnet(as.matrix(X_red),Y, family="multinomial")
plot(alpha_fit)
alpha_est<-0
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
m_coefs<-coef(cfit_assess, cfit_assess$lambda.1se)
write.table(data.frame(as.matrix(m_coefs$`Acute Asthma`)),file="Asthma_Features.csv",sep=",")
write.table(data.frame(as.matrix(m_coefs$`Acute COPD`)),file="COPD_Features.csv",sep=",")
write.table(data.frame(as.matrix(m_coefs$Healthy)),file="Healthy_Features.csv",sep=",")
write.table(data.frame(as.matrix(m_coefs$HF)),file="HF_Features.csv",sep=",")
write.table(data.frame(as.matrix(m_coefs$Pneumonia)),file="Pneumonia_Features.csv",sep=",")
assess.glmnet(cfit_assess,newx=as.matrix(X_red),newy=Y)
assess.glmnet(cfit_assess,newx=as.matrix(X_red_t),newy=Y_t)
dd_disc<-data.frame(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min")),Y)
dd_rep<-data.frame(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min")),Y_t)
names(dd_disc)<-c("Asthma_Score","COPD_Score","Healthy_Score","HF_Score","Pneumonia_Score","Grp")
names(dd_rep)<-c("Asthma_Score","COPD_Score","Healthy_Score","HF_Score","Pneumonia_Score","Grp")
dd_disc$Patient.ID<-X$Patient.ID
dd_rep$Patient.ID<-X_t$Patient.ID
dd_disc$Data_Set<-X$Data_Set
dd_rep$Data_Set<-X_t$Data_Set
dd_scores_disc_rep<-rbind(dd_disc,dd_rep)
write.table(dd_scores_disc_rep,file="Scores_Disc_Rep.csv",row.names=FALSE,sep=",")
dd_disc_prob<-data.frame(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min"),type=c("response")),Y)
dd_rep_prob<-data.frame(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min"),type=c("response")),Y_t)
names(dd_disc_prob)<-c("Asthma_Prob","COPD_Prob","Healthy_Prob","HF_Prob","Pneumonia_Prob","Grp")
names(dd_rep_prob)<-c("Asthma_Prob","COPD_Prob","Healthy_Prob","HF_Prob","Pneumonia_Prob","Grp")

library(bdpv)
confusionMatrix(as.factor(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min"),type=c("class"))),Y)
tt_d<-confusionMatrix(as.factor(predict(cfit_assess,newx=as.matrix(X_red), s = c("lambda.min"),type=c("class"))),Y,mode = "everything")
ci_lst<-list()
for(i in 1:length(levels(Y))){
m_d<-matrix(c(tt_d$table[i,i],sum(tt_d$table[-i,i]),sum(tt_d$table[i,-i]),sum(tt_d$table[-i,-i])),ncol=2)
ci_d<-BDtest(m_d,pr=sum(Y==levels(Y)[i])/length(Y),conf.level = 0.95)
ci_lst[[i]]<-list(ci_d$SESPDAT)
}
names(ci_lst)<-levels(Y)
dd_ci_d<-data.frame(ci_lst)

confusionMatrix(as.factor(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min"),type=c("class"))),Y_t)
tt_r<-confusionMatrix(as.factor(predict(cfit_assess,newx=as.matrix(X_red_t), s = c("lambda.min"),type=c("class"))),Y_t,mode = "everything")
ci_lst_r<-list()
for(i in 1:length(levels(Y_t))){
m_r<-matrix(c(tt_r$table[i,i],sum(tt_r$table[-i,i]),sum(tt_r$table[i,-i]),sum(tt_r$table[-i,-i])),ncol=2)
ci_r<-BDtest(m_r,pr=sum(Y_t==levels(Y_t)[i])/length(Y_t),conf.level = 0.95)
ci_lst_r[[i]]<-list(ci_r$SESPDAT)
}
names(ci_lst_r)<-levels(Y_t)
dd_ci_r<-data.frame(ci_lst_r)
write.table(dd_ci_d,file="Discovery_Ridge_CIs.csv",sep=",")
write.table(dd_ci_r,file="Replication_Ridge_CIs.csv",sep=",")


predictions<-predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.min"),type="class")
predictions<-predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.1se"),type="class")
library(DescTools)
library(jmuOutlier)
my_acc<-c()
my_v2<-c()
my_jac<-c()
for(i in 1:1000){
#set.seed(i)
Y_r<-as.factor(sample(Y,length(Y),replace = FALSE))
my_v2[i]<-CramerV(as.matrix(data.frame(as.numeric(Y),as.numeric(Y_r))),method = c("fisheradj"),correct = FALSE)
my_acc[i]<-mean(predictions==Y_r)
}

tiff(filename = "Model_Accuracy.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(my_v2,my_acc,ylim=c(0:1),xlab="Cramer's V, Y, Y rand",ylab="Model Accuracy")
abline(h=mean(predictions==Y),col="red")
dev.off()

perm.test(mean(predictions==Y),my_acc)



predictions<-predict(cfit_assess,newx=as.matrix(X_red_t),s = c("lambda.1se"),type="class")
library(DescTools)
library(jmuOutlier)
my_acc<-c()
my_v2<-c()
my_jac<-c()
for(i in 1:1000){
#set.seed(i)
Y_r<-as.factor(sample(Y_t,length(Y_t),replace = FALSE))
my_v2[i]<-CramerV(as.matrix(data.frame(as.numeric(Y_t),as.numeric(Y_r))),method = c("fisheradj"),correct = FALSE)
my_acc[i]<-mean(predictions==Y_r)
}

tiff(filename = "Model_Accuracy_Rep.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(my_v2,my_acc,ylim=c(0:1),xlab="Cramer's V, Y, Y rand",ylab="Model Accuracy")
abline(h=mean(predictions==Y_t),col="red")
dev.off()

perm.test(mean(predictions==Y_t),my_acc)

m_err<-c()
for(i in 1:1000){
alpha_est<-0
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red),newy=Y)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.1se"),type="class"))
m_err[i]<-mlr3measures::bacc(Y, predictions)
}

m_err_r<-c()
for(i in 1:1000){
alpha_est<-0
Y_r<-as.factor(sample(Y,length(Y),replace = FALSE))
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err_r[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red),newy=Y_r)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red),s = c("lambda.1se"),type="class"))
m_err_r[i]<-mlr3measures::bacc(Y_r, predictions)
}
tiff(filename = "Model_Accuracy_EN.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(seq(1:1000),m_err,ylim=c(0:1),xlab="Index",ylab="Balanced Accuracy",type="l",col="green")
lines(seq(1:1000),m_err_r,col="red")
dev.off()

m_err<-c()
for(i in 1:1000){
alpha_est<-0
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red_t),newy=Y_t)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red_t),s = c("lambda.1se"),type="class"))
m_err[i]<-mlr3measures::bacc(Y_t, predictions)
}

m_err_r<-c()
for(i in 1:1000){
alpha_est<-0
Y_r<-as.factor(sample(Y_t,length(Y_t),replace = FALSE))
cfit_assess=cv.glmnet(as.matrix(X_red),Y,family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
#m_err_r[i]<-1-assess.glmnet(cfit_assess,newx=as.matrix(X_red_t),newy=Y_r)$class
predictions<-as.factor(predict(cfit_assess,newx=as.matrix(X_red_t),s = c("lambda.1se"),type="class"))
m_err_r[i]<-mlr3measures::bacc(Y_r, predictions)
}
tiff(filename = "Model_Accuracy_EN_Rep.tiff" ,units="in", width=10, height=9, res=300,compression = "lzw")
plot(seq(1:1000),m_err,ylim=c(0:1),xlab="Index",ylab="Balanced Accuracy",type="l",col="green")
lines(seq(1:1000),m_err_r,col="red")
dev.off()


library(ggpubr)
dd_disc_prob$Grp<-gsub("Acute ","",dd_disc_prob$Grp)
p1<-ggviolin(dd_disc_prob, x = "Grp",
          y = c("Asthma_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Asthma Probability", 
          add = "boxplot")

p2<-ggviolin(dd_disc_prob, x = "Grp",
          y = c("COPD_Prob"),
          color = "Grp", palette = "jco",
          ylab = "COPD Probability", 
          add = "boxplot")

p3<-ggviolin(dd_disc_prob, x = "Grp",
          y = c("Healthy_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Healthy Probability", 
          add = "boxplot")

p4<-ggviolin(dd_disc_prob, x = "Grp",
          y = c("HF_Prob"),
          color = "Grp", palette = "jco",
          ylab = "HF Probability", 
          add = "boxplot")

p5<-ggviolin(dd_disc_prob, x = "Grp",
          y = c("Pneumonia_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Pneumonia Probability", 
          add = "boxplot")
tiff(filename = "Disc_Probs.tiff" ,units="in", width=18, height=9, res=300,compression = "lzw")
ggarrange(p1,p2,p3,p4,p5)
dev.off()

dd_rep_prob$Grp<-gsub("Acute ","",dd_rep_prob$Grp)
p1<-ggviolin(dd_rep_prob, x = "Grp",
          y = c("Asthma_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Asthma Probability", 
          add = "boxplot")

p2<-ggviolin(dd_rep_prob, x = "Grp",
          y = c("COPD_Prob"),
          color = "Grp", palette = "jco",
          ylab = "COPD Probability", 
          add = "boxplot")

p3<-ggviolin(dd_rep_prob, x = "Grp",
          y = c("Healthy_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Healthy Probability", 
          add = "boxplot")

p4<-ggviolin(dd_rep_prob, x = "Grp",
          y = c("HF_Prob"),
          color = "Grp", palette = "jco",
          ylab = "HF Probability", 
          add = "boxplot")

p5<-ggviolin(dd_rep_prob, x = "Grp",
          y = c("Pneumonia_Prob"),
          color = "Grp", palette = "jco",
          ylab = "Pneumonia Probability", 
          add = "boxplot")
tiff(filename = "Rep_Probs.tiff" ,units="in", width=18, height=9, res=300,compression = "lzw")
ggarrange(p1,p2,p3,p4,p5)
dev.off()

library(pROC)
tiff(filename = "Green_Amber_Features_Asthma.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
asthma_score_disc<-roc(as.numeric(dd_disc$Grp=="Acute Asthma"),dd_disc$Asthma_Score,ci = TRUE)
asthma_score_rep<-roc(as.numeric(dd_rep$Grp=="Acute Asthma"),dd_rep$Asthma_Score,ci = TRUE)
plot(asthma_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.8,print.auc.y=0.85)
plot(asthma_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.5)
legend("bottomright", legend=c("Asthma Score Discovery", "Asthma Score Replication"),col=c("black", "blue"),lwd=0.7)
dev.off()

tiff(filename = "Green_Amber_Features_COPD.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
copd_score_disc<-roc(as.numeric(dd_disc$Grp=="Acute COPD"),dd_disc$COPD_Score,ci = TRUE)
copd_score_rep<-roc(as.numeric(dd_rep$Grp=="Acute COPD"),dd_rep$COPD_Score,ci = TRUE)
plot(copd_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.9,print.auc.y=0.8)
plot(copd_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.5)
legend("bottomright", legend=c("COPD Score Discovery", "COPD Score Replication"),col=c("black", "blue"),lwd=0.7)
dev.off()

tiff(filename = "Green_Amber_Features_Healthy.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
healthy_score_disc<-roc(as.numeric(dd_disc$Grp=="Healthy"),dd_disc$Healthy_Score,ci = TRUE)
healthy_score_rep<-roc(as.numeric(dd_rep$Grp=="Healthy"),dd_rep$Healthy_Score,ci = TRUE)
plot(healthy_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.87,print.auc.y=0.9)
plot(healthy_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.79, print.auc.y=0.65)
legend("bottomright", legend=c("Healthy Score Discovery", "Healthy Score Replication"),col=c("black", "blue"),lwd=0.7)
dev.off()

tiff(filename = "Green_Amber_Features_HF.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
hf_score_disc<-roc(as.numeric(dd_disc$Grp=="HF"),dd_disc$HF_Score,ci = TRUE)
hf_score_rep<-roc(as.numeric(dd_rep$Grp=="HF"),dd_rep$HF_Score,ci = TRUE)
plot(hf_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5,print.auc.y=0.5)
plot(hf_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.7)
legend("bottomright", legend=c("Heart Failure Score Discovery", "Heart Failure Score Replication"),col=c("black", "blue"),lwd=0.7)
dev.off()

tiff(filename = "Green_Amber_Features_Pneumonia.tiff" ,units="in", width=10, height=7, res=300,compression = "lzw")
pneumonia_score_disc<-roc(as.numeric(dd_disc$Grp=="Pneumonia"),dd_disc$Pneumonia_Score,ci = TRUE)
pneumonia_score_rep<-roc(as.numeric(dd_rep$Grp=="Pneumonia"),dd_rep$Pneumonia_Score,ci = TRUE)
plot(pneumonia_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.9,print.auc.y=0.9)
plot(pneumonia_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.8, print.auc.y=0.45)
legend("bottomright", legend=c("Pneumonia Score Discovery", "Pneumonia Score Replication"),col=c("black", "blue"),lwd=0.7)
dev.off()

asthma_features<-read.table(file="Asthma_Features.csv",header=TRUE,sep=",")
asthma_features$Feature<-rownames(asthma_features)
names(asthma_features)[1]<-"Coef"

copd_features<-read.table(file="COPD_Features.csv",header=TRUE,sep=",")
copd_features$Feature<-rownames(copd_features)
names(copd_features)[1]<-"Coef"

healthy_features<-read.table(file="Healthy_Features.csv",header=TRUE,sep=",")
healthy_features$Feature<-rownames(healthy_features)
names(healthy_features)[1]<-"Coef"

hf_features<-read.table(file="HF_Features.csv",header=TRUE,sep=",")
hf_features$Feature<-rownames(hf_features)
names(hf_features)[1]<-"Coef"

pneumonia_features<-read.table(file="Pneumonia_Features.csv",header=TRUE,sep=",")
pneumonia_features$Feature<-rownames(pneumonia_features)
names(pneumonia_features)[1]<-"Coef"

breath_features<-read.table(file="Breathlessness_Features.csv",header=TRUE,sep=",")
breath_features$Feature<-rownames(breath_features)
names(breath_features)[1]<-"Coef"

asthma_top10<-asthma_features[order(abs(asthma_features$Coef),decreasing=TRUE),][1:10,] 
copd_top10<-copd_features[order(abs(copd_features$Coef),decreasing=TRUE),][1:10,] 
healthy_top10<-healthy_features[order(abs(healthy_features$Coef),decreasing=TRUE),][1:10,] 
hf_top10<-hf_features[order(abs(hf_features$Coef),decreasing=TRUE),][1:10,] 
pneumonia_top10<-pneumonia_features[order(abs(pneumonia_features$Coef),decreasing=TRUE),][1:10,] 
breath_top10<-breath_features[order(abs(breath_features$Coef),decreasing=TRUE),][1:10,] 
tiff(filename = "Green_Amber_Features_Top10_Dists.tiff" ,units="in", width=17, height=10, res=300,compression = "lzw")
ggarrange(p1,p2,p3,p4,p5,p6)
dev.off()



library(ggpubr)

p1<-gghistogram(dd_gcgc[,names(dd_gcgc) %in% asthma_top10$Feature],
       x = asthma_top10$Feature[-grep("Intercept",asthma_top10$Feature)],
       y = "..count..",
       combine = TRUE,                  # Combine the 3 plots
       xlab = "Abundence", 
       add = "median",                  # Add median line. 
       rug = TRUE                       # Add marginal rug
)
p1 <- ggpar(p1, title = "Asthma Top 10 Features Log(x+1) Transformed",subtitle = "All Data")

p2<-gghistogram(dd_gcgc[,names(dd_gcgc) %in% copd_top10$Feature],
       x = copd_top10$Feature[-grep("Intercept",copd_top10$Feature)],
       y = "..count..",
       combine = TRUE,                  # Combine the 3 plots
       xlab = "Abundence", 
       add = "median",                  # Add median line. 
       rug = TRUE                       # Add marginal rug
)
p2 <- ggpar(p2, title = "COPD Top 10 Features Log(x+1) Transformed",subtitle = "All Data")

p3<-gghistogram(dd_gcgc[,names(dd_gcgc) %in% healthy_top10$Feature],
       x = healthy_top10$Feature[-grep("Intercept",healthy_top10$Feature)],
       y = "..count..",
       combine = TRUE,                  # Combine the 3 plots
       xlab = "Abundence", 
       add = "median",                  # Add median line. 
       rug = TRUE                       # Add marginal rug
)
p3 <- ggpar(p3, title = "Healthy Top 10 Features Log(x+1) Transformed",subtitle = "All Data")

p4<-gghistogram(dd_gcgc[,names(dd_gcgc) %in% hf_top10$Feature],
       x = hf_top10$Feature[-grep("Intercept",hf_top10$Feature)],
       y = "..count..",
       combine = TRUE,                  # Combine the 3 plots
       xlab = "Abundence", 
       add = "median",                  # Add median line. 
       rug = TRUE                       # Add marginal rug
)
p4 <- ggpar(p4, title = "HF Top 10 Features Log(x+1) Transformed",subtitle = "All Data")

p5<-gghistogram(dd_gcgc[,names(dd_gcgc) %in% pneumonia_top10$Feature],
       x = pneumonia_top10$Feature[-grep("Intercept",pneumonia_top10$Feature)],
       y = "..count..",
       combine = TRUE,                  # Combine the 3 plots
       xlab = "Abundence", 
       add = "median",                  # Add median line. 
       rug = TRUE                       # Add marginal rug
)
p5 <- ggpar(p5, title = "Pneumonia Top 10 Features Log(x+1) Transformed",subtitle = "All Data")

p6<-gghistogram(dd_gcgc[,names(dd_gcgc) %in% breath_top10$Feature],
       x = breath_top10$Feature[-grep("Intercept",breath_top10$Feature)],
       y = "..count..",
       combine = TRUE,                  # Combine the 3 plots
       xlab = "Abundence", 
       add = "median",                  # Add median line. 
       rug = TRUE                       # Add marginal rug
)
p6 <- ggpar(p6, title = "Breathlessness Top 10 Features Log(x+1) Transformed",subtitle = "All Data")
library(psych)
X_dists<-dd_gcgc[,names(dd_gcgc) %in% dd_feature_space$feature_list]
tiff(filename = "Green_Amber_Features_1_50.tiff" ,units="in", width=17, height=10, res=300,compression = "lzw")
multi.hist(X_dists[,c(1:50)])
dev.off()
tiff(filename = "Green_Amber_Features_51_101.tiff" ,units="in", width=17, height=10, res=300,compression = "lzw")
multi.hist(X_dists[,c(51:101)])
dev.off()
tiff(filename = "Green_Amber_Features_1_101.tiff" ,units="in", width=17, height=13, res=300,compression = "lzw")
multi.hist(X_dists)
dev.off()

#Repeated cross-validation, 5 clinical groups
for(i in 1:100){
set.seed(.Random.seed[i])
cvfit=cv.glmnet(as.matrix(X_red) ,Y, family="multinomial",type.measure="class",alpha=alpha_est,standardize=FALSE,nfolds=10,parallel = TRUE)
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



library(pROC)
dd_disc<-dd_gcgc_clinical_batch[dd_gcgc_clinical_batch$Discovery.1=="Discovery",]
dd_rep<-dd_gcgc_clinical_batch[dd_gcgc_clinical_batch$Discovery.1=="Replication",]

asthma_score_disc<-roc(as.numeric(dd_disc$Grp=="Acute Asthma"),dd_disc$Asthma_Score,ci = TRUE)
asthma_score_rep<-roc(as.numeric(dd_rep$Grp=="Acute Asthma"),dd_rep$Asthma_Score,ci = TRUE)
plot(asthma_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.6,print.auc.y=0.85)
plot(asthma_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.5)
legend("bottomright", legend=c("Asthma Score Discovery", "Asthma Score Replication"),col=c("black", "blue"),lwd=0.7)

copd_score_disc<-roc(as.numeric(dd_disc$Grp=="Acute COPD"),dd_disc$COPD_Score,ci = TRUE)
copd_score_rep<-roc(as.numeric(dd_rep$Grp=="Acute COPD"),dd_rep$COPD_Score,ci = TRUE)
plot(copd_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.79,print.auc.y=0.8)
plot(copd_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.5)
legend("bottomright", legend=c("COPD Score Discovery", "COPD Score Replication"),col=c("black", "blue"),lwd=0.7)

healthy_score_disc<-roc(as.numeric(dd_disc$Grp=="Healthy"),dd_disc$Healthy_Score,ci = TRUE)
healthy_score_rep<-roc(as.numeric(dd_rep$Grp=="Healthy"),dd_rep$Healthy_Score,ci = TRUE)
plot(healthy_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.85,print.auc.y=0.9)
plot(healthy_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.75, print.auc.y=0.65)
legend("bottomright", legend=c("Healthy Score Discovery", "Healthy Score Replication"),col=c("black", "blue"),lwd=0.7)

hf_score_disc<-roc(as.numeric(dd_disc$Grp=="HF"),dd_disc$HF_Score,ci = TRUE)
hf_score_rep<-roc(as.numeric(dd_rep$Grp=="HF"),dd_rep$HF_Score,ci = TRUE)
plot(hf_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.85,print.auc.y=0.5)
plot(hf_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.6, print.auc.y=0.7)
legend("bottomright", legend=c("Heart Failure Score Discovery", "Heart Failure Score Replication"),col=c("black", "blue"),lwd=0.7)

pneumonia_score_disc<-roc(as.numeric(dd_disc$Grp=="Pneumonia"),dd_disc$Pneumonia_Score,ci = TRUE)
pneumonia_score_rep<-roc(as.numeric(dd_rep$Grp=="Pneumonia"),dd_rep$Pneumonia_Score,ci = TRUE)
plot(pneumonia_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.93,print.auc.y=0.9)
plot(pneumonia_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.8, print.auc.y=0.45)
legend("bottomright", legend=c("Pneumonia Score Discovery", "Pneumonia Score Replication"),col=c("black", "blue"),lwd=0.7)

breath_score_disc<-roc(as.numeric(dd_disc$V1a.breathlessness..Y.1..N.0.=="1"),dd_disc$Breathlessness_Score,ci = TRUE)
breath_score_rep<-roc(as.numeric(dd_rep$V1a.breathlessness..Y.1..N.0.=="1"),dd_rep$Breathlessness_Score,ci = TRUE)
plot(breath_score_disc,print.auc.cex=0.8,print.auc=TRUE,print.auc.x=1,print.auc.y=1)
plot(breath_score_rep, add=TRUE,col=c("blue"),print.auc.cex=0.8,print.auc=TRUE,print.auc.x=0.5, print.auc.y=0.8)
legend("bottomright", legend=c("Breathlessness Score Discovery", "Breathlessness Score Replication"),col=c("black", "blue"),lwd=0.7)

#LDA
library(MASS)
library(tidyverse)
library(caret)
library(ggpubr)


X_lda<-data.frame(X_red,Y)
names(X_lda)[ncol(X_lda)]<-c("Y")
X_lda$Y<-as.factor(X_lda$Y)
model_stable <- lda(Y~., data = X_lda)
predictions <- model_stable %>% predict(X_lda)
mean(predictions$class==X_lda$Y)
confusionMatrix(predictions$class,X_lda$Y)
mean(predict(model_stable,X_red_t)$class==Y_t)
confusionMatrix(predict(model_stable,X_red_t)$class,as.factor(Y_t))
lda.data <-cbind(X_lda,predictions$x)


#Boosted classificcation using xgboost
#objective = "multi:softprob"
library(MASS)
library(caret)
library(ggpubr)
library(tidyverse)
library(xgboost)
require(Matrix)
dd_train_sparse <- sparse.model.matrix(Y ~ ., data = X_red_train)[,-1]
dd_test_sparse <- sparse.model.matrix(Y ~ ., data = X_red_test)[,-1]
dimnames(dd_train_sparse)
dimnames(dd_test_sparse)
#Multiclass
mod_xgboost <- xgboost(data = dd_train_sparse, num_class=5,label = (as.integer(X_red_train$Y)-1), max_depth = 4,
               eta = 1, nthread = 2, nrounds = 10,eval_metric = "merror", eval_metric = "mlogloss",objective = "multi:softprob")
predicted.prob <- mod_xgboost %>% predict(dd_train_sparse,reshape=TRUE)
head(predicted.prob)
factor(max.col(predicted.prob)-1)
caret::confusionMatrix(factor(as.integer(X_red_train$Y)-1),factor(max.col(predicted.prob)-1))
importance_matrix <- xgb.importance(model = mod_xgboost)
print(importance_matrix)
head(importance_matrix,10)[,1]
xgb.plot.importance(importance_matrix = importance_matrix,top_n=20,rel_to_first = TRUE)

#Binary
mod_xgboost <- xgboost(data = dd_train_sparse,label = (as.integer(X_red_train$Y)-1), max_depth = 4,
               eta = 1, nthread = 2, nrounds = 10,eval_metric = "error", eval_metric = "logloss",objective = "binary:logistic")
predicted.prob <- mod_xgboost %>% predict(dd_train_sparse,reshape=TRUE)
head(predicted.prob)
as.factor(as.numeric(predicted.prob>0.5))
caret::confusionMatrix(factor(as.integer(X_red_train$Y)-1),as.factor(as.numeric(predicted.prob>0.5)))
importance_matrix <- xgb.importance(model = mod_xgboost)
print(importance_matrix)
head(importance_matrix,10)[,1]
xgb.plot.importance(importance_matrix = importance_matrix,top_n=20,rel_to_first = TRUE)

#Model performance on test data
predicted.prob_test <- mod_xgboost %>% predict(dd_test_sparse,reshape=TRUE)
head(predicted.prob_test)
as.factor(as.numeric(predicted.prob_test>0.5))
caret::confusionMatrix(factor(as.integer(X_red_test$Y)-1),as.factor(as.numeric(predicted.prob_test>0.5)))



#Through caret
set.seed(154)
mod_xgboost <- train(
  Y ~., data =X_red_train , method = "xgbTree",
  trControl = trainControl("cv", number = 15)
  )

mod_xgboost <- train(
  Y ~., data =X_red_train , method = "xgbTree",
  trControl = trainControl(method = "repeatedcv", number = 20,repeats = 10)
  )
predicted.classes <- mod_xgboost %>% predict(X_red_train)
mean(predicted.classes == X_red_train$Y)
caret::confusionMatrix(predicted.classes ,X_red_train$Y)
predicted.classes_test <- mod_xgboost %>% predict(X_red_test)
caret::confusionMatrix(predicted.classes_test ,X_red_test$Y)
varImp(mod_xgboost)
plot(varImp(mod_xgboost),top=20)
mod_xgboost$finalModel
#



#Advanced interface
dd_train = xgb.DMatrix(data = dd_train_sparse , label = (as.integer(X_red_train$Y)-1))
dd_test = xgb.DMatrix(data = dd_test_sparse , label = (as.integer(X_red_test$Y)-1))
watchlist <- list(train=dd_train, test=dd_test)
mod_xgboost <- xgb.train(data=dd_train,num_class=5,max_depth=2, eta=1, nthread = 2, nrounds=10,  eval_metric = "merror", eval_metric = "mlogloss",watchlist=watchlist,objective = "multi:softprob")

predicted.prob <- mod_xgboost %>% predict(dd_train,reshape=TRUE)
head(predicted.prob)
factor(max.col(predicted.prob)-1)
caret::confusionMatrix(factor(as.integer(X_red_train$Y)-1),factor(max.col(predicted.prob)-1))
importance_matrix <- xgb.importance(model = mod_xgboost)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix,top_n=20,rel_to_first = TRUE)

predicted.prob_test <- mod_xgboost %>% predict(dd_test,reshape=TRUE)
head(predicted.prob_test)
factor(max.col(predicted.prob_test)-1)
caret::confusionMatrix(factor(as.integer(X_red_test$Y)-1),factor(max.col(predicted.prob_test)-1))


#Build neural Network
library(neuralnet)
library(nnet)
#Training Data
dd_nnet_train<-dd_train_sparse[,dimnames(dd_train_sparse)[[2]] %in% c("X_Ray1","X_Ray2","Age","BNP","CRP")]
dd_nnet_train<-as.matrix(dd_nnet_train)
dd_nnet_train<-data.frame(dd_nnet_train)
min_max_norm <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
dd_nnet_train<-as.data.frame(lapply(dd_nnet_train,min_max_norm))
dd_nnet_train<-cbind(dd_nnet_train,class.ind(X_red_train$Y))
names(dd_nnet_train)[6]<-"Acute_Asthma"
names(dd_nnet_train)[7]<-"Acute_COPD"
#Test Data
dd_nnet_test<-dd_test_sparse[,dimnames(dd_test_sparse)[[2]] %in% c("X_Ray1","X_Ray2","Age","BNP","CRP")]
dd_nnet_test<-as.matrix(dd_nnet_test)
dd_nnet_test<-data.frame(dd_nnet_test)
dd_nnet_test<-as.data.frame(lapply(dd_nnet_test,min_max_norm))
dd_nnet_test<-cbind(dd_nnet_test,class.ind(X_red_test$Y))
names(dd_nnet_test)[6]<-"Acute_Asthma"
names(dd_nnet_test)[7]<-"Acute_COPD"
#Build Network
n_m <- names(dd_nnet_train) 
f_m <- as.formula(paste("Acute_Asthma+Acute_COPD+Healthy+HF+Pneumonia~", paste(n_m[!n_m %in% c("Acute_Asthma","Acute_COPD","Healthy","HF","Pneumonia")], collapse = " + ")))
f_m
nn <- neuralnet(f_m,
                data = dd_nnet_train,
                hidden = c(0),
                act.fct = "logistic",
                linear.output = FALSE,
                lifesign = "minimal")
nn$result.matrix
plot(nn)
write.table(data.frame(nn$result.matrix),file="Neural_Net_0_layers_Clinical_Vars.csv",sep=",")

#Accuracy on training data
pr.nn <- compute(nn, dd_nnet_train[, 1:5])
prob_train <- pr.nn$net.result
head(prob_train)
y_ob <- max.col(dd_nnet_train[,6:10])
y_pred <- max.col(prob_train)
mean(y_pred == y_ob)
caret::confusionMatrix(as.factor(y_ob),as.factor(y_pred ))

#Accuracy on training data
pr.nn <- compute(nn, dd_nnet_test[, 1:5])
prob_test <- pr.nn$net.result
head(prob_test)
y_ob <- max.col(dd_nnet_test[,6:10])
y_pred <- max.col(prob_test)
mean(y_pred == y_ob)
caret::confusionMatrix(as.factor(y_ob),as.factor(y_pred))
