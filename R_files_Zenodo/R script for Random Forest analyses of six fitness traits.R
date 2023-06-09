# Charlie Waters
# University of Washington School of Aquatic and Fishery Sciences
# Random Forest analysis on CESRF Chinook for weight at Roza

# Genotypes and phenotypes were corrected for all potential covariates

install.packages("randomForest")
library(randomForest)

# Read in data. For all loci, missing values were imputed prior to correction since RF can't take missing values. 
# The 1998 Founders were excluded since it provides the most balanced design between INT and SEG lines. 
# Individuals with missing phenotypic values or covariates needed to be excluded from the analysis, which reduced sample sizes. 

All_data <- read.csv("weight_corrected_all_covariates.csv", header=TRUE)

# Initial RF with all 9108 loci for Weight
rf_weight_250000_1 = randomForest(x = All_data[,3:9110], y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_250000_1,file="rf_weight_250000_1.Rdata")
rf_weight_250000_2 = randomForest(x = All_data[,3:9110], y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_250000_2,file="rf_weight_250000_2.Rdata")

#Check correlation
importance_rf_weight_250000_1<-data.frame(importance(rf_weight_250000_1,type=1))
importance_rf_weight_250000_2<-data.frame(importance(rf_weight_250000_2,type=1))
imp <- cbind(importance_rf_weight_250000_1,importance_rf_weight_250000_2)
colnames(imp)<-c("Imp1","Imp2")
fit<-lm(Imp1~Imp2,data=imp)
summary(fit)

rf_weight_250000_3 = randomForest(x = All_data[,3:9110], y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_250000_3,file="rf_weight_250000_3.Rdata")

#load("rf_weight_250000_1.Rdata")
#load("rf_weight_250000_2.Rdata")
#load("rf_weight_250000_3.Rdata")

mean_rf_weight_250000_1_rsq <- mean(rf_weight_250000_1$rsq)
mean_rf_weight_250000_2_rsq <- mean(rf_weight_250000_2$rsq)
mean_rf_weight_250000_3_rsq <- mean(rf_weight_250000_3$rsq)

#Plot of variable importance for Weight at Roza
varImpPlot(rf_weight_250000_1, main="All Loci weight ntree=250000 1",n.var=50)
varImpPlot(rf_weight_250000_2, main="All Loci weight ntree=250000 2",n.var=50)
varImpPlot(rf_weight_250000_3, main="All Loci weight ntree=250000 3",n.var=50)

#Extract importance (mean decrease in accuracy) measure of each locus for weight
importance_rf_weight_250000_1<-data.frame(importance(rf_weight_250000_1,type=1))
colnames(importance_rf_weight_250000_1)<-c("importance")
importance_rf_weight_250000_2<-data.frame(importance(rf_weight_250000_2,type=1))
colnames(importance_rf_weight_250000_2)<-c("importance")
importance_rf_weight_250000_3<-data.frame(importance(rf_weight_250000_3,type=1))
colnames(importance_rf_weight_250000_3)<-c("importance")

importance_rf_weight_250000_all_loci <-cbind(rownames(importance_rf_weight_250000_1),importance_rf_weight_250000_1,importance_rf_weight_250000_2, importance_rf_weight_250000_3)
colnames(importance_rf_weight_250000_all_loci)<-c("Variable","Importance1","Importance2", "Importance3")
write.csv(importance_rf_weight_250000_all_loci,file="RF_weight_importance_all_loci_250000.csv",row.names=FALSE)

#Plot importance values against each other to check if model converged
plot(importance_rf_weight_250000_all_loci$Importance1,importance_rf_weight_250000_all_loci$Importance2,xlim=c(-15,25),ylim=c(-15,25),main=("Imp Values weight 150k, R2=0.96"))
fit<-lm(Importance2~Importance1, data=importance_rf_weight_250000_all_loci)
summary(fit)
abline(fit)

################################ weight

##### Best 0.5% 
names_best_0.5perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.995))]
names_best_0.5perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.995))]
names_best_0.5perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.995))]
names_best_0.5perc_weight_unique<-unique(c(names_best_0.5perc_weight_1,names_best_0.5perc_weight_2,names_best_0.5perc_weight_3))
genotypes_best0.5perc_weight<-All_data[,colnames(All_data) %in% names_best_0.5perc_weight_unique]
rf_weight_0.5perc_1 = randomForest(x=genotypes_best0.5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_0.5perc_1,file="rf_weight_0.5perc_1.Rdata")
rf_weight_0.5perc_2 = randomForest(x=genotypes_best0.5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_0.5perc_2,file="rf_weight_0.5perc_2.Rdata")
rf_weight_0.5perc_3 = randomForest(x=genotypes_best0.5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_0.5perc_3,file="rf_weight_0.5perc_3.Rdata")

#load("rf_weight_0.5perc_1.Rdata")
#load("rf_weight_0.5perc_2.Rdata")
#load("rf_weight_0.5perc_3.Rdata")

mean_rf_weight_0.5perc_1_rsq <- mean(rf_weight_0.5perc_1$rsq)
mean_rf_weight_0.5perc_2_rsq <- mean(rf_weight_0.5perc_2$rsq)
mean_rf_weight_0.5perc_3_rsq <- mean(rf_weight_0.5perc_3$rsq)

rm(rf_weight_0.5perc_1,rf_weight_0.5perc_2,rf_weight_0.5perc_3)

##### Best 1% 
names_best_1perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.99))]
names_best_1perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.99))]
names_best_1perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.99))]
names_best_1perc_weight_unique<-unique(c(names_best_1perc_weight_1,names_best_1perc_weight_2,names_best_1perc_weight_3))
genotypes_best1perc_weight<-All_data[,colnames(All_data) %in% names_best_1perc_weight_unique]
rf_weight_1perc_1 = randomForest(x=genotypes_best1perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_1perc_1,file="rf_weight_1perc_1.Rdata")
rf_weight_1perc_2 = randomForest(x=genotypes_best1perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_1perc_2,file="rf_weight_1perc_2.Rdata")
rf_weight_1perc_3 = randomForest(x=genotypes_best1perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_1perc_3,file="rf_weight_1perc_3.Rdata")

#load("rf_weight_1perc_1.Rdata")
#load("rf_weight_1perc_2.Rdata")
#load("rf_weight_1perc_3.Rdata")

mean_rf_weight_1perc_1_rsq <- mean(rf_weight_1perc_1$rsq)
mean_rf_weight_1perc_2_rsq <- mean(rf_weight_1perc_2$rsq)
mean_rf_weight_1perc_3_rsq <- mean(rf_weight_1perc_3$rsq)

rm(rf_weight_1perc_1,rf_weight_1perc_2,rf_weight_1perc_3)

##### Best 1.5%
names_best_1.5perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.985))]
names_best_1.5perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.985))]
names_best_1.5perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.985))]
names_best_1.5perc_weight_unique<-unique(c(names_best_1.5perc_weight_1,names_best_1.5perc_weight_2,names_best_1.5perc_weight_3))
genotypes_best1.5perc_weight<-All_data[,colnames(All_data) %in% names_best_1.5perc_weight_unique]
rf_weight_1.5perc_1 = randomForest(x=genotypes_best1.5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_1.5perc_1,file="rf_weight_1.5perc_1.Rdata")
rf_weight_1.5perc_2 = randomForest(x=genotypes_best1.5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_1.5perc_2,file="rf_weight_1.5perc_2.Rdata")
rf_weight_1.5perc_3 = randomForest(x=genotypes_best1.5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_1.5perc_3,file="rf_weight_1.5perc_3.Rdata")

#load("rf_weight_1.5perc_1.Rdata")
#load("rf_weight_1.5perc_2.Rdata")
#load("rf_weight_1.5perc_3.Rdata")

mean_rf_weight_1.5perc_1_rsq <- mean(rf_weight_1.5perc_1$rsq)
mean_rf_weight_1.5perc_2_rsq <- mean(rf_weight_1.5perc_2$rsq)
mean_rf_weight_1.5perc_3_rsq <- mean(rf_weight_1.5perc_3$rsq)

rm(rf_weight_1.5perc_1,rf_weight_1.5perc_2,rf_weight_1.5perc_3)

##### Best 2% 
names_best_2perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.98))]
names_best_2perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.98))]
names_best_2perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.98))]
names_best_2perc_weight_unique<-unique(c(names_best_2perc_weight_1,names_best_2perc_weight_2,names_best_2perc_weight_3))
genotypes_best2perc_weight<-All_data[,colnames(All_data) %in% names_best_2perc_weight_unique]
rf_weight_2perc_1 = randomForest(x=genotypes_best2perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_2perc_1,file="rf_weight_2perc_1.Rdata")
rf_weight_2perc_2 = randomForest(x=genotypes_best2perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_2perc_2,file="rf_weight_2perc_2.Rdata")
rf_weight_2perc_3 = randomForest(x=genotypes_best2perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_2perc_3,file="rf_weight_2perc_3.Rdata")

#load("rf_weight_2perc_1.Rdata")
#load("rf_weight_2perc_2.Rdata")
#load("rf_weight_2perc_3.Rdata")

mean_rf_weight_2perc_1_rsq <- mean(rf_weight_2perc_1$rsq)
mean_rf_weight_2perc_2_rsq <- mean(rf_weight_2perc_2$rsq)
mean_rf_weight_2perc_3_rsq <- mean(rf_weight_2perc_3$rsq)

rm(rf_weight_2perc_1,rf_weight_2perc_2,rf_weight_2perc_3)

##### Best 3% 
names_best_3perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.97))]
names_best_3perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.97))]
names_best_3perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.97))]
names_best_3perc_weight_unique<-unique(c(names_best_3perc_weight_1,names_best_3perc_weight_2,names_best_3perc_weight_3))
genotypes_best3perc_weight<-All_data[,colnames(All_data) %in% names_best_3perc_weight_unique]
rf_weight_3perc_1 = randomForest(x=genotypes_best3perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_3perc_1,file="rf_weight_3perc_1.Rdata")
rf_weight_3perc_2 = randomForest(x=genotypes_best3perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_3perc_2,file="rf_weight_3perc_2.Rdata")
rf_weight_3perc_3 = randomForest(x=genotypes_best3perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_3perc_3,file="rf_weight_3perc_3.Rdata")

#load("rf_weight_3perc_1.Rdata")
#load("rf_weight_3perc_2.Rdata")
#load("rf_weight_3perc_3.Rdata")

mean_rf_weight_3perc_1_rsq <- mean(rf_weight_3perc_1$rsq)
mean_rf_weight_3perc_2_rsq <- mean(rf_weight_3perc_2$rsq)
mean_rf_weight_3perc_3_rsq <- mean(rf_weight_3perc_3$rsq)

rm(rf_weight_3perc_1,rf_weight_3perc_2,rf_weight_3perc_3)

##### Best 4% 
names_best_4perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.96))]
names_best_4perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.96))]
names_best_4perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.96))]
names_best_4perc_weight_unique<-unique(c(names_best_4perc_weight_1,names_best_4perc_weight_2,names_best_4perc_weight_3))
genotypes_best4perc_weight<-All_data[,colnames(All_data) %in% names_best_4perc_weight_unique]
rf_weight_4perc_1 = randomForest(x=genotypes_best4perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_4perc_1,file="rf_weight_4perc_1.Rdata")
rf_weight_4perc_2 = randomForest(x=genotypes_best4perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_4perc_2,file="rf_weight_4perc_2.Rdata")
rf_weight_4perc_3 = randomForest(x=genotypes_best4perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_4perc_3,file="rf_weight_4perc_3.Rdata")

#load("rf_weight_4perc_1.Rdata")
#load("rf_weight_4perc_2.Rdata")
#load("rf_weight_4perc_3.Rdata")

mean_rf_weight_4perc_1_rsq <- mean(rf_weight_4perc_1$rsq)
mean_rf_weight_4perc_2_rsq <- mean(rf_weight_4perc_2$rsq)
mean_rf_weight_4perc_3_rsq <- mean(rf_weight_4perc_3$rsq)

rm(rf_weight_4perc_1,rf_weight_4perc_2,rf_weight_4perc_3)

##### Best 5% 
names_best_5perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.95))]
names_best_5perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.95))]
names_best_5perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.95))]
names_best_5perc_weight_unique<-unique(c(names_best_5perc_weight_1,names_best_5perc_weight_2,names_best_5perc_weight_3))
genotypes_best5perc_weight<-All_data[,colnames(All_data) %in% names_best_5perc_weight_unique]
rf_weight_5perc_1 = randomForest(x=genotypes_best5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_5perc_1,file="rf_weight_5perc_1.Rdata")
rf_weight_5perc_2 = randomForest(x=genotypes_best5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_5perc_2,file="rf_weight_5perc_2.Rdata")
rf_weight_5perc_3 = randomForest(x=genotypes_best5perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_5perc_3,file="rf_weight_5perc_3.Rdata")

#load("rf_weight_5perc_1.Rdata")
#load("rf_weight_5perc_2.Rdata")
#load("rf_weight_5perc_3.Rdata")

mean_rf_weight_5perc_1_rsq <- mean(rf_weight_5perc_1$rsq)
mean_rf_weight_5perc_2_rsq <- mean(rf_weight_5perc_2$rsq)
mean_rf_weight_5perc_3_rsq <- mean(rf_weight_5perc_3$rsq)

rm(rf_weight_5perc_1,rf_weight_5perc_2,rf_weight_5perc_3)

##### Best 10% 
names_best_10perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.90))]
names_best_10perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.90))]
names_best_10perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.90))]
names_best_10perc_weight_unique<-unique(c(names_best_10perc_weight_1,names_best_10perc_weight_2,names_best_10perc_weight_3))
genotypes_best10perc_weight<-All_data[,colnames(All_data) %in% names_best_10perc_weight_unique]
rf_weight_10perc_1 = randomForest(x=genotypes_best10perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_10perc_1,file="rf_weight_10perc_1.Rdata")
rf_weight_10perc_2 = randomForest(x=genotypes_best10perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_10perc_2,file="rf_weight_10perc_2.Rdata")
rf_weight_10perc_3 = randomForest(x=genotypes_best10perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_10perc_3,file="rf_weight_10perc_3.Rdata")

#load("rf_weight_10perc_1.Rdata")
#load("rf_weight_10perc_2.Rdata")
#load("rf_weight_10perc_3.Rdata")

mean_rf_weight_10perc_1_rsq <- mean(rf_weight_10perc_1$rsq)
mean_rf_weight_10perc_2_rsq <- mean(rf_weight_10perc_2$rsq)
mean_rf_weight_10perc_3_rsq <- mean(rf_weight_10perc_3$rsq)

rm(rf_weight_10perc_1,rf_weight_10perc_2,rf_weight_10perc_3)

##### Best 20% 
names_best_20perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.80))]
names_best_20perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.80))]
names_best_20perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.80))]
names_best_20perc_weight_unique<-unique(c(names_best_20perc_weight_1,names_best_20perc_weight_2,names_best_20perc_weight_3))
genotypes_best20perc_weight<-All_data[,colnames(All_data) %in% names_best_20perc_weight_unique]
rf_weight_20perc_1 = randomForest(x=genotypes_best20perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_20perc_1,file="rf_weight_20perc_1.Rdata")
rf_weight_20perc_2 = randomForest(x=genotypes_best20perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_20perc_2,file="rf_weight_20perc_2.Rdata")
rf_weight_20perc_3 = randomForest(x=genotypes_best20perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_20perc_3,file="rf_weight_20perc_3.Rdata")

#load("rf_weight_20perc_1.Rdata")
#load("rf_weight_20perc_2.Rdata")
#load("rf_weight_20perc_3.Rdata")

mean_rf_weight_20perc_1_rsq <- mean(rf_weight_20perc_1$rsq)
mean_rf_weight_20perc_2_rsq <- mean(rf_weight_20perc_2$rsq)
mean_rf_weight_20perc_3_rsq <- mean(rf_weight_20perc_3$rsq)

rm(rf_weight_20perc_1,rf_weight_20perc_2,rf_weight_20perc_3)

##### Best 30% 
names_best_30perc_weight_1<-rownames(importance_rf_weight_250000_1)[which(importance_rf_weight_250000_1$importance > quantile(importance_rf_weight_250000_1$importance, probs=0.70))]
names_best_30perc_weight_2<-rownames(importance_rf_weight_250000_2)[which(importance_rf_weight_250000_2$importance > quantile(importance_rf_weight_250000_2$importance, probs=0.70))]
names_best_30perc_weight_3<-rownames(importance_rf_weight_250000_3)[which(importance_rf_weight_250000_3$importance > quantile(importance_rf_weight_250000_3$importance, probs=0.70))]
names_best_30perc_weight_unique<-unique(c(names_best_30perc_weight_1,names_best_30perc_weight_2,names_best_30perc_weight_3))
genotypes_best30perc_weight<-All_data[,colnames(All_data) %in% names_best_30perc_weight_unique]
rf_weight_30perc_1 = randomForest(x=genotypes_best30perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_30perc_1,file="rf_weight_30perc_1.Rdata")
rf_weight_30perc_2 = randomForest(x=genotypes_best30perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_30perc_2,file="rf_weight_30perc_2.Rdata")
rf_weight_30perc_3 = randomForest(x=genotypes_best30perc_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_weight_30perc_3,file="rf_weight_30perc_3.Rdata")

#load("rf_weight_30perc_1.Rdata")
#load("rf_weight_30perc_2.Rdata")
#load("rf_weight_30perc_3.Rdata")

mean_rf_weight_30perc_1_rsq <- mean(rf_weight_30perc_1$rsq)
mean_rf_weight_30perc_2_rsq <- mean(rf_weight_30perc_2$rsq)
mean_rf_weight_30perc_3_rsq <- mean(rf_weight_30perc_3$rsq)

rm(rf_weight_30perc_1,rf_weight_30perc_2,rf_weight_30perc_3)

## Variance explained-weight
All_initial_Var_explained_weight<-rbind(cbind(mean_rf_weight_250000_1_rsq,mean_rf_weight_250000_2_rsq,mean_rf_weight_250000_3_rsq),
                                        cbind(mean_rf_weight_0.5perc_1_rsq,mean_rf_weight_0.5perc_2_rsq,mean_rf_weight_0.5perc_3_rsq),
                                        cbind(mean_rf_weight_1perc_1_rsq,mean_rf_weight_1perc_2_rsq,mean_rf_weight_1perc_3_rsq),
                                        cbind(mean_rf_weight_1.5perc_1_rsq,mean_rf_weight_1.5perc_2_rsq,mean_rf_weight_1.5perc_3_rsq),
                                        cbind(mean_rf_weight_2perc_1_rsq,mean_rf_weight_2perc_2_rsq,mean_rf_weight_2perc_3_rsq),
                                        cbind(mean_rf_weight_3perc_1_rsq,mean_rf_weight_3perc_2_rsq,mean_rf_weight_3perc_3_rsq),
                                        cbind(mean_rf_weight_4perc_1_rsq,mean_rf_weight_4perc_2_rsq,mean_rf_weight_4perc_3_rsq),
                                        cbind(mean_rf_weight_5perc_1_rsq,mean_rf_weight_5perc_2_rsq,mean_rf_weight_5perc_3_rsq),
                                        cbind(mean_rf_weight_10perc_1_rsq,mean_rf_weight_10perc_2_rsq,mean_rf_weight_10perc_3_rsq),
                                        cbind(mean_rf_weight_20perc_1_rsq,mean_rf_weight_20perc_2_rsq,mean_rf_weight_20perc_3_rsq),
                                        cbind(mean_rf_weight_30perc_1_rsq,mean_rf_weight_30perc_2_rsq,mean_rf_weight_30perc_3_rsq))

#weight
All_initial_Var_explained_weight<-data.frame(All_initial_Var_explained_weight)
All_initial_Var_explained_weight$Number_loci<-c(9108,length(names_best_0.5perc_weight_unique),length(names_best_1perc_weight_unique),length(names_best_1.5perc_weight_unique),length(names_best_2perc_weight_unique),length(names_best_3perc_weight_unique),length(names_best_4perc_weight_unique),length(names_best_5perc_weight_unique),length(names_best_10perc_weight_unique),length(names_best_20perc_weight_unique),length(names_best_30perc_weight_unique))
rownames(All_initial_Var_explained_weight)<-c("All","Best0.5%","Best1%","Best1.5%","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_Var_explained_weight$Average<-apply(All_initial_Var_explained_weight[,1:3],1,mean)
write.csv(All_initial_Var_explained_weight,file="All_initial_Var_explained_weight250k.csv")
par(mar=c(5,5,3,3))
plot(All_initial_Var_explained_weight$Number_loci,All_initial_Var_explained_weight$Average,log="x",
     xlab="Number of Loci",ylab="Proportion Variation Explained",
     main="Initial Variation Explained for Weight at Roza",cex.axis=1.5,cex.lab=1.5,pch=19)


#Based on this table and plot, the best 1% of loci have the largest r-squared value, which approximates to the % variance explained
# The best 1% data set includes 101 unique loci, so I'll run backward purging RF with the best 2% loci (n=207 unique)


#################### Backward purging approach - weight
names_best_weight_unique <- names_best_2perc_weight_unique
length(names_best_weight_unique)

names_best_weight207 <- names_best_weight_unique
genotypes_purging_weight<-All_data[,colnames(All_data) %in% names_best_weight207]
rf_weight_purging_1 = randomForest(x=genotypes_purging_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=100000)
save(rf_weight_purging_1,file="rf_weight_purging_1_100k.Rdata")
rf_weight_purging_2 = randomForest(x=genotypes_purging_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=100000)
save(rf_weight_purging_2,file="rf_weight_purging_2_100k.Rdata")

#Extract importance (mean decrease in accuracy) measure of each locus for weight
importance_rf_weight_purging_1<-data.frame(importance(rf_weight_purging_1,type=1))
colnames(importance_rf_weight_purging_1)<-c("Importance1")
importance_rf_weight_purging_2<-data.frame(importance(rf_weight_purging_2,type=1))
colnames(importance_rf_weight_purging_2)<-c("Importance2")

importance_rf_purging <-cbind(rownames(importance_rf_weight_purging_1),importance_rf_weight_purging_1,importance_rf_weight_purging_2)

#Plot importance values against each other to check if model converged
plot(importance_rf_purging$Importance1,importance_rf_purging$Importance2)
fit<-lm(Importance2~Importance1, data=importance_rf_purging)
summary(fit)
abline(fit)


rf_weight_purging_3 = randomForest(x=genotypes_purging_weight,y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=100000)
save(rf_weight_purging_3,file="rf_weight_purging_3_100k.Rdata")

#load("rf_weight_purging_1.Rdata")
#load("rf_weight_purging_2.Rdata")
#load("rf_weight_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[207]]<-names_best_weight207
all_Var_explained_weight_best<-data.frame(V1=1:207,V2=1:207,V3=1:207)
rownames(all_Var_explained_weight_best)<-1:207
all_Var_explained_weight_best[207,]<-c(mean(rf_weight_purging_1$rsq),mean(rf_weight_purging_2$rsq),mean(rf_weight_purging_3$rsq))
all_mse_weight<-data.frame(V1=1:207,V2=1:207,V3=1:207)
rownames(all_mse_weight)<-1:207
all_mse_weight[207,]<-c(mean(rf_weight_purging_1$mse),mean(rf_weight_purging_2$mse),mean(rf_weight_purging_3$mse))

for (i in 1:206){
  print(i)
  imp_purging_1<-data.frame(importance(rf_weight_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_weight_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_weight_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging_weight)
  colnames(imp_purging_1)<-"X.IncMSE1"
  rownames(imp_purging_2)<-colnames(genotypes_purging_weight)
  colnames(imp_purging_2)<-"X.IncMSE2"
  rownames(imp_purging_3)<-colnames(genotypes_purging_weight)
  colnames(imp_purging_3)<-"X.IncMSE3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[207-i]]<-names_keep
  genotypes_purging_weight<-All_data[,colnames(All_data) %in% names_keep]
  rf_weight_purging_1 = randomForest(x = genotypes_purging_weight, y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=100000)
  rf_weight_purging_2 = randomForest(x = genotypes_purging_weight, y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=100000)
  rf_weight_purging_3 = randomForest(x = genotypes_purging_weight, y = All_data$weight_Roza, importance=TRUE ,proximity=TRUE, ntree=100000)
  all_Var_explained_weight_best[207-i,]<-c(mean(rf_weight_purging_1$rsq),mean(rf_weight_purging_2$rsq),mean(rf_weight_purging_3$rsq))
  all_mse_weight[207-i,]<-c(mean(rf_weight_purging_1$mse),mean(rf_weight_purging_2$mse),mean(rf_weight_purging_3$mse))
  
}

all_Var_explained_weight_best$Average<-apply(all_Var_explained_weight_best,1,mean)

lapply(names_all_iterations, write, "Backward_purging_weight_names_all_iterations.txt", append=TRUE, ncolumns=210)

write.csv(all_Var_explained_weight_best, file="Backward_purging_variance_explained_weight.csv")

par(mar=c(5,5,3,3))
plot(all_Var_explained_weight_best$Average[-c(1)],xlab="Number of Loci", ylab="Proportion Variation Explained",
     main="Backward Purging Weight",cex.lab=1.5,cex.axis=1.5,pch=19)

which(all_Var_explained_weight_best$Average==max(all_Var_explained_weight_best$Average[-c(1)]))

write.csv(names_all_iterations[[37]],file="Purging_optimum_weight_250000.csv")


##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################


# Charlie Waters
# University of Washington School of Aquatic and Fishery Sciences
# Random Forest analysis on CESRF Chinook for daily growth coefficient

# Genotypes and phenotypes were corrected for all potential covariates

install.packages("randomForest")
library(randomForest)

# Read in data. For all loci, missing values were imputed prior to correction since RF can't take missing values. 
# The 1998 Founders were excluded since it provides the most balanced design between INT and SEG lines. 
# Individuals with missing phenotypic values or covariates needed to be excluded from the analysis, which reduced sample sizes. 

All_data <- read.csv("DGC_corrected_all_covariates.csv", header=TRUE)

# Initial RF with all 9108 loci for Daily Growth Coefficient
rf_DGC_250000_1 = randomForest(x = All_data[,3:9110], y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_250000_1,file="rf_DGC_250000_1.Rdata")
rf_DGC_250000_2 = randomForest(x = All_data[,3:9110], y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_250000_2,file="rf_DGC_250000_2.Rdata")

#Check correlation early
importance_rf_DGC_250000_1<-data.frame(importance(rf_DGC_250000_1,type=1))
importance_rf_DGC_250000_2<-data.frame(importance(rf_DGC_250000_2,type=1))
imp<-cbind(importance_rf_DGC_250000_1,importance_rf_DGC_250000_2)
colnames(imp)<-c("Imp1","Imp2")
fit <- lm(Imp1~Imp2,data=imp)
summary(fit)

rf_DGC_250000_3 = randomForest(x = All_data[,3:9110], y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_250000_3,file="rf_DGC_250000_3.Rdata")

#load("rf_DGC_250000_1.Rdata")
#load("rf_DGC_250000_2.Rdata")
#load("rf_DGC_250000_3.Rdata")

mean_rf_DGC_250000_1_rsq <- mean(rf_DGC_250000_1$rsq)
mean_rf_DGC_250000_2_rsq <- mean(rf_DGC_250000_2$rsq)
mean_rf_DGC_250000_3_rsq <- mean(rf_DGC_250000_3$rsq)

#Plot of variable importance for CESRF DGC
varImpPlot(rf_DGC_250000_1, main="All Loci DGC ntree=250000 1",n.var=50)
varImpPlot(rf_DGC_250000_2, main="All Loci DGC ntree=250000 2",n.var=50)
varImpPlot(rf_DGC_250000_3, main="All Loci DGC ntree=250000 3",n.var=50)

#Extract importance (mean decrease in accuracy) measure of each locus for DGC
importance_rf_DGC_250000_1<-data.frame(importance(rf_DGC_250000_1,type=1))
colnames(importance_rf_DGC_250000_1)<-c("importance")
importance_rf_DGC_250000_2<-data.frame(importance(rf_DGC_250000_2,type=1))
colnames(importance_rf_DGC_250000_2)<-c("importance")
importance_rf_DGC_250000_3<-data.frame(importance(rf_DGC_250000_3,type=1))
colnames(importance_rf_DGC_250000_3)<-c("importance")

importance_rf_DGC_250000_all_loci <-cbind(rownames(importance_rf_DGC_250000_1),importance_rf_DGC_250000_1,importance_rf_DGC_250000_2, importance_rf_DGC_250000_3)
colnames(importance_rf_DGC_250000_all_loci)<-c("Variable","Importance1","Importance2", "Importance3")
write.csv(importance_rf_DGC_250000_all_loci,file="RF_DGC_importance_all_loci_250000.csv",row.names=FALSE)

#Plot importance values against each other to check if model converged
plot(importance_rf_DGC_250000_all_loci$Importance1,importance_rf_DGC_250000_all_loci$Importance2,xlim=c(-15,25),ylim=c(-15,25),main=("Imp Values DGC 300k"))
fit<-lm(Importance2~Importance1, data=importance_rf_DGC_250000_all_loci)
summary(fit)
abline(fit)

################################ CESRF DGC

##### Best 0.5% 
names_best_0.5perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.995))]
names_best_0.5perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.995))]
names_best_0.5perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.995))]
names_best_0.5perc_DGC_unique<-unique(c(names_best_0.5perc_DGC_1,names_best_0.5perc_DGC_2,names_best_0.5perc_DGC_3))
genotypes_best0.5perc_DGC<-All_data[,colnames(All_data) %in% names_best_0.5perc_DGC_unique]
rf_DGC_0.5perc_1 = randomForest(x=genotypes_best0.5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_0.5perc_1,file="rf_DGC_0.5perc_1.Rdata")
rf_DGC_0.5perc_2 = randomForest(x=genotypes_best0.5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_0.5perc_2,file="rf_DGC_0.5perc_2.Rdata")
rf_DGC_0.5perc_3 = randomForest(x=genotypes_best0.5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_0.5perc_3,file="rf_DGC_0.5perc_3.Rdata")

#load("rf_DGC_0.5perc_1.Rdata")
#load("rf_DGC_0.5perc_2.Rdata")
#load("rf_DGC_0.5perc_3.Rdata")

mean_rf_DGC_0.5perc_1_rsq <- mean(rf_DGC_0.5perc_1$rsq)
mean_rf_DGC_0.5perc_2_rsq <- mean(rf_DGC_0.5perc_2$rsq)
mean_rf_DGC_0.5perc_3_rsq <- mean(rf_DGC_0.5perc_3$rsq)

rm(rf_DGC_0.5perc_1,rf_DGC_0.5perc_2,rf_DGC_0.5perc_3)

##### Best 1% 
names_best_1perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.99))]
names_best_1perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.99))]
names_best_1perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.99))]
names_best_1perc_DGC_unique<-unique(c(names_best_1perc_DGC_1,names_best_1perc_DGC_2,names_best_1perc_DGC_3))
genotypes_best1perc_DGC<-All_data[,colnames(All_data) %in% names_best_1perc_DGC_unique]
rf_DGC_1perc_1 = randomForest(x=genotypes_best1perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_1perc_1,file="rf_DGC_1perc_1.Rdata")
rf_DGC_1perc_2 = randomForest(x=genotypes_best1perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_1perc_2,file="rf_DGC_1perc_2.Rdata")
rf_DGC_1perc_3 = randomForest(x=genotypes_best1perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_1perc_3,file="rf_DGC_1perc_3.Rdata")

#load("rf_DGC_1perc_1.Rdata")
#load("rf_DGC_1perc_2.Rdata")
#load("rf_DGC_1perc_3.Rdata")

mean_rf_DGC_1perc_1_rsq <- mean(rf_DGC_1perc_1$rsq)
mean_rf_DGC_1perc_2_rsq <- mean(rf_DGC_1perc_2$rsq)
mean_rf_DGC_1perc_3_rsq <- mean(rf_DGC_1perc_3$rsq)

rm(rf_DGC_1perc_1,rf_DGC_1perc_2,rf_DGC_1perc_3)

##### Best 1.5%
names_best_1.5perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.985))]
names_best_1.5perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.985))]
names_best_1.5perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.985))]
names_best_1.5perc_DGC_unique<-unique(c(names_best_1.5perc_DGC_1,names_best_1.5perc_DGC_2,names_best_1.5perc_DGC_3))
genotypes_best1.5perc_DGC<-All_data[,colnames(All_data) %in% names_best_1.5perc_DGC_unique]
rf_DGC_1.5perc_1 = randomForest(x=genotypes_best1.5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_1.5perc_1,file="rf_DGC_1.5perc_1.Rdata")
rf_DGC_1.5perc_2 = randomForest(x=genotypes_best1.5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_1.5perc_2,file="rf_DGC_1.5perc_2.Rdata")
rf_DGC_1.5perc_3 = randomForest(x=genotypes_best1.5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_1.5perc_3,file="rf_DGC_1.5perc_3.Rdata")

#load("rf_DGC_1.5perc_1.Rdata")
#load("rf_DGC_1.5perc_2.Rdata")
#load("rf_DGC_1.5perc_3.Rdata")

mean_rf_DGC_1.5perc_1_rsq <- mean(rf_DGC_1.5perc_1$rsq)
mean_rf_DGC_1.5perc_2_rsq <- mean(rf_DGC_1.5perc_2$rsq)
mean_rf_DGC_1.5perc_3_rsq <- mean(rf_DGC_1.5perc_3$rsq)

rm(rf_DGC_1.5perc_1,rf_DGC_1.5perc_2,rf_DGC_1.5perc_3)

##### Best 2% 
names_best_2perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.98))]
names_best_2perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.98))]
names_best_2perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.98))]
names_best_2perc_DGC_unique<-unique(c(names_best_2perc_DGC_1,names_best_2perc_DGC_2,names_best_2perc_DGC_3))
genotypes_best2perc_DGC<-All_data[,colnames(All_data) %in% names_best_2perc_DGC_unique]
rf_DGC_2perc_1 = randomForest(x=genotypes_best2perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_2perc_1,file="rf_DGC_2perc_1.Rdata")
rf_DGC_2perc_2 = randomForest(x=genotypes_best2perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_2perc_2,file="rf_DGC_2perc_2.Rdata")
rf_DGC_2perc_3 = randomForest(x=genotypes_best2perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_2perc_3,file="rf_DGC_2perc_3.Rdata")

#load("rf_DGC_2perc_1.Rdata")
#load("rf_DGC_2perc_2.Rdata")
#load("rf_DGC_2perc_3.Rdata")

mean_rf_DGC_2perc_1_rsq <- mean(rf_DGC_2perc_1$rsq)
mean_rf_DGC_2perc_2_rsq <- mean(rf_DGC_2perc_2$rsq)
mean_rf_DGC_2perc_3_rsq <- mean(rf_DGC_2perc_3$rsq)

rm(rf_DGC_2perc_1,rf_DGC_2perc_2,rf_DGC_2perc_3)

##### Best 3% 
names_best_3perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.97))]
names_best_3perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.97))]
names_best_3perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.97))]
names_best_3perc_DGC_unique<-unique(c(names_best_3perc_DGC_1,names_best_3perc_DGC_2,names_best_3perc_DGC_3))
genotypes_best3perc_DGC<-All_data[,colnames(All_data) %in% names_best_3perc_DGC_unique]
rf_DGC_3perc_1 = randomForest(x=genotypes_best3perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_3perc_1,file="rf_DGC_3perc_1.Rdata")
rf_DGC_3perc_2 = randomForest(x=genotypes_best3perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_3perc_2,file="rf_DGC_3perc_2.Rdata")
rf_DGC_3perc_3 = randomForest(x=genotypes_best3perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_3perc_3,file="rf_DGC_3perc_3.Rdata")

#load("rf_DGC_3perc_1.Rdata")
#load("rf_DGC_3perc_2.Rdata")
#load("rf_DGC_3perc_3.Rdata")

mean_rf_DGC_3perc_1_rsq <- mean(rf_DGC_3perc_1$rsq)
mean_rf_DGC_3perc_2_rsq <- mean(rf_DGC_3perc_2$rsq)
mean_rf_DGC_3perc_3_rsq <- mean(rf_DGC_3perc_3$rsq)

rm(rf_DGC_3perc_1,rf_DGC_3perc_2,rf_DGC_3perc_3)

##### Best 4% 
names_best_4perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.96))]
names_best_4perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.96))]
names_best_4perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.96))]
names_best_4perc_DGC_unique<-unique(c(names_best_4perc_DGC_1,names_best_4perc_DGC_2,names_best_4perc_DGC_3))
genotypes_best4perc_DGC<-All_data[,colnames(All_data) %in% names_best_4perc_DGC_unique]
rf_DGC_4perc_1 = randomForest(x=genotypes_best4perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_4perc_1,file="rf_DGC_4perc_1.Rdata")
rf_DGC_4perc_2 = randomForest(x=genotypes_best4perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_4perc_2,file="rf_DGC_4perc_2.Rdata")
rf_DGC_4perc_3 = randomForest(x=genotypes_best4perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_4perc_3,file="rf_DGC_4perc_3.Rdata")

#load("rf_DGC_4perc_1.Rdata")
#load("rf_DGC_4perc_2.Rdata")
#load("rf_DGC_4perc_3.Rdata")

mean_rf_DGC_4perc_1_rsq <- mean(rf_DGC_4perc_1$rsq)
mean_rf_DGC_4perc_2_rsq <- mean(rf_DGC_4perc_2$rsq)
mean_rf_DGC_4perc_3_rsq <- mean(rf_DGC_4perc_3$rsq)

rm(rf_DGC_4perc_1,rf_DGC_4perc_2,rf_DGC_4perc_3)

##### Best 5% 
names_best_5perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.95))]
names_best_5perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.95))]
names_best_5perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.95))]
names_best_5perc_DGC_unique<-unique(c(names_best_5perc_DGC_1,names_best_5perc_DGC_2,names_best_5perc_DGC_3))
genotypes_best5perc_DGC<-All_data[,colnames(All_data) %in% names_best_5perc_DGC_unique]
rf_DGC_5perc_1 = randomForest(x=genotypes_best5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_5perc_1,file="rf_DGC_5perc_1.Rdata")
rf_DGC_5perc_2 = randomForest(x=genotypes_best5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_5perc_2,file="rf_DGC_5perc_2.Rdata")
rf_DGC_5perc_3 = randomForest(x=genotypes_best5perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_5perc_3,file="rf_DGC_5perc_3.Rdata")

#load("rf_DGC_5perc_1.Rdata")
#load("rf_DGC_5perc_2.Rdata")
#load("rf_DGC_5perc_3.Rdata")

mean_rf_DGC_5perc_1_rsq <- mean(rf_DGC_5perc_1$rsq)
mean_rf_DGC_5perc_2_rsq <- mean(rf_DGC_5perc_2$rsq)
mean_rf_DGC_5perc_3_rsq <- mean(rf_DGC_5perc_3$rsq)

rm(rf_DGC_5perc_1,rf_DGC_5perc_2,rf_DGC_5perc_3)

##### Best 10% 
names_best_10perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.90))]
names_best_10perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.90))]
names_best_10perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.90))]
names_best_10perc_DGC_unique<-unique(c(names_best_10perc_DGC_1,names_best_10perc_DGC_2,names_best_10perc_DGC_3))
genotypes_best10perc_DGC<-All_data[,colnames(All_data) %in% names_best_10perc_DGC_unique]
rf_DGC_10perc_1 = randomForest(x=genotypes_best10perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_10perc_1,file="rf_DGC_10perc_1.Rdata")
rf_DGC_10perc_2 = randomForest(x=genotypes_best10perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_10perc_2,file="rf_DGC_10perc_2.Rdata")
rf_DGC_10perc_3 = randomForest(x=genotypes_best10perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_10perc_3,file="rf_DGC_10perc_3.Rdata")

#load("rf_DGC_10perc_1.Rdata")
#load("rf_DGC_10perc_2.Rdata")
#load("rf_DGC_10perc_3.Rdata")

mean_rf_DGC_10perc_1_rsq <- mean(rf_DGC_10perc_1$rsq)
mean_rf_DGC_10perc_2_rsq <- mean(rf_DGC_10perc_2$rsq)
mean_rf_DGC_10perc_3_rsq <- mean(rf_DGC_10perc_3$rsq)

rm(rf_DGC_10perc_1,rf_DGC_10perc_2,rf_DGC_10perc_3)

##### Best 20% 
names_best_20perc_DGC_1<-rownames(importance_rf_DGC_250000_1)[which(importance_rf_DGC_250000_1$importance > quantile(importance_rf_DGC_250000_1$importance, probs=0.80))]
names_best_20perc_DGC_2<-rownames(importance_rf_DGC_250000_2)[which(importance_rf_DGC_250000_2$importance > quantile(importance_rf_DGC_250000_2$importance, probs=0.80))]
names_best_20perc_DGC_3<-rownames(importance_rf_DGC_250000_3)[which(importance_rf_DGC_250000_3$importance > quantile(importance_rf_DGC_250000_3$importance, probs=0.80))]
names_best_20perc_DGC_unique<-unique(c(names_best_20perc_DGC_1,names_best_20perc_DGC_2,names_best_20perc_DGC_3))
genotypes_best20perc_DGC<-All_data[,colnames(All_data) %in% names_best_20perc_DGC_unique]
rf_DGC_20perc_1 = randomForest(x=genotypes_best20perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_20perc_1,file="rf_DGC_20perc_1.Rdata")
rf_DGC_20perc_2 = randomForest(x=genotypes_best20perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_20perc_2,file="rf_DGC_20perc_2.Rdata")
rf_DGC_20perc_3 = randomForest(x=genotypes_best20perc_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_DGC_20perc_3,file="rf_DGC_20perc_3.Rdata")

#load("rf_DGC_20perc_1.Rdata")
#load("rf_DGC_20perc_2.Rdata")
#load("rf_DGC_20perc_3.Rdata")

mean_rf_DGC_20perc_1_rsq <- mean(rf_DGC_20perc_1$rsq)
mean_rf_DGC_20perc_2_rsq <- mean(rf_DGC_20perc_2$rsq)
mean_rf_DGC_20perc_3_rsq <- mean(rf_DGC_20perc_3$rsq)

rm(rf_DGC_20perc_1,rf_DGC_20perc_2,rf_DGC_20perc_3)

## Variance explained-DGC
All_initial_Var_explained_DGC<-rbind(cbind(mean_rf_DGC_250000_1_rsq,mean_rf_DGC_250000_2_rsq,mean_rf_DGC_250000_3_rsq),
                                     cbind(mean_rf_DGC_0.5perc_1_rsq,mean_rf_DGC_0.5perc_2_rsq,mean_rf_DGC_0.5perc_3_rsq),
                                     cbind(mean_rf_DGC_1perc_1_rsq,mean_rf_DGC_1perc_2_rsq,mean_rf_DGC_1perc_3_rsq),
                                     cbind(mean_rf_DGC_1.5perc_1_rsq,mean_rf_DGC_1.5perc_2_rsq,mean_rf_DGC_1.5perc_3_rsq),
                                     cbind(mean_rf_DGC_2perc_1_rsq,mean_rf_DGC_2perc_2_rsq,mean_rf_DGC_2perc_3_rsq),
                                     cbind(mean_rf_DGC_3perc_1_rsq,mean_rf_DGC_3perc_2_rsq,mean_rf_DGC_3perc_3_rsq),
                                     cbind(mean_rf_DGC_4perc_1_rsq,mean_rf_DGC_4perc_2_rsq,mean_rf_DGC_4perc_3_rsq),
                                     cbind(mean_rf_DGC_5perc_1_rsq,mean_rf_DGC_5perc_2_rsq,mean_rf_DGC_5perc_3_rsq),
                                     cbind(mean_rf_DGC_10perc_1_rsq,mean_rf_DGC_10perc_2_rsq,mean_rf_DGC_10perc_3_rsq),
                                     cbind(mean_rf_DGC_20perc_1_rsq,mean_rf_DGC_20perc_2_rsq,mean_rf_DGC_20perc_3_rsq))

#DGC
All_initial_Var_explained_DGC<-data.frame(All_initial_Var_explained_DGC)
All_initial_Var_explained_DGC$Number_loci<-c(9108,length(names_best_0.5perc_DGC_unique),length(names_best_1perc_DGC_unique),length(names_best_1.5perc_DGC_unique),length(names_best_2perc_DGC_unique),length(names_best_3perc_DGC_unique),length(names_best_4perc_DGC_unique),length(names_best_5perc_DGC_unique),length(names_best_10perc_DGC_unique),length(names_best_20perc_DGC_unique))
rownames(All_initial_Var_explained_DGC)<-c("All","Best0.5%","Best1%","Best1.5%","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%")
All_initial_Var_explained_DGC$Average<-apply(All_initial_Var_explained_DGC[,1:3],1,mean)
par(mar=c(5,5,3,3))
plot(All_initial_Var_explained_DGC$Number_loci,All_initial_Var_explained_DGC$Average,log="x",
     main="Proportion Variation Explained DGC-RF 250000 Trees",pch=19,xlab="Number of Loci",
     ylab="Proportion of Variation",cex.lab=1.5,cex.axis=1.5)
write.csv(All_initial_Var_explained_DGC,file="Initial_Variance_Explained_DGC_250k.csv")

#Based on this table and plot, the best 0.5% of loci have the largest r-squared value, which approximates to the % variance explained
# The best 0.5% data set includes 51 unique loci, so I'll run backward purging RF with the best 1% loci (n=103 unique)


#################### Backward purging approach - DGC
names_best_DGC_unique <- names_best_1perc_DGC_unique
length(names_best_DGC_unique)

names_best_DGC103 <- names_best_DGC_unique
genotypes_purging_DGC<-All_data[,colnames(All_data) %in% names_best_DGC103]
rf_DGC_purging_1 = randomForest(x=genotypes_purging_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_DGC_purging_1,file="rf_DGC_purging_1_50k.Rdata")
rf_DGC_purging_2 = randomForest(x=genotypes_purging_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_DGC_purging_2,file="rf_DGC_purging_2_50k.Rdata")

#Extract importance (mean decrease in accuracy) measure of each locus for DGC
importance_rf_DGC_purging_1<-data.frame(importance(rf_DGC_purging_1,type=1))
colnames(importance_rf_DGC_purging_1)<-c("Importance1")
importance_rf_DGC_purging_2<-data.frame(importance(rf_DGC_purging_2,type=1))
colnames(importance_rf_DGC_purging_2)<-c("Importance2")

importance_rf_purging <-cbind(rownames(importance_rf_DGC_purging_1),importance_rf_DGC_purging_1,importance_rf_DGC_purging_2)

#Plot importance values against each other to check if model converged
plot(importance_rf_purging$Importance1,importance_rf_purging$Importance2)
fit<-lm(Importance2~Importance1, data=importance_rf_purging)
summary(fit)
abline(fit)

rf_DGC_purging_3 = randomForest(x=genotypes_purging_DGC,y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_DGC_purging_3,file="rf_DGC_purging_3_50k.Rdata")

#load("rf_DGC_purging_1.Rdata")
#load("rf_DGC_purging_2.Rdata")
#load("rf_DGC_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[103]]<-names_best_DGC103
all_Var_explained_DGC_best<-data.frame(V1=1:103,V2=1:103,V3=1:103)
rownames(all_Var_explained_DGC_best)<-1:103
all_Var_explained_DGC_best[103,]<-c(mean(rf_DGC_purging_1$rsq),mean(rf_DGC_purging_2$rsq),mean(rf_DGC_purging_3$rsq))
all_mse_DGC<-data.frame(V1=1:103,V2=1:103,V3=1:103)
rownames(all_mse_DGC)<-1:103
all_mse_DGC[103,]<-c(mean(rf_DGC_purging_1$mse),mean(rf_DGC_purging_2$mse),mean(rf_DGC_purging_3$mse))

for (i in 1:102){
  print(i)
  imp_purging_1<-data.frame(importance(rf_DGC_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_DGC_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_DGC_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging_DGC)
  colnames(imp_purging_1)<-"X.IncMSE1"
  rownames(imp_purging_2)<-colnames(genotypes_purging_DGC)
  colnames(imp_purging_2)<-"X.IncMSE2"
  rownames(imp_purging_3)<-colnames(genotypes_purging_DGC)
  colnames(imp_purging_3)<-"X.IncMSE3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[103-i]]<-names_keep
  genotypes_purging_DGC<-All_data[,colnames(All_data) %in% names_keep]
  rf_DGC_purging_1 = randomForest(x = genotypes_purging_DGC, y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=50000)
  rf_DGC_purging_2 = randomForest(x = genotypes_purging_DGC, y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=50000)
  rf_DGC_purging_3 = randomForest(x = genotypes_purging_DGC, y = All_data$Roza_CESRF_DGC, importance=TRUE ,proximity=TRUE, ntree=50000)
  all_Var_explained_DGC_best[103-i,]<-c(mean(rf_DGC_purging_1$rsq),mean(rf_DGC_purging_2$rsq),mean(rf_DGC_purging_3$rsq))
  all_mse_DGC[103-i,]<-c(mean(rf_DGC_purging_1$mse),mean(rf_DGC_purging_2$mse),mean(rf_DGC_purging_3$mse))
  
  write.csv(all_Var_explained_DGC_best, file="backward_purging_DGC_temp.csv")
  lapply(names_all_iterations, write, "Backward_purging_DGC_names_all_iterations_temp.txt", append=TRUE, ncolumns=110)
}

all_Var_explained_DGC_best$Average<-apply(all_Var_explained_DGC_best,1,mean)

write.csv(all_Var_explained_DGC_best, file="Backward_purging_variance_explained_DGC.csv")

lapply(names_all_iterations, write, "Backward_purging_DGC_names_all_iterations.txt", append=TRUE, ncolumns=110)

par(mar=c(5,5,3,3))
plot(all_Var_explained_DGC_best$Average[-c(1)],xlab="Number of Loci", ylab="Proportion Variation Explained",
     main="Backward Purging DGC Outliers Removed",cex.lab=1.5,cex.axis=1.5,pch=19)

#Find optimum number of loci
which(all_Var_explained_DGC_best$Average==max(all_Var_explained_DGC_best$Average[-c(1)]))

write.csv(names_all_iterations[[35]],file="Purging_optimum_DGC_outliers_removed.csv")



##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################


# Charlie Waters
# University of Washington School of Aquatic and Fishery Sciences
# Random Forest analysis on CESRF Chinook for fork length at Roza Dam

# Genotypes and phenotypes were corrected for all potential covariates

install.packages("randomForest")
library(randomForest)

# Read in data. For all loci, missing values were imputed prior to correction since RF can't take missing values. 
# The 1998 Founders were excluded since it provides the most balanced design between INT and SEG lines. 
# Individuals with missing phenotypic values or covariates needed to be excluded from the analysis, which reduced sample sizes. 

All_data <- read.csv("Forklength_corrected_all_covariates.csv", header=TRUE)

# Initial RF with all 9108 loci for Forklength
rf_forklength_250000_1 = randomForest(x = All_data[,3:9110], y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_250000_1,file="rf_forklength_250000_1.Rdata")
rf_forklength_250000_2 = randomForest(x = All_data[,3:9110], y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_250000_2,file="rf_forklength_250000_2.Rdata")

#Check correlation early
importance_rf_forklength_250000_1<-data.frame(importance(rf_forklength_250000_1,type=1))
importance_rf_forklength_250000_2<-data.frame(importance(rf_forklength_250000_2,type=1))
imp<-cbind(importance_rf_forklength_250000_1,importance_rf_forklength_250000_2)
colnames(imp)<-c("Imp1","Imp2")
fit <- lm(Imp1~Imp2,data=imp)
summary(fit)

rf_forklength_250000_3 = randomForest(x = All_data[,3:9110], y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_250000_3,file="rf_forklength_250000_3.Rdata")

#load("rf_forklength_250000_1.Rdata")
#load("rf_forklength_250000_2.Rdata")
#load("rf_forklength_250000_3.Rdata")

mean_rf_forklength_250000_1_rsq <- mean(rf_forklength_250000_1$rsq)
mean_rf_forklength_250000_2_rsq <- mean(rf_forklength_250000_2$rsq)
mean_rf_forklength_250000_3_rsq <- mean(rf_forklength_250000_3$rsq)

#Plot of variable importance for Roza forklength
varImpPlot(rf_forklength_250000_1, main="All Loci forklength ntree=250000 1",n.var=50)
varImpPlot(rf_forklength_250000_2, main="All Loci forklength ntree=250000 2",n.var=50)
varImpPlot(rf_forklength_250000_3, main="All Loci forklength ntree=250000 3",n.var=50)

#Extract importance (mean decrease in accuracy) measure of each locus for forklength
importance_rf_forklength_250000_1<-data.frame(importance(rf_forklength_250000_1,type=1))
colnames(importance_rf_forklength_250000_1)<-c("importance")
importance_rf_forklength_250000_2<-data.frame(importance(rf_forklength_250000_2,type=1))
colnames(importance_rf_forklength_250000_2)<-c("importance")
importance_rf_forklength_250000_3<-data.frame(importance(rf_forklength_250000_3,type=1))
colnames(importance_rf_forklength_250000_3)<-c("importance")

importance_rf_forklength_250000_all_loci <-cbind(rownames(importance_rf_forklength_250000_1),importance_rf_forklength_250000_1,importance_rf_forklength_250000_2, importance_rf_forklength_250000_3)
colnames(importance_rf_forklength_250000_all_loci)<-c("Variable","Importance1","Importance2", "Importance3")
write.csv(importance_rf_forklength_250000_all_loci,file="RF_forklength_importance_all_loci_250000.csv",row.names=FALSE)

#Plot importance values against each other to check if model converged
plot(importance_rf_forklength_250000_all_loci$Importance1,importance_rf_forklength_250000_all_loci$Importance2,xlim=c(-15,25),ylim=c(-15,25),main=("Imp Values forklength 250k, R2=0.96"))
fit<-lm(Importance2~Importance1, data=importance_rf_forklength_250000_all_loci)
summary(fit)
abline(fit)

################################ CESRF forklength

##### Best 0.5% 
names_best_0.5perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.995))]
names_best_0.5perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.995))]
names_best_0.5perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.995))]
names_best_0.5perc_forklength_unique<-unique(c(names_best_0.5perc_forklength_1,names_best_0.5perc_forklength_2,names_best_0.5perc_forklength_3))
genotypes_best0.5perc_forklength<-All_data[,colnames(All_data) %in% names_best_0.5perc_forklength_unique]
rf_forklength_0.5perc_1 = randomForest(x=genotypes_best0.5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_0.5perc_1,file="rf_forklength_0.5perc_1.Rdata")
rf_forklength_0.5perc_2 = randomForest(x=genotypes_best0.5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_0.5perc_2,file="rf_forklength_0.5perc_2.Rdata")
rf_forklength_0.5perc_3 = randomForest(x=genotypes_best0.5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_0.5perc_3,file="rf_forklength_0.5perc_3.Rdata")

#load("rf_forklength_0.5perc_1.Rdata")
#load("rf_forklength_0.5perc_2.Rdata")
#load("rf_forklength_0.5perc_3.Rdata")

mean_rf_forklength_0.5perc_1_rsq <- mean(rf_forklength_0.5perc_1$rsq)
mean_rf_forklength_0.5perc_2_rsq <- mean(rf_forklength_0.5perc_2$rsq)
mean_rf_forklength_0.5perc_3_rsq <- mean(rf_forklength_0.5perc_3$rsq)

rm(rf_forklength_0.5perc_1,rf_forklength_0.5perc_2,rf_forklength_0.5perc_3)

##### Best 1% 
names_best_1perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.99))]
names_best_1perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.99))]
names_best_1perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.99))]
names_best_1perc_forklength_unique<-unique(c(names_best_1perc_forklength_1,names_best_1perc_forklength_2,names_best_1perc_forklength_3))
genotypes_best1perc_forklength<-All_data[,colnames(All_data) %in% names_best_1perc_forklength_unique]
rf_forklength_1perc_1 = randomForest(x=genotypes_best1perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_1perc_1,file="rf_forklength_1perc_1.Rdata")
rf_forklength_1perc_2 = randomForest(x=genotypes_best1perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_1perc_2,file="rf_forklength_1perc_2.Rdata")
rf_forklength_1perc_3 = randomForest(x=genotypes_best1perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_1perc_3,file="rf_forklength_1perc_3.Rdata")

#load("rf_forklength_1perc_1.Rdata")
#load("rf_forklength_1perc_2.Rdata")
#load("rf_forklength_1perc_3.Rdata")

mean_rf_forklength_1perc_1_rsq <- mean(rf_forklength_1perc_1$rsq)
mean_rf_forklength_1perc_2_rsq <- mean(rf_forklength_1perc_2$rsq)
mean_rf_forklength_1perc_3_rsq <- mean(rf_forklength_1perc_3$rsq)

rm(rf_forklength_1perc_1,rf_forklength_1perc_2,rf_forklength_1perc_3)

##### Best 1.5%
names_best_1.5perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.985))]
names_best_1.5perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.985))]
names_best_1.5perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.985))]
names_best_1.5perc_forklength_unique<-unique(c(names_best_1.5perc_forklength_1,names_best_1.5perc_forklength_2,names_best_1.5perc_forklength_3))
genotypes_best1.5perc_forklength<-All_data[,colnames(All_data) %in% names_best_1.5perc_forklength_unique]
rf_forklength_1.5perc_1 = randomForest(x=genotypes_best1.5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_1.5perc_1,file="rf_forklength_1.5perc_1.Rdata")
rf_forklength_1.5perc_2 = randomForest(x=genotypes_best1.5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_1.5perc_2,file="rf_forklength_1.5perc_2.Rdata")
rf_forklength_1.5perc_3 = randomForest(x=genotypes_best1.5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_1.5perc_3,file="rf_forklength_1.5perc_3.Rdata")

#load("rf_forklength_1.5perc_1.Rdata")
#load("rf_forklength_1.5perc_2.Rdata")
#load("rf_forklength_1.5perc_3.Rdata")

mean_rf_forklength_1.5perc_1_rsq <- mean(rf_forklength_1.5perc_1$rsq)
mean_rf_forklength_1.5perc_2_rsq <- mean(rf_forklength_1.5perc_2$rsq)
mean_rf_forklength_1.5perc_3_rsq <- mean(rf_forklength_1.5perc_3$rsq)

rm(rf_forklength_1.5perc_1,rf_forklength_1.5perc_2,rf_forklength_1.5perc_3)

##### Best 2% 
names_best_2perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.98))]
names_best_2perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.98))]
names_best_2perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.98))]
names_best_2perc_forklength_unique<-unique(c(names_best_2perc_forklength_1,names_best_2perc_forklength_2,names_best_2perc_forklength_3))
genotypes_best2perc_forklength<-All_data[,colnames(All_data) %in% names_best_2perc_forklength_unique]
rf_forklength_2perc_1 = randomForest(x=genotypes_best2perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_2perc_1,file="rf_forklength_2perc_1.Rdata")
rf_forklength_2perc_2 = randomForest(x=genotypes_best2perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_2perc_2,file="rf_forklength_2perc_2.Rdata")
rf_forklength_2perc_3 = randomForest(x=genotypes_best2perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_2perc_3,file="rf_forklength_2perc_3.Rdata")

#load("rf_forklength_2perc_1.Rdata")
#load("rf_forklength_2perc_2.Rdata")
#load("rf_forklength_2perc_3.Rdata")

mean_rf_forklength_2perc_1_rsq <- mean(rf_forklength_2perc_1$rsq)
mean_rf_forklength_2perc_2_rsq <- mean(rf_forklength_2perc_2$rsq)
mean_rf_forklength_2perc_3_rsq <- mean(rf_forklength_2perc_3$rsq)

rm(rf_forklength_2perc_1,rf_forklength_2perc_2,rf_forklength_2perc_3)

##### Best 3% 
names_best_3perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.97))]
names_best_3perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.97))]
names_best_3perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.97))]
names_best_3perc_forklength_unique<-unique(c(names_best_3perc_forklength_1,names_best_3perc_forklength_2,names_best_3perc_forklength_3))
genotypes_best3perc_forklength<-All_data[,colnames(All_data) %in% names_best_3perc_forklength_unique]
rf_forklength_3perc_1 = randomForest(x=genotypes_best3perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_3perc_1,file="rf_forklength_3perc_1.Rdata")
rf_forklength_3perc_2 = randomForest(x=genotypes_best3perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_3perc_2,file="rf_forklength_3perc_2.Rdata")
rf_forklength_3perc_3 = randomForest(x=genotypes_best3perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_3perc_3,file="rf_forklength_3perc_3.Rdata")

#load("rf_forklength_3perc_1.Rdata")
#load("rf_forklength_3perc_2.Rdata")
#load("rf_forklength_3perc_3.Rdata")

mean_rf_forklength_3perc_1_rsq <- mean(rf_forklength_3perc_1$rsq)
mean_rf_forklength_3perc_2_rsq <- mean(rf_forklength_3perc_2$rsq)
mean_rf_forklength_3perc_3_rsq <- mean(rf_forklength_3perc_3$rsq)

rm(rf_forklength_3perc_1,rf_forklength_3perc_2,rf_forklength_3perc_3)

##### Best 4% 
names_best_4perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.96))]
names_best_4perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.96))]
names_best_4perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.96))]
names_best_4perc_forklength_unique<-unique(c(names_best_4perc_forklength_1,names_best_4perc_forklength_2,names_best_4perc_forklength_3))
genotypes_best4perc_forklength<-All_data[,colnames(All_data) %in% names_best_4perc_forklength_unique]
rf_forklength_4perc_1 = randomForest(x=genotypes_best4perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_4perc_1,file="rf_forklength_4perc_1.Rdata")
rf_forklength_4perc_2 = randomForest(x=genotypes_best4perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_4perc_2,file="rf_forklength_4perc_2.Rdata")
rf_forklength_4perc_3 = randomForest(x=genotypes_best4perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_4perc_3,file="rf_forklength_4perc_3.Rdata")

#load("rf_forklength_4perc_1.Rdata")
#load("rf_forklength_4perc_2.Rdata")
#load("rf_forklength_4perc_3.Rdata")

mean_rf_forklength_4perc_1_rsq <- mean(rf_forklength_4perc_1$rsq)
mean_rf_forklength_4perc_2_rsq <- mean(rf_forklength_4perc_2$rsq)
mean_rf_forklength_4perc_3_rsq <- mean(rf_forklength_4perc_3$rsq)

rm(rf_forklength_4perc_1,rf_forklength_4perc_2,rf_forklength_4perc_3)

##### Best 5% 
names_best_5perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.95))]
names_best_5perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.95))]
names_best_5perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.95))]
names_best_5perc_forklength_unique<-unique(c(names_best_5perc_forklength_1,names_best_5perc_forklength_2,names_best_5perc_forklength_3))
genotypes_best5perc_forklength<-All_data[,colnames(All_data) %in% names_best_5perc_forklength_unique]
rf_forklength_5perc_1 = randomForest(x=genotypes_best5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_5perc_1,file="rf_forklength_5perc_1.Rdata")
rf_forklength_5perc_2 = randomForest(x=genotypes_best5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_5perc_2,file="rf_forklength_5perc_2.Rdata")
rf_forklength_5perc_3 = randomForest(x=genotypes_best5perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_5perc_3,file="rf_forklength_5perc_3.Rdata")

#load("rf_forklength_5perc_1.Rdata")
#load("rf_forklength_5perc_2.Rdata")
#load("rf_forklength_5perc_3.Rdata")

mean_rf_forklength_5perc_1_rsq <- mean(rf_forklength_5perc_1$rsq)
mean_rf_forklength_5perc_2_rsq <- mean(rf_forklength_5perc_2$rsq)
mean_rf_forklength_5perc_3_rsq <- mean(rf_forklength_5perc_3$rsq)

rm(rf_forklength_5perc_1,rf_forklength_5perc_2,rf_forklength_5perc_3)

##### Best 10% 
names_best_10perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.90))]
names_best_10perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.90))]
names_best_10perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.90))]
names_best_10perc_forklength_unique<-unique(c(names_best_10perc_forklength_1,names_best_10perc_forklength_2,names_best_10perc_forklength_3))
genotypes_best10perc_forklength<-All_data[,colnames(All_data) %in% names_best_10perc_forklength_unique]
rf_forklength_10perc_1 = randomForest(x=genotypes_best10perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_10perc_1,file="rf_forklength_10perc_1.Rdata")
rf_forklength_10perc_2 = randomForest(x=genotypes_best10perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_10perc_2,file="rf_forklength_10perc_2.Rdata")
rf_forklength_10perc_3 = randomForest(x=genotypes_best10perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_10perc_3,file="rf_forklength_10perc_3.Rdata")

#load("rf_forklength_10perc_1.Rdata")
#load("rf_forklength_10perc_2.Rdata")
#load("rf_forklength_10perc_3.Rdata")

mean_rf_forklength_10perc_1_rsq <- mean(rf_forklength_10perc_1$rsq)
mean_rf_forklength_10perc_2_rsq <- mean(rf_forklength_10perc_2$rsq)
mean_rf_forklength_10perc_3_rsq <- mean(rf_forklength_10perc_3$rsq)

rm(rf_forklength_10perc_1,rf_forklength_10perc_2,rf_forklength_10perc_3)

##### Best 20% 
names_best_20perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.80))]
names_best_20perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.80))]
names_best_20perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.80))]
names_best_20perc_forklength_unique<-unique(c(names_best_20perc_forklength_1,names_best_20perc_forklength_2,names_best_20perc_forklength_3))
genotypes_best20perc_forklength<-All_data[,colnames(All_data) %in% names_best_20perc_forklength_unique]
rf_forklength_20perc_1 = randomForest(x=genotypes_best20perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_20perc_1,file="rf_forklength_20perc_1.Rdata")
rf_forklength_20perc_2 = randomForest(x=genotypes_best20perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_20perc_2,file="rf_forklength_20perc_2.Rdata")
rf_forklength_20perc_3 = randomForest(x=genotypes_best20perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_20perc_3,file="rf_forklength_20perc_3.Rdata")

#load("rf_forklength_20perc_1.Rdata")
#load("rf_forklength_20perc_2.Rdata")
#load("rf_forklength_20perc_3.Rdata")

mean_rf_forklength_20perc_1_rsq <- mean(rf_forklength_20perc_1$rsq)
mean_rf_forklength_20perc_2_rsq <- mean(rf_forklength_20perc_2$rsq)
mean_rf_forklength_20perc_3_rsq <- mean(rf_forklength_20perc_3$rsq)

rm(rf_forklength_20perc_1,rf_forklength_20perc_2,rf_forklength_20perc_3)

##### Best 30% 
names_best_30perc_forklength_1<-rownames(importance_rf_forklength_250000_1)[which(importance_rf_forklength_250000_1$importance > quantile(importance_rf_forklength_250000_1$importance, probs=0.70))]
names_best_30perc_forklength_2<-rownames(importance_rf_forklength_250000_2)[which(importance_rf_forklength_250000_2$importance > quantile(importance_rf_forklength_250000_2$importance, probs=0.70))]
names_best_30perc_forklength_3<-rownames(importance_rf_forklength_250000_3)[which(importance_rf_forklength_250000_3$importance > quantile(importance_rf_forklength_250000_3$importance, probs=0.70))]
names_best_30perc_forklength_unique<-unique(c(names_best_30perc_forklength_1,names_best_30perc_forklength_2,names_best_30perc_forklength_3))
genotypes_best30perc_forklength<-All_data[,colnames(All_data) %in% names_best_30perc_forklength_unique]
rf_forklength_30perc_1 = randomForest(x=genotypes_best30perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_30perc_1,file="rf_forklength_30perc_1.Rdata")
rf_forklength_30perc_2 = randomForest(x=genotypes_best30perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_30perc_2,file="rf_forklength_30perc_2.Rdata")
rf_forklength_30perc_3 = randomForest(x=genotypes_best30perc_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_forklength_30perc_3,file="rf_forklength_30perc_3.Rdata")

#load("rf_forklength_30perc_1.Rdata")
#load("rf_forklength_30perc_2.Rdata")
#load("rf_forklength_30perc_3.Rdata")

mean_rf_forklength_30perc_1_rsq <- mean(rf_forklength_30perc_1$rsq)
mean_rf_forklength_30perc_2_rsq <- mean(rf_forklength_30perc_2$rsq)
mean_rf_forklength_30perc_3_rsq <- mean(rf_forklength_30perc_3$rsq)

rm(rf_forklength_30perc_1,rf_forklength_30perc_2,rf_forklength_30perc_3)

## Variance explained-forklength
All_initial_Var_explained_forklength<-rbind(cbind(mean_rf_forklength_250000_1_rsq,mean_rf_forklength_250000_2_rsq,mean_rf_forklength_250000_3_rsq),
                                            cbind(mean_rf_forklength_0.5perc_1_rsq,mean_rf_forklength_0.5perc_2_rsq,mean_rf_forklength_0.5perc_3_rsq),
                                            cbind(mean_rf_forklength_1perc_1_rsq,mean_rf_forklength_1perc_2_rsq,mean_rf_forklength_1perc_3_rsq),
                                            cbind(mean_rf_forklength_1.5perc_1_rsq,mean_rf_forklength_1.5perc_2_rsq,mean_rf_forklength_1.5perc_3_rsq),
                                            cbind(mean_rf_forklength_2perc_1_rsq,mean_rf_forklength_2perc_2_rsq,mean_rf_forklength_2perc_3_rsq),
                                            cbind(mean_rf_forklength_3perc_1_rsq,mean_rf_forklength_3perc_2_rsq,mean_rf_forklength_3perc_3_rsq),
                                            cbind(mean_rf_forklength_4perc_1_rsq,mean_rf_forklength_4perc_2_rsq,mean_rf_forklength_4perc_3_rsq),
                                            cbind(mean_rf_forklength_5perc_1_rsq,mean_rf_forklength_5perc_2_rsq,mean_rf_forklength_5perc_3_rsq),
                                            cbind(mean_rf_forklength_10perc_1_rsq,mean_rf_forklength_10perc_2_rsq,mean_rf_forklength_10perc_3_rsq),
                                            cbind(mean_rf_forklength_20perc_1_rsq,mean_rf_forklength_20perc_2_rsq,mean_rf_forklength_20perc_3_rsq),
                                            cbind(mean_rf_forklength_30perc_1_rsq,mean_rf_forklength_30perc_2_rsq,mean_rf_forklength_30perc_3_rsq))

#forklength
All_initial_Var_explained_forklength<-data.frame(All_initial_Var_explained_forklength)
All_initial_Var_explained_forklength$Number_loci<-c(9108,length(names_best_0.5perc_forklength_unique),length(names_best_1perc_forklength_unique),length(names_best_1.5perc_forklength_unique),length(names_best_2perc_forklength_unique),length(names_best_3perc_forklength_unique),length(names_best_4perc_forklength_unique),length(names_best_5perc_forklength_unique),length(names_best_10perc_forklength_unique),length(names_best_20perc_forklength_unique),length(names_best_30perc_forklength_unique))
rownames(All_initial_Var_explained_forklength)<-c("All","Best0.5%","Best1%","Best1.5%","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_Var_explained_forklength$Average<-apply(All_initial_Var_explained_forklength[,1:3],1,mean)
write.csv(All_initial_Var_explained_forklength,file="All_initial_Var_explained_forklength_250ktrees.csv")
par(mar=c(5,6,3,3))
plot(All_initial_Var_explained_forklength$Number_loci,All_initial_Var_explained_forklength$Average,log="x",
     main="Variation Explained for Forklength-RF 250000 Trees", pch=19,xlab="Number of Loci", ylab="Proportion of Variation Explained",
     cex.lab=1.5,cex.axis=1.5)

#Based on this table and plot, the best 1% of loci have the largest r-squared value, which approximates to the % variance explained
# The best 1% data set includes 104 unique loci, so I'll run backward purging RF with the best 2% loci (n=215 unique)


#################### Backward purging approach - forklength
names_best_forklength_unique <- names_best_2perc_forklength_unique
length(names_best_forklength_unique)

names_best_forklength215 <- names_best_forklength_unique
genotypes_purging_forklength<-All_data[,colnames(All_data) %in% names_best_forklength215]

rf_forklength_purging_1 = randomForest(x=genotypes_purging_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_forklength_purging_1,file="rf_forklength_purging_1_50k.Rdata")
rf_forklength_purging_2 = randomForest(x=genotypes_purging_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_forklength_purging_2,file="rf_forklength_purging_2_50k.Rdata")

#Extract importance (mean decrease in accuracy) measure of each locus for forklength
importance_rf_forklength_purging_1<-data.frame(importance(rf_forklength_purging_1,type=1))
colnames(importance_rf_forklength_purging_1)<-c("Importance1")
importance_rf_forklength_purging_2<-data.frame(importance(rf_forklength_purging_2,type=1))
colnames(importance_rf_forklength_purging_2)<-c("Importance2")

importance_rf_purging <-cbind(rownames(importance_rf_forklength_purging_1),importance_rf_forklength_purging_1,importance_rf_forklength_purging_2)

#Plot importance values against each other to check if model converged
plot(importance_rf_purging$Importance1,importance_rf_purging$Importance2)
fit<-lm(Importance2~Importance1, data=importance_rf_purging)
summary(fit)
abline(fit)

#Run the third RF once you determine the appropriate number of trees to run per forest
rf_forklength_purging_3 = randomForest(x=genotypes_purging_forklength,y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_forklength_purging_3,file="rf_forklength_purging_3_50k.Rdata")

#load("rf_forklength_purging_1.Rdata")
#load("rf_forklength_purging_2.Rdata")
#load("rf_forklength_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[215]]<-names_best_forklength215
all_Var_explained_forklength_best<-data.frame(V1=1:215,V2=1:215,V3=1:215)
rownames(all_Var_explained_forklength_best)<-1:215
all_Var_explained_forklength_best[215,]<-c(mean(rf_forklength_purging_1$rsq),mean(rf_forklength_purging_2$rsq),mean(rf_forklength_purging_3$rsq))
all_mse_forklength<-data.frame(V1=1:215,V2=1:215,V3=1:215)
rownames(all_mse_forklength)<-1:215
all_mse_forklength[215,]<-c(mean(rf_forklength_purging_1$mse),mean(rf_forklength_purging_2$mse),mean(rf_forklength_purging_3$mse))

for (i in 1:214){
  print(i)
  imp_purging_1<-data.frame(importance(rf_forklength_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_forklength_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_forklength_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging_forklength)
  colnames(imp_purging_1)<-"X.IncMSE1"
  rownames(imp_purging_2)<-colnames(genotypes_purging_forklength)
  colnames(imp_purging_2)<-"X.IncMSE2"
  rownames(imp_purging_3)<-colnames(genotypes_purging_forklength)
  colnames(imp_purging_3)<-"X.IncMSE3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[215-i]]<-names_keep
  genotypes_purging_forklength<-All_data[,colnames(All_data) %in% names_keep]
  rf_forklength_purging_1 = randomForest(x = genotypes_purging_forklength, y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=50000)
  rf_forklength_purging_2 = randomForest(x = genotypes_purging_forklength, y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=50000)
  rf_forklength_purging_3 = randomForest(x = genotypes_purging_forklength, y = All_data$forklgth_Roza, importance=TRUE ,proximity=TRUE, ntree=50000)
  all_Var_explained_forklength_best[215-i,]<-c(mean(rf_forklength_purging_1$rsq),mean(rf_forklength_purging_2$rsq),mean(rf_forklength_purging_3$rsq))
  all_mse_forklength[215-i,]<-c(mean(rf_forklength_purging_1$mse),mean(rf_forklength_purging_2$mse),mean(rf_forklength_purging_3$mse))
  
  write.csv(all_Var_explained_forklength_best, file="All_var_explained_forklength.csv")
  lapply(names_all_iterations, write, "Backward_purging_names_all_iterations.txt", append=TRUE, ncolumns=220)
}

all_Var_explained_forklength_best$Average<-apply(all_Var_explained_forklength_best,1,mean)
write.csv(all_Var_explained_forklength_best, file="All_var_explained_forklength.csv")
lapply(names_all_iterations, write, "Backward_purging_names_all_iterations.txt", append=TRUE, ncolumns=220)

par(mar=c(5,5,3,3))

plot(all_Var_explained_forklength_best$Average[-c(1)],xlab="Number of Loci", ylab="Proportion Variation Explained",
     main="Backward Purging Forklength",cex.lab=1.5,cex.axis=1.5,pch=16)

which(all_Var_explained_forklength_best$Average==max(all_Var_explained_forklength_best$Average[-c(1)]))

write.csv(names_all_iterations[[44]],file="Purging_optimum_forklength_250000.csv")



##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################


# Charlie Waters
# University of Washington School of Aquatic and Fishery Sciences
# Random Forest analysis on CESRF Chinook for spawn time at CESRF

# Genotypes and phenotypes were corrected for all potential covariates

install.packages("randomForest")
library(randomForest)

# Read in data. For all loci, missing values were imputed prior to correction since RF can't take missing values. 
# The 1998 Founders were excluded since it provides the most balanced design between INT and SEG lines. 
# Individuals with missing phenotypic values or covariates needed to be excluded from the analysis, which reduced sample sizes. 

All_data <- read.csv("spawn_corrected_all_covariates.csv", header=TRUE)

# Initial RF with all 9108 loci for spawn timing
rf_spawn_350000_1 = randomForest(x = All_data[,3:9110], y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_350000_1,file="rf_spawn_350000_1.Rdata")
rf_spawn_350000_2 = randomForest(x = All_data[,3:9110], y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_350000_2,file="rf_spawn_350000_2.Rdata")

#Check correlation
importance_rf_spawn_350000_1<-data.frame(importance(rf_spawn_350000_1,type=1))
importance_rf_spawn_350000_2<-data.frame(importance(rf_spawn_350000_2,type=1))
imp <- cbind(importance_rf_spawn_350000_1,importance_rf_spawn_350000_2)
colnames(imp)<-c("Imp1","Imp2")
fit<-lm(Imp1~Imp2,data=imp)
summary(fit)

rf_spawn_350000_3 = randomForest(x = All_data[,3:9110], y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_350000_3,file="rf_spawn_350000_3.Rdata")

#load("rf_spawn_350000_1.Rdata")
#load("rf_spawn_350000_2.Rdata")
#load("rf_spawn_350000_3.Rdata")

mean_rf_spawn_350000_1_rsq <- mean(rf_spawn_350000_1$rsq)
mean_rf_spawn_350000_2_rsq <- mean(rf_spawn_350000_2$rsq)
mean_rf_spawn_350000_3_rsq <- mean(rf_spawn_350000_3$rsq)

#Plot of variable importance for spawn time
varImpPlot(rf_spawn_350000_1, main="All Loci spawn ntree=350000 1",n.var=50)
varImpPlot(rf_spawn_350000_2, main="All Loci spawn ntree=350000 2",n.var=50)
varImpPlot(rf_spawn_350000_3, main="All Loci spawn ntree=350000 3",n.var=50)

#Extract importance (mean decrease in accuracy) measure of each locus for spawn time
importance_rf_spawn_350000_1<-data.frame(importance(rf_spawn_350000_1,type=1))
colnames(importance_rf_spawn_350000_1)<-c("importance")
importance_rf_spawn_350000_2<-data.frame(importance(rf_spawn_350000_2,type=1))
colnames(importance_rf_spawn_350000_2)<-c("importance")
importance_rf_spawn_350000_3<-data.frame(importance(rf_spawn_350000_3,type=1))
colnames(importance_rf_spawn_350000_3)<-c("importance")

importance_rf_spawn_350000_all_loci <-cbind(rownames(importance_rf_spawn_350000_1),importance_rf_spawn_350000_1,importance_rf_spawn_350000_2, importance_rf_spawn_350000_3)
colnames(importance_rf_spawn_350000_all_loci)<-c("Variable","Importance1","Importance2", "Importance3")
write.csv(importance_rf_spawn_350000_all_loci,file="RF_spawn_importance_all_loci_350000.csv",row.names=FALSE)

#Plot importance values against each other to check if model converged
plot(importance_rf_spawn_350000_all_loci$Importance1,importance_rf_spawn_350000_all_loci$Importance2,xlim=c(-15,25),ylim=c(-15,25),main=("Imp Values spawn 250k, R2=0.96"))
fit<-lm(Importance2~Importance1, data=importance_rf_spawn_350000_all_loci)
summary(fit)
abline(fit)

################################ CESRF spawn

##### Best 0.5% 
names_best_0.5perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.995))]
names_best_0.5perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.995))]
names_best_0.5perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.995))]
names_best_0.5perc_spawn_unique<-unique(c(names_best_0.5perc_spawn_1,names_best_0.5perc_spawn_2,names_best_0.5perc_spawn_3))
genotypes_best0.5perc_spawn<-All_data[,colnames(All_data) %in% names_best_0.5perc_spawn_unique]
rf_spawn_0.5perc_1 = randomForest(x=genotypes_best0.5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_0.5perc_1,file="rf_spawn_0.5perc_1.Rdata")
rf_spawn_0.5perc_2 = randomForest(x=genotypes_best0.5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_0.5perc_2,file="rf_spawn_0.5perc_2.Rdata")
rf_spawn_0.5perc_3 = randomForest(x=genotypes_best0.5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_0.5perc_3,file="rf_spawn_0.5perc_3.Rdata")

#load("rf_spawn_0.5perc_1.Rdata")
#load("rf_spawn_0.5perc_2.Rdata")
#load("rf_spawn_0.5perc_3.Rdata")

mean_rf_spawn_0.5perc_1_rsq <- mean(rf_spawn_0.5perc_1$rsq)
mean_rf_spawn_0.5perc_2_rsq <- mean(rf_spawn_0.5perc_2$rsq)
mean_rf_spawn_0.5perc_3_rsq <- mean(rf_spawn_0.5perc_3$rsq)

rm(rf_spawn_0.5perc_1,rf_spawn_0.5perc_2,rf_spawn_0.5perc_3)

##### Best 1% 
names_best_1perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.99))]
names_best_1perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.99))]
names_best_1perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.99))]
names_best_1perc_spawn_unique<-unique(c(names_best_1perc_spawn_1,names_best_1perc_spawn_2,names_best_1perc_spawn_3))
genotypes_best1perc_spawn<-All_data[,colnames(All_data) %in% names_best_1perc_spawn_unique]
rf_spawn_1perc_1 = randomForest(x=genotypes_best1perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_1perc_1,file="rf_spawn_1perc_1.Rdata")
rf_spawn_1perc_2 = randomForest(x=genotypes_best1perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_1perc_2,file="rf_spawn_1perc_2.Rdata")
rf_spawn_1perc_3 = randomForest(x=genotypes_best1perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_1perc_3,file="rf_spawn_1perc_3.Rdata")

#load("rf_spawn_1perc_1.Rdata")
#load("rf_spawn_1perc_2.Rdata")
#load("rf_spawn_1perc_3.Rdata")

mean_rf_spawn_1perc_1_rsq <- mean(rf_spawn_1perc_1$rsq)
mean_rf_spawn_1perc_2_rsq <- mean(rf_spawn_1perc_2$rsq)
mean_rf_spawn_1perc_3_rsq <- mean(rf_spawn_1perc_3$rsq)

rm(rf_spawn_1perc_1,rf_spawn_1perc_2,rf_spawn_1perc_3)

##### Best 1.5%
names_best_1.5perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.985))]
names_best_1.5perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.985))]
names_best_1.5perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.985))]
names_best_1.5perc_spawn_unique<-unique(c(names_best_1.5perc_spawn_1,names_best_1.5perc_spawn_2,names_best_1.5perc_spawn_3))
genotypes_best1.5perc_spawn<-All_data[,colnames(All_data) %in% names_best_1.5perc_spawn_unique]
rf_spawn_1.5perc_1 = randomForest(x=genotypes_best1.5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_1.5perc_1,file="rf_spawn_1.5perc_1.Rdata")
rf_spawn_1.5perc_2 = randomForest(x=genotypes_best1.5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_1.5perc_2,file="rf_spawn_1.5perc_2.Rdata")
rf_spawn_1.5perc_3 = randomForest(x=genotypes_best1.5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_1.5perc_3,file="rf_spawn_1.5perc_3.Rdata")

#load("rf_spawn_1.5perc_1.Rdata")
#load("rf_spawn_1.5perc_2.Rdata")
#load("rf_spawn_1.5perc_3.Rdata")

mean_rf_spawn_1.5perc_1_rsq <- mean(rf_spawn_1.5perc_1$rsq)
mean_rf_spawn_1.5perc_2_rsq <- mean(rf_spawn_1.5perc_2$rsq)
mean_rf_spawn_1.5perc_3_rsq <- mean(rf_spawn_1.5perc_3$rsq)

rm(rf_spawn_1.5perc_1,rf_spawn_1.5perc_2,rf_spawn_1.5perc_3)

##### Best 2% 
names_best_2perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.98))]
names_best_2perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.98))]
names_best_2perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.98))]
names_best_2perc_spawn_unique<-unique(c(names_best_2perc_spawn_1,names_best_2perc_spawn_2,names_best_2perc_spawn_3))
genotypes_best2perc_spawn<-All_data[,colnames(All_data) %in% names_best_2perc_spawn_unique]
rf_spawn_2perc_1 = randomForest(x=genotypes_best2perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_2perc_1,file="rf_spawn_2perc_1.Rdata")
rf_spawn_2perc_2 = randomForest(x=genotypes_best2perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_2perc_2,file="rf_spawn_2perc_2.Rdata")
rf_spawn_2perc_3 = randomForest(x=genotypes_best2perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_2perc_3,file="rf_spawn_2perc_3.Rdata")

#load("rf_spawn_2perc_1.Rdata")
#load("rf_spawn_2perc_2.Rdata")
#load("rf_spawn_2perc_3.Rdata")

mean_rf_spawn_2perc_1_rsq <- mean(rf_spawn_2perc_1$rsq)
mean_rf_spawn_2perc_2_rsq <- mean(rf_spawn_2perc_2$rsq)
mean_rf_spawn_2perc_3_rsq <- mean(rf_spawn_2perc_3$rsq)

rm(rf_spawn_2perc_1,rf_spawn_2perc_2,rf_spawn_2perc_3)

##### Best 3% 
names_best_3perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.97))]
names_best_3perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.97))]
names_best_3perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.97))]
names_best_3perc_spawn_unique<-unique(c(names_best_3perc_spawn_1,names_best_3perc_spawn_2,names_best_3perc_spawn_3))
genotypes_best3perc_spawn<-All_data[,colnames(All_data) %in% names_best_3perc_spawn_unique]
rf_spawn_3perc_1 = randomForest(x=genotypes_best3perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_3perc_1,file="rf_spawn_3perc_1.Rdata")
rf_spawn_3perc_2 = randomForest(x=genotypes_best3perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_3perc_2,file="rf_spawn_3perc_2.Rdata")
rf_spawn_3perc_3 = randomForest(x=genotypes_best3perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_3perc_3,file="rf_spawn_3perc_3.Rdata")

#load("rf_spawn_3perc_1.Rdata")
#load("rf_spawn_3perc_2.Rdata")
#load("rf_spawn_3perc_3.Rdata")

mean_rf_spawn_3perc_1_rsq <- mean(rf_spawn_3perc_1$rsq)
mean_rf_spawn_3perc_2_rsq <- mean(rf_spawn_3perc_2$rsq)
mean_rf_spawn_3perc_3_rsq <- mean(rf_spawn_3perc_3$rsq)

rm(rf_spawn_3perc_1,rf_spawn_3perc_2,rf_spawn_3perc_3)

##### Best 4% 
names_best_4perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.96))]
names_best_4perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.96))]
names_best_4perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.96))]
names_best_4perc_spawn_unique<-unique(c(names_best_4perc_spawn_1,names_best_4perc_spawn_2,names_best_4perc_spawn_3))
genotypes_best4perc_spawn<-All_data[,colnames(All_data) %in% names_best_4perc_spawn_unique]
rf_spawn_4perc_1 = randomForest(x=genotypes_best4perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_4perc_1,file="rf_spawn_4perc_1.Rdata")
rf_spawn_4perc_2 = randomForest(x=genotypes_best4perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_4perc_2,file="rf_spawn_4perc_2.Rdata")
rf_spawn_4perc_3 = randomForest(x=genotypes_best4perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_4perc_3,file="rf_spawn_4perc_3.Rdata")

##load("rf_spawn_4perc_1.Rdata")
##load("rf_spawn_4perc_2.Rdata")
##load("rf_spawn_4perc_3.Rdata")

mean_rf_spawn_4perc_1_rsq <- mean(rf_spawn_4perc_1$rsq)
mean_rf_spawn_4perc_2_rsq <- mean(rf_spawn_4perc_2$rsq)
mean_rf_spawn_4perc_3_rsq <- mean(rf_spawn_4perc_3$rsq)

rm(rf_spawn_4perc_1,rf_spawn_4perc_2,rf_spawn_4perc_3)

##### Best 5% 
names_best_5perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.95))]
names_best_5perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.95))]
names_best_5perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.95))]
names_best_5perc_spawn_unique<-unique(c(names_best_5perc_spawn_1,names_best_5perc_spawn_2,names_best_5perc_spawn_3))
genotypes_best5perc_spawn<-All_data[,colnames(All_data) %in% names_best_5perc_spawn_unique]
rf_spawn_5perc_1 = randomForest(x=genotypes_best5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_5perc_1,file="rf_spawn_5perc_1.Rdata")
rf_spawn_5perc_2 = randomForest(x=genotypes_best5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_5perc_2,file="rf_spawn_5perc_2.Rdata")
rf_spawn_5perc_3 = randomForest(x=genotypes_best5perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_5perc_3,file="rf_spawn_5perc_3.Rdata")

##load("rf_spawn_5perc_1.Rdata")
##load("rf_spawn_5perc_2.Rdata")
##load("rf_spawn_5perc_3.Rdata")

mean_rf_spawn_5perc_1_rsq <- mean(rf_spawn_5perc_1$rsq)
mean_rf_spawn_5perc_2_rsq <- mean(rf_spawn_5perc_2$rsq)
mean_rf_spawn_5perc_3_rsq <- mean(rf_spawn_5perc_3$rsq)

rm(rf_spawn_5perc_1,rf_spawn_5perc_2,rf_spawn_5perc_3)

##### Best 10% 
names_best_10perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.90))]
names_best_10perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.90))]
names_best_10perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.90))]
names_best_10perc_spawn_unique<-unique(c(names_best_10perc_spawn_1,names_best_10perc_spawn_2,names_best_10perc_spawn_3))
genotypes_best10perc_spawn<-All_data[,colnames(All_data) %in% names_best_10perc_spawn_unique]
rf_spawn_10perc_1 = randomForest(x=genotypes_best10perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_10perc_1,file="rf_spawn_10perc_1.Rdata")
rf_spawn_10perc_2 = randomForest(x=genotypes_best10perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_10perc_2,file="rf_spawn_10perc_2.Rdata")
rf_spawn_10perc_3 = randomForest(x=genotypes_best10perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_10perc_3,file="rf_spawn_10perc_3.Rdata")

##load("rf_spawn_10perc_1.Rdata")
##load("rf_spawn_10perc_2.Rdata")
##load("rf_spawn_10perc_3.Rdata")

mean_rf_spawn_10perc_1_rsq <- mean(rf_spawn_10perc_1$rsq)
mean_rf_spawn_10perc_2_rsq <- mean(rf_spawn_10perc_2$rsq)
mean_rf_spawn_10perc_3_rsq <- mean(rf_spawn_10perc_3$rsq)

rm(rf_spawn_10perc_1,rf_spawn_10perc_2,rf_spawn_10perc_3)

##### Best 20% 
names_best_20perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.80))]
names_best_20perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.80))]
names_best_20perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.80))]
names_best_20perc_spawn_unique<-unique(c(names_best_20perc_spawn_1,names_best_20perc_spawn_2,names_best_20perc_spawn_3))
genotypes_best20perc_spawn<-All_data[,colnames(All_data) %in% names_best_20perc_spawn_unique]
rf_spawn_20perc_1 = randomForest(x=genotypes_best20perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_20perc_1,file="rf_spawn_20perc_1.Rdata")
rf_spawn_20perc_2 = randomForest(x=genotypes_best20perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_20perc_2,file="rf_spawn_20perc_2.Rdata")
rf_spawn_20perc_3 = randomForest(x=genotypes_best20perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_20perc_3,file="rf_spawn_20perc_3.Rdata")

##load("rf_spawn_20perc_1.Rdata")
##load("rf_spawn_20perc_2.Rdata")
##load("rf_spawn_20perc_3.Rdata")

mean_rf_spawn_20perc_1_rsq <- mean(rf_spawn_20perc_1$rsq)
mean_rf_spawn_20perc_2_rsq <- mean(rf_spawn_20perc_2$rsq)
mean_rf_spawn_20perc_3_rsq <- mean(rf_spawn_20perc_3$rsq)

rm(rf_spawn_20perc_1,rf_spawn_20perc_2,rf_spawn_20perc_3)

##### Best 30% 
names_best_30perc_spawn_1<-rownames(importance_rf_spawn_350000_1)[which(importance_rf_spawn_350000_1$importance > quantile(importance_rf_spawn_350000_1$importance, probs=0.70))]
names_best_30perc_spawn_2<-rownames(importance_rf_spawn_350000_2)[which(importance_rf_spawn_350000_2$importance > quantile(importance_rf_spawn_350000_2$importance, probs=0.70))]
names_best_30perc_spawn_3<-rownames(importance_rf_spawn_350000_3)[which(importance_rf_spawn_350000_3$importance > quantile(importance_rf_spawn_350000_3$importance, probs=0.70))]
names_best_30perc_spawn_unique<-unique(c(names_best_30perc_spawn_1,names_best_30perc_spawn_2,names_best_30perc_spawn_3))
genotypes_best30perc_spawn<-All_data[,colnames(All_data) %in% names_best_30perc_spawn_unique]
rf_spawn_30perc_1 = randomForest(x=genotypes_best30perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_30perc_1,file="rf_spawn_30perc_1.Rdata")
rf_spawn_30perc_2 = randomForest(x=genotypes_best30perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_30perc_2,file="rf_spawn_30perc_2.Rdata")
rf_spawn_30perc_3 = randomForest(x=genotypes_best30perc_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=350000)
save(rf_spawn_30perc_3,file="rf_spawn_30perc_3.Rdata")

##load("rf_spawn_30perc_1.Rdata")
##load("rf_spawn_30perc_2.Rdata")
##load("rf_spawn_30perc_3.Rdata")

mean_rf_spawn_30perc_1_rsq <- mean(rf_spawn_30perc_1$rsq)
mean_rf_spawn_30perc_2_rsq <- mean(rf_spawn_30perc_2$rsq)
mean_rf_spawn_30perc_3_rsq <- mean(rf_spawn_30perc_3$rsq)

rm(rf_spawn_30perc_1,rf_spawn_30perc_2,rf_spawn_30perc_3)

## Variance explained-spawn
All_initial_Var_explained_spawn<-rbind(cbind(mean_rf_spawn_350000_1_rsq,mean_rf_spawn_350000_2_rsq,mean_rf_spawn_350000_3_rsq),
                                       cbind(mean_rf_spawn_0.5perc_1_rsq,mean_rf_spawn_0.5perc_2_rsq,mean_rf_spawn_0.5perc_3_rsq),
                                       cbind(mean_rf_spawn_1perc_1_rsq,mean_rf_spawn_1perc_2_rsq,mean_rf_spawn_1perc_3_rsq),
                                       cbind(mean_rf_spawn_1.5perc_1_rsq,mean_rf_spawn_1.5perc_2_rsq,mean_rf_spawn_1.5perc_3_rsq),
                                       cbind(mean_rf_spawn_2perc_1_rsq,mean_rf_spawn_2perc_2_rsq,mean_rf_spawn_2perc_3_rsq),
                                       cbind(mean_rf_spawn_3perc_1_rsq,mean_rf_spawn_3perc_2_rsq,mean_rf_spawn_3perc_3_rsq),
                                       cbind(mean_rf_spawn_4perc_1_rsq,mean_rf_spawn_4perc_2_rsq,mean_rf_spawn_4perc_3_rsq),
                                       cbind(mean_rf_spawn_5perc_1_rsq,mean_rf_spawn_5perc_2_rsq,mean_rf_spawn_5perc_3_rsq),
                                       cbind(mean_rf_spawn_10perc_1_rsq,mean_rf_spawn_10perc_2_rsq,mean_rf_spawn_10perc_3_rsq),
                                       cbind(mean_rf_spawn_20perc_1_rsq,mean_rf_spawn_20perc_2_rsq,mean_rf_spawn_20perc_3_rsq),
                                       cbind(mean_rf_spawn_30perc_1_rsq,mean_rf_spawn_30perc_2_rsq,mean_rf_spawn_30perc_3_rsq))

#spawn
All_initial_Var_explained_spawn<-data.frame(All_initial_Var_explained_spawn)
All_initial_Var_explained_spawn$Number_loci<-c(9108,length(names_best_0.5perc_spawn_unique),length(names_best_1perc_spawn_unique),length(names_best_1.5perc_spawn_unique),length(names_best_2perc_spawn_unique),length(names_best_3perc_spawn_unique),length(names_best_4perc_spawn_unique),length(names_best_5perc_spawn_unique),length(names_best_10perc_spawn_unique),length(names_best_20perc_spawn_unique),length(names_best_30perc_spawn_unique))
rownames(All_initial_Var_explained_spawn)<-c("All","Best0.5%","Best1%","Best1.5%","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_Var_explained_spawn$Average<-apply(All_initial_Var_explained_spawn[,1:3],1,mean)
write.csv(All_initial_Var_explained_spawn,file="All_initial_Var_explained_spawn350k.csv")
par(mar=c(5,5,3,3))
plot(All_initial_Var_explained_spawn$Number_loci,All_initial_Var_explained_spawn$Average,log="x",
     xlab="Number of Loci",ylab="Proportion Variation Explained",
     main="Initial Variation Explained for Spawn Time",cex.axis=1.5,cex.lab=1.5,pch=19)


#Based on this table and plot, the best 0.5% of loci have the largest r-squared value, which approximates to the % variance explained
# The best 0.5% data set includes 53 unique loci, so I'll run backward purging RF with the best 1% loci (n=110 unique)


#################### Backward purging approach - spawn
names_best_spawn_unique <- names_best_1perc_spawn_unique
length(names_best_spawn_unique)

names_best_spawn110 <- names_best_spawn_unique
genotypes_purging_spawn<-All_data[,colnames(All_data) %in% names_best_spawn110]
rf_spawn_purging_1 = randomForest(x=genotypes_purging_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=75000)
save(rf_spawn_purging_1,file="rf_spawn_purging_1_75k.Rdata")
rf_spawn_purging_2 = randomForest(x=genotypes_purging_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=75000)
save(rf_spawn_purging_2,file="rf_spawn_purging_2_75k.Rdata")

#Extract importance (mean decrease in accuracy) measure of each locus for spawn timing
importance_rf_spawn_purging_1<-data.frame(importance(rf_spawn_purging_1,type=1))
colnames(importance_rf_spawn_purging_1)<-c("Importance1")
importance_rf_spawn_purging_2<-data.frame(importance(rf_spawn_purging_2,type=1))
colnames(importance_rf_spawn_purging_2)<-c("Importance2")

importance_rf_purging <-cbind(rownames(importance_rf_spawn_purging_1),importance_rf_spawn_purging_1,importance_rf_spawn_purging_2)

#Plot importance values against each other to check if model converged
plot(importance_rf_purging$Importance1,importance_rf_purging$Importance2)
fit<-lm(Importance2~Importance1, data=importance_rf_purging)
summary(fit)
abline(fit)


rf_spawn_purging_3 = randomForest(x=genotypes_purging_spawn,y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=75000)
save(rf_spawn_purging_3,file="rf_spawn_purging_3_75k.Rdata")

##load("rf_spawn_purging_1.Rdata")
##load("rf_spawn_purging_2.Rdata")
##load("rf_spawn_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[110]]<-names_best_spawn110
all_Var_explained_spawn_best<-data.frame(V1=1:110,V2=1:110,V3=1:110)
rownames(all_Var_explained_spawn_best)<-1:110
all_Var_explained_spawn_best[110,]<-c(mean(rf_spawn_purging_1$rsq),mean(rf_spawn_purging_2$rsq),mean(rf_spawn_purging_3$rsq))
all_mse_spawn<-data.frame(V1=1:110,V2=1:110,V3=1:110)
rownames(all_mse_spawn)<-1:110
all_mse_spawn[110,]<-c(mean(rf_spawn_purging_1$mse),mean(rf_spawn_purging_2$mse),mean(rf_spawn_purging_3$mse))

for (i in 1:109){
  print(i)
  imp_purging_1<-data.frame(importance(rf_spawn_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_spawn_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_spawn_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging_spawn)
  colnames(imp_purging_1)<-"X.IncMSE1"
  rownames(imp_purging_2)<-colnames(genotypes_purging_spawn)
  colnames(imp_purging_2)<-"X.IncMSE2"
  rownames(imp_purging_3)<-colnames(genotypes_purging_spawn)
  colnames(imp_purging_3)<-"X.IncMSE3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[110-i]]<-names_keep
  genotypes_purging_spawn<-All_data[,colnames(All_data) %in% names_keep]
  rf_spawn_purging_1 = randomForest(x = genotypes_purging_spawn, y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=75000)
  rf_spawn_purging_2 = randomForest(x = genotypes_purging_spawn, y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=75000)
  rf_spawn_purging_3 = randomForest(x = genotypes_purging_spawn, y = All_data$spawn_day_of_year, importance=TRUE ,proximity=TRUE, ntree=75000)
  all_Var_explained_spawn_best[110-i,]<-c(mean(rf_spawn_purging_1$rsq),mean(rf_spawn_purging_2$rsq),mean(rf_spawn_purging_3$rsq))
  all_mse_spawn[110-i,]<-c(mean(rf_spawn_purging_1$mse),mean(rf_spawn_purging_2$mse),mean(rf_spawn_purging_3$mse))
  
}


write.csv(all_Var_explained_spawn_best, file="Backward_purging_variance_explained_Spawn_75k.csv")

lapply(names_all_iterations, write, "Backward_purging_names_all_iterations.txt", append=TRUE, ncolumns=150)

all_Var_explained_spawn_best$Average<-apply(all_Var_explained_spawn_best,1,mean)

par(mar=c(5,5,3,3))
plot(all_Var_explained_spawn_best$Average[-c(1)],xlab="Number of Loci", ylab="Proportion Variation Explained",
     main="Backward Purging Spawn Timing 75k",cex.lab=1.5,cex.axis=1.5,pch=19)
abline(h=0)

which(all_Var_explained_spawn_best$Average==max(all_Var_explained_spawn_best$Average[-c(1)]))

write.csv(names_all_iterations[[68]],file="Purging_optimum_spawn_350000.csv")


##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################


# Charlie Waters
# University of Washington School of Aquatic and Fishery Sciences
# Random Forest analysis on CESRF Chinook for return time to Roza

# Genotypes and phenotypes were corrected for all potential covariates

install.packages("randomForest")
library(randomForest)

# Read in data. For all loci, missing values were imputed prior to correction since RF can't take missing values. 
# The 1998 Founders were excluded since it provides the most balanced design between INT and SEG lines. 
# Individuals with missing phenotypic values or covariates needed to be excluded from the analysis, which reduced sample sizes. 

All_data <- read.csv("return_corrected_all_covariates.csv", header=TRUE)

# Initial RF with all 9108 loci for Return time
rf_return_250000_1 = randomForest(x = All_data[,3:9110], y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_250000_1,file="rf_return_250000_1.Rdata")
rf_return_250000_2 = randomForest(x = All_data[,3:9110], y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_250000_2,file="rf_return_250000_2.Rdata")

#Check correlation early
importance_rf_return_250000_1<-data.frame(importance(rf_return_250000_1,type=1))
importance_rf_return_250000_2<-data.frame(importance(rf_return_250000_2,type=1))
imp<-cbind(importance_rf_return_250000_1,importance_rf_return_250000_2)
colnames(imp)<-c("Imp1","Imp2")
fit <- lm(Imp1~Imp2,data=imp)
summary(fit)

rf_return_250000_3 = randomForest(x = All_data[,3:9110], y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_250000_3,file="rf_return_250000_3.Rdata")

#load("rf_return_250000_1.Rdata")
#load("rf_return_250000_2.Rdata")
#load("rf_return_250000_3.Rdata")

mean_rf_return_250000_1_rsq <- mean(rf_return_250000_1$rsq)
mean_rf_return_250000_2_rsq <- mean(rf_return_250000_2$rsq)
mean_rf_return_250000_3_rsq <- mean(rf_return_250000_3$rsq)

#Plot of variable importance for return time
varImpPlot(rf_return_250000_1, main="All Loci return ntree=250000 1",n.var=50)
varImpPlot(rf_return_250000_2, main="All Loci return ntree=250000 2",n.var=50)
varImpPlot(rf_return_250000_3, main="All Loci return ntree=250000 3",n.var=50)

#Extract importance (mean decrease in accuracy) measure of each locus for return time
importance_rf_return_250000_1<-data.frame(importance(rf_return_250000_1,type=1))
colnames(importance_rf_return_250000_1)<-c("importance")
importance_rf_return_250000_2<-data.frame(importance(rf_return_250000_2,type=1))
colnames(importance_rf_return_250000_2)<-c("importance")
importance_rf_return_250000_3<-data.frame(importance(rf_return_250000_3,type=1))
colnames(importance_rf_return_250000_3)<-c("importance")

importance_rf_return_250000_all_loci <-cbind(rownames(importance_rf_return_250000_1),importance_rf_return_250000_1,importance_rf_return_250000_2, importance_rf_return_250000_3)
colnames(importance_rf_return_250000_all_loci)<-c("Variable","Importance1","Importance2", "Importance3")
write.csv(importance_rf_return_250000_all_loci,file="RF_return_importance_all_loci_250000.csv",row.names=FALSE)

#Plot importance values against each other to check if model converged
plot(importance_rf_return_250000_all_loci$Importance1,importance_rf_return_250000_all_loci$Importance2,xlim=c(-15,25),ylim=c(-15,25),main=("Imp Values return 150k, R2=0.96"))
fit<-lm(Importance2~Importance1, data=importance_rf_return_250000_all_loci)
summary(fit)
abline(fit)

################################ CESRF return

##### Best 0.5% 
names_best_0.5perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.995))]
names_best_0.5perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.995))]
names_best_0.5perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.995))]
names_best_0.5perc_return_unique<-unique(c(names_best_0.5perc_return_1,names_best_0.5perc_return_2,names_best_0.5perc_return_3))
genotypes_best0.5perc_return<-All_data[,colnames(All_data) %in% names_best_0.5perc_return_unique]
rf_return_0.5perc_1 = randomForest(x=genotypes_best0.5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_0.5perc_1,file="rf_return_0.5perc_1.Rdata")
rf_return_0.5perc_2 = randomForest(x=genotypes_best0.5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_0.5perc_2,file="rf_return_0.5perc_2.Rdata")
rf_return_0.5perc_3 = randomForest(x=genotypes_best0.5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_0.5perc_3,file="rf_return_0.5perc_3.Rdata")

#load("rf_return_0.5perc_1.Rdata")
#load("rf_return_0.5perc_2.Rdata")
#load("rf_return_0.5perc_3.Rdata")

mean_rf_return_0.5perc_1_rsq <- mean(rf_return_0.5perc_1$rsq)
mean_rf_return_0.5perc_2_rsq <- mean(rf_return_0.5perc_2$rsq)
mean_rf_return_0.5perc_3_rsq <- mean(rf_return_0.5perc_3$rsq)

rm(rf_return_0.5perc_1,rf_return_0.5perc_2,rf_return_0.5perc_3)

##### Best 1% 
names_best_1perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.99))]
names_best_1perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.99))]
names_best_1perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.99))]
names_best_1perc_return_unique<-unique(c(names_best_1perc_return_1,names_best_1perc_return_2,names_best_1perc_return_3))
genotypes_best1perc_return<-All_data[,colnames(All_data) %in% names_best_1perc_return_unique]
rf_return_1perc_1 = randomForest(x=genotypes_best1perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_1perc_1,file="rf_return_1perc_1.Rdata")
rf_return_1perc_2 = randomForest(x=genotypes_best1perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_1perc_2,file="rf_return_1perc_2.Rdata")
rf_return_1perc_3 = randomForest(x=genotypes_best1perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_1perc_3,file="rf_return_1perc_3.Rdata")

#load("rf_return_1perc_1.Rdata")
#load("rf_return_1perc_2.Rdata")
#load("rf_return_1perc_3.Rdata")

mean_rf_return_1perc_1_rsq <- mean(rf_return_1perc_1$rsq)
mean_rf_return_1perc_2_rsq <- mean(rf_return_1perc_2$rsq)
mean_rf_return_1perc_3_rsq <- mean(rf_return_1perc_3$rsq)

rm(rf_return_1perc_1,rf_return_1perc_2,rf_return_1perc_3)

##### Best 1.5%
names_best_1.5perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.985))]
names_best_1.5perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.985))]
names_best_1.5perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.985))]
names_best_1.5perc_return_unique<-unique(c(names_best_1.5perc_return_1,names_best_1.5perc_return_2,names_best_1.5perc_return_3))
genotypes_best1.5perc_return<-All_data[,colnames(All_data) %in% names_best_1.5perc_return_unique]
rf_return_1.5perc_1 = randomForest(x=genotypes_best1.5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_1.5perc_1,file="rf_return_1.5perc_1.Rdata")
rf_return_1.5perc_2 = randomForest(x=genotypes_best1.5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_1.5perc_2,file="rf_return_1.5perc_2.Rdata")
rf_return_1.5perc_3 = randomForest(x=genotypes_best1.5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_1.5perc_3,file="rf_return_1.5perc_3.Rdata")

#load("rf_return_1.5perc_1.Rdata")
#load("rf_return_1.5perc_2.Rdata")
#load("rf_return_1.5perc_3.Rdata")

mean_rf_return_1.5perc_1_rsq <- mean(rf_return_1.5perc_1$rsq)
mean_rf_return_1.5perc_2_rsq <- mean(rf_return_1.5perc_2$rsq)
mean_rf_return_1.5perc_3_rsq <- mean(rf_return_1.5perc_3$rsq)

rm(rf_return_1.5perc_1,rf_return_1.5perc_2,rf_return_1.5perc_3)

##### Best 2% 
names_best_2perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.98))]
names_best_2perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.98))]
names_best_2perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.98))]
names_best_2perc_return_unique<-unique(c(names_best_2perc_return_1,names_best_2perc_return_2,names_best_2perc_return_3))
genotypes_best2perc_return<-All_data[,colnames(All_data) %in% names_best_2perc_return_unique]
rf_return_2perc_1 = randomForest(x=genotypes_best2perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_2perc_1,file="rf_return_2perc_1.Rdata")
rf_return_2perc_2 = randomForest(x=genotypes_best2perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_2perc_2,file="rf_return_2perc_2.Rdata")
rf_return_2perc_3 = randomForest(x=genotypes_best2perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_2perc_3,file="rf_return_2perc_3.Rdata")

#load("rf_return_2perc_1.Rdata")
#load("rf_return_2perc_2.Rdata")
#load("rf_return_2perc_3.Rdata")

mean_rf_return_2perc_1_rsq <- mean(rf_return_2perc_1$rsq)
mean_rf_return_2perc_2_rsq <- mean(rf_return_2perc_2$rsq)
mean_rf_return_2perc_3_rsq <- mean(rf_return_2perc_3$rsq)

rm(rf_return_2perc_1,rf_return_2perc_2,rf_return_2perc_3)

##### Best 3% 
names_best_3perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.97))]
names_best_3perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.97))]
names_best_3perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.97))]
names_best_3perc_return_unique<-unique(c(names_best_3perc_return_1,names_best_3perc_return_2,names_best_3perc_return_3))
genotypes_best3perc_return<-All_data[,colnames(All_data) %in% names_best_3perc_return_unique]
rf_return_3perc_1 = randomForest(x=genotypes_best3perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_3perc_1,file="rf_return_3perc_1.Rdata")
rf_return_3perc_2 = randomForest(x=genotypes_best3perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_3perc_2,file="rf_return_3perc_2.Rdata")
rf_return_3perc_3 = randomForest(x=genotypes_best3perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_3perc_3,file="rf_return_3perc_3.Rdata")

#load("rf_return_3perc_1.Rdata")
#load("rf_return_3perc_2.Rdata")
#load("rf_return_3perc_3.Rdata")

mean_rf_return_3perc_1_rsq <- mean(rf_return_3perc_1$rsq)
mean_rf_return_3perc_2_rsq <- mean(rf_return_3perc_2$rsq)
mean_rf_return_3perc_3_rsq <- mean(rf_return_3perc_3$rsq)

rm(rf_return_3perc_1,rf_return_3perc_2,rf_return_3perc_3)

##### Best 4% 
names_best_4perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.96))]
names_best_4perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.96))]
names_best_4perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.96))]
names_best_4perc_return_unique<-unique(c(names_best_4perc_return_1,names_best_4perc_return_2,names_best_4perc_return_3))
genotypes_best4perc_return<-All_data[,colnames(All_data) %in% names_best_4perc_return_unique]
rf_return_4perc_1 = randomForest(x=genotypes_best4perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_4perc_1,file="rf_return_4perc_1.Rdata")
rf_return_4perc_2 = randomForest(x=genotypes_best4perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_4perc_2,file="rf_return_4perc_2.Rdata")
rf_return_4perc_3 = randomForest(x=genotypes_best4perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_4perc_3,file="rf_return_4perc_3.Rdata")

#load("rf_return_4perc_1.Rdata")
#load("rf_return_4perc_2.Rdata")
#load("rf_return_4perc_3.Rdata")

mean_rf_return_4perc_1_rsq <- mean(rf_return_4perc_1$rsq)
mean_rf_return_4perc_2_rsq <- mean(rf_return_4perc_2$rsq)
mean_rf_return_4perc_3_rsq <- mean(rf_return_4perc_3$rsq)

rm(rf_return_4perc_1,rf_return_4perc_2,rf_return_4perc_3)

##### Best 5% 
names_best_5perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.95))]
names_best_5perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.95))]
names_best_5perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.95))]
names_best_5perc_return_unique<-unique(c(names_best_5perc_return_1,names_best_5perc_return_2,names_best_5perc_return_3))
genotypes_best5perc_return<-All_data[,colnames(All_data) %in% names_best_5perc_return_unique]
rf_return_5perc_1 = randomForest(x=genotypes_best5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_5perc_1,file="rf_return_5perc_1.Rdata")
rf_return_5perc_2 = randomForest(x=genotypes_best5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_5perc_2,file="rf_return_5perc_2.Rdata")
rf_return_5perc_3 = randomForest(x=genotypes_best5perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_5perc_3,file="rf_return_5perc_3.Rdata")

#load("rf_return_5perc_1.Rdata")
#load("rf_return_5perc_2.Rdata")
#load("rf_return_5perc_3.Rdata")

mean_rf_return_5perc_1_rsq <- mean(rf_return_5perc_1$rsq)
mean_rf_return_5perc_2_rsq <- mean(rf_return_5perc_2$rsq)
mean_rf_return_5perc_3_rsq <- mean(rf_return_5perc_3$rsq)

rm(rf_return_5perc_1,rf_return_5perc_2,rf_return_5perc_3)

##### Best 10% 
names_best_10perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.90))]
names_best_10perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.90))]
names_best_10perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.90))]
names_best_10perc_return_unique<-unique(c(names_best_10perc_return_1,names_best_10perc_return_2,names_best_10perc_return_3))
genotypes_best10perc_return<-All_data[,colnames(All_data) %in% names_best_10perc_return_unique]
rf_return_10perc_1 = randomForest(x=genotypes_best10perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_10perc_1,file="rf_return_10perc_1.Rdata")
rf_return_10perc_2 = randomForest(x=genotypes_best10perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_10perc_2,file="rf_return_10perc_2.Rdata")
rf_return_10perc_3 = randomForest(x=genotypes_best10perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_10perc_3,file="rf_return_10perc_3.Rdata")

#load("rf_return_10perc_1.Rdata")
#load("rf_return_10perc_2.Rdata")
#load("rf_return_10perc_3.Rdata")

mean_rf_return_10perc_1_rsq <- mean(rf_return_10perc_1$rsq)
mean_rf_return_10perc_2_rsq <- mean(rf_return_10perc_2$rsq)
mean_rf_return_10perc_3_rsq <- mean(rf_return_10perc_3$rsq)

rm(rf_return_10perc_1,rf_return_10perc_2,rf_return_10perc_3)

##### Best 20% 
names_best_20perc_return_1<-rownames(importance_rf_return_250000_1)[which(importance_rf_return_250000_1$importance > quantile(importance_rf_return_250000_1$importance, probs=0.80))]
names_best_20perc_return_2<-rownames(importance_rf_return_250000_2)[which(importance_rf_return_250000_2$importance > quantile(importance_rf_return_250000_2$importance, probs=0.80))]
names_best_20perc_return_3<-rownames(importance_rf_return_250000_3)[which(importance_rf_return_250000_3$importance > quantile(importance_rf_return_250000_3$importance, probs=0.80))]
names_best_20perc_return_unique<-unique(c(names_best_20perc_return_1,names_best_20perc_return_2,names_best_20perc_return_3))
genotypes_best20perc_return<-All_data[,colnames(All_data) %in% names_best_20perc_return_unique]
rf_return_20perc_1 = randomForest(x=genotypes_best20perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_20perc_1,file="rf_return_20perc_1.Rdata")
rf_return_20perc_2 = randomForest(x=genotypes_best20perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_20perc_2,file="rf_return_20perc_2.Rdata")
rf_return_20perc_3 = randomForest(x=genotypes_best20perc_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=250000)
save(rf_return_20perc_3,file="rf_return_20perc_3.Rdata")

#load("rf_return_20perc_1.Rdata")
#load("rf_return_20perc_2.Rdata")
#load("rf_return_20perc_3.Rdata")

mean_rf_return_20perc_1_rsq <- mean(rf_return_20perc_1$rsq)
mean_rf_return_20perc_2_rsq <- mean(rf_return_20perc_2$rsq)
mean_rf_return_20perc_3_rsq <- mean(rf_return_20perc_3$rsq)

rm(rf_return_20perc_1,rf_return_20perc_2,rf_return_20perc_3)

## Variance explained-return
All_initial_Var_explained_return<-rbind(cbind(mean_rf_return_250000_1_rsq,mean_rf_return_250000_2_rsq,mean_rf_return_250000_3_rsq),
                                        cbind(mean_rf_return_0.5perc_1_rsq,mean_rf_return_0.5perc_2_rsq,mean_rf_return_0.5perc_3_rsq),
                                        cbind(mean_rf_return_1perc_1_rsq,mean_rf_return_1perc_2_rsq,mean_rf_return_1perc_3_rsq),
                                        cbind(mean_rf_return_1.5perc_1_rsq,mean_rf_return_1.5perc_2_rsq,mean_rf_return_1.5perc_3_rsq),
                                        cbind(mean_rf_return_2perc_1_rsq,mean_rf_return_2perc_2_rsq,mean_rf_return_2perc_3_rsq),
                                        cbind(mean_rf_return_3perc_1_rsq,mean_rf_return_3perc_2_rsq,mean_rf_return_3perc_3_rsq),
                                        cbind(mean_rf_return_4perc_1_rsq,mean_rf_return_4perc_2_rsq,mean_rf_return_4perc_3_rsq),
                                        cbind(mean_rf_return_5perc_1_rsq,mean_rf_return_5perc_2_rsq,mean_rf_return_5perc_3_rsq),
                                        cbind(mean_rf_return_10perc_1_rsq,mean_rf_return_10perc_2_rsq,mean_rf_return_10perc_3_rsq),
                                        cbind(mean_rf_return_20perc_1_rsq,mean_rf_return_20perc_2_rsq,mean_rf_return_20perc_3_rsq))


#return
All_initial_Var_explained_return<-data.frame(All_initial_Var_explained_return)
All_initial_Var_explained_return$Number_loci<-c(9108,length(names_best_0.5perc_return_unique),length(names_best_1perc_return_unique),length(names_best_1.5perc_return_unique),length(names_best_2perc_return_unique),length(names_best_3perc_return_unique),length(names_best_4perc_return_unique),length(names_best_5perc_return_unique),length(names_best_10perc_return_unique),length(names_best_20perc_return_unique))
rownames(All_initial_Var_explained_return)<-c("All","Best0.5%","Best1%","Best1.5%","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%")
All_initial_Var_explained_return$Average<-apply(All_initial_Var_explained_return[,1:3],1,mean)
write.csv(All_initial_Var_explained_return, file="Intial Variance Explained Return Time.csv")
par(mar=c(5,5,3,3))
plot(All_initial_Var_explained_return$Number_loci,All_initial_Var_explained_return$Average,log="x",
     cex.axis=1.5,cex.lab=1.5,xlab="Number of Loci", ylab="Proportion of Variation Explained",
     main="Variation in Return Time Explained by Loci",pch=19)


#Based on this table and plot, the best 0.5% of loci have the largest r-squared value, which approximates to the % variance explained
# The best 0.5% data set includes 50 unique loci, so I'll run backward purging RF with the best 1% loci (n=101 unique)


#################### Backward purging approach - return time
names_best_return_unique <- names_best_1perc_return_unique
length(names_best_return_unique)

names_best_return101 <- names_best_return_unique
genotypes_purging_return<-All_data[,colnames(All_data) %in% names_best_return101]
rf_return_purging_1 = randomForest(x=genotypes_purging_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_return_purging_1,file="rf_return_purging_1_50k.Rdata")
rf_return_purging_2 = randomForest(x=genotypes_purging_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_return_purging_2,file="rf_return_purging_2_50k.Rdata")

#check correlation between importance values
importance_rf_return_purging_1<-data.frame(importance(rf_return_purging_1,type=1))
importance_rf_return_purging_2<-data.frame(importance(rf_return_purging_2,type=1))
imp<-cbind(importance_rf_return_purging_1,importance_rf_return_purging_2)
colnames(imp)<-c("Imp1","Imp2")
fit <- lm(Imp1~Imp2,data=imp)
summary(fit)

rf_return_purging_3 = randomForest(x=genotypes_purging_return,y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=50000)
save(rf_return_purging_3,file="rf_return_purging_3_50k.Rdata")

#load("rf_return_purging_1.Rdata")
#load("rf_return_purging_2.Rdata")
#load("rf_return_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[101]]<-names_best_return101
all_Var_explained_return_best<-data.frame(V1=1:101,V2=1:101,V3=1:101)
rownames(all_Var_explained_return_best)<-1:101
all_Var_explained_return_best[101,]<-c(mean(rf_return_purging_1$rsq),mean(rf_return_purging_2$rsq),mean(rf_return_purging_3$rsq))
all_mse_return<-data.frame(V1=1:101,V2=1:101,V3=1:101)
rownames(all_mse_return)<-1:101
all_mse_return[101,]<-c(mean(rf_return_purging_1$mse),mean(rf_return_purging_2$mse),mean(rf_return_purging_3$mse))

for (i in 1:100){
  print(i)
  imp_purging_1<-data.frame(importance(rf_return_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_return_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_return_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging_return)
  colnames(imp_purging_1)<-"X.IncMSE1"
  rownames(imp_purging_2)<-colnames(genotypes_purging_return)
  colnames(imp_purging_2)<-"X.IncMSE2"
  rownames(imp_purging_3)<-colnames(genotypes_purging_return)
  colnames(imp_purging_3)<-"X.IncMSE3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[101-i]]<-names_keep
  genotypes_purging_return<-All_data[,colnames(All_data) %in% names_keep]
  rf_return_purging_1 = randomForest(x = genotypes_purging_return, y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=50000)
  rf_return_purging_2 = randomForest(x = genotypes_purging_return, y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=50000)
  rf_return_purging_3 = randomForest(x = genotypes_purging_return, y = All_data$roza_day_of_year, importance=TRUE ,proximity=TRUE, ntree=50000)
  all_Var_explained_return_best[101-i,]<-c(mean(rf_return_purging_1$rsq),mean(rf_return_purging_2$rsq),mean(rf_return_purging_3$rsq))
  all_mse_return[101-i,]<-c(mean(rf_return_purging_1$mse),mean(rf_return_purging_2$mse),mean(rf_return_purging_3$mse))
  
  write.csv(all_Var_explained_return_best, file="All_var_explained_return.csv")
  lapply(names_all_iterations, write, "Backward_purging_return_names_all_iterations_temp.txt", append=TRUE, ncolumns=105)
  
}

all_Var_explained_return_best$Average<-apply(all_Var_explained_return_best,1,mean)

write.csv(all_Var_explained_return_best, file="Backward_purging_variance_explained_return.csv")

lapply(names_all_iterations, write, "Backward_purging_return_names_all_iterations.txt", append=TRUE, ncolumns=105)

par(mar=c(5,5,3,3))
plot(all_Var_explained_return_best$Average[-c(1)],xlab="Number of Loci", ylab="Proportion Variation Explained",
     main="Backward Purging Return Time",cex.lab=1.5,cex.axis=1.5,pch=19)

which(all_Var_explained_return_best$Average==max(all_Var_explained_return_best$Average[-c(1)]))

write.csv(names_all_iterations[[26]],file="Purging_optimum_return_250000.csv") #Insert index of max variation


##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################

# Charlie Waters
# University of Washington School of Aquatic and Fishery Sciences
# Random Forest analysis on CESRF Chinook for age at maturity

# Genotypes were corrected for all potential covariates; phenotypes were not corrected

install.packages("randomForest")
library(randomForest)

# Read in data. For all loci, missing values were imputed prior to correction since RF can't take missing values. 
# The 1998 Founders were excluded since it provides the most balanced design between INT and SEG lines. 
# Individuals with missing phenotypic values or covariates needed to be excluded from the analysis, which reduced sample sizes. 

All_data <- read.csv("Age_genotypes_only_corrected_all_covariates.csv", header=TRUE)

# Need to know how many age 3, 4, and 5 individuals you have, which will inform the random forest 'strata' option

hist(All_data$age)
length(which(All_data$age==3)) #23 fish age 3
length(which(All_data$age==4)) #344 fish age 4
length(which(All_data$age==5)) #16 fish age 5

# Typically, 2/3 of the data is randomly chosen for the training data set, which is then used by RF to grow the trees and forest. 
# To balance the representation of the different age classes, we take 2/3 of the lowest sample size and use it with the strata option 
# Age 5 fish have the lowest sample size, so 16*(2/3), or about 10 fish, will be used in the sampling.
sample_size <- c(10,10,10)

# Initial RF with all 9108 loci for Age at Maturation
rf_age_750000_1 = randomForest(x = All_data[,3:9110], y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_750000_1,file="rf_age_750000_1.Rdata")
rf_age_750000_2 = randomForest(x = All_data[,3:9110], y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_750000_2,file="rf_age_750000_2.Rdata")

#Check correlation early
importance_rf_age_750000_1<-data.frame(importance(rf_age_750000_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy
importance_rf_age_750000_2<-data.frame(importance(rf_age_750000_2,type=1))
imp<-cbind(importance_rf_age_750000_1,importance_rf_age_750000_2)
colnames(imp)<-c("Imp1","Imp2")
fit <- lm(Imp1~Imp2,data=imp)
summary(fit)

rf_age_750000_3 = randomForest(x = All_data[,3:9110], y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_750000_3,file="rf_age_750000_3.Rdata")

#load("rf_age_750000_1.Rdata")
#load("rf_age_750000_2.Rdata")
#load("rf_age_750000_3.Rdata")

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest. 
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_age_750000_1_err_rate_converge <- rf_age_750000_1$err.rate[750000]
rf_age_750000_2_err_rate_converge <- rf_age_750000_2$err.rate[750000]
rf_age_750000_3_err_rate_converge <- rf_age_750000_3$err.rate[750000]

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

#Plot of variable importance for CESRF age
varImpPlot(rf_age_750000_1, main="All Loci Age ntree=750000 1",n.var=50)
varImpPlot(rf_age_750000_2, main="All Loci Age ntree=750000 2",n.var=50)
varImpPlot(rf_age_750000_3, main="All Loci Age ntree=750000 3",n.var=50)

#Extract importance (mean decrease in accuracy) measure of each locus
importance_rf_age_750000_1<-data.frame(importance(rf_age_750000_1,type=1))
colnames(importance_rf_age_750000_1)<-c("importance")
importance_rf_age_750000_2<-data.frame(importance(rf_age_750000_2,type=1))
colnames(importance_rf_age_750000_2)<-c("importance")
importance_rf_age_750000_3<-data.frame(importance(rf_age_750000_3,type=1))
colnames(importance_rf_age_750000_3)<-c("importance")

importance_rf_age_750000_all_loci <-cbind(rownames(importance_rf_age_750000_1),importance_rf_age_750000_1,importance_rf_age_750000_2, importance_rf_age_750000_3)
colnames(importance_rf_age_750000_all_loci)<-c("Variable","Importance1","Importance2", "Importance3")
write.csv(importance_rf_age_750000_all_loci,file="rf_age_importance_all_loci_750000.csv",row.names=FALSE)

#Plot importance values against each other to check if model converged
plot(importance_rf_age_750000_all_loci$Importance1,importance_rf_age_750000_all_loci$Importance2,xlim=c(-15,25),ylim=c(-15,25),main=("Imp Values Age 750k, R2=0.87"))
fit<-lm(Importance2~Importance1, data=importance_rf_age_750000_all_loci)
summary(fit)
abline(fit)


##### Best 0.5% 
names_best_0.5perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.995))]
names_best_0.5perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.995))]
names_best_0.5perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.995))]
names_best_0.5perc_age_unique<-unique(c(names_best_0.5perc_age_1,names_best_0.5perc_age_2,names_best_0.5perc_age_3))
genotypes_best0.5perc_age<-All_data[,colnames(All_data) %in% names_best_0.5perc_age_unique]
rf_age_0.5perc_1 = randomForest(x=genotypes_best0.5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_0.5perc_1,file="rf_age_0.5perc_1.Rdata")
rf_age_0.5perc_2 = randomForest(x=genotypes_best0.5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_0.5perc_2,file="rf_age_0.5perc_2.Rdata")
rf_age_0.5perc_3 = randomForest(x=genotypes_best0.5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_0.5perc_3,file="rf_age_0.5perc_3.Rdata")

#load("rf_age_0.5perc_1.Rdata")
#load("rf_age_0.5perc_2.Rdata")
#load("rf_age_0.5perc_3.Rdata")

rf_age_0.5perc_1_err_rate_converge <- rf_age_0.5perc_1$err.rate[750000]
rf_age_0.5perc_2_err_rate_converge <- rf_age_0.5perc_2$err.rate[750000]
rf_age_0.5perc_3_err_rate_converge <- rf_age_0.5perc_3$err.rate[750000]

rm(rf_age_0.5perc_1,rf_age_0.5perc_2,rf_age_0.5perc_3)

##### Best 1% 
names_best_1perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.99))]
names_best_1perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.99))]
names_best_1perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.99))]
names_best_1perc_age_unique<-unique(c(names_best_1perc_age_1,names_best_1perc_age_2,names_best_1perc_age_3))
genotypes_best1perc_age<-All_data[,colnames(All_data) %in% names_best_1perc_age_unique]
rf_age_1perc_1 = randomForest(x=genotypes_best1perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_1perc_1,file="rf_age_1perc_1.Rdata")
rf_age_1perc_2 = randomForest(x=genotypes_best1perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_1perc_2,file="rf_age_1perc_2.Rdata")
rf_age_1perc_3 = randomForest(x=genotypes_best1perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_1perc_3,file="rf_age_1perc_3.Rdata")

#load("rf_age_1perc_1.Rdata")
#load("rf_age_1perc_2.Rdata")
#load("rf_age_1perc_3.Rdata")

rf_age_1perc_1_err_rate_converge <- rf_age_1perc_1$err.rate[750000]
rf_age_1perc_2_err_rate_converge <- rf_age_1perc_2$err.rate[750000]
rf_age_1perc_3_err_rate_converge <- rf_age_1perc_3$err.rate[750000]

rm(rf_age_1perc_1,rf_age_1perc_2,rf_age_1perc_3)

##### Best 1.5%
names_best_1.5perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.985))]
names_best_1.5perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.985))]
names_best_1.5perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.985))]
names_best_1.5perc_age_unique<-unique(c(names_best_1.5perc_age_1,names_best_1.5perc_age_2,names_best_1.5perc_age_3))
genotypes_best1.5perc_age<-All_data[,colnames(All_data) %in% names_best_1.5perc_age_unique]
rf_age_1.5perc_1 = randomForest(x=genotypes_best1.5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_1.5perc_1,file="rf_age_1.5perc_1.Rdata")
rf_age_1.5perc_2 = randomForest(x=genotypes_best1.5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_1.5perc_2,file="rf_age_1.5perc_2.Rdata")
rf_age_1.5perc_3 = randomForest(x=genotypes_best1.5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_1.5perc_3,file="rf_age_1.5perc_3.Rdata")

#load("rf_age_1.5perc_1.Rdata")
#load("rf_age_1.5perc_2.Rdata")
#load("rf_age_1.5perc_3.Rdata")

rf_age_1.5perc_1_err_rate_converge <- rf_age_1.5perc_1$err.rate[750000]
rf_age_1.5perc_2_err_rate_converge <- rf_age_1.5perc_2$err.rate[750000]
rf_age_1.5perc_3_err_rate_converge <- rf_age_1.5perc_3$err.rate[750000]

rm(rf_age_1.5perc_1,rf_age_1.5perc_2,rf_age_1.5perc_3)

##### Best 2% 
names_best_2perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.98))]
names_best_2perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.98))]
names_best_2perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.98))]
names_best_2perc_age_unique<-unique(c(names_best_2perc_age_1,names_best_2perc_age_2,names_best_2perc_age_3))
genotypes_best2perc_age<-All_data[,colnames(All_data) %in% names_best_2perc_age_unique]
rf_age_2perc_1 = randomForest(x=genotypes_best2perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_2perc_1,file="rf_age_2perc_1.Rdata")
rf_age_2perc_2 = randomForest(x=genotypes_best2perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_2perc_2,file="rf_age_2perc_2.Rdata")
rf_age_2perc_3 = randomForest(x=genotypes_best2perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_2perc_3,file="rf_age_2perc_3.Rdata")

#load("rf_age_2perc_1.Rdata")
#load("rf_age_2perc_2.Rdata")
#load("rf_age_2perc_3.Rdata")

rf_age_2perc_1_err_rate_converge <- rf_age_2perc_1$err.rate[750000]
rf_age_2perc_2_err_rate_converge <- rf_age_2perc_2$err.rate[750000]
rf_age_2perc_3_err_rate_converge <- rf_age_2perc_3$err.rate[750000]

rm(rf_age_2perc_1,rf_age_2perc_2,rf_age_2perc_3)

##### Best 3% 
names_best_3perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.97))]
names_best_3perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.97))]
names_best_3perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.97))]
names_best_3perc_age_unique<-unique(c(names_best_3perc_age_1,names_best_3perc_age_2,names_best_3perc_age_3))
genotypes_best3perc_age<-All_data[,colnames(All_data) %in% names_best_3perc_age_unique]
rf_age_3perc_1 = randomForest(x=genotypes_best3perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_3perc_1,file="rf_age_3perc_1.Rdata")
rf_age_3perc_2 = randomForest(x=genotypes_best3perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_3perc_2,file="rf_age_3perc_2.Rdata")
rf_age_3perc_3 = randomForest(x=genotypes_best3perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_3perc_3,file="rf_age_3perc_3.Rdata")

#load("rf_age_3perc_1.Rdata")
#load("rf_age_3perc_2.Rdata")
#load("rf_age_3perc_3.Rdata")

rf_age_3perc_1_err_rate_converge <- rf_age_3perc_1$err.rate[750000]
rf_age_3perc_2_err_rate_converge <- rf_age_3perc_2$err.rate[750000]
rf_age_3perc_3_err_rate_converge <- rf_age_3perc_3$err.rate[750000]

rm(rf_age_3perc_1,rf_age_3perc_2,rf_age_3perc_3)

##### Best 4% 
names_best_4perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.96))]
names_best_4perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.96))]
names_best_4perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.96))]
names_best_4perc_age_unique<-unique(c(names_best_4perc_age_1,names_best_4perc_age_2,names_best_4perc_age_3))
genotypes_best4perc_age<-All_data[,colnames(All_data) %in% names_best_4perc_age_unique]
rf_age_4perc_1 = randomForest(x=genotypes_best4perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_4perc_1,file="rf_age_4perc_1.Rdata")
rf_age_4perc_2 = randomForest(x=genotypes_best4perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_4perc_2,file="rf_age_4perc_2.Rdata")
rf_age_4perc_3 = randomForest(x=genotypes_best4perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_4perc_3,file="rf_age_4perc_3.Rdata")

#load("rf_age_4perc_1.Rdata")
#load("rf_age_4perc_2.Rdata")
#load("rf_age_4perc_3.Rdata")

rf_age_4perc_1_err_rate_converge <- rf_age_4perc_1$err.rate[750000]
rf_age_4perc_2_err_rate_converge <- rf_age_4perc_2$err.rate[750000]
rf_age_4perc_3_err_rate_converge <- rf_age_4perc_3$err.rate[750000]

rm(rf_age_4perc_1,rf_age_4perc_2,rf_age_4perc_3)

##### Best 5% 
names_best_5perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.95))]
names_best_5perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.95))]
names_best_5perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.95))]
names_best_5perc_age_unique<-unique(c(names_best_5perc_age_1,names_best_5perc_age_2,names_best_5perc_age_3))
genotypes_best5perc_age<-All_data[,colnames(All_data) %in% names_best_5perc_age_unique]
rf_age_5perc_1 = randomForest(x=genotypes_best5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_5perc_1,file="rf_age_5perc_1.Rdata")
rf_age_5perc_2 = randomForest(x=genotypes_best5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_5perc_2,file="rf_age_5perc_2.Rdata")
rf_age_5perc_3 = randomForest(x=genotypes_best5perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_5perc_3,file="rf_age_5perc_3.Rdata")

#load("rf_age_5perc_1.Rdata")
#load("rf_age_5perc_2.Rdata")
#load("rf_age_5perc_3.Rdata")

rf_age_5perc_1_err_rate_converge <- rf_age_5perc_1$err.rate[750000]
rf_age_5perc_2_err_rate_converge <- rf_age_5perc_2$err.rate[750000]
rf_age_5perc_3_err_rate_converge <- rf_age_5perc_3$err.rate[750000]

rm(rf_age_5perc_1,rf_age_5perc_2,rf_age_5perc_3)

##### Best 10% 
names_best_10perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.90))]
names_best_10perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.90))]
names_best_10perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.90))]
names_best_10perc_age_unique<-unique(c(names_best_10perc_age_1,names_best_10perc_age_2,names_best_10perc_age_3))
genotypes_best10perc_age<-All_data[,colnames(All_data) %in% names_best_10perc_age_unique]
rf_age_10perc_1 = randomForest(x=genotypes_best10perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_10perc_1,file="rf_age_10perc_1.Rdata")
rf_age_10perc_2 = randomForest(x=genotypes_best10perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_10perc_2,file="rf_age_10perc_2.Rdata")
rf_age_10perc_3 = randomForest(x=genotypes_best10perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_10perc_3,file="rf_age_10perc_3.Rdata")

#load("rf_age_10perc_1.Rdata")
#load("rf_age_10perc_2.Rdata")
#load("rf_age_10perc_3.Rdata")

rf_age_10perc_1_err_rate_converge <- rf_age_10perc_1$err.rate[750000]
rf_age_10perc_2_err_rate_converge <- rf_age_10perc_2$err.rate[750000]
rf_age_10perc_3_err_rate_converge <- rf_age_10perc_3$err.rate[750000]

rm(rf_age_10perc_1,rf_age_10perc_2,rf_age_10perc_3)

##### Best 20% 
names_best_20perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.80))]
names_best_20perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.80))]
names_best_20perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.80))]
names_best_20perc_age_unique<-unique(c(names_best_20perc_age_1,names_best_20perc_age_2,names_best_20perc_age_3))
genotypes_best20perc_age<-All_data[,colnames(All_data) %in% names_best_20perc_age_unique]
rf_age_20perc_1 = randomForest(x=genotypes_best20perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_20perc_1,file="rf_age_20perc_1.Rdata")
rf_age_20perc_2 = randomForest(x=genotypes_best20perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_20perc_2,file="rf_age_20perc_2.Rdata")
rf_age_20perc_3 = randomForest(x=genotypes_best20perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_20perc_3,file="rf_age_20perc_3.Rdata")

#load("rf_age_20perc_1.Rdata")
#load("rf_age_20perc_2.Rdata")
#load("rf_age_20perc_3.Rdata")

rf_age_20perc_1_err_rate_converge <- rf_age_20perc_1$err.rate[750000]
rf_age_20perc_2_err_rate_converge <- rf_age_20perc_2$err.rate[750000]
rf_age_20perc_3_err_rate_converge <- rf_age_20perc_3$err.rate[750000]

rm(rf_age_20perc_1,rf_age_20perc_2,rf_age_20perc_3)

##### Best 30% 
names_best_30perc_age_1<-rownames(importance_rf_age_750000_1)[which(importance_rf_age_750000_1$importance > quantile(importance_rf_age_750000_1$importance, probs=0.70))]
names_best_30perc_age_2<-rownames(importance_rf_age_750000_2)[which(importance_rf_age_750000_2$importance > quantile(importance_rf_age_750000_2$importance, probs=0.70))]
names_best_30perc_age_3<-rownames(importance_rf_age_750000_3)[which(importance_rf_age_750000_3$importance > quantile(importance_rf_age_750000_3$importance, probs=0.70))]
names_best_30perc_age_unique<-unique(c(names_best_30perc_age_1,names_best_30perc_age_2,names_best_30perc_age_3))
genotypes_best30perc_age<-All_data[,colnames(All_data) %in% names_best_30perc_age_unique]
rf_age_30perc_1 = randomForest(x=genotypes_best30perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_30perc_1,file="rf_age_30perc_1.Rdata")
rf_age_30perc_2 = randomForest(x=genotypes_best30perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_30perc_2,file="rf_age_30perc_2.Rdata")
rf_age_30perc_3 = randomForest(x=genotypes_best30perc_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=750000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_30perc_3,file="rf_age_30perc_3.Rdata")

#load("rf_age_30perc_1.Rdata")
#load("rf_age_30perc_2.Rdata")
#load("rf_age_30perc_3.Rdata")

rf_age_30perc_1_err_rate_converge <- rf_age_30perc_1$err.rate[750000]
rf_age_30perc_2_err_rate_converge <- rf_age_30perc_2$err.rate[750000]
rf_age_30perc_3_err_rate_converge <- rf_age_30perc_3$err.rate[750000]

rm(rf_age_30perc_1,rf_age_30perc_2,rf_age_30perc_3)

## Combine initial out of bag error rates for age
All_initial_err_rate_age<-rbind(cbind(rf_age_750000_1_err_rate_converge,rf_age_750000_2_err_rate_converge,rf_age_750000_3_err_rate_converge),
                                cbind(rf_age_0.5perc_1_err_rate_converge,rf_age_0.5perc_2_err_rate_converge,rf_age_0.5perc_3_err_rate_converge),
                                cbind(rf_age_1perc_1_err_rate_converge,rf_age_1perc_2_err_rate_converge,rf_age_1perc_3_err_rate_converge),
                                cbind(rf_age_1.5perc_1_err_rate_converge,rf_age_1.5perc_2_err_rate_converge,rf_age_1.5perc_3_err_rate_converge),
                                cbind(rf_age_2perc_1_err_rate_converge,rf_age_2perc_2_err_rate_converge,rf_age_2perc_3_err_rate_converge),
                                cbind(rf_age_3perc_1_err_rate_converge,rf_age_3perc_2_err_rate_converge,rf_age_3perc_3_err_rate_converge),
                                cbind(rf_age_4perc_1_err_rate_converge,rf_age_4perc_2_err_rate_converge,rf_age_4perc_3_err_rate_converge),
                                cbind(rf_age_5perc_1_err_rate_converge,rf_age_5perc_2_err_rate_converge,rf_age_5perc_3_err_rate_converge),
                                cbind(rf_age_10perc_1_err_rate_converge,rf_age_10perc_2_err_rate_converge,rf_age_10perc_3_err_rate_converge),
                                cbind(rf_age_20perc_1_err_rate_converge,rf_age_20perc_2_err_rate_converge,rf_age_20perc_3_err_rate_converge),
                                cbind(rf_age_30perc_1_err_rate_converge,rf_age_30perc_2_err_rate_converge,rf_age_30perc_3_err_rate_converge))

#age
All_initial_err_rate_age<-data.frame(All_initial_err_rate_age)
All_initial_err_rate_age$Number_loci<-c(9108,length(names_best_0.5perc_age_unique),length(names_best_1perc_age_unique),length(names_best_1.5perc_age_unique),length(names_best_2perc_age_unique),length(names_best_3perc_age_unique),length(names_best_4perc_age_unique),length(names_best_5perc_age_unique),length(names_best_10perc_age_unique),length(names_best_20perc_age_unique),length(names_best_30perc_age_unique))
rownames(All_initial_err_rate_age)<-c("All","Best0.5%","Best1%","Best1.5%","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_err_rate_age$Average<-apply(All_initial_err_rate_age[,1:3],1,mean)
write.csv(All_initial_err_rate_age,file="All_initial_err_rate_age_geno_corrected_and_balanced.csv")

par(mar=c(5,5,3,3))
plot(All_initial_err_rate_age$Number_loci,All_initial_err_rate_age$Average,log="x",
     main="Out of Bag Error Rate for Age-RF 750000 Trees", xlab="Number of Loci", ylab="OOB Error Rate",
     cex.lab=1.5,cex.axis=1.5,pch=19)

#Based on this table and plot, the best 2% of loci has the lowest out of bag error rate (error is 29.2 %) for age at maturity
# The best 2% data set includes 211 unique loci, so I'll run backward purging RF with the best 3% loci (n=324 unique)

#################### Backward purging approach with top 3% of loci - age

names_best_age_unique <- names_best_3perc_age_unique
length(names_best_age_unique)

names_best_age324 <- names_best_age_unique
genotypes_purging_age<-All_data[,colnames(All_data) %in% names_best_age324]
rf_age_purging_1 = randomForest(x=genotypes_purging_age, y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=250000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_purging_1,file="rf_age_purging_1_250k.Rdata")
rf_age_purging_2 = randomForest(x=genotypes_purging_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=250000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_purging_2,file="rf_age_purging_2_250k.Rdata")

#Extract importance (mean decrease in accuracy) measure of each locus for age
importance_rf_age_purging_1<-data.frame(importance(rf_age_purging_1,type=1))
colnames(importance_rf_age_purging_1)<-c("Importance1")
importance_rf_age_purging_2<-data.frame(importance(rf_age_purging_2,type=1))
colnames(importance_rf_age_purging_2)<-c("Importance2")
importance_rf_purging <-cbind(rownames(importance_rf_age_purging_1),importance_rf_age_purging_1,importance_rf_age_purging_2)

#Plot importance values against each other to check if model converged
plot(importance_rf_purging$Importance1,importance_rf_purging$Importance2)
fit<-lm(Importance2~Importance1, data=importance_rf_purging)
summary(fit)
abline(fit)

rf_age_purging_3 = randomForest(x=genotypes_purging_age,y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=250000, strata=as.factor(All_data$age), sampsize=sample_size)
save(rf_age_purging_3,file="rf_age_purging_3_250k.Rdata")

#load("rf_age_purging_1.Rdata")
#load("rf_age_purging_2.Rdata")
#load("rf_age_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[324]]<-names_best_age324
error_rate_age_best<-data.frame(V1=1:324,V2=1:324,V3=1:324)
rownames(error_rate_age_best)<-1:324
error_rate_age_best[324,]<-c(rf_age_purging_1$err.rate[250000],rf_age_purging_2$err.rate[250000],rf_age_purging_3$err.rate[250000])

for (i in 1:323){
  print(i)
  imp_purging_1<-data.frame(importance(rf_age_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_age_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_age_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging_age)
  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  rownames(imp_purging_2)<-colnames(genotypes_purging_age)
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  rownames(imp_purging_3)<-colnames(genotypes_purging_age)
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[324-i]]<-names_keep
  genotypes_purging_age<-All_data[,colnames(All_data) %in% names_keep]
  rf_age_purging_1 = randomForest(x = genotypes_purging_age, y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=250000, strata=as.factor(All_data$age), sampsize=sample_size)
  rf_age_purging_2 = randomForest(x = genotypes_purging_age, y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=250000, strata=as.factor(All_data$age), sampsize=sample_size)
  rf_age_purging_3 = randomForest(x = genotypes_purging_age, y = as.factor(All_data$age), importance=TRUE ,proximity=TRUE, ntree=250000, strata=as.factor(All_data$age), sampsize=sample_size)
  error_rate_age_best[324-i,]<-c(rf_age_purging_1$err.rate[250000],rf_age_purging_2$err.rate[250000],rf_age_purging_3$err.rate[250000])
  
  write.csv(error_rate_age_best, file="Error_rate_return_best.csv")
  lapply(names_all_iterations, write, "Backward_purging_age_names_all_iterations_temp.txt", append=TRUE, ncolumns=328)
  
}

error_rate_age_best$Average<-apply(error_rate_age_best,1,mean)

write.csv(error_rate_age_best, file="Error_rate_return_best.csv")

lapply(names_all_iterations, write, "Backward_purging_age_names_all_iterations.txt", append=TRUE, ncolumns=328)

par(mar=c(5,5,3,3))
plot(error_rate_age_best$Average[-c(1)],xlab="Number of Loci", ylab="Out of Bag Error Rate",
     main="Backward Purging Age",cex.lab=1.5,cex.axis=1.5,pch=19)

which(error_rate_age_best$Average==min(error_rate_age_best$Average[-c(1)])) # We want the smallest error rate, not the largest

write.csv(names_all_iterations[[30]],file="Purging_optimum_age_750000.csv") #Insert index of max variation




##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
