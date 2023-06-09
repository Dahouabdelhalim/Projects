# set the working directory
setwd("/Users/msmith/Box Sync/MSmith/In Progress - Synced/PredictivePhylgeography_withTraitData/Results/SAL_August08")

# get the data
SALData <- read.table("/Users/msmith/Box Sync/MSmith/In Progress - Synced/PredictivePhylgeography_withTraitData/Data/SAL_MLS_August08_2018_FullDataset_forRF.csv", sep=",", header=T)

# call the library
library(randomForest)

#look at the data
names(SALData)

# get rid of Neotoma and fuscus
SALData <- SALData[which(SALData$complex!= "Neotoma"),]
SALData <- SALData[which(SALData$complex!= "Pipilo fuscus"),]
levels(SALData$complex)
table(SALData$complex)
SALData$complex <- droplevels(SALData$complex, exclude = c("Neotoma", "Pipilo fucus"))

nrow(SALData)
# check out the data
table(SALData$complex, SALData$Group)
table(SALData$complex, SALData$Diet_5Cat)
table(SALData$complex, SALData$Nocturnal)
table(SALData$complex, SALData$BodyMass)
table(SALData$complex, SALData$ClutchSize)
table(SALData$complex, SALData$Reproduction)
table(SALData$complex, SALData$Class)

class(SALData$Group)
class(SALData$Diet_5Cat)
range(SALData$Nocturnal)
SALData$Nocturnal <- as.factor(as.character(SALData$Nocturnal))
class(SALData$Nocturnal)
levels(SALData$Nocturnal)
class(SALData$BodyMass)
class(SALData$ClutchSize)
class(SALData$Reproduction)
class(SALData$Class)

# build the model with only bioclimatic variables and taxonomy
NoTraits_Taxonomy <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                           + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                           + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                           + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                           + Class, data=SALData, ntree=500, importance = TRUE)
notraits_taxonomy <- NoTraits_Taxonomy$err.rate[500,1]

pdf('NoTraits_Taxonomy.pdf')
varImpPlot(NoTraits_Taxonomy)
dev.off()

# look at the model
NoTraits_Taxonomy

# build the model with all traits and taxonomy
AllTraits_Taxonomy <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                           + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                           + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                           + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                           + BodyMass + Diet_5Cat + Nocturnal  + ClutchSize   + Reproduction + Class, data=SALData, ntree=500, importance = TRUE)
alltraits_taxonomy <- AllTraits_Taxonomy$err.rate[500,1]

# look at the model
AllTraits_Taxonomy

# plot MDA and GINI
pdf('AllTraits_Taxonomy.pdf')
varImpPlot(AllTraits_Taxonomy)
dev.off()

# build the model with all traits and no taxonomy
AllTraits_NoTaxonomy <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                           + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                           + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                           + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                           + BodyMass + Diet_5Cat + Nocturnal  + ClutchSize   + Reproduction, data=SALData, ntree=500, importance = TRUE)

alltraits_notaxonomy <- AllTraits_NoTaxonomy$err.rate[500,1]

# look at the model
AllTraits_NoTaxonomy

# plot MDA and GINI
pdf('AllTraits_NoTaxonomy.pdf')
varImpPlot(AllTraits_NoTaxonomy)
dev.off()


# Cross-validation with all information, without taxonomy, and without diet and body mass.
Crossvalresults_1 <- data.frame()
Crossvalresults_2 <- data.frame()
for (item in levels(SALData$complex)){
  SALData_cval_training <- SALData[which(SALData$complex != item),] #create dataframe for training
  SALData_cval_testing <- SALData[which(SALData$complex == item),] #create dataframe for testing
  RF_model_1_cval <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                           + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                           + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                           + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                           + BodyMass + Diet_5Cat + Nocturnal  + ClutchSize   + Reproduction + Class, data=SALData_cval_training, ntree=500)
  RF_model_2_cval <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                                + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                                + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                                + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                                + BodyMass + Diet_5Cat + Nocturnal  + ClutchSize   + Reproduction , data=SALData_cval_training, ntree=500)
  Prediction_1 <- predict(RF_model_1_cval, newdata = SALData_cval_testing, type = "response")
  MyPrediction1 <- cbind(Prediction_1, SALData_cval_testing)
  Prediction_2 <- predict(RF_model_2_cval, newdata = SALData_cval_testing, type = "response")
  MyPrediction2 <- cbind(Prediction_2, SALData_cval_testing)
  Crossvalresults_1 <- rbind(Crossvalresults_1,MyPrediction1)
  Crossvalresults_2 <- rbind(Crossvalresults_2,MyPrediction2)
}

AllTraits_Taxonomy_truecryptic = Crossvalresults_1[Crossvalresults_1$Group == 'C' & Crossvalresults_1$Prediction_1 == 'C',]
AllTraits_Taxonomy_truenoncryptic = Crossvalresults_1[Crossvalresults_1$Group == 'NC' & Crossvalresults_1$Prediction_1 == 'NC',]
AllTraits_Taxonomy_totalcryptic = Crossvalresults_1[Crossvalresults_1$Prediction_1 == 'C',]
AllTraits_Taxonomy_totalnoncryptic = Crossvalresults_1[Crossvalresults_1$Prediction_1 == 'NC',]
AllTraits_Taxonomy_accuracyoverall = (nrow(AllTraits_Taxonomy_truecryptic)+nrow(AllTraits_Taxonomy_truenoncryptic))/nrow(Crossvalresults_1)
AllTraits_Taxonomy_accuracycryptic = (nrow(AllTraits_Taxonomy_truecryptic)/nrow(AllTraits_Taxonomy_totalcryptic))
AllTraits_Taxonomy_accuracynoncryptic = (nrow(AllTraits_Taxonomy_truenoncryptic)/nrow(AllTraits_Taxonomy_totalnoncryptic))
AllTraits_Taxonomy_accuracyoverall
AllTraits_Taxonomy_accuracycryptic
AllTraits_Taxonomy_accuracynoncryptic

AllTraits_NoTaxonomy_truecryptic = Crossvalresults_2[Crossvalresults_2$Group == 'C' & Crossvalresults_2$Prediction_2 == 'C',]
AllTraits_NoTaxonomy_truenoncryptic = Crossvalresults_2[Crossvalresults_2$Group == 'NC' & Crossvalresults_2$Prediction_2 == 'NC',]
AllTraits_NoTaxonomy_totalcryptic = Crossvalresults_2[Crossvalresults_2$Prediction_2 == 'C',]
AllTraits_NoTaxonomy_totalnoncryptic = Crossvalresults_2[Crossvalresults_2$Prediction_2 == 'NC',]
AllTraits_NoTaxonomy_accuracyoverall = (nrow(AllTraits_NoTaxonomy_truecryptic)+nrow(AllTraits_NoTaxonomy_truenoncryptic))/nrow(Crossvalresults_2)
AllTraits_NoTaxonomy_accuracycryptic = (nrow(AllTraits_NoTaxonomy_truecryptic)/nrow(AllTraits_NoTaxonomy_totalcryptic))
AllTraits_NoTaxonomy_accuracynoncryptic = (nrow(AllTraits_NoTaxonomy_truenoncryptic)/nrow(AllTraits_NoTaxonomy_totalnoncryptic))
AllTraits_NoTaxonomy_accuracyoverall
AllTraits_NoTaxonomy_accuracycryptic
AllTraits_NoTaxonomy_accuracynoncryptic


# variable importance

Crossvalresults_5 <- data.frame()
SALdata_forRF <- SALData[,c(2:20, 26:31, 25, 21)]
names(SALdata_forRF)
traitlist <- c("Diet_5Cat", "Nocturnal", "BodyMass" 
               , "ClutchSize", "Reproduction", "Class")


combos_1<- combn(traitlist, m=1,simplyify=F)
combos_2<- combn(traitlist, m=2,simplyify=F)
combos_3<- combn(traitlist, m=3,simplyify=F)
combos_4<- combn(traitlist, m=4,simplyify=F)
combos_5<- combn(traitlist, m=5,simplyify=F)

all_combos <- list(combos_1, combos_2)


for (taxon in levels(SALdata_forRF$complex)){
  for (item in all_combos){
    mycombos <- item
    for (i in 1:range(ncol(item))){
      theitems <- mycombos[,i]
      SAL_corvars <- SALdata_forRF[,!(names(SALdata_forRF) %in% theitems)]
      toremove <- list('complex')
      SALdata_cval_training <- SAL_corvars[which(SAL_corvars$complex != taxon),] #create dataframe for training
      SALdata_cval_testing <- SAL_corvars[which(SAL_corvars$complex == taxon),] #create dataframe for testing
      SALdata_cval_testing <- SALdata_cval_testing[,!(names(SALdata_cval_testing) %in% toremove)]
      SALdata_cval_training <- SALdata_cval_training[,!(names(SALdata_cval_training) %in% toremove)]
      RF_model_5_cval <- randomForest(Group~., data = SALdata_cval_training, 
                                          ntree = 500)
      Prediction_5 <- predict(RF_model_5_cval, newdata = SALdata_cval_testing, type = "response")
      theVar <- rep(paste(theitems, collapse = "_"), length(Prediction_5))
      MyPrediction5 <- cbind(V1 = as.vector(Prediction_5), V2 = as.vector(SALdata_cval_testing$Group),theVar)
      Crossvalresults_5 <- rbind(Crossvalresults_5,MyPrediction5)
    }
  }
}

MyVariableImportance <- data.frame()
for (item in all_combos){
  mycombos <- item
  for (i in 1:range(ncol(item))){
    theitems <- mycombos[,i]
    theitems <- paste(theitems, collapse = "_")
    relevant <- Crossvalresults_5[which(Crossvalresults_5$theVar == theitems),]
    truecryptic = relevant[relevant$V2 == 'C' & relevant$V1 == 'C',]
    truenoncryptic = relevant[relevant$V2 == 'NC' & relevant$V1 == 'NC',]
    totalcryptic = relevant[relevant$V1 == 'C',]
    totalnoncryptic = relevant[relevant$V1 == 'NC',]
    accuracyoverall = (nrow(truecryptic)+nrow(truenoncryptic))/(nrow(totalcryptic)+nrow(totalnoncryptic))
    accuracycryptic = (nrow(truecryptic)/nrow(totalcryptic))
    accuracynoncryptic = (nrow(truenoncryptic)/nrow(totalnoncryptic))
    VariableImportance <- cbind.data.frame(accuracyoverall,accuracycryptic,accuracynoncryptic,theitems)
    MyVariableImportance <- rbind.data.frame(MyVariableImportance, VariableImportance)
  }
}

MyVariableImportance$MDA <- AllTraits_Taxonomy_accuracyoverall - MyVariableImportance$accuracyoverall
par(mai=c(1,2,1,1))
barplot(MyVariableImportance$MDA, horiz = TRUE, names.arg = MyVariableImportance$theitems, axes = TRUE, las=1)

write.table(MyVariableImportance, 'MyVariableImportance_SAL.csv', sep = ',')  



#leaving out body mass and nocturnal, curated traits + taxonomy model

CuratedTraits_Taxonomy <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                           + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                           + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                           + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                           + Diet_5Cat + Reproduction + ClutchSize + Class, data=SALData, ntree=500, importance = TRUE)

curatedtraits_taxonomy <- CuratedTraits_Taxonomy$err.rate[500,1]

# look at the model
CuratedTraits_Taxonomy

# plot MDA and GINI
pdf('CuratedTraits_Taxonomy.pdf')
varImpPlot(CuratedTraits_Taxonomy)
dev.off()


# build the model with the best traits, but without taxonomy
CuratedTraits_NoTaxonomy<- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                                          + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                                          + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                                          + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                                          + Diet_5Cat + Reproduction + ClutchSize, data=SALData, ntree=500, importance = TRUE)

curatedtraits_notaxonomy <- CuratedTraits_NoTaxonomy$err.rate[500,1]

# look at the model
CuratedTraits_NoTaxonomy

# plot MDA and GINI
pdf('CuratedTraits_NoTaxonomy.pdf')
varImpPlot(CuratedTraits_NoTaxonomy)
dev.off()




# Do cross validation for hte model with the best traits + taxonomy
Crossvalresults_6 <- data.frame()
for (item in levels(SALData$complex)){
  SALData_cval_training <- SALData[which(SALData$complex != item),] #create dataframe for training
  SALData_cval_testing <- SALData[which(SALData$complex == item),] #create dataframe for testing
  RF_model_6_cval <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                                  + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                                  + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                                  + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                                  + Diet_5Cat + Reproduction + ClutchSize + Class , data=SALData_cval_training, ntree=500)
  Prediction_6 <- predict(RF_model_6_cval, newdata = SALData_cval_testing, type = "response")
  Prediction_6b <- predict(RF_model_6_cval, newdata = SALData_cval_testing, type = "prob")
  MyPrediction6 <- cbind(Prediction_6, SALData_cval_testing, Prediction_6b)
  Crossvalresults_6 <- rbind(Crossvalresults_6,MyPrediction6)
}

CuratedTraits_Taxonomy_truecryptic = Crossvalresults_6[Crossvalresults_6$Group == 'C' & Crossvalresults_6$Prediction_6 == 'C',]
CuratedTraits_Taxonomy_truenoncryptic = Crossvalresults_6[Crossvalresults_6$Group == 'NC' & Crossvalresults_6$Prediction_6 == 'NC',]
CuratedTraits_Taxonomy_totalcryptic = Crossvalresults_6[Crossvalresults_6$Prediction_6 == 'C',]
CuratedTraits_Taxonomy_totalnoncryptic = Crossvalresults_6[Crossvalresults_6$Prediction_6 == 'NC',]
CuratedTraits_Taxonomy_accuracyoverall = (nrow(CuratedTraits_Taxonomy_truecryptic)+nrow(CuratedTraits_Taxonomy_truenoncryptic))/nrow(Crossvalresults_6)
CuratedTraits_Taxonomy_accuracycryptic = (nrow(CuratedTraits_Taxonomy_truecryptic)/nrow(CuratedTraits_Taxonomy_totalcryptic))
CuratedTraits_Taxonomy_accuracynoncryptic = (nrow(CuratedTraits_Taxonomy_truenoncryptic)/nrow(CuratedTraits_Taxonomy_totalnoncryptic))
CuratedTraits_Taxonomy_accuracyoverall
CuratedTraits_Taxonomy_accuracycryptic
CuratedTraits_Taxonomy_accuracynoncryptic

# Do cross validation for the model with bioclimatic variables and taxonomy only
Crossvalresults_7 <- data.frame()
for (item in levels(SALData$complex)){
  SALData_cval_training <- SALData[which(SALData$complex != item),] #create dataframe for training
  SALData_cval_testing <- SALData[which(SALData$complex == item),] #create dataframe for testing
  RF_model_7_cval <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                                  + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                                  + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                                  + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                                  + Class , data=SALData_cval_training, ntree=500)
  Prediction_7 <- predict(RF_model_7_cval, newdata = SALData_cval_testing, type = "response")
  Prediction_7b <- predict(RF_model_7_cval, newdata = SALData_cval_testing, type = "prob")
  MyPrediction7 <- cbind(Prediction_7, SALData_cval_testing, Prediction_7b)
  Crossvalresults_7 <- rbind(Crossvalresults_7,MyPrediction7)
}


NoTraits_Taxonomy_truecryptic = Crossvalresults_7[Crossvalresults_7$Group == 'C' & Crossvalresults_7$Prediction_7 == 'C',]
NoTraits_Taxonomy_truenoncryptic = Crossvalresults_7[Crossvalresults_7$Group == 'NC' & Crossvalresults_7$Prediction_7 == 'NC',]
NoTraits_Taxonomy_totalcryptic = Crossvalresults_7[Crossvalresults_7$Prediction_7 == 'C',]
NoTraits_Taxonomy_totalnoncryptic = Crossvalresults_7[Crossvalresults_7$Prediction_7 == 'NC',]
NoTraits_Taxonomy_accuracyoverall = (nrow(NoTraits_Taxonomy_truecryptic)+nrow(NoTraits_Taxonomy_truenoncryptic))/nrow(Crossvalresults_7)
NoTraits_Taxonomy_accuracycryptic = (nrow(NoTraits_Taxonomy_truecryptic)/nrow(NoTraits_Taxonomy_totalcryptic))
NoTraits_Taxonomy_accuracynoncryptic = (nrow(NoTraits_Taxonomy_truenoncryptic)/nrow(NoTraits_Taxonomy_totalnoncryptic))
NoTraits_Taxonomy_accuracyoverall
NoTraits_Taxonomy_accuracycryptic
NoTraits_Taxonomy_accuracynoncryptic

# Do cross validation for the model with the best traits and no taxonomy
Crossvalresults_9 <- data.frame()
for (item in levels(SALData$complex)){
  SALData_cval_training <- SALData[which(SALData$complex != item),] #create dataframe for training
  SALData_cval_testing <- SALData[which(SALData$complex == item),] #create dataframe for testing
  RF_model_9_cval <- randomForest(Group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03 + wc2.0_bio_30s_04 
                                  + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 + wc2.0_bio_30s_08 + wc2.0_bio_30s_09
                                  + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14
                                  + wc2.0_bio_30s_15 + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19
                                  + Diet_5Cat + Reproduction + ClutchSize  , data=SALData_cval_training, ntree=500)
  Prediction_9 <- predict(RF_model_9_cval, newdata = SALData_cval_testing, type = "response")
  Prediction_9b <- predict(RF_model_9_cval, newdata = SALData_cval_testing, type = "prob")
  MyPrediction9 <- cbind(Prediction_9, SALData_cval_testing, Prediction_9b)
  Crossvalresults_9 <- rbind(Crossvalresults_9,MyPrediction9)
}

CuratedTraits_NoTaxonomy_truecryptic = Crossvalresults_9[Crossvalresults_9$Group == 'C' & Crossvalresults_9$Prediction_9 == 'C',]
CuratedTraits_NoTaxonomy_truenoncryptic = Crossvalresults_9[Crossvalresults_9$Group == 'NC' & Crossvalresults_9$Prediction_9 == 'NC',]
CuratedTraits_NoTaxonomy_totalcryptic = Crossvalresults_9[Crossvalresults_9$Prediction_9 == 'C',]
CuratedTraits_NoTaxonomy_totalnoncryptic = Crossvalresults_9[Crossvalresults_9$Prediction_9 == 'NC',]
CuratedTraits_NoTaxonomy_accuracyoverall = (nrow(CuratedTraits_NoTaxonomy_truecryptic)+nrow(CuratedTraits_NoTaxonomy_truenoncryptic))/nrow(Crossvalresults_9)
CuratedTraits_NoTaxonomy_accuracycryptic = (nrow(CuratedTraits_NoTaxonomy_truecryptic)/nrow(CuratedTraits_NoTaxonomy_totalcryptic))
CuratedTraits_NoTaxonomy_accuracynoncryptic = (nrow(CuratedTraits_NoTaxonomy_truenoncryptic)/nrow(CuratedTraits_NoTaxonomy_totalnoncryptic))
CuratedTraits_NoTaxonomy_accuracyoverall
CuratedTraits_NoTaxonomy_accuracycryptic
CuratedTraits_NoTaxonomy_accuracynoncryptic


# plot the results
Ammospermophilus_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Ammospermophilus leucurus",c(33:34)]
Ammospermophilus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Ammospermophilus leucurus",c(33:34)]
Ammospermophilus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Ammospermophilus leucurus",c(33:34)]
Ammospermophilus_cv6$RF <- 'both'
Ammospermophilus_cv7$RF <- 'taxonomy'
Ammospermophilus_cv9$RF <- 'traits'
Ammospermophilus_cv <- rbind(Ammospermophilus_cv6,Ammospermophilus_cv7,Ammospermophilus_cv9)

Anaxyrus_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Anaxyrus.punctatus.group",c(33:34)]
Anaxyrus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Anaxyrus.punctatus.group",c(33:34)]
Anaxyrus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Anaxyrus.punctatus.group",c(33:34)]
Anaxyrus_cv6$RF <- 'both'
Anaxyrus_cv7$RF <- 'taxonomy'
Anaxyrus_cv9$RF <- 'traits'
Anaxyrus_cv <- rbind(Anaxyrus_cv6,Anaxyrus_cv7,Anaxyrus_cv9)

Auriparus_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Auriparus flaviceps",c(33:34)]
Auriparus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Auriparus flaviceps",c(33:34)]
Auriparus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Auriparus flaviceps",c(33:34)]
Auriparus_cv6$RF <- 'both'
Auriparus_cv7$RF <- 'taxonomy'
Auriparus_cv9$RF <- 'traits'
Auriparus_cv <- rbind(Auriparus_cv6,Auriparus_cv7,Anaxyrus_cv9)


Callilepla_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Callilepla.californica-gambelli",c(33:34)]
Callilepla_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Callilepla.californica-gambelli",c(33:34)]
Callilepla_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Callilepla.californica-gambelli",c(33:34)]
Callilepla_cv6$RF <- 'both'
Callilepla_cv7$RF <- 'taxonomy'
Callilepla_cv9$RF <- 'traits'
Callilepla_cv <- rbind(Callilepla_cv6,Callilepla_cv7,Callilepla_cv9)

Campylorhynchus_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Campylorhynchus brunneicapillus",c(33:34)]
Campylorhynchus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Campylorhynchus brunneicapillus",c(33:34)]
Campylorhynchus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Campylorhynchus brunneicapillus",c(33:34)]
Campylorhynchus_cv6$RF <- 'both'
Campylorhynchus_cv7$RF <- 'taxonomy'
Campylorhynchus_cv9$RF <- 'traits'
Campylorhynchus_cv <- rbind(Campylorhynchus_cv6,Campylorhynchus_cv7,Campylorhynchus_cv9)

Chaetodipus_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Chaetodipus.baileyi.group",c(33:34)]
Chaetodipus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Chaetodipus.baileyi.group",c(33:34)]
Chaetodipus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Chaetodipus.baileyi.group",c(33:34)]
Chaetodipus_cv6$RF <- 'both'
Chaetodipus_cv7$RF <- 'taxonomy'
Chaetodipus_cv9$RF <- 'traits'
Chaetodipus_cv <- rbind(Chaetodipus_cv6,Chaetodipus_cv7,Chaetodipus_cv9)

Dipodomys_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Dipodomys.merriami.group",c(33:34)]
Dipodomys_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Dipodomys.merriami.group",c(33:34)]
Dipodomys_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Dipodomys.merriami.group",c(33:34)]
Dipodomys_cv6$RF <- 'both'
Dipodomys_cv7$RF <- 'taxonomy'
Dipodomys_cv9$RF <- 'traits'
Dipodomys_cv <- rbind(Dipodomys_cv6,Dipodomys_cv7,Dipodomys_cv9)

Peromyscus_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Peromyscus.group",c(33:34)]
Peromyscus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Peromyscus.group",c(33:34)]
Peromyscus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Peromyscus.group",c(33:34)]
Peromyscus_cv6$RF <- 'both'
Peromyscus_cv7$RF <- 'taxonomy'
Peromyscus_cv9$RF <- 'traits'
Peromyscus_cv <- rbind(Peromyscus_cv6,Peromyscus_cv7,Peromyscus_cv9)

AbCris_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Pi.aberti-crissalis",c(33:34)]
AbCris_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Pi.aberti-crissalis",c(33:34)]
AbCris_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Pi.aberti-crissalis",c(33:34)]
AbCris_cv6$RF <- 'both'
AbCris_cv7$RF <- 'taxonomy'
AbCris_cv9$RF <- 'traits'
AbCris_cv <- rbind(AbCris_cv6,AbCris_cv7,AbCris_cv9)


MelCal_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Po.melanura-californica",c(33:34)]
MelCal_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Po.melanura-californica",c(33:34)]
MelCal_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Po.melanura-californica",c(33:34)]
MelCal_cv6$RF <- 'both'
MelCal_cv7$RF <- 'taxonomy'
MelCal_cv9$RF <- 'traits'
MelCal_cv <- rbind(MelCal_cv6,MelCal_cv7,MelCal_cv9)

CinBin_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "T.cinereum-binderei",c(33:34)]
CinBin_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "T.cinereum-binderei",c(33:34)]
CinBin_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "T.cinereum-binderei",c(33:34)]
CinBin_cv6$RF <- 'both'
CinBin_cv7$RF <- 'taxonomy'
CinBin_cv9$RF <- 'traits'
CinBin_cv <- rbind(CinBin_cv6,CinBin_cv7,CinBin_cv9)

Lecontei_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "Toxostoma lecontei",c(33:34)]
Lecontei_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Toxostoma lecontei",c(33:34)]
Lecontei_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Toxostoma lecontei",c(33:34)]
Lecontei_cv6$RF <- 'both'
Lecontei_cv7$RF <- 'taxonomy'
Lecontei_cv9$RF <- 'traits'
Lecontei_cv <- rbind(Lecontei_cv6,Lecontei_cv7,Lecontei_cv9)

LucRuf_cv6 <- Crossvalresults_6[Crossvalresults_6$complex == "V.luciae-ruficapilla",c(33:34)]
LucRuf_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "V.luciae-ruficapilla",c(33:34)]
LucRuf_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "V.luciae-ruficapilla",c(33:34)]
LucRuf_cv6$RF <- 'both'
LucRuf_cv7$RF <- 'taxonomy'
LucRuf_cv9$RF <- 'traits'
LucRuf_cv <- rbind(LucRuf_cv6,LucRuf_cv7,LucRuf_cv9)

library(ggplot2)
Ammospermophilusplot <-   ggplot(Ammospermophilus_cv,    aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..))  + guides(fill=FALSE) + xlim(0, 1) + labs(title = 'Ammospermophilus') + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Anaxyrusplot <-           ggplot(Anaxyrus_cv,            aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Anaxyrus') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Auriparusplot <-          ggplot(Auriparus_cv,           aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Auriparus') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Callileplaplot <-         ggplot(Callilepla_cv,          aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Callilepla') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Campylorhynchusplot <-    ggplot(Campylorhynchus_cv,     aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Campylorhynchus') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Chaetodipusplot <-        ggplot(Chaetodipus_cv,         aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Chaetodipus') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Dipodomysplot <-          ggplot(Dipodomys_cv,           aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Dipodomys') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Peromyscusplot <-         ggplot(Peromyscus_cv,          aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Peromyscus') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
AbCrisplot <-             ggplot(AbCris_cv,              aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Pipilo aberti/ crissalis') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
MelCalplot <-             ggplot(MelCal_cv,              aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Polioptila melanura/ californica') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
CinBinplot <-             ggplot(CinBin_cv,              aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Toxostoma cinereum/ binderei') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
Leconteiplot <-           ggplot(Lecontei_cv,            aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Toxostoma lecontei') + theme(plot.title = element_text(hjust = 0.5))+ scale_fill_manual( values = c("red","blue", "yellow"))
LucRufplot <-              ggplot(LucRuf_cv,             aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE)+ labs(title = 'Vermivora luciae/ ruficapilla') + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))

library(gridExtra)
pdf("SALComparisons.pdf")
grid.arrange(Ammospermophilusplot, Anaxyrusplot, Auriparusplot, 
             Callileplaplot, Campylorhynchusplot, Chaetodipusplot, 
             Dipodomysplot, Peromyscusplot, 
             AbCrisplot, MelCalplot, 
             CinBinplot, Leconteiplot, LucRufplot, nrow = 5)
dev.off()

Ammospermophilusplotlegend <- ggplot(Ammospermophilus_cv, aes(C, fill = RF)) + geom_density(alpha = 0.2)  + xlim(0, 1) + labs(title = 'Ammospermophilus') + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
pdf("Ammospermophilusforlegend.pdf")
grid.arrange(Ammospermophilusplotlegend)
dev.off()

table(SALData$complex,SALData$Group)

oobrates <- cbind(notraits_taxonomy,
                  alltraits_taxonomy,
                  alltraits_notaxonomy,
                  curatedtraits_taxonomy,
                  curatedtraits_notaxonomy,
                  notraits_taxonomy)

write.csv(oobrates, "OOBRates.csv")

predictionaccuracies <- cbind(AllTraits_Taxonomy_accuracyoverall,
                              AllTraits_Taxonomy_accuracycryptic,
                              AllTraits_Taxonomy_accuracynoncryptic,
                              AllTraits_NoTaxonomy_accuracyoverall,
                              AllTraits_NoTaxonomy_accuracycryptic,
                              AllTraits_NoTaxonomy_accuracynoncryptic,
                              CuratedTraits_Taxonomy_accuracyoverall,
                              CuratedTraits_Taxonomy_accuracycryptic,
                              CuratedTraits_Taxonomy_accuracynoncryptic,
                              NoTraits_Taxonomy_accuracyoverall,
                              NoTraits_Taxonomy_accuracycryptic,
                              NoTraits_Taxonomy_accuracynoncryptic,
                              CuratedTraits_NoTaxonomy_accuracyoverall,
                              CuratedTraits_NoTaxonomy_accuracycryptic,
                              CuratedTraits_NoTaxonomy_accuracynoncryptic)

write.csv(predictionaccuracies, "PredictionAccuracies.csv")


# test how variable importance changes when the two bad traits are omitted.
# variable importance

Crossvalresults_traits2 <- data.frame()
SALdata_forRF_traits2 <- SALData[,c(2:20, 26, 29:31, 25, 21)]
names(SALdata_forRF_traits2)
traitlist <- c("Diet_5Cat" , "ClutchSize", "Reproduction", "Class")


combos_1<- combn(traitlist, m=1,simplyify=F)

all_combos <- list(combos_1)


for (taxon in levels(SALdata_forRF_traits2$complex)){
  for (item in all_combos){
    mycombos <- item
    for (i in 1:range(ncol(item))){
      theitems <- mycombos[,i]
      SAL_corvars <- SALdata_forRF_traits2[,!(names(SALdata_forRF_traits2) %in% theitems)]
      toremove <- list('complex')
      SALdata_cval_training <- SAL_corvars[which(SAL_corvars$complex != taxon),] #create dataframe for training
      SALdata_cval_testing <- SAL_corvars[which(SAL_corvars$complex == taxon),] #create dataframe for testing
      SALdata_cval_testing <- SALdata_cval_testing[,!(names(SALdata_cval_testing) %in% toremove)]
      SALdata_cval_training <- SALdata_cval_training[,!(names(SALdata_cval_training) %in% toremove)]
      RF_model_testraits_cval <- randomForest(Group~., data = SALdata_cval_training, 
                                      ntree = 500)
      Prediction_testtraits <- predict(RF_model_testraits_cval, newdata = SALdata_cval_testing, type = "response")
      theVar <- rep(paste(theitems, collapse = "_"), length(Prediction_testtraits))
      MyPredictiontesttraits <- cbind(V1 = as.vector(Prediction_testtraits), V2 = as.vector(SALdata_cval_testing$Group),theVar)
      Crossvalresults_traits2 <- rbind(Crossvalresults_traits2,MyPredictiontesttraits)
    }
  }
}

MyVariableImportance <- data.frame()
for (item in all_combos){
  mycombos <- item
  for (i in 1:range(ncol(item))){
    theitems <- mycombos[,i]
    theitems <- paste(theitems, collapse = "_")
    relevant <- Crossvalresults_traits2[which(Crossvalresults_traits2$theVar == theitems),]
    truecryptic = relevant[relevant$V2 == 'C' & relevant$V1 == 'C',]
    truenoncryptic = relevant[relevant$V2 == 'NC' & relevant$V1 == 'NC',]
    totalcryptic = relevant[relevant$V1 == 'C',]
    totalnoncryptic = relevant[relevant$V1 == 'NC',]
    accuracyoverall = (nrow(truecryptic)+nrow(truenoncryptic))/(nrow(totalcryptic)+nrow(totalnoncryptic))
    accuracycryptic = (nrow(truecryptic)/nrow(totalcryptic))
    accuracynoncryptic = (nrow(truenoncryptic)/nrow(totalnoncryptic))
    VariableImportance <- cbind.data.frame(accuracyoverall,accuracycryptic,accuracynoncryptic,theitems)
    MyVariableImportance <- rbind.data.frame(MyVariableImportance, VariableImportance)
  }
}

MyVariableImportance$MDA <- 0.928 - MyVariableImportance$accuracyoverall
par(mai=c(1,2,1,1))
barplot(MyVariableImportance$MDA, horiz = TRUE, names.arg = MyVariableImportance$theitems, axes = TRUE, las=1)

write.table(MyVariableImportance, 'MyVariableImportance_SAL.csv', sep = ',')  


