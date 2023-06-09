# set the working directory
setwd("/Users/msmith/Box Sync/MSmith/In Progress - Synced/PredictivePhylgeography_withTraitData/Results/PNW_August08")

# read in the data
PNWdata <- read.table("/Users/msmith/Box Sync/MSmith/In Progress - Synced/PredictivePhylgeography_withTraitData/Data/PNW_MLS_August06_2018_FullDataset_forRF.csv", sep=",", header=T)
head(PNWdata)
names(PNWdata)

# check the data
table(PNWdata$complex, PNWdata$taxon)
table(PNWdata$complex, PNWdata$dispStage)
table(PNWdata$complex, PNWdata$selfOut)
table(PNWdata$complex, PNWdata$dispersion)
table(PNWdata$complex, PNWdata$tropicLevel)
table(PNWdata$complex, PNWdata$maxSize)
class(PNWdata$taxon)
class(PNWdata$dispStage)
class(PNWdata$selfOut)
class(PNWdata$dispersion)
class(PNWdata$tropicLevel)
class(PNWdata$maxSize)

range(PNWdata$x)
range(PNWdata$y)

library(randomForest)

# no traits
AllTaxa_NoTraits_Taxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + taxon , data=PNWdata, ntree=500, importance = T)

oob_alltaxa_notraits_taxonomy <- AllTaxa_NoTraits_Taxonomy$err.rate[500,1]

# build the classifier All Taxa - All Traits - Taxonomy
AllTaxa_AllTraits_Taxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + taxon + dispStage + selfOut + dispersion 
                         + tropicLevel + maxSize, data=PNWdata, ntree=500, importance = T)

oob_alltaxa_alltraits_taxonomy <- AllTaxa_AllTraits_Taxonomy$err.rate[500,1]


# build the classifier All Taxa - All Traits - No Taxonomy
AllTaxa_AllTraits_NoTaxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage + selfOut + dispersion 
                         + tropicLevel + maxSize, data=PNWdata, ntree=500, importance = T)

oob_alltaxa_alltraits_notaxonomy <- AllTaxa_AllTraits_NoTaxonomy$err.rate[500,1]

# Create a dataframe that omits the two problematic taxa
ConfidentPNWdata <- PNWdata[PNWdata$complex!= 'Conaphe.armataEW' & PNWdata$complex!= 'Prophysaon.vanattae.humile',]
levels(ConfidentPNWdata$complex)
ConfidentPNWdata$complex <- droplevels(ConfidentPNWdata$complex, exclude = c("Conaphe.armataEW", "Prophysaon.vanattae.humile"))
table(ConfidentPNWdata$complex)

# build the classifier Curated Taxa -All Traits- Taxonomy
CuratedTaxa_AllTraits_Taxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + taxon + dispStage + selfOut + dispersion 
                         + tropicLevel + maxSize, data=ConfidentPNWdata, ntree=500, importance = T)

oob_curatedtaxa_alltraits_taxonomy <- CuratedTaxa_AllTraits_Taxonomy$err.rate[500,1]


# build the classifier Reduced Taxa - No Taxonomy
CuratedTaxa_AllTraits_NoTaxonomy <- randomForest(group ~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage + selfOut + dispersion 
                         + tropicLevel + maxSize, data=ConfidentPNWdata, ntree=500, importance = T)

oob_curatedtaxa_alltraits_notaxonomy <- CuratedTaxa_AllTraits_NoTaxonomy$err.rate[500,1]


# Plot MDA and mean decrease in GINI for each classifier
pdf('AllTaxa_AllTraits_Taxonomy.pdf')
varImpPlot(AllTaxa_AllTraits_Taxonomy)
dev.off()

pdf('AllTaxa_AllTraits_NoTaxonomy.pdf')
varImpPlot(AllTaxa_AllTraits_NoTaxonomy)
dev.off()

pdf('CuratedTaxa_AllTraits_Taxonomy.pdf')
varImpPlot(CuratedTaxa_AllTraits_Taxonomy)
dev.off()

pdf('CuratedTaxa_AllTraits_NoTaxonomy.pdf')
varImpPlot(CuratedTaxa_AllTraits_NoTaxonomy)
dev.off()

# look at results
AllTaxa_AllTraits_Taxonomy
AllTaxa_AllTraits_NoTaxonomy
CuratedTaxa_AllTraits_Taxonomy
CuratedTaxa_AllTraits_NoTaxonomy

# do cross-validation to access accuracy of the two classifiers with All Taxa included
Crossvalresults_1 <- data.frame() # create empty dataframe for All Taxa - Taxonomy
Crossvalresults_2 <- data.frame() # create empty dataframe for All Taxa - No Taxonomy
for (item in levels(PNWdata$complex)){ # loop through each complex
  PNWdata_cval_training <- PNWdata[which(PNWdata$complex != item),] #create dataframe for training
  PNWdata_cval_testing <- PNWdata[which(PNWdata$complex == item),] #create dataframe for testing
  RF_model_1_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage+ selfOut + dispersion 
                                  + tropicLevel + maxSize+taxon, data = PNWdata_cval_training, 
                                  ntree = 500, importance = T) # build the classifier All Taxa - Taxonomy
  RF_model_2_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage+ selfOut + dispersion 
                                  + tropicLevel + maxSize, data = PNWdata_cval_training, 
                                  ntree = 500, importance = T) # build the classifier All Taxa - No Taxonomy
  Prediction_1 <- predict(RF_model_1_cval, newdata = PNWdata_cval_testing, type = "response") # Predict on the training dataset for classifier 1
  MyPrediction1 <- cbind(Prediction_1, PNWdata_cval_testing) # bind the prediciton to the test data, so we know which taxon was omitted
  Prediction_2 <- predict(RF_model_2_cval, newdata = PNWdata_cval_testing, type = "response") # Predict on the training dataset for classifier 2
  MyPrediction2 <- cbind(Prediction_2, PNWdata_cval_testing) # bind the prediciton to the test data, so we know which taxon was omitted
  Crossvalresults_1 <- rbind(Crossvalresults_1,MyPrediction1) # Bind the previous crossval results and the results from the current taxon
  Crossvalresults_2 <- rbind(Crossvalresults_2, MyPrediction2) # Bind the previous crossval results and the results from the current taxon
}

AllTaxa_AllTraits_Taxonomy_truecryptic = Crossvalresults_1[Crossvalresults_1$group == 'C' & Crossvalresults_1$Prediction_1 == 'C',] # count true cryptic classifications
AllTaxa_AllTraits_Taxonomy_truenoncryptic = Crossvalresults_1[Crossvalresults_1$group == 'NC' & Crossvalresults_1$Prediction_1 == 'NC',] # count true noncryptic classifications
AllTaxa_AllTraits_Taxonomy_totalcryptic = Crossvalresults_1[Crossvalresults_1$Prediction_1 == 'C',] # count all cryptic classifications
AllTaxa_AllTraits_Taxonomy_totalnoncryptic = Crossvalresults_1[Crossvalresults_1$Prediction_1 == 'NC',] # count all non-cryptic classifications
AllTaxa_AllTraits_Taxonomy_accuracyoverall = (nrow(AllTaxa_AllTraits_Taxonomy_truecryptic)+nrow(AllTaxa_AllTraits_Taxonomy_truenoncryptic))/nrow(Crossvalresults_1) # calculate overall accuracy
AllTaxa_AllTraits_Taxonomy_accuracycryptic = (nrow(AllTaxa_AllTraits_Taxonomy_truecryptic)/nrow(AllTaxa_AllTraits_Taxonomy_totalcryptic)) # calculate cryptic accuracy
AllTaxa_AllTraits_Taxonomy_accuracynoncryptic = (nrow(AllTaxa_AllTraits_Taxonomy_truenoncryptic)/nrow(AllTaxa_AllTraits_Taxonomy_totalnoncryptic)) # calculate noncryptic accuracy
AllTaxa_AllTraits_Taxonomy_accuracyoverall
AllTaxa_AllTraits_Taxonomy_accuracycryptic
AllTaxa_AllTraits_Taxonomy_accuracynoncryptic

AllTaxa_AllTraits_NoTaxonomy_truecryptic = Crossvalresults_2[Crossvalresults_2$group == 'C' & Crossvalresults_2$Prediction_2 == 'C',]
AllTaxa_AllTraits_NoTaxonomy_truenoncryptic = Crossvalresults_2[Crossvalresults_2$group == 'NC' & Crossvalresults_2$Prediction_2 == 'NC',]
AllTaxa_AllTraits_NoTaxonomy_totalcryptic = Crossvalresults_2[Crossvalresults_2$Prediction_2 == 'C',]
AllTaxa_AllTraits_NoTaxonomy_totalnoncryptic = Crossvalresults_2[Crossvalresults_2$Prediction_2 == 'NC',]
AllTaxa_AllTraits_NoTaxonomy_accuracyoverall= (nrow(AllTaxa_AllTraits_NoTaxonomy_truecryptic)+nrow(AllTaxa_AllTraits_NoTaxonomy_truenoncryptic))/nrow(Crossvalresults_2)
AllTaxa_AllTraits_NoTaxonomy_accuracycryptic= (nrow(AllTaxa_AllTraits_NoTaxonomy_truecryptic)/nrow(AllTaxa_AllTraits_NoTaxonomy_totalcryptic))
AllTaxa_AllTraits_NoTaxonomy_accuracynoncryptic = (nrow(AllTaxa_AllTraits_NoTaxonomy_truenoncryptic)/nrow(AllTaxa_AllTraits_NoTaxonomy_totalnoncryptic))
AllTaxa_AllTraits_NoTaxonomy_accuracyoverall
AllTaxa_AllTraits_NoTaxonomy_accuracycryptic
AllTaxa_AllTraits_NoTaxonomy_accuracynoncryptic

Crossvalresults_3 <- data.frame()
Crossvalresults_4 <- data.frame()
for (item in levels(ConfidentPNWdata$complex)){
  ConfidentPNWdata_cval_training <- ConfidentPNWdata[which(ConfidentPNWdata$complex != item),] #create dataframe for training
  ConfidentPNWdata_cval_testing <- ConfidentPNWdata[which(ConfidentPNWdata$complex == item),] #create dataframe for testing
  RF_model_3_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage+ selfOut + dispersion 
                                  + tropicLevel + maxSize+taxon, data = ConfidentPNWdata_cval_training, 
                                  ntree = 500, importance = T)
  RF_model_4_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage+ selfOut + dispersion 
                                  + tropicLevel + maxSize, data = ConfidentPNWdata_cval_training, 
                                  ntree = 500, importance = T)
  Prediction_3 <- predict(RF_model_3_cval, newdata = ConfidentPNWdata_cval_testing, type = "response")
  MyPrediction3 <- cbind(Prediction_3, ConfidentPNWdata_cval_testing)
  Prediction_4 <- predict(RF_model_4_cval, newdata = ConfidentPNWdata_cval_testing, type = "response")
  MyPrediction4 <- cbind(Prediction_4, ConfidentPNWdata_cval_testing)
  Crossvalresults_3 <- rbind(Crossvalresults_3,MyPrediction3)
  Crossvalresults_4 <- rbind(Crossvalresults_4, MyPrediction4)

}
ReducedTaxa_AllTraits_Taxonomy_truecryptic = Crossvalresults_3[Crossvalresults_3$group == 'C' & Crossvalresults_3$Prediction_3 == 'C',]
ReducedTaxa_AllTraits_Taxonomy_truenoncryptic = Crossvalresults_3[Crossvalresults_3$group == 'NC' & Crossvalresults_3$Prediction_3 == 'NC',]
ReducedTaxa_AllTraits_Taxonomy_totalcryptic = Crossvalresults_3[Crossvalresults_3$Prediction_3 == 'C',]
ReducedTaxa_AllTraits_Taxonomy_totalnoncryptic = Crossvalresults_3[Crossvalresults_3$Prediction_3 == 'NC',]
ReducedTaxa_AllTraits_Taxonomy_accuracyoverall = (nrow(ReducedTaxa_AllTraits_Taxonomy_truecryptic)+nrow(ReducedTaxa_AllTraits_Taxonomy_truenoncryptic))/nrow(Crossvalresults_3)
ReducedTaxa_AllTraits_Taxonomy_accuracycryptic = (nrow(ReducedTaxa_AllTraits_Taxonomy_truecryptic)/nrow(ReducedTaxa_AllTraits_Taxonomy_totalcryptic))
ReducedTaxa_AllTraits_Taxonomy_accuracynoncryptic = (nrow(ReducedTaxa_AllTraits_Taxonomy_truenoncryptic)/nrow(ReducedTaxa_AllTraits_Taxonomy_totalnoncryptic))
ReducedTaxa_AllTraits_Taxonomy_accuracyoverall
ReducedTaxa_AllTraits_Taxonomy_accuracycryptic
ReducedTaxa_AllTraits_Taxonomy_accuracynoncryptic

ReducedTaxa_AllTraits_NoTaxonomy_truecryptic = Crossvalresults_4[Crossvalresults_4$group == 'C' & Crossvalresults_4$Prediction_4 == 'C',]
ReducedTaxa_AllTraits_NoTaxonomy_truenoncryptic = Crossvalresults_4[Crossvalresults_4$group == 'NC' & Crossvalresults_4$Prediction_4 == 'NC',]
ReducedTaxa_AllTraits_NoTaxonomy_totalcryptic = Crossvalresults_4[Crossvalresults_4$Prediction_4 == 'C',]
ReducedTaxa_AllTraits_NoTaxonomy_totalnoncryptic = Crossvalresults_4[Crossvalresults_4$Prediction_4 == 'NC',]
ReducedTaxa_AllTraits_NoTaxonomy_accuracyoverall = (nrow(ReducedTaxa_AllTraits_NoTaxonomy_truecryptic)+nrow(ReducedTaxa_AllTraits_NoTaxonomy_truenoncryptic))/nrow(Crossvalresults_4)
ReducedTaxa_AllTraits_NoTaxonomy_accuracycryptic = (nrow(ReducedTaxa_AllTraits_NoTaxonomy_truecryptic)/nrow(ReducedTaxa_AllTraits_NoTaxonomy_totalcryptic))
ReducedTaxa_AllTraits_NoTaxonomy_accuracynoncryptic = (nrow(ReducedTaxa_AllTraits_NoTaxonomy_truenoncryptic)/nrow(ReducedTaxa_AllTraits_NoTaxonomy_totalnoncryptic))
ReducedTaxa_AllTraits_NoTaxonomy_accuracyoverall
ReducedTaxa_AllTraits_NoTaxonomy_accuracycryptic
ReducedTaxa_AllTraits_NoTaxonomy_accuracynoncryptic


Crossvalresults_5 <- data.frame()
names(ConfidentPNWdata)
PNWdata_forRF <- ConfidentPNWdata[,c(2:20, 27:31, 25, 26)]
names(PNWdata_forRF)
traitlist <- c("dispStage", "selfOut", "dispersion" 
                          , "tropicLevel", "maxSize")


combos_1<- combn(traitlist, m=1,simplyify=F)
combos_2<- combn(traitlist, m=2,simplyify=F)
combos_3<- combn(traitlist, m=3,simplyify=F)
combos_4<- combn(traitlist, m=4,simplyify=F)
combos_5<- combn(traitlist, m=5,simplyify=F)

all_combos <- list(combos_1)


for (taxon in levels(PNWdata_forRF$complex)){
  for (item in all_combos){
    mycombos <- item
    for (i in 1:range(ncol(item))){
      theitems <- mycombos[,i]
      PNW_corvars <- PNWdata_forRF[,!(names(PNWdata_forRF) %in% theitems)]
      toremove <- list('complex')
      PNWdata_cval_training <- PNW_corvars[which(PNW_corvars$complex != taxon),] #create dataframe for training
      PNWdata_cval_testing <- PNW_corvars[which(PNW_corvars$complex == taxon),] #create dataframe for testing
      PNWdata_cval_testing <- PNWdata_cval_testing[,!(names(PNWdata_cval_testing) %in% toremove)]
      PNWdata_cval_training <- PNWdata_cval_training[,!(names(PNWdata_cval_training) %in% toremove)]
      RF_model_5_cval <- randomForest(group~., data = PNWdata_cval_training, 
                                          ntree = 500)
      Prediction_5 <- predict(RF_model_5_cval, newdata = PNWdata_cval_testing, type = "response")
      theVar <- rep(paste(theitems, collapse = "_"), length(Prediction_5))
      MyPrediction5 <- cbind(V1 = as.vector(Prediction_5), V2 = as.vector(PNWdata_cval_testing$group),theVar)
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


MyVariableImportance$MDA <- ReducedTaxa_AllTraits_NoTaxonomy_accuracyoverall - MyVariableImportance$accuracyoverall
pdf('VariableImportance_PNW.pdf')
par(mai=c(1,2,1,1))
barplot(MyVariableImportance$MDA, horiz = TRUE, names.arg = MyVariableImportance$theitems,  las=1)
dev.off()

write.table(MyVariableImportance, 'MyVariableImportance_PNW.csv', sep = ',')  

# build a new model based on these results
CuratedTaxa_CuratedTraits_NoTaxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                           + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                           + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                           + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                           + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                           + dispStage + selfOut + dispersion
                           + tropicLevel  , data=ConfidentPNWdata, ntree=500, importance = T)

CuratedTaxa_CuratedTraits_NoTaxonomy
oob_curatedtaxa_curatedtraits_notaxonomy <- CuratedTaxa_CuratedTraits_NoTaxonomy$err.rate[500,1]
pdf('CuratedTaxa_CuratedTraits_NoTaxonomy.pdf')
varImpPlot(CuratedTaxa_CuratedTraits_NoTaxonomy)
dev.off()

# build another model with taxonomy
CuratedTaxa_CuratedTraits_Taxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                               + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                               + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                               + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                               + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                               + dispStage + selfOut + dispersion
                               + tropicLevel +taxon , data=ConfidentPNWdata, ntree=500, importance = T)
CuratedTaxa_CuratedTraits_Taxonomy
oob_curatedtaxa_curatedltraits_taxonomy <- CuratedTaxa_CuratedTraits_Taxonomy$err.rate[500,1]
pdf('CuratedTaxa_CuratedTraits_Taxonomy.pdf')
varImpPlot(CuratedTaxa_CuratedTraits_Taxonomy)
dev.off()

# build another model with taxonomy only
CuratedTaxa_NoTraits_Taxonomy <- randomForest(group ~ wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                                                   + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                                                   + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                                                   + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                                                   + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                                                   +taxon , data=ConfidentPNWdata, ntree=500, importance = T)
CuratedTaxa_NoTraits_Taxonomy
oob_curatedtaxa_notraits_taxonomy <- CuratedTaxa_NoTraits_Taxonomy$err.rate[500,1]
pdf('CuratedTaxa_NoTraits_Taxonomy.pdf')
varImpPlot(CuratedTaxa_NoTraits_Taxonomy)
dev.off()

# cross validation for curated traits & no taxonomy
Crossvalresults_7 <- data.frame()
for (item in levels(ConfidentPNWdata$complex)){
  PNWdata_cval_training <- ConfidentPNWdata[which(ConfidentPNWdata$complex != item),] #create dataframe for training
  PNWdata_cval_testing <- ConfidentPNWdata[which(ConfidentPNWdata$complex == item),] #create dataframe for testing
  RF_model_7_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                                  + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                                  + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                                  + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                                  + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                                  + dispStage+ selfOut + dispersion
                                  + tropicLevel, data = PNWdata_cval_training, 
                                  ntree = 500, importance = T)
  Prediction_7 <- predict(RF_model_7_cval, newdata = PNWdata_cval_testing, type = "response")
  Prediction_7b <- predict(RF_model_7_cval, newdata = PNWdata_cval_testing, type = "prob")
  MyPrediction7 <- cbind(Prediction_7, PNWdata_cval_testing, Prediction_7b)
  Crossvalresults_7 <- rbind(Crossvalresults_7,MyPrediction7)
}
CuratedTaxa_CuratedTraits_NoTaxonomy_truecryptic = Crossvalresults_7[Crossvalresults_7$group == 'C' & Crossvalresults_7$Prediction_7 == 'C',]
CuratedTaxa_CuratedTraits_NoTaxonomy_truenoncryptic = Crossvalresults_7[Crossvalresults_7$group == 'NC' & Crossvalresults_7$Prediction_7 == 'NC',]
CuratedTaxa_CuratedTraits_NoTaxonomy_totalcryptic = Crossvalresults_7[Crossvalresults_7$Prediction_7 == 'C',]
CuratedTaxa_CuratedTraits_NoTaxonomy_totalnoncryptic = Crossvalresults_7[Crossvalresults_7$Prediction_7 == 'NC',]
CuratedTaxa_CuratedTraits_NoTaxonomy_accuracyoverall = (nrow(CuratedTaxa_CuratedTraits_NoTaxonomy_truecryptic)+nrow(CuratedTaxa_CuratedTraits_NoTaxonomy_truenoncryptic))/nrow(Crossvalresults_7)
CuratedTaxa_CuratedTraits_NoTaxonomy_accuracycryptic = (nrow(CuratedTaxa_CuratedTraits_NoTaxonomy_truecryptic)/nrow(CuratedTaxa_CuratedTraits_NoTaxonomy_totalcryptic))
CuratedTaxa_CuratedTraits_NoTaxonomy_accuracynoncryptic = (nrow(CuratedTaxa_CuratedTraits_NoTaxonomy_truenoncryptic)/nrow(CuratedTaxa_CuratedTraits_NoTaxonomy_totalnoncryptic))
CuratedTaxa_CuratedTraits_NoTaxonomy_accuracyoverall
CuratedTaxa_CuratedTraits_NoTaxonomy_accuracycryptic
CuratedTaxa_CuratedTraits_NoTaxonomy_accuracynoncryptic


# crossval for taxonomy only
Crossvalresults_8 <- data.frame()
for (item in levels(ConfidentPNWdata$complex)){
  PNWdata_cval_training <- ConfidentPNWdata[which(ConfidentPNWdata$complex != item),] #create dataframe for training
  PNWdata_cval_testing <- ConfidentPNWdata[which(ConfidentPNWdata$complex == item),] #create dataframe for testing
  RF_model_8_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                                  + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                                  + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                                  + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                                  + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                                  + taxon, data = PNWdata_cval_training, 
                                  ntree = 500, importance = T)
  Prediction_8 <- predict(RF_model_8_cval, newdata = PNWdata_cval_testing, type = "response")
  Prediction_8b <- predict(RF_model_8_cval, newdata = PNWdata_cval_testing, type = "prob")
  MyPrediction8 <- cbind(Prediction_8, PNWdata_cval_testing, Prediction_8b)
  Crossvalresults_8 <- rbind(Crossvalresults_8,MyPrediction8)
}

CuratedTaxa_NoTraits_Taxonomy_truecryptic = Crossvalresults_8[Crossvalresults_8$group == 'C' & Crossvalresults_8$Prediction_8 == 'C',]
CuratedTaxa_NoTraits_Taxonomy_truenoncryptic = Crossvalresults_8[Crossvalresults_8$group == 'NC' & Crossvalresults_8$Prediction_8 == 'NC',]
CuratedTaxa_NoTraits_Taxonomy_totalcryptic = Crossvalresults_8[Crossvalresults_8$Prediction_8 == 'C',]
CuratedTaxa_NoTraits_Taxonomy_totalnoncryptic = Crossvalresults_8[Crossvalresults_8$Prediction_8 == 'NC',]
CuratedTaxa_NoTraits_Taxonomy_accuracyoverall = (nrow(CuratedTaxa_NoTraits_Taxonomy_truecryptic)+nrow(CuratedTaxa_NoTraits_Taxonomy_truenoncryptic))/nrow(Crossvalresults_8)
CuratedTaxa_NoTraits_Taxonomy_accuracycryptic = (nrow(CuratedTaxa_NoTraits_Taxonomy_truecryptic)/nrow(CuratedTaxa_NoTraits_Taxonomy_totalcryptic))
CuratedTaxa_NoTraits_Taxonomy_accuracynoncryptic = (nrow(CuratedTaxa_NoTraits_Taxonomy_truenoncryptic)/nrow(CuratedTaxa_NoTraits_Taxonomy_totalnoncryptic))
CuratedTaxa_NoTraits_Taxonomy_accuracyoverall
CuratedTaxa_NoTraits_Taxonomy_accuracycryptic
CuratedTaxa_NoTraits_Taxonomy_accuracynoncryptic

# cross val with curated traits & taxonomy
Crossvalresults_9 <- data.frame()
for (item in levels(ConfidentPNWdata$complex)){
  PNWdata_cval_training <- ConfidentPNWdata[which(ConfidentPNWdata$complex != item),] #create dataframe for training
  PNWdata_cval_testing <- ConfidentPNWdata[which(ConfidentPNWdata$complex == item),] #create dataframe for testing
  RF_model_9_cval <- randomForest(group~wc2.0_bio_30s_01 + wc2.0_bio_30s_02 + wc2.0_bio_30s_03
                                  + wc2.0_bio_30s_04 + wc2.0_bio_30s_05 + wc2.0_bio_30s_06 + wc2.0_bio_30s_07 
                                  + wc2.0_bio_30s_08 + wc2.0_bio_30s_09 + wc2.0_bio_30s_10 + wc2.0_bio_30s_11 
                                  + wc2.0_bio_30s_12 + wc2.0_bio_30s_13 + wc2.0_bio_30s_14 + wc2.0_bio_30s_15 
                                  + wc2.0_bio_30s_16 + wc2.0_bio_30s_17 + wc2.0_bio_30s_18 + wc2.0_bio_30s_19 
                                  + dispStage+ selfOut + dispersion
                                  + tropicLevel+ taxon, data = PNWdata_cval_training, 
                                  ntree = 500, importance = T)
  Prediction_9 <- predict(RF_model_9_cval, newdata = PNWdata_cval_testing, type = "response")
  Prediction_9b <- predict(RF_model_9_cval, newdata = PNWdata_cval_testing, type = "prob")
  MyPrediction9 <- cbind(Prediction_9, PNWdata_cval_testing, Prediction_9b)
  Crossvalresults_9 <- rbind(Crossvalresults_9,MyPrediction9)
}

CuratedTraits_CuratedTaxa_Taxonomy_truecryptic = Crossvalresults_9[Crossvalresults_9$group == 'C' & Crossvalresults_9$Prediction_9 == 'C',]
CuratedTraits_CuratedTaxa_Taxonomy_truenoncryptic = Crossvalresults_9[Crossvalresults_9$group == 'NC' & Crossvalresults_9$Prediction_9 == 'NC',]
CuratedTraits_CuratedTaxa_Taxonomy_totalcryptic = Crossvalresults_9[Crossvalresults_9$Prediction_9 == 'C',]
CuratedTraits_CuratedTaxa_Taxonomy_totalnoncryptic = Crossvalresults_9[Crossvalresults_9$Prediction_9 == 'NC',]
CuratedTraits_CuratedTaxa_Taxonomy_accuracyoverall = (nrow(CuratedTraits_CuratedTaxa_Taxonomy_truecryptic)+nrow(CuratedTraits_CuratedTaxa_Taxonomy_truenoncryptic))/nrow(Crossvalresults_9)
CuratedTraits_CuratedTaxa_Taxonomy_accuracycryptic = (nrow(CuratedTraits_CuratedTaxa_Taxonomy_truecryptic)/nrow(CuratedTraits_CuratedTaxa_Taxonomy_totalcryptic))
CuratedTraits_CuratedTaxa_Taxonomy_accuracynoncryptic = (nrow(CuratedTraits_CuratedTaxa_Taxonomy_truenoncryptic)/nrow(CuratedTraits_CuratedTaxa_Taxonomy_totalnoncryptic))
CuratedTraits_CuratedTaxa_Taxonomy_accuracyoverall
CuratedTraits_CuratedTaxa_Taxonomy_accuracycryptic
CuratedTraits_CuratedTaxa_Taxonomy_accuracynoncryptic

library(ggplot2)

Alnus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Alnus.rubra",c(33:34)]
Alnus_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Alnus.rubra",c(33:34)]
Alnus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Alnus.rubra",c(33:34)]
Alnus_cv7$RF <- 'traits'
Alnus_cv8$RF <- 'taxonomy'
Alnus_cv9$RF <- 'both'
Alnus_cv <- rbind(Alnus_cv7,Alnus_cv8,Alnus_cv9)

Ascaphus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Ascaphus.montanus.truei",c(33:34)]
Ascaphus_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Ascaphus.montanus.truei",c(33:34)]
Ascaphus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Ascaphus.montanus.truei",c(33:34)]
Ascaphus_cv7$RF <-  'traits'
Ascaphus_cv8$RF <-  'taxonomy'
Ascaphus_cv9$RF <-  'both'
Ascaphus_cv <- rbind(Ascaphus_cv7,Ascaphus_cv8,Ascaphus_cv9)


Dicamptadon_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Dicamptodon.ater.tene",c(33:34)]
Dicamptadon_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Dicamptodon.ater.tene",c(33:34)]
Dicamptadon_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Dicamptodon.ater.tene",c(33:34)]
Dicamptadon_cv7$RF <-  'traits'
Dicamptadon_cv8$RF <-  'taxonomy'
Dicamptadon_cv9$RF <-  'both'
Dicamptadon_cv <- rbind(Dicamptadon_cv7,Dicamptadon_cv8,Dicamptadon_cv9)

Haplotrema_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Haplotrema.vancouverense",c(33:34)]
Haplotrema_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Haplotrema.vancouverense",c(33:34)]
Haplotrema_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Haplotrema.vancouverense",c(33:34)]
Haplotrema_cv7$RF <-  'traits'
Haplotrema_cv8$RF <-  'taxonomy'
Haplotrema_cv9$RF <-  'both'
Haplotrema_cv <- rbind(Haplotrema_cv7,Haplotrema_cv8,Haplotrema_cv9)

Microtus_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Microtus.richardsoni",c(33:34)]
Microtus_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Microtus.richardsoni",c(33:34)]
Microtus_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Microtus.richardsoni",c(33:34)]
Microtus_cv7$RF <-  'traits'
Microtus_cv8$RF <-  'taxonomy'
Microtus_cv9$RF <-  'both'
Microtus_cv <- rbind(Microtus_cv7,Microtus_cv8,Microtus_cv9)

Plethodon_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Plethodon.idaho.vandy",c(33:34)]
Plethodon_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Plethodon.idaho.vandy",c(33:34)]
Plethodon_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Plethodon.idaho.vandy",c(33:34)]
Plethodon_cv7$RF <-  'traits'
Plethodon_cv8$RF <-  'taxonomy'
Plethodon_cv9$RF <-  'both'
Plethodon_cv <- rbind(Plethodon_cv7,Plethodon_cv9,Plethodon_cv8)

Andersoni_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Prophysaon.andersoni",c(33:34)]
Andersoni_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Prophysaon.andersoni",c(33:34)]
Andersoni_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Prophysaon.andersoni",c(33:34)]
Andersoni_cv7$RF <-  'traits'
Andersoni_cv8$RF <-  'taxonomy'
Andersoni_cv9$RF <-  'both'
Andersoni_cv <- rbind(Andersoni_cv7,Andersoni_cv9,Andersoni_cv8)

Coeruleum_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Prophysaon.coeruEW",c(33:34)]
Coeruleum_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Prophysaon.coeruEW",c(33:34)]
Coeruleum_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Prophysaon.coeruEW",c(33:34)]
Coeruleum_cv7$RF <-  'traits'
Coeruleum_cv8$RF <-  'taxonomy'
Coeruleum_cv9$RF <-  'both'
Coeruleum_cv <- rbind(Coeruleum_cv7,Coeruleum_cv9,Coeruleum_cv8)

Dubium_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Prophysaon.dubium",c(33:34)]
Dubium_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Prophysaon.dubium",c(33:34)]
Dubium_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Prophysaon.dubium",c(33:34)]
Dubium_cv7$RF <-  'traits'
Dubium_cv8$RF <-  'taxonomy'
Dubium_cv9$RF <-  'both'
Dubium_cv <- rbind(Dubium_cv7,Dubium_cv8,Dubium_cv9)

Salix_cv7 <- Crossvalresults_7[Crossvalresults_7$complex == "Salix.melanopsis",c(33:34)]
Salix_cv8 <- Crossvalresults_8[Crossvalresults_8$complex == "Salix.melanopsis",c(33:34)]
Salix_cv9 <- Crossvalresults_9[Crossvalresults_9$complex == "Salix.melanopsis",c(33:34)]
Salix_cv7$RF <-  'traits'
Salix_cv8$RF <-  'taxonomy'
Salix_cv9$RF <-  'both'
Salix_cv <- rbind(Salix_cv7,Salix_cv8,Salix_cv9)

library(ggplot2)

Alnusplot <-        ggplot(Alnus_cv,        aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Alnus') +                        theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Ascaphusplot <-     ggplot(Ascaphus_cv,     aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Ascaphus') +                     theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Dicamptadonplot <-  ggplot(Dicamptadon_cv,  aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Dicamptadon') +                  theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Haplotremaplot <-   ggplot(Haplotrema_cv,   aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Haplotrema') +                   theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Microtusplot <-     ggplot(Microtus_cv,     aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Microtus') +                     theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Plethodonplot <-    ggplot(Plethodon_cv,    aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Plethodon') +                    theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Andersoniplot <-    ggplot(Andersoni_cv,    aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Prophysaon andersoni') +         theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Coeruleumplot <-    ggplot(Coeruleum_cv,    aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Prophysaon coeruleum') +         theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Dubiumplot <-       ggplot(Dubium_cv,       aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Prophysaon dubium') +            theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
Salixplot <-        ggplot(Salix_cv,        aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + guides(fill=FALSE) + labs(title = 'Salix') +                        theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))

library(gridExtra)
pdf("PNW_comparisons.pdf")
grid.arrange(Alnusplot, Ascaphusplot, Dicamptadonplot, Haplotremaplot, Microtusplot, Plethodonplot, Andersoniplot, Coeruleumplot, Dubiumplot, Salixplot, nrow = 4)
dev.off()

Alnusplotlegend <-        ggplot(Alnus_cv,        aes(C, fill = RF)) + geom_density(alpha = 0.2, aes(y=..scaled..)) + xlim(0, 1) + labs(title = 'Alnus') + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual( values = c("red","blue", "yellow"))
pdf("Alnusforlegend.pdf")
grid.arrange(Alnusplotlegend)
dev.off()


table(PNWdata$complex, PNWdata$group)

# write oob error rates to table
oobrates <- cbind(oob_alltaxa_notraits_taxonomy,
      oob_alltaxa_alltraits_taxonomy,
      oob_alltaxa_alltraits_notaxonomy,
      oob_curatedtaxa_alltraits_taxonomy,
      oob_curatedtaxa_alltraits_notaxonomy,
      oob_curatedtaxa_curatedtraits_notaxonomy,
      oob_curatedtaxa_curatedltraits_taxonomy,
      oob_curatedtaxa_notraits_taxonomy)

write.csv(file = "OOB_Rates.csv", oobrates)

predictionaccuracies <- cbind(AllTaxa_AllTraits_Taxonomy_accuracyoverall,
                              AllTaxa_AllTraits_Taxonomy_accuracycryptic,
                              AllTaxa_AllTraits_Taxonomy_accuracynoncryptic,
                              AllTaxa_AllTraits_NoTaxonomy_accuracyoverall,
                              AllTaxa_AllTraits_NoTaxonomy_accuracycryptic,
                              AllTaxa_AllTraits_NoTaxonomy_accuracynoncryptic,
                              ReducedTaxa_AllTraits_Taxonomy_accuracyoverall,
                              ReducedTaxa_AllTraits_Taxonomy_accuracycryptic,
                              ReducedTaxa_AllTraits_Taxonomy_accuracynoncryptic,
                              ReducedTaxa_AllTraits_NoTaxonomy_accuracyoverall,
                              ReducedTaxa_AllTraits_NoTaxonomy_accuracycryptic,
                              ReducedTaxa_AllTraits_NoTaxonomy_accuracynoncryptic,
                              CuratedTaxa_CuratedTraits_NoTaxonomy_accuracyoverall,
                              CuratedTaxa_CuratedTraits_NoTaxonomy_accuracycryptic,
                              CuratedTaxa_CuratedTraits_NoTaxonomy_accuracynoncryptic,
                              CuratedTaxa_NoTraits_Taxonomy_accuracyoverall,
                              CuratedTaxa_NoTraits_Taxonomy_accuracycryptic,
                              CuratedTaxa_NoTraits_Taxonomy_accuracynoncryptic,
                              CuratedTraits_CuratedTaxa_Taxonomy_accuracyoverall,
                              CuratedTraits_CuratedTaxa_Taxonomy_accuracycryptic,
                              CuratedTraits_CuratedTaxa_Taxonomy_accuracynoncryptic
)

write.csv(file = "Prediction_Accuracies.csv", predictionaccuracies)
