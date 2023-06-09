# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())
graphics.off()

library(ggplot2)
library(cowplot)
library(dplyr)
library(phyloseq)
library(DESeq2)
library(caret)
library(randomForest)
library(pROC)

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
## set your working directory to where files are located
setwd(dir)

#phyloseq object (kmnp2) 
# raw ASV counts including singletons (will be removed)
# marked individuals only
load("Rdata/sifaka_allotus_phyloseq_openref_2016.RData")
head(sample_data(kmnp2))
#######################################################################
## Random Forest Analysis
#######################################################################
kmnp2
kmnp2 <- kmnp2 %>%
  subset_samples(Group %in% c("I","II","III","IV","V","VI")) 

##remove singletons
reduced_kmnp2 <- prune_taxa(taxa_sums(kmnp2)>1, kmnp2)
reduced_kmnp2
sample_data(reduced_kmnp2)
nsamples(reduced_kmnp2)#305
ntaxa(reduced_kmnp2)#3098
dim(otu_table(reduced_kmnp2))#312 3098

red <- subset_taxa(reduced_kmnp2, Phylum!="Unclassified")
red
# otu_table()   OTU Table:         [ 2732 taxa and 305 samples ]
# sample_data() Sample Data:       [ 305 samples by 24 sample variables ]
# tax_table()   Taxonomy Table:    [ 2732 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2732 tips and 2731 internal nodes ]

length(unique(sample_data(red)$Name))#57
## ASVs in 5% of samples
sites.prune <- filter_taxa(red, function(x) sum(x > 1) > (0.05*length(x)), prune=TRUE) 
sites.prune#405
sample_data(sites.prune)$Year  = as.factor(sample_data(sites.prune)$Year)
sample_data(sites.prune)$Group = as.factor(sample_data(sites.prune)$Group)
### normalize sequencing depth
kmnp.deseq <- phyloseq_to_deseq2(sites.prune, ~Group+Year) #convert phyloseq object to deseq2 object
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(kmnp.deseq), 1, gm_mean)
kmnp.deseq= estimateSizeFactors(kmnp.deseq, geoMeans = geoMeans,type="poscounts")
kmnp.deseq = estimateDispersions(kmnp.deseq)
diagvst = getVarianceStabilizedData(kmnp.deseq)
dim(diagvst)
#405 305
sifaka0 <- sites.prune #rename sites.prune to sifaka0 so that we still have the original OTU counts
normalizedCounts <- t( t(counts(kmnp.deseq)) / sizeFactors(kmnp.deseq) )
otu_table(sites.prune) <- otu_table(normalizedCounts, taxa_are_rows = TRUE)

sites.prune <- prune_taxa(taxa_sums(sites.prune)>0,sites.prune)
sites.prune
otu_table(sites.prune)<-round(otu_table(sites.prune))

predictors <- t(otu_table(sites.prune))
dim(predictors)#305 405

######################################
## Year
######################################
response <- as.factor(sample_data(sites.prune)$Year)
rf.data <- data.frame(response, predictors)
head(rf.data)

## stratify training and test sets by year and extraction method
set.seed(998)
df <- sample_data(sites.prune)%>%
  as_tibble()%>%
  mutate(n = row_number()) %>% #create row number if you dont have one
  select(n, everything()) # put 'n' at the front of the dataset
train <- df %>%
  group_by(Year, Extraction) %>% #any number of variables you wish to partition by proportionally
  sample_frac(.8) # '.8' is the proportion of the original df you wish to sample
test <- anti_join(df, train) # creates test dataframe with those observations not in 'train.'
test <- as.matrix(test$n)
train <- as.matrix(train$n)

nrow(rf.data)
nrow(train)
nrow(test)

training <- rf.data[train,]
testing <- rf.data[test,]
table(testing$response) #check that 2015 samples are included
nrow(training)#244
nrow(testing)#61
fit_control <- trainControl( method = "LOOCV",
                             allowParallel = TRUE,
                             savePredictions = TRUE,
                             verboseIter = T)    

##less computationally intensive but also less accurate
# fit_control <- trainControl(method = "repeatedcv",
#                            number = 10,
#                            repeats = 5,
#                            savePredictions = TRUE,
#                            verboseIter = TRUE)

rfFit <- train(response ~ ., data = training, method = "rf", preProc = "center",proximity = TRUE,ntree=1000,importance=TRUE, trControl=fit_control )
save(rfFit,file="Rdata/sifaka_RF_fit_year_classification.Rdata")
# load("Rdata/sifaka_RF_fit_year_classification.Rdata")
rfFit$bestTune #mtry 203
rfFit$finalModel
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$response)
cm = confusionMatrix(rfClasses, testing$response)
cm
postResample(pred = rfClasses, obs = testing$response)
# Accuracy     Kappa 
# 0.9672131 0.9426692

roc.multi <- multiclass.roc(predictor=as.numeric(rfClasses),response=as.numeric(testing$response))
auc(roc.multi)#0.9444
######################################
## Social Group
######################################
predictors <- t(otu_table(sites.prune))
dim(predictors)#305 405
    
response <- as.factor(sample_data(sites.prune)$Group)
rf.data <- data.frame(response, predictors)
head(rf.data)
nrow(rf.data)

## stratify training and test sets by social group and extraction method
set.seed(998)
df <- sample_data(sites.prune)%>%
  as_tibble()%>%
  mutate(n = row_number()) %>% #create row number if you dont have one
  select(n, everything()) # put 'n' at the front of the dataset
train <- df %>%
  group_by(Group, Extraction) %>% #any number of variables you wish to partition by proportionally
  sample_frac(.8) # '.8' is the proportion of the original df you wish to sample
test <- anti_join(df, train) # creates test dataframe with those observations not in 'train.'
test <- as.matrix(test$n)
train <- as.matrix(train$n)

training <- rf.data[train,]
testing <- rf.data[test,]
nrow(training)#243
nrow(testing)#62

fit_control <- trainControl( method = "LOOCV",allowParallel = TRUE,savePredictions = TRUE,verboseIter = T)    

rfFit <- train(response ~ ., data = training, method = "rf", preProc = "center", proximity = TRUE,ntree=1000,importance=TRUE, trControl=fit_control )
save(rfFit,file="Rdata/sifaka_RF_fit_group_classification.Rdata")
load("Rdata/sifaka_RF_fit_group_classification.Rdata")
rfFit$bestTune #mtry 203
rfFit$finalModel
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$response)

cm = confusionMatrix(rfClasses, testing$response)
cm
cm[["byClass"]][ , "F1"]
postResample(pred = rfClasses, obs = testing$response)
# Accuracy     Kappa 
#0.9354839 0.9219635 

roc.multi <- multiclass.roc(predictor=as.numeric(rfClasses),response=as.numeric(testing$response))
auc(roc.multi)#0.9692
