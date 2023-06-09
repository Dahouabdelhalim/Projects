## Date: 12-5-2019
## Personal: Jagdeep Singh Sidhu (jagdeeproots@gmail.com)
## Purpose: This code:
##          1. is used to conduct discriminant analysis for different root phenotypes taking each location*depth as separate variable
##
## Data needed: Data file named "alllocs_allphenotypes_data.csv" is provided along with 
##              An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" manuscript.

## Figures/tables in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)"
##             from this code: Figure 9.


######### Needed Packages ###########
#install.packages("finalfit")
library(finalfit)
library(dplyr)
library(plotrix)
library("tripack")
library(tidyverse)
library(ggplot2)
library(reshape2)
library(Rmisc) ## for summarySE
library(ggrepel) ## for avoiding overlapping labels
library(magrittr) ## for changing class of mulitple columns
library(ggthemes)
library(ggpubr) # for multiple comps on ggplot
library(rlang)
library(caret)
library(MASS) # QDA qdrtc discrimant

############# SET directory function ##############
  directory <- paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Comparingdifferenephenotypes/D95/Discrminant analysis")

############# read file ###############
########## By using RLD ##############
RVW_6locs <- read.csv("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Comparingdifferenephenotypes/D95/alllocs_allphenotypes_data.csv") %>%
mutate (Location = case_when( Location == "1" ~ "Loc1",
                              Location == "2" ~ "Loc2",
                              Location == "3" ~ "Loc3", 
                              Location == "4" ~ "Loc4",
                              Location == "5" ~ "Loc5",
                              Location == "6" ~ "Loc6",
                              # Location == "rustic" ~ "rustic",
                              # Location == "voronoi" ~ "voronoi",
                              # Location == "Wholeplot" ~"Wholeplot"
                              ))


RVW_6locs$Location <- as.factor(RVW_6locs$Location)
summary(RVW_6locs)

RVW_6locs_dcast <-  dcast(RVW_6locs, Rep + phenotype + crop ~ Location + Depth, value.var = "Length.cm")
summary(RVW_6locs_dcast)
write.csv(RVW_6locs_dcast, "RVW_6locs_dcast.csv")

########################################################################
############### Creating training and testing datasets #################

crops <- c("maize", "bean")
reps <- c(11,16,21,26,31,36,41,46,51)
random_rep <- c(1,2,3,4,5,6,7,8,9,10)
phenotypes <- c("maize_", "maize_shallow_", "maize_deep_", "bean_", "bean_shallow_", "bean_deep_")
train.data.1 <-  data.frame()

# Split the data into training (80%) and test set (20%)
for (crp in crops){
  
  ##### Uncomment this section for running QDA on 3 phenotypes
  # if (crp == "maize"){
  # data <- subset(RVW_6locs_dcast, crop == crp )}else{
  #   data <- subset(RVW_6locs_dcast, crop == crp )
  # }
  
  #### Run this section for running QDA on 2 phenotypes (Deep vs intermediate)
  if (crp == "maize"){
  data <- subset(RVW_6locs_dcast, crop == crp & phenotype != "maize_")}else{
   data <- subset(RVW_6locs_dcast, crop == crp & phenotype != "bean_" )
  }



for (rep in reps) {
  for (rand_rep in random_rep){
    for (pheno in phenotypes){
  # Randomly pick number of rows #
  data_training <- subset(data, Rep < 51 & phenotype == pheno, drop = TRUE)

  randomRows <- sample(1:length(data_training[,1]), rep, replace=T)
  train.data <- data_training %>% slice(randomRows)
  
  train.data.1 <- rbind (train.data.1, train.data)
  }
  train.data.1$random_rep <- rand_rep  
  #train.data <- subset(data, Rep < rep, drop = TRUE)
  test.data <- subset(data, Rep >= rep, drop = TRUE)

  write.csv(train.data.1, paste(crp,"_",rep,"_",rand_rep,"_train_data",".csv", sep = ""))
  write.csv(test.data, paste(crp,"_",rep,"_test_data",".csv", sep = ""))
  train.data.1 = data.frame()
  }
}
}

#############################################
################ QDA ########################

##### data.frame to store error rates #######
Error.rates.train <- data.frame()
Error.rates.test <- data.frame()
confusion_matrix <- data.frame()

crops <- c("maize", "bean")
reps <- c(51)
Locations <- c("Loc1", "Loc2", "Loc3", "Loc6")
Best_locs <- c("Loc3", "Loc6")

for (crp in crops){
  for(rep in reps){
   # for(rand_rep in random_rep){
            for(loc in Locations){
              for (selected_loc in Best_locs){
              train.data1 <- read.csv(paste(crp,"_",rep,"_",rand_rep,"_train_data",".csv", sep = ""))
              test.data <- read.csv(paste(crp,"_",51,"_test_data",".csv", sep = ""))

########## QDA model  ################

train.data <- dplyr::select(train.data1, "phenotype" , c( starts_with(paste(selected_loc)), starts_with(paste(loc))))
#train.data <- dplyr::select(train.data1, "phenotype" , starts_with(paste(loc)))
              
              QDA_model <- qda(phenotype ~ ., data = train.data)

              QDA_model
              predicted.classes <- QDA_model %>% predict(test.data)
             

train.qda.pred <- predict(QDA_model)$class ### errors in training set
test.qda.pred <- predict(QDA_model, newdata = test.data) ###### prediction error


table_test <- table(test.data$phenotype, test.qda.pred$class, dnn = c('Actual Group','Predicted Group'))
table_train <- table(train.data$phenotype, train.qda.pred, dnn = c('Actual Group','Predicted Group'))



print(c("Crop ",crp))
print(c("Location",loc, selected_loc))
print(c("Training set reps used",rep))
print("Accuracy")
print(mean(test.qda.pred$class == test.data$phenotype))
print(table_test)

#print(table_train)

###################################################
########### calculate error rates #################
## training set 
error_train = (train.data %>% 
  mutate(qda_pred = (train.qda.pred)) %>%
  summarise(qda.error = mean(phenotype != qda_pred)))
 #print(error_train)

## test set  
 error_test = (test.data %>% 
                  mutate(qda_pred = (test.qda.pred$class)) %>%
                  summarise(qda.error = mean(phenotype != qda_pred)))
 
 #print(error_test)
 #print(error_train)

 
 #### store values in dataframe 
 df <- data.frame(crp,rep,random_rep,selected_loc,error_train$qda.error, loc)
Error.rates.train <- rbind(df, Error.rates.train)

df <- data.frame(crp,rep,random_rep,selected_loc, error_test$qda.error, loc)
Error.rates.test <- rbind(df, Error.rates.test)
}
}
}
}

## Depending on your run (3 or 2 phenotypes comparision), choose from following to save files: 
#write.csv(Error.rates.train, "Error_rates_train_combined_deepvsshallow.csv")
#write.csv(Error.rates.test, "Error_rates_test_combined_deepvsshallow.csv")
#write.csv(Error.rates.train, "Error_rates_train_combined_deepvsshallowvsmid.csv")
#write.csv(Error.rates.test, "Error_rates_test_combined_deepvsshallowvsmid.csv")


################################
########### Main figure 9 ######

Error.rates.test$loc_comb <- paste(Error.rates.test$loc,"_",Error.rates.test$selected_loc)
summarized_data <- summarySE(Error.rates.test, measurevar="error_test.qda.error", groupvars=c("loc_comb", "crp", "rep"))
summary(summarized_data)

summarized_data$loc_comb = as.factor(summarized_data$loc_comb)
#write.csv(summarized_data,"summarized_error_rates_three_phenotypes.csv")
#write.csv(summarized_data,"summarized_error_rates_two_phenotypes.csv")

summarized_error_test_set <- ggplot(summarized_data, aes(rep, error_test.qda.error)) + geom_point(aes(colour = loc_comb), size = 1) + 
  geom_line(aes(colour = loc_comb), size = 1) +
    ggtitle("Misclassification rates for test set of 50 reps") + 
  geom_errorbar(data = summarized_data, aes(ymin=error_test.qda.error-se, ymax = error_test.qda.error+se),  width=.1) +
  facet_wrap(~crp, scales = "free")+
   xlab("No. of reps used in training set") + ylab ("Misclassification Error rate (0-1)") +
  # theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 15))
  theme(axis.text = element_text(colour = "black", size = 20, hjust = 1, vjust=0.2), axis.title = element_text(colour = "black", size = 20), axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.key = element_blank(), legend.position = c(0.9, 0.8)) + #+scale_y_continuous(breaks=0:60*4)
  scale_y_continuous(limits = c(0,0.7),expand = c(0,0))+
  scale_x_continuous(limits = c(0,60),expand = c(0,0))

summarized_error_test_set 

### save your plot accordingly (comparing two phenotypes or three)
#ggsave("2.error_test_set_Wholeplot_witherror_rate.eps", summarized_error_test_set)
#ggsave("3.error_test_set_Wholeplot_witherror_rate.eps", summarized_error_test_set)

#############################################################################################
## After running this, for each crop you should have:
## 1. Training data sets and test data sets for three phenotypes QDA and two phenotype QDA. 
## 2. Error rate files for three phenotypes QDA and two phenotype QDA. Data used in Table 1.
## 3. Plots for each QDA (two and three phenotypes). Which are then shown in Figure 9








############################################
###########################################
############## Supplementary ##############
error_training_set <- ggplot(Error.rates.train, aes(rep, error_train.qda.error)) + geom_point(aes(colour = loc), size = 2) + ggtitle("Misclassification rates for training set of 50 reps") + 
  facet_wrap(~crp,scales = "free")+
  xlab("No. of reps used") + ylab ("Misclassification Error rate (0-1)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 15))

error_test_set <- ggplot(Error.rates.test, aes(rep, error_test.qda.error)) + geom_point(aes(colour = loc), size = 2) + ggtitle("Misclassification rates for test set of 50 reps") + 
  facet_wrap(~crp, scales = "free")+
  xlab("No. of reps used in training set") + ylab ("Misclassification Error rate (0-1)") +
  # theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 15))
  theme(axis.text = element_text(colour = "black", size = 14, hjust = 1, vjust=0.2), axis.title = element_text(colour = "black", size = 20), axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.key = element_blank(), legend.position = c(0.9, 0.8)) + #+scale_y_continuous(breaks=0:60*4)
  scale_y_continuous(limits = c(0,0.7),expand = c(0,0))+
  scale_x_continuous(limits = c(0,60),expand = c(0,0))

error_training_set
error_test_set

## Choose from following to save your file accordingly
#ggsave("3.error_training_set_Wholeplot.jpg", error_training_set )
#ggsave("3.error_test_set_Wholeplot.eps", error_test_set )

#ggsave("2.error_training_set_Wholeplot.jpg", error_training_set )
#ggsave("2.error_test_set_Wholeplot.eps", error_test_set )

