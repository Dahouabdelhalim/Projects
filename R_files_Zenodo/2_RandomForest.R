# This R script uses random forest models to determine which 
# floral traits differentiate primary floral visitors (Figure 3). 
# Method follows that of Dellinger et al. 2019 and 2021. It uses 
# the floral trait data "imputed.csv" and the data file 
# "Pollinator_observations.csv"

library(randomForest)
require(caTools)
library(sciplot) #bargraph with confidence intervals

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

data <- read.csv("imputed.csv")

#pollinator observations (for sub-setting data)
polobs <- read.csv("Pollinator_observations.csv")
polobs <- polobs[,c(1,3)]

setdiff(polobs$tip_label, data$sp_tip_label)
setdiff(data$sp_tip_label,polobs$tip_label)

data <- merge(data, polobs, by.x="sp_tip_label", by.y="tip_label", all.x=T, all.y=F)

#subset traits to remove those highly correlated (determined in previous analysis)
data.sub <- data[,c(1,2,5,6,7,8,11,12,13,14,15,16,17,18,19,20)]

#taxa with pollinator obs
data.obs <- subset(data, observed=="yes") #28 species with obs
traits <- data.obs[,2:19]

#taxa missing pollinator obs (note that NA obs are excluded: these are african taxa)
data.missing <- subset(data,  is.na(observed)) #21 missing pollinator obs + 3 outgroups
traits.missing <- data.missing[,2:18]

#get gini coef., robustness, and predicted syndromes for 100 random forests
robust <- c()
predict <- c()
gini <- c()

for (i in 1:1000) {
rf <- randomForest(
  as.factor(Syndrome) ~ .,
  data=traits,
  ntree=500,
  mtry=4, #this is equal to sqrt(number traits), see help file
)
myrobust <- data.obs$Syndrome == rf$predicted #to test robustness, examine whether a species was correctly classified in the training models
robust <- cbind(robust, myrobust)
mygini <- as.character(rf$importance)
gini <- cbind(gini, mygini)
mysyndrome <- as.character(predict(rf, newdata=traits.missing)) #predict syndrome for missing taxa
mypredict <- mysyndrome[1:21] == data.missing$Syndrome
predict <- cbind(predict, mypredict)
}

# proportion of 1000 RFs where syndrome was correctly determined
prop <- c()
for (i in 1:nrow(robust)) {
  myprop <- table(robust[i,])["TRUE"]/1000
  prop <- c(prop, myprop)
}
sum.robust <- as.data.frame(cbind(data.obs$sp_tip_label,as.numeric(prop)))
colnames(sum.robust) <- c("sp_tip_label", "prop_syn_correct")
sum.robust #proportion RF for taxa with pollinator obs where model correctly predicted pollinator

# proportion of 1000 RFs where predicted syndrome matches apriori predictions for taxa missing pollinator observations
prop <- c()
for (i in 1:nrow(predict)) {
  myprop <- table(predict[i,])["TRUE"]/1000 #add line here so if NA set value to 0
  prop <- c(prop, myprop)
}
sum.predict <- as.data.frame(cbind(data.missing$sp_tip_label,as.numeric(prop)))
colnames(sum.predict) <- c("sp_tip_label", "prop_syn_matched_apriori")
sum.predict[c(19,14,8),2] <- 0 #NA value is 0 (was NA because 0 matches with apriori 0/1000 = NA)
sum.predict #predicted pollinators for 21 taxa missing observations

# get mean and SD gini coefficients
mean <- c()
sd <- c()
for (i in 1:nrow(gini)) {
  mymean <- mean(as.numeric(gini[i,]))
  mysd <- sd(as.numeric(gini[i,]))
  mean <- c(mean, mymean)
  sd <- c(sd, mysd)
}
sum.gini <- as.data.frame(cbind(names(traits[1:17]),mean,sd))
colnames(sum.gini) <- c("trait", "gini.mean", "gini.sd")
sum.gini #ranked importance of each floral trait to predicting syndrome

## examine key traits by syndrome for 28 taxa with pollinator obs
table(data.obs[,c("Syndrome","bract_color")])
table(data.obs[,c("Syndrome","yellow_labstripe")])
table(data.obs[,c("Syndrome","red_labstripe")])
b <- subset(data.obs, Syndrome=="bee")
h <- subset(data.obs, Syndrome=="hummingbird")
mean(b$Labellum_Length)/mean(h$Labellum_Length) #bee taxa 2.03 times size of hummer taxa on average
mean(b$Labellum_Width)/mean(h$Labellum_Width) #bee taxa 3.27 times size of hummer taxa on average
mean(b$Stamen_Width)/mean(h$Stamen_Width) #bee taxa 1.6 times size of hummer taxa on average

#number predicted to be bee vs bird
tmp<-merge(sum.predict,data.missing[,c(1,19)], by="sp_tip_label")
tmp$Syndrome[c(19,14,8)] = c("hummingbird","hummingbird","hummingbird")
table(tmp$Syndrome)

#######
## Figure gini
#######
#relabel traits and order by mean gini value
sum.gini$trait <- c("corolla length", "corolla lobe length","corolla tube length","stamen exsertion","labellum length","labellum width","stamen length","stamen width","anther length", "style length","bract color","corolla color","labellum color","yellow labellum stripe", "red labellum stripe","stamen color","stamen tip color")
sum.gini <- sum.gini[order(sum.gini$gini.mean,decreasing = TRUE),]

pdf(file="Fig3.pdf", width=4,height=4) # specifications for your 
par (mfrow=c(1,1), mar=c(9,5,1,1))
base_r_barplot <- barplot(as.numeric(sum.gini$gini.mean),names.arg = sum.gini$trait, ylim = c(0, 3),xaxt="none", cex.lab=0.8,
                          ylab="Mean gini index")
axis(1, at=base_r_barplot,labels=sum.gini$trait, las=2, font=1, cex.axis=0.8, tick=F, hadj=1)
arrows(x0 = base_r_barplot, # Add error bars
         y0 = as.numeric(sum.gini$gini.mean) + as.numeric(sum.gini$gini.sd),
         y1 = as.numeric(sum.gini$gini.mean) - as.numeric(sum.gini$gini.sd),
         angle = 90,
         code = 3,
         length = 0)
dev.off()