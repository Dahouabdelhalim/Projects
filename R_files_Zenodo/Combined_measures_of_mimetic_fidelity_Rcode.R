library(MASS)
library(dplyr)
library(tidyr)
library(tidyverse)
library(PerformanceAnalytics)
library(hda)
library(irr)
library(ggplot2)


setwd("C:/Users/tdixi/Documents/Granularity")
AllEggs <- read.csv("Prinia_CF_DataSheet_alleggs.csv")
names(AllEggs)
AllEggs <- AllEggs[complete.cases(AllEggs), ] 
#select one row per Female ID 
AllEggs <- AllEggs %>% mutate(Year_Nest_Species = paste(Year, Nest, Species, sep = "_"))  %>%  distinct(Year_Nest_Species, .keep_all = TRUE) 

#make Year_Nest_Species the row label
AllEggs <- AllEggs %>% remove_rownames %>% column_to_rownames(var="Year_Nest_Species")
#correlation between variables
Corr <- AllEggs %>% dplyr::select(u, s, m, l, lum, FN, FS, FSD, PC_a, PD_a, OldEmax_ac,OldEmax1_ac,OldEprop_ac,OldEtot_ac,OldEsd_ac,
                                  MF_Area,MF_Perim,MF_Euler,MF_Psqovera,MaxMF_Area,MaxMF_Perim,MaxMF_Euler,MaxMF_Psqovera)
chart.Correlation(Corr, method = "spearman")

#PCA of colour and luminance
allpredictors<-AllEggs[c("u","s","m","l","lum")]
str(allpredictors) #confirm that all are numeric
allpredictors<-na.omit(allpredictors)
prin_comp <- prcomp(allpredictors, center = TRUE, scale. = T) #run pca on all data, ensuring it is scaled
names(prin_comp) #returns "sdev"     "rotation" "center"   "scale"    "x"    
#outputs the mean of variables
prin_comp$center
#outputs the standard deviation of variables
prin_comp$scale
#PC loading
prin_comp$rotation
dim(prin_comp$x) # dimensions should be (n-1)x(p) - i.e. the dimensions of the dataframe
biplot(prin_comp, scale = 0) 

#standard deviation of each principal component
pc_std_dev <- prin_comp$sdev
#variance of each component
pc_var <- pc_std_dev^2 # nb these should obviously descend
#variance of first 5 components
pc_var[1:5]
#find proportion variance explained by first x components
#prop_varexplained <- pc_var/sum(pc_var)
prop_varex <- pc_var/sum(pc_var)
#e.g. first 5 components
prop_varex[1:5]
prop_varex[1]
#use a scree plot to see which components explain most variation in data
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")
#first PC explains 60% of variation. First two explain >90%
#and PCs to dataframe
AllEggs$PC1col<-prin_comp$x[,1]
AllEggs$PC2col<-prin_comp$x[,2]

#PCA of NPM variables
allpredictors<-AllEggs[c("FN","FS","FSD")]
str(allpredictors) #confirm that all are numeric
allpredictors<-na.omit(allpredictors)
prin_comp <- prcomp(allpredictors, center = TRUE, scale. = T) #run pca on all data, ensuring it is scaled
names(prin_comp) #returns "sdev"     "rotation" "center"   "scale"    "x"    
#outputs the mean of variables
prin_comp$center
#outputs the standard deviation of variables
prin_comp$scale
#PC loading
prin_comp$rotation
dim(prin_comp$x) # dimensions should be (n-1)x(p) - i.e. the dimensions of the dataframe
biplot(prin_comp, scale = 0) 

#standard deviation of each principal component
pc_std_dev <- prin_comp$sdev
#variance of each component
pc_var <- pc_std_dev^2 # nb these should obviously descend
#variance of first 5 components
pc_var[1:5]
#find proportion variance explained by first x components
#prop_varexplained <- pc_var/sum(pc_var)
prop_varex <- pc_var/sum(pc_var)
#e.g. first 5 components
prop_varex[1:5]
prop_varex[1]
#use a scree plot to see which components explain most variation in data
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")
#1st dim is 72% of var
#and PCs to dataframe
AllEggs$PC1NPM<-prin_comp$x[,1] #FS is mostly PCA1
AllEggs$PC2NPM<-prin_comp$x[,2]

Corr <- AllEggs %>% dplyr::select(u, s, m, l, lum, PC1NPM,PC2NPM, PC_a, PD_a, OldEprop_ac,OldEtot_ac,
                                  MF_Area,MF_Perim,MF_Euler,MF_Psqovera,MaxMF_Area,MaxMF_Perim,MaxMF_Euler,MaxMF_Psqovera)
chart.Correlation(Corr, method = "spearman")
Corr <- AllEggs %>% dplyr::select(PC1col,PC2col, PC1NPM,PC2NPM, PC_a, PD_a, OldEprop_ac,OldEtot_ac,
                                  MF_Area,MF_Perim,MF_Euler,MF_Psqovera,MaxMF_Area,MaxMF_Perim,MaxMF_Euler,MaxMF_Psqovera)
chart.Correlation(Corr, method = "spearman")
#MF_PSqovera highly correlated with PC1NPM(which itself is mostly FS), makes sense as P^2/A and FS both measure size in some way?



#check normality for all variables to determine which discriminant analysis to use
library(ggpubr)
ggqqplot(AllEggs$PC1col)
ggqqplot(AllEggs$PC2col)
ggqqplot(AllEggs$PC1NPM)
ggqqplot(AllEggs$PC2NPM)
ggqqplot(AllEggs$MF_Euler)
ggqqplot(AllEggs$MF_Psqovera)
ggqqplot(AllEggs$OldEprop_ac)
ggqqplot(AllEggs$OldEtot_ac)
ggqqplot(AllEggs$PC_a)
ggqqplot(AllEggs$PD_a)


#can also do bartlett tests e.g.:
bartlett.test(logMF_Psqovera~Species,data=AllEggs)
#data is heteroschedastic - so use FDA

#### fda ####
library(mda)
#scale variables
AllEggs$PC1col<-scale(AllEggs$PC1col)
AllEggs$PC2col<-scale(AllEggs$PC2col)
AllEggs$OldEmax1_ac<-scale(AllEggs$OldEmax1_ac)
AllEggs$OldEprop_ac<-scale(AllEggs$OldEprop_ac)
AllEggs$OldEtot_ac<-scale(AllEggs$OldEtot_ac)
AllEggs$PC_a<-scale(AllEggs$PC_a)
AllEggs$PD_a<-scale(AllEggs$PD_a)
AllEggs$PC1NPM<-scale(AllEggs$PC1NPM)
AllEggs$PC2NPM<-scale(AllEggs$PC2NPM)
AllEggs$MF_Euler<-scale(AllEggs$MF_Euler)
AllEggs$MF_Psqovera<-scale(AllEggs$MF_Psqovera)

#fda
ffit4 <- fda(as.factor(Species) ~ PC1col + PC2col+ OldEmax1_ac + 
               OldEprop_ac + OldEtot_ac + PC_a + PD_a + 
               PC1NPM +PC2NPM + MF_Euler +MF_Psqovera
             , prior = c(1,1)/2, data=AllEggs, CV = TRUE)
#predict
predict(ffit4)
pred<-predict(ffit4,prior = c(1,1)/2)
predarranged<-data.frame(AllEggs$Species,pred)
predarranged$correct <- ifelse(predarranged$AllEggs.Species == predarranged$pred, 0, 1)
sum(predarranged$correct) #37 errors
predarrangedCF<-subset(predarranged,predarranged$AllEggs.Species=="CF")
sum(predarrangedCF$correct) #11 errors. Therefore prinia errors = 26
#therefore confusion matrix is 112,26,11,189

#confusion
ffit4$confusion
#get 106,20,17,195

#coefficients
ffit4$fit
#overall result is the MFs do add more info - Euler and PSqovera are actually most important out of all variables!
round(ffit4$fit$coefficients,2)


ffit4$percent.explained #100 because only 2 factors

## for 1,1 predict function
ct4 <- table(predarranged[,2:1])
colnames(ct4) <- c("CF","Prinia")
rownames(ct4) <- c("CF","Prinia")
ct4
diag(prop.table(ct4, 1)) #     CF    Prinia 
#                         0.8115942  0.9450000  
sum(ct4[row(ct4) == col(ct4)]) / sum(ct4) #0.8846154
plot(ct4, main = "", xlab = "Observation species", ylab = "Species assigned", col = c("tomato", "darkgrey"))

smoke <- matrix(c(112,26,11,189),ncol=2,byrow=TRUE) #for 1,1 predict function
colnames(smoke) <- c("CF","Prinia")
rownames(smoke) <- c("CF","Prinia")
smoke
fisher.test(smoke)


#boxplot MF_Psqovera and Euler#
ggplot(AllEggs,aes(Species,MF_Psqovera))+
  geom_boxplot(fill="grey",outlier.shape = NA)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()+
          theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14)))+
  geom_point(position="jitter")+
  labs(x="Species",
       y="MF P^2/A")+
  scale_x_discrete(limits=c("CF","P"),labels=c("Cuckoo finch","Prinia"))

ggplot(AllEggs,aes(Species,MF_Euler))+
  geom_boxplot(fill="grey",outlier.shape = NA)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()+
          theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14)))+
  geom_point(position="jitter")+
  labs(x="Species",
       y="MF Euler")+
  scale_x_discrete(limits=c("CF","P"),labels=c("Cuckoo finch","Prinia"))


#### cuckoo finch experiments - test if MFs or treatment predicts rejection ####
setwd("C:/Users/tdixi/Documents/Granularity")
CFPSexp<-read.csv("CFPSexperiments_final.csv",header=T,sep=",")

#test whether addition of squiggles affects rejection, and if MFs predict rejection
CFPSexp$Treatment<-as.factor(CFPSexp$Treatment)
CFPSexp$Treatment <- relevel(CFPSexp$Treatment, ref = "W")
tc1<-glm(Egg_rejected ~ Treatment, data=CFPSexp, family = binomial(link = "logit"))
summary(tc1)
#Treatment does not predict rejection

##number of rejections plot for squiggles, blotches, water ##
CFPSexpS<-subset(CFPSexp,CFPSexp$Treatment=="S") #rejections = 11/17
CFPSexpB<-subset(CFPSexp,CFPSexp$Treatment=="B") #rejections = 8/13
CFPSexpW<-subset(CFPSexp,CFPSexp$Treatment=="W") #rejections = 8/14
ybar<-as.numeric(c(11/17,8/14,8/13))
xbar<-c("Scribbles","Blotches","Water")
bardata<-as.data.frame(cbind(xbar,ybar))
colnames(bardata)<-c("Treatment","Rejection_frequency")
bardata[,2]<-as.numeric(bardata[,2])
bardata$Treatment <- factor(bardata$Treatment,                                    # Change ordering manually
                            levels = c("Scribbles","Blotches","Water"))


ggplot(data=bardata, aes(Treatment,Rejection_frequency))+
  geom_bar(stat="identity",width = 0.8, position = position_dodge())+
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()+
          theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14)))+
  ylim(0,1)+
  ylab("Rejection Frequency")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

#### check if squiggles increase intMFs ####
setwd("C:/Users/tdixi/Documents/Granularity")
CFmanip<-read.csv("IntMFs_CFs_pre_manip_comparison.csv")

CFmanipS<-subset(CFmanip,CFmanip$SpeciesSBW=="S")
CFmanipB<-subset(CFmanip,CFmanip$SpeciesSBW=="B")
CFmanipW<-subset(CFmanip,CFmanip$SpeciesSBW=="W")

t.test(CFmanipS$Sum_Psqovera,CFmanipS$Sum_Psqoveramanip, paired = TRUE, alternative = "two.sided")

#plot psqovera for pre- and post-painted eggs
CFmanipscribbles<-read.csv("IntMFs_CFs_pre_manip_comparison_scribbles.csv")

ggplot(CFmanipscribbles,aes(BeforePost,Psqovera))+
  geom_boxplot(fill="grey",outlier.shape = NA)+theme_bw()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  geom_point(position="jitter")+
  labs(x="Before or after painting with scribbles",
       y="P^2/A")+
  scale_x_discrete(limits=c("Before","Post"),labels=c("Before","After"))

