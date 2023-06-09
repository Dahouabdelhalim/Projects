

### R analyses for Euglossa CHC manuscript. Accompanying data is found in the supplemental files. ####


#load necessary packages
library(ecodist)
library(vegan)
library(gplots)
library(MASS)
library(ggplot2)
library(corrplot)
library(car)
library(Steel.Dwass.test)
library(ggpubr)
library(RVAideMemoire)

############# NMDS for all populations, with 4 outliers removed (data in supplemental files) ###########################


#read in and rename data
CHC = read.csv("All_Bees_CR_FL_MX.csv")

#check data import
dim(CHC)
head(CHC)

#find relative abundances (i.e. percents)
rsums <- rowSums (CHC[,c(5:15)])
CHC.norm <- CHC[,c(5:15)]/rsums
rowSums(CHC.norm)
head(CHC.norm)

#add back in descriptive columns

CHC <- cbind(CHC[,1:4], CHC.norm)
head(CHC)


#heatmap
heatmap.2(as.matrix(CHC[,5:15]), Colv=NA, col=terrain.colors(250), labRow=paste (CHC$Sample,sep="_"), main="Female CHCs", cexRow=2, cexCol=.5, trace="none")


#nMDS plots

#create triangular distance matrix
bc.dist <- bcdist(CHC[,5:15]) 
bc.dist

# nMDS,note: may take a long time to run. Running at 10 iterations is faster and provides very similar grouping with comparable (low) stress value.
CHC
CHC.nmds <- nmds (bc.dist, mindim=2, maxdim=2, nits=100)
CHC.nmin <- nmds.min (CHC.nmds)

CHC.nmin <- cbind(CHC.nmin, CHC$Sample, CHC$Behavior, CHC$Species_Pop_Sex, CHC$chemotype)
CHC.nmin
head(CHC.nmin)
colnames (CHC.nmin)[3]<- "Sample"
colnames (CHC.nmin)[4]<- "Behavior"
colnames(CHC.nmin)[5] <- "Species_Pop_Sex"
colnames(CHC.nmin) [6] <- "chemotype"
head(CHC.nmin)

#quick plot of  nMDS

colors <- c("blue","red","purple","green", "orange", "olivedrab","cyan","firebrick","gray","magenta","lightgreen", "yellow", "tan", "thistle", "deepskyblue4","limegreen","wheat", "sienna","orchid4","black")

plot(CHC.nmin$X1, CHC.nmin$X2, cex=1.2, main="Female CHCs",col=colors[factor(CHC.nmin$Behavior)],pch=c(8,17,16,15,18,14,19,5,6,7)[factor(CHC.nmin$Behavior)], xlab="", ylab="")

#put sample name labels on points
text(CHC.nmin$X1, CHC.nmin$X2, pos=1, labels=CHC.nmin$Sample, cex=0.2)

# output NMDS coordinates for ggplot2 plotting
#write.table(CHC.nmin,"NMDS_all_ggplot.txt")

# Simper Analysis of chemotypes

dim(CHC)
comm = CHC[ ,5:15]
simpcomm= simper(comm, group = CHC$chemotype)
summary(simpcomm)
lapply(simpcomm, FUN=function(x){x$overall})


# ggplot for nicer NMDS visualization based on previous NMDS output

# rename NMDS matrix for ggplot2
NMDS_all_ggplot = CHC.nmin

NMDS_all_ggplot$Behavior = as.factor(NMDS_all_ggplot$Behavior)
gp =ggplot(NMDS_all_ggplot, aes(x=NMDS_all_ggplot$X1, y=NMDS_all_ggplot$X2)) +
  geom_point(aes(fill = NMDS_all_ggplot$Behavior ,shape = NMDS_all_ggplot$Behavior, color= NMDS_all_ggplot$Behavior), size = 3) +
  scale_shape_manual(values = c(21,22,10,7,24,23,9,2))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")
gp

gp + scale_fill_manual(breaks = c("1", "2", "3", "4","5", "6","7","8"),
                       values=c("red", "blue", "red3", "midnightblue", "palegreen4", "purple", "darkorchid4", "darkgreen"))+
  scale_color_manual(breaks = c("1", "2", "3", "4","5", "6","7","8"),
                     values=c("black", "black", "red3", "midnightblue", "black", "black", "darkorchid4", "darkgreen"))

######################## NMDS For SOcial Behavior #####################################



# read in and rename data
CHC = read.csv("Social_FL_NewlyEmerged.csv")

# note: in data file behavior_num 1 = dominant, 2 = subordinate, 3 = foundress, 4 = guard, 5 = newly emerged

#check that data read in
dim(CHC)
head(CHC)

#find relative abundances (i.e. percents)
rsums <- rowSums (CHC[,c(4:20)])
CHC.norm <- CHC[,c(4:20)]/rsums
rowSums(CHC.norm)
head(CHC.norm)

#add back in descriptive columns
CHC <- cbind(CHC[,1:3], CHC.norm)
head(CHC)


#heatmap
heatmap.2(as.matrix(CHC[,4:20]), Colv=NA, col=terrain.colors(250), labRow=paste (CHC$Sample,sep="_"), main="Female CHCs", cexRow=2, cexCol=.5, trace="none")


#nMDS plots

#create triangular distance matrix
bc.dist <- bcdist(CHC[,4:20]) 
bc.dist

# nMDS, note: may take a long time to run. Running at 10 iterations is faster and provides very similar grouping with comparable (low) stress value.
CHC
CHC.nmds <- nmds (bc.dist, mindim=2, maxdim=2, nits=100)
CHC.nmin <- nmds.min (CHC.nmds)

CHC.nmin <- cbind(CHC.nmin, CHC$Sample, CHC$Behavior_num, CHC$Behavior_cat)
CHC.nmin
head(CHC.nmin)
colnames (CHC.nmin)[3]<- "Sample"
colnames (CHC.nmin)[4]<- "Behavior_num"
colnames(CHC.nmin)[5] <- "Behavior_cat"
head(CHC.nmin)

# quick plot of nMDS

#symbols <- c(pch=16,15,17,18,3)
colors <- c("blue","red","purple","green", "orange", "olivedrab","cyan","firebrick","gray","magenta","lightgreen", "yellow", "tan", "thistle", "deepskyblue4","limegreen","wheat", "sienna","orchid4","black")
plot(CHC.nmin$X1, CHC.nmin$X2, cex=1, main="Female CHCs",col=colors[factor(CHC.nmin$Behavior_num)],pch=c(8,17,16,15,18,14,19,5,6,7)[factor(CHC.nmin$Behavior_num)], xlab="", ylab="")

# put text labels of sample ID on points
text(CHC.nmin$X1, CHC.nmin$X2, pos=1, labels=CHC.nmin$Sample, cex=0.2)

# output NMDS coordinates for ggplot2 plotting
#write.table(CHC.nmin,"Social_NMDS.txt")

#SimperAnalysis for contribution of peaks to behavior

dim(CHC)
comm = CHC[ ,4:20]
simpcomm= simper(comm, group = CHC$Behavior_num)
summary(simpcomm)
lapply(simpcomm, FUN=function(x){x$overall})

# ggplot for nicer NMDS plot of social individuals, based on output of above NMDS 

#rename NMDS matrix for ggplot2
Social_NMDS = CHC.nmin

Social_NMDS$Behavior_cat = as.factor(Social_NMDS$Behavior_cat)
gp =ggplot(Social_NMDS, aes(x=Social_NMDS$X1, y=Social_NMDS$X2)) +
  geom_point(aes(fill = Social_NMDS$Behavior_cat ,shape = Social_NMDS$Behavior_cat), size = 3) +
  scale_shape_manual(values = c(21,24,23,22,25))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")
gp

gp + scale_fill_manual(breaks = c("1", "2", "3", "4","5"),
                       values=c("red", "blue", "palegreen4", "purple", "light blue"))



################# Correlation Plots made with FL females (foundress, guard, dominant, subordinate) ###########

#read in and rename data
CHC_social_correlation_matrix = read.csv("CHC_social_correlation_matrix.csv")

#first check the best fit for # of clusters (3 is cleary the best fit here)

# make correlation plot with three possible clusters

m = cor(CHC_social_correlation_matrix, method = "spearman")
corrplot(m, method = 'circle', order = "hclust", addrect = 3, hclust.method = "ward.D2")

# make correlation plot with two possible clusters

m = cor(CHC_social_correlation_matrix, method = "spearman")
corrplot(m, method = 'circle', order = "hclust", addrect = 2, hclust.method = "ward.D2")

# make correlation plot with 4 possible clusters

m = cor(CHC_social_correlation_matrix, method = "spearman")
corrplot(m, method = 'circle', order = "hclust", addrect = 4, hclust.method = "ward.D2")


############################CHC Peak Anovas, boxplots, regression #####################



OvData = read.csv("Social_chemotypes_ovarysize.csv")

# lettered variables in data sheet: a = foundress, b = guard, c = dominant, d = subordinate ##
# numbered variables in data sheet: 1 = foundress, 2 = guard, 3 = dominant, 4 = subordinate


OI = OvData$OvaryIndex
Behav = OvData$Behavior_lettered
test = OvData$module_three

## homogeneity of variances and normality ##
leveneTest(test~Behav, data =OvData)
shapiro.test(test)

## Fails normality so test transformed ##
peak = sqrt(OvData$module_three)
shapiro.test(peak)


#Peaks vs Behavior
fit <- aov(peak ~ Behav , data=OvData)
summary(fit)
TukeyHSD(fit)


#module-three Boxplots 

peak = OvData$module_three
Behav = OvData$Behavior_lettered

gp = ggplot(OvData, aes(x = Behav, y = peak))+
  geom_boxplot(width=0.25)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp + geom_jitter(shape=16, position=position_jitter(0.05)) + ylim(0.1,0.35)
# note: outlier points are duplicated on this figure due to jittering and the duplicate values were removed in final figure editing

#### module-three vs Ovary Size ####

gp =ggplot(OvData, aes(x=OvData$OvaryIndex, y=OvData$module_three)) +
  geom_point(aes(fill = OvData$Behavior_lettered ,shape = OvData$Behavior_lettered), size = 4) +
  scale_shape_manual(values = c(23,22,21,24))+
  geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(method = "spearman", color = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("a", "b", "c", "d"),
                       values=c("palegreen4", "purple", "red", "blue"))

# add data labels
gp +geom_text(aes(label=OvData$Sample),hjust=0, vjust=0)


########################### Repeat analysis With only Chemotype A individuals ##############################


OvData = read.csv("Chemotype_A_Social_ovaries.csv")
# module-three check w/in chemotype A individuals 
# lettered variables in data sheet: a = foundress, b = guard, c = dominant, d = subordinate
# numbered variables in data sheet: 1 = foundress, 2 = guard, 3 = dominant, 4 = subordinate

OI = OvData$OvaryIndex
Behav = OvData$Behavior_lettered
peak = OvData$module_three


## homogeneity of variances and normality ##
leveneTest(peak~Behav, data =OvData)
shapiro.test(OvData$module_three)

## Passes both, so continue

#Peaks vs Behavior
fit <- aov(peak ~ Behav , data=OvData)
summary(fit)
TukeyHSD(fit)

## boxplot of within chemotype A only, signal peaks. Replicates finding of larger dataset ###

gp = ggplot(OvData, aes(x = Behav, y = peak))+
  geom_boxplot(width=0.25)+
 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp + geom_jitter(shape=16, position=position_jitter(0.05))

# module-three vs Ovary Size for only chemotype A individuals #

gp =ggplot(OvData, aes(x=OvData$OvaryIndex, y=OvData$module_three)) +
  geom_point(aes(fill = OvData$Behavior_lettered ,shape = OvData$Behavior_lettered), size = 4) +
  scale_shape_manual(values = c(23,22,21,24))+
  geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(method = "pearson", color = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("a", "b", "c", "d"),
                       values=c("palegreen4", "purple", "red", "blue"))

# add sample labels
gp +geom_text(aes(label=OvData$Sample),hjust=0, vjust=0)

# chemotype A check w/in chemotype A individuals

OI = OvData$OvaryIndex
Behav = OvData$Behavior_lettered
peak = OvData$chemotype_A_peaks

## homogeneity of variances and normality ##
leveneTest(peak~Behav, data =OvData)
shapiro.test(OvData$chemotype_A_peaks)

# fails normality so sqrt transform for stats 

sqrtChemoA = sqrt(OvData$chemotype_A_peaks)
leveneTest(sqrtChemoA~Behav, data =OvData)
shapiro.test(sqrtChemoA)

#chemotype A still failed normality after transformation, so move forward with kruskal-wallis on untransformed data

kruskal.test(peak~Behav, data = OvData)
Steel.Dwass(OvData$chemotype_A_peaks, group = OvData$Behavior_numbered )

#ggplot boxplot of chemotype A peaks with chemotype A only individuals

gp = ggplot(OvData, aes(x = Behav, y = peak))+
  geom_boxplot(width=0.25)+
  #geom_point(aes())+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp + geom_jitter(shape=16, position=position_jitter(0.05))
# note: outlier points are duplicated on this figure due to jittering and the duplicate values were removed in final figure editing

# CHC Module 1 (diagnostic of chemotype-A) vs Ovary Size ###

gp =ggplot(OvData, aes(x=OvData$OvaryIndex, y=OvData$chemotype_A_peaks)) +
  geom_point(aes(fill = OvData$Behavior_lettered ,shape = OvData$Behavior_lettered), size = 4) +
  scale_shape_manual(values = c(23,22,21,24))+
  geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(method = "spearman", color = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("a", "b", "c", "d"),
                       values=c("palegreen4", "purple", "red", "blue"))


# chemotype B check w/in chemotype A individuals

OI = OvData$OvaryIndex
Behav = OvData$Behavior_lettered
peak = OvData$chemotype_B_peaks

## homogeneity of variances and normality ##
leveneTest(peak~Behav, data =OvData)
shapiro.test(peak)

# fails normality so proceed with sqrt transformation

sqrtChemoB = sqrt(peak)
leveneTest(sqrtChemoB~Behav, data = OvData)
shapiro.test(sqrtChemoB)

# fails homogeneity of variances, proceed with kruskal wallis. 

kruskal.test(peak~Behav, data = OvData)
Steel.Dwass(OvData$chemotype_B_peaks, group = OvData$Behavior_numbered )

# ggplot boxplot of chemotype B peaks in chemotype A individuals

gp = ggplot(OvData, aes(x = Behav, y = peak))+
  geom_boxplot(width=0.25)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp + geom_jitter(shape=16, position=position_jitter(0.05))
# note: outlier points are duplicated on this figure due to jittering and the duplicate values were removed in final figure editing

# CHC Module 2 (diagnostic of chemotype-B) vs Ovary Size ###

gp =ggplot(OvData, aes(x=OvData$OvaryIndex, y=OvData$chemotype_B_peaks)) +
  geom_point(aes(fill = OvData$Behavior_lettered ,shape = OvData$Behavior_lettered), size = 4) +
  scale_shape_manual(values = c(23,22,21,24))+
  geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(method = "spearman", color = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("a", "b", "c", "d"),
                       values=c("palegreen4", "purple", "red", "blue"))


################# Gene Expression CHC Correlations ############################


####BRAIN ##########

BrainData = read.csv("Brains_CHCs_Genes.csv")
  
#individual correlations can be explored by changing the XY columns in the code below

gp =ggplot(BrainData, aes(x=BrainData$Edil_04295, y=BrainData$module_three)) +
  geom_point(aes( fill=BrainData$Behavior, shape = BrainData$Behavior),size = 4) +
  scale_shape_manual(values = c(23,22,21,24))+
  geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(method = "spearman", color = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("a", "b", "c", "d"),
                       values=c("palegreen4", "purple", "red", "blue"))



############# Ovaries ##################################


OvaryData = read.csv("Ovaries_CHCs_Genes.csv")

#individual correlations can be explored by changing the XY columns in the code below

  gp =ggplot(OvaryData, aes(x=OvaryData$Edil_04108, y=OvaryData$module_three)) +
  geom_point(aes( fill = OvaryData$Behavior, shape = OvaryData$Behavior), size = 4) +
  scale_shape_manual(values = c(23,22,21,24))+
  geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(method = "spearman", color = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("a", "b", "c", "d"),
                       values=c("palegreen4", "purple", "red", "blue"))
gp +geom_text(aes(label=OvaryData$Sample),hjust=0, vjust=0)


########### Relatedness Histograms for each nest, relative to mean relatedness of all samples #############################


RelatednessEstimates = read.csv("RelatednessEstimates.csv")


#### Nest 100 control ######
hist(RelatednessEstimates$Wang_Statistic, main = "Control Nest ID:100 \\n Mother and Daughters \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nesmate relatedness
abline(v=0.371,col="red", lty =3)   
abline(v=0.404,col="red", lty =3)
abline(v=0.349,col="red", lty =5)
abline(v=0.374,col="red", lty =3)
abline(v=0.354,col="red", lty =5)
abline(v=0.293,col="red", lty =5)


#### Nest 13 control ######
hist(RelatednessEstimates$Wang_Statistic, main = "Control Nest ID:13 \\n Mother and Daughters \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.51,col="red", lty =3)
abline(v=0.499,col="red", lty =3)
abline(v=0.472,col="red", lty =3)
abline(v=0.331,col="red", lty =5)
abline(v=0.443,col="red", lty =3)
abline(v=0.449,col="red", lty =3)
abline(v=0.453,col="red", lty =3)
abline(v=0.35,col="red", lty =5)
abline(v=0.324,col="red", lty =5)
abline(v=0.342,col="red", lty =5)

#### Nest 28B control ######
hist(RelatednessEstimates$Wang_Statistic, main = "Control Nest ID:28B \\n Non-kin \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.214,col="red", lty =6)

###### Nest 51B control #####
hist(RelatednessEstimates$Wang_Statistic, main = "Control Nest ID:51B \\n Mother and Daughter \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.318,col="red", lty =5)


###### Nest 51A ############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID: 51B \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.351,col="red", lty =3)
abline(v=0.355,col="red", lty =3)
abline(v=0.47,col="red", lty =3)
abline(v=0.292,col="red", lty =3)
abline(v=0.334,col="red", lty =3)
abline(v=0.352,col="red", lty =3)

######## Nest 59A ###############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:59A \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")

# nestmate relatedness
abline(v=0.351,col="red", lty =6)

########## Nest 59B ###########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:59B \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")

# nestmate relatedness
abline(v=0.349,col="red", lty =6)

########### Nest 75 ############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:75 \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

#nestmate relatedness
abline(v=0.26,col="red", lty =6)

########## Nest 82 ##############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:82 \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.409,col="red", lty =6)

########## Nest 5 ##############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:5 \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.175,col="red", lty =6)

####### Nest 47 ############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:47 \\n All Chemotype B", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

#nestmate relatedness
abline(v=0.376,col="red", lty =6)

####### Nest 44 ############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:44 \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$V6),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

#nestmate relatedness
abline(v=0.297,col="red", lty =6)

###### Nest 33 ############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:33 \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

#nestmate relatedness
abline(v=0.282,col="red", lty =6)
abline(v=0.385,col="red", lty =6)  ## chemotype Mix ##
abline(v=0.337,col="red", lty =6)  ## chemotype Mix ##

####### Nest 29 ############
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:29 \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.449,col="red", lty =6)
abline(v=0.392,col="red", lty =6)
abline(v=0.385,col="red", lty =6)

###### Nest 28A ##########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:28A \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.36,col="red", lty =6)
abline(v=0.43,col="red", lty =6)
abline(v=0.363,col="red", lty =6)
abline(v=0.379,col="red", lty =6)
abline(v=0.454,col="red", lty =6)

###### Nest 27 ##########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:27 \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.409,col="red", lty =6)
abline(v=0.464,col="red", lty =6)
abline(v=0.376,col="red", lty =6)

####### Nest 2 #################
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:2 \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))


#nestmate relatedness
abline(v=0.357,col="red", lty =6)

###### Nest 26 ##########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:26 \\n All Chemotype B", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")


# nestmate relatedness
abline(v=0.236,col="red", lty =6)

###### Nest 20 ##########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:20 \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.381,col="red", lty =6)
abline(v=0.311,col="red", lty =6)
abline(v=0.348,col="red", lty =6)

###### Nest 16A ##########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:16A \\n All Chemotype A", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.335,col="red", lty =6)

###### Nest 16B ##########
hist(RelatednessEstimates$Wang_Statistic, main = "Nest ID:16B \\n Mixed Chemotypes", xlab = "Relatedness Estimate", ylab = "Dyad Frequency")
abline(v=mean(RelatednessEstimates$Wang_Statistic),col="blue")
text(0.01, 300, paste( "Mean =",0.218, "\\n St. Deviation = 0.067", "\\n 90th Percentile = 0.296"))

# nestmate relatedness
abline(v=0.284,col="red", lty =6)





################# R code for supplemental comparison of body, abdomen,and wing CHC extractions #########################


# rename data
CHC = read.csv("body_abdomen_wing_extractions.csv")
#CHC = All_Bees_CR_FL
dim(CHC)
head(CHC)

#find relative abundances (i.e. percents)
rsums <- rowSums (CHC[,c(5:21)])
CHC.norm <- CHC[,c(5:21)]/rsums
rowSums(CHC.norm)
head(CHC.norm)

#add back in descriptive columns
CHC <- cbind(CHC[,1:4], CHC.norm)
head(CHC)
#write.csv(CHC, file= "body_abdomen_wing_histogram.csv")

#create triangular distance matrix
bc.dist <- bcdist(CHC[,5:21]) 
bc.dist

# nMDS
CHC
CHC.nmds <- nmds (bc.dist, mindim=2, maxdim=2, nits=10)
CHC.nmin <- nmds.min (CHC.nmds)

CHC.nmin <- cbind(CHC.nmin, CHC$Sample, CHC$method, CHC$individual, CHC$method_numbered)
CHC.nmin
head(CHC.nmin)
colnames (CHC.nmin)[3]<- "Sample"
colnames (CHC.nmin)[4]<- "method"
colnames(CHC.nmin)[5] <- "individual"
colnames(CHC.nmin) [6] <- "method_numbered"
head(CHC.nmin)

#plot nMDS


colors <- c("blue","red","purple","green", "orange", "olivedrab","cyan","firebrick","gray","magenta","lightgreen", "yellow", "tan", "thistle", "deepskyblue4","limegreen","wheat", "sienna","orchid4","black")

plot(CHC.nmin$X1, CHC.nmin$X2, cex=1, main="Female CHCs",col=colors[factor(CHC.nmin$method)],pch=c(8,17,16,15,18,14,19,5,6,7)[factor(CHC.nmin$method)], xlab="", ylab="")

# add sample names to points
text(CHC.nmin$X1, CHC.nmin$X2, pos=1, labels=CHC.nmin$Sample, cex=0.2)


# Simper Analysis  #

dim(CHC)
comm = CHC[ ,5:21]
simpcomm= simper(comm, group = CHC$method)
summary(simpcomm)
lapply(simpcomm, FUN=function(x){x$overall})


#### statistical and NMDS graphing approach based on Bruckner and Heethoff, 2016 in Chemoecology ########

my.data = read.csv("body_abdomen_wing_permanova.csv")
# Calculating the distance/dissimilarity matrix -> Euclidean/Bray-Curtis
distm1 <- vegdist(my.data[,-1], method="bray")
mean(distm1)  ### Mean distance of all data in the distance matrix 


### Permutational multivariate analysis of variance  -> PERMANOVA (Anderson 2001) ###
PERMANOVA <- adonis(distm1 ~ my.data$method, permutations=10000)
PERMANOVA

pairwise.perm.manova(distm1,my.data$method,nperm=10000, p.method = "fdr" )

mds1<-metaMDS(my.data[,-1],"bray",2)
mds1$stress

plot(mds1, type="n")
#mds1 = ordiplot3d(mds1, angle=45)
cols<-c("red", "blue", "darkorchid", "darkorange3", "chartreuse3","darkgreen","lightblue","darkred")
points(mds1, col=cols[my.data$method],pch=c(8,17,16,15,18,20,19,22)[unclass(my.data$method)], cex=1.5)

ef <- envfit(mds1, my.data[,-1], display="sites", permu=10000)
plot(ef, col="black", p.max=0.001, cex=0.75) 


######### run statistical tests on specific CHC differences

CHC = read.csv("body_abdomen_wing_histogram.csv")

#C21
peak= CHC$C21 
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C23:1a
peak= CHC$C23.1a 
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )


#C23:1b
peak= CHC$C231.b 
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C23
peak= CHC$C23
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C24.1
peak= CHC$C24.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C24
peak= CHC$C24
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C25.1
peak= CHC$C25.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C25
peak= CHC$C25
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C26.1
peak= CHC$C26.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C26
peak= CHC$C26
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C27.1
peak= CHC$C27.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C27
peak= CHC$C27
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C29.1
peak= CHC$C29.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C29
peak= CHC$C29
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C31.1
peak= CHC$C31.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C31
peak= CHC$C31
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

#C33.1
peak= CHC$C33.1
method = CHC$method_numbered

kruskal.test(peak~method, data = CHC)
Steel.Dwass(peak, group = method )

