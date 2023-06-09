#This file contains the commands for the processing of the raw data and the statistical analyses.
#The code for the final version of Figure 2 and Figure 2S are in a separate file.
Quantity<-read.csv(file ='QuantityVQuality_MainDataset.csv')

Quantity[,4]<-as.character(Quantity[,4])
Quantity[1468:1551,4]<-Quantity[1468:1551,44]

Quantity.a<-Quantity[which(Quantity$Adult.Weight>0),]

# Principal Component Analysis 

Caste.Quantity<-Quantity[1488:1550,]

PCA<-Caste.Quantity[,21:28]
PCA1<-prcomp(PCA[,1:7], scale=FALSE)

PCA2<-cbind(PCA1$x[,1:2],PCA[,8])
PCA2<-data.frame(PCA2)

Sup.Quantity<-Quantity[1:1487,]
Sup.Quantity<-Sup.Quantity[which(Sup.Quantity$Adult.Weight>0),]
Sup.Quantity<-Sup.Quantity[,21:28]

res.pca<-predict(PCA1, newdata=Sup.Quantity[,1:7])

res.pca2<-rbind(res.pca[,1:2], PCA2[,1:2])
res.pca2<-cbind(res.pca2,Quantity.a$Quantity,Quantity.a$Treatment,Quantity.a$Diet,
                Quantity.a$ProPerc,Quantity.a$CarbPerc,Quantity.a$WaterPerc, Quantity.a$p.c.ratio,Quantity.a$Hive)
res.pca2<-data.frame(res.pca2)
colnames(res.pca2)<-c("PC1","PC2","Quantity","color","diet","properc","carbperc","waterperc","pcratio","hive")

res.pca2$Quantity<-as.character(res.pca2$Quantity)
res.pca2$PC1<-as.numeric(res.pca2$PC1)
res.pca2$PC2<-as.numeric(res.pca2$PC2)


res.pca2$PC1<-res.pca2$PC1*-1
res.pca2$PC2<-res.pac2$PC2*-1

#to view the importance of components:
summary(PCA1)

#To view the eignvalues from the PCA. 

library(factoextra)

eig.val <- get_eigenvalue(PCA1)
eig.val
#for graph versions of eigenvalues
fviz_eig(PCA1)
#for graph versions individual/variables
fviz_pca(PCA1)

# Results for Variables
res.var <- get_pca_var(PCA1)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(PCA1)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

summary(PCA1)
library(psych)

pcacor<-cor(PCA[,1:7])
pcmodell<-principal(pcacor, rotate="none")

#Clustering Analysis
library(dplyr)
pca.clustering<-Quantity[,21:28]
pca.clustering<-pca.clustering[which(pca.clustering$Adult.Weight>0),]

clusters <- hclust(dist(pca.clustering), "complete")


plot(clusters)
summary(clusters)

#How we came up with 3 clusters as the optimum
k.max <- 15 # Maximal number of clusters
data <- pca.clustering[,-8]
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=10 )$tot.withinss})
summary(wss)

plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

#The elbow is the optimal number of clusters, at which within cluster variation 
#or within cluster sum of squares is minimized.
#The elbow method is implemented in factoextra package and can be computed using the function fviz_nbclust()

library(factoextra)
fviz_nbclust(data, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)

# You can use base graphics. A negative "hang" makes labels at the bottom even.
plot(clusters, labels = NULL, hang = -1,
     main = "Clustering of in vitro bees", sub = NULL,
     xlab = NULL, ylab = "Height")

#We used a package called sparcl to make color the dendrogram used in the paper. 
#But, in the newer versions of R, sparcl is no longer supported. 
#Here we have a similar verion using dendroextras

library(dendroextras)
#The 3 labels the dendrogram by three groups
colorden<-color_clusters(clusters, 3,groupLabels=TRUE)
dencolorleaf<-slice(colorden, h=1.9)
plot(colorden)
plot(dencolorleaf)


ColorDendrogram(clusters,y=pca.clustering$Treatment,branchlength=.8,ylab="Linkage Distance",xlab=NULL)
abline(h=1.9, col="red")
legend(322,2.7,legend=c("370ul","340ul","310ul","280ul","250ul","220ul","190ul","160ul",
                    "Ad-libitum","Q-Control","W-Control"), 
       col=c("red3","orange","yellow","green4","cyan","cadetblue","dodgerblue4","darkviolet",
             "deeppink","black","dimgray"), pch=19, cex=.5, box.col="black")
text(30,2.76,"Intercaste Cluster",cex=.7)
text(180,2,"Worker Cluster",cex=.7)
text(280,2,"Queen Cluster",cex=.7)

#Here is the K-means clustering analysis
k.clustering<-pca.clustering[1:347,1:8]

irisCluster <- kmeans(k.clustering[,1:7], 3)
irisCluster
table(irisCluster$cluster, k.clustering$Treatment)
summary(irisCluster)

#Here is the table organized by diet quality, not quantity. This required me to go back to
#the original data frama and remake a new version of pca.clustering that included the 'diet' column
pca.clustering.diet<-Quantity %>% select(3,21,22,23,24,25,26,27,28)
pca.clustering.diet<-pca.clustering.diet[which(pca.clustering.diet$Adult.Weight>0),]

k.clustering.diet<-pca.clustering.diet[1:347,1:9]
irisCluster.diet <- kmeans(k.clustering.diet[,2:8], 3)
irisCluster.diet
table(irisCluster.diet$cluster, k.clustering.diet$Diet)
summary(irisCluster.diet)

#What happens when you tell it to do 4 clusters?
k.clustering<-pca.clustering[1:347,1:8]

irisCluster <- kmeans(k.clustering[,1:7], 4)
irisCluster
table(irisCluster$cluster, k.clustering$Treatment)
summary(irisCluster)

#Quantity Line graph
res.pca2[,3]<-as.numeric(res.pca2[,3])
res.pca2<-res.pca2[1:269,]

ggplot(res.pca2, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=res.pca2$color,cex=2) +
  geom_smooth(method = "lm", se = TRUE,size=.5,color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30)) +
  theme_bw() + xlab("Total Consumption (ul)") + ylab("Principal Component 1 (64.42%)") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10))

linearModelVar <- lm(PC1 ~ Quantity, res.pca2)
summary(linearModelVar)

#Quality Line Graph

res.pca2[,6:9]<-as.numeric(res.pca2[,6:9])
res.pca2<-res.pca2[which(res.pca2$properc>.02),]
#Pro Per
p1<-ggplot(res.pca2, aes(x=properc, y=PC1, colour=Quantity)) +
  geom_point(color="black", cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5,color="black") +
  theme_bw() + xlab("Protein proportion") + ylab("Principal Component 1") +
  annotate("text", x=.035, y=.52, label= "A", size=5) +
  scale_x_continuous(breaks=seq(.02, .08, .01)) +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10))
print(p1)
linearModelVar <- lm(PC1 ~ properc, res.pca2)
summary(linearModelVar)

#Carb Perc
p2<-ggplot(res.pca2, aes(x=carbperc, y=PC1, colour=Quantity)) +
  geom_point(color="black",cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5,color="black") +
  theme_bw() + xlab("Carbohydrate Proportion") + ylab("Principal Component 1")+
  annotate("text", x=.135, y=.55, label= "B", size=5) +
  scale_x_continuous(breaks=seq(.1, .3, .05))+
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10))

print(p2)
linearModelVar <- lm(PC1 ~ carbperc, res.pca2)
summary(linearModelVar)

#water Perc
p3<-ggplot(res.pca2, aes(x=waterperc, y=PC1, colour=Quantity)) +
  geom_point(color="black",cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5, color="black") +
  theme_bw() + xlab("Water proportion") + ylab("Principal Component 1")+
  annotate("text", x=.6, y=.52, label= "C", size=5)+
  scale_x_continuous(breaks=seq(.6, .8, .05)) +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10))

print(p3)
linearModelVar <- lm(PC1 ~ waterperc, res.pca2)
summary(linearModelVar)

#P:C ratio
p4<-ggplot(res.pca2, aes(x=pcratio, y=PC1, colour=Quantity)) +
  geom_point(color="black",cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5,color="black") +
  theme_bw() + xlab("Protein:Carbohydrate Ratio") + ylab("Principal Component 1")+
  annotate("text", x=.16, y=.52, label= "D", size=5) +
  scale_x_continuous(breaks=seq(.1, .5, .05)) +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10))

print(p4)
linearModelVar <- lm(PC1 ~ pcratio, res.pca2)
summary(linearModelVar)

library(gridExtra)
grid.arrange(p1, p2,p3,p4, nrow=2, ncol=2)

#Diet Quality and Total Consumption

#HPHS
HPHS.pca<-res.pca2[which(res.pca2$diet=="HPHS"),]
HPHS<-ggplot(HPHS.pca, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=HPHS.pca$color,cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5, color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30)) +
  theme_bw() + xlab("")+ ylab("Principal Component 1") +
  ggtitle("A. High Protein-High Carb") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),plot.title=element_text(size=10))

print(HPHS)
linearModelVar <- lm(PC1 ~ Quantity, HPHS.pca)
summary(linearModelVar)

#HPMS
HPMS.pca<-res.pca2[which(res.pca2$diet=="HPMS"),]
HPMS<-ggplot(HPMS.pca, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=HPMS.pca$color,cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5,color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30))+
  theme_bw() + xlab("") + ylab("") +
  ggtitle("B. High Protein-Medium Carb") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),plot.title=element_text(size=10))

print(HPMS)
linearModelVar <- lm(PC1 ~ Quantity, HPMS.pca)
summary(linearModelVar)

#HPLS
HPLS.pca<-res.pca2[which(res.pca2$diet=="HPLS"),]
HPLS<-ggplot(HPLS.pca, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=HPLS.pca$color,cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5, color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30))+
  theme_bw() + xlab("") + ylab("") +
  ggtitle("C.High Protein-Low Carb") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),plot.title=element_text(size=10))

print(HPLS)
linearModelVar <- lm(PC1 ~ Quantity, HPLS.pca)
summary(linearModelVar)

#MPHS
MPHS.pca<-res.pca2[which(res.pca2$diet=="MPHS"),]
MPHS<-ggplot(MPHS.pca, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=MPHS.pca$color,cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5, color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30)) +
  theme_bw()+ xlab("Total Consumption (ul)") + ylab("Principal Component 1") +
  ggtitle("D. Medium Protein-High Carb") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),plot.title=element_text(size=10))

print(MPHS)
linearModelVar <- lm(PC1 ~ Quantity, MPHS.pca)
summary(linearModelVar)

#MPMS
MPMS.pca<-res.pca2[which(res.pca2$diet=="MPMS"),]
MPMS<-ggplot(MPMS.pca, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=MPMS.pca$color,cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5, color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30)) +
  theme_bw()+ xlab("Total Consumption (ul)") + ylab("") +
  ggtitle("E. Medium Protein-Medium Carb") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),plot.title=element_text(size=10))

print(MPMS)
linearModelVar <- lm(PC1 ~ Quantity, MPMS.pca)
summary(linearModelVar)

#MPLS
MPLS.pca<-res.pca2[which(res.pca2$diet=="MPLS"),]
MPLS<-ggplot(MPLS.pca, aes(x=Quantity, y=PC1, colour=Quantity)) +
  geom_point(color=MPLS.pca$color,cex=1) +
  geom_smooth(method = "lm", se = TRUE,size=.5, color="black") +
  scale_x_continuous(breaks=seq(160, 370, 30))+
  theme_bw() + xlab("Total Consumption (ul)") + ylab("") +
  ggtitle("F. Medium Protein-Low Carb") +
  theme(legend.text = element_text(size=10),legend.title=element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),plot.title=element_text(size=10))

print(MPLS)
linearModelVar <- lm(PC1 ~ Quantity, MPLS.pca)
summary(linearModelVar)

library(gridExtra)
grid.arrange(HPHS,HPMS,HPLS,MPHS,MPMS,MPLS, nrow=2, ncol=3)

#Below is the code for performing the GLMM

res.pca2$Quantity<-as.numeric(res.pca2$Quantity)
res.pca2$properc<-as.numeric(res.pca2$properc)
res.pca2$carbperc<-as.numeric(res.pca2$carbperc)
res.pca2$waterperc<-as.numeric(res.pca2$waterperc)



library(lme4)
library(lmerTest)
H<-lmer(PC1  ~ properc + carbperc + waterperc + Quantity + (1|hive),data=res.pca2)

summary(H)
library(sjPlot)
sjt.lmer(H,show.aic=TRUE,show.icc=TRUE, show.dev=TRUE,
         pred.labels = c("Dietary Protein Proportion",
                               "Dietary Carbohydrate Proportion",
                               "Dietary Water Proportion",
                               "Total Quantity"),
           depvar.labels = c("Principal Component Analaysis (59.02%)"))


summary(H)
anova(H)

library(MuMIn)
r.squaredGLMM(H)
