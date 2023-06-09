#Franklin, Janet et al. Journal of Biogeography. "Geographical Ecology of Dry Forest Tree Communities in the West Indies" 
#Script includes ordination and clustering of tree species presence/absence data per main published results

#install packages if not installed already
install.packages(pkgs=c("BiodiversityR", "vegan", "Rcmdr", "MASS", "mgcv", "cluster", "RODBC", "rpart", 
                        "effects", "multcomp", "ellipse", "maptree", "sp", "splancs", "spatial", "akima", 
                        "nnet", "dismo", "raster", "rgdal", "gbm", "randomForest", "gam", "earth", "mda", 
                        "kernlab", "e1071", "sem", "rgl", "relimp", "lmtest", "leaps", "Hmisc", 
                        "colorspace", "aplpack", "abind", "XLConnect", "car", "markdown", "knitr"), 
                        dependencies=c("Depends", "Imports"))
install.packages("labdsv")
install.packages("dplyr")
install.packages("plyr")
library(plyr)
library(dplyr)
library(vegan)
library(labdsv)
library(BiodiversityR) 
library(ellipse)

#row.names = 1 replaces plot_id as the row label in first row and now it is row.names
com.data <- read.csv("WestIndiesSDTF_tree_species_PA_572sites_649taxa.csv", header=T, row.names=1)
#the archived environmental data WestIndiesSDTF_environment_572sites_with_descriptions.csv includes variable description in the first 28 records
#Delete these records and rename WestIndiesSDTF_environment_572sites_with_descriptions.csv
env <- read.csv("WestIndiesSDTF_environment_572sites.csv", header=T)

#calculate distaince matrix, Jaccard distance, from community presence-absence data
dist<-vegdist(com.data, binary=TRUE, method="jaccard")

#data summary plots by island as in Table S2
com.env.data <- com.data
com.env.data$island <- env$island
com.env.data$archipelago <- env$archipelago
#number of records per island and archipelago as in Table S2
table(com.env.data$island)
table(com.env.data$archipelago)

#Non-metric multidimensional scaling indirect ordination of community distance matrix for Figures 2, 3
nms.out=nmds(dist, k=2,maxit=50)

env.sub = env[,c(5:8,11,17:21)]
#subset to include long, lat (rounded to 2nd decimal place to protect plot locations
#and also note FIA plot locations are fuzzed),
#...alt, MAT, CV_T, Twarm, Tcold, MAP, Pwet, Pdry 

###FIGURE 3
#apply envfit to nms - WARNING - controversial. Oksanen says "You should not correlate environmental variables with axes in any ordination method."
(fit <- envfit(nms.out,env.sub,perm=999))
scores(fit,"vectors")
#Plot ordination showing sites on first 2 NMDS axes, label by archipelago
nms.plot<-ordiplot(nms.out,type="none",xlim=c(-2,1.5),ylim=c(-1,1),xlab="NMS1",ylab="NMS2")
plot(fit)
ordisymbol(nms.plot,y=env,factor="archipelago",rainbow=T, col=env,legend=T, with(legend="bottomright"))

# cluster community presence-absence distance matrix using Ward's linkage, cut to 14 clusters
iv.clust.ward<-hclust(dist,method = "ward.D2")
rect.hclust(iv.clust.ward, k=14)
iv.ward.cut<-cutree(iv.clust.ward,k=14)
# the cluster numbers is added to the env.sub frame
env.sub$ward = iv.ward.cut 
#plot the cluster dendrogram with sites labeled by their group number as shown in Figure 1b
plot(iv.clust.ward,label=env.sub$ward)
rect.hclust(iv.clust.ward, k=14)

###FIGURE 2 (BUT WITH MORE VECTORS SHOWN THAN IN RESUBMITTED MS)
#Plot NMDS ordination showing Ward's clusters and environmental vectors
nms.plot<-ordiplot(nms.out,type="none",xlim=c(-2,1.5),ylim=c(-1,1),xlab="NMS1",ylab="NMS2")
ordisymbol(nms.plot,y=env.sub,factor="ward",rainbow=T, col=env,legend=T)
(fit <- envfit(nms.out,env.sub[,1:10],perm=999))
scores(fit,"vectors")
plot(fit, p.max=0.096)

#Indicator Species Analysis applied to groups of sites formed by clustering 
#Information used in Tables 3 and S3
clusward_IS <- indval(com.data, iv.ward.cut)
summary(clusward_IS, p=0.05, digits=2, show=p,sort=FALSE, too.many=350)

#create crosstabluations of clusters by island and archipelago and write externally for Table S2.
ward14_arch <- table(env.sub$ward,env$archipelago)
cluster_by_island <- table(env.sub$ward,env$island)
#write as tables if needed
#write.csv(cluster_by_island, "cluster_by_island.csv", row.names=TRUE)
#write.csv(ward14_arch, "ward14_arch.csv", row.names=TRUE)

#distribution most frequent species across groups, used in Supplemental Information Table S4
table(env.sub$ward,com.data$Bursera.simaruba)
table(env.sub$ward,com.data$Coccoloba.diversifolia)
table(env.sub$ward,com.data$Metopium.toxiferum)
table(env.sub$ward,com.data$Bourreria.baccata)
table(env.sub$ward,com.data$Krugiodendron.ferreum)
table(env.sub$ward,com.data$Guapira.discolor)
table(env.sub$ward,com.data$Amyris.elemifera)
table(env.sub$ward,com.data$Eugenia.foetida)
table(env.sub$ward,com.data$Eugenia.axillaris)
table(env.sub$ward,com.data$Sideroxylon.salicifolium)
table(env.sub$ward,com.data$Vachellia.choriophylla)
table(env.sub$ward,com.data$Exothea.paniculata)
table(env.sub$ward,com.data$Gymnanthes.lucida)
table(env.sub$ward,com.data$Sideroxylon.foetidissimum)
table(env.sub$ward,com.data$Leucaena.leucocephala)
table(env.sub$ward,com.data$Exostema.caribaeum)
table(env.sub$ward,com.data$Reynosia.septentrionalis)
table(env.sub$ward,com.data$Tabebuia.heterophylla)
table(env.sub$ward,com.data$Swietenia.mahagoni)
table(env.sub$ward,com.data$Lysiloma.latisiliquum)

#Non-parametric analysis of variance examines differences among groups of sites defined by clustering,
#as well as groups defined by archipelagos
#MUST DEFINE GROUP NUMBER AS A FACTOR
env.sub$ward <- factor(env.sub$ward)
adonis(com.data ~ ward, data = env.sub, permutations=999)     #PERMANOVA examines differences in group means
adonis(com.data ~ archipelago, data = env, permutations=999)     #PERMANOVA examines differences in group means

#Also examine DCA (indirect ordination) for comparison of species turnoven on indirect gradient with genus level analysis
tdf.dca<-decorana(com.data)
tdf.dca

###FIGURE S3
#PCA for climate data 
env.pca <-pca(env.sub[,3:10], cor = TRUE)

install.packages(c("ggfortify"))
#Climate PCA figure with environmental vectors
library(ggfortify)
summary(prcomp(env.sub[,3:10],scale=TRUE))
autoplot(prcomp(env.sub[,3:10],scale=TRUE), data = env, colour = 'archipelago',
         loadings=TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3)

#non-parametric ANOVA of climate variables versus archipelago
adonis(env.sub[,3:10] ~ archipelago, data=env, permutations=999)     #examines differences in group mean

#partitioning sources of variation in beta diversity (turnover) using package betapart()
install.packages("betapart")
library(betapart)

#make sure the correct com.data is in the R Environment
tdf_pa.core <- betapart.core(com.data)
#yields 3 multi-site dissimilarities reported in Results

tdf_pa.multi.j <- beta.multi(tdf_pa.core,index.family="jac")
tdf_pa.multi.j

#yield 3 dissimilarity matrices pairwise disimilaries for all sites
tdf.dist<-beta.pair(tdf_pa.core, index.family="jac")

#prepare input for gdm()
dist.turn<-tdf.dist$beta.jtu

# sampling across equal sites
tdf.samp <- beta.sample(tdf_pa.core, sites=10, samples=100)
# plotting the distributions of components dist.s <- ceram.s.samp$sampled.values
distrib <- tdf.samp$sampled.values 
plot(density(distrib$beta.SOR), xlim=c(0,1), ylim=c(0, 25), xlab='Beta diversity', main='', lwd=3) 
lines(density(distrib$beta.SNE), lty=3, lwd=2)
lines(density(distrib$beta.SIM), lty=2, lwd=2) 

#GDM
install.packages("gdm")
library(gdm)

attach(env)
Climate <- data.frame(TDF_raw.plot_id,lat,long, MAT,CV_T,Twarm,Tcold,MAP,Pdry) 
Distance<-data.frame(TDF_raw.plot_id,lat,long)

##...Climate + Geographical distance dataframe / "B_S_sim_100" is a dissimilarity matrix / 
## Note that bioFormat have to be "3" here
#Jaccard turnover from betapart
DistMat<-as.matrix(dist.turn)
TDF_raw.plot_id<-rownames(DistMat)
bioData=data.frame(TDF_raw.plot_id, DistMat)

data_GDM_S_sim_Clim<-formatsitepair(bioData,bioFormat=3,siteColumn="TDF_raw.plot_id",XColumn="long",YColumn="lat",predData=Climate)

data_GDM_S_sim_Clim<-formatsitepair(bioData=data.frame(TDF_raw.plot_id,as.matrix(dist)),bioFormat=3,siteColumn="TDF_raw.plot_id",XColumn="long",YColumn="lat",predData=Climate)

##...Geographical distance dataframe
data_GDM_S_sim_Dist<-formatsitepair(bioData,bioFormat=3,siteColumn="TDF_raw.plot_id",XColumn="long",YColumn="lat",predData=Distance)
#summary(data_GDM_S_sim_Dist)
data_GDM_S_sim_Dist<-formatsitepair(bioData=data.frame(TDF_raw.plot_id,as.matrix(dist)),bioFormat=3,siteColumn="TDF_raw.plot_id",XColumn="long",YColumn="lat",predData=Distance)

detach(env)

##Fit GDM with only Climate (geo=FALSE), only Distance and both Climate and Distance variables (geo=TRUE)
tab_GDM_S_sim_all_Clim<-gdm(data_GDM_S_sim_Clim,geo=FALSE)
tab_GDM_S_sim_all_Dist<-gdm(data_GDM_S_sim_Dist,geo=TRUE)
tab_GDM_S_sim_all_Both<-gdm(data_GDM_S_sim_Clim,geo=TRUE)

##...A function to generate "R2" (variance explained) from the GDM results from T. Ibanez
toR2a<-function(x){
  R2<-1-(x$gdmdeviance/x$nulldeviance)
  total_df<-x$sample-1
  residual_df<-x$sample-length(x$predictors)-1
  R2A<-1-(1-R2)*(total_df/residual_df)
  return(R2A)
}

##...Here the partitioning of the variance (following Legendre 2008)
A_B<-toR2a(tab_GDM_S_sim_all_Clim) ##...Variance explained by Climate (A) + the interaction between Climate and Distance (B).
B_C<-toR2a(tab_GDM_S_sim_all_Dist) ##...Variance explained by Distance (C) + the interaction between Distance and Climate (B).
A_B_C<-toR2a(tab_GDM_S_sim_all_Both) ##...Variance explained by Climate (A) + Distance (C) + the interaction between Distance and Climate (B).

B<-round(A_B+B_C-A_B_C,2) ##...Variance explained by the interaction between Climate and Distance (B)
A<-round(A_B-B,2) ##...Variance explained by Climate (A)
C<-round(B_C-B,2) ##...Variance explained by Distance (C)
D<-round(1-A_B_C,2) ##...Unexplained variance (residuals)
##These values used to create Figure 4.
C ##...Variance explained by Distance (C)
A ##...Variance explained by Climate (A)
B ##...Variance explained by the combined Climate and Distance (B)
D ##...Unexplained variance (residuals)

#Predictor relative importance used in Table 4
VAR_IMP_SPECIES_ALL<-numeric()
for (i in 1:length(tab_GDM_S_sim_all_Both$predictor)) VAR_IMP_SPECIES_ALL[i]<-sum(tab_GDM_S_sim_all_Both$coefficients[(1+(3*(i-1))):(i*3)])

EXP<-tab_GDM_S_sim_all_Both$explained/100

VAR_IMP_SPECIES_ALL<-(VAR_IMP_SPECIES_ALL*EXP)/sum(VAR_IMP_SPECIES_ALL)
names(VAR_IMP_SPECIES_ALL)<-tab_GDM_S_sim_all_Both$predictor
VAR_IMP_SPECIES_ALL
