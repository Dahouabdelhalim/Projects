## Messerman, A. and M. Leal 2020
## Ambystoma juvenile salamanders (5 species)
## Respiratory Surface Area Water Loss (RSAWL) principle component analysis (PCA) 
## RSAWL figures
## Abiotic covariate PCA

data<-read.csv("WVP_RSAWL_Data_Means.csv", header=TRUE)
head(data)
data$round<-as.factor(data$round)
data$chamber<-as.factor(data$chamber)
data$channel<-as.factor(data$channel)
str(data)

data$log.rwl<-log10(data$mean.rwl)
data$log.svl<-log10(data$svl.test)
data$log.tl<-log10(data$tl.test)
data$log.wide<-log10(data$wide.test)
data$log.mass<-log10(data$premass)
str(data)

##Add column categorizing site latitude
data$pop [100:150]
data$lat<-ifelse(data$pop=="FLW", "2mid", "1north")
data$lat<-ifelse(data$pop=="Mingo", "3south", data$lat)
data$lat<-as.factor(data$lat)
data$lat [100:150]

##Column of continuous latitudinal value
data$lat.vec<-ifelse(data$pop=="FLW", 37.7415, data$pop)
data$lat.vec<-ifelse(data$pop=="Mingo", 36.9615, data$lat.vec)
data$lat.vec<-ifelse(data$pop=="Forum", 38.9222, data$lat.vec)
data$lat.vec<-ifelse(data$pop=="Danville", 38.8717, data$lat.vec)
data$lat.vec<-ifelse(data$pop=="DBCA", 38.7861, data$lat.vec)
data$lat.vec<-as.numeric(data$lat.vec)

##Standardize latitude covariate
for (i in 1:length(data$lat.vec)) {
  data$lat.vec.std[i] <- (data$lat.vec[i]-mean(data$lat.vec[]))/sd(data$lat.vec[])
}

##Column of spp*lat interaction
data$sppvlat<-interaction(data$spp,data$lat)
str(data)


##Remove first 2 rounds as flush, retain rounds 3-5 for analyses
data$round<-as.integer(data$round)
unique(data$round)
data3<-subset(data, data$round>=3)
str(data3)
data3$round<-as.factor(data3$round)


##Build dataframe of averaged RSAWL rounds 3-5
names(data3)
mean.df<-aggregate(data3, list(data3$ID), mean)
str(mean.df)
cleaned<-unique.data.frame(data3[,c(1:6,33)], incomparables = FALSE)
cleaned<-cleaned[order(cleaned$ID),]
str(cleaned)
rwl<-cbind(cleaned, mean.df$mean.rwl)
rwl<-cbind(rwl, mean.df$days)
rwl<-cbind(rwl, mean.df$premass)
rwl<-cbind(rwl, mean.df$log.mass)
rwl<-cbind(rwl, mean.df$svl.test)
rwl<-cbind(rwl, mean.df$tl.test)
rwl<-cbind(rwl, mean.df$wide.test)
rwl<-cbind(rwl, mean.df$log.svl)
rwl<-cbind(rwl, mean.df$log.tl)
rwl<-cbind(rwl, mean.df$log.wide)
rwl<-cbind(rwl, mean.df$lat.vec)
rwl<-cbind(rwl, mean.df$lat.vec.std)
names(rwl)
names(rwl)[8]<-paste("mean.rwl")
names(rwl)[9]<-paste("days")
names(rwl)[10]<-paste("premass")
names(rwl)[11]<-paste("log.mass")
names(rwl)[12]<-paste("svl")
names(rwl)[13]<-paste("tl")
names(rwl)[14]<-paste("widest")
names(rwl)[15]<-paste("log.svl")
names(rwl)[16]<-paste("log.tl")
names(rwl)[17]<-paste("log.wide")
names(rwl)[18]<-paste("lat.vec")
names(rwl)[19]<-paste("lat.vec.std")
str(rwl)
names(rwl)

#PCA Function and Summary
pca <- prcomp(~ tl + svl + log.mass + widest, data = rwl, center = TRUE, scale = TRUE)
pca #Scores

#extract eigenvalues and examine Kaiser's Rule
eigen <- pca$sdev^2

plot(eigen)
abline(h = 1, col = "blue")

summary(pca) #0.78 proportion variance explained on PC1

#Format data for PCA
pc1and2 <- as.data.frame(pca$x[,c(1,2)])

pc1<-c(pca$x[,1])#Coordinates per individual on PC1
rwl$pc1<-pc1
str(rwl)

#Visualize loadings of individual morphometrics 
biplot(pca, choices = c(1,2), xlabs=rep(".", nrow(rwl)), cex=1.5)

require(ggplot2)
library(ggfortify)
library(tidyverse)
library(cowplot)
library(dplyr)

species<-c(as.expression(bquote(italic(.("A. annulatum")))),as.expression(bquote(italic(.("A. maculatum")))), 
           as.expression(bquote(italic(.("A. opacum")))),as.expression(bquote(italic(.("A. talpoideum")))), 
           as.expression(bquote(italic(.("A. texanum")))))
spp<-rwl$spp
spp<-ifelse(spp=="AMAN", "purple", spp)
spp<-ifelse(spp=="AMMA", "blue", spp)
spp<-ifelse(spp=="AMOP", "red", spp)
spp<-ifelse(spp=="AMTA", "orange", spp)
spp<-ifelse(spp=="AMTE", "darkgreen", spp)
unique(spp)
tab<-table(rwl$spp)

a1 <- autoplot(pca,data=rwl,loadings=TRUE,loadings.colour=1,loadings.size=3)
rwl$spp <- as.factor(rwl$spp)
legends <- tibble(
  nums = c(unique(rwl$spp)),
  names = species,
  clrs = c("purple", "blue","red","orange", "darkgreen")
)

gg<-a1 + geom_point(aes(colour=spp)) +
  labs(x= "PC1",y="PC2") +
  theme(legend.position = c(0.15,0.9),
        legend.title=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        text = element_text(size=20))+
  annotate("text", x=0.1, y=0.17, label="DM", size=6)+
  annotate("text", x=0.137, y=0.041, label=expression(log[10] ~ (mass)), size=6)+
  annotate("text", x=0.11, y=-0.1, label="SVL/TL", size=6)

tiff("PCA_averaged_3-5_300dpi.tiff", width = 8.5, height = 8, units = 'in', res=300, compression = 'none')
par(mfrow=c(1,1), las=1)
gg + scale_color_manual(values = legends$clrs, labels = legends$names)
dev.off()

mod.pc<-lm(rwl$mean.rwl~rwl$pc1)
summary(mod.pc)
resid.pc<-resid(mod.pc)
rwl$resid.pc<-resid.pc
str(rwl)

##################################################
###Build figures for residual RSAWL rounds 3-5
#################################################

par(mfrow=c(1,1), las=1)
tiff("RSAWL_PC1-300dpi.tiff", width = 8, height = 6.5, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(2.2,1,0))
plot(rwl$pc1,rwl$mean.rwl, bty="l", pch=c(rwl$spp), col=c(spp),ylab=expression(italic(RSAWL) ~ (mg ~ cm^{-2} ~ h^{-1})), xlab="Coordinates on PC1", 
     cex=1.3, cex.lab=1.2, cex.axis=1.2)
legend(1.55, 10.35,bty="n", legend=c(as.expression(bquote(italic(.("A. annulatum")))),as.expression(bquote(italic(.("A. maculatum")))), 
                           as.expression(bquote(italic(.("A. opacum")))),as.expression(bquote(italic(.("A. talpoideum")))), 
                           as.expression(bquote(italic(.("A. texanum"))))), col=c("purple", "blue", "red", "orange", "dark green"), pch=c(1:5), cex=1.3)
abline(lm(rwl$mean.rwl~rwl$pc1), lwd=3)
dev.off()

summary(lm(rwl$mean.rwl~rwl$pc1))

aman3<-subset(rwl,rwl$spp=="AMAN")
amop3<-subset(rwl,rwl$spp=="AMOP")
amma3<-subset(rwl,rwl$spp=="AMMA")
amte3<-subset(rwl,rwl$spp=="AMTE")
amta3<-subset(rwl,rwl$spp=="AMTA")

x.aman<-mean(aman3$resid.pc)
x.amma<-mean(amma3$resid.pc)
x.amop<-mean(amop3$resid.pc)
x.amta<-mean(amta3$resid.pc)
x.amte<-mean(amte3$resid.pc)

resid.means<-c(x.aman, x.amma, x.amop, x.amta, x.amte)

#Average of differences in mean residual RSAWL by species
((x.aman-x.amma)+(x.aman-x.amop)+(x.aman-x.amta)+(x.aman-x.amte)+
  (x.amma-x.amop)+(x.amma-x.amta)+(x.amma-x.amte)+(x.amop-x.amta)+
  (x.amop-x.amte)+(x.amta-x.amte))/10 #0.16 (mg ~ cm^{-2} ~ h^{-1}) ~ by ~ PC1 ~ coordinates

sd.aman<-sd(aman3$resid.pc)
sd.amma<-sd(amma3$resid.pc)
sd.amop<-sd(amop3$resid.pc)
sd.amta<-sd(amta3$resid.pc)
sd.amte<-sd(amte3$resid.pc)

resid.sd<-c(sd.aman, sd.amma, sd.amop, sd.amta, sd.amte)

se.aman<-sd(aman3$residm) /  sqrt(length(aman3$residm)) 
se.amma<-sd(amma3$residm) /  sqrt(length(amma3$residm)) 
se.amop<-sd(amop3$residm) /  sqrt(length(amop3$residm)) 
se.amta<-sd(amta3$residm) /  sqrt(length(amta3$residm)) 
se.amte<-sd(amte3$residm) /  sqrt(length(amte3$residm)) 

upper.aman<-x.aman+se.aman
lower.aman<-x.aman-se.aman
sd.means<-c(sd.aman, sd.amma, sd.amop, sd.amta, sd.amte)
se.means<-c(se.aman, se.amma, se.amop, se.amta, se.amte)

par(mfrow=c(1,1), las=1)
tiff("RSAWL_SPP-300dpi.tiff", width = 8, height = 6, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(3,1,0))
plot(resid.means, cex=3, cex.axis=1.1, cex.lab=1.1, bg=c("purple", "blue", "red", "orange", "dark green"), pch=21, bty="l",
     ylab=expression(Mean ~ residual ~ italic(RSAWL) ~ (mg ~ cm^{-2} ~ h^{-1}) ~ by ~ PC1 ~ coordinates), xlab="Species", 
     ylim = c(-1.5,1), xaxt="n")
axis(1, 1:5, labels=species, cex.axis=1.1)
mtext(side=3,at = 1:5,c("A", "B", "BC", "D", "AC"), cex=1.1)
segments(1, x.aman+sd.aman,1,x.aman-sd.aman, col="purple", lwd=2)
segments(2, x.amma+sd.amma,2,x.amma-sd.amma, col="blue", lwd=2)
segments(3, x.amop+sd.amop,3,x.amop-sd.amop, col="red", lwd=2)
segments(4, x.amta+sd.amta,4,x.amta-sd.amta, col="orange", lwd=2)
segments(5, x.amte+sd.amte,5,x.amte-sd.amte, col="dark green", lwd=2)
#legend(3, 0.3,legend=species, col=c("purple", "blue", "red", "orange", "dark green"), pch=1, cex=3)
dev.off()


######################################################################################
##Abiotic variables by location
##Data obtained from https://www.ncdc.noaa.gov/cdo-web/datatools/normals
######################################################################################
covs<-read.csv("normals-metric.csv", header=TRUE)
str(covs)

#PCA Function and Summary
pca.abio <- prcomp(~ lat.vec + precip.cm + min.temp + x.temp + max.temp,
              data = covs, center = TRUE, scale = TRUE)
pca.abio #Scores

#extract eigenvalues and examine Kaiser's Rule
eigen.abio <- pca.abio$sdev^2

plot(eigen.abio)
abline(h = 1, col = "blue")

summary(pca.abio) #0.804 proportion variance explained on PC1

#Format data for PCA
pc1and2.abio <- as.data.frame(pca.abio$x[,c(1,2)])

pc1.abio<-c(pca.abio$x[,1])#Coordinates per location on PC1
covs$pc.abio<-pc1.abio

#Visualize loadings of locality environmental covariates
tiff("Covariate_PCA_300dpi_nolab.tiff", width = 6, height = 6, units = 'in', res=300, compression = 'none')
biplot(pca.abio, choices = c(1,2), xlabs=covs$co, cex=1, xlim=c(-1, 0.8), ylim=c(-1,0.8))
dev.off()

library(corrplot)
str(covs)
A<-cor(covs[3:7])

par(mfrow=c(1,1))
tiff("Corr_matrix_300dpi.tiff", width = 6, height = 6, units = 'in', res=300, compression = 'none')
corrplot(A, method="number")#min temp and mean temp <|0.7| correlated with latitude
dev.off()
#max temp has greatest correlation with all other variables (>=0.79)