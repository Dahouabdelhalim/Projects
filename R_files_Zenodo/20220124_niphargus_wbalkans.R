# R code for manuscript entitled "A hotspot of groundwater amphipod diversity on a crossroad of evolutionary radiations"
# 24th January 2022

#Set wd (command is commented)
##setwd(insert the path to "Analysis" folder)

# Dependables ------------------------------------------------------------------
library(rgdal)
library(raster)
library(fuzzySim)
library(phytools)
library(picante)
library(vegan)
library(iNEXT)

# Data import and transformation -----------------------------------------------
# . Lambert projection, W Balkans####
cs_lambert_proj4=c("+proj=lcc +lat_0=12 +lon_0=18 +lat_1=42 +lat_2=46 +x_0=600000 +y_0=-3000000 +datum=WGS84 +units=m +no_defs")

# . Map, W Balkans ####
dinMapLambert=readOGR(dsn="./data/maps", layer="JugaSkupaj_lamb")
projection(dinMapLambert) <- CRS(cs_lambert_proj4)

# . Species occurrence table ####
# As spatial data frame
speciesCoordinates <- read.csv("./data/species_occurence.csv", encoding="UTF-8-BOM", fileEncoding="UTF-8-BOM")
coordinates(speciesCoordinates) <- ~x + y
proj4string(speciesCoordinates) <- CRS(cs_lambert_proj4)

# . Raster 20 km x 20 km ####
r20=raster(extent(matrix(c(220000, 460000, 820000,  1020000), nrow=2)),res=20000,crs = cs_lambert_proj4)
r20[] <- 1:ncell(r20)

# .. Species incidence table ####
# Translate SDF to incidence table in given raster
ext<-raster::extract(r20, speciesCoordinates) #assign raster id to each species occurrence
ext<-matrix(ext, ncol=1) # translate to matrix
colnames(ext)<-c("rasterID") #name cols with IDs of raster cells
voucher<-speciesCoordinates$voucher
ext<-cbind(ext,voucher)

speciesCoordinates<-cbind(speciesCoordinates,ext) #merge SDF and raster ID
speciesCoordinates<-speciesCoordinates[,-6] #remove redundant column

species10 <- as.data.frame(speciesCoordinates) # make new df, for occurence
species10=species10[,c(4,5)] # keep only species and rasterID
species10=splist2presabs(species10, 2, 1, keep.n = FALSE) # translate to incidence
rownames(species10)=species10$rasterID
species10=species10[,-1]

# . Tree ####
beastTree <- read.nexus("./data/beast.tree")
beastTree <- drop.tip(beastTree, "Nl_nolli") # drop outgroup


# Species Richness -------------------------------------------------------------
AlphaSR <- rasterize(speciesCoordinates, r20, 'species',
                     function(x, ...) length(unique(na.omit(x)))) # rasterize SR df
AlphaSR[is.na(AlphaSR[])] <- 0 

#plot species richness
cutsSR=c(2, 5, 8, 11, 16)
palleteSR <- c("white","#FFFFBE","#F4B800","#F57B00", "#F53D00", "#720000")
#pdf(file="SR.pdf") # un-comment for export to pdf
plot(AlphaSR, col = palleteSR, breaks = c(0,0.1,cutsSR), axes=F, main="Species richness")
plot(dinMapLambert, add=TRUE,col='darkgray', lwd=0.7)
plot(rasterToPolygons(AlphaSR, na.rm=TRUE, fun=function(x){x>0}), add=TRUE, border='black', lwd=0.7) 
#dev.off() # un-comment for export to pdf

# Phylogenetic Diversity (PD) -------------------------------------------------------

# calculate PD
PDbeastTree <- pd(species10,beastTree, include.root = F)
PDbeastTree[is.na(PDbeastTree)] <- 0 # assign null value to single species cells

# assign Raster ID
PDbeastTree$rasterID=rownames(PDbeastTree)
# make raster
rasPD <- raster(extent(matrix(c(220000, 460000, 820000,  1020000), nrow=2)),res=20000, 
                crs = cs_lambert_proj4)
rasPD[] <- -0.1

# fill raster with PD
for (i in 1:nrow(PDbeastTree)) {
  x=as.double(PDbeastTree[i,]$rasterID)
  rasPD[x]=as.double(PDbeastTree[i,]$PD)
}
writeRaster(rasPD,filename = "alphaPD")

# plot PD
cutsPD=seq(50,200,by=50)
palletePD <- c("white","lightgray","#FFFFBE","#F4B800","#F57B00", "#F53D00", "#720000")
#pdf(file="PD.pdf") # un-comment for export to pdf
plot(rasPD, col = palletePD, breaks = c(-0.1,-0.01,0.1,cutsPD,253), axes=F, main="Phylogenetic Diversity")
plot(dinMapLambert, add=TRUE,col='darkgray', lwd=0.7)
plot(rasterToPolygons(rasPD, na.rm=TRUE, fun=function(x){x>=0}), add=TRUE, border='black', lwd=0.7) 
scalebar(200000, divs=4, type="bar")
#dev.off() # un-comment for export to pdf

# Standardised Phylogenetic Diversity: sesPD ------------------------------------------

# Randomize community data matrix with the independent swap algorithm (Gotelli 2000)
# maintaining species occurrence frequency and sample species richness
sesPDbeastTreeIS <- ses.pd(species10,beastTree,runs = 9999, iterations = 10000,
                           null.model = "independentswap", include.root=FALSE)
sesPDbeastTreeIS <- sesPDbeastTreeIS[complete.cases(sesPDbeastTreeIS), ] # remove single species cells
sesPDbeastTreeIS$rasterID=rownames(sesPDbeastTreeIS)

# make raster
rasPDses <- raster(extent(matrix(c(220000, 460000, 820000,  1020000), nrow=2)),res=20000, 
                   crs = cs_lambert_proj4)
rasPDses[] <- -3.6

# fill raster with sesPD
for (i in 1:nrow(sesPDbeastTreeIS)) {
  x=as.double(sesPDbeastTreeIS[i,]$rasterID)
  rasPDses[x]=as.double(sesPDbeastTreeIS[i,]$pd.obs.z)
}

# make raster
rasPDp <- raster(extent(matrix(c(220000, 460000, 820000,  1020000), nrow=2)),res=20000, 
                 crs = cs_lambert_proj4)
rasPDp[] <- NA
# fill raster with p-value of sesPD
for (i in 1:nrow(sesPDbeastTreeIS)) {
  x=as.double(sesPDbeastTreeIS[i,]$rasterID)
  rasPDp[x]=as.double(sesPDbeastTreeIS[i,]$pd.obs.p)
}
writeRaster(rasPDses,filename = "alphasesPD")

# plot sesPD
cutsPDses=seq(-3.5,2.5,by=1)
palletePDses <- c("white","#0C4C96","#3985D1","#93CCED","#FFFFBE","#F4B800","#F53D00")
#pdf(file="PDses.pdf") # un-comment for export to pdf
plot(rasPDses, col = palletePDses, breaks = c(-3.6,cutsPDses), axes=F, main="standardised PD all Dinaric species")
plot(dinMapLambert, add=TRUE)
plot(rasterToPolygons(rasPDp, na.rm=TRUE, fun=function(x){x<=0.025 | x>=0.975}), lwd=2, add=TRUE) 
#dev.off() # un-comment for export to pdf

# Standardised Phylogenetic Diversity: Residuals from linear regression ------------------------------------------
# Filter only cells with more than one species
PDbeastTree_filtered <- PDbeastTree[PDbeastTree$SR > 1,]
#make linear regression
lm_pdsr <- lm(PD~SR, data = PDbeastTree_filtered)
#print summary of linear regression (supplementary material )
summary(lm_pdsr)
#save residuals from lm
lm_pdsr_resid <- lm_pdsr$residuals

#make raster
rasPDSRlm <- raster(extent(matrix(c(220000, 460000, 820000,  1020000), nrow=2)),res=20000, 
                    crs = cs_lambert_proj4)
rasPDSRlm[] <- -57
# fill raster with residuals
for (i in 1:length(lm_pdsr_resid)) {
  x=as.double(names(lm_pdsr_resid[i]))
  rasPDSRlm[x]=as.double(lm_pdsr_resid[i])
}

# plot residuals on a map, use same colour scale as for sesPD
cutsPDSRres=c(-57,-56,-40,-20,0,20,40,46)
palletePDSRres <- c("white","#0C4C96","#3985D1","#93CCED","#FFFFBE","#F4B800","#F53D00")
#pdf(file="PDresiduals.pdf") # un-comment for export to pdf
plot(rasPDSRlm, col = palletePDSRres, breaks = c(cutsPDSRres), axes=F, main="Residuals PD~SR")
plot(dinMapLambert, add=TRUE,col='darkgray', lwd=0.7)
plot(rasterToPolygons(rasPDSRlm, na.rm=TRUE, fun=function(x){x>-57}), add=TRUE, border='black', lwd=0.7) 
#dev.off() # un-comment for export to pdf

# Species Richness versus Standardised Phylogenetic Diversity ------------------
PDSRratio=as.data.frame(sesPDbeastTreeIS)
# Make categories: species richness < or >=7, sesPD low (>-0.5), around zero, or high (>0.5)
PDSRratio$conservation[PDSRratio$ntaxa<7 | (PDSRratio$pd.obs.z<0.5&PDSRratio$pd.obs.z>-0.5)] <- 0
PDSRratio$conservation[PDSRratio$ntaxa>=7 & (PDSRratio$pd.obs.z<0.5&PDSRratio$pd.obs.z>0)]<- 3
PDSRratio$conservation[PDSRratio$ntaxa>=7 & PDSRratio$pd.obs.z>0.5] <- 1
PDSRratio$conservation[PDSRratio$ntaxa>=7 & PDSRratio$pd.obs.z< -0.5] <- 2

#make raster
rasPDSR <- raster(extent(matrix(c(220000, 460000, 820000,  1020000), nrow=2)),res=20000, 
                        crs = cs_lambert_proj4)
rasPDSR[] <- NA

# fill raster with values of SR vs sesPD
for (i in 1:nrow(PDSRratio)) {
  x=as.double(PDSRratio[i,]$rasterID)
  rasPDSR[x]=as.double(PDSRratio[i,]$conservation)
}

# plot SR vs sesPD
cutsPDSR=c(-0.1,0.1,1.1,2.1,3.1)
palletePDSR <- c("#FFFFBE","#F53D00","#2166AC","green")
#pdf(file="PDSR.pdf") # un-comment for export to pdf
plot(rasPDSR, legend = FALSE, col = palletePDSR, breaks = cutsPDSR, axes=F, main="SR versus sesPD")
plot(dinMapLambert, add=TRUE, col='darkgray', lwd=0.7)
plot(rasterToPolygons(rasPDSR, na.rm=TRUE, fun=function(x){x>=0}), add=TRUE, border='black', lwd=0.7) 
legend("bottomleft", 
       legend = c("SR<7 |(-0.5<sesPD<0.5)", "SR>=7 & sesPD>0.5", "SR>=7 & sesPD<-0.5"), 
       fill = palletePDSR)
#dev.off() # un-comment for export to pdf

# Plot dotplot of SR versus sesPD, coloured by discrete groups
PDSRratio$col[PDSRratio$ntaxa<7 | (PDSRratio$pd.obs.z<0.5&PDSRratio$pd.obs.z>-0.5)]<- "#FFFFBE"
PDSRratio$col[PDSRratio$ntaxa>=7 & (PDSRratio$pd.obs.z<0.5&PDSRratio$pd.obs.z>-0.5)]<- "green"
PDSRratio$col[PDSRratio$ntaxa>=7 & PDSRratio$pd.obs.z>0.5]<- "#F53D00"
PDSRratio$col[PDSRratio$ntaxa>=7 & PDSRratio$pd.obs.z< -0.5]<- "#2166AC"
#pdf(file="PDSRplot.pdf") # un-comment for export to pdf
plot(PDSRratio$ntaxa,PDSRratio$pd.obs.z,pch=21, bg=PDSRratio$col)
#dev.off() # un-comment for export to pdf

# Appendix S5, Fig.S5.1a: Species Accumulation Curve  ------------------
species10NW <- subset(species10,as.numeric(rownames(species10))%in%c(1:311, 336:340))
species10SE <- subset(species10,!(as.numeric(rownames(species10))%in%c(1:311, 336:340)))
specaccum_NW <- specaccum(species10NW, method = "random")
specaccum_SE <- specaccum(species10SE, method = "random")

lightblue50 <- rgb(173,216,230, max = 255, alpha = 125)
pink50 <- rgb(255,192,203, max = 255, alpha = 125)
plot(specaccum_NW, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col=lightblue50, ylab="SAC (random order)", xlab="Cells")
plot(specaccum_SE, add=T, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col=pink50)
legend('bottomright', c('northwest', 'southeast'), col=c("blue", "red"), lty=1, 
       bty='n', inset=0.025)

# Appendix S5, Fig.S5.1b: Rarefaction and Extrapolation Sampling Curve  ------------------
dataNW <- as.data.frame(t(species10NW))
dataSE <- as.data.frame(t(species10SE))
#inext needs a named list
data_list <- list(dataNW,dataSE)
names(data_list) <- c("northwest", "southeast")
#perform iNEXT
inext_curve = iNEXT(data_list, datatype = "incidence_raw", endpoint = 600)
ggiNEXT(inext_curve, facet.var = "none") + labs(x = "nr. of cells")+theme_bw()

# --------------------end of file -----------------------------------