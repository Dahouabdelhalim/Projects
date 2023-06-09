

library(GenABEL)
library(plyr)
library(dplyr)


# load pheno file
pheno_month <- read.table("pheno_month_final.txt", header = T)
head(pheno_month)

pheno_month <- "pheno_month_final.txt"

# Estimate phenotypic variance
var(pheno_month$new_wing_mean_month_as_factor, na.rm = T)

#use the ped file and map file to make the .raw file which is an internal genotypic data for genabel ####
convert.snp.ped(ped="../Original\\ Files/pdosnp_6400_n1898.ped", map="../Original\\ Files/pdosnp_6400_n1898.map", out="pass_data.raw", mapHasHeaderLine=FALSE)

geno <- "pass_data.raw" #assign the .raw file to genotypes
gwaa <- load.gwaa.data(pheno_month,geno,force=TRUE,makemap=F)

nids(gwaa)


# quality control
qc1 <- check.marker(gwaa, call=0.95, perid.call=0.95, maf=0.01, p.lev =0, ibs.threshold=0.9)
summary(qc1)

gwaa0.1<-gwaa[qc1$idok, qc1$snpok]
# gwaa0.1.1<- Xfix(gwaa0.1)

nsnps(gwaa0.1)
nids(gwaa0.1)

#### PCA  ####

# Create ibs matrix
data.qc1.gkin <- ibs(gwaa0.1, weight="freq") #Note: this line may take a long time, depending on # of inds
# where gwaa0.1.1 is a genABEL gwaa.data.class object that has gone through preliminary quality control for call rates, ibs, and maf

data.qc1.dist <- as.dist(0.5-data.qc1.gkin)

#Create PC's, where k= # of PC's
data.qc1.mds <- cmdscale(data.qc1.dist, k=4) 

#Generate plot of within group sums of squares in order to decide how many clusters are proper
wss <- NA
#can let the for loop run for longer and get the same answer...
for (i in 1:15) wss[i] <- sum(kmeans(data.qc1.mds, centers=i, nstart=nids(gwaa0.1))$withinss)
plot(1:15, wss, type="b", xlab="Number of clusters", ylab="WSS")
# Curve flattens out at ~6, 6 clusters should be good

km <- kmeans(data.qc1.mds, centers=6, nstart=nids(gwaa0.1)) #Assigns each individual to a cluster

# Plot subpopulations #

# Create vectors for plots
cl1 <- which(km$cluster==1)
cl1num <- as.numeric(cl1)
cl1x <- data.qc1.mds[cl1num,1] #PC1 for individs in cluster 1
cl1y <- data.qc1.mds[cl1num,2] #PC2 for individs in cluster 1
cl1z <- data.qc1.mds[cl1num,3]
cl1w <- data.qc1.mds[cl1num,4]

cl2 <- which(km$cluster==2)
cl2num <- as.numeric(cl2)
cl2x <- data.qc1.mds[cl2num,1] #PC1 for individs in cluster 2
cl2y <- data.qc1.mds[cl2num,2]
cl2z <- data.qc1.mds[cl2num,3]
cl2w <- data.qc1.mds[cl2num,4]

cl3 <- which(km$cluster==3)
cl3num <- as.numeric(cl3)
cl3x <- data.qc1.mds[cl3num,1]
cl3y <- data.qc1.mds[cl3num,2]
cl3z <- data.qc1.mds[cl3num,3]
cl3w <- data.qc1.mds[cl3num,4]

cl4 <- which(km$cluster==4)
cl4num <- as.numeric(cl4)
cl4x <- data.qc1.mds[cl4num,1]
cl4y <- data.qc1.mds[cl4num,2]
cl4z <- data.qc1.mds[cl4num,3]
cl4w <- data.qc1.mds[cl4num,4]

cl5 <- which(km$cluster==5)
cl5num <- as.numeric(cl5)
cl5x <- data.qc1.mds[cl5num,1]
cl5y <- data.qc1.mds[cl5num,2]
cl5z <- data.qc1.mds[cl5num,3]
cl5w <- data.qc1.mds[cl5num,4]

cl6 <- which(km$cluster==6)
cl6num <- as.numeric(cl6)
cl6x <- data.qc1.mds[cl6num,1]
cl6y <- data.qc1.mds[cl6num,2]
cl6z <- data.qc1.mds[cl6num,3]
cl6w <- data.qc1.mds[cl6num,4]


par(mfrow=c(2,3))

#PCA 1 and 2
plot(data.qc1.mds, type="n", xlab="MDS1", ylab="MDS2", main="1 vs 2 (4 PC's, 6 clusters)")
points(cl1x,cl1y, pch=16, cex=0.7, col="dodgerblue")
points(cl2x,cl2y,pch=19,cex=.5,col="red")
points(cl3x,cl3y,pch=19,cex=.5,col="green")
points(cl4x,cl4y,pch=19,cex=.5,col="black")
points(cl5x,cl5y,pch=19,cex=.5,col="orange")
points(cl6x,cl6y,pch=19,cex=.5,col="grey")

#PCA 1 and 3
plot(data.qc1.mds, type="n", xlab="MDS1", ylab="MDS3", main="1 vs 3 (4 PC's, 6 clusters)")
points(cl1x,cl1z, pch=16, cex=0.7, col="dodgerblue")
points(cl2x,cl2z,pch=19,cex=.5,col="red")
points(cl3x,cl3z,pch=19,cex=.5,col="green")
points(cl4x,cl4z,pch=19,cex=.5,col="black")
points(cl5x,cl5z,pch=19,cex=.5,col="orange")
points(cl6x,cl6z,pch=19,cex=.5,col="grey")

#PCA 1 and 4
plot(data.qc1.mds, type="n", xlab="MDS1", ylab="MDS4", main="1 vs 4 (4 PC's, 6 clusters)")
points(cl1x,cl1w, pch=16, cex=0.7, col="dodgerblue")
points(cl2x,cl2w,pch=19,cex=.5,col="red")
points(cl3x,cl3w,pch=19,cex=.5,col="green")
points(cl4x,cl4w,pch=19,cex=.5,col="black")
points(cl5x,cl5w,pch=19,cex=.5,col="orange")
points(cl6x,cl6w,pch=19,cex=.5,col="grey")

#PCA 2 and 3
plot(data.qc1.mds, type="n", xlab="MDS2", ylab="MDS3", main="2 vs 3 (4 PC's, 6 clusters)")
points(cl1y,cl1z, pch=16, cex=0.7, col="dodgerblue")
points(cl2y,cl2z,pch=19,cex=.5,col="red")
points(cl3y,cl3z,pch=19,cex=.5,col="green")
points(cl4y,cl4z,pch=19,cex=.5,col="black")
points(cl5y,cl5z,pch=19,cex=.5,col="orange")
points(cl6y,cl6z,pch=19,cex=.5,col="grey")

#PCA 2 and 4
plot(data.qc1.mds, type="n", xlab="MDS2", ylab="MDS4", main="2 vs 4 (4 PC's, 6 clusters)")
points(cl1y,cl1w, pch=16, cex=0.7, col="dodgerblue")
points(cl2y,cl2w,pch=19,cex=.5,col="red")
points(cl3y,cl3w,pch=19,cex=.5,col="green")
points(cl4y,cl4w,pch=19,cex=.5,col="black")
points(cl5y,cl5w,pch=19,cex=.5,col="orange")
points(cl6y,cl6w,pch=19,cex=.5,col="grey")

#PCA 3 and 4
plot(data.qc1.mds, type="n", xlab="MDS3", ylab="MDS4", main="3 vs 4 (4 PC's, 6 clusters)")
points(cl1z,cl1w, pch=16, cex=0.7, col="dodgerblue")
points(cl2z,cl2w,pch=19,cex=.5,col="red")
points(cl3z,cl3w,pch=19,cex=.5,col="green")
points(cl4z,cl4w,pch=19,cex=.5,col="black")
points(cl5z,cl5w,pch=19,cex=.5,col="orange")
points(cl6z,cl6w,pch=19,cex=.5,col="grey")

par(mfrow=(c(1,1)))

cl1names <- names(which(km$cluster==1)) 
cl1_phen <- subset(pheno_month, pheno_month$ringnr %in% cl1names)
table(cl1_phen$hatchisland)

cl2names <- names(which(km$cluster==2)) 
cl2_phen <- subset(pheno_month, pheno_month$ringnr %in% cl2names)
table(cl2_phen$hatchisland)

cl3names <- names(which(km$cluster==3)) 
cl3_phen <- subset(pheno_month, pheno_month$ringnr %in% cl3names)
table(cl3_phen$hatchisland)

cl4names <- names(which(km$cluster==4)) 
cl4_phen <- subset(pheno_month, pheno_month$ringnr %in% cl4names)
table(cl4_phen$hatchisland)

cl5names <- names(which(km$cluster==5)) 
cl5_phen <- subset(pheno_month, pheno_month$ringnr %in% cl5names)
table(cl5_phen$hatchisland)

cl6names <- names(which(km$cluster==6)) 
cl6_phen <- subset(pheno_month, pheno_month$ringnr %in% cl6names)
table(cl6_phen$hatchisland)

# Note that cluster names are arbitrarily assigned and individs in cluster "1" this time
# may be in a different cluster if you were to rerun the script

# Quality Control for HWE ####
# assess HWE for each cluster separately, then only remove SNP's out of HWE in all cluster

qc_cl1<- check.marker(gwaa0.1, call=0.95, perid.call=0.95, maf=0, p.lev =0.001, hweidsubset=cl1names, ibs.threshold=1)
qc_cl2<- check.marker(gwaa0.1, call=0.95, perid.call=0.95, maf=0, p.lev =0.001, hweidsubset=cl2names, ibs.threshold=1)
qc_cl3<- check.marker(gwaa0.1, call=0.95, perid.call=0.95, maf=0, p.lev =0.001, hweidsubset=cl3names, ibs.threshold=1)
qc_cl4<- check.marker(gwaa0.1, call=0.95, perid.call=0.95, maf=0, p.lev =0.001, hweidsubset=cl4names, ibs.threshold=1)
qc_cl5<- check.marker(gwaa0.1, call=0.95, perid.call=0.95, maf=0, p.lev =0.001, hweidsubset=cl5names, ibs.threshold=1)
qc_cl6<- check.marker(gwaa0.1, call=0.95, perid.call=0.95, maf=0, p.lev =0.001, hweidsubset=cl6names, ibs.threshold=1)

# Subset of snp's where only those that failed qc in all clusters are removed
gwaa0.2 <- gwaa0.1[,gwaa0.1@gtdata@snpnames %in% qc_cl1$snpok | 
                     gwaa0.1@gtdata@snpnames %in% qc_cl2$snpok | 
                     gwaa0.1@gtdata@snpnames %in% qc_cl3$snpok | 
                     gwaa0.1@gtdata@snpnames %in% qc_cl4$snpok |
                     gwaa0.1@gtdata@snpnames %in% qc_cl5$snpok |
                     gwaa0.1@gtdata@snpnames %in% qc_cl6$snpok]

nsnps(gwaa0.2)
nids(gwaa0.2)

# data descriptives

descriptives.trait(gwaa0.2) #prints out nr records, mean and SD

descriptives.marker(gwaa0.2) #prints out descriptive stats for the markers

excluded_birds_gwaa0.2 <- setdiff(pheno_month$id,idnames(gwaa0.2)) #find excluded ind
excluded_birds_gwaa0.2

#make kinship matrix
data.agemonth.gkin <- ibs(gwaa0.2, weight="freq")

#examine distribution of kinship coefficients

hist(data.agemonth.gkin)




#### 	GWAS USING GRAMMAR    ####
# GAMMA APPROACH FROM SVICEHVA ET AL. 2012 NATURE GENETICS


#### TARSUS  ####

tarsus <- polygenic(new_tars_mean_month_as_factor ~ hatchyear + sex, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

tarsus$esth2
tarsus$h2an

var(pheno_month$new_tars_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_tarsus <- grammar(tarsus, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_tarsus, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_tarsus)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_tarsus),chromosome=chromosome(grgamma_tarsus), position=map(grgamma_tarsus),refallele=refallele(grgamma_tarsus), codeallele=effallele(grgamma_tarsus), refallelefreq=refallfreq, beta=grgamma_tarsus[,"effB"], se_beta=grgamma_tarsus[,"se_effB"], P_val=grgamma_tarsus[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_tars_mean.csv",row.names=FALSE)


# Figs

png("Manhattan tarsus mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_tarsus, df="Pc1df", sub="", main="Tarsus", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot tarsus mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_tarsus)
estlambda(grgamma_tarsus[,"P1df"], plot=T, main = "Q-Q Tarsus")
dev.off()




#### WING ####
wing <- polygenic(new_wing_mean_month_as_factor ~ sex, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

wing$esth2
wing$h2an

var(pheno_month$new_wing_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_wing <- grammar(wing, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_wing, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_wing)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_wing),chromosome=chromosome(grgamma_wing), position=map(grgamma_wing),refallele=refallele(grgamma_wing), codeallele=effallele(grgamma_wing), refallelefreq=refallfreq, beta=grgamma_wing[,"effB"], se_beta=grgamma_wing[,"se_effB"], P_val=grgamma_wing[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_wing_mean.csv",row.names=FALSE)

# Figs

png("Manhattan wing mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_wing, df="Pc1df", sub="", main="Wing", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot wing mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_wing)
estlambda(grgamma_wing[,"P1df"], plot=T, main = "Q-Q Wing")
dev.off()


#### MASS ####
mass <- polygenic(new_mass_mean_month_as_factor ~ sex, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

mass$esth2
mass$h2an

var(pheno_month$new_mass_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_mass <- grammar(mass, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_mass, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_mass)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_mass),chromosome=chromosome(grgamma_mass), position=map(grgamma_mass),refallele=refallele(grgamma_mass), codeallele=effallele(grgamma_mass), refallelefreq=refallfreq, beta=grgamma_mass[,"effB"], se_beta=grgamma_mass[,"se_effB"], P_val=grgamma_mass[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_mass_mean.csv",row.names=FALSE)

# Figs

png("Manhattan mass mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_mass, df="Pc1df", sub="", main="Mass", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot mass mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_mass)
estlambda(grgamma_mass[,"P1df"], plot=T, main = "Q-Q Mass")
dev.off()



# BILL DEPTH ####
billD <- polygenic(new_bill_d_mean_month_as_factor ~ sex+ hatchisland, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

billD$esth2
billD$h2an

var(pheno_month$new_bill_d_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_billD <- grammar(billD, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_billD, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_billD)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_billD),chromosome=chromosome(grgamma_billD), position=map(grgamma_billD),refallele=refallele(grgamma_billD), codeallele=effallele(grgamma_billD), refallelefreq=refallfreq, beta=grgamma_billD[,"effB"], se_beta=grgamma_billD[,"se_effB"], P_val=grgamma_billD[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_billD_mean.csv",row.names=FALSE)


# Figs

png("Manhattan billD mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_billD, df="Pc1df", sub="", main="Bill Depth", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot billD mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_billD)
estlambda(grgamma_billD[,"P1df"], plot=T, main = "Q-Q Bill Depth")
dev.off()



# BILL LENGTH ####
billL <- polygenic(new_bill_l_mean_month_as_factor ~ hatchisland+ sex, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

billL$esth2
billL$h2an

var(pheno_month$new_bill_l_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_billL <- grammar(billL, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_billL, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_billL)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_billL),chromosome=chromosome(grgamma_billL), position=map(grgamma_billL),refallele=refallele(grgamma_billL), codeallele=effallele(grgamma_billL), refallelefreq=refallfreq, beta=grgamma_billL[,"effB"], se_beta=grgamma_billL[,"se_effB"], P_val=grgamma_billL[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_billL_mean.csv",row.names=FALSE)


# Figs

png("Manhattan billL mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_billL, df="Pc1df", sub="", main="Bill Length", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot billL mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_billL)
estlambda(grgamma_billL[,"P1df"], plot=T, main = "Q-Q Bill Length")
dev.off()



# TOTAL BADGE ####
totbadge <- polygenic(new_TotalBadge_mean_month_as_factor ~ hatchyear, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

totbadge$esth2

var(pheno_month$new_TotalBadge_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_totbadge <- grammar(totbadge, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_totbadge, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_totbadge)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_totbadge),chromosome=chromosome(grgamma_totbadge), position=map(grgamma_totbadge),refallele=refallele(grgamma_totbadge), codeallele=effallele(grgamma_totbadge), refallelefreq=refallfreq, beta=grgamma_totbadge[,"effB"], se_beta=grgamma_totbadge[,"se_effB"], P_val=grgamma_totbadge[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_totbadge_mean.csv",row.names=FALSE)


# Figs

png("Manhattan totbadge mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_totbadge, df="Pc1df", sub="", main="Total Badge", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot totbadge mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_totbadge)
estlambda(grgamma_totbadge[,"P1df"], plot=T, main = "Q-Q Total Badge")
dev.off()




# VIS BADGE ####
visbadge <- polygenic(new_VisibleBadge_mean_month_as_factor ~ hatchisland, kin=data.agemonth.gkin, gwaa0.2, trait.type= "gaussian", quiet=TRUE)

visbadge$esth2

var(pheno_month$new_VisibleBadge_mean_month_as_factor, na.rm = T)

#Fit GWAS using Gramma gamma approach
grgamma_visbadge <- grammar(visbadge, data=gwaa0.2, method="gamma")

#Pull out top 20 SNPs
descriptives.scan(grgamma_visbadge, top=20)[,c(1,2,4,5,6,7,8,9,15)]

gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(grgamma_visbadge)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grgamma_visbadge),chromosome=chromosome(grgamma_visbadge), position=map(grgamma_visbadge),refallele=refallele(grgamma_visbadge), codeallele=effallele(grgamma_visbadge), refallelefreq=refallfreq, beta=grgamma_visbadge[,"effB"], se_beta=grgamma_visbadge[,"se_effB"], P_val=grgamma_visbadge[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication

write.csv(tabOutSorted,"TopSNPs GWAS_visbadge_mean.csv",row.names=FALSE)


# Figs

png("Manhattan visbadge mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grgamma_visbadge, df="Pc1df", sub="", main="Visible Badge", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot visbadge mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grgamma_visbadge)
estlambda(grgamma_visbadge[,"P1df"], plot=T, main = "Q-Q Visible Badge")
dev.off()


# Export data for Plink ####

export.merlin(gwaa0.2, pedfile = "batch1_10k.ped", datafile = "batch1_10k.dat", mapfile = "batch1_10k.map", format = "plink", fixstrand = "no", extendedmap = TRUE, traits = 1, order = TRUE, stepids = 100)
