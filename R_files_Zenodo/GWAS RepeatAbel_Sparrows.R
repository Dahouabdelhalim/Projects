
install.packages("GenABEL")
install.packages("hglm")
install.packages("gdata")
install.packages("RepeatABEL")

library(GenABEL)
library(hglm)
library(gdata)
library(RepeatABEL)


# DATA IMPORT ####

# reorganise repeat file - note: gwaa object is obtained using mean values!

pheno_rep <- read.table("pheno_rep_final.txt", header = T)
head(pheno_rep)

# Phenotypic variance
var(pheno_rep$tarsus, na.rm = T)


# load genotype and phenotye files into a gwaa.data object useable in genable. makemap =T orders the SNP markers by total genome.


nids(gwaa0.2)
nsnps(gwaa0.2)


# data descriptives

descriptives.trait(gwaa0.2) #prints out nr records, mean and SD

descriptives.marker(gwaa0.2) #prints out descriptive stats for the markers



# KINSHIP MATRIX ####

data.agemonth.gkin <- ibs(gwaa0.2, weight="freq")


# examine distribution of kinship coefficients

hist(data.agemonth.gkin)

#plot genetic distance between individuals
data.dist <- as.dist(0.5-data.gkin)
data.mds <- cmdscale(data.dist)
plot(data.mds, pch=19, cex=0.5, xlab="MDS1", ylab="MDS2")



# REPEAT GWAS	####


# TARSUS ####
GWAS_tars <- rGLS(formula.FixedEffects= tarsus ~sex + age, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_tars@call$hglm$varRanef[1]-> va
GWAS_tars@call$hglm$varRanef[2]-> vpe
GWAS_tars@call$hglm$varFix-> vr

 va/(va+vpe+vr)
 
# Pull out top 20 SNPs

descriptives.scan(GWAS_tars, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_tars)[1:2,]

# table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_tars),chromosome=chromosome(GWAS_tars), position=map(GWAS_tars),refallele=refallele(GWAS_tars), codeallele=effallele(GWAS_tars), refallelefreq=refallfreq, beta=GWAS_tars[,"effB"], se_beta=GWAS_tars[,"se_effB"], P_val=GWAS_tars[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_tars_repeat.csv",row.names=FALSE)

# Figs

png("Tarsus repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_tars, sub="", main="Tarsus", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot tarsus repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_tars)
estlambda(GWAS_tars[,"P1df"], plot=T, main = "Q-Q Tarsus")
dev.off()


# WING ####
GWAS_wing <- rGLS(formula.FixedEffects=wing~ sex+ age+ hatchyear, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_wing@call$hglm$varRanef[1]-> va
GWAS_wing@call$hglm$varRanef[2]-> vpe
GWAS_wing@call$hglm$varFix-> vr

va/(va+vpe+vr)


#Pull out top 20 SNPs

descriptives.scan(GWAS_wing, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_wing)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_wing),chromosome=chromosome(GWAS_wing), position=map(GWAS_wing),refallele=refallele(GWAS_wing), codeallele=effallele(GWAS_wing), refallelefreq=refallfreq, beta=GWAS_wing[,"effB"], se_beta=GWAS_wing[,"se_effB"], P_val=GWAS_wing[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_wing_repeat.csv",row.names=FALSE)

# Figs

png("Wing repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_wing, sub="", main="Wing", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot wing repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_wing)
estlambda(GWAS_wing[,"P1df"], plot=T, main = "Q-Q Wing")
dev.off()


# MASS ####
GWAS_mass <- rGLS(formula.FixedEffects=mass~ age+ hatchyear+ hatchisland+ month, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_mass@call$hglm$varRanef[1]-> va
GWAS_mass@call$hglm$varRanef[2]-> vpe
GWAS_mass@call$hglm$varFix-> vr

va/(va+vpe+vr)


#Pull out top 20 SNPs

descriptives.scan(GWAS_mass, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_mass)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_mass),chromosome=chromosome(GWAS_mass), position=map(GWAS_mass),refallele=refallele(GWAS_mass), codeallele=effallele(GWAS_mass), refallelefreq=refallfreq, beta=GWAS_mass[,"effB"], se_beta=GWAS_mass[,"se_effB"], P_val=GWAS_mass[,"P1df"], stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_MASS_repeat.csv",row.names=FALSE)

# Figs

png("Mass repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_mass, sub="", main="Mass", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot Mass repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_mass)
estlambda(GWAS_mass[,"P1df"], plot=T, main = "Q-Q Mass")
dev.off()




# billD ####
GWAS_billD <- rGLS(formula.FixedEffects=billD~age+ hatchisland+ hatchyear+ sex+ month, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_billD@call$hglm$varRanef[1]-> va
GWAS_billD@call$hglm$varRanef[2]-> vpe
GWAS_billD@call$hglm$varFix-> vr

va/(va+vpe+vr)


#Pull out top 20 SNPs

descriptives.scan(GWAS_billD, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_billD)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_billD),chromosome=chromosome(GWAS_billD), position=map(GWAS_billD),refallele=refallele(GWAS_billD), codeallele=effallele(GWAS_billD), refallelefreq=refallfreq, beta=GWAS_billD[,"effB"], se_beta=GWAS_billD[,"se_effB"], P_val=GWAS_billD[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_billD_repeat.csv",row.names=FALSE)

# Figs

png("BillD repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_billD, sub="", main="Bill Depth", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot BillD repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_billD)
estlambda(GWAS_billD[,"P1df"], plot=T, main = "Q-Q Bill Depth")
dev.off()




# billL ####
GWAS_billL <- rGLS(formula.FixedEffects=billL~age+ hatchisland+ hatchyear+ month, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_billL@call$hglm$varRanef[1]-> va
GWAS_billL@call$hglm$varRanef[2]-> vpe
GWAS_billL@call$hglm$varFix-> vr

va/(va+vpe+vr)


#Pull out top 20 SNPs

descriptives.scan(GWAS_billL, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_billL)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_billL),chromosome=chromosome(GWAS_billL), position=map(GWAS_billL),refallele=refallele(GWAS_billL), codeallele=effallele(GWAS_billL), refallelefreq=refallfreq, beta=GWAS_billL[,"effB"], se_beta=GWAS_billL[,"se_effB"], P_val=GWAS_billL[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_billL_repeat.csv",row.names=FALSE)

# Figs

png("BillL repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_billL, sub="", main="Bill Length", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot billL repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_billL)
estlambda(GWAS_billL[,"P1df"], plot=T, main = "Q-Q Bill Length")
dev.off()



# totbadge ####
GWAS_totbadge <- rGLS(formula.FixedEffects=totalbadge~age+ hatchisland+ hatchyear+ month, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_totbadge@call$hglm$varRanef[1]-> va
GWAS_totbadge@call$hglm$varRanef[2]-> vpe
GWAS_totbadge@call$hglm$varFix-> vr

va/(va+vpe+vr)


#Pull out top 20 SNPs

descriptives.scan(GWAS_totbadge, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_totbadge)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_totbadge),chromosome=chromosome(GWAS_totbadge), position=map(GWAS_totbadge),refallele=refallele(GWAS_totbadge), codeallele=effallele(GWAS_totbadge), refallelefreq=refallfreq, beta=GWAS_totbadge[,"effB"], se_beta=GWAS_totbadge[,"se_effB"], P_val=GWAS_totbadge[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_totbadge_repeat.csv",row.names=FALSE)

# Figs

png("Totbadge repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_totbadge, sub="", main="Total Badge", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot totbadge repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_totbadge)
estlambda(GWAS_totbadge[,"P1df"], plot=T, main = "Q-Q Total Badge")
dev.off()



# visbadge ####
GWAS_visbadge <- rGLS(formula.FixedEffects=visiblebadge~age+ hatchisland+ hatchyear+ month, genabel.data=gwaa0.2, phenotype.data=pheno_rep, id="id")

GWAS_visbadge@call$hglm$varRanef[1]-> va
GWAS_visbadge@call$hglm$varRanef[2]-> vpe
GWAS_visbadge@call$hglm$varFix-> vr

va/(va+vpe+vr)


#Pull out top 20 SNPs

descriptives.scan(GWAS_visbadge, top=20)
gtdatasum <- summary(gtdata(gwaa0.2))
gtdatasum[1:2,]
refallfreq <- gtdatasum[,"Q.2"]
results(GWAS_visbadge)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(GWAS_visbadge),chromosome=chromosome(GWAS_visbadge), position=map(GWAS_visbadge),refallele=refallele(GWAS_visbadge), codeallele=effallele(GWAS_visbadge), refallelefreq=refallfreq, beta=GWAS_visbadge[,"effB"], se_beta=GWAS_visbadge[,"se_effB"], P_val=GWAS_visbadge[,"P1df"], stringsAsFactors=FALSE)
tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs GWAS_visbadge_repeat.csv",row.names=FALSE)

# Figs

png("Visbadge repeat.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_visbadge, sub="", main="Visible Badge", ylim=c(0,6),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa0.2)), lty=2)
dev.off()

# QQplot

png("QQ plot visbadge repeat.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(GWAS_visbadge)
estlambda(GWAS_visbadge[,"P1df"], plot=T, main = "Q-Q Visible Badge")
dev.off()
