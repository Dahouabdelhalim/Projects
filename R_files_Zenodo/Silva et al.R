####These are the GenABEL and RepeatABEL analyses for the flycatcher results in Silva et al. ####
### in this case, the code refers to the 'thinned' data sets (i.e. Table 4), but the same code was used for analyses for Tables 2 and 4. 

library(RepeatABEL)

###read in repeated measures phenotypes
read.csv("repeatedmeasuresphenotypes.csv")->phenos

###read in mean phenotypes
pheno <- "meanphenos_flycatchers.txt" 

###read in plink files, code for which is described in "Silva et al. Flycatcher GCTA code"
convert.snp.ped("collareds_pruned_R.ped","collareds_pruned_R.map","CFdata.raw", mapHasHeaderLine=FALSE)
geno <- "CFdata.raw"
gwaa <- load.gwaa.data(pheno,geno,makemap=F)

##rename the chromosomes from PLINK so that it works in GenABEL
levels(gwaa@gtdata@chromosome) <- list('0'='0','1'='1','2'='2','3'='3','4'='4','5'='5','6'='6','7'='7','8'='8',
'9'='9','10'='10','11'='11','12'='12','13'='13','14'='14','15'='15',
'17'='17','18'='18','19'='19','20'='20','21'='21',
'22'='22','23'='23','24'='24','25'='25','26'='26', '27'='27','28'='28','29'='29','30'='30','31'='31','32'='32','33'='33','X'='34')

###First round quality control
qc1 <- check.marker(gwaa, call=0.90, perid.call=0.90, maf=0.01, p.lev =1e-03)

data.qc1 <- gwaa[qc1$idok, qc1$snpok]

excluded_birds_qc1 <- setdiff(phdata(gwaa)$id,idnames(data.qc1)) #find excluded ind

excluded_birds_qc1 #prints the list of ringnrs of birds excluded  due to low call rate or out of HWE

### CHECK FOR POPULATION STRUCTURE ###

autosomalMarkers <- which(chromosome(data.qc1)!= "34")
length(autosomalMarkers)
autosomalMarkerNames <- snpnames(data.qc1)[autosomalMarkers]

data.qc1.gkin <- ibs(data.qc1[,autosomalMarkerNames], weight="freq")
data.qc1.dist <- as.dist(0.5-data.qc1.gkin)
data.qc1.mds <- cmdscale(data.qc1.dist,k=4)


# examine the distribution of kinship coeficients

summary(data.qc1.gkin[lower.tri(data.qc1.gkin)])

#MAKE A MDS PLOT OF KINSHIP ACROSS ALL AUTOSOMAL MARKERS
png("MDSplot.png", width=30, height=17.5, units="cm",res=155)
plot(data.qc1.mds[,1],data.qc1.mds[,2], pch=19, xlab="MDS1", ylab="MDS2", cex=0.5)
dev.off()

#outliers <- identify(data.qc1.mds[,1],data.qc1.mds[,2]) #clicking on the outliers puts them into this object

outliers <- c(351, 298, 430, 350,432,431 ) #manually written down the outliers

outlierNames <- phdata(data.qc1)[outliers,]$id #ringnr of the outliers

nonOutlierNames <- phdata(data.qc1)[-outliers,]$id #ringnr of the nonoutliers

data.qc1<-del.phdata(data.qc1, outlierNames)

nids(data.qc1)


###Second round of quality control
qc2 <- check.marker(data.qc1, call=0.95, perid.call=0.95, maf=0.01, p.lev =1e-05, het.fdr=0.1)

data.qc2 <- gwaa[qc2$idok, qc2$snpok]

excluded_birds_qc2 <- setdiff(phdata(data.qc1)$id,idnames(data.qc2)) #find excluded ind

excluded_birds_qc2 #prints the list of ringnrs of birds excluded  due to low call rate or out of HWE

#fix X-linked errors (female heterozygosity)
data.qc2 <- Xfix(data.qc2)


### CHECK FOR POPULATION STRUCTURE ###

autosomalMarkers <- which(chromosome(data.qc2)!= "34")
length(autosomalMarkers)
autosomalMarkerNames <- snpnames(data.qc2)[autosomalMarkers]

data.qc2.gkin <- ibs(data.qc2[,autosomalMarkerNames], weight="freq")
data.qc2.dist <- as.dist(0.5-data.qc2.gkin)
data.qc2.mds <- cmdscale(data.qc2.dist,k=4)

##get minor allele frequencies

sumgt<-summary(gtdata(data.qc2))
afr<-sumgt[, "Q.2"]
maf<-pmin(afr, (1.-afr))



####Grammar Analyses start here
gwaa2<-data.qc2

####Tarsus####
tarsus_h2 <- polygenic(tarsus, kinship.matrix=data.qc2.gkin, gwaa2)

tarsus_h2$esth2

grammar<- grammar(tarsus_h2, data=gwaa2)

lambda(grammar)

descriptives.scan(grammar, top=20)

Vp<-var(phdata(gwaa2)$tarsus, na.rm=T)
Va<-tarsus_h2$esth2*Vp
Ve = Vp - Va
K = compute.GRM(gwaa2)

h2_SE <- get.SEh2_genabel(tarsus ~1, gwaa2, Va, Ve, GRM=K)
cat("Estimated heritability", h2_SE[1], "with SE", h2_SE[2], "\\n")

#Q_Q plot 
png("QQ plot tarsus mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grammar)
estlambda(grammar[,"P1df"], plot=T, main = "Q-Q Tarsus")
dev.off()

png("tarsus mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grammar, df="Pc1df", sub="", main="Tarsus", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

gtdatasum <- summary(gtdata(gwaa2))

gtdatasum[1:2,]

refallfreq <- gtdatasum[,"Q.2"]

results(grammar)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grammar),chromosome=chromosome(grammar), position=map(grammar),
refallele=refallele(grammar), codeallele=effallele(grammar), refallelefreq=refallfreq, n= grammar[,"N"],
beta= grammar[,"effB"], se_beta= grammar[,"se_effB"],maf=maf, P_val= grammar[,"Pc1df"], VSNP=2*(1-maf)*maf*grammar[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs_tarsus_dec27.csv",row.names=FALSE)


####Body Mass#####
bodymass_h2 <- polygenic(body.mass~sex+area, kinship.matrix=data.qc2.gkin, gwaa2)

bodymass_h2$esth2

Vp<-var(phdata(gwaa2)$body.mass, na.rm=T)
Va<-bodymass_h2$esth2*Vp
Ve = Vp - Va
#K = compute.GRM(gwaa2)## I don't need to do this again, this only needs to be done once for the whole data set

h2_SE <- get.SEh2_genabel(body.mass ~sex+area, gwaa2, Va, Ve, GRM=K)
cat("Estimated heritability", h2_SE[1], "with SE", h2_SE[2], "\\n")

grammar<- grammar(bodymass_h2, data=gwaa2)

lambda(grammar)

descriptives.scan(grammar, top=5)

#Q_Q plot 
png("QQ plot mass mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grammar)
estlambda(grammar[,"P1df"], plot=T, main = "Q-Q Body Mass")
dev.off()


png("mass mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grammar, df="Pc1df", sub="", main="Body Mass", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

gtdatasum <- summary(gtdata(gwaa2))

gtdatasum[1:2,]

refallfreq <- gtdatasum[,"Q.2"]

results(grammar)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grammar),chromosome=chromosome(grammar), position=map(grammar),
refallele=refallele(grammar), codeallele=effallele(grammar), refallelefreq=refallfreq, n= grammar[,"N"],
beta= grammar[,"effB"], se_beta= grammar[,"se_effB"],maf=maf, P_val= grammar[,"Pc1df"], VSNP=2*(1-maf)*maf*grammar[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs_bodymass_dec27.csv",row.names=FALSE)


####Wing Length####

winglength_h2 <- polygenic(wing.length~sex, kinship.matrix=data.qc2.gkin, gwaa2)

winglength_h2$esth2

Vp<-var(phdata(gwaa2)$wing.length, na.rm=T)
Va<-winglength_h2$esth2*Vp
Ve = Vp - Va
#K = compute.GRM(gwaa2)## I don't need to do this again, this only needs to be done once for the whole data set

h2_SE <- get.SEh2_genabel(wing.length ~sex, gwaa2, Va, Ve, GRM=K)
cat("Estimated heritability", h2_SE[1], "with SE", h2_SE[2], "\\n")

grammar<- grammar(winglength_h2, data=gwaa2)

	lambda(grammar)

descriptives.scan(grammar, top=5)

#Q_Q plot 
png("QQ plot Wing Length mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grammar)
estlambda(grammar[,"P1df"], plot=T, main = "Q-Q Wing Length")
dev.off()

png("wing mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grammar, df="Pc1df", sub="", main="Wing length", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

gtdatasum <- summary(gtdata(gwaa2))

gtdatasum[1:2,]

refallfreq <- gtdatasum[,"Q.2"]

results(grammar)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grammar),chromosome=chromosome(grammar), position=map(grammar),
refallele=refallele(grammar), codeallele=effallele(grammar), refallelefreq=refallfreq, n= grammar[,"N"],
beta= grammar[,"effB"], se_beta= grammar[,"se_effB"],maf=maf, P_val= grammar[,"Pc1df"], VSNP=2*(1-maf)*maf*grammar[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs_winglength_dec27.csv",row.names=FALSE)

####White on the wings####

wingwhite_h2 <- polygenic(wing.white~sex, kinship.matrix=data.qc1.gkin, gwaa2)

wingwhite_h2$esth2

Vp<-var(phdata(gwaa2)$wing.white, na.rm=T)
Va<-wingwhite_h2$esth2*Vp
Ve = Vp - Va
#K = compute.GRM(gwaa2)## I don't need to do this again, this only needs to be done once for the whole data set

h2_SE <- get.SEh2_genabel(wing.white ~sex, gwaa2, Va, Ve, GRM=K)
cat("Estimated heritability", h2_SE[1], "with SE", h2_SE[2], "\\n")

grammar<- grammar(wingwhite_h2, data=gwaa2)

lambda(grammar)

descriptives.scan(grammar, top=5)

png("QQ plot White mean.png", width=18, height=18, units="cm",res=155,pointsize = 22)
lambda(grammar)
estlambda(grammar[,"P1df"], plot=T, main = "Q-Q White on the Wings")
dev.off()

png("white mean.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(grammar, df="Pc1df", sub="", main="White on the wings", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

gtdatasum <- summary(gtdata(gwaa2))

gtdatasum[1:2,]

refallfreq <- gtdatasum[,"Q.2"]

results(grammar)[1:2,]

#table with effect size and MAF for all SNPs

tabOut <- data.frame(name=snpnames(grammar),chromosome=chromosome(grammar), position=map(grammar),
refallele=refallele(grammar), codeallele=effallele(grammar), refallelefreq=refallfreq, n= grammar[,"N"],
beta= grammar[,"effB"], se_beta= grammar[,"se_effB"],maf=maf, P_val= grammar[,"Pc1df"], VSNP=2*(1-maf)*maf*grammar[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPs_white_dec27.csv",row.names=FALSE)


####RepeatABEL analyses start here#####


library(RepeatABEL)

 subset(phenos, area2!="NA")->phenos2

GWAS_tarsus<- rGLS(tarsus~area2, genabel.data=gwaa2, phenotype.data=phenos2, id="ring")
GWAS_tarsus@call$hglm$varRanef[1]->va
GWAS_tarsus@call$hglm$varRanef[2]->vpe
GWAS_tarsus@call$hglm$varFix->vr
 va/(va+vpe+vr)
 
SE_tarsus<-get.SEh2(tarsus~area2, genabel.data=gwaa2, phenotype.data=phenos2, GWAS.output=GWAS_tarsus, id.name="ring")
 
 
 summary(GWAS_tarsus)
 
 
 tabOut <- data.frame(name=snpnames(GWAS_tarsus),chromosome=chromosome(GWAS_tarsus), position=map(GWAS_tarsus),
refallele=refallele(GWAS_tarsus), codeallele=effallele(GWAS_tarsus), refallelefreq=refallfreq, 
beta= GWAS_tarsus[,"effB"], se_beta= GWAS_tarsus[,"se_effB"],maf=maf,  P_val= GWAS_tarsus[,"P1df"], VSNP=2*(1-maf)*maf*GWAS_tarsus[,"effB"]^2,
stringsAsFactors=FALSE)


tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPstarsus_repeated_dec27.csv",row.names=FALSE)
 
png("tarsus repeated.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_tarsus, sub="", main="Tarsus", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()


png("QQ plot tarsus repeated.png", width=18, height=18, units="cm",res=155,pointsize = 22)
estlambda(GWAS_tarsus[,"P1df"], plot=T, main= "Q-Q Tarsus")
dev.off()


 
GWAS_wing<- rGLS(formula.FixedEffects=wing.length~sex+area2, genabel.data= gwaa2, phenotype.data=phenos2, id="ring")
GWAS_wing@call$hglm$varRanef[1]->va
GWAS_wing@call$hglm$varRanef[2]->vpe
GWAS_wing@call$hglm$varFix->vr
 va/(va+vpe+vr)
 
 SE_wing<-get.SEh2(wing.length~sex+area2, genabel.data= gwaa2, phenotype.data=phenos2, id="ring", GWAS.output=GWAS_wing)
 
 
 summary(GWAS_wing)
 
 tabOut <- data.frame(name=snpnames(GWAS_wing),chromosome=chromosome(GWAS_wing), position=map(GWAS_wing),
refallele=refallele(GWAS_wing), codeallele=effallele(GWAS_wing), refallelefreq=refallfreq, 
beta= GWAS_wing[,"effB"], se_beta= GWAS_wing[,"se_effB"],maf=maf,  P_val= GWAS_wing[,"P1df"], VSNP=2*(1-maf)*maf*GWAS_wing[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPswing.length_repeated_July7.csv",row.names=FALSE)

png("wing repeated.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_wing, sub="", main="Wing Length", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

png("QQ plot wing repeated.png", width=18, height=18, units="cm",res=155,pointsize = 22)
estlambda(GWAS_wing[,"P1df"], plot=T, main= "Q-Q Wing")
dev.off()
 

 ##beak length
  GWAS_mass<- rGLS(formula.FixedEffects=body.mass~sex, genabel.data=gwaa2, phenotype.data=phenos2, id="ring")
GWAS_mass@call$hglm$varRanef[1]->va
GWAS_mass@call$hglm$varRanef[2]->vpe
GWAS_mass@call$hglm$varFix->vr
 va/(va+vpe+vr)
 
  SE_mass<-get.SEh2(wing.length~sex, genabel.data= gwaa2, phenotype.data=phenos2, id="ring", GWAS.output=GWAS_mass)
 
 
 summary(GWAS_mass)
 
 tabOut <- data.frame(name=snpnames(GWAS_mass),chromosome=chromosome(GWAS_mass), position=map(GWAS_mass),
refallele=refallele(GWAS_mass), codeallele=effallele(GWAS_mass), refallelefreq=refallfreq, 
beta= GWAS_mass[,"effB"], se_beta= GWAS_mass[,"se_effB"],maf=maf,  P_val= GWAS_mass[,"P1df"], VSNP=2*(1-maf)*maf* GWAS_mass[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPsbody.mass_repeated_July7.csv",row.names=FALSE)

png("mass repeated.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_mass, sub="", main="Body Mass", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

png("QQ plot mass repeated.png", width=18, height=18, units="cm",res=155,pointsize = 22)
estlambda(GWAS_mass[,"P1df"], plot=T, main= "Q-Q Mass")
dev.off()

 ##tail length
  GWAS_white<- rGLS(formula.FixedEffects=wing.white~sex+area2, genabel.data=gwaa2, phenotype.data=phenos2, id="ring")
GWAS_white@call$hglm$varRanef[1]->va
GWAS_white@call$hglm$varRanef[2]->vpe
GWAS_white@call$hglm$varFix->vr
 va/(va+vpe+vr)
 
   SE_white<-get.SEh2(wing.white~sex+area2, genabel.data= gwaa2, phenotype.data=phenos2, id="ring", GWAS.output=GWAS_white)
   
 summary(GWAS_white)
 
 tabOut <- data.frame(name=snpnames(GWAS_white),chromosome=chromosome(GWAS_white), position=map(GWAS_white),
refallele=refallele(GWAS_white), codeallele=effallele(GWAS_white), refallelefreq=refallfreq, 
beta= GWAS_white[,"effB"], se_beta= GWAS_white[,"se_effB"],maf=maf,  P_val= GWAS_white[,"P1df"], VSNP=2*(1-maf)*maf* GWAS_white[,"effB"]^2,
stringsAsFactors=FALSE)

tabOutSorted <- tabOut[order(tabOut$P_val),]

#write a CVS table with all variants, their effect size, allele frequency and P-value for publication
write.csv(tabOutSorted,"TopSNPswing.white_repeated_July7.csv",row.names=FALSE)

 png("white repeated.png", width=25, height=18, units="cm",res=155,pointsize = 22)
plot(GWAS_white, sub="", main="White on the wings", ylim=c(0,8),col=c("darkblue","darkred"))
abline(h=-log10(0.05/nsnps(gwaa2)), lty=2)
dev.off()

png("QQ plot white repeated.png", width=18, height=18, units="cm",res=155,pointsize = 22)
estlambda(GWAS_white[,"P1df"], plot=T, main= "Q-Q White")
dev.off()

summary(GWAS_tarsus@call$hglm)$SummVC2
summary(GWAS_mass@call$hglm)$SummVC2
summary(GWAS_white@call$hglm)$SummVC2
summary(GWAS_wing@call$hglm)$SummVC2

###to get the standard errors for the RepeatABEL models

source(file="RepeatABEL_ses.R")

h2.SE.tarsus<- get.SEh2(formula.FixedEffects = tarsus ~ area2, genabel.data=gwaa2, phenotype.data=phenos2, GWAS.output=GWAS_tarsus, id.name="ring")

h2.SE.wing<- get.SEh2(formula.FixedEffects = wing.length ~ sex+area2, genabel.data=gwaa2, phenotype.data=phenos2, GWAS.output=GWAS_wing, id.name="ring")

h2.SE.mass<- get.SEh2(formula.FixedEffects = body.mass ~ sex, genabel.data=gwaa2, phenotype.data=phenos, GWAS.output=GWAS_mass, id.name="ring")

h2.SE.white<- get.SEh2(formula.FixedEffects = wing.white ~ sex+area2, genabel.data=gwaa2, phenotype.data=phenos, GWAS.output=GWAS_white, id.name="ring")