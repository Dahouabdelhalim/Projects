# 1 - LIBRARIES AND FUNCTIONS #####

library(vegan)
library(adegenet)
library(tidyverse)
library(readxl)
library(hierfstat)
library(pegas)
library(genetics)
library(genepop)
library(zvau)

six.digit.genotypes <- function(x) {
  ifelse(x == 0,"000000", ifelse(x < 100000, paste("0",as.character(x),sep=""),as.character(x)))
}

# 2 - LOAD DATA ####

setwd("C:/Users/06023338/OneDrive - Nord universitet/Documents/COLLABORATIONS/11 - Sara Vandamme/EVA/Comparatieve Popgen Flatfish")

brill0 <- read_xlsx("Vandamme_et_al_EVA_2020_ALL.xlsx", sheet = "brill") # brill
turbot0 <- read_xlsx("Vandamme_et_al_EVA_2020_ALL.xlsx", sheet = "turbot") # turbot
sole0 <- read_xlsx("Vandamme_et_al_EVA_2020_ALL.xlsx", sheet = "sole") # sole

setwd("C:/Users/06023338/OneDrive - Nord universitet/Documents/COLLABORATIONS/11 - Sara Vandamme/EVA/Comparatieve Popgen Flatfish/popgen")

# 3 - IDENTIFIERS: AREA, LONG/LAT, YEAR, LOCI, ENVIRONMENTAL VARIABLES ####

popinfo <- c("NewID","Area","Long","Lat","Year","Area.Year")

loci.brill <-c("scor26","scor28","sma7","sma6","scor12","smae32","scor27","scor16","SmaE41","ScoR5","ScoR2","ScoR11","ScoR4","ScoR6")
loci.turbot <- c("Sma3","Sma4","SmaE28","SmaE32","SmaE36","SmaE41","Sma5","SmaE26","SmaE40","Sma2","Sma8","SmaE8","SmaE10","SmaE21")
loci.sole <- c("SolCA13","AC20","AC45","F13","F8_ica9","F8_itg11","Solga12","Sosac6","Sseca28","Ssegata26")

order.brill <- c("BEL10","BEL09","KAT09","SKR09","ENS10","ENS09","CNS07","SNS09","SNS10","EEC07","EEC09","EEC10","WEC10","BCH07","BCH09","SEI09","IRS07","IRS09","WSC09","WIR09","BOB06","BOB07","NWS00")
order.turbot <- c("BEL10","BEL09","KAT09","ENS10","CNS10","CNS07","SNS07","SNS09","EEC07","EEC09","WEC10","BCH07","BCH09","BCH10","SEI09","IRS06","IRS07","IRS09","WIR09","BOB07","BOB09","NWS00","POR00")
order.sole <- c("BEL07","KAT07","SKR07","ENS07","CNS07","CNS08","SNS07","SNS08","EEC08","WEC09","BCH08","IRS08","BOB07")

Order.brill <- factor(order.brill, levels = order.brill)
Order.turbot <- factor(order.turbot, levels = order.turbot)
Order.sole <- factor(order.sole, levels = order.sole)

# 4 - DATA CONVERSIONS AND SUBSETS ####

# 4.1 - Brill ####

# Area:Year identifier

brill0$Area.Year <- factor(paste(brill0$Area,substr(brill0$Year,3,4),sep=""))
brill0$Area.Year

# convert to dataframe

brill1 <- data.frame(brill0)

# convert genotypes to 6 digit format

brill1.loci <- apply(brill1[,loci.brill], 2, six.digit.genotypes)
brill2 <- cbind(brill1[,popinfo],brill1.loci)
head(brill2)
colnames(brill2)

# spatiotemporal extent of the data for RDA

table(factor(brill2$Year))

# 2000 2006 2007 2009 2010 
# 30   18  240  438  153 

table(factor(brill2$Area))

# BCH BEL BOB CNS EEC ENS IRS KAT NWS SEI SKR SNS WEC WIR WSC 
# 46  54  67  66 147  40 146  30  30  78  17  75  35  29  19 

table(factor(brill2$Area.Year))

# BEL09 BEL10 CNS07 EEC07 EEC09 EEC10 ENS09 ENS10 KAT09 SKR09 SNS09 SNS10 
# 38    16    66    37    66    44    15    25    30    17    42    33 

brill <- brill2
dim(brill) # 879 x 20

# exclusion of missing genotypes

gen.brill.all <- df2genind(X=brill[,loci.brill],ncode=3,ploidy=2,NA.char="000")
gen.brill.all$loc.n.all
gen.brill.all@all.names # shows all alleles
dim(gen.brill.all@tab) # dimensions of allelic matrix: 879 x 240
rownames(gen.brill.all@tab) <- 1:nrow(gen.brill.all@tab)
complete.genotypes.brill <- as.numeric(rownames(na.omit(gen.brill.all$tab)))
length(complete.genotypes.brill) # 737
brill.complete <- brill[complete.genotypes.brill,]

# population identifier for popgen

gen.brill.all@pop <- brill$Area.Year
popNames(gen.brill.all)

# write to genepop file

pop(gen.brill.all)
writeGenPop(gen.brill.all, file.name = "brill_genepop", comment = "gen.brill.all")

# 4.2 - Turbot ####

# Area:Year identifier
  
turbot0$Area.Year <- factor(paste(turbot0$Area,substr(turbot0$Year,3,4),sep=""))
turbot0$Area.Year

# convert to dataframe

turbot1 <- data.frame(turbot0)

# convert genotypes to 6 digit format

turbot1.loci <- apply(turbot1[,loci.turbot], 2, six.digit.genotypes)
turbot2 <- cbind(turbot1[,popinfo],turbot1.loci)
head(turbot2)
colnames(turbot2)

# spatiotemporal extent of the data for RDA

table(factor(turbot2$Year))

# 2000 2006 2007 2009 2010 
# 46   21  156  360  165

table(factor(turbot2$Area))

# BCH BEL BOB CNS EEC ENS IRS KAT NWS POR SEI SNS WEC WIR 
# 79  65  43  62  80  53 123  15  27  19  90  50  16  26 

table(factor(turbot2$Area.Year))

# BCH07 BCH09 BCH10 BEL09 BEL10 BOB07 BOB09 CNS07 CNS10 EEC07 EEC09 ENS10 IRS06 IRS07 IRS09 KAT09 NWS00 POR00 SEI09 SNS07 SNS09 WEC10 WIR09 
# 16    20    43    26    39    25    18    48    14    29    51    53    21    20    82    15    27    19    90    18    32    16    26 

turbot <- turbot2

dim(turbot) # 748 x 20

# exclusion of missing genotypes

gen.turbot.all <- df2genind(X=turbot[,loci.turbot],ncode=3,ploidy=2,NA.char="000")
gen.turbot.all$loc.n.all
gen.turbot.all@all.names # shows all alleles
dim(gen.turbot.all@tab) # dimensions of allelic matrix: 911 x 164
rownames(gen.turbot.all@tab) <- 1:nrow(gen.turbot.all@tab)
complete.genotypes.turbot <- as.numeric(rownames(na.omit(gen.turbot.all$tab)))
length(complete.genotypes.turbot) # 841
turbot.complete <- turbot[complete.genotypes.turbot,]

# population identifier for popgen

gen.turbot.all@pop <- turbot$Area.Year
popNames(gen.turbot.all)

# write to genepop file

pop(gen.turbot.all)
writeGenPop(gen.turbot.all, file.name = "turbot_genepop", comment = "gen.turbot.all")

# 4.3 - Sole ####

# Area:Year identifier

sole0$Area.Year <- factor(paste(sole0$Area,substr(sole0$Year,3,4),sep=""))
sole0$Area.Year

# convert to dataframe

sole1 <- data.frame(sole0)

# convert genotypes to 6 digit format

sole1.loci <- apply(sole1[,loci.sole], 2, six.digit.genotypes)
sole2 <- cbind(sole1[,popinfo],sole1.loci)
head(sole2)
colnames(sole2)

# spatiotemporal extent of the data for RDA

table(factor(sole2$Year))

# 2007 2008 2009 
# 636  415   74 

table(factor(sole2$Area))

# BCH BEL BOB CNS EEC ENS IRS KAT SKR SNS WEC 
# 72  40 171  59  45  33  88  71  24 448  74 

table(factor(sole2$Area.Year))

# BCH08 BEL07 BOB07 CNS07 CNS08 EEC08 ENS07 IRS08 KAT07 SKR07 SNS07 SNS08 WEC09 
# 72    40   171    20    39    45    33    88    71    24   277   171    74 

sole <- sole2
dim(sole) # 1125 x 16

# exclusion of missing genotypes

gen.sole.all <- df2genind(X=sole[,loci.sole],ncode=3,ploidy=2,NA.char="000")
gen.sole.all$loc.n.all
gen.sole.all@all.names # shows all alleles
dim(gen.sole.all@tab) # dimensions of allelic matrix: 1125 x 228
rownames(gen.sole.all@tab) <- 1:nrow(gen.sole.all@tab)
complete.genotypes.sole <- as.numeric(rownames(na.omit(gen.sole.all$tab)))
length(complete.genotypes.sole) # 1025
sole.complete <- sole[complete.genotypes.sole,]

# population identifier for popgen

gen.sole.all@pop <- sole$Area.Year
popNames(gen.sole.all)

# write to genepop file

pop(gen.sole.all)
writeGenPop(gen.sole.all, file.name = "sole_genepop", comment = "gen.sole.all")

# save(brill, turbot, sole, 
#     gen.brill.all, gen.turbot.all, gen.sole.all,
#     loci.brill, loci.turbot, loci.sole,
#     popinfo,
#     Order.brill, Order.turbot, Order.sole,
#     six.digit.genotypes, file = "flatfish.popgen.rdata")

# 5 - POPGEN #####

load("flatfish.popgen.rdata")

# 5.1 - Brill ####

# Hardy-Weinberg (pegas)

hw.test.global.brill <- hw.test(gen.brill.all) # global test 
hw.test.global.brill

# We can also perform the test in each population separately
brillbypop <- seppop(gen.brill.all)
names(brillbypop)
OUT <- do.call(rbind,lapply(brillbypop,hw.test));Locus <- rownames(OUT);rownames(OUT) <- 1: nrow(OUT)
Populations.brill <- rep(brill$Area.Year,each=nLoc(gen.brill.all))
HWresults.brill <- data.frame(popNames(gen.brill.all),Locus,OUT)
HWresults.brill
tests <- length(levels(pop(gen.brill.all))) * nLoc(gen.brill.all) # number of pops x number of loci
HWresults.brill.significant <- subset(HWresults.brill,Pr.chi.2... < 0.05/tests)
HWresults.brill.significant
HWresults.brill.significant.exact <- subset(HWresults.brill,Pr.exact < 0.05/tests)
HWresults.brill.significant.exact

# GENVAR (hierfstat/adegenet)

h.gen.brill.all <- genind2hierfstat(gen.brill.all)
sum.final.brill <- summary(gen.brill.all)
Popgen.brill <- basic.stats(h.gen.brill.all);Popgen.brill
AR.brill <- allelic.richness(h.gen.brill.all)
AR.brill.frame <- as.data.frame(AR.brill)
colnames(AR.brill.frame) <- c("min.all",popNames(gen.brill.all))
AR.brill.frame
VAR.brill <- cbind(sum.final.brill$n.by.pop, sum.final.brill$pop.n.all, colMeans(Popgen.brill$Ho), colMeans(Popgen.brill$Hs), colMeans(Popgen.brill$Fis, na.rm=TRUE))
temp1 <- colMeans(AR.brill.frame)
VAR.brill <- cbind(VAR.brill,temp1[-1])
overall <- c(sum(VAR.brill[,1]),mean(VAR.brill[,2]),Popgen.brill$overall[1],Popgen.brill$overall[2],Popgen.brill$overall[9],mean(VAR.brill[,6]))
VAR.brill <- rbind(VAR.brill,overall)
rownames(VAR.brill) <- c(popNames(gen.brill.all),"overall")
colnames(VAR.brill) <- c("N","N alleles", "Hobs","Hexp","Fis","AR")
varfinal.brill <- as.data.frame(VAR.brill)
varfinal.brill

# FIS (hierfstat)

FIS.brill <- boot.ppfis(dat=gen.brill.all,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
varfinal.brill.FIS.CI <- cbind(varfinal.brill[-nrow(varfinal.brill),],FIS.brill$fis.ci)
colnames(varfinal.brill.FIS.CI) <- c("N","N alleles", "Hobs","Hexp","Fis","AR","Fis.ll","Fis.hl")
varfinal.brill.FIS.CI

# FST (hierfstat)

brill.fst <- fstat(gen.brill.all) # Has to be a genind object
brill.fst
brill.pwfst0 <- data.frame(pairwise.fst(gen.brill.all,res.type='matrix'))
brill.pwfst <- brill.pwfst0[order(factor(rownames(brill.pwfst0), levels = Order.brill)), order(factor(rownames(brill.pwfst0), levels = Order.brill))]
brill.pwfst

# FIGURES

pdf("brill.pdf", 12, 8.5)

# CMDSCALE

par(mfrow=c(1,1),pty="s")
plot(cmdscale(as.dist(brill.pwfst)), type="n", xlab = "Dimension 1", ylab= "Dimension 2", main = "CMDS brill")
text(cmdscale(as.dist(brill.pwfst)), rownames(brill.pwfst), cex=0.8)
par(mfrow=c(1,1),pty="m")

# DAPC(adegenet)

dapc.brill <- dapc(gen.brill.all,n.pca=10,n.da = 2)
scatter(dapc.brill)

# PCA (hierfstat)

labels <- gen.brill.all@pop
head(gen.brill.all)

par(mfrow=c(1,1),pty="s")
x <- indpca(gen.brill.all)
#plot(x, cex = 0.7)
#head(x$ipca$li)
#rownames(x$ipca$li)
plot(x$ipca$li[,1], x$ipca$li[,2], type="n", main = "PCA brill", xlab = "PC1", ylab= "PC2", cex=0.7)
text(x$ipca$li[,1], x$ipca$li[,2], labels, cex=0.7)
par(mfrow=c(1,1),pty="m")

dev.off()

# Hardy-Weinberg, LD and Fst (genepop package)

test_HW("brill_genepop", # HW
        which = "Proba",
        outputFile = "brill_genepop_HW", 
        enumeration = FALSE,
        dememorization = 10000,
        batches = 20,
        iterations = 5000,
        verbose = interactive()
)

genepop::Fst("brill_genepop",outputFile = "brill_genepop_Fst") # Fst

test_LD("brill_genepop", # LD,
        outputFile = "brill_genepop_LD",
        dememorization = 10000,
        batches = 100,
        iterations = 5000,
        verbose = interactive())


# 5.2 - Turbot ####

# Hardy-Weinberg (pegas)

hw.test.global.turbot <- hw.test(gen.turbot.all) # global test 
hw.test.global.turbot

# We can also perform the test in each population separately
turbotbypop <- seppop(gen.turbot.all)
names(turbotbypop)
OUT <- do.call(rbind,lapply(turbotbypop,hw.test));Locus <- rownames(OUT);rownames(OUT) <- 1: nrow(OUT)
Populations.turbot <- rep(turbot$Area.Year,each=nLoc(gen.turbot.all))
HWresults.turbot <- data.frame(popNames(gen.turbot.all),Locus,OUT)
HWresults.turbot
tests <- length(levels(pop(gen.turbot.all))) * nLoc(gen.turbot.all) # number of pops x number of loci
HWresults.turbot.significant <- subset(HWresults.turbot,Pr.chi.2... < 0.05/tests)
HWresults.turbot.significant
HWresults.turbot.significant.exact <- subset(HWresults.turbot,Pr.exact < 0.05/tests)
HWresults.turbot.significant.exact

# GENVAR (hierfstat/adegenet)

h.gen.turbot.all <- genind2hierfstat(gen.turbot.all)
sum.final.turbot <- summary(gen.turbot.all)
Popgen.turbot <- basic.stats(h.gen.turbot.all);Popgen.turbot
AR.turbot <- allelic.richness(h.gen.turbot.all)
AR.turbot.frame <- as.data.frame(AR.turbot)
colnames(AR.turbot.frame) <- c("min.all",popNames(gen.turbot.all))
AR.turbot.frame
VAR.turbot <- cbind(sum.final.turbot$n.by.pop, sum.final.turbot$pop.n.all, colMeans(Popgen.turbot$Ho), colMeans(Popgen.turbot$Hs), colMeans(Popgen.turbot$Fis, na.rm=TRUE))
temp1 <- colMeans(AR.turbot.frame)
VAR.turbot <- cbind(VAR.turbot,temp1[-1])
overall <- c(sum(VAR.turbot[,1]),mean(VAR.turbot[,2]),Popgen.turbot$overall[1],Popgen.turbot$overall[2],Popgen.turbot$overall[9],mean(VAR.turbot[,6]))
VAR.turbot <- rbind(VAR.turbot,overall)
rownames(VAR.turbot) <- c(popNames(gen.turbot.all),"overall")
colnames(VAR.turbot) <- c("N","N alleles", "Hobs","Hexp","Fis","AR")
varfinal.turbot <- as.data.frame(VAR.turbot)
varfinal.turbot

# FIS (hierfstat)

FIS.turbot <- boot.ppfis(dat=gen.turbot.all,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
varfinal.turbot.FIS.CI <- cbind(varfinal.turbot[-nrow(varfinal.turbot),],FIS.turbot$fis.ci)
colnames(varfinal.turbot.FIS.CI) <- c("N","N alleles", "Hobs","Hexp","Fis","AR","Fis.ll","Fis.hl")
varfinal.turbot.FIS.CI

# FST (hierfstat)

turbot.fst <- fstat(gen.turbot.all) # Has to be a genind object
turbot.fst
turbot.pwfst0 <- data.frame(pairwise.fst(gen.turbot.all,res.type='matrix'))
turbot.pwfst <- turbot.pwfst0[order(factor(rownames(turbot.pwfst0), levels = Order.turbot)), order(factor(rownames(turbot.pwfst0), levels = Order.turbot))]
turbot.pwfst

# FIGURES

pdf("turbot.pdf", 12, 8.5)

# CMDSCALE

par(mfrow=c(1,1),pty="s")
plot(cmdscale(as.dist(turbot.pwfst)), type="n", xlab = "Dimension 1", ylab= "Dimension 2", main = "CMDS turbot")
text(cmdscale(as.dist(turbot.pwfst)), rownames(turbot.pwfst), cex=0.8)
par(mfrow=c(1,1),pty="m")

# DAPC(adegenet)

dapc.turbot <- dapc(gen.turbot.all,n.pca=10,n.da = 2)
scatter(dapc.turbot)

# PCA (hierfstat)

labels <- gen.turbot.all@pop
head(gen.turbot.all)

par(mfrow=c(1,1),pty="s")
x <- indpca(gen.turbot.all)
#plot(x, cex = 0.7)
#head(x$ipca$li)
#rownames(x$ipca$li)
plot(x$ipca$li[,1], x$ipca$li[,2], type="n", main = "PCA turbot", xlab = "PC1", ylab= "PC2", cex=0.7)
text(x$ipca$li[,1], x$ipca$li[,2], labels, cex=0.7)
par(mfrow=c(1,1),pty="m")

dev.off()

# Hardy-Weinberg, LD and Fst (genepop package)

test_HW("turbot_genepop", # HW
        which = "Proba",
        outputFile = "turbot_genepop_HW", 
        enumeration = FALSE,
        dememorization = 10000,
        batches = 20,
        iterations = 5000,
        verbose = interactive()
)

genepop::Fst("turbot_genepop",outputFile = "turbot_genepop_Fst") # Fst

test_LD("turbot_genepop", # LD,
        outputFile = "turbot_genepop_LD",
        dememorization = 10000,
        batches = 100,
        iterations = 5000,
        verbose = interactive())

# 5.3 - Sole ####

# Hardy-Weinberg (pegas)

hw.test.global.sole <- hw.test(gen.sole.all) # global test 
hw.test.global.sole

# We can also perform the test in each population separately
solebypop <- seppop(gen.sole.all)
names(solebypop)
OUT <- do.call(rbind,lapply(solebypop,hw.test));Locus <- rownames(OUT);rownames(OUT) <- 1: nrow(OUT)
Populations.sole <- rep(sole$Area.Year,each=nLoc(gen.sole.all))
HWresults.sole <- data.frame(popNames(gen.sole.all),Locus,OUT)
HWresults.sole
tests <- length(levels(pop(gen.sole.all))) * nLoc(gen.sole.all) # number of pops x number of loci
HWresults.sole.significant <- subset(HWresults.sole,Pr.chi.2... < 0.05/tests)
HWresults.sole.significant
HWresults.sole.significant.exact <- subset(HWresults.sole,Pr.exact < 0.05/tests)
HWresults.sole.significant.exact

# GENVAR (hierfstat/adegenet)

h.gen.sole.all <- genind2hierfstat(gen.sole.all)
sum.final.sole <- summary(gen.sole.all)
Popgen.sole <- basic.stats(h.gen.sole.all);Popgen.sole
AR.sole <- allelic.richness(h.gen.sole.all)
AR.sole.frame <- as.data.frame(AR.sole)
colnames(AR.sole.frame) <- c("min.all",popNames(gen.sole.all))
AR.sole.frame
VAR.sole <- cbind(sum.final.sole$n.by.pop, sum.final.sole$pop.n.all, colMeans(Popgen.sole$Ho), colMeans(Popgen.sole$Hs), colMeans(Popgen.sole$Fis, na.rm=TRUE))
temp1 <- colMeans(AR.sole.frame)
VAR.sole <- cbind(VAR.sole,temp1[-1])
overall <- c(sum(VAR.sole[,1]),mean(VAR.sole[,2]),Popgen.sole$overall[1],Popgen.sole$overall[2],Popgen.sole$overall[9],mean(VAR.sole[,6]))
VAR.sole <- rbind(VAR.sole,overall)
rownames(VAR.sole) <- c(popNames(gen.sole.all),"overall")
colnames(VAR.sole) <- c("N","N alleles", "Hobs","Hexp","Fis","AR")
varfinal.sole <- as.data.frame(VAR.sole)
varfinal.sole

# FIS (hierfstat)

FIS.sole <- boot.ppfis(dat=gen.sole.all,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
varfinal.sole.FIS.CI <- cbind(varfinal.sole[-nrow(varfinal.sole),],FIS.sole$fis.ci)
colnames(varfinal.sole.FIS.CI) <- c("N","N alleles", "Hobs","Hexp","Fis","AR","Fis.ll","Fis.hl")
varfinal.sole.FIS.CI

# FST (hierfstat)

sole.fst <- fstat(gen.sole.all) # Has to be a genind object
sole.fst
sole.pwfst0 <- data.frame(pairwise.fst(gen.sole.all,res.type='matrix'))
sole.pwfst <- sole.pwfst0[order(factor(rownames(sole.pwfst0), levels = Order.sole)), order(factor(rownames(sole.pwfst0), levels = Order.sole))]
sole.pwfst

# FIGURES

pdf("sole.pdf", 12, 8.5)

# CMDSCALE

par(mfrow=c(1,1),pty="s")
plot(cmdscale(as.dist(sole.pwfst)), type="n", xlab = "Dimension 1", ylab= "Dimension 2", main = "CMDS sole")
text(cmdscale(as.dist(sole.pwfst)), rownames(sole.pwfst), cex=0.8)
par(mfrow=c(1,1),pty="m")

# DAPC(adegenet)

dapc.sole <- dapc(gen.sole.all,n.pca=10,n.da = 2)
scatter(dapc.sole)

# PCA (hierfstat)

labels <- gen.sole.all@pop
head(gen.sole.all)

par(mfrow=c(1,1),pty="s")
x <- indpca(gen.sole.all)
#plot(x, cex = 0.7)
#head(x$ipca$li)
#rownames(x$ipca$li)
plot(x$ipca$li[,1], x$ipca$li[,2], type="n", main = "PCA sole", xlab = "PC1", ylab= "PC2", cex=0.7)
text(x$ipca$li[,1], x$ipca$li[,2], labels, cex=0.7)
par(mfrow=c(1,1),pty="m")

dev.off()

# Hardy-Weinberg, LD and Fst (genepop package)

test_HW("sole_genepop", # HW
        which = "Proba",
        outputFile = "sole_genepop_HW", 
        enumeration = FALSE,
        dememorization = 10000,
        batches = 20,
        iterations = 5000,
        verbose = interactive()
)

genepop::Fst("sole_genepop",outputFile = "sole_genepop_Fst") # Fst

test_LD("sole_genepop", # LD,
        outputFile = "sole_genepop_LD",
        dememorization = 10000,
        batches = 100,
        iterations = 5000,
        verbose = interactive())

