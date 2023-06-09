################################################################################
### Methods Section "Analyses of parallel trait evolution and parallel QTL" ####
################################################################################
#File "re_sized_adjusted_final_qtl_pops_traits.csv" contains the size-corrected data for 64 traits measured on wild-caught lake and stream fish from Boot, Misty, Pye, and Roberts watersheds. Raw data were taken from Stuart et al. (2017) Nature Ecology and Evolution 1:0158 (doi:10.1038/s41559-017-0158).
#Do trait-by-trait linear models to test for effects of habitat (lake vs stream), watershed, and habitat x watershed interaction, and calculate effect sizes (eta squared)

morpho<-read.csv("~/Desktop/LakeStreamQTL/re_sized_adjusted_final_qtl_pops_traits.csv",header=TRUE,sep=",",na.strings="NA",stringsAsFactors=FALSE)
morpho.lm.data <- data.frame(matrix(NA, nrow = ncol(morpho) - 1, ncol = 7))
colnames(morpho.lm.data) <- c("morpho.trait", "etaSq.habitat", "etaSq.watershed", "etaSq.habitat.x.watershed","partial.etaSq.habitat", "partial.etaSq.watershed", "partial.etaSq.habitat.x.watershed")
head(morpho)
habitat <- substr(morpho$fishID.univ,4,4)
watershed <- substr(morpho$fishID.univ,1,3)
trait.names <- colnames(morpho)
cbind(c(1:length(trait.names)), trait.names)

# === f.colname.number ===
# function which tkes a data frame and return a data frame with  
# number of column and its column name
f.colname.number <- function(dataframe){
  #dataframe <- morpho
  namescol <- colnames(dataframe)
  numberscol <- c(1:length(namescol))
  return(data.frame(numberscol, namescol))
}
# === end) f.colname.number ===

f.colname.number(morpho)
m.morpho <- c(2:65) #The morph traits to use for etaSq

# === f.EtaSq from BaylorEdPysch===
EtaSq<-function (x) 
{
    anovaResults <- summary.aov(x)[[1]]
    anovaResultsNames <- rownames(anovaResults)
    SS <- anovaResults[,2] #SS effects and residuals
    k <- length(SS) - 1  # Number of factors 
    ssResid <- SS[k + 1]  # Sum of Squares Residual
    ssTot <- sum(SS)  # Sum of Squares Total
    SS <- SS[1:k] # takes only the effect SS
    anovaResultsNames <- anovaResultsNames[1:k]
    etaSquared <- SS/ssTot # Should be the same as R^2 values
    partialEtaSquared <- SS/(SS + ssResid)
    res <- cbind(etaSquared, partialEtaSquared)
    colnames(res) <- c("Eta^2", "Partial Eta^2")
    rownames(res) <- anovaResultsNames
    return(res)
}
# === end) f.EtaSq ===

for(ctr.m in 1:(length(trait.names) - 1)){
  # ctr.m <- 1
  print(ctr.m)
  morpho.lm.data[ctr.m, 1] <- trait.names[ctr.m + 1]
  m.t <- aov(morpho[, ctr.m + 1] ~ habitat * watershed)
  #summary(m.t)
  eta.sq <- EtaSq(m.t)
  # eta squared for the habitat, watershed and habitat * watershed
  morpho.lm.data[ ctr.m,2:4] <- eta.sq[,1]
  # partial eta squared for the habitat, watershed and habitat * watershed
  morpho.lm.data[ ctr.m,5:7] <- eta.sq[,2]
}

#calculate p-values from anova
morpho.lm.data.p <- data.frame(matrix(NA, nrow = ncol(morpho) - 1, ncol = 4))
colnames(morpho.lm.data.p) <- c("morpho.trait", "habitat.p", "watershed.p", "habitatxwatershed.p")
head(morpho)
habitat <- substr(morpho$fishID.univ,4,4)
watershed <- substr(morpho$fishID.univ,1,3)
trait.names <- colnames(morpho)

cbind(c(1:length(trait.names)), trait.names)
f.colname.number(morpho)
m.morpho <- c(2:65)

for(ctr.m in 1:(length(trait.names) - 1)){
  # ctr.m <- 1
  print(ctr.m)
  morpho.lm.data.p[ctr.m, 1] <- trait.names[ctr.m + 1]
  m.t <- aov(morpho[, ctr.m + 1] ~ habitat * watershed)
  #summary(m.t)[1,1]
  morpho.lm.data.p[ctr.m, 2:4] <- round(summary(m.t)[[1]][1:3, 5], 3)
  
}

morpho.lm.data.all<-cbind(morpho.lm.data,morpho.lm.data.p[,2:4])

#write to file
write.csv(morpho.lm.data.all,file= "~/Desktop/LakeStreamQTL/etas_morpho_linear_model.csv", row.names = FALSE)

# Make Figure 2, following Figure 1 of Stuart et al. (2017)
# Figure 2A
par(mar=c(5,7,2,2))
plot(morpho.lm.data.all$etaSq.habitat.x.watershed,morpho.lm.data.all$etaSq.habitat,pch=16,xlab="habitat x watershed effect size",ylab="habitat effect size
",abline(0,1,lty=2),cex.lab=1.25,las=1)

# make graph for each lake stream pair for each trait in m.morpho
# use habitat and watershed
d.subset <- data.frame(watershed, habitat, morpho[, m.morpho])
head(d.subset)

#Figure 2B
#body_depth_adj, column 47
trait.e <- colnames(d.subset)[47]
etas <- morpho.lm.data[morpho.lm.data$morpho.trait == trait.e,]
title.e <- paste(round(etas$etaSq.habitat.x.watershed,2), ".wsxh_", round(etas$etaSq.habitat,2), ".h_",trait.e, sep = "")
ws.mean <- tapply(d.subset[,47], list(watershed, habitat), mean, na.rm = TRUE)
range.t <- range(ws.mean[, 1:2], na.rm = TRUE)
pdf(paste("~/Desktop/LakeStreamQTL", "/", title.e, ".pdf", sep = "" ), width = 5, height = 5 )
plot(NA, xlim = c(0.5, 2.5), ylim = range.t, xaxt = 'n', xlab = "", ylab = "body depth", las = 1, cex.lab=1.25)
axis(1,c(1,2), c(colnames(ws.mean)))
for(j in 1:nrow(ws.mean)){
    #j <- 1
    points(ws.mean[j,], pch = 19)
    lines(ws.mean[j,])
  }
  tt <- is.na(range(ws.mean[ , 1:2])) == TRUE 
  if( length( tt[tt == TRUE]) > 0 ) {
    mtext("missing value in means", 1, 2)
  }
  dev.off()

#Figure 2C
#right_gill_raker_count, column 65
trait.e <- colnames(d.subset)[65]
etas <- morpho.lm.data[morpho.lm.data$morpho.trait == trait.e,]
title.e <- paste(round(etas$etaSq.habitat.x.watershed,2), ".wsxh_", round(etas$etaSq.habitat,2), ".h_",trait.e, sep = "")
ws.mean <- tapply(d.subset[,65], list(watershed, habitat), mean, na.rm = TRUE)
range.t <- range(ws.mean[, 1:2], na.rm = TRUE)
pdf(paste("~/Desktop/LakeStreamQTL", "/", title.e, ".pdf", sep = "" ), width = 5, height = 5 )
plot(NA, xlim = c(0.5, 2.5), ylim = range.t, xaxt = 'n', xlab = "", ylab = "log right gill raker number", las = 1, cex.lab=1.25)
axis(1,c(1,2), c(colnames(ws.mean)))
 for(j in 1:nrow(ws.mean)){
    #j <- 1
    points(ws.mean[j,], pch = 19)
    lines(ws.mean[j,])
  }

  tt <- is.na(range(ws.mean[ , 1:2])) == TRUE 
  if( length( tt[tt == TRUE]) > 0 ) {
    mtext("missing value in means", 1, 2)
  }
  dev.off()

#Figure 2D
#mean_pelvic_spine_length_adj, column 66
trait.e <- colnames(d.subset)[66]
etas <- morpho.lm.data[morpho.lm.data$morpho.trait == trait.e,]
title.e <- paste(round(etas$etaSq.habitat.x.watershed,2), ".wsxh_", round(etas$etaSq.habitat,2), ".h_",trait.e, sep = "")
ws.mean <- tapply(d.subset[,66], list(watershed, habitat), mean, na.rm = TRUE)
range.t <- range(ws.mean[, 1:2], na.rm = TRUE)
pdf(paste("~/Desktop/LakeStreamQTL", "/", title.e, ".pdf", sep = "" ), width = 5, height = 5 )
plot(NA, xlim = c(0.5, 2.5), ylim = range.t, xaxt = 'n', xlab = "", ylab = "mean pelvic spine length", las = 1, cex.lab=1.25)
axis(1,c(1,2), c(colnames(ws.mean)))
for(j in 1:nrow(ws.mean)){
    #j <- 1
    points(ws.mean[j,], pch = 19)
    lines(ws.mean[j,])
  }
  tt <- is.na(range(ws.mean[ , 1:2])) == TRUE 
  if( length( tt[tt == TRUE]) > 0 ) {
    mtext("missing value in means", 1, 2)
  }
  dev.off()

#Figure 2E
#snout_length_adj, col 49
trait.e <- colnames(d.subset)[49]
etas <- morpho.lm.data[morpho.lm.data$morpho.trait == trait.e,]
title.e <- paste(round(etas$etaSq.habitat.x.watershed,2), ".wsxh_", round(etas$etaSq.habitat,2), ".h_",trait.e, sep = "")
ws.mean <- tapply(d.subset[,49], list(watershed, habitat), mean, na.rm = TRUE)
range.t <- range(ws.mean[, 1:2], na.rm = TRUE)
pdf(paste("~/Desktop/LakeStreamQTL", "/", title.e, ".pdf", sep = "" ), width = 5, height = 5 )
plot(NA, xlim = c(0.5, 2.5), ylim = range.t, xaxt = 'n', xlab = "", ylab = "snout length", las = 1, cex.lab=1.25)
axis(1,c(1,2), c(colnames(ws.mean)))
for(j in 1:nrow(ws.mean)){
    #j <- 1
    points(ws.mean[j,], pch = 19)
    lines(ws.mean[j,])
  }
  tt <- is.na(range(ws.mean[ , 1:2])) == TRUE 
  if( length( tt[tt == TRUE]) > 0 ) {
    mtext("missing value in means", 1, 2)
  }
  dev.off()


#examine correlation between habitat effect size and parallel QTL
#to the file with the etas, I added the number of LGs with QTL (number.LG), the number of parallel QTL(number.parallel.LG), and the proportion of parallel QTL (prop.parallel.LG) by hand based on data from Supplemental Table S4 in the manuscript: "Lake_stream_QTL_Table S4.csv"
etas<-read.csv("~/Desktop/LakeStreamQTL/etas_morpho_linear_model.csv",header=TRUE,sep=",",na.strings="NA",stringsAsFactors=FALSE)

#Figure 3
par(mar=c(5,7,2,2))
plot(etas$etaSq.habitat,etas$prop.parallel.LG,xlab="habitat effect size",ylab="proportion shared LGs",pch=16, cex.lab=1.25)
mean(etas$prop.parallel.LG,na.rm=TRUE)
#0.1544419

library(Hmisc)
corr<-rcorr(etas$etaSq.habitat,etas$prop.parallel.LG,type="spearman")
corr$r
#0.091
corr$P
#0.561


################################################################################
################## Results Section "Caveats of our analyses" ###################
################################################################################

#examine correlation between marker number and QTL number (data taken from Supplemental Table S2)
marker<-c(592,353,230,353,420,715,526,635,440)
QTL<-c(86,37,19,32,7,24,34,6,12)
corr<-rcorr(marker,QTL,type="spearman")
corr$r
#-0.0837
corr$P
#0.831

#examine correlation between F2 number and QTL number (data taken from Supplemental Table S2)
F2<-c(259,274,198,214,72,91,166,70,141)
QTL<-c(86,37,19,32,7,24,34,6,12)
corr<-rcorr(F2,QTL,type="spearman")
corr$r
#0.867
corr$P
#0.0025

#examine correlation between difference in F2 number between crosses and proportion of parallel QTL (data taken from and analyses in Supplemental Table S6)
diff.F2.number<-c(15,16,19,75,94,71,121,322,204,201,83,118)
prop.parallel.QTL<-c(0.18,0.25,0,0.45,0,0.25,0.154,0.15,0.154,0.111,0.026,0.063)
corr<-rcorr(diff.F2.number,prop.parallel.QTL,type="spearman")
corr$r
#-0.267
corr$P
#0.401


########################################################################################
### Results Section "Differences in the extent of gene flow between the populations" ###
########################################################################################

#examine distribution of QTL across linkage groups relative to length of chromosome and gene number on chromosome (data taken from and analyses presented in Supplemental Table S7: "Lake_stream_QTL_TableS7.csv")
#length of chromosome and gene numbers are from original Jones et al. 2012 assembly because we mapped to this assembly
#LG19 is excluded as it was not included in the QTL mapping study
QTL<-c(29,12,9,41,18,8,14,9,13,4,17,21,6,4,2,14,11,4,13,8)
length<-c(0.074066678,0.061216094,0.044142955,0.085752552,0.032194105,0.04489232,0.073413748,0.050896897,0.053211389,0.041144473,0.043900005,0.04835415,0.052774259,0.040064506,0.042566959,0.047604496,0.038373996,0.042787567,0.051851749,0.030791101)

z<-chisq.test(QTL,y=NULL,p=length,simulate.p.value=TRUE)
z
#based on 2000 replicates
#X-squared = 71.099, df = NA, p-value = 0.0004998

z$stdres
#[1]  2.37357803 -0.97122944 -0.71203420  4.22428178  3.43708704 -1.06560527 -1.16410698 -1.15809648 -0.18768032 -2.06461633
#[11]  1.74088772  2.49293450 -2.11003099 -2.00279781 -2.76227095  0.51725424  0.36949575 -2.15647983 -0.09168479  0.03130158

genes<-c(0.069612893,0.047626959,0.051614332,0.073267985,0.040538295,0.039929113,0.073101844,0.048789943,0.056044747,0.045134851,0.058592236,0.055546326,0.053718779,0.040759816,0.043085784,0.044359528,0.03887689,0.042199701,0.051558952,0.025641026)

z<-chisq.test(QTL,y=NULL,p=genes,simulate.p.value=TRUE)
z
#based on 2000 replicates
#X-squared = 72.318, df = NA, p-value = 0.0004998

z$stdres
#[1]  2.72301943 -0.07033105 -1.20243870  5.30722413  2.39801048 -0.72058913 -1.14718238 -1.02473722 -0.38062971 -2.28349668
#[11]  0.51573662  1.83139288 -2.15960127 -2.04273351 -2.78729587  0.78758867  0.32548730 -2.12390566 -0.07070418  0.55655014


#examine correlation between number of traits that map to a LG and number of parallel QTL on that LG (data taken from and analyses presented in Supplemental Table S8)
traits<-c(22,11,8,28,17,8,10,6,12,4,14,17,6,4,2,12,11,4,12,7)
parallel<-c(0.483,0.167,0.222,0.512,0.111,0.000,0.429,0.667,0.154,0.000,0.353,0.286,0.000,0.000,0.000,0.286,0.000,0.000,0.154,0.250)
corr<-rcorr(traits,parallel,type="spearman")
corr$r
#0.583
corr$P
#0.007



