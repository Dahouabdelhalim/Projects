# Electronic supplementary material
#
# Agiadi K., Quillevere F., Nawrot R., Sommeville T., Coll M., Koskeridou E., Fietzke J. 
# and Zuschin M., 2021. Community-level decrease in mesopelagic fish size during past 
# climate warming.
# 
# Corresponding author: Konstantina Agiadi
# Department of Paleontology, University of Vienna
# Josef-Holaubek-Platz 2, 1090 Vienna, Austria
# E-mail: konstantina.agiadi@univie.ac.at
#
# --------------------------------------------------------------------------------------------
# R code for estimating fish size from otolith measurements
# Rafal Nawrot, University of Vienna, rafal.nawrot@univie.ac.at 
# Last updated: 11 June 2022
#
#---------------------------------------------------------------------------------------------
library(scales)

##### Importing data ######################################################################## 
species.data<-read.csv("Taxon data.csv", stringsAsFactors=T, header=T, skip=5)
otoliths<-read.csv("Otolith data.csv", stringsAsFactors=T, header=T, skip=4)

str(otoliths)
#add information on the age
otoliths$Age<-factor(otoliths$Sample, levels=levels(otoliths$Sample), 
                     labels=c("MIS 20", "MIS 19", "MIS 18"))
#define colors for ages for plotting
agecols<-c("skyblue", "lightcoral", "skyblue4")

#add names of genera
species.data$Genus<-sub(" .*", "", species.data$Taxon)
table(species.data$Genus)

#add climatic affinity
species.data$Climate<-NA
wm<-(species.data$Tr==1 | species.data$ST==1) & (species.data$Te==0 & species.data$SP==0)
species.data$Climate[wm]<-"Warm-water" 
cl<-(species.data$Te==1 | species.data$SP==1) & (species.data$Tr==0 & species.data$ST==0)
species.data$Climate[cl]<-"Cold-water"
species.data$Climate<-as.factor(species.data$Climate)
species.data$Climate<-factor(species.data$Climate, levels=c("Cold-water", "Warm-water"))

#add information on families and regression coefficients to the otolith dataset
otoliths<-cbind(otoliths, species.data[match(otoliths$Taxon, species.data$Taxon), -2])

#total number of otoliths per species and the number with length (L) with width (W) measured
n.otoliths.sp<-cbind(N.total=tapply(otoliths$L, otoliths$Taxon, length),
                     with.L=tapply(otoliths$L, otoliths$Taxon, function(X) sum(!is.na(X))),
                     with.W=tapply(otoliths$W, otoliths$Taxon, function(X) sum(!is.na(X))),
                     with.LW=tapply(complete.cases(otoliths[,c("L", "W")]), 
                                    otoliths$Taxon, sum),
                     without.L=tapply(otoliths$L, otoliths$Taxon, function(X) sum(is.na(X))))


##### Estimating fish length from otolith measurements ###################################### 
otoliths$Fish.length<-NA

##based on otolith length
#power function: fish length = a * otolith length^b 
x<-otoliths$Type.ofl %in% "Power"
otoliths$Fish.length[x]<-otoliths$a.ofl[x] * otoliths$L[x]^otoliths$b.ofl[x]

#linear function: fish length = a * otolith length + b
x<-otoliths$Type.ofl %in% "Linear"
otoliths$Fish.length[x]<-otoliths$a.ofl[x] * otoliths$L[x] + otoliths$b.ofl[x]

##when length is missing for incomplete specimens use width 
x<-is.na(otoliths$L) & otoliths$Type.ofw %in% "Power" 
otoliths$Fish.length[x]<-otoliths$a.ofw[x] * otoliths$W[x]^otoliths$b.ofw[x]

x<-is.na(otoliths$L) & otoliths$Type.ofw %in% "Linear"
otoliths$Fish.length[x]<-otoliths$a.ofw[x] * otoliths$W[x] + otoliths$b.ofw[x]

#the original fish length estimates are in mm so we convert them into cm
otoliths$Fish.length<-otoliths$Fish.length/10
#add log10-transformed length
otoliths$Log.length<-log10(otoliths$Fish.length)
summary(otoliths$Fish.length)
summary(otoliths$Log.length)


###### Estimating fish weight from fish length ###############################################
otoliths$Fish.weight<-NA

#power function: fish weight = a * fish length^b 
x<-otoliths$Type.flw %in% "Power"
otoliths$Fish.weight[x]<-otoliths$a.flw[x] * otoliths$Fish.length[x]^otoliths$b.flw[x]

#linear function: fish weight = a * fish length + b
x<-otoliths$Type.flw %in% "Linear"
otoliths$Fish.weight[x]<-otoliths$a.flw[x] * otoliths$Fish.length[x] + otoliths$b.flw[x]

#add log10-transformed weight (to make the code simpler later)
otoliths$Log.weight<-log10(otoliths$Fish.weight)
summary(otoliths$Fish.weight)
summary(otoliths$Log.weight)

#save size estimates
size.est<-otoliths[,c("Age", "Specimen", "Taxon", "L", "W", "Fish.length", "Fish.weight")]
write.csv(size.est, "Size estimates.csv", row.names=F)

#median weight each taxon across all samples
x<-tapply(otoliths$Log.weight, otoliths$Taxon, median, na.rm=T)
species.data$Median.log.weight <- x[match(species.data$Taxon, names(x))]

#keep only specimens with estimated weight and exclude indetermined specimens
#this data is used in all size analyses
otolithsR<-droplevels(otoliths[!is.na(otoliths$Fish.weight) & otoliths$Taxon!="indet.",])


##### Basic summary of the results ########################################################## 
#number of all specimens per time interval
n.tot<-c(table(otoliths$Age))
#number of indeterminable specimens
n.indet<-c(table(otoliths$Age[otoliths$Taxon=="indet."]))
#number determined specimens with fish weight estimated
n.fish.weight<-c(table(otolithsR$Age))
#proportion of determined specimens with fish weight estimated per sample
prop.fish.weight<-n.fish.weight/(n.tot-n.indet)

#number of all taxa (excluding indet.)
S.tot<-tapply(otolithsR$Taxon, otolithsR$Age, function(X) length(unique(X)))
#number of identified species
ids<-otoliths$Species.identified=="YES"
S.id<-tapply(otoliths$Taxon[ids], otoliths$Age[ids], function(X) length(unique(X)))

#data on length 
mean.length<-tapply(otolithsR$Fish.length, otolithsR$Age, mean, na.rm=T)
median.length<-tapply(otolithsR$Fish.length, otolithsR$Age, median, na.rm=T)
min.length<-tapply(otolithsR$Fish.length, otolithsR$Age, min, na.rm=T)
max.length<-tapply(otolithsR$Fish.length, otolithsR$Age, max, na.rm=T)
log.median.length<-tapply(otolithsR$Log.length, otolithsR$Age, median, na.rm=T)

#data on weight 
mean.weight<-tapply(otolithsR$Fish.weight, otolithsR$Age, mean, na.rm=T)
median.weight<-tapply(otolithsR$Fish.weight, otolithsR$Age, median, na.rm=T)
min.weight<-tapply(otolithsR$Fish.weight, otolithsR$Age, min, na.rm=T)
max.weight<-tapply(otolithsR$Fish.weight, otolithsR$Age, max, na.rm=T)
log.median.weight<-tapply(otolithsR$Log.weight, otolithsR$Age, median, na.rm=T)


#summary table for the total assemblage
tot.ass<-data.frame(Number.of.all.otoliths=n.tot, 
                    Number.of.indetermined.otoliths=n.indet,
                    Number.of.taxa=S.tot, 
                    Number.of.species=S.id,
                    Number.with.size.estimate=n.fish.weight, 
                    Proportion.with.size.estimate=prop.fish.weight*100,
                    Mean.length=mean.length, Median.length=median.length,
                    Min.length=min.length, Max.length=max.length,
                    Mean.weight=mean.weight, Median.weight=median.weight, 
                    Min.weight=min.weight, Max.weight=max.weight)

#replace dots with spaces and add units
colnames(tot.ass)<-gsub("\\\\.", " ", colnames(tot.ass))
colnames(tot.ass)<-gsub("length", "length (cm)", colnames(tot.ass))
colnames(tot.ass)<-gsub("weight", "weight (g)", colnames(tot.ass))

#save the table (rounded to 3 digits)
write.csv(t(round(tot.ass,3)), file="Table S3 Summary of the otolith sizeâ€“fish size results.csv")


##### Assemblage composition ########## 
#number of specimens with weight estimate per taxon and sample
n.sp.sample<-as.matrix(table(otolithsR$Taxon, otolithsR$Age))
#order according to decreasing total abundance
n.sp.sample.ord<-n.sp.sample[order(rowSums(n.sp.sample), decreasing=T),]
#top taxa
n.sp.sample.ord[1:5,]
#save the results
write.csv(n.sp.sample.ord, file="Table S4 abundances of the taxa.csv")
#relate abundances
prop.sp.sample.ord<-prop.table(n.sp.sample.ord, margin=2)

#number of specimens in each genus per sample
n.gen.sample<-as.matrix(table(otolithsR$Genus, otolithsR$Age))
n.gen.sample.ord<-n.gen.sample[order(rowSums(n.gen.sample), decreasing=T),]
prop.gen.sample.ord<-prop.table(n.gen.sample.ord, margin=2)

#number of specimens in each family per sample
n.fam.sample<-as.matrix(table(otolithsR$Family, otolithsR$Age))
n.fam.sample.ord<-n.fam.sample[order(rowSums(n.fam.sample), decreasing=T),]
prop.fam.sample.ord<-prop.table(n.fam.sample.ord, margin=2)

#relative abundance of families and dominant species 
cairo_pdf("Fig. S3 Assemblage composition.pdf", height=7, width=7)
op<-par(mfrow=c(2,1), mar=c(5.5,4,0,1), oma=c(1,1,1,0), cex.axis=0.9, cex.lab=0.9)
bp<-barplot(t(prop.fam.sample.ord), beside=T, las=2, ylim=c(0,1), 
            ylab="Relative abundance", col=agecols, xaxt="n", legend=T)
text(x=bp[2,], y=-0.06, labels=rownames(prop.fam.sample.ord),
     srt=40, cex=0.7, adj=1, xpd=T)
mtext("A", side=2, line=3.5, at=1, cex=1.6, las=2)
bp<-barplot(t(prop.sp.sample.ord[1:12,]), beside=T, las=2, ylim=c(0,0.4), 
            ylab="Relative abundance", col=agecols, xaxt="n")
text(x=bp[2,], y=-0.02, labels=rownames(prop.sp.sample.ord)[1:12],
     srt=40, cex=0.7, adj=1, xpd=T)
mtext("B", side=2, line=3.5, at=0.4, cex=1.6, las=2)
par(op)
dev.off()


##### Size distribution in the total assemblages ########## 
#Kruskal-Wallis test and pairwise Wilcoxon (= Mann-Whitney) test
kruskal.test(Log.weight ~ Age, data=otolithsR)
pairwise.wilcox.test(otolithsR$Log.weight, g=otolithsR$Age, p.adjust.method="bonferroni")

# % change in median length from MIS 20 to MIS 19 and from MIS 19 to MIS 18
round((median.length["MIS 19"]-median.length["MIS 20"])/median.length["MIS 20"]*100, 2)
round((median.length["MIS 18"]-median.length["MIS 19"])/median.length["MIS 19"]*100, 2)

# % change in median weight from MIS 20 to MIS 19 and from MIS 19 to MIS 18
round((median.weight["MIS 19"]-median.weight["MIS 20"])/median.weight["MIS 20"]*100, 2)
round((median.weight["MIS 18"]-median.weight["MIS 19"])/median.weight["MIS 19"]*100, 2)


#Figure 2 - histograms of fish length and weight (log-transformed)
cairo_pdf("Fig. 2 Fish length and weight distributions (log-transformed).pdf", height=6, width=9)
op<-par(mfcol=c(3,2), mar=c(2,4,2,1), oma=c(3,1,1,0))
for(i in 1:3){
  ok<-otoliths$Age==levels(otoliths$Age)[i]
  hist(otoliths$Log.length[ok], breaks=seq(-0.1, 1.6, by=0.05), main="", 
       las=1, ylab="Number of individuals", col=agecols[i])
  if(i==1) mtext("A. Length", side=3, line=0.5, at=-0.1, cex=1.2)
  if(i==3) mtext(bquote(paste(Log[10], " length (cm)")), side=1, line=3, cex=0.9)
  #add median size and basic information
  abline(v=log.median.length[i], lty=2,  lwd=1.5)
  mtext(levels(otoliths$Age)[i], side=3, at=-0.1, cex=0.8, line=-1.5, adj=0)
  mtext(side=3, at=-0.1, cex=0.7, line=-3.5, adj=0,
        text=paste("N =", n.fish.weight[i], "\\nMedian =", round(median.length[i], 3), "cm"))
}
for(i in 1:3){
  ok<-otoliths$Age==levels(otoliths$Age)[i]
  hist(otoliths$Log.weight[ok], breaks=seq(-3,3, by=0.2), main="", 
       las=1, ylab="Number of individuals", col=agecols[i])
  if(i==1) mtext("B. Weight", side=3, line=0.5, at=-3, cex=1.2)
  if(i==3) mtext(bquote(paste(Log[10], " weight (g)")), side=1, line=3, cex=0.9)
  #add median size and basic information
  abline(v=log.median.weight[i], lty=2,  lwd=1.5)
  mtext(levels(otoliths$Age)[i], side=3, at=-3, cex=0.8, line=-1.5, adj=0)
  mtext(side=3, at=-3, cex=0.7, line=-3.5, adj=0,
        text=paste("N =", n.fish.weight[i], "\\nMedian =", round(median.weight[i], 3), "g"))
}
par(op)
dev.off()

#Figure S2 - histograms of fish length and weight 
cairo_pdf("Fig. S2 Fish length and weight distributions.pdf", height=6, width=9)
op<-par(mfcol=c(3,2), mar=c(2,4,2,1), oma=c(3,1,1,0))
for(i in 1:3){
  ok<-otoliths$Age==levels(otoliths$Age)[i]
  hist(otoliths$Fish.length[ok], breaks=seq(0, 25, by=1), main="", 
       las=1, ylab="Number of individuals", col=agecols[i])
  if(i==1) mtext("A. Length", side=3, line=0.5, at=0, cex=1.2)
  if(i==3) mtext("Length (cm)", side=1, line=3, cex=0.9)
  #add median size and basic information
  abline(v=median.length[i], lty=2,  lwd=1.5)
  mtext(levels(otoliths$Age)[i], side=3, at=15, cex=0.8, line=-1.5, adj=0)
  mtext(side=3, at=15, cex=0.7, line=-3.5, adj=0,
        text=paste("N =", n.fish.weight[i], "\\nMedian =", round(median.length[i], 3), "cm"))
}
for(i in 1:3){
  ok<-otoliths$Age==levels(otoliths$Age)[i]
  hist(otoliths$Fish.weight[ok], breaks=seq(0, 110, by=2), main="", 
       las=1, ylab="Number of individuals", col=agecols[i])
  if(i==1) mtext("B. Weight", side=3, line=0.5, at=-3, cex=1.2)
  if(i==3) mtext("Weight (g)", side=1, line=3, cex=0.9)
  #add median size and basic information
  abline(v=median.weight[i], lty=2,  lwd=1.5)
  mtext(levels(otoliths$Age)[i], side=3, at=66, cex=0.8, line=-1.5, adj=0)
  mtext(side=3, at=66, cex=0.7, line=-3.5, adj=0,
        text=paste("N =", n.fish.weight[i], "\\nMedian =", round(median.weight[i], 3), "g"))
}
par(op)
dev.off()


##### Size trends in dominant taxa ###########################################################
#selected top taxa (all Diaphus species combined together)
top.sp<-c("Hygophum benoiti", "Ceratoscopelus maderensis", "Lobianchia dofleini",
          "Diaphus spp.","Myctophidae indet.")

###boostrapped 95% confidence intervals (CIs) around median weight 
boot.med.fun<-function(x, nrep=10000){
  if(length(x)>=5){ #ignore samples with < 5 specimens
    res<-c()
    for(i in 1:nrep){
      s<-sample(x, replace=T)
      res<-c(res, median(s))
    }
    ci<-quantile(res, c(0.025, 0.975))
  } else {
    ci<-c(NA, NA)
  }
  out<-data.frame(N=length(x), Median=median(x), lCI=ci[1], uCI=ci[2])
  out
}

#median weight +/- CIs in the total assemblage
temp<-split(otolithsR$Log.weight, otolithsR$Age)
boot.weight<-do.call(rbind.data.frame, lapply(temp, boot.med.fun))

#median weight +/- CIs the for top taxa
ts<-droplevels(otolithsR[otolithsR$Taxon %in% top.sp,])
temp<-split(ts$Log.weight, list(ts$Age, ts$Taxon))
boot.weight.sp<-do.call(rbind.data.frame, lapply(temp, boot.med.fun))
boot.weight.sp$Age<-substr(rownames(boot.weight.sp), 1, 6)
boot.weight.sp$Taxon<-substr(rownames(boot.weight.sp), 8, nchar(rownames(boot.weight.sp)))
#for all Diaphus species combined
boot.weight.diaphus<-tapply(otolithsR$Log.weight[otolithsR$Genus=="Diaphus"], 
                            otolithsR$Age[otolithsR$Genus=="Diaphus"], boot.med.fun)
boot.weight.diaphus<-do.call(rbind.data.frame, boot.weight.diaphus)
boot.weight.diaphus$Age<-rownames(boot.weight.diaphus)
boot.weight.diaphus$Taxon<-"Diaphus spp."
boot.weight.sp<-rbind(boot.weight.sp, boot.weight.diaphus)

###plot size trends in the most abundant taxa
#colors for top taxa
spcols<-c("royalblue", "darkorchid2", "darkorange", "brown2", "darkolivegreen")

#scaling of the points according to the relative abundance
prop.sp.top<-prop.sp.sample.ord[top.sp[-4],]
prop.sp.top<-rbind(prop.sp.top, prop.gen.sample.ord["Diaphus",])
rownames(prop.sp.top)[5]<-"Diaphus spp."
spcex<-log(prop.sp.top+1)*10

cairo_pdf("Fig. 3 Median weight in the dominant taxa.pdf", height=5, width=10)
op<-par(mfrow=c(1,2), oma=c(0,5,0,5), mar=c(5, 0.5, 2, 0))
#off<-seq(-0.06, 0.06, by=0.02)
off<-seq(-0.26, 0.26, by=0.13)

labs<-paste0(levels(otolithsR$Age), "\\n",c("glacial", "interglacial", "glacial"))   
plot(boot.weight.sp$Median[1:3], 1:3, type="n", ylim=c(0.7, 3.3), xlim=c(-1.3, 0), las=1,
     ylab="", xlab="", yaxt="n")
axis(2, at=1:3, labels=labs, las=1)   
mtext(bquote(paste(Log[10], " weight (g)")), side=1, line=2.5)
rect(par("usr")[1], par("usr")[3], par("usr")[2], 1.5, col=alpha(agecols[1], 0.3), lwd=0) 
rect(par("usr")[1], 1.5, par("usr")[2], 2.5, col=alpha(agecols[2], 0.3), lwd=0)
rect(par("usr")[1], 2.5, par("usr")[2],  par("usr")[4], col=alpha(agecols[3], 0.3), lwd=0) 

for (i in 1:5){
  ok<-boot.weight.sp[boot.weight.sp$Taxon==top.sp[i],]
  points(ok$Median, 1:3+off[i], col=spcols[i], type="o", cex=1.5, #cex=spcex[top.sp[i],],
         pch=19, lwd=1.5)
  segments(ok$lCI, 1:3+off[i], ok$uCI, 1:3+off[i],col=spcols[i], lwd=1.5)
}
#add patterns in total assemblage
points(log.median.weight, 1:3, col=1, pch=18, type="o", cex=2, lwd=1.2)
segments(boot.weight[,"lCI"], 1:3, boot.weight[,"uCI"], 1:3, col=1, lwd=1.5)

barplot(prop.sp.top*100, beside=T, horiz=T, col=alpha(spcols,0.8), names.arg=rep(NA, 3),
        xlim=c(0, 60), xlab="", xaxt="n")
axis(1, at=seq(0,40, by=10))
mtext("Relative abundance (%)", side=1, line=2.5, at=20)
legend("topright", col=c("black", spcols), leg=c("Total assemblage", top.sp), 
       pch=c(18, rep(19, 5)), cex=0.9)

par(op)
dev.off()


###### Size trends in families ############################################################## 

cairo_pdf("Fig. S4 Fish weight in families.pdf", height=5, width=10)
pa<-par(mfrow=c(2,7), mar=c(5,3,2,1), mgp=c(1.5, 0.7, 0))
for (i in 1:nrow(n.fam.sample.ord)){
  fam<-rownames(n.fam.sample.ord)[i]
  boxplot(Log.weight ~ Age, data=otolithsR[otolithsR$Family %in% fam, ], las=2, xlab="", 
          ylim=c(-2.9,3), main=fam, ylab=bquote(paste(Log[10], " weight (g)")))
  text(lab=n.fam.sample.ord[i,], y=2.8, x=1:3)
  abline(h=median(otolithsR$Log.weight), lty=3)
}
par(op)
dev.off()

kruskal.test(Log.weight ~ Age, data=otolithsR, subset=Family=="Myctophidae")
pairwise.wilcox.test(otolithsR$Log.weight[otolithsR$Family=="Myctophidae"], 
                     g=otolithsR$Age[otolithsR$Family=="Myctophidae"] )
kruskal.test(Log.weight ~ Age, data=otolithsR, subset=Family=="Gadidae")
pairwise.wilcox.test(otolithsR$Log.weight[otolithsR$Family=="Gadidae"], 
                     g=otolithsR$Age[otolithsR$Family=="Gadidae"] )


###### Climatic affinity #################################################################### 
#number of specimens per taxon in each category
table(otolithsR$Taxon, otolithsR$Climate, useNA="ifany")

#number of specimens per climatic affinity in each time interval
n.climate.sample<-table(otolithsR$Climate, otolithsR$Age)
prop.climate.sample<-prop.table(n.climate.sample, margin=2)

median.weight.climate<-tapply(otolithsR$Fish.weight, 
                              list(otolithsR$Age, otolithsR$Climate), median)

#basic statistical tests
kruskal.test(Log.weight ~ Age, data=otolithsR, subset=Climate=="Cold-water")
kruskal.test(Log.weight ~ Age, data=otolithsR, subset=Climate=="Warm-water")
pairwise.wilcox.test(otolithsR$Log.weight[otolithsR$Climate=="Warm-water"], 
                     g=otolithsR$Age[otolithsR$Climate=="Warm-water"] )

#plot relative abundance and size of the climatic groups
climcol<-c("snow1","coral")

cairo_pdf("Fig. S5 climatic affinity.pdf", height=4, width=8)
op<-par(mfrow=c(1,2), mar=c(3,3.1,2,0), oma=c(0,1,0,1), 
        mgp=c(2, 0.7, 0), cex.axis=0.8, cex=0.9)
bp<-barplot(prop.climate.sample*100, beside=F, las=1, ylim=c(0,100), ylab="Relative abundance (%)", 
            col=climcol, xlim=c(0,6), legend=T, args.legend=list(cex=0.9, bty="n"))
mtext("A", side=2, line=2.1, cex=1.4, at=105, las=1)

boxplot(Log.weight ~ Age + Climate, data=otolithsR, xaxt="n", ylim=c(-3, 3.5), xlab="",
        las=1,
        ylab=bquote(paste(Log[10], " weight (g)")), col=c(agecols, agecols))
axis(1, at=1:9, lab=rep(levels(otolithsR$Age), 3))
abline(v=3.5)
abline(h=median(otolithsR$Log.weight), lty=3)
mtext(c("N = ", c(n.climate.sample)), side=3, line=-1.5, cex=0.7, at=c(0.5, 1:6))
mtext(levels(otolithsR$Climate), side=3, line=0.5, cex=1, at=c(2,5))
mtext("B", side=2, line=2.1, cex=1.4, at=4.1, las=1)

par(op)
dev.off()


##### Shifts in relative abundances ##########################################################

#relative abundances of taxa
relab.ord<-as.data.frame.matrix(prop.sp.sample.ord) #convert to data frame
ok<-match(rownames(relab.ord), species.data$Taxo)
#add climate affinity
relab.ord$Climate<-as.character(species.data$Climate[ok])
#set category "Undetermined" for "Myctophidae indet."
relab.ord$Climate[is.na(relab.ord$Climate)]<-"Undet."
relab.ord$Climate<-factor(relab.ord$Climate, levels=c("Undet.","Cold-water","Warm-water"))
#add median weight across all samples
relab.ord$Log.weight<-species.data$Median.log.weight[ok]
#add median weight in each time interval
median.weight.sp<-tapply(otolithsR$Log.weight, list(otolithsR$Taxon, otolithsR$Age), median)
median.weight.sp<-median.weight.sp[rownames(relab.ord),]
colnames(median.weight.sp)<-paste("Log.weight", colnames(median.weight.sp))
relab.ord<-cbind(relab.ord, median.weight.sp)

#plot shift in abundance vs. specie median size
cairo_pdf("Fig. S6 shifts in relative abundance.pdf", height=6, width=12)
op<-par(mfrow=c(1,2), oma=c(0,1,0,0), las=1)
#between MIS 20 and MIS 19
relab.ord1<-relab.ord[rowSums(relab.ord[,1:2])>0,] #omit taxa absent in both intervals
abund.shift1<-relab.ord1[,2]-relab.ord1[,1] #shift in abundance
s<-abund.shift1 > 0.05 | abund.shift1 < -0.05 #which change of more than 5%
#new appearances - green/up triangle, disappearance - red/down triangle
colb<-ifelse(relab.ord1[,1]==0, 4, ifelse(relab.ord1[,2]==0, 2, 1))
pchb<-ifelse(relab.ord1[,2]==0, 24, ifelse(relab.ord1[,3]==0, 25, 21))
#shift in all Myctophidae combined
abs.myc1<-prop.fam.sample.ord[1,2]-prop.fam.sample.ord[1,1]
med.weigth.myct<-median(otolithsR$Log.weight[otolithsR$Family=="Myctophidae"])

plot(relab.ord1$Log.weight, abund.shift1, xlim=c(-3, 3), ylim=c(-0.19, 0.19),
     xlab=bquote(paste("Median ", log[10], " weight (g) of a taxon")), 
     ylab="Change in relative abundance", main="MIS 20 to MIS 19",
     pch=pchb, bg=c("black", climcol)[relab.ord1$Climate], cex=1.1, col=1)
points(med.weigth.myct, abs.myc1, pch=15, cex=1.2) #Myctophidae
text(med.weigth.myct, abs.myc1, lab="Myctophidae", pos=1, cex=0.8)
abline(h=0, lty=2)
text(relab.ord1$Log.weight[s], abund.shift1[s], lab=rownames(relab.ord1)[s], pos=4, cex=0.8)
legend("topleft", leg=levels(species.data$Climate), pch=21, pt.bg=climcol, cex=0.8, bty="n")
legend("bottomleft", leg=c("Present in both intervals","Disappearnce", "Reappearance"), 
       pch=c(21, 25, 24), col=1, cex=0.8, bty="n")

#between MIS 19 and MIS 18
relab.ord2<-relab.ord[rowSums(relab.ord[,2:3])>0,] #omit taxa absent in both intervals
abund.shift2<-relab.ord2[,3]-relab.ord2[,2] #shift in abundance
#which change of more than 5%
s<-abund.shift2 > 0.05 | abund.shift2 < -0.05 #which change of more than 5%
#new appearances - green/up triangle, disappearance - red/down triangle
colb<-ifelse(relab.ord2[,2]==0, 4, ifelse(relab.ord2[,3]==0, 2, 1))
pchb<-ifelse(relab.ord2[,2]==0, 24, ifelse(relab.ord2[,3]==0, 25, 21))
#shift in all Myctophidae combined
abs.myc2<-prop.fam.sample.ord[1,3]-prop.fam.sample.ord[1,2]

plot(relab.ord2$Log.weight, abund.shift2, xlim=c(-3, 3), ylim=c(-0.19, 0.19),
     xlab=bquote(paste("Median ", log[10], " weight (g) of a taxon")), 
     ylab="Change in relative abundance", main="MIS 19 to MIS 18",
     pch=pchb, bg=c("black", climcol)[relab.ord2$Climate], cex=1.1, col=1)
points(med.weigth.myct, abs.myc2, pch=15, cex=1.2) #Myctophidae
text(med.weigth.myct, abs.myc2, lab="Myctophidae", pos=1, cex=0.8)
abline(h=0, lty=2)
text(relab.ord2$Log.weight[s], abund.shift2[s], lab=rownames(relab.ord2)[s], pos=4, cex=0.8)

par(op)
dev.off()


##### Otolith preservation ################################################################## 
#for taphonomic analyses all specimens (including indet.) are used

#names of taphonomic signatures
taphs<-c("Completeness", "Translucency", "Bioerosion", "Edge.preservation", 
         "Dissolution", "Ornamentation.loss")

###patterns in average taphonomic signature
otoliths$Mean.taph.score<-apply(otoliths[,taphs], 1, mean, na.rm=T)

#kruskal-wallis test and pairwise Wilcoxon (= Mann-Whitney) test
kruskal.test(Mean.taph.score ~ Age, data=otoliths)

cairo_pdf("Fig. S7 average taphonomic signature.pdf", height=6, width=6)
boxplot(Mean.taph.score ~ Age, data=otoliths, ylim=c(0,3), las=1, col=agecols, 
        xlab="", ylab="Taphonomic score")
dev.off()

###patterns in each taphonomic signature separately
taphs2<-gsub("\\\\.", " ", taphs) #remove "." for plotting

#average score for each signature per age
for (i in 1:6){
  print(taphs[i])
  print(round(tapply(otoliths[,taphs[i]], otoliths$Age, mean, na.rm=T), 2))
}

#Number and proportion of otolith specimens in each score category
taph.sig.n<-list()
taph.sig.prop<-list()
for (i in 1:6){
  taph.sig.n[[i]]<-table(otoliths$Age, otoliths[,taphs[i]])
  taph.sig.prop[[i]]<-prop.table(taph.sig.n[[i]], margin=1)
}
names(taph.sig.n)<-taphs
names(taph.sig.prop)<-taphs


cairo_pdf("Fig. S8 taphonomic signatures.pdf", height=7, width=7)
op<-par(mfrow=c(2,3))
for(i in 1:6) {
  if(i==1){
    barplot(t(taph.sig.prop[[i]]), main=taphs2[i], las=1, ylab="Proportion of specimens",
            legend=T, args.legend=list(x=3.5, y=0.4, cex=1.1))
  } else {
    barplot(t(taph.sig.prop[[i]]), main=taphs2[i], las=1, ylab="Proportion of specimens")
  }
}
par(op)
dev.off()

