#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# analyses for Addis & Lowe "Environmentally associated variation in dispersal distance affects inbreeding risk in a stream salamander." The American Naturalist.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(lme4)
library(nlme)
library(plotrix)
library(moments)
library(psych)
library(emmeans)
library(plyr)
library(ggplot2)
library(tidyr)
library(effects)
library(adegenet)
library(ecodist)
library(h2o)
library(ggpubr)
library(ggpattern)  ## 

#install 'related' package (not on CRAN)
install.packages("related", repos="http://R-Forge.R-project.org")
library(related)

######################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#dispersal distance analyses, using all recaptured individuals
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data=read.table("Dispersal.csv", sep=",", header=T)
nrow(data) #663 individuals

#subset to movements within reaches (< 500m), then by dispersal status and stream reach
disp=data[which(data$absolute_distance<500 & data$absolute_distance >=10),]
nodisp=data[which(data$absolute_distance<500 & data$absolute_distance < 10),]
fishdisp=disp[which(disp$reach=="downstream"),] 
fishnodisp=nodisp[which(nodisp$reach=="downstream"),]
fishlessdisp=disp[which(disp$reach=="upstream"),]
fishlessnodisp=nodisp[which(nodisp$reach=="upstream"),]
nrow(fishdisp) #43 individuals
nrow(fishlessdisp) #50 individuals

#test for difference in proportion dispersing between reaches
prop.test(x=c(43,50), n=c(43+184,50+384), conf.level=0.95)
43/184 #raw proportion dispersing downstream = 23.4%
50/384 #raw proportion dispersing upstream= 13%

# compute mean distances of dispersers in upstream and downstream reaches
mean(fishdisp$absolute_distance)
std.error(fishdisp$absolute_distance)
mean(fishlessdisp$absolute_distance)
std.error(fishlessdisp$absolute_distance)

# test for diffs in dispersal between upstream and downstream reaches
t.test(fishdisp$absolute_distance, fishlessdisp$absolute_distance)
wilcox.test(fishdisp$absolute_distance, fishlessdisp$absolute_distance)

#use log-transformed distance
log.fishdisp<-log(fishdisp$absolute_distance)
log.fishlessdisp<-log(fishlessdisp$absolute_distance)
t.test(log.fishdisp, log.fishlessdisp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test for sex-biased dispersal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#subset by sex and dispersal status
males<-data[which(data$sex=="M" & data$absolute_distance <=500),]
nrow(males)#45
male.disp<-males[which(males$absolute_distance>=10 & males$absolute_distance<=500),]
nrow(male.disp)#11
females<-data[which(data$sex=="F" & data$absolute_distance <=500),]
nrow(females)#80
female.disp<-females[which(females$absolute_distance>=10 & females$absolute_distance<=500),]
nrow(female.disp)#16

#subset by sex and reach
male.up<-data[which(data$sex=="M" & data$reach=="upstream"),]
nrow(male.up)#33
male.down<-data[which(data$sex=="M" & data$reach=="downstream"),]
nrow(male.down)#13
female.up<-data[which(data$sex=="F" & data$reach=="upstream"),]
nrow(female.up)#55
female.down<-data[which(data$sex=="F" & data$reach=="downstream"),]
nrow(female.down)#25

#test for difference in proportion dispersing >= 10m
prop.test(x=c(11,16), n=c(11+34,16+64), conf.level = 0.95)

#means
mean(males$absolute_distance)
std.error(males$absolute_distance)
mean(females$absolute_distance)
std.error(females$absolute_distance)

#log-transformed distances
log.male<-log(males$absolute_distance+1)
log.female<-log(females$absolute_distance+1)
t.test(log.male, log.female)
wilcox.test(log.male, log.female)

#test for differences in proportion dispersing between reaches
#subset males by dispersal status and reach
male.disp.up<-male.up[which(male.up$absolute_distance>9),]
nrow(male.disp.up)#9
male.nodisp.up<-male.up[which(male.up$absolute_distance<=9),]
nrow(male.nodisp.up)#24
male.disp.down<-male.down[which(male.down$absolute_distance>9),]
nrow(male.disp.down)#3
male.nodisp.down<-male.down[which(male.down$absolute_distance<=9),]
nrow(male.nodisp.down)#10
prop.test(x=c(9,3), n =c(9+24, 3+10),conf.level=0.95)

#subset females by dispersal status and reach
female.disp.up<-female.up[which(female.up$absolute_distance>9),]
nrow(female.disp.up)#9
female.nodisp.up<-female.up[which(female.up$absolute_distance<=9),]
nrow(female.nodisp.up)#46
female.disp.down<-female.down[which(female.down$absolute_distance>9),]
nrow(female.disp.down)#7
female.nodisp.down<-female.down[which(female.down$absolute_distance<=9),]
nrow(female.nodisp.down)#18
prop.test(x=c(9,7), n =c(9+46, 7+18),conf.level=0.95)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plots of dispersal distance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#calculate group means
mu<-ddply(disp, "reach", summarise, grp.mean=mean(absolute_distance))

ggplot(disp, aes(x=absolute_distance,color=reach, fill=reach))+
  geom_histogram(alpha=0.7)+
  facet_grid(reach ~ .)+
  geom_vline(data=mu,aes(xintercept=grp.mean, color=reach),linetype="dashed")+
  geom_density(alpha=0.6)+
  scale_color_manual(values=c("#000000","#000000"))+
  scale_fill_manual(values=c("#000000","#999999"))+
  labs(x="Dispersal distance (m)", y ="Number of individuals")+
  theme_classic()+
  theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14), 
        axis.text.x=element_text(size=14),strip.background = element_blank(),
        strip.text.y = element_blank(),legend.title=element_blank(),
        legend.position="top")


#####################################################################3
#Mantel test

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate geographic distances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xycoords=read.table("XYcoords.csv", sep=",", header=T)
str(xycoords)

#subset by stream
parxy=xycoords[which(xycoords$stream=="paradise"),]
bearxy=xycoords[which(xycoords$stream=="bear"),]
zigxy=xycoords[which(xycoords$stream=="zigzag"),]
casxy=xycoords[which(xycoords$stream=="cascade"),]
canxy=xycoords[which(xycoords$stream=="canyon"),]

#remove extra columns
parxy$ID<-NULL
parxy$stream<-NULL
bearxy$ID<-NULL
bearxy$stream<-NULL
zigxy$ID<-NULL
zigxy$stream<-NULL
casxy$ID<-NULL
casxy$stream<-NULL
canxy$ID<-NULL
canxy$stream<-NULL


as.matrix(parxy)
as.matrix(bearxy)
as.matrix(zigxy)
as.matrix(casxy)
as.matrix(canxy)

par.geo<-dist(parxy) 
bear.geo<-dist(bearxy)
zig.geo<-dist(zigxy)
cas.geo<-dist(casxy)
can.geo<-dist(canxy)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate genetic distances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in structure files for each stream
bear.str<-read.structure("bear_297loci_structure.str") #150 genotypes, 297 loci, 1,2,none,1,n
cascade.str<-read.structure("cascade_297loci_structure.str")#22 genotypes, 297 loci, 1,2,none,1,n
canyon.str<-read.structure("canyon_297loci_structure.str")#36 genotypes, 297 loci, 1,2,none,1,n
paradise.str<-read.structure("paradise_297loci_structure.str")#112 genotypes, 297 loci, 1,2,none,1,n
zigzag.str<-read.structure("zigzag_297loci_structure.str")#62 genotypes, 297 loci, 1,2,none,1,n


#calculate genetic distances for each stream
bear.distgenEUCL <- dist(bear.str, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
cascade.distgenEUCL <- dist(cascade.str, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
canyon.distgenEUCL <- dist(canyon.str, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
paradise.distgenEUCL <- dist(paradise.str, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
zigzag.distgenEUCL <- dist(zigzag.str, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

#mantel Rs
mantel.bear<-mantel(bear.distgenEUCL~bear.geo,  nperm = 999,nboot = 500, pboot = 0.9, cboot = 0.95)
mantel.cascade<-mantel(cascade.distgenEUCL~cas.geo,  nperm = 999,nboot = 500, pboot = 0.9, cboot = 0.95)
mantel.canyon<-mantel(canyon.distgenEUCL~can.geo,  nperm = 999,nboot = 500, pboot = 0.9, cboot = 0.95)
mantel.paradise<-mantel(paradise.distgenEUCL~par.geo,  nperm = 999,nboot = 500, pboot = 0.9, cboot = 0.95)
mantel.zigzag<-mantel(zigzag.distgenEUCL~zig.geo,  nperm = 999,nboot = 500, pboot = 0.9, cboot = 0.95)

#mantel plots for each stream
breaks<-c(-50,50,150,250,350,450,550,650,750,850,950,1050,1150,1250,1350,1450)
par.gram <- mgram(paradise.distgenEUCL, par.geo,breaks=breaks,nperm=999,nboot = 500, pboot = 0.9, cboot = 0.95)
plot(par.gram)
bear.gram<-mgram(bear.distgenEUCL,bear.geo,breaks=breaks,nperm=999,nboot = 500, pboot = 0.9, cboot = 0.95)
plot(bear.gram)
zig.gram<-mgram(zigzag.distgenEUCL,zig.geo,breaks=breaks,nperm=999,nboot = 500, pboot = 0.9, cboot = 0.95)
plot(zig.gram)
can.gram<-mgram(canyon.distgenEUCL,can.geo, breaks=breaks,nperm=999,nboot = 500, pboot = 0.9, cboot = 0.95)
plot(can.gram)
cas.gram<-mgram(cascade.distgenEUCL,cas.geo, breaks=breaks,nperm=999,nboot = 500, pboot = 0.9, cboot = 0.95)
plot(cas.gram)

#save data for paradise and bear for plotting (because these streams have significant IBD)
par.gram<-as.data.frame(par.gram$mgram)
par.gram<-par.gram[par.gram$ngroup != 0,] # remove rows at distances with no data
bear.gram<-as.data.frame(bear.gram$mgram)

# combine data frames and add identifiers for stream and significance level
bear.par<-rbind(par.gram, bear.gram)
bear.par$stream<-c(rep("paradise", 13), rep("bear", 15))
bear.par$sig<-c(rep("y",1), rep("n",1), rep("y", 1), rep("n",10), 
                rep("y",1), rep("n",1), rep("y", 1), rep("n", 6), rep("y", 2), rep("n",4 ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot mantel correlograms for bear and paradise
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(bear.par, aes(x=lag, y=mantelr, group=stream)) +
  geom_hline(aes(yintercept=0), colour="grey")+
  geom_line(aes(linetype=stream), size=1,show.legend = FALSE)+
  geom_point(aes(shape=sig),size=3,show.legend = FALSE)+
  scale_shape_manual(values=c(1,16))+
  ylim(-.15,.15)+
  xlab("Distance (m)")+
  ylab("Mantel r correlation")+
  scale_x_continuous(breaks=seq(0,1400,200))+
  theme_classic()+
  theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14), axis.text.x=element_text(size=14))


#############################################################################
### get estimates of relatedness using coancestry############
input <- readgenotypedata("297loci_genotypes_related.txt")

#look at data frames
input$gdata
input$nloci
input$nalleles
input$ninds
input$freqs

#compare relatedness estimators (only moments based estimators (but not ritland), not likelihood)
out<-compareestimators(input, 100)

#Correlation Coefficients Between Observed & Expected Values:
#wang		0.955558
#lynchli		0.955932
#lynchrd		0.962566  # so this is best
#quellergt	0.961201

# now do custom comparisons including two likelihood methods and lynchrd estimator

##simulations to determine which simulator is best
sim<-familysim(input$freqs, 100)
outputsim<-coancestry(sim, dyadml=1, trioml=1,lynchrd=1) ## this takes a long time
simrel<-cleanuprvals(outputsim$relatedness, 100)
write.table(as.data.frame(simrel),file="simrel.csv",quote=F,sep=",")
simrel=read.table("simrel297.csv", sep=",",header=T)

#parse data based on relatedness type and estimator
triomlpo<-simrel[1:100, 5]
triomlfs<-simrel[(100+1):(2*100), 5]
triomlhs<-simrel[((2*100)+1):(3*100),5]
triomlur<-simrel[((3*100)+1):(4*100),5]

dyadmlpo<-simrel[1:100, 11]
dyadmlfs<-simrel[(100+1):(2*100), 11]
dyadmlhs<-simrel[((2*100)+1):(3*100),11]
dyadmlur<-simrel[((3*100)+1):(4*100),11]

lynchrdpo<-simrel[1:100, 8]
lynchrdfs<-simrel[(100+1):(2*100), 8]
lynchrdhs<-simrel[((2*100)+1):(3*100),8]
lynchrdur<-simrel[((3*100)+1):(4*100),8]

#create labels for the different estimators
trioml<- rep('tri', 100)
dyadml<- rep('di', 100)
lynchrd<- rep('LR', 100)
estimator2<-c(trioml,dyadml,lynchrd)
Estimator<- rep(estimator2, 4)

#labels for relatedness types
po<-rep("Parent-Offspring", (3*100))
fs<-rep("Full-Sibs", (3*100))
hs<-rep("Half-Sibs", (3*100))
ur<-rep("Unrelated", (3*100))
relationship<-c(po,fs,hs,ur)

#combine different values for each estimator based on relatedness type, as lists
relatednesspo<-c(triomlpo, dyadmlpo, lynchrdpo)
relatednessfs<-c(triomlfs, dyadmlfs, lynchrdfs)
relatednesshs<-c(triomlhs, dyadmlhs, lynchrdhs)
relatednessur<-c(triomlur, dyadmlur, lynchrdur)
Relatedness_Value<-c(relatednesspo, relatednessfs, relatednesshs, relatednessur)

#combine the data
combineddata<-as.data.frame(cbind(Estimator, relationship, Relatedness_Value))
combineddata$Relatedness_Value<-as.numeric(as.character(combineddata$Relatedness_Value))

#plot the data
ggplot(combineddata, aes(x=Estimator, y=Relatedness_Value), ylim=c(-0.5,1.0))+
  geom_boxplot() +
  facet_wrap (~ relationship)

#calculate correlation coefficient between the observed values and expected values
urval<-rep(0, 100)
hsval<-rep(0.25, 100)
fsval<-rep(0.5, 100)
poval<-rep(0.5, 100)
relvals<-c(poval,fsval,hsval,urval)

cor(relvals, simrel[,5])
# 0.9760151, highest is trioml
cor(relvals, simrel[,8])
# 0.9685344
cor(relvals, simrel[,11])
# 0.9754507

# determine cutoffs for relatedness categories for trioml estimator
relvalues<-simrel[, 5]
label1<-rep("PO", 100)
label2<-rep("Full",100)
label3<-rep("Half",100)
label4<-rep("Unrelated",100)
labels<-c(label1,label2,label3,label4)
plot(as.factor(labels), relvalues, ylab="Relatedness Value", xlab="Relatedness")

pdf("FigS1.pdf", useDingbats = FALSE)
qplot(as.factor(labels), relvalues, geom ="boxplot", ylab ="Relatedness Values", xlab ="Relatedness")
dev.off()

PO=simrel[which(simrel$group=="POPO"),]
SB=simrel[which(simrel$group=="SBSB"),]
HS=simrel[which(simrel$group=="HSHS"),]
UR=simrel[which(simrel$group=="URUR"),]

#these are the cutoffs for relationship categories
quantile(PO$trioml,probs=c(.025,.975)) #0.500-0.609
quantile(SB$trioml,probs=c(.025,.975)) #0.422-0.613
quantile(HS$trioml,probs=c(.025,.975))#0.131-0.384
quantile(UR$trioml,probs=c(.025,.975))#0-0.123

hist(PO$trioml, breaks=seq(0,1,.01), freq=T, xlim=c(0,1),ylim=c(0,20),col="grey", border = "black", xlab="Relatedness", ylab="Frequency", cex.axis=1.8, cex.lab=1.8,main = "Parent-offspring")
qts<-quantile(PO$trioml,probs=c(.025,.975))
abline(v=qts[1],col="red")
abline(v=qts[2],col="red")
hist(SB$trioml, breaks=seq(0,1,.01), freq=T, xlim=c(0,1),ylim=c(0,20),col="grey", border = "black", xlab="Relatedness", ylab="Frequency", cex.axis=1.8, cex.lab=1.8,main = "Full-siblings")
qts<-quantile(SB$trioml,probs=c(.025,.975))
abline(v=qts[1],col="red")
abline(v=qts[2],col="red")
hist(HS$trioml, breaks=seq(0,1,.01), freq=T, xlim=c(0,1),ylim=c(0,20),col="grey", border = "black", xlab="Relatedness", ylab="Frequency", cex.axis=1.8, cex.lab=1.8,main = "Half-siblings")
qts<-quantile(HS$trioml,probs=c(.025,.975))
abline(v=qts[1],col="red")
abline(v=qts[2],col="red")
hist(UR$trioml, breaks=seq(-.2,1,.01), freq=T, xlim=c(-.2,1),ylim=c(0,40),col="grey", border = "black", xlab="Relatedness", ylab="Frequency", cex.axis=1.8, cex.lab=1.8,main = "Unrelated")
qts<-quantile(UR$trioml,probs=c(.025,.975))
abline(v=qts[1],col="red")
abline(v=qts[2],col="red")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## now calculate relatedness using 297 snps and trioml estimator
# multithread the process to take less time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h2o.init(nthreads=8)
output<-coancestry("297loci_genotypes_related.txt", trioml=2)
write.table(as.data.frame(output$relatedness),file="5pops_relatedness_trioml.csv",quote=F,sep=",")
write.table(as.data.frame(output$relatedness.ci95),file="5pops_relatedness_ci95_trioml.csv",quote=F,sep=",")
# this yields 72,771 estimates of pairwise relatedness (pairwise for 382 individuals)


###################################################################################
# calculate proportion of individuals within 50m that are relatives
relatedness<-read.csv("5pops_relatedness_trioml.csv", stringsAsFactors=FALSE, strip.white=TRUE)
str(relatedness)
nrow(relatedness)

#delete extra columns
relatedness$lynchli<-NULL
relatedness$lynchrd<-NULL
relatedness$ritland<-NULL
relatedness$quellergt<-NULL
relatedness$dyadml<-NULL
relatedness$wang<-NULL

# only keep relatedness estimates among individuals within same stream
#add column designating within-stream comparisons
relatedness$stream<-relatedness$group

relatedness$stream[relatedness$stream==1111]="bear"
relatedness$stream[relatedness$stream==2222]="cascade"
relatedness$stream[relatedness$stream==3333]="canyon"
relatedness$stream[relatedness$stream==4444]="paradise"
relatedness$stream[relatedness$stream==5555]="zigzag"

#new dataframe subset to only include within stream relatedness estimates
instream<-relatedness[relatedness$stream=="bear" | relatedness$stream=="cascade" | relatedness$stream=="canyon"
                      | relatedness$stream=="paradise" | relatedness$stream=="zigzag",]

str(instream)
nrow(instream) #20143 pairwise comparisons
unique(instream$stream) #check to make sure only within-stream comparisons included

# need to make IDs match dispersal dataframe
instream<-separate(instream, ind1.id, c("group2", "ind1"), sep="_", remove=FALSE)
instream<-separate(instream, ind2.id, c("group3", "ind2"), sep="_", remove=FALSE)

#get rid of extra columns
instream$group2<-NULL
instream$group3<-NULL
instream$ind1.id<-NULL
instream$ind2.id<-NULL

# add location data for each individual from dispersal dataframe to figure out how far apart individuals are
# bring in dispersal data
dispdata=read.table("Dispersal.csv", sep=",", header=T)
str(dispdata)

#need to make ids character to match with instream df
dispdata$id<-as.character(dispdata$id)
str(dispdata)
str(instream)

#add final location data for each individual
instream$ind1.loc <- dispdata$final_location[match(instream$ind1, dispdata$id)]
instream$ind2.loc <- dispdata$final_location[match(instream$ind2, dispdata$id)]
str(instream)
tail(instream)

# calculate distance apart between individuals
instream$dist.apart<-abs(instream$ind1.loc-instream$ind2.loc)

#new dataframe with individuals 50m apart or less
instream.50<-instream[instream$dist.apart < 51,]
instream.50
tail(instream.50)

# get rid of rows with NAs
instream.50<-na.omit(instream.50)
nrow(instream.50)
table(instream.50$dist.apart)

# count number of individuals within 50m of focal individual
ind1.count<-as.data.frame(table(instream.50$ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"
ind2.count<-as.data.frame(table(instream.50$ind2))
ind2.count

#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

setdiff(ind1.count$Id, ind2.count$Id) # 40 individuals in ind1 but not in ind2
setdiff(ind2.count$Id, ind1.count$Id) # 38 individuals in ind2 but not in ind1


#merge count dataframes to make mastercount dataframe
mastercount<-merge(ind1.count, ind2.count, 
                   by = "Id",
                   all = TRUE)
nrow(mastercount)

#replace NAs with zero
mastercount$Ind1.Count[is.na(mastercount$Ind1.Count)] = 0
mastercount$Ind2.Count[is.na(mastercount$Ind2.Count)] = 0

# sum counts to get total number of individuals within 50m for each focal individual
mastercount$total.count<-mastercount$Ind1.Count+mastercount$Ind2.Count
max(mastercount$total.count) #30
min(mastercount$total.count) #1

##new dataframe with individuals 50m apartt or less AND relatedness higher than .132
instream.related<-instream.50[instream.50$trioml > 0.132,]
instream.related
tail(instream.related)
nrow(instream.related)
str(instream.related)

# count number of individuals that are related
ind1.count<-as.data.frame(table(instream.related$ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"
ind2.count<-as.data.frame(table(instream.related$ind2))
ind2.count

#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

#merge count dataframes to make mastercount dataframe
relatedcount<-merge(ind1.count, ind2.count, 
                    by = "Id",
                    all = TRUE)
nrow(relatedcount)

#replace NAs with zero
relatedcount$Ind1.Count[is.na(relatedcount$Ind1.Count)] = 0
relatedcount$Ind2.Count[is.na(relatedcount$Ind2.Count)] = 0

# sum counts to get total number of individuals within 50m for each focal individual
relatedcount$related<-relatedcount$Ind1.Count+relatedcount$Ind2.Count
max(relatedcount$related) #10
min(relatedcount$related) #1

#now merged related counts with mastercounts
all.counts<-merge(mastercount, relatedcount, 
                  by = "Id",
                  all = TRUE)

#replace NAs with zero
all.counts$related[is.na(all.counts$related)] = 0
all.counts

#calculate proportion relatives
all.counts$prop.related<-all.counts$related/all.counts$total.count
all.counts
max(all.counts$prop.related) #1
min(all.counts$prop.related) #0
median(all.counts$prop.related)# 0.111

#get rid of unnecessary columns
all.counts$Ind1.Count.x<-NULL
all.counts$Ind2.Count.x<-NULL
all.counts$Ind1.Count.y<-NULL
all.counts$Ind2.Count.y<-NULL

# add dispersal data to all.counts df
str(dispdata) # add final_location,absolute_distance, stream, reach
all.counts$stream <- dispdata$stream[match(all.counts$Id, dispdata$id)]
all.counts$final_location <- dispdata$final_location[match(all.counts$Id, dispdata$id)]
all.counts$absolute_distance <- dispdata$absolute_distance[match(all.counts$Id, dispdata$id)]
all.counts$reach <- dispdata$reach[match(all.counts$Id, dispdata$id)]

# add column designating disperser status
all.counts$disperser10<- "no"
all.counts$disperser10[all.counts$absolute_distance >=10]<- "yes"
str(all.counts)
unique(all.counts$stream)
unique(all.counts$reach)
unique(all.counts$disperser10)
table(all.counts$disperser10)

# need to exclude GBO4 and GZY1-7-8 because they dispersed between stream reaches
all.counts<-all.counts[all.counts$Id!="GBO4",]
all.counts<-all.counts[all.counts$Id!="GZY1-7-8",]
nrow(all.counts)

# test for correlation between salamander density and proportion related
ranges<-cbind(all.counts$total.count, all.counts$prop.related)
cor(ranges, use="pairwise.complete.obs", method="pearson")
corr.test(ranges)


#################################################################################
#### glmm testing for effects of dispersal on inbreeding risk

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# logistic regression with binomial distribution, dispersal treated categorically
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod<-glmer(cbind(related, total.count)~reach+disperser10+reach*disperser10 +(1|stream),na.action = na.exclude,data=all.counts, family="binomial")

summary(mod)
print(mod)
fixef(mod) #get fixed effects estimates
ranef(mod) #get random effects estimates

#function from Ben Bolker to test for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(mod)


# post-hoc tukey tests for dispersal x reach groups
em<-emmeans(mod, specs=pairwise~reach:disperser10, adjust="tukey")
em.response<-summary(em, type="response")


#save estimates for plotting
em.estimates<-as.data.frame(em.response$emmeans)
em.estimates$upper.se<-em.estimates$prob + em.estimates$SE
em.estimates$lower.se<-em.estimates$prob - em.estimates$SE

#rename column 
names(em.estimates)[names(em.estimates)=="disperser10"] <-"Dispersal"
em.estimates

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### plot estimated marginal means
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#greyscale plot
a<-ggplot(em.estimates, aes(x=reach, y=prob))+
  geom_point(size=4, aes(color=reach, shape=Dispersal),position=position_dodge(width=0.4))+
  geom_errorbar(aes(ymin=lower.se, ymax=upper.se, color=reach, group=Dispersal),width = .1,position=position_dodge(width=0.4))+
  scale_shape_manual(values=c(17,16), guide=guide_legend(override.aes=list(shape=c(2,1), size=4)))+
  scale_color_manual(values=c('#000000', '#999999'), guide="none")+
  scale_fill_discrete()+
  labs(x="Reach", y= "Proportion relatives within 50m (EMM \\u00B1 SE)")+
  theme_classic()+
  theme(axis.title=element_text(size=9),axis.text.y=element_text(size=8), axis.text.x=element_text(size=8),
        legend.position="top", legend.text=element_text(size=8),legend.title=element_text(size=8))
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# logistic regression with binomial distribution, dispersal treated continuously
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mod2<-glmer(cbind(related, total.count)~reach*absolute_distance +(1|stream),na.action = na.exclude,data=all.counts, family="binomial")
summary(mod2)

em<-as.data.frame(effect("reach*absolute_distance", mod2, xlevels=list(absolute_distance=seq(0,351,1))))
str(em)

#greyscale
b<-ggplot(em, aes(x=absolute_distance, y=fit, fill=reach, color=reach))+
  geom_line(size=1)+
  geom_ribbon_pattern(aes(ymin=lower, ymax=upper, pattern=reach), 
                      alpha=0.3, size=0.1, 
                      pattern_fill="black",
                      pattern_alpha=0.4)+
  scale_x_continuous(name="Dispersal distance (m)")+
  scale_y_continuous(name="Proportion relatives within 50m")+
  scale_color_manual(values=c("#000000","#999999"))+
  scale_fill_manual(values=c("#000000","#999999"))+
  theme_bw()+
  theme(legend.position="top",legend.text=element_text(size=8),
        legend.title=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.line=element_line(colour="black"),
        axis.title=element_text(size=9),axis.text.y=element_text(size=8), axis.text.x=element_text(size=8))


##################################################################################################
#calculate number of conspecifics within 50m each individual -- test of intraspecific competition

alldata<-read.table("CaptureData.csv", sep=",", header=T)

#subset data by stream
bear<-alldata[which(alldata$Stream=="Bear"),]
cascade<-alldata[which(alldata$Stream=="Cascade"),]
canyon<-alldata[which(alldata$Stream=="Canyon"),]
paradise<-alldata[which(alldata$Stream=="Paradise"),]
zigzag<-alldata[which(alldata$Stream=="Zig Zag"),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate pairwise distances apart for zigzag
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zigzag.pair<-data.frame()
for (i in 1:nrow(zigzag)){
  for (j in 1:nrow(zigzag)){
    output=c(zigzag$ID[i], zigzag$ID[j], zigzag$EndLoc[i], zigzag$EndLoc[j])
    zigzag.pair=rbind(zigzag.pair,output)
  }
}
zigzag.pair
colnames(zigzag.pair)<-c("Ind1", "Ind2", "EndLoc1", "EndLoc2")
nrow(zigzag.pair)#306916, =554*554. so need to get rid of duplicate comparisons

#get rid of rows that compare individuals to themselves
zigzag.pair<-zigzag.pair[which(zigzag.pair$Ind1 != zigzag.pair$Ind2),]
nrow(zigzag.pair) #306362, = 306916-554, good

#find distance between individuals
zigzag.pair$EndLoc1<-as.numeric(zigzag.pair$EndLoc1)
zigzag.pair$EndLoc2<-as.numeric(zigzag.pair$EndLoc2)
zigzag.pair$diff<-abs(zigzag.pair$EndLoc1-zigzag.pair$EndLoc2)

#delete end locs
zigzag.pair$EndLoc1<-NULL
zigzag.pair$EndLoc2<-NULL
zigzag.pair

#sort to identify duplicate comparisons
zigzag.pair<-data.frame(t(apply(zigzag.pair,1,sort)))
zigzag.pair

zigzag.pair<-unique(zigzag.pair)
nrow(zigzag.pair)# 153181, yay!

#rename columns
colnames(zigzag.pair)<-c("DistApart", "Ind1", "Ind2")

#slim down to only individuals that are 50m apart or less
zigzag.pair$DistApart<-as.numeric(zigzag.pair$DistApart)
zigzag.50<-zigzag.pair[which(zigzag.pair$DistApart < 51),]
max(zigzag.50$DistApart)# good
table(zigzag.50$DistApart)

# count number of individuals within 50m of focal individual
ind1.count<-as.data.frame(table(zigzag.50$Ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"

ind2.count<-as.data.frame(table(zigzag.50$Ind2))
ind2.count

#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

#merge count dataframes to make mastercount dataframe
mastercount<-merge(ind1.count, ind2.count, 
                   by = "Id",
                   all = TRUE)
nrow(mastercount)

#replace NAs with zero
mastercount$Ind1.Count[is.na(mastercount$Ind1.Count)] = 0
mastercount$Ind2.Count[is.na(mastercount$Ind2.Count)] = 0

#sum counts to get total number of individuals within 50m for each focal individual
mastercount$total.count<-mastercount$Ind1.Count+mastercount$Ind2.Count
max(mastercount$total.count) #117
min(mastercount$total.count) #1

#add # of inds within 50m  to Dispersal dataframe
dispersal=read.table("Dispersal.csv", sep=",", header=T)
dispersal
names(dispersal)
names(mastercount)

dispersal$NumIndsWithin50m.zigzag<- mastercount$total.count[match(dispersal$id, mastercount$Id)]
dispersal

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate pairwise distances apart for canyon
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
canyon.pair<-data.frame()
for (i in 1:nrow(canyon)){
  for (j in 1:nrow(canyon)){
    output=c(canyon$ID[i], canyon$ID[j], canyon$EndLoc[i], canyon$EndLoc[j])
    canyon.pair=rbind(canyon.pair,output)
  }
}
canyon.pair
colnames(canyon.pair)<-c("Ind1", "Ind2", "EndLoc1", "EndLoc2")
nrow(canyon.pair)#152100, =390*390. so need to get rid of duplicate comparisons

#get rid of rows that compare individuals to themselves
canyon.pair<-canyon.pair[which(canyon.pair$Ind1 != canyon.pair$Ind2),]
nrow(canyon.pair) #151710, = 152100-390, good

#find distance between individuals
canyon.pair$EndLoc1<-as.numeric(canyon.pair$EndLoc1)
canyon.pair$EndLoc2<-as.numeric(canyon.pair$EndLoc2)
canyon.pair$diff<-abs(canyon.pair$EndLoc1-canyon.pair$EndLoc2)

#delete end locs
canyon.pair$EndLoc1<-NULL
canyon.pair$EndLoc2<-NULL
canyon.pair

#sort to identify duplicate comparisons
canyon.pair<-data.frame(t(apply(canyon.pair,1,sort)))
canyon.pair

canyon.pair<-unique(canyon.pair)
nrow(canyon.pair)# 75855, yay!

#rename columns 
colnames(canyon.pair)<-c("DistApart", "Ind1", "Ind2")

#slim down to only inds that are 50m apart or less
canyon.pair$DistApart<-as.numeric(canyon.pair$DistApart)
canyon.50<-canyon.pair[which(canyon.pair$DistApart < 51),]
max(canyon.50$DistApart)# good
table(canyon.50$DistApart)

# count number of individuals within 50m of focal individual
ind1.count<-as.data.frame(table(canyon.50$Ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"

ind2.count<-as.data.frame(table(canyon.50$Ind2))
ind2.count
#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

#merge count dataframes to make mastercount dataframe
mastercount<-merge(ind1.count, ind2.count, 
                   by = "Id",
                   all = TRUE)
nrow(mastercount)

#replace NAs with zero
mastercount$Ind1.Count[is.na(mastercount$Ind1.Count)] = 0
mastercount$Ind2.Count[is.na(mastercount$Ind2.Count)] = 0

# sum counts to get total number of individuals within 50m for each focal individual
mastercount$total.count<-mastercount$Ind1.Count+mastercount$Ind2.Count
max(mastercount$total.count) #100
min(mastercount$total.count) #21

#add # of inds within 50m  to Dispersal dataframe
setwd("C:/Users/addis/OneDrive/manuscripts/genomics/to submit")
dispersal=read.table("Dispersal.csv", sep=",", header=T)
dispersal
names(dispersal)
names(mastercount)

dispersal$NumIndsWithin50m.canyon<- mastercount$total.count[match(dispersal$id, mastercount$Id)]
dispersal


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate pairwise distances apart for paradise
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
paradise.pair<-data.frame()
for (i in 1:nrow(paradise)){
  for (j in 1:nrow(paradise)){
    output=c(paradise$ID[i], paradise$ID[j], paradise$EndLoc[i], paradise$EndLoc[j])
    paradise.pair=rbind(paradise.pair,output)
  }
}
paradise.pair
colnames(paradise.pair)<-c("Ind1", "Ind2", "EndLoc1", "EndLoc2")
nrow(paradise.pair)#753424, =868*868. so need to get rid of duplicate comparisons

#get rid of rows that compare individuals to themselves
paradise.pair<-paradise.pair[which(paradise.pair$Ind1 != paradise.pair$Ind2),]
nrow(paradise.pair) #752556, = 753424-868, good

#find distance between individuals
paradise.pair$EndLoc1<-as.numeric(paradise.pair$EndLoc1)
paradise.pair$EndLoc2<-as.numeric(paradise.pair$EndLoc2)
paradise.pair$diff<-abs(paradise.pair$EndLoc1-paradise.pair$EndLoc2)

#delete end locs
paradise.pair$EndLoc1<-NULL
paradise.pair$EndLoc2<-NULL
paradise.pair

#sort to identify duplicate comparisons
paradise.pair<-data.frame(t(apply(paradise.pair,1,sort)))
paradise.pair

paradise.pair<-unique(paradise.pair)
nrow(paradise.pair)# 376278, yay!

#rename columns 
colnames(paradise.pair)<-c("DistApart", "Ind1", "Ind2")

#slim down to only inds that are 50m apart or less
paradise.pair$DistApart<-as.numeric(paradise.pair$DistApart)
paradise.50<-paradise.pair[which(paradise.pair$DistApart < 51),]
max(paradise.50$DistApart)# good
table(paradise.50$DistApart)

# count number of individuals within 50m of focal individual
ind1.count<-as.data.frame(table(paradise.50$Ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"

ind2.count<-as.data.frame(table(paradise.50$Ind2))
ind2.count

#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

#merge count dataframes to make mastercount dataframe
mastercount<-merge(ind1.count, ind2.count, 
                   by = "Id",
                   all = TRUE)
nrow(mastercount)

#replace NAs with zero
mastercount$Ind1.Count[is.na(mastercount$Ind1.Count)] = 0
mastercount$Ind2.Count[is.na(mastercount$Ind2.Count)] = 0

# sum counts to get total number of individuals within 50m for each focal individual
mastercount$total.count<-mastercount$Ind1.Count+mastercount$Ind2.Count
max(mastercount$total.count) #158
min(mastercount$total.count) #22

#add # of inds within 50m  to Dispersal dataframe
setwd("C:/Users/addis/OneDrive/manuscripts/genomics/to submit")
dispersal=read.table("Dispersal.csv", sep=",", header=T)
dispersal
names(dispersal)
names(mastercount)

dispersal$NumIndsWithin50m.paradise<- mastercount$total.count[match(dispersal$id, mastercount$Id)]
dispersal

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate pairwise distances apart for bear
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bear.pair<-data.frame()
for (i in 1:nrow(bear)){
  for (j in 1:nrow(bear)){
    output=c(bear$ID[i], bear$ID[j], bear$EndLoc[i], bear$EndLoc[j])
    bear.pair=rbind(bear.pair,output)
  }
}
bear.pair
colnames(bear.pair)<-c("Ind1", "Ind2", "EndLoc1", "EndLoc2")
nrow(bear.pair)#863041, =929*929. so need to get rid of duplicate comparisons

#get rid of rows that compare individuals to themselves
bear.pair<-bear.pair[which(bear.pair$Ind1 != bear.pair$Ind2),]
nrow(bear.pair) #862110, = 863041-929, good

#find distance between individuals
bear.pair$EndLoc1<-as.numeric(bear.pair$EndLoc1)
bear.pair$EndLoc2<-as.numeric(bear.pair$EndLoc2)
bear.pair$diff<-abs(bear.pair$EndLoc1-bear.pair$EndLoc2)

#delete end locs
bear.pair$EndLoc1<-NULL
bear.pair$EndLoc2<-NULL
bear.pair

#sort to identify duplicate comparisons
bear.pair<-data.frame(t(apply(bear.pair,1,sort)))
bear.pair

bear.pair<-unique(bear.pair)
nrow(bear.pair)# 431055, yay!

#rename columns 
colnames(bear.pair)<-c("DistApart", "Ind1", "Ind2")

#slim down to only inds that are 50m apart or less
bear.pair$DistApart<-as.numeric(bear.pair$DistApart)
bear.50<-bear.pair[which(bear.pair$DistApart < 51),]
max(bear.50$DistApart)# good
table(bear.50$DistApart)

# count number of individuals within 50m of focal individual
ind1.count<-as.data.frame(table(bear.50$Ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"

ind2.count<-as.data.frame(table(bear.50$Ind2))
ind2.count

#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

#merge count dataframes to make mastercount dataframe
mastercount<-merge(ind1.count, ind2.count, 
                   by = "Id",
                   all = TRUE)
nrow(mastercount)

#replace NAs with zero
mastercount$Ind1.Count[is.na(mastercount$Ind1.Count)] = 0
mastercount$Ind2.Count[is.na(mastercount$Ind2.Count)] = 0

# sum counts to get total number of individuals within 50m for each focal individual
mastercount$total.count<-mastercount$Ind1.Count+mastercount$Ind2.Count
max(mastercount$total.count) #221
min(mastercount$total.count) #17

#add # of inds within 50m  to Dispersal dataframe
setwd("C:/Users/addis/OneDrive/manuscripts/genomics/to submit")
dispersal=read.table("Dispersal.csv", sep=",", header=T)
dispersal
names(dispersal)
names(mastercount)

dispersal$NumIndsWithin50m.bear<- mastercount$total.count[match(dispersal$id, mastercount$Id)]
dispersal

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate pairwise distances apart for cascade
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cas.pair<-data.frame()
for (i in 1:nrow(cascade)){
  for (j in 1:nrow(cascade)){
    output=c(cascade$ID[i], cascade$ID[j], cascade$EndLoc[i], cascade$EndLoc[j])
    cas.pair=rbind(cas.pair,output)
  }
}
cas.pair
colnames(cas.pair)<-c("Ind1", "Ind2", "EndLoc1", "EndLoc2")
nrow(cas.pair)#57121, =239*239. so need to get rid of duplicate comparisons

#get rid of rows that compare individuals to themselves
cas.pair<-cas.pair[which(cas.pair$Ind1 != cas.pair$Ind2),]
nrow(cas.pair) #56882, = 57121-239, good

#find distance between individuals
cas.pair$EndLoc1<-as.numeric(cas.pair$EndLoc1)
cas.pair$EndLoc2<-as.numeric(cas.pair$EndLoc2)
cas.pair$diff<-abs(cas.pair$EndLoc1-cas.pair$EndLoc2)

#delete end locs
cas.pair$EndLoc1<-NULL
cas.pair$EndLoc2<-NULL
cas.pair

#sort to identify duplicate comparisions
cas.pair<-data.frame(t(apply(cas.pair,1,sort)))
cas.pair

cas.pair<-unique(cas.pair)
nrow(cas.pair)# 28441, yay!

#rename columns 
colnames(cas.pair)<-c("DistApart", "Ind1", "Ind2")

#slim down to only inds that are 50m apart or less
cas.pair$DistApart<-as.numeric(cas.pair$DistApart)
cas.50<-cas.pair[which(cas.pair$DistApart < 51),]
max(cas.50$DistApart)# good
table(cas.50$DistApart)

# count number of individuals within 50m of focal individual
ind1.count<-as.data.frame(table(cas.50$Ind1))
ind1.count
str(ind1.count)
names(ind1.count)

#rename columns
names(ind1.count)[names(ind1.count)=="Var1"] <-"Id"
names(ind1.count)[names(ind1.count)=="Freq"] <-"Ind1.Count"

ind2.count<-as.data.frame(table(cas.50$Ind2))
ind2.count

#rename columns
names(ind2.count)[names(ind2.count)=="Var1"] <-"Id"
names(ind2.count)[names(ind2.count)=="Freq"] <-"Ind2.Count"

#merge count dataframes to make mastercount dataframe
mastercount<-merge(ind1.count, ind2.count, 
                   by = "Id",
                   all = TRUE)
nrow(mastercount)

#replace NAs with zero
mastercount$Ind1.Count[is.na(mastercount$Ind1.Count)] = 0
mastercount$Ind2.Count[is.na(mastercount$Ind2.Count)] = 0

# sum counts to get total number of individuals within 50m for each focal individual
mastercount$total.count<-mastercount$Ind1.Count+mastercount$Ind2.Count
max(mastercount$total.count) #61
min(mastercount$total.count) #3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add # of inds within 50m  to Dispersal dataframe

dispersal=read.table("Dispersal.csv", sep=",", header=T)
dispersal
names(dispersal)
names(mastercount)

dispersal$NumIndsWithin50m.cascade<- mastercount$total.count[match(dispersal$id, mastercount$Id)]
dispersal


#replace NAs with 0
dispersal$NumIndsWithin50m.bear[is.na(dispersal$NumIndsWithin50m.bear)]<-0
dispersal$NumIndsWithin50m.canyon[is.na(dispersal$NumIndsWithin50m.canyon)]<-0
dispersal$NumIndsWithin50m.cascade[is.na(dispersal$NumIndsWithin50m.cascade)]<-0
dispersal$NumIndsWithin50m.paradise[is.na(dispersal$NumIndsWithin50m.paradise)]<-0
dispersal$NumIndsWithin50m.zigzag[is.na(dispersal$NumIndsWithin50m.zigzag)]<-0

#create one column
dispersal$NumIndsWithin50m<-dispersal$NumIndsWithin50m.bear+dispersal$NumIndsWithin50m.canyon+
  dispersal$NumIndsWithin50m.cascade+dispersal$NumIndsWithin50m.paradise+dispersal$NumIndsWithin50m.zigzag

#delete extra columns
dispersal$NumIndsWithin50m.bear<-NULL
dispersal$NumIndsWithin50m.canyon<-NULL
dispersal$NumIndsWithin50m.cascade<-NULL
dispersal$NumIndsWithin50m.paradise<-NULL
dispersal$NumIndsWithin50m.zigzag<-NULL
dispersal

#calculate mean number of conspecifics within 50m
dispersal<-dispersal[!is.na(dispersal$NumIndsWithin50m),]# get rid of rows with NA
min(dispersal$NumIndsWithin50m)#6
max(dispersal$NumIndsWithin50m)#221
median(dispersal$NumIndsWithin50m)#89

#add dispersal key
dispersal$disp<-"NULL"
dispersal$disp[dispersal$absolute_distance>=10]<-"yes"
dispersal$disp[dispersal$absolute_distance<10]<-"no"
dispersal[,c("id", "absolute_distance", "disp")]
table(dispersal$disp)

# test for correlation between salamander density and dispersal distance
#exclude distance > 500m

dispersal2<-dispersal[which(dispersal$absolute_distance<501),]
dispersal2$log.dist<-log(dispersal2$absolute_distance+1)
str(dispersal2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#linear mixed model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LMM<-lme(NumIndsWithin50m~disp*reach, random = ~1|stream, 
          data=dispersal2, na.action=na.exclude)

summary(LMM)

#get estimated marginal means
em.lmm<-emmeans(LMM, specs=pairwise~reach:disp, adjust="tukey")
em.response<-summary(em.lmm, type="response")
em.response

#save estimates for plotting
em.estimates<-as.data.frame(em.response$emmeans)
em.estimates$upper.se<-em.estimates$emmean + em.estimates$SE
em.estimates$lower.se<-em.estimates$emmean - em.estimates$SE

#rename column 
names(em.estimates)[names(em.estimates)=="disp"] <-"Dispersal"
em.estimates


#greyscale
c<-ggplot(em.estimates, aes(x=reach, y=emmean))+
  geom_point(size=4, aes(color=reach, shape=Dispersal),position=position_dodge(width=0.4))+
  geom_errorbar(aes(ymin=lower.se, ymax=upper.se, color=reach, group=Dispersal),width = .1,position=position_dodge(width=0.4))+
  scale_shape_manual(values=c(17,16), guide=guide_legend(override.aes=list(shape=c(2,1), size=4)))+
  scale_color_manual(values=c('#000000', '#999999'), guide="none")+
  scale_fill_discrete()+
  labs(x="Reach", y= "Number of conspecifics within 50m (EMM \\u00B1 SE)")+
  theme_classic()+
  theme(axis.title=element_text(size=9),axis.text.y=element_text(size=8), axis.text.x=element_text(size=8),
        legend.position="top",legend.text=element_text(size=8), legend.title=element_text(size=8))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create multiplanel plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggarrange(a,b,c, nrow=2, ncol=3)



