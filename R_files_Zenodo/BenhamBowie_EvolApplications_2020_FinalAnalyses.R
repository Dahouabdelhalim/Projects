#Primary code used to perform analyses and generate most figures for:
#Benham, P.M. & Bowie, R.C.K. (2020). The influence of spatially heterogenous anthropogenic change #on bill #size evolution in a coastal songbird. Evolutionary Applications.

#Pathways to datasets and directories will need to be modified prior to running on a different computer.

#last edited: 22 September 2020 by Phred M. Benham.

#packages needed to run code
library(ggplot2)
library(gridExtra)
library(car)
library(AICcmodavg)
library(grDevices)
library(MASS)
library(lme4)
library(stargazer)
library(dotwhisker)
library(arsenal)
library(RColorBrewer)

#####################################################################################################
#############################Spatial variation in bill morphology####################################

#initial filtering to remove adults, retain only resident/breeding populations of the 4 CA subspecies.
SAVSbills<-read.csv("~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/FinalData_Code_BenhamBowie_2020_EvolApplications/MorphologicalData_BenhamBowie_EvolApplications_2020.csv",header=TRUE)

SAVSbills<-SAVSbills[!(SAVSbills$Subspecies== "alaudinus" & SAVSbills$tarsus > 70.0),]

SAVSbills2<-subset(SAVSbills, Country == "USA" & Subspecies != "rostratus" & Subspecies != "migrant" & age == "adult" & sex != "undetermined")

SAVS.CA.final<-droplevels(SAVSbills2)

#Make table with final numbers for each subspecies, males, females, etc. 
SAVS.CA.final$TimeSex<-paste(SAVS.CA.final$Time,SAVS.CA.final$sex,sep="_")
SumTable<-tableby(Subspecies~TimeSex, data = SAVS.CA.final)
print(summary(SumTable))

traits.z<-scale(SAVS.CA.final[,c(16:48)],center=TRUE,scale=TRUE)
SAVS.CA.final2<-as.data.frame(cbind(SAVS.CA.final[,c(1:15)],traits.z))

#Perform AIC analyses on several mixed effects models where temporal variation
#and subspecific variation is controlled for (i.e. are random effects)
Cand.models<-list()
Cand.models[[1]]<-lmer(SurfaceArea~sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[2]]<-lmer(SurfaceArea~Tmax_year + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[3]]<-lmer(SurfaceArea~VPDmax_max + sex + tarsus  + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[4]]<-lmer(SurfaceArea~Prec_year + sex + tarsus  +  (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[5]]<-lmer(SurfaceArea~Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[6]]<-lmer(SurfaceArea~Tmax_year*VPDmax_max + sex + tarsus  + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[7]]<-lmer(SurfaceArea~Tmax_year*Prec_year + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[8]]<-lmer(SurfaceArea~Tmax_year*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[9]]<-lmer(SurfaceArea~VPDmax_max*Prec_year + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[10]]<-lmer(SurfaceArea~VPDmax_max*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[11]]<-lmer(SurfaceArea~Prec_year*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[12]]<-lmer(SurfaceArea~Tmax_year*Prec_year*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[13]]<-lmer(SurfaceArea~Tmax_year*VPDmax_max*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[14]]<-lmer(SurfaceArea~Prec_year*VPDmax_max*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[15]]<-lmer(SurfaceArea~Prec_year*VPDmax_max*Tmax_year + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[16]]<-lmer(SurfaceArea~Prec_year*VPDmax_max*Tmax_year*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
Cand.models[[17]]<-lmer(SurfaceArea~Ecoregion*Habitat + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)

Modnames<-c("base","Tmax","VPDmax","Prec","Salinity","Tmax:VPDmax","Tmax:Prec","Tmax:Salinity","VPDmax:Prec","VPDmax:Salinity","Prec:Salinity","Tmax:Prec:Salinity","Tmax:VPDmax:Salinity","Prec:VPDmax:Salinity","Prec:VPDmax:Tmax","Prec:VPDmax:Tmax:Salinity","Ecoregion")

aic.table<-aictab(cand.set=Cand.models, modnames=Modnames, sort=TRUE)
print(aic.table, digits=3, LL=TRUE)

write.table(aic.table, file="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/Clim.mixedmodel.AIC.tsv", quote=FALSE, sep='\\t',col.names=NA)

#output final mixed effects model results from top model selected by AIC results. These results are reported in Table 2 (AIC results) and Table 3 (LFMM results) of the main paper. 
lmm1<-lmer(SurfaceArea~Prec_year*VPDmax_max*Tmax_year*Mean_Salinity + sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
#lmm2<-lmer(SurfaceArea~ sex + tarsus + (1|Subspecies) + (1|decade), data=SAVS.CA.final2, REML=FALSE)
print(summary(lmm1))
print(Anova(lmm1))

##################################################################################################################

#correct for body size for each bill character using tarsus length
savs.sa<-SAVS.CA.final[,c("SurfaceArea", "tarsus")]
SAVS.new<-SAVS.CA.final[complete.cases(savs.sa),]

tarsusSA.lm<-lm(SAVS.new$SurfaceArea~SAVS.new$tarsus)
Blengthtarsus.lm<-lm(SAVS.new$BillLength~SAVS.new$tarsus)
Bwidthtarsus.lm<-lm(SAVS.new$BillWidth~SAVS.new$tarsus)
Bdepthtarsus.lm<-lm(SAVS.new$BillDepth~SAVS.new$tarsus)

SAtarCorr<-residuals(tarsusSA.lm)
BlengthCorr<-residuals(Blengthtarsus.lm)
BwidthCorr<-residuals(Bwidthtarsus.lm)
BdepthCorr<-residuals(Bdepthtarsus.lm)

SAVS.temp<-as.data.frame(cbind(SAVS.new[,-49],SAtarCorr,BlengthCorr,BwidthCorr,BdepthCorr))
traits2.z<-scale(SAVS.temp[,c(16:52)],center=TRUE,scale=TRUE)
SAVS.CA.final3<-as.data.frame(cbind(SAVS.temp[,c(1:15)], traits2.z))

#PCA on statewide and temporal climate data from the PRISM dataset.
SAVS.CA.final3<-na.omit(SAVS.CA.final3)

Clim.pca<-prcomp(SAVS.CA.final3[,c(17,18,20,21)])
print(summary(Clim.pca))
print(Clim.pca)
PC1<-Clim.pca$x[,1]
PC2<-Clim.pca$x[,2]
SAVS.CA.final4<-as.data.frame(cbind(SAVS.CA.final3,PC1,PC2))

#Plot of PC1 and PC2 of climate variables in model, color by subspecies. This creates Figure 2 of main paper. 
cbPalette <- c("#009E73","#0072B2", "#CC79A7", "#D55E00")
#
nevadensis.col = rgb(136/255,204/255,238/255)
brooksi.col = rgb(51/255,34/255,136/255)
alaudinus.col = rgb(17/255,119/255,51/255)
belding.col = rgb(170/255,68/255,153/255)

CA.bill.clim.plot<-ggplot(SAVS.CA.final4, aes(x=PC1,y=PC2, group=Subspecies, fill=Subspecies))  + theme_bw() + geom_point(aes(shape=Time, size=Time),color='black') + scale_shape_manual(values=c(21,24)) + scale_size_manual(values=c(3,6)) + labs(x="PC1 (55.8%)", y="PC2 (34.2%)") + scale_fill_manual(values=c(alaudinus.col, belding.col,brooksi.col,nevadensis.col)) + theme(axis.text=element_text(size=16),axis.title=element_text(size=18), legend.title=element_text(size=16, color="black"), legend.text=element_text(color="black", size=14))

print(CA.bill.clim.plot)

ggsave("California_Subspecies_ClimatePCplot.pdf", plot=CA.bill.clim.plot, width=9, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#####################################################################################################
#############################Temporal variation in bill morphology###################################

#plot relationship between subspecies and Time. This code generates Figure 3 in the main paper.
Subspp.var.Time.plot<-ggplot(SAVS.temp, aes(x=year, y=SAtarCorr, group=Subspecies, color=Subspecies)) + facet_wrap(~Subspecies) + theme_bw()+ scale_x_continuous(limits=c(1850,2020), expand=c(0,0), breaks=seq(1880,2000,40)) + geom_point(size=1.5) + geom_smooth(method=lm) + scale_color_manual(values=c(alaudinus.col, belding.col,brooksi.col,nevadensis.col)) + labs(y="Bill surface area (size corrected)", x="Year") + theme(axis.text=element_text(size=16),axis.title=element_text(size=18), legend.title=element_text(size=16, color="black"), legend.text=element_text(color="black", size=14))

print(Subspp.var.Time.plot)

ggsave("Subspp.var.Time.plot.pdf", plot=Subspp.var.Time.plot, width=9, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#BillChar_Time.plot is Supplemental Figure S2.
Subspp.blength.Time.plot<-ggplot(SAVS.temp, aes(x=year, y=BlengthCorr, group=Subspecies, color=Subspecies))  + theme_bw() + geom_point(size=1) + geom_smooth(method=lm) + scale_color_manual(values=c(alaudinus.col, belding.col,brooksi.col,nevadensis.col)) + labs(y="Bill length (size corrected)", x="Year") + theme(axis.text=element_text(size=12),axis.title=element_text(size=14), legend.title=element_text(size=12, color="black"), legend.text=element_text(color="black", size=12))

Subspp.bwidth.Time.plot<-ggplot(SAVS.temp, aes(x=year, y=BwidthCorr, group=Subspecies, color=Subspecies))  + theme_bw() + geom_point(size=1) + geom_smooth(method=lm) + scale_color_manual(values=c(alaudinus.col, belding.col,brooksi.col,nevadensis.col)) + labs(y="Bill width (size corrected)", x="Year") + theme(axis.text=element_text(size=12),axis.title=element_text(size=14), legend.title=element_text(size=12, color="black"), legend.text=element_text(color="black", size=12))

Subspp.bdepth.Time.plot<-ggplot(SAVS.temp, aes(x=year, y=BdepthCorr, group=Subspecies, color=Subspecies))  + theme_bw() + geom_point(size=1) + geom_smooth(method=lm) + scale_color_manual(values=c(alaudinus.col, belding.col,brooksi.col,nevadensis.col)) + labs(y="Bill depth (size corrected)", x="Year") + theme(axis.text=element_text(size=12),axis.title=element_text(size=14), legend.title=element_text(size=12, color="black"), legend.text=element_text(color="black", size=12))

BillChar_Time.plot<-grid.arrange(Subspp.blength.Time.plot, Subspp.bwidth.Time.plot, Subspp.bdepth.Time.plot, ncol=1,nrow=3)

print(BillChar_Time.plot)

ggsave("BillChar_Time.plot.pdf", plot=BillChar_Time.plot, width=6, height=11, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")
############################################################################################################
#regression analyses comparing relationship between year and and different bill dimensions.
#across all four subspecies. Results can be found in Fig. 3 of main paper for bill surface area and
#in Supplemental Table S4 for other characters.

subspp.Time.aov<-aov(SAtarCorr~year*Subspecies, data=SAVS.temp)
print(summary(subspp.Time.aov))
bill.char<-c( "SAtarCorr", "BlengthCorr", "BwidthCorr", "BdepthCorr")
for (i in 1:4){
	print(bill.char[i])
	alaud.bill<-droplevels(subset(SAVS.temp,Subspecies=="alaudinus"))
	belding.bill<-droplevels(subset(SAVS.temp,Subspecies=="beldingi"))
	brooks.bill<-droplevels(subset(SAVS.temp,Subspecies=="brooksi"))
	nevadensis.bill<-droplevels(subset(SAVS.temp,Subspecies=="nevadensis"))
	
	alaud.lm<-lm(alaud.bill[,i+47]~alaud.bill$year, data=alaud.bill)
	belding.lm<-lm(belding.bill[,i+47]~belding.bill$year, data=belding.bill)
	brooksi.lm<-lm(brooks.bill[,i+47]~brooks.bill$year, data=brooks.bill)
	nevadensis.lm<-lm(nevadensis.bill[,i+47]~nevadensis.bill$year, data=nevadensis.bill)
	
	print(summary(alaud.lm))
	print(summary(belding.lm))
	print(summary(brooksi.lm))
	print(summary(nevadensis.lm))
	
	stargazer(alaud.lm, belding.lm, brooksi.lm,nevadensis.lm, type="text",style="all", digits=3, 					star.cutoffs=c(0.05,0.01,0.001),digit.separator="", out="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/Year.subspp.txt")
}

#####################################################################################################
#############################Temporal variation in bill morphology for alaudinus#####################
#subset for just alaudinus
SAVS.alaudinus<-subset(SAVS.temp, Subspecies == "alaudinus")

#This code creates Figure 4 of main paper.
BillSurfaceArea.plot<-ggplot(SAVS.alaudinus, aes(y=SAtarCorr, x=latitude, group=Time, color=Time)) + theme_bw() + geom_point(size=1) + stat_smooth(method="lm") + labs(title="Bill surface area", y="Bill surface area (residuals)", x="Latitude") + theme(axis.text=element_text(size=14),axis.title=element_text(size=16), legend.position="none",plot.title=element_text(size=16,face="bold",hjust=0.5)) 

BillLength.plot<-ggplot(SAVS.alaudinus, aes(y=BlengthCorr, x=latitude, group=Time, color=Time)) + theme_bw() + geom_point(size=1) + stat_smooth(method="lm") + labs(title="Bill length",y="Bill length (residuals)", x="Latitude") + theme(axis.text=element_text(size=14),axis.title=element_text(size=16), legend.position="none",plot.title=element_text(size=16,face="bold",hjust=0.5)) 

Billwidth.plot<-ggplot(SAVS.alaudinus, aes(y=BwidthCorr, x=latitude, group=Time, color=Time)) + theme_bw() + geom_point(size=1) + stat_smooth(method="lm") + labs(title="Bill width", y="Bill width (residuals)", x="Latitude") + theme(axis.text=element_text(size=14),axis.title=element_text(size=16), legend.position="none",plot.title=element_text(size=16,face="bold",hjust=0.5)) 

Billdepth.plot<-ggplot(SAVS.alaudinus, aes(y=BdepthCorr, x=latitude, group=Time, color=Time)) + theme_bw() + geom_point(size=1) + stat_smooth(method="lm") + labs(title="Bill depth", y="Bill depth (residuals)", x="Latitude",fill="Time period:") + theme(axis.text=element_text(size=14),axis.title=element_text(size=16), legend.position="bottom",plot.title=element_text(size=16,face="bold",hjust=0.5)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Billdepth.plot)


Billsplot<-grid.arrange(arrangeGrob(BillSurfaceArea.plot + theme(legend.position="none"),BillLength.plot + theme(legend.position="none"), Billwidth.plot + theme(legend.position="none"), Billdepth.plot + theme(legend.position="none"), ncol=2,nrow=2), mylegend, nrow=2,heights=c(10,1))


ggsave("Billcharacters_Timebylat.pdf", plot=Billsplot, width=9, height=8, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#ANCOVA analyses assessing influence of time period (historic vs. modern), latitude and time by latitude
#interactions on each bill character. Results reported in Table 4 of the main paper.
billSA.Time.lm<-aov(SAtarCorr~Time*latitude, data=SAVS.alaudinus)
billLength.Time.lm<-aov(BlengthCorr~Time*latitude, data=SAVS.alaudinus)
billWidth.Time.lm<-aov(BwidthCorr~Time*latitude, data=SAVS.alaudinus)
billDepth.Time.lm<-aov(BdepthCorr~Time*latitude, data=SAVS.alaudinus)

print(summary(billSA.Time.lm))
print(summary(billLength.Time.lm))
print(summary(billWidth.Time.lm))
print(summary(billDepth.Time.lm))

#####################################################################################################
##########################   Accounting for genetic drift in shaping temporal trends  ###############
#Code to generate Figure 5 of the main paper. Estimates of amount of drift generated using the Lande1976.py 
#python code. 
traits.z<-scale(SAVS.alaudinus[,c(15:51)],scale=TRUE, center=TRUE)
SAVS.alaudinus2<-as.data.frame(cbind(SAVS.alaudinus[,c(1:14)],traits.z))

mor.t<-t.test(SurfaceArea~Time, data=subset(SAVS.alaudinus2, region=="Morro" & sex == "male"))
sfbay.t<-t.test(SurfaceArea~Time, data=subset(SAVS.alaudinus2, region=="SFBay" & sex == "male"))
sanpab.t<-t.test(SurfaceArea~Time, data=subset(SAVS.alaudinus2, region=="SanPablo" & sex == "male"))
suisun.t<-t.test(SurfaceArea~Time, data=subset(SAVS.alaudinus2, region=="Suisun" & sex == "male"))
humbol.t<-t.test(SurfaceArea~Time, data=subset(SAVS.alaudinus2, region=="HumboldtBay" & sex == "male"))

id <-as.numeric(c(1:5))
region.names<-c("MorroBay","SFBay", "SanPabloBay", "SuisunBay", "HumboldtBay")

diff.means<-round(as.numeric(c(mor.t$estimate[2] - mor.t$estimate[1], sfbay.t$estimate[2] - sfbay.t$estimate[1], sanpab.t$estimate[2] - sanpab.t$estimate[1], suisun.t$estimate[2] - suisun.t$estimate[1], humbol.t$estimate[2] - humbol.t$estimate[1])),digits=3)

Up.conf.int<-round(as.numeric(c(-1*(mor.t$conf.int[1]),-1*(sfbay.t$conf.int[1]),-1*(sanpab.t$conf.int[1]),-1*(suisun.t$conf.int[1]),-1*(humbol.t$conf.int[1]))),digits=3)

low.conf.int<-as.numeric(c(-1*(mor.t$conf.int[2]),-1*(sfbay.t$conf.int[2]),-1*(sanpab.t$conf.int[2]),-1*(suisun.t$conf.int[2]),-1*(humbol.t$conf.int[2]))
)

low.conf.int<-round(low.conf.int,digits=3)

t.test.results<-cbind.data.frame(id,region.names, diff.means,low.conf.int,Up.conf.int)

print(str(t.test.results))

#plots difference in z-transformed means between Modern and Historical savannah sparrows across 5 coastal populations.

morro.col = rgb(136/255,204/255,238/255)
sanpablo.col = rgb(51/255,34/255,136/255)
sfbay.col = rgb(17/255,119/255,51/255)
suisun.col = rgb(170/255,68/255,153/255)
humboldt.col = rgb(221/255,204/255,119/255)
 
BillSizeChange.plot<-ggplot(data=t.test.results, mapping=aes(x = reorder(region.names,diff.means), y=diff.means, color=region.names)) + geom_hline(yintercept=c(-0.028,0.028), color="red", size=0.5, linetype="dashed") + scale_fill_manual(values=c(morro.col, sanpablo.col, sfbay.col, suisun.col, humboldt.col)) + theme_bw() +
geom_pointrange(mapping=aes(ymin=low.conf.int, ymax=Up.conf.int)) + labs(x="", y="Temporal change in bill surface area (z-transformed)") + coord_flip() + theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

print(BillSizeChange.plot)

ggsave("BillSizeChange.drift.pdf", plot=BillSizeChange.plot, width=9, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#####################################################################################################
############## How does climate and habitat change influence variation in bill change?###############

#Upload file with differences in bill morphology, climate, habitat, etc. between historic and modern periods
#Also includes estimates of Fst divergence between P. s. beldingi vs. alaudinus and P.s. nevadensis/brooksi vs. #alaudinus.

SAVSdiff<-read.csv("~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/FinalData_Code_BenhamBowie_2020_EvolApplications/BillDifferences_BenhamBowie_EvolApplications_2020.csv",header=TRUE)

Region<-SAVSdiff$Region
SAVSdiff2.z<-scale(SAVSdiff[,c(2:25)], center=TRUE, scale=TRUE)
SAVSdiff2<-cbind.data.frame(Region,SAVSdiff2.z)

#Principal components analysis of temporal change in different climatic, human/urbanization, and habitat
#parameters. See Supplemental Table S5 for loadings. 
clim.pca<-prcomp(SAVSdiff2[,6:13])
human.pca<-prcomp(SAVSdiff2[,16:19])
habitat.pca<-prcomp(SAVSdiff2[,20:25])

print(summary(clim.pca))
print(clim.pca)

print(summary(human.pca))
print(human.pca)

print(summary(habitat.pca))
print(habitat.pca)

clim.PC1<-clim.pca$x[,1]
clim.PC2<-clim.pca$x[,2]

human.PC1<-human.pca$x[,1]
human.PC2<-human.pca$x[,2]

habitat.PC1<-habitat.pca$x[,1]
habitat.PC2<-habitat.pca$x[,2]

SAVSdiff.final<-cbind.data.frame(SAVSdiff2, clim.PC1, clim.PC2, human.PC1, human.PC2, habitat.PC1, habitat.PC2)
Bill.character<-c("blank","Surface area", "Bill Length", "Bill Width", "Bill Depth")
Modnames<-c("Clim PC1", "Clim PC2", "Human PC1", "Human PC2", "Habitat PC1", "Habitat PC2", "Fst interior", "Fst beldingi", "Tmax", "Tmin", "Prec","VPDmax", "VPDmin", "TDmean", "SampleSalinity","PctPopulation", "CurrentPop","HumanPopDensity","HousingDensity", "TidalMarshAcresCurrent","PercentOriginalTidalHabitat","PctTidalMarsh","PctUrbanHabitat","PctCropland","PctPasture")
#AIC analysis comparing fit of various linear models to spatial variation in temporal change of bill characters.
#AIC results are reported in supplemental tables S6 & S7.
for (i in 2:5) {
print(Bill.character[i])	
Cand.models<-list()
Cand.models[[1]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$clim.PC1)
Cand.models[[2]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$clim.PC2)
Cand.models[[3]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$human.PC1)
Cand.models[[4]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$human.PC2)
Cand.models[[5]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$habitat.PC1)
Cand.models[[6]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$habitat.PC2)
Cand.models[[7]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$Fst_int)
Cand.models[[8]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$Fst_bel)
Cand.models[[9]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$Tmax)
Cand.models[[10]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$Tmin)
Cand.models[[11]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$Prec)
Cand.models[[12]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$VPDmax)
Cand.models[[13]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$VPDmin)
Cand.models[[14]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$TDmean)
Cand.models[[15]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$SampleSalinity)
Cand.models[[16]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$PctPopulation)
Cand.models[[17]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$CurrentPop)
Cand.models[[18]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$HumanPopDensity)
Cand.models[[19]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$HousingDensity)
Cand.models[[20]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$TidalMarshAcresCurrent)
Cand.models[[21]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$PercentOriginalTidalHabitat)
Cand.models[[22]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$PctTidalMarsh)
Cand.models[[23]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$PctUrbanHabitat)
Cand.models[[24]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$PctCropland)
Cand.models[[25]]<-lm(SAVSdiff.final[,i]~SAVSdiff.final$PctPasture)


aic.table<-aictab(cand.set=Cand.models, modnames=Modnames, sort=TRUE)

print(aic.table, digits=2, LL=TRUE)
}

#plot of best model from AIC results above. Compares variation in temporal change of bill surface area across
#different estuaries and variation in climate change across same five estuaries.
#Plot shown as Figure 6 in main paper and regression results can be found in Table 5.
#surface area top model
SA.lm<-lm(SAVSdiff.final$SAeffect~SAVSdiff.final$clim.PC2)
print(summary(SA.lm))

SA_clim.plot<-ggplot(SAVSdiff.final, aes(y=SAeffect, x=clim.PC2, color=Region)) + theme_bw() + 			
	stat_smooth(method="lm", color="black") +
	geom_point(size=4) +
	labs(y="Bill surface area change (effect size)", x="Environment PC2 (35.3)") + 
	scale_color_manual(values=c(humboldt.col, morro.col,SanPablo.col,sfbay.col,suisun.col)) +
	theme(axis.text=element_text(size=16),axis.title=element_text(size=18)) 

print(SA_clim.plot)

ggsave("SA_clim.plot.pdf", plot=SA_clim.plot, width=7, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#Linear regression analyses and figures comparing spatial variation in temporal change between
#different environmental parameters and other bill morphology dimensions.
#This code generates supplemental figures S4 & S5 and regression results in Table 5.

#bill length top model
BL.lm<-lm(SAVSdiff.final$BLeffect~SAVSdiff.final$clim.PC2)
BL2.lm<-lm(SAVSdiff.final$BLeffect~SAVSdiff.final$habitat.PC1)
print(summary(BL.lm))
print(summary(BL2.lm))

BL_clim.plot<-ggplot(SAVSdiff.final, aes(y=BLeffect, x=clim.PC2, color=Region)) + theme_bw() + 						 
	stat_smooth(method="lm", color="black") +
	geom_point(size=4) +
	labs(y="Bill length change (effect size)", x="Environment PC2 (35.3)") + 
	scale_color_manual(values=c(humboldt.col, morro.col,SanPablo.col,sfbay.col,suisun.col)) +
	theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) 

BL_hab.plot<-ggplot(SAVSdiff.final, aes(y=BLeffect, x=habitat.PC1, color=Region)) + theme_bw() + 					
	stat_smooth(method="lm", color="black") +
	geom_point(size=4) +
	labs(y="Bill length change (effect size)", x="Habitat PC1 (41.4)") + 
	scale_color_manual(values=c(humboldt.col, morro.col,SanPablo.col,sfbay.col,suisun.col)) +
	theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) 

BillLength.plot<-grid.arrange(BL_clim.plot, BL_hab.plot, ncol=2)
print(BillLength.plot)

ggsave("BillLength.plot.pdf", plot=BillLength.plot, width=12, height=4, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#bill width top model
BW.lm<-lm(SAVSdiff.final$BWeffect~SAVSdiff.final$Fst_bel)
print(summary(BW.lm))

BW_fst.plot<-ggplot(SAVSdiff.final, aes(y=BWeffect, x=Fst_bel, color=Region)) + theme_bw() + 	
	stat_smooth(method="lm", color="black") +
	geom_point(size=4) + 
	labs(y="Bill width change (effect size)", x="beldingi - alaudinus Fst") + 
	scale_color_manual(values=c(humboldt.col, morro.col, SanPablo.col, sfbay.col, suisun.col)) +
	theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) 

#print(BW_fst.plot)

ggsave("BW_fst.plot.pdf", plot=BW_fst.plot, width=7, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#bill depth top model
BD.lm<-lm(SAVSdiff.final$BDeffect~SAVSdiff.final$Fst_int)
print(summary(BD.lm))


BD_fst.plot<-ggplot(SAVSdiff.final, aes(y=BDeffect, x=Fst_int, color=Region)) + theme_bw() + 				
	stat_smooth(method="lm", color="black") +
	geom_point(size=4) + 
	labs(y="Bill depth change (effect size)", x="Interior-coastal Fst") + 
	scale_color_manual(values=c(humboldt.col, morro.col,SanPablo.col,sfbay.col,suisun.col)) +
	theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) 

#print(BD_fst.plot)

ggsave("BD_fst.plot.pdf", plot=BD_fst.plot, width=7, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

#bill Fst plot
BillWidthDepth_Fst<-grid.arrange(BW_fst.plot,BD_fst.plot,ncol=2)

print(BillWidthDepth_Fst)

ggsave("BillWidthDepth_Fst.pdf", plot=BillWidthDepth_Fst, width=12, height=4, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

BD.clim.lm<-lm(SAVSdiff.final$BDeffect~SAVSdiff.final$clim.PC1)
print(summary(BD.clim.lm))

###########################################################################################################
################### Potential functional significance of change in bill surface area change ###############

#code to create Figure 7 of the main paper and predict the amount of water savings possible
#given magnitude of bill surface area change in different population. Results in supplemental Table S8.
TEWL.data<-read.csv("~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/FinalData_Code_BenhamBowie_2020_EvolApplications/TEWL_BillSize_BenhamBowie_EvolApplications_2020.csv",header=TRUE)

TEWL.data <- subset(TEWL.data, Subspp == "alaudinus")
SAtarsus.lm<-lm(BillSurfaceArea~Tarsus, data=TEWL.data)
print(summary(SAtarsus.lm))
TEWL.data$SAtarCorr<-residuals(SAtarsus.lm)

humboldt.col = rgb(221/255,204/255,119/255)
morro.col = rgb(136/255,204/255,238/255)
sfbay.col = rgb(17/255,119/255,51/255)
suisun.col = rgb(170/255,68/255,153/255)
print(tail(TEWL.data))

tewl.bill.plot<-ggplot(TEWL.data, aes(x=BillSurfaceArea,y=TEWL_g_h, group=Region, color=Region)) +
theme_bw() + geom_point(size=4)  + geom_smooth(method="lm", color="black", group=1) +  labs(y="Total evaporative water loss (mg/g/h)", x="Bill surface area (mm^2)") + scale_color_manual(values=c(humboldt.col, morro.col,sfbay.col,suisun.col)) + theme(axis.text=element_text(size=16),axis.title=element_text(size=18), legend.title=element_text(size=16, color="black"), legend.text=element_text(color="black", size=14)) +
geom_vline(xintercept=c(61.77,66.32), col=c("red","blue"), linetype="dotted", size=0.5) +
geom_hline(yintercept=c(6.214,5.326), col=c("red","blue"), linetype="dotted", size=0.5) 

ggsave("tewl.bill.relationship.pdf", plot=tewl.bill.plot, width=9, height=7, units=c("in"), device="pdf", dpi=500, path="~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/BillMorph_manuscript/FinalTables_Figures/")

print(tewl.bill.plot)

SA.tewl.lm<-lm(TEWL_g_h~BillSurfaceArea, data=TEWL.data)
print(summary(SA.tewl.lm))

humboldt.bill<-c(61.09, 58.73)
SanPablo.bill<-c(56.66, 59.73)
SFbay.bill<-c(58.17,59.23)
Morro.bill<-c(61.77,66.32)
Suisun.bill<-c(57.95, 58.78)

bill.predict<-data.frame(BillSurfaceArea=c(61.09, 58.73,56.66, 59.73,58.17,59.23,61.77,66.32,57.95, 58.78))
new.tewl<-predict(SA.tewl.lm, bill.predict)
print(new.tewl)

hum.diff<-new.tewl[2] - new.tewl[1]
spb.diff<-new.tewl[4] - new.tewl[3]
sfb.diff<-new.tewl[6] - new.tewl[5]
morro.diff<-new.tewl[8] - new.tewl[7]
suisun.diff<-new.tewl[10] - new.tewl[9]

locality<-c("Humboldt", "Suisun", "San Pablo","San Francisco", "Morro")
TewlChange<-c(hum.diff, suisun.diff, spb.diff,sfb.diff,morro.diff)

PredictTempChangeTEWL<-cbind.data.frame(locality,TewlChange)
print(PredictTempChangeTEWL)