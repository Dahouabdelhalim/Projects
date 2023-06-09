library(ape)
library(geiger)
library(phytools)
library(nlme)
library(MuMIn)
library(effects)
#trait information
data.file<-read.table(file='Lopez-Uribe_etal_data_FINAL.txt',sep="\\t",header=T)
head(data.file)
#data transformation of some variables
data.file$log.Body_size<-log(data.file$Body_size)
data.file$log.Body_size_Q<-log(data.file$Body_size_Q)
data.file$asin.sqrt_He<-asin(sqrt(data.file$He))
data.file$asin.sqrt_Gst<-asin(sqrt(data.file$Gst))
data.file$asin.sqrt_Gst_H<-asin(sqrt(data.file$Gst_H))
data.file$asin.sqrt_Dest<-asin(sqrt(data.file$Dest))
data.file$Z.Mean_Geo=as.vector(scale(data.file$Mean_Geo)) 
data.file$Z.Mean_Geo
data.frame(data.file[,2:23])


#tree information
phylo_raw<-read.tree(file='Lopez-Uribe_RAxML_Hedtke_bestTree.result.tre')
plot(phylo_raw, cex=0.5)
#check overlap between traits and phylogenetic tree
phylo<-drop.tip(phylo_raw,c("Bembix_amoena",
                         "Bembix_americana",
                         "Andrena_cineraria",
                         "Andrena_nasonii",
                         "Andrena_vaga4",
                         "Anthophora_montana",
                         "Apis_andreniformis",
                         "Apis_cerana1",
                         "Apis_cerana2",
                         "Apis_cerana3",
                         "Apis_dorsata",
                         "Apis_florea",
                         "Apis_mellifera2",
                         "Apis_mellifera3",
                         "Apis_mellifera4",
                         "Bombus_bifarius2",
                         "Bombus_bifarius3",
                         "Bombus_pascuorum3",
                         "Colletes_pascoensis",
                         "Colletes_seminitidus",
                         "Dasypoda_argentata",
                         "Dasypoda_hirtipes",
                         "Dasypoda_visnaga",
                         "Euglossa_piliventris",
                         "Eufriesea_surinamensis",
                         "Exoneura_bicolor",
                         "Hylaeus_proximus",
                         "Lasioglossum_umbripenne",
                         "Megachile_Chelostomoda_sp.",
                         "Megachile_ericetorum",
                         "Megachile_mandibularis",
                         "Melipona_quadrifasciata",
                         "Melipona_quinquefasciata",
                         "Nannotrigona_perilampoides",
                         "Osmia_apicata",
                         "Plebeia_minima",
                         "Scaptotrigona_mexicana",
                         "Tetragonisca_weyrauchi",
                         "Tetraloniella_glauca",
                         "Trigona_corvina",
                         "Osmia_rufa","Bombus_flavifrons","Bombus_mixtus",
                         "Bombus_vosnesenskii",
                         "Megachile_Chelostomoda_sp._CJP",
                         "Osmia_cornuta",
                         "Bombus_melanopygus",
                         "Bombus_sylvicola",
                         "Euglossa_imperialis",
                         "Megachile_rotundata",
                         "Megachile_rotundata",
                         "Scaptotrigona_xanthotricha"),trim.internal=TRUE)
plot(phylo, cex=0.5)

#library(geiger)
#check if phylogenetic tree matches the observations in the summary table
rownames(data.file)<-(data.file[,1])
name.check(phylo, data.file)


##############PGLS Nei's Gst####################
#library(nlme)
gls.phylo_Gst<-drop.tip(phylo,c("Euglossa_championi2","Eulaema_bombiformis"))
name.check(gls.phylo_Gst, data.file)

#Phylogenetically uncorrected model
options(na.action="na.fail")
gls8_Gst<-gls(Gst~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,data=data.file,na.action=na.exclude, method="ML")
summary(gls8_Gst)
plot(Gst~log.Body_size_Q,data=data.file,pch=19, ylab="Nei's Gst", xlab="Body Size (log)")
abline(gls8_Gst)

#Brownian motion
bm.data_Gst<-corBrownian(phy=gls.phylo_Gst)
bm.gls8_Gst<-gls(Gst~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,correlation=bm.data_Gst,data=data.file,na.action=na.exclude,method="ML")
summary(bm.gls8_Gst) 

#OU Model
ou.data_Gst<-corMartins(1,phy=gls.phylo_Gst)
ou.gls8_Gst<-gls(Gst~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,correlation=ou.data_Gst,data=data.file,na.action=na.exclude,method="ML")
summary(ou.gls8_Gst)

#model selection
#library(MuMIn)
best_full_Gst_model<- model.sel(gls8_Gst,bm.gls8_Gst,ou.gls8_Gst, rank=AIC)
best_full_Gst_model

#best model (combination of variables) within OU models
ou.data_Gst<-corMartins(1,phy=gls.phylo_Gst)
ou.gls8_Gst<-gls(Gst~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,correlation=ou.data_Gst,data=data.file,na.action=na.exclude,method="ML")
summary(ou.gls8_Gst)

ou.gls1_Gst<-gls(Gst~log.Body_size+Z.Mean_Geo+He,correlation=ou.data_Gst, data=data.file,na.action=na.exclude, method="ML")
summary(ou.gls1_Gst)
##########FIGURE##########
plot(Gst~log.Body_size,data=data.file,pch=19, ylab="Nei's Gst", xlab="Body Size (log)")
abline(ou.gls1_Gst)
##########END FIGURE##########
ou.gls2_Gst<-gls(Gst~Diet+Z.Mean_Geo+He,correlation=ou.data_Gst, data=data.file,na.action=na.exclude, method="ML")
summary(ou.gls2_Gst)
##########FIGURE##########
boxplot(Gst~Diet,data=data.file,pch=19, ylab="Nei's Gst", xlab="Diet Breadth)")
##########END FIGURE##########
ou.gls3_Gst<-gls(Gst~Behavior+Z.Mean_Geo+He,correlation=ou.data_Gst, data=data.file,na.action=na.exclude, method="ML")
summary(ou.gls3_Gst)
##########FIGURE##########
boxplot(Gst~Behavior,data=data.file,pch=19, ylab="Nei's Gst", xlab="Behavior)")
##########END FIGURE##########
ou.gls4_Gst<-gls(Gst~log.Body_size+Diet+Z.Mean_Geo+He,correlation=ou.data_Gst, data=data.file,na.action=na.exclude, method="ML")
summary(ou.gls4_Gst)
ou.gls5_Gst<-gls(Gst~log.Body_size+Behavior+Z.Mean_Geo+He,correlation=ou.data_Gst, data=data.file,na.action=na.exclude, method="ML")
summary(ou.gls5_Gst)
ou.gls6_Gst<-gls(Gst~Diet+Behavior+Z.Mean_Geo+He,correlation=ou.data_Gst, data=data.file,na.action=na.exclude, method="ML")
summary(ou.gls6_Gst)


best_gls_modelGst<- model.sel(ou.gls8_Gst, ou.gls1_Gst,ou.gls2_Gst,ou.gls3_Gst,ou.gls4_Gst,ou.gls5_Gst,ou.gls6_Gst, rank = AIC)
best_gls_modelGst

options(na.action="na.omit")
plot(allEffects(ou.gls5_Gst))


##############PGLS Hedricks's Gst####################
#library(nlme)
gls.phylo_GstH<-drop.tip(phylo,c("Bombus_pascuorum1","Bombus_pascuorum2","Bombus_sylvarum","Bombus_terrestris1","Bombus_terrestris2","Bombus_vosnesenskii1","Bombus_vosnesenskii2","Macrotera_portalis","Melipona_scutellaris","Eulaema_atleticana"))
name.check(gls.phylo_GstH, data.file)

#Phylogenetically uncorrected model
options(na.action="na.fail")
gls8_GstH<-gls(Gst_H~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,data=data.file,na.action=na.exclude, method="ML")
summary(gls8_GstH)

#Brownian motion
bm.data_GstH<-corBrownian(phy=gls.phylo_GstH)
bm.gls8_GstH<-gls(Gst_H~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,correlation=bm.data_GstH,data=data.file,na.action=na.exclude,method="ML")
summary(bm.gls8_GstH)

#OU Model
ou.data_GstH<-corMartins(1,phy=gls.phylo_GstH)
ou.gls8_GstH<-gls(Gst_H~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,correlation=ou.data_GstH,data=data.file,na.action=na.exclude,method="ML")
summary(ou.gls8_GstH)

#model selection
best_full_Gst_H_model<- model.sel(gls8_GstH,bm.gls8_GstH,ou.gls8_GstH, rank=AIC)
best_full_Gst_H_model

#best model (combination of variables) within non phylogenetically corrected (GLS) for Gst_H
gls1_GstH<-gls(Gst_H~log.Body_size+He+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls1_GstH)
gls2_GstH<-gls(Gst_H~Diet+He+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls2_GstH)
gls3_GstH<-gls(Gst_H~Behavior+He+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls3_GstH)
gls4_GstH<-gls(Gst_H~log.Body_size+Diet+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls4_GstH)
gls5_GstH<-gls(Gst_H~log.Body_size+Behavior+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls5_GstH)
gls6_GstH<-gls(Gst_H~Diet+Behavior+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls6_GstH)

best_gls_modelGst_H<- model.sel(gls8_GstH,gls1_GstH,gls2_GstH,gls3_GstH,gls4_GstH,gls5_GstH,gls6_GstH, rank=AIC)
best_gls_modelGst_H

options(na.action="na.omit")
plot(allEffects(gls5_GstH))

##############PGLS Dest####################
#library(nlme)
gls.phylo_Dest<-drop.tip(phylo,c("Bombus_impatiens2","Bombus_jonellus","Bombus_muscorum","Bombus_pascuorum1","Bombus_pascuorum2","Bombus_pennsylvanicus2","Bombus_sylvarum","Bombus_terrestris1","Bombus_terrestris2","Melipona_scutellaris","Scaptotrigona_hellwegeri", "Bombus_ruderatus"))
name.check(gls.phylo_Dest, data.file)

#Phylogenetically uncorrected model
options(na.action="na.fail")
gls8_Dest<-gls(Dest~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,data=data.file,na.action=na.exclude, method="ML")
summary(gls8_Dest)

#Brownian motion for Dest
bm.data_Dest<-corBrownian(phy=gls.phylo_Dest)
bm.gls8_Dest<-gls(Dest~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,correlation=bm.data_Dest,data=data.file,na.action=na.exclude,method="ML")
summary(bm.gls8_Dest)

#OU Model
ou.data_Dest<-corMartins(1,phy=gls.phylo_Dest)
ou.gls8_Dest<-gls(Dest~log.Body_size+Diet+Behavior+Z.Mean_Geo+He, correlation=ou.data_Dest,data=data.file,na.action=na.exclude,method="ML")
summary(ou.gls8_Dest)


#model selection
best_full_Dest_model<- model.sel(gls8_Dest,bm.gls8_Dest,ou.gls8_Dest, rank=AIC)
best_full_Dest_model

#best model (combination of variables) within GLS models for Dest
gls8_Dest<-gls(Dest~log.Body_size+Diet+Behavior+Z.Mean_Geo+He,data=data.file,na.action=na.exclude, method="ML")
summary(gls8_Dest)
gls1_Dest<-gls(Dest~log.Body_size+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls1_Dest)
gls2_Dest<-gls(Dest~Diet+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls2_Dest)
gls3_Dest<-gls(Dest~Behavior+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls3_Dest)
gls4_Dest<-gls(Dest~log.Body_size+Diet+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls4_Dest)
gls5_Dest<-gls(Dest~log.Body_size+Behavior+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls5_Dest)
gls6_Dest<-gls(Dest~Diet+Behavior+Z.Mean_Geo+He, data=data.file,na.action=na.exclude, method="ML")
summary(gls6_Dest)

best_gls_modelDest<- model.sel(gls8_Dest,gls1_Dest,gls2_Dest,gls3_Dest,gls4_Dest,gls5_Dest,gls6_Dest, rank = AIC)
best_gls_modelDest

options(na.action="na.omit")
plot(allEffects(gls5_Dest))

################### FIGURE 1 ###################
GLM_Dest<-gls(Dest~log.Body_size,data=data.file,na.action=na.exclude, method="ML")
summary(GLM_Dest)
plot(Dest~log.Body_size, data=data.file, cex=1.2, pch=19, xlab="Body Size (log)", ylab="Jost D")
abline(GLM_Dest,lty=1)
#####################
pdf("Figure1.pdf")
par(mfrow=c(2,2))
plot(Dest~log.Body_size, data=data.file, cex=1, pch=20, xlab="Body Size (log)", ylab="Jost's D")
abline(GLM_Dest,lty=1)
boxplot(Dest~Behavior,col=c("white","lightgrey"),xlab="Sociality", ylab="Jost D", outline=FALSE, 
        ylim=c(-0.01,0.5), data=data.file, na.action=na.exclude)
stripchart(Dest~Behavior, vertical = TRUE, data = data.file, 
           method = "jitter", add = TRUE, pch = 20, col = 'black', na.action=na.exclude)
dev.off()
#####################
GLM_Gst<-gls(Gst~log.Body_size,data=data.file,na.action=na.exclude, method="ML")
summary(GLM_Gst)
GLM_Gst_H<-gls(Gst_H~log.Body_size,data=data.file,na.action=na.exclude, method="ML")
summary(GLM_Gst_H)

pdf("FigureS2.pdf")
par(mfrow=c(2,2))
plot(Gst~log.Body_size, data=data.file, cex=1, pch=20, xlab="Body Size (log)", ylab="Nei's Gst")
abline(GLM_Gst,lty=1)
boxplot(Gst~Behavior,col=c("white","lightgrey"),xlab="Sociality", ylab="Nei's Gst", 
        outline=FALSE, ylim=c(-0.01,0.4), data=data.file, na.action=na.exclude)
stripchart(Gst~Behavior, vertical = TRUE, data = data.file, 
           method = "jitter", add = TRUE, pch = 20, col = 'black', na.action=na.exclude)
plot(Gst_H~log.Body_size, data=data.file, cex=1, pch=20, xlab="Body Size (log)", ylab="Hedricks's Gst")
abline(GLM_Gst_H,lty=1)
boxplot(Gst_H~Behavior,col=c("white","lightgrey"),xlab="Sociality", ylab="Hedricks's Gst", outline=FALSE, ylim=c(-0.01,0.7), data=data.file, na.action=na.exclude)
stripchart(Gst_H~Behavior, vertical = TRUE, data = data.file, 
           method = "jitter", add = TRUE, pch = 20, col = 'black', na.action=na.exclude)
dev.off()

pdf("FigureS3.pdf")
par(mfrow=c(3,3))

boxplot(Gst~Diet,col=c("white","lightgrey"),xlab="Diet Breadth", ylab="Nei's Gst", outline=FALSE, ylim=c(-0.01,0.4), data=data.file, na.action=na.exclude)
stripchart(Gst~Diet, vertical = TRUE, data = data.file, 
           method = "jitter", add = TRUE, pch = 20, col = 'black',na.action=na.exclude)

boxplot(Gst_H~Diet,col=c("white","lightgrey"),xlab="Diet Breadth", ylab="Hedricks's Gst", outline=FALSE, ylim=c(-0.01,0.7), data=data.file, na.action=na.exclude)
stripchart(Gst_H~Diet, vertical = TRUE, data = data.file, 
           method = "jitter", add = TRUE, pch = 20, col = 'black',na.action=na.exclude)

boxplot(Dest~Diet,col=c("white","lightgrey"),xlab="Diet Breadth", ylab="Jost's D", outline=FALSE, ylim=c(-0.01,0.7), data=data.file, na.action=na.exclude)
stripchart(Gst_H~Diet, vertical = TRUE, data = data.file, 
           method = "jitter", add = TRUE, pch = 20, col = 'black',na.action=na.exclude)

dev.off()

#####################