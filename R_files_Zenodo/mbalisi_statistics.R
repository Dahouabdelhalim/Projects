# Code for Balisi et al. (submitted)
# "Dietary specialization is linked to reduced species durations in North American fossil canids"

# This file contains scripts for 
# principal component analysis, linear regressions, quantile regressions
# and the associated plots



# initialize packages
library(ggplot2)
library(quantreg)
library(plot3D)
library(cowplot)



# Van Valkenburgh et al (2004)'s trait list: JD/DL, RBL, and RUGA


# first: create comparative dataset

# read in extant comparative dataset
# these data are available in the Supplementary Information PDF
# comp.data = read.csv("comp-data_45spp-excl-feliforms_3vars.csv", header=1)
# comp.data.num.3 = comp.data[, 4:6]  # extract the three traits
# comp.data.diet.fam = comp.data[ , c(1, 8)]  # extract family and diet categories
# rownames(comp.data.num.3) = rownames(comp.data.diet.fam) = 
#   paste(comp.data$GENUS, sep="_", comp.data$SPECIES)  # assign each row to Genus_species

################################################################################
# FIGURE S5
# create morphospace based on extant comparative dataset
corrPCA3 = prcomp(comp.data.num.3, scale.=T)  # run correlation-matrix-based PCA
sdevs.corr3 = corrPCA3$sdev  # extract standard deviations per component
importance.corr3 = sdevs.corr3^2 / sum(sdevs.corr3^2)  # variance explained per component
biplot(corrPCA3, xlab=paste("PC 1 (", round(importance.corr3[1]*100, 2), "%)"), 
       ylab=paste("PC 2 (", round(importance.corr3[2]*100, 2), "%)"), cex=c(0.75, 1))
################################################################################

# create dataframes with extant PC scores for later plotting
extant.scores3 = as.data.frame(corrPCA3$x)
extant.scores3 = merge(extant.scores3, comp.data.diet.fam, by=0)
  # match family and diet categories to extant scores
rownames(extant.scores3) = extant.scores3$Row.names
extant.scores3 = extant.scores3[, 2:6]  # clean up


# PCA with fossils

# load the fossil unknowns 
fossils = read.csv("mbalisi_fossilData.R1_.csv")

fossils.num.3 = fossils[, 4:6]
  # extract VV et al (2004)'s trait list: JD/DL, RBL, and RUGA 
rownames(fossils.num.3) = paste(fossils$GENUS, sep="_", fossils$SPECIES)
  # assign each row to Genus_species
fossils.num.3 = fossils.num.3[complete.cases(fossils.num.3), ]
  # we can run PCA only on species with no NAs
  # 93 spp with no NAs 

# generate PC scores for fossil taxa based on extant-taxon morphospace
pc.fossil3 = predict(corrPCA3, fossils.num.3)

# create fossil-score dataframes for plotting
fossil.scores3 = as.data.frame(pc.fossil3)
fossil.scores3$FAMILY = fossils$SUBFAMILY[paste(fossils$GENUS, sep="_", 
                                                fossils$SPECIES) %in% 
                                            rownames(fossils.num.3)]
fossil.scores3$DIET = rep("unknown")


# create data frames with all scores
all.scores3 = rbind(fossil.scores3, extant.scores3)


# added Oct 6, 2017 - fossil / extant distinctions among "family" categories
colnames(all.scores3)[4] = "(SUB)FAMILY"
all.scores3$`(SUB)FAMILY` = as.character(all.scores3$`(SUB)FAMILY`)
  # change to character to alter value
all.scores3$`(SUB)FAMILY`[all.scores3$`(SUB)FAMILY`=="Canidae"] = "extant Canidae"
all.scores3$`(SUB)FAMILY`[all.scores3$`(SUB)FAMILY`=="Caninae"] = "fossil Caninae"
all.scores3$`(SUB)FAMILY` = as.factor(all.scores3$`(SUB)FAMILY`)
class(all.scores3$`(SUB)FAMILY`)  # ok now back to factor


# principal component ggplot of extant comparative database and fossils superimposed

################################################################################
# FIGURE 1
# for paper - modified 06/29/2017 with symbols
# and again modified 10/05/2017 with coding for diet AND family
# and modified another time 02/2018 with hollow triangles changed to squares / circles
all.scores3$`(SUB)FAMILY` = factor(all.scores3$`(SUB)FAMILY`, 
                                   levels=c("extant Canidae", "Hyaenidae", "Mephitidae", 
                                            "Mustelidae", "Procyonidae", "Hesperocyoninae", 
                                            "Borophaginae", "fossil Caninae"))

ggplot(data=all.scores3, aes(x=PC1, y=PC2, color=DIET, shape=`(SUB)FAMILY`)) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_point(data=all.scores3, aes(x=PC1, y=PC2), alpha=0.8, size=10) + 
  scale_shape_manual(values=c(15, 16, 17, 18, 20, 21, 22, 23)) +
  xlab(paste("PC 1 (", round(importance.corr3[1]*100, 2), "%)\\nincreasing carnivory")) +
  ylab(paste("PC 2 (", round(importance.corr3[2]*100, 2), "%)\\ndecreasing durophagy")) +
  theme_bw(base_size=25) +
  coord_fixed(ratio=1) +
  theme(legend.position="bottom") +
  scale_color_manual(values=rainbow(6, 0.6)) +
  guides(shape=guide_legend(title.position="top"), colour=guide_legend(title.position="top"))
################################################################################



# Setting up principal-component and other data for regression analyses 

# extract body mass, duration, locality coverage, and taxonomic data from fossils
MaMass = fossils[, c(7, 10, 11)]
rownames(MaMass) = paste(fossils$GENUS, sep="_", fossils$SPECIES)

# merge PC scores data with body mass, duration, locality coverage, and taxonomy
maMaMass = merge(fossil.scores3, MaMass, by=0)
rownames(maMaMass) = maMaMass$Row.names
maMa = maMaMass[, 2:9]
colnames(maMa)[4] = "SUBFAMILY"



################################################################################
# FIGURE 2

par(mar=c(6, 6, 0, 0)+0.1, oma=c(0, 0, 0, 0))

m = rbind(c(1, 1, 1), c(1, 1, 1), c(2, 3, 4))
layout(m)

# normalize PC1 so that median lies at 0
maMa$PC1trans = maMa$PC1 - median(maMa$PC1, na.rm=T) 
# positive values for hypercarnivores, negative values for hypocarnivores

plot(maMa$PC1trans, maMa$PyRate.Duration, pch=c(15:17)[maMa$SUBFAMILY], cex.lab=2,
     col=alpha(c("red", "green", "blue")[maMa$SUBFAMILY], 0.6), cex=3, cex.axis=2,
     xlab=expression(paste("increasing carnivory ", symbol('\\256'))), 
     ylab="species duration (Ma)")
abline(v=median(maMa$PC1trans), lty=2, lwd=2)
legend(1, 13, legend=levels(maMa$SUBFAMILY), pch=c(15:17)[unique(maMa$SUBFAMILY)],
       col=alpha(c("red", "green", "blue")[unique(maMa$SUBFAMILY)], 0.6), # bty="n",
       cex=2, y.intersp=0.5, text.width=1)

plot(10^maMa$LOGMASS, maMa$PyRate.Duration, pch=c(15:17)[maMa$SUBFAMILY], cex.lab=2,
     col=alpha(c("red", "green", "blue")[maMa$SUBFAMILY], 0.6), cex=2, cex.axis=2,
     xlab="body mass (kg)", ylab="species duration (Ma)", log="x")

plot(maMa$PC1, maMa$maxLocCover, pch=c(15:17)[maMa$SUBFAMILY], cex.lab=2,
     col=alpha(c("red", "green", "blue")[maMa$SUBFAMILY], 0.6), cex=2, cex.axis=2,
     xlab=expression(paste("increasing carnivory ", symbol('\\256'))), 
     ylab="maximum locality coverage")

plot(10^maMa$LOGMASS, maMa$maxLocCover, pch=c(15:17)[maMa$SUBFAMILY], cex.lab=2,
     col=alpha(c("red", "green", "blue")[maMa$SUBFAMILY], 0.6), cex=2, cex.axis=2,
     xlab="body mass (kg)", ylab="maximum locality coverage", log="x")
################################################################################



# check for relationships

# duration v carnivory
cor(maMa$PC1, maMa$PyRate.Duration, method="spearman", use="pairwise.complete.obs")
mDietDur = lm(PyRate.Duration ~ PC1, data=maMa, na.action=na.omit)
summary(mDietDur, se="boot", R=10000, bsmethod="xy")

# duration v mass
cor(maMa$LOGMASS, maMa$PyRate.Duration, method="spearman", use="pairwise.complete.obs")
mMassDur = lm(PyRate.Duration ~ LOGMASS, data=maMa, na.action=na.omit)
summary(mMassDur, se="boot", R=10000, bsmethod="xy")

# maxLocCover v carnivory
cor(maMa$PC1, maMa$maxLocCover, method="spearman", use="pairwise.complete.obs")
mMassLoc = lm(maxLocCover ~ PC1, data=maMa, na.action=na.omit)
summary(mMassLoc, se="boot", R=10000, bsmethod="xy")

# maxLocCover v mass
cor(maMa$LOGMASS, maMa$maxLocCover, method="spearman", use="pairwise.complete.obs")
mDietLoc = lm(maxLocCover ~ LOGMASS, data=maMa, na.action=na.omit)
summary(mDietLoc, se="boot", R=10000, bsmethod="xy")



# Quantile regression

# Before quantile regression... first need to DIVIDE the distribution into two: 
# above the median and below the median. 
# Keep these distributions separate.

# QR to be carried out only for duration v PC1
# so weed out singletons, i.e. species for which duration cannot be reliably calculated
maMaDur = maMa[is.na(maMa$PyRate.Duration)==FALSE, ]
# 87 resulting species with duration

# Divide 87-spp dataset between less and more carnivorous
maMaHyper = maMaDur[maMaDur$PC1 > median(maMaDur$PC1, na.rm=T), ]  # 43 spp
maMaHypo = maMaDur[maMaDur$PC1 <= median(maMaDur$PC1, na.rm=T), ]  # 44 spp
# there is one species almost right on the median (Eucyon davisi)
##  assigned it to be with "less carnivorous" because its PC1trans value is -0.000...
# Microtomarctus conferta's PC1trans = 0.000 (actual)
##  included automatically in "more carnivorous"

# Set hypos to positive for plotting
maMaHypo$PC1trans = -maMaHypo$PC1trans


# hypocarnivores: traditional linear regression
m1hypo = lm(PyRate.Duration ~ PC1trans, data=maMaHypo)
summary(m1hypo, se="boot", R=10000, bsmethod="xy")

# quantreg for hypocarnivores 
q1hypo = rq(PyRate.Duration ~ PC1trans, tau=seq(0.1, 0.9, 0.1), data=maMaHypo)
summary(q1hypo, se="boot", R=10000, bsmethod="xy")


# hypercarnivores: traditional linear regression
m1hyper = lm(PyRate.Duration ~ PC1trans, data=maMaHyper)
summary(m1hyper, se="boot", R=10000, bsmethod="xy")

# quantreg for hypercarnivores
q1hyper = rq(PyRate.Duration ~ PC1trans, tau=seq(0.1, 0.9, 0.1), data=maMaHyper)
summary(q1hyper, se="boot", R=10000, bsmethod="xy")

# remove Enhydrocyon + Epicyon haydeni as outliers. 
excl = c("Enhydrocyon_basilatus", "Enhydrocyon_crassidens", "Epicyon_haydeni")
maMaHyperNoE = maMaHyper[!rownames(maMaHyper) %in% excl, ]

# rerun linear regression without outliers
m1hyperNoE = lm(PyRate.Duration ~ PC1trans, data=maMaHyperNoE)
summary(m1hyperNoE, se="boot", R=10000, bsmethod="xy")

# rerun quantreg without hyper outliers
q1hyperNoE = rq(PyRate.Duration ~ PC1trans, tau=seq(0.1, 0.9, 0.1), data=maMaHyperNoE)
summary(q1hyperNoE, se="boot", R=10000, bsmethod="xy")


################################################################################
# FIGURE 3
# plotting quantreg

par(mfrow=c(1, 2), tcl=-0.5, mar=c(2, 2, 2, 0), oma=c(2, 3, 0, 0))

plot(PyRate.Duration ~ PC1trans, data=maMaHypo, col=alpha("black", 0.6), ylab="", 
     cex.axis=1.5, xaxp=c(0, 3, 3), cex.main=2, cex=2,
     main="less carnivorous canids (n=44)", xlim=range(maMaHyper$PC1trans))
grid()
abline(m1hypo, col=alpha("red", 0.6), lwd=2)
taus = c(0.6, 0.7, 0.8, 0.9)
for(i in taus){
  abline(rq(PyRate.Duration ~ PC1trans, data=maMaHypo, tau=i, na.action=na.omit),
         col=alpha("blue", 0.6), lwd=2)
}
text(2.5, c(0.25, 1.25, 4.5, 2), c(expression(paste(tau, "=0.9")), 
                                   expression(paste(tau, "=0.8")), 
                                   expression(paste(tau, "=0.7")), 
                                   expression(paste(tau, "=0.6"))))
     
plot(PyRate.Duration ~ PC1trans, data=maMaHyper, col=alpha("black", 0.6), yaxt="n", 
     xaxt="n", cex.axis=1.5, cex.main=2, cex=2,
     main="more carnivorous canids (n=43)", xlim=range(maMaHyper$PC1trans),
     ylim=range(maMaHypo$PyRate.Duration, na.rm=T), xlab="", ylab="")
axis(2, labels=F)
axis(1, at=c(0, 1, 2), cex.axis=1.5)
grid()
abline(m1hyperNoE, col=alpha("red", 0.6), lwd=2)
taus = c(0.6, 0.7, 0.8, 0.9)
for(i in taus){
  abline(rq(PyRate.Duration ~ PC1trans, data=maMaHyperNoE, tau=i), 
         col=alpha("blue", 0.6), lwd=2)
}
text(2.5, c(0.25, 0.75, 1.75, 2.5), c(expression(paste(tau, "=0.9")), 
                                expression(paste(tau, "=0.8")), 
                                expression(paste(tau, "=0.7")), 
                                expression(paste(tau, "=0.6"))))

title(ylab="species duration (Ma)", outer=T, line=1, cex.lab=2,
      xlab=expression(paste("increasing specialization ", symbol('\\256'))))
################################################################################



# return to maMaMass because that has ALL values including NAs for duration
##  because it would be good to show how many species on the specialized ends
##  have no duration and/or locality coverage data

# normalize PC1 so that median lies at 0
maMaMass$PC1trans = maMaMass$PC1 - median(maMaMass$PC1, na.rm=T) 
# positive values for hypercarnivores, negative values for hypocarnivores

# divide maMaMass into less- and more-carnivorous
maMaMassHyper = maMaMass[maMaMass$PC1 > median(maMaMass$PC1, na.rm=T), ]  # 46 spp
maMaMassHypo = maMaMass[maMaMass$PC1 <= median(maMaMass$PC1, na.rm=T), ]  # 47 spp

# Set hypos to positive for plotting
maMaMassHypo$PC1trans = -maMaMassHypo$PC1trans


# exploratory plot of all species
cc = ggplot(data=maMaMass, aes(x=PC1trans, y=10^LOGMASS, color=PyRate.Duration, 
                                 shape=FAMILY)) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_point(data=maMaMass, aes(x=PC1trans, y=10^LOGMASS), alpha=0.8, size=10) + 
  xlab("increasing carnivory >") +
  ylab("body mass (kg)") +
  theme_bw(base_size=25) +
  scale_x_continuous(breaks=c(-2, 1, 0, 1, 2)) +
  scale_y_continuous(breaks=c(5, 10, 20, 40), trans="log10") +
  scale_shape_manual(values=c(15, 18, 16)) +
  scale_color_gradientn(colors=rev(jet.col((64))), 
                        limits=range(maMa$PyRate.Duration, na.rm=T)) +
  theme(plot.margin=unit(c(6,0,6,0), "pt")) +
  labs(shape="Subfamily", colour="Duration (Ma)") +
  theme(legend.position="none")


# less-carnivorous species only
aa = ggplot(data=maMaMassHypo, aes(x=PC1trans, y=10^LOGMASS, color=PyRate.Duration, 
                           shape=FAMILY)) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_point(data=maMaMassHypo, aes(x=PC1trans, y=10^LOGMASS), alpha=0.8, size=10) + 
  xlab(expression(paste("increasing specialization (hypocarn) ", symbol('\\256')))) +
  ylab("body mass (kg)") +
  theme_bw(base_size=25) +
  scale_x_continuous(breaks=c(0, 1, 2)) +
  scale_y_continuous(limits=range(10^maMaMass$LOGMASS, na.rm=T), breaks=c(5, 10, 20, 40),
                     trans="log10") +
  scale_shape_manual(values=c(15, 18, 16)) +
  scale_color_gradientn(colors=rev(jet.col((64))), 
                        limits=range(maMaMass$PyRate.Duration, na.rm=T)) +
  theme(plot.margin=unit(c(6,0,6,0), "pt"))

# more carnivorous species only
bb = ggplot(data=maMaMassHyper, aes(x=PC1trans, y=10^LOGMASS, color=PyRate.Duration, 
                              shape=FAMILY)) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_point(data=maMaMassHyper, aes(x=PC1trans, y=10^LOGMASS), alpha=0.8, size=10) + 
  xlab(expression(paste("increasing specialization (hypercarn) ", symbol('\\256')))) +
  ylab("body mass (kg)") +
  theme_bw(base_size=25) +
  scale_y_continuous(limits=range(10^maMaMass$LOGMASS, na.rm=T), breaks=c(5, 10, 20, 40),
                     trans="log10") +
  scale_shape_manual(values=c(15, 18, 16)) +
  scale_color_gradientn(colors=rev(jet.col((64))), 
                        limits=range(maMaMass$PyRate.Duration, na.rm=T)) +
  theme(plot.margin=unit(c(6, 0, 6, 0), "pt"), legend.position=c(0.75, 0.25)) + 
  ylab("") +
  labs(shape="Subfamily", colour="Duration (Ma)")

################################################################################
# FIGURE 4
# plot in a grid
plot_grid(aa + theme(legend.position="none"), bb,
                 align="h", labels=c("A", "B"), hjust=-1, nrow=1)
################################################################################



################################################################################
# FIGURE S6
# Relationships between dietary measures and body mass

plot(fossils[ , 4:7])
################################################################################