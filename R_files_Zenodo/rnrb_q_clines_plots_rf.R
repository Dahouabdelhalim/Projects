### rnrb_q_clines_plots.R
### Libby Natola
### Started 22 July 2021, major revisions 9 December 2021
### Plot admixture q values per individual along a distance transect

# set working directory
setwd("~/Documents/UBC/Bioinformatics/rnrb/clines/")

# libraries
library(ggplot2)
library(geosphere)
library(hzar)
library(viridis)
library(dplyr)
library(forcats) 
library(sp)

# load the admixture data, turn it into a df, lad the id data, label the columns
nbc_admix <- read.table("../admixture/nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.pruned.2.Q")
nbc_admix <- as.data.frame(nbc_admix)
nbc_fam <- read.table("../admixture/nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.pruned.fam")
nbc_fam <- as.data.frame(nbc_fam)
nbc_admix <- cbind(nbc_fam$V2, nbc_admix)
colnames(nbc_admix) <- c("ID", "p", "q")

#load the species, location, and isocline distance data
nbc_info <- read.table("../sample_info_nbc.txt", header=TRUE)
nbc_admix <- left_join(nbc_admix, nbc_info)
nbc_dist <- read.csv("../isoclines/nbc_dist.csv")
nbc_admix$dist <- nbc_dist$distance


# load the admixture data, turn it into a df, load the id data, label the columns
caor_admix <- read.table("../admixture/caor_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.2.Q")
caor_admix <- as.data.frame(caor_admix)
caor_fam <- read.table("../admixture/caor_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.fam")
caor_fam <- as.data.frame(caor_fam)
caor_admix <- cbind(caor_fam$V2, caor_admix)
colnames(caor_admix) <- c("ID", "p", "q")

#load the species, location, and isocline distance data
caor_info <- read.table("../sample_info_caor.txt", header=TRUE)
caor_admix <- left_join(caor_admix, caor_info)
caor_dist <- read.csv("../isoclines/caor_dist.csv")
caor_admix$dist <- caor_dist$distance


#load the species and location data into a df, load the id data, label the columns
man_admix <- read.table("../admixture/man_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.2.Q")
man_admix <- as.data.frame(man_admix)
man_fam <- read.table("../admixture/man_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.fam")
man_fam <- as.data.frame(man_fam)
man_admix <- cbind(man_fam$V2, man_admix)
colnames(man_admix) <- c("ID", "q", "p")

#load the species, location, and isocline distance data
man_info <- read.table("../sample_info_man.txt", header=TRUE)
man_admix <- left_join(man_admix, man_info)
man_dist <- read.csv("../isoclines/man_dist.csv")
man_admix$dist <- man_dist$distance


# just taking a look-see
ggplot(man_admix, aes(x=long, y=lat, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + ggtitle("man lat-long")
ggplot(nbc_admix, aes(x=long, y=lat, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + ggtitle("nbc lat-long")
ggplot(caor_admix, aes(x=long, y=lat, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + ggtitle("caor lat-long")


ggplot(nbc_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient")

ggplot(man_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient")

ggplot(caor_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient")

### can I find the distance of overlap between the westernmost sample with q ~ 0 and the easternmost sample with q ~ 1
nbc_rns <- filter(nbc_admix, q < 0.05)
nbc_rbs <- filter(nbc_admix, q > 0.95)
nbc_rn_min <- min(nbc_rns$dist) 
nbc_rn_min
nbc_rb_max <- max(nbc_rbs$dist)
nbc_rb_max
nbc_rb_max-nbc_rn_min
#[1] -14.32969s


pdf(file="nbc_symptatry.pdf")
theme_set(theme_bw())
theme <- theme_update(text = element_text(size=18), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_rect(color=rgb(253/255, 231/255, 37/255), size = 2))
ggplot(nbc_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient") + annotate("rect", xmin = nbc_rn_min, xmax = nbc_rb_max, ymin = 0.0, ymax = 1, alpha = .5) + annotate(geom="text", x=nbc_rb_max-((nbc_rb_max-nbc_rn_min)/2), y=0.65, label="14.3 KM",color="black")
dev.off()

man_rns <- filter(man_admix, q < 0.05)
man_rbs <- filter(man_admix, q > 0.95)
man_rn_min <- min(man_rns$dist) 
man_rn_min
man_rb_max <- max(man_rbs$dist)
man_rb_max
man_rb_max-man_rn_min
# [1] 4.960015

pdf(file="man_sympatry.pdf")
theme_set(theme_bw())
theme <- theme_update(text = element_text(size=18), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_rect(color=rgb(68/255, 1/255, 84/255), size = 2))
ggplot(man_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient") + annotate("rect", xmin = man_rn_min, xmax = man_rb_max, ymin = 0.0, ymax = 1, alpha = .5) + annotate(geom="text", x=man_rb_max-((man_rb_max-man_rn_min)/2), y=0.25, label="5.0 km",color="black")
dev.off()

caor_rns <- filter(caor_admix, q < 0.05)
caor_rbs <- filter(caor_admix, q > 0.95)
caor_rn_min <- min(caor_rns$dist) 
caor_rn_min
caor_rb_max <- max(caor_rbs$dist)
caor_rb_max
caor_rb_max-caor_rn_min
# [1] 53.93073

pdf(file="caor_sympatry.pdf")
theme_set(theme_bw())
theme <- theme_update(text = element_text(size=18), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_rect(color=rgb(33/255, 145/255, 140/255), size = 2))
ggplot(caor_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient") + annotate("rect", xmin = caor_rn_min, xmax = caor_rb_max, ymin = 0.0, ymax = 1, alpha = .5) + annotate(geom="text", x=caor_rb_max-((caor_rb_max-caor_rn_min)/2), y=0.25, label="53.9 km",color="black")
dev.off()

### can I find the distance from p > 0.05 to p < 0.95

nbc_hybs <- filter(nbc_admix, q >= 0.05 & q <= 0.95)
nbc_hybs_min <- min(nbc_hybs$dist) 
nbc_hybs_min
nbc_hybs_max <- max(nbc_hybs$dist)
nbc_hybs_max
nbc_hybs_max-nbc_hybs_min
#[1] 260.4347

pdf(file="nbc_admix_range.pdf")
theme_set(theme_bw())
theme <- theme_update(text = element_text(size=18), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_rect(color=rgb(253/255, 231/255, 37/255), size = 2))
ggplot(nbc_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient") + annotate("rect", xmin = nbc_hybs_min, xmax = nbc_hybs_max, ymin = 0.0, ymax = 1, alpha = .5) + annotate(geom="text", x=nbc_hybs_max-((nbc_hybs_max-nbc_hybs_min)/2), y=0.6, label="<          260.4 KM         >",color="black")
dev.off()

man_hybs <- filter(man_admix, p >= 0.05 & p <= 0.95)
man_hybs_min <- min(man_hybs$dist) 
man_hybs_max <- max(man_hybs$dist)
man_hybs_max
man_hybs_max-man_hybs_min
# [1] 13.2786

pdf(file="man_admix_range.pdf")
theme_set(theme_bw())
theme <- theme_update(text = element_text(size=18), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_rect(color=rgb(68/255, 1/255, 84/255), size = 2))
ggplot(man_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient") + annotate("rect", xmin = man_hybs_min, xmax = man_hybs_max, ymin = 0.0, ymax = 1, alpha = .5) + annotate(geom="text", x=man_hybs_max-((man_hybs_max-man_hybs_min)/2), y=0.25, label="13.3 KM",color="black")
dev.off()

caor_hybs <- filter(caor_admix, q > 0.05 & q < 0.95)
caor_hybs_min <- min(caor_hybs$dist) 
caor_hybs_max <- max(caor_hybs$dist)
caor_hybs_max
caor_hybs_max-caor_hybs_min
# [1] 212.9006

pdf(file="caor_admix_range.pdf")
theme_set(theme_bw())
theme <- theme_update(text = element_text(size=18), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_rect(color=rgb(33/255, 145/255, 140/255), size = 2))
ggplot(caor_admix, aes(x=dist, y=q, color=spp)) + geom_point() + scale_color_manual(values=c("#DE73FF", "#FF4343", "#4343FF")) + xlab("Distance in KM") + ylab("Ancestry Coefficient") + annotate("rect", xmin = caor_hybs_min, xmax = caor_hybs_max, ymin = 0.0, ymax = 1, alpha = .5) + annotate(geom="text", x=caor_hybs_max-((caor_hybs_max-caor_hybs_min)/2), y=0.25, label="<     212.9 KM    >",color="black")
dev.off()


###################### make hzar clines ######################
# commented out all the models except the one chosen by AICc for easier plug n chug

### add an nSamples column with 1 for each because we are using samples not pops
man_admix$nSamples <- rep(1, nrow(man_admix))

### make data object of allele frequency data
man_Q <-
  hzar.doMolecularData1DPops(man_admix$dist,
                             man_admix$q,
                             man_admix$nSamples)

### Plot the associated observed frequency versus distance
hzar.plot.obsData(man_Q)


### Construct a clineMetaModel object for use with hzar.first.fitRequest.old.ML
#man_Q_model_free_both <- hzar.makeCline1DFreq(man_Q, scaling="free", tails="both")
#man_Q_model_free_none <- hzar.makeCline1DFreq(man_Q, scaling="free", tails="none")
#man_Q_model_free_right <- hzar.makeCline1DFreq(man_Q, scaling="free", tails="right")
# man_Q_model_free_left <- hzar.makeCline1DFreq(man_Q, scaling="free", tails="left")
# man_Q_model_free_mirror <- hzar.makeCline1DFreq(man_Q, scaling="free", tails="mirror")
# # 
# man_Q_model_fixed_both <- hzar.makeCline1DFreq(man_Q, scaling="fixed", tails="both")
# man_Q_model_fixed_none <- hzar.makeCline1DFreq(man_Q, scaling="fixed", tails="none")
# man_Q_model_fixed_right <- hzar.makeCline1DFreq(man_Q, scaling="fixed", tails="right")
# man_Q_model_fixed_left <- hzar.makeCline1DFreq(man_Q, scaling="fixed", tails="left")
# man_Q_model_fixed_mirror <- hzar.makeCline1DFreq(man_Q, scaling="fixed", tails="mirror")
# # 
# man_Q_model_none_both <- hzar.makeCline1DFreq(man_Q, scaling="none", tails="both")
man_Q_model_none_none <- hzar.makeCline1DFreq(man_Q, scaling="none", tails="none")
# man_Q_model_none_right <- hzar.makeCline1DFreq(man_Q, scaling="none", tails="right")
# man_Q_model_none_left <- hzar.makeCline1DFreq(man_Q, scaling="none", tails="left")
# man_Q_model_none_mirror <- hzar.makeCline1DFreq(man_Q, scaling="none", tails="mirror")


### The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.
# man_Q_model_free_both <- hzar.model.addBoxReq(man_Q_model_free_both,-30,160)
# man_Q_model_free_left <- hzar.model.addBoxReq(man_Q_model_free_left,-30,160)
# man_S1_22a_model_free_right <- hzar.model.addBoxReq(man_Q_model_free_right,-30,160)
# man_S1_22a_model_free_mirror <- hzar.model.addBoxReq(man_Q_model_free_mirror,-30,160)
# man_S1_22a_model_free_none <- hzar.model.addBoxReq(man_Q_model_free_none,-30,160)
# 
# man_Q_model_fixed_both <- hzar.model.addBoxReq(man_Q_model_fixed_both,-30,160)
# man_Q_model_fixed_left <- hzar.model.addBoxReq(man_Q_model_fixed_left,-30,160)
# man_Q_model_fixed_right <- hzar.model.addBoxReq(man_Q_model_fixed_right,-30,160)
# man_Q_model_fixed_mirror <- hzar.model.addBoxReq(man_Q_model_fixed_mirror,-30,160)
# man_Q_model_fixed_none <- hzar.model.addBoxReq(man_Q_model_fixed_none,-30,160)
# 
# man_Q_model_none_both <- hzar.model.addBoxReq(man_Q_model_none_both,-30,160)
# man_Q_model_none_left <- hzar.model.addBoxReq(man_Q_model_none_left,-30,160)
# man_Q_model_none_right <- hzar.model.addBoxReq(man_Q_model_none_right,-30,160)
# man_Q_model_none_mirror <- hzar.model.addBoxReq(man_Q_model_none_mirror,-30,160)
man_Q_model_none_none <- hzar.model.addBoxReq(man_Q_model_none_none,-30,160)


###cline model fitting
##generate an hzar.fitRequest object suitable for hzar.doFit
# man_Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_free_both, man_Q, verbose=FALSE)
# man_Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_free_left, man_Q, verbose=FALSE)
# man_Q_model_free_rightFitR <-hzar.first.fitRequest.old.ML(model=man_Q_model_free_right, man_Q, verbose=FALSE)
# man_Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_free_mirror, man_Q, verbose=FALSE)
# man_Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_free_none, man_Q, verbose=FALSE)

# man_Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_fixed_both, man_Q, verbose=FALSE)
# man_Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_fixed_left, man_Q, verbose=FALSE)
# man_Q_model_fixed_rightFitR <-hzar.first.fitRequest.old.ML(model=man_Q_model_fixed_right, man_Q, verbose=FALSE)
# man_Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_fixed_mirror, man_Q, verbose=FALSE)
# man_Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_fixed_none, man_Q, verbose=FALSE)

# man_Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_none_both, man_Q, verbose=FALSE)
# man_Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_none_left, man_Q, verbose=FALSE)
# man_Q_model_none_rightFitR <-hzar.first.fitRequest.old.ML(model=man_Q_model_none_right, man_Q, verbose=FALSE)
# man_Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_none_mirror, man_Q, verbose=FALSE)
man_Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=man_Q_model_none_none, man_Q, verbose=FALSE)

##set mcmc chain length and burn in
# man_Q_model_free_bothFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_free_bothFitR$mcmcParam$burnin <- 5e5
# man_Q_model_free_leftFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_free_leftFitR$mcmcParam$burnin <- 5e5
# man_Q_model_free_rightFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_free_rightFitR$mcmcParam$burnin <- 5e5
# man_Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_free_mirrorFitR$mcmcParam$burnin <- 5e5
# man_Q_model_free_noneFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_free_noneFitR$mcmcParam$burnin <- 5e5

# man_Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_fixed_bothFitR$mcmcParam$burnin <- 5e5
# man_Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_fixed_leftFitR$mcmcParam$burnin <- 5e5
# man_Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_fixed_rightFitR$mcmcParam$burnin <- 5e5
# man_Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_fixed_mirrorFitR$mcmcParam$burnin <- 5e5
# man_Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

# man_Q_model_none_bothFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_none_bothFitR$mcmcParam$burnin <- 5e5
# man_Q_model_none_leftFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_none_leftFitR$mcmcParam$burnin <- 5e5
# man_Q_model_none_rightFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_none_rightFitR$mcmcParam$burnin <- 5e5
# man_Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e5
# man_Q_model_none_mirrorFitR$mcmcParam$burnin <- 5e5
man_Q_model_none_noneFitR$mcmcParam$chainLength <- 1e5
man_Q_model_none_noneFitR$mcmcParam$burnin <- 5e5

##Run the optimizer using the parameters listed in the hzar.fitRequest given.
# man_Q_model_free_bothFit <- hzar.doFit(man_Q_model_free_bothFitR)
# man_Q_model_free_leftFit <- hzar.doFit(man_Q_model_free_leftFitR)
# man_Q_model_free_rightFit <- hzar.doFit(man_Q_model_free_rightFitR)
# man_Q_model_free_mirrorFit <- hzar.doFit(man_Q_model_free_mirrorFitR)
# man_Q_model_free_noneFit <- hzar.doFit(man_Q_model_free_noneFitR)

# man_Q_model_fixed_bothFit <- hzar.doFit(man_Q_model_fixed_bothFitR)
# man_Q_model_fixed_leftFit <- hzar.doFit(man_Q_model_fixed_leftFitR)
# man_Q_model_fixed_rightFit <- hzar.doFit(man_Q_model_fixed_rightFitR)
# man_Q_model_fixed_mirrorFit <- hzar.doFit(man_Q_model_fixed_mirrorFitR)
# man_Q_model_fixed_noneFit <- hzar.doFit(man_Q_model_fixed_noneFitR)

# man_Q_model_none_bothFit <- hzar.doFit(man_Q_model_none_bothFitR)
# man_Q_model_none_leftFit <- hzar.doFit(man_Q_model_none_leftFitR)
# man_Q_model_none_rightFit <- hzar.doFit(man_Q_model_none_rightFitR)
# man_Q_model_none_mirrorFit <- hzar.doFit(man_Q_model_none_mirrorFitR)
man_Q_model_none_noneFit <- hzar.doFit(man_Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and return the mcmc data with an added a log likelihood column.
# par(mar=c(1,1,1,1))
# plot(hzar.mcmc.bindLL(man_Q_model_free_bothFit))
# plot(hzar.mcmc.bindLL(man_Q_model_free_leftFit))
# plot(hzar.mcmc.bindLL(man_Q_model_free_rightFit))
# plot(hzar.mcmc.bindLL(man_Q_model_free_mirrorFit))
# plot(hzar.mcmc.bindLL(man_Q_model_free_noneFit))

### group multiple fits of the same model and the same observation data into a single object
# man_Q_model_free_bothData <- hzar.dataGroup.add(man_Q_model_free_bothFit)
# man_Q_model_free_leftData <- hzar.dataGroup.add(man_Q_model_free_leftFit) 
# man_Q_model_free_rightData <- hzar.dataGroup.add(man_Q_model_free_rightFit)
# man_Q_model_free_mirrorData <- hzar.dataGroup.add(man_Q_model_free_mirrorFit)
# man_Q_model_free_noneData <- hzar.dataGroup.add(man_Q_model_free_noneFit)

# man_Q_model_fixed_bothData <- hzar.dataGroup.add(man_Q_model_fixed_bothFit)
# man_Q_model_fixed_leftData <- hzar.dataGroup.add(man_Q_model_fixed_leftFit) 
# man_Q_model_fixed_rightData <- hzar.dataGroup.add(man_Q_model_fixed_rightFit)
# man_Q_model_fixed_mirrorData <- hzar.dataGroup.add(man_Q_model_fixed_mirrorFit)
# man_Q_model_fixed_noneData <- hzar.dataGroup.add(man_Q_model_fixed_noneFit)

# man_Q_model_none_bothData <- hzar.dataGroup.add(man_Q_model_none_bothFit)
# man_Q_model_none_leftData <- hzar.dataGroup.add(man_Q_model_none_leftFit) 
# man_Q_model_none_rightData <- hzar.dataGroup.add(man_Q_model_none_rightFit)
# man_Q_model_none_mirrorData <- hzar.dataGroup.add(man_Q_model_none_mirrorFit)
man_Q_model_none_noneData <- hzar.dataGroup.add(man_Q_model_none_noneFit)


### Generate a hzar.dataGroup object representing a fit of the null model to a hzar.obsData object
man_Q_modelNull <- hzar.dataGroup.null(man_Q)

##make list of cline models and null models
man_Q_dGs <- list(cline_free_bothModel = man_Q_model_free_bothData,
                  cline_free_leftModel = man_Q_model_free_leftData,
                  cline_free_rightModel = man_Q_model_free_rightData,
                  cline_free_mirrorModel = man_Q_model_free_mirrorData,
                  cline_free_noneModel = man_Q_model_free_noneData,
                  cline_fixed_bothModel = man_Q_model_fixed_bothData,
                  cline_fixed_leftModel = man_Q_model_fixed_leftData,
                  cline_fixed_rightModel = man_Q_model_fixed_rightData,
                  cline_fixed_mirrorModel = man_Q_model_fixed_mirrorData,
                  cline_fixed_noneModel = man_Q_model_fixed_noneData,
                  cline_none_bothModel = man_Q_model_none_bothData,
                  cline_none_leftModel = man_Q_model_none_leftData,
                  cline_none_rightModel = man_Q_model_none_rightData,
                  cline_none_mirrorModel = man_Q_model_none_mirrorData,
                  cline_none_noneModel = man_Q_model_none_noneData,
                           nullModel = man_Q_modelNull)

### Collect optimizer output based on the same hzar.obsData object
man_Q_oDG <- hzar.make.obsDataGroup(man_Q_dGs)
 
### Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
man_Q_oDG <- hzar.copyModelLabels(man_Q_dGs,man_Q_oDG)
 
### Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(man_Q_oDG)
 
### Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(man_Q_oDG))
#                             AICc
# cline_free_bothModel    38.74405
# cline_free_leftModel    32.61257
# cline_free_rightModel   32.76373
# cline_free_mirrorModel  32.74017
# cline_free_noneModel    27.38734
# cline_fixed_bothModel   32.19472
# cline_fixed_leftModel   26.99131
# cline_fixed_rightModel  26.47339
# cline_fixed_mirrorModel 26.98859
# cline_fixed_noneModel   22.23681
# cline_none_bothModel    32.27241
# cline_none_leftModel    26.98895
# cline_none_rightModel   26.88876
# cline_none_mirrorModel  26.99258
# cline_none_noneModel    22.23561
# nullModel               56.57411

########## nbc ##########

### add an nSamples column with 1 for each because we are using samples not pops
nbc_admix$nSamples <- rep(1, nrow(nbc_admix))

### make data object of allele frequency data
nbc_Q <-
  hzar.doMolecularData1DPops(nbc_admix$dist,
                             nbc_admix$q,
                             nbc_admix$nSamples)

### Plot the associated observed frequency versus distance
hzar.plot.obsData(nbc_Q)


### Construct a clineMetaModel object for use with hzar.first.fitRequest.old.ML
# nbc_Q_model_free_both <- hzar.makeCline1DFreq(nbc_Q, scaling="free", tails="both")
# nbc_Q_model_free_none <- hzar.makeCline1DFreq(nbc_Q, scaling="free", tails="none")
# nbc_Q_model_free_right <- hzar.makeCline1DFreq(nbc_Q, scaling="free", tails="right")
# nbc_Q_model_free_left <- hzar.makeCline1DFreq(nbc_Q, scaling="free", tails="left")
# nbc_Q_model_free_mirror <- hzar.makeCline1DFreq(nbc_Q, scaling="free", tails="mirror")

# nbc_Q_model_fixed_both <- hzar.makeCline1DFreq(nbc_Q, scaling="fixed", tails="both")
# nbc_Q_model_fixed_none <- hzar.makeCline1DFreq(nbc_Q, scaling="fixed", tails="none")
# nbc_Q_model_fixed_right <- hzar.makeCline1DFreq(nbc_Q, scaling="fixed", tails="right")
# nbc_Q_model_fixed_left <- hzar.makeCline1DFreq(nbc_Q, scaling="fixed", tails="left")
# nbc_Q_model_fixed_mirror <- hzar.makeCline1DFreq(nbc_Q, scaling="fixed", tails="mirror")

# nbc_Q_model_none_both <- hzar.makeCline1DFreq(nbc_Q, scaling="none", tails="both")
# nbc_Q_model_none_none <- hzar.makeCline1DFreq(nbc_Q, scaling="none", tails="none")
nbc_Q_model_none_right <- hzar.makeCline1DFreq(nbc_Q, scaling="none", tails="right")
# nbc_Q_model_none_left <- hzar.makeCline1DFreq(nbc_Q, scaling="none", tails="left")
# nbc_Q_model_none_mirror <- hzar.makeCline1DFreq(nbc_Q, scaling="none", tails="mirror")


##The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.
# nbc_Q_model_free_both <- hzar.model.addBoxReq(nbc_Q_model_free_both,-100,650)
# nbc_Q_model_free_left <- hzar.model.addBoxReq(nbc_Q_model_free_left,-100,650)
# nbc_Q_model_free_right <- hzar.model.addBoxReq(nbc_Q_model_free_right,-100,650)
# nbc_Q_model_free_mirror <- hzar.model.addBoxReq(nbc_Q_model_free_mirror,-100,650)
# nbc_Q_model_free_none <- hzar.model.addBoxReq(nbc_Q_model_free_none,-100,650)

# nbc_Q_model_fixed_both <- hzar.model.addBoxReq(nbc_Q_model_fixed_both,-100,650)
# nbc_Q_model_fixed_left <- hzar.model.addBoxReq(nbc_Q_model_fixed_left,-100,650)
# nbc_Q_model_fixed_right <- hzar.model.addBoxReq(nbc_Q_model_fixed_right,-100,650)
# nbc_Q_model_fixed_mirror <- hzar.model.addBoxReq(nbc_Q_model_fixed_mirror,-100,650)
# nbc_Q_model_fixed_none <- hzar.model.addBoxReq(nbc_Q_model_fixed_none,-100,650)

# nbc_Q_model_none_both <- hzar.model.addBoxReq(nbc_Q_model_none_both,-100,650)
# nbc_Q_model_none_left <- hzar.model.addBoxReq(nbc_Q_model_none_left,-100,650)
nbc_Q_model_none_right <- hzar.model.addBoxReq(nbc_Q_model_none_right,-100,650)
# nbc_Q_model_none_mirror <- hzar.model.addBoxReq(nbc_Q_model_none_mirror,-100,650)
# nbc_Q_model_none_none <- hzar.model.addBoxReq(nbc_Q_model_none_none,-100,650)


#### cline model fitting
### generate an hzar.fitRequest object suitable for hzar.doFit
# nbc_Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_free_both, nbc_Q, verbose=FALSE)
# nbc_Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_free_left, nbc_Q, verbose=FALSE)
# nbc_Q_model_free_rightFitR <-hzar.first.fitRequest.old.ML(model=nbc_Q_model_free_right, nbc_Q, verbose=FALSE)
# nbc_Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_free_mirror, nbc_Q, verbose=FALSE)
# nbc_Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_free_none, nbc_Q, verbose=FALSE)

# nbc_Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_fixed_both, nbc_Q, verbose=FALSE)
# nbc_Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_fixed_left, nbc_Q, verbose=FALSE)
# nbc_Q_model_fixed_rightFitR <-hzar.first.fitRequest.old.ML(model=nbc_Q_model_fixed_right, nbc_Q, verbose=FALSE)
# nbc_Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_fixed_mirror, nbc_Q, verbose=FALSE)
# nbc_Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_fixed_none, nbc_Q, verbose=FALSE)

# nbc_Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_none_both, nbc_Q, verbose=FALSE)
# nbc_Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_none_left, nbc_Q, verbose=FALSE)
nbc_Q_model_none_rightFitR <-hzar.first.fitRequest.old.ML(model=nbc_Q_model_none_right, nbc_Q, verbose=FALSE)
# nbc_Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_none_mirror, nbc_Q, verbose=FALSE)
# nbc_Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_model_none_none, nbc_Q, verbose=FALSE)

##set mcmc chain length and burn in
# nbc_Q_model_free_bothFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_free_bothFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_free_leftFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_free_leftFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_free_rightFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_free_rightFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_free_mirrorFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_free_noneFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_free_noneFitR$mcmcParam$burnin <- 5e5

# nbc_Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_fixed_bothFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_fixed_leftFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_fixed_rightFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_fixed_mirrorFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

# nbc_Q_model_none_bothFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_none_bothFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_none_leftFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_none_leftFitR$mcmcParam$burnin <- 5e5
nbc_Q_model_none_rightFitR$mcmcParam$chainLength <- 1e5
nbc_Q_model_none_rightFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_none_mirrorFitR$mcmcParam$burnin <- 5e5
# nbc_Q_model_none_noneFitR$mcmcParam$chainLength <- 1e5
# nbc_Q_model_none_noneFitR$mcmcParam$burnin <- 5e5

##Run the optimizer using the parameters listed in the hzar.fitRequest given.
# nbc_Q_model_free_bothFit <- hzar.doFit(nbc_Q_model_free_bothFitR)
# nbc_Q_model_free_leftFit <- hzar.doFit(nbc_Q_model_free_leftFitR)
# nbc_Q_model_free_rightFit <- hzar.doFit(nbc_Q_model_free_rightFitR)
# nbc_Q_model_free_mirrorFit <- hzar.doFit(nbc_Q_model_free_mirrorFitR)
# nbc_Q_model_free_noneFit <- hzar.doFit(nbc_Q_model_free_noneFitR)

# nbc_Q_model_fixed_bothFit <- hzar.doFit(nbc_Q_model_fixed_bothFitR)
# nbc_Q_model_fixed_leftFit <- hzar.doFit(nbc_Q_model_fixed_leftFitR)
# nbc_Q_model_fixed_rightFit <- hzar.doFit(nbc_Q_model_fixed_rightFitR)
# nbc_Q_model_fixed_mirrorFit <- hzar.doFit(nbc_Q_model_fixed_mirrorFitR)
# nbc_Q_model_fixed_noneFit <- hzar.doFit(nbc_Q_model_fixed_noneFitR)

# nbc_Q_model_none_bothFit <- hzar.doFit(nbc_Q_model_none_bothFitR)
# nbc_Q_model_none_leftFit <- hzar.doFit(nbc_Q_model_none_leftFitR)
nbc_Q_model_none_rightFit <- hzar.doFit(nbc_Q_model_none_rightFitR)
# nbc_Q_model_none_mirrorFit <- hzar.doFit(nbc_Q_model_none_mirrorFitR)
# nbc_Q_model_none_noneFit <- hzar.doFit(nbc_Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and return the mcmc data with an added a log likelihood column.
# par(mar=c(1,1,1,1))
# plot(hzar.mcmc.bindLL(nbc_Q_model_free_bothFit))
# plot(hzar.mcmc.bindLL(nbc_Q_model_free_leftFit))
# plot(hzar.mcmc.bindLL(nbc_Q_model_free_rightFit))
# plot(hzar.mcmc.bindLL(nbc_Q_model_free_mirrorFit))
# plot(hzar.mcmc.bindLL(nbc_Q_model_free_noneFit))
# 
##group multiple fits of the same model and the same observation data into a single object
# nbc_Q_model_free_bothData <- hzar.dataGroup.add(nbc_Q_model_free_bothFit)
# nbc_Q_model_free_leftData <- hzar.dataGroup.add(nbc_Q_model_free_leftFit) 
# nbc_Q_model_free_rightData <- hzar.dataGroup.add(nbc_Q_model_free_rightFit)
# nbc_Q_model_free_mirrorData <- hzar.dataGroup.add(nbc_Q_model_free_mirrorFit)
# nbc_Q_model_free_noneData <- hzar.dataGroup.add(nbc_Q_model_free_noneFit)

# nbc_Q_model_fixed_bothData <- hzar.dataGroup.add(nbc_Q_model_fixed_bothFit)
# nbc_Q_model_fixed_leftData <- hzar.dataGroup.add(nbc_Q_model_fixed_leftFit) 
# nbc_Q_model_fixed_rightData <- hzar.dataGroup.add(nbc_Q_model_fixed_rightFit)
# nbc_Q_model_fixed_mirrorData <- hzar.dataGroup.add(nbc_Q_model_fixed_mirrorFit)
# nbc_Q_model_fixed_noneData <- hzar.dataGroup.add(nbc_Q_model_fixed_noneFit)

# nbc_Q_model_none_bothData <- hzar.dataGroup.add(nbc_Q_model_none_bothFit)
# nbc_Q_model_none_leftData <- hzar.dataGroup.add(nbc_Q_model_none_leftFit) 
nbc_Q_model_none_rightData <- hzar.dataGroup.add(nbc_Q_model_none_rightFit)
# nbc_Q_model_none_mirrorData <- hzar.dataGroup.add(nbc_Q_model_none_mirrorFit)
# nbc_Q_model_none_noneData <- hzar.dataGroup.add(nbc_Q_model_none_noneFit)


### Generate a hzar.dataGroup object representing a fit of the null model to a hzar.obsData object
nbc_Q_modelNull <- hzar.dataGroup.null(nbc_Q)

##make list of cline models and null models
nbc_Q_dGs <- list(cline_free_bothModel = nbc_Q_model_free_bothData,
                  cline_free_leftModel = nbc_Q_model_free_leftData,
                  cline_free_rightModel = nbc_Q_model_free_rightData,
                  cline_free_mirrorModel = nbc_Q_model_free_mirrorData,
                  cline_free_noneModel = nbc_Q_model_free_noneData,
                  cline_fixed_bothModel = nbc_Q_model_fixed_bothData,
                  cline_fixed_leftModel = nbc_Q_model_fixed_leftData,
                  cline_fixed_rightModel = nbc_Q_model_fixed_rightData,
                  cline_fixed_mirrorModel = nbc_Q_model_fixed_mirrorData,
                  cline_fixed_noneModel = nbc_Q_model_fixed_noneData,
                  cline_none_bothModel = nbc_Q_model_none_bothData,
                  cline_none_leftModel = nbc_Q_model_none_leftData,
                  cline_none_rightModel = nbc_Q_model_none_rightData,
                  cline_none_mirrorModel = nbc_Q_model_none_mirrorData,
                  cline_none_noneModel = nbc_Q_model_none_noneData,
                  nullModel = nbc_Q_modelNull)

##Collect optimizer output based on the same hzar.obsData object
nbc_Q_oDG <- hzar.make.obsDataGroup(nbc_Q_dGs)

##Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
nbc_Q_oDG <- hzar.copyModelLabels(nbc_Q_dGs,nbc_Q_oDG)

##Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(nbc_Q_oDG)

##Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(nbc_Q_oDG))
# AICc
# cline_free_bothModel    37.84759
# cline_free_leftModel    32.93323
# cline_free_rightModel   32.77541
# cline_free_mirrorModel  33.01367
# cline_free_noneModel    28.41212
# cline_fixed_bothModel   32.92816
# cline_fixed_leftModel   37.27875
# cline_fixed_rightModel  28.50004
# cline_fixed_mirrorModel 29.98919
# cline_fixed_noneModel   32.87428
# cline_none_bothModel    32.53676
# cline_none_leftModel    43.87432
# cline_none_rightModel   28.25540
# cline_none_mirrorModel  29.36315
# cline_none_noneModel    39.47224
# nullModel               78.41916





########## caor ##########

### add an nSamples column with 1 for each because we are using samples not pops
caor_admix$nSamples <- rep(1, nrow(caor_admix))

### make data object of allele frequency data
caor_Q <-
  hzar.doMolecularData1DPops(caor_admix$dist,
                             caor_admix$q,
                             caor_admix$nSamples)

### Plot the associated observed frequency versus distance
hzar.plot.obsData(caor_Q)


### Construct a clineMetaModel object for use with hzar.first.fitRequest.old.ML
 # caor_Q_model_free_both <- hzar.makeCline1DFreq(caor_Q, scaling="free", tails="both")
 # caor_Q_model_free_none <- hzar.makeCline1DFreq(caor_Q, scaling="free", tails="none")
 # caor_Q_model_free_right <- hzar.makeCline1DFreq(caor_Q, scaling="free", tails="right")
 # caor_Q_model_free_left <- hzar.makeCline1DFreq(caor_Q, scaling="free", tails="left")
 # caor_Q_model_free_mirror <- hzar.makeCline1DFreq(caor_Q, scaling="free", tails="mirror")
 
#  caor_Q_model_fixed_both <- hzar.makeCline1DFreq(caor_Q, scaling="fixed", tails="both")
# caor_Q_model_fixed_none <- hzar.makeCline1DFreq(caor_Q, scaling="fixed", tails="none")
#  caor_Q_model_fixed_right <- hzar.makeCline1DFreq(caor_Q, scaling="fixed", tails="right")
#  caor_Q_model_fixed_left <- hzar.makeCline1DFreq(caor_Q, scaling="fixed", tails="left")
#  caor_Q_model_fixed_mirror <- hzar.makeCline1DFreq(caor_Q, scaling="fixed", tails="mirror")
 
 # caor_Q_model_none_both <- hzar.makeCline1DFreq(caor_Q, scaling="none", tails="both")
 caor_Q_model_none_none <- hzar.makeCline1DFreq(caor_Q, scaling="none", tails="none")
 # caor_Q_model_none_right <- hzar.makeCline1DFreq(caor_Q, scaling="none", tails="right")
 # caor_Q_model_none_left <- hzar.makeCline1DFreq(caor_Q, scaling="none", tails="left")
 # caor_Q_model_none_mirror <- hzar.makeCline1DFreq(caor_Q, scaling="none", tails="mirror")
 

## The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.
 # caor_Q_model_free_both <- hzar.model.addBoxReq(caor_Q_model_free_both,-100,750)
 # caor_Q_model_free_left <- hzar.model.addBoxReq(caor_Q_model_free_left,-100,750)
 # caor_Q_model_free_right <- hzar.model.addBoxReq(caor_Q_model_free_right,-100,750)
 # caor_Q_model_free_mirror <- hzar.model.addBoxReq(caor_Q_model_free_mirror,-100,750)
 # caor_Q_model_free_none <- hzar.model.addBoxReq(caor_Q_model_free_none,-100,750)
 
#  caor_Q_model_fixed_both <- hzar.model.addBoxReq(caor_Q_model_fixed_both,-100,750)
#  caor_Q_model_fixed_left <- hzar.model.addBoxReq(caor_Q_model_fixed_left,-100,750)
#  caor_Q_model_fixed_right <- hzar.model.addBoxReq(caor_Q_model_fixed_right,-100,750)
#  caor_Q_model_fixed_mirror <- hzar.model.addBoxReq(caor_Q_model_fixed_mirror,-100,750)
# caor_Q_model_fixed_none <- hzar.model.addBoxReq(caor_Q_model_fixed_none,-100,750)

 # caor_Q_model_none_both <- hzar.model.addBoxReq(caor_Q_model_none_both,-100,750)
 # caor_Q_model_none_left <- hzar.model.addBoxReq(caor_Q_model_none_left,-100,750)
 # caor_Q_model_none_right <- hzar.model.addBoxReq(caor_Q_model_none_right,-100,750)
 # caor_Q_model_none_mirror <- hzar.model.addBoxReq(caor_Q_model_none_mirror,-100,750)
caor_Q_model_none_none <- hzar.model.addBoxReq(caor_Q_model_none_none,-100,750)


#### cline model fitting
## generate an hzar.fitRequest object suitable for hzar.doFit
# caor_Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_free_both, caor_Q, verbose=FALSE)
# caor_Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_free_left, caor_Q, verbose=FALSE)
# caor_Q_model_free_rightFitR <-hzar.first.fitRequest.old.ML(model=caor_Q_model_free_right, caor_Q, verbose=FALSE)
# caor_Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_free_mirror, caor_Q, verbose=FALSE)
# caor_Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_free_none, caor_Q, verbose=FALSE)
 
# caor_Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_fixed_both, caor_Q, verbose=FALSE)
# caor_Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_fixed_left, caor_Q, verbose=FALSE)
# caor_Q_model_fixed_rightFitR <-hzar.first.fitRequest.old.ML(model=caor_Q_model_fixed_right, caor_Q, verbose=FALSE)
# caor_Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_fixed_mirror, caor_Q, verbose=FALSE)
# caor_Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_fixed_none, caor_Q, verbose=FALSE)

# caor_Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_none_both, caor_Q, verbose=FALSE)
# caor_Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_none_left, caor_Q, verbose=FALSE)
# caor_Q_model_none_rightFitR <-hzar.first.fitRequest.old.ML(model=caor_Q_model_none_right, caor_Q, verbose=FALSE)
# caor_Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_none_mirror, caor_Q, verbose=FALSE)
caor_Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_model_none_none, caor_Q, verbose=FALSE)
 
 ## set mcmc chain length and burn in
# caor_Q_model_free_bothFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_free_bothFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_free_leftFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_free_leftFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_free_rightFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_free_rightFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_free_mirrorFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_free_noneFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_free_noneFitR$mcmcParam$burnin <- 5e5

# caor_Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_fixed_bothFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_fixed_leftFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_fixed_rightFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_fixed_mirrorFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

# caor_Q_model_none_bothFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_none_bothFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_none_leftFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_none_leftFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_none_rightFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_none_rightFitR$mcmcParam$burnin <- 5e5
# caor_Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e5
# caor_Q_model_none_mirrorFitR$mcmcParam$burnin <- 5e5
caor_Q_model_none_noneFitR$mcmcParam$chainLength <- 1e5
caor_Q_model_none_noneFitR$mcmcParam$burnin <- 5e5
 
# ### Run the optimizer using the parameters listed in the hzar.fitRequest given.
# caor_Q_model_free_bothFit <- hzar.doFit(caor_Q_model_free_bothFitR)
# caor_Q_model_free_leftFit <- hzar.doFit(caor_Q_model_free_leftFitR)
# caor_Q_model_free_rightFit <- hzar.doFit(caor_Q_model_free_rightFitR)
# caor_Q_model_free_mirrorFit <- hzar.doFit(caor_Q_model_free_mirrorFitR)
# caor_Q_model_free_noneFit <- hzar.doFit(caor_Q_model_free_noneFitR)

# caor_Q_model_fixed_bothFit <- hzar.doFit(caor_Q_model_fixed_bothFitR)
# caor_Q_model_fixed_leftFit <- hzar.doFit(caor_Q_model_fixed_leftFitR)
# caor_Q_model_fixed_rightFit <- hzar.doFit(caor_Q_model_fixed_rightFitR)
# caor_Q_model_fixed_mirrorFit <- hzar.doFit(caor_Q_model_fixed_mirrorFitR)
# caor_Q_model_fixed_noneFit <- hzar.doFit(caor_Q_model_fixed_noneFitR)

# caor_Q_model_none_bothFit <- hzar.doFit(caor_Q_model_none_bothFitR)
# caor_Q_model_none_leftFit <- hzar.doFit(caor_Q_model_none_leftFitR)
# caor_Q_model_none_rightFit <- hzar.doFit(caor_Q_model_none_rightFitR)
# caor_Q_model_none_mirrorFit <- hzar.doFit(caor_Q_model_none_mirrorFitR)
caor_Q_model_none_noneFit <- hzar.doFit(caor_Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and return the mcmc data with an added a log likelihood column.
# par(mar=c(1,1,1,1))
# plot(hzar.mcmc.bindLL(caor_Q_model_free_bothFit))
# plot(hzar.mcmc.bindLL(caor_Q_model_free_leftFit))
# plot(hzar.mcmc.bindLL(caor_Q_model_free_rightFit))
# plot(hzar.mcmc.bindLL(caor_Q_model_free_mirrorFit))
# plot(hzar.mcmc.bindLL(caor_Q_model_free_noneFit))
# 
##group multiple fits of the same model and the same observation data into a single object
# caor_Q_model_free_bothData <- hzar.dataGroup.add(caor_Q_model_free_bothFit)
# caor_Q_model_free_leftData <- hzar.dataGroup.add(caor_Q_model_free_leftFit) 
# caor_Q_model_free_rightData <- hzar.dataGroup.add(caor_Q_model_free_rightFit)
# caor_Q_model_free_mirrorData <- hzar.dataGroup.add(caor_Q_model_free_mirrorFit)
# caor_Q_model_free_noneData <- hzar.dataGroup.add(caor_Q_model_free_noneFit)

# caor_Q_model_fixed_bothData <- hzar.dataGroup.add(caor_Q_model_fixed_bothFit)
# caor_Q_model_fixed_leftData <- hzar.dataGroup.add(caor_Q_model_fixed_leftFit) 
# caor_Q_model_fixed_rightData <- hzar.dataGroup.add(caor_Q_model_fixed_rightFit)
# caor_Q_model_fixed_mirrorData <- hzar.dataGroup.add(caor_Q_model_fixed_mirrorFit)
# caor_Q_model_fixed_noneData <- hzar.dataGroup.add(caor_Q_model_fixed_noneFit)

# caor_Q_model_none_bothData <- hzar.dataGroup.add(caor_Q_model_none_bothFit)
# caor_Q_model_none_leftData <- hzar.dataGroup.add(caor_Q_model_none_leftFit) 
# caor_Q_model_none_rightData <- hzar.dataGroup.add(caor_Q_model_none_rightFit)
# caor_Q_model_none_mirrorData <- hzar.dataGroup.add(caor_Q_model_none_mirrorFit)
caor_Q_model_none_noneData <- hzar.dataGroup.add(caor_Q_model_none_noneFit)

### Generate a hzar.dataGroup object representing a fit of the null model to a hzar.obsData object
caor_Q_modelNull <- hzar.dataGroup.null(caor_Q)

##make list of cline models and null models
caor_Q_dGs <- list(cline_free_bothModel = caor_Q_model_free_bothData,
                  cline_free_leftModel = caor_Q_model_free_leftData,
                  cline_free_rightModel = caor_Q_model_free_rightData,
                  cline_free_mirrorModel = caor_Q_model_free_mirrorData,
                  cline_free_noneModel = caor_Q_model_free_noneData,
                  cline_fixed_bothModel = caor_Q_model_fixed_bothData,
                  cline_fixed_leftModel = caor_Q_model_fixed_leftData,
                  cline_fixed_rightModel = caor_Q_model_fixed_rightData,
                  cline_fixed_mirrorModel = caor_Q_model_fixed_mirrorData,
                  cline_fixed_noneModel = caor_Q_model_fixed_noneData,
                  cline_none_bothModel = caor_Q_model_none_bothData,
                  cline_none_leftModel = caor_Q_model_none_leftData,
                  cline_none_rightModel = caor_Q_model_none_rightData,
                  cline_none_mirrorModel = caor_Q_model_none_mirrorData,
                  cline_none_noneModel = caor_Q_model_none_noneData,
                  nullModel = caor_Q_modelNull)

##Collect optimizer output based on the same hzar.obsData object
caor_Q_oDG <- hzar.make.obsDataGroup(caor_Q_dGs)

##Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
caor_Q_oDG <- hzar.copyModelLabels(caor_Q_dGs,caor_Q_oDG)

##Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(caor_Q_oDG)

##Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(caor_Q_oDG))
# AICc
# cline_free_bothModel    36.72147
# cline_free_leftModel    31.79695
# cline_free_rightModel   32.07694
# cline_free_mirrorModel  31.86228
# cline_free_noneModel    27.35920
# cline_fixed_bothModel   31.57778
# cline_fixed_leftModel   27.14017
# cline_fixed_rightModel  27.14300
# cline_fixed_mirrorModel 27.13850
# cline_fixed_noneModel   22.86554
# cline_none_bothModel    31.56564
# cline_none_leftModel    27.14127
# cline_none_rightModel   27.13960
# cline_none_mirrorModel  27.13947
# cline_none_noneModel    22.86538
# nullModel               79.89871




### make em purty, get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width
min(hzar.AICc.hzar.obsDataGroup(nbc_Q_oDG))
# nbc = cline_none_rightModel   28.25540

min(hzar.AICc.hzar.obsDataGroup(man_Q_oDG))
# man = cline_none_noneModel    22.23561 

min(hzar.AICc.hzar.obsDataGroup(caor_Q_oDG))
# caor = # cline_none_noneModel    22.86538

pdf(file="man_q_cline.pdf")
hzar.plot.fzCline(man_Q_model_none_noneData, fzCol = rgb(68/255, 1/255, 84/255, 0.3), pch = 20, col=rgb(68/255, 1/255, 84/255, 0.7), xlab = "Distance (km)", ylab = "Q value")
dev.off()
print(man_Q_model_none_noneData$ML.cline$param.free$center)
#[1] 0.8665297
print(man_Q_model_none_noneData$ML.cline$param.free$width)
#[1] 20.0455
print(hzar.getLLCutParam(man_Q_model_none_noneData,c("center","width")))
# center2LLLow center2LLHigh width2LLLow width2LLHigh
#   -8.183703      5.566193    10.53658     47.68713

pdf(file="caor_q_cline.pdf")
hzar.plot.fzCline(caor_Q_model_none_noneData,  fzCol = rgb(33/255, 145/255, 140/255, 0.3), pch = 20, col=rgb(33/255, 145/255, 140/255, 0.7), xlab = "Distance (km)", ylab = "Q value")
dev.off()
print(caor_Q_model_none_noneData$ML.cline$param.free$center)
#[1] -41.10436
print(caor_Q_model_none_noneData$ML.cline$param.free$width)
#[1] 138.0059
print(hzar.getLLCutParam(caor_Q_model_none_noneData,c("center","width")))
# center2LLLow center2LLHigh width2LLLow width2LLHigh
# 1    -69.69375      2.401899    79.31715     228.5547


pdf(file="nbc_q_cline.pdf")
hzar.plot.fzCline(nbc_Q_model_none_rightData, fzCol = rgb(253/255, 231/255, 37/255, 0.3), pch = 20, col=rgb(253/255, 231/255, 37/255, 0.7), xlab = "Distance (km)", ylab = "Q value" )
dev.off()
print(nbc_Q_model_none_rightData$ML.cline$param.free$center)
#[1]  -0.2135945
print(nbc_Q_model_none_rightData$ML.cline$param.free$width)
#[1] 27.59858
print(hzar.getLLCutParam(nbc_Q_model_none_rightData,c("center","width")))
# center2LLLow center2LLHigh width2LLLow width2LLHigh
#   -5.608157      13.69607    10.18367     63.09187


############## Plot all the clines together ##############
### remake all the models with the points subtracted around the centre
man_centered <- man_admix
man_centered$distance <- man_centered$dist-(man_Q_model_none_noneData$ML.cline$param.free$center)
nbc_centered <- nbc_admix
nbc_centered$distance <- nbc_admix$dist-(nbc_Q_model_none_rightData$ML.cline$param.free$center)
caor_centered <- caor_admix
caor_centered$distance <- caor_admix$dist-(caor_Q_model_none_noneData$ML.cline$param.free$center)

man_Q_centered <-
  hzar.doMolecularData1DPops(man_centered$distance,
                             man_centered$q,
                             man_centered$nSamples)
nbc_Q_centered <-
  hzar.doMolecularData1DPops(nbc_centered$distance,
                             nbc_centered$q,
                             nbc_centered$nSamples)
caor_Q_centered <-
  hzar.doMolecularData1DPops(caor_centered$distance,
                             caor_centered$q,
                             caor_centered$nSamples)

hzar.plot.obsData(man_Q_centered)
hzar.plot.obsData(nbc_Q_centered)
hzar.plot.obsData(caor_Q_centered)

### make ClineMetaModel object
man_Q_centered_model_none_none<- hzar.makeCline1DFreq(man_Q_centered, scaling="none", tails="none")
nbc_Q_centered_model_none_right <- hzar.makeCline1DFreq(nbc_Q_centered, scaling="none", tails="right")
caor_Q_centered_model_none_none <- hzar.makeCline1DFreq(caor_Q_centered, scaling="none", tails="none")

### set BoxReq
man_Q_centered_model_none_none <- hzar.model.addBoxReq(man_Q_centered_model_none_none,-100,60)
nbc_Q_centered_model_none_right <- hzar.model.addBoxReq(nbc_Q_centered_model_none_right,-100,550)
caor_Q_centered_model_none_none <- hzar.model.addBoxReq(caor_Q_centered_model_none_none,-325,425)

### generate hzar fitRequests
man_Q_centered_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=man_Q_centered_model_none_none, man_Q_centered, verbose=FALSE)
nbc_Q_centered_model_none_rightFitR <- hzar.first.fitRequest.old.ML(model=nbc_Q_centered_model_none_right, nbc_Q_centered, verbose=FALSE)
caor_Q_centered_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=caor_Q_centered_model_none_none, caor_Q_centered, verbose=FALSE)

### set mcmc chainLength and burnin
man_Q_centered_model_none_noneFitR$mcmcParam$chainLength <- 1e5
man_Q_centered_model_none_noneFitR$mcmcParam$burnin <- 5e5
nbc_Q_centered_model_none_rightFitR$mcmcParam$chainLength <- 1e5
nbc_Q_centered_model_none_rightFitR$mcmcParam$burnin <- 5e5
caor_Q_centered_model_none_noneFitR$mcmcParam$chainLength <- 1e5
caor_Q_centered_model_none_noneFitR$mcmcParam$burnin <- 5e5

### Run the optimizer using the parameters listed in the hzar.fitRequest given.
man_Q_centered_model_none_noneFit <- hzar.doFit(man_Q_centered_model_none_noneFitR)
nbc_Q_centered_model_none_rightFit <- hzar.doFit(nbc_Q_centered_model_none_rightFitR)
caor_Q_centered_model_none_noneFit <- hzar.doFit(caor_Q_centered_model_none_noneFitR)

### group multiple fits of the same model and the same observation data into a single object
man_Q_centered_model_none_noneData <- hzar.dataGroup.add(man_Q_centered_model_none_noneFit)
nbc_Q_centered_model_none_rightData <- hzar.dataGroup.add(nbc_Q_centered_model_none_rightFit)
caor_Q_centered_model_none_noneData <- hzar.dataGroup.add(caor_Q_centered_model_none_noneFit)

### get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width
print(hzar.getLLCutParam(man_Q_centered_model_none_noneData,c("center","width")))
print(hzar.getLLCutParam(nbc_Q_centered_model_none_rightData,c("center","width")))
print(hzar.getLLCutParam(caor_Q_centered_model_none_noneData,c("center","width")))


### MANNY BOEHM CODE
# hzar.overPlot.fzCline doesn't contain any code for creating plots :(
# so, let's extract info from 'hzar.dataGroup' objects in the same way
# that hzar.plot.fzCline does

# --------------------
# 95% CI

# hzar.plot.fzCline() fits a line (and the 95% CI) to the points by first 
# generating a series of 107 x coordinates bounded by the current plot window
# later, the y coordinates are estimated by accessing 
# the function describing the cline in the dataGroup object
## sometimes this dies and makes the polygons real small but if I rerun this code a couple times they go backto normal
# the numbers are from hzar, no relation to my observations being 108 in caor
xSeries <- seq(from = par("usr")[1], to = par("usr")[2], 
               length.out = 109)
if (par("xaxs") == "r") 
  xSeries <- xSeries[2:108]

# hzar.plot.fzCline() plots the 95% CI using graphics::polygon() 
# you can customize the colours here
man_fzCline = hzar.getCredParamRed(man_Q_centered_model_none_noneData)
man_fzCoor <- man_fzCline$fzCline(xSeries)
nbc_fzCline = hzar.getCredParamRed(nbc_Q_centered_model_none_rightData)
nbc_fzCoor <- nbc_fzCline$fzCline(xSeries)
caor_fzCline = hzar.getCredParamRed(caor_Q_centered_model_none_noneData)
caor_fzCoor <- caor_fzCline$fzCline(xSeries)

### plot the polygons and lines
# the stats::line() function uses the x and y coordinates from above to draw the line

pdf(file = "centered_clines.pdf")

# first create an empty plot bounded by the biggest dataset
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent")

#add polygons and lines
polygon(x = c(caor_fzCoor$x, rev(caor_fzCoor$x)), y = c(caor_fzCoor$yMin, rev(caor_fzCoor$yMax)), border = rgb(33/255, 145/255, 140/255, 0.9), col=rgb(33/255, 145/255, 140/255, 0.8))
lines(x = xSeries, y = caor_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries), col = rgb(33/255, 145/255, 140/255), lwd=2)
polygon(x = c(nbc_fzCoor$x, rev(nbc_fzCoor$x)), y = c(nbc_fzCoor$yMin, rev(nbc_fzCoor$yMax)), border = rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.7))
lines(x = xSeries, y = nbc_Q_centered_model_none_rightData$ML.cline$clineFunc(xSeries), col=rgb(253/255, 231/255, 37/255), lwd=2)
polygon(x = c(man_fzCoor$x, rev(man_fzCoor$x)), y = c(man_fzCoor$yMin, rev(man_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.75), col=rgb(68/255, 1/255, 84/255, 0.3))
lines(x = xSeries, y = man_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries),  col=rgb(68/255, 1/255, 84/255), lwd=2)
dev.off()   

# do a zoomed in version on just the centers
pdf(file = "centered_clines_zoom.pdf")

# set x axis limits at -200 and 200
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent", xlim=c(-200,200))

#add polygons and lines
polygon(x = c(caor_fzCoor$x, rev(caor_fzCoor$x)), y = c(caor_fzCoor$yMin, rev(caor_fzCoor$yMax)), border = rgb(33/255, 145/255, 140/255, 0.9), col=rgb(33/255, 145/255, 140/255, 0.6))
lines(x = xSeries, y = caor_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries), col = rgb(33/255, 145/255, 140/255), lwd=4)
polygon(x = c(nbc_fzCoor$x, rev(nbc_fzCoor$x)), y = c(nbc_fzCoor$yMin, rev(nbc_fzCoor$yMax)), border = rgb(253/255, 231/255, 37/255, 0.9), col=rgb(253/255, 231/255, 37/255, 0.5))
lines(x = xSeries, y = nbc_Q_centered_model_none_rightData$ML.cline$clineFunc(xSeries), col=rgb(253/255, 231/255, 37/255), lwd=4)
polygon(x = c(man_fzCoor$x, rev(man_fzCoor$x)), y = c(man_fzCoor$yMin, rev(man_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.75), col=rgb(68/255, 1/255, 84/255, 0.25))
lines(x = xSeries, y = man_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries),  col=rgb(68/255, 1/255, 84/255), lwd=4)
dev.off()

# --------------------
# points

# hzar.plot.fzCline() finds points to plot via hzar.plot.obsData() which extracts
# points from dataGroup objects using hzar.extract.obsData()
pts_man <- 
  hzar.extract.obsData(man_Q_centered_model_none_noneData)$frame %>% 
  dplyr::select(1:2)
pts_man$site <- rep("man", nrow(pts_man))
pts_nbc <- 
  hzar.extract.obsData(nbc_Q_centered_model_none_rightData)$frame %>% 
  dplyr::select(1:2)
pts_nbc$site <- rep("nbc", nrow(pts_nbc))
pts_caor <- 
  hzar.extract.obsData(caor_Q_centered_model_none_noneData)$frame %>% 
  dplyr::select(1:2)
pts_caor$site <- rep("caor", nrow(pts_caor))
all_pts <- full_join(pts_man, pts_nbc) %>%
  full_join(.,pts_caor)

# save plot
pdf(file = "centered_samples.pdf")

# first create an empty plot bounded by the biggest dataset
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent")

points(pts_caor, col=rgb(33/255, 145/255, 140/255, 0.7), pch=20, cex = 1.25)
points(pts_nbc, col=rgb(253/255, 231/255, 37/255, 0.7), pch=20, cex = 1.25)
points(pts_man, col=rgb(68/255, 1/255, 84/255, 0.3), pch=20, cex = 1.25)

dev.off()


#plot man
pdf(file="man_cline.pdf")
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent", xlim=c(-200,200))points(pts_man, col=rgb(68/255, 1/255, 84/255, 0.3), pch=20, cex = 1.25)
polygon(x = c(man_fzCoor$x, rev(man_fzCoor$x)), y = c(man_fzCoor$yMin, rev(man_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.7), col=rgb(68/255, 1/255, 84/255, 0.3))
lines(x = xSeries, y = man_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries),  col=rgb(68/255, 1/255, 84/255), lwd=2)
dev.off()

pdf(file="man_cline_stdscale.pdf")
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent", xlim=c(-300,400))
points(pts_man, col=rgb(68/255, 1/255, 84/255, 0.3), pch=20, cex = 1.25)
polygon(x = c(man_fzCoor$x, rev(man_fzCoor$x)), y = c(man_fzCoor$yMin, rev(man_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.7), col=rgb(68/255, 1/255, 84/255, 0.3))
lines(x = xSeries, y = man_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries),  col=rgb(68/255, 1/255, 84/255), lwd=2)
dev.off()

#plot nbc
pdf(file = "nbc_cline.pdf")
hzar.plot.obsData(nbc_Q_centered_model_none_rightData, col = "transparent")
points(pts_nbc, col=rgb(253/255, 231/255, 37/255, 0.5), pch=20, cex = 1.25)
polygon(x = c(nbc_fzCoor$x, rev(nbc_fzCoor$x)), y = c(nbc_fzCoor$yMin, rev(nbc_fzCoor$yMax)), border = rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.5))
lines(x = xSeries, y = nbc_Q_centered_model_none_rightData$ML.cline$clineFunc(xSeries), col=rgb(253/255, 231/255, 37/255), lwd=2)
dev.off()

#plot nbc
pdf(file = "nbc_cline_stdscale.pdf")
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent", xlim=c(-300,400))
points(pts_nbc, col=rgb(253/255, 231/255, 37/255, 0.5), pch=20, cex = 1.25)
polygon(x = c(nbc_fzCoor$x, rev(nbc_fzCoor$x)), y = c(nbc_fzCoor$yMin, rev(nbc_fzCoor$yMax)), border = rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.5))
lines(x = xSeries, y = nbc_Q_centered_model_none_rightData$ML.cline$clineFunc(xSeries), col=rgb(253/255, 231/255, 37/255), lwd=2)
dev.off()

#plot caor
pdf(file = "caor_cline.pdf")
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent")
points(pts_caor, col=rgb(33/255, 145/255, 140/255, 0.5), pch=20, cex = 1.25)
polygon(x = c(caor_fzCoor$x, rev(caor_fzCoor$x)), y = c(caor_fzCoor$yMin, rev(caor_fzCoor$yMax)), border = rgb(33/255, 145/255, 140/255, 0.7), col=rgb(33/255, 145/255, 140/255, 0.5))
lines(x = xSeries, y = caor_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries), col=rgb(33/255, 145/255, 140/255, 0.3), lwd=2)
dev.off()

pdf(file = "caor_cline_std_scale.pdf")
hzar.plot.obsData(caor_Q_centered_model_none_noneData, col = "transparent", xlim=c(-300,400))
points(pts_caor, col=rgb(33/255, 145/255, 140/255, 0.5), pch=20, cex = 1.25)
polygon(x = c(caor_fzCoor$x, rev(caor_fzCoor$x)), y = c(caor_fzCoor$yMin, rev(caor_fzCoor$yMax)), border = rgb(33/255, 145/255, 140/255, 0.7), col=rgb(33/255, 145/255, 140/255, 0.5))
lines(x = xSeries, y = caor_Q_centered_model_none_noneData$ML.cline$clineFunc(xSeries), col=rgb(33/255, 145/255, 140/255, 0.3), lwd=2)
dev.off()


### plot widths and CIs for all
#make df first
widths <- c(man_Q_model_none_noneData$ML.cline$param.free$width, nbc_Q_model_none_rightData$ML.cline$param.free$width,caor_Q_model_none_noneData$ML.cline$param.free$width)
man_CI <- hzar.getLLCutParam(man_Q_model_none_noneData,c("width"))
nbc_CI <- hzar.getLLCutParam(nbc_Q_model_none_rightData,c("width")) 
caor_CI <- hzar.getLLCutParam(caor_Q_model_none_noneData,c("width"))
CIs <- rbind(man_CI, nbc_CI, caor_CI)
width.df <- cbind(widths, CIs)
width.df <- cbind(c("man", "nbc", "caor"), width.df)
colnames(width.df) <- c("transect", "width", "low", "high")
width.df

colorses <- c(rgb(253, 231, 37, maxColorValue = 255), rgb(68, 1, 84, maxColorValue = 255),rgb(33, 145, 140, maxColorValue = 255))

pdf(file = "cline_widths_chart.pdf")
width.df %>%
  mutate(transect = fct_relevel(transect, "nbc", "man", "caor")) %>%
  ggplot(aes(x=transect, y=width, color=transect)) + geom_point(stat="identity") + geom_errorbar(ymin = width.df$low, ymax = width.df$high) + ylim(0,250) + scale_color_manual(values=colorses) + annotate("text", x = as.factor(unique(width.df$transect)), y=c(77.82, 105.61, 231.62), label = c("a","a","b"))
dev.off()

### put elevation behind 

### load in the dfs with the coordinates, values, and distances to isocline
caor_transect <- read.csv("../worldclim/caor_trans_worldclim.csv")
man_transect <- read.csv("../worldclim/man_trans_worldclim.csv")
nbc_transect <- read.csv("../worldclim/nbc_trans_worldclim.csv")

#Please note that the temperature data are in °C * 10. This means that a value of 231 represents 23.1 °C. I only need to change the ones in ˚C
caor_transect$Mean_Diurnal_Range <- caor_transect$Mean_Diurnal_Range/10
man_transect$Mean_Diurnal_Range <- man_transect$Mean_Diurnal_Range/10
man_transect$Mean_Temperature_of_Driest_Quarter <- man_transect$Mean_Temperature_of_Driest_Quarter/10
nbc_transect$Mean_Temperature_of_Driest_Quarter <- nbc_transect$Mean_Temperature_of_Driest_Quarter/10
nbc_transect$Temperature_Seasonality <- nbc_transect$Temperature_Seasonality/10

#plot caor which were 1. Precipitation Seasonality, 2. Annual Precipitation, 3. Mean Diurnal Range

pdf(file="caor_cline_elevation.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(caor_Q_model_none_noneData, pch = 20, fzCol=rgb(33/255, 145/255, 140/255, 0.5), col=rgb(33/255, 145/255, 140/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-307.4367,311.7423) ) 
par(new = TRUE)                             # Add new plot
plot(caor_transect$caor_trans_dist, caor_transect$caor_trans_pts_elev_values, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-307.4367,311.7423))  
axis(side = 4, at = pretty(range(caor_transect$caor_trans_pts_elev_values)))      # Add second axis
mtext("Elevation (m)", side = 4, line = 3)            
dev.off()

pdf(file="caor_cline_precip_driest_month.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(caor_Q_model_none_noneData, fzCol=rgb(33/255, 145/255, 140/255, 0.5), pch = 20, col=rgb(33/255, 145/255, 140/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-307.4367,311.7423) )
par(new = TRUE)                             # Add new plot
plot(caor_transect$caor_trans_dist, caor_transect$Precipitation_of_Driest_Month, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-307.4367,311.7423)) 
axis(side = 4, at = pretty(range(caor_transect$Precipitation_of_Driest_Month)))      # Add second axis
mtext("Precipitation of Driest Month (mm)", side = 4, line = 3)            
dev.off()

pdf(file="caor_cline_precip_warmest_quarter.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(caor_Q_model_none_noneData,fzCol=rgb(33/255, 145/255, 140/255, 0.5), pch = 20, col=rgb(33/255, 145/255, 140/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-307.4367,311.7423) )
par(new = TRUE)                             # Add new plot
plot(caor_transect$caor_trans_dist, caor_transect$Precipitation_of_Warmest_Quarter, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-307.4367,311.7423)) 
axis(side = 4, at = pretty(range(caor_transect$Precipitation_of_Warmest_Quarter)))      # Add second axis
mtext("Precipitation of Warmest quarter (mm)", side = 4, line = 3)            
dev.off()

pdf(file="caor_cline_mean_temp_wettest_qtr.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(caor_Q_model_none_noneData, pch = 20,fzCol=rgb(33/255, 145/255, 140/255, 0.5), col=rgb(33/255, 145/255, 140/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-307.4367,311.7423) )
par(new = TRUE)                             # Add new plot
plot(caor_transect$caor_trans_dist, caor_transect$Mean_Temperature_of_Wettest_Quarter, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-307.4367,311.7423)) 
axis(side = 4, at = pretty(range(caor_transect$Mean_Temperature_of_Wettest_Quarter)))     # Add second axis
mtext("Mean Temperature of Wettest Quarter (\\u00B0 C)", side = 4, line = 3)         
dev.off()

#plot man

### elevation
pdf(file="man_cline_elevation.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(man_Q_model_none_noneData, pch = 20,fzCol=rgb(68/255, 1/255, 84/255, 0.5), col=rgb(68/255, 1/255, 84/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-70.80133,56.29867))
par(new = TRUE)                             # Add new plot
plot(man_transect$man_trans_dist, man_transect$man_trans_pts_elev_values, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-70.80133,56.29867)) 
axis(side = 4, at = pretty(range(man_transect$man_trans_pts_elev_values)))      # Add second axis
mtext("Elevation (m)", side = 4, line = 3)            
dev.off()

pdf(file="man_cline_mean_diurnal_range.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(man_Q_model_none_noneData, pch = 20,fzCol=rgb(68/255, 1/255, 84/255, 0.5), col=rgb(68/255, 1/255, 84/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-70.80133,56.29867))
par(new = TRUE)                             # Add new plot
plot(man_transect$man_trans_dist, man_transect$Mean_Diurnal_Range, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-70.80133,56.29867)) 
axis(side = 4, at = pretty(range(man_transect$Mean_Diurnal_Range)))    # Add second axis
mtext("Mean Diurnal Range (\\u00B0 C)", side = 4, line = 3)            
dev.off()

pdf(file="man_cline_annual_precip.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(man_Q_model_none_noneData, pch = 20,fzCol=rgb(68/255, 1/255, 84/255, 0.5), col=rgb(68/255, 1/255, 84/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-70.80133,56.29867))
par(new = TRUE)                             # Add new plot
plot(man_transect$man_trans_dist, man_transect$Annual_Precipitation, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-70.80133,56.29867)) 
axis(side = 4, at = pretty(range(man_transect$Annual_Precipitation)))      # Add second axis
mtext("Annual Precipitation (mm)", side = 4, line = 3)   
dev.off()

pdf(file="man_cline_precip_seasonality.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(man_Q_model_none_noneData, pch = 20,fzCol=rgb(68/255, 1/255, 84/255, 0.5), col=rgb(68/255, 1/255, 84/255, 0.0 ), xlab = "Distance (km)", ylab = "Q value", xlim = c(-70.80133,56.29867))
par(new = TRUE)                             # Add new plot
plot(man_transect$man_trans_dist, man_transect$Precipitation_Seasonality, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-70.80133,56.29867)) 
axis(side = 4, at = pretty(range(man_transect$Precipitation_Seasonality)))      # Add second axis
mtext("Precipitation Seasonality (\\u00B0 C)", side = 4, line = 3)            
dev.off()

pdf(file="man_cline_annual_mean_temp.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(man_Q_model_none_noneData, pch = 20,fzCol=rgb(68/255, 1/255, 84/255, 0.5), col=rgb(68/255, 1/255, 84/255, 0.0 ), xlab = "Distance (km)", ylab = "Q value", xlim = c(-70.80133,56.29867))
par(new = TRUE)                             # Add new plot
plot(man_transect$man_trans_dist, man_transect$Annual_Mean_Temperature, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-70.80133,56.29867)) 
axis(side = 4, at = pretty(range(man_transect$Annual_Mean_Temperature)))      # Add second axis
mtext("Annual Mean Temperature (\\u00B0 C)", side = 4, line = 3)            
dev.off()

#plot nbc

### elevation
pdf(file="nbc_cline_elevation.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(nbc_Q_model_none_rightData, pch = 20,fzCol=rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.0), xlab = "Distance (km)", ylab = "Q value", xlim = c(-41.83419, 489.36941) )
par(new = TRUE)                             # Add new plot
plot(nbc_transect$nbc_trans_dist, nbc_transect$nbc_trans_pts_elev_values, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-41.83419, 489.36941)) 
axis(side = 4, at = pretty(range(nbc_transect$nbc_trans_pts_elev_values)))      # Add second axis
mtext("Elevation (m)", side = 4, line = 3)            
dev.off()

pdf(file="nbc_cline_temp_seasonality.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(nbc_Q_model_none_rightData, pch = 20,fzCol=rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.0), xlab = "Distance (km)", ylab = "Q value" , xlim = c(-41.83419, 489.36941))
par(new = TRUE)                             # Add new plot
plot(nbc_transect$nbc_trans_dist, nbc_transect$Temperature_Seasonality, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-41.83419, 489.36941)) 
axis(side = 4, at = pretty(range(nbc_transect$Temperature_Seasonality)))      # Add second axis
mtext("Temperature Seasonality", side = 4, line = 3)            
dev.off()

pdf(file="nbc_cline_annual_precip.pdf")
par(mar = c(5, 4, 4, 4) + 0.3) 
hzar.plot.fzCline(nbc_Q_model_none_rightData, pch = 20,fzCol=rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.0), xlab = "Distance (km)", ylab = "Q value" , xlim = c(-41.83419, 489.36941))
par(new = TRUE)                             # Add new plot
plot(nbc_transect$nbc_trans_dist, nbc_transect$Annual_Precipitation, type = "l", lty = 3, lwd=6, axes = FALSE, xlab = "", ylab = "", xlim = c(-41.83419, 489.36941)) 
axis(side = 4, at = pretty(range(nbc_transect$Annual_Precipitation)))      # Add second axis
mtext("Annual Precipitation (mm)", side = 4, line = 3)            
dev.off()



