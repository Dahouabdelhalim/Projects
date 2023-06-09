############################
# Code for Parasites have variable effects on the outcomes of host species interactions
############################

#############
# Ensure the following data files are in the working directory to enable the code to run smoothly
#############

# es_db.csv
# host_tax_list.csv
# host_tree.tre
# host_tree_no_out.tre
# outlier.list.no.phy.with.virus.csv
# outlier.list.phy
# outlier.list.no.phy.csv
# parasite_tax_list.csv
# parasite_tree.tre
# parasite_tree_no_out.tre


setwd("C:/Users/azhas/Dropbox/Meta_Analysis/AmNat Submission")

# clear memory
rm(list=ls())

# load required packages
library(devtools)
library(orchaRd)
require(metafor)
require(Gmisc)
require(forestplot)
require(ggplot2)
require(ggpubr)
library(ggfortify)
library(gridExtra)
library(cowplot)
require(grid)
require(gridExtra)
require(EnvStats)
require(rcompanion)
require(RColorBrewer)
require(plyr)
require(dplyr)
require(tidyr)
library(rotl)
library(ape)
library(phytools)
library(reshape2)

# Load two color blind palettes
cbbPalette <- c( "#D55E00","#56B4E9", "#009E73","#CC79A7", "#999999")
cbbPalette.par <- c( "#D55E00","#D55E00","#56B4E9","#56B4E9", "#009E73", "#009E73")

# load database of studies
es_calculator <- read.csv("es_db.csv")

# check for heterscedascity of population variances before calculating effect sizes
var.test(es_calculator$control_mean_transformed,
         es_calculator$treatment_mean_transformed)
# variances are heteroscedastic, use SMDH measure in escalc()

# calculate effect sizes
dat <- escalc(measure="SMDH",
              n1i=n_treatment,
              n2i=n_control,
              m1i=treatment_mean_transformed,
              m2i=control_mean_transformed,
              sd1i=sd_treatment,
              sd2i=sd_control,
              data=es_calculator,
              append=FALSE)

es_calculator$yi<-dat$yi #add effect size back into dataframe
es_calculator$vi<-dat$vi #add effect size variance back into dataframe

# make a dataframe without viruses for phylogeny
es_calculator_no_virus<-subset(es_calculator,es_calculator$Parasite_Narrow!="Virus"&
                                 es_calculator$Parasite_Narrow!="Unclassified")
# add host species names to dataframe
list<-read.csv("host_tax_list.csv") #load in list of host species for phylogeny
es_calculator_no_virus$names<-na.omit(list$name) #add host names to dataframe
# add parasite species names to dataframe
p.list<-read.csv("parasite_tax_list.csv") #list of parasite species for phylogeny
es_calculator_no_virus$p.names<-p.list$name #add parasite names to dataframe

# update the five host taxa that need it such that their names match the names on the tips of the tree
es_calculator_no_virus$names<-gsub("blattella germanica","Blattella germanica",es_calculator_no_virus$names)
es_calculator_no_virus$names<-gsub("Bombus terrestris","Bombus terrestris -species in domain Eukaryota",es_calculator_no_virus$names)
es_calculator_no_virus$names<-gsub("Apis mellifera","Apis mellifera -species in domain Eukaryota",es_calculator_no_virus$names)
es_calculator_no_virus$names<-gsub("Oncorhynchus mykiss","Oncorhynchus mykiss -species in domain Eukaryota",es_calculator_no_virus$names)
es_calculator_no_virus$names<-gsub("Rhynchophorus ferrugineus","Rhynchophorus ferrugineus -species in domain Eukaryota",es_calculator_no_virus$names)
es_calculator_no_virus$names<-gsub(" ","_",es_calculator_no_virus$names)
es_calculator_no_virus$names<-as.factor(es_calculator_no_virus$names)

# update the parasite taxa that need it such that their names match the names on the tips of the tree
es_calculator_no_virus$p.names<-gsub("Trypanosoma","Trypanosoma -genus in infrakingdom Excavata",es_calculator_no_virus$p.names)
es_calculator_no_virus$p.names<-gsub("Pseudacteon","Pseudacteon -genus in Ecdysozoa",es_calculator_no_virus$p.names)
es_calculator_no_virus$p.names<-gsub(" ","_",es_calculator_no_virus$p.names)
es_calculator_no_virus$p.names<-as.factor(es_calculator_no_virus$p.names)

# create additional factors needed below for host x parasite correlation matrix interaction
es_calculator_no_virus$hostxparasite <- paste(es_calculator_no_virus$names, es_calculator_no_virus$p.names, sep="XX")
es_calculator_no_virus$hostxparasite <- factor(es_calculator_no_virus$hostxparasite)

##########################################################
# read in plant and fungal phylogenetic trees
##########################################################

tree.h <- read.tree(file="host_tree.tre")
tree.h<-compute.brlen(tree.h) #calculate branch lengths
R.H <- vcv(phy=tree.h, corr=TRUE) # calculate correlation matrix
round(R.H[1:10,1:10],3) #view part of the R.H object

tree.p <- read.tree(file="parasite_tree.tre")
tree.p<-compute.brlen(tree.p) #calculate branch lengths
R.P <- vcv(phy=tree.p, corr=TRUE) # calculate correlation matrix
round(R.P[1:10,1:10],3) # view part of the R.P object

# check that all of the hosts/parasites in the data are actually in the trees
is.element(unique(es_calculator_no_virus$names), tree.h$tip) #check hosts
is.element(unique(es_calculator_no_virus$p.names), tree.p$tip) #check parasites

# Calculate dataframe for tensor products of phylogenetic vcv matrices, for observed species in dataset
# for use with host x parasite interaction effects
hostparphy <- data.frame(levels(es_calculator_no_virus$hostxparasite))
colnames(hostparphy) <- "HPcombo"
hostparphy_split<-colsplit(hostparphy$HPcombo, "XX", c("host","parasite"))
RO.HP <- matrix(NA, nrow=length(hostparphy$HPcombo), ncol=length(hostparphy$HPcombo))
rownames(RO.HP) <- levels(es_calculator_no_virus$hostxparasite) #rownames to interaction matrix
colnames(RO.HP) <- levels(es_calculator_no_virus$hostxparasite)

# calculate the tensor products for each level of host x parasite interaction
for (i in 1:length(hostparphy_split$host)) {
  for (j in 1:length(hostparphy_split$parasite)) {
    if (i <= j)
      next 
         RO.HP[i,j] <- R.H[hostparphy_split$host[i], hostparphy_split$host[j]] * R.P[hostparphy_split$parasite[i], hostparphy_split$parasite[j]]
  }
}
diag(RO.HP) <- 1
RO.HP[upper.tri(RO.HP)] <- t(RO.HP)[upper.tri(RO.HP)]

# calculation of average sample size in database
n.all<-(es_calculator_no_virus$n_control+es_calculator_no_virus$n_treatment)/2 #for dataframe that does not include viruses
mean(n.all)
# average sample size is ~ 40, REML estimator is appropriate to use

n.all.virus<-(es_calculator$n_control+es_calculator$n_treatment)/2 #for dataframe that includes viruses
mean(n.all.virus)
# average sample size is ~ 40, REML estimator is appropriate to use

###############################################
## OVERALL INTERCEPT ONLY MODEL RUNS
###############################################

# fit basic multi-level models

# with phylogeny
res.phy<- rma.mv(yi, vi, data=es_calculator_no_virus,
             random= list(~ 1 | Study_ID/Effect_Size_ID,
                          ~ 1 | names,
                          ~ 1 | p.names,
                          ~ 1 | hostxparasite),
             R=list(names=R.H,p.names=R.P,hostxparasite=RO.HP),method = "REML")
# without phylogeny
res.no.phy<- rma.mv(yi, vi, data=es_calculator_no_virus,
                 random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
# without phylogeny and with viruses
res<- rma.mv(yi, vi, data=es_calculator,
              random= list(~ 1 | Study_ID/Effect_Size_ID),
             method = "REML")
res.phy #overall effect of parasites with phylogeny -0.61 [-0.96,-0.26]
res.no.phy #overall effect of parasites without phylogeny -0.63 [-0.85,-0.41]
res #overall effect of parasites without phylogeny and with viruses -0.63 [-0.84,-0.41]

######################
# CREATE ORCHARD PLOTS FOR FIG. 1
######################

# with phylogeny
phy.orchard<-orchard_plot(res.phy, mod = "Int",
                          xlab = "Standardized mean difference",
                          transfm = "none",cb=TRUE,k=FALSE,alpha=.2)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15))
# without phylogeny
no.phy.orchard<-orchard_plot(res.no.phy, mod = "Int",xlab = "Standardized mean difference",
                             transfm = "none",cb=TRUE,k=FALSE,alpha=.2)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15))
# without phylogeny and with viruses
viruses.orchard<-orchard_plot(res, mod = "Int",xlab = "Standardized mean difference",
                                 transfm = "none",cb=TRUE,k=FALSE,alpha=.2)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15)) 
fig.1a <-
  ggdraw() +
  draw_plot(phy.orchard+xlim(-5,5)+labs(tag = "(a)")+
              annotate(geom="text", x=-5, y=1.5,label="Controlling for phylogeny",hjust=0,
                       size=4.5,fontface="italic")+
              theme(legend.position = "none",
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"))) +
  draw_plot(phy.orchard+theme(axis.text.y = element_blank(),
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              legend.position = "none",
                              panel.border = element_blank(),
                              axis.line = element_line(color="black")),
            x = .7, y = .58,
            width = .3, height = .35)
fig.1b <-
  ggdraw() +
  draw_plot(no.phy.orchard+xlim(-5,5)+labs(tag = "(b)")+
              annotate(geom="text", x=-5, y=1.5,hjust=0,
                       label="Without controlling for phylogeny",size=4.5,fontface="italic")+
              theme(legend.position = "none",
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"))) +
  draw_plot(no.phy.orchard+theme(axis.text.y = element_blank(),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 legend.position = "none",
                                 panel.border = element_blank(),
                                 axis.line = element_line(color="black")),
            x = .7, y = .58,
            width = .3, height = .35)
fig.1c <-
  ggdraw() +
  draw_plot(viruses.orchard+xlim(-5,5)+labs(tag = "(c)")+
              annotate(geom="text", x=-5, y=1.5,hjust=0,
                       label="Without controlling for phylogeny, with viruses",size=4.5,fontface="italic")+
              theme(axis.text.y = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"),
                    legend.position = c(.94,0))) +
  draw_plot(viruses.orchard+theme(axis.text.y = element_blank(),
                                     axis.title.x=element_blank(),
                                     legend.position = "none",
                                  panel.border = element_blank(),
                                  axis.line = element_line(color="black")),
            x = .7, y = .58,
            width = .3, height = .35)
fig.1<-plot_grid(fig.1a,
                 fig.1b,
                 fig.1c,align = "hv",hjust=-1,
                 ncol = 1,rel_heights = c(1,1,1))
fig.1
# ggsave(fig.1, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure 1.jpeg",
#        height=9, width=9, units="in", dpi=1500)

####################
# CREATE PLOTS FOR FIG. S3
####################

####################
# STUDY TYPE x HOST ORGANISM
####################
# set up a dataframe that is grouped by host_broad and study_type for all data (including viruses)
dfsum <- group_by(es_calculator,Host_Broad,Study_Type)
# get summary counts for the each level of study type crossed with host broad
num_stud<-summarize(dfsum,study_label=n_distinct(study_label))
num_stud$Host_Broad<-as.factor(num_stud$Host_Broad)

# graph data
host.study.plot = ggplot(data=num_stud,
                         aes(x = Host_Broad,y = study_label,
                             color=Study_Type,shape=Study_Type,
                             fill=Study_Type,group=Study_Type))+
  geom_col(position = position_dodge2(preserve = "single", width = 1))+
  xlab('Host Organism')+ ylab("Number of Studies")+
  theme_classic()+scale_color_manual(values=cbbPalette,
                                     name="Study Type")+
  scale_fill_manual(values=cbbPalette,
                    name="Study Type")+
  scale_x_discrete(breaks=c("Fungi","Invertebrate","Plant","Vertebrate"))+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90))+
  theme(text = element_text(size=20),axis.text.x = element_text(angle = -45,hjust = .1))+labs(tag="(a)")
host.study.plot
####################
# STUDY TYPE x PARASITE ORGANISM
####################

# set up a dataframe that is grouped by host_broad and study_type for all data (including viruses)
dfsum.p <- group_by(es_calculator,Parasite_Narrow,Study_Type)
# get summary counts for the each level of study type crossed with host broad
num_stud.p<-summarize(dfsum.p,study_label=n_distinct(study_label))
num_stud.p

# graph data
para.study.plot = ggplot(data=num_stud.p,
                         aes(x = Parasite_Narrow,y = study_label,
                             color=Study_Type,shape=Study_Type,
                             fill=Study_Type,group=Study_Type))+
  geom_col(position = position_dodge2(preserve = "single", width = 1))+
  xlab('Parasite Organism')+ ylab("Number of Studies")+
  theme_classic()+scale_color_manual(values=cbbPalette,
                                     name="Study Type")+
  scale_fill_manual(values=cbbPalette,
                    name="Study Type")+
  scale_x_discrete(breaks=c("Acanthocephalan","Apicomplexan","Arthropod","Bacteria",
                            "Fungus","Nematode","Plant","Platyhelminth",
                            "Protist","Protozoan","Unclassified","Virus"))+
  scale_y_continuous(limits = c(0,80),
                     breaks=c(0,10,20,30,40,50,60,70,80))+
  theme(text = element_text(size=20),
        # aspect.ratio = 1,
        axis.text.x = element_text(angle = -45,hjust = .1))+labs(tag="(b)")

# Make Figure S3
legend<-get_legend(para.study.plot+theme(legend.position = "right", #extract legend
                                         text = element_text(size=20),
                                         plot.margin=grid::unit(c(-1,0,0,0), "mm")))
fig.s3<-plot_grid(host.study.plot+theme(legend.position = "none",
                                        text = element_text(size=20),
                                        plot.margin=grid::unit(c(0,0,0,0), "mm")),
                  para.study.plot+theme(legend.position = "none",
                                        text = element_text(size=20),
                                        plot.margin=grid::unit(c(0,0,0,0), "mm")),
                  align="hv",
                  hjust=-1,
                  ncol=1,
                  rel_heights=c(1,1))
fig.s3<-plot_grid(fig.s3,legend,align='hv',hjust=-1,ncol=2,rel_widths = c(1,.2))
fig.s3
# ggsave(fig.s3, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S3.jpeg",
#        height=9, width=11, units="in", dpi=1500)

##############################################
#check heterogeneity and sampling variance of the intercept only models
##############################################

###
# WITH PHYLOGENY
###

# first check total I2 (total variance due to heterogeneity)
W <- diag(1/es_calculator_no_virus$vi)
X <- model.matrix(res.phy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res.phy$sigma2) / (sum(res.phy$sigma2) + (res.phy$k-res.phy$p)/sum(diag(P)))
# 97.39% of the variance is due to heterogeneity

# now check the between and within-cluster heterogeneity
100 * res.phy$sigma2 / (sum(res.phy$sigma2) + (res.phy$k-res.phy$p)/sum(diag(P)))
# 62.60% is due to between-study, 32.31% due to within-study
# 2.48% is due to host phylogeny, 3.86e-7% is due to parasite phylogeny
# 2.13e-9 is due to the interaction of host and parasite phylogeny

###
# WITHOUT PHYLOGENY
###

# first check total I2 (total variance due to heterogeneity)
W <- diag(1/es_calculator_no_virus$vi)
X <- model.matrix(res.no.phy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res.no.phy$sigma2) / (sum(res.no.phy$sigma2) + (res.no.phy$k-res.no.phy$p)/sum(diag(P)))
# 97.37% of the variance is due to heterogeneity

# now check the between and within-cluster heterogeneity
100 * res.no.phy$sigma2 / (sum(res.no.phy$sigma2) + (res.no.phy$k-res.no.phy$p)/sum(diag(P)))
# 64.71% is due to between-study, 32.65% due to within-study

###
# WITHOUT PHYLOGENY AND WITH VIRUSES
###

# first check total I2 (total variance due to heterogeneity)
W <- diag(1/es_calculator$vi)
X <- model.matrix(res)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))
# 97.31% of the variance is due to heterogeneity

# now check the between and within-cluster heterogeneity
100 * res$sigma2 / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))
# 63.77% is due to between-study, 33.54% due to within-study

##############################################
## CHECK FOR PUBLICATION BIAS AND OUTLIERS
##############################################

####
# WITH PHYLOGENY
####

# publication bias
funnel(res.phy)
# the plot is mostly symmetrical, simply scattered more horizontally

####
# WITHOUT PHYLOGENY
####

# publication bias
funnel(res.no.phy)
# the plot is mostly symmetrical, simply scattered more horizontally

####
# WITHOUT PHYLOGENY AND WITH VIRUSES
####

# publication bias
funnel(res)
# the plot is mostly symmetrical, simply scattered more horizontally

#######
# CREATE FUNNEL PLOTS FOR FIG. S2
#######
# jpeg("plots/Figure S2.jpeg",width=4,height=8,units="in",res=1500)
par(mfcol=c(3,1),oma=c(3,2,0,0),mai=c(.15,.7,.01,.7))
xlmts <- c(-22,22)
ylmts <- c(0,7)
funnel(res.phy,xlab = "",ylab="",xlim=xlmts,ylim=ylmts,xaxt='n')
funnel(res.no.phy,ylab="",
       xlim=xlmts,ylim=ylmts,xlab = "",xaxt='n')
funnel(res,xlab="",ylab="",xlim=xlmts,ylim=ylmts,)
mtext("Standardized mean difference",side=1,line=2,outer=TRUE,cex=1.3)
mtext("Standard error",side=2,line=0,outer=TRUE,cex=1.3,las=0)
# dev.off()
par(mfrow=c(1,1))

####
# TRIM AND FILL ANALYSIS
### 

# build a less complex model
res.less<- rma(yi, vi, data=es_calculator,measure = "SMDH")
funnel(res.less)
# carry out trim-and-fill analysis
taf <- trimfill(res.less)
# draw funnel plot with missing studies filled in
funnel(taf, legend=TRUE)
# nothing appears to be missing

#######################
#CHECK FOR OUTLIERS
#######################

# WARNING: this code to calculate cook's distance may take hours to run

####
# WITH PHYLOGENY
####

# # Cook's distances for each effect size
# x <- cooks.distance(res.phy, cluster=es_calculator_no_virus$Accession_number,progbar = TRUE)
# mean(x)
# mean distance is 0.0003077945
# rule of thumb is any distance > three times the mean distance
# is a potential outlier and influential on the regression results
# mean(x)*3
# all distances above 0.0009233835 may be outliers
# create list of influential effect sizes
# x_no_out<-subset(x,x>(mean(x)*3))
# names(x_no_out)
# x_no_out<-as.numeric(names(x_no_out)) #create vector of influential outliers
# # write.csv(x_no_out,"outlier.list.csv")

x_no_out<-read.csv("outlier.list.phy.csv") #load list of influential outliers for model with phylogeny
x_no_out<-x_no_out$x
# remove effect sizes that are outliers
es_calculator_no_virus_no_out<-es_calculator_no_virus %>%
  filter(!Accession_number %in% x_no_out)

# #subset the host phylogeny and correlation matrix to exclude influential outliers
# list.no.out<-cbind(list,es_calculator_no_virus$Accession_number)
# list.no.out<-list.no.out %>%
#   filter(!`es_calculator_no_virus$Accession_number` %in% x_no_out)
# resolved_names.no.out<-tnrs_match_names(names=(list.no.out$name))
# my_tree.no.out <- tol_induced_subtree(ott_ids = resolved_names.no.out$ott_id, label_format ="name")
# 
# plot(my_tree.no.out,no.margin = TRUE)
# write.tree(my_tree.no.out,"host_tree_no_out.tre")

# #subset the parasite phylogeny and correlation matrix
# plist.no.out<-cbind(p.list,es_calculator_no_virus$Accession_number)
# plist.no.out<-plist.no.out %>%
#   filter(!`es_calculator_no_virus$Accession_number` %in% x_no_out)
# plist.no.out<-na.omit(plist.no.out)
# p.resolved_names.no.out<-tnrs_match_names(names=(plist.no.out$name))
# p.resolved_names.no.out<-update(p.resolved_names.no.out,taxon_name="trypanosoma",new_ott_id=779558) #update trypanosoma info
# p_tree.no.out <- tol_induced_subtree(ott_ids = p.resolved_names.no.out$ott_id, label_format ="name")
# 
# plot(p_tree.no.out,no.margin = TRUE)
# write.tree(p_tree.no.out,"parasite_tree_no_out.tre")

##########################################################
# read in phylogenetic trees without outliers
tree.h.no.out <- read.tree(file="host_tree_no_out.tre")
tree.h.no.out<-compute.brlen(tree.h.no.out)
R.H.no.out <- vcv(phy=tree.h.no.out, corr=TRUE) # corr=TRUE gives a correlation matrix
round(R.H.no.out[1:10,1:10],3) #view part of the R.H object

tree.p.no.out <- read.tree(file="parasite_tree_no_out.tre")
tree.p.no.out<-compute.brlen(tree.p.no.out)
R.P.no.out <- vcv(phy=tree.p.no.out, corr=TRUE)
round(R.P.no.out[1:10,1:10],3) # view part of the R.P object

# check that all of the hosts/parasites in the data are actually in the trees
is.element(unique(es_calculator_no_virus_no_out$names), tree.h.no.out$tip)
is.element(unique(es_calculator_no_virus_no_out$p.names), tree.p.no.out$tip)

# create additional factors needed below for host x parasite interaction
es_calculator_no_virus_no_out$hostxparasite <- paste(es_calculator_no_virus_no_out$names, es_calculator_no_virus_no_out$p.names, sep="XX")
es_calculator_no_virus_no_out$hostxparasite <- factor(es_calculator_no_virus_no_out$hostxparasite)

# Calculate tensor product of phylogenetic vcv matrices, for observed species in dataset
# for use with host x parasite interaction effects without influential outliers
hostparphy.no.out <- data.frame(levels(es_calculator_no_virus_no_out$hostxparasite))
colnames(hostparphy.no.out) <- "HPcombo"
hostparphy_split<-colsplit(hostparphy.no.out$HPcombo, "XX", c("host","parasite"))
RO.HP.no.out <- matrix(NA, nrow=length(hostparphy.no.out$HPcombo), ncol=length(hostparphy.no.out$HPcombo))
rownames(RO.HP.no.out) <- levels(es_calculator_no_virus_no_out$hostxparasite) #rownames to interaction matrix
colnames(RO.HP.no.out) <- levels(es_calculator_no_virus_no_out$hostxparasite)

#calculate the tensor products for each level of host x parasite interaction without influential outliers
for (i in 1:length(hostparphy_split$host)) {
  for (j in 1:length(hostparphy_split$parasite)) {
    if (i <= j)
      next 
    RO.HP.no.out[i,j] <- R.H.no.out[hostparphy_split$host[i], hostparphy_split$host[j]] * R.P.no.out[hostparphy_split$parasite[i], hostparphy_split$parasite[j]]
  }
}
diag(RO.HP.no.out) <- 1
RO.HP.no.out[upper.tri(RO.HP.no.out)] <- t(RO.HP.no.out)[upper.tri(RO.HP.no.out)]

# run overall model with phylogeny and without influential outliers
res.phy.no.out<- rma.mv(yi, vi, data=es_calculator_no_virus_no_out,
                 random= list(~ 1 | Study_ID/Effect_Size_ID,
                              ~ 1 | names,
                              ~ 1 | p.names,
                              ~ 1 | hostxparasite),
                 R=list(names=R.H.no.out,p.names=R.P.no.out,hostxparasite=RO.HP.no.out),method = "REML")
res.phy.no.out #overall theta is -0.54 [-1.11, 0.03]

####
# WITHOUT PHYLOGENY
####

# ## Cook's distances for each effect size
# x <- cooks.distance(res.no.phy, cluster=es_calculator_no_virus$Accession_number,progbar = TRUE)
# mean(x)
# #mean distance is 0.0006010538
# #rule of thumb is any distance > three times the mean distance
# #is a potential outlier and influential on the regression results
# mean(x)*3
# #all distances above 0.0009233835 may be outliers
# #list of influential effect sizes
# x_no_out<-subset(x,x>(mean(x)*3))
# names(x_no_out)
# x_no_out<-as.numeric(names(x_no_out)) #create vector of influential outliers
# write.csv(x_no_out,"outlier.list.no.phy.csv")

x_no_out<-read.csv("outlier.list.no.phy.csv") #load list of influential outliers for model without phylogeny
x_no_out<-x_no_out$x
#remove effect sizes that are outliers
es_calculator_no_phy_no_out<-es_calculator_no_virus %>%
  filter(!Accession_number %in% x_no_out)

# run overall model without phylogeny and without influential outliers
res.no.phy.no.out<- rma.mv(yi, vi, data=es_calculator_no_phy_no_out,
                           random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.no.phy.no.out #overall theta is -0.48 [-0.68, -0.29]

####
# WITHOUT PHYLOGENY AND WITH VIRUSES
####

# ## Cook's distances for each effect size
# x <- cooks.distance(res, cluster=es_calculator$Accession_number,progbar = TRUE)
# mean(x)
# #mean distance is 0.0005420208
# #rule of thumb is any distance > three times the mean distance
# #is a potential outlier and influential on the regression results
# mean(x)*3
# #all distances above 0.001626062 may be outliers
# #list of influential effect sizes
# x_no_out<-subset(x,x>(mean(x)*3))
# # names(x_no_out)
# # x_no_out<-as.numeric(names(x_no_out)) #create vector of influential outliers
# # write.csv(x_no_out,"outlier.list_no_phy.csv")

x_no_out<-read.csv("outlier.list.no.phy.with.virus.csv") #load list of influential outliers for model without phylogeny
x_no_out<-x_no_out$x                                     #and with viruses
# remove effect sizes that are outliers
es_calculator_no_out<-es_calculator %>%
  filter(!Accession_number %in% x_no_out)

# run overall model without phylogeny, with viruses and without influential outliers
res.no.out<- rma.mv(yi, vi, data=es_calculator_no_out,
                    random= list(~ 1 | Study_ID/Effect_Size_ID),
                    method = "REML")
res.no.out #overall theta is -0.49 [-0.68, -0.30]

######################
# CREATE ORCHARD PLOTS FOR FIGURE S4
######################

# with phylogeny
phy.no.out.orchard<-orchard_plot(res.phy.no.out, mod = "Int",
                          xlab = "Standardized mean difference",
                          transfm = "none",cb=TRUE,k=FALSE,alpha=.2)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15))

# without phylogeny
no.phy.no.out.orchard<-orchard_plot(res.no.phy.no.out, mod = "Int",xlab = "Standardized mean difference",
                             transfm = "none",cb=TRUE,k=FALSE,alpha=.2)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.position = "none")

# without phylogeny and with viruses
viruses.no.out.orchard<-orchard_plot(res.no.out, mod = "Int",xlab = "Standardized mean difference",
                                 transfm = "none",cb=TRUE,k=FALSE,alpha=.2)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15)) 

fig.s4a <-
  ggdraw() +
  draw_plot(phy.no.out.orchard+labs(tag = "(a)")+
              annotate(geom="text", x=-20, y=1.5,hjust=0,
                       label="Controlling for phylogeny",size=4.5,fontface="italic")+
              theme(legend.position = "none",
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"))+expand_limits(y=c(1,1)))
fig.s4b <-
  ggdraw() +
  draw_plot(no.phy.no.out.orchard+labs(tag = "(b)")+
              annotate(geom="text", x=-20, y=1.5,hjust=0,
                       label="Without controlling for phylogeny",size=4.5,fontface="italic")+
              theme(legend.position = "none",
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"))+expand_limits(y=c(1,1)))
fig.s4c <-
  ggdraw() +
  draw_plot(viruses.no.out.orchard+labs(tag = "(c)")+
              annotate(geom="text", x=-20, y=1.5,hjust=0,
                       label="Without controlling for phylogeny, with viruses",size=4.5,fontface="italic")+
              theme(axis.text.y = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"),
                    legend.position = c(1,0))+expand_limits(y=c(1,1)))

fig.s4<-plot_grid(fig.s4a,
                 fig.s4b,
                 fig.s4c,align = "hv",hjust=-1,
                 ncol = 1,rel_heights = c(1,1,1))
fig.s4

# ggsave(fig.s4, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S4.jpeg",
#        height=9, width=9, units="in", dpi=1500)


########
# SUBSET DATAFRAMES TO INCLUDE ONLY PREDATION, COMPEITION, AND REPRODUCTION STUDIES FOR MODEL COMPARISONS
########

# subset database without viruses
es_calculator_no_virus_subset<-subset(es_calculator_no_virus,es_calculator_no_virus$Study_Type=="Predation"|
                                        es_calculator_no_virus$Study_Type=="Competition"|
                                        es_calculator_no_virus$Study_Type=="Reproduction")

# subset database without viruses for models with phylogeny without outliers
es_calculator_no_virus_no_out_subset<-subset(es_calculator_no_virus_no_out,
                                             es_calculator_no_virus_no_out$Study_Type=="Predation"|
                                               es_calculator_no_virus_no_out$Study_Type=="Competition"|
                                               es_calculator_no_virus_no_out$Study_Type=="Reproduction")

# subset database without viruses for models without phylogeny and without outliers
es_calculator_no_phy_no_out_subset<-subset(es_calculator_no_phy_no_out,
                                           es_calculator_no_phy_no_out$Study_Type=="Predation"|
                                             es_calculator_no_phy_no_out$Study_Type=="Competition"|
                                             es_calculator_no_phy_no_out$Study_Type=="Reproduction")

# subset database with viruses
es_calculator_subset<-subset(es_calculator,es_calculator$Study_Type=="Predation"|
                              es_calculator$Study_Type=="Competition"|
                              es_calculator$Study_Type=="Reproduction")

# subset database with viruses and without outliers
es_calculator_subset_no_out<-subset(es_calculator_no_out,
                                     es_calculator_no_out$Study_Type=="Predation"|
                                       es_calculator_no_out$Study_Type=="Competition"|
                                       es_calculator_no_out$Study_Type=="Reproduction")
####################
# EFFECTS OF PARASITES BY STUDY TYPE
####################

# subset phylogeny correlation matrices
subset.names<-as.vector(unique(es_calculator_no_virus_subset$names)) #subset host cor matrix
R.H.subset<-R.H[subset.names,subset.names]
subset.p.names<-as.vector(unique(es_calculator_no_virus_subset$p.names)) #subset par cor matrix
R.P.subset<-R.P[subset.p.names,subset.p.names]
subset.hostxparasite<-as.vector(unique(es_calculator_no_virus_subset$hostxparasite)) #subset host x par cor matrix
RO.HP.subset<-RO.HP[subset.hostxparasite,subset.hostxparasite]

# model with phylogeny
res.phy.study<- rma.mv(yi, vi,mods=~Study_Type-1,data=es_calculator_no_virus_subset,
                        random= list(~ 1 | Study_ID/Effect_Size_ID,
                                     ~ 1 | names,
                                     ~ 1 | p.names,
                                     ~ 1 | hostxparasite),
                        R=list(names=R.H.subset,
                               p.names=R.P.subset,hostxparasite=RO.HP.subset),method = "REML",
                       control=list(optimizer="optim", optmethod="BFGS"))
res.phy.study

# model without phylogeny
res.no.phy.study<- rma.mv(yi, vi,mods=~Study_Type-1,data=es_calculator_no_virus_subset,
                       random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.no.phy.study

# model without phylogeny and with viruses
res.study<- rma.mv(yi, vi,mods=~Study_Type-1,data=es_calculator_subset,
                   random= list(~ 1 | Study_ID/Effect_Size_ID))
res.study

####################
# CREATE ORCHARD PLOTS FOR FIG. 2
####################

fig.2a <- ggdraw()+draw_plot(orchard_plot(res.phy.study, mod = "Study_Type", xlab = "Standardized mean difference",
                       angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color="black")) +labs(tag = "(a)")+
    annotate(geom="text", x=-5, y=3.6,label="Controlling for phylogeny",hjust=0,
             size=4.5,fontface="italic")+
  scale_color_manual(values=cbbPalette[1:3])+
  scale_fill_manual(values=cbbPalette[1:3])+expand_limits(y=c(0,4)))

fig.2b <- ggdraw()+draw_plot(orchard_plot(res.no.phy.study, mod = "Study_Type", xlab = "Standardized mean difference",
                       angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color="black")) +labs(tag = "(b)")+
    annotate(geom="text", x=-5, y=3.6,label="Without controlling for phylogeny",hjust=0,
             size=4.5,fontface="italic")+
  scale_color_manual(values=cbbPalette[1:3])+
  scale_fill_manual(values=cbbPalette[1:3])+expand_limits(y=c(0,4)))

fig.2c <- ggdraw()+draw_plot(orchard_plot(res.study, mod = "Study_Type", xlab = "Standardized mean difference",
                       angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = c(1,0)) +labs(tag = "(c)")+
    annotate(geom="text", x=-5, y=3.6,label="Without controlling for phylogeny, with viruses",hjust=0,
             size=4.5,fontface="italic")+
  scale_color_manual(values=cbbPalette[1:3])+
  scale_fill_manual(values=cbbPalette[1:3])+expand_limits(y=c(0,4)))

fig.2<-plot_grid(fig.2a,fig.2b,fig.2c,
                 align = "hv",hjust=-1,
                 ncol = 1,rel_heights = c(1,1,1))
fig.2

# ggsave(fig.2, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure 2.jpeg",
#        height=9, width=10, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY STUDY TYPE WITHOUT OUTLIERS
####################

# subset phylogenetic correlation matrices
subset.names.no.out<-as.vector(unique(es_calculator_no_virus_no_out_subset$names)) #subset host cor matrix
R.H.subset.no.out<-R.H[subset.names.no.out,subset.names.no.out]
subset.p.names.no.out<-as.vector(unique(es_calculator_no_virus_no_out_subset$p.names)) #subset par cor matrix
R.P.subset.no.out<-R.P[subset.p.names.no.out,subset.p.names.no.out]
subset.hostxparasite.no.out<-as.vector(unique(es_calculator_no_virus_no_out_subset$hostxparasite)) #subset host x par cor matrix
RO.HP.subset.no.out<-RO.HP[subset.hostxparasite.no.out,subset.hostxparasite.no.out]

# model with phylogeny
res.phy.study.no.out<- rma.mv(yi, vi,mods=~Study_Type-1,data=es_calculator_no_virus_no_out_subset,
                              random= list(~ 1 | Study_ID/Effect_Size_ID,
                                           ~ 1 | names,
                                           ~ 1 | p.names,
                                           ~ 1 | hostxparasite),
                              R=list(names=R.H.subset.no.out,
                                     p.names=R.P.subset.no.out,hostxparasite=RO.HP.subset.no.out),method = "REML",
                              control=list(optimizer="optim", optmethod="BFGS"))
res.phy.study.no.out

# model without phylogeny
res.no.phy.study.no.out<- rma.mv(yi, vi,mods=~Study_Type-1,data=es_calculator_no_phy_no_out_subset,
                              random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.no.phy.study.no.out

# model without phylogeny and with viruses
res.study.no.out<- rma.mv(yi, vi,mods=~Study_Type-1,data=es_calculator_subset_no_out,
                          random= list(~ 1 | Study_ID/Effect_Size_ID))
res.study.no.out

####################
# CREATE ORCHARD PLOTS FOR FIG. S5
####################

fig.s5a <-ggdraw()+
  draw_plot(orchard_plot(res.phy.study.no.out, 
                         mod = "Study_Type",
                         xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black")) +labs(tag = "(a)")+
                               annotate(geom="text", x=-5, y=3.6,label="Controlling for phylogeny",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette[1:3])+
                               scale_fill_manual(values=cbbPalette[1:3])+expand_limits(y=c(0,4)))

fig.s5b <- ggdraw()+
  draw_plot(orchard_plot(res.no.phy.study.no.out, 
                         mod = "Study_Type",
                         xlab = "Standardized mean difference",
                         angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
              theme(text = element_text(size=15),
                    axis.text.y = element_text(size=15),
                    legend.position = "none",
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black")) +labs(tag = "(b)")+
              annotate(geom="text", x=-5, y=3.6,label="Without controlling for phylogeny",hjust=0,
                       size=4.5,fontface="italic")+
              scale_color_manual(values=cbbPalette[1:3])+
              scale_fill_manual(values=cbbPalette[1:3])+expand_limits(y=c(0,4)))

fig.s5c <-ggdraw()+
  draw_plot(orchard_plot(res.study.no.out, 
                         mod = "Study_Type",
                         xlab = "Standardized mean difference",
                         angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
              theme(text = element_text(size=15),
                    axis.text.y = element_text(size=15),
                    panel.border = element_blank(),
                    axis.line = element_line(color="black"),
                    legend.position = c(1,0)) +labs(tag = "(b)")+
              annotate(geom="text", x=-5, y=3.6,label="Without controlling for phylogeny, with viruses",hjust=0,
                       size=4.5,fontface="italic")+
              scale_color_manual(values=cbbPalette[1:3])+
              scale_fill_manual(values=cbbPalette[1:3])+expand_limits(y=c(0,4)))
fig.s5<-plot_grid(fig.s5a,fig.s5b,fig.s5c,
                 align = "hv",hjust=-1,ncol = 1,rel_heights = c(1,1,1))
fig.s5
# ggsave(fig.s5, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S5.jpeg",
#        height=9, width=10, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY STUDY TYPE AND PARASITE GROUP
####################

# make new moderator for the parasite type and study type interaction
es_calculator_no_virus_subset$p.stud<-paste(es_calculator_no_virus_subset$Study_Type,
                                            es_calculator_no_virus_subset$Parasite_Broad,sep = " - ")
# with phylogeny
res.study.par.phy<- rma.mv(yi, vi,mods=~p.stud-1,
                       data=es_calculator_no_virus_subset,
                       random= list(~ 1 | Study_ID/Effect_Size_ID,
                                    ~ 1 | names,
                                    ~ 1 | p.names,
                                    ~ 1 | hostxparasite),
                       R=list(names=R.H.subset,
                              p.names=R.P.subset,hostxparasite=RO.HP.subset),method = "REML",
                       control=list(optimizer="optim", optmethod="BFGS"))
summary.rma(res.study.par.phy)

# without phylogeny
res.study.par.no.phy<- rma.mv(yi, vi,mods=~p.stud-1,
                           data=es_calculator_no_virus_subset,
                           random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary.rma(res.study.par.no.phy)

# without phylogeny with viruses
# make new moderator for the parasite type and study type combo
es_calculator_subset$p.stud<-paste(es_calculator_subset$Study_Type,
                                   es_calculator_subset$Parasite_Broad,sep = " - ")
res.study.par<- rma.mv(yi, vi, mods=~p.stud-1, 
                              data=es_calculator_subset,
                              random= list(~ 1 | Study_ID/Effect_Size_ID),
                              method = "REML")
summary.rma(res.study.par)

####################
# CREATE ORCHARD PLOTS FOR FIG. S6
####################

fig.s6a <-ggdraw()+draw_plot(orchard_plot(res.study.par.phy, mod = "p.stud", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black")) +labs(tag = "(a)")+
                               annotate(geom="text", x=-5, y=7,label="Controlling for phylogeny",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s6b <-ggdraw()+draw_plot(orchard_plot(res.study.par.no.phy, mod = "p.stud", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black")) +labs(tag = "(b)")+
                               annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s6c <-ggdraw()+draw_plot(orchard_plot(res.study.par, mod = "p.stud", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black"),
                                     legend.position = c(1,-0.015)) +labs(tag = "(c)")+
                               annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny, with viruses",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s6<-plot_grid(fig.s6a,fig.s6b,fig.s6c,align='hv',hjust=-1,ncol=1,rel_heights = c(1,1,1))
fig.s6
# ggsave(fig.s6, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S6.jpeg",
#        height=9, width=11.5, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY STUDY TYPE AND PARASITE GROUP WITHOUT OUTLIERS
####################

# with phylogeny
# make new moderator for the parasite type and study type interaction
es_calculator_no_virus_no_out_subset$p.stud<-paste(es_calculator_no_virus_no_out_subset$Study_Type,
                                                   es_calculator_no_virus_no_out_subset$Parasite_Broad,sep = " - ")
# fit multi-level model and check for an effect of study type and parasite type
res.study.par.phy.no.out<- rma.mv(yi, vi,mods=~p.stud-1,
                           data=es_calculator_no_virus_no_out_subset,
                           random= list(~ 1 | Study_ID/Effect_Size_ID,
                                        ~ 1 | names,
                                        ~ 1 | p.names,
                                        ~ 1 | hostxparasite),
                           R=list(names=R.H.subset.no.out,
                                  p.names=R.P.subset.no.out,hostxparasite=RO.HP.subset.no.out),method = "REML",
                           control=list(optimizer="optim", optmethod="BFGS"))
summary.rma(res.study.par.phy.no.out)

# without phylogeny
# make new moderator for the parasite type and study type interaction
es_calculator_no_phy_no_out_subset$p.stud<-paste(es_calculator_no_phy_no_out_subset$Study_Type,
                                                 es_calculator_no_phy_no_out_subset$Parasite_Broad,sep = " - ")
# fit multi-level model and check for an effect of study type and parasite type
res.study.par.no.phy.no.out<- rma.mv(yi, vi,mods=~p.stud-1,
                                  data=es_calculator_no_phy_no_out_subset,
                                  random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary.rma(res.study.par.no.phy.no.out)

# without phylogeny and with viruses
# make new moderator for the parasite type and study type interaction
es_calculator_subset_no_out$p.stud<-paste(es_calculator_subset_no_out$Study_Type,
                                          es_calculator_subset_no_out$Parasite_Broad,sep = " - ")
res.study.par.no.out<- rma.mv(yi, vi, mods=~p.stud-1, 
                       data=es_calculator_subset_no_out,
                       random= list(~ 1 | Study_ID/Effect_Size_ID),
                       method = "REML")
summary.rma(res.study.par.no.out)

####################
# CREATE ORCHARD PLOTS FOR FIG. S7
####################

fig.s7a <-ggdraw()+draw_plot(orchard_plot(res.study.par.phy.no.out, mod = "p.stud", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black")) +labs(tag = "(a)")+
                               annotate(geom="text", x=-5, y=7,label="Controlling for phylogeny",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s7b <-ggdraw()+draw_plot(orchard_plot(res.study.par.no.phy.no.out, mod = "p.stud", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black")) +labs(tag = "(b)")+
                               annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s7c <-ggdraw()+draw_plot(orchard_plot(res.study.par.no.out, mod = "p.stud", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black"),
                                     legend.position = c(1,-0.015)) +labs(tag = "(c)")+
                               annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny, with viruses",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s7<-plot_grid(fig.s7a,fig.s7b,fig.s7c,align='hv',hjust=-1,ncol=1,rel_heights = c(1,1,1))
fig.s7
# ggsave(fig.s7, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S7.jpeg",
#        height=9, width=11.5, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY FITNESS COMPONENT
####################

# with phylogeny
res.fitness.phy<- rma.mv(yi, vi,mods=~Fitness_components-1,
                           data=es_calculator_no_virus_subset,
                           random= list(~ 1 | Study_ID/Effect_Size_ID,
                                        ~ 1 | names,
                                        ~ 1 | p.names,
                                        ~ 1 | hostxparasite),
                           R=list(names=R.H.subset,
                                  p.names=R.P.subset,hostxparasite=RO.HP.subset),method = "REML",
                           control=list(optimizer="optim", optmethod="BFGS"))
summary.rma(res.fitness.phy)

# without phylogeny
res.fitness.no.phy<- rma.mv(yi, vi,mods=~Fitness_components-1,
                         data=es_calculator_no_virus_subset,
                         random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary.rma(res.fitness.no.phy)

# without phylogeny and with viruses
res.fitness<- rma.mv(yi, vi, mods=~Fitness_components-1,
                     data=es_calculator_subset,
                     random= list(~ 1 | Study_ID/Effect_Size_ID),
                     method = "REML")
res.fitness

####################
# CREATE ORCHARD PLOTS FOR FIG. 3
####################

fig.3a <- ggdraw()+draw_plot(
  orchard_plot(res.fitness.phy, mod = "Fitness_components", xlab = "Standardized mean difference",
                        angle = 0,cb=TRUE,alpha = .2)+xlim(-5,5)+
    theme(text = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color="black")) +labs(tag = "(a)")+
    annotate(geom="text", x=-5, y=3.6,hjust=0,
             label="Controlling for phylogeny",size=4.5,fontface="italic")+
    scale_color_manual(values=cbbPalette)+
  scale_fill_manual(values=cbbPalette)+expand_limits(y=c(0,4)))

fig.3b <- ggdraw()+draw_plot(
  orchard_plot(res.fitness.no.phy, mod = "Fitness_components", xlab = "Standardized mean difference",
                        angle = 0,cb=TRUE,alpha = .2)+xlim(-5,5)+
    theme(text = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color="black")) +labs(tag = "(b)")+
    annotate(geom="text", x=-5, y=3.6,hjust=0,
             label="Without controlling for phylogeny",size=4.5,fontface="italic")+
    scale_color_manual(values=cbbPalette)+
  scale_fill_manual(values=cbbPalette)+expand_limits(y=c(0,4)))

fig.3c <- ggdraw()+draw_plot(
  orchard_plot(res.fitness, mod = "Fitness_components", xlab = "Standardized mean difference",
               angle = 0,cb=TRUE,alpha = .2)+xlim(-5,5)+
    theme(text = element_text(size=15),
          axis.text.y = element_text(size=15),
          panel.border = element_blank(),
          axis.line = element_line(color="black"),
          legend.position = c(1,0)) +labs(tag = "(c)")+
  annotate(geom="text", x=-5, y=3.6,hjust=0,
           label="Without controlling for phylogeny, with viruses",size=4.5,fontface="italic")+
    scale_color_manual(values=cbbPalette)+
  scale_fill_manual(values=cbbPalette)+expand_limits(y=c(0,4)))

fig.3<-plot_grid(fig.3a,fig.3b,fig.3c,align='hv',hjust=-1,ncol=1,rel_heights = c(1,1,1))
fig.3
# ggsave(fig.3, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure 3.jpeg",
#        height=9, width=10, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY FITNESS COMPONENT WITHOUT OUTLIERS
####################

# with phylogeny
res.fitness.phy.no.out<- rma.mv(yi, vi,mods=~Fitness_components-1,
                                data=es_calculator_no_virus_no_out_subset,
                                random= list(~ 1 | Study_ID/Effect_Size_ID,
                                             ~ 1 | names,
                                             ~ 1 | p.names,
                                             ~ 1 | hostxparasite),
                                R=list(names=R.H.subset,
                                       p.names=R.P.subset,hostxparasite=RO.HP.subset),method = "REML",
                                control=list(optimizer="optim", optmethod="BFGS"))
summary.rma(res.fitness.phy.no.out)

# without phylogeny
res.fitness.no.phy.no.out<- rma.mv(yi, vi,mods=~Fitness_components-1,
                                   data=es_calculator_no_phy_no_out_subset,
                                   random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary.rma(res.fitness.no.phy.no.out)

# without phylogeny and with viruses
res.fitness.no.out<- rma.mv(yi, vi, mods=~Fitness_components-1,
                            data=es_calculator_subset_no_out,
                            random= list(~ 1 | Study_ID/Effect_Size_ID),
                            method = "REML")
res.fitness.no.out


####################
# CREATE ORCHARD PLOTS FOR FIG. S8
####################

fig.s8a <- ggdraw()+draw_plot(
  orchard_plot(res.fitness.phy.no.out, mod = "Fitness_components", xlab = "Standardized mean difference",
               angle = 0,cb=TRUE,alpha = .2)+xlim(-5,5)+
    theme(text = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color="black")) +labs(tag = "(a)")+
    annotate(geom="text", x=-5, y=3.6,hjust=0,
             label="Controlling for phylogeny",size=4.5,fontface="italic")+
    scale_color_manual(values=cbbPalette)+
    scale_fill_manual(values=cbbPalette)+expand_limits(y=c(0,4)))

fig.s8b <- ggdraw()+draw_plot(
  orchard_plot(res.fitness.no.phy.no.out, mod = "Fitness_components", xlab = "Standardized mean difference",
               angle = 0,cb=TRUE,alpha = .2)+xlim(-5,5)+
    theme(text = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color="black")) +labs(tag = "(b)")+
    annotate(geom="text", x=-5, y=3.6,hjust=0,
             label="Without controlling for phylogeny",size=4.5,fontface="italic")+
    scale_color_manual(values=cbbPalette)+
    scale_fill_manual(values=cbbPalette)+expand_limits(y=c(0,4)))

fig.s8c <- ggdraw()+draw_plot(
  orchard_plot(res.fitness.no.out, mod = "Fitness_components", xlab = "Standardized mean difference",
               angle = 0,cb=TRUE,alpha = .2)+xlim(-5,5)+
    theme(text = element_text(size=15),
          axis.text.y = element_text(size=15),
          panel.border = element_blank(),
          axis.line = element_line(color="black"),
          legend.position = c(1,0)) +labs(tag = "(c)")+
    annotate(geom="text", x=-5, y=3.6,hjust=0,
             label="Without controlling for phylogeny, with viruses",size=4.5,fontface="italic")+
    scale_color_manual(values=cbbPalette)+
    scale_fill_manual(values=cbbPalette)+expand_limits(y=c(0,4)))
fig.s8<-plot_grid(fig.s8a,fig.s8b,fig.s8c,align='hv',hjust=-1,ncol=1,rel_heights = c(1,1,1))
fig.s8
# ggsave(fig.s8, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S8.jpeg",
#        height=9, width=10, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY FITNESS COMPONENT AND PARASITE GROUP
####################

# with phylogeny
# make new moderator for the parasite type and fitness component interaction
es_calculator_no_virus_subset$p.fit<-paste(es_calculator_no_virus_subset$Fitness_components,
                                           es_calculator_no_virus_subset$Parasite_Broad,sep = " - ")
res.fit.par.phy<- rma.mv(yi, vi,mods=~p.fit-1,
                         data=es_calculator_no_virus_subset,
                         random= list(~ 1 | Study_ID/Effect_Size_ID,
                                      ~ 1 | names,
                                      ~ 1 | p.names,
                                      ~ 1 | hostxparasite),
                         R=list(names=R.H.subset,
                                p.names=R.P.subset,hostxparasite=RO.HP.subset),method = "REML",
                         control=list(optimizer="optim", optmethod="BFGS"))
summary.rma(res.fit.par.phy)

# without phylogeny
res.fit.par.no.phy<- rma.mv(yi, vi,mods=~p.fit-1,
                            data=es_calculator_no_virus_subset,
                            random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary.rma(res.fit.par.no.phy)

# without phylogeny and with viruses
# make new moderator for the parasite type and fitness component interaction
es_calculator_subset$p.fit<-paste(es_calculator_subset$Fitness_components,
                                  es_calculator_subset$Parasite_Broad,sep = " - ")
res.fit.par<- rma.mv(yi, vi, mods=~p.fit-1,
                     data=es_calculator_subset,
                     random= list(~ 1 | Study_ID/Effect_Size_ID),
                     method = "REML")
res.fit.par

####################
# CREATE ORCHARD PLOTS FOR FIG. S9
####################

fig.s9a <- ggdraw()+draw_plot(orchard_plot(res.fit.par.phy, mod = "p.fit", xlab = "Standardized mean difference",
                                           angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                                theme(text = element_text(size=15),
                                      axis.text.y = element_text(size=15),
                                      legend.position = "none",
                                      axis.text.x = element_blank(),
                                      axis.title.x = element_blank(),
                                      panel.border = element_blank(),
                                      axis.line = element_line(color="black")) +labs(tag = "(a)")+
                                annotate(geom="text", x=-5, y=7,label="Controlling for phylogeny",hjust=0,
                                         size=4.5,fontface="italic")+
                                scale_color_manual(values=cbbPalette.par)+
                                scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s9b <-ggdraw()+draw_plot(orchard_plot(res.fit.par.no.phy, mod = "p.fit", xlab = "Standardized mean difference",
                                           angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                                theme(text = element_text(size=15),
                                      axis.text.y = element_text(size=15),
                                      legend.position = "none",
                                      axis.text.x = element_blank(),
                                      axis.title.x = element_blank(),
                                      panel.border = element_blank(),
                                      axis.line = element_line(color="black")) +labs(tag = "(b)")+
                                annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny",hjust=0,
                                         size=4.5,fontface="italic")+
                                scale_color_manual(values=cbbPalette.par)+
                                scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s9c <-ggdraw()+draw_plot(orchard_plot(res.fit.par, mod = "p.fit", xlab = "Standardized mean difference",
                                           angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                                theme(text = element_text(size=15),
                                      axis.text.y = element_text(size=15),
                                      panel.border = element_blank(),
                                      axis.line = element_line(color="black"),
                                      legend.position = c(1,-0.015)) +labs(tag = "(c)")+
                                annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny, with viruses",hjust=0,
                                         size=4.5,fontface="italic")+
                                scale_color_manual(values=cbbPalette.par)+
                                scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))
fig.s9<-plot_grid(fig.s9a,fig.s9b,fig.s9c,align='hv',hjust=-1,ncol=1,rel_heights = c(1,1,1))
fig.s9
# ggsave(fig.s9, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S9.jpeg",
#        height=9, width=11.5, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES BY FITNESS COMPONENT AND PARASITE GROUP WITHOUT OUTLIERS
####################

# with phylogeny
# make new moderator for the parasite type and study type interaction
es_calculator_no_virus_no_out_subset$p.fit<-paste(es_calculator_no_virus_no_out_subset$Fitness_components,
                                                  es_calculator_no_virus_no_out_subset$Parasite_Broad,sep = " - ")
# fit multi-level model and check for an effect of fitness component and parasite type
res.fit.par.phy.no.out<- rma.mv(yi, vi,mods=~p.fit-1,
                         data=es_calculator_no_virus_no_out_subset,
                         random= list(~ 1 | Study_ID/Effect_Size_ID,
                                      ~ 1 | names,
                                      ~ 1 | p.names,
                                      ~ 1 | hostxparasite),
                         R=list(names=R.H.subset.no.out,
                                p.names=R.P.subset.no.out,hostxparasite=RO.HP.subset.no.out),method = "REML",
                         control=list(optimizer="optim", optmethod="BFGS"))
summary.rma(res.fit.par.phy.no.out)

# without phylogeny
# make new moderator for the parasite type and study type interaction
es_calculator_no_phy_no_out_subset$p.fit<-paste(es_calculator_no_phy_no_out_subset$Fitness_components,
                                                es_calculator_no_phy_no_out_subset$Parasite_Broad,sep = " - ")
res.fit.par.no.phy.no.out<- rma.mv(yi, vi,mods=~p.fit-1,
                            data=es_calculator_no_phy_no_out_subset,
                            random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary.rma(res.fit.par.no.phy.no.out)

# without phylogeny with viruses
# make new moderator for the parasite type and study type interaction
es_calculator_subset_no_out$p.fit<-paste(es_calculator_subset_no_out$Fitness_components,
                                         es_calculator_subset_no_out$Parasite_Broad,sep = " - ")
res.fit.par.no.out<- rma.mv(yi, vi, mods=~p.fit-1,
                     data=es_calculator_subset_no_out,
                     random= list(~ 1 | Study_ID/Effect_Size_ID),
                     method = "REML")
res.fit.par.no.out

####################
# CREATE ORCHARD PLOTS FOR FIG. S10
####################

fig.s10a <- ggdraw()+draw_plot(orchard_plot(res.fit.par.phy.no.out, mod = "p.fit", xlab = "Standardized mean difference",
                                            angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                                 theme(text = element_text(size=15),
                                       axis.text.y = element_text(size=15),
                                       legend.position = "none",
                                       axis.text.x = element_blank(),
                                       axis.title.x = element_blank(),
                                       panel.border = element_blank(),
                                       axis.line = element_line(color="black")) +labs(tag = "(a)")+
                                 annotate(geom="text", x=-5, y=7,label="Controlling for phylogeny",hjust=0,
                                          size=4.5,fontface="italic")+
                                 scale_color_manual(values=cbbPalette.par)+
                                 scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s10b <-ggdraw()+draw_plot(orchard_plot(res.fit.par.no.phy.no.out, mod = "p.fit", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black")) +labs(tag = "(b)")+
                               annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))

fig.s10c <-ggdraw()+draw_plot(orchard_plot(res.fit.par.no.out, mod = "p.fit", xlab = "Standardized mean difference",
                                          angle = 0,cb=TRUE,alpha=.2)+xlim(-5,5)+
                               theme(text = element_text(size=15),
                                     axis.text.y = element_text(size=15),
                                     panel.border = element_blank(),
                                     axis.line = element_line(color="black"),
                                     legend.position = c(1,-0.015)) +labs(tag = "(c)")+
                               annotate(geom="text", x=-5, y=7,label="Without controlling for phylogeny, with viruses",hjust=0,
                                        size=4.5,fontface="italic")+
                               scale_color_manual(values=cbbPalette.par)+
                               scale_fill_manual(values=cbbPalette.par)+expand_limits(y=c(0,7.5)))
fig.s10<-plot_grid(fig.s10a,fig.s10b,fig.s10c,align='hv',hjust=-1,ncol=1,rel_heights = c(1,1,1))
fig.s10
# ggsave(fig.s10, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S10.jpeg",
#        height=9, width=11.5, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES ACROSS LATITUDE
####################

# create dataframe with only field studies without viruses
es_calculator_no_virus_lat<-subset(es_calculator_no_virus_subset,es_calculator_no_virus_subset$lab.field=="field")
# create dataframe with only field studies with viruses
es_calculator_subset_lat<-subset(es_calculator_subset,es_calculator_subset$lab.field=="field")

# create dataframe with only field studies without viruses and without outliers for analysis with phylogeny
es_calculator_no_virus_lat_no_out<-subset(es_calculator_no_virus_no_out_subset,es_calculator_no_virus_no_out_subset$lab.field=="field")
# create dataframe with only field studies without viruses and without outliers for analysis without phylogeny
es_calculator_no_phy_lat_no_out<-subset(es_calculator_no_phy_no_out_subset,es_calculator_no_phy_no_out_subset$lab.field=="field")
# create dataframe with only field studies with viruses and without outliers
es_calculator_subset_lat_no_out<-subset(es_calculator_subset_no_out,es_calculator_subset_no_out$lab.field=="field")

####################
# subset corr matrices
####################
# phylogeny with outliers
subset.names.lat<-as.vector(unique(es_calculator_no_virus_lat$names)) #subset host cor matrix
R.H.subset.lat<-R.H[subset.names.lat,subset.names.lat]
subset.p.names.lat<-as.vector(unique(es_calculator_no_virus_lat$p.names)) #subset par cor matrix
R.P.subset.lat<-R.P[subset.p.names.lat,subset.p.names.lat]
subset.hostxparasite.lat<-as.vector(unique(es_calculator_no_virus_lat$hostxparasite)) #subset host x par cor matrix
RO.HP.subset.lat<-RO.HP[subset.hostxparasite.lat,subset.hostxparasite.lat]
# phylogeny without outliers
subset.names.lat.no.out<-as.vector(unique(es_calculator_no_virus_lat_no_out$names)) #subset host cor matrix
R.H.subset.lat.no.out<-R.H[subset.names.lat.no.out,subset.names.lat.no.out]
subset.p.names.lat.no.out<-as.vector(unique(es_calculator_no_virus_lat_no_out$p.names)) #subset par cor matrix
R.P.subset.lat.no.out<-R.P[subset.p.names.lat.no.out,subset.p.names.lat.no.out]
subset.hostxparasite.lat.no.out<-as.vector(unique(es_calculator_no_virus_lat_no_out$hostxparasite)) #subset host x par cor matrix
RO.HP.subset.lat.no.out<-RO.HP[subset.hostxparasite.lat.no.out,subset.hostxparasite.lat.no.out]

#####
# run models to test for linear and non-linear effects
#####

#############################################################
# with phylogeny, linear model
res.lat.phy<- rma.mv(yi, vi,mods=~abs(Latitude):Study_Type,
                            data=es_calculator_no_virus_lat,
                            random= list(~ 1 | Study_ID/Effect_Size_ID,
                                         ~ 1 | names,
                                         ~ 1 | p.names,
                                         ~ 1 | hostxparasite),
                            R=list(names=R.H.subset.lat,
                                   p.names=R.P.subset.lat,hostxparasite=RO.HP.subset.lat),method = "REML",
                            control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy)

# with phylogeny, non-linear model
res.lat.phy.quad<- rma.mv(yi, vi,mods=~abs(Latitude):Study_Type+
                            I(abs(Latitude)^2):Study_Type,
                          data=es_calculator_no_virus_lat,
                          random= list(~ 1 | Study_ID/Effect_Size_ID,
                                       ~ 1 | names,
                                       ~ 1 | p.names,
                                       ~ 1 | hostxparasite),
                          R=list(names=R.H.subset.lat,
                                 p.names=R.P.subset.lat,hostxparasite=RO.HP.subset.lat),method = "REML",
                          control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.quad)
AIC(res.lat.phy,res.lat.phy.quad) #model with quadratic term is better
anova(res.lat.phy.quad,btt=5:7) #interaction is significant, split up data
                                
# create dataframe for each study type
es_calculator_no_virus_lat_reprod<-subset(es_calculator_no_virus_lat,
                                          es_calculator_no_virus_lat$Study_Type=="Reproduction")
es_calculator_no_virus_lat_pred<-subset(es_calculator_no_virus_lat,
                                          es_calculator_no_virus_lat$Study_Type=="Predation")
es_calculator_no_virus_lat_comp<-subset(es_calculator_no_virus_lat,
                                          es_calculator_no_virus_lat$Study_Type=="Competition")
###########################
# subset corr matrices

# competition
subset.names.lat.comp<-as.vector(unique(es_calculator_no_virus_lat_comp$names)) #subset host cor matrix
R.H.subset.lat.comp<-R.H[subset.names.lat.comp,subset.names.lat.comp]
subset.p.names.lat.comp<-as.vector(unique(es_calculator_no_virus_lat_comp$p.names)) #subset par cor matrix
R.P.subset.lat.comp<-R.P[subset.p.names.lat.comp,subset.p.names.lat.comp]
subset.hostxparasite.lat.comp<-as.vector(unique(es_calculator_no_virus_lat_comp$hostxparasite)) #subset host x par cor matrix
RO.HP.subset.lat.comp<-RO.HP[subset.hostxparasite.lat.comp,subset.hostxparasite.lat.comp]

# predation
subset.names.lat.pred<-as.vector(unique(es_calculator_no_virus_lat_pred$names)) #subset host cor matrix
R.H.subset.lat.pred<-R.H[subset.names.lat.pred,subset.names.lat.pred]
subset.p.names.lat.pred<-as.vector(unique(es_calculator_no_virus_lat_pred$p.names)) #subset par cor matrix
R.P.subset.lat.pred<-R.P[subset.p.names.lat.pred,subset.p.names.lat.pred]
subset.hostxparasite.lat.pred<-as.vector(unique(es_calculator_no_virus_lat_pred$hostxparasite)) #subset host x par cor matrix
RO.HP.subset.lat.pred<-RO.HP[subset.hostxparasite.lat.pred,subset.hostxparasite.lat.pred]

# reproduction
subset.names.lat.reprod<-as.vector(unique(es_calculator_no_virus_lat_reprod$names)) #subset host cor matrix
R.H.subset.lat.reprod<-R.H[subset.names.lat.reprod,subset.names.lat.reprod]
subset.p.names.lat.reprod<-as.vector(unique(es_calculator_no_virus_lat_reprod$p.names)) #subset par cor matrix
R.P.subset.lat.reprod<-R.P[subset.p.names.lat.reprod,subset.p.names.lat.reprod]
subset.hostxparasite.lat.reprod<-as.vector(unique(es_calculator_no_virus_lat_reprod$hostxparasite)) #subset host x par cor matrix
RO.HP.subset.lat.reprod<-RO.HP[subset.hostxparasite.lat.reprod,subset.hostxparasite.lat.reprod]

######
# model runs

#####
#with phylogeny
#####

# linear model for competition
res.lat.phy.comp<- rma.mv(yi, vi,mods=~abs(Latitude),
                          data=es_calculator_no_virus_lat_comp,
                          random= list(~ 1 | Study_ID/Effect_Size_ID,
                                       ~ 1 | names,
                                       ~ 1 | p.names,
                                       ~ 1 | hostxparasite),
                          R=list(names=R.H.subset.lat.comp,
                                 p.names=R.P.subset.lat.comp,
                                 hostxparasite=RO.HP.subset.lat.comp),method = "REML",
                          control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.comp) 

# non-linear model for competition
res.lat.phy.comp.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                               data=es_calculator_no_virus_lat_comp,
                               random= list(~ 1 | Study_ID/Effect_Size_ID,
                                            ~ 1 | names,
                                            ~ 1 | p.names,
                                            ~ 1 | hostxparasite),
                               R=list(names=R.H.subset.lat.comp,
                                      p.names=R.P.subset.lat.comp,
                                      hostxparasite=RO.HP.subset.lat.comp),method = "REML",
                               control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.comp.quad)
AIC(res.lat.phy.comp,res.lat.phy.comp.quad) #quadratic model isn't better, keep simpler linear model
# set up predicted values for plot
lat.new<-seq(min(es_calculator_no_virus_lat_comp$Latitude),
             max(es_calculator_no_virus_lat_comp$Latitude),length=10)
preds.phy.comp<-predict.rma(res.lat.phy.comp,newmods=lat.new, addx = TRUE, levels = 0) 

#####
# linear model for predation
res.lat.phy.pred<- rma.mv(yi, vi,mods=~abs(Latitude),
                            data=es_calculator_no_virus_lat_pred,
                            random= list(~ 1 | Study_ID/Effect_Size_ID,
                                         ~ 1 | names,
                                         ~ 1 | p.names,
                                         ~ 1 | hostxparasite),
                            R=list(names=R.H.subset.lat.pred,
                                   p.names=R.P.subset.lat.pred,
                                   hostxparasite=RO.HP.subset.lat.pred),method = "REML",
                            control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.pred) 

# non-linear model for predation
res.lat.phy.pred.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                                 data=es_calculator_no_virus_lat_pred,
                                 random= list(~ 1 | Study_ID/Effect_Size_ID,
                                              ~ 1 | names,
                                              ~ 1 | p.names,
                                              ~ 1 | hostxparasite),
                                 R=list(names=R.H.subset.lat.pred,
                                        p.names=R.P.subset.lat.pred,
                                        hostxparasite=RO.HP.subset.lat.pred),method = "REML",
                                 control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.pred.quad)
AIC(res.lat.phy.pred,res.lat.phy.pred.quad) #quadratic model is better
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_pred$Latitude)),
             max(abs(es_calculator_no_virus_lat_pred$Latitude)),length=10)
preds.phy.pred.quad<-predict.rma(res.lat.phy.pred.quad,newmods=cbind(lat.new,lat.new^2), addx = TRUE, levels = 0)

#####
# linear model for reproduction
res.lat.phy.reprod<- rma.mv(yi, vi,mods=~abs(Latitude),
                          data=es_calculator_no_virus_lat_reprod,
                          random= list(~ 1 | Study_ID/Effect_Size_ID,
                                       ~ 1 | names,
                                       ~ 1 | p.names,
                                       ~ 1 | hostxparasite),
                          R=list(names=R.H.subset.lat.reprod,
                                 p.names=R.P.subset.lat.reprod,
                                 hostxparasite=RO.HP.subset.lat.reprod),method = "REML",
                          control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.reprod) 

# non-linear model for reproduction
res.lat.phy.reprod.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                               data=es_calculator_no_virus_lat_reprod,
                               random= list(~ 1 | Study_ID/Effect_Size_ID,
                                            ~ 1 | names,
                                            ~ 1 | p.names,
                                            ~ 1 | hostxparasite),
                               R=list(names=R.H.subset.lat.reprod,
                                      p.names=R.P.subset.lat.reprod,
                                      hostxparasite=RO.HP.subset.lat.reprod),method = "REML",
                               control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.reprod.quad)
AIC(res.lat.phy.reprod,res.lat.phy.reprod.quad) #quadratic model is better
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_reprod$Latitude)),
             max(abs(es_calculator_no_virus_lat_reprod$Latitude)),length=10)
preds.phy.reprod.quad<-predict.rma(res.lat.phy.reprod.quad,newmods=cbind(lat.new,lat.new^2), addx = TRUE, levels = 0)


#####
#without phylogeny
#####

# linear model for competition
res.lat.no.phy.comp<- rma.mv(yi, vi,mods=~abs(Latitude),
                             data=es_calculator_no_virus_lat_comp,
                             random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML",
                             control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.no.phy.comp) 

# non-linear model for competition
res.lat.no.phy.comp.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                                  data=es_calculator_no_virus_lat_comp,
                                  random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML",
                                  control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.no.phy.comp.quad)
AIC(res.lat.no.phy.comp,res.lat.no.phy.comp.quad) #quadratic model is barely better, stick with linear
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_comp$Latitude)),
             max(abs(es_calculator_no_virus_lat_comp$Latitude)),length=10)
preds.no.phy.comp<-predict.rma(res.lat.no.phy.comp,
                               newmods=cbind(lat.new),
                               addx = TRUE, levels = 0)

#####
# linear model for predation
res.lat.no.phy.pred<- rma.mv(yi, vi,mods=~abs(Latitude),
                          data=es_calculator_no_virus_lat_pred,
                          random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML",
                          control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.no.phy.pred) 

# non-linear model for predation
res.lat.no.phy.pred.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                               data=es_calculator_no_virus_lat_pred,
                               random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML",
                               control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.no.phy.pred.quad)
AIC(res.lat.no.phy.pred,res.lat.no.phy.pred.quad) #quadratic model is better
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_pred$Latitude)),
             max(abs(es_calculator_no_virus_lat_pred$Latitude)),length=10)
preds.no.phy.pred.quad<-predict.rma(res.lat.no.phy.pred.quad,
                                      newmods=cbind(lat.new,lat.new^2),
                                      addx = TRUE, levels = 0) 

#####
# linear model for reproduction
res.lat.no.phy.reprod<- rma.mv(yi, vi,mods=~abs(Latitude),
                               data=es_calculator_no_virus_lat_reprod,
                               random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML",
                               control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.no.phy.reprod) 

# non-linear model for reproduction
res.lat.no.phy.reprod.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                                    data=es_calculator_no_virus_lat_reprod,
                                    random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML",
                                    control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.no.phy.reprod.quad)
AIC(res.lat.no.phy.reprod,res.lat.no.phy.reprod.quad) #quadratic model is better
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_reprod$Latitude)),
             max(abs(es_calculator_no_virus_lat_reprod$Latitude)),length=10)
preds.no.phy.reprod.quad<-predict.rma(res.lat.no.phy.reprod.quad,
                                      newmods=cbind(lat.new,lat.new^2),
                                      addx = TRUE, levels = 0)

########################
# test for effect of latitude with data including viruses

# without phylogeny and with viruses, linear model
res.lat<- rma.mv(yi, vi, mods=~abs(Latitude):Study_Type,data=es_calculator_subset_lat,
                 random= list(~ 1 | Study_ID/Effect_Size_ID))
res.lat

# without phylogeny and with viruses, non-linear model
res.lat.quad<- rma.mv(yi, vi, mods=~abs(Latitude):Study_Type+
                        I(abs(Latitude)^2):Study_Type,data=es_calculator_subset_lat,
                 random= list(~ 1 | Study_ID/Effect_Size_ID))
res.lat.quad
AIC(res.lat,res.lat.quad) #model with quadratic term is better
anova(res.lat.quad,btt=5:7) #interaction is significant, split up data

# subset data by study type
es_calculator_subset_lat_reprod<-subset(es_calculator_subset_lat,
                                        es_calculator_subset_lat$Study_Type=="Reproduction")
es_calculator_subset_lat_pred<-subset(es_calculator_subset_lat,
                                        es_calculator_subset_lat$Study_Type=="Predation")
es_calculator_subset_lat_comp<-subset(es_calculator_subset_lat,
                                        es_calculator_subset_lat$Study_Type=="Competition")


#####
#without phylogeny and with viruses
#####
# linear model for competition
res.lat.comp<- rma.mv(yi, vi,mods=~abs(Latitude),
                        data=es_calculator_subset_lat_comp,
                        random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.comp) 

# non-linear model for competition
res.lat.comp.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                             data=es_calculator_subset_lat_comp,
                             random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.comp.quad)
AIC(res.lat.comp,res.lat.comp.quad) #no real difference, stick with simpler model
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_subset_lat_comp$Latitude)),
             max(abs(es_calculator_subset_lat_comp$Latitude)),length=10)
preds.virus.comp<-predict.rma(res.lat.comp,newmods=lat.new,addx = TRUE, levels = 0)

#####
# linear model for predation
res.lat.pred<- rma.mv(yi, vi,mods=~abs(Latitude),
                      data=es_calculator_subset_lat_pred,
                      random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.pred) 

# non-linear model for predation
res.lat.pred.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                           data=es_calculator_subset_lat_pred,
                           random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.pred.quad)
AIC(res.lat.pred,res.lat.pred.quad) #no real difference, stick with simpler model
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_subset_lat_pred$Latitude)),
             max(abs(es_calculator_subset_lat_pred$Latitude)),length=10)
preds.virus.pred<-predict.rma(res.lat.pred,newmods=lat.new,addx = TRUE, levels = 0)

#####
# linear model for reproduction
res.lat.reprod<- rma.mv(yi, vi,mods=~abs(Latitude),
                        data=es_calculator_subset_lat_reprod,
                        random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.reprod) 

# non-linear model for reproduction
res.lat.reprod.quad<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                             data=es_calculator_subset_lat_reprod,
                             random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.reprod.quad)
AIC(res.lat.reprod,res.lat.reprod.quad) #quadratic model is better
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_subset_lat_reprod$Latitude)),
             max(abs(es_calculator_subset_lat_reprod$Latitude)),length=10)
preds.virus.reprod.quad<-predict.rma(res.lat.reprod.quad,
                                     newmods=cbind(lat.new,lat.new^2),
                                     addx = TRUE, levels = 0) #set up predicted values for plot

####################
# CREATE PLOTS FOR FIG. 4
####################

# model with phylogeny

# create dataframe for predicted values for all three species interactions
# add interaction type to individual dataframes
preds.phy.reprod.quad<-as.data.frame(preds.phy.reprod.quad)
preds.phy.pred.quad<-as.data.frame(preds.phy.pred.quad)
preds.phy.comp<-as.data.frame(preds.phy.comp)

preds.phy.reprod.quad$Study_Type<-"Reproduction"
preds.phy.pred.quad$Study_Type<-"Predation"
preds.phy.comp$Study_Type<-"Competition"
preds.phy<-rbind(preds.phy.reprod.quad[,-9],preds.phy.pred.quad[,-9],preds.phy.comp)
colnames(preds.phy)[8]<-"Latitude"
es_calculator_no_virus_lat_all<-rbind(es_calculator_no_virus_lat_comp,
                                      es_calculator_no_virus_lat_reprod,
                                      es_calculator_no_virus_lat_pred)

# plot for non-linear model with phylogeny
mod.plot.phy.quad<-ggplot(es_calculator_no_virus_lat, aes(x = abs(Latitude), y = yi,col=Study_Type,fill=Study_Type)) +
  geom_ribbon(data = preds.phy,aes(x = abs(Latitude), ymax=ci.ub,ymin=ci.lb,col=Study_Type,
                                   fill=Study_Type),alpha=2/5,
              linetype=1,inherit.aes = FALSE)+
  scale_x_continuous(name="Latitude",
                     limits=c(0,80),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80))+
  geom_line(data = preds.phy, aes(x = abs(Latitude), y = pred,col=Study_Type),inherit.aes = FALSE)+
  geom_point(data = es_calculator_no_virus_lat, aes(x = abs(Latitude), y = yi, color=Study_Type, fill = Study_Type),
             shape=21,size=2.75,alpha=.2,color="black")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  ylim(-6,6)+
  xlab(expression(Latitude))+
  ylab("SMD") +
  theme_classic()+
  annotate(geom="text", x=1, y=5,
           label="Controlling for phylogeny",
           size=3,fontface="italic",hjust=0)+
  scale_color_manual(values=cbbPalette,name="Study Type")+
  scale_fill_manual(values=cbbPalette, name="Study Type")+
  theme(text = element_text(size=15),legend.position = "right",aspect.ratio = 1)+labs(tag = "(b)")

#######
# model without phylogeny
#######

# create dataframe for predicted values for all three species interactions
# add interaction type to individual dataframes
preds.no.phy.reprod.quad<-as.data.frame(preds.no.phy.reprod.quad)
preds.no.phy.pred.quad<-as.data.frame(preds.no.phy.pred.quad)
preds.no.phy.comp<-as.data.frame(preds.no.phy.comp)

preds.no.phy.reprod.quad$Study_Type<-"Reproduction"
preds.no.phy.pred.quad$Study_Type<-"Predation"
preds.no.phy.comp$Study_Type<-"Competition"
preds.no.phy<-rbind(preds.no.phy.reprod.quad[,-9],preds.no.phy.pred.quad[,-9],preds.no.phy.comp)
colnames(preds.no.phy)[8]<-"Latitude"

# plot for non-linear model without phylogeny
mod.no.phy.plot.quad<-ggplot(es_calculator_no_virus_lat, aes(x = abs(Latitude), y = yi,col=Study_Type,fill=Study_Type)) +
  geom_ribbon(data = preds.no.phy,aes(x = abs(Latitude), ymax=ci.ub,ymin=ci.lb,col=Study_Type,
                                   fill=Study_Type),alpha=2/5,
              linetype=1,inherit.aes = FALSE)+
  scale_x_continuous(name="Latitude",
                     limits=c(0,80),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80))+
  geom_line(data = preds.no.phy, aes(x = abs(Latitude), y = pred,col=Study_Type),inherit.aes = FALSE)+
  geom_point(data = es_calculator_no_virus_lat, aes(x = abs(Latitude), y = yi, color=Study_Type, fill = Study_Type),
             shape=21,size=2.75,alpha=.2,color="black")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  ylim(-6,6)+
  xlab(expression(Latitude))+
  ylab("SMD") +
  theme_classic()+
  annotate(geom="text", x=1, y=5,
           label="Without controlling for phylogeny",
           size=3,fontface="italic",hjust=0)+
  scale_color_manual(values=cbbPalette,name="Study Type")+
  scale_fill_manual(values=cbbPalette, name="Study Type")+
  theme(text = element_text(size=15),legend.position = "right",aspect.ratio = 1)+labs(tag = "(c)")

#######
# model without phylogeny and with viruses
#######

# create dataframe for predicted values for all three species interactions
# add interaction type to individual dataframes
preds.virus.reprod.quad<-as.data.frame(preds.virus.reprod.quad)
preds.virus.pred<-as.data.frame(preds.virus.pred)
preds.virus.comp<-as.data.frame(preds.virus.comp)

preds.virus.reprod.quad$Study_Type<-"Reproduction"
preds.virus.pred$Study_Type<-"Predation"
preds.virus.comp$Study_Type<-"Competition"
preds.virus<-rbind(preds.virus.reprod.quad[,-9],preds.virus.pred,preds.virus.comp)
colnames(preds.virus)[8]<-"Latitude"
es_calculator_subset_lat_all<-rbind(es_calculator_subset_lat_comp,
                                    es_calculator_subset_lat_reprod,
                                    es_calculator_subset_lat_pred)

# plot for non-linear model without phylogeny and with viruses
mod.plot.quad<-ggplot(es_calculator_subset_lat, aes(x = abs(Latitude), y = yi,col=Study_Type,fill=Study_Type)) +
  geom_ribbon(data = preds.virus,aes(x = abs(Latitude), ymax=ci.ub,ymin=ci.lb,col=Study_Type,
                                      fill=Study_Type),alpha=2/5,
              linetype=1,inherit.aes = FALSE)+
  scale_x_continuous(name="Latitude",
                     limits=c(0,80),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80))+
  geom_line(data = preds.virus, aes(x = abs(Latitude), y = pred,col=Study_Type),inherit.aes = FALSE)+
  geom_point(data = es_calculator_subset_lat, aes(x = abs(Latitude), y = yi, color=Study_Type, fill = Study_Type),
             shape=21,size=2.75,alpha=.2,color="black")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  ylim(-6,6)+
  xlab(expression(Latitude))+
  ylab("SMD") +
  theme_classic()+
  annotate(geom="text", x=1, y=5,
           label="Without controlling for phylogeny, \\nwith viruses",
           size=3,fontface="italic",hjust=0)+
  scale_color_manual(values=cbbPalette,name="Study Type")+
  scale_fill_manual(values=cbbPalette, name="Study Type")+
  theme(text = element_text(size=15),legend.position = "right",aspect.ratio = 1)+labs(tag = "(d)")

# create world map
map.world <-borders("world", colour = "grey", fill = "grey")

mp<- ggplot() + 
  map.world +
  xlab(expression(Longitude))+
  scale_y_continuous(name="Latitude",
                     limits=c(-80,80),
                     breaks=c(-80,-60,-40,-20,0,20,40,60,80),
                     labels=c(-80,-60,-40,-20,0,20,40,60,80))+
  ylab(expression(Latitude))+
  geom_point(data= es_calculator_subset_lat, aes(x = Longitude,
                                                   y = Latitude,
                                                   color = Study_Type),
             size=2.75) +
  theme_classic()+scale_color_manual(values=cbbPalette,
                                     name="Study Type")+
  scale_fill_manual(values=cbbPalette)+theme(text = element_text(size=15),
                                             axis.title.y = element_blank(),
                                             legend.position = "none",plot.margin = unit(c(0,0,0,0),"cm"))+labs(tag = "(a)")
mp

# extract legend

get_legend<-function(myggplot){ #function to extract legend for use across plots
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend<-get_legend(mp+theme(legend.position = "top",
                            text = element_text(size=15),
                            plot.margin=grid::unit(c(0,0,0,0), "mm")))
# create common y label
y.grob <- textGrob("Latitude",gp=gpar(fontsize=15), rot=90)
# add to plot
mp<-grid.arrange(arrangeGrob(mp,left = y.grob))

lat.fig<-plot_grid(mod.plot.phy.quad+theme(legend.position = "none", axis.title.y = element_blank()),
                   NULL,
                   mod.no.phy.plot.quad+theme(legend.position = "none",axis.title.y = element_blank()),
                   NULL,
                   mod.plot.quad+theme(legend.position = "none",axis.title.y = element_blank()),
                   align = c("h"), axis='l',nrow = 1,rel_widths = c(1,0,1,0,1))
# create common y label
y.grob <- textGrob("SMD",gp=gpar(fontsize=15), rot=90)
# add to plot
lat.fig<-grid.arrange(arrangeGrob(lat.fig, left = y.grob))

map.fig<-plot_grid(mp,lat.fig,legend,align = "hv",
                          rel_heights = c(1,1,.09),axis='l',rel_widths = c(1.5,.85,.85),ncol=1)
map.fig
fig.4<-map.fig
# ggsave(fig.4, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure 4.jpeg",
#        height=7, width=9, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES ACROSS LATITUDE WITHOUT OUTLIERS
####################

#######
# with phylogeny linear model
res.lat.phy.no.out<- rma.mv(yi, vi,mods=~abs(Latitude):Study_Type,
                     data=es_calculator_no_virus_lat_no_out,
                     random= list(~ 1 | Study_ID/Effect_Size_ID,
                                  ~ 1 | names,
                                  ~ 1 | p.names,
                                  ~ 1 | hostxparasite),
                     R=list(names=R.H.subset.lat.no.out,
                            p.names=R.P.subset.lat.no.out,hostxparasite=RO.HP.subset.lat.no.out),method = "REML",
                     control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.no.out)

# with phylogeny linear model without interaction
res.lat.phy.no.int.no.out<- rma.mv(yi, vi,mods=~abs(Latitude),
                            data=es_calculator_no_virus_lat_no_out,
                            random= list(~ 1 | Study_ID/Effect_Size_ID,
                                         ~ 1 | names,
                                         ~ 1 | p.names,
                                         ~ 1 | hostxparasite),
                            R=list(names=R.H.subset.lat.no.out,
                                   p.names=R.P.subset.lat.no.out,hostxparasite=RO.HP.subset.lat.no.out),method = "REML",
                            control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.no.int.no.out)
AIC(res.lat.phy.no.out,res.lat.phy.no.int.no.out) #no real difference, keep simpler model without interaction

# with phylogeny, non-linear model without interaction
res.lat.phy.no.int.quad.no.out<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                                 data=es_calculator_no_virus_lat_no_out,
                                 random= list(~ 1 | Study_ID/Effect_Size_ID,
                                              ~ 1 | names,
                                              ~ 1 | p.names,
                                              ~ 1 | hostxparasite),
                                 R=list(names=R.H.subset.lat.no.out,
                                        p.names=R.P.subset.lat.no.out,hostxparasite=RO.HP.subset.lat.no.out),method = "REML",
                                 control=list(optimizer="optim", optmethod="BFGS"))
summary(res.lat.phy.no.int.quad.no.out)

AIC(res.lat.phy.no.int.no.out,res.lat.phy.no.int.quad.no.out) #no difference in AIC, keep simpler linear model
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_no_out$Latitude)),
             max(abs(es_calculator_no_virus_lat_no_out$Latitude)),length=10)
preds.phy.lat.no.out<-predict.rma(res.lat.phy.no.int.no.out,
                                  newmods=lat.new,
                                  addx = TRUE, levels = 0)

#######
# without phylogeny, linear model without interactions
res.lat.no.phy.no.int.no.out<- rma.mv(yi, vi,mods=~abs(Latitude),
                               data=es_calculator_no_virus_lat_no_out,
                               random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.no.phy.no.int.no.out)
# without phylogeny, non-linear without interactions 
res.lat.no.phy.no.int.quad.no.out<- rma.mv(yi, vi,mods=~abs(Latitude)+I(abs(Latitude)^2),
                                    data=es_calculator_no_virus_lat_no_out,
                                    random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
summary(res.lat.no.phy.no.int.quad.no.out)
AIC(res.lat.no.phy.no.int.no.out,res.lat.no.phy.no.int.quad.no.out) #no difference in models, keep simpler linear model
# set up predicted values for plot
lat.new<-seq(min(abs(es_calculator_no_virus_lat_no_out$Latitude)),
             max(abs(es_calculator_no_virus_lat_no_out$Latitude)),length=10)
preds.no.phy.lat.no.out<-predict.rma(res.lat.no.phy.no.int.no.out,
                                     newmods=lat.new,
                                     addx = TRUE, levels = 0) 

# without phylogeny, with viruses, linear model without interactions
res.lat.virus.no.out<- rma.mv(yi, vi, mods=~abs(Latitude),data=es_calculator_subset_lat_no_out,
                         random= list(~ 1 | Study_ID/Effect_Size_ID))
res.lat.virus.no.out
# without phylogeny, with viruses, non-linear model without interactions
res.lat.virus.quad.no.out<- rma.mv(yi, vi, mods=~abs(Latitude)+I(abs(Latitude)^2),data=es_calculator_subset_lat_no_out,
                              random= list(~ 1 | Study_ID/Effect_Size_ID))
res.lat.virus.quad.no.out
AIC(res.lat.virus.no.out,res.lat.virus.quad.no.out) #no difference between models, keep simpler linear model
lat.new<-seq(min(abs(es_calculator_subset_lat_no_out$Latitude)),
             max(abs(es_calculator_subset_lat_no_out$Latitude)),length=10)
preds.lat.virus.no.out<-predict.rma(res.lat.virus.no.out,
                                    newmods=lat.new,
                                    addx = TRUE, levels = 0) #set up predicted values for plot
####################
# CREATE PLOTS FOR FIG. S11
####################
preds.phy.lat.no.out<-as.data.frame(preds.phy.lat.no.out)
colnames(preds.phy.lat.no.out)[8]<-"Latitude"
#######
# model with phylogeny
mod.plot.phy.quad.no.out<-ggplot(es_calculator_no_virus_lat_no_out, aes(x = abs(Latitude), y = yi,col=Study_Type,fill=Study_Type)) +
  geom_ribbon(data = preds.phy.lat.no.out,aes(x = abs(Latitude), ymax=ci.ub,ymin=ci.lb),alpha=2/5,
              linetype=1,inherit.aes = FALSE)+
  scale_x_continuous(name="Latitude",
                     limits=c(0,80),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80))+
  geom_line(data = preds.phy.lat.no.out, aes(x = abs(Latitude),
                                                y = pred),inherit.aes = FALSE)+
  geom_point(data = es_calculator_no_virus_lat_no_out, aes(x = abs(Latitude), y = yi),
             shape=21,size=2.75,alpha=.2,color="black")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  ylim(-6,6)+
  xlab(expression(Latitude))+
  ylab("SMD") +
  theme_classic()+
  annotate(geom="text", x=1, y=5,
           label="Controlling for phylogeny",
           size=3,fontface="italic",hjust=0)+
  scale_color_manual(values=cbbPalette,name="Study Type")+
  scale_fill_manual(values=cbbPalette, name="Study Type")+
  theme(text = element_text(size=15),legend.position = "right",aspect.ratio = 1)+labs(tag = "(b)")

#######
# model without phylogeny
preds.no.phy.lat.no.out<-as.data.frame(preds.no.phy.lat.no.out)
colnames(preds.no.phy.lat.no.out)[8]<-"Latitude"

mod.no.phy.plot.quad.no.out<-ggplot(es_calculator_no_virus_lat_no_out, aes(x = abs(Latitude), y = yi,col=Study_Type,fill=Study_Type)) +
  geom_ribbon(data = preds.no.phy.lat.no.out,aes(x = abs(Latitude), ymax=ci.ub,ymin=ci.lb),alpha=2/5,
              linetype=1,inherit.aes = FALSE)+
  scale_x_continuous(name="Latitude",
                     limits=c(0,80),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80))+
  geom_line(data = preds.no.phy.lat.no.out, aes(x = abs(Latitude),
                                                y = pred),inherit.aes = FALSE)+
  geom_point(data = es_calculator_no_virus_lat_no_out, aes(x = abs(Latitude), y = yi),
             shape=21,size=2.75,alpha=.2,color="black")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  ylim(-6,6)+
  xlab(expression(Latitude))+
  ylab("SMD") +
  theme_classic()+
  annotate(geom="text", x=1, y=5,
           label="Without controlling for phylogeny",
           size=3,fontface="italic",hjust=0)+
  scale_color_manual(values=cbbPalette,name="Study Type")+
  scale_fill_manual(values=cbbPalette, name="Study Type")+
  theme(text = element_text(size=15),legend.position = "right",aspect.ratio = 1)+labs(tag = "(c)")

#######
# model without phylogeny and with viruses
preds.lat.virus.no.out<-as.data.frame(preds.lat.virus.no.out)
colnames(preds.lat.virus.no.out)[8]<-"Latitude"

mod.plot.quad.no.out<-ggplot(es_calculator_subset_lat_no_out, aes(x = abs(Latitude), y = yi,col=Study_Type,fill=Study_Type)) +
  geom_ribbon(data = preds.lat.virus.no.out,aes(x = abs(Latitude), ymax=ci.ub,ymin=ci.lb),alpha=2/5,
              linetype=1,inherit.aes = FALSE)+
  scale_x_continuous(name="Latitude",
                     limits=c(0,80),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80))+
  geom_line(data = preds.lat.virus.no.out, aes(x = abs(Latitude),
                                             y = pred),inherit.aes = FALSE)+
  geom_point(data = es_calculator_subset_lat_no_out, aes(x = abs(Latitude), y = yi),
             shape=21,size=2.75,alpha=.2,color="black")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  ylim(-6,6)+
  xlab(expression(Latitude))+
  ylab("SMD") +
  theme_classic()+
  annotate(geom="text", x=1, y=5,
           label="Without controlling for phylogeny, \\nwith viruses",
           size=3,fontface="italic",hjust=0)+
  scale_color_manual(values=cbbPalette,name="Study Type")+
  scale_fill_manual(values=cbbPalette, name="Study Type")+
  theme(text = element_text(size=15),legend.position = "right",aspect.ratio = 1)+labs(tag = "(d)")


# create world map
mp.no.out<- ggplot() + 
  map.world +
  xlab(expression(Longitude)) +
  scale_y_continuous(name="Latitude",
                     limits=c(-80,80),
                     breaks=c(-80,-60,-40,-20,0,20,40,60,80),
                     labels=c(-80,-60,-40,-20,0,20,40,60,80))+
  ylab(expression(Latitude))+
  geom_point(data= es_calculator_subset_lat_no_out, aes(x = Longitude,
                                                 y = Latitude,
                                                 color = Study_Type),
             size=2.75) +
  theme_classic()+scale_color_manual(values=cbbPalette,
                                     name="Study Type")+
  scale_fill_manual(values=cbbPalette)+theme(text = element_text(size=15),
                                             axis.title.y = element_blank(),
                                             legend.position = "none",plot.margin = unit(c(0,0,0,0),"cm"))+labs(tag = "(a)")
mp.no.out
# extract legend
legend.no.out<-get_legend(mp.no.out+theme(legend.position = "top",
                            text = element_text(size=15),
                            plot.margin=grid::unit(c(0,0,0,0), "mm")))
# create common y label
y.grob <- textGrob("Latitude",gp=gpar(fontsize=15), rot=90)
# add to plot
mp.no.out<-grid.arrange(arrangeGrob(mp.no.out,left = y.grob))

lat.fig.no.out<-plot_grid(mod.plot.phy.quad.no.out+theme(legend.position = "none", axis.title.y = element_blank()),
                   NULL,
                   mod.no.phy.plot.quad.no.out+theme(legend.position = "none",axis.title.y = element_blank()),
                   NULL,
                   mod.plot.quad.no.out+theme(legend.position = "none",axis.title.y = element_blank()),
                   align = c("h"), axis='l',nrow = 1,rel_widths = c(1,0,1,0,1))
# create common y label
y.grob <- textGrob("SMD",gp=gpar(fontsize=15), rot=90)
# add to plot
lat.fig.no.out<-grid.arrange(arrangeGrob(lat.fig.no.out, left = y.grob))

map.fig.no.out<-plot_grid(mp.no.out,lat.fig.no.out,legend,align = "hv",
                   rel_heights = c(1,1,.09),axis='l',rel_widths = c(1.5,.85,.85),ncol=1)
map.fig.no.out
fig.s11<-map.fig.no.out
# ggsave(fig.s11, file="C:/Users/azhas/Dropbox/Meta_Analysis/plots/Figure S11.jpeg",
#        height=7, width=9, units="in", dpi=1500)

####################
# EFFECTS OF PARASITES IN THE LAB VS FIELD
####################

# model with phylogeny
res.loc.phy<-rma.mv(yi, vi,mods=~factor(lab.field)-1,
                data=es_calculator_no_virus,
                random= list(~ 1 | Study_ID/Effect_Size_ID,
                                               ~ 1 | names,
                                               ~ 1 | p.names,
                                               ~ 1 | hostxparasite),
                                  R=list(names=R.H,p.names=R.P,hostxparasite=RO.HP),method = "REML")
res.loc.phy
anova(res.loc.phy, L=c(1,-1))

# model without phylogeny
res.loc.no.phy<-rma.mv(yi, vi,mods=~factor(lab.field)-1,
                    data=es_calculator_no_virus,
                    random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.loc.no.phy
anova(res.loc.no.phy, L=c(1,-1))

# model without phylogeny and with viruses
res.loc<-rma.mv(yi, vi,mods=~factor(lab.field)-1,
                    data=es_calculator,
                    random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.loc
anova(res.loc, L=c(1,-1))

####################
# EFFECTS OF PARASITES IN THE LAB VS FIELD WITHOUT OUTLIERS
####################

# with phylogeny
res.loc.phy.no.out<-rma.mv(yi, vi,mods=~factor(lab.field)-1,
                    data=es_calculator_no_virus_no_out,
                    random= list(~ 1 | Study_ID/Effect_Size_ID,
                                 ~ 1 | names,
                                 ~ 1 | p.names,
                                 ~ 1 | hostxparasite),
                    R=list(names=R.H.no.out,p.names=R.P.no.out,hostxparasite=RO.HP.no.out),method = "REML")
res.loc.phy.no.out

# without phylogeny
res.loc.no.phy.no.out<-rma.mv(yi, vi,mods=~factor(lab.field)-1,
                           data=es_calculator_no_virus_no_out,
                           random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.loc.no.phy.no.out

# without phylogeny and with viruses
res.loc.no.out<-rma.mv(yi, vi,mods=~factor(lab.field)-1,
                data=es_calculator_no_out,
                random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.loc.no.out

# calculate difference in the magnitude of the effects with and without outliers

# with phylogeny
res.loc.phy
res.loc.phy.no.out
(.6837-.6362)/.6837 #difference in the lab
(.3621-.2904)/.3621 #difference in the field

# without phylogeny
res.loc.no.phy
res.loc.no.phy.no.out
(.6964-.6384)/.6964 #difference in the lab
(.3637-.2483)/.3637 #difference in the field

# without phylogeny and with viruses
res.loc
res.loc.no.out
(.6851-.5356)/.6851 #difference in the lab
(.3984-.3082)/.3984 #difference in the field

####################
# EFFECTS OF PARASITES IN EXPERIMENTAL VS OBSERVATIONAL STUDIES
####################

# with phylogeny
res.exob.phy<- rma.mv(yi, vi,mods=~factor(experimental.observational)-1,
                  data=es_calculator_no_virus,
                  random= list(~ 1 | Study_ID/Effect_Size_ID,
                               ~ 1 | names,
                               ~ 1 | p.names,
                               ~ 1 | hostxparasite),
                  R=list(names=R.H,p.names=R.P,hostxparasite=RO.HP),method = "REML")
res.exob.phy
anova(res.exob.phy, L=c(1,-1))

# without phylogeny
res.exob.no.phy<- rma.mv(yi, vi,mods=~factor(experimental.observational)-1,
                      data=es_calculator_no_virus,
                      random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.exob.no.phy
anova(res.exob.no.phy, L=c(1,-1))

# without phylogeny and with viruses
res.exob<-rma.mv(yi, vi,mods=~factor(experimental.observational)-1,
                       data=es_calculator,
                       random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.exob
anova(res.exob, L=c(1,-1))

####################
# EFFECTS OF PARASITES IN EXPERIMENTAL VS OBSERVATIONAL STUDIES WITHOUT OUTLIERS
####################

# with phylogeny
res.exob.phy.no.out<- rma.mv(yi, vi,mods=~factor(experimental.observational)-1,
                      data=es_calculator_no_virus_no_out,
                      random= list(~ 1 | Study_ID/Effect_Size_ID,
                                   ~ 1 | names,
                                   ~ 1 | p.names,
                                   ~ 1 | hostxparasite),
                      R=list(names=R.H.no.out,p.names=R.P.no.out,hostxparasite=RO.HP.no.out),method = "REML")
res.exob.phy.no.out

# without phylogeny
res.exob.no.phy.no.out<- rma.mv(yi, vi,mods=~factor(experimental.observational)-1,
                             data=es_calculator_no_virus_no_out,
                             random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.exob.no.phy.no.out

# without phylogeny and with viruses
res.exob.no.out<-rma.mv(yi, vi,mods=~factor(experimental.observational)-1,
                 data=es_calculator_no_out,
                 random= list(~ 1 | Study_ID/Effect_Size_ID),method = "REML")
res.exob.no.out

# calculate difference in the magnitude of the effects with and without outliers

# with phylogeny
res.exob.phy
res.exob.phy.no.out
(.6374-.5869)/.6374 #difference in experimental studies
(.3767-.2562)/.3767 #difference in observational studies

# without phylogeny
res.exob.no.phy
res.exob.no.phy.no.out
(.6580-.6172)/.6580 #difference in experimental studies
(.3744-.1417)/.3744 #difference in obsevational studies

# without phylogeny and with viruses
res.exob
res.exob.no.out
(.6556-.5291)/.6556 #difference in experimental studies
(.3530-.1650)/.3530 #difference in observational studies