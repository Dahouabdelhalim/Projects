# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())
graphics.off()

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

library(dplyr)
library(phyloseq) #Bioconductor
library(vegan)

library(ggpubfigs) #devtools::install_github("JLSteenwyk/ggpubfigs")

#plotting
library(ggplot2)
library(cowplot)
library(RColorBrewer)

## for rooting tree
library(magrittr)
library(data.table)
library(ape) 

## import phyloseq object
## normalized sequence counts, limited to ASVs shared across at least two samples 
## marked individuals only 
load("Rdata/sifaka_trimmed_phyloseq_normalized_dada2.RData") #phyloseq object: kmnp_trimmed (beta diversity, shared OTUs)
head(sample_data(kmnp_trimmed))
length(unique(sample_data(kmnp_trimmed)$Name))#58
#######################################################################
## Principal Coordinates Analysis
#######################################################################
## Figure S5: Effect of DNA extraction on 2012 samples
## Figure 1: Bray-Curtis PCOA colored by year
## Figure 2: Bray-Curtis PCOA colored by social group
## Figure S6: weighted Unifrac PCOA colored by year
## Figure S9: weighted Unifrac PCOA colored by social group

#######################################################################
## Figure S5: check if 2012 samples cluster by extraction method
#######################################################################
subset_2012 = kmnp_trimmed %>%
  subset_samples(Year=="2012"&Sex!="unknown") 
subset_2012 = prune_taxa(taxa_sums(subset_2012)>0, subset_2012)
length(unique(sample_data(subset_2012)$Name))#35
phyloseq::nsamples(subset_2012)#140
braydist <- phyloseq::distance(subset_2012, method="bray")
ord = ordinate(subset_2012, method = "PCoA", distance = braydist)
p = plot_ordination(subset_2012,ord,color="Extraction") 
p
df = p$data[, 1:2]
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, p$data)

r = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  geom_point(size=3,aes(fill=Extraction,color=Extraction,shape=Group),alpha=0.7) + 
  scale_shape_manual(values = c(21,22,23,24,25,22))+
  scale_fill_manual(values=brewer.pal(3,"Paired")[2:3])+
  scale_color_manual(values=brewer.pal(3,"Paired")[2:3])+
  theme_bw(base_size = 15)+
  background_grid(major = "xy", minor = "none")+
  xlab("PC 1 (21.7%)")+
  ylab("PC 2 (16.2%)")+
  theme(text = element_text(size = 15),axis.text=element_text(size=15), axis.title=element_text(size=16),
        strip.text = element_text(size=17),
        strip.background = element_blank(),
        legend.background = element_blank())
r
save_plot(r,filename = "figures/SFig5_2012_pcoa_extraction_method.pdf",base_width = 7,base_height = 5)
#######################################################################
## Figure 1: BC PCoA by year
#######################################################################

sample_data(kmnp_trimmed)$Date <- as.Date(sample_data(kmnp_trimmed)$Date,"%m/%d/%y")
sample_data(kmnp_trimmed)$Year <- as.factor(sample_data(kmnp_trimmed)$Year)

## limit to longitudinally sampled animals
unique(sample_data(kmnp_trimmed)$Group)

## limit to Groups I to VI
subset_marked = kmnp_trimmed %>%
  subset_samples(Group != "Bachelor_1" & Group!="Bachelor_2" & Group != "VII") 
length(unique(sample_data(subset_marked)$Name))#57
phyloseq::nsamples(subset_marked)#305
unique(sample_data(subset_marked)$Group)
subset_marked = prune_taxa(taxa_sums(subset_marked)>0, subset_marked)

sample_data(subset_marked)$Group = as.factor(sample_data(subset_marked)$Group)
sample_data(subset_marked)$Group <- droplevels(sample_data(subset_marked)$Group)

braydist <- phyloseq::distance(subset_marked, method="bray")
ord = ordinate(subset_marked, method = "PCoA", distance = braydist)
nsamples(subset_marked)
length(unique(sample_data(subset_marked)$Name))

p = plot_ordination(subset_marked,ord,color="Year",shape="Group") 
p
df = p$data[, 1:2]
head(df)
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, p$data)

yearColors = c("#BB5566", "#004488", "#DDAA33")
names(yearColors) <- c("2012","2015","2016")
yearColorScale <- scale_colour_manual(name = "Year", values = yearColors)
yearFillScale <- scale_fill_manual(name = "Year", values = yearColors)

## figure 1A
p = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  stat_ellipse(geom = "polygon",type="t",alpha=0.1, aes(group=Year,fill=Year),lty="dashed")+
  stat_ellipse(type="t",alpha=0.8, aes(group=Year,color=Year),lty="dashed",lwd=0.7)+
  geom_point(size=3,alpha=0.6,aes(fill=Year,shape=Year))+
  scale_shape_manual(values=c(21,22,24))+
  yearFillScale +
  yearColorScale+
  theme_bw(base_size = 15)+
  background_grid(major = "xy", minor = "none")+
  xlab("PC 1 (20.1%)")+
  ylab("PC 2 (12.4%)")+
  theme(text = element_text(size = 15),axis.text=element_text(size=15), axis.title=element_text(size=16),
        strip.text = element_text(size=17),strip.background = element_blank(),
        legend.background = element_blank())+
  theme(legend.position = c(0.1,0.85))
p

braydist <- phyloseq::distance(subset_marked, method="bray")
ord = ordinate(subset_marked, method = "PCoA", distance = braydist)
### Figure 1B
q <- plot_ordination(subset_marked, ord, color="Year") + 
  yearColorScale + 
  facet_wrap("Group") + 
  geom_density2d(contour_var = "ndensity",alpha=0.7) + 
  theme_bw(base_size = 15)+
  background_grid(major = "none", minor = "xy")+
  xlab("PC 1 (20.1%)")+
  ylab("PC 2 (12.4%)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
        strip.text = element_text(size=17),strip.background = element_blank(),
        legend.background = element_blank())+
  theme(legend.position = "none")
q

combined = plot_grid(p,q,nrow=1,labels=c("a","b"))
combined
save_plot(combined,filename = "figures/Fig1_bray_curtis_pcoa_year.pdf",base_width = 10,base_height = 5)
save_plot(combined,filename = "figures/Fig1_bray_curtis_pcoa_year.png",dpi=600,base_width = 10,base_height = 5)
#######################################################################
## Figure 2 BC PCoA by social group (Groups I to VI only)
#######################################################################
p = plot_ordination(subset_marked, ordinate(subset_marked, distance="bray", "PCoA"),color="Group",shape="Group") 
p
df = p$data[, 1:2]
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, p$data)

p = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  geom_point(size=3,aes(fill=Group,shape=Group),alpha=0.7) + 
  scale_shape_manual(values = c(21,22,23,24,25,22))+
  background_grid(major = "xy", minor = "none")+
  scale_fill_manual(values=friendly_pal("bright_seven"))+
  guides(fill="legend")+
  facet_wrap("Year")+
  xlab("PC 1 (20.1%)")+
  ylab("PC 2 (12.4%)")+
  theme(legend.position = "none")+
  theme_bw(base_size = 15)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=16),
        strip.text = element_text(size=17),strip.background = element_blank(),
        legend.background = element_blank())
p


# Propithecus verreauxi is listed as a critically endangered species. 
# Publishing exact geographic coordinates could adversely impact conservation efforts or 
# introduce unintended risk to the species. 
# home range boundaries in Figure 2B are available upon reasonable request to the corresponding author

##combines PCoA and home range figures
# source("Fig2_home_range_maps.R")
# source("Fig2_home_range_maps_updated.R") 
# multi # combined pcoa and HR maps
# save_plot(multi,filename = "figures/Fig2_group_pcoa_home_range_plots_2012_2016.pdf",base_width = 15,base_height = 9)
# save_plot(multi,filename = "figures/Fig2_group_pcoa_home_range_plots_2012_2016.png",dpi=600,base_width = 15,base_height = 9)

#######################################################################
## Figure S6: WU PCoA by year
#######################################################################
## root tree
# https://john-quensen.com/r/unifrac-and-tree-roots/
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }

my.tree <- phy_tree(subset_marked)
out.group <- pick_new_outgroup(my.tree)
out.group 
new.tree <- ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
phy_tree(subset_marked) <- new.tree
phy_tree(subset_marked) 

## Figure S6A
## if tree is not rooted, phyloseq::distance will randomly assign an outgroup
wunifracdist <- phyloseq::distance(subset_marked, method="wunifrac")
ord = ordinate(subset_marked, method = "PCoA", distance = wunifracdist)
p <- plot_ordination(subset_marked, ord,shape="Year",label ="Name") 
p
df = p$data[, 1:2]
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, p$data)

p = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  stat_ellipse(geom = "polygon",type="t",alpha=0.1, aes(group=Year,fill=Year),lty="dashed")+
  stat_ellipse(type="t",alpha=0.8, aes(group=Year,color=Year),lty="dashed",lwd=0.7)+
  geom_point(size=3,alpha=0.7,aes(fill=Year,shape=Year))+
  scale_shape_manual(values=c(21,22,24))+
  yearFillScale +
  yearColorScale+
  theme_bw(base_size = 15)+
  background_grid(major = "xy", minor = "none")+
  xlab("PC 1 (33.3%)")+
  ylab("PC 2 (18.1%)")+
  theme(text = element_text(size = 15),axis.text=element_text(size=15), axis.title=element_text(size=16),
        strip.text = element_text(size=17),strip.background = element_blank(),
        legend.background = element_blank())+
  theme(legend.position = c(0.8,0.2))
p

### Figure S6B
q <- plot_ordination(subset_marked, ord, color="Year") + 
  yearColorScale + 
  facet_wrap("Group") + 
  geom_density2d(contour_var = "ndensity",alpha=0.7) + 
  theme_bw(base_size = 15)+
  background_grid(major = "none", minor = "xy")+
  xlab("PC 1 (33.3%)")+
  ylab("PC 2 (18.1%)")+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=16),
        strip.text = element_text(size=17),strip.background = element_blank(),
        legend.background = element_blank())+
  theme(legend.position = "none")
q
combined = plot_grid(p,q,nrow=1,labels=c("a","b"),rel_widths = c(1,1))
save_plot(combined,filename = "figures/SFig6_wunifrac_pcoa_year.pdf",base_width = 15,base_height = 6)##supplementary figure

#######################################################################
## Figure S9: WU PCoA by social group
#######################################################################
# set.seed(200)## setting the seed is not necessary b/c tree is already rooted
p = plot_ordination(subset_marked, ordinate(subset_marked, distance="wunifrac", "PCoA"),
                    color="Group",shape="Group") 
p
df = p$data[, 1:2]
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, p$data)
nsamples(subset_marked)#305
length(unique(sample_data(subset_marked)$Name))#57

p = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  scale_fill_manual(values=friendly_pal("bright_seven"))+
  scale_shape_manual(values = c(21,22,23,24,25,22))+
  geom_point(size=3,aes(fill=Group,shape=Group),alpha=0.8) + 
  background_grid(major = "xy", minor = "none")+
  guides(fill="legend")+
  facet_wrap("Year")+
  xlab("PC 1 (33.3%)")+
  ylab("PC 2 (18.1%)")+
  theme(legend.position = "none")+
  theme_bw(base_size = 15)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=16),
        strip.text = element_text(size=17),strip.background = element_blank(),
        legend.background = element_blank()) 
p
save_plot("figures/SFig9_wunifrac_group_pcoa_all_years.pdf",p,base_height = 5, base_width = 12)
