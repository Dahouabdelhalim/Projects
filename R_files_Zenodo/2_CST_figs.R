# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())
graphics.off()

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

library(phyloseq)
library(vegan) 
library(dplyr)
library(data.table)

## clustering
library(cluster)
library(factoextra)
library(RColorBrewer)

## plotting
library(ggplot2) 
library(cowplot)
library(ggpubfigs) #devtools::install_github("JLSteenwyk/ggpubfigs")
library(scales)
library(pals)

## phyloseq object (kmnp_trimmed)
# normalized sequence counts, limited to ASVs shared across at least two samples 
## marked individuals only 
load("Rdata/sifaka_trimmed_phyloseq_normalized_dada2.RData") #phyloseq object: kmnp_trimmed (beta diversity, shared OTUs)
kmnp_trimmed
head(sample_data(kmnp_trimmed))
length(unique(sample_data(kmnp_trimmed)$Name))
#######################################################################
## Community state type (CST) analysis
#######################################################################
## Figure 3: CST PCoA and time course
## Figure S10: CSTs by individual animal
## Figure S12: CST relative abundance
#######################################################################

#######################################################################
#### CST clustering: weighted Unifrac
#######################################################################
## weighted unifrac distances do not produce distinct clusters

# wudist <- phyloseq::distance(kmnp_trimmed, method = "wunifrac")
# ord = ordinate(kmnp_trimmed, method = "PCoA", distance = wudist)
# plot_scree(ord) + xlim(as.character(seq(1,20))) + ggtitle("MDS-WU ordination eigenvalues")
# evs <- ord$value$Eigenvalues
# print(evs[1:20])## First 4 eigenvalues
# 
# pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
# x = phyloseq:::scores.pcoa(ord, display="sites")
# gskmn = clusGap(x[, 1:4], FUN=pam1, K.max = 20, B = 100)
# plot_clusgap(gskmn)
# 
# # visualize optimal number of clusters
# p = fviz_gap_stat(gskmn,
#                   maxSE = list(method = "Tibs2001SEmax"))+ggtitle("") #7 clusters
# p

#######################################################################
# CST clustering: Bray Curtis
#######################################################################
braydist <- phyloseq::distance(kmnp_trimmed, method = "bray")
ord = ordinate(kmnp_trimmed, method = "PCoA", distance = braydist)

## decide which eigenvalues to include
plot_scree(ord) + xlim(as.character(seq(1,20))) + ggtitle("MDS-bray ordination eigenvalues")
evs <- ord$value$Eigenvalues
print(evs[1:20])
h_sub5 <- hist(evs[6:length(evs)], 100)

## keep first 7 eigenvalues
# Gap Statistic
pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
x = phyloseq:::scores.pcoa(ord, display="sites")
gskmn = clusGap(x[, 1:7], FUN=pam1, K.max = 20, B = 100)
plot_clusgap(gskmn)

### Figure S11: visualize optimal number of clusters 
p = fviz_gap_stat(gskmn, 
                  maxSE = list(method = "Tibs2001SEmax"))+ggtitle("") #7 clusters
p
save_plot(p,filename = "figures/SFig11_Bray_Curtis_PAM_clusters.pdf",base_width = 5, base_height = 4)##supplementary figure

## pam clustering with k = 7
set.seed(123)
km.res <- pam(ord$vectors[,1:7], k=7)
head(km.res$cluster, 20)
sample_data(kmnp_trimmed)$CST <- as.factor(km.res$cluster)  # add cst column to sample data 
levels(sample_data(kmnp_trimmed)$CST)<- c("1","2","3","5","4","6","7") # make sure cluster order puts more similar clusters next to each other


#save normalized phyloseq object with CST assignments added
save(kmnp_trimmed, file= "Rdata/sifaka_trimmed_phyloseq_normalized_dada2_CSTs.RData")
# load("Rdata/sifaka_trimmed_phyloseq_normalized_dada2_CSTs.RData")
CSTs <- unique(sample_data(kmnp_trimmed)$CST)
sample_data(kmnp_trimmed)$Date <- as.Date(sample_data(kmnp_trimmed)$Date,"%m/%d/%y")
sample_data(kmnp_trimmed)$Year <- as.factor(sample_data(kmnp_trimmed)$Year)

## limit to Groups I to VI
subset_marked = kmnp_trimmed %>%
  subset_samples(Group != "Bachelor_1" & Group!="Bachelor_2" & Group != "solitary" & Group != "VII") 
subset_marked = prune_taxa(taxa_sums(subset_marked)>0, subset_marked)
unique(sample_data(subset_marked)$Group)
sample_data(subset_marked)$Group = as.factor(sample_data(subset_marked)$Group)
sample_data(subset_marked)$Group <- droplevels(sample_data(subset_marked)$Group)

length(unique(sample_data(subset_marked)$Name))#57
phyloseq::nsamples(subset_marked)#305

#######################################################################
#Figure 3A: stacked bar plot of CST proportions by social group and year

df = data.frame(sample_data(subset_marked))

brewer = brewer.pal(n = 8, name = "RdBu")

#date frame of CST proportions
dfc <- df %>% filter(!(Group %in% c("Bachelor_1","Bachelor_2","VII","solitary"))) %>%
  bind_rows() %>%
  dplyr::group_by(Group,Year,CST)  %>%
  dplyr::tally() %>%
  dplyr::mutate(prop = n / sum(n)) 

p <- ggplot(dfc, aes(x = Year, y = prop, fill = CST)) +
  scale_fill_manual(values=brewer,labels=c("1","2","3","4","5","6","7"))+
  geom_bar(stat = "identity",colour="black") +
  theme_cowplot(font_size = 15)+
  ylab("Proportion")+
  facet_wrap("Group",ncol = 1)+
  background_grid(major="none",minor="none")
p
#######################################################################
##  Fig 3B plot sampling time course for entire period
#######################################################################
samdf <- data.frame(sample_data(subset_marked))
samdf$Date <- as.Date(samdf$Date,"%m/%d/%y")
nrow(data.frame(sample_data(subset_marked)))#305

set.seed(452)
samdf <- samdf %>% 
  arrange(Name,Date)%>%
  group_by(Name,Sex) %>%
  mutate(og_group = dplyr::first(Group, order_by = Date)) %>%
  dplyr::select(Name,Date,Year,CST,og_group,Group) %>%
  group_by(Name,Date)%>%
  sample_n(1)%>%##randomly select one sample per individual per date for timeline
  ungroup()

brewer = brewer.pal(n = 8, name = "RdBu")

p2 <- ggplot(samdf, aes(x = Date, y = Name))+ 
  facet_grid(Group~ Year,scales="free",space = "free")

p3 <- p2 + geom_point(size=4) + 
  ylab("Individual") + 
  xlab("Sampling Date") +
  aes(fill=CST,shape=CST) + 
  scale_fill_manual(values=brewer,labels=c("1","2","3","4","5","6","7"))+
  scale_shape_manual(values=c(21,22,23,24,25,21,22),labels=c("1","2","3","4","5","6","7"))+
  scale_x_date(expand=c(0.2,0.2),date_labels = "%m-%d")+
  theme_cowplot(font_size = 15) 
p3

p4 <- p3 + 
  theme(
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title=element_text(size=14),
    legend.key.height = unit(0.15, "in"),
    legend.key.width = unit(0.15, "in"),
    legend.background = element_blank(),
    plot.margin= unit(c(0.5,0.5,0.5,0.5),"mm"),
    strip.text = element_text(size=12),
    panel.grid.minor = element_line(size = 0.05), 
    panel.grid.major = element_line(size = 0.1,color="black"))+
  panel_border()+
  scale_y_discrete(breaks=samdf$Name)
p4
combined_cst = plot_grid(p+theme(legend.position = "none"),p4,labels = c("a","b"),rel_widths = c(1,3))
combined_cst

CSTColors <-  brewer.pal(n = 8, name = "RdBu")
names(CSTColors) <- unique(CSTs)
CSTColorScale <- scale_colour_manual(name = "CST", values = CSTColors)
CSTFillScale <- scale_fill_manual(name = "CST", values = CSTColors)

braydist <- phyloseq::distance(subset_marked, method="bray")
ord = ordinate(subset_marked, method = "PCoA", distance = braydist)

q = plot_ordination(subset_marked,ord,color="CST") 
q
df = q$data[, 1:2]
head(df)
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, q$data)
names(df)

cst_ord = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  geom_point(size=4,aes(fill=CST,shape=CST),alpha=0.8) + 
  theme_cowplot(font_size = 15)+
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15)) + 
  panel_border() +
  aes(fill=CST,shape=CST) + 
  scale_fill_manual(values=brewer,labels=c("1","2","3","4","5","6","7"))+
  scale_shape_manual(values=c(21,22,23,24,25,21,22),labels=c("1","2","3","4","5","6","7"))+
  background_grid(major = "xy", minor = "none")+
  guides(fill="legend")+
  facet_wrap("Year")+
  xlab("PC 1 (20.1%)")+
  ylab("PC 2 (12.4%)")+
  theme(legend.position = "right")
cst_ord

##### Figure 3
combined_cst = plot_grid(p+theme(legend.position = "none"),p4+theme(legend.position = "none"),labels = c("b","c"),rel_widths = c(1,3))
combined_cst_pcoa = plot_grid(cst_ord,combined_cst,nrow=2,rel_heights = c(1.4,4),labels = c("a",""))
combined_cst_pcoa
save_plot(combined_cst_pcoa,filename = "figures/Fig3_combined_cst_plots_pcoa.pdf",base_width = 16,base_height = 18)
save_plot(combined_cst_pcoa,filename = "figures/Fig3_combined_cst_plots_pcoa.png",dpi=600,base_width = 16,base_height = 18)

#######################################################################
##  Fig S10: CST assignments for individual animals
#######################################################################
q = plot_ordination(subset_marked,ord,color="CST") 
q
df = q$data[, 1:2]
head(df)
colnames(df) = c("Axis_1", "Axis_2")
df = cbind(df, q$data)
names(df)

## limit to longitudinally sampled sifaka
keep = df %>% filter(Sex!="unknown")%>% group_by(Year,Name)%>%
       tally()%>%filter(n>1)%>%distinct(Year,Name)%>% pull(Name)
df = df %>% filter(Name %in% keep)

df$Sex2 = ifelse(df$Sex=="Female","F","M")
f_labels <- data.frame(Name=df$Name, label = paste(df$Sex2,df$Group,sep="-")) %>% distinct()
f_labels2 <- f_labels %>% filter(!(Name %in% c("GJW","WXM"))) ## individuals that switched groups

## immigrants
df %>% filter(Name %in% c("GJW","WXM"))%>%distinct(Name,Group,Sex)
extra = data.frame(Name=c("WXM","GJW"),label=c("M-III,IV","F-II,III"))
f_labels3 = bind_rows(f_labels2,extra)

yearColors = c("#BB5566", "#004488", "#DDAA33")
names(yearColors) <- c("2012","2015","2016")
yearColorScale <- scale_colour_manual(name = "Year", values = yearColors)
yearFillScale <- scale_fill_manual(name = "Year", values = yearColors)

order = df%>%arrange(Group,Name)%>%pull(Name)%>%unique()
length(unique(df$X.SampleID))#303
length(unique(df$Name))#55
ind = ggplot(df,aes(x=Axis.1,y=Axis.2))+
  geom_point(size=4,aes(shape=Year,fill=Year),alpha=0.8) + 
  geom_path(alpha = 0.5, size = 0.5) + # , fill=NA
  geom_text(x = 0.2, y = 0.2, aes(label = label), data = f_labels3)+
  theme_cowplot(font_size = 15)+
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15)) + 
  panel_border() +
  scale_shape_manual(values=c(21,22,24))+
  facet_wrap(~factor(Name,levels=order))+
  background_grid()+
  yearFillScale+
  xlab("PC 1 (20.1%)")+
  ylab("PC 2 (12.4%)")
ind
save_plot(ind,filename = "figures/SFig10_ind_pcoa_all_years.pdf",base_width = 15,base_height = 12)## supplementary figure

#######################################################################
##  Figure S12: CST relative abundance
#######################################################################

### Relative abundance plot of bacterial order
df <- data.frame(sample_data(kmnp_trimmed))
samples <- unique(df$X.SampleID)
samples <- as.vector(samples)
length(samples)#315

kmnp_cst <- subset_taxa(kmnp_trimmed, Phylum!= "Unclassified")
# Prune taxa from the OTU table that are in zero samples
kmnp_cst <- phyloseq::prune_taxa(taxa_sums(kmnp_cst)>0,kmnp_cst) 

data.frame(tax_table(kmnp_cst))%>%distinct(Order)

kmnp_df_cst <- as.data.frame(sample_data(kmnp_cst))
kmnp_df <- setorder(kmnp_df_cst,CST)
sample_data(kmnp_cst) <- kmnp_df
merged_kmnp = merge_samples(kmnp_cst, "CST")
kmnp_order_cst <- merged_kmnp %>% # agglomerate at phylum level
  tax_glom(taxrank = "Order")  %>%                     # Transform to rel. abundance
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%  # Melt to long format
  arrange(Order) %>%
  droplevels()

order_order = kmnp_order_cst %>% group_by(Order) %>%
  dplyr::summarize(mean_Abund = mean(Abundance)) %>%
  dplyr::arrange(desc(mean_Abund))

order_order %>% print(n=80)
nrow(order_order)

low_abund = order_order %>% filter(mean_Abund< 0.005)

rem = low_abund$Order
kmnp_order_cst$Order[kmnp_order_cst$Order %in% rem] <- "< 0.5% abundance"
keep = unique(kmnp_order_cst$Order)[unique(kmnp_order_cst$Order)!="< 0.5% abundance"]
kmnp_order_cst$Order = as.factor(kmnp_order_cst$Order)
kmnp_order_cst$Order = factor(kmnp_order_cst$Order,levels=c(c(keep,"< 0.5% abundance")))
levels(kmnp_order_cst$Order)

col <- levels(as.factor(kmnp_order_cst$Order))
col <- as.factor(col)
length(col)
order_colors = c(
  "#79a5ff",
  "#83d8a5",
  "#01bfbc",
  "#0170a2",
  "#012a94",
  "#b78980",
  "#ff7b6a",
  "#006670",
  "#57003f",
  "#b4c6f4",
  "#ff8dac",
  "#1d224c",
  "#8f8000",
  "#a04fff",
  "#e3c53b",
  "#59001b",
  "#0022ae",
  "#c48300",
  "#a400c9",
  "#009044",
  "#7c7eff",
  "#ff50cf",
  "#006307",
  "#f1001e",
  "#af0031",
  "#d9c67d",
  "#aa006e",
  "#455400")

order_colors <- order_colors[1:length(col)]
names(order_colors) <- levels(col)
colScale <- scale_fill_manual(name="Order",values=order_colors)

labs = c("1","2","3","4","5","6","7")
kmnp_order_cst$CST <- factor(kmnp_order_cst$CST, levels=c("1","2","3","4","5","6","7"))
p <- ggplot(kmnp_order_cst, aes(x = as.factor(CST), y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity",color="black",size=0.5) +
  scale_y_continuous(labels = percent_format())+
  scale_x_discrete(labels=labs,name="CST")+
  colScale + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(colour = "black", size=15),
    axis.title.y = element_text(size=17),
    axis.title.x = element_text(size=17),
    # axis.ticks = element_blank(),
    axis.text.x = element_text(size=14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size=14),
    legend.position="bottom",
    #Manipulating the facet features
    strip.text.x = element_text(size=15, face="bold"),
    strip.text.y = element_text(size=15, face="bold"),
    strip.background = element_rect(colour="black", fill="snow2")) +
  ylab("Relative Abundance")
p
p2 <- p + labs(fill="")   
p2

#phylum plot
kmnp_phylum_cst <- merged_kmnp %>% # agglomerate at phylum level
  tax_glom(taxrank = "Phylum")  %>%                     # Transform to rel. abundance
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%  # Melt to long format
  arrange(Phylum)

kmnp_phylum_cst$Phylum <- as.character(kmnp_phylum_cst$Phylum) #convert to character

phylum_order = kmnp_phylum_cst %>% group_by(Phylum) %>%
  dplyr::summarize(mean_Abund = mean(Abundance)) %>%
  dplyr::arrange(desc(mean_Abund))

phylum_order

Count = length(unique(kmnp_phylum_cst$Phylum))
Count
kmnp_phylum_cst$Phylum = as.factor(kmnp_phylum_cst$Phylum)

col <- levels(as.factor(kmnp_phylum_cst$Phylum))
col <- as.factor(col)
unique(col)

labs = c("1","2","3","4","5","6","7")
kmnp_phylum_cst$CST <- factor(kmnp_phylum_cst$CST, levels=c("1","2","3","4","5","6","7"))
p <- ggplot(kmnp_phylum_cst, aes(x = as.factor(CST), y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity",color="black",size=0.5) +
  scale_y_continuous(labels = percent_format())+
  scale_x_discrete(labels=labs,name="CST")+
  scale_fill_manual(values=tol())+
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(colour = "black", size=15),
    axis.title.y = element_text(size=17),
    axis.title.x = element_text(size=17),
    # axis.ticks = element_blank(),
    axis.text.x = element_text(size=14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size=14),
    legend.position="bottom",
    strip.text.x = element_text(size=15, face="bold"),
    strip.text.y = element_text(size=15, face="bold"),
    strip.background = element_rect(colour="black", fill="snow2")) +
  ylab("Relative Abundance")  # Black rectangle around facet title
p
p4 <- p + labs(fill="")   
p4

combined_rel_abund = plot_grid(p4,p2,nrow=2,labels=c("a","b"))
combined_rel_abund
save_plot(combined_rel_abund,base_width = 12,base_height = 18, filename = "figures/FigS12_combined_rel_abund_plots.png",dpi=300)## Figure S11
