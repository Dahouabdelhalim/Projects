# ----- Load libraries ----- 
# install.packages("readr")
library(readr) # Loading csv files
# install.packages("vegan")
library(vegan) # PERMANOVA and CCA analyses
# install.packages("iNEXT")
library(iNEXT) # Rarefaction and Extrapolation methods
# install.packages("ggplot2")
library(ggplot2) # Create figures
# install.packages("ggpubr")
library(ggpubr) # Arrange multiple figures
# install.packages("ade4")
library(ade4) # Correspondence analysis
# install.packages("factoextra")
library(factoextra) # Plot ordination output
# install.packages("cluster")
library(cluster) # Hierarchical agglomerative clustering
# install.packages("ape")
library(ape) # Load and edit Newick trees
# install.packages("phylobase")
library(phylobase) # Create phylo4d objects
# install.packages("phylosignal")
library(phylosignal) # Calculate phylogenetic signal


####################################################################################################
######    Pooling interception records across regions: partial Constrained Correspondence    #######
######    Analysis (CCA) of species' commodity profiles in each region                       #######
####################################################################################################

# ------ Import data ------
# The proportion of interceptions per HS-2 commodity in each region, for the 59 species intercepted more than 20 times in two or more regions
spe_profiles_allregions <- read.csv("spe_profiles_allregions.csv", row.names = "species")

# The species and region each species-by-region commodity profile belongs to
species_region<- read.csv("species_region.csv")


# ------ Format species and interception region data ------
### Table of the species-by-region combinations per species
species_profiles <- xtabs(formula = intercepted ~ species_region + species, data = species_region)

### Table of the species-by-region combinations per interception region 
region_profiles <- xtabs(formula = intercepted ~ species_region + region, data = species_region)


# ------ Partial CCA ------ 
# Estimate the variance in commodity associations explained by species, once the effect of interception region is removed 

cca<- cca(spe_profiles_allregions ~ species_profiles + Condition(region_profiles)) # The interception region of each commodity profile as a conditioning variable

# ANOVA like permutation test for CCA
anova.cca(cca) 



####################################################################################################
###### Estimate species richness and Shannon diversity transported by each commodity class  ########
####################################################################################################

# ------ Import data ------
# The number of interceptions per species on each commodity class
spe_comclass <- read.csv("spe_comclass.csv", row.names = "species")


# ------ Calculate richness and diversity estimates ------
# Estimate species richness per commodity class
chao_rich<- ChaoRichness(spe_comclass, datatype = "abundance", conf = 0.95)
chao_rich$class<- as.character(rownames(chao_rich))
names(chao_rich)[4:5] <- c("lower", "upper") # Rename columns for plotting
chao_rich[order(chao_rich$Estimator),] # Order by richness

# Estimate Shannon diversity per commodity class
chao_SH<- ChaoShannon(spe_comclass, datatype = "abundance", conf = 0.95)
chao_SH$class<- as.character(rownames(chao_SH))
names(chao_SH)[4:5] <- c("lower", "upper") # Rename columns for plotting
chao_SH[order(chao_SH$Estimator),] # Order by diversity


# ------ Plot data ------
# Species richness per commodity class: separate columns into observed and estimated richness
chao_rich$obs_est<- chao_rich$Estimator - chao_rich$Observed 
obs<- as.data.frame(chao_rich[c(6,1,4,5)]) 
obs$species<- rep("Observed", length(obs$class)) 
obs_est<- as.data.frame(chao_rich[c(6,7,4,5)])
obs_est$species<- rep("Estimated", length(obs_est$class)) 
obs_est[c(3,4)] <- NA
names(obs_est) <- c("class","Observed","lower","upper","species")
richness<- rbind(obs,obs_est) 

# Create barplot
Richness_plot<- ggplot(data = richness, aes(x = reorder(class, -Observed), y = Observed, fill = species)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y ="Number of species", x = "") + scale_y_continuous(breaks = seq(0, 13000, by = 1500)) + ggtitle("a)") + geom_errorbar(aes(ymin= lower, ymax= upper), width=.2) + theme(legend.position = "none", legend.title = element_blank()) + theme(text=element_text(size=14))

# Shannon diversity per commodity class: separate columns into observed and estimated diversity
chao_SH$obs_est<- chao_SH$Estimator - chao_SH$Observed
obs<- as.data.frame(chao_SH[c(6,1,4,5)]) 
obs$diversity<- rep("Observed", length(obs$class)) 
obs_est<- as.data.frame(chao_SH[c(6,7,4,5)])
obs_est$diversity<- rep("Estimated", length(obs_est$class)) 
obs_est[c(3,4)] <- NA
names(obs_est) <- c("class","Observed","lower","upper","diversity")
diversity<- rbind(obs,obs_est)

# Create barplot
Diversity_plot<- ggplot(data = diversity, aes(x = reorder(class, -Observed), y = Observed, fill = diversity)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) + labs(y ="Shannon diversity", x = "") + ggtitle("b)") + geom_errorbar(aes(ymin= lower, ymax= upper), width=.2) + theme(text=element_text(size=14))


# Figure 2 #
Figure_2 <- ggpubr::ggarrange(Richness_plot, Diversity_plot, common.legend = TRUE, legend = c("right"))
plot(Figure_2)


####################################################################################################
######  Estimate the species richness and Shannon diversity transported by HS-2 commodities   ######
######  classed as plant products and wood products                                           ######
####################################################################################################

# ------ Import data ------
# The number of interceptions per species on each HS-2 commodity group within plant products and wood products
spe_HS2 <-  read.csv("spe_HS2.csv", row.names = "species") 

colnames(spe_HS2) <- c("10 Cereals", "09 Coffee/tea/spices","45 Cork", "11 Flours", "08 Fruit/nuts","13 Gum/resin","06 Live plants/cut flowers","48 Paper", "46 Plaiting materials", "49 Printed matter", "12 Seeds/grains/medicinal","53 Vegetable fibres", "14 Vegetable products/bamboo", "07 Vegetables", "44 Wood/wood articles") # Change column names for figure labels

# ------ Calculate richness and diversity estimates ------
# Estimate species richness per HS-2 commodity group
chao_rich2<- ChaoRichness(spe_HS2, datatype = "abundance", conf = 0.95)
chao_rich2$name<- as.character(rownames(chao_rich2))
names(chao_rich2)[4:5] <- c("lower", "upper") # Rename columns for plotting
chao_rich2[order(chao_rich2$Estimator),] # Order by diversity

# Estimate Shannon diversity per HS-2 commodity group
chao_SH2<- ChaoShannon(spe_HS2, datatype = "abundance", conf = 0.95)
chao_SH2$name<- as.character(rownames(chao_SH2))
names(chao_SH2)[4:5] <- c("lower", "upper") # Rename columns for plotting
chao_SH2[order(chao_SH2$Estimator),] # Order by diversity


# ------ Plot data ------
# Species richness per HS-2 commodity: separate columns into observed and estimated richness
chao_rich2$obs_est<- chao_rich2$Estimator - chao_rich2$Observed
obs<- as.data.frame(chao_rich2[c(6,1,4,5)])
obs$species<- rep("Observed", length(obs$name))
obs_est<- as.data.frame(chao_rich2[c(6,7,4,5)])
obs_est$species<- rep("Estimated", length(obs_est$name))
obs_est[c(3,4)] <- NA
names(obs_est) <- c("name","Observed","lower","upper","species")
richness2<- rbind(obs,obs_est)

# Create barplot
Richness_plot2<- ggplot(data = richness2, aes(x = reorder(name, -Observed), y = Observed, fill = species)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y ="Number of species", x = "") + scale_y_continuous(breaks = seq(0, 13000, by = 1500)) + ggtitle("a)") + geom_errorbar(aes(ymin= lower, ymax= upper), width=.2) + theme(legend.position = "none", legend.title = element_blank()) + theme(text=element_text(size=14)) + theme(plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))

# Shannon diversity per HS-2 commodity: separate columns into observed and estimated richness
chao_SH2$obs_est<- chao_SH2$Estimator - chao_SH2$Observed
obs<- as.data.frame(chao_SH2[c(6,1,4,5)]) 
obs$diversity<- rep("Observed", length(obs$name))
obs_est<- as.data.frame(chao_SH2[c(6,7,4,5)])
obs_est$diversity<- rep("Estimated", length(obs_est$name))
obs_est[c(3,4)] <- NA
names(obs_est) <- c("name","Observed","lower","upper","diversity")
diversity2<- rbind(obs,obs_est)

# Create barplot
Diversity_plot2 <- ggplot(data = diversity2, aes(x = reorder(name, -Observed), y = Observed, fill = diversity)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) + labs(y ="Shannon diversity", x = "") + ggtitle("b)") + geom_errorbar(aes(ymin= lower, ymax= upper), width=.2) + theme(text=element_text(size=14), legend.title = element_blank()) + theme(plot.margin = unit(c(0.2,0.2,0.2,1), "cm"))


# Figure 3 #
Figure_3 <- ggpubr::ggarrange(Richness_plot2, Diversity_plot2, common.legend = TRUE, legend = "right")
Figure_3 + theme(plot.margin=margin(10,10,10,20)) # Larger margins for figure labels



####################################################################################################
######     The taxonomic composition of insect orders transported with the same commodity     ######
######     classes in the cargo, baggage and mail pathways                                    ######
####################################################################################################

# ------ Import data ------
# Interceptions per insect order on the five commodity classes found in all three pathways (animal products, foodstuffs, machinery/electrical, plant products and wood products)
orders_pathway_comclass<- read.csv("orders_pathway_comclass.csv", row.names = 1)


# ------ Plot data ------
# Set colours for each insect order
order_cols<- c("Blattodea" = "#093186", "Coleoptera" = "#005230", "Dermaptera" = "#67A300", "Diptera" = "#AAFF99", "Hemiptera" = "#7A0002", "Hymenoptera" =  "#D33D22", "Lepidoptera" = "#FF999B", "Orthoptera" = "#85008F", "Psocodea" = "#D80085", "Thysanoptera" = "#FFA400", "Zygentoma" = "#8C7F62", "Mantodea" = "#D4AFB9", "Phasmida" = "#D1CFE2", "Neuroptera" = "#9CADCE", "Odonata" = "#7EC4CF", "Siphonaptera" = "#52B2CF", "Embioptera" = "#4F86C6", "Trichoptera" = "#744FC6", "Plecoptera" = "#379392", "Mecoptera" = "#17301C") 


# Figure 1 #
Figure_1 <- ggplot(data = orders_pathway_comclass, mapping = aes(x = pathway, y = interceptions, fill = Order)) + geom_bar(stat="identity", position = "fill") + facet_grid(~ class) + theme_minimal() + scale_fill_manual(values = order_cols) + ylab("Proportion of interceptions") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10, colour = "black")) + xlab("") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.text.x = element_text(size = 10)) + scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + theme(panel.spacing = unit(2, "lines"))
plot(Figure_1)



# ------ Calculate Permutational Multivariate Analysis of Variance (PERMANOVA) ------
# Compare the differences in taxonomic composition between commodity classes, and between pathways 

# Interceptions per insect order on each commodity class & pathway combination
tab<- as.data.frame.matrix(xtabs(formula = interceptions ~ com_path + Order, data = orders_pathway_comclass)) 

# Independent variables: the pathway and the commodity class of each commodity class & pathway combination
tab2<- unique(orders_pathway_comclass[,c(1,4,5)][order(orders_pathway_comclass$com_path),]) 
rownames(tab2) <- rownames(tab) # Set rownames to match contingency table (tab)

permanova<- adonis2(formula = tab ~ class + pathway, data = tab2, method = "bray", by = "terms", permutations = 9999) # method = "bray" uses Brayâ€“Curtis dissimilarity for pairwise distances, by = "terms" assesses the significance for each term sequentially 



####################################################################################################
######  Correspondence Analysis (CA) and hierarchical agglomerative clustering of species'    ######
######  commodity associations                                                                ######
####################################################################################################

# ------ Import  data ------
# Commodity profiles of species with > 20 interceptions: the proportion of interceptions per HS-2 commodity group
spe_commodity_profiles<- read.csv("spe_commodity_profiles.csv", row.names = "species")


# ------ Correspondence analysis of species' commodity profiles ------
coa_std <- dudi.coa(spe_commodity_profiles, scannf = FALSE, nf = 8) # Retain the eight first CA axes


# ------ Hierarchical agglomerative clustering of commodity profiles ------
# Species' commodity profiles based on CA coordinates 
species_profile<- coa_std$li 

hc <- agnes(species_profile, method = "ward")
hc$ac # Agglomerative coefficient - closer to 1 means stronger clustering structure

# Permutation test to determine whether non-random levels of clustering are present - Greenacre, M. (2011). A Simple Permutation Test for Clusteredness. Economics Working Papers, (April)
"rtest.hclust" <- function(li, method = "ward.D", nsim = 999){
  d <- dist(li)
  hc <- hclust(d, method = "ward.D")
  obs <- hc$height
  
  sim <- matrix(0, nrow = nsim, ncol = length(obs))
  for(i in 1:nsim){
    li_sim <- apply(li, 2, function(x) sample(x, replace = FALSE))
    hc_sim <- hclust(dist(li_sim), method = "ward.D")
    sim[i,] <- hc_sim$height
  }
  
  pval <- numeric(0)
  for(i in 1:length(obs)){
    p <- sum(sim[,i] < obs[i]) / nsim
    pval <- c(pval, p)
  }
  s <- 0.05 / (length(obs))
  n <- sum(pval < s)
  print(paste("% of significative nodes :", n / length(obs) * 100))
  
  index <- which(pval < s)
  if(length(index) == 0){
    cut <- 0
  }	else{
    m <- max(index)
    cut <- mean(obs[m:(m+1)])				
  }
  return(cut)
} 

# Define the cutting point according to the permutation test introduced by Greenacre (2011)
cut_r <- rtest.hclust(coa_std$li)
plot(hc)
abline(h = cut_r, col = "red", lty = 2) # 11 clusters optimal

clusters <- as.data.frame(cutree(hc, k = 11)) # Cut tree 
table(clusters) # The number of species per cluster


# ------ Plot data ------
# Load colours for HS-2 commodity groups - coloured by commodity class
class_cols <- read.csv("class_cols.csv", row.names = "species") 

# Load colours for species - coloured by order
order_cols<- read.csv("order_cols.csv", row.names = "species")

# Load colours for species - coloured by cluster
cluster_cols <- read.csv("cluster_cols.csv", row.names = "species")


# Figure 4 #
options(ggrepel.max.overlaps = Inf) # Allow overlapping labels

# a) Ordination output with commodities coloured by commodity class:
Figure_4a <- fviz_ca_col(coa_std, axes = c(1, 2), label ="col", arrows = c(FALSE, FALSE),  repel = TRUE, col.col = as.factor(rownames(class_cols)), col.row = "#3C3C3C", title = "", pointsize = "contrib", alpha.row = 0.25, labelsize = 6) + theme(legend.position = "none") + scale_colour_manual(values = class_cols$cols_class) + ggtitle("a)")
plot(Figure_4a)

# b) Ordination output with species coloured by cluster:
rownames(coa_std$co)[c(41,77,16,31,79,42,44,76,9)] <- c("06 Live plants/cut flowers","07 Vegetables","09 Coffee/tea/spices","08 Fruit/nuts","44 Wood/wood articles","84 Machinery","16 Meat/crustacean preparations","14 Vegetable products","69 Ceramics") # Rename commodities for plot labels

Figure_4b <- fviz_ca_biplot(coa_std, label ="col", arrows = c(FALSE, FALSE),  repel = TRUE, col.col = "#3C3C3C", col.row = as.factor(rownames(cluster_cols)), title = "", pointsize = "contrib", labelsize = 4, select.col = list(name = c("06 Live plants/cut flowers","07 Vegetables","09 Coffee/tea/spices","08 Fruit/nuts","44 Wood/wood articles","84 Machinery","16 Meat/crustacean preparations", "14 Vegetable products", "69 Ceramics"))) + theme(legend.position = "none") + scale_colour_manual(values = as.character(cluster_cols$cluster_cols), aes(alpha = 0.1)) + ggtitle("b)") 
plot(Figure_4b)

# c) Ordination output with species coloured by order:
Figure_4c <- fviz_ca_row(coa_std,axes = c(1, 2), label ="none", arrows = c(FALSE, FALSE),  repel = TRUE, col.col = "#3C3C3C", col.row = as.factor(rownames(order_cols)), title = "", pointsize = 0.7) + theme(legend.position = "none") + scale_colour_manual(values = as.character(order_cols$order_cols)) + ggtitle("c)")
plot(Figure_4c)
add.scatter.eig(coa_std$eig, nf = 7, posi = c("bottomleft"), ratio = 0.10, xax = 1, yax = 2) # Add eigenvalues to plot



####################################################################################################
######           Test for phylogenetic signal to species' commodity associations              ######
####################################################################################################

# ------ Import data ------
# Tree of species based on taxonomy
taxo_tree<- read.tree("taxo_tree.txt")
taxo_tree$tip.label <- gsub("_", " ",taxo_tree$tip.label) # Replace underscore to match species names from CA

# Species' commodity profiles from CA coordinates
species_profile<- coa_std$li

# Create phylo4 object
p4d <- phylo4d(taxo_tree, species_profile)


# Tree of 150 species based on phylogeny from Timetree.org
phylo_tree<- read.tree("phylo_tree.nwk")
phylo_tree$tip.label <- gsub("_", " ",phylo_tree$tip.label) # replace underscore to match species names from CA

# Species' commodity profiles from CA coordinates
species_profile_150<- species_profile[rownames(species_profile) %in% phylo_tree$tip.label,] 

# phylo4 object
p4d_timetree <- phylo4d(phylo_tree, species_profile_150)


# ------ Calculate phylogenetic signal  ------
# Abouheif's Cmean statistic for each CA axis

# Based on the taxonomic tree
phylo<- phyloSignal(p4d = p4d, method = "Cmean")

# Based on the phylogenetic tree
phylo2<- phyloSignal(p4d = p4d_timetree, method = "Cmean") 



####################################################################################################
######      Similarities in commodity associations within higher taxa: partial CCA of the     ######
######     variance in species' commodity associations explained by order, family and genus   ######
####################################################################################################

# ------ Import data ------
# Which order each species belongs to 
order_profiles<- read.csv("order_profiles.csv", row.names = "species")

# Which family each species belongs to - families where commodity associations are available for  only one species are excluded
family_profiles<- read.csv("family_profiles.csv", row.names = "species")

# Which genus each species belongs to - genera where commodity associations are available for only one species are excluded
genus_profiles<- read.csv("genus_profiles.csv", row.names = "species")

# The contingency table  used earlier (spe_commodity_profiles) includes additional genera with no species identified to species level representing at least one additional species, exclude these here
index <- grep(" ", rownames(spe_commodity_profiles))
spe_commodity_profiles2 <- spe_commodity_profiles[index,]


# ------ CCA of species' commodity profiles ------
# Variance explained by order
cca_orders <- cca(spe_commodity_profiles2, order_profiles)

# Variance explained by family
spe_commodity_profiles3 <- spe_commodity_profiles2[rownames(spe_commodity_profiles2) %in% rownames(family_profiles),] # Exclude species with no other members in the same family

cca_families <- cca(spe_commodity_profiles3, family_profiles) 

# Variance explained by genus
spe_commodity_profiles4<- spe_commodity_profiles2[rownames(spe_commodity_profiles2) %in% rownames(genus_profiles),] # Exclude species with no other members in the same genus

cca_genera <- cca(spe_commodity_profiles4, genus_profiles) 