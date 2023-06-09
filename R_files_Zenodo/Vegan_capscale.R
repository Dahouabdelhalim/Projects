if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(vegan)
library(qiime2R)
#library(MASS)

feature_table<-read_qza("C:/Gresham/Project_Gravimondo/BGS_final/qiime_results/table-dada2.qza")
info_data <-feature_table$data

Environmental_variables <- read.delim("C:/Gresham/Project_Gravimondo/Boris_resubmit/Supplemental_code/COSantana_2020/env_variables.csv")

my_varechem = env_variables
my_varespec = taxa_counts_transposed
rownames(my_varespec) <- my_varespec$OTU_ID
my_varespec <- my_varespec[-c(1)]

#Using Mantel

Salinity_DM <- vegdist(env_variables[,3])
Temperature_DM <- vegdist(env_variables[,4])
OrganicMatter_DM <- vegdist(env_variables[,5])
taxa_dist <- vegdist(my_varespec) # Creates distance matrix of the OTUs

taxa_salinity_Mantel <- mantel(taxa_dist, Salinity_DM, permutations = 999)
taxa_salinity_Mantel
taxa_temp_Mantel <- mantel(taxa_dist, Temperature_DM, permutations = 999)
taxa_temp_Mantel
taxa_OM_Mantel <- mantel(taxa_dist, OrganicMatter_DM, permutations = 999)
taxa_OM_Mantel

#Using NMDS, envfit
feature_table<-read_qza("C:/Gresham/Project_Gravimondo/Project_Impact/table-dada2.qza")
info_data <-feature_table$data

Environmental_variables <- read.delim("C:/Gresham/Project_Gravimondo/Boris_resubmit/Supplemental_code/COSantana_2020/env_variables.csv")

Otu_veg <- t(info_data) # Samples become rows and OTUs are columns
Otus.pca <- rda(Otu_veg)
OTUs_dist <- vegdist(Otu_veg) # Creates distance matrix of the OTUs#Plot PCA
PCA_plot<-ordiplot(Otus.pca)
PCA_plot# Plotting environmental variables
Envdt <- Environmental_variables[,2:12]
Tax_env <- envfit(OTUs_dist, Envdt, permu=999)
Tax_env
plot(Tax_env, p.max = 0.1)
# add labels to the samples in the ordiplot
orditorp(PCA_plot, "site", pch="+", pcol="grey")

#pdf('C:/Gresham/Project_Gravimondo/Boris_resubmit/figure_5/NMDS_test_taxa.pdf')
ord <- metaMDS(Otu_veg)
#ord <- metaMDS((info_data), try=1000, k = 3)
(fit <- envfit(ord, Envdt, perm = 999))

priSite <- diversity(my_varespec, index = "invsimpson", MARGIN = 1)

plot(ord)
orditorp(ord, display = "sites", priority = priSite, scaling = 3,
         col = "blue", cex = 1, pch = 19)

scores(fit, "vectors")

plot(fit)
#plot(fit, p.max = 0.9, col = "blue")
#plot(fit, p.max = 0.05, col = "red")
dev.off()

print(fit)

Otu_veg <- t(info_data) # Samples become rows and OTUs are columns
Envdt <- Environmental_variables[,2:12]
ord <- metaMDS(Otu_veg)
(fit <- envfit(ord, Envdt, perm = 999))
scores(fit, "vectors")
plot(ord)
plot(fit)
plot(fit, p.max = 0.05, col = "red")

cc_info_data <- complete.cases(info_data)
df1_info_data <- info_data[cc_info_data, ]

df1_info_data <- info_data[1:10, ]
###
feature_table<-read_qza("C:/Gresham/Project_Gravimondo/quarantine/table-dada2.qza")
info_data <-feature_table$data

Environmental_variables <- read.delim("C:/Gresham/Project_Gravimondo/quarantine/env_variables.csv")

Otu_veg <- t(info_data) # Samples become rows and OTUs are columns
Envdt <- Environmental_variables[,2:5]

vare.cap <- capscale(Otu_veg ~ (Salinity + Organic_Matter + Temperature + Water_content) , Envdt,
                     dist="bray", sqrt.dist= TRUE,  metaMDS = TRUE)
vare.cap
pdf('C:/Gresham/Project_Gravimondo/Boris_resubmit/figure_5/_capscale.pdf')
plot(vare.cap)
dev.off()
anova(vare.cap)


library(BiocManager)
install("RCM")

library(phyloseq)
data(Zeller)