library('tidyverse')
library('phyloseq')
library('cowplot')
library('wesanderson')
library('metagMisc')
library('vegan')
library('ggvegan')
library('directlabels')

annadelHood <- readRDS(file = "phyloseq.rds")
sampleData <- read.csv(file = "sampleData.csv")

annadelHood <- phyloseq_standardize_otu_abundance(annadelHood, method = "total") %>%
  subset_samples(is.na(cPerc) == FALSE) %>%
  subset_samples(is.na(pH) == FALSE)

#### fixing taxonomy

taxa <- annadelHood@tax_table@.Data %>% as.data.frame()

#cryptococcus shows up more than once because it's classified as both
#Cryptococcaceae and Tremellaceae so we have to manually merge them
jointGenera <- subset_taxa(annadelHood, is.na(Genus) == FALSE) %>%
  subset_taxa(Genus != "g__Cryptococcus") %>%
  tax_glom(taxrank = "Genus")
cryptococcus <- subset_taxa(annadelHood, Genus == "g__Cryptococcus") %>%
  tax_glom(taxrank = "Order")
cryptococcusTaxa <- cryptococcus@tax_table@.Data %>% as.data.frame()

rownames(cryptococcusTaxa) <- taxa_names(cryptococcus)
cryptococcusTaxa$Genus <- "g__Cryptococcus"
cryptococcusTaxa$Family <- "g__Tremellaceae"

tax_table(cryptococcus) <- cryptococcusTaxa %>% as.matrix()
jointGenera <- merge_phyloseq(jointGenera, cryptococcus)
jointGeneraTaxa <- jointGenera@tax_table@.Data %>% as.data.frame() 

taxa_names(jointGenera) <- jointGeneraTaxa$Genus %>%
  as.character()

#extract the top 25 taxa
jointCCA <- fantaxtic::get_top_taxa(jointGenera, 25, discard_other = TRUE)

#generate a species table
sppTab <- jointCCA@otu_table@.Data %>% as.data.frame()

#fix rownames
sppTab$sample <- rownames(sppTab)

#check DCA axis lengths
dplyr::select(sppTab, -sample) %>% as.matrix() %>%
  decorana()
#greatest axis length is over 4 so CCA is appropriate

#prepare the matrix by adding the predictor columns
#include factors found important in the GDM
ccaData <- dplyr::select(sampleData, sample, burned, ecosystem, ABS_COVER, cPerc, pH) %>%
  left_join(sppTab, ., by = "sample")

#turn ecosystem and burned into binary categories
ccaData$ecosystem <- (ccaData$ecosystem == "DF") %>% as.numeric() %>% as.factor()
ccaData$burned <- (ccaData$burned == "B") %>% as.numeric() %>% as.factor()

#make a distance matrix
distMatrix <- filter(sampleData, is.na(cPerc) == FALSE, is.na(pH) == FALSE) %>%
  dplyr::select(sample, lon, lat)
rownames(distMatrix) <- distMatrix$sample

distMatrix <- dplyr::select(distMatrix, -sample) %>%
  dist() %>%
  pcnm()

ccaTab <- filter(sppTab, sample %in% sampleData$sample)

rownames(ccaTab) <- ccaTab$sample

ccaTab <- ccaTab %>% dplyr::select(-sample)

#screen the scores to see how many to include
set.seed(900)
distAnova <- cca(ccaTab ~ (scores(distMatrix)))

#check significance of the distance axes
set.seed(900)
distAnova <- anova(distAnova, by = 'axis', permutations = 200)

#only axis 1 is significant at alpha = 0.05

#run the partial cca including axis 1
set.seed(900)
partialCCA <- cca(ccaTab ~ ccaData$burned + ccaData$ecosystem + ccaData$ABS_COVER + ccaData$cPerc + ccaData$pH +
                    Condition(scores(distMatrix, choices = 1)))

#check significance
set.seed(900)
ccaAnova <- partialCCA %>% anova(model = "reduced", by = "mar")

#make separate dataframes
specCCA <- fortify(partialCCA) %>% filter(Score == "species")
biplotCCA <- fortify(partialCCA) %>% filter(Score == "biplot")

#add guild data to species dataframe
specCCA <- annadelHood@tax_table@.Data %>% as.data.frame() %>% dplyr::select(Genus, Phylum, funcGroup) %>%
  distinct() %>%
  left_join(specCCA, ., by = c("Label" = "Genus"))

#fix labels
specCCA$Label <- specCCA$Label %>% str_remove_all("g__")
biplotCCA$Label <- biplotCCA$Label %>% str_remove_all("ccaData\\\\$")
biplotCCA$Label <- biplotCCA$Label %>% str_replace("burned1", "Fire")
biplotCCA$Label <- biplotCCA$Label %>% str_replace("ecosystem1", "Evergreen forest")
biplotCCA$Label <- biplotCCA$Label %>% str_replace("pH", "Soil pH")
biplotCCA$Label <- biplotCCA$Label %>% str_replace("cPerc", "Soil C (%)")

#include only biplot labels that are significant at alpha = 0.05
biplotCCA <- filter(biplotCCA, Label != "ABS_COVER")

#fix the phylum labels
specCCA$Phylum <- specCCA$Phylum %>%
  str_remove_all("p__")

specCCA$funcGroup <- specCCA$funcGroup %>% str_replace_all("mycorrhiza", "Mycorrhiza") %>%
  str_replace_all("pathogen", "Pathogen") %>%
  str_replace_all("saprotroph", "Saprotroph") %>%
  str_replace_all("unclassified", "Other")

ccaPlot_spec <- ggplot() +
  geom_point(aes(x = CCA1, y = CCA2, fill = funcGroup, color = funcGroup, shape = Phylum), size = 5, data = specCCA,
             alpha = 0.7, stroke = 1) +
  geom_dl(aes(x = CCA1, y = CCA2, label=Label), 
          alpha = 1, 
          method = list('top.bumpup', cex=0.75), 
          data = specCCA) +
  geom_segment(aes(x = CA1, y = 0, xend = CCA1, yend = CCA2), 
               data = biplotCCA, alpha = 0.3, 
               color = wes_palette("Rushmore1")[5],
               arrow = arrow(angle = 20, length = unit(0.2, "cm"),
                             ends = "last", type = "open")) +
  scale_fill_manual(values = c(wes_palette("Rushmore1")[3], 
                                "#bababa",
                                wes_palette("Rushmore1")[4], 
                                wes_palette("Rushmore1")[1])) +
  scale_color_manual(values = c(wes_palette("Rushmore1")[3], 
                               "#bababa",
                               wes_palette("Rushmore1")[4], 
                               wes_palette("Rushmore1")[1])) +
  scale_shape_manual(values = c(22, 23, 25)) +
  directlabels::geom_dl(aes(x = CCA1, y = CCA2, label=Label), color = wes_palette("Rushmore1")[5], alpha = 1, method= list(directlabels::last.bumpup, cex = 0.75), data = biplotCCA) +
  theme_minimal() +
  theme(legend.title = element_blank())

#save the cca species plot so it can be knitted together with the nmds
saveRDS(ccaPlot_spec, file = "ccaPlot.rds")
