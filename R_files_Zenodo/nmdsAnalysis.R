library('tidyverse')
library('phyloseq')
library('cowplot')
library('wesanderson')
library('metagMisc')
library('vegan')
library('ggvegan')
library('directlabels')

annadelHood <- readRDS(file = "phyloseq.rds")

#standardize asv abundance
annadelHood <- phyloseq_standardize_otu_abundance(annadelHood, method = "total")

annadelHood <- merge_samples(annadelHood, group = "specSite")
  
set.seed(900)
ord <- annadelHood %>%
  ordinate(method="NMDS", distance="bray") %>%
  scores() %>% data.frame()

ord$specSite <- rownames(ord)

#remake metadata
annadelHood@sam_data$specSite <- rownames(data.frame(annadelHood@sam_data))
annadelHood@sam_data$ecosystem <- rownames(data.frame(annadelHood@sam_data)) %>% str_split(" ") %>%
  sapply(function(x) x[1])
annadelHood@sam_data$burned <- rownames(data.frame(annadelHood@sam_data)) %>% str_split(" ") %>%
  sapply(function(x) x[2])

ord <- full_join(annadelHood@sam_data %>% data.frame(), ord, by = "specSite")

ord$ecosystem <- ord$ecosystem %>% str_replace_all("DF", "Evergreen forest") %>%
  str_replace_all("OW", "Oak woodland")
ord$burned <- ord$burned %>% str_replace_all("B", "Burned") %>%
  str_replace_all("U", "Unburned")

ordPlot <- ggplot(ord, aes(x = NMDS1, y = NMDS2, shape = ecosystem, color = burned)) +
  geom_point(size = 5, alpha = 0.7) +
  scale_shape_manual(values = c(17, 16)) +
  scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
  theme_minimal()

#save the nmds plot so it can be knitted together with the cca plot
saveRDS(ordPlot, file = "ordPlot.rds")

set.seed(900)
adonis(annadelHood@otu_table  ~ burned*ecosystem, data = data.frame(annadelHood@sam_data))

       