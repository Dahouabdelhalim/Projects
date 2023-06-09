##### PACKAGES + DATA #####

library(tidyverse)
library(ggplot2)
library(iNEXT)
library(vegan)
library(randomcoloR)
library(colorRamps)
library(stringr)
library(sp)
library(cowplot)
library(patchwork)
library(rstan)
library(brms)
library(tidybayes)
library(ggridges)
library(ggstance)


bali_COI <- read.csv("data/COI_OTU_Taxonomy.csv")
bali_18S <- read.csv("data/18S_OTU_Taxonomy.csv")
coralnet <- read_csv("data/coralnet.csv")
bali_meta <- read.csv("data/Bali_Metadata.csv")


########*#########
###### COI #######
########*#########

### GET RID OF SINGLETONS

COI.nosingleton <- bali_COI %>% 
  mutate(rowsum = rowSums(.[2:28])) %>%
  filter(rowsum > 1) %>%
  select(-rowsum)

# number of OTUs annotated with Phylum "Unknown" - 15,770

length(which(COI.nosingleton$Phylum == "Unknown"))

# total number of sequence reads per ARMS fraction and overarching sum

COI.total <- COI.nosingleton %>%
  select(-c(OTUID, Phylum)) %>%
  summarize_all(funs(sum)) %>%
  mutate(sum = rowSums(.))


##### RAREFACTION CURVES #####

### RAREFACTION BY 500 FRACTION

# select 500 fraction samples

COI.500 <- COI.nosingleton %>%
  select(contains("F500"))

COI.500[c(1:9),c(1:9)]

rownames(COI.500) <- COI.500$OTUID

COI.list500 <- as.data.frame(t(COI.500))
str(COI.list500)

COI.list500 <- as.data.frame(COI.list500)

COI.list500[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

COI.list500 <- split(COI.list500, rownames(COI.list500))

COI.list500$INDO12S1AF500 <- t(COI.list500$INDO12S1AF500)
COI.list500$INDO12S1BF500 <- t(COI.list500$INDO12S1BF500)
COI.list500$INDO12S1CF500 <- t(COI.list500$INDO12S1CF500)

COI.list500$INDO12S2AF500 <- t(COI.list500$INDO12S2AF500)
COI.list500$INDO12S2BF500 <- t(COI.list500$INDO12S2BF500)
COI.list500$INDO12S2CF500 <- t(COI.list500$INDO12S2CF500)

COI.list500$INDO13S1AF500 <- t(COI.list500$INDO13S1AF500)
COI.list500$INDO13S1BF500 <- t(COI.list500$INDO13S1BF500)
COI.list500$INDO13S1CF500 <- t(COI.list500$INDO13S1CF500)

# rarefaction curves by samples

COI.rf.500 <- iNEXT(COI.list500, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
COI.rf.500

COI.rare.500 <- fortify(COI.rf.500, type=1)
head(COI.rare.500)

COI.rare.pt.500 <-COI.rare.500[which(COI.rare.500$method=="observed"),]
COI.rare.ln.500 <-COI.rare.500[which(COI.rare.500$method!="observed"),]
COI.rare.ln.500$method <- factor(COI.rare.ln.500$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
COI.500.points <- COI.rare.500 %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
COI.500.rarefaction <- ggplot(COI.rare.500, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=COI.rare.ln.500) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = COI.500.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO12S1AF500="coral1", INDO12S1BF500="coral1", INDO12S1CF500="coral1", 
                                 INDO12S2AF500="gold", INDO12S2BF500="gold", INDO12S2CF500="gold", 
                                 INDO13S1AF500="deepskyblue", INDO13S1BF500="deepskyblue", INDO13S1CF500="deepskyblue")) +
  scale_fill_manual(values = c(INDO12S1AF500="coral1", INDO12S1BF500="coral1", INDO12S1CF500="coral1", 
                               INDO12S2AF500="gold", INDO12S2BF500="gold", INDO12S2CF500="gold", 
                               INDO13S1AF500="deepskyblue", INDO13S1BF500="deepskyblue", INDO13S1CF500="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
COI.500.rarefaction

# ggsave("plots/COI_Rare_500.pdf", COI.500.rarefaction, width = 8, height = 6, useDingbats = FALSE)


### RAREFACTION BY 100 FRACTION

# select 100 fraction samples

COI.100 <- COI.nosingleton %>%
  select(contains("F100"))

COI.100[c(1:9),c(1:9)]

rownames(COI.100) <- COI.100$OTUID

COI.list100 <- as.data.frame(t(COI.100))
str(COI.list100)

COI.list100 <- as.data.frame(COI.list100)

COI.list100[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

COI.list100 <- split(COI.list100, rownames(COI.list100))

COI.list100$INDO12S1AF100 <- t(COI.list100$INDO12S1AF100)
COI.list100$INDO12S1BF100 <- t(COI.list100$INDO12S1BF100)
COI.list100$INDO12S1CF100 <- t(COI.list100$INDO12S1CF100)

COI.list100$INDO12S2AF100 <- t(COI.list100$INDO12S2AF100)
COI.list100$INDO12S2BF100 <- t(COI.list100$INDO12S2BF100)
COI.list100$INDO12S2CF100 <- t(COI.list100$INDO12S2CF100)

COI.list100$INDO13S1AF100 <- t(COI.list100$INDO13S1AF100)
COI.list100$INDO13S1BF100 <- t(COI.list100$INDO13S1BF100)
COI.list100$INDO13S1CF100 <- t(COI.list100$INDO13S1CF100)

# rarefaction curves by samples

COI.rf.100 <- iNEXT(COI.list100, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
COI.rf.100

COI.rare.100 <- fortify(COI.rf.100, type=1)
head(COI.rare.100)

COI.rare.pt.100 <-COI.rare.100[which(COI.rare.100$method=="observed"),]
COI.rare.ln.100 <-COI.rare.100[which(COI.rare.100$method!="observed"),]
COI.rare.ln.100$method <- factor(COI.rare.ln.100$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
COI.100.points <- COI.rare.100 %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
COI.100.rarefaction <- ggplot(COI.rare.100, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=COI.rare.ln.100) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = COI.100.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO12S1AF100="coral1", INDO12S1BF100="coral1", INDO12S1CF100="coral1", 
                                 INDO12S2AF100="gold", INDO12S2BF100="gold", INDO12S2CF100="gold", 
                                 INDO13S1AF100="deepskyblue", INDO13S1BF100="deepskyblue", INDO13S1CF100="deepskyblue")) +
  scale_fill_manual(values = c(INDO12S1AF100="coral1", INDO12S1BF100="coral1", INDO12S1CF100="coral1", 
                               INDO12S2AF100="gold", INDO12S2BF100="gold", INDO12S2CF100="gold", 
                               INDO13S1AF100="deepskyblue", INDO13S1BF100="deepskyblue", INDO13S1CF100="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
COI.100.rarefaction

# ggsave("plots/COI_Rare_100.pdf", COI.100.rarefaction, width = 8, height = 6, useDingbats = FALSE)


### RAREFACTION BY SESSILE FRACTION

# select sessile fraction samples

COI.SES <- COI.nosingleton %>%
  select(contains("SES"))

COI.SES[c(1:9),c(1:9)]

rownames(COI.SES) <- COI.SES$OTUID

COI.listSES <- as.data.frame(t(COI.SES))
str(COI.listSES)

COI.listSES <- as.data.frame(COI.listSES)

COI.listSES[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

COI.listSES <- split(COI.listSES, rownames(COI.listSES))

COI.listSES$INDO12S1ASES <- t(COI.listSES$INDO12S1ASES)
COI.listSES$INDO12S1BSES <- t(COI.listSES$INDO12S1BSES)
COI.listSES$INDO12S1CSES <- t(COI.listSES$INDO12S1CSES)

COI.listSES$INDO12S2ASES <- t(COI.listSES$INDO12S2ASES)
COI.listSES$INDO12S2BSES <- t(COI.listSES$INDO12S2BSES)
COI.listSES$INDO12S2CSES <- t(COI.listSES$INDO12S2CSES)

COI.listSES$INDO13S1ASES <- t(COI.listSES$INDO13S1ASES)
COI.listSES$INDO13S1BSES <- t(COI.listSES$INDO13S1BSES)
COI.listSES$INDO13S1CSES <- t(COI.listSES$INDO13S1CSES)

# rarefaction curves by samples

COI.rf.SES <- iNEXT(COI.listSES, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
COI.rf.SES

COI.rare.SES <- fortify(COI.rf.SES, type=1)
head(COI.rare.SES)

COI.rare.pt.SES <-COI.rare.SES[which(COI.rare.SES$method=="observed"),]
COI.rare.ln.SES <-COI.rare.SES[which(COI.rare.SES$method!="observed"),]
COI.rare.ln.SES$method <- factor(COI.rare.ln.SES$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
COI.SES.points <- COI.rare.SES %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
COI.SES.rarefaction <- ggplot(COI.rare.SES, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=COI.rare.ln.SES) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = COI.SES.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO12S1ASES="coral1", INDO12S1BSES="coral1", INDO12S1CSES="coral1", 
                                 INDO12S2ASES="gold", INDO12S2BSES="gold", INDO12S2CSES="gold", 
                                 INDO13S1ASES="deepskyblue", INDO13S1BSES="deepskyblue", INDO13S1CSES="deepskyblue")) +
  scale_fill_manual(values = c(INDO12S1ASES="coral1", INDO12S1BSES="coral1", INDO12S1CSES="coral1", 
                               INDO12S2ASES="gold", INDO12S2BSES="gold", INDO12S2CSES="gold", 
                               INDO13S1ASES="deepskyblue", INDO13S1BSES="deepskyblue", INDO13S1CSES="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
COI.SES.rarefaction

# ggsave("plots/COI_Rare_SES.pdf", COI.SES.rarefaction, width = 8, height = 6, useDingbats = FALSE)


### RAREFACTION BY ARMS

COI.nosingleton[c(1:20),c(1:20)]

# merge COI data with metadata

COI.long <- COI.nosingleton %>%
  select(-Phylum) %>%
  pivot_longer(names_to = "COI", values_to = "value", -OTUID) %>%
  inner_join(bali_meta)

# merge fractions

COI.merged <- COI.long %>% group_by(Year, Site, Unit, OTUID) %>% summarize(summed.value = sum(value)) %>%
  mutate(ARMSID = paste("INDO",Year,Site,Unit, sep = "")) %>%
  ungroup() %>%
  select(-Year, -Site, -Unit) %>%
  pivot_wider(names_from = OTUID, values_from = summed.value)

# total sequence reads for each ARMS

COI.count <- COI.merged %>%
  select(-ARMSID) %>%
  mutate(sums = rowSums(.))
COI.count$sums
  
# total OTUs for each ARMS

COI.otu <- COI.long %>% group_by(Year, Site, Unit, OTUID) %>% summarize(summed.value = sum(value)) %>%
  mutate(ARMSID = paste("INDO",Year,Site,Unit, sep = "")) %>%
  ungroup() %>%
  select(-Year, -Site, -Unit) %>%
  mutate(pa.otu = case_when(summed.value > 0 ~ 1,
                            TRUE ~ 0)) %>%
  group_by(ARMSID) %>%
  summarize(otu.count = sum(pa.otu))

# prepare data for rarefaction

COI.list <- as.data.frame(COI.merged)
rownames(COI.list) <- COI.list$ARMSID
COI.listARMS <- COI.list[-1]

COI.listARMS[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

COI.listARMS <- split(COI.listARMS, rownames(COI.listARMS))

COI.listARMS$INDO2012S1A <- t(COI.listARMS$INDO2012S1A)
COI.listARMS$INDO2012S1B <- t(COI.listARMS$INDO2012S1B)
COI.listARMS$INDO2012S1C <- t(COI.listARMS$INDO2012S1C)

COI.listARMS$INDO2012S2A <- t(COI.listARMS$INDO2012S2A)
COI.listARMS$INDO2012S2B <- t(COI.listARMS$INDO2012S2B)
COI.listARMS$INDO2012S2C <- t(COI.listARMS$INDO2012S2C)

COI.listARMS$INDO2013S1A <- t(COI.listARMS$INDO2013S1A)
COI.listARMS$INDO2013S1B <- t(COI.listARMS$INDO2013S1B)
COI.listARMS$INDO2013S1C <- t(COI.listARMS$INDO2013S1C)

# rarefaction curves by samples

COI.rf.ARMS <- iNEXT(COI.listARMS, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
COI.rf.ARMS

COI.rare.ARMS <- fortify(COI.rf.ARMS, type=1)
head(COI.rare.ARMS)

COI.rare.pt.ARMS <-COI.rare.ARMS[which(COI.rare.ARMS$method=="observed"),]
COI.rare.ln.ARMS <-COI.rare.ARMS[which(COI.rare.ARMS$method!="observed"),]
COI.rare.ln.ARMS$method <- factor(COI.rare.ln.ARMS$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
COI.ARMS.points <- COI.rare.ARMS %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
COI.ARMS.rarefaction <- ggplot(COI.rare.ARMS, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=COI.rare.ln.ARMS) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = COI.ARMS.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO2012S1A="coral1", INDO2012S1B="coral1", INDO2012S1C="coral1", 
                                 INDO2012S2A="gold", INDO2012S2B="gold", INDO2012S2C="gold", 
                                 INDO2013S1A="deepskyblue", INDO2013S1B="deepskyblue", INDO2013S1C="deepskyblue")) +
  scale_fill_manual(values = c(INDO2012S1A="coral1", INDO2012S1B="coral1", INDO2012S1C="coral1", 
                               INDO2012S2A="gold", INDO2012S2B="gold", INDO2012S2C="gold", 
                               INDO2013S1A="deepskyblue", INDO2013S1B="deepskyblue", INDO2013S1C="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
COI.ARMS.rarefaction

# ggsave("plots/COI_Rare_ARMS.pdf", COI.ARMS.rarefaction, width = 8, height = 6, useDingbats = FALSE)


##### RELATIVE ABUNDANCES ##### 

head(COI.nosingleton)

COI.rra <- COI.nosingleton %>%
  select(-Phylum) %>%
  pivot_longer(names_to = "ARMS", values_to = "seq", -OTUID) %>%
  group_by(ARMS) %>%
  mutate(rra = seq/sum(seq)) %>%
  pivot_wider(id_cols = ARMS, names_from = OTUID, values_from = rra, values_fill = list(rra = 0))

head(COI.rra)


##### MDS #####

# run nmds ordination

COI.nmds <- metaMDS(COI.rra[-1], k = 3, metric = "bray", trymax = 1000)
stressplot(COI.nmds)
COI.nmds

# plot ordination with ggplot

COI.scores <- as.data.frame(scores(COI.nmds))
COI.scores$fraction <- COI.rra$ARMS
COI.scores1 <- as.data.frame(COI.scores %>% 
                               mutate(ARMS = substr(fraction, 1, 9)))
COI.scores1

# convex hulls

COI.convex.hull <- function(COI.scores1) COI.scores1[chull(COI.scores1$NMDS1, COI.scores1$NMDS2),]
COI.convex.hull(COI.scores1)

# ddply to get hulls based on species

COI.species.hulls <- plyr::ddply(COI.scores1, "ARMS", COI.convex.hull)
COI.species.hulls

COI.scores2 <- as.data.frame(COI.scores1[order(COI.scores1$ARMS),])
COI.spp <- as.data.frame(unique(COI.scores2$ARMS)) %>%
  rename(ARMS = "unique(COI.scores2$ARMS)")

# mds plot with color by species

COI.mdsplot <- ggplot(COI.scores2, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = COI.species.hulls, aes(x = NMDS1, y = NMDS2,
                                             fill = ARMS, group = ARMS), alpha = 0.4, lty = 1, lwd = 0.3, color = "black") +
  geom_point(size = 2, aes(color = ARMS, group = ARMS, alpha = 0.6)) +
  geom_text(data = COI.scores2, aes(x = NMDS1, y = NMDS2, label = fraction), size = 3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.title = element_text(size=14),
                     axis.text.x = element_text(size = 14),
                     axis.text.y = element_text(size = 14),
                     legend.position="none") +
  scale_color_manual(values = c(rep("coral1",3), rep("gold",3), rep("deepskyblue2",3))) + 
  scale_fill_manual(values=c(rep("coral1",3), rep("gold",3), rep("deepskyblue2",3)))
COI.mdsplot

# ggsave("plots/COI_MDS.pdf", COI.mdsplot, width = 10, height = 10, useDingbats = FALSE)


##### PHYLA SUMMARY #####

# summarize by phyla

COI.rra.phyla <- COI.nosingleton %>%
  select(-OTUID) %>%
  pivot_longer(names_to = "ARMS", values_to = "seq", -Phylum) %>%
  group_by(ARMS, Phylum) %>%
  summarize(seq.phylum = sum(seq)) %>%
  ungroup() %>%
  group_by(ARMS) %>%
  mutate(rra = seq.phylum/sum(seq.phylum)) %>%
  pivot_wider(id_cols = Phylum, names_from = ARMS, values_from = rra, values_fill = list(rra = 0))
head(COI.rra.phyla)

COI.rra.phyla %>% 
  tbl_df %>% 
  print(n=39)

# calculate means for each phlya

COI.rra.mean <- COI.rra.phyla %>% 
  mutate(mean_all = rowMeans(.[2:28]))

# fetch average unknown taxonomic assignments across samples - 41.1%

COI.rra.mean[38,"mean_all"]

# sum fractions from each ARMS

COI.phylum.helper <- COI.nosingleton %>%
  select(OTUID, Phylum)

COI.group <- COI.nosingleton %>%
  select(-Phylum) %>%
  pivot_longer(names_to = "COI", values_to = "value", -OTUID) %>%
  mutate(ARMS = substr(COI, 1, 9)) %>%
  group_by(ARMS, OTUID) %>%
  summarize(sum = sum (value)) %>%
  pivot_wider(id_cols = OTUID , names_from = ARMS, values_from = sum, values_fill = list(sum = 0)) %>%
  left_join(COI.phylum.helper)

head(COI.group)

# remove unknowns & summarize by phyla

COI.phyla.known <- COI.group %>%
  select(-OTUID) %>%
  pivot_longer(names_to = "ARMS", values_to = "seq", -Phylum) %>%
  filter(Phylum != "Unknown") %>%
  group_by(ARMS, Phylum) %>%
  summarize(seq.phylum = sum(seq)) %>%
  ungroup() %>%
  group_by(ARMS) %>%
  mutate(rra = seq.phylum/sum(seq.phylum)) %>%
  pivot_wider(id_cols = Phylum, names_from = ARMS, values_from = rra, values_fill = list(rra = 0))
head(COI.phyla.known)

COI.phyla.known %>% 
  tbl_df %>% 
  print(n=38)

# rank phlya by row means

COI.rank <- COI.phyla.known %>% 
  mutate(mean_all = rowMeans(.[2:10])) %>%
  arrange(-mean_all)

head(COI.rank)

COI.rank %>% 
  tbl_df %>% 
  pull(Phylum)

# get rid of top 10 phyla, combine all other phyla into "Other" category

COI.other <- COI.rank %>%
  filter(mean_all<0.01) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Phylum="Other") %>%
  select(Phylum,everything())

# add top 10 phyla back to "Other" category

COI.phyla.final <- bind_rows(COI.rank,COI.other) %>%
filter(mean_all>0.01)

# unite Site and Year from metadata

COI.meta <- bali_meta %>%
  unite(Treatment, c("Site","Year"), remove = FALSE) %>%
  select(ARMS, Treatment)

COI.join <- COI.phyla.final %>%
  select(-mean_all) %>%
  pivot_longer(names_to = "ARMS", values_to = "value", -Phylum) %>% 
  pivot_wider(id_cols = ARMS, names_from = Phylum, values_from = value, values_fill = list(value = 0)) %>%
  left_join(COI.meta) %>%
  distinct()

COI.phyla.sum <- COI.join %>%
  pivot_longer(names_to = "Phylum", values_to = "value", -c(ARMS, Treatment)) %>%
  mutate(percent = value*100) %>%
  group_by(Treatment, Phylum) %>%
  summarize(mean = mean(percent), sd = sd(percent), n = n(), se = sd/sqrt(n)) %>%
  arrange(-mean) %>%
  mutate(ordervar = case_when(Treatment == "S1_2012" ~ 1,
                              Treatment == "S2_2012" ~ 2,
                              TRUE ~ 3))

COI.phyla.sum 

# phyla to be plotted

COI.mean.rank <- COI.phyla.final %>%
  select(Phylum)

# caterpillar plot of percent composition of top phylum

COI.caterpillar <- ggplot(COI.phyla.sum, aes(y = reorder(Treatment, -ordervar), 
                                             x = mean, color = reorder(Treatment, -ordervar))) +
  geom_point() +
  geom_errorbarh(aes(xmin = mean-se, xmax = mean+se), height = 0) +
  facet_grid(reorder(Phylum, -mean)~., switch = "both") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text.y = element_text(color = "black", angle = 180),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "top", legend.title = element_blank()) +
  scale_colour_manual(values = c("deepskyblue", "gold", "coral1"), 
                      labels = c("2013_Site1", "2012_Site2", "2012_Site1"),
                      guide = guide_legend(reverse = TRUE)) +
  ylab("Phylum") +
  xlab("Relative abundance (mean ± SE)") + 
  ggtitle("(a) COI") +
  theme(plot.title = element_text(size = 14))
COI.caterpillar

# ggsave("plots/COI_Caterpillar.pdf", COI.caterpillar, width = 6, height = 8, useDingbats = FALSE)



########*#########
###### 18S #######
########*#########

### GET RID OF SINGLETONS

x18S.nosingleton <- bali_18S %>% 
  mutate(rowsum = rowSums(.[2:28])) %>%
  filter(rowsum > 1) %>%
  select(-rowsum)

# number of OTUs annotated with Phylum "Unknown" - 209

length(which(x18S.nosingleton$Phylum == "Unknown"))

# total number of sequence reads per ARMS fraction and overarching sum

x18S.total <- x18S.nosingleton %>%
  select(-c(OTUID, Phylum)) %>%
  summarize_all(funs(sum)) %>%
  mutate(sum = rowSums(.))


##### RAREFACTION CURVES #####

### RAREFACTION BY 500 FRACTION

# select 500 fraction samples

x18S.500 <- x18S.nosingleton %>%
  select(contains("F500"))

x18S.500[c(1:9),c(1:9)]

rownames(x18S.500) <- x18S.500$OTUID

x18S.list500 <- as.data.frame(t(x18S.500))
str(x18S.list500)

x18S.list500 <- as.data.frame(x18S.list500)

x18S.list500[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

x18S.list500 <- split(x18S.list500, rownames(x18S.list500))

x18S.list500$INDO18S12S1AF500 <- t(x18S.list500$INDO18S12S1AF500)
x18S.list500$INDO18S12S1BF500 <- t(x18S.list500$INDO18S12S1BF500)
x18S.list500$INDO18S12S1CF500 <- t(x18S.list500$INDO18S12S1CF500)

x18S.list500$INDO18S12S2AF500 <- t(x18S.list500$INDO18S12S2AF500)
x18S.list500$INDO18S12S2BF500 <- t(x18S.list500$INDO18S12S2BF500)
x18S.list500$INDO18S12S2CF500 <- t(x18S.list500$INDO18S12S2CF500)

x18S.list500$INDO18S13S1AF500 <- t(x18S.list500$INDO18S13S1AF500)
x18S.list500$INDO18S13S1BF500 <- t(x18S.list500$INDO18S13S1BF500)
x18S.list500$INDO18S13S1CF500 <- t(x18S.list500$INDO18S13S1CF500)

# rarefaction curves by samples

x18S.rf.500 <- iNEXT(x18S.list500, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
x18S.rf.500

x18S.rare.500 <- fortify(x18S.rf.500, type=1)
head(x18S.rare.500)

x18S.rare.pt.500 <-x18S.rare.500[which(x18S.rare.500$method=="observed"),]
x18S.rare.ln.500 <-x18S.rare.500[which(x18S.rare.500$method!="observed"),]
x18S.rare.ln.500$method <- factor(x18S.rare.ln.500$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
x18S.500.points <- x18S.rare.500 %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
x18S.500.rarefaction <- ggplot(x18S.rare.500, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=x18S.rare.ln.500) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = x18S.500.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO18S12S1AF500="coral1", INDO18S12S1BF500="coral1", INDO18S12S1CF500="coral1", 
                                 INDO18S12S2AF500="gold", INDO18S12S2BF500="gold", INDO18S12S2CF500="gold", 
                                 INDO18S13S1AF500="deepskyblue", INDO18S13S1BF500="deepskyblue", INDO18S13S1CF500="deepskyblue")) +
  scale_fill_manual(values = c(INDO18S12S1AF500="coral1", INDO18S12S1BF500="coral1", INDO18S12S1CF500="coral1", 
                               INDO18S12S2AF500="gold", INDO18S12S2BF500="gold", INDO18S12S2CF500="gold", 
                               INDO18S13S1AF500="deepskyblue", INDO18S13S1BF500="deepskyblue", INDO18S13S1CF500="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
x18S.500.rarefaction

# ggsave("plots/18S_Rare_500.pdf", x18S.500.rarefaction, width = 8, height = 6, useDingbats = FALSE)


### RAREFACTION BY 100 FRACTION

# select 100 fraction samples

x18S.100 <- x18S.nosingleton %>%
  select(contains("F100"))

x18S.100[c(1:9),c(1:9)]

rownames(x18S.100) <- x18S.100$OTUID

x18S.list100 <- as.data.frame(t(x18S.100))
str(x18S.list100)

x18S.list100 <- as.data.frame(x18S.list100)

x18S.list100[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

x18S.list100 <- split(x18S.list100, rownames(x18S.list100))

x18S.list100$INDO18S12S1AF100 <- t(x18S.list100$INDO18S12S1AF100)
x18S.list100$INDO18S12S1BF100 <- t(x18S.list100$INDO18S12S1BF100)
x18S.list100$INDO18S12S1CF100 <- t(x18S.list100$INDO18S12S1CF100)

x18S.list100$INDO18S12S2AF100 <- t(x18S.list100$INDO18S12S2AF100)
x18S.list100$INDO18S12S2BF100 <- t(x18S.list100$INDO18S12S2BF100)
x18S.list100$INDO18S12S2CF100 <- t(x18S.list100$INDO18S12S2CF100)

x18S.list100$INDO18S13S1AF100 <- t(x18S.list100$INDO18S13S1AF100)
x18S.list100$INDO18S13S1BF100 <- t(x18S.list100$INDO18S13S1BF100)
x18S.list100$INDO18S13S1CF100 <- t(x18S.list100$INDO18S13S1CF100)

# rarefaction curves by samples

x18S.rf.100 <- iNEXT(x18S.list100, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
x18S.rf.100

x18S.rare.100 <- fortify(x18S.rf.100, type=1)
head(x18S.rare.100)

x18S.rare.pt.100 <-x18S.rare.100[which(x18S.rare.100$method=="observed"),]
x18S.rare.ln.100 <-x18S.rare.100[which(x18S.rare.100$method!="observed"),]
x18S.rare.ln.100$method <- factor(x18S.rare.ln.100$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
x18S.100.points <- x18S.rare.100 %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
x18S.100.rarefaction <- ggplot(x18S.rare.100, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=x18S.rare.ln.100) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = x18S.100.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO18S12S1AF100="coral1", INDO18S12S1BF100="coral1", INDO18S12S1CF100="coral1", 
                                 INDO18S12S2AF100="gold", INDO18S12S2BF100="gold", INDO18S12S2CF100="gold", 
                                 INDO18S13S1AF100="deepskyblue", INDO18S13S1BF100="deepskyblue", INDO18S13S1CF100="deepskyblue")) +
  scale_fill_manual(values = c(INDO18S12S1AF100="coral1", INDO18S12S1BF100="coral1", INDO18S12S1CF100="coral1", 
                               INDO18S12S2AF100="gold", INDO18S12S2BF100="gold", INDO18S12S2CF100="gold", 
                               INDO18S13S1AF100="deepskyblue", INDO18S13S1BF100="deepskyblue", INDO18S13S1CF100="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
x18S.100.rarefaction

# ggsave("plots/18S_Rare_100.pdf", x18S.100.rarefaction, width = 8, height = 6, useDingbats = FALSE)


### RAREFACTION BY SESSILE FRACTION

# select sessile fraction samples

x18S.SES <- x18S.nosingleton %>%
  select(contains("SES"))

x18S.SES[c(1:9),c(1:9)]

rownames(x18S.SES) <- x18S.SES$OTUID

x18S.listSES <- as.data.frame(t(x18S.SES))
str(x18S.listSES)

x18S.listSES <- as.data.frame(x18S.listSES)

x18S.listSES[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

x18S.listSES <- split(x18S.listSES, rownames(x18S.listSES))

x18S.listSES$INDO18S12S1ASES <- t(x18S.listSES$INDO18S12S1ASES)
x18S.listSES$INDO18S12S1BSES <- t(x18S.listSES$INDO18S12S1BSES)
x18S.listSES$INDO18S12S1CSES <- t(x18S.listSES$INDO18S12S1CSES)

x18S.listSES$INDO18S12S2ASES <- t(x18S.listSES$INDO18S12S2ASES)
x18S.listSES$INDO18S12S2BSES <- t(x18S.listSES$INDO18S12S2BSES)
x18S.listSES$INDO18S12S2CSES <- t(x18S.listSES$INDO18S12S2CSES)

x18S.listSES$INDO18S13S1ASES <- t(x18S.listSES$INDO18S13S1ASES)
x18S.listSES$INDO18S13S1BSES <- t(x18S.listSES$INDO18S13S1BSES)
x18S.listSES$INDO18S13S1CSES <- t(x18S.listSES$INDO18S13S1CSES)

# rarefaction curves by samples

x18S.rf.SES <- iNEXT(x18S.listSES, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
x18S.rf.SES

x18S.rare.SES <- fortify(x18S.rf.SES, type=1)
head(x18S.rare.SES)

x18S.rare.pt.SES <-x18S.rare.SES[which(x18S.rare.SES$method=="observed"),]
x18S.rare.ln.SES <-x18S.rare.SES[which(x18S.rare.SES$method!="observed"),]
x18S.rare.ln.SES$method <- factor(x18S.rare.ln.SES$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
x18S.SES.points <- x18S.rare.SES %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
x18S.SES.rarefaction <- ggplot(x18S.rare.SES, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=x18S.rare.ln.SES) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = x18S.SES.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO18S12S1ASES="coral1", INDO18S12S1BSES="coral1", INDO18S12S1CSES="coral1", 
                                 INDO18S12S2ASES="gold", INDO18S12S2BSES="gold", INDO18S12S2CSES="gold", 
                                 INDO18S13S1ASES="deepskyblue", INDO18S13S1BSES="deepskyblue", INDO18S13S1CSES="deepskyblue")) +
  scale_fill_manual(values = c(INDO18S12S1ASES="coral1", INDO18S12S1BSES="coral1", INDO18S12S1CSES="coral1", 
                               INDO18S12S2ASES="gold", INDO18S12S2BSES="gold", INDO18S12S2CSES="gold", 
                               INDO18S13S1ASES="deepskyblue", INDO18S13S1BSES="deepskyblue", INDO18S13S1CSES="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
x18S.SES.rarefaction

# ggsave("plots/18S_Rare_SES.pdf", x18S.SES.rarefaction, width = 8, height = 6, useDingbats = FALSE)


### RAREFACTION BY ARMS

x18S.nosingleton[c(1:20),c(1:20)]

# merge x18S data with metadata

x18S.long <- x18S.nosingleton %>%
  select(-Phylum) %>%
  pivot_longer(names_to = "x18S", values_to = "value", -OTUID) %>%
  inner_join(bali_meta)

# merge fractions

x18S.merged <- x18S.long %>% group_by(Year, Site, Unit, OTUID) %>% summarize(summed.value = sum(value)) %>%
  mutate(ARMSID = paste("INDO18S",Year,Site,Unit, sep = "")) %>%
  ungroup() %>%
  select(-Year, -Site, -Unit) %>%
  pivot_wider(names_from = OTUID, values_from = summed.value)

# total sequence reads for each ARMS

x18S.count <- x18S.merged %>%
  select(-ARMSID) %>%
  mutate(sums = rowSums(.))
x18S.count$sums

# total OTUs for each ARMS

x18S.otu <- x18S.long %>% group_by(Year, Site, Unit, OTUID) %>% summarize(summed.value = sum(value)) %>%
  mutate(ARMSID = paste("INDO18S",Year,Site,Unit, sep = "")) %>%
  ungroup() %>%
  select(-Year, -Site, -Unit) %>%
  mutate(pa.otu = case_when(summed.value > 0 ~ 1,
                   TRUE ~ 0)) %>%
  group_by(ARMSID) %>%
  summarize(otu.count = sum(pa.otu))

# prepare data for rarefaction

x18S.list <- as.data.frame(x18S.merged)
rownames(x18S.list) <- x18S.list$ARMSID
x18S.listARMS <- x18S.list[-1]

x18S.listARMS[c(1:9), c(1:9)]

# split dataframe into list, with each ARMS as an element

x18S.listARMS <- split(x18S.listARMS, rownames(x18S.listARMS))

x18S.listARMS$INDO18S2012S1A <- t(x18S.listARMS$INDO18S2012S1A)
x18S.listARMS$INDO18S2012S1B <- t(x18S.listARMS$INDO18S2012S1B)
x18S.listARMS$INDO18S2012S1C <- t(x18S.listARMS$INDO18S2012S1C)

x18S.listARMS$INDO18S2012S2A <- t(x18S.listARMS$INDO18S2012S2A)
x18S.listARMS$INDO18S2012S2B <- t(x18S.listARMS$INDO18S2012S2B)
x18S.listARMS$INDO18S2012S2C <- t(x18S.listARMS$INDO18S2012S2C)

x18S.listARMS$INDO18S2013S1A <- t(x18S.listARMS$INDO18S2013S1A)
x18S.listARMS$INDO18S2013S1B <- t(x18S.listARMS$INDO18S2013S1B)
x18S.listARMS$INDO18S2013S1C <- t(x18S.listARMS$INDO18S2013S1C)

# rarefaction curves by samples

x18S.rf.ARMS <- iNEXT(x18S.listARMS, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=20, se=TRUE, conf=0.95, nboot=50)
x18S.rf.ARMS

x18S.rare.ARMS <- fortify(x18S.rf.ARMS, type=1)
head(x18S.rare.ARMS)

x18S.rare.pt.ARMS <-x18S.rare.ARMS[which(x18S.rare.ARMS$method=="observed"),]
x18S.rare.ln.ARMS <-x18S.rare.ARMS[which(x18S.rare.ARMS$method!="observed"),]
x18S.rare.ln.ARMS$method <- factor(x18S.rare.ln.ARMS$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
x18S.ARMS.points <- x18S.rare.ARMS %>% filter(method == "observed") %>% add_column(shape.vector = rep(1:3,3))
x18S.ARMS.rarefaction <- ggplot(x18S.rare.ARMS, aes(x=x, y=y, color=site)) +
  geom_line(aes(linetype=method), lwd=0.5, data=x18S.rare.ln.ARMS) + 
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, color=NULL), alpha=0.2) + 
  geom_point(data = x18S.ARMS.points, aes(x = x, y = y, color = site, shape = as.factor(shape.vector)), size = 2) + 
  labs(x="Sequences", y="Number of OTUs") +
  scale_colour_manual(values = c(INDO18S2012S1A="coral1", INDO18S2012S1B="coral1", INDO18S2012S1C="coral1", 
                                 INDO18S2012S2A="gold", INDO18S2012S2B="gold", INDO18S2012S2C="gold", 
                                 INDO18S2013S1A="deepskyblue", INDO18S2013S1B="deepskyblue", INDO18S2013S1C="deepskyblue")) +
  scale_fill_manual(values = c(INDO18S2012S1A="coral1", INDO18S2012S1B="coral1", INDO18S2012S1C="coral1", 
                               INDO18S2012S2A="gold", INDO18S2012S2B="gold", INDO18S2012S2C="gold", 
                               INDO18S2013S1A="deepskyblue", INDO18S2013S1B="deepskyblue", INDO18S2013S1C="deepskyblue")) +
  scale_shape_manual(values = 15:17, name = "ARMS", labels = c("ARMSA", "ARMSB", "ARMSC")) +
  theme (legend.position = "right", legend.title = element_blank(), text=element_text(size=18)) + theme_classic()
x18S.ARMS.rarefaction

# ggsave("plots/18S_Rare_ARMS.pdf", x18S.ARMS.rarefaction, width = 8, height = 6, useDingbats = FALSE)


##### RELATIVE ABUNDANCES ##### 

head(x18S.nosingleton)

x18S.rra <- x18S.nosingleton %>%
  select(-Phylum) %>%
  pivot_longer(names_to = "ARMS", values_to = "seq", -OTUID) %>%
  group_by(ARMS) %>%
  mutate(rra = seq/sum(seq)) %>%
  pivot_wider(id_cols = ARMS, names_from = OTUID, values_from = rra, values_fill = list(rra = 0))

head(x18S.rra)


##### MDS #####

# run nmds ordination

x18S.nmds <- metaMDS(x18S.rra[-1], k = 3, metric = "bray", trymax = 1000)
stressplot(x18S.nmds)
x18S.nmds

# plot ordination with ggplot

x18S.scores <- as.data.frame(scores(x18S.nmds))
x18S.scores$fraction <- x18S.rra$ARMS
x18S.scores1 <- as.data.frame(x18S.scores %>% 
                               mutate(ARMS = substr(fraction, 1, 12)))
x18S.scores1

# convex hulls

x18S.convex.hull <- function(x18S.scores1) x18S.scores1[chull(x18S.scores1$NMDS1, x18S.scores1$NMDS2),]
x18S.convex.hull(x18S.scores1)

# ddply to get hulls based on species

x18S.species.hulls <- plyr::ddply(x18S.scores1, "ARMS", x18S.convex.hull)
x18S.species.hulls

x18S.scores2 <- as.data.frame(x18S.scores1[order(x18S.scores1$ARMS),])
x18S.spp <- as.data.frame(unique(x18S.scores2$ARMS)) %>%
  rename(ARMS = "unique(x18S.scores2$ARMS)")

# mds plot with color by species

x18S.mdsplot <- ggplot(x18S.scores2, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = x18S.species.hulls, aes(x = NMDS1, y = NMDS2,
                                             fill = ARMS, group = ARMS), alpha = 0.4, lty = 1, lwd = 0.3, color = "black") +
  geom_point(size = 2, aes(color = ARMS, group = ARMS, alpha = 0.6)) +
  geom_text(data = x18S.scores2, aes(x = NMDS1, y = NMDS2, label = fraction), size = 3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.title = element_text(size=14),
                     axis.text.x = element_text(size = 14),
                     axis.text.y = element_text(size = 14),
                     legend.position="none") +
  scale_color_manual(values = c(rep("coral1",3), rep("gold",3), rep("deepskyblue2",3))) + 
  scale_fill_manual(values=c(rep("coral1",3), rep("gold",3), rep("deepskyblue2",3)))
x18S.mdsplot

# ggsave("plots/18S_MDS.pdf", x18S.mdsplot, width = 10, height = 10, useDingbats = FALSE)


##### PHYLA SUMMARY #####

# summarize by phyla

x18S.rra.phyla <- x18S.nosingleton %>%
  select(-OTUID) %>%
  pivot_longer(names_to = "ARMS", values_to = "seq", -Phylum) %>%
  group_by(ARMS, Phylum) %>%
  summarize(seq.phylum = sum(seq)) %>%
  ungroup() %>%
  group_by(ARMS) %>%
  mutate(rra = seq.phylum/sum(seq.phylum)) %>%
  pivot_wider(id_cols = Phylum, names_from = ARMS, values_from = rra, values_fill = list(rra = 0))
head(x18S.rra.phyla)

x18S.rra.phyla %>% 
  tbl_df %>% 
  print(52)

# calculate means for each phlya

x18S.rra.mean <- x18S.rra.phyla %>% 
  mutate(mean_all = rowMeans(.[2:28]))

# fetch average unknown taxonomic assignments across samples - 0.62%

x18S.rra.mean[50,"mean_all"]

# sum fractions from each ARMS

x18S.phylum.helper <- x18S.nosingleton %>%
  select(OTUID, Phylum)

x18S.group <- x18S.nosingleton %>%
  select(-Phylum) %>%
  pivot_longer(names_to = "x18S", values_to = "value", -OTUID) %>%
  mutate(ARMS = substr(x18S, 1, 12)) %>%
  group_by(ARMS, OTUID) %>%
  summarize(sum = sum (value)) %>%
  pivot_wider(id_cols = OTUID , names_from = ARMS, values_from = sum, values_fill = list(sum = 0)) %>%
  left_join(x18S.phylum.helper)

head(x18S.group)

# remove unknowns & summarize by phyla

x18S.phyla.known <- x18S.group %>%
  select(-OTUID) %>%
  pivot_longer(names_to = "ARMS", values_to = "seq", -Phylum) %>%
  filter(Phylum != "Unknown") %>%
  group_by(ARMS, Phylum) %>%
  summarize(seq.phylum = sum(seq)) %>%
  ungroup() %>%
  group_by(ARMS) %>%
  mutate(rra = seq.phylum/sum(seq.phylum)) %>%
  pivot_wider(id_cols = Phylum, names_from = ARMS, values_from = rra, values_fill = list(rra = 0))
head(x18S.phyla.known)

x18S.phyla.known %>% 
  tbl_df %>% 
  print(n=38)

# rank phlya by row means

x18S.rank <- x18S.phyla.known %>% 
  mutate(mean_all = rowMeans(.[2:10])) %>%
  arrange(-mean_all)

head(x18S.rank)

x18S.rank %>% 
  tbl_df %>% 
  pull(Phylum)

# get rid of top 10 phyla, combine all other phyla into "Other" category

x18S.other <- x18S.rank %>%
  filter(mean_all<0.01471) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Phylum="Other") %>%
  select(Phylum,everything())

# add top 10 phyla back to "Other" category

x18S.phyla.final <- bind_rows(x18S.rank,x18S.other) %>%
  filter(mean_all>0.01471)

# unite Site and Year from metadata

x18S.meta <- bali_meta %>%
  unite(Treatment, c("Site","Year"), remove = FALSE) %>%
  select(ARMS.18S, Treatment)

x18S.join <- x18S.phyla.final %>%
  select(-mean_all) %>%
  pivot_longer(names_to = "ARMS.18S", values_to = "value", -Phylum) %>% 
  pivot_wider(id_cols = ARMS.18S, names_from = Phylum, values_from = value, values_fill = list(value = 0)) %>%
  left_join(x18S.meta) %>%
  distinct()

x18S.phyla.sum <- x18S.join %>%
  pivot_longer(names_to = "Phylum", values_to = "value", -c(ARMS.18S, Treatment)) %>%
  mutate(percent = value*100) %>%
  group_by(Treatment, Phylum) %>%
  summarize(mean = mean(percent), sd = sd(percent), n = n(), se = sd/sqrt(n)) %>%
  arrange(-mean) %>%
  mutate(ordervar = case_when(Treatment == "S1_2012" ~ 1,
                              Treatment == "S2_2012" ~ 2,
                              TRUE ~ 3))

x18S.phyla.sum 

# phyla to be plotted

x18S.mean.rank <- x18S.phyla.final %>%
  select(Phylum)

# caterpillar plot of percent composition of top phylum

x18S.caterpillar <- ggplot(x18S.phyla.sum, aes(y = reorder(Treatment, -ordervar), 
                                               x = mean, color = reorder(Treatment, -ordervar))) +
  geom_point() +
  geom_errorbarh(aes(xmin = mean-se, xmax = mean+se), height = 0) +
  facet_grid(reorder(Phylum, -mean)~., switch = "both") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text.y = element_text(color = "black", angle = 180),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "top", legend.title = element_blank()) +
  scale_colour_manual(values = c("deepskyblue", "gold", "coral1"), 
                      labels = c("2013_Site1", "2012_Site2", "2012_Site1"),
                      guide = guide_legend(reverse = TRUE)) +
  ylab("Phylum") +
  xlab("Relative abundance (mean ± SE)") +
  ggtitle("(b) 18S") +
  theme(plot.title = element_text(size = 14))
x18S.caterpillar

# ggsave("plots/18S_Caterpillar.pdf", x18S.caterpillar, width = 6, height = 8, useDingbats = FALSE)



##############*##############
###### COI + 18S PLOTS ######
##############*##############

# rarefaction curves by fraction

plot_grid(COI.500.rarefaction, x18S.500.rarefaction, 
          COI.100.rarefaction, x18S.100.rarefaction, 
          COI.SES.rarefaction, x18S.SES.rarefaction, 
          align="hv", nrow=3, ncol=2)

# ggsave("plots/COI18S_Rare_Fractions.pdf", plot_grid(COI.500.rarefaction, x18S.500.rarefaction, 
#                                              COI.100.rarefaction, x18S.100.rarefaction, 
#                                              COI.SES.rarefaction, x18S.SES.rarefaction,
#                                              align="hv", nrow=3, ncol=2), width = 16, height = 20, useDingbats = FALSE)

# rarefaction curves by ARMS

plot_grid(COI.ARMS.rarefaction, x18S.ARMS.rarefaction, align="h")

# ggsave("plots/COI18S_Rare_ALL.pdf", plot_grid(COI.ARMS.rarefaction, x18S.ARMS.rarefaction, 
#                                        align="h"), width = 16, height = 6, useDingbats = FALSE)

# MDS plots

plot_grid(COI.mdsplot, x18S.mdsplot, align="h")

# ggsave("plots/COI18S_MDS.pdf", plot_grid(COI.mdsplot, x18S.mdsplot, 
#                                   align="h"), width = 15, height = 7, useDingbats = FALSE)

# caterpillar plots

plot_grid(COI.caterpillar, x18S.caterpillar, align="h")

#ggsave("plots/COI18S_Caterpillar.pdf", plot_grid(COI.caterpillar, x18S.caterpillar, 
#                                           align="h"), width = 12, height = 8, useDingbats = FALSE)

patch <- COI.caterpillar + x18S.caterpillar
patch



#############*###############
###### COI v 18S MODEL ######
#############*###############

### data wrangling

# create metadata for model

meta.siteyear <- bali_meta %>%
  unite(Treatment, c("Site","Year"), remove = FALSE) %>%
  mutate(Treatment = paste(Year, Site, sep = "_")) %>%
  select(ARMS, ARMS.18S, Treatment, Site, Year) %>%
  unique()

# COI data to long format

COI.long <- COI.phyla.known %>%
  filter_at(vars(-Phylum), all_vars(.>0)) %>%
  pivot_longer(names_to = "ARMS", values_to = "rra", -Phylum) %>%
  left_join(meta.siteyear) %>%
  mutate(marker = "COI")

# 18S data to long format

x18S.long <- x18S.phyla.known %>%
  filter_at(vars(-Phylum), all_vars(.>0)) %>%
  pivot_longer(names_to = "ARMS.18S", values_to = "rra", -Phylum) %>%
  left_join(meta.siteyear) %>%
  mutate(marker = "x18S")

# phyla that co-occur in COI and 18S datasets

match.phyla <- unique(COI.long$Phylum)[unique(COI.long$Phylum) %in% unique(x18S.long$Phylum)]
length(match.phyla)

# combine COI and 18S datasets

combined.bali <- full_join(COI.long, x18S.long) %>%
  filter(Phylum %in% match.phyla) %>%
  ungroup()

# plot summary of combined datasets

combined.plot <- ggplot(combined.bali)+
  geom_boxplot(aes(x = marker, y = log(rra), fill = Treatment)) +
  facet_wrap(~Phylum, scales = "free") +
  theme_bw()


##### RUN MODEL #####

# primer and phylum as fixed effects, random effect of treatment (year and site)

fit <-  brms::brm(log(rra) ~ 0 + Phylum:marker + (1|Phylum:Treatment), 
            data = combined.bali, family = "student",
            control = list(adapt_delta = 0.9), 
            prior = brms::set_prior("uniform(-100,0)", lb = -100, ub = 0, class = "b"))
summary(fit)

plot(fit) # trace plots look good
brms::pp_check(fit) # posterior prediction plot looks good

# dataframe for predictions

prr <- unique(select(combined.bali, Phylum, Year, Site, Treatment, marker))

# arrange model to create prediction plot

pred <- fitted(fit, prr, summary = FALSE)
pred <- t(pred) 
pred <- cbind(prr, pred)
pred <- tidyr::gather(pred, key = "key", value = "value", -Phylum, -Year, -Site, -marker, -Treatment) 


##### MODEL PLOTS #####

### prediction plots of 20 co-occuring phyla

COI.18S.plot <- ggplot(pred) + 
  geom_density_ridges(aes(x = exp(value), fill = marker, y = Treatment), alpha = 0.7, 
                       quantile_lines = TRUE, quantiles = c(0.025, 0.975)) +
  scale_fill_manual(values = c("royalblue3", "aquamarine2"),
                    labels = c("COI", "18S"), name = "Marker") +
  theme_bw() +
  labs(x = "Relative abundance", y = "") +
  scale_y_discrete(limits = rev(levels(as.factor(pred$Treatment))), labels = c("2013_Site1", "2012_Site2", "2012_Site1")) +
  facet_wrap(~Phylum, scales = "free", ncol = 3) +
  theme(strip.background =element_rect(fill = "white"),
        panel.grid = element_blank(),
        legend.position = c(0.83, 0.06))

ggsave("plots/COIvs18S_Model.png", COI.18S.plot, width = 10, height = 14)


### random effect plots

# concatenate phylum, site, and marker

random.effects <- fit %>% spread_draws(`r_Phylum:Treatment`[level,])  %>% 
  separate(level, into = c("Phylum", "Year", "Site"), sep = "_", remove = FALSE) %>%
  mutate(Treatment = paste(Year, Site, sep = "_"))
  

# random effect plot

random.plot <- random.effects %>%
  ggplot(aes(y = Phylum, x = `r_Phylum:Treatment`, fill = Treatment), alpha = 0.5) +
  geom_density_ridges(alpha = 0.5, quantile_lines = TRUE, quantiles = c(0.025, 0.975)) +
  scale_fill_manual(values = c( "coral1", "gold", "deepskyblue"), 
                      labels = c("2012_Site1", "2012_Site2", "2013_Site1"),
                    name = "") +  
  theme_bw() +
  scale_y_discrete(limits = rev(levels(as.factor(pred$Phylum)))) +
  labs(y = "", x = "random effect") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0), 
        legend.position = "top")

 ggsave("plots/COIv18S_RandomEffects.png", random.plot, width = 5, height = 8)


##### MODEL SUMMARY #####

# prepare data frame

prr.sum <- unique(select(combined.bali, Phylum, marker)) %>% 
  as.data.frame()

# predict only to phyla and marker, without random effects

phyla.sum <- 
  cbind(prr.sum,
        as.data.frame(fitted(fit, newdata = prr.sum, 
                             re_formula = "log(rra) ~ 0 + Phylum:marker"))) %>%
  mutate(rra.median = exp(Estimate)) %>%
  mutate(Q2.5.bt = exp(Q2.5)) %>%
  mutate(Q97.5.bt = exp(Q97.5))



#############*##############
###### CORALNET MODEL ######
#############*##############

# subset COI sessile data

COI.SES <- COI.rra.phyla %>%
  pivot_longer(names_to = "COI", values_to = "rra", -Phylum) %>%
  left_join(bali_meta) %>%
  filter(Fraction == "SES") %>%
  mutate(marker = "COI")

# determine distinct number of identified phyla in sessile fraction - 32

COI.SES.distinct <- COI.rra.phyla %>%
  pivot_longer(names_to = "COI", values_to = "rra", -Phylum) %>%
  left_join(bali_meta) %>%
  filter(Fraction == "SES") %>%
  pivot_wider(id_cols = Phylum, names_from = ARMS, values_from = rra, values_fill = list(rra = 0)) %>%
  mutate(rowsum = rowSums(.[2:10])) %>%
  filter(rowsum > 0) %>%
  select(-rowsum)

# subset 18S sessile data

x18S.SES <- x18S.rra.phyla %>%
  pivot_longer(names_to = "x18S", values_to = "rra", -Phylum) %>%
  left_join(bali_meta) %>%
  filter(Fraction == "SES")%>%
  mutate(marker = "x18S")

# determine distinct number of identified phyla in sessile fraction - 44

x18S.SES.distinct <- x18S.rra.phyla %>%
  pivot_longer(names_to = "x18S", values_to = "rra", -Phylum) %>%
  left_join(bali_meta) %>%
  filter(Fraction == "SES") %>%
  pivot_wider(id_cols = Phylum, names_from = ARMS, values_from = rra, values_fill = list(rra = 0)) %>%
  mutate(rowsum = rowSums(.[2:10])) %>%
  filter(rowsum > 0) %>%
  select(-rowsum)

# COI sessile data to long format

COI.long.ses <- COI.SES.distinct %>%
  filter_at(vars(-Phylum), all_vars(.>0)) %>%
  pivot_longer(names_to = "ARMS", values_to = "rra", -Phylum) %>%
  left_join(meta.siteyear) %>%
  mutate(marker = "COI")

# 18S sessile data to long format

x18S.long.ses <- x18S.SES.distinct %>%
  filter_at(vars(-Phylum), all_vars(.>0)) %>%
  pivot_longer(names_to = "ARMS.18S", values_to = "rra", -Phylum) %>%
  left_join(meta.siteyear) %>%
  mutate(marker = "x18S")

# phyla that co-occur in sessile COI and 18S datasets -14

match.phyla.ses <- unique(COI.long.ses$Phylum)[unique(COI.long.ses$Phylum) %in% unique(x18S.long.ses$Phylum)]
length(match.phyla.ses)


# combine COI and 18S sessile datasets

combined.SES <- full_join(COI.SES, x18S.SES) %>%
  mutate(ARMS = str_remove(ARMS, "COI")) %>%
  select(Phylum, ARMS, Year, Site, marker, rra)

# format and subset coralnet dataset (remove unknown, unavailable, and none)
  
coralnet.subset <- coralnet %>%
  select(-unavailable, -unknown, -none) %>%
  pivot_longer(names_to = "Phylum", values_to = "n", 3:11) %>%
  group_by(ARMS, Phylum) %>%
  summarise(n = sum(n)) %>%
  group_by(ARMS) %>%
  mutate(ntot = sum(n)) %>%
  ungroup() %>%
  mutate(rra = n/ntot) %>%
  mutate(Phylum = str_to_title(Phylum)) %>%
  select(ARMS, Phylum, rra)

# phyla to remove - any phyla with zero occurrences for any sammple - only Cnidaria

phyla.remove <- coralnet.subset[coralnet.subset$rra==0, "Phylum"] %>%
  unique()

# remove Cnidaria from coralnet dataset

coralnet.subset <- coralnet.subset %>%
  filter(!Phylum %in% phyla.remove) %>%
  mutate(marker = "coralnet") 

# combine COI, 18S, and coralnet datasets

combined.ses.cn <- coralnet.subset %>%
  left_join(unique(select(combined.SES, ARMS, Site, Year))) %>%
  full_join(combined.SES) %>%
  # only keep Phyla that are present in all datasets
  filter(Phylum %in% coralnet.subset$Phylum,
         Phylum %in% combined.SES$Phylum) %>%
  # rescale rra
  group_by(ARMS, Year, Site, marker) %>%
  mutate(tot = sum(rra)) %>%
  mutate(rra = rra/tot) %>%
  select(-tot)  %>%
  mutate(treatment = paste(Site, Year)) %>%
  ungroup()

# plot summary of combined datasets

combined.plot.cn <- ggplot(combined.ses.cn) +
  geom_boxplot((aes(x = treatment, y = log(rra), color = marker))) +
  facet_wrap(~Phylum)


###### RUN MODEL ######

# primer and phylum as fixed effects, random effect of treatment (year and site)

fit.cn <- brm(log(rra) ~ 0 + Phylum:marker + (1|Phylum:treatment), 
              data = combined.ses.cn, family = "student",
              control = list(adapt_delta = 0.9), 
              prior = set_prior("uniform(-100,0)", lb = -100, ub = 0, class = "b"))

summary(fit.cn)
plot(fit.cn) # trace plots look good
brms::pp_check(fit.cn) # posterior prediction plot looks good too

# dataframe for predictions

prr.cn <- unique(select(combined.ses.cn, Phylum, Year, Site, marker, treatment)) %>% 
  as.data.frame()

# arrange model to create prediction plot

pred.cn <- fitted(fit.cn, newdata = prr.cn, summary = FALSE)
pred.cn <- t(pred.cn) %>% as.data.frame()
pred.cn <- cbind(prr.cn, pred.cn)
pred.cn <- tidyr::gather(pred.cn, key = "key", value = "value", 
                         -Phylum, -marker, -Site, -Year, -treatment) 


##### MODEL PLOTS #####

# reorder factors for figure legend

pred.cn$marker <- factor(pred.cn$marker, levels=c("coralnet", "COI", "x18S"))

# reorder treatment factors - year and site - for y-axis

pred.cn <- pred.cn %>%
  mutate(treat.order = recode(treatment, 
                              "S1 2012" = "2012_Site1", 
                              "S2 2012" = "2012_Site2", 
                              "S1 2013" = "2013_Site1"))

### prediction plots of 7 co-occuring phyla

COI.18S.cn.plot <- ggplot(pred.cn) + 
  geom_density_ridges(aes(x = exp(value), fill = marker,  y = treat.order), 
                       quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
  scale_fill_manual(values = c("grey80", "royalblue3", "aquamarine2"),
                    labels = c("CoralNet", "COI", "18S"),
                    name = "Technique") +
  labs(x = "Relative abundance", y = "") +
  scale_y_discrete(limits = rev(levels(as.factor(pred.cn$treat.order)))) +
  facet_wrap( ~ Phylum, scales = "free", ncol = 2, strip.position = "top") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        legend.position = c(0.75, 0.12))

 ggsave("plots/COI_18S_CoralNet.png", COI.18S.cn.plot, width = 8, height = 10)


##### MODEL SUMMARY #####

# prepare data frame

prr.sum.cn <- unique(select(combined.ses.cn, Phylum, marker)) %>% 
  as.data.frame()

# predict only to phyla and marker, without random effects

phyla.sum.cn <- 
  cbind(prr.sum.cn,
        as.data.frame(fitted(fit.cn, newdata = prr.sum.cn, 
               re_formula = "log(rra) ~ 0 + Phylum:marker"))) %>%
  mutate(rra.median = exp(Estimate)) %>%
  mutate(Q2.5.bt = exp(Q2.5)) %>%
  mutate(Q97.5.bt = exp(Q97.5))
  
