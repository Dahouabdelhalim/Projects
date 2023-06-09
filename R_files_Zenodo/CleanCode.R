### Microbiomes associated with avian malaria survival differ between susceptible Hawaiian honeycreepers and sympatric malaria-resistant introduced birds ####
## Statistical analysis and data visualization 
## Prepared by Amanda K. Navine
## Based on code by Elin Videvall and Kristina L. Paxton
## Last edited 2022-09-22

## Install specific packages with bioconductor installer
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'decontam', 'phyloseq'))

## Load libraries
library('ggplot2')
library('dplyr')
library('vegan')
library('DESeq2')
library('decontam')
library('phyloseq')
library('MASS')
library('scales')
library('ggpubr')
library('tidyr')

## Define colors
yr.col = c('#9C179EFF','#ED7953FF') # From option 'C'
reg.col = c('#310A5DFF', '#801E6CFF', '#CA404AFF', '#F78410FF', '#FBB71CFF', '#FCFFA4FF') # From option 'B'
inf.col = c('#F8E025FF', '#0D0887FF') # From option 'C'
sp.col = c('#38286FFF', '#62FC6BFF') # From option 'H'
age.col = c('#440154FF', '#277F8EFF', '#B2DD2DFF') # From option 'D'

## Load processed data
load('ps.clean.RData')

### Processing ##############################################################################
#### Import to phyloseq ####
## Read in ASV table
otu = read.table(file = 'otu_table.tsv', header = TRUE, sep = '\\t', fill=T)
rownames(otu) = otu$OTU
otu$OTU = NULL # Changes ASV names to row names instead of giving them a separate column
otu = as.matrix(otu)

## Read in EDITED taxonomy (separated by kingdom, phylum, class, order, family, genus, species)
## Note: make sure this is fixed manually prior to R
taxonomy = read.table('edited-tax-Run1.tsv', sep = '\\t', row.names = 1, header=T)
taxonomy = as.matrix(taxonomy)

## Read in metadata
metadata = read.table('Metadata.tsv', row.names = 1, header=T, sep = '\\t')
metadata = metadata[!(row.names(metadata) %in% 'NEGCAMP09'),] # Removed due to low read count

## Make sure metadata columns are formatted correctly
library('plyr')
metadata$Elevation = as.character(metadata$Elevation)
metadata$Elevation = as.numeric(metadata$Elevation)
metadata$RelInf = as.character(metadata$RelInf)
metadata$RelInf = as.numeric(metadata$RelInf)
metadata$Date = as.Date(metadata$Date, '%m-%d-%Y')
metadata[,'Year'] = as.factor(format(metadata[,'Date'], '%Y'))
metadata[,'Month'] = as.factor(format(metadata[,'Date'], '%m'))
metadata$Host = revalue(metadata$Host, c('Amakihi' = "`Amakihi"))
metadata$Region = revalue(metadata$Region, c('Kau' = "Ka`\\u016B", 'Puu Waawaa' = "Pu`u Wa`awa`a", 'Kilauea' = 'K\\u012Blauea'))

## Read in tree
phy_tree = read_tree('tree.nwk')

## Import all as phyloseq objects
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

## Merge all four objects to one phyloseq object
ps = phyloseq(OTU, TAX, META, phy_tree)

## Return some metadata sampleID to a column for easier use later in Adonis
sample_data(ps)$SampleID = rownames(sample_data(ps))

#### Decontamination ####
## Use Decontam to remove contaminants based on negative controls (threshold = 0.1, default)
sample_data(ps)$is.neg = sample_data(ps)$SampleorControl == 'C'
contamdf.prev = isContaminant(ps, method='prevalence', neg='is.neg')

## Make phyloseq object of presence-absence in negative controls and true samples
ps.pa = transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg = prune_samples(sample_data(ps.pa)$SampleorControl == 'C', ps.pa)
ps.pa.pos = prune_samples(sample_data(ps.pa)$SampleorControl == 'S', ps.pa)
ps.nocontam = prune_taxa(!contamdf.prev$contaminant, ps)

## Remove blank samples + ASVs with 0 remaining reads
ps.sample = prune_samples(sample_data(ps.nocontam)$SampleorControl == 'S', ps.nocontam)
ps.sample = prune_taxa(taxa_sums(ps.sample) > 0, ps.sample)

## Remove mitochondria, chloroplast, eukaryote, and unassigned reads
mitochondria = subset_taxa(ps.sample, Family == 'Mitochondria')
chloroplast = subset_taxa(ps.sample, Order == 'Chloroplast')
eukaryote = subset_taxa(ps.sample, Kingdom == 'Eukaryota')
unassigned = subset_taxa(ps.sample, Kingdom == 'Unassigned')

badTaxa = row.names(tax_table(mitochondria)[,1])
allTaxa = taxa_names(ps.sample)
goodTaxa = allTaxa[!(allTaxa %in% badTaxa)]
ps.sample = prune_taxa(goodTaxa, ps.sample)

badTaxa = row.names(tax_table(chloroplast)[,1])
allTaxa = taxa_names(ps.sample)
goodTaxa = allTaxa[!(allTaxa %in% badTaxa)]
ps.sample = prune_taxa(goodTaxa, ps.sample)

badTaxa = row.names(tax_table(eukaryote)[,1])
allTaxa = taxa_names(ps.sample)
goodTaxa = allTaxa[!(allTaxa %in% badTaxa)]
ps.sample = prune_taxa(goodTaxa, ps.sample)

badTaxa = row.names(tax_table(unassigned)[,1])
allTaxa = taxa_names(ps.sample)
goodTaxa = allTaxa[!(allTaxa %in% badTaxa)]

#### Filtering ####
## Filter samples based on prevalence
ps.rel = transform_sample_counts(ps.sample, function(x) x / sum(x))
ps.rel.F = filter_taxa(ps.rel, function(x) mean(x) > 1e-5, T)
psF = prune_taxa(rownames(otu_table(ps.sample)) %in% rownames(otu_table(ps.rel.F)), ps.sample)
ps.clean = prune_samples(sample_sums(psF) >= 500, psF)
ps.clean = prune_taxa(taxa_sums(ps.clean) > 0, ps.clean)

## Remove samples with unknown infection status
ps.known = subset_samples(ps.clean, !(Status %in% 'Unknown'))
summary(sample_data(ps.known))

## Create subset for each host
ps.HAAM = subset_samples(ps.known, (Host %in% '`Amakihi'))
summary(sample_data(ps.HAAM))
ps.WAWE = subset_samples(ps.known, (Host %in% 'White-eye'))
summary(sample_data(ps.WAWE))
count(ps.WAWE@sam_data$Year == 2020)
## Transform to RELATIVE ABUNDANCES for weighted Unifrac and Bray-Curtis dissimilarity
ps.known.rel = transform_sample_counts(ps.known, function(x) x / sum(x))
ps.HAAM.rel = transform_sample_counts(ps.HAAM, function(x) x / sum(x))
ps.WAWE.rel = transform_sample_counts(ps.WAWE, function(x) x / sum(x))

rm(OTU, TAX, META, otu, taxonomy, phy_tree,
   ps.pa, ps.pa.neg, ps.pa.pos, contamdf.prev, ps.nocontam,
   mitochondria, chloroplast, eukaryote, unassigned, allTaxa, goodTaxa, badTaxa,
   ps.rel, ps.rel.F, psF)

save(ps, ps.known, ps.known.rel, ps.HAAM, ps.HAAM.rel, ps.WAWE, ps.WAWE.rel, file = 'ps.clean.RData')

### Statistics and Figures ##############################################################################
#### Assessing cloacal microbiome alpha diversity ####
## Negative binomial generalized linear model for malaria infection status
p1 = plot_richness(ps.known, measures = 'Shannon')
p1$data$value = round(p1$data$value, digits = 0)
glm.H = glm.nb(value ~ PosOrNeg + Host + Age + Region + Year + Month, data = p1$data)
res.H = car::Anova(glm.H)
res.H
capture.output(res.H, file = 'StatusGLM.tsv')

## Check model for overdispersion and fit
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type='pearson')
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(glm.H) # Not overdispersed (p > 0.05)

## nbGLM for malaria relative infection intensity
glm.H.Ct = glm.nb(value ~ RelInf + Host + Age + Region + Year + Month, data = p1$data)
res.H.Ct = car::Anova(glm.H.Ct)
res.H.Ct
capture.output(res.H.Ct, file = 'RelInfGLM.tsv')

## Check model for overdispersion and fit
overdisp_fun(glm.H.Ct) # Not overdispersed (p > 0.05)

##### Figure 2 ####
## Boxplots for alpha diversity
library('plyr')
p1.all = plot_richness(ps.known, measures = 'Shannon')

fig2a = ggplot(p1.all$data, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(color = 'black', width = .75,
               outlier.color = 'white', outlier.size = 0, lwd=0.36,
               position = position_dodge2(preserve = 'single')) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.8, color = 'black', pch = 21, show.legend = F) +
  scale_fill_manual(values = reg.col) +
  ylab('Shannon diversity index') + xlab('') + ylim(0,6) + theme_bw() +
  facet_wrap(~ Host, nrow = 2, strip.position = 'top', scales = 'free') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))

fig2a

fig2b = ggplot(p1.all$data, aes(x = Year, y = value, fill = Year)) +
  geom_boxplot(color = 'black', width = .75,
               outlier.color = 'white', outlier.size = 0, lwd=0.36,
               position = position_dodge2(preserve = 'single')) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.8, color = 'black', pch = 21, show.legend = F) +
  scale_fill_manual(values = yr.col) +
  ylab('Shannon diversity index') + xlab('') + ylim(0,6) + theme_bw() +
  facet_wrap(~ Host, nrow = 1, strip.position = 'top', scales = 'free') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))

fig2b

fig2c = ggplot(p1.all$data, aes(x = PosOrNeg, y = value, fill = PosOrNeg)) +
  geom_boxplot(color = 'black', width = .75,
               outlier.color = 'white', outlier.size = 0, lwd=0.36,
               position = position_dodge2(preserve = 'single')) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             alpha = 1, size = 1.8, color = 'black', pch = 21, show.legend = F) +
  scale_fill_manual(values = inf.col) +
  ylab('Shannon diversity index') + xlab('') + ylim(0,6) + 
  labs(fill = 'Infection status') + theme_bw() +
  facet_wrap(~ Host, nrow = 1, strip.position = 'top', scales = 'free') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))

fig2c

fig2bc = ggarrange(fig2b, fig2c, nrow = 1, ncol = 2, labels = c('(b)', '(c)'))
Fig2.plots = ggarrange(fig2a, fig2bc, nrow = 2, labels = '(a)', heights = c(3, 2))
Fig2.cap = ggarrange(ggparagraph(text = paste('FIGURE 2. Shannon diversity indices of Hawai`i `amakihi (N = 174) and warbling white-eye (N = 172) cloacal microbiomes sampled from (a) different regions on Hawai`i Island and categorized by (b) sampling year and (c) malaria infection status. Boxplots show the median and interquartile range for each population, and whiskers represent the 25th and 75th percentile.'), size = 16, color = 'black'))
Fig2 = ggarrange(Fig2.plots, Fig2.cap, nrow = 2, heights = c(7,1))
Fig2
ggsave(Fig2, filename = 'Figure2.pdf', device = 'pdf')
ggsave(Fig2, filename = 'Figure2.svg', device = 'svg')

##### Figure S1 ####
## Linear regression for alpha diversity vs relative infection intensity
figS1 = ggplot(p1.all$data, aes(x = RelInf, y = value, fill = Year, color = Year)) +
  geom_point(size = 4) + geom_smooth(aes(), method = 'lm', alpha = 0.25) +
  facet_wrap(~ Host, nrow = 2, strip.position = 'top') +
  ylab('Shannon diversity index') + xlab('Relative parasitemia intensity') + theme_bw() +
  ylim(0,5.1) + scale_x_log10(expand = c(0.01, 0), breaks = trans_breaks("log10", function(x) 10^x)) +
  theme(plot.tag.position = 'bottom',
        plot.tag = element_text(size = 16, color = 'black', hjust = 0, vjust = -1)) +
  scale_fill_manual(values = yr.col) +
  scale_color_manual(values = yr.col) +
  guides(fill = guide_legend(nrow = 2, byrow= T)) +
  theme(axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black'),
        strip.text = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = 'bottom',
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))
figS1

FigS1.cap = ggarrange(ggparagraph(text = paste('FIGURE S1. Shannon diversity indices versus malaria parasitemia intensity relative to the sampled population of Hawai`i `amakihi (N = 174) and warbling white-eye (N = 172) cloacal microbiomes in 2019 and 2020. Curves represent linear regressions and shaded bands denote the 95% confidence interval.'), size = 16, color = 'black'))
FigS1 = ggarrange(figS1, FigS1.cap, nrow = 2, heights = c(7,1))
FigS1
ggsave(FigS1, filename = 'FigureS1.pdf', device = 'pdf')
ggsave(FigS1, filename = 'FigureS1.svg', device = 'svg')

#### Assessing cloacal microbiome beta diversity ####
## Permutational multivariate analysis of variance model for malaria infection status
set.seed(350)
df = as(sample_data(ps.known.rel), 'data.frame')
ad = adonis(phyloseq::distance(ps.known.rel, method='bray') ~ PosOrNeg + Host + Age + Region + Year + Month, data = df, permutations = 1000)
ad 
capture.output(ad, file = 'StatusPERMANOVA.tsv')

## PERMANOVA for relative parasitemia intensity
set.seed(350)
df = as(sample_data(ps.known.rel), 'data.frame')
ad.Ct = adonis(phyloseq::distance(ps.known.rel, method='bray') ~ RelInf + Host + Age + Region + Year + Month, data = df, permutations = 1000)
ad.Ct
capture.output(ad.Ct, file = 'RelInfPERMANOVA.tsv')

##### Figure S2 ####
## Principal coordinate analyses for beta diversity
pS2 = ordinate(ps.known.rel, 'PCoA', 'bray')
gS2a = plot_ordination(ps.known.rel, pS2, color = 'Region', shape = 'Host')
gS2b = plot_ordination(ps.known.rel, pS2, color = 'Year', shape = 'Host')
gS2c = plot_ordination(ps.known.rel, pS2, color = 'PosOrNeg', shape = 'Host')
gS2d = plot_ordination(ps.known.rel, pS2, color = 'Age', shape = 'Host')

figS2a = gS2a +  geom_point(size = 4, aes(shape = Host)) +
  theme_bw() + ylim(-0.4, 0.4) + xlim(-0.4, 0.4) +
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = reg.col) +
  scale_color_manual(values = reg.col) +
  theme(axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black'),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = 'right',
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))
figS2a

figS2b = gS2b +  geom_point(size = 4, aes(shape = Host)) +
  theme_bw() + ylim(-0.4, 0.4) + xlim(-0.4, 0.4) +
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = yr.col) +
  scale_color_manual(values = yr.col) +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
  theme(axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black'),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = 'right',
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))
figS2b

figS2c = gS2c +  geom_point(size = 4, aes(shape = Host)) +
  theme_bw() + ylim(-0.4, 0.4) + xlim(-0.4, 0.4) +
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = inf.col) +
  scale_color_manual(values = inf.col) +
  labs(color = 'Infection Status') +
  theme(axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black'),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = 'right',
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))
figS2c

figS2d = gS2d +  geom_point(size = 4, aes(shape = Host)) +
  theme_bw() + ylim(-0.4, 0.4) + xlim(-0.4, 0.4) +
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = age.col) +
  scale_color_manual(values = age.col) +
  labs(color = 'Age') +
  theme(axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=8)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black'),
        strip.text.x = element_text(size=16),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = 'right',
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1,'cm'))
figS2d

FigS2.plots = ggarrange(figS2a, figS2b, figS2c, figS2d, nrow = 2, ncol = 2, labels = c('(a)', '(b)', '(c)', '(d)'))
FigS2.cap = ggarrange(ggparagraph(text = paste('FIGURE S2. Principal coordinate analysis plots of the beta diversity of Hawai`i `amakihi (N = 174) and warbling white-eye (N = 172) cloacal microbiomes based on Bray-Curtis dissimilarity matrices colored by (a) different regions on Hawai`i Island, (b) sampling year, (c) malaria infection status, and (d) age. Sample counts were transformed to relative abundances.'), size = 16, color = 'black'))
FigS2 = ggarrange(FigS2.plots, FigS2.cap, nrow = 2, heights = c(7,1))
FigS2
ggsave(FigS2, filename = 'FigureS2.pdf', device = 'pdf')
ggsave(FigS2, filename = 'FigureS2.svg', device = 'svg')

#### Top taxa ####
##### Top genera ####
gen = psmelt(ps.known.rel)
gen_summary = aggregate(gen$Abundance, by = list(Category = gen$Genus), FUN = sum) # calculate sum of abundance for each genera
gen_summary = as.data.frame(gen_summary[order(-gen_summary$x),])
capture.output(gen_summary, file = 'TopGenera.tsv')
print(tail(gen_WAWE_summary, 10))
gen.HAAM = psmelt(ps.HAAM.rel)
gen_HAAM_summary = aggregate(gen.HAAM$Abundance, by = list(Category = gen.HAAM$Genus), FUN = sum) # calculate sum of abundance for each genera
gen_HAAM_summary = as.data.frame(gen_HAAM_summary[order(-gen_HAAM_summary$x),])
capture.output(gen_HAAM_summary, file = 'TopGeneraHAAM.tsv')

gen.WAWE = psmelt(ps.WAWE.rel)
gen_WAWE_summary = aggregate(gen.WAWE$Abundance, by = list(Category = gen.WAWE$Genus), FUN = sum) # calculate sum of abundance for each genera
gen_WAWE_summary = as.data.frame(gen_WAWE_summary[order(-gen_WAWE_summary$x),])
capture.output(gen_WAWE_summary, file = 'TopGeneraWAWE.tsv')

##### Top phyla ####
phy = psmelt(ps.known.rel)
phy_summary = aggregate(phy$Abundance, by = list(Category = phy$Phylum), FUN = sum) # calculate sum of abundance for each phyla
phy_summary = as.data.frame(phy_summary[order(-phy_summary$x),])
capture.output(phy_summary, file = 'TopPhyla.tsv')

phy.HAAM = psmelt(ps.HAAM.rel)
phy_HAAM_summary = aggregate(phy.HAAM$Abundance, by = list(Category = phy.HAAM$Phylum), FUN = sum) # calculate sum of abundance for each phyla
phy_HAAM_summary = as.data.frame(phy_HAAM_summary[order(-phy_HAAM_summary$x),])
capture.output(phy_HAAM_summary, file = 'TopPhylaHAAM.tsv')

phy.WAWE = psmelt(ps.WAWE.rel)
phy_WAWE_summary = aggregate(phy.WAWE$Abundance, by = list(Category = phy.WAWE$Phylum), FUN = sum) # calculate sum of abundance for each phyla
phy_WAWE_summary = as.data.frame(phy_WAWE_summary[order(-phy_WAWE_summary$x),])
capture.output(phy_WAWE_summary, file = 'TopPhylaWAWE.t')

##### Figure 3 ####
## Top phyla barplots
phy = psmelt(ps.known.rel)
taxa_summary = aggregate(phy$Abundance, by = list(Category = phy$Phylum), FUN = sum) # calculate sum of abundance for each phyla
taxa_summary = taxa_summary[order(-taxa_summary$x),]
## Reduce number of categories by combining rare ASVs
table(phy$Phylum, useNA ='ifany') # table of phyla
phy$Phylum[phy$Phylum == ''] = 'Unassigned' # change if there are empty names to NA
table(phy$Phylum, useNA = 'ifany') # check
list(as.character(taxa_summary[-c(1:11), 1])) # taxa NOT in the top 10 for Abundance, using 11 to account for unassigned ASVs
library('plyr')
phy$Phylum = revalue(phy$Phylum, c('Chloroflexi' = 'Other', 'WPS-2' = 'Other',
                                   'Patescibacteria' = 'Other', 'Fusobacteriota' = 'Other',
                                   'Planctomycetota' = 'Other', 'Deinococcota' = 'Other',
                                   'Gemmatimonadota' = 'Other', 'Bdellovibrionota' = 'Other',
                                   'Armatimonadota' = 'Other', 'Synergistota' = 'Other',
                                   'Methylomirabilota' = 'Other', 'Entotheonellaeota' = 'Other',
                                   'Nitrospirota' = 'Other'))
table(phy$Phylum, useNA = 'ifany')
table(is.na(phy$Phylum))
phy$Phylum = factor(phy$Phylum, levels = c('Acidobacteriota', 'Actinobacteriota',
                                           'Bacteroidota', 'Campilobacterota',
                                           'Cyanobacteria', 'Firmicutes',
                                           'Myxococcota','Proteobacteria',
                                           'Spirochaetota','Verrucomicrobiota',
                                           'Other', 'Unassigned'))
levels(phy$Phylum)

fig3a = ggplot(phy, aes(x = Host, y = Abundance, fill = Phylum)) +
  geom_bar(stat = 'identity', position = 'fill', show.legend = F) +
  facet_grid(~Region, scale = 'free', space = 'free_x') +
  scale_y_continuous(expand = c(0,0), labels = percent) +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  labs(x = '', y = 'Relative abundance') +
  guides(fill = guide_legend(ncol = 1, bycol = TRUE)) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 7)),
        axis.title.y = element_text(size = 16, margin = margin(r = 2)),
        axis.text.y  = element_text(size = 14, color= 'black'),
        axis.text.x  = element_text(size = 14, color= 'black'),
        strip.text.x = element_text(size = 14, color = 'black'),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size = 14),
        legend.key.size = unit(0.75, 'cm'),
        legend.text  = element_text(size = 14))

fig3a

fig3b = ggplot(phy, aes(x = Host, y = Abundance, fill = Phylum)) +
  geom_bar(stat = 'identity', position = 'fill', show.legend = F) +
  facet_grid(~Year, scale = 'free', space = 'free_x') +
  scale_y_continuous(expand = c(0,0), labels = percent) +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  labs(x = '', y = 'Relative abundance') +
  guides(fill = guide_legend(ncol = 1, bycol = TRUE)) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 7)),
        axis.title.y = element_text(size = 16, margin = margin(r = 2)),
        axis.text.y  = element_text(size = 14, color= 'black'),
        axis.text.x  = element_text(size = 14, color= 'black'),
        strip.text.x = element_text(size = 14, color = 'black'),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size = 14),
        legend.key.size = unit(0.75, 'cm'),
        legend.text  = element_text(size = 14))

fig3b

fig3c = ggplot(phy, aes(x = Host, y = Abundance, fill = Phylum)) +
  geom_bar(stat = 'identity', position = 'fill', show.legend = T) +
  facet_grid(~PosOrNeg, scale = 'free', space = 'free_x') +
  scale_y_continuous(expand = c(0,0), labels = percent) +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  labs(x = '', y = 'Relative abundance') +
  guides(fill = guide_legend(ncol = 1, bycol = TRUE)) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 7)),
        axis.title.y = element_text(size = 16, margin = margin(r = 2)),
        axis.text.y  = element_text(size = 14, color= 'black'),
        axis.text.x  = element_text(size = 14, color= 'black'),
        strip.text.x = element_text(size = 14, color = 'black'),
        strip.background = element_blank()) +
  theme(panel.border = element_rect(color = 'grey29'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.title = element_text(size = 14),
        legend.key.size = unit(0.75, 'cm'),
        legend.text  = element_text(size = 14))

fig3c

fig3bc = ggarrange(fig3b, fig3c, nrow = 1, labels = c('(b)', '(c)'), widths = c(1,1.5))
Fig3.plots = ggarrange(fig3a, fig3bc, nrow = 2, labels = c('(a)', ''))
Fig3.cap = ggarrange(ggparagraph(text = paste('FIGURE 3. Relative abundances of the 10 most common phyla in the cloacal microbiomes of Hawai`i `amakihi (174) and warbling white-eyes (172) organized by (a) sampling region, (b) sampling year, and (c) malaria infection status.'), size = 16, color = 'black'))
Fig3 = ggarrange(Fig3.plots, Fig3.cap, nrow = 2, heights = c(7,1))
Fig3
ggsave(Fig3, filename = 'Figure3.pdf', device = 'pdf', units = 'cm', height = 50, width = 30)
# ggsave(Fig3, filename = 'Figure3.svg', device = 'svg', units = 'cm', height = 50, width = 30) # Very large file

#### Investigating differential ASV abundance ####
## Comparing infected 'amakihi vs uninfected 'amakihi (remove actively infected birds)
## Create subsets to remove 'amakihi with active infection
ps.noactive = subset_samples(ps.HAAM, !(Status == 'Active'))
summary(sample_data(ps.noactive))

## Convert phyloseq object to deseq2 object
des.obj.noactive = phyloseq_to_deseq2(ps.noactive, ~ Region + Year + Month + PosOrNeg)
des.obj.noactive = DESeq(des.obj.noactive, sfType = 'poscounts', test = 'Wald', fitType = 'local')
noactive.norm.tab = as.data.frame(counts(des.obj.noactive, normalized = T))
write.table(noactive.norm.tab, file = 'DeSeq2_normtab.tsv', row.names = T, sep ='\\t') # Table used in identifying ASVs positively associated with disease pressure
des.res.noactive = results(des.obj.noactive, cooksCutoff = F, 
                           contrast = c('PosOrNeg', 'Infected', 'Uninfected'), alpha = 0.05, independentFiltering = T)

## Find significantly more abundant ASVs (padj < 0.5)
sigtab.noactive = des.res.noactive[which(des.res.noactive$padj < 0.05), ]
sigtab.noactive = cbind(as(sigtab.noactive, 'data.frame'), as(tax_table(ps.noactive)[rownames(sigtab.noactive), ], 'matrix'))
write.table(sigtab.noactive, file = 'Deseq2_A.tsv')

##### Figure 4 ####
## Differential ASV abundance
## ASVs positively associated with 'amakihi infected with avian malaria
library('plyr')
library('dplyr')

sigtab.noactive = sigtab.noactive[which(sigtab.noactive$log2FoldChange > 0), ]
x = tapply(sigtab.noactive$log2FoldChange, sigtab.noactive$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.noactive = subset(sigtab.noactive, !(Phylum %in% ''))
sigtab.noactive$Phylum = sub('^$', 'Unassigned', sigtab.noactive$Phylum)
sigtab.noactive$Phylum = factor(as.character(sigtab.noactive$Phylum), levels=names(x))
x = tapply(sigtab.noactive$log2FoldChange, sigtab.noactive$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.noactive$Genus = sub('^$', 'Unassigned', sigtab.noactive$Genus)
sigtab.noactive$Genus = factor(as.character(sigtab.noactive$Genus), levels=names(x))
sigtab.noactive$Genus = revalue(sigtab.noactive$Genus, c('Methylobacterium-Methylorubrum' = 'Methylobacterium'))
sigtab.noactive$Genus = revalue(sigtab.noactive$Genus, c('Candidatus_Rhabdochlamydia' = 'Rhabdochlamydia'))
sigtab.noactive$Genus = revalue(sigtab.noactive$Genus, c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = 'Rhizobium sensu lato'))
sigtab.noactive$Genus = revalue(sigtab.noactive$Genus, c('Nostoc_PCC-7107' = 'Nostoc'))

fig4a = ggplot(sigtab.noactive, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(size=6) +
  theme_bw() + labs(title = "Chronically infected versus uninfected `amakihi") +
  ylab('Log2 fold change') + xlab('') + ylim(0, 31) +
  scale_color_viridis_d() +
  theme(title  = element_text(size=14, color='black'),
        axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=10)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black', angle = -60, hjust = 0, vjust=0.5)) +
  theme(legend.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16),
        legend.key.size = unit(0.9,'cm'))
fig4a

## Comparing uninfected white-eyes and uninfected 'amakihi
## Create subsets to remove infected birds
ps.uninf = subset_samples(ps.known, (PosOrNeg %in% 'Uninfected'))
summary(sample_data(ps.uninf))

## Convert phyloseq object to deseq2 object
des.obj.uninf = phyloseq_to_deseq2(ps.uninf, ~ Region + Year + Month + Host)
des.obj.uninf = DESeq(des.obj.uninf, sfType = 'poscounts', test = 'Wald', fitType = 'local')
des.res.uninf = results(des.obj.uninf, cooksCutoff = F, 
                        contrast = c('Host', 'White-eye', '`Amakihi'), alpha = 0.05, independentFiltering = T)

## Find significantly differentially abundant ASVs
sigtab.uninf = des.res.uninf[which(des.res.uninf$padj < 0.05), ]
sigtab.uninf = cbind(as(sigtab.uninf, 'data.frame'), as(tax_table(ps.uninf)[rownames(sigtab.uninf), ], 'matrix'))
write.table(sigtab.uninf, file = 'Deseq2_B.tsv')

## Figure 4b - ASVs positively associated with white-eyes compared to uninfected 'amakihi
library('plyr')
library('dplyr')

sigtab.uninf = sigtab.uninf[which(sigtab.uninf$log2FoldChange > 0), ]
x = tapply(sigtab.uninf$log2FoldChange, sigtab.uninf$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.uninf = subset(sigtab.uninf, !(Phylum %in% ''))
sigtab.uninf$Phylum = sub('^$', 'Unassigned', sigtab.uninf$Phylum)
sigtab.uninf$Phylum = factor(as.character(sigtab.uninf$Phylum), levels=names(x))
x = tapply(sigtab.uninf$log2FoldChange, sigtab.uninf$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.uninf$Genus = sub('^$', 'Unassigned', sigtab.uninf$Genus)
sigtab.uninf$Genus = factor(as.character(sigtab.uninf$Genus), levels=names(x))
sigtab.uninf$Genus = revalue(sigtab.uninf$Genus, c('Methylobacterium-Methylorubrum' = 'Methylobacterium'))
sigtab.uninf$Genus = revalue(sigtab.uninf$Genus, c('Clostridium_sensu_stricto_1' = 'Clostridium sensu stricto'))

fig4b = ggplot(sigtab.uninf, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(size=6) +
  theme_bw() + labs(title = "Uninfected white-eyes versus uninfected `amakihi") +
  ylab('Log2 fold change') + xlab('') + ylim(0, 31) +
  scale_color_viridis_d() +
  theme(title  = element_text(size=14, color='black'),
        axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=10)),
        axis.text.y  = element_text(size=14, color='black'),
        axis.text.x  = element_text(size=14, color='black', angle = -60, hjust = 0, vjust=0.5)) +
  theme(legend.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16),
        legend.key.size = unit(0.9,'cm'))
fig4b

## Comparing uninfected white-eyes and chronically infected 'amakihi
## Create subsets to remove uninfected 'amakihi and infected white-eyes
ps.inf = subset_samples(ps.known, !(Host == '`Amakihi' & PosOrNeg == 'Uninfected'))
ps.inf = subset_samples(ps.inf, !(Host == '`Amakihi' & Status == 'Active'))
ps.inf = subset_samples(ps.inf, !(Host == 'White-eye' & PosOrNeg == 'Infected'))
summary(sample_data(ps.inf))

## Convert phyloseq object to deseq2 object
des.obj.inf = phyloseq_to_deseq2(ps.inf, ~ Region + Year + Host)
des.obj.inf = DESeq(des.obj.inf, sfType = 'poscounts', test = 'Wald', fitType = 'local')
des.res.inf = results(des.obj.inf, cooksCutoff = F, 
                      contrast = c('Host', 'White-eye', '`Amakihi'), alpha = 0.05, independentFiltering = T)

## Find significantly differentially abundant ASVs
sigtab.inf = des.res.inf[which(des.res.inf$padj < 0.05), ]
sigtab.inf = cbind(as(sigtab.inf, 'data.frame'), as(tax_table(ps.inf)[rownames(sigtab.inf), ], 'matrix'))
write.table(sigtab.inf, file = 'Deseq2_C.tsv')

## Figure 4c - ASVs positively associated with white-eyes compared to infected 'amakihi
library('plyr')
library('dplyr')

sigtab.inf = sigtab.inf[which(sigtab.inf$log2FoldChange > 0), ]
x = tapply(sigtab.inf$log2FoldChange, sigtab.inf$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.inf = subset(sigtab.inf, !(Phylum %in% ''))
sigtab.inf$Phylum = sub('^$', 'Unassigned', sigtab.inf$Phylum)
sigtab.inf$Phylum = factor(as.character(sigtab.inf$Phylum), levels=names(x))
x = tapply(sigtab.inf$log2FoldChange, sigtab.inf$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.inf$Genus = sub('^$', 'Unassigned', sigtab.inf$Genus)
sigtab.inf$Genus = factor(as.character(sigtab.inf$Genus), levels=names(x))
sigtab.inf$Genus = revalue(sigtab.inf$Genus, c('Methylobacterium-Methylorubrum' = 'Methylobacterium'))
sigtab.inf$Genus = revalue(sigtab.inf$Genus, c('Burkholderia-Caballeronia-Paraburkholderia' = 'Burkholderia sensu lato'))

fig4c = ggplot(sigtab.inf, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(size=6) +
  theme_bw()+ labs(title = "Uninfected white-eyes versus chronically infected `amakihi") +
  ylab('Log2 fold change') + xlab('Genus') + ylim(0, 31) +
  scale_color_viridis_d() +
  theme(title  = element_text(size=14, color='black'),
        axis.title.x = element_text(size=16, margin=margin(t=5)),
        axis.title.y = element_text(size=16, margin=margin(r=10)),
        axis.text.y  = element_text(size=16, color='black'),
        axis.text.x  = element_text(size=14, color='black', angle = -60, hjust = 0)) +
  theme(legend.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16),
        legend.key.size = unit(0.9,'cm'))

fig4c

Fig4.plots = ggarrange(fig4a, fig4b, fig4c, nrow = 3, labels = c('(a)', '(b)', '(c)'))
Fig4.cap = ggarrange(ggparagraph(text = paste('FIGURE 4. Log2 fold increase of ASV abundance (p-adj < 0.05) in the cloacal microbiomes of (a) chronically malaria infected (N = 40) versus uninfected (N = 125) Hawai`i `amakihi, (b) uninfected warbling white-eyes (N = 133) versus uninfected `amakihi (N = 125), and (c) uninfected white-eyes (N = 133) versus chronically infected `amakihi (N = 40). Points indicate ASVs more abundant in infected `amakihi in (a), and in white-eyes in (b, c). NA indicates ASVs without taxonomic genus classification.'), size = 16, color = 'black'))
Fig4 = ggarrange(Fig4.plots, Fig4.cap, nrow = 2, heights = c(7,1))
Fig4
ggsave(Fig4, filename = 'Figure4.pdf', device = 'pdf', units = 'cm', height = 45, width = 35)
ggsave(Fig4, filename = 'Figure4.svg', device = 'svg', units = 'cm', height = 45, width = 35)

#### Identify ASVs positively associated with disease pressure ####
prev.dat = read.table(file = 'Prevalence.tsv', header = T, sep = '\\t')

asv.dat = read.table(file = 'DeSeq2_normtab.tsv', header = T, sep = '\\t')
asv.dat = t(asv.dat)
prev.dat.all = subset(prev.dat, SampleID %in% rownames(asv.dat))

asv_prev_lm = function(abs, prev){  
  corr_results = data.frame(asv=character(),  # empty table to put output of each linear model
                            coeff=numeric(), SE=numeric(),adjR2=numeric(), 
                            Fvalue=numeric(), df=numeric(), resid=numeric(), P=numeric())
  
  abs = abs[,colSums(abs) > 1] # removes ASVs that are absent
  for(i in 1:ncol(abs)){ # for each ASV (which is a column) do the following
    prev$abs_i = abs[,i] # adds a column of abundance value for each ASV to prevalence dataframe
    rm.model = lm(log10(abs_i+0.1) ~ Prevalence, data = prev) # model
    Ftest = anova(rm.model, test='F') # F-test of model
    res = data.frame(asv=colnames(abs)[i], coeff=coef(summary(rm.model))[2,1], # fill output table
                     SE=coef(summary(rm.model))[2,2], adjR2=summary(rm.model)$adj.r.squared,
                     Fvalue=Ftest$'F value'[1], df=Ftest$'Df'[1], resid=Ftest$'Df'[2], P=Ftest$'Pr(>F)'[1]) 
    corr_results = rbind(corr_results, res) # correlation test of rests
  }
  return(corr_results)
}

lm_results = asv_prev_lm(asv.dat, prev.dat.all)
lm_results$Padjust = p.adjust(lm_results$P, method = 'BH') # Get adjusted p-values
lm_sig = lm_results[lm_results$Padjust < 0.05, ] # Isolate just significant factors
lm_sig
lm_sig_pos = lm_sig %>% filter(coeff >= 0)
lm_sig_pos # 3 ASVs associated with disease pressure
capture.output(lm_sig_pos, file = 'PressureASVs.tsv')
