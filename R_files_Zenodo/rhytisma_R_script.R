############################################################################################
#decomposition
df <- read.csv("rhytisma_AFDM_spor_data.csv")
head(df)

str(df)

model <- lm(ln_AFDM_remaining ~ Treatment * Harvest_days, data=df)
model

summary(model)

anova(model)

R <- subset(df, df$Treatment == "InR")
B <- subset(df, df$Treatment == "InB")
U <- subset(df, df$Treatment == "CtrlU")

lmR <- lm(R$ln_AFDM_remaining ~ R$Harvest_days)
lmR

lmB <- lm(B$ln_AFDM_remaining ~ B$Harvest_days)
lmB

lmU <- lm(U$ln_AFDM_remaining ~ U$Harvest_days)
lmU

plot(AFDM_remaining ~ Harvest_days, data=df, xlab="Harvest (days)",ylab="% AFDM Remaining",type='n')

points(R$AFDM_remaining~R$Harvest_days, pch=15,col="black")
points(U$AFDM_remaining~U$Harvest_days, pch=2,col="black")
points(B$AFDM_remaining~B$Harvest_days, pch=16,col="darkgrey")

lmB
lmR
lmU

abline(lm(B$AFDM_remaining ~ B$Harvest_days), lty=2, col="darkgrey")
abline(lm(R$AFDM_remaining ~ R$Harvest_days), lty=1,col="black")
abline(lm(U$AFDM_remaining ~ U$Harvest_days),lty=3,col="black")
legend("topright", c("Rhytisma-infected","bullseye-infected","uninfected"), lty=c(1,2,3), pch=c(15,16,2), col=c("black","darkgrey","black"))

dev.off()
############################################################################################

#sporulation rate
dfsum
948.96425+618.58957 #sporulation_rate+se
38.99219+13.62726
63.52872+42.67612

limits <- aes(ymax = dfsum$sporulation_rate + dfsum$se, ymin = dfsum$sporulation_rate)
p <- ggplot(data = dfsum, aes(x =Treatment, y =sporulation_rate, fill = Treatment))
p
Pgraph <- p + geom_bar(stat = "identity",position = position_dodge(0.9),
                       color="#000000")  + geom_errorbar(limits, position = position_dodge(0.9), 
                                                         width = 0.25) + labs(x = "Tissue infection type", 
                                                                              y = "Sporulation rate (conidia/mg AFDM * day)") + 
  ggtitle("") + scale_fill_manual(name = "Treatment",
                                  values=c("grey", "black","white")) + theme(legend.position="none",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                                                                             panel.background=element_rect(fill="white"),
                                                                             panel.border=element_rect(fill=NA,color="black"),axis.text.x = element_text (color="black"),
                                                                             axis.text.y = element_text (color="black")) + annotate("text",x=3,y=1600,label="B") + 
  annotate("text",x=1,y=85,label="A")+annotate("text",x=2,y=140,label="A")
Pgraph
############################################################################################

#community analyses and figures were completed using PC-ORD
############################################################################################

#supplementary phyloseq figures
library(phyloseq)
library(ggplot2)
fungi <- read.csv("fungi.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
taxa <- read.csv("taxa.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
env <- read.csv("env.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

fungi <- as.matrix(fungi)
taxa <- as.matrix(taxa)
#env should be df not matrix

OTU = otu_table(fungi, taxa_are_rows = TRUE)
TAX = tax_table(taxa)

physeq = phyloseq(OTU, TAX)

envdata = sample_data(env)
sample_names(physeq) #make sure names match
sample_names(envdata) #make sure names match

merged = merge_phyloseq(physeq, envdata)

#Convert to relative abundance and remove low abundance taxa (<0.1%)
relabund = transform_sample_counts(merged, function(x) x / sum (x))
no_low_taxa = filter_taxa(relabund, function(x) mean (x) > 1e-3, TRUE)

no_low_taxa_fungi <- otu_table(no_low_taxa, taxa_are_rows = FALSE)
write.csv(no_low_taxa_fungi, file = 'no_low_taxa_fungi.csv')

relabund_fungi <- otu_table(relabund, taxa_are_rows = FALSE)
write.csv(relabund_fungi, file = 'relabund_fungi2.csv')

#per Brett's suggestion, removing any OTUs not resolved past order, which will remove 27 #OTUS:
#removes rows 2, 8 ,28, 44, and 46
#remaining taxa saved as no_low_taxa_fungi-resolved.csv in the usual dir

fungi <- read.csv("fungi-resolved2.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
taxa <- read.csv("taxa_fungi-resolved.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
env <- read.csv("env.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

fungi <- as.matrix(fungi)
taxa <- as.matrix(taxa)
#env should be df not matrix

OTU = otu_table(fungi, taxa_are_rows = TRUE)
TAX = tax_table(taxa)

physeq = phyloseq(OTU, TAX)

envdata = sample_data(env)
sample_names(physeq) #make sure names match
sample_names(envdata) #make sure names match

merged = merge_phyloseq(physeq, envdata)

#Convert to relative abundance and remove low abundance taxa (<0.1%)
relabund = transform_sample_counts(merged, function(x) x / sum (x))
no_low_taxa = filter_taxa(relabund, function(x) mean (x) > 1e-3, TRUE)

otus <- otu_table(no_low_taxa, taxa_are_rows = FALSE)
write.csv(otus, file = 'fungi-resolved_OTUs2.csv')

fungi <- read.csv("fungi-resolved_OTUs2.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
taxa <- read.csv("taxa_fungi-resolved.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
env <- read.csv("env.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)

fungi <- as.matrix(fungi)
taxa <- as.matrix(taxa)
#env should be df not matrix

OTU = otu_table(fungi, taxa_are_rows = TRUE)
TAX = tax_table(taxa)

physeq = phyloseq(OTU, TAX)

envdata = sample_data(env)
sample_names(physeq) #make sure names match
sample_names(envdata) #make sure names match

merged = merge_phyloseq(physeq, envdata)
relabund = transform_sample_counts(merged, function(x) x / sum (x))

otus <- otu_table(relabund, taxa_are_rows = FALSE)
write.csv(otus, file = 'fungi-resolved_rel2.csv')
#this was renamed to a .xlsx for calculations
#new sheet with just values given same name as fungi-resolved_rel2_final.csv to
#avoid confusion


#after relativizing the relative abundance matrix (2 punches for A0, so divide all counts in #each column by 2, etc. so that the rel abund on y=1)
fungi <- read.csv("fungi-resolved_rel2_final.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
taxa <- read.csv("taxa_fungi-resolved_OTUs.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
env <- read.csv("env.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
fungi <- as.matrix(fungi)
taxa <- as.matrix(taxa)

levels(env$Treatment) <- c("bullseye-infected","Rhytisma-infected","lesion-free")

env$Harvest <- as.factor(env$Harvest)
#env should be df not matrix

OTU = otu_table(fungi, taxa_are_rows = TRUE)
TAX = tax_table(taxa)

physeq = phyloseq(OTU, TAX)
envdata = sample_data(env)
sample_names(physeq) #make sure names match
sample_names(envdata) #make sure names match
merged = merge_phyloseq(physeq, envdata)
#rolled 43-color palette w/ soft k-means on IWantHue
f <- plot_bar(merged, "Harvest", fill = "Family", facet_grid = ~Treatment) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom") + scale_fill_manual(values=c("#dc375f", "#9d4e55",
                                                                                                                                                                                                                               "#d96e76","#b63931","#e86a55", "#d75623", "#dd9a78", "#9c5b2f", "#dd9653",
                                                                                                                                                                                                                               "#df9529", "#9e7423", "#c6ae31", "#766f2f", "#bcb365", "#899332", "#a8bd33",
                                                                                                                                                                                                                               "#567531", "#6eb729", "#86be59", "#538d30", "#339e36", "#53cb5b", "#79b779",
                                                                                                                                                                                                                               "#2f7c52", "#38c481", "#56bca0", "#3abec8", "#53a4d6", "#5b8cdf", "#5765a5",
                                                                                                                                                                                                                               "#6666d4", "#b098dc", "#8555a3", "#a44ec9", "#cf7fe1", "#ac3d96", "#9e5d8d",
                                                                                                                                                                                                                               "#df51b7", "#e080bf", "#d06792", "#e0458b", "#a13462", "#e796b0")) +geom_bar(colour="transparent", stat="identity")
f
ggsave(f, filename = "fig_fun.jpeg", dpi = 300)

#for bacteria
library(phyloseq)
bac <- read.csv("bac-resolved_OTUs_rel.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
taxa <- read.csv("taxa_bac_resolved2.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
env <- read.csv("env_b.csv", header=TRUE, row.names=1, stringsAsFactors=TRUE)
levels(env$Treatment) <- c("bullseye-infected","Rhytisma-infected","lesion-free")

bac <- as.matrix(bac)
taxa <- as.matrix(taxa)
env$Harvest <- as.factor(env$Harvest)
#env should be df not matrix

OTU = otu_table(bac, taxa_are_rows = TRUE)
TAX = tax_table(taxa)

physeq = phyloseq(OTU, TAX)
envdata = sample_data(env)
sample_names(physeq) #make sure names match
sample_names(envdata) #make sure names match
merged = merge_phyloseq(physeq, envdata)


#now for family as the fill
f <- plot_bar(merged, "Harvest", fill = "Family", facet_grid = ~Treatment) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position="bottom") + scale_fill_manual(values=c("#df3666","#b04e57","#e3897a","#d14832","#a0542d","#d97c3f","#dd9a35",
f
ggsave(f, filename = "fig_bac.jpeg", dpi = 300)