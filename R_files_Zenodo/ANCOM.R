########## packges load ###################################################################################
library(vegan)
library(car)
library (DESeq2); packageVersion("DESeq2")
library(Rfast) #to take row maximums
library (ggthemes)
library(data.table); packageVersion("data.table")
library(adaptiveGPCA)
library(ggrepel)
library(phyloseq)
library(exactRankTests)
library(nlme)
library(dplyr)
library(compositions)
library(readr)
library(tidyverse)
library(ggpubr)
library(rstatix)

source("ANCOM_master/scripts/ancom_v2.1.R")
# https://github.com/FrederickHuangLin/ANCOM

########### load clean and rerified data for analysis ###############################################################
data_in <- readRDS(file="data_in_rarefied_8000.rds")
taxa_names(data_in) <- paste0("SV", seq(ntaxa(data_in))) # rename to short names

OTU = as(otu_table(data_in), "matrix")
OTUdf = data.frame(OTU)
OTUdf = t(OTUdf)

#rownames(OTUdf) are the taxa
meta<-data.frame(sample_data(data_in))
meta$SampleID<-rownames(meta)

Tax<-as(tax_table(data_in), "matrix")

##--- making the analysis on Genera-----------------------------------------------------
data_in_Genus <- tax_glom(data_in, taxrank = 'Genus')
taxa_names(data_in_Genus) <- paste0("Gen", seq(ntaxa(data_in_Genus))) # rename to short names

OTUGen = as(otu_table(data_in_Genus), "matrix")
OTUGendf = data.frame(OTUGen)
OTUGendf = t(OTUGendf)

#rownames(OTUdf) are the taxa
metaGen<-data.frame(sample_data(data_in_Genus))
metaGen$SampleID<-rownames(metaGen)

TaxGen<-as(tax_table(data_in_Genus), "matrix")

#----- Step 1: Data preprocessing---------------------------------------------------------------------
feature_table = OTUGendf; sample_var = "SampleID"; group_var = NULL
out_cut = 0.05 
zero_cut = 0.90 # Taxa with proportion of zeroes greater than zero_cut are not included in the analysis
#zero_cut = 0.75
lib_cut = 1000 # Nor relevant for rarefied data==>Shoud we use non-rarefied
neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, metaGen, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

#----- Step 2: ANCOM----------------------------------------------------------------------------------

# with  arandom factor
main_var = "ControlGroupWeek"; p_adj_method = "BH"; alpha = 0.01
adj_formula = NULL; rand_formula = "~ 1 | WeekSinceBreeding"
t_start = Sys.time()
resGen = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
             alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # 19.56198 secs

ResultsAncomGen <- data.frame(resGen["out"])
#SigDifGen <-  ResultsAncomGen$out.taxa_id[ResultsAncomGen$out.detected_0.7]
SigDifGen091 <-  ResultsAncomGen$out.taxa_id[ResultsAncomGen$out.detected_0.9]
W <-ResultsAncomGen$out.W[ResultsAncomGen$out.detected_0.9]
# take cut off of 85%
Wper<-ResultsAncomGen$out.W/max(ResultsAncomGen$out.W)
SigDifGen091 <-  ResultsAncomGen$out.taxa_id[Wper>0.85]
W <-ResultsAncomGen$out.W[Wper>0.85]
# make data frame of the different genera
AllGen<-rownames(TaxGen)
INDKeep<-c()
for(i in 1:length(SigDifGen091)){
  IND <- which(AllGen==SigDifGen091[i])
  INDKeep[i]<-IND
}

TaxDif<-TaxGen[INDKeep, ]
TaxDif<-as.data.frame(TaxDif)
TaxDif$W<-W
TaxDif<-TaxDif[with(TaxDif, order(Genus)),] # order by genus, the same way the levels are ordered bt defult
TaxDif$ord<-seq(1, nrow(TaxDif), by=1) 
# reaorder the levels of the sample to create x ordered axis ordered 
TaxDif<-TaxDif[with(TaxDif, order(Phylum, Class,Genus)),]

#----- Step 3: Plot and add statistics----------------------------------------------------------------------------------

#--save the significant genera table--------
saveRDS (TaxDif, file="DifGeneraAncom.rds")
TaxDif <- readRDS(file="DifGeneraAncom.rds")
#--plot the significant ones counts --------
dat_gen <- psmelt(data_in_Genus)

# take only the relevant genera for all our samples
dat_DifGen <- as.data.frame(matrix(data=NA,ncol=ncol(dat_gen))) #@ create empty dataframe
colnames(dat_DifGen)<-colnames(dat_gen)
SpOerder<-c()
for(i in 1:nrow(TaxDif)){
  dat_DifGen <- rbind(dat_DifGen,dat_gen[as.character(dat_gen$Genus)==TaxDif$Genus[i],])
}
dat_DifGen<-dat_DifGen[-1,]
# we want to order it by Phylum, class, genus
dat_DifGen$Genus<-as.factor(dat_DifGen$Genus)
dat_DifGen$Genus<- factor(dat_DifGen$Genus,levels(dat_DifGen$Genus)[TaxDif$ord])
dat_DifGen$ControlGroupWeek<- factor(dat_DifGen$ControlGroupWeek,
                                 levels = c("Pre_migration", "Fall", "Winter_fields","Winter_feeding_station"))

Col<-c("#ffa321","#48b823","#559bab","red4")

#add 1 to make log
dat_DifGen$LogAbundance<-log(dat_DifGen$Abundance+1)

#---MAKE FACET PLOT-------------------------   
# reorder again for facet
p <- ggplot(dat_DifGen, aes(x = factor(ControlGroupWeek), y = LogAbundance, fill=factor(ControlGroupWeek))) +
  geom_boxplot(alpha = 1,width = 0.9, fatten = 4) +
  facet_wrap(. ~ Genus,ncol=3)+
  scale_color_manual(values =Col)+
  scale_fill_manual(values =Col)+
  ylab("log(abundance)")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                               size = 0.5, linetype = "solid"),axis.text.x=element_blank(),
                                strip.background =element_rect(fill="#e8e8e8"),strip.text = element_text(size = 10))
######################################
# Add significance
######################################

try <-dat_DifGen[c(55,47,54)]

##### use significance from Pairwise comparisons using Wilcoxon's test ##########
pwc2 <- try %>% 
  group_by(Genus) %>%
  wilcox_test(LogAbundance ~ ControlGroupWeek, p.adjust.method = "bonferroni")
pwc2<- pwc2 %>% add_xy_position(x = "ControlGroupWeek")

pwc2$p.adj.signif<-"ns"
pwc2$p.adj.signif[pwc2$p.adj*6<0.05] <- "*"
pwc2$p.adj.signif[pwc2$p.adj*6<0.001] <- "**"

# plot with stats
ggboxplot(try, x = "ControlGroupWeek", y = "LogAbundance", facet.by = "Genus",
          fill = "ControlGroupWeek", palette = Col, alpha = 1,width = 0.9, fatten = 4,ncol=6) +
  stat_pvalue_manual(pwc2, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),axis.text.x=element_blank(),
        strip.background =element_rect(fill="#e8e8e8"),strip.text = element_text(size = 10))


################################################
###############################################

#---MAKE FACET PLOT OF ONLY THE MOST ABYNDANT/PREVELENT GENERA-------------------------
ps2_rel_gen = transform_sample_counts(data_in_Genus, function(x){x / sum(x)}) # convert to relative obundance
# create dataframe from phyloseq object
dat_gen <- psmelt(ps2_rel_gen)
# convert Genus to a character vector from a factor
dat_gen$Genus <- as.character(dat_gen$Genus)

aa_Gen<-plyr::ddply(dat_gen, ~Genus, summarize, 
                    mean.abundance = round(mean(Abundance), digits = 3),
                    sd.abundance=round(sd(Abundance), digits = 3),
                    meadian.abundance=round(median(Abundance), digits = 3))
GeneraMain <- aa_Gen[aa_Gen$mean.abundance>0.002,] # more than 0.1%

#-- Only most prevelent Genera------------------------------
prev0_genera  = apply(X = otu_table(data_in_Genus),
                      MARGIN = ifelse(taxa_are_rows(data_in_Genus), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})
Tax<-tax_table(data_in_Genus)
prevdf_genera  = data.frame(Prevalence = prev0_genera ,
                            Tax@.Data)

NumberOfSamples<-nrow(sample_data(data_in_Genus))
prevdf_genera$Relative.Prevelence  = prevdf_genera$Prevalence/NumberOfSamples
GeneraPrevelent <- prevdf_genera[prevdf_genera$Relative.Prevelence>0.5,] # more than 50% prevelence

## take indexes of the abundant/prevelent different genera

INDKeep<-c()
for(i in 1:nrow(GeneraMain)){
  IND <- which(dat_DifGen$Genus==GeneraMain$Genus[i])
  if (!isEmpty(IND)) {
  INDKeep<-rbind(INDKeep,IND)
  }
}

INDKeepPR<-c()
for(i in 1:nrow(GeneraPrevelent)){
  IND <- which(dat_DifGen$Genus==GeneraPrevelent$Genus[i])
  if (!isEmpty(IND)) {
    INDKeepPR<-rbind(INDKeepPR,IND)
  }
}


dat_DifGen1<-dat_DifGen[INDKeep, ]
dat_DifGen2<-dat_DifGen[INDKeepPR, ]
# drop unused levels and reorder for plotting
dat_DifGen1$Genus<-droplevels(dat_DifGen1)$Genus
dat_DifGen1$Genus <- factor(dat_DifGen1$Genus, levels = c("Pseudarthrobacter", "Cetobacterium", "Catellicoccus",
                                              "Lactobacillus", "Streptococcus", "Clostridium_sensu_stricto_1", 
                                              "Romboutsia", "Terrisporobacter", 
                                              "Turicibacter","Pseudomonas"))

######################################
# Add significance 
######################################
INDKeep<-c()
for(i in 1:nrow(GeneraMain)){
  IND <- which(pwc2$Genus==GeneraMain$Genus[i])
  if (!isEmpty(IND)) {
    INDKeep<-rbind(INDKeep,IND)
  }
}

INDKeep<-c(INDKeep)   
pwc2_main<-pwc2[INDKeep,]
pwc2_main$Genus<-droplevels(pwc2_main)$Genus
pwc2_main$Genus <- factor(pwc2_main$Genus, levels = c("Pseudarthrobacter", "Cetobacterium", "Catellicoccus",
                                                          "Lactobacillus", "Streptococcus", "Clostridium_sensu_stricto_1", 
                                                          "Romboutsia", "Terrisporobacter", 
                                                          "Turicibacter","Pseudomonas"))

# plot with stats (Wilcoxon's test)
ggboxplot(dat_DifGen1, x = "ControlGroupWeek", y = "LogAbundance", facet.by = "Genus",
          fill = "ControlGroupWeek", palette = Col, alpha = 1,width = 0.9, fatten = 4,ncol=5) +
  stat_pvalue_manual(pwc2_main, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
         size = 0.5, linetype = "solid"),axis.text.x=element_blank(),
        strip.background =element_rect(fill="#e8e8e8"),strip.text = element_text(size = 10))

# plot with no stats
ggboxplot(dat_DifGen1, x = "ControlGroupWeek", y = "LogAbundance", facet.by = "Genus",
          fill = "ControlGroupWeek", palette = Col, alpha = 1,width = 0.9, fatten = 4,ncol=5) +
  theme(panel.border = element_blank(),text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
          size = 0.5, linetype = "solid"),axis.text.x=element_blank(),
        strip.background =element_rect(fill="#e8e8e8",colour="#e8e8e8"),strip.text = element_text(size = 10))+
        ylim(0, 10)


##--------Add the abundance and make a table of different taxa-------------

TaxDif[,"Mean_abundance"] <- NA
TaxDif[,"Relative_prevelence"] <- NA
for(i in 1:nrow(TaxDif)){
  IND_AB <- which(aa_Gen$Genus==TaxDif$Genus[i])
  IND_PR <- which(prevdf_genera$Genus==TaxDif$Genus[i])
  if (!isEmpty(IND_AB)) {
    TaxDif$Mean_abundance[i]<-aa_Gen$mean.abundance[IND_AB]
    TaxDif$Relative_prevelence[i]<-prevdf_genera$Relative.Prevelence[IND_PR]
  }
}
write.csv(TaxDif, file = "ANCOMCutOff0.85.csv")
