#### INFORMATION ABOUT THE DATA AND SCRIPT ####
### This script contains the codes for generating Figure 2 in the main article and Figure S1, Figure S2, and Figure S3 in the supporting information of the publication: Krishnan, A. and Osuri, A.M. (2022), Beyond the passive-active dichotomy: aligning research with the intervention continuum framework of ecological restoration. Restor Ecol. Accepted Author Manuscript e13828. https://doi.org/10.1111/rec.13828

# Authors: Aparna Krishnan and Anand Osuri*
# *aosuri@ncf-india.org

# The original data set on forest restoration used for this analysis is from:
# Crouzeilles R, Ferreira MS, Curran M (2016) Forest restoration: a global data set for biodiversity and vegetation structure. Ecology 97:2167. https://doi.org/ 10.1002/ecy.1474

# For more information on the methods please refer to Appendix 1 in the supporting information of the publication: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Frec.13828&file=rec13828-sup-0001-supinfo.docx

#### LOADING DATA AND CLEANING ####

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

#### Data from Crouzeilles et al. (2016) ###
# downloaded from https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.1474
forstdat <- read.csv(file = "Forest_restoration.csv", header = T, stringsAsFactors = F)

# Cleaning the Crouzeilles et al. 2016 data sheet

#subseting relevant columns from the restoration data frame
forstdat_sub <- forstdat[,c("ID", "Study", "Site", "Reference", "Restoration_activity", "Past_disturbance_type", "Time_restored", "Cover_10",  "Vegetation_structure", "Taxonomic_group", "Ecological_metric", "Restored", "Reference.1")]

# Selecting only active and passive restoration studies
forstdat_sub <- forstdat_sub %>% filter(forstdat_sub$Restoration_activity == "passive"| forstdat_sub$Restoration_activity == "active_managed" | forstdat_sub$Restoration_activity == "active_plantation") 

forstdat_sub$PA_activity <- ifelse(forstdat_sub$Restoration_activity == "active_plantation" | forstdat_sub$Restoration_activity == "active_managed",  "active" , "passive")

# removing all rows with restoration activity "-" and  with reference vaue 0
forstdat_sub <- forstdat_sub %>% filter(Reference.1 != 0.0, Restored != "-")

forstdat_sub$Reference.1 <- as.numeric(forstdat_sub$Reference.1) 
forstdat_sub$Restored <- as.numeric(forstdat_sub$Restored)

# Calculating Response Ratios
forstdat_sub$RR <- log(forstdat_sub$Restored/forstdat_sub$Reference.1)
forstdat_sub <- forstdat_sub[-which(forstdat_sub$RR == -Inf),]

#### Loading data collated for this manuscript ###
mot <- read.csv(file = "Restoration_Methods_Goals.csv", stringsAsFactors = F, header = T)

# Information on the data collated for this manuscript: We revisited the active restoration studies in the Crouzeilles et al. (2016) data set and ascertained details on the goals for restoration and on methods used and reported. This data set contains the goals, methods, and reporting details associated with each row of the Crouzeilles et al. (2016) data set.

# Merging the two data sheets
dat.mot <- merge(forstdat_sub, mot, by = c("ID", "Study", "Site", "Reference"), all.x = T)

#### DATA SUMMARIES ####

# Total number of comparisons - 3352
length((dat.mot$ID))

# Total number of studies - 235
length(unique(dat.mot$Reference))

# Total number of sites - 204
length(unique(dat.mot$Site))

# Number of active and passive restoration studies 
# Active - Studies 73, Rows 919, Sites 72
# Passive - Studies 188, Rows 2433, Sites 159
dat.mot %>% group_by(PA_activity) %>% summarise(
  No_Studies = length(unique(Reference)),
  No_rows = length(Restored),
  No_Sites = length(unique(Site)))


#### MAIN ARTICLE ####
#### FIGURE 2 ####

# The workflow and code for Figure 2 was modified from:
# Crouzeilles R, Ferreira MS, Chazdon RL, Lindenmayer DB, Sansevero JBB, Monteiro L, Iribarrem A, Latawiec AE, Strassburg BBN (2017) Ecological restoration success is higher for natural regeneration than for active restora- tion in tropical forests. Science Advances 3:e1701345. https://doi.org/10. 1126/sciadv.1701345

## A function which samples one reponse per study landscape
randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}

# Only active ecological (ER) and active economic (EP) studies
tree_active <- dat.mot %>% filter(Restoration_activity == "active_plantation" )
tree_active$Goal <- as.character(tree_active$Goal)
tree_active <- tree_active %>% filter(Goal == "EP" | Goal == "ER")

# Sample Sizes of biodiversity indicators
ss <- tree_active %>% group_by(Taxonomic_group, Goal, Ecological_metric) %>% summarise(
  nosStudies = length(unique(Reference)),
  nosSite = length(unique(Site)),
  nosrows = length(Goal))

bp_ss <- ss %>% filter(Taxonomic_group == "Birds" | Taxonomic_group == "Plants") %>% filter(Ecological_metric == "richness" | Ecological_metric == "abundance")
sum(bp_ss$nosrows)

#### Plants Richness ###

plants_ric <- subset(tree_active, tree_active$Taxonomic_group == "Plants" & tree_active$Ecological_metric == "richness")

Plants<-plants_ric

b.Plants_ER <- c()
b.Plants_EP <- c()

for (i in 1:2000){
  bsample <- ddply(Plants,.(Site),randomRows,1)
  t.stat <-as.numeric(summary(aov(bsample$Time_restored ~ bsample$Goal))[[1]][1,5])
  c.stat <- as.numeric(summary(aov(bsample$Cover_10 ~ bsample$Goal))[[1]][1,5])
  chisquare<-ftable(bsample$Past_disturbance ~ bsample$Goal)
  p.stat <- chisq.test(chisquare)[[3]]
  
  if ((t.stat > 0.05) & (c.stat > 0.05) & (p.stat > 0.05)){
    Plants_ER<-subset(bsample, bsample$Goal=="ER") # This is modified for ER vs. EP
    bestimate_ER <- mean(Plants_ER$RR)
    b.Plants_ER <- c(b.Plants_ER,bestimate_ER)
    
    Plants_EP<-subset(bsample, bsample$Goal=="EP")
    bestimate_EP <- mean(Plants_EP$RR)
    b.Plants_EP <- c(b.Plants_EP,bestimate_EP)}}

overall.Plants.Ric <- data.frame(Taxa = rep("Plants", 1000), Metric = rep("Richness", 1000), EP = b.Plants_EP[1:1000], ER =  b.Plants_ER[1:1000])

#### Plants Abundance ###

plants_abd <- subset(tree_active, tree_active$Taxonomic_group == "Plants" & tree_active$Ecological_metric == "abundance")

Plants<- plants_abd

b.Plants_ER <- c()
b.Plants_EP <- c()

for (i in 1:2000){
  bsample <- ddply(Plants,.(Site),randomRows,1)
  t.stat <-as.numeric(summary(aov(bsample$Time_restored ~ bsample$Goal))[[1]][1,5])
  c.stat <- as.numeric(summary(aov(bsample$Cover_10 ~ bsample$Goal))[[1]][1,5])
  chisquare<-ftable(bsample$Past_disturbance ~ bsample$Goal)
  p.stat <- chisq.test(chisquare)[[3]]
  
  if ((t.stat > 0.05) & (c.stat > 0.05) & (p.stat > 0.05)){
    Plants_ER<-subset(bsample, bsample$Goal=="ER") # This is modified for ER vs. EP
    bestimate_ER <- mean(Plants_ER$RR)
    b.Plants_ER <- c(b.Plants_ER,bestimate_ER)
    
    Plants_EP<-subset(bsample, bsample$Goal=="EP")
    bestimate_EP <- mean(Plants_EP$RR)
    b.Plants_EP <- c(b.Plants_EP,bestimate_EP)}}

overall.Plants.Abd <- data.frame(Taxa = rep("Plants", 1000), Metric = rep("Abundance", 1000),  EP = b.Plants_EP[1:1000], ER =  b.Plants_ER[1:1000])


#### Birds Richness ###

Birds_ric <- subset(tree_active, tree_active$Taxonomic_group == "Birds" & tree_active$Ecological_metric == "richness")

Birds<-Birds_ric

b.Birds_ER <- c()
b.Birds_EP <- c()

for (i in 1:2000){
  bsample <- ddply(Birds,.(Site),randomRows,1)
  t.stat <-as.numeric(summary(aov(bsample$Time_restored ~ bsample$Goal))[[1]][1,5])
  c.stat <- as.numeric(summary(aov(bsample$Cover_10 ~ bsample$Goal))[[1]][1,5])
  chisquare<-ftable(bsample$Past_disturbance ~ bsample$Goal)
  p.stat <- chisq.test(chisquare)[[3]]
  
  if ((t.stat > 0.05) & (c.stat > 0.05) & (p.stat > 0.05)){
    Birds_ER<-subset(bsample, bsample$Goal=="ER") # This is modified for ER vs. EP
    bestimate_ER <- mean(Birds_ER$RR)
    b.Birds_ER <- c(b.Birds_ER,bestimate_ER)
    
    Birds_EP<-subset(bsample, bsample$Goal=="EP")
    bestimate_EP <- mean(Birds_EP$RR)
    b.Birds_EP <- c(b.Birds_EP,bestimate_EP)}}

overall.Birds.Ric <- data.frame(Taxa = rep("Birds", 1000), Metric = rep("Richness", 1000), EP = b.Birds_EP[1:1000], ER =  b.Birds_ER[1:1000])

#### Birds Abundance ###

Birds_abd <- subset(tree_active, tree_active$Taxonomic_group == "Birds" & tree_active$Ecological_metric == "abundance")

Birds<- Birds_abd

b.Birds_ER <- c()
b.Birds_EP <- c()

for (i in 1:2000){
  bsample <- ddply(Birds,.(Site),randomRows,1)
  t.stat <-as.numeric(summary(aov(bsample$Time_restored ~ bsample$Goal))[[1]][1,5])
  c.stat <- as.numeric(summary(aov(bsample$Cover_10 ~ bsample$Goal))[[1]][1,5])
  chisquare<-ftable(bsample$Past_disturbance ~ bsample$Goal)
  p.stat <- chisq.test(chisquare)[[3]]
  
  if ((t.stat > 0.05) & (c.stat > 0.05) & (p.stat > 0.05)){
    Birds_ER<-subset(bsample, bsample$Goal=="ER") # This is modified for ER vs. EP
    bestimate_ER <- mean(Birds_ER$RR)
    b.Birds_ER <- c(b.Birds_ER,bestimate_ER)
    
    Birds_EP<-subset(bsample, bsample$Goal=="EP")
    bestimate_EP <- mean(Birds_EP$RR)
    b.Birds_EP <- c(b.Birds_EP,bestimate_EP)}}

overall.Birds.Abd <- data.frame(Taxa = rep("Birds", 1000), Metric = rep("Abundance", 1000),  EP = b.Birds_EP[1:1000], ER =  b.Birds_ER[1:1000])


#### Plotting Plants and Birds ###
All_Birds_Plants <- rbind(overall.Plants.Ric, overall.Plants.Abd, overall.Birds.Ric, overall.Birds.Abd)
All_Birds_Plants_long <- pivot_longer(All_Birds_Plants, cols = c(3,4), names_to = "Goals", values_to = "Outcome")
All_Birds_Plants_long$taxa_metric <- paste(All_Birds_Plants_long$Taxa, All_Birds_Plants_long$Metric, sep = " ")

richness_plant_bird <- All_Birds_Plants_long %>% filter(taxa_metric == "Plants Richness" | taxa_metric == "Birds Richness")
abundance_plant_bird <- All_Birds_Plants_long %>% filter(taxa_metric == "Plants Abundance" | taxa_metric == "Birds Abundance")

# Richness
richness <- ggplot(richness_plant_bird, aes(x = Goals, y = Outcome))+
  geom_abline(slope = 0, intercept = 0,linetype = "dashed",  col = "grey60")+
  geom_boxplot(width = 0.5, outlier.shape = NA,
               position = position_dodge(0.0002),
               notch = T,
               coef = 0,
               aes(fill = Goals))+
  ylab("Richness LRR")+
  labs(fill = "")+
  scale_fill_manual(label = c( "Economic", "Ecological"), 
                    values=c("#E67E22", "#2E86C1"))+
  scale_y_continuous(n.breaks = 5)+
  theme(text=element_text(size=14),
        strip.background = element_rect(fill=NA),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA, color = "white", size=0.2),
        strip.placement = "outside",
        legend.position= "none",
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10),
        axis.ticks.x = element_blank())+
  scale_x_discrete(name = c(""), labels = c("", ""))+
  facet_wrap(.~reorder(Taxa, rep(c(1,1,2,2), each = 1000 )), ncol = 2)

# Abundance
abundance <- ggplot(abundance_plant_bird, aes(x = Goals, y = Outcome))+
  geom_abline(slope = 0, intercept = 0,linetype = "dashed",  col = "grey60")+
  geom_boxplot(width = 0.5, outlier.shape = NA,
               position = position_dodge(0.0002),
               notch = T,
               coef = 0,
               aes(fill = Goals))+
  labs(fill = "")+
  ylab("Abundance LRR")+
  scale_fill_manual(label = c( "Economic", "Ecological"), 
                    values=c("#E67E22", "#2E86C1"))+
  theme(text=element_text(size=12),
        strip.background = element_rect(fill=NA),
        strip.text = element_blank(),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA, color = "white", size=0.2),
        strip.placement = "outside",
        legend.position= "none",
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10),
        axis.ticks.x = element_blank())+
  scale_x_discrete(name = c("Goal"), labels = c("Economic", "Ecological"))+
  scale_y_continuous(n.breaks = 5, limits = c(-1, 1))+
  facet_wrap(.~reorder(Taxa, rep(c(1,1,2,2), each = 1000 )), ncol = 2)

Fig2 <- richness / abundance
Fig2

#ggsave(filename = "Fig2",  width = 4.5, height = 4.5, device='tiff', dpi=700)

#### SUPPORTING INFORMATION ####
#### FIGURE S1 ####

# Summarizing intervention comparisons for active studies from the 'mot' data sheet
active <- mot %>% group_by(Reference) %>% summarise(
  rest = "active",
  int_comp = unique(Interventions_Compared),
  type = unique(Passive_Active)) %>% 
  filter(!is.na(int_comp))

# Summarizing intervention comparisons for passive studies from the 'forstdat_sub' data sheet
passive <- forstdat_sub %>% filter(Restoration_activity =="passive" ) %>% group_by(Reference) %>% summarise(rest_act = "passive")

# Merging the active and passive studies
comparisons <- merge(active, passive, by = "Reference", all.y = T, all.x = T)

# All passive studies which are not compared to active have "One" comparison
comparisons$int_comp[which(is.na(comparisons$int_comp) & comparisons$rest_act == "passive")] <- "One"

# Total studies = 235
tot_studies <- length(comparisons$Reference)

# summarizing the proportion of studies in each intervention comparison category
fig3b <- comparisons %>% group_by(int_comp) %>% summarise(nos = length(int_comp),
                                                          prop = 100*length(int_comp)/tot_studies)

ggplot(fig3b, aes(x = reorder(int_comp, c(1,3,2)), y = nos))+
  geom_bar(stat = "identity", fill = "#2E86C1", width = 0.7)+
  theme_classic(base_size = 17)+
  xlab("Restoration method comprisons")+ ylab("Number of studies")+
  theme(axis.line = element_line(size = 0.7, colour = "black"),
        legend.position = "none",
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0,0.1), limits = c(0,220), n.breaks = 7)+
  scale_x_discrete(labels = c("One", "Two", "   Three or more"))

# ggsave(filename = "Figure2b",  width = 6, height = 5.4, device='tiff', dpi=700)

#### FIGURE S2 ####
active_goals <- dat.mot %>% filter(PA_activity == "active" )

# ER,EP is when interventions/goals associated with the rows are unclear, hence dropped from analysis. This was the case primarily with Armstrong_&_Nichols_2000 (and a few row IDs in Kanowski_et_al_2003).
active_goals <- active_goals[-which(active_goals$Goal == "ER,EP"),] 

interventions <- active_goals %>% group_by(Reference, Restoration_activity, Goal) %>% summarise(nos = unique(No_Active_Interventions))

# Total number of interventions
tot.interventions <- sum(interventions$nos, na.rm = T)

FigS2 <- interventions %>% group_by(Goal, Restoration_activity) %>% summarise(
  prop = 100*sum(nos, na.rm = T)/tot.interventions,
  totals = sum(nos, na.rm = T))

# Making NA as Unknown
FigS2$Goal[which(is.na(FigS2$Goal))] <- "Unknown"

ggplot(FigS2, aes(x = reorder(Goal, c(2,2,1,1,3,4,4)), y = totals))+
  geom_bar(stat = "identity", aes(fill = Restoration_activity), width = 0.7)+
  theme_classic(base_size = 17)+
  xlab("Restoration goal")+ ylab("Number of treatments")+
  theme(legend.position = c(0.8, 0.8),
        axis.line = element_line(size = 0.7, colour = "black"),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0,0.1), limits = c(0,60), n.breaks = 7)+
  scale_fill_manual(name = "", labels = c("Active management", "Active plantation"), values = c("#A9CCE3", "#2E86C1"))+
  scale_x_discrete(labels = c( "Ecological","Economic", "Ecological \\nand Economic", "Unknown"))

# ggsave(filename = "FigureS2", width = 6, height = 5.4, device='tiff', dpi=700)

# Only tree planting studies summary
tree_planting_goals <- FigS2 %>% filter(Restoration_activity == "active_plantation") %>% select("Goal", "Restoration_activity", "totals")

# total number of tree planting cases 678
all_tp <- active_goals %>% filter(Restoration_activity=="active_plantation")
length(all_tp$ID)

# total number of tree planting interventions - 73
tot_tp <- sum(tree_planting_goals$totals)

# proportion of interventions with different goals
tree_planting_goals$pct_goals <- 100*tree_planting_goals$totals/tot_tp
tree_planting_goals

#### FIGURE S3 ####
# Reporting data sheet
dat <- read.csv(file = "Restoration_Methods_Goals.csv", header = T, stringsAsFactors = F)

# subsetting ecological restoration studies
report <- dat[,c("Reference", "Site_Preparation_Reported", "Main_Intervention_Reported", "Maintenance_Reported")]
report <- na.omit(report)
colnames(report) <- c("Reference","Site_Preparation", "Main_Intervention", "Maintenance")
tot.studies<- length(report$Reference)

# Number and percentage of studies reporting 
report_nos <- pivot_longer(report, cols = 2:4, names_to = "Method", values_to = "Reporting")
report_nos %>% group_by(Method) %>% summarise(
  No_Studies_Not_Reported = sum(Reporting == "N", na.rm = T),
  Pct_Studies_Not_Reported = 100*sum(Reporting == "N", na.rm = T)/tot.studies)

FigS3 <- report_nos %>% group_by(Method, Reporting) %>% summarise(
  precent = 100*sum(length(Reporting))/tot.studies,
  nums = length(Reporting))

ggplot(FigS3, aes(x = reorder(Method, c(2,3,3,1,1)), y = precent))+
  geom_bar(stat = "identity",  aes(fill = Reporting),colour = "grey30",size = 0.2)+
  xlab("")+ ylab("Percent studies")+
  scale_fill_manual(label = c( "Not Reported", "Reported"), 
                    values=c("white", "grey30"))+
  theme_classic(base_size = 17)+
  theme(axis.line = element_line(size = 0.7, colour = "black"),
        legend.title = element_blank(),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0,0.1), limits = c(0,100), n.breaks = 7)+
  scale_x_discrete(labels = c("Site\\npreparation", "Main\\nintervention", "Maintenance"))

# ggsave(filename = "FigureS3", width = 6.3, height = 5.4, device='tiff', dpi=700)