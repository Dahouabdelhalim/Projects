### clear working environment
rm(list = ls())

### open necessary libraries
library(car)
library(tidyr)
library(ggplot2)
library(dplyr)
library(wesanderson)
library(mosaic)
library(cowplot)
library(lme4)
library(multcompView)
library(multcomp)
library(ggrepel)
library(grid)
library(dglm)
library(vegan)
library(tidyverse)
library(lsmeans)
library(magrittr)
library(mctoolsr)
library(broom)

### set contrasts
options(contrasts=c('contr.sum', 'contr.poly'))

### read in siderophore data
SID<-read.csv("siderophore.csv", header=TRUE)
attach(SID)
names(SID)

### normalize siderophore data (CAS_cult) by setting minimum value to zero
SID <- SID %>% mutate(CAS_norm = CAS_cult+(-1*min(CAS_cult, na.rm=T))) %>% filter(Time!="3")

# -------------------------- #
# Code accompanying Figure 1 #
# -------------------------- #

### select final time point
SID_T6 <- SID %>% filter (Time==6) %>%
  dplyr::mutate(., Culture = factor(Culture, levels = c("Mono","Co"))) %>% 
  dplyr::mutate(., Sample = factor(Sample, levels = c("NMC","Wt"))) %>% 
  dplyr::mutate(., Cu = factor(Cu, levels = c("0","1")))

### select compost community data (NMC)
SID_NMC <- SID_T6 %>% filter(Sample=="NMC") 

### select data without focal species P. fluorescens
SID_NMC_mono <- SID_NMC %>% filter(Culture =="Mono")

# ----- #
# STATS #
# ----- #

### test whether variance differs between copper treatments in absence of P. fluorescens
model_disp_mono <- dglm(CAS_norm ~ Cu, ~Cu, data=SID_NMC_mono, family=gaussian)

summary(model_disp_mono)
anova(model_disp_mono)

# ---- #
# PLOT #
# ---- #

### set up colour scheme 
control <- wes_palette(n=5, name="Darjeeling1")
control_col <- control[5]

copper <- wes_palette(n=1, name="Darjeeling1")
copper_col <- copper[1]

# ------- #
# Panel A #
# ------- #

Fig_1A <- ggplot(SID_NMC_mono, aes(x=CAS_norm, colour=Cu, fill=Cu)) +
  scale_fill_manual(values=c(control_col,copper_col), name="Treatment", breaks=c("0","1"), labels=c("Control", "Copper")) +
  geom_density(alpha=0.8, adjust=3, size=1) +
  scale_color_manual(values=c(control_col,copper_col), name="Treatment", breaks=c("0","1"), labels=c("Control", "Copper")) +
  geom_rug(size=0.8) +
  theme_bw()

Fig_1A <- Fig_1A + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_1A <- Fig_1A + xlab("\\nSiderophore production") 
Fig_1A <- Fig_1A + ylab("Frequency\\n")

Fig_1A <- Fig_1A + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_1A <- Fig_1A + theme(axis.title.x = element_text(size=16))
Fig_1A <- Fig_1A + theme(axis.text.x = element_text(size=14))
Fig_1A <- Fig_1A + theme(axis.text.y = element_text(size=14))

Fig_1A <- Fig_1A + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_1A <- Fig_1A + theme(legend.position="none")
Fig_1A <- Fig_1A + theme(legend.title = element_text(size=14))
Fig_1A <- Fig_1A + theme(legend.text =  element_text(size=12))

Fig_1A <- Fig_1A + theme(axis.title.y = element_text(size=16))

Fig_1A


# ------- #
# Panel B #
# ------- #

### Read in data
reg <- read.csv("upregulation.csv", header=T)
attach(reg)
names(reg)

# ----- #
# STATS #
# ----- #

# Normalize data by setting minimum value to zero
reg <- reg %>% mutate(cas_norm = CAS+0.53)
reg$Cu<-factor(reg$Cu)
reg$Genus<-plyr::revalue(reg$Genus, c("Arthrobacter"="A", "Bacillus"="B", "Microbacterium"= "M", "Pseudomonas" ="P", "Stenotrophomonas"="S", "Cupriavidus"="C"))
reg$Genus<-factor(reg$Genus, levels=c("M", "C","A","P","S", "B"))

# 2-way ANOVA
model <- glm(cas_norm ~ Genus * Cu, data=reg)
summary(model)
anova(model, test="F")

model_b <- glm(cas_norm ~ Genus + Cu, data=reg)
summary(model_b)
anova(model, model_b, test="F")

# ---- #
# PLOT #
# ---- #

# create data frame
df<-data.frame(reg$Cu, reg$Genus, reg$cas_norm)
inter<-interaction(df$reg.Genus,df$reg.Cu)

# set up colours
my_palette<-c(rep("red", 6), rep("skyblue", 6))

# ggplot
Fig_1B <- ggplot(df, aes(as.factor(reg.Genus), reg.cas_norm, fill=inter))
Fig_1B <- Fig_1B + geom_boxplot(position="dodge", size=1)
Fig_1B <- Fig_1B + scale_x_discrete(position="top")

Fig_1B <- Fig_1B + coord_flip()
Fig_1B <- Fig_1B + scale_colour_manual(values=my_palette)

Fig_1B <- Fig_1B + scale_fill_manual(values=my_palette)

Fig_1B <- Fig_1B + scale_alpha_manual(values=c(1,1))

Fig_1B <- Fig_1B + ylab("Siderophore production")

Fig_1B <- Fig_1B + xlab("")

Fig_1B <- Fig_1B + theme(axis.line = element_line(colour = "black", size=1), panel.grid.major = element_blank(), panel.border = element_rect(colour='black', size=1, linetype="solid"))

Fig_1B <- Fig_1B + theme(axis.title.x = element_text(size=16))
Fig_1B <- Fig_1B + theme(axis.text.x = element_text(size=14))

Fig_1B <- Fig_1B + theme(axis.line.x=element_blank())
Fig_1B <- Fig_1B + theme(axis.line.y=element_blank())

Fig_1B <- Fig_1B + theme(axis.text.y = element_text(face="italic"))
Fig_1B <- Fig_1B + theme(axis.text.y = element_text(size=14))
Fig_1B <- Fig_1B + theme(legend.position="none")

Fig_1B <- Fig_1B + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
Fig_1B <- Fig_1B + theme(axis.title.x=element_text(margin=margin(20,0,0,0)))

Fig_1B <- Fig_1B + scale_x_discrete(position="top", labels=c("B" = "Bacillus", "S" = "Stenotrophomonas", "P" = "Pseudomonas", "A"="Arthrobacter", "C" = "Cupriavidus", "M"="Microbacterium"))

Fig_1B

# combine panels and save as pdf <- 12 x 6.85 inches
Fig.1 <- ggdraw(xlim=c(-0.1,3.1),ylim=c(-0.1,1.6)) +
  draw_plot(Fig_1A, 0, 0, 1.45, 1.5) +
  draw_plot(Fig_1B, 1.6, -0.006, 1.5, 1.508) +
  draw_plot_label(c("A", "B"), c(0.156, 1.596), c(1.586, 1.586), size = 24)

Fig.1

# -------------------------- #
# Code accompanying Figure 2 #
# -------------------------- #

### read in pairwise competition data compost isolates
competition <- read.csv("Comp_densities.csv", header=T)
attach(competition)
names(competition)

### normalize data and set minimum observed threshold
competition <- competition %>% mutate(cas_kb_norm = cas_kb_cfu + (-1*min(cas_kb_cfu))) %>%
  mutate(., replicate=as.factor(replicate)) %>%
  mutate(., treatment=as.factor(treatment)) %>%
  mutate(., m_growth_threshold = ifelse(is.na(m) == TRUE, min(m, na.rm=T), m))

### calculate mean for each pairwise interaction
competition_new <-  competition %>% group_by(pair, treatment, species) %>% 
  dplyr::summarise(m=mean(m_growth_threshold, na.rm=T)) %>% ungroup()

mean_competition <- competition_new %>% dplyr::mutate(row=1:n()) %>% group_by(row) %>% 
  dplyr::mutate(species1=gsub(as.character(species),'', as.character(pair)))

mean_competition$treatment <- recode_factor(mean_competition$treatment, `C` = "control", `CU` = "copper")

### re-order species levels based on levels of siderophore production (Table S1)
mean_competition <- mean_competition %>% mutate(new_species = fct_relevel(species, "A", "C", "E", "F","J","I","G","H","K","D"))

### read in monoculture data
mono<- read.csv("mono_compost_growth.csv", header=T)
attach(mono)
names(mono)

mono <- mono %>% mutate(., m_threshold = ifelse(is.na(m) == TRUE, min(m, na.rm=T), m))
mono_mean <- group_by(mono, species, treatment) %>% dplyr::summarise(mean_m = mean(m_threshold)) %>% ungroup()

mono_combine <- merge(mono,mono_mean, all.x=TRUE)

# define pairs
pairs <- data.frame(crossing(mono_combine, species, species))
pairs<-pairs[!(pairs$species==pairs$species1),]
pairs$combine <- paste(pairs$species,pairs$species1, sep="_")

mono_tidy <- merge(mono_combine, pairs, all.x=TRUE)

mean_comps <- dplyr::select(mono_tidy, species, treatment, mean_m) %>% 
  distinct() %>%
  dplyr::rename(., species1 = species,
                mean_m_spp1 = mean_m)

mono_tidy2 <- merge(mono_tidy, mean_comps, by = c('species1', 'treatment'), all.x = TRUE) %>%
  arrange(species) %>%
  mutate(sel_coef=m_threshold-mean_m_spp1)

# create data frame with means for each unique treatment combination
# reorder levels of species

mono_tidy2 <- mono_tidy2 %>% mutate(new_species = fct_relevel(species, "A", "C", "E", "F","J","I","G","H","K","D"))

means_species <-  mono_tidy2 %>% group_by(new_species, treatment) %>% dplyr::summarise(m=mean(m_threshold, na.rm=T)) %>% 
  ungroup()

means_species$treatment <- recode_factor(means_species$treatment, `C` = "control", `CU` = "copper")

# calculate m competition - m mono for each species combination
# create a variable summarizing mean growth rate for each species growing in monoculture

mean_mono <- means_species %>% dplyr::rename(mean_mono_m=m)

growth_difference <- merge(mean_competition, mean_mono, by = c('new_species', 'treatment'), all.x = TRUE) %>% 
  arrange(new_species) %>%
  dplyr::mutate(difference=m-mean_mono_m) %>%
  arrange(pair)

# calculate mean cas assay for each species
mean_kb <- competition %>% group_by(species) %>% dplyr::summarise(cas_kb=mean(cas_kb_cfu, na.rm=T)) %>% ungroup()
data_cas <- mean_kb %>% mutate(cas_norm = cas_kb+(-1*min(cas_kb)))

# set up colour palette
col1 <- "Gold"
col2 <- wes_palette(n=1, name="Moonrise3")
col3 <- wes_palette(n=1, name="FantasticFox1")
col4 <- wes_palette(n=5, name="FantasticFox1")

col_A <- colorRampPalette(c(col1, col2))(6)
col_B <- colorRampPalette(c(col3, col4[5]))(4)
col_C <- c(col_A, col_B)

# label names
strain_names <- c(
  `A` = "Arthrobacter",
  `C` = "Devosia",
  `E` = "Staphylococcus",
  `F` = "Oerskovia",
  `J` = "Flavobacterium",
  `I` = "Bordetella",
  `G` = "Rhodococcus",
  `H` = "Pedobacter",
  `K` = "Variovorax",
  `D` = "Brevundimonas")

# plot difference in Malthusian growth rate in co-culture between copper and control and plot against siderophore production
competition_co_summary <- merge(mean_competition, data_cas, by = 'species', all.x = TRUE)

# calculate difference
competition_co_difference <- competition_co_summary %>% dplyr::select(new_species, species1, treatment, m, cas_norm) %>%
  group_by(new_species, species1) %>%
  spread(., treatment, m) %>%
  dplyr::mutate(difference_co = copper-control) %>%
  mutate(species1 =fct_relevel(species1, "A", "C", "E", "F","J","I","G","H","K","D"))

# subset for copper only
copper_mono <- mean_mono %>% filter(treatment=="copper")
copper_co <- competition_co_summary %>% filter(treatment=="copper")

copper_total <- merge(copper_co, copper_mono, by = 'new_species', all.x = TRUE) %>% mutate(species1 =fct_relevel(species1, "A", "C", "E", "F","J","I","G","H","K","D"))

# culculate difference between co and monoculture for copper treatment
copper_total <- copper_total %>% mutate(difference = m-mean_mono_m)

# ------- #
# PANEL A #
# ------- #

Fig_2A <- ggplot(data=copper_total, aes(x=cas_norm, y=difference, colour=species1))
Fig_2A  <- Fig_2A  + geom_abline(intercept = 0, slope = 0, colour="grey")

Fig_2A  <- Fig_2A + geom_smooth(color="black", method="lm", se=TRUE, linetype="solid", fullrange=TRUE, size=0.5, alpha=0.1)

Fig_2A  <- Fig_2A + ylab(expression(italic(m) [competition] * " - " * italic(m) [monoculture]))

Fig_2A  <- Fig_2A + xlab("Per capita siderophore production\\n")

Fig_2A  <- Fig_2A + theme_bw()
Fig_2A  <- Fig_2A  + theme(axis.line = element_line(colour="black"))
Fig_2A  <- Fig_2A + theme(axis.title.x = element_text(size=16))
Fig_2A  <- Fig_2A  + theme(axis.title.y = element_text(size=18))
Fig_2A  <- Fig_2A  + theme(axis.text.x = element_text(size=14))

Fig_2A  <- Fig_2A  + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          strip.background = element_blank(),
                          panel.border = element_rect(colour = "black"))

Fig_2A  <- Fig_2A  + scale_x_continuous(expand=c(0,0))
Fig_2A  <- Fig_2A  + scale_y_continuous()
Fig_2A  <- Fig_2A  + coord_cartesian(xlim=c(-0.001, 0.0275), ylim=c(-3,2))


Fig_2A  <- Fig_2A + theme(axis.line.x=element_blank())
Fig_2A  <- Fig_2A + theme(axis.text.y = element_text(size=14))
Fig_2A  <- Fig_2A + theme(legend.position = "none")

Fig_2A  <- Fig_2A + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
Fig_2A  <- Fig_2A  + theme(axis.title.x=element_text(margin=margin(10,0,0,0)))

Fig_2A  <- Fig_2A + geom_jitter(colour="white", width=0.0001, alpha=0.8, size=3.2, shape=21,  stroke=1, aes(fill=species1))+scale_fill_manual(values=col_C, name="", labels=c("Arthrobacter", "Devosia", "Staphylococcus", "Oerskovia", "Flavobacterium", "Bordetella", "Rhodococcus", "Pedobacter", "Variovorax", "Brevundimonas"))

Fig_2A  <- Fig_2A + stat_summary(fun.y = mean, geom = "point", aes(group = new_species), shape=21, fill = col_C, size = 5, stroke=1)

Fig_2A 

# ----- # 
# STATS #
# ----- #

model <- lm(data=copper_total, difference~cas_norm)
summary(model)

model_b <- lm(data=copper_total, difference~1)
summary(model_b)

anova(model, model_b)
confint(model,level=0.95)


# ------- #
# PANEL B #
# ------- #

control_mono <- mean_mono %>% filter(treatment=="control")
control_co <- competition_co_summary %>% filter(treatment=="control")

control_total <- merge(control_co, control_mono, by = 'new_species', all.x = TRUE) %>% mutate(species1 =fct_relevel(species1, "A", "C", "E", "F","J","I","G","H","K","D"))

# culculate difference between co and monoculture
control_total <- control_total %>% mutate(difference = m-mean_mono_m)

Fig_2B <- ggplot(data=control_total, aes(x=cas_norm, y=difference, colour=species1))
Fig_2B <- Fig_2B + geom_abline(intercept = 0, slope = 0, colour="grey")

Fig_2B <- Fig_2B + geom_smooth(color="black", method="lm", se=TRUE, linetype="solid", fullrange=TRUE, size=0.5, alpha=0.1)

Fig_2B <- Fig_2B + ylab("")

Fig_2B <- Fig_2B + xlab("Per capita siderophore production\\n")

Fig_2B <- Fig_2B + theme_bw()
Fig_2B <- Fig_2B  + theme(axis.line = element_line(colour="black"))
Fig_2B <- Fig_2B + theme(axis.title.x = element_text(size=16))
Fig_2B <- Fig_2B  + theme(axis.title.y = element_text(size=16))
Fig_2B <- Fig_2B  + theme(axis.text.x = element_text(size=14))

Fig_2B <- Fig_2B  + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          strip.background = element_blank(),
                          panel.border = element_rect(colour = "black"))

Fig_2B <- Fig_2B  + scale_x_continuous(expand=c(0,0))
Fig_2B <- Fig_2B  + scale_y_continuous()
Fig_2B <- Fig_2B  + coord_cartesian(xlim=c(-0.001, 0.0275), ylim=c(-3,2))

Fig_2B <- Fig_2B + theme(axis.line.x=element_blank())
Fig_2B <- Fig_2B + theme(axis.text.y = element_text(size=14))

Fig_2B <- Fig_2B + theme(legend.text=element_text(size=14))
Fig_2B <- Fig_2B + theme(legend.key.size = unit(1, "cm"))

Fig_2B <- Fig_2B + theme(axis.title.y=element_text(margin=margin(0,10,0,0)))
Fig_2B <- Fig_2B  + theme(axis.title.x=element_text(margin=margin(10,0,0,0)))

Fig_2B <- Fig_2B + geom_jitter(colour="white", width=0.0001, alpha=0.8, size=3.2, shape=21,  stroke=1, aes(fill=species1))+scale_fill_manual(values=col_C, name="", labels=c("Arthrobacter", "Devosia", "Staphylococcus", "Oerskovia", "Flavobacterium", "Bordetella", "Rhodococcus", "Pedobacter", "Variovorax", "Brevundimonas"))

Fig_2B <- Fig_2B + stat_summary(fun.y = mean, geom = "point", aes(group = new_species), shape=21, fill = col_C, size = 5, stroke=1)

Fig_2B

# ----- #
# STATS #
# ----- #

model_control <- lm(data=control_total, difference~cas_norm)
summary(model_control)

model_control_b <- lm(data=control_total, difference~1)
summary(model_control_b)

anova(model_control, model_control_b, test="F")

confint(model_control,level=0.95)

# Combine two panels
Fig.2 <- ggdraw(xlim=c(-0.1,3.5),ylim=c(-0.1,1.6)) +
  draw_plot(Fig_2A, 0, 0, 1.45, 1.5) +
  draw_plot(Fig_2B, 1.6, 0, 1.85, 1.5) +
  draw_plot_label(c("A", "B"), c(0.2, 1.74), c(1.586, 1.586), size = 24)

Fig.2

# -------------------------- #
# Code accompanying Figure 3 #
# -------------------------- #

SID_NMC_co <- SID_NMC %>% filter(Culture =="Co")

# ------- #
# STATS A #
# ------- #

# test whether mean levels of siderophore production differ between copper treatments when NMC is evolving in presence of P. fluorescens
model_lmer <- lmer(CAS_norm~Cu + (1|Rep), data=SID_NMC_co)
print(summary(model_lmer))

model_lmer_b <- lmer(CAS_norm ~ 1 + (1|Rep), data=SID_NMC_co)
anova(model_lmer, model_lmer_b, test = "LRT")

# ------- #
# PANEL A #
# ------- #

# create data frame with means per community/strain
sid_df <- SID %>%
  dplyr::group_by(Sample, Culture, Rep, Time, Cu) %>%
  dplyr::summarise(mean_sids=mean(CAS_norm, na.rm=T)) %>%
  ungroup()

# collapse copper treatment for T = 0
sid_df <- sid_df %>%
  mutate(., new = case_when(Time == 0 ~ "I",
                            TRUE ~ as.character(.$Cu))) %>%
  mutate(., new=as.factor(new)) %>%
  mutate(., new=relevel(new, ref="I")) %>%
  mutate(., new_cult = case_when(Time == 0 ~ "No",
                                 TRUE ~ as.character(.$Culture))) %>%
  mutate(., new_cult=as.factor(new_cult)) %>%
  mutate(., new_cult=relevel(new_cult, ref="No")) %>%
  mutate(., Sample=relevel(Sample, ref="Wt"))

# collapse copper treatment for T=0
SID_co <-   sid_df %>% filter(Sample=="NMC" & Culture=="Co")  %>%
  dplyr::mutate(., Cu = factor(Cu, levels = c("1","0"))) %>%
  dplyr::mutate(., Time = factor(Time, levels = c("0","6"))) %>%
  dplyr::mutate(., new = factor(new, levels = c("I", "0","1"))) %>%
  arrange(., Time) %>%
  dplyr::mutate(., treatment = as.factor(c(rep("I",12), rep("con",6), rep("cu",6)))) %>%
  dplyr::mutate(., treatment = factor(treatment, levels = c("I", "con", "cu"), labels = c("ancestral", "control", "copper")))

# set up colour scheme
control <- wes_palette(n=5, name="Darjeeling1")
control_col <- control[5]

copper <- wes_palette(n=1, name="Darjeeling1")
copper_col <- copper[1]

dodge = position_dodge(width=0.8)

Fig_3A <- ggplot(data=SID_co, aes(x=treatment, y=mean_sids))
Fig_3A <- Fig_3A + geom_boxplot(aes(fill=treatment), outlier.size = 0) + scale_fill_manual(values=c("black", control_col, copper_col))

dat <- ggplot_build(Fig_3A)$data[[1]]

xmin=dat$xmin[1]
xmax=dat$xmax[1]
middle=dat$middle[1]

Fig_3A <- Fig_3A + geom_segment(aes(x=xmin, y=middle, xend=xmax, yend=middle), col="white", size=0.5)

Fig_3A <- Fig_3A + coord_cartesian(ylim=c(0.4, 0.7)) + scale_y_continuous(breaks=seq(0.4, 0.7, 0.1))

Fig_3A <- Fig_3A + geom_jitter(width=0.15, size=3, shape=21,  stroke=1, col="black", fill="white") 

Fig_3A <- Fig_3A + theme_bw()

Fig_3A <- Fig_3A + xlab("\\nSelection regime") 

Fig_3A <- Fig_3A + ylab("Siderophore production\\n")

Fig_3A <- Fig_3A + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_3A <- Fig_3A + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=18)
) 

Fig_3A <- Fig_3A + theme(axis.title.x = element_text(size=16))
Fig_3A <- Fig_3A + theme(axis.text.x = element_text(size=14))
Fig_3A <- Fig_3A + theme(axis.text.y = element_text(size=14))
Fig_3A <- Fig_3A + theme(legend.position="none")
Fig_3A <- Fig_3A + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_3A <- Fig_3A + theme(axis.title.y = element_text(size=16))

Fig_3A

# add text to graph
stat_data <- data.frame(text = c("","a","b"),
                        y = dat$ymin-0.02,
                        x = dat$xid, stringsAsFactors = FALSE)

Fig_3A <- Fig_3A + geom_text(data = stat_data, aes(label = text, x = x, y = y), colour = "black", size=6, fontface = "bold")

Fig_3A

# ------- #
# STATS B #
# ------- #

# test whether variance differs between copper treatments in presence of P. fluorescens, using dglm
model_disp_co <- dglm(CAS_norm ~ Cu, ~Cu, data=SID_NMC_co, family=gaussian)
summary(model_disp_co)

anova(model_disp_co)

# ------- #
# PANEL B #
# ------- #

# set up colour scheme
control <- wes_palette(n=5, name="Darjeeling1")
control_col <- control[5]

copper <- wes_palette(n=1, name="Darjeeling1")
copper_col <- copper[1]

# plot
Fig_3B <- ggplot(SID_NMC_co, aes(x=CAS_norm, colour=Cu, fill=Cu)) +
  scale_fill_manual(values=c(control_col,copper_col), name="Treatment", breaks=c("0","1"), labels=c("Control", "Copper")) +
  geom_density(alpha=0.8, adjust=3, size=1) +
  scale_color_manual(values=c(control_col,copper_col), name="Treatment", breaks=c("0","1"), labels=c("Control", "Copper")) +
  geom_rug(size=0.8) +
  theme_bw()

Fig_3B <- Fig_3B + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_3B <- Fig_3B + xlab("\\nSiderophore production") 
Fig_3B <- Fig_3B + ylab("Frequency\\n")

Fig_3B <- Fig_3B + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_3B <- Fig_3B + theme(axis.title.x = element_text(size=16))
Fig_3B <- Fig_3B + theme(axis.text.x = element_text(size=14))
Fig_3B <- Fig_3B + theme(axis.text.y = element_text(size=14))

Fig_3B <- Fig_3B + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_3B <- Fig_3B + theme(legend.position="none")
Fig_3B <- Fig_3B + theme(legend.title = element_text(size=14))
Fig_3B <- Fig_3B + theme(legend.text =  element_text(size=12))

Fig_3B <- Fig_3B + theme(axis.title.y = element_text(size=16))

Fig_3B

# ------- #
# STATS C #
# ------- #

Bray <- read.csv("Bray_curtis.csv", header=T)
attach(Bray)
names(Bray)

# first 6 rows is control, last 6 rows is copper

d_vars <- data.frame(groups = as.character(c(rep('control',6), rep('copper',6))),
                     row = row.names(Bray))

# create a distance matrix - want samples as rows, species as columns
d_dist <- vegdist(Bray, method = 'bray')

# do permanova and test whether variance is equally distributed across treatments (YES)
mod1 <- adonis(d_dist ~ groups, data = d_vars, method = "bray")
anova(betadisper(d_dist, d_vars$groups))

# ------- #
# PANEL C #
# ------- #

# extract all the data using broom
# get PCA summary
d_PCA_sum <- broom::tidy(PCA, matrix = 'pcs')

# get PCA loadings
d_PCA_load <- broom::tidy(PCA, matrix = 'v') %>%
  mutate(., PC = paste('PC_', PC, sep = '')) %>%
  spread(., PC, value)

# subset for just top 15 loadings
# if you have a lot of genes this is advisable otherwise you will have so much text!
# I calculate the distance of each point in the PC1 and PC2 from 0,0 and subset those

# function to get distance from 00
dist_from_00 <- function(x, y){
  return(sqrt((0 - x)^2+(0-y)^2))
}

top_15_loadings <- mutate(d_PCA_load, distance = dist_from_00(PC_1, PC_2)) %>%
  top_n(., 15, distance)

# get PCA points
d_PCA_points <- broom::tidy(PCA, matrix = 'samples') %>%
  mutate(., row = as.character(row),
         PC = paste('PC_', PC, sep = '')) %>%
  merge(., d_vars, by = 'row') %>%
  spread(., PC, value)

# do adonis and multiple contrasts
# distance matrices
Euclid_matrix <- dist(Bray)

# make some character vectors factors
d_vars <- mutate_at(d_vars, c('groups'), as.factor)
row.names(d_vars) <- d_vars$row

# adonis
mod1_euclid <- adonis(Euclid_matrix ~ groups, data = d_vars)

# do pairwise contrasts - look at different ways of doing this! ####
# using mctoolsr from GitHub - this looks good
# NOTE - need to make sure that variable dataset has rownames equivalent to the distance matrix

euclid_contrasts <- mctoolsr::calc_pairwise_permanovas(Euclid_matrix, d_vars, 'groups')

# obviously a bit redundant when there are only two levels
# make a column of all unique combinations of factors and insert instead of "groups"

# betadisper() analysis and plotting in ggplot2 ####

# Looking at variance across groups ####

# can only do one factor would use every combination of id
mod1_dispers <- betadisper(Euclid_matrix, d_vars$groups)

# plot of model
plot(mod1_dispers)
boxplot(mod1_dispers)

# anova
anova(mod1_dispers)

# Permutation test for F
pmod <- permutest(mod1_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod1_dispers)

# get betadisper dataframes ####
# have written functions to grab the necessary data from the betadisper object

# functions ####
# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

# get betadisper data ####
betadisper_dat <- get_betadisper_data(mod1_dispers)

# do some transformations on the data
betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig))

# add convex hull points ####
# this could be put in a function
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(dplyr::select(betadisper_dat$centroids, group, PCoA1, PCoA2), dplyr::select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# Now the dataframes are all ready to be completely customisable in ggplot
# plot betadispersion plot

Fig_3C <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 8, alpha = 0.9) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$eigenvector, size = 2, alpha=1) +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull, alpha=0.5) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = group, alpha=0.5), betadisper_lines) +scale_color_manual('', values = c(copper_col, control_col), labels = c("Copper", 'Control'))

Fig_3C <- Fig_3C + theme_bw() + theme(legend.position = 'none')

Fig_3C <- Fig_3C + ylab("Axis 2 (14.7%)\\n")
Fig_3C <- Fig_3C + xlab("\\nAxis 1 (58.5%)")

Fig_3C <- Fig_3C + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_3C <- Fig_3C + theme(axis.title.x = element_text(size=16))
Fig_3C <- Fig_3C + theme(axis.text.x = element_text(size=14))
Fig_3C <- Fig_3C + theme(axis.text.y = element_text(size=14))

Fig_3C <- Fig_3C + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_3C <- Fig_3C + theme(axis.title.y = element_text(size=16))

Fig_3C 

# ------- #
# PANEL D #
# ------- #

#####################################################################################
# relative abundance of ten most common culturable genera in compost microbial      # 
# community after 6 weeks pf incubation in copper and control microcosms            #
# with high siderophore producing labstrain of Pseudomonas fluorescens              #
#####################################################################################

# relative abundance with colour for mean siderophores
library(plyr)
library(data.table)

source("https://raw.githubusercontent.com/janhove/janhove.github.io/master/RCode/sortLvls.R")

Gen<-read.csv("Genera_abundance.csv", header=T)
attach(Gen)
names(Gen)

Gen <- Gen %>% mutate(CAS_mean)

Gen_sub <- group_by(Gen, Treatment) %>%
  subset(., !Genus=="Unknown")

# subset to only have co-culture NMC!!!!!!

Gen_sum_co <- Gen_sub[Gen_sub$Treatment=="PF" | Gen_sub$Treatment=="B", ] %>%
  subset(., !Genus=="Unknown")

Gen_sum_co$Genus <- factor(Gen_sum_co$Genus)

# calculate and rank abundance
summed <- as.vector(tapply(Gen_sum_co$n, Gen_sum_co$Genus, sum))
gen_summed <- unique(Gen_sum_co$Genus)

abun <- data.frame(gen_summed, summed)
abun2 <- dplyr::mutate(abun, rank = rank(desc(summed), ties.method = 'first'))

colnames(abun2) <- c("Genus", "Sum", "Rank")

# merge two data frames  
df <- merge(Gen_sum_co, abun2, by="Genus", all.x=T)

# subset for ten most common genera  
df_final <- df %>%  filter(Rank<11)

# calculate relative frequencies and remove unknown genus from data frame
rel_freq <- group_by(df_final, Treatment) %>%
  dplyr::mutate(., prop = n/sum(n))

# now calculate mean CAS production across treatments -> does not work!!
rel_freq_2 <- rel_freq  %>%
  group_by(Genus) %>%
  dplyr::mutate(m=mean(CAS_mean, na.rm=TRUE)) %>%
  ungroup()

levels(rel_freq_2$Genus)

# Sort factor levels by the factor level mean of another covariate
sortLvlsByVar.fnc <- function(oldFactor, sortingVariable, ascending = TRUE) {
  
  require("dplyr")
  require("magrittr")
  
  # Combine into data frame
  df <- data.frame(oldFactor, sortingVariable)
  
  # Compute average of sortingVariable and arrange (ascending)
  if (ascending == TRUE) {
    df_av <- df %>% group_by(oldFactor) %>% dplyr::summarise(meanSortingVariable = mean(sortingVariable)) %>% 
      arrange(meanSortingVariable)
  }
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = df_av$oldFactor)
  return(newFactor)
}

rel_freq_2$Genus <- sortLvlsByVar.fnc(rel_freq_2$Genus, rel_freq_2$m)

# colour palette
col1 <- "Gold"
col2 <- wes_palette(n=1, name="Moonrise3")
col3 <- wes_palette(n=1, name="FantasticFox1")
col4 <- wes_palette(n=5, name="FantasticFox1")

col_A <- colorRampPalette(c(col1, col2))(6)
col_B <- colorRampPalette(c(col3, col4[5]))(4)
col_C <- c(col_A, col_B)

# select data needed for plotting
df_plot <-  dplyr::select(rel_freq_2, Genus,Treatment, m, prop)  %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("PF","B"), labels=c("Control", "Copper")))

# actual plot & arringing fill according to levels of Genus (sorted by increasing mean siderophore production)
Fig_3D <- ggplot(arrange(df_plot, m), aes(x=Treatment, y=prop, fill=Genus))
Fig_3D <- Fig_3D + geom_bar(stat="identity", colour="black", size=0.4)
Fig_3D <- Fig_3D + scale_fill_manual(values=col_C, name="Genera\\n")
Fig_3D <- Fig_3D + xlab("\\nSelection regime") 
Fig_3D <- Fig_3D + ylab("Relative abundance\\n")
Fig_3D <- Fig_3D + scale_y_continuous(expand = c(0,0.02))
Fig_3D <- Fig_3D + theme_bw()

Fig_3D <- Fig_3D + theme(axis.line = element_line(colour = "black"), panel.grid.minor =  element_blank(), panel.grid.major = element_blank(), panel.border=element_rect(colour="black", size=0.95))

Fig_3D <- Fig_3D + theme(legend.title = element_text(size=12, face="bold"))
Fig_3D <- Fig_3D + theme(legend.text = element_text(size=12, face="italic"))
Fig_3D <- Fig_3D + theme(legend.key.size = unit(0.7, "cm"))
Fig_3D <- Fig_3D + theme (legend.position = "right")

Fig_3D <- Fig_3D + theme(axis.title.x = element_text(size=16))
Fig_3D <- Fig_3D + theme(axis.title.y = element_text(size=16))
Fig_3D <- Fig_3D + theme(axis.text.y = element_text(size=14))
Fig_3D <- Fig_3D + theme(axis.text.x = element_text(size=14))

Fig_3D <- Fig_3D + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_3D 

# combine two altered plots -> 14 inch by 10 inch

Fig.3 <- ggdraw(xlim=c(-0.2,3.3),ylim=c(-0.6,4.2)) +
  draw_plot(Fig_3A, 0, 1.9, 1.4, 2) +
  draw_plot(Fig_3B, 1.63, 1.9, 1.4, 2) +
  draw_plot(Fig_3C, 0.02, -0.405, 1.38, 2) +
  draw_plot(Fig_3D, 1.576, -0.405, 1.42, 2) +
  draw_plot_label(c("A", "B", "C", "D"), c(0.2, 1.8, 0.2, 1.8), c(4.1, 4.1,1.8, 1.8), size = 24)

Fig.3


# -------------------------- # 
# Code accompanying Figure 4 #
# -------------------------- # 

# ----- #
# STATS #
# ----- #

SID_PF_social <- SID_T6 %>% filter (Sample=="Wt") %>%
  dplyr::mutate(., Cu = factor(Cu, levels = c("1","0")))

model_lmer_PF_social <- lmer(CAS_norm ~ Cu * Culture + (1|Rep), data=SID_PF_social)
print(summary(model_lmer_PF_social))

model_lmer_PF_social_b <- lmer(CAS_norm ~ Cu + Culture + (1|Rep), data=SID_PF_social)
print(summary(model_lmer_PF_social_b))

# test two-way interaction
anova(model_lmer_PF_social, model_lmer_PF_social_b)

# test main effect of community
model_lmer_PF_social_c <- lmer(CAS_norm ~ Cu + (1|Rep), data=SID_PF_social)
print(summary(model_lmer_PF_social_c))

anova(model_lmer_PF_social_b, model_lmer_PF_social_c)

model_lmer_PF_social_d <- lmer(CAS_norm ~ 1 + (1|Rep), data=SID_PF_social)
print(summary(model_lmer_PF_social_d))

anova(model_lmer_PF_social_c, model_lmer_PF_social_d)

# package lsmeans
ls <-lsmeans(model_lmer_PF_social, ~Cu)

contrast(ls, alpha=0.05, method="pairwise", adjust='bonferroni')

# ------- #
# PANEL A #
# ------- #

# tidy up data set
SID_co_PF <- sid_df %>% filter(Sample=="Wt" & Culture=="Co")  %>%
  dplyr::mutate(., Cu = factor(Cu, levels = c("1","0"))) %>%
  dplyr::mutate(., Time = factor(Time, levels = c("0","6"))) %>%
  dplyr::mutate(., new = factor(new, levels = c("I", "0","1"))) %>%
  arrange(., Time) %>%
  dplyr::mutate(., treatment = as.factor(c(rep("I",12), rep("con",6), rep("cu",6)))) %>%
  dplyr::mutate(., treatment = factor(treatment, levels = c("I", "con", "cu"), labels = c("ancestral", "control", "copper")))

# these are the data for mono-cultures needed for Fig 4B
SID_mo_PF <- sid_df %>% filter(Sample=="Wt" & Culture=="Mono")  %>%
  dplyr::mutate(., Cu = factor(Cu, levels = c("1","0"))) %>%
  dplyr::mutate(., Time = factor(Time, levels = c("0","6"))) %>%
  dplyr::mutate(., new = factor(new, levels = c("I", "0","1"))) %>%
  arrange(., Time) %>%
  dplyr::mutate(., treatment = as.factor(c(rep("I",12), rep("con",6), rep("cu",6)))) %>%
  dplyr::mutate(., treatment = factor(treatment, levels = c("I", "con", "cu"), labels = c("ancestral", "control", "copper")))

# set up colour scheme
control <- wes_palette(n=5, name="Darjeeling1")
control_col <- control[5]

copper <- wes_palette(n=1, name="Darjeeling1")
copper_col <- copper[1]

# plot raw data using different colours 
dodge = position_dodge(width=0.8)

Fig_4A <- ggplot(data=SID_co_PF, aes(x=treatment, y=mean_sids))
Fig_4A <- Fig_4A + geom_boxplot(aes(fill=treatment), outlier.size = 0) + scale_fill_manual(values=c("black", control_col, copper_col))

dat_co <- ggplot_build(Fig_4A)$data[[1]]

xmin_co = dat_co$xmin[1]
xmax_co = dat_co$xmax[1]
middle_co = dat_co$middle[1]

Fig_4A <- Fig_4A + geom_segment(aes(x=xmin_co, y=middle_co, xend=xmax_co, yend=middle_co), col="white", size=0.6)

Fig_4A <- Fig_4A + coord_cartesian(ylim=c(0.6, 0.8)) + scale_y_continuous(breaks=seq(0.6, 0.8, 0.05))

Fig_4A <- Fig_4A + geom_jitter(width=0.15, size=3, shape=21,  stroke=1, col="black", fill="white") 

Fig_4A <- Fig_4A + theme_bw()

Fig_4A <- Fig_4A + xlab("\\nSelection regime") 

Fig_4A <- Fig_4A + ylab("Siderophore production\\n")

Fig_4A <- Fig_4A + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_4A <- Fig_4A + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=18)
) 

Fig_4A <- Fig_4A + theme(axis.title.x = element_text(size=16))
Fig_4A <- Fig_4A + theme(axis.text.x = element_text(size=14))
Fig_4A <- Fig_4A + theme(axis.text.y = element_text(size=14))
Fig_4A <- Fig_4A + theme(legend.position="none")
Fig_4A <- Fig_4A + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_4A <- Fig_4A + theme(axis.title.y = element_text(size=16))

# add text to graph
# create summarizing data frame
stat_data <- data.frame(text = c("","a","b"),
                        y = dat_co$ymin-0.025,
                        x = dat_co$xid, stringsAsFactors = FALSE)

Fig_4A <- Fig_4A + geom_text(data = stat_data, aes(label = text, x = x, y = y), colour = "black", size=6, fontface = "bold")

Fig_4A


# ------- #
# PANEL B #
# ------- #

dodge = position_dodge(width=0.8)

Fig_4B <- ggplot(data=SID_mo_PF, aes(x=treatment, y=mean_sids))
Fig_4B  <- Fig_4B  + geom_boxplot(aes(fill=treatment), outlier.size = 0) + scale_fill_manual(values=c("black", control_col, copper_col))

dat_mono <- ggplot_build(Fig_4B )$data[[1]]

xmin=dat_mono$xmin[1]
xmax=dat_mono$xmax[1]
middle=dat_mono$middle[1]

Fig_4B  <- Fig_4B + geom_segment(aes(x=xmin, y=middle, xend=xmax, yend=middle), col="white", size=0.6)

Fig_4B  <- Fig_4B + coord_cartesian(ylim=c(0.6, 0.8)) + scale_y_continuous(breaks=seq(0.6, 0.8, 0.05))

Fig_4B  <- Fig_4B + geom_jitter(width=0.15, size=3, shape=21,  stroke=1, col="black", fill="white") 

Fig_4B  <- Fig_4B + theme_bw()

Fig_4B  <- Fig_4B + xlab("\\nSelection regime") 

Fig_4B  <- Fig_4B + ylab("Siderophore production\\n")

Fig_4B  <- Fig_4B + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_4B  <- Fig_4B + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=18)
) 

Fig_4B  <- Fig_4B + theme(axis.title.x = element_text(size=16))
Fig_4B  <- Fig_4B + theme(axis.text.x = element_text(size=14))
Fig_4B  <- Fig_4B + theme(axis.text.y = element_text(size=14))
Fig_4B  <- Fig_4B + theme(legend.position="none")
Fig_4B  <- Fig_4B + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_4B  <- Fig_4B + theme(axis.title.y = element_text(size=16))

# add text to graph
# create summarizing data frame
stat_data_mono <- data.frame(text = c("","a","b"),
                        y = dat_mono$ymin-0.025,
                        x = dat_mono$xid, stringsAsFactors = FALSE)

Fig_4B  <- Fig_4B + geom_text(data = stat_data_mono, aes(label = text, x = x, y = y), colour = "black", size=6, fontface = "bold")

Fig_4B


# combine graphs and save pdf <- 12 x 6.85 inches
Fig.4 <- ggdraw(xlim=c(-0.1,3.2),ylim=c(-0.1,1.62)) +
  draw_plot(Fig_4A, 0, 0, 1.45, 1.5) +
  draw_plot(Fig_4B, 1.6, -0.006, 1.5, 1.508) +
  draw_plot_label(c("A", "B"), c(0.25, 1.85), c(1.6, 1.6), size = 24)

Fig.4

# -------------------------- #
# Code accompanying Figure 5 #
# -------------------------- #

growth_FDS <- read.csv("PF_m_FDS.csv", header=T)
attach(growth_FDS)
names(growth_FDS)

growth_FDS <- growth_FDS %>% dplyr::mutate(frequency = factor(frequency, levels = c("mono", "common", "rare"), labels=c("mono", "common","rare")), focal = factor(focal, levels=c("wt", "ch"), labels = c("producer", "cheat")))

# ---- #
# PLOT #
# ---- #

dodge = position_dodge(width=0.8)

Fig_5A <- ggplot(growth_FDS, aes(x=frequency, y=m,  colour=frequency, frquency))

Fig_5A <- Fig_5A + geom_hline(yintercept=0, col="grey")


Fig_5A <- Fig_5A + geom_boxplot(outlier.alpha = 0, fill = c("black", "gold", "skyblue","black", "gold", "skyblue"), colour=c("black", "gold", "skyblue","black", "gold", "skyblue")) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, colour="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

Fig_5A <- Fig_5A + geom_jitter(aes(x=frequency, y=m, fill = frequency, colour=frequency), width=0.15, size=3, shape=21,  stroke=0.5, col="black", fill="white") 

Fig_5A <- Fig_5A + facet_grid(~focal)

Fig_5A <- Fig_5A + coord_cartesian(ylim=c(-1, 2)) + scale_y_continuous(breaks=seq(-1, 2, 0.5))

Fig_5A <- Fig_5A + theme_bw()

Fig_5A <- Fig_5A + xlab("\\nSocial background") 

Fig_5A <- Fig_5A  + ylab("Malthusian growth rate (m)\\n")

Fig_5A <- Fig_5A + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_5A <- Fig_5A + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_5A <- Fig_5A+ theme(axis.title.x = element_text(size=16))
Fig_5A <- Fig_5A + theme(axis.text.x = element_text(size=14))
Fig_5A <- Fig_5A  + theme(axis.text.y = element_text(size=14))
Fig_5A <- Fig_5A + theme(legend.position="none")
Fig_5A <- Fig_5A + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_5A <- Fig_5A + theme(axis.title.y = element_text(size=16))

Fig_5A

# ----- #
# STATS #
# ----- #

model_FDS <- glm(m~frequency*focal, data=growth_FDS)
anova(model_FDS, test="F")

model_FDS_b <- glm(m~frequency+focal, data=growth_FDS)
anova(model_FDS, model_FDS_b, test="F")

ls_FDS<-lsmeans(model_FDS, ~focal*frequency)

cld(ls_FDS,Letters=letters,adjust='bonferroni')
contrast(ls_FDS, alpha=0.05, method="pairwise", adjust='bonferroni')

# 
dat_5A <- ggplot_build(Fig_5A)$data[[2]]

# add text to graph
stat_data_5A <- data.frame(text = c("ab","c","d","ab", "a", "b"),
                           y = dat_5A$ymax_final + 0.2,
                           focal = factor(c("producer", "producer", "producer", "cheat", "cheat", "cheat"), levels = c("producer", "cheat")),
                           x = dat_5A$xid, stringsAsFactors = FALSE)

Fig_5A <- Fig_5A + geom_text(data = stat_data_5A, aes(label = text, x = x, y = y, fill="white"), colour=c("black","black", "black", "black", "black", "black"), size=5, fontface="bold")

Fig_5A

# ------- #
# PANEL C #
# ------- #

growth_tr_mean_FDS <- growth_FDS %>% group_by(treatment, frequency, focal) %>% dplyr::summarise(mean_tr = mean(m, na.rm=T)) %>% ungroup() %>% filter(focal!="producer")

selection_wt_FDS <-  merge(growth_FDS, growth_tr_mean_FDS, by='treatment', all=TRUE) %>% filter(focal.x!="cheat") %>% rowwise() %>% mutate(sel_wt = m-mean_tr) %>% dplyr::select(-frequency.y, -focal.y, -mean_tr)

# ----- #
# STATS #
# ----- #

model_sel_FDS <-aov(sel_wt ~ frequency.x, data = selection_wt_FDS)
summary(model_sel_FDS)

model_sel_FDS_2 <-aov(sel_wt ~ 1, data = selection_wt_FDS)

anova(model_sel_FDS, model_sel_FDS_2)

ls_FDS_sel<-lsmeans(model_sel_FDS, ~frequency.x)

cld(ls_FDS_sel,Letters=letters,adjust='bonferroni')
contrast(ls_FDS_sel, alpha=0.05, method="pairwise", adjust='bonferroni')

# ---- #
# PLOT #
# ---- #
Fig_5B <- ggplot(selection_wt_FDS, aes(x=frequency.x, y = sel_wt, fill = frequency.x))

Fig_5B <- Fig_5B + geom_hline(yintercept=0, col="grey")

Fig_5B  <- Fig_5B + geom_boxplot(outlier.alpha = 0, fill=c("black", "gold", "skyblue"), colour=c("black", "gold", "skyblue")) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, colour="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

Fig_5B  <- Fig_5B + geom_jitter(width=0.15, size=3, stroke = 0.5, shape=21,  stroke=1, col="black", fill="white") 

Fig_5B  <- Fig_5B + coord_cartesian(ylim=c(-2, 0.75)) + scale_y_continuous(breaks=seq(-2, 0.75, 0.5))

Fig_5B  <- Fig_5B + theme_bw()

Fig_5B  <- Fig_5B + xlab("\\nSocial background") 

Fig_5B  <- Fig_5B  + ylab("Selection coefficient (s)\\n")

Fig_5B  <- Fig_5B + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_5B  <- Fig_5B + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_5B  <- Fig_5B + theme(axis.title.x = element_text(size=16))
Fig_5B  <- Fig_5B + theme(axis.text.x = element_text(size=14))
Fig_5B  <- Fig_5B  + theme(axis.text.y = element_text(size=14))
Fig_5B  <- Fig_5B + theme(legend.position="none")
Fig_5B  <- Fig_5B + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_5B  <- Fig_5B  + theme(axis.title.y = element_text(size=16))

Fig_5B

dat_Fig_5B <- ggplot_build(Fig_5B)$data[[2]]

# add text to graph
stat_data_Fig_5B <- data.frame(text = c("a","b", "c"),
                               y = dat_Fig_5B$ymax_final + 0.2,
                               x = dat_Fig_5B$xid, stringsAsFactors = FALSE)

Fig_5B <- Fig_5B + geom_text(data = stat_data_Fig_5B, mapping = aes(label = text, x = x, y = y, fill="white"), colour = "black", size=5, fontface = "bold")

Fig_5B


##########################################################################################################
# growth assays SWB25 producer and cheat in presence and absence of community in copper-polluted compost #
##########################################################################################################

# ------------ #
# PANELS B & D #
# ------------ #

growth <- read.csv("growth_PF_copper.csv", header=T)
attach(growth)
names(growth)

growth <- growth %>% dplyr::mutate(frequency = factor(frequency, levels = c("mono", "poly"), labels=c("mono", "co")), focal = factor(focal, levels=c("wt", "ch"), labels = c("producer", "cheat")))

# ---- #
# PLOT #
# ---- #

dodge = position_dodge(width=0.8)

Fig_5C <- ggplot(growth, aes(x=frequency, y=m, fill = frequency, colour=frequency))

Fig_5C <- Fig_5C + geom_hline(yintercept=0, col="grey")

Fig_5C = Fig_5C + geom_boxplot(outlier.alpha = 0, colour=c("black", "tomato1", "black", "tomato1"), fill = c("black", "tomato1", "black", "tomato1")) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, colour="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

Fig_5C = Fig_5C + geom_jitter(aes(x=frequency, y=m, fill = frequency, colour=frequency), width=0.15, size=3, shape=21,  stroke=0.5, col="black", fill="white")

Fig_5C = Fig_5C + facet_grid(~focal)

Fig_5C = Fig_5C + coord_cartesian(ylim=c(-1, 1)) + scale_y_continuous(breaks=seq(-1, 1, 0.5))

Fig_5C = Fig_5C + theme_bw()

Fig_5C = Fig_5C + xlab("\\nSocial background") 

Fig_5C = Fig_5C  + ylab("Malthusian growth rate (m)\\n")

Fig_5C = Fig_5C + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_5C = Fig_5C + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_5C = Fig_5C + theme(axis.title.x = element_text(size=16))
Fig_5C = Fig_5C + theme(axis.text.x = element_text(size=14))
Fig_5C = Fig_5C + theme(axis.text.y = element_text(size=14))
Fig_5C = Fig_5C + theme(legend.position="none")
Fig_5C = Fig_5C + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_5C = Fig_5C + theme(axis.title.y = element_text(size=16))

Fig_5C

# ----- #
# STATS #
# ----- #

model <- glm(m~frequency*focal, data=growth)
anova(model, test="F")

model_b <- glm(m~frequency+focal, data=growth)
anova(model, model_b, test="F")

ls<-lsmeans(model, ~focal*frequency)

print(ls)
cld(ls,Letters=letters,adjust='bonferroni')
contrast(ls, alpha=0.05, method="pairwise", adjust='bonferroni')

# 
dat_Fig_5C <- ggplot_build(Fig_5C)$data[[2]]

# add text to graph
stat_data_Fig_5C <- data.frame(text = c("a","b","a","a"),
                               y = dat_Fig_5C$ymax_final + 0.2,
                               focal = factor(c("producer", "producer", "cheat", "cheat"), levels = c("producer", "cheat")),
                               x = dat_Fig_5C$xid, stringsAsFactors = FALSE)

Fig_5C <- Fig_5C + geom_text(data = stat_data_Fig_5C, aes(label = text, x = x, y = y, fill="white"), colour=c("black","black", "black", "black"), size=5, fontface="bold")

Fig_5C

# ------- #
# PANEL D #
# ------- #

growth_tr_mean <- growth %>% group_by(treatment, frequency, focal) %>% dplyr::summarise(mean_tr = mean(m, na.rm=T)) %>% ungroup() %>% filter(focal!="producer")

selection_wt <-  merge(growth, growth_tr_mean, by='treatment', all=TRUE) %>% filter(focal.x!="cheat") %>% rowwise() %>% mutate(sel_wt = m-mean_tr) %>% dplyr::select(-frequency.y, -focal.y, -mean_tr)

# first model this
model <-aov(sel_wt ~ frequency.x, data = selection_wt)
summary(model)

model_2 <-aov(sel_wt ~ 1, data = selection_wt)
summary(model_2)

anova(model, model_2)

ls_model <- ls(model, ~frequency.x)
cld(ls_model,Letters=letters,adjust='bonferroni')
contrast(ls, alpha=0.05, method="pairwise", adjust='bonferroni')

# ---- #
# PLOT #
# ---- #

Fig_5D <- ggplot(selection_wt, aes(x=frequency.x, y = sel_wt, fill = frequency.x))

Fig_5D <- Fig_5D + geom_hline(yintercept=0, col="grey")

Fig_5D <- Fig_5D + geom_boxplot(outlier.alpha = 0, fill= c("black", "tomato1"), colour=c("black", "tomato1")) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, colour="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

Fig_5D <- Fig_5D + geom_jitter(width=0.15, size=3, shape=21,  stroke=0.5, col="black", fill="white") 

Fig_5D <- Fig_5D + coord_cartesian(ylim=c(-1, 0.75)) + scale_y_continuous(breaks=seq(-1, 0.75, 0.5))

Fig_5D <- Fig_5D + theme_bw()

Fig_5D <- Fig_5D + xlab("\\nSocial background") 

Fig_5D <- Fig_5D  + ylab("Selection coefficient (s)\\n")

Fig_5D <- Fig_5D + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_5D <- Fig_5D + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_5D <- Fig_5D + theme(axis.title.x = element_text(size=16))
Fig_5D <- Fig_5D + theme(axis.text.x = element_text(size=14))
Fig_5D <- Fig_5D  + theme(axis.text.y = element_text(size=14))
Fig_5D <- Fig_5D + theme(legend.position="none")
Fig_5D <- Fig_5D + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_5D <- Fig_5D  + theme(axis.title.y = element_text(size=16))

Fig_5D 

dat_Fig_5D  <- ggplot_build(Fig_5D )$data[[2]]

# add text to graph
stat_data_Fig_5D  <- data.frame(text = c("a","b"),
                                y = dat_Fig_5D $ymax_final + 0.2,
                                x = dat_Fig_5D $xid, stringsAsFactors = FALSE)

Fig_5D  <- Fig_5D  + geom_text(data = stat_data_Fig_5D , mapping = aes(label = text, x = x, y = y, fill="white"), colour = "black", size=5, fontface = "bold")

Fig_5D 

# combine these into one plot with four panels <- Fig. 5 -> 14 inch by 10 inch

Fig.5 <- ggdraw(xlim=c(-0.2,3.1),ylim=c(-0.6,4.3)) +
  draw_plot(Fig_5A, 0, 1.9, 1.6, 2.132) +
  draw_plot(Fig_5B, 1.83, 1.9, 0.8, 2) +
  draw_plot(Fig_5C, 0, -0.5, 1.6, 2.132) +
  draw_plot(Fig_5D, 1.83, -0.5, 0.8, 2) +
  draw_plot_label(c("A", "C", "B", "D"), c(0.25, 2.07, 0.25, 2.07), c(4.2, 4.2, 1.8, 1.8), size = 24)

Fig.5

##########
# Fig_S1 #
##########

den<-read.csv("Densities.csv", header=TRUE)
attach(den)
names(den)

den_PF <- den %>% filter (Strain =="Wt", Time=="6") %>% dplyr::mutate(Cu=factor(Cu, levels=c("0","1"), labels=c("control", "copper"))) %>% dplyr::mutate(Culture = factor(Culture, levels = c("Mono", "Co"), labels=c("monoculture", "co-culture")))

# set up colour scheme
control <- wes_palette(n=5, name="Darjeeling1")
control_col <- control[5]

copper <- wes_palette(n=1, name="Darjeeling1")
copper_col <- copper[1]

# ------- # 
# PLOT S1 #
# ------- # 

Fig_S1 <- ggplot(den_PF, aes(x=Cu, y=Density, colour = Cu))

Fig_S1  = Fig_S1 + facet_grid(~ Culture)

Fig_S1  = Fig_S1  + geom_boxplot(outlier.alpha = 0, colour=c(control_col, copper_col,control_col, copper_col), fill = c(control_col, copper_col,control_col, copper_col)) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, colour="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

Fig_S1  = Fig_S1 + geom_jitter(aes(x=Cu, y=Density, colour=Cu, fill = Cu), width=0.15, size=3, shape=21,  stroke=0.5, col="black", fill="white")

# add text then change to log10 scale
dat_Fig_S1 <- ggplot_build(Fig_S1)$data[[1]]

# add text to graph
S1 <- data.frame(text = c("a","b","c","d"),
                 y = c(10^8.2, 10^7.6, 10^5.4, 10^4.4),
                 Culture = factor(c("monoculture", "monoculture", "co-culture", "co-culture")),
                 Cu = factor(c("control", "copper", "control", "copper")),
                 x = dat_Fig_S1$xid, stringsAsFactors = FALSE)

Fig_S1 <- Fig_S1 + geom_text(data = S1, aes(label = text, x = x, y = y, colour=Cu), colour=c("black","black", "black", "black"), size=5, fontface="bold")

Fig_S1  = Fig_S1 + scale_y_log10(limits = c(1,1e9), expand = c(0, 0),
                                 breaks = scales::trans_breaks("log10", function(x) 10^x),
                                 labels = scales::trans_format("log10", scales::math_format(10^.x)))

Fig_S1  = Fig_S1 + annotation_logticks(base = 10, sides = "l", short = unit(0.1, "cm"))

Fig_S1  = Fig_S1 + theme_bw()

Fig_S1  = Fig_S1 + xlab("\\nTreatment") 

Fig_S1  = Fig_S1  + ylab("Density (log scale)\\n")

Fig_S1  = Fig_S1 + theme(panel.border = element_rect(color = "black", fill = NA, size = 1), panel.grid.minor=element_blank(), panel.grid.major = element_blank(),  panel.background = element_rect(fill = 'white'))

Fig_S1  = Fig_S1 + theme(
  strip.background = element_blank(),
  strip.text.x = element_text(size=16)
) 

Fig_S1  = Fig_S1 + theme(axis.title.x = element_text(size=16))
Fig_S1  = Fig_S1 + theme(axis.text.x = element_text(size=14))
Fig_S1  = Fig_S1 + theme(axis.text.y = element_text(size=14))
Fig_S1  = Fig_S1 + theme(legend.position="none")
Fig_S1  = Fig_S1 + theme(axis.title.y=element_text(margin=margin(0,2,0,0)))

Fig_S1  = Fig_S1 + theme(axis.title.y = element_text(size=16))

Fig_S1

# ----- #
# STATS #
# ----- #

model <- glm(log10(Density)~Cu*Culture, data=den_PF)
anova(model, test="F")

model_b <- glm(log10(Density)~Cu+Culture, data=den_PF)
summary(model_b)

anova(model, model_b, test="F")

model_c <- glm(log10(Density)~Cu, data=den_PF)
anova(model_b, model_c, test="F")

model_d <- glm(log10(Density)~Culture, data=den_PF)
anova(model_b, model_d, test="F")

ls<-lsmeans(model_b, ~Cu+Culture)

cld(ls,Letters=letters,adjust='bonferroni')
contrast(ls, alpha=0.05, method="pairwise", adjust='bonferroni')