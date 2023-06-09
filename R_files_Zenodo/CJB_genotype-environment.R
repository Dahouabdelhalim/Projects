#R script used in Bracic, Bohn et al (2022)
#author: Marko Bracic

###################################################################################################
# Model for analyzing influence of genotype and environment on CJB
###################################################################################################
# This code generates values from the statistical analysis of cognitive judgment bias test used in:
  # In-text values in the manuscript Results section
  # Figure 3: Cognitive judgment bias
  # Supplementary Figure 2: Choice score for each cue 
  # Supplementary Table 3: Summary statistic and pairwise comparison of choice scores 
  # Supplementary Table 4: Statistical analysis of CJB test and behavioral test battery


### Load packages ###
library(lme4)#mixed models
library(lmerTest)#gives p-values from lme4 models
library(tibble) #data wrangling 
library(dplyr) #data wrangling 
library(ggplot2)#plotting
library(ggpubr)#plotting
library(ggbeeswarm)#plots that spread outs dots instead of jittering it
library(RVAideMemoire)# for Wilcox paired test


### Loading the data ###

input <- read.table("CJB_genotype-environment.txt", header = TRUE, stringsAsFactors = FALSE)


### Preparing the data ###

# Using just ambiguous cues, set them as factors and preparing coding types
md <- input
md = subset(md, md$cue== "NP" | md$cue== "M" | md$cue== "NN")
md = droplevels (md)
md$cue = factor(md$cue, c("NP", "M", "NN"))
md$cueC = relevel(md$cue, ref="M")
  # M cue set as reference cue -> not relevant for statistical model

# Sum-contrast coding of cue three-level factor
md$cue_CS = md$cue
contrasts(md$cue_CS) <- contr.sum(3) #three-level factor

# Genotype as factor with 2 levels: C57BL/6 (C57) and B6D2F1 (F1)
md$genotype = as.factor(md$genotype)
md$genotype = relevel(md$genotype, ref="C57BL/6")
  # C57BL/6 strain as reference -> not relevant for statistical model

# Environment as factor with 2 levels: scarce and complex environment
md$environment = as.factor(md$environment)
md$environment = relevel(md$environment, ref="scarce")
  # scarce environment as reference -> not relevant for statistical model

# Centering factors for extracting model estimates  -> increases model interpretation 
md$genotypeC = as.numeric(md$genotype=="B6D2F1") - 0.5
# == B6D2F1" gives B6D2F1=1 and C57BL/6=0
  # this makes C57BL/6 reference level and - 0.5 centers it (B6D2F1=0.5, C57BL/6=-0.5)
md$environmentC = as.numeric(md$environment=="complex") - 0.5
#centering and making st reference level, see logic above


### Linear mixed-effect model  ###

#Model with sum-contrast coded cues and centered G and E
modC = lmer(score ~ cue_CS * genotypeC * environmentC + (1|ind) + (1|cage), data=md)

# ANOVA type III tables: used to report model output
anova(modC) #anova function from lmerTest: sum-contrasts codes all the factors 

# ANOVA table for in-text stats and Supplementary Table 4 "Statistical analysis"
options(pillar.sigfig = 5) #regulates number of shown digits in tibble table
ANOVA_table= modC %>% 
  anova() %>% 
  rownames_to_column("Factor") %>% #from tibble::
  as_tibble() %>% #
  mutate(across(where(is.numeric), round,3)) %>% 
  mutate(`Pr(>F)`= format.pval(.$`Pr(>F)`, 2, 0.001,0)) #rounding bigger p-values to 2 digits
  #provides all values except model estimates (+SE)


#Extracting model estimates and SE (reported in-text)
summary(modC)

model_estimates <- summary(modC)$coefficients
model_estimates <- model_estimates %>% 
  subset( select = c("Estimate", "Std. Error")) %>% 
  round(3) 
#table with model estimates and SE values reported in-text


##########################################################################################
#                           Figure Cognitive judgment bias             
##########################################################################################
# This code generates figure titled Cognitive judgment bias

#extracting p-values from ANOVA table for ploting
stats = paste("Genotype: p = ", 
      ANOVA_table$`Pr(>F)`[ANOVA_table$Factor == "genotypeC"], 
      "\\nEnvironment: p = ", 
      ANOVA_table$`Pr(>F)`[ANOVA_table$Factor == "environmentC"], 
      "\\nGxE: p = ", 
      ANOVA_table$`Pr(>F)`[ANOVA_table$Factor == "genotypeC:environmentC"],
      "\\nCue: p < 0.001",#manual entry!
      sep="")

#preparing data frame for plot annotation (p values used from The ANOVA_table)
plot_text1 <- data.frame(
  cue = as.factor(c("NP", "M", "NN")),
  xstats = c(1, 1.5, 1),
  ystats = c(1, 0.9, 1),
  lab = c("", stats, "")
)


#beeswarm plot 
F3=ggplot(data= md, 
       aes(x=genotype, y=score)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,
                size=2,
                alpha=0.5,
                cex = 2.2) +
  stat_summary(aes(colour = environment),
               fun=median,geom="crossbar", 
               position=position_dodge(width=0.5), 
               width=0.3, 
               size=0.3,
               fatten=3)+
  stat_summary(aes(colour = environment),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  geom_text(data=plot_text1,  
            aes(x = xstats,  y = ystats, label = lab), 
            size = 2.3,
            hjust = 0.5) +
  ylim(-1.1,1.1) +
  ylab("Choice score") +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",
                      breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),
                      values = c("#fc7e4c", "#799df2")) + 
  scale_shape_manual(name="Environment",
                     breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),
                     values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6" = "C57", "B6D2F1" = "F1"),
                   expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(
          hjust = .5,
          vjust = .0,
          face = "plain",
          size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle ("Ambiguous cues")  +
  facet_wrap(vars(cue), 
             ncol=3, 
             strip.position = "top",
             labeller = labeller(cue = c("NP" = "Near Positive",
                                         "M" = "Middle", 
                                         "NN" = "Near Negative")))  

#tiff("Figure 3.tiff", res=350,width = 17.7, height = 12.8, units = "cm")
#F3
#dev.off()
#run the line above to save the figure


##########################################################################################
#       Exploring potential influence of design effects (reported in Material and Methods)          
##########################################################################################

#Influence of trial_type, batch and training duration
anova(lmer(score ~ trial_type + batch + training_sessions + 
             cue + genotype * environment + (1|ind) + (1|cage), data=md))



# Checking a model with all cues (not just ambiguous ones) 

# Preparing data set with all cues
md_all <- input
md_all$cue = factor(md_all$cue, c("P", "NP", "M", "NN", "N"))
#cue, genotype and environment can again be coded as factors/centered /contrast coded
  #but this will not influence the anova result 

# Model with all cues
anova(lmer(score ~ cue * genotype * environment + (1|ind) + (1|cage), 
                 data=md_all))
#similar result as the model with just ambiguous cues


##########################################################################################
#             Additional analysis reported in Supplementary Data          
##########################################################################################

# Preparing data set with all cues (again)
md_all <- input
md_all$cue = factor(md_all$cue, c("P", "NP", "M", "NN", "N"))


# Supp. Figure 2. Response curve 
ggplot(md_all, aes(x=cue, y=score)) +
  stat_summary(aes(y = score,group=1), fun=median,geom="line",group=1) +
  stat_summary(aes(y = score,group=1), fun=median,geom="point",group=1, size=2) +
  stat_summary(aes(y = score),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               geom = "errorbar", width=0.1) +
  ylab("Choice score") +
  xlab("") +
  scale_x_discrete(labels=c("P" = "Positive", "NP" = "Near Positive",
                            "M" = "Middle", "NN" = "Near Negative", "N" = "Negative")) +
  theme_classic() +
  theme(axis.title=element_text(size=20), axis.text  = element_text(size=16))
#ggsave("Supplementary Figure 2.tiff", width = 24.1, height = 15, units = "cm", dpi=600)



### Summary statistic and pairwise comparison of cues (reported in Sup.Table 3 )###

# Summary statistics
md_all$Cue = factor(md_all$cue, 
                    levels = c("P", "NP", "M", "NN", "N"),
                    labels = c("Positive", "Near Positive", "Middle","Near Negative", "Negative"))

# Preparing a table
raw_table = md_all %>% group_by(Cue) %>% 
  summarise(Mean = mean(score),
            SD = sd(score),
            Median=median(score),
            "Range (min, max)" = paste (round(min(score),2), 
                                        round(max(score), 2), sep=", ")) %>% 
  mutate(across(where(is.numeric), round, digits=2))

raw_table
#final summary statistic table 


# Pairwise cue comparison 

# Wilcoxon signed-rank test
cue_diff_wilcox = wilcox.paired.multcomp( score ~ cue | ind, 
                                          md_all,
                                          p.method = "holm")

cue_diff_wilcox
#final cue comparison



### Full model test (reported in Supplementary Table 4 "Statistical analysis" )###

mod_final = lmer(score ~ cue_CS * genotypeC * environmentC + (1|ind) + (1|cage), 
                 data=md, REML = FALSE)
mod_null = lmer(score ~  1 + (1|ind) + (1|cage), data=md, REML = FALSE)
anova(mod_final, mod_null)