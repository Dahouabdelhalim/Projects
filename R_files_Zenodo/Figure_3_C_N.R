#Ryan Stephens; Finalized November 9, 2020
rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/TDF_Lit_Review")

library(dplyr)
################################################################################################
#Get data ready
################################################################################################
C_N<- read.csv("MetaAnalysis_TDF_data.csv",header=T)#data
head(C_N)

library(forcats)#use forcats to rename levels
C_N<-mutate(C_N, Diet_Source = fct_recode(Diet_Source, "C[3]" = "C3", "C[4]" = "C4"))#rename diet source for graphing

#remove outliers which strongly influence results
library(ggplot2)
ggplot(C_N, aes(x = Diet_Source, y = Diet_C.N, fill = Diet_Source))+
  geom_boxplot()

C_N<-filter(C_N, Diet_C.N < 40)#remove two outliers (one for C3 and one for mixed diets)

#summary range of C/N for the different diet sources
C_N %>% 
  group_by(Diet_Source )%>%
  summarise(min= min(Diet_C.N), max=max(Diet_C.N))
################################################################################################




################################################################################################
#C/N diet vs TDF
################################################################################################


##############################
#Regressions for C/N vs TDF for C3 and mixed diets
##############################
C_N_C3_Mixed<-filter(C_N, Diet_Source == "C[3]" | Diet_Source == "Mixed" )#filter C3 and mixed diets

#Use broom to do regression for combinations of both isotopes and diets
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
regressions <- C_N_C3_Mixed %>%
  nest(data = -Diet_Source, -Isotope) %>% 
  mutate(
    fit = map(data, ~ lm(Diet_C.N ~ TDF, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

regressions %>% 
  unnest(tidied)

CN_TDF<-regressions%>%#regressions with p-values and r squared values
  unnest(glanced)

#R squared
CN_TDF$r.squared<- sprintf('"%.2f"',round(CN_TDF$r.squared, 2))#round to 2 sig digits and keep trailing zero
CN_TDF$r.squared_graph<- paste("italic(R)^2 ==~", CN_TDF$r.squared)#use for graphing r squared
#P-value
CN_TDF$p.value<- ifelse(CN_TDF$p.value < .001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(CN_TDF$p.value,digits=3)))#adjust P-value for graphing
CN_TDF$p.value_graph<- ifelse(CN_TDF$p.value=="\\"0.0001\\"",paste("italic(P)<",CN_TDF$p.value), paste("italic(P)==",CN_TDF$p.value))#P-value for graphing
head(CN_TDF)  
#Sample size
CN_TDF$n<- paste("italic(n)==",CN_TDF$nobs)#use n-observations for sample size
##############################


##############################
#graphing
library(lemon)
library(ggplot2)
library(ggpmisc)

label_parse <- function(breaks) {#parse text for legend with subscript
  parse(text = breaks)
}

#add letters for the figure
C_N_C3_Mixed<-mutate(C_N_C3_Mixed, Letters = ifelse(Isotope == "d13C", "(b)", "(d)"))

TDF_graph<-ggplot(C_N_C3_Mixed, aes(x = Diet_C.N, y = TDF, colour = Diet_Source, fill=Diet_Source))+
#Regression lines
  stat_smooth(aes(color = Diet_Source), method = "lm", alpha=.25)+ #regression lines for C3 and mixed diets
#point shape and colors
  geom_point(aes(shape=Consumer_type,colour = Diet_Source, fill=Diet_Source),size=1.7, alpha=.65)+#shape by consumer type and color by diet source
  scale_color_manual(values =c("chartreuse4","mediumorchid3"),label = label_parse)+#custom outline color colors
  scale_fill_manual(values =c("chartreuse4","mediumorchid3"))+#custom fill colors
  scale_shape_manual(values =c(24,21,22))+#shapes for diet class
#background and facets
  facet_rep_wrap(.~ interaction(Isotope), labeller = "label_parsed", scales='free_y')+#"lemon" package for facets (parsed for isotopes labels of facets)
  coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+#remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#removes gridlines, but keeps border
  theme(strip.text.x = element_blank())+#remove facet labels for isotopes
  theme(panel.grid.major.x = element_blank())+#remove facets
#labels for regression results, protein/energy rich, letters
  #n
  geom_text(data=CN_TDF,aes(x=c(21,.2,5,20),y= c(6,.6,7.6,7.6), hjust = 0, label=(n)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
              size=2.7,inherit.aes=F, color =c("chartreuse4","mediumorchid3","chartreuse4","mediumorchid3"), alpha = .65)+ #size and colors
  #p-value
  geom_text(data=CN_TDF,aes(x=c(21,.2,5,20),y= c(5.4,0,7,7), hjust = 0, label=(p.value_graph)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("chartreuse4","mediumorchid3","chartreuse4","mediumorchid3"), alpha = .65)+ #size and colors
  #r-squared
  geom_text(data=CN_TDF,aes(x=c(21,.2,5,20),y= c(4.8,-.6,6.4,6.4), hjust = 0, label=(r.squared_graph)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
              size=2.7,inherit.aes=F, color =c("chartreuse4","mediumorchid3","chartreuse4","mediumorchid3"), alpha = .65)+ #size and colors
  #Energy rich and protein rich
  annotate("text", label = "High-protein", x = 0, y = -2.5, size = 2.7, hjust = 0, colour = "black", alpha = .65)+#protein-rich
  annotate("text", label = "Low-protein", x = 24, y = -2.5, size = 2.7, hjust = 0, colour = "black", alpha = .65)+#energy-rich
  annotate("segment", x = 12, xend = 23, y = -2.5, yend = -2.5, colour = "gray65", lwd=.6,  linejoin = "round", arrow=arrow(length=unit(2, "mm")))+
  #letters
  geom_text(aes(x = -Inf, y = Inf, vjust=1.2, hjust=-.2,label = Letters), color = "black", size=3.7)+#letters
#Legend
  guides(color=guide_legend(override.aes=list(fill=NA),#remove fill behind color legend for Diet Source
         title = "Diet"),#rename title
         fill=FALSE,#remove the fill legend
         shape = guide_legend(title = "Consumer"))+#rename title
  theme(legend.box.background = element_rect(colour = "gray90"),#add rectangle around legend
    legend.box = "horizontal",#put legends side by side
    legend.position = c(0.83, 0.21),#legend position
    legend.text.align = 0,#make text align on left for Diet source
    legend.background=element_blank(),#remove legend background
    legend.text=element_text(size=7),#text size
    legend.title = element_text(size=8),#legend title text size
    legend.margin = margin(.1,.1,.1,.1, unit="cm"),
    legend.key.size = unit(.15, "cm"))+#size of legend
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-3, 9),breaks=c(-2,0,2,4,6,8))+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 35),breaks=c(0, 10, 20, 30))+
  ylab("Trophic discrimination factor (TDF)")+#y axis title
  xlab("C:N (Diet)")+#x axis title
  theme(axis.title.x = element_text(colour="black", size=10),
        axis.text.x  = element_text(angle=0, hjust = .5,size=8,colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=10),
        axis.text.y  = element_text(angle=0, hjust = .5,size=8,colour="black"))

TDF_graph#print to check
##############################

#regression with Nyctalus_noctula removed (bat with low C/N but high d13C TDF of 5.9 removed)
C_new<-filter(C_N, Reference != "Roswag et al. (2015)" & Isotope == "d13C" & Diet_Source == "C[3]")
C_lm<-lm(Diet_C.N ~ TDF, data = C_new)
summary(C_lm)
################################################################################################




################################################################################################
#C/N diet vs stable isotope ratio of diet
################################################################################################
C_N<-mutate(C_N, Isotope = fct_recode(Isotope, "paste(delta^13 * \\"C (\\u2030)\\")" = "d13C",#rename for graphing
                                               "paste(delta^15 * \\"N (\\u2030)\\")" = "d15N"))

##############################
#Regressions for C/N vs isotopic ratio for C3 and mixed diets
##############################
C_N_C3_Mixed<-filter(C_N, Diet_Source == "C[3]" | Diet_Source == "Mixed" )#filter C3 and mixed diets

#Use broom to do regression for combinations of both isotopes and diets
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
regressions <- C_N_C3_Mixed %>%
  nest(data = -Diet_Source, -Isotope) %>% 
  mutate(
    fit = map(data, ~ lm(Diet_C.N ~ Diet_isotope, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

regressions %>% 
  unnest(tidied)

CN_Diet_isotope<-regressions%>% #regressions with p-values and r squared values
  unnest(glanced)

#R squared
CN_Diet_isotope$r.squared<- sprintf('"%.2f"',round(CN_Diet_isotope$r.squared, 2))#round to 2 sig digits and keep trailing zero
CN_Diet_isotope$r.squared_graph<- paste("italic(R)^2 ==~", CN_Diet_isotope$r.squared)#use for graphing r squared
#P-value
CN_Diet_isotope$p.value<- ifelse(CN_Diet_isotope$p.value < .001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(CN_Diet_isotope$p.value,digits=3)))#adjust P-value for figure
CN_Diet_isotope$p.value_graph<- ifelse(CN_Diet_isotope$p.value=="\\"0.0001\\"",paste("italic(P)<",CN_Diet_isotope$p.value), paste("italic(P)==",CN_Diet_isotope$p.value))#P-value for graphing
head(CN_Diet_isotope)  
#Sample size
CN_Diet_isotope$n<- paste("italic(n)==",CN_Diet_isotope$nobs)# use n-observations for sample size
##############################


##############################
#graphing
library(lemon)
library(ggplot2)
library(ggpmisc)

label_parse <- function(breaks) {#parse text for legend
  parse(text = breaks)
}

#add letters for the figure
C_N_C3_Mixed<-mutate(C_N_C3_Mixed, Letters = ifelse(Isotope == "paste(delta^13 * \\"C (\\u2030)\\")", "(a)", "(c)"))

Diet_isotope_graph<- 
  ggplot(C_N_C3_Mixed, aes(x = Diet_C.N, y = Diet_isotope, colour = Diet_Source, fill=Diet_Source))+
#Regression lines
  stat_smooth(aes(color = Diet_Source), method = "lm", alpha=.25)+#regression lines for C3 and mixed diets
#point shape and colors
  geom_point(aes(shape=Consumer_type,colour = Diet_Source, fill=Diet_Source),size=1.7, alpha=.65)+#shape by consumer type and color by diet source
  scale_color_manual(values =c("chartreuse4","mediumorchid3"),label = label_parse)+#custom outline color colors
  scale_fill_manual(values =c("chartreuse4","mediumorchid3"))+#custom fill colors
  scale_shape_manual(values =c(24,21,22))+#shapes for diet class
#background and facets
  facet_rep_wrap(.~ interaction(Isotope), labeller = "label_parsed", scales='free_y')+#"lemon" package for facets (parsed for isotopes labels of facets)
  coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+#remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#removes gridlines, but keeps border
  
  theme(strip.background = element_rect(color=NA))+
  theme(strip.text.x = element_text(size =10,margin = margin(.05,0,.05,0, "cm")))+#white text and smaller box height
#labels for regression resutls, letters
  geom_text(data=CN_Diet_isotope,aes(x=c(1,20,1,20),y= c(-29,-24,2,13), hjust = 0, label=(n)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("chartreuse4","mediumorchid3","chartreuse4","mediumorchid3"), alpha = .65)+#size and colors
  geom_text(data=CN_Diet_isotope,aes(x=c(1,20,1,20),y= c(-30,-25,1,12), hjust = 0, label=(p.value_graph)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("chartreuse4","mediumorchid3","chartreuse4","mediumorchid3"), alpha = .65)+#size and colors
  geom_text(data=CN_Diet_isotope,aes(x=c(1,20,1,20),y= c(-31,-26,0,11), hjust = 0, label=(r.squared_graph)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("chartreuse4","mediumorchid3","chartreuse4","mediumorchid3"), alpha = .65)+#size and colors
  #letters
  geom_text(aes(x = -Inf, y = Inf, vjust=1.9, hjust=-.2,label = Letters), color = "black", size=3.7)+#letters
#Legend
    theme(legend.position = "none")+
#Axes
  scale_x_continuous(expand = c(0, 0),limits = c(0, 35),breaks=c(0, 10, 20, 30))+
  ylab("Stable isotopic ratio (Diet)")+#y axis title
  theme(axis.title.x= element_blank(),
        axis.text.x  = element_blank())+
  theme(axis.title.y = element_text(colour="black", size=10),
        axis.text.y  = element_text(angle=0, hjust = .5,size=8,colour="black"))

Diet_isotope_graph#print to check
##############################
################################################################################################




################################################################################################
#Final Plot
library(cowplot)

ggdraw() +
  draw_plot(Diet_isotope_graph, x = -.015, y = .5, width = 1.01, height = .5)+
  draw_plot(TDF_graph, x = 0, y = 0, width = 1, height = .5) 

ggsave("Figure_3_CN_Mixed_C3.tiff",
       plot = last_plot(), width = 4.5, height = 5.5, units = "in",
       dpi = 600) 
################################################################################################



