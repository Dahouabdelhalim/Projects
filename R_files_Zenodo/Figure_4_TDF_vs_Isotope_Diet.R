#Ryan Stephens; finalized November 18, 2020
rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/TDF_Lit_Review")

TDF_C_N<- read.csv("MetaAnalysis_TDF_data.csv",header=T)#data
head(TDF_C_N)

################################################################################################
#Get data ready
################################################################################################
library(dplyr)
TDF_C_N$Consumer_type <- factor(TDF_C_N$Consumer_type, levels = c("Herbivore","Omnivore","Carnivore"))#change order of consumers

#add letters for the figure
TDF_C_N<-mutate(TDF_C_N, Letters = ifelse(Isotope == "d13C" & Diet_Source == "C3","(a)",
                                   ifelse(Isotope == "d13C" & Diet_Source == "C4","(b)",
                                   ifelse(Isotope == "d13C" & Diet_Source == "Marine","(c)",     
                                   ifelse(Isotope == "d13C" & Diet_Source == "Mixed","(d)",  
                                   ifelse(Isotope == "d15N" & Diet_Source == "C3","(e)",
                                   ifelse(Isotope == "d15N" & Diet_Source == "C4","(f)",
                                   ifelse(Isotope == "d15N" & Diet_Source == "Marine","(g)","(h)"))))))))

library(forcats)#use forcats to rename levels
TDF_C_N<-mutate(TDF_C_N, Diet_Source = fct_recode(Diet_Source, "C[3]" = "C3", "C[4]" = "C4"))#rename diet sources

TDF_C_N<-mutate(TDF_C_N, Diet_Source_label = ifelse(Isotope == "d15N",paste(Diet_Source), NA))#use this column for adding diet source labels

TDF_C_N<-mutate(TDF_C_N, Average_line = ifelse(Isotope == "d15N",3, 1))#lines of commonly used TDFs for d13C (1) and d15N (3)

TDF_C_N<-mutate(TDF_C_N, Isotope = fct_recode(Isotope, "paste(delta^13 * \\"C (\\u2030)\\")" = "d13C",#rename d13C and d15N for graphing
                                      "paste(delta^15 * \\"N (\\u2030)\\")" = "d15N"))

#make background to graph all points
TDF_C_N$Diet_Source1<-(TDF_C_N$Diet_Source)
Background<-select(TDF_C_N,-Diet_Source)
head(Background)

TDF_C_N_field<-filter(TDF_C_N,Study=="Field")#field based studies
################################################################################################




################################################################################################
#Regressions
################################################################################################
#Use broom to do regression for combinations of both isotopes and diets
#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

#overall regression for text in manuscript
Overall <- TDF_C_N %>%
  nest(data = -Isotope) %>% 
  mutate(
    fit = map(data, ~ lm(Diet_isotope ~ TDF, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

CN_TDF_overall<-Overall%>% #regressions with p-values and r squared values
  unnest(glanced)

#Each diet source separately for figure
regressions <- TDF_C_N %>%
  nest(data = -Diet_Source, -Isotope) %>% 
  mutate(
    fit = map(data, ~ lm(Diet_isotope ~ TDF, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

CN_TDF<-regressions%>% #regressions with p-values and r squared values
  unnest(glanced)

#R squared
CN_TDF$r.squared<- sprintf('"%.2f"',round(CN_TDF$r.squared, 2)) #round to 2 sig digits and keep trailing zero
CN_TDF$r.squared_graph<- paste("italic(R)^2 ==~", CN_TDF$r.squared) #use for graphing r squared
#P-value
CN_TDF$p.value<- ifelse(CN_TDF$p.value < .0001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(CN_TDF$p.value,digits=3)))#adjust P-value for figure
CN_TDF$p.value_graph<- ifelse(CN_TDF$p.value=="\\"0.0001\\"",paste("italic(P)<",CN_TDF$p.value), paste("italic(P)==",CN_TDF$p.value))#P-value for graphing
head(CN_TDF)  
#Sample size
CN_TDF$n<- paste("italic(n)==",CN_TDF$nobs)# use n-observations for sample size
################################################################################################




################################################################################################
#Graph
################################################################################################
library(ggplot2)
library(lemon)
  ggplot(TDF_C_N, aes(x = Diet_isotope, y = TDF,fill = Diet_Source)) +
#points
  geom_point(data = Background, aes(shape=Consumer_type),size=2.3,alpha = .35,color = "grey",fill = "grey")+#background
  scale_shape_manual(values =c(21,22,24))+#shapes for consumer
  geom_point(aes(shape=Consumer_type,colour = Diet_Source, fill=Diet_Source),size=2.3, alpha=.65)+#colors by diet source for each panel
  geom_point(data = TDF_C_N_field, aes(shape=Consumer_type),size=2.3, fill="NA", stroke = 1)+#put bolded shape over field studies
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3","mediumorchid3"))+#custom colors for outline
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3","mediumorchid3"))+#custom color for fill
#regression line
  geom_smooth(method=lm,color="black",alpha=.25,size=.75)+#Add linear regression line
#Line for commonly used TDF values
  geom_hline(aes(yintercept=Average_line), color='grey50',linetype = "dotted")+
#Error bars
  geom_errorbar(aes(ymin=TDF-SD_TDF, ymax=TDF+SD_TDF), colour="black", width=0,alpha=.35)+#sd for TDF
  geom_errorbarh(aes(xmin=Diet_isotope-SD_Diet_isotope, xmax=Diet_isotope+SD_Diet_isotope), colour="black", height=0,alpha=.35)+#sd for diet  
#Facets
  facet_rep_grid(Diet_Source ~ Isotope, labeller = "label_parsed", scales='free_x') + coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+#remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#removes gridlines, but keeps border
  theme(panel.grid.major.x = element_blank())+#remove facets
  theme(panel.spacing.x = unit(1.5, "lines"))+#add space between panels
  theme(strip.text.y = element_blank())+
  theme(strip.text.x = element_text(size=10,margin = margin(.04,0,.04,0, "cm")))+ 
  theme(strip.background = element_rect(color=NA))+#remove box around text
#labels for regression results
  #n
  geom_text(data=CN_TDF,aes(x=-Inf, y= 0, hjust = -.4, label=(n)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color = "gray50", alpha = .65)+#size and colors
  #p-value
  geom_text(data=CN_TDF,aes(x=-Inf,y = -1, hjust = -.2, label=(p.value_graph)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color = "gray50", alpha = .65)+#size and colors
  #r-squared
  geom_text(data=CN_TDF,aes(x=-Inf,y= -2, hjust = -.2, label=(r.squared_graph)), parse=T,#location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color = "gray50", alpha = .65)+#size and colors
#legend
  guides(color=FALSE,fill=FALSE)+
  theme(legend.title=element_blank(),
        legend.position = c(0.92, 0.05),
        legend.background=element_blank(),
        legend.text=element_text(size=7),
        legend.key.size = unit(.35, "cm"))+
#Labels
  geom_text(aes(x = Inf, y = Inf, vjust=1.3, hjust="inward",#diet source labels on d15N panel
                label = paste(Diet_Source_label), color=Diet_Source), parse = TRUE, size=3.7)+
  geom_text(aes(x = -Inf, y = Inf, vjust=1.2, hjust=-.2,label = Letters), color = "black", size=3.7)+#letters
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-4, 9),breaks=c(-4,0,4,8))+#custom scale for y-axis
  xlab("Isotope value of diet")+#x axis label
  ylab("Trophic discrimination factor (TDF)")+
  theme(axis.title.x= element_text(size=12),
        axis.text.x  = element_text(angle=0, hjust = 1,size=9, colour="black"))+
  theme(axis.title.y = element_text(size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=9, colour="black"))

ggsave("Figure_4_TDF_Diet.tiff",#save figure
       plot = last_plot(), width =5, height = 6, units = "in",
       dpi = 600) 
################################################################################################






