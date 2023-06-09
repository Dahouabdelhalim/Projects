#Ryan Stephens March 7, 2021
rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/Field_data")#set working directory


##################################################################################
#Get data ready
##################################################################################
Isotopes<- read.csv("Field_study_data.csv",header=T)
Isotopes$Date<-as.Date(Isotopes$Date)#make date a date
head(Isotopes)
str(Isotopes)

library(dplyr)
Isotopes$Species_Abr<-Isotopes$Species#make species abbreviation column
library(forcats)#use forcats to rename species levels
Isotopes<-mutate(Isotopes, Species = fct_recode(Species,"B. brevicauda" = "BLBR",
                                                            "M. gapperi" = "MYGA",
                                                            "N. insignis" = "NAIN",  
                                                            "P. maniculatus" = "PEMA"))
#Calculate Diet from stomach contents
Diet<-filter(Isotopes, Type=="Stomach")#select stomach data
Diet <- Diet %>%  mutate(d15N=ifelse(Species_Abr== "MYGA",d15N-1.85,d15N-1))#subtract 1 to account digestive enzymes of stomach (1.85 for MGYA)
Diet<-mutate(Diet, Type = "Diet")#code 'type' as diet

Isotopes_Diet<-rbind(Isotopes,Diet)#append diet to isotope data
head(Isotopes_Diet)

#Use only hair and muscle samples from September to avoid a diet switch from red maple to ECM fungi in M. gapperi 
Isotopes_Diet<-mutate(Isotopes_Diet, Remove = (ifelse(Species_Abr == "MYGA" & Type == "Hair" & DOY < 244, "remove",
                                               ifelse(Species_Abr == "MYGA" & Type == "Muscle" & DOY < 244, "remove", ""))))  
##################################################################################                                                     
                                                 



##################################################################################
#Field-derived TDFs using a resampling approach
##################################################################################
library(dplyr)
library(tidyr)
library(data.table)

##############################
#Hair
##############################
Diet_Hair<-filter(Isotopes_Diet, Type == "Hair" | Type == "Diet")#filter hair and Diet
Diet_Hair<-filter(Diet_Hair, Remove == "")
#resample three replicates 1500 times
Resampled_Hair<-replicate(1500, Diet_Hair %>%#Specify number of times to resample a given dataset
                       group_by(Species, Species_Abr, Type) %>%#Resample within group(s) 
                       sample_n(3), simplify = FALSE) %>%#Specify number of samples to take
                       bind_rows(.id = "ID" )%>%as.data.frame#bind resample data and name with a replicate number using data.table
head(Resampled_Hair)

#get mean value for the three selected replicates
Resampled_long_Hair<-Resampled_Hair%>%
  group_by(Species, Species_Abr, ID, Type)%>%
  summarise(d15N=mean(d15N), d13C=mean(d13C))
head(Resampled_long_Hair)

#Rearrange data
Resampled_group_Hair <- gather(Resampled_long_Hair, Isotope, Value, d15N:d13C, factor_key=TRUE)#put data in long format
head(Resampled_group_Hair)

Resampled_wide_Hair<-spread(Resampled_group_Hair, Type, Value)#put back in wide with a diet and hair column
head(Resampled_wide_Hair)

Resampled_wide_Hair<-mutate(Resampled_wide_Hair, TDF = Hair - Diet)#calculate offset for each iteration

#get mean and sd
Resampled_Hair_TDF<-Resampled_wide_Hair%>%
  group_by(Species, Species_Abr, Isotope)%>%
  summarise(TDF_Mean=round(mean(TDF),2),TDF_SD=round(sd(TDF),2))%>%
  mutate(Method = "Field", Tissue = "Hair")


Resampled_Hair_TDF<-select(Resampled_Hair_TDF, Species,  Species_Abr, Method, Isotope, TDF_Mean, TDF_SD)
##############################


##############################
#Muscle
##############################
Diet_Muscle<-filter(Isotopes_Diet, Type == "Muscle" | Type == "Diet")#filter Muscle and Diet
Diet_Muscle<-filter(Diet_Muscle, Remove == "")

#Resample three replicates 1500 times
Resampled_Muscle<-replicate(1500, Diet_Muscle %>%#Specify number of times to resample a given dataset
                            group_by(Species, Species_Abr, Type) %>%#Resample within group(s) 
                            sample_n(3), simplify = FALSE) %>%#Specify number of samples to take
                            bind_rows(.id = "ID" )%>%as.data.frame#bind resample data and name with a replicate number using data.table
head(Resampled_Muscle)

#get mean value for the three selected replicates
Resampled_long_Muscle<-Resampled_Muscle%>%
  group_by(Species, Species_Abr, ID, Type)%>%
  summarise(d15N=mean(d15N), d13C=mean(d13C))
head(Resampled_long_Muscle)

#Rearange data
Resampled_group_Muscle <- gather(Resampled_long_Muscle, Isotope, Value, d15N:d13C, factor_key=TRUE)#put data in long format
head(Resampled_group_Muscle)

Resampled_wide_Muscle<-spread(Resampled_group_Muscle, Type, Value)#put back in wide with a diet and hair colum
head(Resampled_wide_Muscle)

Resampled_wide_Muscle<-mutate(Resampled_wide_Muscle, TDF = Muscle - Diet)#calculate offset for each iteration
#get mean and sd
Resampled_Muscle_TDF<-Resampled_wide_Muscle%>%
  group_by(Species, Species_Abr, Isotope)%>%
  summarise(TDF_Mean=round(mean(TDF),2),TDF_SD=round(sd(TDF),2))%>%
  mutate(Method = "Field", Tissue = "Muscle")


Resampled_Muscle_TDF<-select(Resampled_Muscle_TDF, Species,  Species_Abr, Method, Isotope, TDF_Mean, TDF_SD)
##############################


##############################
#Merge hair and muscle
##############################
Hair_Muslce_TDF<-rbind(Resampled_Hair_TDF, Resampled_Muscle_TDF)
Hair_Muslce_TDF$Species<- factor(Hair_Muslce_TDF$Species, levels = c("M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#make stomach come first for graphing
Hair_Muslce_TDF<-arrange(Hair_Muslce_TDF, Isotope, Tissue, Species)

write.csv(Hair_Muslce_TDF, "Field_derived_Hair_Muscle_TDF_resampled.csv", row.names = F)
##############################

##################################################################################




##################################################################################
#Hair & diet d13C vs time graph for Appendix
##################################################################################

##############################
#d13C vs time graph
##############################
#Use broom to do regression for combinations of both species and tissue type
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

Stomach_Hair<-filter(Isotopes_Diet, Type == "Hair" | Type == "Stomach")#filter Stomach and Hair
Stomach_Hair$Type <- factor(Stomach_Hair$Type, levels = c("Stomach", "Hair"))#make stomach come first for graphing
#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
regressions_C_Hair <- Stomach_Hair %>%
  nest(data = -Type, -Species) %>% 
  mutate(
    fit = map(data, ~ lm(DOY ~ d13C, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

regressions_C_Hair %>% 
  unnest(tidied)

d13C_regressions_Hair<-regressions_C_Hair%>% #regressions with p-values and r squared values
  unnest(glanced)

#R squared
d13C_regressions_Hair$r.squared<- sprintf('"%.2f"',round(d13C_regressions_Hair$r.squared, 2)) #round to 2 sig digits and keep trailing zero
d13C_regressions_Hair$r.squared_graph<- paste("italic(R)^2 ==~", d13C_regressions_Hair$r.squared) #use for graphing r squared
#P-value
d13C_regressions_Hair$p.value<- ifelse(d13C_regressions_Hair$p.value < .001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(d13C_regressions_Hair$p.value,digits=3)))#adjust P-value for table
d13C_regressions_Hair$p.value_graph<- ifelse(d13C_regressions_Hair$p.value=="\\"0.0001\\"",paste("italic(P)<",d13C_regressions_Hair$p.value), paste("italic(P)==",d13C_regressions_Hair$p.value))#P-value for graphing
head(d13C_regressions_Hair)  


#Graph
head(Stomach_Hair)
library(lemon)
library(ggplot2)

d13C_plot_Hair<-ggplot(Stomach_Hair, aes(x=Date, y=d13C, color = Type, fill = Type)) +
#Points and regressions
  geom_smooth(aes(fill = Type), method='lm', alpha=.25)+#regression line and se
  geom_point(shape= 21, size =2, alpha = .65)+#points
  scale_color_manual(values =c("gray30","orange"),name="Tissue type")+ #custom outline color colors
  scale_fill_manual(values =c("gray30","orange"),name="Tissue type")+ #custom outline color colors
#background and facets
  facet_rep_grid(Species~.)+ #"lemon" package for facets (parsed for isotopes labels of facets)
  coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+ #remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #removes gridlines, but keeps border
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "italic", size = 12))+#strip text for facets
#labels for regression resutls
  #p-value
  geom_text(data=d13C_regressions_Hair,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(-21, -25, -22, -26, -22, -28, -21, -26), 
                hjust = 0, label=(p.value_graph)), parse=T,
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
  #r-squared
  geom_text(data=d13C_regressions_Hair,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(-22.5, -26.5, -23.5, -27.5, -23.5, -29.5, -22.5, -27.5), 
                hjust = 0, label=(r.squared_graph)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
#legend
  theme(legend.position="bottom")+
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-33, -19),breaks=c(-32, -28, -24, -20))+#set y-axis scale and breaks
  theme(panel.spacing.y = unit(0, "lines"))+#decrease space between y-axis facets
  scale_x_date(limits = as.Date(c("2015-07-15","2015-09-17")))+#set x-axis scale
  ylab(expression(delta^13 * "C (\\211)"))+
  xlab("Collection date")+ #x axis title
  theme(axis.title.x = element_text(colour="black", size=12),
        axis.text.x  = element_text(angle=0, hjust = .5,size=8,colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=12),
        axis.text.y  = element_text(angle=0, hjust = .5,size=8,colour="black"))
d13C_plot_Hair
##############################


##############################
#d15N vs time graph
##############################
#Use broom to do regression for combinations of both species and tissue type
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

Stomach_Hair<-filter(Isotopes_Diet, Type == "Hair" | Type == "Stomach")#filter Stomach and Hair
Stomach_Hair$Type <- factor(Stomach_Hair$Type, levels = c("Stomach", "Hair"))#make stomach come first for graphing
#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
regressions_N_Hair <- Stomach_Hair %>%
  nest(data = -Type, -Species) %>% 
  mutate(
    fit = map(data, ~ lm(DOY ~ d15N, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

regressions_N_Hair %>% 
  unnest(tidied)

d15N_regressions_Hair<-regressions_N_Hair%>% #regressions with p-values and r squared values
  unnest(glanced)

#R squared
d15N_regressions_Hair$r.squared<- sprintf('"%.2f"',round(d15N_regressions_Hair$r.squared, 2)) #round to 2 sig digits and keep trailing zero
d15N_regressions_Hair$r.squared_graph<- paste("italic(R)^2 ==~", d15N_regressions_Hair$r.squared) #use for graphing r squared
#P-value
d15N_regressions_Hair$p.value<- ifelse(d15N_regressions_Hair$p.value < .001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(d15N_regressions_Hair$p.value,digits=3)))#adjust P-value for table
d15N_regressions_Hair$p.value_graph<- ifelse(d15N_regressions_Hair$p.value=="\\"0.0001\\"",paste("italic(P)<",d15N_regressions_Hair$p.value), paste("italic(P)==",d15N_regressions_Hair$p.value))#P-value for graphing
head(d15N_regressions_Hair)  


#Graph
head(Stomach_hair)
library(lemon)
library(ggplot2)

d15N_plot_Hair<-ggplot(Stomach_Hair, aes(x=Date, y=d15N, color = Type, fill = Type)) +
#Points and regressions
  geom_smooth(aes(fill = Type), method='lm', alpha=.25)+#regression line and se
  geom_point(shape= 21, size =2, alpha = .65)+#points
  scale_color_manual(values =c("gray30","orange"))+ #custom outline color colors
  scale_fill_manual(values =c("gray30","orange"))+ #custom outline color colors
#background and facets
  facet_rep_grid(Species~.)+ #"lemon" package for facets (parsed for isotopes labels of facets)
  coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+ #remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #removes gridlines, but keeps border
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "italic", size = 12))+#strip text for facets
#labels for regression results
  #p-value
  geom_text(data=d15N_regressions_Hair,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(8, 3, 9, 4, 5, 0, 4, -1), 
                hjust = 0, label=(p.value_graph)), parse=T,
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
#r-squared
  geom_text(data=d15N_regressions_Hair,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(6, 1, 7, 2, 3, -2, 2, -3), 
                hjust = 0, label=(r.squared_graph)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
#legend
  theme(legend.position="none")+
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-4, 12),breaks=c(-4, 0, 4, 8, 12))+#set y-axis scale and breaks
  theme(panel.spacing.y = unit(-0, "lines"))+#decrease space between y-axis facets
  scale_x_date(limits = as.Date(c("2015-07-15","2015-09-17")))+#set x-axis scale
  ylab(expression(delta^15 * "N (\\211)"))+
  xlab("Collection date")+ #x axis title
  theme(axis.title.x = element_text(colour="black", size=12),
        axis.text.x  = element_text(angle=0, hjust = .5,size=8,colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=12),
        axis.text.y  = element_text(angle=0, hjust = .5,size=8,colour="black"))
d15N_plot_Hair
##############################


##############################
#Merge d13C & d15N
##############################
library(cowplot)
library(rlang)
# extract a legend that is laid out horizontally
Legend <- get_legend(
  d13C_plot_Hair + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

ggdraw() +
  draw_plot(d13C_plot_Hair+theme(legend.position="none"), x = 0, y = .1, width = .5, height = .9)+#placement and remove legend
  draw_label("(a)", x = 0.03, y = .98, size = 12, fontface = "bold", angle =0, vjust = .5, color="black")+
  draw_plot(d15N_plot_Hair, x = .5, y = .1, width = .5, height = .9)+#placement and remove legend
  draw_label("(b)", x = .53, y = .98, size = 12, fontface = "bold", angle =0, vjust = .5, color="black")+
  draw_plot(Legend, x = .25, y = 0, width = .5, height = .1)#placement of legend

ggsave("Hair_Stomach_d15N_d13C_time.jpeg", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =6, height = 6, units = "in",
       dpi = 600) 
##############################

##################################################################################




##################################################################################
#Muscle & diet d13C vs time graph
##################################################################################

##############################
#d13C vs time graph
##############################
#Use broom to do regression for combinations of both species and tissue type
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

Stomach_Muscle<-filter(Isotopes_Diet, Type == "Muscle" | Type == "Stomach")#filter Stomach and Muscle
Stomach_Muscle$Type <- factor(Stomach_Muscle$Type, levels = c("Stomach", "Muscle"))#make stomach come first for graphing
#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
regressions_C_Muscle <- Stomach_Muscle %>%
  nest(data = -Type, -Species) %>% 
  mutate(
    fit = map(data, ~ lm(DOY ~ d13C, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

regressions_C_Muscle %>% 
  unnest(tidied)

d13C_regressions_Muscle<-regressions_C_Muscle%>% #regressions with p-values and r squared values
  unnest(glanced)

#R squared
d13C_regressions_Muscle$r.squared<- sprintf('"%.2f"',round(d13C_regressions_Muscle$r.squared, 2)) #round to 2 sig digits and keep trailing zero
d13C_regressions_Muscle$r.squared_graph<- paste("italic(R)^2 ==~", d13C_regressions_Muscle$r.squared) #use for graphing r squared
#P-value
d13C_regressions_Muscle$p.value<- ifelse(d13C_regressions_Muscle$p.value < .001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(d13C_regressions_Muscle$p.value,digits=3)))#adjust P-value for table
d13C_regressions_Muscle$p.value_graph<- ifelse(d13C_regressions_Muscle$p.value=="\\"0.0001\\"",paste("italic(P)<",d13C_regressions_Muscle$p.value), paste("italic(P)==",d13C_regressions_Muscle$p.value))#P-value for graphing
head(d13C_regressions_Muscle)  


#Graph
head(Stomach_Muscle)
library(lemon)
library(ggplot2)

d13C_plot_Muscle<-ggplot(Stomach_Muscle, aes(x=Date, y=d13C, color = Type, fill = Type)) +
#Points and regressions
  geom_smooth(aes(fill = Type), method='lm', alpha=.25)+#regression line and se
  geom_point(shape= 21, size =2, alpha = .65)+#points
  scale_color_manual(values =c("gray30","orange"),name="Tissue type")+ #custom outline color colors
  scale_fill_manual(values =c("gray30","orange"),name="Tissue type")+ #custom outline color colors
  #background and facets
  facet_rep_grid(Species~.)+ #"lemon" package for facets (parsed for isotopes labels of facets)
  coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+ #remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #removes gridlines, but keeps border
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "italic", size = 12))+#strip text for facets
#labels for regression results
  #p-value
  geom_text(data=d13C_regressions_Muscle,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(-21, -25, -22, -26, -22, -28, -21, -26), 
                hjust = 0, label=(p.value_graph)), parse=T,
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
  #r-squared
  geom_text(data=d13C_regressions_Muscle,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(-22.5, -26.5, -23.5, -27.5, -23.5, -29.5, -22.5, -27.5), 
                hjust = 0, label=(r.squared_graph)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
#legend
  theme(legend.position="bottom")+
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-33, -19),breaks=c(-32, -28, -24, -20))+#set y-axis scale and breaks
  theme(panel.spacing.y = unit(0, "lines"))+#decrease space between y-axis facets
  scale_x_date(limits = as.Date(c("2015-07-15","2015-09-17")))+#set x-axis scale
  ylab(expression(delta^13 * "C (\\211)"))+
  xlab("Collection date")+ #x axis title
  theme(axis.title.x = element_text(colour="black", size=12),
        axis.text.x  = element_text(angle=0, hjust = .5,size=8,colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=12),
        axis.text.y  = element_text(angle=0, hjust = .5,size=8,colour="black"))
d13C_plot_Muscle
##############################


##############################
#d15N vs time graph
##############################
#Use broom to do regression for combinations of both species and tissue type
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

Stomach_Muscle<-filter(Isotopes_Diet, Type == "Muscle" | Type == "Stomach")#filter Stomach and Muscle
Stomach_Muscle$Type <- factor(Stomach_Muscle$Type, levels = c("Stomach", "Muscle"))#make stomach come first for graphing
#https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
regressions_N_Muscle <- Stomach_Muscle %>%
  nest(data = -Type, -Species) %>% 
  mutate(
    fit = map(data, ~ lm(DOY ~ d15N, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

regressions_N_Muscle %>% 
  unnest(tidied)

d15N_regressions_Muscle<-regressions_N_Muscle%>% #regressions with p-values and r squared values
  unnest(glanced)

#R squared
d15N_regressions_Muscle$r.squared<- sprintf('"%.2f"',round(d15N_regressions_Muscle$r.squared, 2)) #round to 2 sig digits and keep trailing zero
d15N_regressions_Muscle$r.squared_graph<- paste("italic(R)^2 ==~", d15N_regressions_Muscle$r.squared) #use for graphing r squared
#P-value
d15N_regressions_Muscle$p.value<- ifelse(d15N_regressions_Muscle$p.value < .001, sprintf('"%s"',"0.0001"),sprintf('"%.3f"',round(d15N_regressions_Muscle$p.value,digits=3)))#adjust P-value for table
d15N_regressions_Muscle$p.value_graph<- ifelse(d15N_regressions_Muscle$p.value=="\\"0.0001\\"",paste("italic(P)<",d15N_regressions_Muscle$p.value), paste("italic(P)==",d15N_regressions_Muscle$p.value))#P-value for graphing
head(d15N_regressions_Muscle)  


#Graph
head(Stomach_Muscle)
library(lemon)
library(ggplot2)

d15N_plot_Muscle<-ggplot(Stomach_Muscle, aes(x=Date, y=d15N, color = Type, fill = Type)) +
#Points and regressions
  geom_smooth(aes(fill = Type), method='lm', alpha=.25)+#regression line and se
  geom_point(shape= 21, size =2, alpha = .65)+#points
  scale_color_manual(values =c("gray30","orange"))+ #custom outline color colors
  scale_fill_manual(values =c("gray30","orange"))+ #custom outline color colors
#background and facets
  facet_rep_grid(Species~.)+ #"lemon" package for facets (parsed for isotopes labels of facets)
  coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+ #remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #removes gridlines, but keeps border
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "italic", size = 12))+#strip text for facets
#labels for regression results
  #p-value
  geom_text(data=d15N_regressions_Muscle,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(8, 3, 9, 4, 5, 0, 4, -1), 
                hjust = 0, label=(p.value_graph)), parse=T,
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
  #r-squared
  geom_text(data=d15N_regressions_Muscle,
            aes(x=as.Date(c("2015-08-24","2015-08-24","2015-07-17","2015-07-17","2015-08-24","2015-08-24","2015-07-17","2015-07-17")),
                y= c(6, 1, 7, 2, 3, -2, 2, -3), 
                hjust = 0, label=(r.squared_graph)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=2.7,inherit.aes=F, color =c("orange","gray30","orange","gray30","orange","gray30","orange","gray30"))+ #size and colors
#legend
  theme(legend.position="none")+
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-4, 12),breaks=c(-4, 0, 4, 8, 12))+#set y-axis scale and breaks
  theme(panel.spacing.y = unit(-0, "lines"))+#decrease space between y-axis facets
  scale_x_date(limits = as.Date(c("2015-07-15","2015-09-17")))+#set x-axis scale
  ylab(expression(delta^15 * "N (\\211)"))+
  xlab("Collection date")+ #x axis title
  theme(axis.title.x = element_text(colour="black", size=12),
        axis.text.x  = element_text(angle=0, hjust = .5,size=8,colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=12),
        axis.text.y  = element_text(angle=0, hjust = .5,size=8,colour="black"))
d15N_plot_Muscle
##############################


##############################
#Merge d13C & d15N
##############################
library(cowplot)
library(rlang)
# extract a legend that is laid out horizontally
Legend <- get_legend(
  d13C_plot_Muscle + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

ggdraw() +
  draw_plot(d13C_plot_Muscle+theme(legend.position="none"), x = 0, y = .1, width = .5, height = .9)+#placement and remove legend
  draw_label("(a)", x = 0.03, y = .98, size = 12, fontface = "bold", angle =0, vjust = .5, color="black")+
  draw_plot(d15N_plot_Muscle, x = .5, y = .1, width = .5, height = .9)+#placement and remove legend
  draw_label("(b)", x = .53, y = .98, size = 12, fontface = "bold", angle =0, vjust = .5, color="black")+
  draw_plot(Legend, x = .25, y = 0, width = .5, height = .1)#placement of legend

ggsave("Muscle_Stomach_d15N_d13C_time.jpeg", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =6, height = 6, units = "in",
       dpi = 600) 
##############################

##################################################################################






