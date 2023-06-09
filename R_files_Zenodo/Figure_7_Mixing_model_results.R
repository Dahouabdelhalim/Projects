#Ryan Stephens June 9, 2021
rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure")

##################################################################################################
#Diet, food sources, and hair adjusted by SIBER TDF and Hair TDF
##################################################################################################

##################################################
#Potential food sources
##################################################
Food<- read.csv("Bartlett_diet_items.csv",header=T)#Potential dietary food items at Bartlett Experimental Forest

library(dplyr)
#calculate standard deviations of food sources
Food_SD<-Food %>%
  group_by (Type) %>%
  summarise (Md13C = mean(d13C), Md15N = mean(d15N),  SD_d13C = sd(d13C), SD_d15N = sd(d15N))%>%
  as.data.frame       
Food_SD<-select(Food_SD, Type,d13C=Md13C, d15N=Md15N, SD_d13C, SD_d15N)#change column names

#Locations for food source labels
Food_text <- tibble(
  Species = c("B. brevicauda","B. brevicauda","B. brevicauda","B. brevicauda","B. brevicauda"),
  Type = c("Red Maple", "EM Fungi", "AM Fungi", "Berries", "Arthropods"),
  x     = c(-22.8,-26,-30,-33,-23.8),
  y     = c(-4,11,.7,2.5,7))
Diet_Space = c(" Isospace of food sources & diet", " Isospace of food sources & diet", " Isospace of food sources & diet", " Isospace of food sources & diet"," Isospace of food sources & diet")
Food_text$Species<-as.factor(Food_text$Species)
str(Food_text)
##################################################


##################################################
#small mammal diets (from stomach samples)
##################################################
#Calculate Diet from stomach contents
Field_data<- read.csv("Field_study_data.csv",header=T)#Stomach samples from Bartlett small mammals
Diet<-filter(Field_data, Type=="Stomach")#select stomach data
Diet <- Diet %>% mutate(d15N=ifelse(Species == "MYGA",d15N-1.85,d15N-1))#subtract 1 to account digestive enzymes of stomach (1.85 for MGYA)
Diet<-select(Diet, Species, d15N, d13C)

library(forcats)#use forcats to rename levels
Diet<-mutate(Diet, Species = fct_recode(Species,"B. brevicauda" = "BLBR",
                                                "M. gapperi" = "MYGA",
                                                "N. insignis" = "NAIN",  
                                                "P. maniculatus" = "PEMA"))
Diet<-mutate(Diet, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder
Diet$TDF<-"Diet"
Diet$Diet_Space<-" Isospace of food sources & diet"
head(Diet)

###############
#Diet contour
###############
#function to get contour of confidence interval
#https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in
Interval_95<-function(n)
{xy <- select(n, x = d13C, y = d15N)
 kd <- ks::kde(xy, compute.cont=TRUE)
 contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
 contour_95 <- data.frame(contour_95)               
 contour_95 <- select(contour_95, d13C=x, d15N=y) 
}

BLBR_Diet<-filter(Diet, Species == "B. brevicauda")#B. brevicauda
BLBR_Diet_Interval<-Interval_95(BLBR_Diet)
BLBR_Diet_Interval$Species<-"B. brevicauda"
MYGA_Diet<-filter(Diet, Species == "M. gapperi")#M. gapperi
MYGA_Diet_Interval<-Interval_95(MYGA_Diet)
MYGA_Diet_Interval$Species<-"M. gapperi"
NAIN_Diet<-filter(Diet, Species == "N. insignis")#N. insignis
NAIN_Diet_Interval<-Interval_95(NAIN_Diet)
NAIN_Diet_Interval$Species<-"N. insignis"
PEMA_Diet<-filter(Diet, Species == "P. maniculatus")#P. maniculatus
PEMA_Diet_Interval<-Interval_95(PEMA_Diet)
PEMA_Diet_Interval$Species<-"P. maniculatus"
Diet_Contour<-rbind(BLBR_Diet_Interval, MYGA_Diet_Interval, NAIN_Diet_Interval, PEMA_Diet_Interval)
Diet_Contour$Diet_Space<-" Isospace of food sources & diet"
Diet_Contour$Species<-as.factor(Diet_Contour$Species)
Diet_Contour<-mutate(Diet_Contour, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder
levels(Diet_Contour$Species)
###############
##################################################


##################################################
#Hair samples
##################################################
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)#hair samples from live trapping
Hair<-mutate(Hair, Species = fct_recode(Species,"B. brevicauda" = "BLBR",
                                        "M. gapperi" = "MYGA",
                                        "N. insignis" = "NAIN",  
                                        "P. maniculatus" = "PEMA"))
head(Hair)
str(Hair)

###############
#SIDER adjusted TDF
###############
TDF<- read.csv("Bartlett_estimated_TDF.csv",header=T)#estimated TDF dataset
TDF_SIDER<-filter(TDF, Method == "SIDER")#select SIDER TDF
TDF_SIDER<-select(TDF_SIDER, Species, TDF_Mean, Isotope)
library(tidyr)
TDF_SIDER_wide <- spread(TDF_SIDER, Isotope, TDF_Mean)#put data in wide format
TDF_SIDER_wide<-rename(TDF_SIDER_wide, TDF_d13C = d13C, TDF_d15N = d15N)#change column names
head(TDF_SIDER_wide)
SIDER_Hair<-left_join(Hair, TDF_SIDER_wide, by = "Species")#Join TDF data to hair data
SIDER_Hair <- SIDER_Hair %>%  mutate(adj_d13C=d13C-TDF_d13C)#calculate adjusted d13C value for hair
SIDER_Hair <- SIDER_Hair %>%  mutate(adj_d15N=d15N-TDF_d15N)#calculate adjusted d15N value for hair
SIDER<-select(SIDER_Hair, Species, d15N=adj_d15N, d13C=adj_d13C, ID)
SIDER$TDF<-"SIDER"
SIDER$Diet_Space<-" Isospace of food sources & diet"
head(SIDER)
SIDER<-mutate(SIDER, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder

#SIDER contours
BLBR_SIDER<-filter(SIDER, Species == "B. brevicauda" & ID != "1892")#B. brevicauda (exclude sample 1892 from interval for one contour line)
BLBR_SIDER_Interval<-Interval_95(BLBR_SIDER)
BLBR_SIDER_Interval$Species<-"B. brevicauda"
MYGA_SIDER<-filter(SIDER, Species == "M. gapperi")#M. gapperi
MYGA_SIDER_Interval<-Interval_95(MYGA_SIDER)
MYGA_SIDER_Interval$Species<-"M. gapperi"
NAIN_SIDER<-filter(SIDER, Species == "N. insignis")#N. insignis
NAIN_SIDER_Interval<-Interval_95(NAIN_SIDER)
NAIN_SIDER_Interval$Species<-"N. insignis"
PEMA_SIDER<-filter(SIDER, Species == "P. maniculatus")#P. maniculatus
PEMA_SIDER_Interval<-Interval_95(PEMA_SIDER)
PEMA_SIDER_Interval$Species<-"P. maniculatus"
SIDER_Contour<-rbind(BLBR_SIDER_Interval, MYGA_SIDER_Interval, NAIN_SIDER_Interval, PEMA_SIDER_Interval)
SIDER_Contour$SIDER_Space<-" Isospace of food sources & SIDER"
SIDER_Contour$Species<-as.factor(SIDER_Contour$Species)
SIDER_Contour<-mutate(SIDER_Contour, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder
###############


###############
#Field adjusted TDF
###############
TDF<- read.csv("Bartlett_estimated_TDF.csv",header=T)#estimated TDF dataset
TDF<-mutate(TDF, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder

TDF_Field<-filter(TDF, Method == "Field")#select Field TDF
TDF_Field<-select(TDF_Field, Species, TDF_Mean, Isotope)
library(tidyr)
TDF_Field_wide <- spread(TDF_Field, Isotope, TDF_Mean)#put data in wide format
TDF_Field_wide<-rename(TDF_Field_wide, TDF_d13C = d13C, TDF_d15N = d15N)#change column names
head(TDF_Field_wide)
Field_Hair<-left_join(Hair, TDF_Field_wide, by = "Species")#Join TDF data to hair data
Field_Hair <- Field_Hair %>%  mutate(adj_d13C=d13C-TDF_d13C)#calculate adjusted d13C value for hair
Field_Hair <- Field_Hair %>%  mutate(adj_d15N=d15N-TDF_d15N)#calculate adjusted d15N value for hair
Field<-select(Field_Hair, Species, d15N=adj_d15N, d13C=adj_d13C, ID)
Field$TDF<-"Field"
Field$Diet_Space<-" Isospace of food sources & diet"
head(Field)

#Field contours
BLBR_Field<-filter(Field, Species == "B. brevicauda" & ID != "1892")#B. brevicauda (exclude sample 1892 from interval for one contour line)
BLBR_Field_Interval<-Interval_95(BLBR_Field)
BLBR_Field_Interval$Species<-"B. brevicauda"
MYGA_Field<-filter(Field, Species == "M. gapperi")#M. gapperi
MYGA_Field_Interval<-Interval_95(MYGA_Field)
MYGA_Field_Interval$Species<-"M. gapperi"
NAIN_Field<-filter(Field, Species == "N. insignis")#N. insignis
NAIN_Field_Interval<-Interval_95(NAIN_Field)
NAIN_Field_Interval$Species<-"N. insignis"
PEMA_Field<-filter(Field, Species == "P. maniculatus")#P. maniculatus
PEMA_Field_Interval<-Interval_95(PEMA_Field)
PEMA_Field_Interval$Species<-"P. maniculatus"
Field_Contour<-rbind(BLBR_Field_Interval, MYGA_Field_Interval, NAIN_Field_Interval, PEMA_Field_Interval)
Field_Contour$Field_Space<-" Isospace of food sources & Field"
Field_Contour$Species<-as.factor(Field_Contour$Species)
Field_Contour<-mutate(Field_Contour, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder
###############
##################################################


##################################################
#Combine data sets for graphing
##################################################
SIDER<-select(SIDER, -ID)
Field<-select(Field, -ID)
Diet_SIDER_Field<-rbind(Diet, SIDER, Field)
head(Diet_SIDER_Field)
str(Diet_SIDER_Field)
Diet_SIDER_Field<-mutate(Diet_SIDER_Field, Diet =ifelse(Species == "M. gapperi", "Herbivore", #text for diet class
                                                 ifelse(Species == "B. brevicauda", "Carnivore", "Omnivore")))
Diet_SIDER_Field$TDF<-as.factor(Diet_SIDER_Field$TDF)
##################################################


##################################################
#Plot
##################################################
library(ggplot2)
library(lemon)
library(cowplot)

Diet_Plot<-
  ggplot(Diet_SIDER_Field, aes(x=d13C, y=d15N))+
#Base plot
  geom_polygon(data = Diet_Contour, size=.6, color="gray85", fill = "gray90")+#diet polygons
  facet_rep_grid(Species~Diet_Space, switch = "y") + coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw()+
  theme_cowplot(12)+#make background transparent
  theme(panel.border=element_blank(), axis.line=element_line(size=.3))+
  scale_fill_manual(values = c("gray30", "orange","mediumorchid3"))+#fill of TDfs
#Food sources
  geom_point(data=Food,aes(x=d13C, y=d15N,color=Type,shape=Type),alpha = 0.4,size=.7)+#add food sources
  scale_shape_manual(values = c(15, 17, 18, 15, 19))+#shape of food sources
  scale_color_manual(values = c("cyan4","gray40","coral3","Tan4","olivedrab4"))+#outline of food sources 
  guides(color=FALSE)+
  geom_errorbarh(data=Food_SD,aes(xmin = d13C - SD_d13C,xmax = d13C + SD_d13C,y = d15N,color = Type, height = 0), size=.6,alpha = 1)+ #error bars
  geom_errorbar(data=Food_SD,aes(ymin = d15N - SD_d15N, ymax = d15N + SD_d15N,x = d13C,color = Type, width = 0), size=.6,alpha = 1)+ #error bars
  geom_text(data = Food_text,mapping = aes(x = x, y = y, label = Type,color=Type),size=2,alpha = 0.7)+#text for food sources
  geom_text(aes(x = -33, y = 6,label = Diet), color="gray50", size=2)+#text for diet class
#Overlay hair and stomach stomach samples points
  geom_point(aes(fill =TDF), color="white",alpha = .45, size=1.7, shape = 21)+#diet samples, color = "white"
  geom_point(data= Diet, fill="gray30",color="gray30", size=1.7, shape = 21)+#diet samples, color = "white"
#overlay polygons
  #Field TDF adjusted hair samples
  geom_polygon(data = Field_Contour, size=.6, fill=NA, color="orange")+
  geom_point(data=Field,aes(x=d13C, y=d15N), alpha = .1, size=1.7, colour= "orange")+  
  #SIDER TDF adjusted hair samples  
  geom_polygon(data = SIDER_Contour, size=.6, fill=NA, color="mediumorchid3")+
  geom_point(data=SIDER,aes(x=d13C, y=d15N), alpha = .1, size=1.7, colour= "mediumorchid3")+ 
#legend
  guides(shape = FALSE)+#remove food sources
  theme(legend.title=element_blank(),#remove title
        legend.box.margin = margin(.5, .5, .5, .5, "mm"),#set margin limits
        legend.box.background = element_rect(size = .2, colour = "gray35"),#Size and color of border
        legend.position = c(0.78, 0.97),#legend position
        legend.background=element_blank(),#background color of legend
        legend.text=element_text(size=6, color ="gray35"),#text color
        legend.spacing.x = unit(0, 'cm'),#horizontal spacing
        legend.spacing.y = unit(0, "mm"), 
        legend.key.size = unit(.25, "cm"))+#vertical spacing
#Formats the text of the labels
  theme(strip.background = element_rect(fill="gray85"))+
  theme(panel.spacing.x=unit(8, "lines"))+
  theme(strip.placement = "outside")+
  theme(strip.text.y = element_blank())+#remove strip text
  theme(strip.text.x = element_text(size =6,margin = margin(.07,0,.07,0, "cm")))+#set margins
#Axes portion
  scale_x_continuous(limits = c(-36.3, -20),breaks=c(-36,-32,-28,-24,-20), expand = c(0, 0))+#manually put in scales on y axis
  scale_y_continuous(limits = c(-6, 13),breaks=c(-4,0,4,8,12), expand = c(0, 0))+#manually put in scales on y axis
  xlab(expression(delta^13 * "C (\\211)"))+
  ylab(expression(delta^15 * "N (\\211)"))+
  theme(axis.title.x = element_text(colour="black", size=8),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0, size=6,colour="black"))+#adjusts angles, size, etc of labels for x axis
  theme(axis.title.y = element_text(colour="black", vjust=0, size=8),#adjusts size and color of y axis title
        axis.text.y  = element_text(angle=0, vjust=0.5, size=6, colour="black"))+
    theme(panel.spacing.y=unit(.6, "lines"))

Diet_Plot
##################################################################################################




##################################################################################################
#Results of mixing models for diet and different TDF estimates
##################################################################################################


##################################################
#Mixing model quantiles 
##################################################
#Load all results of mixing models and combine them into one dataframe (note that individual output
#files were saves separately for each model). Data were then combined into one dataframe.
library(dplyr)
library(data.table)
#files <- list.files(path = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure", pattern = "Posterior", full.names = TRUE)
#Mixing_model_output <- bind_rows(lapply(files, fread))
#write.csv(Mixing_model_output, "Mixing_model_output.csv", row.names = F)#data were saved to file 
Mixing_Model_output<- read.csv("Mixing_model_output.csv", header=T)#Load mixing model output
head(Mixing_Model_output)

#Quantile functions
Quant_025 <- function(x) {r <- quantile(x, probs=c(0.025))}
Quant_125 <- function(x) {r <- quantile(x, probs=c(0.125))}
Quant_250 <- function(x) {r <- quantile(x, probs=c(0.250))}
Quant_750 <- function(x) {r <- quantile(x, probs=c(0.750))}
Quant_875 <- function(x) {r <- quantile(x, probs=c(0.875))}
Quant_975 <- function(x) {r <- quantile(x, probs=c(0.975))}

#summarize quantiles for each food Source
Diets_Quantile<-Mixing_Model_output%>%
  group_by(Species,Source, Model)%>%
  summarize(Median = median(value), Mean = mean(value),Q_975=Quant_975(value),Q_875=Quant_875(value),Q_75=Quant_750(value), 
            Q_25=Quant_250(value),Q_125=Quant_125(value),Q_025=Quant_025(value))%>%
  as.data.frame

#Select needed levels of food sources and reorder
Diets_Quantile$Source<-as.factor(Diets_Quantile$Source)
Diets_Quantile<-filter(Diets_Quantile, Source == "Red_Maple" |Source == "EM_Fungi" |Source == "AM_Fungi" |Source == "Berries" |Source == "Arthropods")
Diets_Quantile<-droplevels(Diets_Quantile)
Diets_Quantile$Source<-as.factor(Diets_Quantile$Source)

library(forcats)
Diets_Quantile<-mutate(Diets_Quantile, Source = fct_relevel(Source,"Red_Maple", "EM_Fungi", "AM_Fungi", "Berries", "Arthropods"))#reorder

#Reorder levels of models
Diets_Quantile$Model<-as.factor(Diets_Quantile$Model)
levels(Diets_Quantile$Model)
Diets_Quantile<-mutate(Diets_Quantile, Model = fct_relevel(Model, "Diet", "SIDER", "SIDER_update", "Meta_analysis", "Field"))#reorder       
Diets_Quantile<-mutate(Diets_Quantile, Species_Source_Model = paste(Species, Source, Model, sep="_"))
head(Diets_Quantile)
##################################################


##################################################
#Significant differences between diet and TDF models
##################################################
#Determine if diet is significantly different than model output                                                      
Mixing_Model_output<-mutate(Mixing_Model_output, ID = paste(Number, Source, Species,sep="_"))#make matching ID
Models<-filter(Mixing_Model_output, Model != "Diet")#Models
Diet<-filter(Mixing_Model_output, Model == "Diet")#Diet
Diet_value_ID<-select(Diet, ID, Diet_value = value)#select ID and posterior values
head(Diet_value_ID)
Models_Diet<-left_join(Models, Diet_value_ID, by = "ID")#join diet posterior values to model posterior values
head(Models_Diet)

Models_Diet<-mutate(Models_Diet, Model_Diet = value-Diet_value, Diet_Model = Diet_value-value)#subtract model from diet
Models_Diet<-filter(Models_Diet, Source == "Beech" | Source == "Red_Maple" |Source == "EM_Fungi" |Source == "AM_Fungi" |Source == "Berries" |Source == "Arthropods")

#Proportion of posterior values that overlap (at 95%)
Probability<-Models_Diet%>%
  group_by(Species,Source, Model)%>%
  summarize(P_Model_Diet = sum(Model_Diet < 0)/n())%>%
  as.data.frame()
Probability<-mutate(Probability, Sig = ifelse(P_Model_Diet> 0.95 | P_Model_Diet< 0.05, "**", ""))#determine which ones are significant 
head(Probability)
Probability<-mutate(Probability, Species_Source_Model = paste(Species, Source, Model, sep="_"))#make matching ID
Prob<-select(Probability,Species_Source_Model, Sig)#select ID and Sig column

Diets_Quantile<-left_join(Diets_Quantile, Prob, by = "Species_Source_Model")#join significant column to quantile dataset
Diets_Quantile[is.na(Diets_Quantile)] <- ""#Make NA for diet to blank space

#Determine if model output falls within 75% quanitles of diet estimates  
Diets_Quantile<-mutate(Diets_Quantile, Species_Source = paste(Species, Source, sep="_"))#make ID
Diet_Q<-filter(Diets_Quantile, Model == "Diet")#select only diet
Diet_Q<-select(Diet_Q, Species_Source,  Diet_Q_875 = Q_875, Diet_Q_125 = Q_125)#75% quantiles
Diets_Quantile<-left_join (Diets_Quantile, Diet_Q, by = "Species_Source")#join diet quantiles back to dataset

Diets_Quantile<-mutate(Diets_Quantile, Model_Diet_Range = ifelse(Diet_Q_875 > Mean & Diet_Q_125 < Mean, "", "*"))#determine if model means fall between diet quantiles

#for models with adjusted alpha this comparison doesn't work
Diets_Quantile<-mutate(Diets_Quantile, Model_Diet_Range = ifelse(Species == "MYGA" & Source == "Arthropods", "", Model_Diet_Range))
Diets_Quantile<-mutate(Diets_Quantile, Model_Diet_Range = ifelse(Species == "BLBR" & Source == "EM_Fungi", "", Model_Diet_Range))

head(Diets_Quantile)
Diets_Quantile<-mutate(Diets_Quantile, asterisk = ifelse(Sig == "**", Sig, Model_Diet_Range))#if significant put that otherwise use quantile overlap
##################################################


##################################################
#Plotting 
##################################################
#rename levels
library(forcats)
levels(Diets_Quantile$Model)
Diets_Quantile<-mutate(Diets_Quantile, Model = fct_recode(Model,"Meta-analysis" = "Meta_analysis",
                                                                "SIDER+Meta" = "SIDER_update"))

Diets_Quantile<-mutate(Diets_Quantile, Source = fct_recode(Source,"AM Fungi" = "AM_Fungi",
                                                                  "EM Fungi" = "EM_Fungi",
                                                                  "Red Maple" = "Red_Maple"))

Diets_Quantile<-mutate(Diets_Quantile, Species = fct_recode(Species,"B. brevicauda" = "BLBR",
                                                                    "M. gapperi" = "MYGA",
                                                                    "N. insignis" = "NAIN",  
                                                                    "P. maniculatus" = "PEMA"))

Diets_Quantile<-mutate(Diets_Quantile, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus", "B. brevicauda"))#reorder
head(Diets_Quantile)
levels(Diets_Quantile$Species)

#make data frame for labeling facets without models
Species<-c("P. maniculatus", "P. maniculatus")
Model<-c("SIDER", "SIDER+Meta")
Source<-c("AM Fungi", "AM Fungi")
Median<-c(.5, .5)
No_data_label<-data.frame(Species, Model, Source, Median)
No_data_label$Label<- "TDF ouside of\\ndiatary mixing\\nspace"
No_data_label<-mutate(No_data_label, Model = fct_relevel(Model,"SIDER", "SIDER+Meta"))#reorder
No_data_label$Species<-as.factor(No_data_label$Species)
str(No_data_label)

library(lemon)
library(ggplot2)
library(grid)

Mixing_Model_Plot<- 
  ggplot(Diets_Quantile, aes(x=Source, y=Median, fill=Source)) + 
#crossbars to represent posterior quantiles
  geom_crossbar(aes(ymin = Q_025, ymax = Q_975, color= Source), size=.3, width = 0.05)+# 95% probability of data distribution
  geom_crossbar(aes(ymin = Q_125, ymax = Q_875, color= Source), size=.3, width = 0.3)+# 75% probability of data distribution
  geom_crossbar(aes(ymin = Q_25, ymax = Q_75, color= Source), size=.3, fatten = 2, width = 0.6)+# 50% probability of data distribution
  scale_color_manual(values = c("olivedrab4","Tan4","cyan4","coral3","gray40"))+#colors of food source outline  
  scale_fill_manual(values = c("#C2D2A4","#DAC3AB","#9ED3D7","#FFA7A2","#AAAAAA"))+#colors of food source fill  
#means
  geom_point(aes(y=Mean),shape=21,fill = "white",size=.9)+#Means plotted as white dot
#facets
  facet_rep_grid(Species~Model) + coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+
#text
  geom_text(aes(x = Source, y = Q_975, label = asterisk), vjust=.2,size=3.5, color="black")+#add text for significance
  geom_text(data = No_data_label, aes(x = Source, y = Median, label = Label),inherit.aes = FALSE, size = 2, color = "gray65")+
#Legend
  theme(legend.position = "none")+
#background and axes
  theme(panel.spacing.y =unit(-.7, "lines"))+
  theme(panel.spacing.x =unit(-.2, "lines"))+
  scale_y_continuous(limits = c(0, 1),breaks=c(0,.2,.4,.6,.8, 1), expand = c(0, 0))+#manually put in scales on y axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#removes gridlines, but keeps border
  theme(panel.grid.major.x = element_blank())+
  theme(strip.background.y = element_blank())+#remove background of y axis text
  theme(strip.text.y = element_blank())+
  theme(strip.text.x = element_text(color = "white", size =6,margin = margin(.05,0,.05,0, "cm")))+#white text and smaller box hight
  ylab("Proportion of dietary food sources")+
  theme(axis.title.x = element_blank(),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=45, hjust = 1,size=6, face="italic",colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=8),#adjusts size and color of y axis title
        axis.text.y  = element_text(angle=0, vjust=0.5, size=6, colour="black"))

Mixing_Model_Plot

#change color of top panel
#https://github.com/tidyverse/ggplot2/issues/2096
library(tidyverse)
library(grid)

Mixing_Model_Plot#plot
g <- ggplot_gtable(ggplot_build(Mixing_Model_Plot))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("gray30","mediumorchid3","darkmagenta", "chartreuse4","orange")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
##################################################
##################################################################################################




##################################################################################################
#Combine diet and mixing model plots
##################################################################################################
library(cowplot)
ggdraw() +
#Small mammal illustrations
  draw_image("~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/Southern_Redbacked_Vole_White.jpg",  x = -.36, y = 0.425, scale = .15) +
  draw_image("~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/Woodland_Jumping_Mouse_White.jpg", x = -.35, y = 0.20, scale = .15) +
  draw_image("~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/Woodland_Deer_Mouse_White.jpg",   x = -.35, y = -.032, scale = .15) +
  draw_image("~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/Northern_short-tailed_shrew_white.jpg",  x = -.35, y = -0.27, scale = .16) +
#Diet, food sources, and hair adjusted by SIBER TDF and Hair TDF
  draw_plot(Diet_Plot, x = 0, y = .008,#X and y coordinates to start plot (range from 0 to 1)
            width = .4, height = .995)+
#Results of mixing models for diet and different TDF estimates
  draw_plot(g, x = .4, y = 0,#X and y coordinates to start plot (range from 0 to 1)
            width = .6, height = 1)+
#letters
  draw_label("(a)",x = .09, y = .96, size = 10, hjust = 0, color="black")+
  draw_label("(b)",x = .09, y = .73, size = 10, hjust = 0, color="black")+ 
  draw_label("(c)",x = .09, y = .50, size = 10, hjust = 0, color="black")+
  draw_label("(d)",x = .09, y = .27, size = 10, hjust = 0, color="black")
 
ggsave("Diet_Mixing_Model.jpeg", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =6, height = 6.5, units = "in",
       dpi = 600) 
##################################################################################################






