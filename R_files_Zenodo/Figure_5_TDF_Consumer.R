#Ryan Stephens; Finalized Nov. 21, 2020
rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/TDF_Lit_Review")


#######################################################################################
#Get data ready
#######################################################################################
TDF_data<- read.csv("MetaAnalysis_TDF_data.csv",header=T)#data

library(dplyr)
TDF_data<-TDF_data %>% filter(!(Consumer_type=="Carnivore" & Diet_Source=='C4'))#Filter out C4 carnivore

TDF_data$Consumer_type<- factor(TDF_data$Consumer_type, levels = c("Herbivore","Omnivore","Carnivore"))#reorder levels

TDF_data$Diet_Source<- factor(TDF_data$Diet_Source, levels = c("C3", "C4", "Marine", "Mixed"))#reorder levels

TDF_data<-TDF_data%>%mutate(Diet_Source=ifelse(Diet_Source == "C3","C[3]",#rename levels for graphing
                                 ifelse(Diet_Source == "C4","C[4]",paste(Diet_Source))))
head(TDF_data)
str(TDF_data)
#######################################################################################




#######################################################################################
#d13C
#######################################################################################
C<-filter(TDF_data, Isotope == "d13C")
head(C)
##############################
#d13C anova
##############################
#two way anova with interaction
ANOVA_interactions_C <- aov(TDF ~ Diet_Source*Consumer_type, data=C)
library(car)
Anova(ANOVA_interactions_C, type="II")#two-way ANOVA for unbalanced group sample sizes
#kept interaction to do pairwise comparison

#Model check 
plot(ANOVA_interactions_C, 1)#homogeneity of variances
leveneTest(TDF ~ Diet_Source*Consumer_type, data=C)
plot(ANOVA_interactions_C, 2)#normality
hist(residuals(ANOVA_interactions_C),col="darkgray")

#Post hoc comparisons
library(emmeans)
lsmeans(ANOVA_interactions_C, pairwise ~  Consumer_type|Diet_Source,adjust = "Tukey")
#For mixed diets carnivores are significantly higher than herbivores or omnivores

lsmeans(ANOVA_interactions_C, pairwise ~ Diet_Source |Consumer_type,adjust = "Tukey")
#For herbivores, those feeding on a C3 diet have significantly higher TDFs than those feeding on a mixed diet

#Summary by groups
C%>%
  group_by(Diet_Source, Consumer_type)%>%
  summarise(Mean_TDF = mean(TDF), sd_TDF = sd(TDF), TDF = n())

C%>%
  group_by(Consumer_type)%>%
  summarise(Mean_TDF = mean(TDF), sd_TDF = sd(TDF), TDF = n())
##############################


##############################
#d13C plot
##############################
#Extract point of top whisker of boxplots (location to place letters)
library(ggplot2)
test<-ggplot(C, aes(x=Diet_Source, y=TDF))+#plot for data
  geom_boxplot()+
  facet_grid(~Consumer_type)
Plot_Data<-as.data.frame(layer_data(test, 1))
Ymax<-select(Plot_Data,ymax)#get ymax

Max<-C %>%#get max value (not used) along with Diet_Source,Consumer_type, TDF_13C
  group_by (Diet_Source,Consumer_type) %>%
  summarise(max = max(TDF))
Max<-as.data.frame(Max)

Max<-Max[with(Max, order(Consumer_type,Diet_Source)),]#reorder variables so letters graph correctly  
MaxY<-cbind(Max,Ymax)#bind data to ymax         

MaxY
#letters from ANOVAs
MaxY$letters<-c("","","a","","a","","","b")  

MaxY$TDF_13C<-(MaxY$ymax)#change TDF_13C so matches C data frame

#Significant bars
library(ggsignif)
Annotate_sig <- data.frame(Consumer_type=c("Herbivore", "Carnivore", "Omnivore"), 
                           start=c("C[3]", NA,NA), 
                           end=c("Mixed", NA,NA),
                           y=c(8.2,0,0),
                           label=c('italic(P)*""=="0.0001"', "",""))#"italic(P) == 0.0001"

Annotate_sig$Consumer_type<- factor(Annotate_sig$Consumer_type, levels = c("Herbivore","Omnivore","Carnivore"))#reorder levels

library(EnvStats)#sample size labels
library(ggplot2)
Plot_C<-ggplot(C, aes(x =Diet_Source , y = TDF))+
#Points
  geom_point(aes(colour=Diet_Source,shape = Consumer_type, fill=Diet_Source), size=1.5,alpha=0.3)+
  scale_shape_manual(values =c(21,22,24))+
#boxplots
  geom_boxplot(width = 0.8, alpha = 0.5, lwd=.5, aes(fill=factor(Diet_Source)),outlier.size=-2)+#boxplot fill
  geom_boxplot(width = 0.8, alpha = 0, lwd=.5,aes(color=factor(Diet_Source)),outlier.size=-2)+#boarder lines
  stat_summary(geom = "crossbar", width=0.72, fatten=2, color="black",#black median bar
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })+
  stat_summary(fun=mean, geom="point",size=2,fill="white",aes(shape = Consumer_type))+#white mean with shape indicating consumer type
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3","mediumorchid3"))+#line colors
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3","mediumorchid3"))+#fill colors
#sample size
  stat_n_text(y.pos = -2.2,size=2.5)+
#Facets and background
  facet_grid(~Consumer_type,scales = "free", space = "free")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"),#remove background lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  theme(strip.background = element_blank())+#removes gray background from titles
  theme(strip.text.x = element_text(size = 9.5))+
#letters and significance bars
  geom_text(data=MaxY,aes(label = letters, x=Diet_Source,y=(TDF_13C),vjust=-.7, color= Diet_Source),size=3)+#letters above the boxplots  
  geom_signif(data=Annotate_sig,#bars to denote significant differences for herbivores
              aes(xmin=start, xmax=end, annotations=label, y_position=y),parse = T, 
              textsize = 2.5, vjust = 1.5, color="gray65",tip_length=.02,
              manual=TRUE)+
#legend
  theme(legend.position="none")+
#Line for commonly used TDFs
  geom_hline(yintercept=1, color='grey40',linetype = "dotted")+
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-2.7, 9),breaks=c( -2, 0, 2, 4, 6, 8))+#custom scale for y-axis
  ylab(expression("TDF " * delta^13 * "C"))+#Y-axis label
  theme(axis.title.x= element_blank(),axis.text.x  = element_blank())+#remove x axis labels
  theme(axis.title.y = element_text(size=11),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=9.5, colour="black"))

Plot_C
##############################
#######################################################################################




#######################################################################################
#d15N
#######################################################################################
library(dplyr)
N<-filter(TDF_data, Isotope == "d15N")
head(N)

##############################
#d15N anova
##############################
#Two way anova with interaction
ANOVA_interactions_N<- aov(TDF ~ Diet_Source*Consumer_type, data=N)
library(car)
Anova(ANOVA_interactions_N, type="II")#two-way ANOVA for unbalanced group sample sizes

#interaction is not significant so it is dropped
ANOVA_N <- aov(TDF ~ Diet_Source+Consumer_type, data=N)
Anova(ANOVA_N, type="II")#two-way ANOVA for unbalanced group sample sizes

#Model check
plot(ANOVA_N, 1)#homogeneity of variances
leveneTest(TDF ~ Diet_Source*Consumer_type, data=N)
plot(ANOVA_N, 2)#normality
hist(residuals(ANOVA_N),col="darkgray")

library(emmeans)
#post hoc test for consumer type
lsmeans(ANOVA_N, pairwise ~  Consumer_type,adjust = "Tukey")

#summary by groups
N%>%
  group_by(Consumer_type)%>%
  summarise(Mean_TDF = mean(TDF), sd_TDF = sd(TDF), TDF = n())
##############################



##############################
#d15N plot
##############################
#Make axis labels have subscripts
N$Diet_Source<-as.factor(N$Diet_Source)
x_breaks <- levels(N$Diet_Source)
x_labels <- parse(text = x_breaks)

#Significant bars
library(ggsignif)
Annotate_sig_N <- data.frame(Consumer_type=c("Herbivore", "Omnivore", "Carnivore"), 
                           start=c("C[3]","C[3]", "C[3]"), 
                           end=c("Mixed", "Mixed", "Mixed"),
                           y=c(8, 8, 8),
                           label=c("a", "b", "b"))

Annotate_sig_N$Consumer_type<- factor(Annotate_sig_N$Consumer_type, levels = c("Herbivore","Omnivore","Carnivore"))#reorder levels



Plot_N<-
  ggplot(N, aes(x =Diet_Source , y = TDF))+
#Points
  geom_point(aes(colour=Diet_Source,shape = Consumer_type,fill=Diet_Source), size=1.5,alpha=0.3) +
  scale_shape_manual(values =c(21,22,24))+
#boxplots
  geom_boxplot(width = 0.8, alpha = 0.5, lwd=.5, aes(fill=factor(Diet_Source)),outlier.size=-2)+#boxplot fill
  geom_boxplot(width = 0.8, alpha = 0, lwd=.5,aes(color=factor(Diet_Source)),outlier.size=-2)+#boarder lines
  stat_summary(geom = "crossbar", width=0.72, fatten=2, color="black",
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })+#white median bar
  stat_summary(fun=mean, geom="point",size=2,fill="white",aes(shape = Consumer_type))+#white mean with shape indicating consumer type
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3","mediumorchid3"))+#line colors
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3","mediumorchid3"))+#fill colors
#sample size
  stat_n_text(y.pos = -2.2,size=2.5) +
#Facets and background
  facet_grid(~Consumer_type,scales = "free",space = "free")+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),#remove background lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+#remove facets
#letters and significance bars
  geom_signif(data=Annotate_sig_N,#bars to denote significant differences for herbivores
              aes(xmin=start, xmax=end, annotations=label, y_position=y),parse = TRUE, 
              textsize = 2.5, vjust = 0, color="gray65",tip_length=0,
              manual=TRUE)+
#legend
  theme(legend.position="none")+
#Line for commonly used TDFs
  geom_hline(yintercept=3, color='grey50',linetype = "dotted")+
#Axes
  scale_y_continuous(expand = c(0, 0),limits = c(-2.7, 9.5),breaks=c(-2, 0, 2, 4, 6, 8))+#custom scale
  scale_x_discrete(breaks = x_breaks, label = x_labels)+#Make x-axis labels subscripts
  xlab("Diet source")+#x-axis label
  ylab(expression("TDF " * delta^15 * "N"))+#y-axis label
  theme(axis.title.x= element_text(size=11),
        axis.text.x  = element_text(angle=40, hjust = 1,size=9.5, colour="black"))+
  theme(axis.title.y = element_text(size=11),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=9.5, colour="black"))

Plot_N
##############################
#######################################################################################





#######################################################################################
#Final plot
#######################################################################################
#Silhouettes for consumers
library(rphylopic)#pull silhouette from phylopic
#used
cow <- image_data("415714b4-859c-4d1c-9ce0-9e1081613df7", size = "128")[[1]]#http://phylopic.org/image/415714b4-859c-4d1c-9ce0-9e1081613df7/
mouse <- image_data("3b89954e-b28e-4062-89cd-1ad5df0eb431", size = "128")[[1]]#http://phylopic.org/image/3b89954e-b28e-4062-89cd-1ad5df0eb431/
tiger <- image_data("e148eabb-f138-43c6-b1e4-5cda2180485a", size = "128")[[1]]#http://phylopic.org/image/e148eabb-f138-43c6-b1e4-5cda2180485a/

#not used
#deer <- image_data("68ce0b27-c1af-40df-bf89-12e048ea893e", size = "128")[[1]]
#cow2 <- image_data("aae1d17f-149f-41fd-b14b-72f1c95e3f28", size = "128")[[1]]#http://phylopic.org/image/aae1d17f-149f-41fd-b14b-72f1c95e3f28/
#hare <- image_data("8e61e166-11f4-4377-a923-9b5b597b6eba", size = "128")[[1]]#http://phylopic.org/image/8e61e166-11f4-4377-a923-9b5b597b6eba/
#caribou <- image_data("e6e864fd-8e3d-435f-9db3-dc6869c589f1", size = "128")[[1]]#http://phylopic.org/image/e6e864fd-8e3d-435f-9db3-dc6869c589f1/
#mouse2 <- image_data("570c7d9e-e6d1-46f5-b165-988981bfc5f6", size = "128")[[1]]#http://phylopic.org/image/570c7d9e-e6d1-46f5-b165-988981bfc5f6/
#wolf <- image_data("8cad2b22-30d3-4cbd-86a3-a6d2d004b201", size = "128")[[1]]#http://phylopic.org/image/8cad2b22-30d3-4cbd-86a3-a6d2d004b201/


library(cowplot)


#graph using cowplot
ggdraw() +
#plots
  draw_plot(Plot_N, x = 0, y = 0, width = 1, height = .53)+#d15N plot
  draw_plot(Plot_C, x = 0, y = .51, width = 1, height = .4625)+#d13C plot
#text for anova results
  #draw_label(expression("Diet source: F"["("][3][","][120][")"] *" = 7.19; " *italic(P) == "0.0002"),#Diet d13C
  #           x = 1, y = .88, size = 7, hjust = 1, color="gray65")+
  #draw_label(expression("Consumer: F"["("][2][","][120][")"] *" = 6.36; " *italic(P) == 0.002),#Consumer d13C
  #           x = 1, y = .86, size = 7, hjust = 1, color="gray65")+
  #
  #draw_label(expression("Diet source: F"["("][3][","][99][")"] *" = 0.37; " *italic(P) == "0.330"),#Diet d15N
  #           x = 1, y = .44, size = 7, hjust = 1, color="gray65")+
  #draw_label(expression("Consumer: F"["("][2][","][99][")"] *" = 3.70; " *italic(P) == 0.005),#Consumer d15N
  #           x = 1, y = .42, size = 7, hjust = 1, color="gray65")+
#letters
  draw_label("(a)",x = 0, y = .90, size = 11, hjust = 0, color="black")+#Diet d13C
  draw_label("(b)",x = 0, y = .50, size = 11, hjust = 0, color="black")+#Diet d13C          
#silhouettes 
  add_phylopic(cow, alpha=0.50, x=.32, y=.97, ysize=.08, color="black")+
  add_phylopic(mouse, alpha=0.50, x=.58, y=.97, ysize=.08, color="black")+
  add_phylopic(tiger, alpha=0.50, x=.84, y=.97, ysize=.08, color="black")

ggsave("Figure_5_TDF_boxplot.tiff", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =3, height = 5.5, units = "in",
       dpi = 600) 
#######################################################################################
 



#######################################################################################
#Table of average TDFs for Supplemental Information (Table S2)
#######################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/TDF_Lit_Review")
TDF_data<- read.csv("MetaAnalysis_TDF_data.csv",header=T)#data
TDF_data<-TDF_data %>% filter(!(Consumer_type=="Carnivore" & Diet_Source=='C4'))#Filter out C4 carnivore

##############################
#TDF averages (unique combinations of Diet and consumer for both isotopes)
##############################
TDF_Diet_consumer<-TDF_data%>%
  group_by(Isotope, Consumer_type, Diet_Source)%>% #Grouping variables
  summarise(TDF_Mean=round(mean(TDF),2),TDF_SD=round(sd(TDF),2), n = n())
head(TDF_Diet_consumer)
##############################


##############################
#d15N-TDF averages (consumer groups)
##############################
TDF_d15N<-filter(TDF_data, Isotope == "d15N")#filter nitrogen

TDF_d15N_consumer<-TDF_d15N%>%
  group_by(Isotope, Consumer_type)%>% #Grouping variables
  summarise(TDF_Mean=round(mean(TDF),2),TDF_SD=round(sd(TDF),2), n = n())

TDF_d15N_consumer<-mutate(TDF_d15N_consumer, Diet_Source = ifelse(Consumer_type == "Herbivore", "C3 + C4 + Mixed",#herbivores are from C3, C4, and mixed systems
                                                         ifelse(Consumer_type == "Omnivore", "C3 + Mixed",#Omnivores are from C3, C4, and mixed systems
                                                                "C3 + Marine + Mixed")))#carnivores are from C3, marine, and mixed systems

TDF_d15N_consumer<-select(TDF_d15N_consumer, Isotope, Consumer_type, Diet_Source, TDF_Mean, TDF_SD, n)#reorder
head(TDF_d15N_consumer)
##############################


##############################
#d15N-TDF averages (terrestrial carnivore)
TDF_d15N_Ter_Carn<-filter(TDF_d15N, Consumer_type == "Carnivore" & Diet_Source == "C3" |#filter terrestrial carnivore
                                    Consumer_type == "Carnivore" & Diet_Source == "Mixed")

TDF_d15N_Ter_Carn_mean<-TDF_d15N_Ter_Carn%>%
  group_by(Isotope, Consumer_type)%>% #Grouping variables
  summarise(TDF_Mean=round(mean(TDF),2),TDF_SD=round(sd(TDF),2), n = n())

TDF_d15N_Ter_Carn_mean<-mutate(TDF_d15N_Ter_Carn_mean, Diet_Source =  "C3 + Mixed")

TDF_d15N_Ter_Carn_mean<-select(TDF_d15N_Ter_Carn_mean, Isotope, Consumer_type, Diet_Source, TDF_Mean, TDF_SD, n)#reorder
head(TDF_d15N_Ter_Carn_mean)
##############################


##############################
#Merge TDF data
##############################
TDF_d13C_d15N<-rbind(TDF_Diet_consumer, TDF_d15N_consumer, TDF_d15N_Ter_Carn_mean)
TDF_d13C_d15N<-filter(TDF_d13C_d15N, n > 2)#select TDFs derived from at least three samples
##############################


##############################
#Adjust hair to other tissue types based on values from Table S2
##############################
TDF_d13C_d15N<-rename(TDF_d13C_d15N, Hair = TDF_Mean)#rename TDF_Mean to Hair
TDF_d13C_d15N<- mutate(TDF_d13C_d15N, Collagen = ifelse(Isotope == "d13C", Hair+1.21, Hair-0.22))#Collagen
TDF_d13C_d15N<- mutate(TDF_d13C_d15N, Muscle = ifelse(Isotope == "d13C", Hair-1.31, Hair-0.15))#Muscle
TDF_d13C_d15N<- mutate(TDF_d13C_d15N, Liver = ifelse(Isotope == "d13C", Hair-1.65, Hair+0.32))#Liver
TDF_d13C_d15N<- mutate(TDF_d13C_d15N, Blood = ifelse(Isotope == "d13C", Hair-1.50, Muscle-0.26))#Blood (Adjusted from muscle since hair and muscle were not correlated)

TDF_d13C_d15N<-select(TDF_d13C_d15N, Isotope, Consumer = Consumer_type, Diet_Source, 
                                     Hair, Collagen, Muscle, Liver, Blood, SD = TDF_SD, n)
TDF_d13C_d15N<-TDF_d13C_d15N %>%#arrange table
  arrange(Isotope, match(Consumer, c("Herbivore","Omnivore","Carnivore")))
                        



write.csv(TDF_d13C_d15N, "Table_TDF_Averages.csv", row.names = F)
#######################################################################################

