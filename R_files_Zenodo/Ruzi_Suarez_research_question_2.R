## Describing the chamber depths of Ectatomma ruidum nests using wax casts.
## Comparing chamber depth to estimated maximum emergence depth, and
## conducting an LMM to assess whether seeds could survive various burial depths.

## RQ2: Can seeds survive and emerge from the deposition location?

# abbreviations: Zan = Zanthoxylum ekmanii; ER = Ectatomma ruidum

# -- load libraries ####
library(here)
library(tidyverse) # includes ggplot2
library(plyr)
library(nlme)
library(emmeans)

# -- set paths ####
data_path <- here::here("data")
figure_path <- here::here("figs")

# -- read in the data ####

depth_raw <- read_csv(paste(data_path, "Ruzi_Suarez_data_chamber_depths.csv", sep = "/"),
         col_types = 
           cols(
             .default = col_double(),
             Chamber = col_character(),
             Colony = col_character()))
nrow(depth_raw)

seedling_data <- read_csv(paste(data_path, "Ruzi_Suarez_data_seedling_emergence.csv", sep = "/"),
         col_types = 
           cols(
             Depth_cm = col_double(),
             Replicate = col_double(),
             Num_seeds_buried = col_double(),
             Num_seeds_emerged = col_double(),
             Prop_emerged = col_double(),
             Mass_soil = col_double(),
             Soil_autoclaved = col_character(),
             Holes = col_character()))

nrow(seedling_data)

## - to determine chamber depths and if different from the estiamted maximum emergence depth ####
# resuls summarized in Table 1 of the main text
# maximum emergence depth was calculated using the formula from Bond, Honig, and Maze 1999

depth_summary <- ddply(depth_raw, "Chamber", summarise,
                       n=sum(!is.na(Depth)),
                       Mean=round(mean(Depth, na.rm=TRUE),2),
                       Median=round(median(Depth, na.rm=TRUE),2),
                       SD=round(sd(Depth, na.rm=TRUE),2),
                       IQR=round(IQR(Depth, na.rm=TRUE),2),
                       min=round(min(Depth, na.rm=TRUE),2),
                       max=round(max(Depth, na.rm=TRUE),2))
depth_summary
#   Chamber  n  Mean Median    SD  IQR  min  max
#1 chamber1 13  7.72   8.50  2.29  2.5  3.8 13.2
#2 chamber2 13  9.92   9.80  2.39  1.7  5.7 15.7
#3 chamber3 10 13.58  12.05  3.94  5.6  6.7 19.2
#4 chamber4  5 15.88  15.30  6.63  8.5  6.5 22.7
#5 chamber5  3 26.37  22.70 15.33 15.0 13.2 43.2
#6 chamber6  1 25.20  25.20    NA  0.0 25.2 25.2

### -- chamber 1 data
depth_chamber1 <- depth_raw %>%
  filter(Chamber == "chamber1")
depth_cham1_v <- depth_chamber1$Depth
shapiro.test(depth_cham1_v)
#Shapiro-Wilk normality test
#data:  depth_cham1_v
#W = 0.90184, p-value = 0.1417

t.test(depth_cham1_v, mu=7.026, alternative="two.sided") 
#One Sample t-test
#data:  depth_cham1_v
#t = 1.0981, df = 12, p-value = 0.2937
#alternative hypothesis: true mean is not equal to 7.026
#95 percent confidence interval:
#  6.339971 9.106183
#sample estimates:
#  mean of x 
#7.723077 

### -- chamber 2 data
depth_chamber2 <- depth_raw %>%
  filter(Chamber == "chamber2")

depth_cham2_v <- depth_chamber2$Depth
shapiro.test(depth_cham2_v)
#Shapiro-Wilk normality test
#data:  depth_cham2_v
#W = 0.92551, p-value = 0.2974


t.test(depth_cham2_v, mu=7.026, alternative="two.sided") 
#One Sample t-test
#data:  depth_cham2_v
#t = 4.3789, df = 12, p-value = 0.000898
#alternative hypothesis: true mean is not equal to 7.026
#95 percent confidence interval:
#  8.481579 11.364575
#sample estimates:
#  mean of x 
#9.923077 



### -- chamber 3 data
depth_chamber3 <- depth_raw %>%
  filter(Chamber == "chamber3")

depth_cham3_v <- depth_chamber3$Depth
shapiro.test(depth_cham3_v)
#Shapiro-Wilk normality test
#data:  depth_cham3_v
#W = 0.91152, p-value = 0.2916

t.test(depth_cham3_v, mu=7.026, alternative="two.sided") 
#One Sample t-test
#data:  depth_cham3_v
#t = 5.2546, df = 9, p-value = 0.0005245
#alternative hypothesis: true mean is not equal to 7.026
#95 percent confidence interval:
#  10.75843 16.40157
#sample estimates:
#  mean of x 
#13.58 


### -- chamber 4 data
depth_chamber4 <- depth_raw %>%
  filter(Chamber == "chamber4")

depth_cham4_v <- depth_chamber4$Depth
shapiro.test(depth_cham4_v)

#Shapiro-Wilk normality test
#data:  depth_cham4_v
#W = 0.93423, p-value = 0.6255

t.test(depth_cham4_v, mu=7.026, alternative="two.sided") 

#One Sample t-test
#data:  depth_cham4_v
#t = 2.9856, df = 4, p-value = 0.04051
#alternative hypothesis: true mean is not equal to 7.026
#95 percent confidence interval:
#  7.646354 24.113646
#sample estimates:
#  mean of x 
#15.88 


### -- chamber 5 data
depth_chamber5 <- depth_raw %>%
  filter(Chamber == "chamber5")

depth_cham5_v <- depth_chamber5$Depth
shapiro.test(depth_cham5_v)
#Shapiro-Wilk normality test
#data:  depth_cham5_v
#W = 0.95711, p-value = 0.6016

t.test(depth_cham5_v, mu=7.026, alternative="two.sided") 
#One Sample t-test
#data:  depth_cham5_v
#t = 2.1848, df = 2, p-value = 0.1605
#alternative hypothesis: true mean is not equal to 7.026
#95 percent confidence interval:
#  -11.72119  64.45453
#sample estimates:
#  mean of x 
#26.36667 


### -- chamber 6 data
depth_chamber6 <- depth_raw %>%
  filter(Chamber == "chamber6")

depth_cham6_v <- depth_chamber6$Depth
shapiro.test(depth_cham6_v)
#Error in shapiro.test(depth_cham6_v) : 
#  sample size must be between 3 and 5000

t.test(depth_cham6_v, mu=7.026, alternative="two.sided") 
# not enough observations to run which is expected since only one colony has this chamber

## - to describe where larvae, workers, and Zan seeds are found within the colony ####


larvae_num_summary<-ddply(depth_raw, "Chamber", summarise,
                              n=sum(!is.na(Larvae.num)),
                              Mean=round(mean(Larvae.num, na.rm=TRUE),2),
                              Median=round(median(Larvae.num, na.rm=TRUE),2),
                              SD=round(sd(Larvae.num, na.rm=TRUE),2),
                              IQR=round(IQR(Larvae.num, na.rm=TRUE),2))
larvae_num_summary
#   Chamber  n  Mean Median    SD   IQR
#1 chamber1 13 30.46     27 25.08 33.00
#2 chamber2 13 13.92      6 20.56  8.00
#3 chamber3 10  3.70      2  4.88  5.25
#4 chamber4  6  1.50      0  3.21  0.75
#5 chamber5  3  0.33      0  0.58  0.50
#6 chamber6  1  0.00      0    NA  0.00

larvae_percent_summary<-ddply(depth_raw, "Chamber", summarise,
                              n=sum(!is.na(Larvae.percent)),
                              Mean=round(mean(Larvae.percent, na.rm=TRUE),2),
                              Median=round(median(Larvae.percent, na.rm=TRUE),2),
                              SD=round(sd(Larvae.percent, na.rm=TRUE),2),
                              IQR=round(IQR(Larvae.percent, na.rm=TRUE),2))
larvae_percent_summary
#   Chamber  n  Mean Median    SD   IQR
#1 chamber1 13 63.00  72.97 35.41 46.87
#2 chamber2 13 24.92  16.22 27.17 28.16
#3 chamber3 10  8.09   6.08  9.36 12.50
#4 chamber4  6 11.31   0.00 27.12  0.89
#5 chamber5  3  2.78   0.00  4.81  4.17
#6 chamber6  1  0.00   0.00    NA  0.00

worker_percent_summary<-ddply(depth_raw, "Chamber", summarise,
                              n=sum(!is.na(Worker.percent)),
                              Mean=round(mean(Worker.percent, na.rm=TRUE),2),
                              Median=round(median(Worker.percent, na.rm=TRUE),2),
                              SD=round(sd(Worker.percent, na.rm=TRUE),2),
                              IQR=round(IQR(Worker.percent, na.rm=TRUE),2))
worker_percent_summary
#  Chamber  n  Mean Median    SD   IQR
#1 chamber1 13 28.87  33.33 17.20 22.57
#2 chamber2 13 35.48  32.93 23.19 26.01
#3 chamber3 10 26.60  25.30 12.79 19.22
#4 chamber4  6  9.68   6.85  7.19  5.12
#5 chamber5  3 28.69  18.82 22.36 20.66
#6 chamber6  1  3.70   3.70    NA  0.00


ZAEK_percent_summary<- ddply(depth_raw, "Chamber", summarise,
                             n=sum(!is.na(ZAEKseed.percent)),
                             Mean=round(mean(ZAEKseed.percent, na.rm=TRUE),2),
                             Median=round(median(ZAEKseed.percent, na.rm=TRUE),2),
                             SD=round(sd(ZAEKseed.percent, na.rm=TRUE),2),
                             IQR=round(IQR(ZAEKseed.percent, na.rm=TRUE),2))
ZAEK_percent_summary
#   Chamber n  Mean Median    SD  IQR
#1 chamber1 7 14.69      0 19.69 30.0
#2 chamber2 7 48.57     50 41.80 65.0
#3 chamber3 4 37.50     25 47.87 62.5
#4 chamber4 1  0.00      0    NA  0.0
#5 chamber5 1  0.00      0    NA  0.0
#6 chamber6 0   NaN     NA    NA   NA


## - to graph chamber depths and where Zan seeds are found ####

depth_raw2 <- depth_raw %>%
  # create a new column so can have an example of depths for one colony
  mutate(Colonyalt = if_else(Colony != "SAR843", "two", "one"))

depth <- ggplot(depth_raw2, aes(x=Chamber, y=Depth, color=Colonyalt))+ #shape=Colony
  #but if i change color to fill then the stat summary goes wonky and i get the 
  #black average bars on every dot
  geom_point(size=3, position=position_jitter(width=.3, height=0))+ #shape = 21
  stat_summary(fun="mean", geom="point", shape=95, size=10, colour="black")+
  ylab("Chamber depth (cm)\\n")+
  xlab("Chamber Number")+
  scale_colour_manual(values=c("grey40", "red"), breaks=c("two","one"), labels=c(" ", "SAR843"), guide="none")+
  scale_y_reverse(limits=c(45,-4))+
  geom_hline(yintercept=7.026, colour="grey20")+
  geom_hline(yintercept=0, colour="brown")+
  annotate("text", x=5.5, size=3, y=-1.5, colour="brown", label="Ground Surface")+
  annotate("text", x=5.5, size=3, y=4, colour="grey20", label="Calculated Maximum\\nEmergence Depth" )+
  #annotate("text", x=.1,family="Times",size=5, y=(600/ratiosample), label=paste("mu(mean) is ",rhizobium_sample_mean))
  scale_x_discrete(limit=c("chamber1", "chamber2", "chamber3", "chamber4", "chamber5", "chamber6"), 
                   labels=c("1", "2", "3", "4", "5", "6"))+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_blank(), #gets rid of square going around the entire graph
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        legend.position=c(0.10,0.25), # moves the location of the legend
        legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(face="bold", size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"),
        axis.line.x	=element_line(colour="black"),
        axis.line.y	=element_line(colour="black"))#size of y-axis text
depth

Zan_seeds <- ggplot(depth_raw, aes(x=Chamber, y=ZAEKseed.percent)) +
  geom_boxplot()+
  ylab("Relative Abundance\\nof Seeds (%)")+
  #ggtitle("Larvae")+
  xlab(" ")+
  #scale_y_reverse(limits=c(45,0))+
  #geom_hline(yintercept=7, colour="grey54")+
  #geom_hline(yintercept=0, colour="brown")+
  #annotate("text", x=5, family="Times", size=4, y=1, colour="brown", label="Ground Surface")+
  #annotate("text", x=5, family="Times", size=3, y=8, colour="grey54", label="Calculated Maximum Emergence Depth" )+
  #annotate("text", x=.1,family="Times",size=5, y=(600/ratiosample), label=paste("mu(mean) is ",rhizobium_sample_mean))
  scale_x_discrete(limit=c("chamber1", "chamber2", "chamber3", "chamber4", "chamber5", "chamber6"), 
                   labels=c(" "," "," "," "," "," "))+
  #labels=c("Chamber 1", "Chamber 2", "Chamber 3", "Chamber 4", "Chamber 5", "Chamber 6"))+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_blank(), #gets rid of square going around the entire graph
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(face="bold", size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"))#size of y-axis text
Zan_seeds

## - seedling emergence analyses ####
seedling_data

# to check whether the density of the soil compares to what has been published in the literature
# for the region (from Cavelier 1992)

rad <- 1.25 # the radius of the clear tube used to plant the seeds in
rad_squared <- rad*rad
seedling_data <- seedling_data %>%
  # update the height as all tubes had 5 cm of soil below were the seeds were placed
  mutate(height = 5 + Depth_cm) %>%
  mutate(volume = pi*rad_squared*height) %>%
  mutate(density = Mass_soil/volume) %>%
  # get a percentage of seedlings that emerged
  mutate(percent = Prop_emerged*100) %>%
  mutate(Depth2 = as_factor(Depth_cm))

# formula for standard error
se <- function(x) sqrt(var(x)/length(x))

density_summary <- ddply(seedling_data, "Depth2", summarise,
                         n=sum(!is.na(density)),
                         Mean=mean(density,na.rm=TRUE),
                         SD=sd(density,na.rm=TRUE),
                         se=se(density))                        
density_summary
#  Depth2  n      Mean         SD          se
#1        0  4 0.8043176 0.05825266 0.029126331
#2        2 14 0.8267212 0.03395842 0.009075769
#3        5 12 0.8255668 0.03889498 0.011228013
#4        7  8 0.8204096 0.04718091 0.016680971

mean(seedling_data$density) # 0.8226696 is the mean soil density used
sd(seedling_data$density) # 0.03997068
se(seedling_data$density) # 0.0064841

#to compare to the published density of 0.74

shapiro.test(seedling_data$density)
#Shapiro-Wilk normality test
#data:  seedling_data$density
#W = 0.92313, p-value = 0.01226

wilcox.test(seedling_data$density, mu=0.74, alternative="two.sided")
# Wilcoxon signed rank test
#data:  seedling_data2$density
#V = 731, p-value = 3.129e-10
#alternative hypothesis: true location is not equal to 0.74

3.129e-10 < 0.001 # TRUE


# for seedling emergence running a LMM

Model3 <- lme(Prop_emerged~Depth2+Soil_autoclaved, random=~1|Replicate,
              method="ML", data=seedling_data)
summary(Model3)

plot(fitted(Model3), residuals(Model3), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Model3), residuals(Model3)))
h<-qqnorm(residuals((Model3)))
qqline(residuals(Model3))
anova(Model3)
#                 numDF denDF   F-value p-value
#(Intercept)         1    20 15.523701  0.0008
#Depth2              3    20  7.682207  0.0013
#Soil_autoclaved     1    20  0.980156  0.3340

emmeans(Model3, pairwise~Depth2, adjust="tukey")
#$emmeans
#Depth2 emmean     SE df lower.CL upper.CL
#0      0.4116 0.1067 13    0.181    0.642
#2      0.3077 0.0629 13    0.172    0.444
#5      0.0305 0.0682 13   -0.117    0.178
#7      0.0359 0.0791 13   -0.135    0.207
#
#Results are averaged over the levels of: Soil_autoclaved 
#Degrees-of-freedom method: containment 
#Confidence level used: 0.95 
#
#$contrasts
#contrast estimate     SE df t.ratio p.value
#0 - 2     0.10389 0.1139 20  0.912  0.7988 
#0 - 5     0.38103 0.1160 20  3.283  0.0179 
#0 - 7     0.37564 0.1191 20  3.153  0.0238 
#2 - 5     0.27715 0.0740 20  3.746  0.0064 
#2 - 7     0.27176 0.0852 20  3.189  0.0220 
#5 - 7    -0.00539 0.0870 20 -0.062  0.9999 
#
#Results are averaged over the levels of: Soil_autoclaved 
#Degrees-of-freedom method: containment 
#P value adjustment: tukey method for comparing a family of 4 estimates 

## to plot the figure for the supplement

seedling_boxplot <- ggplot(seedling_data, aes(x=Depth2, y=percent))+
  geom_boxplot(fill="grey40")+
  ylab("Seedlings Emerged (%)")+
  scale_y_continuous(limits=c(0,110), expand=c(0,0))+
  xlab("Depth")+
  annotate("text", x="0", size=4, 
           y=105, colour="black", label="A")+
  annotate("text", x="2", size=4, 
           y=105, colour="black", label="A")+
  annotate("text", x="5", size=4, 
           y=105, colour="black", label="B")+
  annotate("text", x="7", size=4, 
           y=105, colour="black", label="B")+
  theme(axis.line=element_line(colour="black"),
        strip.text.y=element_text(size=16, face="italic", angle=0),
        strip.background=element_blank(),
        strip.text.x=element_text(face="bold", colour="black", angle=45, hjust=1, vjust=1, size=16),
        axis.text.y=element_text(colour="black", size=12),
        axis.title.x=element_text(face="bold", colour="black", size=14,margin = margin(20,0,0,0)),
        axis.title.y=element_text(face="bold", colour="black", size=14, margin=margin(0,20,0,0)),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.title=element_text(colour="black", face="bold", size=12),
        legend.text=element_text(colour="black", size=12),
        #plot.title=element_text(family="Arial", face="bold", colour="black", size=rel(1.0)),
        panel.border=element_blank(),
        #legend.key=element_rect(colour="black", size=1),
        legend.position="top",
        panel.background=element_blank(),
        axis.line.y	=element_line(colour="black"),
        axis.line.x=element_line(colour="black"))
seedling_boxplot

## - to save all the figures ####
ggsave("Ruzi_Suarez_Figure2b.png", width = 6, height = 4.5,
       units = "in", dpi = 600, plot = depth,
       path = figure_path)


ggsave("Ruzi_Suarez_Figure2a.png", width = 6, height = 3,
       units = "in", dpi = 600, plot = Zan_seeds,
       path = figure_path)


ggsave("Ruzi_Suarez_Figure2c.png", width = 6, height = 4,
       units = "in", dpi = 600, plot = seedling_boxplot,
       path = figure_path)
