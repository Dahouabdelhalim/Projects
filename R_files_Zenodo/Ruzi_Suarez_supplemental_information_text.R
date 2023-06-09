## Assessing seed removal by year and whether baiting influenced seed removal amounts

# This uses data collected under sampling schemes from Ruzi et al. 2017. Plant Ecology. doi:10.1007/s11258-017-0745-7
# and Ruzi et al. 2021. Biotropica. doi:10.1111/btp.12904
# This time however, seed removal is analyzed between 1000-1200 hours

# This analysis appears in the Supplemental Information text

# abbreviations: Zan = Zanthoxylum ekmanii; ER = Ectatomma ruidum

# -- load libraries ####
library(here)
library(tidyverse) # includes ggplot2
library(nlme)
library(emmeans)
library(plyr)

# -- set paths ####
data_path <- here::here("data")
figure_path <- here::here("figs")

# -- read in the data ####

removal_data <- read_csv(paste(data_path, "Ruzi_Suarez_data_seed_removal_2hours.csv", sep = "/"),
         col_types = 
           cols(
             Plot = col_character(),
             Year = col_character(),
             Rep = col_double(),
             Side = col_character(),
             Seeds_placed = col_double(),
             Seeds_removed = col_double(),
             Seeds_remaining = col_double(),
             Date = col_character()))
removal_data

names(removal_data)


# for seed removal by year running a LMM ####

# need to first set up the response variable into a proportion
removal_data2 <- removal_data %>%
  mutate(Prop_removed = Seeds_removed/Seeds_placed)

removal_data2

Model1 <- lme(Prop_removed~Year, random = ~1|Plot/Rep,
              method="ML", data=removal_data2)

summary(Model1)

plot(fitted(Model1), residuals(Model1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Model1), residuals(Model1)))
h<-qqnorm(residuals((Model1)))
qqline(residuals(Model1))
anova(Model1)
#            numDF denDF   F-value p-value
#(Intercept)     1    18 20.240006  0.0003
#Year            2    18  2.114003  0.1497

# summary of seed removal by year ####

removal_summary <- ddply(removal_data2, "Year", summarise,
                       n=sum(!is.na(Prop_removed)),
                       Mean=round(mean(Prop_removed, na.rm=TRUE),2),
                       Median=round(median(Prop_removed, na.rm=TRUE),2),
                       SD=round(sd(Prop_removed, na.rm=TRUE),2),
                       IQR=round(IQR(Prop_removed, na.rm=TRUE),2),
                       min=round(min(Prop_removed, na.rm=TRUE),2),
                       max=round(max(Prop_removed, na.rm=TRUE),2))
removal_summary
#   Year  n Mean Median   SD  IQR min max
#1 2013 10 0.31    0.2 0.34 0.28   0   1
#2 2015 15 0.59    0.8 0.44 0.80   0   1
#3 2016 10 0.58    0.7 0.44 0.80   0   1

# making the figure for removal by year and saving it to include in the supplement ####


removal_by_year_boxplot <- ggplot(removal_data2, aes(x=Year, y=Prop_removed))+
  geom_boxplot(fill="grey40")+
  ylab("Proportion of\\nSeeds Removed")+
  scale_y_continuous(limits=c(0,1.2), expand=c(0,0))+
  xlab("Year")+
  theme(axis.line=element_line(colour="black"),
        strip.text.y=element_text(size=16, face="italic", angle=0),
        strip.background=element_blank(),
        strip.text.x=element_text(face="bold", colour="black", angle=45, hjust=1, vjust=1, size=16),
        axis.text.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=12),
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
removal_by_year_boxplot

# to save the figure ####


ggsave("Ruzi_Suarez_FigureS1.png", width = 4, height = 3,
       units = "in", dpi = 300, plot = removal_by_year_boxplot,
       path = figure_path)
