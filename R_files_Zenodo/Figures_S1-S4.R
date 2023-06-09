###--------------------------Figure S1----------------------------###
# loading packages
library(ggplot2)
library(ggthemes)
library(gridExtra)

# reading data of publication_year in csv file
pub_stemflow <- read.table("publication_year.csv", header = TRUE, sep = ",")

# 1. plotting publication ~ year
pub_number <- ggplot(data = pub_stemflow)+
  geom_bar(aes(x = as.factor(Year), y = Publications, fill = Year), stat = "identity")+ 
  scale_y_continuous(limits = c(0,80))+
  geom_text(aes(x = as.factor(Year), y = Publications + 2, label = Publications), size = 3, fontface = "bold")+
  labs(x = "Year", y = "Number of Publications", fontface = "bold")+
  theme_light()+
  theme(axis.text.x = element_text(size = 12, angle = -90, vjust = 0.5), axis.text.y = element_text(size = 12), axis.title= element_text(size=15), legend.position = "none") + 
  theme(panel.grid.major = element_line(linetype = "dashed"), panel.grid.minor = element_line(linetype = "blank"))

# 2. plotting selected publication ~ year
pub_selected <- ggplot(data = pub_stemflow)+
  geom_bar(aes(x = as.factor(Year), y = Selected, fill = Year), stat = "identity")+ 
  scale_y_continuous(limits = c(0,80))+
  geom_text(aes(x = as.factor(Year), y = Selected + 2, label = Selected), size = 3, fontface = "bold")+
  labs(x = "Year", y = "Number of Publications", fontface = "bold")+
  theme_light()+
  theme(axis.text.x = element_text(size = 12, angle = -90, vjust = 0.5), axis.text.y = element_text(size = 12), axis.title= element_text(size=15), legend.position = "none") + 
  theme(panel.grid.major = element_line(linetype = "dashed"), panel.grid.minor = element_line(linetype = "blank"))

# Printing Figure S1
Figure_S1 <- grid.arrange(pub_number, pub_selected, ncol=1, nrow =2)

# ---subfigure in Figure S1a---
library(tmap)
data("World")
library(sp)

df1 <- read.csv("publications_world.csv", header = TRUE, sep=",")
World <- sp::merge(World, df1, by = "name")

tm_shape(World, projection = 4326) +
  tm_polygons("Publications",
              breaks = c(0, 1, 5, 10, 20, 50, 100, 200, 300), 
              palette="-Set2") +
  tm_style("white")+
  tm_format("World")+
  tm_legend(position=c("left", "bottom"), scale = 1.2, bg.color="grey95", bg.alpha = 0.2, frame=FALSE)

###--------------------------Figure S2----------------------------###
library(PerformanceAnalytics)

#reading data of correlation_matrix in csv file
sf_variables <- read.table("correlation_matrix.csv", header = TRUE, sep = ",") 
names(sf_variables)
# testing normality of input continuous variables (MAP, MAT, Height, Density, Age, Height, DBH, Basal area)
# 1. MAP
shapiro.test(sf_variables$MAP)
# 2. MAT
shapiro.test(sf_variables$MAT)
# 3. Height
shapiro.test(sf_variables$Height)
# 4. Density
shapiro.test(sf_variables$Density)
# 5. Age
shapiro.test(sf_variables$Age)
# 6. LAI
shapiro.test(sf_variables$LAI)
# 7. DBH
shapiro.test(sf_variables$DBH)
# 8. Basal area
shapiro.test(sf_variables$Basal_area)

# According to above normality test, the data were not normally distributed, therefore we choose spearman test 
# Printing Figure S2
chart.Correlation(sf_variables, histogram = TRUE, pch = 20, method = c("spearman"))

###--------------------------Figure S3----------------------------###
library(tidyverse)
library(ggbeeswarm) 
library(ggthemes)
library(ggpubr) 

# reading source_data
GSF<-read.table("Source_data.csv", header = T, sep=",") 
GSF_Bark <- GSF %>% drop_na(Stemflow, Bark)

SF_Bark <- ggplot(data = GSF_Bark, aes(x = Bark, y = Stemflow, fill = Bark)) +
  geom_violin(alpha = 0.3) +
  geom_quasirandom(alpha = 0.8, size = 3, 
                   varwidth = TRUE,
                   aes(color = Bark)) +
  geom_boxplot(width = .15,
               notch = TRUE,
               fill = "orange",
               outlier.shape = NA,
               alpha = 0.5)+
  scale_y_continuous(limits = c(0, 52))+
  stat_compare_means(comparisons = list(c("M", "R"), c("M", "S"), c("R", "S")), step.increase = 0.05, hide.ns = FALSE, vjust = 0.5, size = 6, tip.length = 0.02, bracket.size = 0.5, label = "p.signif", label.y = c(50, 52, 51))+
  labs(x = "Bark", y = "Stemflow percentage (%)") +
  theme_base() +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.position = "none")

# Printing Figure S3
print(SF_Bark)

###--------------------------Figure S4----------------------------###
# loading packages
library(ggthemes)
library(ggplot2)

# reading Source_data in csv file
GSF<-read.table("Source_data.csv", header = T, sep = ",") 

Figure_S4 <- ggplot(GSF, aes(MAP, fill = Community_type)) + 
  geom_density(alpha = 0.6, outline.type = "lower")+
  labs(x = "MAP", y = "Probability density")+
  theme_base()+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.9), legend.text = element_text(size = 20))+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30))

# Printing Figure S4
print(Figure_S4)
