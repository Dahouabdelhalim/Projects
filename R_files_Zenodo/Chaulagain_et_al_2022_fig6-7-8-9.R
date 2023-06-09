###############################################################
###### Disease Gradient plots #######
############ BY Bhim Chaulagain #################
############### Oregon State University ######################

# install required packages
install.packages(c("readxl", "reshape2", "cowplot", "scales", "grid", "tidyverse", "RColorBrewer", "viridis"))

# Load required packages
library(readxl)
library(RColorBrewer)
library(grid)
library(scales)
library(cowplot)
library(reshape2)
library(tidyverse)
library(viridis)

# read in data
data<-read_excel("Chaulagain_et_al_2022_Data_fig_6_7_8_9", sheet = "data")
head(data)

#reshape data
test_data_long <- melt(data, id=c("Distance", "Size", 'Treatment', 'Time'))%>%
  select(Distance:value)%>%
  filter(Treatment == "90 FE")### you can filter based on treatment name as in the data for different figures
head(test_data_long) 

test_data_long$Time<-factor(test_data_long$Time, levels = c('1 day', "3 day", "5 day", "10 day", "18 day", "Control" ), labels = c ("1.071 LP (1 day)", "1.214 LP (3 days)", "1.357 LP (5 days)", "1.714 LP (10 days)", "2.285 LP (18 days)", "Inoculated control"))
test_data_long$Size = factor(test_data_long$Size, levels=c('1','2','3', 
                                                          '4','5','6', '7', '8'),
                                  labels = c("1 x 1 focus width","3 x 3 focus width","3 x 7 focus width", "3 x 11 focus width", "3 x 15 focus width","3 x 19 focus width","3 x 23 focus width", "3 x 45 focus width (entire plot)"))

cbp1 <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', "#293352")

plot_hyslp<-ggplot(test_data_long, aes(Distance, value))+
  geom_line(aes(color=Time), size=0.3, alpha=1)+
  geom_point(aes(shape = Time, color = Time), size=0.3)+
  scale_colour_manual(values = cbp1)+
  scale_shape_manual(values = c( 0, 1, 2, 3, 4, 5))+
  xlab("Distance (m)")+
  ylab("Disease prevalence (%)")+
  theme_bw()+
  facet_wrap(.~Size)+
  theme(text = element_text(face='plain', size = 9))+ 
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title.y= element_text(size = 9))+
  theme(axis.title.x = element_text(size = 9))+
  theme(legend.position=c (0.75, 0.15), legend.direction = 'vertical', legend.title = element_blank(), legend.text = element_text (size = 8), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits=c(0, 105), breaks = c(0, 20, 40, 60, 80, 100))+
  scale_x_continuous(limits=c(-15, 15), breaks = c(-15,  -10,  -5,  0,   5,  10, 15))+
  theme(strip.text.x = element_text(size = 7))+ 
  theme(strip.text.x = element_text(margin = margin(2, 2, 2, 2)))

plot_hyslp

ggsave(filename = "Fig.png", width = , height = , units = "in",
       dpi = 600, limitsize = TRUE)
