###############################################################
###### Disease Gradient plots #######
############ BY Bhim Chaulagain #################
############### Oregon State University ######################

# install required packages
install.packages(c("readxl", "reshape2", "cowplot", "scales", "grid", "tidyverse"))

# Load required packages
library(readxl)
library(grid)
library(scales)
library(cowplot)
library(reshape2)
library(tidyverse)

# read in data
data<-read_excel("Chaulagain_et_al_2022_Data_fig_2_3_4_5", sheet = "data")
head(data)

#reshape data
test_data_long <- melt(data, id=c("Treatment", "Distance", "Location", 'Type'))%>%
  select(Treatment:value)%>%
  filter(Location == "M" & Type == 'F') # filter the data based on location and field/simulation studies as in data file for different figures
head(test_data_long) 


test_data_long$variable<-factor(test_data_long$variable, levels = c('2019', "2020"), labels = c ("Year 2019", "Year 2020"))
test_data_long$Treatment = factor(test_data_long$Treatment, levels=c('Inoculated control','CP cull at 3 days after first sporulation','CP cull ring 3', 
                                                                     'CP cull at 1 day after first sporulation','Fungicide protection of CP at 3 days', 'CP cull at 3 days_FP of non-CP', 
                                                                     'Fungicide protection of entire plot at 3 days'),
                                  labels = c("T1: Inoculated Control","T2: CP cull at 3 days","T3: large CP cull at 3 days","T4: CP cull at 1 day","T5: FP of CP at 3 days",
                                              "T6: CP cull + FP of non-CP at 3 days", "T7: FP of entire plot at 3 days"))



cbp2<-c("red", "blue")

plot_hyslp<-ggplot(test_data_long, aes(Distance, value))+
  geom_line(aes(color=variable), size=0.5, alpha=1)+
  geom_point(aes(shape = variable), size=0.6)+
  scale_colour_manual(values = cbp2)+
  scale_shape_manual(values = c(1, 2))+
  xlab("Distance (m)")+
  ylab("Disease prevalence (%)")+
  theme_bw()+
  facet_wrap(.~Treatment)+
  theme(text = element_text(face='plain', size = 9))+ 
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title.y= element_text(size = 9))+
  theme(axis.title.x = element_text(size = 9))+
  theme(legend.position=c (0.5, 0.25), legend.direction = 'horizontal', legend.title = element_blank(), legend.text = element_text (size = 9), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits=c(0, 105), breaks = c(0, 20, 40, 60, 80, 100))+
  scale_x_continuous(limits=c(-15, 15), breaks = c(-15,  -10,  -5,  0,   5,  10, 15))+
  theme(strip.text.x = element_text(size = 7))+
  theme(strip.text.x = element_text(margin = margin(2, 2, 2, 2)))

plot_hyslp


#Figure used in Inset 1
newdata <- test_data_long%>%
  select(Treatment, Distance, Location, Type, variable, value)%>%
  filter(Treatment =='T6: CP cull + FP of non-CP at 3 days')
head(newdata)

inset1 <- ggplot(newdata, aes(Distance, value))+
  geom_line(aes(color=variable), size=0.5, alpha=1)+
  geom_point(aes(shape = variable), size=0.6)+
  scale_colour_manual(values = cbp2)+
  scale_shape_manual(values = c(1, 2))+
  xlab("Distance (m)")+
  ylab("Log (Disease prevalence (%))")+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title.y= element_text(size = 6))+
  theme(axis.title.x = element_text(size = 6))+
  theme(legend.position='none', legend.title = element_blank(), legend.text = element_text (size = 9), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous (trans='log10', labels = scientific)+
  scale_x_continuous(limits=c(-15, 15), breaks = c(-15,  -10,  -5,  0,   5,  10, 15))+
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))

inset1

#Figure used in Inset 2
newdata2 <- test_data_long%>%
  select(Treatment, Distance, Location, Type, variable, value)%>%
  filter(Treatment =='T7: FP of entire plot at 3 days')
head(newdata2)
inset2 <- ggplot(newdata2, aes(Distance, value))+
  geom_line(aes(color=variable), size=0.5, alpha=1)+
  geom_point(aes(shape = variable), size=0.6)+
  scale_colour_manual(values = cbp2)+
  scale_shape_manual(values = c(1, 2))+
  xlab("Distance (m)")+
  ylab("Log (Disease prevalence (%))")+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title.y= element_text(size = 6))+
  theme(axis.title.x = element_text(size = 6))+
  theme(legend.position='none', legend.title = element_blank(), legend.text = element_text (size = 9), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous (trans='log10', labels = scientific)+
  scale_x_continuous(limits=c(-15, 15), breaks = c(-15,  -10,  -5,  0,   5,  10, 15))+
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))
inset2

####Draw rectangles (to reflect data points not used for AUDPC estimation)
#### you may need to adjust the position of the shape for each location
rect1 <- rectGrob(
  x = unit(1.6, "in"),
  y = unit(1, "npc") - unit(0.23, "in"),
  width = unit(0.51, "in"),
  height = unit(1.525, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

rect2 <- rectGrob(
  x = unit(4.437, "in"),
  y = unit(1, "npc") - unit(0.23, "in"),
  width = unit(0.51, "in"),
  height = unit(1.525, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

rect3 <- rectGrob(
  x = unit(7.27, "in"),
  y = unit(1, "npc") - unit(0.23, "in"),
  width = unit(0.51, "in"),
  height = unit(1.525, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

rect4 <- rectGrob(
  x = unit(1.6, "in"),
  y = unit(1, "npc") - unit(1.98, "in"),
  width = unit(0.51, "in"),
  height = unit(1.525, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

rect5 <- rectGrob(
  x = unit(4.437, "in"),
  y = unit(1, "npc") - unit(1.98, "in"),
  width = unit(0.51, "in"),
  height = unit(1.525, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

rect6 <- rectGrob(
  x = unit(7.4, "in"),
  y = unit(1, "npc") - unit(2.169, "in"),
  width = unit(0.29, "in"),
  height = unit(0.867, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

rect7 <- rectGrob(
  x = unit(1.81, "in"),
  y = unit(1, "npc") - unit(3.935, "in"),
  width = unit(0.289, "in"),
  height = unit(0.857, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "skyblue2", alpha = 0.5)
)

ggdraw(plot_hyslp)+ 
  draw_plot(inset1, .7, .4, .23, .23)+
  draw_plot(inset2, .08, .09, .23, .23)+
  draw_grob(rect1)+
  draw_grob(rect2)+
  draw_grob(rect3)+
  draw_grob(rect4)+
  draw_grob(rect5)+
  draw_grob(rect6)+
  draw_grob(rect7)

ggsave(filename = "Fig.png", width = , height = , units = "in",
       dpi = 600, limitsize = TRUE)

