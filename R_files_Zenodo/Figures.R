
#Combinding all plots for figure 2 and supplementary files

Reponsecurves_Maintext <- ggarrange(Snowdaysplot, Humanfootprintplot, Coniferousplot, Openlandscapeplot,
                                    labels = c("a)", "b)", "c)", "d)"),
                                    ncol = 2, nrow = 2)
Reponsecurves_Maintext

ggsave(filename = "Responsecurves Maintext.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 14.6,
       dpi = 450)


Reponsecurves_Suplementary_files <- ggarrange(Waterplot, Soilmoistureplot,
                                    labels = c("a)", "b)"))
                                    
Reponsecurves_Suplementary_files

ggsave(filename = "Responsecurves Supplementary files.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 7.7,
       dpi = 450)


#Plots
Coniferousplot <- ggplot(Coniferoustotal, aes(x))+
  geom_line(aes(y=y1, color = "Individual runs"))+
  geom_line(aes(y=y2, color = "Individual runs"))+
  geom_line(aes(y=y3, color = "Individual runs"))+
  geom_line(aes(y=y4, color = "Individual runs"))+
  geom_line(aes(y=y5, color = "Individual runs"))+
  geom_line(aes(y=y6, color = "Individual runs"))+
  geom_line(aes(y=y7, color = "Individual runs"))+
  geom_line(aes(y=y8, color = "Individual runs"))+
  geom_line(aes(y=y9, color = "Individual runs"))+
  geom_line(aes(y=y10, color = "Individual runs"))+
  geom_line(aes(y=Mean, color = "Mean"),size=2)+
  xlab('Proportion coniferous forest')+
  ylab('Predicted suitability')+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(family = "Helvetica", size = (15)),
          axis.text = element_text(family = "Courier", size = (15)))+


  scale_color_manual(values= c('Grey','Blue', 'Red'))

Coniferousplot
print(Coniferousplot + mynamestheme)


ggsave(filename = "Coniferousplot.jpeg",
  path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 12.6,
       dpi = 300)




Humanfootprintplot <- ggplot(Humanfootprinttotal, aes(x))+
  geom_line(aes(y=y1, color = "Individual runs"))+
  geom_line(aes(y=y2, color = "Individual runs"))+
  geom_line(aes(y=y3, color = "Individual runs"))+
  geom_line(aes(y=y4, color = "Individual runs"))+
  geom_line(aes(y=y5, color = "Individual runs"))+
  geom_line(aes(y=y6, color = "Individual runs"))+
  geom_line(aes(y=y7, color = "Individual runs"))+
  geom_line(aes(y=y8, color = "Individual runs"))+
  geom_line(aes(y=y9, color = "Individual runs"))+
  geom_line(aes(y=y10, color = "Individual runs"))+
  geom_line(aes(y=Mean, color = "Mean"),size=2)+
  xlab('Human footprint index')+
  ylab('Predicted suitability')+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Courier", size = (15)))+
  
  scale_color_manual(values= c('Grey','Blue', 'Red'))

Humanfootprintplot
print(Humanfootprintplot + mynamestheme)

ggsave(filename = "Humanfootprint.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 12.6,
       dpi = 300)

Openlandscapeplot <- ggplot(Openlandscapetotal, aes(x))+
  geom_line(aes(y=y1, color = "Individual runs"))+
  geom_line(aes(y=y2, color = "Individual runs"))+
  geom_line(aes(y=y3, color = "Individual runs"))+
  geom_line(aes(y=y4, color = "Individual runs"))+
  geom_line(aes(y=y5, color = "Individual runs"))+
  geom_line(aes(y=y6, color = "Individual runs"))+
  geom_line(aes(y=y7, color = "Individual runs"))+
  geom_line(aes(y=y8, color = "Individual runs"))+
  geom_line(aes(y=y9, color = "Individual runs"))+
  geom_line(aes(y=y10, color = "Individual runs"))+
  geom_line(aes(y=Mean, color = "Mean"),size=2)+
  xlab('Proportion open landscape')+
  ylab('Predicted suitability')+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Courier", size = (15)))+
  
  scale_color_manual(values= c('Grey','Blue', 'Red'))

Openlandscapeplot
print(Openlandscapeplot + mynamestheme)

ggsave(filename = "Openlandscape.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 12.6,
       dpi = 300)


Snowdaysplot <- ggplot(Snowdaystotal, aes(x))+
  geom_line(aes(y=y1, color = "Individual runs"))+
  geom_line(aes(y=y2, color = "Individual runs"))+
  geom_line(aes(y=y3, color = "Individual runs"))+
  geom_line(aes(y=y4, color = "Individual runs"))+
  geom_line(aes(y=y5, color = "Individual runs"))+
  geom_line(aes(y=y6, color = "Individual runs"))+
  geom_line(aes(y=y7, color = "Individual runs"))+
  geom_line(aes(y=y8, color = "Individual runs"))+
  geom_line(aes(y=y9, color = "Individual runs"))+
  geom_line(aes(y=y10, color = "Individual runs"))+
  geom_line(aes(y=Mean, color = "Mean"),size=2)+
  xlab('Number of snow days')+
  ylab('Predicted suitability')+
  theme_classic()+

  theme(legend.position="none", axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Courier", size = (15)))+
  
  scale_color_manual(values= c('Grey','Blue', 'Red'))

Snowdaysplot
print(Snowdaysplot + mynamestheme)

ggsave(filename = "Snowdays.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 12.6,
       dpi = 300)


Waterplot <- ggplot(Watertotal, aes(x))+
  geom_line(aes(y=y1, color = "Individual runs"))+
  geom_line(aes(y=y2, color = "Individual runs"))+
  geom_line(aes(y=y3, color = "Individual runs"))+
  geom_line(aes(y=y4, color = "Individual runs"))+
  geom_line(aes(y=y5, color = "Individual runs"))+
  geom_line(aes(y=y6, color = "Individual runs"))+
  geom_line(aes(y=y7, color = "Individual runs"))+
  geom_line(aes(y=y8, color = "Individual runs"))+
  geom_line(aes(y=y9, color = "Individual runs"))+
  geom_line(aes(y=y10, color = "Individual runs"))+
  geom_line(aes(y=Mean, color = "Mean"),size=2)+
  xlab('Water availability')+
  ylab('Predicted suitability')+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Courier", size = (15)))+
  
  scale_color_manual(values= c('Grey','Blue', 'Red'))

Waterplot
print(Waterplot + mynamestheme)

ggsave(filename = "Water availability2.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 12.6,
       dpi = 300)

Soilmoistureplot <- ggplot(Soilmoisturetotal, aes(x))+
  geom_line(aes(y=y1, color = "Individual runs"))+
  geom_line(aes(y=y2, color = "Individual runs"))+
  geom_line(aes(y=y3, color = "Individual runs"))+
  geom_line(aes(y=y4, color = "Individual runs"))+
  geom_line(aes(y=y5, color = "Individual runs"))+
  geom_line(aes(y=y6, color = "Individual runs"))+
  geom_line(aes(y=y7, color = "Individual runs"))+
  geom_line(aes(y=y8, color = "Individual runs"))+
  geom_line(aes(y=y9, color = "Individual runs"))+
  geom_line(aes(y=y10, color = "Individual runs"))+
  geom_line(aes(y=Mean, color = "Mean"),size=2)+
  xlab('Proportion Soilmoisture')+
  ylab('Predicted suitability')+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Courier", size = (15)))+
  

  scale_color_manual(values= c('Grey','Blue', 'Red'))

Soilmoistureplot
print(Soilmoistureplot + mynamestheme)


ggsave(filename = "Soilmoisture.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 20,
       height = 12.6,
       dpi = 300)

