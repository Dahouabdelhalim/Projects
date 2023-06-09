library(ggplot2)

remove(list = ls())
set.seed(3)
Dating_For_R <- read.delim("Fig3_data_for_R.txt")
### Plot for Dating Runs ###

Dating_For_R$Dataset <- factor(Dating_For_R$Dataset, levels=c("PB1", "PB2", "PB3","PB4", "PB5", "PB6", "PB7", "PB8", "PB9","PB10","PRIOR"),labels=c("Dataset 1", "Dataset 2", "Dataset 3", "Dataset 4", "Dataset 5", "Dataset 6", "Dataset 7", "Dataset 8", "Dataset 9", "Dataset 10", "Prior, No Data"))

Date_plot <- ggplot(Dating_For_R, #setup plot object
                 aes(x = as.factor(Number), y = Age)) #, group = as.factor(Age)
                 
#Add boxplot
Date_plot <- Date_plot + geom_boxplot(outlier.shape=NA,
              notch=FALSE, width=.75, position=position_dodge(3))

#Add points
Date_plot <- Date_plot + geom_point(position = position_jitter(0.2), size=1.1, aes(color=Dataset, shape = Dataset), show.legend = TRUE, na.rm = TRUE) + 
				  scale_shape_manual( values=c(3, 15, 17, 3, 16, 18, 16, 15, 18,  17,  8)) 

#Remove grid background and add borders

Date_plot <- Date_plot + theme_bw(base_size = 19) + theme(panel.border = element_blank(),
				   panel.background = element_rect(fill = "grey98", colour = "grey50"), 
                   panel.grid.major.y = element_line(colour = "grey",size=0.3),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   axis.text = element_text(size=8),
                   axis.ticks.y = element_line(),
                   axis.title = element_text(size=12),
                   legend.title = element_text(size=11),
                   legend.text = element_text(size=10),
                   legend.key.size = unit(1.0, "cm"),
                   plot.margin=unit(c(.5,.5,.5,.5),"cm")) +
                   scale_y_continuous(limit = c(0,75), breaks = seq(0, 80, 10)) +
                   scale_x_discrete( expand = c(.025, .025))

#Adjust axis labels

Date_plot <- Date_plot + labs(x = expression(paste(
  'Node Number')),
  y = expression('Millions of Years'))

#Change the size of the Legend symbol sizes
Date_plot <- Date_plot + guides(color = guide_legend(override.aes = list(size = 3.5)))

Date_plot <- Date_plot + annotate('text', label = 'Prior\\n=105', x = 7, y = 70.5, size = 2.2)
Date_plot <- Date_plot + annotate('text', label = 'Prior\\n=149', x = 8, y = 74.5, size = 2.2)
Date_plot
