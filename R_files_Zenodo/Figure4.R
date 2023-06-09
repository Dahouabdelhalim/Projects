### Script to  analyze chlamy timelapse data ###

library(ggplot2)
library(dplyr)
library(reshape2)


data <- read.csv(file=("/Users/matt/Documents/Dropbox/LifeCycleRevision/ScientificReports/Dryad/Figure4Data.csv"), header=T)

# Originally (in 2017-04-12 script) I used unweighted means for ordering strains
#  but used weighted means for plotting size distributions
# In this current script, I use weighted means for both plotting and calculating mean progeny sizes
# In other words, I use weighted means for everything in this analysis.


# The following weights cell counts of progeny size data by proportion biomass
# This script repeats a given row by a frequency equal to the progeny size 
#   E.g. num.cells.progeny = 8 will have 8 rows, therefore weighting the
#        count towards these larger but fewer progeny

# First I create a new, separate dataframe from the original
#   containing just strain and number of cells in the progeny
data.new <- data.frame(data[1:4])
#   and I create a last column which is the desired frequency (or biomass) which
#   is equal to the number of cells in the progeny
data.new$freq <- data.new$Num.Cells.Progeny
# I omit NAs which cannot be handled by the next line of code
data.new <- na.omit(data.new)
# This line of code repeats rows by times=$freq
#   which has now created extra data points and weighted the data by biomass
#   E.g. a row of progeny = 8 cells is now repeated 8 times
data.weight <- data.new[rep(row.names(data.new), data.new$freq), 1:4]

# Log-base 2 transforms the weighted data
data.weight$log.progeny <- log(data.weight$Num.Cells.Progeny, 2)

# Calculate WEIGHTED mean progeny sizes and standard deviations for each strain at reproduction
#   aggregate will go through a particular vector and perform function
means.progeny <- aggregate(data.weight$Num.Cells.Progeny ~ data.weight$Strain, FUN = mean)
sd.progeny <- aggregate(data.weight$Num.Cells.Progeny ~ data.weight$Strain, FUN = sd)
colnames(means.progeny) <- c("Strain","Mean.Num.Progeny")
colnames(sd.progeny) <- c("Strain","SD.Num.Progeny")

# Assign weighted mean progeny size and SD to weighted counts by strain
data.weight <- inner_join(data.weight,means.progeny, by = "Strain")
data.weight <- inner_join(data.weight,sd.progeny, by = "Strain")

# Going back to working with the data.new dataframe to plot actual data points
# over boxplots

#   Log transform the unweighted data
data.new$log.progeny <- log(data.new$Num.Cells.Progeny, 2)

# PLOT PROGENY SIZE ACROSS STRAINS 
# note that the default geom_boxplot command plots:
#   lower hinge marks 25th quantile
#   middle hash marks median = 50% quantile
#   upper hinge marks 75th quantile
#   whiskers 

# Manually obtain and input sample sizes above boxplots
samplesizes <- summary(data.new$Strain)
samplesizes

# This plot orders the strains along the x-axis by mean progeny size EXCEPT
#  it has been manually ordered for K101, K106, B506, and B210 because all have a mean of 1

data.weight$Strain <- factor(data.weight$Strain, levels = c('K101','K106',
                                                            'B506','B210',
                                                            'B201','B505',
                                                            'B211','B203',
                                                            'B204'))

data.new$Strain <- factor(data.new$Strain, levels = c('K101','K106',
                                                      'B506','B210',
                                                      'B201','B505',
                                                      'B211','B203',
                                                      'B204'))

myplot <- ggplot(data = data.weight,
                  aes(x=Strain,
                      y=log.progeny)) +
  geom_boxplot(size = 1.2) +
  geom_jitter(data = data.new,
              aes(x=factor(Strain),
                  y=log.progeny),
              size = 0.7, alpha = 0.3, color = 'green4', fill = 'green4') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = 'black'),
        text = element_text(size=20),
        axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black')) +
  scale_y_continuous(limits=c(-0.5,6)) +
  scale_x_discrete(labels = c('K101','K106',
                              'B506','B210',
                              'B201','B505',
                              'B211','B203',
                              'B204')) +
  labs(x="Strain", y=expression("Cells per Propagule (Log"[2]*")")) +
  annotate('text', x=1, y=-0.5, label="N = 668") +
  annotate('text', x=2, y=-0.5, label="162") +
  annotate('text', x=3, y=-0.5, label="162") +
  annotate('text', x=4, y=-0.5, label="158") +
  annotate('text', x=5, y=-0.5, label="119") +
  annotate('text', x=6, y=-0.5, label="704") +
  annotate('text', x=7, y=-0.5, label="51") +
  annotate('text', x=8, y=-0.5, label="115") +
  annotate('text', x=9, y=-0.5, label="103")

ggsave(filename=paste("/Users/matt/Documents/Dropbox/LifeCycleRevision/ScientificReports/Dryad/",
                      paste("Figure4", ".png", sep=""), sep=""), plot=myplot,
       width = 10, height = 4, bg = 'transparent')




