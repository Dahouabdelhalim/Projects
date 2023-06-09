##########################################################
#   Gawel et al. RSOS Contrasting Ecological Roles      ##
#   of non-native ungulates in a novel ecosystem.       ##
#   abundance in scats vs. species abundance in forests ##  
##########################################################

###load packages
library(ggplot2)
library(tidyr)
library(dplyr)

###load data
pooportions <- read.csv("pooportions.csv")

##create separate "native" and "exotic" dataframes for figures
native <- subset(pooportions, native_exotic == "native")
exotic <- subset (pooportions, native_exotic == "exotic")


###   CREATE FIGURES   ###

##change facet labels###
facet_names <- c(
  `prop_nature` = "Forest",
  `prop_pig` = "Pig scat",
  `prop_deer` = "Deer scat")

nativeplot <-ggplot(native,aes(plant_spp,proportion))+
  geom_bar(position="dodge",stat="identity")+
  coord_flip()+
  facet_grid(. ~ treatment,labeller = as_labeller(facet_names)) +
  xlab("Species")+
  ylab("Proportional abundance")+
  theme_minimal()+
  theme(axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Native plants")


nonnatplot <-ggplot(exotic,aes(plant_spp,proportion))+
  geom_bar(position="dodge",stat="identity")+
  coord_flip()+
  facet_grid(. ~ treatment,labeller = as_labeller(facet_names)) +
  xlab("Species")+
  ylab("Proportional abundance")+
  theme_minimal()+
  theme(axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Non-native plants")


###Multiplot function###
#credit: Cookbook for R, retrieved at http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist = NULL, file, cols=1, rows=1, layout = NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(nativeplot, nonnatplot, rows=2, cols=1)


