##########################################################
#   Gawel et al. RSOS Contrasting Ecological Roles      ##
#   of non-native ungulates in a novel ecosystem.       ##
#   Scat abundance and forest community abundances      ##
#   analysis and ggplot                                 ##
##########################################################

###Load data######
vegandsign<-read.csv("vegandsign.csv")

###load packages###
library(ggplot2)
library(lme4)
library(AICcmodavg)
library(MuMIn)
library(tidyr)

#############################################
#             LINEAR MODELS                 #
#############################################

##AICc values reported in supplementary materials##

##First, test lm with different response variables [total seedlings, native and exotic seedlings, and vines] with both pig and deer scat densities:

##LOOK AT MODELS PREDICTING TOTAL SEEDLING ABUNDANCE FIRST##
bothtot <-lm(totalsdl~pig*deer, vegandsign)
nointxntot<-lm(totalsdl~pig+deer, vegandsign)
pigtot <-lm(totalsdl~pig, vegandsign)
deertot <-lm(totalsdl~deer, vegandsign)
logpigtot<-lm(totalsdl~log(pig), vegandsign)
logdeertot<-lm(totalsdl~log(deer), vegandsign)

aictab(list(bothtot, nointxntot, pigtot, deertot,logpigtot,logdeertot))

###TEST NATIVE AND NON-NATIVE SEEDLINGS SEPARATELY###
##NATIVE##
bothnat <- lm(nativsdl~pig+deer, vegandsign)
pignat <-lm(nativsdl~pig, vegandsign)
deernat <-lm(nativsdl~deer, vegandsign)
logpignat<-lm(nativsdl~log(pig), vegandsign)
logdeernat<-lm(nativsdl~log(deer), vegandsign)

aictab(list(bothnat, pignat, deernat, logpignat, logdeernat))

##NON-NATIVE##
bothexo <- lm(nonnativesdl~pig+deer, vegandsign)
pigexo <-lm(nonnativesdl~pig, vegandsign)
deerexo <-lm(nonnativesdl~deer, vegandsign)
logpigexo<-lm(nonnativesdl~log(pig), vegandsign)
logdeerexo<-lm(nonnativesdl~log(deer), vegandsign)

aictab(list(bothexo, pigexo, deerexo, logpigexo, logdeerexo))

###TEST MODELS WITH VINE ABUNDANCE###
bothvin <- lm(vines~pig+deer, vegandsign)
pigvin <-lm(vines~pig, vegandsign)
deervin <-lm(vines~deer, vegandsign)
logpigvin<-lm(vines~log(pig), vegandsign)
logdeervin<-lm(vines~log(deer), vegandsign)
aictab(list(bothvin, pigvin, deervin, logpigvin, logdeervin))

##########################################################
#                     FIGURES                            #
##########################################################

###r-squared values for input###
summary(lm(totalsdl~pig, vegandsign))$r.squared

summary(lm(nativsdl~pig, vegandsign))$r.squared

summary(lm(nonnativesdl~pig, vegandsign))$r.squared

summary(lm(vines~pig, vegandsign))$r.squared

r01 <- paste("r^2 <", "0.01")
r02 <- paste("r^2 ==", "0.02")

###PIGS and TOTAL SEEDLINGS### 
p1<- ggplot(data=vegandsign, aes(y=totalsdl, x=pig))+
  geom_point()+
  xlab(expression(paste("Pig scats per 100", m^2)))+
  ylab("Total seedling abundance")+
  theme_bw() +
  theme(axis.line = element_line())+
  theme (panel.background = element_rect(fill = "transparent", colour = NA),
         panel.grid.minor = element_line(colour = NA),
         panel.grid.major = element_line(colour = NA), legend.position="none")+
  theme(axis.title.x=element_text(size=9))+
  theme(axis.title.y=element_text(size=9))+
  annotate("text", x = 0.75, y = 800, label = "r^2<0.01", parse=T)

####PIGS AND TWO TRENDLINES FOR NATIVE/EXOTIC SEEDLINGS###
#First, create longform so that one column indicates whether a value comes from native or exotic seedlings#
vegsignlong <- gather(vegandsign, native_exotic, sdls, totalsdl:nonnativesdl, factor_key=TRUE)
vegsignlong1 <-subset(vegsignlong, native_exotic == "nativsdl"| native_exotic == "nonnativesdl")

p2<-ggplot(vegsignlong1, aes(y=sdls, x=pig, color= native_exotic))+
  geom_point()+
  xlab(expression(paste("Pig scats per 100", m^2)))+
  ylab("Seedling abundance")+
  theme_bw() +
  theme(axis.line = element_line())+
  theme (panel.background = element_rect(fill = "transparent", colour = NA),
         panel.grid.minor = element_line(colour = NA),
         panel.grid.major = element_line(colour = NA), legend.position="none")+
  scale_colour_manual(name = "Nat/Exo",
                      labels = c("Exotic", "Native"),
                      values = c("gray", "black")) +
  scale_shape_manual(name = "Nat/Exo",
                     labels = c("Exotic", "Native"),
                     values = c(19, 17))+
  theme(axis.title.x=element_text(size=9))+
  theme(axis.title.y=element_text(size=9))+
  annotate("text", x = 0.7, y = 700, label = "native~r^2<0.01", parse=T)+
  annotate("text", x = 0.7, y = 600, label = "nonnative~r^2==0.02", parse=T)

###PIGS AND VINES####
p3<-ggplot(data=vegandsign, aes(y=vines, x=pig))+
  geom_point()+
  xlab(expression(paste("Pig scats per 100", m^2)))+
  ylab("Vine abundance")+
  theme_bw() +
  theme(axis.line = element_line())+
  theme (panel.background = element_rect(fill = "transparent", colour = NA),
         panel.grid.minor = element_line(colour = NA),
         panel.grid.major = element_line(colour = NA), legend.position="none")+
  theme(axis.title.x=element_text(size=9))+
  theme(axis.title.y=element_text(size=9))+
  annotate("text", x = 0.75, y = 200, label = "r^2<0.01", parse=T)

###DEER AND TOTAL SEEDLINGS###
###r-squared values for input###
summary(lm(totalsdl~log(deer), vegandsign))$r.squared

summary(lm(nativsdl~log(deer), vegandsign))$r.squared

summary(lm(nonnativesdl~log(deer), vegandsign))$r.squared

summary(lm(vines~log(deer), vegandsign))$r.squared

p4<- ggplot(data=vegandsign, aes(y=totalsdl, x=deer))+
  geom_smooth(se=FALSE, colour="black",method="lm",formula=y~log(x))+
  geom_point()+
  xlab(expression(paste("Deer scats per 100", m^2)))+
  ylab("Total seedling abundance")+
  theme_bw() +
  theme(axis.line = element_line())+
  theme (panel.background = element_rect(fill = "transparent", colour = NA),
         panel.grid.minor = element_line(colour = NA),
         panel.grid.major = element_line(colour = NA), legend.position="none")+
  theme(axis.title.x=element_text(size=9))+
  theme(axis.title.y=element_text(size=9))+
  annotate("text", x = 7.2, y = 800, label = "r^2==0.71", parse=T)

####DEER AND TWO TRENDLINES FOR NATIVE/EXOTIC SEEDLINGS###
#Use long dataset created above on lines 95-96#
p5<-ggplot(vegsignlong1, aes(y=sdls, x=deer, color= native_exotic))+
  geom_smooth(se=FALSE, method="lm",formula=y~log(x))+
  geom_point()+
  xlab(expression(paste("Deer scats per 100", m^2)))+
  ylab("Seedling abundance")+
  theme_bw() +
  theme(axis.line = element_line())+
  theme (panel.background = element_rect(fill = "transparent", colour = NA),
         panel.grid.minor = element_line(colour = NA),
         panel.grid.major = element_line(colour = NA), legend.position="none")+
  scale_colour_manual(name = "Nat/Exo",
                      labels = c("Exotic", "Native"),
                      values = c("gray", "black")) +
  scale_shape_manual(name = "Nat/Exo",
                     labels = c("Exotic", "Native"),
                     values = c(19, 17))+
  theme(axis.title.x=element_text(size=9))+
  theme(axis.title.y=element_text(size=9))+
  annotate("text", x = 7.2, y = 700, label = "native~r^2==0.64", parse=T)+
  annotate("text", x = 7.2, y = 600, label = "nonnative~r^2==0.79", parse=T)

###DEER AND VINES###
p6<- ggplot(data=vegandsign, aes(y=vines, x=deer))+
  geom_smooth(se=FALSE, colour="black",method="lm",formula=y~log(x))+
  geom_point()+
  xlab(expression(paste("Deer scats per 100", m^2)))+
  ylab("Vine abundance")+
  theme_bw() +
  theme(axis.line = element_line())+
  theme (panel.background = element_rect(fill = "transparent", colour = NA),
         panel.grid.minor = element_line(colour = NA),
         panel.grid.major = element_line(colour = NA), legend.position="none")+
  theme(axis.title.x=element_text(size=9))+
  theme(axis.title.y=element_text(size=9))+
  annotate("text", x = 7.2, y = 220, label = "r^2==0.79", parse=T)

###Multiplot function###
#credit: Cookbook for R, retrieved at http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist = NULL, file, cols = 2, layout = NULL) {
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

###MULTIPLOT###
multiplot(p1,p2,p3,p4,p5,p6,cols=2)
