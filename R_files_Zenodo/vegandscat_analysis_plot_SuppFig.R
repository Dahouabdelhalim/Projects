##########################################################
#   Gawel et al. RSOS Contrasting Ecological Roles      ##
#   of non-native ungulates in a novel ecosystem.       ##
#   Supp Information - sans highest ung values          ##
#   analysis and ggplot                                 ##
##########################################################

###load library###
library(dplyr)
library(ggplot2)
library(lme4)
library(AICcmodavg)
library(MuMIn)

###upload dataset###
vegandsign<-read.csv("vegandsign.csv")


#######################################################
#   Rerun analysis without highest                    # 
#   value for deer and without highest                #
#   value for pigs (scat counts).                     #
#######################################################

vegandsign2<-vegandsign[vegandsign$pig<1,] #without high pig site
vegandsign3<-vegandsign[vegandsign$deer<8,] #without high deer site

##########################################
#         LINEAR MODELS                  #
##########################################

###### Just pigs #########################

linpigtot2<-lm(totalsdl~pig, vegandsign2)
logpigtot2<-lm(totalsdl~log(pig), vegandsign2)
summary(logpigtot2)
summary(linpigtot2)

linpignat2<-lm(nativsdl~(pig), vegandsign2)
logpignat2<-lm(nativsdl~log(pig), vegandsign2)
summary(linpignat2)
summary(logpiglnat2)

linpigexo2<-lm(nonnativesdl~pig, vegandsign2)
logpigexo2<-lm(nonnativesdl~log(pig), vegandsign2)
summary(linpigexo2)
summary(logpigexo2)

linpigvin2<-lm(vines~pig, vegandsign2)
logpigvin2<-lm(vines~log(pig), vegandsign2)
summary(logpigvin2)
summary(logpiglvin2)

###### Just deer ##############
lindeertot3<-lm(totalsdl~deer, vegandsign3)
logdeertot3<-lm(totalsdl~log(deer), vegandsign3)
summary(lindeertot3)
summary(logdeertot3)

lindeernat3<-lm(nativsdl~deer, vegandsign3)
logdeernat3<-lm(nativsdl~log(deer), vegandsign3)
summary(lindeernat3)
summary(logdeernat3)

lindeernnat3<-lm(nonnativesdl~deer, vegandsign3)
logdeernnat3<-lm(nonnativesdl~log(deer), vegandsign3)
summary(lindeernnat3)
summary(logdeernnat3)

lindeervin3<-lm(vines~deer, vegandsign3)
logdeervin3<-lm(vines~log(deer), vegandsign3)
summary(lindeervin3)
summary(logdeervin3)

###PIG PLOTS, LEFT-HAND PANELS###
###r-squared values for input###

summary(lm(totalsdl~pig, vegandsign2))$r.squared

summary(lm(nativsdl~pig, vegandsign2))$r.squared

summary(lm(nonnativesdl~pig, vegandsign2))$r.squared

summary(lm(vines~pig, vegandsign2))$r.squared


###PIG and TOTAL SEEDLINGS###  bquote('x axis'~(Ã…^2))
p1<- ggplot(data=vegandsign2, aes(y=totalsdl, x=pig))+
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
  annotate("text", x = 0.5, y = 800, label = "r^2==0.21", parse=T)


####PIGS AND TWO TRENDLINES FOR NATIVE/EXOTIC SEEDLINGS###
#First, create longform so that one column indicates whether a value comes from native or exotic seedlings#
vegsignlong_p <- gather(vegandsign2, native_exotic, sdls, totalsdl:nonnativesdl, factor_key=TRUE)
vegsignlong_p1 <-subset(vegsignlong_p, native_exotic == "nativsdl"| native_exotic == "nonnativesdl")

p2<-ggplot(vegsignlong_p1, aes(y=sdls, x=pig, color= native_exotic))+
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
  annotate("text", x = 0.4, y = 700, label = "native~r^2==0.21", parse=T)+
  annotate("text", x = 0.4, y = 550, label = "nonnative~r^2==0.09", parse=T)

###PIGS AND VINES####
p3<-ggplot(data=vegandsign2, aes(y=vines, x=pig))+
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
  annotate("text", x = 0.5, y = 200, label = "r^2==0.25", parse=T)

####DEER PLOTS - RIGHT-HAND PANELS#################

###r-squared values for input###
summary(lm(totalsdl~log(deer), vegandsign3))$r.squared

summary(lm(nativsdl~log(deer), vegandsign3))$r.squared

summary(lm(nonnativesdl~log(deer), vegandsign3))$r.squared

summary(lm(vines~log(deer), vegandsign3))$r.squared


###deer and total seedlings####
p4<- ggplot(data=vegandsign3, aes(y=totalsdl, x=deer))+
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
  annotate("text", x = 5, y = 800, label = "r^2==0.69", parse=T)

###deer and separate trendlines for native and non-native## 
#First, create longform so that one column indicates whether a value comes from native or exotic seedlings#
vegsignlong_d <- gather(vegandsign3, native_exotic, sdls, totalsdl:nonnativesdl, factor_key=TRUE)
vegsignlong_d1 <-subset(vegsignlong_d, native_exotic == "nativsdl"| native_exotic == "nonnativesdl")
p5<-ggplot(vegsignlong_d1, aes(y=sdls, x=deer, color= native_exotic))+
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
  annotate("text", x = 5, y = 700, label = "native~r^2==0.62", parse=T)+
  annotate("text", x = 5, y = 550, label = "nonnative~r^2==0.78", parse=T)


###DEER and VINES###
p6<- ggplot(data=vegandsign3, aes(y=vines, x=deer))+
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
  annotate("text", x = 5, y = 220, label = "r^2==0.78", parse=T)


###now that you have individual plots, combine into multiplot using function below###

####MULTI-PANEL####
###create multiplot function###
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
