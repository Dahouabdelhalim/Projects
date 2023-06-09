rm(list=ls())
# Load packages into session 
packages.vec <- c("sp", "fields", "rgeos", "maps", "ggmap", "rgdal", "gstat", "dplyr",
                  "ggplot2", "scales", "magrittr", 'gridExtra', 'splancs', 'parallel', 'maptools', 'automap'
                  , 'plyr', 'plotrix', 'readr')
lapply(packages.vec, require, character.only=TRUE)

# read in data
dat <- read_csv("All_Ancestry_K17_and_IndRefSet_Updated.csv")
# remove points where coordinates are NA
dat2 <- dat[!is.na(dat$Long) & ! is.na(dat$Lat), ] 
# remove points outside USA 
top = 49.3457868 # north lat
left = -124.7844079 # west long
right = -66.9513812 # east long
bottom = 24.7433195 # south lat
dat3 <- dat2[dat2$Lat < top & dat2$Lat > bottom & dat2$Long > left & dat2$Long < right,]
dat2 <- dat3
# create commercial breed column
dat2$commercial_breed <- dat2$Berkshire + dat2$Hampshire + 
  dat2$Duroc_Other.Breeds + dat2$Landrace + dat2$Yorkshires
dat2[dat2$State=="NJ", "commercial_breed"] <- 0.01

#- convert dataframe to spatial points datafrmae (SPDF)
# specify which points are the coordinates
coordinates(dat2) <- ~ Long + Lat
# set projection of the data
dat2@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
# display the four counts that denote the spatial extent of the data
bbox(dat2); dat2@bbox

# remove excess points
usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +ellps=WGS84  +datum=WGS84"))
staying <- dat2[!is.na(over(dat2, as(usa, "SpatialPolygons"))),]
dat2 <- staying

# dealing with colors
num.colors <- 20
blue.end.pnts <- c("#c6dbef","#08306b")
red.end.pnts <- c("#ffeda0","#730022")
# function for color
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return(colorRampPalette(colors,interpolate = c("spline"))(colsteps)[ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] 
  )
}

# Blue Palette (cut point at .25 so 0 to .25)
blue.palette <- color.gradient(1:5,colors=blue.end.pnts)

# Red Palette (.25 to 1)
red.palette <- color.gradient(1:15,colors=red.end.pnts)

# Combine Palettes into one
red.blue.palette <- c(rev(blue.palette),red.palette)

#----- mapping Western Heritage, Wild Boar, and Commercial Breeds
#-- Western Heritage:: Calling Western Heritage "NE" for convenience
x <- dat2@coords
y <- dat2$Western.Heritage
#tps1 <- Tps(dircos(x), y)
tps1 <- Tps(x, y)

# extracting values and setting up graphical parameters
dfNE <- data.frame(tps1$x, tps1$fitted.values)
colnames(dfNE) <- c('x', 'y', 'z')
state <- map_data("state")
state$group2 <- rep(1, nrow(state))
pts_sz <- 4
dfNE <- dfNE[order(dfNE$z),]
# plot with ggplot
plotNE <- dfNE %>% as.data.frame %>%
  ggplot(aes(x=x, y=y))+
  geom_polygon(data=state, aes(x=long, y=lat, group=group, fill="white", colour='white'), fill='white',col="black") + 
  geom_point(aes(colour=z), size=pts_sz)+ 
  coord_equal() +
  scale_colour_gradientn(colors=red.blue.palette, breaks=seq(0,1,by=.05), guide=F, limits=c(0,1))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  guides(fill=FALSE, point=FALSE) + 
  ggtitle('A. Western heritage breeds') +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())
plotNE

#-- Wild Boar Heritage
x <- dat2@coords
y <- dat2$European.Wild.Boar
tps3 <- Tps(x, y)

# extracting values and setting up graphical parameters
dfWB <- data.frame(tps3$x, tps3$fitted.values)
colnames(dfWB) <- c('x', 'y', 'z')
dfWB <- dfWB[dfWB$z > 0.001,]
state <- map_data("state")
state$group2 <- rep(1, nrow(state))
dfWB <- dfWB[order(dfWB$z),]

# plot with legend
plotWB <- dfWB %>% as.data.frame %>%
  ggplot(aes(x=x, y=y))+
  geom_polygon(data=state, aes(x=long, y=lat, group=group, fill="white", colour='white'), fill='white',col="black") + 
  geom_point(aes(colour=z), size=pts_sz)+ 
  coord_equal() +
  scale_colour_gradientn(colors=red.blue.palette, breaks=seq(0,1,by=.2), 
                         #*** if making legend, use next line ***#
                         guide="colourbar", 
                         limits=c(0,1))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  labs(colour='Proportion') + 
  guides(fill=FALSE, point=FALSE) + 
  ggtitle('B. European wild boar') +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())
#pdf("WB_withLegend.pdf")
plotWB
#dev.off()
theme_get()$legend.title$hjust

# plot WB without legend
plotWB <- dfWB %>% as.data.frame %>%
  ggplot(aes(x=x, y=y))+
  geom_polygon(data=state, aes(x=long, y=lat, group=group, fill="white", colour='white'), fill='white',col="black") + 
  geom_point(aes(colour=z), size=pts_sz)+ 
  coord_equal() +
  scale_colour_gradientn(colors=red.blue.palette, breaks=seq(0,1,by=.05), guide=F, limits=c(0,1))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  labs(colour='') + 
  guides(fill=FALSE, point=FALSE) + 
  ggtitle('B. European wild boar') +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())
plotWB

#-- commercial breed
x <- dat2@coords
y <- dat2$commercial_breed
tps4 <- Tps(x, y)

# extracting values and setting up graphical parameters
dfCB <- data.frame(tps4$x, tps4$fitted.values)
colnames(dfCB) <- c('x', 'y', 'z')
state <- map_data("state")
state$group2 <- rep(1, nrow(state))
dfCB <- dfCB[order(dfCB$z),]

# plot CB
plotCB <- dfCB %>% as.data.frame %>%
  ggplot(aes(x=x, y=y))+
  geom_polygon(data=state, aes(x=long, y=lat, group=group, fill="white", colour='white'), fill='white',col="black") + 
  geom_point(aes(colour=z), size=pts_sz)+ 
  coord_equal() +
  scale_colour_gradientn(colors=red.blue.palette, breaks=seq(0,1,by=.05), guide=F, limits=c(0,1))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  labs(colour='') + 
  guides(fill=FALSE, point=FALSE) + 
  ggtitle('C. Commercial breeds') +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())
plotCB

# plot all panels
#pdf('heritage_3panels.pdf')
grid.arrange(plotNE, plotWB, plotCB, ncol=1)
#dev.off()
