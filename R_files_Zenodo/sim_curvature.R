#!/usr/bin/env Rscript

args = commandArgs(TRUE)

if( length(args) != 2 ){stop("Usage: Rscript simhair.R <radius> <output_dir_path> \\n or \\n for i in {1..1000}; do Rscript simhair.R $i data/simArcs/; done")}

#specify radius of the arc in mm and output file name
radius = as.numeric(args[1])
output_dir = args[2]

#load packages that are necessary
require(ggplot2)
require(dplyr)
require(data.table)
require(tidyr)


#no. of hair to simulate - 25 for now
nhair = 25

#pick different lengths of hair between 1.5 and 3mm
hair.lengths = runif(nhair, min = 1.0, max = 3.5)

#pick a starting angles for hair segments
start.theta = runif(nhair, min = 0, max = pi)

#define angle of the arc (in radians)
#arc.angle = pi/(2*radius)
arc.angle = hair.lengths/radius

#set end value of the angle
end.theta = start.theta + arc.angle


#create labels for each arc
dat = data.table(start.theta = start.theta,
                 end.theta = end.theta,
                 arc.name = paste("arc_",seq(1,25),sep=""))

#function to generate arc given the start and end angles
apoints = function(stheta,etheta){
  rthetas = seq(stheta,etheta,
                length.out = 25)
  x = radius*cos(rthetas)
  y = radius*sin(rthetas)
  return(data.table(x = x, y = y))}

#generate arcs
dat = dat%>%
  group_by(arc.name)%>%
  mutate(dats = list(apoints(start.theta,
                             end.theta)))%>%
  unnest(dats)

#center the arcs so they appear at the center of each 'window'
dat = dat %>%
  mutate(x2 = x - mean(x),
         y2 = y - mean(y))

#set limits for plots - the length of the arc/line
#the length of each subpanel is 6 units
#the length of entire plot is 6*5 = 30 units
#since the plot will be 30mm, each unit = 1 mm
xlims = ylims = c(-3,3)

#plot size in 1 dimension (mm)
plt.size = 3960/132

#plot resolution.
#we need the plot to be 3960 x 3960 pixels with a resolution of 132 pixels per mm
plt.resolution = 132 * 25.4 # ppi

#plot the windows altogether
line.width = (1/3900) * plt.size * .pt

output_filename = paste(output_dir,
                        "/simArc_",
                        radius, ".tiff", sep = "")

tiff(output_filename, units="mm",
     width = plt.size, height = plt.size,
     res = plt.resolution)

ggplot(data = dat)+
  geom_path(aes(x2, y2),
            size = line.width)+
  theme_bw()+
  facet_wrap(.~arc.name,
             scale = "fixed")+
  coord_cartesian(xlim = xlims, ylim = ylims)+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    panel.spacing = unit(0,"lines"),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0,0), "line"))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

dev.off()



