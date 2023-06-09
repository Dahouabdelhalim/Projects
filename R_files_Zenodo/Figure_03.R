rm(list=ls())

########
# Xing et al., Figure 03
########


###goal: 

##compare sleral ring and orbit morphology of Oculudentavis to modern birds, lizards, and other Mesozoic birds (data from Schmitz and Motani 2011, Science)

#differentiate between
#diurnal birds, nocturnal birds, cathemeral birds, => circles
#fossil birds => open squares 
#diurnal lizards, nocturnal lizards => triangles
#HPG-15-3 left and right => filled squares


##### 1 Preliminaries

#please set working directory
setwd(path/to/your/working/directory)

#load data
dat0 <- read.csv("morphospace.csv", header=T)
head(dat0)
dat <- dat0

#call libraries
library(ggplot2)
library(ggpubr)
library(cowplot) 
library(data.table)

##### 2 Preparing plot

#log10-transform the geometric mean
dat$geomm <- log10(dat$geomm); head(dat)

amber <- dat[1:2,]
fossil <- dat[(dat$dap == "unknown") & (dat$group == "fossil"), ]
fossilr <- as.numeric(rownames(fossil))
dat <- dat[-fossilr,]
dat <- dat[-c(1,2),]

###### 3 final plot without sensitivity analysis

pdf("Fig_03.pdf", 8, 6)
sp <- ggscatter(dat, x="geomm", y="opt",
                color = "dap", shape = "group",
                palette = "Dark2",
                ellipse = TRUE, #mean.point = TRUE,
                star.plot = TRUE,
                size = 4, alpha = 0.4) +
  border(size=1.2)   


spplus <- sp   +
  geom_point(data=as.data.frame(fossil[5:6]), shape=15, size=4) +
  geom_point(data=as.data.frame(amber[5:6]), shape=18, size=5)

pmain <- spplus

xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = dat, aes(x = geomm, fill = dap),
               alpha = 0.4, size = 0.2)+
  ggpubr::fill_palette("Dark2")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = dat, aes(x = opt, fill = dap),
               alpha = 0.4, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("Dark2")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
dev.off()


############################


###sensitivity analysis
oculudentavis <- dat0[1:2,]
avg.oculudentavis <- apply(oculudentavis[,5:6], 2, mean)

#range for geomm
from <- 0.7*avg.oculudentavis[1]
to <- 1.3*avg.oculudentavis[1]
length.out <- 10
geomm <- seq(from, to, by = ((to - from)/(length.out - 1)))

#range for opt
from_opt <- 0.7*avg.oculudentavis[2]
to_opt <- 1.3*avg.oculudentavis[2]
length.out <- 10
opt <- seq(from_opt, to_opt, by = ((to_opt - from_opt)/(length.out - 1)))

#distribution
simulated <- data.frame(geomm=geomm, opt=opt)
simulated[,1] <- log10(simulated[,1])

all.simulated <- CJ(geomm=simulated[,1], opt=simulated[,2], unique=T)


#plot everything, including sensitivity analysis
pdf("Fig_03_sensitivity.pdf", 8, 6)
sp <- ggscatter(dat, x="geomm", y="opt",
                color = "dap", shape = "group",
                palette = "Dark2",
                ellipse = TRUE, #mean.point = TRUE,
                star.plot = TRUE,
                size = 4, alpha = 0.4) +
  border(size=1.2)   


spplus <- sp   +
  geom_point(data=as.data.frame(fossil[5:6]), shape=15, size=4) +
  geom_point(data=as.data.frame(amber[5:6]), shape=18, size=5) +
  #geom_point(data=all.simulated, shape=15, size=1) +
  stat_ellipse(data=all.simulated)

pmain <- spplus

xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = dat, aes(x = geomm, fill = dap),
               alpha = 0.4, size = 0.2)+
  ggpubr::fill_palette("Dark2")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = dat, aes(x = opt, fill = dap),
               alpha = 0.4, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("Dark2")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
dev.off()
