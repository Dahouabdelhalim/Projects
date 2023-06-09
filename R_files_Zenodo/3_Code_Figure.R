# ****************************************************************************
# MODEL CODE FOR:
#   "General rules for environmental management to prioritise social-ecological
#   systems research based on a value of information approach"
# AUTHORS:
#   Katrina J Davis, Iadine Chadès, Jonathan R Rhodes & Michael Bode
# DESCRIPTION:
#   Code recreates Figure 4 in the manuscript (but at lower discretisation)
# NOTE: 
#   This code should be run after Matlab file 'Code_EVPXI'
# ****************************************************************************

# Update with personal file path
workDir <- ""

# Set working directory
setwd(workDir) 

# libraries
library(grDevices) # colours

# clear environment
rm(list=ls())

# Color settings ---------------------------------------------------------

# Colours
DBlue = "#1D7BE5" # Social
LBlue = "#00B0F0" # Social-ecological
LGreen = "#64C39B" # Ecological
White = "White"

# Order for figure: Soc, Soc, Ecol, Ecol, Soc-Ecol, space, Soc, Ecol, Soc-Ecol
color <- c(DBlue,DBlue,LGreen,LGreen,LBlue,White,DBlue,LGreen,LBlue)


# Aggregate systems data ------------------------------------------------------------------


# Read in systems data
S1 <- read.csv("S1.csv", header = F) # System 1
S2 <- read.csv("S2.csv", header = F) # System 2
S3 <- read.csv("S3.csv", header = F) # System 3
S4 <- read.csv("S4.csv", header = F) # System 4

# Combine data into single dataframe
d <- rbind(S1, S2, S3, S4)

# Specify column names
colnames(d) <- c("I","q","r","C","H","Soc","Eco","Soc-Eco")
# Add an empty column
d$No <- NA
# Re-order
d <- d[,c("I","q","r","C","H","No","Soc","Eco","Soc-Eco")]

# Data to percentages
d=d*100

# Create Manuscript Figure 4 as PDF ---------------------------------------------------------------------

# Individual plot titles
mains <- c("1. Territorial use rights fishery", 
           "2. Salmon fishery", 
           "3. Agricultural system", 
           "4. Non-timber forest products")

pdf(file = "EVPXI.pdf", width = 3.5, height = 8) # call pdf (units in inches)
par(mfrow=c(4,1), # multipane plot
    mar=c(0,5,2,0.5), # inner plot margins (space btw plots)
    oma=c(5,0,0,0)) # outer margins (space outside plots)

j <- 0 # initialise counter j
for (i in c(2,5,8,11)){
  j=j+1
  barCenters <- barplot(as.numeric(d[i,]), col=color,
                        names.arg = colnames(d),
                        xaxt="n",
                        axes=T,
                        ylim = c(0, 100),
                        main = mains[j],
                        ylab = "", cex.lab=1, cex.main=1.25, cex.names=0.5, cex.axis=1)
  segments(6.7, 0, 6.7, 85, lwd=1.5, lty=2, col="gray")
  
  if (j==4) {
    par(xpd=NA)
    text(x=barCenters[1:5],y = par("usr")[3] - 9,labels=colnames(d[1:5]), 
         cex = 1, pos = 1, offset = 0.01)
    text(x=barCenters[7:9],y = par("usr")[3] - 9,labels=colnames(d[7:9]), 
         cex = 1, pos = 1, offset = 0.3, srt=45)
  }
  
  segments(barCenters,as.matrix(d[i-1,]), barCenters, as.matrix(d[i+1,]), lwd = 1)
  arrows(barCenters, as.matrix(d[i-1,]), barCenters, as.matrix(d[i+1,]), lwd = 1, 
         angle = 90,code = 3, length = 0.05)
}
par(xpd=NA)

# y axis label
text(x = -2, y = 240, labels = "Expected value of partial information (%)", srt=90, cex=1.25)

dev.off() # close graph