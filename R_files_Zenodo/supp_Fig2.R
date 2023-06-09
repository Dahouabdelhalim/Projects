# load packages
library("PerformanceAnalytics")

# load data 
d <- read.csv("survival_dataset.csv", stringsAsFactors = F)

# subset data
dd <- d[c("gg_mean_s", "gr_mean_s", "sg_mean_s", "sr_mean_s", "f_mean_s")]

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  r <- round(cor(x, y), digits=2)
  txt <- r 
  cex.cor <-  2 
  text(0.05, 0.8, txt, cex = cex.cor) 
}

# Customize upper panel
upper.panel.D<-function(x, y){
  points(x,y, pch = 20) 
}

# Figure 2
pairs(dd, 
      lower.panel = panel.cor,
      upper.panel = upper.panel.D,
      labels = c("Grooming \\ngiving", "Grooming \\nreceiving", 
                 "Support \\ngiving", "Support \\nreceiving",
                 "Foraging \\nin prox"))


