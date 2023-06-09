rm(list=ls())

library(plyr)
library(sciplot)


# Set new working directory to save plots
setwd(".../Data_and_Code")

rcampdata <- read.csv("rhodo_frequency_balanced.csv")
abiesdata <- read.csv("abies_frequency_balanced.csv")

head(rcampdata)
head(abiesdata)


# ======================================================================
#  BARPLOT - COMPOSITE PLOT OF ABUNDANCE AND MORTALITY OF TWO SPECIES
# ======================================================================

# statistical significance of above vs below treeline differences is evaluated against the variable "AbiesTL"
# statistical analysis is shown below this plot

pdf("Figure_4_density_mortality_across_treeline.pdf", width=4.75, height=4.25)

    par(mfrow = c(2, 2))
    par(oma = c(1.5,1.5, 1.5, 0))   # make room  for the overall x and y axis titles
    par(mar = c(0, 2, 1.5, 1))      # make the plots be closer together
    par(ps = 8.5, cex = 1, cex.main = 1)
    
    head(rcampdata)
    
    
    bargraph.CI(x.factor = sizeclass, response = freq_live, group = AbiesTL1, data = rcampdata, na.rm=TRUE, 
                cex = 1.2, xlab = "", ylab="", main="", xaxt='n',
                cex.lab = 0.925, ylim=c(0,40), xlim=c(0.5,9),
                legend=TRUE, x.leg = 2.5, y.leg = 43, cex.leg=0.85, ncol=1, cex.axis=1,   
                col = c("gray80", "gray20"), border=NA, err.col = "gray50", err.width = 0.035)
    text(x = 1, y = 39, labels = "(a)", font=2, cex=1)
    text(x = 2, y = 23.5, labels = "*", font=2, cex=1.5)
    text(x = 5, y = 23.5, labels = "#", font=2, cex=0.85)
    text(x = 8, y = 29, labels = "#", font=2, cex=0.85)
    
    bargraph.CI(x.factor = sizeclass, response = freq_live, group = AbiesTL, data = abiesdata, na.rm=TRUE, 
                cex = 1.2, xlab = "", ylab = "", main="", cex.lab = 0.925, ylim=c(0,30), xlim=c(0.5,9), xaxt='n',
                legend=FALSE, x.leg = 0.775, y.leg = 15, cex.leg=0.85, ncol=1, cex.axis=1,   
                col = c("gray80", "gray20"), border=NA, err.col = "gray50", err.width = 0.035)
    text(x = 1, y = 29, labels = "(c)", font=2, cex=1)
    text(x = 2, y = 26, labels = "***", font=2, cex=1.5)
    text(x = 5, y = 4.5, labels = "**", font=2, cex=1.5)
    text(x = 8, y = 7.5, labels = "***", font=2, cex=1.5)
    
    bargraph.CI(x.factor = sizeclass, response = mortality, group = AbiesTL, data = rcampdata, na.rm=TRUE, 
                cex = 1.0, xlab = "", ylab = "", main="", cex.lab = 0.925, ylim=c(0,0.6), xlim=c(0.5,9), xaxt='n',
                legend=FALSE, x.leg = 0.775, y.leg = 0.5, cex.leg=0.85, ncol=1, cex.axis=1,   
                col = c("gray80", "gray20"), border=NA, err.col = "gray50", err.width = 0.035)
    text(x = 1, y = 0.58, labels = "(b)", font=2, cex=1)
    text(x = 2, y = 0.32, labels = "NS", font=2, cex=0.75)
    text(x = 5, y = 0.5, labels = "***", font=2, cex=1.5)
    text(x = 8, y = 0.25, labels = "*", font=2, cex=1.5)
    axis(1, tck=FALSE, line = -0.75, lwd = 0, at=c(2,5,8), cex.axis=1, labels=c("0-1m", "1-2m", ">2m"))
    
    bargraph.CI(x.factor = sizeclass, response = mortality, group = AbiesTL, data = abiesdata, na.rm=TRUE, 
                cex = 1.0, xlab = "", ylab = "", main="", cex.lab = 0.925, ylim=c(0,0.5), xlim=c(0.5,9),  xaxt='n',
                legend=FALSE, x.leg = 0.775, y.leg = 0.5, cex.leg=0.85, ncol=1, cex.axis=1,   
                col = c("gray80", "gray20"), border=NA, err.col = "gray50", err.width = 0.035)
    text(x = 1, y = 0.475, labels = "(d)", font=2, cex=1)
    text(x = 2, y = 0.195, labels = "NS", font=2, cex=0.75)
    text(x = 5, y = 0.375, labels = "NS", font=2, cex=0.75)
    text(x = 8, y = 0.175, labels = "NS", font=2, cex=0.75)
    axis(1, tck=FALSE, line = -0.75, lwd = 0, at=c(2,5,8), cex.axis=1, labels=c("0-1m", "1-2m", ">2m"))
    
    
    mtext(side=2, line=0.25, cex=1, adj=0.85, expression("Density (number per 100" ~ m^{2} ~ ")"), outer=TRUE) 
    mtext(side=2, line=0.35, cex=1, adj=0.16, "Mortality", outer=TRUE) 
    mtext(side=3, line=0, cex=1, adj=0.1, "Rhododendron campanulatum", outer=TRUE, font=4) 
    mtext(side=3, line=0, cex=1, adj=0.85, "Abies spectabilis", outer=TRUE, font=4) 


dev.off()




# ANOVA to determine treeline effect on each of the sizeclasses
# ==================================================================

# Rhododendron density
for (i in c("0-1m", "1-2m", "2m+")) {
    print(cat(paste("\\n----------- size class:", i, "-----------")))
    rm(sub)
    sub <- subset(rcampdata, sizeclass==i)
    print(summary(sub$sizeclass))
    rm(anova)
    anova <- aov(freq_live ~ AbiesTL*site*distance*transect, data=sub)
    print(summary(anova))
}

# Abies density
for (i in c("0-1m", "1-2m", "2m+")) {
    print(cat(paste("\\n----------- size class:", i, "-----------")))
    rm(sub)
    sub <- subset(abiesdata, sizeclass==i)
    print(summary(sub$sizeclass))
    rm(anova)
    anova <- aov(freq_live ~ AbiesTL*site*distance*transect, data=sub)
    print(summary(anova))
}

# Rhododendron mortality
for (i in c("0-1m", "1-2m", "2m+")) {
    print(cat(paste("\\n----------- size class:", i, "-----------")))
    rm(sub)
    sub <- subset(rcampdata, sizeclass==i)
    print(summary(sub$sizeclass))
    rm(anova)
    anova <- aov(mortality ~ AbiesTL*site*distance*transect, data=sub)
    print(summary(anova))
}

# Abies mortality
for (i in c("0-1m", "1-2m", "2m+")) {
    print(cat(paste("\\n----------- size class:", i, "-----------")))
    rm(sub)
    sub <- subset(abiesdata, sizeclass==i)
    print(summary(sub$sizeclass))
    rm(anova)
    anova <- aov(mortality ~ AbiesTL*site*distance*transect, data=sub)
    print(summary(anova))
}



# -------------------

# Does treeline interact with sizeclass for determining density?
# yes, it does as shown below
head(rcampdata)
intdata <- subset(rcampdata, sizeclass != "1-2m")
head(intdata)
rm(anova)
anova <- aov(freq_live ~ AbiesTL*sizeclass, data=intdata)
print(summary(anova))









