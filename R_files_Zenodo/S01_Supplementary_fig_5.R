### R SCRIPT TO CREATE SI APPENDIX FIGURE S1###

#read in coefficients from RMark CMR analysis#
regimes <- read.delim(file="CMR_coeffs_bodysize.txt", sep="\\t")

regimes$Phi.coef <- -1*regimes$Phi.coef #convert Phi from survival to extinction 
regimes$p.coef <- -1*regimes$p.coef ##convert p from survival to extinction 

taxa <- c("Trilobita", "Rhynchonellata", "Strophomenata", "Crinoidea", "Ostracoda", "Cephalopoda", "Echinoidea", "Bony fish",
          "Bivalvia", "Gastropoda")

#add numerical values to taxon field for plotting#
regimes$tax.num[regimes$class=="Trilobita"] <- 1
regimes$tax.num[regimes$class=="Rhynchonellata"] <- 2
regimes$tax.num[regimes$class=="Strophomenata"] <- 3
regimes$tax.num[regimes$class=="Ostracoda"] <- 4
regimes$tax.num[regimes$class=="Crinoidea"] <- 5
regimes$tax.num[regimes$class=="Cephalopoda"] <- 6
regimes$tax.num[regimes$class=="Echinoidea"] <- 7
regimes$tax.num[regimes$class=="bony fish"] <- 8
regimes$tax.num[regimes$class=="Bivalvia"] <- 9
regimes$tax.num[regimes$class=="Gastropoda"] <- 10 

#add colors and symbols for eras for plotting#
regimes$color <- "black"

regimes$tax.num[regimes$regime=="mass.extinction/recovery"] <- NA

pdf("Supplementary Sampling Figure.pdf", height=3.62, width=4.5)
par(oma=c(7,2.5,0,1), mar=c(1,2,1,2))

y.top <- 1
y.bot <- -1

plot(regimes$p.coef~regimes$tax.num, pch=16,
     xaxt="n", xlab="", ylab="", las=1,
     ylim=c(y.bot,y.top), col="white",
xlim=c(1,10))
abline(h=0, col="gray")
points(regimes$p.coef~regimes$tax.num,  pch=16,
     xaxt="n", xlab="", ylab="", las=1, col=regimes$color)
arrows(regimes$tax.num, regimes$p.coef+1.96*regimes$p.se,
      regimes$tax.num, regimes$p.coef-1.96*regimes$p.se, 
       angle=90, length=0.001, col=regimes$color, lwd=1.5)
axis(side=1, at=c(1:10), labels=taxa, las=2)
text(x=5, y=1, pos=1, labels="Sampling", cex=1, font=2)
mtext(side=2, "coefficient", line=2.5)
dev.off()


