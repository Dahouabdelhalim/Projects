###R SCRIPT TO CREATE FIGURE 2###

#read in coefficients from RMark CMR analysis#
regimes <- read.delim(file="CMR_coeffs_bodysize_final.txt", sep="\\t")

regimes$Phi.coef <- -1*regimes$Phi.coef #convert Phi from survival to regimes 
regimes$p.coef <- -1*regimes$p.coef ##convert p from survival to regimes 

taxa <- c("Trilobita", "Rhynchonellata", "Strophomenata", "Ostracoda", "Crinoidea",  "Cephalopoda", "Echinoidea", "bony fish", "Bivalvia", "Gastropoda")

regimes$tax.num <- NA
#add numerical values to class field for plotting#
regimes$tax.num[regimes$class=="Trilobita"] <- 0.5
regimes$tax.num[regimes$class=="Rhynchonellata"] <- 1.5
regimes$tax.num[regimes$class=="Strophomenata"] <- 2.5
regimes$tax.num[regimes$class=="Ostracoda"] <- 3.5
regimes$tax.num[regimes$class=="Crinoidea"] <- 4.5
regimes$tax.num[regimes$class=="Cephalopoda"] <- 5.5
regimes$tax.num[regimes$class=="Echinoidea"] <- 6.5
regimes$tax.num[regimes$class=="bony fish"] <- 7.5
regimes$tax.num[regimes$class=="Bivalvia"] <- 8.5
regimes$tax.num[regimes$class=="Gastropoda"] <- 9.5 
 
#add colors for plotting#
regimes$color1 <- "#d95f02"
regimes$color2 <- "#7570b3"

#add symbols for plotting
regimes$sym[regimes$regime=="background"] <- 15
regimes$sym[regimes$regime=="mass.extinction/recovery"] <- 16

regimes$sym[regimes$regime=="background" & regimes$class=="Strophomenata"] <- 17
regimes$sym[regimes$regime=="background" & regimes$class=="Echinoidea"] <- 17
regimes$sym[regimes$regime=="background" & regimes$class=="bony fish"] <- 17
regimes$sym[regimes$regime=="background" & regimes$class=="Trilobita"] <- 17
regimes$sym[regimes$regime=="background" & regimes$class=="Ostracoda"] <- 17
regimes$sym[regimes$regime=="background" & regimes$class=="Bivalvia"] <- 17

#add space for plotting points
regimes$tax.num[regimes$regime=="background"] <- regimes$tax.num[regimes$regime=="background"]-0.1
regimes$tax.num[regimes$regime=="mass.extinction/recovery"] <- regimes$tax.num[regimes$regime=="mass.extinction/recovery"]+0.1



pdf("Figure 2_ProcB.pdf",  h=3.63, w=6.82)
par(mfrow=c(1,2), oma=c(7.5,2.5,0,1), mar=c(1,2,1,2))

y.top <- 3.5
y.bot <- -3.5

plot(regimes$Phi.coef~regimes$tax.num, pch=regimes$sym,
     xaxt="n", xlab="", ylab="", las=1,
     ylim=c(y.bot,y.top), col="white",
xlim=c(0,10))
abline(h=0, col="gray")
points(regimes$Phi.coef~regimes$tax.num,  pch=regimes$sym,
     xaxt="n", xlab="", ylab="", las=1, col=regimes$color1)
arrows(regimes$tax.num, regimes$Phi.coef+1.96*regimes$Phi.se,
      regimes$tax.num, regimes$Phi.coef-1.96*regimes$Phi.se, 
       angle=90, length=0.001, col=regimes$color1, lwd=1.5)
axis(side=1, at=c(0.5:9.5), labels=taxa, las=2, cex.axis=0.85)
text(x=5, y=3.2, pos=1, labels="Extinction", cex=1, font=1)
text(0.1,2.9 , "A")
mtext(side=2, "coefficient", line=2.5)
legend("bottomright", pch=c(17, 15, 16), col="#d95f02",
       legend=c("single extinction regime","background extinction", "mass extinction"),
bty="n", cex=0.6)

plot(regimes$Gamma.coef~regimes$tax.num, pch=regimes$sym,
     xaxt="n", xlab="", ylab="", las=1,
     ylim=c(y.bot,y.top), col="white",
xlim=c(0,10))
abline(h=0, col="gray")
points(regimes$Gamma.coef~regimes$tax.num, pch=regimes$sym,
     xaxt="n", xlab="", ylab="", las=1, col=regimes$color2)
arrows(regimes$tax.num, regimes$Gamma.coef+1.96*regimes$Gamma.se,
      regimes$tax.num, regimes$Gamma.coef-1.96*regimes$Gamma.se, 
       angle=90, length=0.001, col=regimes$color2, lwd=1.5)
axis(side=1, at=c(0.5:9.5), labels=taxa, las=2, cex.axis=0.85)
text(x=5, y=3.2, pos=1, labels="Origination", cex=1, font=1)
text(0.1,2.9 , "B")
mtext(side=2, "coefficient", line=2.5)
legend("bottomright", pch=c(17, 15, 16), col="#7570b3",
       legend=c("single origination regime", "background origination", "recovery origination"),
bty="n", cex=0.6)
dev.off()
