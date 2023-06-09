###############################################################################
###############################################################################
# Publication Bias in the Social Sciences: Unlocking the File Drawer
# (Annie Franco, Neil Malhotra, Gabor Simonovits)
###############################################################################
###############################################################################
# REPLICATION CODE (SUPPLEMENTARY MATERIAL)
###############################################################################
# Fig. S1. Sensitivity of Pearson chi-squared test of independence in Table 3 to
# misclassification of TESS studies. 
###############################################################################

library(foreign)


# load dataset
pubbias <- read.dta("filedrawer.dta")
names(pubbias)

bounds <- as.matrix(table(pubbias$anyresults, pubbias$written))
rownames(bounds)<-c("Null", "Significant")
colnames(bounds)<-c("Unwritten", "Written")

chisq.test(bounds)

x <- seq(0,20,1)
y <- seq(0,70,1)

grid <- as.matrix(expand.grid(x, y))

chisq_pval <- function(q) {
  
  a<-q[1]
  b<-q[2]
  x<- cbind(c(a,-a), c(-b,b))
  chisq.test(bounds-x)$p.value
  
}

pvalues <- apply(grid, 1, FUN = chisq_pval)

sensitivity <- data.frame(unwritten = grid[,1], 
                          written = grid[,2], 
                          pvalue = pvalues)

sensitivity$sig <- as.numeric(sensitivity$pvalue<0.05)

pdf(file = "bounds.pdf", width=6.8, height=4.2)
par(mfrow=c(1,1), bg="white", mgp=c(2,.5,0), mar=c(4,4,2,2))
with(sensitivity[sensitivity$pvalue<0.05,], 
     plot(written, unwritten, 
          xlim=c(0,70), 
          cex=0.7,
          ylim=c(0,20),
          pch=16, 
          col="black", 
          axes=FALSE,
          cex.axis=0.9,
          cex.lab=1.1,
          ylab="",
          xlab="Written studies recoded as null (out of 159)"))
title(ylab="Unwritten studies recoded\\nas significant (out of 31)", 
      line=1.5, cex.lab=1.1)
axis(side=1, tick=T, at=seq(0,100,10), cex.axis=0.8)
axis(side=2, tick=T, at=seq(0,100,10), cex.axis=0.8)
with(sensitivity[sensitivity$pvalue>0.05,], 
     points(written, unwritten, col="black", pch=1, cex=0.7))
box()
legend("topright", 
       c("P-value < 0.05", "P-value > 0.05"),
       pch=c(16,1), cex=0.95)
dev.off()
