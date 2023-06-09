install.packages("qtl")
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

setwd("working directory")

library(qtl)
library(snow)
library(qtl2)

#read in the data
Costus <- read.cross("csvs", ".", genfile="F2gen.csv", phefile="F2phe.csv", genotypes = c("PP", "PS","SS","D","C"), alleles = c("P", "S"))

#with PRF and POLMX log transformed
CostusT <- read.cross("csvs", ".", genfile="F2gen.csv", phefile="F2phe_trans.csv", genotypes = c("PP", "PS","SS","D","C"), alleles = c("P", "S"))

Costus <- calc.genoprob(Costus, step=2)
CostusT <- calc.genoprob(CostusT, step=2)

#scantwo to do permutations for likelihood penalties

operm2.hk <- scantwo(Costus, method = "hk", pheno.col = 2:16, n.cluster= 3, n.perm=1000, verbose = TRUE)
summary(operm2.hk)
summary(out2.hk, perms=operm2.hk, pvalues=TRUE,
        alphas=c(0.05, 0.05, 0.05, 0.05, 0.05), pheno.col = 2)

operm2.PRFPOLMX<- scantwo(CostusT, method = "hk", pheno.col = c("PRF","POLMX"), n.cluster= 3, n.perm=1000, verbose = TRUE)
summary(operm2.PRFPOLMX)

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
Costus_penalties <- calc.penalties(operm2.hk, alpha = 0.05)
Costus_penaltiesT <- calc.penalties(operm2.PRFPOLMX, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(Costus_penalties, file = "Costus_penalties.csv")
write.csv(Costus_penaltiesT, file = "Costus_penaltiesT.csv")

Costus_penalties <- read.csv(file = "Costus_penalties.csv", header = TRUE, row.names = 1)
Costus_penaltiesT <- read.csv(file = "Costus_penaltiesT.csv", header = TRUE, row.names = 1 )

#APH
APHqtl <- stepwiseqtl(Costus, pheno.col = "APH", method = "hk", model = "normal", 
            penalties = as.vector(Costus_penalties["APH",1:3]), max.qtl = 8, keeplodprofile = TRUE,
            verbose = TRUE, refine.locations = TRUE)
summary(APHqtl)
plot(APHqtl)
APH.fit <- fitqtl(Costus, pheno.col = "APH", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = APHqtl, formula = attr(APHqtl,"formula"), dropone = TRUE)
summary(APH.fit)
#get bayesint for location
bayesint(APHqtl, qtl.index = 1)
bayesint(APHqtl, qtl.index = 2)
bayesint(APHqtl, qtl.index = 3)

pdf(file="APH_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(APHqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="APH", frame.plot=FALSE)
dev.off()

#APW
APWqtl <- stepwiseqtl(Costus, pheno.col = "APW", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["APW",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(APWqtl)
plot(APWqtl)
APW.fit <- fitqtl(Costus, pheno.col = "APW", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = APWqtl, formula = attr(APWqtl,"formula"))
summary(APW.fit)
#get bayesint for location
bayesint(APWqtl, qtl.index = 1)

pdf(file="APW_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(APWqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="APW", frame.plot=FALSE)
dev.off()

#ANL
ANLqtl <- stepwiseqtl(Costus, pheno.col = "ANL", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["ANL",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(ANLqtl)
plot(ANLqtl)
ANL.fit <- fitqtl(Costus, pheno.col = "ANL", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = ANLqtl, formula = attr(ANLqtl,"formula"))
summary(ANL.fit)
#get bayesint for location
bayesint(ANLqtl, qtl.index = 1)
bayesint(ANLqtl, qtl.index = 2)
bayesint(ANLqtl, qtl.index = 3)

pdf(file="ANL_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(ANLqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="ANL", frame.plot=FALSE)
dev.off()

#plot overdominant QTL
pdf(file = "ANL2@106.25.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=2, pos=106.25)
plotPXG(Costus, marker = marker, pheno.col="ANL")
dev.off()

#ANW
ANWqtl <- stepwiseqtl(Costus, pheno.col = "ANW", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["ANW",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(ANWqtl)
plot(ANWqtl)
ANW.fit <- fitqtl(Costus, pheno.col = "ANW", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = ANWqtl, formula = attr(ANWqtl,"formula"))
summary(ANW.fit)
#get bayesint for location
bayesint(ANWqtl, qtl.index = 1)
bayesint(ANWqtl, qtl.index = 2)
bayesint(ANWqtl, qtl.index = 3)

pdf(file="ANW_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(ANWqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="ANW", frame.plot=FALSE)
dev.off()

#plot overdominant QTL
pdf(file = "ANW1@32.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=1, pos=32)
plotPXG(Costus, marker = marker, pheno.col="ANW")
dev.off()

#STIH
STIHqtl <- stepwiseqtl(Costus, pheno.col = "STIH", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["STIH",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(STIHqtl)
plot(STIHqtl)
STIH.fit <- fitqtl(Costus, pheno.col = "STIH", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = STIHqtl, formula = attr(STIHqtl,"formula"))
summary(STIH.fit)
#get bayesint for location
bayesint(STIHqtl, qtl.index = 1)
bayesint(STIHqtl, qtl.index = 2)
bayesint(STIHqtl, qtl.index = 3)

pdf(file="STIH_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(STIHqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="STIH", frame.plot=FALSE)
dev.off()
#consider plotting interactions with qtlcharts iplotScantwo

#STIW
STIWqtl <- stepwiseqtl(Costus, pheno.col = "STIW", method = "hk", model = "normal", 
                       penalties = as.vector(Costus_penalties["STIW",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                       verbose = TRUE, refine.locations = TRUE)
summary(STIWqtl)
plot(STIWqtl)
STIW.fit <- fitqtl(Costus, pheno.col = "STIW", method = "hk", model = "normal", get.ests = TRUE,
                   qtl = STIWqtl, formula = attr(STIWqtl,"formula"))
summary(STIW.fit)
#get bayesint for location
bayesint(STIWqtl, qtl.index = 1)
bayesint(STIWqtl, qtl.index = 2)
bayesint(STIWqtl, qtl.index = 3)
bayesint(STIWqtl, qtl.index = 4)
bayesint(STIWqtl, qtl.index = 5)
bayesint(STIWqtl, qtl.index = 6)
bayesint(STIWqtl, qtl.index = 7)
bayesint(STIWqtl, qtl.index = 8)

pdf(file="STIW_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(STIWqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="STIW", frame.plot=FALSE)
dev.off()

#plot overdominant QTL
pdf(file = "STIW2@52.97.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=2, pos=52.97)
plotPXG(Costus, marker = marker, pheno.col="STIW")
dev.off()

#not overdominant but antagonistic
pdf(file = "STIW3@42.4.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=3, pos=42.4)
plotPXG(Costus, marker = marker, pheno.col="STIW")
dev.off()

pdf(file = "STIW4@111.32.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=4, pos=111.32)
plotPXG(Costus, marker = marker, pheno.col="STIW")
dev.off()

#LABL
LABLqtl <- stepwiseqtl(Costus, pheno.col = "LABL", method = "hk", model = "normal", 
                       penalties = as.vector(Costus_penalties["LABL",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                       verbose = TRUE, refine.locations = TRUE)
summary(LABLqtl)
plot(LABLqtl)
LABL.fit <- fitqtl(Costus, pheno.col = "LABL", method = "hk", model = "normal", get.ests = TRUE,
                   qtl = LABLqtl, formula = attr(LABLqtl,"formula"))
summary(LABL.fit)
#get bayesint for location
bayesint(LABLqtl, qtl.index = 1)
bayesint(LABLqtl, qtl.index = 2)
bayesint(LABLqtl, qtl.index = 3)
bayesint(LABLqtl, qtl.index = 4)
bayesint(LABLqtl, qtl.index = 5)
bayesint(LABLqtl, qtl.index = 6)
bayesint(LABLqtl, qtl.index = 7)
bayesint(LABLqtl, qtl.index = 8)

pdf(file="LABL_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(LABLqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="LABL", frame.plot=FALSE)
dev.off()

#CLL 
CLLqtl <- stepwiseqtl(Costus, pheno.col = "CLL", method = "hk", model = "normal",
                      penalties = as.vector(Costus_penalties["CLL",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(CLLqtl)
plot(CLLqtl)
CLL.fit <- fitqtl(Costus, pheno.col = "CLL", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = CLLqtl, formula = attr(CLLqtl,"formula"))
summary(CLL.fit)
#get bayesint for location
bayesint(CLLqtl, qtl.index = 1)
bayesint(CLLqtl, qtl.index = 2)
bayesint(CLLqtl, qtl.index = 3)
bayesint(CLLqtl, qtl.index = 4)
bayesint(CLLqtl, qtl.index = 5)
bayesint(CLLqtl, qtl.index = 6)

pdf(file="CLL_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(CLLqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="CLL", frame.plot=FALSE)
dev.off()


#STYL
STYLqtl <- stepwiseqtl(Costus, pheno.col = "STYL", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["STYL",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(STYLqtl)
plot(STYLqtl)
STYL.fit <- fitqtl(Costus, pheno.col = "STYL", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = STYLqtl, formula = attr(STYLqtl,"formula"))
summary(STYL.fit)
#get bayesint for location
bayesint(STYLqtl, qtl.index = 1)
bayesint(STYLqtl, qtl.index = 2)
bayesint(STYLqtl, qtl.index = 3)
bayesint(STYLqtl, qtl.index = 4)
bayesint(STYLqtl, qtl.index = 5)
bayesint(STYLqtl, qtl.index = 6)
bayesint(STYLqtl, qtl.index = 7)

pdf(file="STYL_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(STYLqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="STYL", frame.plot=FALSE)
dev.off()

pdf(file = "STYL2@106.37.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=2, pos=106.37)
plotPXG(Costus, marker = marker, pheno.col="STYL")
dev.off()

#STAL
STALqtl <- stepwiseqtl(Costus, pheno.col = "STAL", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["STAL",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(STALqtl)
plot(STALqtl)
STAL.fit <- fitqtl(Costus, pheno.col = "STAL", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = STALqtl, formula = attr(STALqtl,"formula"))
summary(STAL.fit)
#get bayesint for location
bayesint(STALqtl, qtl.index = 1)
bayesint(STALqtl, qtl.index = 2)
bayesint(STALqtl, qtl.index = 3)
bayesint(STALqtl, qtl.index = 4)
bayesint(STALqtl, qtl.index = 5)
bayesint(STALqtl, qtl.index = 6)
bayesint(STALqtl, qtl.index = 7)

pdf(file="STAL_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(STALqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="STAL", frame.plot=FALSE)
dev.off()

pdf(file = "STAL2@82.27.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=2, pos=82.27)
plotPXG(Costus, marker = marker, pheno.col="STAL")
dev.off()

pdf(file = "STAL2@106.37.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=2, pos=106.37)
plotPXG(Costus, marker = marker, pheno.col="STAL")
dev.off()

pdf(file = "STAL3@92.11.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=3, pos=92.11)
plotPXG(Costus, marker = marker, pheno.col="STAL")
dev.off()

pdf(file = "STAL2@82.27x5@84.37.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=2, pos=82.27)
marker2 <- find.marker(Costus, chr=5, pos=84.37)
plotPXG(Costus, marker = c(marker, marker2), pheno.col="STAL")
dev.off()


#LLBS
LLBSqtl <- stepwiseqtl(Costus, pheno.col = "LLBS", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["LLBS",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(LLBSqtl)
plot(LLBSqtl)
LLBS.fit <- fitqtl(Costus, pheno.col = "LLBS", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = LLBSqtl, formula = attr(LLBSqtl,"formula"))
summary(LLBS.fit)
#get bayesint for location
bayesint(LLBSqtl, qtl.index = 1)
bayesint(LLBSqtl, qtl.index = 2)
bayesint(LLBSqtl, qtl.index = 3)
bayesint(LLBSqtl, qtl.index = 4)
bayesint(LLBSqtl, qtl.index = 5)
bayesint(LLBSqtl, qtl.index = 6)
bayesint(LLBSqtl, qtl.index = 7)
bayesint(LLBSqtl, qtl.index = 8)

pdf(file="LLBS_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(LLBSqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="LLBS", frame.plot=FALSE)
dev.off()

pdf(file = "LLBS1@33.7.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=1, pos=33.7)
plotPXG(Costus, marker = marker, pheno.col="LLBS")
dev.off()

#LANG
LANGqtl <- stepwiseqtl(Costus, pheno.col = "LANG", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["LANG",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(LANGqtl)
plot(LANGqtl)
LANG.fit <- fitqtl(Costus, pheno.col = "LANG", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = LANGqtl, formula = attr(LANGqtl,"formula"))
summary(LANG.fit)
#get bayesint for location
bayesint(LANGqtl, qtl.index = 1)
bayesint(LANGqtl, qtl.index = 2)
bayesint(LANGqtl, qtl.index = 3)
bayesint(LANGqtl, qtl.index = 4)

pdf(file="LANG_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(LANGqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="LANG", frame.plot=FALSE)
dev.off()

#PRF
PRFqtl <- stepwiseqtl(CostusT, pheno.col = "PRF", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penaltiesT["PRF",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(PRFqtl)
plot(PRFqtl)
#fitting the model to the transformed values first to get LOD scores, then to the raw values to get effects
PRF.fitT <- fitqtl(CostusT, pheno.col = "PRF", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = PRFqtl, formula = attr(PRFqtl,"formula"))
summary(PRF.fitT)

PRF.fit <- fitqtl(Costus, pheno.col = "PRF", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = PRFqtl, formula = attr(PRFqtl,"formula"))
summary(PRF.fit)
#get bayesint for location
bayesint(PRFqtl, qtl.index = 1)

pdf(file="PRF_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(PRFqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="PRF", frame.plot=FALSE)
dev.off()

#CST
CSTqtl <- stepwiseqtl(Costus, pheno.col = "CST", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties["CST",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(CSTqtl)
plot(CSTqtl)
CST.fit <- fitqtl(Costus, pheno.col = "CST", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = CSTqtl, formula = attr(CSTqtl,"formula"))
summary(CST.fit)

#get bayesint for location
bayesint(CSTqtl, qtl.index = 1)

pdf(file="CST_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(CSTqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="CST", frame.plot=FALSE)
dev.off()

#POLMX
POLMXqtl <- stepwiseqtl(CostusT, pheno.col = "POLMX", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penaltiesT["POLMX",1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(POLMXqtl)
plot(POLMXqtl)
#fit model with transformed data to get LODs, then fit to raw phenotypes to get effects
POLMX.fitT <- fitqtl(CostusT, pheno.col = "POLMX", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = POLMXqtl, formula = attr(POLMXqtl,"formula"))
summary(POLMX.fitT)

POLMX.fit <- fitqtl(Costus, pheno.col = "POLMX", method = "hk", model = "normal", get.ests = TRUE,
                    qtl = POLMXqtl, formula = attr(POLMXqtl,"formula"))
summary(POLMX.fit)

#plot underdominant QTL
pdf(file = "POLMX7@80.pdf", width=3.14, height = 3.14, pointsize=9)
marker <- find.marker(Costus, chr=7, pos=80)
plotPXG(Costus, marker = marker, pheno.col="POLMX")
dev.off()

#get bayesint for location
bayesint(POLMXqtl, qtl.index = 1)

pdf(file="POLMX_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
plotLodProfile(POLMXqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="POLMX", frame.plot=FALSE)
dev.off()

##########################################################################################
#use the plot_peaks function in RQTL2 to plot all the QTL

setwd("working directory")


#load the control file, which loads all the data
costusPS <- read_cross2("control.yaml")

#isolate the linkage map
map<- costusPS$gmap

#import summary file of QTL found above
df <- read.csv("QTL_summary.csv")

pdf(file="QTL_summary.pdf", width=7.08, height= 7.08 )
plot_peaks(df, map,
           gap = NULL, tick_height = 0.3, lod_labels = FALSE,
           bgcolor = "gray90", altbgcolor = "gray75")

dev.off()


