## Load packages
library (phytools)
library(geiger)

##Set the working directory

## Read tree
tree <- read.tree ("C_tree.nwk")

## Plot tree
plot(tree)

## Read the CSV file
X <- read.csv("Discrete_13.csv", row.names=1, header=TRUE)

## Define the characters
slender <-setNames(X$slender,rownames(X)) 
appear <-setNames(X$appear,rownames(X))
rachis <-setNames(X$rachis,rownames(X))
bract <-setNames(X$bract,rownames(X))
density <-setNames(X$density,rownames(X))
cincinnus <-setNames(X$cincinnus,rownames(X)) 
flowers <-setNames(X$flowers,rownames(X))
calyx <-setNames(X$calyx,rownames(X))
bending <-setNames(X$bending,rownames(X))
claw <-setNames(X$claw,rownames(X))
labellum <-setNames(X$labellum,rownames(X))
blotch <-setNames(X$blotch,rownames(X))
stigma <-setNames(X$stigma,rownames(X))

##29. slender vs appear
##30. slender vs rachis
##31. slender vs bract
##32. slender vs density
##33. slender vs cincinnus
##34. slender vs flowers
##35. slender vs claw
##36. slender vs labellum
##37. slender vs blotch
##38. slender vs stigma
##39. slender vs bending
##40. slender vs calyx
##41. appear vs rachis
##42. appear vs bract
##43. appear vs density
##44. appear vs cincinnus
##45. appear vs flowers
##46. appear vs claw
##47. appear vs labellum
##48. appear vs blotch
##49. appear vs stigma
##50. appear vs bending
##51. appear vs calyx
##52. rachis vs bract
##53. rachis vs density 
##54. rachis vs cincinnus
##55. rachis vs flowers
##56. rachis vs claw
##57. rachis vs labellum
##58. rachis vs blotch
##59. rachis vs stigma
##60. rachis vs bending
##61. rachis vs calyx
##62. calyx vs bract
##63. calyx vs density
##64. calyx vs cincinnus
##65. calyx vs flowers
##66. calyx vs claw
##67. calyx vs labellum
##68. calyx vs blotch
##69. calyx vs stigma
##70. calyx vs bending
##71. bending vs bract
##72. bending vs density
##73. bending vs cincinnus
##74. bending vs flowers
##75. bending vs claw
##76. bending vs labellum
##77. bending vs blotch
##78. bending vs stigma


##29. slender vs appear
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_appear <-fitPagel(tree,slender,appear, model = "ARD", method = "fitDiscrete")
fit.xy_slender_appear
plot(fit.xy_slender_appear)
plot(fit.xy_slender_appear, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_appear_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_appear_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_appear)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_appear_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_appear, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##30. slender vs rachis
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_rachis <-fitPagel(tree,slender,rachis, model = "ARD", method = "fitDiscrete")
fit.xy_slender_rachis
plot(fit.xy_slender_rachis)
plot(fit.xy_slender_rachis, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_rachis_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_rachis_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_rachis)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_rachis_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_rachis, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##31. slender vs bract
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_bract <-fitPagel(tree,slender,bract, model = "ARD", method = "fitDiscrete")
fit.xy_slender_bract
plot(fit.xy_slender_bract)
plot(fit.xy_slender_bract, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_bract_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_bract_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_bract)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_bract_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_bract, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##32. slender vs density
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_density <-fitPagel(tree,slender,density, model = "ARD", method = "fitDiscrete")
fit.xy_slender_density
plot(fit.xy_slender_density)
plot(fit.xy_slender_density, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_density_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_density_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_density)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_density_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_density, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##33. slender vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_cincinnus <-fitPagel(tree,slender,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_slender_cincinnus
plot(fit.xy_slender_cincinnus)
plot(fit.xy_slender_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##34. slender vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_flowers <-fitPagel(tree,slender,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_slender_flowers
plot(fit.xy_slender_flowers)
plot(fit.xy_slender_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##35. slender vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_claw <-fitPagel(tree,slender,claw, model = "ARD", method = "fitDiscrete")
fit.xy_slender_claw
plot(fit.xy_slender_claw)
plot(fit.xy_slender_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##36. slender vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_labellum <-fitPagel(tree,slender,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_slender_labellum
plot(fit.xy_slender_labellum)
plot(fit.xy_slender_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##37. slender vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_blotch <-fitPagel(tree,slender,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_slender_blotch
plot(fit.xy_slender_blotch)
plot(fit.xy_slender_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##38. slender vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_stigma <-fitPagel(tree,slender,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_slender_stigma
plot(fit.xy_slender_stigma)
plot(fit.xy_slender_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##39. slender vs bending
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_bending <-fitPagel(tree,slender,bending, model = "ARD", method = "fitDiscrete")
fit.xy_slender_bending
plot(fit.xy_slender_bending)
plot(fit.xy_slender_bending, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_bending_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_bending_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_bending)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_bending_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_bending, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##40. slender vs calyx
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_slender_calyx <-fitPagel(tree,slender,calyx, model = "ARD", method = "fitDiscrete")
fit.xy_slender_calyx
plot(fit.xy_slender_calyx)
plot(fit.xy_slender_calyx, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("slender_calyx_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(slender[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_calyx_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_calyx)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("slender_calyx_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_slender_calyx, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##41. appear vs rachis
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_rachis <-fitPagel(tree,appear,rachis, model = "ARD", method = "fitDiscrete")
fit.xy_appear_rachis
plot(fit.xy_appear_rachis)
plot(fit.xy_appear_rachis, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_rachis_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_rachis_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_rachis)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_rachis_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_rachis, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##42. appear vs bract
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_bract <-fitPagel(tree,appear,bract, model = "ARD", method = "fitDiscrete")
fit.xy_appear_bract
plot(fit.xy_appear_bract)
plot(fit.xy_appear_bract, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_bract_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_bract_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_bract)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_bract_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_bract, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##43. appear vs density
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_density <-fitPagel(tree,appear,density, model = "ARD", method = "fitDiscrete")
fit.xy_appear_density
plot(fit.xy_appear_density)
plot(fit.xy_appear_density, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_density_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_density_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_density)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_density_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_density, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()

##44. appear vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_cincinnus <-fitPagel(tree,appear,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_appear_cincinnus
plot(fit.xy_appear_cincinnus)
plot(fit.xy_appear_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##45. appear vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_flowers <-fitPagel(tree,appear,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_appear_flowers
plot(fit.xy_appear_flowers)
plot(fit.xy_appear_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##46. appear vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_claw <-fitPagel(tree,appear,claw, model = "ARD", method = "fitDiscrete")
fit.xy_appear_claw
plot(fit.xy_appear_claw)
plot(fit.xy_appear_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##47. appear vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_labellum <-fitPagel(tree,appear,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_appear_labellum
plot(fit.xy_appear_labellum)
plot(fit.xy_appear_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##48. appear vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_blotch <-fitPagel(tree,appear,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_appear_blotch
plot(fit.xy_appear_blotch)
plot(fit.xy_appear_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##49. appear vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_stigma <-fitPagel(tree,appear,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_appear_stigma
plot(fit.xy_appear_stigma)
plot(fit.xy_appear_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##50. appear vs bending
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_bending <-fitPagel(tree,appear,bending, model = "ARD", method = "fitDiscrete")
fit.xy_appear_bending
plot(fit.xy_appear_bending)
plot(fit.xy_appear_bending, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_bending_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_bending_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_bending)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_bending_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_bending, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##51. appear vs calyx
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_appear_calyx <-fitPagel(tree,appear,calyx, model = "ARD", method = "fitDiscrete")
fit.xy_appear_calyx
plot(fit.xy_appear_calyx)
plot(fit.xy_appear_calyx, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("appear_calyx_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(appear[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_calyx_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_calyx)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("appear_calyx_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_appear_calyx, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##52. rachis vs bract
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_bract <-fitPagel(tree,rachis,bract, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_bract
plot(fit.xy_rachis_bract)
plot(fit.xy_rachis_bract, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_bract_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_bract_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_bract)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_bract_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_bract, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##53. rachis vs density 
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_density <-fitPagel(tree,rachis,density, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_density
plot(fit.xy_rachis_density)
plot(fit.xy_rachis_density, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_density_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_density_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_density)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_density_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_density, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##54. rachis vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_cincinnus <-fitPagel(tree,rachis,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_cincinnus
plot(fit.xy_rachis_cincinnus)
plot(fit.xy_rachis_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##55. rachis vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_flowers <-fitPagel(tree,rachis,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_flowers
plot(fit.xy_rachis_flowers)
plot(fit.xy_rachis_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##56. rachis vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_claw <-fitPagel(tree,rachis,claw, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_claw
plot(fit.xy_rachis_claw)
plot(fit.xy_rachis_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##57. rachis vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_labellum <-fitPagel(tree,rachis,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_labellum
plot(fit.xy_rachis_labellum)
plot(fit.xy_rachis_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##58. rachis vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_blotch <-fitPagel(tree,rachis,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_blotch
plot(fit.xy_rachis_blotch)
plot(fit.xy_rachis_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##59. rachis vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_stigma <-fitPagel(tree,rachis,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_stigma
plot(fit.xy_rachis_stigma)
plot(fit.xy_rachis_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##60. rachis vs bending
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_bending <-fitPagel(tree,rachis,bending, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_bending
plot(fit.xy_rachis_bending)
plot(fit.xy_rachis_bending, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_bending_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_bending_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_bending)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_bending_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_bending, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##61. rachis vs calyx
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_rachis_calyx <-fitPagel(tree,rachis,calyx, model = "ARD", method = "fitDiscrete")
fit.xy_rachis_calyx
plot(fit.xy_rachis_calyx)
plot(fit.xy_rachis_calyx, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("rachis_calyx_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(rachis[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_calyx_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_calyx)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("rachis_calyx_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_rachis_calyx, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##62. calyx vs bract
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_bract <-fitPagel(tree,calyx,bract, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_bract
plot(fit.xy_calyx_bract)
plot(fit.xy_calyx_bract, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_bract_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_bract_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_bract)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_bract_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_bract, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##63. calyx vs density
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_density <-fitPagel(tree,calyx,density, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_density
plot(fit.xy_calyx_density)
plot(fit.xy_calyx_density, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_density_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_density_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_density)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_density_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_density, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##64. calyx vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_cincinnus <-fitPagel(tree,calyx,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_cincinnus
plot(fit.xy_calyx_cincinnus)
plot(fit.xy_calyx_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##65. calyx vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_flowers <-fitPagel(tree,calyx,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_flowers
plot(fit.xy_calyx_flowers)
plot(fit.xy_calyx_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##66. calyx vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_claw <-fitPagel(tree,calyx,claw, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_claw
plot(fit.xy_calyx_claw)
plot(fit.xy_calyx_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##67. calyx vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_labellum <-fitPagel(tree,calyx,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_labellum
plot(fit.xy_calyx_labellum)
plot(fit.xy_calyx_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##68. calyx vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_blotch <-fitPagel(tree,calyx,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_blotch
plot(fit.xy_calyx_blotch)
plot(fit.xy_calyx_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##69. calyx vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_stigma <-fitPagel(tree,calyx,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_stigma
plot(fit.xy_calyx_stigma)
plot(fit.xy_calyx_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##70. calyx vs bending
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_calyx_bending <-fitPagel(tree,calyx,bending, model = "ARD", method = "fitDiscrete")
fit.xy_calyx_bending
plot(fit.xy_calyx_bending)
plot(fit.xy_calyx_bending, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("calyx_bending_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(calyx[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_bending_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_bending)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("calyx_bending_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_calyx_bending, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##71. bending vs bract
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_bract <-fitPagel(tree,bending,bract, model = "ARD", method = "fitDiscrete")
fit.xy_bending_bract
plot(fit.xy_bending_bract)
plot(fit.xy_bending_bract, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_bract_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_bract_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_bract)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_bract_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_bract, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##72. bending vs density
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_density <-fitPagel(tree,bending,density, model = "ARD", method = "fitDiscrete")
fit.xy_bending_density
plot(fit.xy_bending_density)
plot(fit.xy_bending_density, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_density_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_density_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_density)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_density_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_density, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##73. bending vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_cincinnus <-fitPagel(tree,bending,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_bending_cincinnus
plot(fit.xy_bending_cincinnus)
plot(fit.xy_bending_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##74. bending vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_flowers <-fitPagel(tree,bending,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_bending_flowers
plot(fit.xy_bending_flowers)
plot(fit.xy_bending_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##75. bending vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_claw <-fitPagel(tree,bending,claw, model = "ARD", method = "fitDiscrete")
fit.xy_bending_claw
plot(fit.xy_bending_claw)
plot(fit.xy_bending_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##76. bending vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_labellum <-fitPagel(tree,bending,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_bending_labellum
plot(fit.xy_bending_labellum)
plot(fit.xy_bending_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##77. bending vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_blotch <-fitPagel(tree,bending,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_bending_blotch
plot(fit.xy_bending_blotch)
plot(fit.xy_bending_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


##78. bending vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bending_stigma <-fitPagel(tree,bending,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_bending_stigma
plot(fit.xy_bending_stigma)
plot(fit.xy_bending_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bending_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bending[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bending_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bending_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()