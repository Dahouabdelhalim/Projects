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

##1. bract vs density
##2. bract vs cincinnus
##3. bract vs flowers
##4. bract vs claw
##5. bract vs labellum
##6. bract vs blotch
##7. bract vs stigma
##8. density vs cincinnus
##9. density vs flowers
##10. density vs claw
##11. density vs labellum
##12. density vs blotch
##13. density vs stigma
##14. cincinnus vs flowers
##15. cincinnus vs claw
##16. cincinnus vs labellum
##17. cincinnus vs blotch
##18. cincinnus vs stigma
##19. flowers vs claw
##20. flowers vs labellum
##21. flowers vs blotch
##22. flowers vs stigma
##23. claw vs labellum
##24. claw vs blotch
##25. claw vs stigma 
##26. labellum vs blotch
##27. labellum vs stigma
##28. blotch vs stigma

## (1) bract vs density
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_density <-fitPagel(tree,bract,density, model = "ARD", method = "fitDiscrete")
fit.xy_bract_density
plot(fit.xy_bract_density)
plot(fit.xy_bract_density, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_density_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_density_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_density)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_density_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_density, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (2) bract vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_cincinnus <-fitPagel(tree,bract,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_bract_cincinnus
plot(fit.xy_bract_cincinnus)
plot(fit.xy_bract_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (3) bract vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_flowers <-fitPagel(tree,bract,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_bract_flowers
plot(fit.xy_bract_flowers)
plot(fit.xy_bract_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (4) bract vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_claw <-fitPagel(tree,bract,claw, model = "ARD", method = "fitDiscrete")
fit.xy_bract_claw
plot(fit.xy_bract_claw)
plot(fit.xy_bract_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (5) bract vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_labellum <-fitPagel(tree,bract,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_bract_labellum
plot(fit.xy_bract_labellum)
plot(fit.xy_bract_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (6) bract vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_blotch <-fitPagel(tree,bract,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_bract_blotch
plot(fit.xy_bract_blotch)
plot(fit.xy_bract_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (7) bract vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_bract_stigma <-fitPagel(tree,bract,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_bract_stigma
plot(fit.xy_bract_stigma)
plot(fit.xy_bract_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("bract_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(bract[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("bract_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("bract_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_bract_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (8) density vs cincinnus
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_density_cincinnus <-fitPagel(tree,density,cincinnus, model = "ARD", method = "fitDiscrete")
fit.xy_density_cincinnus
plot(fit.xy_density_cincinnus)
plot(fit.xy_density_cincinnus, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("density_cincinnus_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("density_cincinnus_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_cincinnus)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("density_cincinnus_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_cincinnus, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (9) density vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_density_flowers <-fitPagel(tree,density,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_density_flowers
plot(fit.xy_density_flowers)
plot(fit.xy_density_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("density_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("density_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("density_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (10) density vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_density_claw <-fitPagel(tree,density,claw, model = "ARD", method = "fitDiscrete")
fit.xy_density_claw
plot(fit.xy_density_claw)
plot(fit.xy_density_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("density_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("density_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("density_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (11) density vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_density_labellum <-fitPagel(tree,density,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_density_labellum
plot(fit.xy_density_labellum)
plot(fit.xy_density_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("density_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("density_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("density_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (12) density vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_density_blotch <-fitPagel(tree,density,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_density_blotch
plot(fit.xy_density_blotch)
plot(fit.xy_density_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("density_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("density_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("density_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (13) density vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_density_stigma <-fitPagel(tree,density,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_density_stigma
plot(fit.xy_density_stigma)
plot(fit.xy_density_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("density_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(density[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("density_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("density_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_density_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (14) cincinnus vs flowers
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_cincinnus_flowers <-fitPagel(tree,cincinnus,flowers, model = "ARD", method = "fitDiscrete")
fit.xy_cincinnus_flowers
plot(fit.xy_cincinnus_flowers)
plot(fit.xy_cincinnus_flowers, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("cincinnus_flowers_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("cincinnus_flowers_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_flowers)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("cincinnus_flowers_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_flowers, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (15) cincinnus vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_cincinnus_claw <-fitPagel(tree,cincinnus,claw, model = "ARD", method = "fitDiscrete")
fit.xy_cincinnus_claw
plot(fit.xy_cincinnus_claw)
plot(fit.xy_cincinnus_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("cincinnus_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("cincinnus_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("cincinnus_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (16) cincinnus vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_cincinnus_labellum <-fitPagel(tree,cincinnus,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_cincinnus_labellum
plot(fit.xy_cincinnus_labellum)
plot(fit.xy_cincinnus_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("cincinnus_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("cincinnus_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("cincinnus_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (17) cincinnus vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_cincinnus_blotch <-fitPagel(tree,cincinnus,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_cincinnus_blotch
plot(fit.xy_cincinnus_blotch)
plot(fit.xy_cincinnus_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("cincinnus_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("cincinnus_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("cincinnus_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (18) cincinnus vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_cincinnus_stigma <-fitPagel(tree,cincinnus,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_cincinnus_stigma
plot(fit.xy_cincinnus_stigma)
plot(fit.xy_cincinnus_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("cincinnus_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(cincinnus[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("cincinnus_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("cincinnus_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_cincinnus_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (19) flowers vs claw
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_flowers_claw <-fitPagel(tree,flowers,claw, model = "ARD", method = "fitDiscrete")
fit.xy_flowers_claw
plot(fit.xy_flowers_claw)
plot(fit.xy_flowers_claw, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("flowers_claw_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("flowers_claw_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_claw)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("flowers_claw_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_claw, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (20) flowers vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_flowers_labellum <-fitPagel(tree,flowers,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_flowers_labellum
plot(fit.xy_flowers_labellum)
plot(fit.xy_flowers_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("flowers_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("flowers_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("flowers_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (21) flowers vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_flowers_blotch <-fitPagel(tree,flowers,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_flowers_blotch
plot(fit.xy_flowers_blotch)
plot(fit.xy_flowers_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("flowers_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("flowers_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("flowers_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (22) flowers vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_flowers_stigma <-fitPagel(tree,flowers,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_flowers_stigma
plot(fit.xy_flowers_stigma)
plot(fit.xy_flowers_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("flowers_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(flowers[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("flowers_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("flowers_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_flowers_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (23) claw vs labellum
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_claw_labellum <-fitPagel(tree,claw,labellum, model = "ARD", method = "fitDiscrete")
fit.xy_claw_labellum
plot(fit.xy_claw_labellum)
plot(fit.xy_claw_labellum, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("claw_labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("claw_labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_claw_labellum)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("claw_labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_claw_labellum, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (24) claw vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_claw_blotch <-fitPagel(tree,claw,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_claw_blotch
plot(fit.xy_claw_blotch)
plot(fit.xy_claw_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("claw_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("claw_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_claw_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("claw_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_claw_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (25) claw vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_claw_stigma <-fitPagel(tree,claw,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_claw_stigma
plot(fit.xy_claw_stigma)
plot(fit.xy_claw_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("claw_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(claw[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("claw_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_claw_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("claw_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_claw_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (26) labellum vs blotch
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_labellum_blotch <-fitPagel(tree,labellum,blotch, model = "ARD", method = "fitDiscrete")
fit.xy_labellum_blotch
plot(fit.xy_labellum_blotch)
plot(fit.xy_labellum_blotch, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("labellum_blotch_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("labellum_blotch_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_labellum_blotch)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("labellum_blotch_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_labellum_blotch, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (27) labellum vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_labellum_stigma <-fitPagel(tree,labellum,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_labellum_stigma
plot(fit.xy_labellum_stigma)
plot(fit.xy_labellum_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("labellum_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(labellum[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("labellum_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_labellum_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("labellum_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_labellum_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()


## (28) blotch vs stigma
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("0","1")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(stigma[tree$tip.label],c("0","1")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
fit.xy_blotch_stigma <-fitPagel(tree,blotch,stigma, model = "ARD", method = "fitDiscrete")
fit.xy_blotch_stigma
plot(fit.xy_blotch_stigma)
plot(fit.xy_blotch_stigma, lwd.by.rate=TRUE)

## Saving files in PDF format
# 1. Open pdf file
pdf("blotch_stigma_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(blotch[tree$tip.label],c("0","1")),piecol=c("blue","red"),
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
pdf("blotch_stigma_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_blotch_stigma)
# 3. Close the file
dev.off()

# 1. Open pdf file
pdf("blotch_stigma_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(fit.xy_blotch_stigma, lwd.by.rate=TRUE)
# 3. Close the file
dev.off()