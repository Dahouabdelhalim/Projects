# R script for the analyses in 
# "Rates of morphological evolution, asymmetry and morphological integration of shell shape in scallops"
# Sherratt, Adams & Serb
# BMC Evolutionary Biology
# Accepted November 2017
# Code written April 30 2014; Modified October 18 2016 & April 18 2017

# Written in R v.3.3.1 
library(geomorph) # v3.0.4
# loads library(ape) # v4.1
library(geiger) # v.2.0.6
library(gplots)

# Read in specimen data
indivL.y <- as.matrix(read.csv("avspecimen_coords_L.csv", header =T, row.names = 1))
indivR.y <- as.matrix(read.csv("avspecimen_coords_R.csv", header =T, row.names = 1))

# Read in classifier
indiv.class <- read.csv("ID classifier all.csv", header = T, row.names = 1)
  indivL.class <- indiv.class[rownames(indivL.y),] #subsets for lefts in this dataset
  #933 individuals
  indivR.class <- indiv.class[rownames(indivR.y),] #subsets for rights in this dataset
  #859

# Remove 'L' or 'R' from names, and replace dimnames with new name
IDs.L <- row.names(indivL.class)
  IDs.L <- sub("L", "", IDs.L) ; IDs.L <- sub("X", "", IDs.L)
  row.names(indivL.y)<- IDs.L 
IDs.R <- row.names(indivR.class)
  IDs.R <- sub("R", "", IDs.R) ; IDs.R <- sub("X", "", IDs.R)
  row.names(indivR.y)<- IDs.R

# Find specimens common to both datasets and remove specimens that are unpaired
symm <- IDs.L[IDs.L %in% IDs.R] # which species in right are also in left
indivL.y <- indivL.y[symm,] 
indivR.y <- indivR.y[symm,]
# now 711 individuals with BOTH left and right valves
indiv.class <- indivL.class[symm,] 
  row.names(indiv.class) <- row.names(indivL.y)

# Subset specimen data and classifier to have only those species in the Phylogeny
  symm <- symm[which(indiv.class$phylo_ID != "")]
  indivL.y <- indivL.y[which(indiv.class$phylo_ID != ""),]
  indivR.y <- indivR.y[which(indiv.class$phylo_ID != ""),]
  indiv.class <- indiv.class[rownames(indivL.y),]
  # now 669 specimens; still in alphabetical order
  indiv.class$habit <- factor(indiv.class$habit, levels= c("cement", "nestle", "byssal", "recess", "free", "glide"))
  indiv.class$phylo_ID <- factor(indiv.class$phylo_ID)
  
# Read in phylogenetic tree, already pruned to match dataset
  tree <- read.tree("86sp scallop tree.tre")
  tree.tips <- tree$tip.label
  # 86 species
  
# Average specimen data by species & reorder by tree
  specL.y <- aggregate(indivL.y ~ indiv.class$phylo_ID,FUN= mean)
    rownames(specL.y) <- specL.y[,1] ; specL.y <- as.matrix(specL.y[,-1]) #remove first column and use as rownames
    specL.y <- specL.y[tree$tip.label,] # reorder to match tree
  specR.y <- aggregate(indivR.y ~ indiv.class$phylo_ID,FUN= mean)
    rownames(specR.y) <- specR.y[,1] ; specR.y <- as.matrix(specR.y[,-1])
    specR.y <- specR.y[tree$tip.label,]
    specL.Y <- arrayspecs(specL.y, 202, 3)
    specR.Y <- arrayspecs(specR.y, 202, 3)
    
# Get species IDs and subset to the 86 species & reorder by tree
spec.class <- read.csv("Spec classifier_ALL.csv", header = T, row.names = 1)
spec.class <- spec.class[match(row.names(specL.y), spec.class$phylo_ID),]
spec.class$phylo_ID <- factor(spec.class$phylo_ID, levels=tree$tip.label) # reset levels to remove empty
    
# Habit grouping variable
habit <- factor(spec.class$habit, levels=c("cement", "nestle", "byssal", "recess", "free", "glide")) 
names(habit) <- spec.class$phylo_ID
  # how many?
summary(habit)
# cement nestle byssal recess   free  glide 
#   2      1     48     11     18      6  
# Species
species <- spec.class$phylo_ID

################################################################################################
# First, examine integration patterns between left and right valve
ind.PLS <- two.b.pls(indivL.y, indivR.y, iter=1)
ind.PLS.shapes <- plot(ind.PLS, shapes = TRUE)

colours <- c( "#00CD66", "red" ,"#4876FF", "purple","black","orange")
names(colours) <- levels(indiv.class$habit )
indiv.habit.col <- colours[match(indiv.class$habit, names(colours))]

plot(ind.PLS$XScores[,1], ind.PLS$YScores[,1], pch=21, col="black", cex=1, bg=indiv.habit.col, 
     xlab="PLS1 left valve", ylab="PLS1 right valve")

# Test for phylogenetic Integration between left and right valve
pcax <- prcomp(specL.y)
d <- which(zapsmall(pcax$sdev) > 0.0001)
x <- pcax$x[,d]

pcay <- prcomp(specR.y)
d <- which(zapsmall(pcay$sdev) > 0.0001)
y <- pcay$x[,d]

phylo.integration(x,y, tree, iter = 999)
# r-PLS: 0.803
# P-value: 0.001
# Based on 1000 random permutations

# Subset by group
L.coords.gp <- coords.subset(specL.Y[,,-which(habit=="nestle")], factor(habit[-which(habit=="nestle")]))
R.coords.gp <- coords.subset(specR.Y[,,-which(habit=="nestle")], factor(habit[-which(habit=="nestle")]))

integ.tests <- Map(function(x,y) integration.test(x, y, iter=999), L.coords.gp, R.coords.gp)
# $byssal
# 
# Call:
#   integration.test(A = x, A2 = y, iter = 999) 

# r-PLS: 0.882
# 
# P-value: 0.001
# 
# Based on 1000 random permutations
# $recess
# 
# Call:
#   integration.test(A = x, A2 = y, iter = 999) 

# r-PLS: 0.882
# 
# P-value: 0.006
# 
# Based on 1000 random permutations
# $free
# 
# Call:
#   integration.test(A = x, A2 = y, iter = 999) 

# r-PLS: 0.946
# 
# P-value: 0.001
# 
# Based on 1000 random permutations
# $glide
# 
# Call:
#   integration.test(A = x, A2 = y, iter = 999) 

# r-PLS: 0.943
# 
# P-value: 0.048
# 
# Based on 1000 random permutations

# Use Adams & Collyer 2016 Comparisons of Effect Sizes from Partial Least Squares
compare.pls(integ.tests)
  # Effect sizes
  # 
  # cement   byssal   recess     free    glide 
  # 0.9870831 6.4578693 2.6144598 4.4156495 1.4685588 
  # 
  # Effect sizes for pairwise differences in PLS effect size
  # 
  #         cement    byssal    recess      free     glide
  # cement 0.0000000 0.5861513 0.6919401 0.5570565 0.7988822
  # byssal 0.5861513 0.0000000 0.8584409 0.2426747 1.5678478
  # recess 0.6919401 0.8584409 0.0000000 0.9280193 0.6502578
  # free   0.5570565 0.2426747 0.9280193 0.0000000 1.5570038
  # glide  0.7988822 1.5678478 0.6502578 1.5570038 0.0000000
  # 
  # P-values
  # 
  #         cement     byssal    recess       free      glide
  # cement 1.0000000 0.27888691 0.2444875 0.28874442 0.21217936
  # byssal 0.2788869 1.00000000 0.1953245 0.40412872 0.05845833
  # recess 0.2444875 0.19532452 1.0000000 0.17669878 0.25776285
  # free   0.2887444 0.40412872 0.1766988 1.00000000 0.05973479
  # glide  0.2121794 0.05845833 0.2577629 0.05973479 1.00000000
  # No significant differences

compare.pls(integ.tests[2:5])
  # Effect sizes
  # 
  # byssal   recess     free    glide
  # 6.457869 2.614460 4.415650 1.468559 
  # 
  # Effect sizes for pairwise differences in PLS effect size
  # 
  #           byssal    recess      free     glide
  # byssal 0.0000000 0.8584409 0.2426747 1.5678478
  # recess 0.8584409 0.0000000 0.9280193 0.6502578
  # free   0.2426747 0.9280193 0.0000000 1.5570038
  # glide  1.5678478 0.6502578 1.5570038 0.0000000
  # 
  # P-values
  # 
  #           byssal    recess       free      glide
  # byssal 1.00000000 0.1953245 0.40412872 0.05845833
  # recess 0.19532452 1.0000000 0.17669878 0.25776285
  # free   0.40412872 0.1766988 1.00000000 0.05973479
  # glide  0.05845833 0.2577629 0.05973479 1.00000000
  # No significant differences

############################################################################################   
## RATES Q 1: Does the rate of left-right asymmetry evolution differ among life habits? 
# Measure of asymmetry between left and right valves |R-L| = Procrustes distance between left and right
# Import Procrustes aligned data from a Procrustes fit of BOTH lefts and flipped rights 
L.Y <- readland.tps("avspecimen_coords_L_ProcLR.tps", specID = "ID")
dimnames(L.Y)[[3]] <- IDs.L # gives new names, without L or R
R.Y <- readland.tps("avspecimen_coords_R_ProcLR.tps", specID = "ID")
dimnames(R.Y)[[3]] <- IDs.R # gives new names, without L or R
# alphabetical datasets

# Subset to include only those with both left and right and in the phylogeny (i.e. symm)
L.Y <- L.Y[,,symm] ; L.y <- two.d.array(L.Y)
R.Y <- R.Y[,,symm] ; R.y <- two.d.array(R.Y)

# Measure asymmetry as Procrustes distance between L and R
asymm <- matrix(NA, ncol=1, nrow=nrow(L.y))
for (i in 1:nrow(asymm)){ asymm[i,] <- dist(rbind(L.y[i,], R.y[i,])) }
rownames(asymm)<- rownames(L.y)

phylo.asymm <- as.matrix(asymm[which(indiv.class$phylo_ID !=""),])

# Plot a histogram of the amount of asymmetry across individuals
hist(asymm, col="grey")
mean(asymm) # 0.09270106
range(asymm) # 0.02823582 0.31282148

# Is the amount of asymmetry in individuals different among habits?
# Plot save as 660 x 660px
boxplot(asymm~indiv.class$habit, ylab = "Procrustes distance")

# Plot supp materials graph save as 643 x 800px
par(las = 1, mar=c(4,9,0,2), cex=0.5)
boxplot(asymm~factor(indiv.class$phylo_ID, levels=tree$tip.label), xlab = "Procrustes distance", col= habit.col, horizontal = TRUE, names = spec.class$phylo_ID2)

# How does the amount of asymmetry vary over a phylogeny?
# Average asymmetry measure by species
phylo.asymm.AV <- tapply(phylo.asymm, factor(indiv.class$phylo_ID), mean)

procD.pgls(as.matrix(phylo.asymm.AV) ~ habit, tree, iter=999)
# Type I (Sequential) Sums of Squares and Cross-products
# Randomized Residual Permutation Procedure Used
# 1000 Permutations
# 
#           Df        SS         MS      Rsq      F       Z Pr(>F)
# habit      5 0.0002858 5.7150e-05 0.069486 1.1948 0.16101  0.805
# Residuals 80 0.0038266 4.7832e-05                               
# Total     85 0.0041123  
boxplot(phylo.asymm.AV~habit, ylab = "Procrustes distance")

# How does rate of asymmetry vary across habits?
aymm.shape.rate <- compare.evol.rates(as.matrix(phylo.asymm.AV), tree, habit, iter=999)
# Observed Rate Ratio: 28.117
# 
# P-value: 0.573
# 
# Based on 1000 random permutations
# 
# The rate for group cement is 2.47445332905639e-05  
# 
# The rate for group nestle is 4.7456644933818e-05  
# 
# The rate for group byssal is 9.14523434393924e-06 
# 
# The rate for group recess is 0.000148275014520253  
# 
# The rate for group free is 2.23653600394923e-05  
# 
# The rate for group glide is 0.000257136450017144 

#Bootstrap individuals and recalculate species averages
iter=100
sigmad.rand = NULL
for(j in 1:iter){
  asymm.rand =NULL
  for(i in 1:length(species)){
    x <- phylo.asymm[which(indiv.class$phylo_ID == species[i]),]
    if(!is.null(length(x))){ x <- mean(x[sample(length(x), replace=T)]) }
    asymm.rand <- c(asymm.rand, x)
  }
  names(asymm.rand) <- species
  sigmad <- compare.evol.rates(as.matrix(asymm.rand), tree, habit, iter=1, print.progress = F)$sigma.d.gp #calculate sigmas
  sigmad.rand <- rbind(sigmad.rand, sigmad)
}
remove(x, sigmad, asymm.rand)

#Calculate confidence intervals per group
STD.L <- apply(sigmad.rand, MARGIN=2, sd)
CI.u.L <- aymm.shape.rate$sigma.d.gp + STD.L
CI.l.L <- aymm.shape.rate$sigma.d.gp - STD.L

# Plot save as 660 x 660px
plotCI(aymm.shape.rate$sigma.d.gp, li=CI.l.L, ui=CI.u.L, xlab="Life Habit",  
       ylab=expression(paste(sigma^2)), pch=21,cex=1.5, pt.bg = "black")
abline(h = aymm.shape.rate$sigma.d.all, lty=2) # plots the overal rate

############################################################################################
## RATES Q 2: Compare rates between lefts and rights pooled (no habit info)
gp = c(rep("L", 606), rep("R", 606)) # variables left, variables right valve
# Compare Lefts versus Rights
shape.rateLvR <- compare.multi.evol.rates(cbind(specL.y, specR.y), gp, tree, Subset = FALSE, iter=999) 
  # Observed Rate Ratio: 1.1728
  # 
  # P-value: 0.001
  # 
  # Based on 1001 random permutations
  # 
  # The rate for group L is 0.000104007584391734  
  # 
  # The rate for group R is 0.000121978541994984 
# must correct for number of variables by dividing by 606
shape.rateLvR$sigma.d.gp / 606
# L 1.716297e-07  R 2.012847e-07

# Nestler and cementers can't be analysed here because of N, but
sigmad.R/sigmad.L  
# Results not used, because pooling does not make sense given the difference between L&R


############################################################################################
## RATES Q 3: Compare evolutionary rates of valve SHAPE among habit groups for left and right
  ### a) LEFT VALVE ###
  habits.shape.rateL <- compare.evol.rates(specL.y, tree, habit, iter=999)  
  # Observed Rate Ratio: 2.9215
  # 
  # P-value: 0.028
  # 
  # Based on 1000 random permutations
  # 
  # The rate for group cement is 1.58848136687151e-07  
  # The rate for group nestle is 1.04972094086336e-07  
  # The rate for group byssal is 1.16966416392516e-07  
  # The rate for group recess is 3.03109804142868e-07  
  # The rate for group free is 1.97156153479448e-07  
  # The rate for group glide is 3.06679546318017e-07
  
  ### b) RIGHT VALVE ###
  habits.shape.rateR <- compare.evol.rates(specR.y, tree, habit, iter=999)
  # Observed Rate Ratio: 5.1553
  # 
  # P-value: 0.001
  # 
  # Based on 1000 random permutations
  # 
  # The rate for group cement is 2.74906903376711e-07  
  # The rate for group nestle is 6.11331086808484e-07  
  # The rate for group byssal is 1.18582966700041e-07  
  # The rate for group recess is 4.33915870094619e-07  
  # The rate for group free is 1.62656159568322e-07  
  # The rate for group glide is 4.59412236544915e-07 
   
  ### Exporting pairwsie difference P-value as a table:
  pairwise.pvalue <- as.matrix(habits.shape.rateL$pairwise.pvalue) # lefts lower tri
  pairwise.pvalue[upper.tri(pairwise.pvalue)] <- as.matrix(habits.shape.rateR$pairwise.pvalue)[upper.tri(habits.shape.rateR$pairwise.pvalue)] # right upper tri
  # write.csv(pairwise.pvalue, "pairwise Pvalue LR habits.csv")
  
  #Bootstrap individuals and recalculate species averages for left and right
  iter = 100
  sigmad.randL = NULL ; sigmad.randR = NULL
  for(j in 1:iter){
    specL.y.rand =NULL ; specR.y.rand =NULL
    for(i in 1:length(species)){
      x <- indivL.y[which(indiv.class$phylo_ID == species[i]),]
      if(!is.null(nrow(x))){ x <- colMeans(x[sample(nrow(x), replace=T),]) }
      specL.y.rand <- rbind(specL.y.rand, x) ; specR.y.rand <- rbind(specR.y.rand, x)
    }
    rownames(specL.y.rand) <- species ; rownames(specR.y.rand) <- species
    sigmadL <- compare.evol.rates(specL.y.rand, tree, habit, iter=1, print.progress = FALSE)$sigma.d.gp #calculate sigmas
    sigmadR <- compare.evol.rates(specR.y.rand, tree, habit, iter=1, print.progress = FALSE)$sigma.d.gp #calculate sigmas
    sigmad.randL <- rbind(sigmad.randL, sigmadL) ; sigmad.randR <- rbind(sigmad.randR, sigmadR)
  } 
  remove(specL.y.rand, sigmadL, specR.y.rand, sigmadR)
  
  # Set up confidence intervals
  sigmad.L <- habits.shape.rateL$sigma.d.gp; sigmad.R <- habits.shape.rateR$sigma.d.gp
  
  #Calculate confidence intervals per group
  STD.L <- apply(sigmad.randL, MARGIN=2, sd)
  CI.u.L <- sigmad.L + STD.L ; CI.l.L <- sigmad.L - STD.L
  STD.R <- apply(sigmad.randR, MARGIN=2, sd)
  CI.u.R <- sigmad.R + STD.R ; CI.l.R <- sigmad.R - STD.R
  
  # Plot on single graph
  sigmad <- as.vector(rbind(sigmad.L, sigmad.R)) # binds as pairs
  names(sigmad) <- as.vector(rbind(names(sigmad.L), names(sigmad.R))) # check
  CI.l <- as.vector(rbind(CI.l.L, CI.l.R)) 
  CI.u <- as.vector(rbind(CI.u.L, CI.u.R)) 
  pch <- rep(c(21,22),12)
  options(warn=-1)
  # save as a 660 x 660 px
  plotCI(sigmad, li= CI.l, ui= CI.u ,xlab="Life Habit", 
         ylab=expression(paste(sigma^2)), pch=pch ,cex=1.2, pt.bg = c("black", "grey"))
  abline(h = habits.shape.rateL$sigma.d.all, lty=2)
  abline(h = habits.shape.rateR$sigma.d.all, lty=3)

############################################################################################  
## RATES Q 4: Rates in Left versus right valve for each life habit  
  # choose only habits with multiple species, so nestle is out
  # using gp, made above, which is a vector denoting L & R for a 2D matrix
  for(i in levels(habit)){
    if(length(which(habit == i)) > 2){
      phy <- drop.tip(tree, tree.tips[which(habit != i)])
      y <- cbind(specL.y[which(habit == i),],specR.y[which(habit == i),])
      res <- compare.multi.evol.rates(y,gp,phy,iter=999, print.progress = F)
      assign(paste(i,"res",sep="."),res)
    } }
  byssal.res
  recess.res
  free.res
  glide.res

  # Cementers
  # 2.74906903376711e-07  / 1.58848136687151e-07 
  # 1.730627
  # Nestlers
  # 6.11331086808484e-07  / 1.04972094086336e-07  
  # 5.823749
  
  # > byssal.res
  # 
  # Call:
  #   
  #   
  # Observed Rate Ratio: 1.0024
  # 
  # P-value: 0.965
  # 
  # Based on 1000 random permutations
  # 
  # The rate for group L is 1.04580293187419e-07  
  # 
  # The rate for group R is 1.04326562342066e-07 
  # >   free.res
  # 
  # Call:
  #   
  #   
  #   Observed Rate Ratio: 1.0685 (L/R) (0.9359306 is R/L)
  # 
  # P-value: 0.5944
  # 
  # Based on 1001 random permutations
  # 
  # The rate for group L is 1.56322719404212e-07  
  # 
  # The rate for group R is 1.46307211590742e-07 
  # >   glide.res
  # 
  # Call:
  #   
  #   
  #   Observed Rate Ratio: 1.9632
  # 
  # P-value: 0.031
  # 
  # Based on 1001 random permutations
  # 
  # The rate for group L is 2.41409920059941e-08  
  # 
  # The rate for group R is 1.22969272755927e-08 
  # >   recess.res
  # 
  # Call:
  #   
  #   
  #   Observed Rate Ratio: 1.6747
  # 
  # P-value: 0.014
  # 
  # Based on 1001 random permutations
  # 
  # The rate for group L is 1.97918319344385e-07  
  # 
  # The rate for group R is 3.31455674388108e-07 

################################################################################################
### Matching symmetry
# Read in raw data replicates & Replace dimnames with new
L.1 <- readland.tps("rawspecimens_202_L_1.tps", specID= "ID"); dimnames(L.1)[[3]]<- IDs.L 
L.2 <- readland.tps("rawspecimens_202_L_2.tps", specID= "ID"); dimnames(L.2)[[3]]<- IDs.L 
R.1 <- readland.tps("rawspecimens_202_R_1.tps", specID= "ID"); dimnames(R.1)[[3]]<- IDs.R
R.2 <- readland.tps("rawspecimens_202_R_2.tps", specID= "ID"); dimnames(R.2)[[3]]<- IDs.R

#subset left and right so only symm dataset individuals are used.
L.1 <- L.1[,,symm] ; L.2 <- L.2[,,symm]
R.1 <- R.1[,,symm] ; R.2 <- R.2[,,symm]

#reflect rights - flip x axis to REFLECT
flip <- array(1, dim=dim(R.1))
flip[,1,] <- rep(-1, 202)
R.1 <- R.1*flip ; R.2 <- R.2*flip

# make vectors to use in bilat.symmetry()
side <- c(rep("L", dim(L.1)[3]*2), rep("R", dim(R.1)[3]*2))
ind <- c(dimnames(L.1)[[3]], dimnames(L.2)[[3]], dimnames(R.1)[[3]], dimnames(R.2)[[3]])
replicate <- c(rep(1, dim(L.1)[3]), rep(2, dim(L.2)[3]), 
               rep(1, dim(R.1)[3]), rep(2, dim(R.2)[3]))
A <- abind::abind(L.1,L.2,R.1,R.2)
match.symm <- bilat.symmetry(A, ind, side, replicate, object.sym=FALSE)

write.csv(match.symm$ANOVA.Shape, "matching.symm.ANOVA.res.csv")

#END



