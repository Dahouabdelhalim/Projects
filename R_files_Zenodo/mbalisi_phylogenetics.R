# Code for Balisi et al. (submitted)
# "Dietary specialization is linked to reduced species durations in North American fossil canids"

# This file contains scripts for phylogenetic analysis of continuous trait data.



# load libraries
library(phytools)
library(geiger)



# read in Slater (2015)'s 500 trees.
# load("Slater2015_data/canid_tree_data.rdata")

# from Slater (2015)'s code in "continuous trait analyses.R": 
# a function to compute AIC and AICc, which are not reported for phytools diversity models 
# aka fitDiversityModel() in phytools ##
aic.comp = function(lnl, k, n = NULL, small.sample = F) {
  AIC <- -2*lnl + 2*k
  if(small.sample ==  F){
    return(AIC)
  } else{
    if(is.null(n)==T) stop("you must specify the sample size to use AIC.c")
    weight <- (2*k * (k+1)) / (n-k-1)
    return(AIC + weight)	
  }
} 

# Taxonomic spelling corrections
best.tree$tip.label[best.tree$tip.label=="Cuon_javanicus"] = "Cuon_alpinus"
best.tree$tip.label[best.tree$tip.label=="Cynarctoides_gawanae"] = "Cynarctoides_gawnae"
best.tree$tip.label[best.tree$tip.label=="Cynarctoides_accridens"] = "Cynarctoides_acridens"
best.tree$tip.label[best.tree$tip.label=="Paracynarctus_sinclari"] = "Paracynarctus_sinclairi"
best.tree$tip.label[best.tree$tip.label=="Phlaocyon_marshlandensis"] = "Phlaocyon_marslandensis"
best.tree$tip.label[best.tree$tip.label=="Protomarctus_opatus"] = "Protomarctus_optatus"
best.tree$tip.label[best.tree$tip.label=="Rhizocyon_oreganensis"] = "Rhizocyon_oregonensis"
best.tree$tip.label[best.tree$tip.label=="Urocyon_galushi"] = "Urocyon_galushai"
best.tree$tip.label[best.tree$tip.label=="Urocyon_citronus"] = "Urocyon_citrinus"

# More corrections through all 500 trees
t.file = .uncompressTipLabel(t.file1)
for(i in 1:length(t.file)) {
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Cuon_javanicus"] = "Cuon_alpinus"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Cynarctoides_gawanae"] = "Cynarctoides_gawnae"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Cynarctoides_accridens"] = "Cynarctoides_acridens"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Paracynarctus_sinclari"] = "Paracynarctus_sinclairi"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Phlaocyon_marshlandensis"] = "Phlaocyon_marslandensis"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Protomarctus_opatus"] = "Protomarctus_optatus"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Rhizocyon_oreganensis"] = "Rhizocyon_oregonensis"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Urocyon_galushi"] = "Urocyon_galushai"
  t.file[[i]]$tip.label[t.file[[i]]$tip.label=="Urocyon_citronus"] = "Urocyon_citrinus"
}
t.file1 = .compressTipLabel(t.file, ref=NULL)



# Added Oct 7, 2017
# Prune trees by subfamily

# read in file containing body mass and saved principal-component trait data
dd = read.csv("mbalisi_massPCsubfam_fossilsOnly.csv", stringsAsFactors=F, row.names=1)

# start with best.tree
Hesperocyoninae = treedata(best.tree, dd[dd$subfamily=="Hesperocyoninae", ], sort=T)$phy
Borophaginae = treedata(best.tree, dd[dd$subfamily=="Borophaginae", ], sort=T)$phy
Caninae = treedata(best.tree, dd[dd$subfamily=="Caninae", ], sort=T)$phy

# now try with all 500 trees
t.file = .uncompressTipLabel(t.file1)
Hesp.t.file = Boro.t.file = Cani.t.file = t.file
for(i in 1:length(t.file)) {
  Hesp.t.file[[i]] = treedata(t.file[[i]], dd[dd$subfamily=="Hesperocyoninae", ], sort=T)$phy
  Boro.t.file[[i]] = treedata(t.file[[i]], dd[dd$subfamily=="Borophaginae", ], sort=T)$phy
  Cani.t.file[[i]] = treedata(t.file[[i]], dd[dd$subfamily=="Caninae", ], sort=T)$phy
}
t.file1 = .compressTipLabel(t.file, ref=NULL)
Hesp.t.file1 = .compressTipLabel(Hesp.t.file, ref=NULL)
Boro.t.file1 = .compressTipLabel(Boro.t.file, ref=NULL)
Cani.t.file1 = .compressTipLabel(Cani.t.file, ref=NULL)



# set up results matrices

bmModelResMass = bmModelResPC1 = acdcModelResMass = acdcModelResPC1 = 
  trendResMass = trendResPC1 = driftResMass = driftResPC1 = divModelResMass = 
  divModelResPC1 = ouModelResMass = ouModelResPC1 = 
  bmModelResMassHesp = bmModelResPC1Hesp = acdcModelResMassHesp = 
  acdcModelResPC1Hesp = trendResMassHesp = trendResPC1Hesp = driftResMassHesp = 
  driftResPC1Hesp = divModelResMassHesp = divModelResPC1Hesp = 
  ouModelResMassHesp = ouModelResPC1Hesp = 
  bmModelResMassBoro = bmModelResPC1Boro = acdcModelResMassBoro = 
  acdcModelResPC1Boro = trendResMassBoro = trendResPC1Boro = driftResMassBoro = 
  driftResPC1Boro = divModelResMassBoro = divModelResPC1Boro = 
  ouModelResMassBoro = ouModelResPC1Boro = 
  bmModelResMassCani = bmModelResPC1Cani = acdcModelResMassCani = 
  acdcModelResPC1Cani = trendResMassCani = trendResPC1Cani = driftResMassCani = 
  driftResPC1Cani = divModelResMassCani = divModelResPC1Cani = ouModelResMassCani = 
  ouModelResPC1Cani = 
  matrix(data=NA, nrow=500, ncol=5, 
         dimnames=list(1:500, c("sample", "LnL", "aicc", "Sig2", "param")))

weightsMass = weightsPC1 = weightsMassHesp = weightsPC1Hesp = weightsMassBoro = 
  weightsPC1Boro = weightsMassCani = weightsPC1Cani = 
  matrix(data=NA, nrow=500, ncol=6, 
         dimnames=list(1:500, c("BM", "ACDC", "trend", "drift", "diversity", "OU")))



# prepare body mass and dietary data using "best tree"
d = dd[, 1:4]



# create two sets of data: one for mass and one for PC1
mass = treedata(best.tree, d[is.na(d[, 1])==F, ], sort=TRUE)$data
mass = mass[ , 1]
nMass = length(mass)
PC1 = treedata(best.tree, d[is.na(d[, 2])==F, ], sort=TRUE)$data
PC1 = PC1[ , 2]
nPC1 = length(PC1)



# read in the same randomly generated 500 samples from the posterior distribution
# as in Slater (2015)
samples = read.csv("Slater2015_data/posteriorSample.csv")[, 1];



# loop trait-evolution model analysis through all 500 trees as in Slater (2015)
for(i in 1:500) {
  
  tree = samples[i]
  phy.tmp = t.file1[[tree]]
  phy.tmp$edge.length = phy.tmp$edge.length / p.file1[tree, "Clockrate"]
  
  # prepare trees for mass and PC1 separately
  # because not all species have both mass and PC1
  phyMass = treedata(phy.tmp, mass)$phy
  phyPC1 = treedata(phy.tmp, PC1)$phy

  # fit models
  # mass first, then PC1
  
  # Bm 
  bmMass = fitContinuous(phyMass, mass, model="BM", control=list(niter=10))$opt
  bmModelResMass[i, ] = c(tree, bmMass$lnL, bmMass$aicc, bmMass$sigsq, NA)
  
  bmPC1 = fitContinuous(phyPC1, PC1, model="BM", control=list(niter=10))$opt
  bmModelResPC1[i, ] = c(tree, bmPC1$lnL, bmPC1$aicc, bmPC1$sigsq, NA)
  
  # ACDC 
  acdcMass = fitContinuous(phyMass, mass, model="EB", bounds=list(a=c(-0.2, 0.2)), 
                           control=list(niter=10))$opt
  acdcModelResMass[i, ] = c(tree, acdcMass$lnL, acdcMass$aicc, acdcMass$sigsq, 
                            acdcMass$a)
  
  acdcPC1 = fitContinuous(phyPC1, PC1, model="EB", bounds=list(a=c(-0.2, 0.2)), 
                          control=list(niter=10))$opt
  acdcModelResPC1[i, ] = c(tree, acdcPC1$lnL, acdcPC1$aicc, acdcPC1$sigsq, 
                           acdcPC1$a)
  
  # Trend
  trendMass = fitContinuous(phyMass, mass, model="trend", control=list(niter=10))$opt
  trendResMass[i, ] = c(tree, trendMass$lnL, trendMass$aicc, trendMass$sigsq, 
                        trendMass$slope)
  
  trendPC1 = fitContinuous(phyPC1, PC1, model="trend", control=list(niter=10))$opt
  trendResPC1[i, ] = c(tree, trendPC1$lnL, trendPC1$aicc, trendPC1$sigsq, 
                       trendMass$slope)
  
  # Drift
  driftMass = fitContinuous(phyMass, mass, model="drift", control=list(niter=10))$opt
  driftResMass[i, ] = c(tree, driftMass$lnL, driftMass$aicc, driftMass$sigsq, 
                        driftMass$drift)
  
  driftPC1 = fitContinuous(phyPC1, PC1, model="drift", control=list(niter=10))$opt
  driftResPC1[i, ] = c(tree, driftPC1$lnL, driftPC1$aicc, driftPC1$sigsq, 
                       driftPC1$drift)
  
  # diversity NOT by diet
  divMass = fitDiversityModel(phyMass, mass, showTree=F)
  divMass$aicc = aic.comp(divMass$logL, k=3, n=nMass, small.sample=T)
  divModelResMass[i, ] = c(tree, divMass$logL, divMass$aicc, divMass$sig0, 
                           divMass$psi)
  
  divPC1 = fitDiversityModel(phyPC1, PC1, showTree=F)
  divPC1$aicc = aic.comp(divPC1$logL, k=3, n=nPC1, small.sample=T)
  divModelResPC1[i, ] = c(tree, divPC1$logL, divPC1$aicc, divPC1$sig0, divPC1$psi)
  
  # OU
  ouMass = fitContinuous(phyMass, mass, model="OU", control=list(niter=10))$opt
  ouModelResMass[i, ] = c(tree, ouMass$lnL, ouMass$aicc, ouMass$sigsq, ouMass$alpha)
  
  ouPC1 = fitContinuous(phyPC1, PC1, model="OU", control=list(niter=10))$opt
  ouModelResPC1[i, ] = c(tree, ouPC1$lnL, ouPC1$aicc, ouPC1$sigsq, ouMass$alpha)
  
  aicMass = setNames(c(bmMass$aicc, acdcMass$aicc, trendMass$aicc, driftMass$aicc, 
                       divMass$aicc, ouMass$aicc), 
                     c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsMass[i, ] = aicw(aicMass)[,3]
  
  aicPC1 = setNames(c(bmPC1$aicc, acdcPC1$aicc, trendPC1$aicc, driftPC1$aicc, 
                      divPC1$aicc, ouPC1$aicc), 
                    c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsPC1[i, ] = aicw(aicPC1)[,3];
  
  if(i%%10 == 0) {
    print(paste("Samples complete =", i))
  }	
}



# prepare body mass and dietary data using per-subfamily trees

# create three sets of data: Hesperocyoninae, Borophaginae, and Caninae

massHesp = treedata(Hesperocyoninae, d[is.na(d[, 1])==F, ], sort=TRUE)$data
massHesp = massHesp[ , 1]
nMassHesp = length(massHesp)
PC1Hesp = treedata(Hesperocyoninae, d[is.na(d[, 2])==F, ], sort=TRUE)$data
PC1Hesp = PC1Hesp[ , 2]
nPC1Hesp = length(PC1Hesp)

massBoro = treedata(Borophaginae, d[is.na(d[, 1])==F, ], sort=TRUE)$data
massBoro = massBoro[ , 1]
nMassBoro = length(massBoro)
PC1Boro = treedata(Borophaginae, d[is.na(d[, 2])==F, ], sort=TRUE)$data
PC1Boro = PC1Boro[ , 2]
nPC1Boro = length(PC1Boro)

massCani = treedata(Caninae, d[is.na(d[, 1])==F, ], sort=TRUE)$data
massCani = massCani[ , 1]
nMassCani = length(massCani)
PC1Cani = treedata(Caninae, d[is.na(d[, 2])==F, ], sort=TRUE)$data
PC1Cani = PC1Cani[ , 2]
nPC1Cani = length(PC1Cani)



# HESPEROCYONINAE

samples = read.csv("Slater2015_data/posteriorSample.csv")[, 1]
# same 500 samples as above

for(i in 1:500) {
  
  tree = samples[i]
  phy.tmp = Hesp.t.file1[[tree]]
  phy.tmp$edge.length = phy.tmp$edge.length / p.file1[tree, "Clockrate"]
  
  # prep trees for mass and PC1 separately
  phyMass = treedata(phy.tmp, massHesp)$phy
  phyPC1 = treedata(phy.tmp, PC1Hesp)$phy

  # fit models
  # first mass, then PC1
  
  # Bm 
  bmMassHesp = fitContinuous(phyMass, massHesp, model="BM", control=list(niter=10))$opt
  bmModelResMassHesp[i, ] = c(tree, bmMassHesp$lnL, bmMassHesp$aicc, bmMassHesp$sigsq, 
                              NA)
  
  bmPC1Hesp = fitContinuous(phyPC1, PC1Hesp, model="BM", control=list(niter=10))$opt
  bmModelResPC1Hesp[i, ] = c(tree, bmPC1Hesp$lnL, bmPC1Hesp$aicc, bmPC1Hesp$sigsq, 
                             NA)
  
  # ACDC 
  acdcMassHesp = fitContinuous(phyMass, massHesp, model="EB", 
                               bounds=list(a=c(-0.2, 0.2)),
                               control=list(niter=10))$opt
  acdcModelResMassHesp[i, ] = c(tree, acdcMassHesp$lnL, acdcMassHesp$aicc, 
                                acdcMassHesp$sigsq, acdcMassHesp$a)
  
  acdcPC1Hesp = fitContinuous(phyPC1, PC1Hesp, model="EB", 
                              bounds=list(a=c(-0.2, 0.2)), 
                              control=list(niter=10))$opt
  acdcModelResPC1Hesp[i, ] = c(tree, acdcPC1Hesp$lnL, acdcPC1Hesp$aicc, 
                               acdcPC1Hesp$sigsq, acdcPC1Hesp$a)
  
  # Trend
  trendMassHesp = fitContinuous(phyMass, massHesp, model="trend", 
                                control=list(niter=10))$opt
  trendResMassHesp[i, ] = c(tree, trendMassHesp$lnL, trendMassHesp$aicc, 
                            trendMassHesp$sigsq, trendMassHesp$slope)
  
  trendPC1Hesp = fitContinuous(phyPC1, PC1Hesp, model="trend", 
                               control=list(niter=10))$opt
  trendResPC1Hesp[i, ] = c(tree, trendPC1Hesp$lnL, trendPC1Hesp$aicc, 
                           trendPC1Hesp$sigsq, trendMassHesp$slope)
  
  # Drift
  driftMassHesp = fitContinuous(phyMass, massHesp, model="drift", 
                                control=list(niter=10))$opt
  driftResMassHesp[i, ] = c(tree, driftMassHesp$lnL, driftMassHesp$aicc, 
                            driftMassHesp$sigsq, driftMassHesp$drift)
  
  driftPC1Hesp = fitContinuous(phyPC1, PC1Hesp, model="drift", 
                               control=list(niter=10))$opt
  driftResPC1Hesp[i, ] = c(tree, driftPC1Hesp$lnL, driftPC1Hesp$aicc, 
                           driftPC1Hesp$sigsq, driftPC1Hesp$drift)
  
  # diversity NOT by diet
  divMassHesp = fitDiversityModel(phyMass, massHesp, showTree=F)
  divMassHesp$aicc = aic.comp(divMassHesp$logL, k=3, n=nMassHesp, small.sample=T)
  divModelResMassHesp[i, ] = c(tree, divMassHesp$logL, divMassHesp$aicc, 
                               divMassHesp$sig0, divMassHesp$psi)
  
  divPC1Hesp = fitDiversityModel(phyPC1, PC1Hesp, showTree=F)
  divPC1Hesp$aicc = aic.comp(divPC1Hesp$logL, k=3, n=nPC1Hesp, small.sample=T)
  divModelResPC1Hesp[i, ] = c(tree, divPC1Hesp$logL, divPC1Hesp$aicc, 
                              divPC1Hesp$sig0, divPC1Hesp$psi)
  
  # OU
  ouMassHesp = fitContinuous(phyMass, massHesp, model="OU", control=list(niter=10))$opt
  ouModelResMassHesp[i, ] = c(tree, ouMassHesp$lnL, ouMassHesp$aicc, 
                              ouMassHesp$sigsq, ouMassHesp$alpha)
  
  ouPC1Hesp = fitContinuous(phyPC1, PC1Hesp, model="OU", control=list(niter=10))$opt
  ouModelResPC1Hesp[i, ] = c(tree, ouPC1Hesp$lnL, ouPC1Hesp$aicc, ouPC1Hesp$sigsq, 
                             ouMassHesp$alpha)
  
  aicMassHesp = setNames(c(bmMassHesp$aicc, acdcMassHesp$aicc, trendMassHesp$aicc, 
                           driftMassHesp$aicc, divMassHesp$aicc, ouMassHesp$aicc), 
                     c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsMassHesp[i, ] = aicw(aicMassHesp)[,3]
  
  aicPC1Hesp = setNames(c(bmPC1Hesp$aicc, acdcPC1Hesp$aicc, trendPC1Hesp$aicc, 
                          driftPC1Hesp$aicc, divPC1Hesp$aicc, ouPC1Hesp$aicc), 
                    c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsPC1Hesp[i, ] = aicw(aicPC1Hesp)[,3]
  
  if(i%%10 == 0) {
    print(paste("Samples complete =", i))
  }	
}



# BOROPHAGINAE

samples = read.csv("Slater2015_data/posteriorSample.csv")[, 1]
# same 500 samples as above

for(i in 1:500) {
  
  tree = samples[i]
  phy.tmp = Boro.t.file1[[tree]]
  phy.tmp$edge.length = phy.tmp$edge.length / p.file1[tree, "Clockrate"]
  
  # prep trees
  phyMass = treedata(phy.tmp, massBoro)$phy
  phyPC1 = treedata(phy.tmp, PC1Boro)$phy

  # fit models
  # first mass, then PC1
  
  # Bm 
  bmMassBoro = fitContinuous(phyMass, massBoro, model="BM", 
                             control=list(niter=10))$opt
  bmModelResMassBoro[i, ] = c(tree, bmMassBoro$lnL, bmMassBoro$aicc, 
                              bmMassBoro$sigsq, NA)
  
  bmPC1Boro = fitContinuous(phyPC1, PC1Boro, model="BM", 
                            control=list(niter=10))$opt
  bmModelResPC1Boro[i, ] = c(tree, bmPC1Boro$lnL, bmPC1Boro$aicc, bmPC1Boro$sigsq, 
                             NA)
  
  # ACDC 
  acdcMassBoro = fitContinuous(phyMass, massBoro, model="EB", 
                               bounds=list(a=c(-0.2, 0.2)), 
                               control=list(niter=10))$opt
  acdcModelResMassBoro[i, ] = c(tree, acdcMassBoro$lnL, acdcMassBoro$aicc, 
                                acdcMassBoro$sigsq, acdcMassBoro$a)
  
  acdcPC1Boro = fitContinuous(phyPC1, PC1Boro, model="EB", 
                              bounds=list(a=c(-0.2, 0.2)), 
                              control=list(niter=10))$opt
  acdcModelResPC1Boro[i, ] = c(tree, acdcPC1Boro$lnL, acdcPC1Boro$aicc, 
                               acdcPC1Boro$sigsq, acdcPC1Boro$a)
  
  # Trend
  trendMassBoro = fitContinuous(phyMass, massBoro, model="trend", 
                                control=list(niter=10))$opt
  trendResMassBoro[i, ] = c(tree, trendMassBoro$lnL, trendMassBoro$aicc, 
                            trendMassBoro$sigsq, trendMassBoro$slope)
  
  trendPC1Boro = fitContinuous(phyPC1, PC1Boro, model="trend", 
                               control=list(niter=10))$opt
  trendResPC1Boro[i, ] = c(tree, trendPC1Boro$lnL, trendPC1Boro$aicc, 
                           trendPC1Boro$sigsq, trendMassBoro$slope)
  
  # Drift
  driftMassBoro = fitContinuous(phyMass, massBoro, model="drift", 
                                control=list(niter=10))$opt
  driftResMassBoro[i, ] = c(tree, driftMassBoro$lnL, driftMassBoro$aicc, 
                            driftMassBoro$sigsq, driftMassBoro$drift)
  
  driftPC1Boro = fitContinuous(phyPC1, PC1Boro, model="drift", 
                               control=list(niter=10))$opt
  driftResPC1Boro[i, ] = c(tree, driftPC1Boro$lnL, driftPC1Boro$aicc, 
                           driftPC1Boro$sigsq, driftPC1Boro$drift)
  
  # diversity NOT by diet
  divMassBoro = fitDiversityModel(phyMass, massBoro, showTree=F)
  divMassBoro$aicc = aic.comp(divMassBoro$logL, k=3, n=nMassBoro, small.sample=T)
  divModelResMassBoro[i, ] = c(tree, divMassBoro$logL, divMassBoro$aicc, 
                               divMassBoro$sig0, divMassBoro$psi)
  
  divPC1Boro = fitDiversityModel(phyPC1, PC1Boro, showTree=F)
  divPC1Boro$aicc = aic.comp(divPC1Boro$logL, k=3, n=nPC1Boro, small.sample=T)
  divModelResPC1Boro[i, ] = c(tree, divPC1Boro$logL, divPC1Boro$aicc, 
                              divPC1Boro$sig0, divPC1Boro$psi)
  
  # OU
  ouMassBoro = fitContinuous(phyMass, massBoro, model="OU", 
                             control=list(niter=10))$opt
  ouModelResMassBoro[i, ] = c(tree, ouMassBoro$lnL, ouMassBoro$aicc, 
                              ouMassBoro$sigsq, ouMassBoro$alpha)
  
  ouPC1Boro = fitContinuous(phyPC1, PC1Boro, model="OU", 
                            control=list(niter=10))$opt
  ouModelResPC1Boro[i, ] = c(tree, ouPC1Boro$lnL, ouPC1Boro$aicc, ouPC1Boro$sigsq, 
                             ouMassBoro$alpha)
  
  aicMassBoro = setNames(c(bmMassBoro$aicc, acdcMassBoro$aicc, trendMassBoro$aicc, 
                           driftMassBoro$aicc, divMassBoro$aicc, ouMassBoro$aicc), 
                         c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsMassBoro[i, ] = aicw(aicMassBoro)[,3]
  
  aicPC1Boro = setNames(c(bmPC1Boro$aicc, acdcPC1Boro$aicc, trendPC1Boro$aicc, 
                          driftPC1Boro$aicc, divPC1Boro$aicc, ouPC1Boro$aicc), 
                        c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsPC1Boro[i, ] = aicw(aicPC1Boro)[,3]
  
  if(i%%10 == 0) {
    print(paste("Samples complete =", i))
  }	
}



# CANINAE

samples = read.csv("Slater2015_data/posteriorSample.csv")[, 1]
# same 500 samples as above

for(i in 1:500) {
  
  tree = samples[i]
  phy.tmp = Cani.t.file1[[tree]]
  phy.tmp$edge.length = phy.tmp$edge.length / p.file1[tree, "Clockrate"]
  
  # prep trees
  phyMass = treedata(phy.tmp, massCani)$phy
  phyPC1 = treedata(phy.tmp, PC1Cani)$phy
  
  # fit models
  # first mass, then PC1
  
  # Bm 
  bmMassCani = fitContinuous(phyMass, massCani, model="BM", 
                             control=list(niter=10))$opt
  bmModelResMassCani[i, ] = c(tree, bmMassCani$lnL, bmMassCani$aicc, 
                              bmMassCani$sigsq, NA)
  
  bmPC1Cani = fitContinuous(phyPC1, PC1Cani, model="BM", 
                            control=list(niter=10))$opt
  bmModelResPC1Cani[i, ] = c(tree, bmPC1Cani$lnL, bmPC1Cani$aicc, bmPC1Cani$sigsq, 
                             NA)
  
  # ACDC 
  acdcMassCani = fitContinuous(phyMass, massCani, model="EB", 
                               bounds=list(a=c(-0.2, 0.2)), 
                               control=list(niter=10))$opt
  acdcModelResMassCani[i, ] = c(tree, acdcMassCani$lnL, acdcMassCani$aicc, 
                                acdcMassCani$sigsq, acdcMassCani$a)
  
  acdcPC1Cani = fitContinuous(phyPC1, PC1Cani, model="EB", 
                              bounds=list(a=c(-0.2, 0.2)), 
                              control=list(niter=10))$opt
  acdcModelResPC1Cani[i, ] = c(tree, acdcPC1Cani$lnL, acdcPC1Cani$aicc, 
                               acdcPC1Cani$sigsq, acdcPC1Cani$a)
  
  # Trend
  trendMassCani = fitContinuous(phyMass, massCani, model="trend", 
                                control=list(niter=10))$opt
  trendResMassCani[i, ] = c(tree, trendMassCani$lnL, trendMassCani$aicc, 
                            trendMassCani$sigsq, trendMassCani$slope)
  
  trendPC1Cani = fitContinuous(phyPC1, PC1Cani, model="trend", 
                               control=list(niter=10))$opt
  trendResPC1Cani[i, ] = c(tree, trendPC1Cani$lnL, trendPC1Cani$aicc, 
                           trendPC1Cani$sigsq, trendMassCani$slope)
  
  # Drift
  driftMassCani = fitContinuous(phyMass, massCani, model="drift", 
                                control=list(niter=10))$opt
  driftResMassCani[i, ] = c(tree, driftMassCani$lnL, driftMassCani$aicc, 
                            driftMassCani$sigsq, driftMassCani$drift)
  
  driftPC1Cani = fitContinuous(phyPC1, PC1Cani, model="drift", 
                               control=list(niter=10))$opt
  driftResPC1Cani[i, ] = c(tree, driftPC1Cani$lnL, driftPC1Cani$aicc, 
                           driftPC1Cani$sigsq, driftPC1Cani$drift)
  
  # diversity NOT by diet
  divMassCani = fitDiversityModel(phyMass, massCani, showTree=F)
  divMassCani$aicc = aic.comp(divMassCani$logL, k=3, n=nMassCani, small.sample=T)
  divModelResMassCani[i, ] = c(tree, divMassCani$logL, divMassCani$aicc, 
                               divMassCani$sig0, divMassCani$psi)
  
  divPC1Cani = fitDiversityModel(phyPC1, PC1Cani, showTree=F)
  divPC1Cani$aicc = aic.comp(divPC1Cani$logL, k=3, n=nPC1Cani, small.sample=T)
  divModelResPC1Cani[i, ] = c(tree, divPC1Cani$logL, divPC1Cani$aicc, 
                              divPC1Cani$sig0, divPC1Cani$psi)
  
  # OU
  ouMassCani = fitContinuous(phyMass, massCani, model="OU", 
                             control=list(niter=10))$opt
  ouModelResMassCani[i, ] = c(tree, ouMassCani$lnL, ouMassCani$aicc, 
                              ouMassCani$sigsq, ouMassCani$alpha)
  
  ouPC1Cani = fitContinuous(phyPC1, PC1Cani, model="OU", 
                            control=list(niter=10))$opt
  ouModelResPC1Cani[i, ] = c(tree, ouPC1Cani$lnL, ouPC1Cani$aicc, ouPC1Cani$sigsq, 
                             ouMassCani$alpha)
  
  aicMassCani = setNames(c(bmMassCani$aicc, acdcMassCani$aicc, trendMassCani$aicc, 
                           driftMassCani$aicc, divMassCani$aicc, ouMassCani$aicc), 
                         c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsMassCani[i, ] = aicw(aicMassCani)[,3]
  
  aicPC1Cani = setNames(c(bmPC1Cani$aicc, acdcPC1Cani$aicc, trendPC1Cani$aicc, 
                          driftPC1Cani$aicc, divPC1Cani$aicc, ouPC1Cani$aicc), 
                        c("bm", "acdc", "trend", "drift", "div", "ou"))
  weightsPC1Cani[i, ] = aicw(aicPC1Cani)[,3]
  
  if(i%%10 == 0) {
    print(paste("Samples complete =", i))
  }	
}



# process posterior sample output


# function to calculate median values
ff = function(x)
  return(apply(x, 2, median, na.rm=T))


# mass

resMass = list(bm=bmModelResMass, acdc=acdcModelResMass, trend=trendResMass, 
               drift=driftResMass, div=divModelResMass, ou=ouModelResMass)
resMassHesp = list(bm=bmModelResMassHesp, acdc=acdcModelResMassHesp, 
                   trend=trendResMassHesp, drift=driftResMassHesp, 
                   div=divModelResMassHesp, ou=ouModelResMassHesp)
resMassBoro = list(bm=bmModelResMassBoro, acdc=acdcModelResMassBoro, 
                   trend=trendResMassBoro, drift=driftResMassBoro, 
                   div=divModelResMassBoro, ou=ouModelResMassBoro)
resMassCani = list(bm=bmModelResMassCani, acdc=acdcModelResMassCani, 
                   trend=trendResMassCani, drift=driftResMassCani, 
                   div=divModelResMassCani, ou=ouModelResMassCani)

massRes = matrix(unlist(lapply(resMass, ff)), nrow=length(resMass), ncol=5, 
                 byrow=T)[, -1]
massResHesp = matrix(unlist(lapply(resMassHesp, ff)), nrow=length(resMassHesp), 
                     ncol=5, byrow=T)[, -1]
massResBoro = matrix(unlist(lapply(resMassBoro, ff)), nrow=length(resMassBoro), 
                     ncol=5, byrow=T)[, -1]
massResCani = matrix(unlist(lapply(resMassCani, ff)), nrow=length(resMassCani), 
                     ncol=5, byrow=T)[, -1]

colnames(massRes) = colnames(massResHesp) = colnames(massResBoro) = 
  colnames(massResCani) = c("loglk", "AICc", "sigmasq", "PARAM")
rownames(massRes) = rownames(massResHesp) = rownames(massResBoro) = 
  rownames(massResCani) = names(resMass)

avgWeightsMass = apply(weightsMass, 2, median)
avgWeightsMassHesp = apply(weightsMassHesp, 2, median)
avgWeightsMassBoro = apply(weightsMassBoro, 2, median)
avgWeightsMassCani = apply(weightsMassCani, 2, median)


# repeat for PC1

resPC1 = list(bm=bmModelResPC1, acdc=acdcModelResPC1, trend=trendResPC1, 
              drift=driftResPC1, div=divModelResPC1, ou=ouModelResPC1)
resPC1Hesp = list(bm=bmModelResPC1Hesp, acdc=acdcModelResPC1Hesp, trend=trendResPC1Hesp, 
              drift=driftResPC1Hesp, div=divModelResPC1Hesp, ou=ouModelResPC1Hesp)
resPC1Boro = list(bm=bmModelResPC1Boro, acdc=acdcModelResPC1Boro, trend=trendResPC1Boro, 
              drift=driftResPC1Boro, div=divModelResPC1Boro, ou=ouModelResPC1Boro)
resPC1Cani = list(bm=bmModelResPC1Cani, acdc=acdcModelResPC1Cani, trend=trendResPC1Cani, 
              drift=driftResPC1Cani, div=divModelResPC1Cani, ou=ouModelResPC1Cani)

PC1res = matrix(unlist(lapply(resPC1, ff)), nrow=length(resPC1), ncol=5, 
                byrow=T)[, -1]
PC1resHesp = matrix(unlist(lapply(resPC1Hesp, ff)), nrow=length(resPC1), ncol=5, 
                    byrow=T)[, -1]
PC1resBoro = matrix(unlist(lapply(resPC1Boro, ff)), nrow=length(resPC1), ncol=5, 
                    byrow=T)[, -1]
PC1resCani = matrix(unlist(lapply(resPC1Cani, ff)), nrow=length(resPC1), ncol=5, 
                    byrow=T)[, -1]

colnames(PC1res) = colnames(PC1resHesp) = colnames(PC1resBoro) = colnames(PC1resCani) = 
  c("loglk", "AICc", "sigmasq", "PARAM")
rownames(PC1res) = rownames(PC1resHesp) = rownames(PC1resBoro) = rownames(PC1resCani) = 
  names(resPC1)

avgWeightsPC1 = apply(weightsPC1, 2, median)
avgWeightsPC1Hesp = apply(weightsPC1Hesp, 2, median)
avgWeightsPC1Boro = apply(weightsPC1Boro, 2, median)
avgWeightsPC1Cani = apply(weightsPC1Cani, 2, median)



################################################################################
# FIGURE 5
# plot Akaike weights

par(oma=c(2, 4, 0, 0), mar=c(3.5, 3, 5, 0), mfrow=c(2, 4))

barplot(avgWeightsMass, ylim=c(0, 1), ylab="", cex.axis=1.5, cex.lab=1.5, 
        cex.names=1.5, yaxt="n")
title("1) all canids", adj=0, cex.main=1.5, line=1)
mtext(text="A. Body mass", adj=0, line=3, cex=1.5)
barplot(avgWeightsMassHesp, ylim=c(0, 1), ylab="", cex.axis=1.5, cex.lab=1.5, 
        cex.names=1.5, yaxt="n")
title("2) Hesperocyoninae", adj=0, cex.main=1.5, line=1)
barplot(avgWeightsMassBoro, ylim=c(0, 1), ylab="", cex.axis=1.5, cex.lab=1.5, 
        cex.names=1.5, yaxt="n")
title("3) Borophaginae", adj=0, cex.main=1.5, line=1)
barplot(avgWeightsMassCani, ylim=c(0, 1), ylab="", cex.axis=1.5, cex.lab=1.5, 
        cex.names=1.5, yaxt="n")
title("4) fossil Caninae", adj=0, cex.main=1.5, line=1)

axis(side=2, outer=T, line=-1.5, cex.axis=1.5)
mtext(text="Akaike weight", side=2, line=2, outer=TRUE, cex=1.5)

barplot(avgWeightsPC1, ylim=c(0, 1), ylab="", cex.axis=1.5, 
        cex.lab=1.5, cex.names=1.5, yaxt="n")
title("1) all canids", adj=0, cex.main=1.5, line=1)
mtext(text="B. Carnivory", adj=0, line=3, cex=1.5)
barplot(avgWeightsPC1Hesp, ylim=c(0, 1), ylab="", cex.axis=1.5, 
        cex.lab=1.5, cex.names=1.5, yaxt="n")
title("2) Hesperocyoninae", adj=0, cex.main=1.5, line=1)
barplot(avgWeightsPC1Boro, ylim=c(0, 1), ylab="", cex.axis=1.5, 
        cex.lab=1.5, cex.names=1.5, yaxt="n")
title("3) Borophaginae", adj=0, cex.main=1.5, line=1)
barplot(avgWeightsPC1Cani, ylim=c(0, 1), ylab="", cex.axis=1.5, 
        cex.lab=1.5, cex.names=1.5, yaxt="n")
title("4) fossil Caninae", adj=0, cex.main=1.5, line=1)

axis(side=2, outer=T, line=-1.5, cex.axis=1.5)

mtext(text="model", side=1, line=0, outer=TRUE, cex=1.5)
################################################################################