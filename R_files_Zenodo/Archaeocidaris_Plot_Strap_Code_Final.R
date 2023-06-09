#The subsequent lines of code replicate analyses of stratigraphic fit of archaeocidarid phylogenies from Thompson et al. 20__
#The packages Ape, Strap, and Paleotree should be installed. Directory names are also specific to my computer, and should be updated as appropriate. This code
#was directly inspired by the tutorial for STRAP, by Mark Bell and Graham Lloyd.

#Plot The Strict Conensus from Parsimony Analysis Against Stratigraphy
Archaeocidaris_Tree<-read.nexus("/Volumes/Seagate_Backup_Plus_Drive/Archaeocidarid_Phylogeny/Revision/Updated_Figures/Reweighted_Contree_Parsimony.tre")
Archaeocidaris_Ages<-read.table("/Volumes/Seagate_Backup_Plus_Drive/Archaeocidarid_Phylogeny/Revision/Manuscript_Revisions/For_Dryad/Archaeocidaris_Ages.txt")

Archaeocidaris_TS<-DatePhylo(Archaeocidaris_Tree, Archaeocidaris_Ages, 2, "equal")

geoscalePhylo(Archaeocidaris_TS, Archaeocidaris_Ages, 
              direction="rightwards", units=c("Period", "Age"), vers="ICS2012", cex.tip=.5, cex.ts=.5)

#Examine how these parsimony trees hold up against stratigraphy
Archaeocidaris_All_Trees<-read.nexus("/Volumes/Seagate_Backup_Plus_Drive/Archaeocidarid_Phylogeny/Revision/Manuscript_Revisions/For_Dryad/MPTs_Reweighted.tre")

Strat_Fit_Archaeocidaris<-StratPhyloCongruence(Archaeocidaris_All_Trees,Archaeocidaris_Ages, 
                                               hard=TRUE, randomly.sample.ages=TRUE, fix.topology=TRUE, fix.outgroup=TRUE)

write.table(Strat_Fit_Archaeocidaris$input.tree.results, 
            "/Volumes/Seagate_Backup_Plus_Drive/Archaeocidarid_Phylogeny/Revision/Manuscript_Revisions/For_Dryad/Strat_Fit_Archaeocidaris.csv", sep=",")

#SCI Plots and Analysis
Alpha_Value_SCI<- qnorm(0.95, mean(Strat_Fit_Archaeocidaris$rand.permutations[, "SCI"]), sd(Strat_Fit_Archaeocidaris$rand.permutations[, "SCI"]))
Alpha_Value_GER<- qnorm(0.95, mean(Strat_Fit_Archaeocidaris$rand.permutations[, "GER"]), sd(Strat_Fit_Archaeocidaris$rand.permutations[, "GER"]))

pdf("/Volumes/Seagate_Backup_Plus_Drive/Archaeocidarid_Phylogeny/Revision/Manuscript_Revisions/For_Dryad/Archaeocidairs_Strat_Metrics_Plot.pdf")
par(mfrow=c(1,2))

Hist<-hist(Strat_Fit_Archaeocidaris$rand.permutations[, "SCI"], xlim=c(0, .6), xlab="SCI",
           main=NULL, col=rgb(1,0,0,0.5))
hist(Strat_Fit_Archaeocidaris$input.tree.results[, "SCI"], xlim=c(0, .6), xlab="SCI", col=rgb(0,0,1,0.5), add=T)
lines(x=c(Alpha_Value_SCI, Alpha_Value_SCI), y=c(0, 1000) , lty=2)

Hist<-hist(Strat_Fit_Archaeocidaris$rand.permutations[, "GER"], xlim=c(0, 1), xlab="GER",
           main=NULL, col=rgb(1,0,0,0.5))
hist(Strat_Fit_Archaeocidaris$input.tree.results[, "GER"], xlim=c(0, 1), xlab="GER", col=rgb(0,0,1,0.5), add=T, breaks=c(.05, .1, 1.5, .2, 2.5, .3, 3.5, .4, 4.5, .5, 5.5, .6, 6.5, .7, 7.5, .8, 8.5, .9, 9.5))
lines(x=c(Alpha_Value_GER, Alpha_Value_GER), y=c(0, 1000) , lty=2)

dev.off()


