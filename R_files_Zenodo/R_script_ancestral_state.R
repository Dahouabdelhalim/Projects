library(ape); library(caper); library(geiger); library(phytools); library(phangorn); library(xlsx)
citation("ape"); citation("caper"); citation("geiger"); citation("phytools"); citation("phangorn"); citation('xlsx')

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\\n\\n")
  tree
}

ReorderData <- function(tree, data, taxa.names="row names"){
  new.data <- data
  if(is.vector(data)){
    for(i in 1:length(tree$tip.label)){
      new.data[i] <- data[which(names(data) == tree$tip.label[i])]
      names(new.data)[i] <- names(data)[which(names(data) == tree$tip.label[i])]
    }
  }
  if(is.data.frame(data) || is.matrix(data)){
    if(taxa.names == "row names"){
      row.names(new.data) <- 1:length(tree$tip.label)
      for(i in 1:length(tree$tip.label)){
        new.data[i,] <- data[which(row.names(data) == tree$tip.label[i]),]
        row.names(new.data)[i] <- row.names(data)[which(row.names(data) == tree$tip.label[i])]
      }
    }
    if(is.numeric(taxa.names)){
      for(i in 1:length(tree$tip.label)){
        new.data[i,] <- data[which(data[,taxa.names] == tree$tip.label[i]),]
      }
    }
  }
  return(new.data)  
}





####  ARTHROPODS ####

Ancestral_state_tree <- read.tree("arthro1.tre") # just using a small one for ease
plot(Ancestral_state_tree)

is.binary.tree(Ancestral_state_tree) #Test is the tree is truly dichotomous

#This line transforms the branches of our tree with 0 in length to be 1e-6 of the maximum of the tree length (to avoid errors later)
Ancestral_state_tree$edge.length[Ancestral_state_tree$edge.length == 0] <- max(nodeHeights(Ancestral_state_tree))*1e-6
plot(Ancestral_state_tree,cex=0.5)
Ancestral_state_tree <- force.ultrametric(Ancestral_state_tree)
is.ultrametric(Ancestral_state_tree) #TRUE
plot(Ancestral_state_tree,cex=0.5)


#Import the homemade data set with the 2011 tree tips and the corresponding values for regen and autotomy based on cladogram I mapped
#from Regier et al. 2010
Regen2 <- read.table("Regen_Data.txt")
Regen2$Real_Regen <- as.factor(Regen2$Real_Regen)
Regen2$Real_Autotomy <- as.factor(Regen2$Real_Autotomy)
summary(Regen2)
#Lots of missing data so two ways about it, either get rid of the missing data or consider them as uncertain states (0.5)

#Getting rid of the missing data and matching the tree to the new data set
Regen2_naexclude <- na.exclude(Regen2); summary(Regen2_naexclude)
Ancestral_state_tree_naexclude <- drop.tip(Ancestral_state_tree, setdiff(Ancestral_state_tree$tip.label, Regen2_naexclude$animal))
Ancestral_state_tree_naexclude$edge.length[Ancestral_state_tree_naexclude$edge.length == 0] <-
  max(nodeHeights(Ancestral_state_tree_naexclude))*1e-6
Ancestral_state_tree_naexclude <- force.ultrametric(Ancestral_state_tree_naexclude)
Ancestral_state_tree_naexclude$edge.length[Ancestral_state_tree_naexclude$edge.length == 0] <-
  max(nodeHeights(Ancestral_state_tree_naexclude))*1e-6 #I know, we do this twice, but that's what it takes
Ancestral_state_tree_naexclude <- force.ultrametric(Ancestral_state_tree_naexclude)
is.ultrametric(Ancestral_state_tree_naexclude) #TRUE
is.binary.tree(Ancestral_state_tree_naexclude) #TRUE OK still good
plot(Ancestral_state_tree_naexclude,cex=0.5)

#Now let's change those tip names to be more elegant
Regen <- read.table("Regen_Data_2.txt")
Regen$Real_Regen <- as.factor(Regen$Real_Regen)
Regen$Real_Autotomy <- as.factor(Regen$Real_Autotomy)
summary(Regen)
Regen_naexclude <- na.exclude(Regen); summary(Regen_naexclude)
Ancestral_state_tree_naexclude$tip.label <- as.character(Regen_naexclude$animal)
plot(Ancestral_state_tree_naexclude,cex=0.5)

#Let's make sure the data set is ordered as the tip labels of the tree
Regen_naexclude <- ReorderData(Ancestral_state_tree_naexclude, Regen_naexclude, taxa.names = 1)

#These commands issue lots of output, of which you are likely to be most interested in the log likelihoods and the inferred transition rates.
ERreconstruction <- ace(Regen_naexclude$Real_Autotomy, Ancestral_state_tree_naexclude, type="discrete", model="ER")
SYMreconstruction <- ace(Regen_naexclude$Real_Autotomy, Ancestral_state_tree_naexclude, type="discrete", model="SYM")
ARDreconstruction <- ace(Regen_naexclude$Real_Autotomy, Ancestral_state_tree_naexclude, type="discrete", model="ARD")

#You can access the log likelihoods for each model as follows
ERreconstruction$loglik
SYMreconstruction$loglik
ARDreconstruction$loglik

#The transition rates are accessed with
ERreconstruction$rates
SYMreconstruction$rates
ARDreconstruction$rates

#You may have noticed that the output for the ER model matches that of the SYM model. This is not surprising. For a binary character these models are
#identical, though they are distinct for characters with three or more states. For a three-state character, ER is a one parameter model, SYM a three parameter
#model, and ARD a six parameter model.
#Which model should you prefer?
#In this example, the ARD model gives the highest likelihood. However, it also includes more parameters that the ER and SYM models, and it is well known that
#adding parameters to a model generally increases its likelihood. To determine whether the use of the more heavily parameterized model is appropriate, you
#should execute a likelihood test.
#Twice the difference in log likelihood between two models (known as the G statistic) is distributed as chi square, with degrees of freedom equal to the
#number of parameters that differ between the models. Thus a likelihood test of these models asks whether the difference in likelihoods is large enough to
#lie in the rightmost tail of the chisquare distribution, typically considered to be the largest 5% of values.
#The command pchisq(value, df) gives the percentage of the cumulative distribution function for chisquare that lies to the left of the given value for a
#desired degree of freedom(df). Thus, the following command gives the percentage of the chisquare distribution that lies to the right of the observed
#likelihood difference for the ER and ARD models (which differ by one parameter), given the Geospiza data.
1-pchisq(2*abs(ERreconstruction$loglik - ARDreconstruction$loglik), 1)

#Evaluation of the equation above yields the p-value for the likelihood test, in this case significant. This is a significant result at a 5% error rate
#, so we do accept the less heavily parameterized ER model.

#Are there any caveats about the reconstruction of discrete characters that I should be aware of?
#Yes. Ace deals with this analysis by parameterizing each node within the phylogeny and reconstructing its ancestral state separately. This is a
#computationally intensive problem, and ace sometime returns spurious results when attempting a solution. If the returned log likelihood is tiny or enormous,
#or some of the reconstructed ancestral states have greater than 100% probability (or negative probability) at certain nodes, it is likely that ape is not
#reconstructing an accurate solution for the character in question. You are most likely to run into this issue when using complicated transition models or
#trees that imply many character state changes.
#You can see the probability of each of the possible ancestral states at each of the nodes by typing the following:
head(ERreconstruction$lik.anc)

#In other cases, the likelihood surface for the rates will be essentially flat, and you can place little confidence in the specific rates being reconstructed
#for each node. If the rate values are important in your application, check their standard errors.
ERreconstruction$se

#If the standard errors are very large, or represented by NaN (not a number), then you are dealing with a flat or nearly flat likelihood surface, and should
#not place any confidence in the reconstructed rates. It's not the case here.

#Good old classic representation
plotTree(Ancestral_state_tree_naexclude, type = "fan", fsize = 0.8, ftype = "i")
cols <- setNames(palette()[1:length(unique(Regen_naexclude$Real_Autotomy))], sort(unique(Regen_naexclude$Real_Autotomy)))
nodelabels(node = 1:Ancestral_state_tree_naexclude$Nnode + Ntip(Ancestral_state_tree_naexclude),
           pie = ERreconstruction$lik.anc, piecol = cols, cex = 0.5)
tiplabels(pie = to.matrix(Regen_naexclude$Real_Regen, sort(unique(Regen_naexclude$Real_Regen))), piecol = "green", cex = 0.1)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1],
                  y = -max(nodeHeights(Ancestral_state_tree_naexclude)), fsize = 0.8)

#An alternative approach to the one outline above is to use an MCMC approach to sample character histories from their posterior probability
#distribution. This is called stochastic character mapping (Huelsenbeck et al. 2003). The model is the same but in this case we get a sample of
#unambiguous histories for our discrete character's evolution on the tree - rather than a probability distribution for the character at nodes.
#For instance, given the data simulated above - we can generate the stochastic character map as follows:
# simulate single stochastic character map using empirical Bayes method

#Trying something more classy
M <- Regen_naexclude$Real_Autotomy; names(M) <- Regen_naexclude$animal
mtree<-make.simmap(Ancestral_state_tree_naexclude,M,model="ER")
mtree
plot(mtree,cols,type="fan",fsize=0.8,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(Ancestral_state_tree_naexclude)),fsize=0.8)

#A single stochastic character map does not mean a whole lot in isolation - we need to look at the whole distribution from a sample of
#stochastic maps. This can be a bit overwhelming. For instance, the following code generates 100 stochastic character maps from our dataset
#and plots them in a grid:
mtrees<-make.simmap(Ancestral_state_tree_naexclude,M,model="ER",nsim=100)
mtrees
par(mfrow=c(10,10)); null<-sapply(mtrees,plot,colors=cols,lwd=1,ftype="off")

#It's possible to summarize a set of stochastic maps in a much more meaningful way. For instance, we can estimate the number of changes of
#each type, the proportion of time spent in each state, and the posterior probabilities that each internal node is in each state, under our
#model. For example:
pd<-summary(mtrees,plot=FALSE); pd
dev.off(); plot(pd,fsize=0.6,ftype="i")

## now let's plot a random map, and overlay the posterior probabilities
plot(mtrees[[1]],cols,type="fan",fsize=0.8,ftype="i")
##SetEnv=TRUE for this type is experimental. please be patient with bugs
nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(mtrees)),fsize=0.8)

#Finally, since we obtained these inferences under exactly the same model, let's compare the posterior probabilities from stochastic
#mapping with our marginal ancestral states. In the former case, our probabilities were obtained by sampling from the joint (rather than
#marginal) probability distribution for the ancestral states.
plot(ERreconstruction$lik.anc,pd$ace,xlab="marginal ancestral states", ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)

obj <- densityMap(mtrees, res = 1000, lwd = 3, outline = FALSE, type = "fan")
tiplabels(pie = to.matrix(Regen_naexclude$Real_Autotomy, sort(unique(Regen_naexclude$Real_Autotomy))), piecol = cols, cex = 0.13)




#### REPTILES ####
Ancestral_state_tree2 <- read.nexus("GEA_Combined_MkA_Calibrated.tre") #Figure 6 from 'Pyron, R. (2017) Systematic Biology
plot(Ancestral_state_tree2)

is.binary.tree(Ancestral_state_tree2) #Test is the tree is truly dichotomous

#This line transforms the branches of our tree with 0 in length to be 1e-6 of the maximum of the tree length (to avoid errors later)
Ancestral_state_tree2$edge.length[Ancestral_state_tree2$edge.length == 0] <- max(nodeHeights(Ancestral_state_tree2))*1e-6
plot(Ancestral_state_tree2,cex=0.5)
Ancestral_state_tree2 <- force.ultrametric(Ancestral_state_tree2)
is.ultrametric(Ancestral_state_tree2) #TRUE
plot(Ancestral_state_tree2,cex=0.5)

#Extract the tip labels to fill in with the regeneration/autotomy information I've got
write.table(Ancestral_state_tree2$tip.label, "reptile_data.txt", sep="\\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

  ############################################################

#Middle step here where I import back the file that I just created, but with the regen/autotomy data added
#Additionally, taxa without information were dropped
Reptile_data <- read.table("reptile_data.txt", header = TRUE)
Reptile_data$Regen <- as.factor(Reptile_data$Regen)
Reptile_data$Autotomy <- as.factor(Reptile_data$Autotomy)
summary(Reptile_data)

  ############################################################

#Let's get rid off the NAs
Reptile_data <- na.exclude(Reptile_data)
summary(Reptile_data)

#Let's make sure our tree correspond to our data set
Ancestral_state_tree2_naexclude <- drop.tip(Ancestral_state_tree2, setdiff(Ancestral_state_tree2$tip.label, Reptile_data$animal))
Ancestral_state_tree2_naexclude$edge.length[Ancestral_state_tree2_naexclude$edge.length == 0] <- max(nodeHeights(Ancestral_state_tree2_naexclude))*1e-6
is.binary.tree(Ancestral_state_tree2_naexclude) #Test is the tree is truly dichotomous
is.ultrametric(Ancestral_state_tree2_naexclude) #TRUE
Ancestral_state_tree2_naexclude <- force.ultrametric(Ancestral_state_tree2_naexclude)
is.ultrametric(Ancestral_state_tree2_naexclude) #TRUE
plot(Ancestral_state_tree2_naexclude,cex=0.5)

#Let's make sure the data set is ordered as the tip labels of the tree
Reptile_data <- ReorderData(Ancestral_state_tree2_naexclude, Reptile_data, taxa.names = 1)





### AUTOTOMY ###
#These commands issue lots of output, of which you are likely to be most interested in the log likelihoods and the inferred transition rates.
ERreconstruction2.5 <- ace(Reptile_data$Autotomy, Ancestral_state_tree2_naexclude, type="discrete", model="ER")
SYMreconstruction2.5 <- ace(Reptile_data$Autotomy, Ancestral_state_tree2_naexclude, type="discrete", model="SYM")
ARDreconstruction2.5 <- ace(Reptile_data$Autotomy, Ancestral_state_tree2_naexclude, type="discrete", model="ARD")

#You can access the log likelihoods for each model as follows
ERreconstruction2.5$loglik
SYMreconstruction2.5$loglik
ARDreconstruction2.5$loglik

#The transition rates are accessed with
ERreconstruction2.5$rates
SYMreconstruction2.5$rates
ARDreconstruction2.5$rates

1-pchisq(2*abs(ERreconstruction2.5$loglik - ARDreconstruction2.5$loglik), 1)

head(ERreconstruction2.5$lik.anc)

ERreconstruction2.5$se

#Good old classic representation
plotTree(Ancestral_state_tree2_naexclude, type = "phylogram", fsize = 0.8, ftype = "i")
cols <- setNames(palette()[1:length(unique(Reptile_data$Autotomy))], sort(unique(Reptile_data$Autotomy)))
nodelabels(node = 1:Ancestral_state_tree2_naexclude$Nnode + Ntip(Ancestral_state_tree2_naexclude),
           pie = ERreconstruction2.5$lik.anc, piecol = cols, cex = 0.4)
tiplabels(pie = to.matrix(Reptile_data$Autotomy, sort(unique(Reptile_data$Autotomy))), piecol = cols, cex = 0.3)

#Trying something more classy
M2.5 <- Reptile_data$Autotomy; names(M2.5) <- Reptile_data$animal
mtree2.5 <- make.simmap(Ancestral_state_tree2_naexclude, M2.5, model = "ER")
mtree2.5
plot(mtree2.5, cols, type = "fan", fsize = 0.8, ftype = "i")
##SetEnv=TRUE for this type is experimental. please be patient with bugs
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(Ancestral_state_tree2_naexclude)), fsize = 0.78)

#A single stochastic character map does not mean a whole lot in isolation - we need to look at the whole distribution from a sample of
#stochastic maps. This can be a bit overwhelming. For instance, the following code generates 100 stochastic character maps from our dataset
#and plots them in a grid:
mtrees2.5 <- make.simmap(Ancestral_state_tree2_naexclude, M2.5, model = "ER", nsim = 100)
mtrees2.5
par(mfrow = c(10, 10)); null <- sapply(mtrees2.5, plot, colors = cols, lwd = 1, ftype = "off")

#It's possible to summarize a set of stochastic maps in a much more meaningful way. For instance, we can estimate the number of changes of
#each type, the proportion of time spent in each state, and the posterior probabilities that each internal node is in each state, under our
#model. For example:
pd2.5 <- summary(mtrees2.5, plot = FALSE); pd2.5
dev.off(); plot(pd2.5, fsize = 0.6, ftype = "i", cex = 0.3)

## now let's plot a random map, and overlay the posterior probabilities
plot(mtrees2.5[[1]], cols, type = "fan", fsize = 0.8, ftype = "i")
## setEnv=TRUE for this type is experimental. please be patient with bugs
nodelabels(pie = pd2.5$ace, piecol = cols, cex = 0.5)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(Ancestral_state_tree2_naexclude)), fsize = 0.8)

#Finally, since we obtained these inferences under exactly the same model, let's compare the posterior probabilities from stochastic
#mapping with our marginal ancestral states. In the former case, our probabilities were obtained by sampling from the joint (rather than
#marginal) probability distribution for the ancestral states.
plot(ERreconstruction2.5$lik.anc, pd2.5$ace, xlab = "marginal ancestral states", ylab = "posterior probabilities from stochastic mapping")
lines(c(0, 1), c(0, 1), lty = "dashed", col = "red", lwd = 2)

obj2 <- densityMap(mtrees2.5, lwd = 3, res = 1000, outline = TRUE, type = "fan")
tiplabels(pie = to.matrix(Reptile_data$Autotomy, sort(unique(Reptile_data$Autotomy))), piecol = cols, cex = 0.1)


  ### REGENERATION ###
ERreconstruction3.5 <- ace(Reptile_data$Regen, Ancestral_state_tree2_naexclude, type="discrete", model="ER")
SYMreconstruction3.5 <- ace(Reptile_data$Regen, Ancestral_state_tree2_naexclude, type="discrete", model="SYM")
ARDreconstruction3.5 <- ace(Reptile_data$Regen, Ancestral_state_tree2_naexclude, type="discrete", model="ARD")

ERreconstruction3.5$loglik
SYMreconstruction3.5$loglik
ARDreconstruction3.5$loglik

ERreconstruction3.5$rates
SYMreconstruction3.5$rates
ARDreconstruction3.5$rates

1-pchisq(2*abs(ERreconstruction3.5$loglik - ARDreconstruction3.5$loglik), 1)

head(ERreconstruction3.5$lik.anc)

ERreconstruction3.5$se

plotTree(Ancestral_state_tree2_naexclude, type = "phylogram", fsize = 0.8, ftype = "i")
cols <- setNames(palette()[1:length(unique(Reptile_data$Regen))], sort(unique(Reptile_data$Regen)))
nodelabels(node = 1:Ancestral_state_tree2_naexclude$Nnode + Ntip(Ancestral_state_tree2_naexclude),
           pie = ERreconstruction3.5$lik.anc, piecol = cols, cex = 0.4)
tiplabels(pie = to.matrix(Reptile_data$Regen, sort(unique(Reptile_data$Regen))), piecol = cols, cex = 0.3)

#Trying something more classy
M3.5 <- Reptile_data$Regen; names(M3.5) <- Reptile_data$animal
mtree3.5 <- make.simmap(Ancestral_state_tree2_naexclude, M3.5, model = "ER")
mtree3.5
plot(mtree3.5, cols, type = "fan", fsize = 0.8, ftype = "i")
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(Ancestral_state_tree2_naexclude)), fsize = 0.78)

mtrees3.5 <- make.simmap(Ancestral_state_tree2_naexclude, M3.5, model = "ER", nsim = 100)
mtrees3.5
par(mfrow = c(10, 10)); null <- sapply(mtrees3.5, plot, colors = cols, lwd = 1, ftype = "off")

pd3.5 <- summary(mtrees3.5, plot = FALSE); pd3.5
dev.off(); plot(pd3.5, fsize = 0.6, ftype = "i", cex = 0.3)

## now let's plot a random map, and overlay the posterior probabilities
plot(mtrees3.5[[1]], cols, type = "fan", fsize = 0.8, ftype = "i")
## setEnv=TRUE for this type is experimental. please be patient with bugs
nodelabels(pie = pd3.5$ace, piecol = cols, cex = 0.5)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(Ancestral_state_tree2_naexclude)), fsize = 0.8)

#Finally, since we obtained these inferences under exactly the same model, let's compare the posterior probabilities from stochastic
#mapping with our marginal ancestral states. In the former case, our probabilities were obtained by sampling from the joint (rather than
#marginal) probability distribution for the ancestral states.
plot(ERreconstruction3.5$lik.anc, pd3.5$ace, xlab = "marginal ancestral states", ylab = "posterior probabilities from stochastic mapping")
lines(c(0, 1), c(0, 1), lty = "dashed", col = "red", lwd = 2)

obj4 <- densityMap(mtrees3.5, res = 1000, lwd = 3, outline = TRUE, type = "fan")
tiplabels(pie = to.matrix(Reptile_data$Regen, sort(unique(Reptile_data$Regen))), piecol = cols, cex = 0.1)
