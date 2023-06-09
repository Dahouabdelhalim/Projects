#Full script for replicating Winchell, Schliep, Mahler, Revell (2020) â€” Evolution

#Finalized, annotated, and checked 2/8/2020

#packages used
library(phytools)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(coda)
library(reshape)
library(MASS)
library(emmeans)
library(caper)



##########################################################################
######                                                              ######
######      (1)  Determining Probability of Being Urban             ######
######                                                              ######
##########################################################################

#read in data
scores<-read.csv("UrbanProbability_FINAL.csv", header=T)

# model: 
# We used a binomial glm where successes are urban observations and failures are the number of additional 
# observations we would need to be confident the species is urban.

X = scores[(scores$n_wild!=0) | (scores$n_urban!=0),] #exclude any that have no data
species <- X$species
failures <- ifelse((X$n_urban)>12, 0 ,(12)-( X$n_urban))  #If more than 12 urban observations probability of being urban will be 1.
fm = cbind(n_urban, failures)  ~   I(species) - 1
fit = glm(fm, X, family = "binomial") #model
res = data.frame(species = as.character(X$species), score =predict(fit, type="response"))
X$score <- round(res$score,3) #urban score, rounded to 3 digits
X$confidence <- ((X$n_urban + X$n_wild)/50)^(1/4.8)  #assign a confidence value based on number of total observations

#adjust urban scores based on the sampling uncertainty
X$urban <- ifelse((X$n_urban + X$n_wild) >= (50), #if number of total observations is at least 50
                        round(X$score,3), #model score is the estimate for urban and remainder goes to avoid
                        round((X$score * X$confidence) + 0.5*(1-X$confidence) ,3) ) #otherwise need to account for uncertainty and split remainder 
X$avoid <-  round((1 - X$urban),3)

#Set equal probability for both states for any species with no information.
NAgroup <- scores[(scores$n_wild==0) & (scores$n_urban==0),]
NAgroup$score <- NA
NAgroup$confidence <- NA
NAgroup$urban <- rep((1/2), nrow(NAgroup))
NAgroup$avoid <- rep((1/2), nrow(NAgroup))

#Combine both groups
scoresALL <- rbind(X, NAgroup)

#Reset the following to equal probability for both states because no urbanization exists on within the range of these species
scoresALL[scoresALL$species=="desechensis",c("urban", "avoid")] <- c(1/2, 1/2)
scoresALL[scoresALL$species=="ernestwilliamsi",c("urban",  "avoid")] <- c(1/2, 1/2)
scoresALL[scoresALL$species=="monensis",c("urban",  "avoid")] <- c(1/2, 1/2)
scoresALL[scoresALL$species=="nubilus",c("urban", "avoid")] <- c(1/2, 1/2)

scoresALL[,c("species", "island", "n_urban", "n_total","score","confidence", "urban", "avoid")]  # -> TipStates_Final


##########################################################################
######                                                              ######
######      (2)  Add and Drop Tips from Tree                        ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######                                                              ######
##########################################################################


#read in Gamble et al. 2014 tree
tree <- read.tree("MCC_Gamble_BEAST1.tre") #dowloaded at: https://doi.org/10.5061/dryad.dp848 
#NOTE: some species names may be spelled incorrectly

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2) 
species <- as.character(row.names(treespecies))  #vector of 131 species names

#Add missing tips
tipadds <- c("agueroi", "fairchildi", "litoralis", "roosevelti", "terraealtae")
sisteradds <- c("chamaeleonides", "smaragdinus", "argillaceus", "cuvieri", "oculatus") #specify sister species for each

tree2 <- bind.tip(tree, tipadds[1], where=which(tree$tip.label==sisteradds[1]),
             position=0.5*tree$edge.length[which(tree$edge[,2]==which(tree$tip.label==sisteradds[1]))])

tree2b <- bind.tip(tree2, tipadds[2], where=which(tree2$tip.label==sisteradds[2]), 
             position=0.5*tree2$edge.length[which(tree2$edge[,2]==which(tree2$tip.label==sisteradds[2]))])  

tree2c <- bind.tip(tree2b, tipadds[3], where=which(tree2b$tip.label==sisteradds[3]),
             position=0.5*tree2b$edge.length[which(tree2b$edge[,2]==which(tree2b$tip.label==sisteradds[3]))])

tree2d <- bind.tip(tree2c, tipadds[4], where=which(tree2c$tip.label==sisteradds[4]), 
             position=0.5*tree2c$edge.length[which(tree2c$edge[,2]==which(tree2c$tip.label==sisteradds[4]))])

tree2e <- bind.tip(tree2d, tipadds[5], where=which(tree2d$tip.label==sisteradds[5]),
             position=0.5*tree2d$edge.length[which(tree2d$edge[,2]==which(tree2d$tip.label==sisteradds[5]))])

#Drop tips from the tree that are not in our dataset
edited_tree <- drop.tip(tree2e,tree2e$tip.label[-match(species,tree2e$tip.label)])


##########################################################################
######                                                              ######
######      (3)  Threshold Model Acestral State Reconstruction      ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######                                                              ######
##########################################################################

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2) 
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states:

#read in edited tree generated in (2)
tree<- read.tree("edited_tree.tre")  


################ RUNNING ANCTHRESH #################
### WARNING: This will take approximately 2 hours per chain!

#Specs
ngen <- 10000000  #number of generations to run MCMC
sample <- 1000 #sampling interval
burnin <- .2*ngen #burnin
nspecies <- length(tree$tip.label) 
n <- nrow(states)
seq <- c("avoid", "urban") #order of state transitions

###
### RUN EACH CHAIN IN NEW R SESSION 
###

#chain 1:
set.seed(8675309)
mcmcBMchain1 <- ancThresh(tree=tree, x=states, ngen = ngen, sequence=seq, model="BM", control = list(sample = sample, burnin=burnin))

#save output !!!
#save(mcmcBMchain1, file="mcmcBMchain1.Rdata")

###

#chain 2:
set.seed(666) 
mcmcBMchain2 <- ancThresh(tree=tree, x=states, ngen = ngen, sequence=seq, model="BM", control = list(sample = sample, burnin=burnin))

#save output !!!
#save(mcmcBMchain2, file="mcmcBMchain2.Rdata")

###

#chain 3:
set.seed(1111111) 
mcmcBMchain3 <- ancThresh(tree=tree, x=states, ngen = ngen, sequence=seq, model="BM", control = list(sample = sample, burnin=burnin))

#save output !!!
#save(mcmcBMchain3, file="mcmcBMchain3.Rdata")

###

#chain 4:
set.seed(61484) 
mcmcBMchain4 <- ancThresh(tree=tree, x=states, ngen = ngen, sequence=seq, model="BM", control = list(sample = sample, burnin=burnin))

#save output !!!
#save(mcmcBMchain4, file="mcmcBMchain4.Rdata")

##########################################################


##########################################################################
######                                                              ######
######      (4)  Threshold Model Diagnostics and Combine Chains     ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######           uses "mcmcBMchain1.Rdata" generated in (3)         ######
######           uses "mcmcBMchain2.Rdata" generated in (3)         ######
######           uses "mcmcBMchain3.Rdata" generated in (3)         ######
######           uses "mcmcBMchain4.Rdata" generated in (3)         ######
######                                                              ######
##########################################################################


#read in after running
load("mcmcBMchain1.Rdata")
load("mcmcBMchain2.Rdata")
load("mcmcBMchain3.Rdata")
load("mcmcBMchain4.Rdata")


#Trim to post burnin chains (burnin 0.2)
rownum_parliab <- nrow(mcmcBMchain1$par)

mcmcBMchain1_trim <- list(par=mcmcBMchain1$par[(0.2*rownum_parliab):rownum_parliab,],liab=mcmcBMchain1$liab[(0.2*rownum_parliab):rownum_parliab,],
                     ace=mcmcBMchain1$ace, mcmc=mcmcBMchain1$mcmc[(0.2*rownum_parliab):rownum_parliab,]) 
mcmcBMchain2_trim <- list(par=mcmcBMchain2$par[(0.2*rownum_parliab):rownum_parliab,],liab=mcmcBMchain2$liab[(0.2*rownum_parliab):rownum_parliab,],
                     ace=mcmcBMchain2$ace, mcmc=mcmcBMchain2$mcmc[(0.2*rownum_parliab):rownum_parliab,])
mcmcBMchain3_trim <- list(par=mcmcBMchain3$par[(0.2*rownum_parliab):rownum_parliab,],liab=mcmcBMchain3$liab[(0.2*rownum_parliab):rownum_parliab,],
                     ace=mcmcBMchain3$ace, mcmc=mcmcBMchain3$mcmc[(0.2*rownum_parliab):rownum_parliab,])
mcmcBMchain4_trim <- list(par=mcmcBMchain4$par[(0.2*rownum_parliab):rownum_parliab,],liab=mcmcBMchain4$liab[(0.2*rownum_parliab):rownum_parliab,],
                     ace=mcmcBMchain4$ace, mcmc=mcmcBMchain4$mcmc[(0.2*rownum_parliab):rownum_parliab,])


#Geweke's convergence diagnostic - run on individual models post burn in
gewekeBMchain1 <- geweke.diag(mcmcBMchain1_trim$liab); summary(gewekeBMchain1$z)
pgew<- pnorm(abs(gewekeBMchain1$z), lower.tail=F)*2 #test significance
length(which(pgew>0.05))/length(pgew)  #fraction that converged

gewekeBMchain2 <- geweke.diag(mcmcBMchain2_trim$liab); summary(gewekeBMchain2$z)
pgew<- pnorm(abs(gewekeBMchain2$z), lower.tail=F)*2 #test significance
length(which(pgew>0.05))/length(pgew)  #fraction that converged

gewekeBMchain3 <- geweke.diag(mcmcBMchain3_trim$liab); summary(gewekeBMchain3$z)
pgew<- pnorm(abs(gewekeBMchain3$z), lower.tail=F)*2 #test significance
length(which(pgew>0.05))/length(pgew)  #fraction that converged

gewekeBMchain4 <- geweke.diag(mcmcBMchain4_trim$liab); summary(gewekeBMchain4$z)
pgew<- pnorm(abs(gewekeBMchain4$z), lower.tail=F)*2 #test significance
length(which(pgew>0.05))/length(pgew)  #fraction that converged


#Gelman convergence diagnostic - run on post-burnin chains as a list
mcmcBM.list <- mcmc.list(as.mcmc(mcmcBMchain1_trim$liab),  as.mcmc(mcmcBMchain2_trim$liab),  as.mcmc(mcmcBMchain3_trim$liab), as.mcmc(mcmcBMchain4_trim$liab))
BMgelman <- gelman.diag(mcmcBM.list)
summary(BMgelman$psrf[,2])  #convergence


#Combine post-burnin chains
mcmcBM <-list(par=rbind(mcmcBMchain1_trim$par, mcmcBMchain2_trim$par, mcmcBMchain3_trim$par, mcmcBMchain4_trim$par), 
              liab=rbind(mcmcBMchain1_trim$liab, mcmcBMchain2_trim$liab, mcmcBMchain3_trim$liab, mcmcBMchain4_trim$liab), 
              ace=(mcmcBMchain1$ace + mcmcBMchain2$ace + mcmcBMchain3$ace + mcmcBMchain4$ace)/4, #don't use post-burnin mcmc
              mcmc=rbind(mcmcBMchain1_trim$mcmc, mcmcBMchain2_trim$mcmc, mcmcBMchain3_trim$mcmc, mcmcBMchain4_trim$mcmc)) 

summary(effectiveSize(mcmcBM$liab)) #effective size of combined chains


#save output !!!
#save(mcmcBM, file="mcmcBM_combined_FINAL.Rdata")


#Posterior probabilities from combined mcmc
L<-as.matrix(mcmcBM$liab)
Th<-as.matrix(mcmcBM$par[,2:3])
STATES<-matrix(NA,nrow(L),ncol(L),dimnames=list(rownames(L),colnames(L)))
for(i in 1:nrow(L)) STATES[i,]<-threshState(L[i,],Th[i,])
foo<-function(x){ y<-summary(factor(x,levels=c("avoid","urban"))); y/sum(y) }
PP<-apply(STATES,2,foo)
PP <- t(PP)

PP_combined_FINAL <- PP

#save output !!!
#saveRDS(PP, "PP_combined_FINAL.rds")



##########################################################################
######                                                              ######
######      (5)  Fit Mk Models & Phylogenetic Signal                ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######                                                              ######
##########################################################################

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2) 
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states:

#read in edited tree generated in (2)
tree<- read.tree("edited_tree.tre")  


#Models:

#asymmetric ordered reversible
model_asymOR<-matrix(c(0,1,2,0),2,2)
rownames(model_asymOR)<-colnames(model_asymOR)<-colnames(states); model_asymOR
fit1 <- fitMk(tree, states, model=model_asymOR)

#asymmetric ordered irreversible
model_IRR <- matrix(c(0,1,0,0),2,2)
rownames(model_IRR)<-colnames(model_IRR)<-colnames(states); model_IRR
fit2 <- fitMk(tree, states, model=model_IRR)

#equal-rates ordered  reversible 
model_equalrate<-matrix(c(0,1,1,0),2,2)
rownames(model_equalrate)<-colnames(model_equalrate)<-colnames(states); model_equalrate
fit3 <- fitMk(tree,states,model=model_equalrate)

#asymmetric ordered irreversible 2
model_IRR_2 <- matrix(c(0,0,1,0),2,2)
rownames(model_IRR_2)<-colnames(model_IRR_2)<-colnames(states); model_IRR_2
fit4 <- fitMk(tree, states, model=model_IRR_2)


#compare models
aic<-setNames(sapply(list(fit1, fit2, fit3, fit4),AIC), c("model_asymOR","model_IRR", "model_equalrate", "model_IRR_2"))
aic
aic.w(aic)


#ancestral reconstruction of best model 
mtrees_equalrate<-make.simmap(tree,states,model=model_equalrate,nsim=500)
MkMOD<-summary(mtrees_equalrate)   

#save output !!!
#save(MkMOD, file="MkMOD_FINAL.Rdata")


#Test for phylogenetic signal in best model
lk.lambda<-function(lambda,tree,x,...) 
  -logLik(fitMk(phytools:::lambdaTree(tree,lambda),
                x,...))
opt<-optimize(lk.lambda,c(0,phytools:::maxLambda(tree)),tree=tree,
              x=states,model=model_equalrate)
opt
lam<-opt$minimum
lam

fit.lambda<-fitMk(phytools:::lambdaTree(tree,lam),states,model=model_equalrate)
fit.h0<-fitMk(phytools:::lambdaTree(tree,0),states,model=model_equalrate)
LR<--2*(logLik(fit.h0)-logLik(fit.lambda))
P.chisq<-pchisq(LR,df=1,lower.tail=FALSE)
P.chisq


##########################################################################
######                                                              ######
######      (6)  Figures from Threshold and Mk Models               ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######           uses "mcmcBM_combined_FINAL.Rdata" generated in (4)######
######           uses "PP_combined_FINAL.rds" generated in (4)      ######
######           uses "MkMOD_FINAL.Rdata" generated in (5)          ######
######                                                              ######
##########################################################################

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2)
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states:

#read in edited tree generated in (2), combined mcmc chains and posterior probabilities from (4), and Mk model from (5)
tree<- read.tree("edited_tree.tre")  
load("mcmcBM_combined_FINAL.Rdata")
PP <- readRDS("PP_combined_FINAL.rds")
load("MkMOD_FINAL.Rdata")


#Tree with original tip state probabilities
plotTree(tree,type="fan",ftype="i",fsize=0.6)
tiplabels(pie=states[tree$tip.label,],piecol=c("white","black"),cex=0.3) 

#Tree with AncThresh estimated tips and nodes
plotTree(tree,type="fan",ftype="i",fsize=0.6)
tiplabels(pie=PP[tree$tip.label,],piecol = c("black",  "white"), cex = .3) 
nodelabels(node = as.numeric(rownames(mcmcBM$ace)), pie = as.matrix(mcmcBM$ace), piecol = c("black","white"), cex = .4) 

#Tree with Mk estimated tips and nodes
plotTree(tree,type="fan",ftype="i",fsize=0.6)
tiplabels(pie=MkMOD$tips,piecol = c("black",  "white"), cex = .3) #need to add third color if using three state model
nodelabels(node = as.numeric(rownames(MkMOD$ace)), pie = as.matrix(MkMOD$ace), piecol = c("black",  "white"), cex = .4) 


#Phenogram from Threshold Model
x<-colMeans(mcmcBM$liab)
phenogram(tree,x,spread.labels=TRUE,spread.cost=c(1,0), link=15, offset=.4, fsize=0.5,
          axes=list(trait=c(min(x),max(x)), time=c(1.5,300)), xlab="", ylab="Urban Tolerance Liability")
th<-colMeans(mcmcBM$par[,c("avoid", "urban")])
lines(c(0,100),c(th[1],th[1]),lty="dashed",col="red") #line where they split



#compare threshold model and best Mk model

#internal nodes
plot(as.matrix(mcmcBM$ace), as.matrix(MkMOD$ace), xlab="BM ancthresh model", ylab="Mk model")

#tips only  
thresh.tips1<-PP[rownames(MkMOD$tips),colnames(MkMOD$tips)]
plot(MkMOD$tips,thresh.tips1, xlab="BM ancthresh model", ylab="Mk model")

#combine ace and tips for both models and color in plot
MkcomboPP <- (rbind(MkMOD$tips, MkMOD$ace))
MkPP_melt <- melt(MkcomboPP) 
colnames(MkPP_melt) <- c("species", "Mkstate", "Mkvalue")
PPmelt <- melt(PP)
colnames(PPmelt) <- c("species", "Thstate", "Thvalue")
combinedMkBM <- data.frame(MkPP_melt, PPmelt, node = c(rep("tip", 131), rep("node",nrow(MkPP_melt)-131 ))) 

ggplot(combinedMkBM, aes(x=Thvalue, y=Mkvalue)) + geom_point(cex=2, aes(color=node)) + theme_bw() + geom_abline(lty=2) + 
  scale_colour_grey(start=.6, end=0) + xlab("Threshold Model") + ylab("Mk Model") +
  geom_point(data = subset(combinedMkBM, node == 'tip'), aes(x = Thvalue, y = Mkvalue, color = node)) + 
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color="gray",linetype="solid"))

#plot facing trees
layout(matrix(1:3,1,3),widths=c(0.42,0.1,0.42))
plotTree(tree,ftype="off")
tiplabels(pie=PP[tree$tip.label,],piecol = c("black", "white"), cex = .5) 
nodelabels(node = as.numeric(rownames(mcmcBM$ace)), pie = as.matrix(mcmcBM$ace), piecol = c("black",  "white"), cex = .5) 
plot.new()
plot.window(xlim=c(-0.1,0.1),ylim=c(1, length(tree$tip.label)))
par(cex=.5)
text(rep(0,length(tree$tip.label)), 1:length(tree$tip.label),tree$tip.label)
plotTree(tree,direction="leftwards", ftype="off")
tiplabels(pie=MkMOD$tips,piecol = c("black",  "white"), cex = .5) 
nodelabels(node = as.numeric(rownames(MkMOD$ace)), pie = as.matrix(MkMOD$ace), piecol = c("black", "white"), cex = .5) 


#Plot threshold model with custom color gradient
mcmcliab<-colMeans(mcmcBM$liab)
CPancthresh<- contMap(tree, mcmcliab[1:131],plot=FALSE, method="user",anc.states=mcmcliab[132:length(mcmcliab)])

th<-colMeans(mcmcBM$par[,c("avoid","urban")]) #our thresholds
threshold<-th[1]

colors<-list(c("yellow","lightgreen","royalblue2"),c("royalblue2","darkorchid1", "red")) #color palette for the segments


#generate color gradient
lims<-CPancthresh$lims
colvals<-seq(lims[1],lims[2],by=diff(lims)/(length(CPancthresh$cols)-1))
ii<-which((colvals)^2==min((colvals)^2))
newcols<-colorRampPalette(colors[[1]])(ii)
jj<-which((colvals-lims[2])^2==min((colvals-lims[2])^2))
newcols<-c(newcols,colorRampPalette(colors[[2]])(jj-ii))
CPancthresh$cols[]<-newcols

leg.width<-8 #desired lengend width
h<-max(nodeHeights(CPancthresh$tree)) #to adjust space at tips
offset.factor<-1.02 ## increase this for greater offset


plotTree((CPancthresh$tree), xlim=c(-1.6*h,1.1*h),
         color="transparent",ftype="i",type="fan",fsize=1, add=F)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)

h<-max(nodeHeights(CPancthresh$tree))
plot(CPancthresh,type="fan",legend=FALSE,xlim=obj$x.lim,ylim=obj$y.lim, fsize=0.01, ftype="i", add=T)

plotTree((CPancthresh$tree), xlim=c(-1.6*h,1.1*h),color="transparent",ftype="i",type="fan",fsize=1, add=T)

add.color.bar(2*h,CPancthresh$cols,title="liability",
              lims=NULL,digits=3,direction="upwards",
              subtitle="",lwd=leg.width,x=-1.6*h,y=-h,prompt=FALSE)
LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
lines(x=rep(-1.6*h+LWD*leg.width/2,2),y=c(-h,h))
ticks<-c(lims[1],threshold,lims[2])
nticks<-length(ticks)
Y<-cbind((ticks-min(ticks))/diff(range(ticks))*2*h-h,
         (ticks-min(ticks))/diff(range(ticks))*2*h-h)
X<-cbind(rep(-1.6*h+LWD*leg.width/2,nticks),
         rep(-1.6*h+LWD*leg.width/2+0.04*h,nticks))
for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
text(x=X[,2],y=Y[,2],c("avoid", "", "urban"),pos=4,cex=1) 



##########################################################################
######                                                              ######
######      (7)  Leave one out analysis: Mk Model                   ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######                                                              ######
##########################################################################


#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2)
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states:

#read in edited tree generated in (2) and Mk model from (5)
tree<- read.tree("edited_tree.tre")  

#Specify which columns are tip states:
states <- as.matrix(treespecies[,c(7,8)]) 
X<-states
X<-X[tree$tip.label,]

model_equalrate<-matrix(c(0,1,1,0),2,2)
rownames(model_equalrate)<-colnames(model_equalrate)<-colnames(states); model_equalrate

fitQ<-fitMk(tree,X,model=model_equalrate)

Q<-matrix(NA,length(fitQ$states),length(fitQ$states))
Q[]<-c(0,fitQ$rates)[fitQ$index.matrix+1]
diag(Q)<-0
diag(Q)<--rowSums(Q)
colnames(Q)<-rownames(Q)<-fitQ$states

PP<-X

### save Q - needed for (8) !!!
#write.csv(Q, "Q.csv")

################ RUNNING LEAVE ONE OUT MK LOOP #################

### WARNING: This will take 36 hours to run! 

for(i in 1:nrow(X)){
Xp<-X
Xp[i,]<-rep(1/2,2)
fit<-rerootingMethod(tree,Xp,fixedQ=Q)
PP[i,]<-fit$marginal.anc[i,]
cat(paste("Done loop",i,"\\n"))
flush.console()
}

#save output !!!
#saveRDS(PP, "Mkloo_FINAL.rds")



##########################################################################
######                                                              ######
######      (8)  Mk Leave One Out Comparison to Original Model      ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######           uses "MkMOD_FINAL.Rdata" generated in (5)          ######
######           uses "Q.csv" generated in (7)                      ######
######           uses "Mkloo_FINAL.rds" generated in (7)            ######
######                                                              ######
##########################################################################


#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2)
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states

#read in edited tree generated in (2), Mk model from (5), and Q from (7)
tree<- read.tree("edited_tree.tre")
PP <- readRDS("Mkloo_FINAL.rds")
Q <- as.matrix(read.csv("Q.csv", row.names=1))

#Specify which columns are tip states:
states <- as.matrix(treespecies[,c(7,8)])
X<-states
X<-X[tree$tip.label,]


#How well does it recover the known state:
obj<-rerootingMethod(tree,X,fixedQ=Q)
obj$marginal.anc[rownames(PP),]->PP.full

X <- PP.full 
Xe <- PP

PPscale<- mean(apply((states),1,max)) #scale by mean row max
right<-sum(Xe*X)/nrow(X)
right/PPscale 

#significance
Xexp<-matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow=TRUE)
by.chance<-sum(Xexp*X)/nrow(X)

p<-vector()
for(i in 1:100000){
  Xp<-X[sample(nrow(X)),]
  p[i]<-sum(Xp*X)/nrow(X)
}
P.value<-mean(p>=right)
P.value


##########################################################################
######                                                              ######
######      (9)  Leave One Out Analysis: Threshold Model            ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######           uses "mcmcBM_combined_FINAL.Rdata" generated in (4)######
######           uses "PP_combined_FINAL.rds" generated in (4)      ######
######                                                              ######
##########################################################################

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2)
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states

#read in edited tree generated in (2)
tree<- read.tree("edited_tree.tre")  

set.seed(61484) #use seed from our model that had highest convergence


X<-states
X<-X[tree$tip.label,]

#Specs for ancthresh
ngen <- 10000000  #number of generations to run MCMC
sample <- 1000 #sampling interval
burnin <- .2*ngen #burnin
nspecies <- length(tree$tip.label)
n <- nrow(states)
seq <- c("avoid", "urban")

PP<-X  

################ RUNNING LEAVE ONE OUT THRESHOLD LOOP #################
### WARNING: This will take 14 days to run! 

# for(i in 1:nrow(X)){
#  Xp<-X
#  Xp[i,]<-rep(1/2,2)
#  fit<-ancThresh(tree=tree, x=Xp, ngen = ngen, sequence=seq, model="BM", control = list(sample = sample, burnin=burnin))
#
#  L<-as.matrix(fit$liab)
#  Th<-as.matrix(fit$par[,2:3])
#  STATES<-matrix(NA,nrow(L),ncol(L),dimnames=list(rownames(L),colnames(L)))
#  for(k in 1:nrow(L)) STATES[k,]<-threshState(L[k,],Th[k,])
#  foo<-function(x){ y<-summary(factor(x,levels=c("avoid","urban"))); y/sum(y) }
#  PP_temp<-apply(STATES,2,foo)
#  PP_temp <- t(PP_temp)
#
#  PP[i,]<-PP_temp[i,]
#  cat(paste("Done loop",i,"\\n"))
#  flush.console()
# }

#save output !!! 
#saveRDS(PP, "PPthreshloop_FINAL.rds")


##########################################################################
######                                                              ######
######      (10) Threshold Model Leave One Out Compare to Original  ######
######           uses "TipStates_Final.csv" generated in (1)        ######
######           uses "edited_tree.tre" generated in (2)            ######
######           uses "PP_combined_FINAL.rds" generated in (4)      ######
######           uses "PPthreshloop_FINAL.rds" generated in (9)     ######
######                                                              ######
##########################################################################

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2)
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states

#read in edited tree generated in (2), combined mcmc chains and posterior probabilities from (4)
tree<- read.tree("edited_tree.tre")  
fullPP_all <- readRDS("PP_combined_FINAL.rds")

#read in leave one out analysis output from (9)
LOO_PP <- readRDS("PPthreshloop_FINAL.rds")


fullPP <- fullPP_all[1:131,] 
fullPP[rownames(LOO_PP),]->PP.full

X <- PP.full 
Xe <- LOO_PP

#How well does it recover the known state:
PPscale<- mean(apply((states),1,max))
right<-sum(Xe*X)/nrow(X)
right/PPscale

Xexp<-matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow=TRUE)
by.chance<-sum(Xexp*X)/nrow(X)

p<-vector()
for(i in 1:100000){
  Xp<-X[sample(nrow(X)),]
  p[i]<-sum(Xp*X)/nrow(X)
}
P.value<-mean(p>=right)
P.value



##########################################################################
######                                                              ######
######     (11) PGLS                                                ######
######         uses "TipStates_Final.csv" generated in (1)          ######
######         uses "edited_tree.tre" generated in (2)              ######
######         uses "mcmcBM_combined_FINAL.Rdata" generated in (4)  ######
######                                                              ######
##########################################################################

#read in tip states generated in (1)
treespecies <- read.csv("TipStates_Final.csv", header=T, row.names=2)
species <- as.character(row.names(treespecies))  #vector of 131 species names
states <- as.matrix(treespecies[,c("urban",  "avoid")])  #Specify which columns are tip states

#read in edited tree generated in (2), combined mcmc chains and posterior probabilities from (4)
tree<- read.tree("edited_tree.tre")  

#read in new datasets
bioclim <- read.csv("bioclim_FINAL.csv", row.names = 1) #climate data
traits <- read.csv("Traits_FINAL.csv") #phenotypic data
rownames(traits) <- traits$species

#load data from combined mcmcBM 
load("mcmcBM_combined_FINAL.Rdata") 


#need mean and SE for liability from posterior distribution
obj <- lapply(mcmcBM$liab[,1:131], function(x) {class(x) <- "mcmc" ; summary(x) })
stats <- t(sapply(obj, function(x) setNames(x$statistics[c(1,2)], c("estimate", "sd"))))
liabilities <- stats[,"estimate"]
sd <- stats[,"sd"]

#custom color palette for plotting
colors<-(c("yellow","lightgreen","royalblue2","darkorchid1", "red")) #color palette for the segments to match earlier plots


#Bioclim phylogenetic PCA
bioclim_phylopca <- phyl.pca(tree, bioclim[,c(2:20)], method="BM", mode="corr") #phylogenetic PCA

summary(bioclim_phylopca)
plot(bioclim_phylopca)
biplot(bioclim_phylopca)

diag(bioclim_phylopca$Eval ) #eiegenvalues
bioclim_phylopca$L[,c(1:4)] #loadings
bcpca <- data.frame(species =rownames(bioclim_phylopca$S), bioclim_phylopca$S[,c(1:4)]) #scores 

#Combine traits with liability from threshold model
correl <- data.frame(liability = liabilities[rownames(traits)], 
                     traits[,c(2:17)], 
                     se = sd[rownames(traits)], observations = treespecies$n_total  ) #liabilities and PP
correl$species <- row.names(correl)

levels(correl$island)[2] <- "Cuba" #count any species found in Bahamas and Cuba as "Cuba" for this analysis


correl_climate <- merge(correl, bcpca, by="species") #make a second dataframe for climate because 3 species no climate data ("ernestwilliamsi" "desechensis" "fairchildi" )
rownames(correl_climate) <- correl_climate$species



##### PGLS MODELS #####

### 1. Island
fit.island <- pgls.SEy(liability ~ island , data=correl, corClass=corBrownian, tree=tree, se=sd,  method="ML")
summary(fit.island)
anova(fit.island) 

isl.emm <- emmeans(fit.island, c("island")) #compare marginal means and group
isl.CLD <- CLD(isl.emm, Letters=LETTERS)

ggplot(correl, aes(x=island, y=liability)) + geom_boxplot() + theme_bw() + geom_hline(yintercept = 0, lty=2, col="red") + 
  geom_point(position=position_jitter(w = 0.25, h = 0), aes(col=ecomorph))  +
  annotate("text", x=c( 1:7), y=rep(10.5,7), label= c(isl.CLD[4,7], isl.CLD[3,7], isl.CLD[5,7], isl.CLD[2,7], isl.CLD[6,7], isl.CLD[7,7], isl.CLD[1,7]),  size=5) +
  scale_x_discrete(name ="Island",  labels=c("BAH","CUBA","CAYM","HISP", "JAM", "LA", "PR"))


### 2. Ecomorph
coreco <- correl[correl$ecomorph!="0",] #drop non ecomorphs
tree_temp <- drop.tip(tree, tree$tip.label[-match(rownames(coreco),tree$tip.label)]) #drop tips from tree not in our trait dataset

fit.eco.island <- pgls.SEy(liability ~   island + ecomorph , data=coreco, corClass=corBrownian, tree=tree_temp, se=sd,  method="ML")
#fit model with island

summary(fit.eco.island)
anova(fit.eco.island) 

ggplot(correl, aes(x=ecomorph, y=liability)) + geom_boxplot() + theme_bw() + geom_hline(yintercept = 0, lty=2, col="red") + 
  geom_point(position="jitter", aes(col=island)) 


### 3. Body temperature 
corBT <- correl[is.na(correl$Tb.C)!=T,] #subset dataset because not all species have Tb data
tree_temp <- drop.tip(tree, tree$tip.label[-match(rownames(corBT),tree$tip.label)]) #drop tips from tree not in our dataset

fit.BT.island <- pgls.SEy(liability ~ island + Tb.C , data=corBT, corClass=corBrownian, tree=tree_temp, se=sd,  method="ML")
summary(fit.BT.island)  
anova(fit.BT.island)


### 4. Morphology 
cormorph_full <- correl[is.na(correl$dorsal.scale)!=T & is.na(correl$ventral.scale)!=T & is.na(correl$Flimb)!=T,]  #subset to shared dataset
tree_temp <- drop.tip(tree, tree$tip.label[-match(rownames(cormorph_full),tree$tip.label)]) #drop tips not in our trait dataset

fit.morph.island <- pgls.SEy(liability ~ island +  
                               log(dorsal.scale) +  log(ventral.scale) +  
                               log(RLAM)+log(FLAM) + log(Hlimb)  +log(Flimb) +log(SVL),  
                      data=cormorph_full, corClass=corBrownian, tree=tree_temp, se=sd, method="ML")
summary(fit.morph.island) #full model
anova(fit.morph.island)

#stepwise model selection 
fit.morph.step<-stepAIC(fit.morph.island,direction="both", trace=0, scope = list(lower = ~ log(SVL)+ island)) #keep SVL
summary(fit.morph.step) #reduced model
anova(fit.morph.step)
#drop1(fit.morph.step, test="Chisq") #drop method to check term order doesn't affect result


### 5. Ecology
tree_temp <- drop.tip(tree, tree$tip.label[-match(rownames(correl_climate),tree$tip.label)]) #drop tips from tree not in our trait dataset

fit.clim.island <- pgls.SEy(liability ~ island + PC1 + PC2+PC3 + PC4 + log(range_area_km2)+  log(Lat_range) + Nsym_MAP ,  
                     data=correl_climate, corClass=corBrownian, tree=tree_temp, se=sd,  method="ML")
summary(fit.clim.island) #full model
anova(fit.clim.island)


fit.clim.step<- stepAIC(fit.clim.island, direction="both", trace=0, scope = list(lower = ~ island))
summary(fit.clim.step) #reduced model
anova(fit.clim.step)
#drop1(fit.clim.step, test="Chisq") #drop method to check term order doesn't affect result


