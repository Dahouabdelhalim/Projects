rm(list=ls())

#please carefully review the README file first!

###########
#geckoEyes#
###########

#Lars Schmitz and Tim Higham, ecomorph@gmail.com

#####Preliminaries

#set working directory, please make sure to change it to match your directory

#PC
setwd("C:/Users/lschmitz/Dropbox/geckoEyes") #adjust!

#MAC
#setwd("/Users/larsschmitz/dropbox/geckoEyes") #adjust!


#load libraries
require(phytools)
require(ape)
require(geiger)
require(bayou)
require(nlme)
require(MASS)
require(ggplot2)
require(OUwie)
require(mvMORPH)


#load data and trees
wiens <- read.tree("./phylogenies/wiens.time.tre")
geckos <- read.csv("./data/tableS1.csv", header=TRUE, stringsAsFactors=FALSE)

#swap in species names for Lepidoblepharis sp (=> festae) and Pesudogonatodes sp (=> guianensis)
lepido <- which(geckos$taxon == "Lepidoblepharis_sp")
geckos$taxon[lepido] <- "Lepidoblepharis_festae"

pseudo <- which(geckos$taxon == "Pseudogonatodes_sp")
geckos$taxon[pseudo] <- "Pseudogonatodes_guianensis"

#align data with phylogeny
rownames(geckos) <- geckos$taxon
compare.wiens <- treedata(wiens, geckos, sort=T, warnings=T)
wiens.tre <- compare.wiens$phy
wiens.data <- data.frame(taxon=as.character(compare.wiens$data[,1]), 
                         dap=as.character(compare.wiens$data[,2]),
                         cat=as.character(compare.wiens$data[,3]),
                         svl=as.numeric(compare.wiens$data[,4]),
                         ed=as.numeric(compare.wiens$data[,5]))
rownames(wiens.data) <- wiens.data$taxon
wiens.data <- wiens.data[wiens.tre$tip.label,]

#basic descriptors
length(which(wiens.data$dap=="cathemeral")) #that's cathemeral/crepuscular combined
length(which(wiens.data$dap=="diurnal"))
length(which(wiens.data$dap=="nocturnal"))



########
#Step 1#
########



#reconstruct history of activity patterns

#make DAP vector and colors that go with it
dap <- as.vector(wiens.data$dap)
names(dap) <- wiens.tre$tip.label
dap.colors <- rep("black", length(dap)) 
diurnals <- which(dap=="diurnal")
cathemerals <- which(dap=="cathemeral")
dap.colors[diurnals] <- "gold"
dap.colors[cathemerals] <- "grey"

#plot tree and add tip labels according to DAP
pdf(file="./figures/labeled_tips.pdf", 6, 8, useDingbats=F)
plot(wiens.tre, cex=0.3, label.offset=2)
tiplabels(pch=21, bg=dap.colors, adj=1)
dev.off()

#initial exploration with maximum likelihood
ER <- ace(dap, wiens.tre, type="discrete", model="ER")
SYM <- ace(dap, wiens.tre, type="discrete", model="SYM")
ARD <- ace(dap, wiens.tre, type="discrete", model="ARD")

#symmetric-rates and all-rates-different models throw errors: 
#warning message:In sqrt(diag(solve(h))) : NaNs produced
#probably related to insufficeintly large tip number for 3 states
# => only use ER model

#root state estimates
ER$lik.anc[1,]

#stochastic character mapping
dap.simmap <- make.simmap(wiens.tre, dap, nsim=1000, model="ER") 

#plot one example
regcols <- setNames(c("grey", "gold", "black"),c("cathemeral", "diurnal", "nocturnal"))

#ER
pdf(file="./figures/example.mapping.pdf", 6, 8, useDingbats=F)
plotSimmap(dap.simmap[[1]], regcols, ftype="i", fsize=0.6)
dev.off()

pdf(file="./figures/4.example.mappings.pdf", 6, 8, useDingbats=F)
par(mfrow=c(2,2))
plotSimmap(dap.simmap[[1]], regcols, ftype="i", fsize=0.3, lwd=1)
plotSimmap(dap.simmap[[2]], regcols, ftype="i", fsize=0.3, lwd=1)
plotSimmap(dap.simmap[[3]], regcols, ftype="i", fsize=0.3, lwd=1)
plotSimmap(dap.simmap[[4]], regcols, ftype="i", fsize=0.3, lwd=1)
dev.off()

#summary of the 1000 trees
dap.evol <- summary(dap.simmap)

#node estimates (needed for later)
node_estimates_ER_basic <- dap.evol$ace
scores_ER_basic <- apply(node_estimates_ER_basic, 1, which.max)

#plot
pdf(file="./figures/dap.evolution.pdf", 6, 8, useDingbats=F)
plot(dap.evol, colors=regcols, ftype="i", fsize=0.5, lwd=1)
add.simmap.legend(leg=c("c", "d", "n"), colors=regcols,x=0,y=78,prompt=FALSE)
dev.off()



###FIGURE 1a###
#plot for manuscript
pdf(file="./figures/dap.evolution_FIG.pdf", 3, 8, useDingbats=F)
plot(dap.evol, colors=regcols, ftype="i", fsize=0.4, lwd=0.5)
add.simmap.legend(leg=c("c", "d", "n"), colors=regcols,x=0,y=78,prompt=FALSE)
dev.off()



#descriptive summary
dap.evol.sum <- describe.simmap(dap.simmap)
dap.evol.sum

#get the minimum age of some diurnal radiations mentioned in the text
max(nodeHeights(wiens.tre)) - fastHeight(wiens.tre, "Phelsuma_parkeri", "Lygodactylus_pictus")
max(nodeHeights(wiens.tre)) - fastHeight(wiens.tre, "Sphaerodactylus_argus", "Gonatodes_annularis")
max(nodeHeights(wiens.tre)) - fastHeight(wiens.tre, "Naultinus_manukanus", "Naultinus_elegans")


#are ASR biased by small sample size? what if we used all available information on DAP?
#let's try it out

#load full dataset on DAP
dap.geckos <- read.csv("./data/tableS2.csv")

#align these data with phylogeny
rownames(dap.geckos) <- dap.geckos$taxon
compare.wiens.again <- treedata(wiens, dap.geckos, sort=T, warnings=T) # the treedata function compares data with tree and prunes the tree automatically

wiens.tre.dap <- compare.wiens.again$phy
wiens.data.dap <- data.frame(taxon=as.character(compare.wiens.again$data[,1]), 
                         dap=as.character(compare.wiens.again$data[,2]))
rownames(wiens.data.dap) <- wiens.data.dap$taxon
wiens.data.dap <- wiens.data.dap[wiens.tre.dap$tip.label,]

droppers <- setdiff(as.vector(wiens.data.dap$taxon), as.vector(wiens.data$taxon)) #need this for pruning later!

dap <- as.vector(wiens.data.dap$dap)
names(dap) <- wiens.tre.dap$tip.label

#stochastic character mapping (stick to equal rates)
dap.simmap.full <- make.simmap(wiens.tre.dap, dap, nsim=1000, model="ER") 

#summary of the 1000 trees
dap.evol.full <- summary(dap.simmap.full)

#plot
pdf(file="./figures/dap.full.pdf", 6, 20)
plot(dap.evol.full, colors=regcols, ftype="i", fsize=0.5, lwd=1)
add.simmap.legend(leg=c("c", "d", "n"), colors=regcols,x=0,y=240, fsize=2, prompt=FALSE)
dev.off()

png(file="./figures/dap.full.figure.png", 8, 11, res=600, units="in")
plot(dap.evol.full, colors=regcols, ftype="i", fsize=0.4, lwd=1)
add.simmap.legend(leg=c("c", "d", "n"), colors=regcols,x=0,y=240, fsize=2, prompt=FALSE)
dev.off()

#get discrete node estimates
node_estimates_ER_full <- dap.evol.full$ace
scores_ER_full <- apply(node_estimates_ER_full, 1, which.max)


# 1: cathemeral/crepuscular
# 2: diurnal
# 3: nocturnal

#assign node estimates to tree
wiens.tre.dap$node.label <- scores_ER_full

plot(wiens.tre.dap); nodelabels(wiens.tre.dap$node.label)


#prune tree to match species with eye data
pruned.tre <- drop.tip(wiens.tre.dap, droppers)
scores_ER_reduced <- pruned.tre$node.label

#check for discrepancies among node estimates => these should all be equal to zero
plot(pruned.tre); nodelabels(scores_ER_basic-scores_ER_reduced, frame="circle", bg="light green")

#this seems fine


########
#Step 2#
########



#PGLS

#preparing data
groups <- wiens.data$dap
names(groups) <- wiens.tre$tip.label
svl <- signif(log10(wiens.data$svl), 3)
names(svl) <- wiens.tre$tip.label
ed <- signif(log10(wiens.data$ed),3)
names(ed) <- wiens.tre$tip.label

#putting it all in a dataframe
eyes<-data.frame(taxon=as.character(wiens.tre$tip.label), svl=as.numeric(svl),ed=as.numeric(ed), groups=as.character(groups))
rownames(eyes) <- wiens.tre$tip.label

#defining different subgroups
cat1 <- eyes[eyes$groups=="cathemeral",]
cat2 <- eyes[eyes$groups=="diurnal",]
cat3 <- eyes[eyes$groups=="nocturnal",]



#PGLS with BM correlation matrix

#simple PGLS model: ED ~ SVL
BM1 <- gls(ed ~ svl, data=eyes, correlation=corPagel(0.1, wiens.tre, fixed=FALSE), method="ML")
BM1_sum <- summary(BM1)

#check for heteroscedasticity
plot(BM1, resid(., type="n")~fitted(.), col="blue", main="Normalized Residuals vs. Fitted Values",
     abline=c(0,0))

#check for departures from normal distribution of residuals
res <- resid(BM1, type="n")
qqnorm(res, col="blue")
qqline(res, col="blue")

#define the range for lambda and at which intervals we are sampling
lambda <- seq(0, 1, length.out = 100)

#Calculate likelihood over the entire range[modified from Symonds&Blombergs tutorial]
lik <- sapply(lambda, function(lambda) logLik(gls(ed ~ svl, data=eyes, correlation = corPagel(value=lambda, phy=wiens.tre, fixed = TRUE))))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#adding our ML estimate as a red line
abline(v=BM1$modelStruct, col = "red")

#intercept PGLS model: ED ~ SVL + DAP
BM2 <- gls(ed ~ svl+groups, data=eyes, correlation=corPagel(0.1,wiens.tre,fixed=FALSE),method="ML")
BM2_sum <- summary(BM2)

#check for heteroscedasticity
plot(BM2, resid(., type="n")~fitted(.), col="blue", main="Normalized Residuals vs. Fitted Values",
     abline=c(0,0))

#check for departures from normal distribution of residuals
res <- resid(BM2, type="n")
qqnorm(res, col="blue")
qqline(res, col="blue")

#define the range for lambda and at which intervals we are sampling
lambda <- seq(0, 1, length.out = 100)

#Calculate likelihood over the entire range[modified from Symonds&Blombergs tutorial]
lik <- sapply(lambda, function(lambda) logLik(gls(ed ~ svl+groups, data=eyes, correlation = corPagel(value=lambda, phy=wiens.tre, fixed = TRUE))))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#adding our ML estimate as a red line
abline(v=BM2$modelStruct, col = "red")

#evaluating models
BM1_sum$AIC
BM2_sum$AIC

#let's calculate AICc scores
aic<-function(logL, k) -2*(logL)+2*k #testing if the function works... (compare to results above)
aic(BM1_sum$logLik, 4)
aic(BM2_sum$logLik, 6)

#AICc
aicc<-function(logL, k, n) -2*(logL)+2*k*(n/(n-k-1))

aicc(BM1_sum$logLik, 99, 4)
aicc(BM2_sum$logLik, 99, 6)

#deltaAICc of BM1
delta1 <- aicc(BM1_sum$logLik, 99, 4)-aicc(BM2_sum$logLik, 99, 6)
delta1

#Akaike weights
akwBM1 <- exp(-0.5 * delta1) / (exp(-0.5 * delta1) + exp(-0.5 * 0))
akwBM2 <- exp(-0.5 * 0) / (exp(-0.5 * delta1) + exp(-0.5 * 0))



#PGLS with OU correlation matrix

#simple PGLS; having computational issues with alpha <1.6
OU1 <- gls(ed~svl, data=eyes, correlation=corMartins(2, wiens.tre, fixed=FALSE), method="ML")
OU1_sum <- summary(OU1)
OU1_sum$AIC #larger than the corresponding BM model

#given the computational problems for small alpha, let's calculate this model over a range of different alpha value before doing anything else
#define the range for alpha and at which intervals we are sampling
alpha <- seq(0, 10, length.out = 100)

#Calculate likelihood over the entire range[modified from Symonds&Blombergs tutorial]
lik <- sapply(alpha, function(alpha) logLik(gls(ed ~ svl, data=eyes, correlation = corMartins(value=alpha, phy=wiens.tre, fixed = TRUE))))
plot(lik ~ alpha, type = "l", main = expression(paste("Likelihood Plot for ",alpha)), ylab = "Log Likelihood", xlab = expression(alpha))

#adding our ML estimate as a red line
abline(v=OU1$modelStruct, col = "red")

#likelihood optimization is failing
#likelihood goes up with very small alpha
#stick to BM correlation matrix

#plot results!

#make color vector for cathemeral/crepuscular, diurnal, nocturnal (alphabetical order!)
color <- c("grey", "gold", "black")

#creating a dataframe for x coordinates
xnew <- data.frame(range(cat1[,2]), range(cat2[,2]), range(cat3[,2]))

#next you want to extract the coefficients from the regression:
cof <- BM2$coef
slope <- cof[2]
int <- c(cof[1], cof[1]+cof[3], cof[1]+cof[4])

# plotting
pdf(file="./figures/PGLS_BM_groups.pdf", 4, 4, useDingbats=F)
plot(cat2$svl, cat2$ed, 
     pch=21, col="black", bg="gold", 
     xlim=range(eyes$svl), ylim=range(eyes$ed),
     xlab="log10 SVL [mm]", ylab="log10 Eye Diameter [mm]", cex.lab=1, cex=1)
points(cat3$svl, cat3$ed, 
       pch=21, col="black", bg="black", cex=1)
points(cat1$svl, cat1$ed, 
       pch=21, col="black", bg="grey", cex=1)
for(i in 1:3) {
  yhat <- slope*xnew[,i] + int[i]
  lines(xnew[,i], yhat, col=color[i], lwd=2)}
box(lwd=2)
dev.off()




###FIGURE 1b###
# plotting for manuscript figure
pdf(file="./figures/PGLS_BM_groups_FIG.pdf", 4, 4, useDingbats=F)
par(mar=c(3,3,0.5,0.5), mgp=c(2,0.5,0))
plot(cat2$svl, cat2$ed, las=1,
     pch=21, col="black", bg="gold", 
     xlim=range(eyes$svl), ylim=range(eyes$ed),
     xlab="log10 SVL [mm]", ylab="log10 Eye Diameter [mm]", cex.lab=1, cex.axis=0.75, cex=1, tck=-0.01)
points(cat3$svl, cat3$ed, 
       pch=21, col="black", bg="black", cex=1)
points(cat1$svl, cat1$ed, 
       pch=21, col="black", bg="grey", cex=1)
for(i in 1:3) {
  yhat <- slope*xnew[,i] + int[i]
  lines(xnew[,i], yhat, col=color[i], lwd=2)}
box(lwd=2)
dev.off()



#stripchart (residual dot plot)

#residuals from simple model (all combined)
res <- as.vector(BM1$residuals)
names(res) <- wiens.tre$tip.label

cat <- as.vector(wiens.data$cat)
names(cat) <- wiens.tre$tip.label

#putting it all in a dataframe
strip.data <- data.frame(taxon=as.character(wiens.tre$tip.label), 
                         res=as.numeric(res), cat=as.factor(cat))
                                           
rownames(strip.data) <- wiens.tre$tip.label

#colors
fill <- c(rep("grey", 4), rep("gold", 8), rep("black", 3))

#basic stripchart
ggplot(strip.data, aes(x=cat, y=res, color=cat)) + geom_jitter()

# Change the position
# 0.2 : degree of jitter in x direction
# horizontal red line, dashed
p <- ggplot(strip.data, aes(x=cat, y=res, color=cat)) + geom_jitter(position=position_jitter(0.1), size=2) + geom_hline(aes(yintercept=0))
p

#use custom color palettes
p2 <- p + scale_color_manual(values=fill) + theme_classic() + theme(legend.position="none") + scale_x_discrete(limits=c("ORIG", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "N1", "N2"))

final_plot <- p2 + labs(x="Evolutionary Bins of DAP", y="Residual Eye Size") + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=8), axis.title=element_text(size=14))
final_plot



###FIGURE 1c###
pdf(file="./figures/stripchart_WIENSv2.pdf", 4, 4, useDingbats=F)
final_plot
dev.off()



########
#Step 3#
########

#model fitting with mvMORPH
data <- data.frame(svl, ed)

#set up empty vectors for storing results
all.res <- matrix(NA, nrow = 1000, ncol = 18)

#start the loop
for (i in 1:1000) try({
OU1<- mvOU(dap.simmap[[i]], data, model="OU1")
OUM<- mvOU(dap.simmap[[i]], data, model="OUM")

all.res[i,][1] <- OU1$AICc
all.res[i,][2] <- OUM$AICc

all.res[i,][3] <- OUM$theta[1]
all.res[i,][4] <- OUM$theta[2]
all.res[i,][5] <- OUM$theta[3]
all.res[i,][6] <- OUM$theta[4]
all.res[i,][7] <- OUM$theta[5]
all.res[i,][8] <- OUM$theta[6]
all.res[i,][9] <- OUM$theta[7]
all.res[i,][10] <- OUM$theta[8]

all.res[i,][11] <- OUM$alpha[1]
all.res[i,][12] <- OUM$alpha[2]
all.res[i,][13] <- OUM$alpha[3]
all.res[i,][14] <- OUM$alpha[4]

all.res[i,][15] <- OUM$sigma[1]
all.res[i,][16] <- OUM$sigma[2]
all.res[i,][17] <- OUM$sigma[3]
all.res[i,][18] <- OUM$sigma[4]
})

write.csv(all.res, file="./figures/mvmorph_res.csv")

pdf(file="./figures/mvmorph.deltaAICc.pdf", 4, 4, useDingbats=F)
hist(all.res[,1]-all.res[,2], main="deltaAICc", xlab="AICc(OU1) - AICc(OUM)", col="light blue")
dev.off()

mean(all.res[,1]-all.res[,2])
range(all.res[,1]-all.res[,2])
sd(all.res[,1]-all.res[,2])

#estimate type I and II error rates through simulation

#WARNING: this will run a couple hours!

#number of simulations
nsim=1000

#dataset simulated with the OUM maximum likelihood estimates
data1<-simulate(OUM, nsim=nsim, tree=dap.simmap[[1]])

#dataset simulated with the MLE of the OU1 model
data2<-simulate(OU1, nsim=nsim, tree=dap.simmap[[1]])

#model fitting for simulated data #1 (OUM)
oum_data1<- lapply(1:nsim, function(x){
  mvOU(dap.simmap[[1]], data1[[x]], model="OUM", method="sparse", diagnostic=F, echo=F)
})

ou1_data1<- lapply(1:nsim, function(x){
  mvOU(dap.simmap[[1]], data1[[x]], model="OU1", method="sparse", diagnostic=F, echo=F)
})


#same for the other simulated dataset (with OU1)
oum_data2<- lapply(1:nsim, function(x){
  mvOU(dap.simmap[[1]], data2[[x]], model="OUM", method="sparse", diagnostic=F, echo=F)
})

ou1_data2<- lapply(1:nsim, function(x){
  mvOU(dap.simmap[[1]], data2[[x]], model="OU1", method="sparse", diagnostic=F, echo=F)
})


# get results (AICc) from simulations
OUM_simul<-sapply(1:nsim, function(x){
  c(oum_data1[[x]]$AICc,ou1_data1[[x]]$AICc)
})

OU1_simul<-sapply(1:nsim, function(x){
  c(oum_data2[[x]]$AICc,ou1_data2[[x]]$AICc)
})

#now compute type I error and power (type II)
sum(OU1_simul[1,]<OU1_simul[2,])/nsim
sum(OUM_simul[1,]<OUM_simul[2,])/nsim


#let's run this over a subset of SIMMAP trees

#number of simulations
nsim=100

#set up empty vectors for storing results
sim.res <- matrix(NA, nrow = 100, ncol = 2)

for (i in 1:2){

  #dataset simulated with the OUM maximum likelihood estimates
  data1<-simulate(OUM, nsim=nsim, tree=dap.simmap[[1]])

  #dataset simulated with the MLE of the OU1 model
  data2<-simulate(OU1, nsim=nsim, tree=dap.simmap[[1]])

  #model fitting for simulated data #1 (OUM)
  oum_data1<- lapply(1:nsim, function(x){
      mvOU(dap.simmap[[1]], data1[[x]], model="OUM", method="sparse", diagnostic=F, echo=F)
    })

  ou1_data1<- lapply(1:nsim, function(x){
      mvOU(dap.simmap[[1]], data1[[x]], model="OU1", method="sparse", diagnostic=F, echo=F)
    })

  #same for the other simulated dataset (with OU1)
  oum_data2<- lapply(1:nsim, function(x){
      mvOU(dap.simmap[[1]], data2[[x]], model="OUM", method="sparse", diagnostic=F, echo=F)
    })

  ou1_data2<- lapply(1:nsim, function(x){
      mvOU(dap.simmap[[1]], data2[[x]], model="OU1", method="sparse", diagnostic=F, echo=F)
    })

  #get results (AICc) from simulations
  OUM_simul<-sapply(1:nsim, function(x){
      c(oum_data1[[x]]$AICc,ou1_data1[[x]]$AICc)
    })

  OU1_simul<-sapply(1:nsim, function(x){
      c(oum_data2[[x]]$AICc,ou1_data2[[x]]$AICc)
    })

  #now compute type I error and power (type II)
  out1 <- sum(OU1_simul[1,]<OU1_simul[2,])/nsim
  out2 <- sum(OUM_simul[1,]<OUM_simul[2,])/nsim

  sim.res[i,][1] <- out1
  sim.res[i,][2] <- out2

  cat("loop", i, "\\n")
  
}

write.csv(sim.res, "sim.res.csv")



########
#Step 4#
########

#WARNING: this will run several hours!

#bayou
#residuals for bayOU
res <- as.vector(BM1$residuals)
names(res) <- wiens.tre$tip.label
res <- res*10
res



#prior
prior <- make.prior(wiens.tre, dists=list(dalpha="dlnorm", 
                                          dsig2="dlnorm",
                                          dsb="dsb", 
                                          dk="cdpois", 
                                          dtheta="dnorm",
                                          dloc="dunif"),
                    param=list(dk=list(lambda=20, kmax=196),
                               dtheta=list(mean=mean(res), sd=5), 
                               dsig2=list(meanlog=2, sdlog=1),
                               dalpha=list(meanlog=-1, sdlog=1)
                    )
)




fit1 <- bayou.mcmc(wiens.tre, res, SE=0.001, model="OU", prior, ngen=2000000, new.dir=getwd(), plot.freq=NULL, ticker.freq=10000)
chain <- load.bayou(fit1, save.Rdata=T, cleanup=T)
chain <- set.burnin(chain, 0.3)
out <- summary(chain)

write.csv(out$statistics, "./bayou_results/bayOU_a_stats.csv")
write.csv(out$branch.posteriors, "./bayou_results/bayOU_a_posteriors.csv")

pdf("./bayou_results/bayOUpar_a.pdf")
par(mfrow=c(1,2))
truehist(chain$alpha, col="blue", xlab="alpha", ylab="density")
curve(dlnorm(x, -1, 1), col ="red", add=T)
truehist(chain$sig2, col="blue", xlab="sig2", ylab="density")
curve(dlnorm(x, 2, 1), col ="red", add=T)
dev.off()

pdf("./bayou_results/bayOUpar_a2.pdf")
par(mfrow=c(1,2))
truehist(chain$alpha, col="blue", xlab="alpha", ylab="density", xlim=c(0,10), ylim=c(0,1))
curve(dlnorm(x, -1, 1), col ="red", add=T)
truehist(chain$sig2, col="blue", xlab="sig2", ylab="density", xlim=c(0,6), ylim=c(0,1))
curve(dlnorm(x, 2, 1), col ="red", add=T)
dev.off()

pdf("./bayou_results/bayOU_a.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain, burnin=0.3, lwd=2, edge.type="theta", pal=colorRampPalette(c("blue", "red")), show.tip.label=T, circle.col="black", cex=0.7)
dev.off()



fit2 <- bayou.mcmc(wiens.tre, res, SE=0.001, model="OU", prior, ngen=2000000, new.dir=getwd(), plot.freq=NULL, ticker.freq=10000)
chain2 <- load.bayou(fit2, save.Rdata=T, cleanup=T)
chain2 <- set.burnin(chain2, 0.3)
out2 <- summary(chain2)

write.csv(out2$statistics, "./bayou_results/bayOU_b_stats.csv")
write.csv(out2$branch.posteriors, "./bayou_results/bayOU_b_posteriors.csv")

pdf("./bayou_results/bayOUpar_b.pdf")
par(mfrow=c(1,2))
truehist(chain2$alpha, col="blue", xlab="alpha", ylab="density")
curve(dlnorm(x, -1, 1), col ="red", add=T)
truehist(chain2$sig2, col="blue", xlab="sig2", ylab="density")
curve(dlnorm(x, 2, 1), col ="red", add=T)
dev.off()

pdf("./bayou_results/bayOUpar_b2.pdf")
par(mfrow=c(1,2))
truehist(chain2$alpha, col="blue", xlab="alpha", ylab="density", xlim=c(0,10), ylim=c(0,1))
curve(dlnorm(x, -1, 1), col ="red", add=T)
truehist(chain2$sig2, col="blue", xlab="sig2", ylab="density", xlim=c(0,10), ylim=c(0,1))
curve(dlnorm(x, 2, 1), col ="red", add=T)
dev.off()

pdf("./bayou_results/bayOU_b.pdf")
par(mfrow=c(1,1))
plotSimmap.mcmc(chain2, burnin=0.3, lwd=2, edge.type="theta", pal=colorRampPalette(c("blue", "red")), show.tip.label=T, circle.col="black", cex=0.7)
dev.off()



###FIGURE S1a###
pdf("./bayou_results/bayOU_gelman.pdf")
par(mfrow=c(3,1))
RlnL <- gelman.R("lnL", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=666666.67, lty=2, col="red")
Ralpha <- gelman.R("alpha", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=666666.67, lty=2, col="red")
Rsig2 <- gelman.R("sig2", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=666666.67, lty=2, col="red")
dev.off()



L1 <- Lposterior(chain, wiens.tre, burnin=0.3)
L2 <- Lposterior(chain2, wiens.tre, burnin=0.3)



###FIGURE S1b###
pdf("./bayou_results/bayOU_posterior_diagnostics0120.pdf", 4, 4, useDingbats=F)
par(mfrow=c(1,1))
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)
dev.off()




pdf("./bayou_results/bayOU_pp_prior_a0120.pdf")
par(mfrow=c(1,1))
plot(density(L1[,1]), main="PP compared to prior")
abline(v=dsb(1, ntips = length(wiens.tre$tip.label), bmax = 1, prob = 1, log = FALSE), col="red", lty=2)
dev.off()

pdf("./bayou_results/bayOU_pp_prior_b0120.pdf")
par(mfrow=c(1,1))
plot(density(L2[,1]), main="PP compared to prior")
abline(v=dsb(1, ntips = length(wiens.tre$tip.label), bmax = 1, prob = 1, log = FALSE), col="red", lty=2)
dev.off()

#plotting the distribution of posterior probabilities (mean from both chains)
PP1 <- read.csv("./bayou_results/bayOU_a_posteriors.csv")
PP2 <- read.csv("./bayou_results/bayOU_b_posteriors.csv")
PP.all <- cbind(PP1[,2], PP2[,2])
PP <- apply(PP.all, 1, mean)

which(PP>0.8)

#77 141 154 178

topPP <- PP[which(PP>0.8)]
#only draw one arrow for the average of the top4



###FIGURE S1c###
pdf("./bayou_results/PP_distribution.pdf", 4, 4)
den <- density(PP, adjust=5); plot(den, xlab="Posterior Probability", main="")
cols <- rgb(200, 0,0, 150, maxColorValue=255)
polygon(den, col=cols)
abline(v=1/196, col="black", lty=2)
arrows(mean(topPP), 3, mean(topPP), 1, col="steel blue", lwd=3, length=.06)
dev.off()





