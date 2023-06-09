library(phytools) # package with tools for reconstruction and Pagel's lambda
library(ggtree) #package with tools for creating and displaying dendrograms
library(ggplot2) #package for graph manipulation utilities
library(ggimage)
library(coda) #package for MCMC diagnostics
library(matrixStats) #package for calculating row means (not in base)
library(dplyr)#for a short filter used to remove outgroup
library(geiger)

#Reading in the species data#
Rallid_CSV<-read.csv("Data/RallidLogCertain.csv",header=TRUE,row.names=1)
Rallid_CSV
head(Rallid_CSV)

Rallids_Named<-read.csv("Data/RallidLogCertain.csv")
str(Rallids_Named)

#Reading in the rooted tree#
RallidTree<-read.tree("Data/Rallid_Tree.tre")
RallidTree

#Making a slightly different tree to remove species without duet information
CertainTree <- drop.tip(RallidTree, c("Rallina_canningi", "Laterallus_levraudi",
                                      "Laterallus_ruber", "Habroptila_wallacii",
                                      "Hypotaenidia_owstoni", 
                                      "Pardirallus_maculatus", 
                                      "Hapalocrex_flaviventer", 
                                      "Laterallus_spilopterus", "Porzana_fluminea",
                                      "Zapornia_akool", "Zapornia_bicolor",
                                      "Paragallinula_angulata",
                                      "Tribonyx_ventralis", "Fulica_cornuta"))
CertainTree
plotTree(CertainTree,ftype="i",fsize=0.6,lwd=1)
Ntip(CertainTree)
print(CertainTree$edge.length)

###################### make vector data for computational reasons #####################################
duetDataV <- as.data.frame(Rallid_CSV$Duets) #defines vector of states for the MCMC
rownames(duetDataV) <- Rallids_Named$Species
duetDataV
name.check(CertainTree,duetDataV)
duetDataV <-treedata(CertainTree, duetDataV)
duetDataV$phy
duetDataV$data
duetData_phy <- duetDataV$phy #assign phylogeny to object
duetData_dat <- as.character(duetDataV$data) #assign data to character vector
names(duetData_dat) <- rownames(duetDataV$data) #give names to entries from original data
duetData_dat
duetData_matrix <- as.matrix(duetData_dat)
duetData_matrix

############################functions for getting the calculated tip proportions#############################
estimateTip<-function(index)
{length(tip_liab[index][tip_liab[index]>0])/nrow(tip_liab[index])}
estimateTip2<-function(index)
{1-length(tip_liab[index][tip_liab[index]>0])/nrow(tip_liab[index])}
#not a great solution but I needed a quick way to get the proportions for the tips, 
#basically liability >0 / total liability

########################### Reconstruct notes (the MCMC)#######################################
ngen <- 50000000 #number of generations to run
burnin <- ngen*.2 #approximate number of initial generations to eliminate 
sample <- 100 #sample posterior distribution every 1000 samples

AncDuet <- ancThresh(duetData_phy, duetData_dat[1:length(duetData_phy$tip.label)], 
              ngen = ngen, burnin = burnin, model = "BM", sequence = c("0", "1"),
              control = list(sample = sample, tipcol = "estimated"), 
              label.offset = 0.01, plot = TRUE)

########################### Metrics function & calculation#####################################
metricsBot <- function(ancObj){
  test <- as.mcmc(ancObj$liab[-(1:burnin/sample),])#gets the post-burn in liability
  size = effectiveSize(test) #gets effective sample size 
  gew = geweke.diag(test, frac1 = 0.1, frac2 = 0.5) #performs Geweke's diagnostic to get z scores
  gewz <- gew$z[!is.na(gew$z)]
  gew_p = pnorm(-abs(gewz)) #creates p values for the z scores
  
  #this next bit is supposed to calculate diff in means over std, see below
  sig = names(gew_p[which(gew_p>0.05)]) #indexes the significant nodes
  test_sig <- test[,sig] 
  means1 <- colMeans(test_sig[1:(length(test_sig[,1])*0.1),]) 
  means2 <- colMeans(test_sig[(length(test_sig[,1])*0.5):(length(test_sig[,1])),])
  std_all <- colSds(test_sig)
  x <- as.data.frame((means1-means2)/std_all)
  colnames(x) <- 'stdiff'
  p <- ggplot(data = x, aes(x$stdiff)) + xlim(-0.5, 0.5) + geom_histogram(binwidth=.01, color="black", fill="grey") + 
    xlab("Standardized mean difference") +
    ylab("Number of variables")
  p
  ret <- list(sig,size)
  return(ret)
}

#Duet metrics 
test <- as.mcmc(AncDuet$liab[-(1:burnin/sample),])#gets the post-burn in liability
summary(test)
plot(test)
size <- effectiveSize(test) #gets effective sample size 
size
geweke.diag(test, frac1 = 0.1, frac2 = 0.5)#performs Geweke's diagnostic to get z scores
geweke.plot(test, frac1 = 0.1, frac2 = 0.5, pvalue = 0.05)

met1 <- metricsBot(AncDuet)
met1
summary(met1)
ggsave("Output/Duet_Means_Certain.pdf", width = 7, height = 7, limitsize=FALSE)

########################################Correlation################################################

#important that we don't have unknown tips in the tree, because they are inferred from other tips: pseudoreplication issue

#defines known duet priors
Kprior_duet <- as.numeric(Rallids_Named$Duet)
names(Kprior_duet) <- Rallids_Named$Species
#gets phylogenetically independent values
cont_duet <- pic(Kprior_duet, CertainTree) 

#####################################Dendrograms###############################################
plot(AncDuet$par[, "logLik"], type = "l", xlab = "generation", ylab = "logLik")
#Converged, shows a fuzzy caterpillar trace
#Get the estimates of thresholds of our data
colMeans(AncDuet$par[(0.2 * ngen/sample):(ngen/sample) + 1, c("0", "1")])

duetData <- cbind(Rallids_Named[3])
numtip <- length(rownames(duetData))
numnode <- CertainTree$Nnode
maxnode <- numtip+numnode

#Duet tree
if(TRUE){
  scoreData <- duetData_phy
  char <- "Duetting"
  nonchar <- "Non-Duetting"
  tip_liab <- AncDuet$liab[1:numtip]
  tip_purp <- lapply(1:numtip,FUN=estimateTip)
  Chart <- as.data.frame(as.matrix(tip_purp))
  Chart$V2 <- 0
  rownames(Chart) <- rownames(scoreData)
  colnames(Chart) <- c(nonchar, char)
  Chart[,char] <- 1-as.numeric(unlist(Chart[,nonchar]))
  Chart[,nonchar] <- as.numeric(Chart[,nonchar])
  nodelist <- data.frame(rbind(AncDuet$ace))
  nodelist$node <- c(1:numtip, (numtip+1):maxnode)
  write.csv(nodelist, "Output/Duet_Nodes_Certain.csv") #for referencing specific node values
  #begin creating ggobject for tree
  dendroplot <- ggtree(duetData_phy, size=0.4) + geom_tiplab(size=1.7, align=TRUE, linesize=0.5, offset=2) + xlim(0,110)
  Npies <- nodepie(nodelist, cols=1:2, alpha=0.8,  color = c("red", "blue") ) 
  dendroplot <- (dendroplot + ggplot2::scale_y_continuous(expand=c(0, 0.3))) %>%
    gheatmap(scoreData[1], offset=0, width=0.015, low="blue", high="red", colnames=FALSE, color="black")
  dendroplot <- dendroplot + theme(legend.position="left") + ggtitle("Duet Tree")
  dendroplot <- inset(dendroplot, Npies, width = 0.01, height = 0.01)
}
plot(dendroplot)

#node reference, need to correspond to the proportions in the CSV files above
if(TRUE){
  nodetree <- RallidTree
  nodetree$node.label <- (numtip+1):maxnode
  dendroplot <- ggtree(nodetree, size=0.4) + geom_tiplab(size=3, align=TRUE, linesize=0.5) + xlim(0,110)
  dendroplot <- dendroplot + geom_nodelab(hjust = 1.2)
  dendroplot <- dendroplot + theme(legend.position="left") + ggtitle("Duet Tree")
  ggsave(plot = dendroplot, "Output/Node_Reference_Certain.pdf", width = 30, height = 30, units = "cm", limitsize = FALSE)
}

#make a histogram of Geweke's test diagnostics of stability and convergence
Geweke_Z<-read.csv("Output/Duet_Geweke_Diagnostics.csv")
Hist<-ggplot(Geweke_Z, aes(x=GewekeZ)) + xlim(-6, 6) +
  geom_histogram(color="black", fill="gray", bins=15, binwidth=0.15, size=1.25)
Hist