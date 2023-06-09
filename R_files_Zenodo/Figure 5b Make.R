## R code to produce Figure 5b in:
## Johnston EC, Cunning R, Burgess SC. Cophylogeny and specificity between cryptic coral species (Pocillopora spp.) at Moâ€™orea and their symbionts (Symbiodiniaceae).
## Code written by Erika Johnston. February 2022. Send comments or corrections to ejohnston@bio.fsu.edu
## R version 4.1.2 (2022-03-01)

# initialization ----------------------------------------------------------
rm(list=ls()) #clear all variables

# Load libraries
library(vegan)
library(ape)
library(paco)
library(ggplot2)
library(dplyr)
library(stringr)

### Establish PACo function: adjustment prior to procustes analysis
PACo <- function(H.dist, P.dist, HP.bin)
{ HP.bin <- which(HP.bin > 0, arr.in=TRUE)
  H.PCo <- pcoa(H.dist, correction="cailliez") #Performs PCo of Host distances 
  P.PCo <- pcoa(P.dist, correction="cailliez") #Performs PCo of Parasite distances
  if (is.null(H.PCo$vectors.cor)==TRUE) H.PCo <- H.PCo$vectors else
    H.PCo <- H.PCo$vectors.cor      # returns corrected pcoord 
  if (is.null(P.PCo$vectors.cor)==TRUE) P.PCo <- P.PCo$vectors else
    P.PCo <- P.PCo$vectors.cor
  H.PCo <- H.PCo[HP.bin[,1],]  #adjust Host PCo vectors 
  P.PCo <- P.PCo[HP.bin[,2],]  #adjust Symbiont PCo vectors
  list (H.PCo = H.PCo, P.PCo = P.PCo)}

#Read in association data matrix
HP <- read.csv("Figure 5b.csv", header = FALSE)
HP <- data.frame(HP[,-1], row.names = HP[,1])
names(HP) <- as.matrix(HP[1, ])
HP <- HP[-1, ]
HP[] <- lapply(HP, function(x) type.convert(as.character(x)))

### 2. DATA INPUT
#2.1 Phylogenetic trees:
#Read in Pocillopora host tree
TreeH <- read.nexus("Figure 5b_Pocillopora.nex")

#Read in Cladocopium host tree
TreeP <- read.nexus("Figure 5b_Cladocopium.nex")

host.D <- cophenetic(TreeH)
para.D <- cophenetic(TreeP)

#Sort host and parasite taxa in distance matrices to match the HP matrix:
host.D <- host.D[rownames(HP), rownames(HP)]
para.D <- para.D[colnames(HP), colnames(HP)]

### 3. APPLY PACo FUNCTION  
PACo.fit <- PACo(host.D, para.D, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) #Procrustes Ordination 
NLinks = sum(HP) #Number of H-P links; needed for further computations

#3.1 Goodness-of-fit-test --- ## Takes a long time to run
m2.obs <- HP.proc$ss #observed sum of squares
N.perm = 100000 #set number of permutations for testing
P.value = 0
set.seed(8765) ### use this option to obtain reproducible randomizations
for (n in c(1:N.perm))
{ 
  if (NLinks <= nrow(HP) | NLinks <= ncol(HP)) 	#control statement to avoid all symbionts being associated to a single host 
  {	flag2 <- TRUE 
  while (flag2 == TRUE)	{ 
    HP.perm <- t(apply(HP,1,sample))
    if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
  }  
  } else { HP.perm <- t(apply(HP,1,sample))} #permutes each HP row independently
  PACo.perm <- PACo(host.D, para.D, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
  if (m2.perm <= m2.obs)
  {P.value = P.value + 1} 
}
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\\n", "P-value = ", P.value, " based on ", N.perm," permutations.")

#3.2 Contribution of individual links
HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)# empty matrix of jackknifed squared residuals
colnames (SQres.jackn) <- paste(rownames(HP.proc$X),rownames(HP.proc$Yrot), sep="-") #colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) #Needed to compute 95% confidence intervals.
for(i in c(1:NLinks)) #PACo setting the ith link = 0
{HP.ind <- HP
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(host.D, para.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo) 
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind	#Append residuals to matrix of jackknifed squared residuals
} 
SQres.jackn <- SQres.jackn**2 #Jackknifed residuals are squared
SQres <- (residuals (HP.proc))**2 # Vector of original square residuals
#jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) #apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) #mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) #standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) #upper 95% confidence interval

#convert vectors into dataframes
phi.UCI.df <- as.data.frame(phi.UCI)
phi.mean.df <- as.data.frame(phi.mean)
phi.UCI.df <- cbind(UCI_links = rownames(phi.UCI.df), phi.UCI.df)
rownames(phi.UCI.df) <- 1:nrow(phi.UCI.df)
phi.mean.df <- cbind(mean_links = rownames(phi.mean.df), phi.mean.df)
rownames(phi.mean.df) <- 1:nrow(phi.mean.df)

#Merge dataframes
phi.UCI.mean = merge(phi.UCI.df, phi.mean.df, by.x=c("UCI_links"), by.y=c("mean_links"))
#Add new column for Pocillopora species/haplotype
phi.UCI.mean$Pocillopora <- phi.UCI.mean$UCI_links
phi.UCI.mean$Pocillopora <- gsub("\\\\-.*","",phi.UCI.mean$Pocillopora)

#Merge with Cladocopium clades
clades <- read.csv("Figure 5b - Clad clades.csv")
PACo_clades = merge(phi.UCI.mean, clades, by.x=c("UCI_links"), by.y=c("UCI_links"))
#Add new column for Cladocopium clades
PACo_clades$Clade_assign <- PACo_clades$Clad_clades
PACo_clades$Clade_assign <- gsub(".*-","",PACo_clades$Clade_assign)

#Find mean of phi.mean
median(phi.UCI.mean$phi.mean) # median = 1.195144e-05

PACo_clades$Pocillopora <- factor(PACo_clades$Pocillopora,levels=c("Haplotype_10",
                                                                   "P_verrucosa",
                                                                   "Haplotype_8a",
                                                                   "P_meandrina",
                                                                   "P_eydouxi",
                                                                   "Haplotype_2",
                                                                   "Haplotype_11"))

## Plot
ggplot(PACo_clades, aes(x=reorder(Clad_clades, phi.mean), y=phi.mean, fill=Pocillopora))+
  geom_bar(stat='identity')+
  geom_hline(yintercept=1.195144e-05, linetype='dashed', col = 'red')+
  theme_minimal(base_size=15)+
  geom_errorbar(aes(ymin=phi.mean, ymax=phi.UCI), width=.1)+ 
  scale_fill_manual(name="Pocillopora species/haplotype",
                    values = c("Haplotype_10"= "#D55E00",
                               "P_verrucosa"= "#E69F00",
                               "P_meandrina"= "#0072B2",
                               "Haplotype_8a"= "#CC79A7",
                               "P_eydouxi"= "#56B4E9",
                               "Haplotype_11"= "#009E73",
                               "Haplotype_2"= "#B2DF8A"))+
  coord_flip()+
  theme(axis.text=element_text(size=8))+
  labs(title="",
       x ="Pocillopora - Cladocopium links", y = "Squared residuals")
