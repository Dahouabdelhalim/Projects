setwd("~/Documents/Gray lab/Negation/2019_03_29_plotting copy 2021 01 16 for zenodo")

dat <- read.csv("NegEx_IEdata.txt", sep = "\\t", header = T)

library(ape)
library(phytools)
library(corHMM)

IE_tree <- read.nexus("IE.txt6.nex")

IE_tree <- IE_tree[[3]]
plot(IE_tree)

IE_tree2 <- compute.brlen(IE_tree)
IE_tree2 <- compute.brlen(IE_tree, 1)
plot(IE_tree2)

# match order of languages so tiplabels doesn't screw up

match <- match(IE_tree$tip.label, dat$NAME)
dat2 <- dat[,][match,]

#data file
# for ASE and node labels
dat_mat <- matrix(ncol = 2, nrow = 55, NA)
rownames(dat_mat) <- as.character(dat2[,1])
dat_mat[,1] <- as.character(dat2[,1])
dat_mat[,2] <- as.character(dat2[,11])
unique(dat_mat[,2])

# for tip labels
dat_mat_tips <- matrix(ncol = 6, nrow = 55, NA)
rownames(dat_mat_tips) <- as.character(dat2[,1])
colnames(dat_mat_tips) <- c("V1","V2","V3","V4","V5","V6")

dat_mat_vect <- dat2[,11]

for(i in 1:nrow(dat_mat_tips)){
  if(dat_mat_vect[i] == "1"){
    dat_mat_tips[i,1] <- 1
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 0
    dat_mat_tips[i,4] <- 0
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 0
  }
  if(dat_mat_vect[i] == "1&3"){
    dat_mat_tips[i,1] <- 1
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 1
    dat_mat_tips[i,4] <- 0
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 0
  }
  if(dat_mat_vect[i] == "3"){
    dat_mat_tips[i,1] <- 0
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 1
    dat_mat_tips[i,4] <- 0
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 0
  }
  if(dat_mat_vect[i] == "4"){
    dat_mat_tips[i,1] <- 0
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 0
    dat_mat_tips[i,4] <- 1
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 0
  }
  if(dat_mat_vect[i] == "6"){
    dat_mat_tips[i,1] <- 0
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 0
    dat_mat_tips[i,4] <- 0
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 1
  }
  if(dat_mat_vect[i] == "1&4"){
    dat_mat_tips[i,1] <- 1
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 0
    dat_mat_tips[i,4] <- 1
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 0
  }
  if(dat_mat_vect[i] == "2"){
    dat_mat_tips[i,1] <- 0
    dat_mat_tips[i,2] <- 1
    dat_mat_tips[i,3] <- 0
    dat_mat_tips[i,4] <- 0
    dat_mat_tips[i,5] <- 0
    dat_mat_tips[i,6] <- 0
  }
  if(dat_mat_vect[i] == "3&5"){
    dat_mat_tips[i,1] <- 0
    dat_mat_tips[i,2] <- 0
    dat_mat_tips[i,3] <- 1
    dat_mat_tips[i,4] <- 0
    dat_mat_tips[i,5] <- 1
    dat_mat_tips[i,6] <- 0
  }
}



ASE <- rayDISC(data = dat_mat, phy = IE_tree2, ntraits = 1, charnum = 1, model = c("ER"))

col5 <- c("blue", "purple","hotpink","red","orange","yellow")

plotTree(IE_tree2, setEnv = TRUE, offset = 0.5, fsize= 0.8)
tiplabels(pie = dat_mat_tips, piecol = col5, cex = 0.4)
nodelabels(pie = ASE$states, piecol = col5, cex = 0.5)

legend("topleft", c("A","AB","B","BC","C","CA"), fill = col5)



