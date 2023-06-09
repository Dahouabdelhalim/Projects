

# This code allows to find out whether a set of samples form well-defined groups in 
# a multidimensional space and how well-defined these groups are. The script also includes 
# a small code (1) with which to generate artificial groups in a multidimensional space
# with which to compare the results of applying the AMD index (2) to the target dataset.


#  1. ------------ Generating artificial groups in a multidimensional space ----------

# with which to compare the results of applying the AMDi (below) to the target dataset


#----- Initiate Variables

ng = 6 # Number of groups to be generated
nd = 9 # Number of dimensions of the space
ns = 8000 # Total Nº of samples
dg = 100 # Side of the multidimensional space 
sd = 5 # Standard deviation of the groups, which represents their degree of compactness.

fg<-c()

for(v in 1:ng){ 
  
  cld <-c()
  
  for (l in 1:nd){ 
    
    co<-runif(1,0,dg) # random coordinates of the centroid of each cluster
    kk<-abs(round(rnorm(ns/ng, co, sd),0)) # generation of each of the groups
    cld <- cbind(cld,kk) # generating as many columns/variables as dimensions of the space (nd)
  }
  
  fg <- rbind(fg, cld) # union of all groups in the same data base
}

Data<-as.data.frame(fg) 



# 2. ------------- AMD index  ------------


#----- Loading Libraries

library(e1071())

# For applying the AMD index to artificial groups: 

Data<-as.data.frame(fg) 

# For applying the AMD index to the the species-level trophic space, open Data_I and copy all excluding TG

# For applying the AMD index to the the community-level trophic space, open Data_II and copy from SIn to Grn


read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\\t",header=header,...)
}

Data=read.excel()


#----- Initiate Variables

its <- 10 # number of times the analyses are repeated, to avoid the effect of chance on clustering
nin <- 2 # Minimun number of user-defined de clusters
nsp <- 15 # Míaximun number of user-defined de clusters
r   <- nsp-nin
nc  <- ncol(Data)
nr  <- nrow(Data)
Maxpms <- c()

#----- Model
# to repeat the analyses several times (to avoid the effect of randomness on the clustering).
for (k in 1:its) { 
  Maxpm <- c()
  
  # to apply the analyses to each number of user-defined clusters.
  for (i in nin:nsp) { 
    # fuzzy clustering with different No of groups, from nin to nsp.
    cl <- cmeans(Data, i, 20, verbose = FALSE, method = "cmeans", m = 2)
    # matrix with the probabilities of each sample belonging to each group
    Prob <- cl$membership 
    
    # Calculate the probability of each sample belonging to the group to which it was assigned in a cluster.
    Maxprob <- apply(Prob, 1, max)
    
    # calculate the average of the maximum probabilities obtained in the previous step.
    Mpm<-mean(Maxprob)-1/i # and the inverse of the number of clusters is subtracted.
    # create a df with the Mpm value of each of the analyses, with different No of groups, nin-to nsp     
    Maxpm <- cbind(Maxpm, Mpm)
    
    N <- ncol(Maxpm)
    # names the number of clusters to which each value of CRc corresponds.
    colnames(Maxpm)[N] <- paste(i, "Cls", sep = "") 
  }
  
  # creates a df with the Maxpm values of each itr
  Maxpms <- rbind(Maxpms, Maxpm) 
  # to indicate, while the script is running, which itr number it is going through.
  print(k)
}

colMax <- function(colData) { apply(colData, MARGIN = c(2), max) }

# calculates the maximum of all CRs values obtained for each number of clusters.
Maxpmean              <- data.frame(colMax(Maxpms)) 
colnames(Maxpmean)[1] <- "Maxpmean"
clusters              <- c(nin:nsp)
Results               <- cbind(clusters,Maxpmean)

#----- Plot
library(ggplot2)
res_df <-  data.frame(
  x = nin:nsp,
  y = Results$Maxpmean
)
windows(); ggplot(res_df, aes(x, y )) +
  geom_point() +
  labs(x = "N. of Clusters", y = "Mean of max probabilities") +
  theme_bw()
  
  