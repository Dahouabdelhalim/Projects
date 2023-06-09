library(TESS)
library(ape)

source("SelectModel.R")

# get the data set
tree <- read.nexus("data/Snakes.tre")

# some parameters
REPS <- 10
POSTERIOR_SIMS <- 10000
sampledSpecies <- tree$Nnode + 1
total <- 3500
missingSpecies <- total - sampledSpecies

# for each data set
name <- "Snakes"

cat("Running analysis for",name,".\\n")
    
# print the latex table header
cat("Model & log-Likelihood & AICc & BIC & p-value(gamma) & p-value(taxa) & p-value(treeheight)\\\\\\\\\\n",file=sprintf("../results/%s.txt",name),append=FALSE)
    
rho <- (sampledSpecies/(missingSpecies+sampledSpecies))

analysis <- selectModel(tree, rho, REPS, PosteriorPreditiveTesting = TRUE, POSTERIOR_SIMS = POSTERIOR_SIMS)
aicc <- analysis[[1]]
bic <- analysis[[2]]
logL <- analysis[[3]]
p.gamma <- analysis[[4]]
p.taxa <- analysis[[5]]
p.height <- analysis[[6]]


###
best_aicc <- which.min(aicc)
best_bic <- which.min(bic)

cat("Best AICc:\\t",best_aicc,"\\n")
cat("Best BIC:\\t",best_bic,"\\n")

for (k in 1:12) {
  cat(k,"&",logL[k],"&",aicc[k],"&",bic[k],"&",p.gamma[k],"&",p.taxa[k],"&",p.height[k],"\\\\\\\\\\n",file=sprintf("../results/%s.txt",name),append=TRUE)
}
    
