rm(list=ls())

########
# Xing et al., Extended Data Figure 1
########


##goal: determine whether phylogenetic covariance affects the scaling relation between orbit length and postnasal skull length
##PGLS (BM) with two different sets of 1000 trees
  #Stage2_MayrParSho_Ericson.nex
  #Stage2_MayrParSho_Hackett.nex

  #data: "skull.orbit.avg.csv"

#####Outline of the phylogenetic comparative analysis

## 1 Preliminaries: load libraries, load data and phylogenies
## 2 Fit conventional linear model
## 3 PGLS with Ericson trees
        #Perform PGLS with a BM correlation structure
        #Extract slopes and their standard errors
## 4 PGLS with Hackett trees
        #Perform PGLS with a BM correlation structure
        #Extract slopes and their standard errors
## 5 Summarize results over each tree set and superimpose the slope from the conventional linear fit



##### 1 Preliminaries

#please set working directory
setwd(path/to/your/working/directory)

#call libraries
require(ape)
require(geiger)
require(nlme)
require(phytools)

#load tree and data
phy.all.1 <- read.nexus("Stage2_MayrParSho_Ericson.nex") #reading in tree samples
phy.all.2 <- read.nexus("Stage2_MayrParSho_Hackett.nex") #reading in tree samples
sockets <- read.csv("skull.orbit.avg.csv", header=T)
rownames(sockets) <- sockets$taxon #assigning rownames (used for matching purposes)


##### 2 Fit conventional linear model
fit0 <- lm(log10(orbit_length)~log10(skull_length), data=sockets)
fit0.sum <- summary(fit0)
conv.slope <- fit0.sum$coefficients[2,1]
conv.slope
conv.slope.se <- fit0.sum$coefficients[2,2]
conv.slope.se

##### 3 PGLS with Ericson trees

#loop PGLS over trees and save slopes and their standard errors

#set up empty vectors for storing results

result_1_slope <- as.vector(rep(NA, 1000))
result_1_se <- as.vector(rep(NA, 1000))
result_1_lambda <- as.vector(rep(NA, 1000))

#start the loop
for (i in 1:1000) try({

  print(i)
  
  #match data with tree
  phy0 <- phy.all.1[[i]] 
  data <- sockets[phy0$tip.label,]
  compare.data <- treedata(phy0, data, sort=T, warnings=F) # the treedata function compares data with tree and prunes the tree automatically
  phy <- compare.data$phy
  
  ordered.sockets <- data.frame(taxon=as.character(compare.data$data[,1]), 
                                sk=as.numeric(compare.data$data[,3]),
                                ol=as.numeric(compare.data$data[,4]))
  rownames(ordered.sockets) <- ordered.sockets$taxon
  ordered.sockets <- ordered.sockets[phy$tip.label,]
  #head(ordered.sockets)
  
  #sort data
  variables <- ordered.sockets[,2:3]
  sk <- variables[,1]
  names(sk) <- phy$tip.label
  orb <- variables[,2]
  names(orb) <- phy$tip.label
  
  #round values
  sk <- signif(log10(sk), 3)
  orb <- signif(log10(orb),3)
  
  #putting it all in a dataframe
  M <- as.data.frame(cbind(sk, orb))
  M <- cbind(phy$tip.label, M)
  colnames(M)<-c("taxon","sk","orb")
  
  #PGLS
  fit1 <- try(gls(orb~sk, data=M, correlation=corPagel(0.5, phy, fixed=FALSE), method="ML"))
  
  sum.fit1 <- summary(fit1)
  result_1_slope[[i]] <- sum.fit1$tTable[[2]]
  result_1_se[[i]] <- sum.fit1$tTable[[4]]
  result_1_lambda[[i]] <- sum.fit1$modelStruct[[1]]


})



#function to visualize the slopes and their standard errors over the entire tree sample
slopesum <- function (slopes, se, minslope, maxslope, title, ...){
  llim <- slopes-se
  ulim <- slopes+se
  
  plot(slopes, as.vector(1:length(slopes)), type="n", main=title, xlab="slope estimates", ylab="tree sample #", xlim=(c(minslope, maxslope)))
  
  for (i in 1:length(slopes)){
    lines(x=c(llim[[i]], ulim[[i]]), y=c(i,i), col="light grey")
  }
  
  points(slopes, as.vector(1:length(slopes)), pch=4, cex=0.1)
  
}

mean(result_1_slope, na.rm=T)
range(result_1_slope, na.rm=T)

#plot the slope estimates in comparison to isometery and conventional slope
slopesum(result_1_slope, result_1_se, 0.7, 1.1, title="Slope estimates")
abline(v=1, lty=2, col="red") #isometry
abline(v=conv.slope, lty=2, col="blue") #slope estimate from conventional statistics



##### 4 PGLS with Hackett trees

#loop PGLS over trees and save slopes and their standard errors

#set up empty vectors for storing results

result_2_slope <- as.vector(rep(NA, 1000))
result_2_se <- as.vector(rep(NA, 1000))
result_2_lambda <- as.vector(rep(NA, 1000))

#start the loop
for (i in 1:1000) try({
  
  print(i)
  
  #match data with tree
  phy0 <- phy.all.2[[i]] 
  data <- sockets[phy0$tip.label,]
  compare.data <- treedata(phy0, data, sort=T, warnings=F) # the treedata function compares data with tree and prunes the tree automatically
  phy <- compare.data$phy
  
  ordered.sockets <- data.frame(taxon=as.character(compare.data$data[,1]), 
                                sk=as.numeric(compare.data$data[,3]),
                                ol=as.numeric(compare.data$data[,4]))
  rownames(ordered.sockets) <- ordered.sockets$taxon
  ordered.sockets <- ordered.sockets[phy$tip.label,]
  #head(ordered.sockets)  
  
  #sort data
  variables <- ordered.sockets[,2:3]
  sk <- variables[,1]
  names(sk) <- phy$tip.label
  orb <- variables[,2]
  names(orb) <- phy$tip.label
  
  #round values
  sk <- signif(log10(sk), 3)
  orb <- signif(log10(orb),3)
  
  #putting it all in a dataframe
  M <- as.data.frame(cbind(sk, orb))
  M <- cbind(phy$tip.label, M)
  colnames(M)<-c("taxon","sk","orb")
  
  #PGLS
  fit2 <- try(gls(orb~sk, data=M, correlation=corPagel(0.5, phy, fixed=FALSE), method="ML"))
  
  sum.fit2 <- summary(fit2)
  result_2_slope[[i]] <- sum.fit2$tTable[[2]]
  result_2_se[[i]] <- sum.fit2$tTable[[4]]
  result_2_lambda[[i]] <- sum.fit2$modelStruct[[1]]
})


mean(result_1_slope, na.rm=T)
range(result_1_slope, na.rm=T)

mean(result_2_slope, na.rm=T)
range(result_2_slope, na.rm=T)

mean(result_1_lambda, na.rm=T)
mean(result_2_lambda, na.rm=T)

#plot the slope estimates in comparison to isometery and conventional slope
slopesum(result_2_slope, result_2_se, 0.7, 1.1, title="Slope estimates")
abline(v=1, lty=2, col="red") #isometry
abline(v=conv.slope, lty=2, col="blue") #slope estimate from conventional statistics




##### 5 Visualize all results in one combined figure

pdf("Extended_Data_Fig_09.pdf", useDingbats = FALSE, 4, 8)

par(mfrow=c(2,1))


#plot the slope estimates in comparison to isometery and conventional slope
slopesum(result_1_slope, result_1_se, 0.7, 1.1, title="Slope estimates, Ericson trees")
abline(v=1, lty=2, col="red") #isometry
abline(v=conv.slope, lty=2, col="blue") #slope estimate from conventional statistics

#plot the slope estimates in comparison to isometery and conventional slope
slopesum(result_2_slope, result_2_se, 0.7, 1.1, title="Slope estimates, Hackett trees")
abline(v=1, lty=2, col="red") #isometry
abline(v=conv.slope, lty=2, col="blue") #slope estimate from conventional statistics


dev.off()
