#bee vac network data

library(bipartite)
library(plyr)
library(vegan)

#Note:  This table is then imported into R and aggregated into one or more webs (or an array or a list of webs) using frame2webs. 

#change file name depending on what year you are analyzing. Note: columns must be "lower", "higher", "webID" and "freq"
webData <- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/2015NetworkFreq.csv",header = T)
webData <- as.data.frame(webData)
names(webData) <- c("Lower", "Higher", "webID","Freq")

a <- frame2webs(webData,varnames = c("Lower", "Higher", "webID", "Freq"), type.out = "list")

bee.indices <- list()
outData <- list()
for (i in c(1:length(a))){
  bee.indices[[i]] <- grouplevel(a[[i]],index = c("generality", "mean number of shared partners", "niche overlap"), level = "higher")
  
  singleWeb <- a[[i]]
  singleWeb_groupIndices <- grouplevel(a[[i]],index = c("generality", "mean number of shared partners", "niche overlap"), level = "higher")
  
  # generate null models and calculate their index values:
  singleWeb.nulls <- nullmodel(singleWeb, "r2dtable")
  
  #calc indices for null models -- time intensive (save)
  singleWeb.null.res.group <- lapply(singleWeb.nulls, function(x) grouplevel(x, index = c("generality", "mean number of shared partners", "niche overlap"), level = "higher"))
  
  #calc z scores for group level
  z.mat.group <- matrix(0, 1, 3)
  colnames(z.mat.group) <- names(singleWeb.null.res.group[[1]])
  
  for (j in 1:3){
    mean.n <- apply(as.matrix(sapply(singleWeb.null.res.group, function(x) x[j])), 2, mean, na.rm=T)
    sd.n <- apply(as.matrix(sapply(singleWeb.null.res.group, function(x) x[j])), 2, sd, na.rm=T)
    z <- (singleWeb_groupIndices[j] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are constant (i.e. mean.n==index);then the z-score should be 0
    z.mat.group[1,j] <- z
  }

  z.mat.group2 <- as.matrix(z.mat.group)
  z.all <- as.data.frame(cbind(Site = i, Metric = "z score", z.mat.group2))
  
  # calculate p-values: for group level
  p.mat.group <- matrix(0, 1, 3)
  colnames(p.mat.group) <- names(singleWeb.null.res.group[[1]])
  
  #rownames(p.mat) <- colnames(exa)
  for (k in 1:3){
    null.vals <- sapply(singleWeb.null.res.group, function(x) x[k])
    p2 <- 1:3
    p.lower <- null.vals < singleWeb_groupIndices[k]
    p.higher <- null.vals > singleWeb_groupIndices[k]
    p.equal <- null.vals == singleWeb_groupIndices[k]
    sample.size <- NROW(null.vals) - max(sum(is.na(null.vals)), sum(is.na(p.lower)),
                                         sum(is.na(p.higher)), sum(is.na(p.equal)))
    #p2[j] <- sum(p.higher)/sample.size
    p2 <- (min(sum(p.lower, na.rm=T), sum(p.higher, na.rm=T))+sum(p.equal, na.rm=T))/sample.size #turn this into probabilities, acknowleding possible NAs (in NSI!)
    p.mat.group[k] <- p2 #ifelse(p2>0.5, 1-p2, p2)*2 #two-tailed test
  }
  p.mat.group2 <- as.matrix(p.mat.group)
  p.all <- as.data.frame(cbind(Site = i, Metric = "p value", p.mat.group2))
  
  z.all2 <- as.data.frame(z.mat.group2)
  p.all2 <- as.data.frame(p.mat.group2)
  
  outData[[i]] <- rbind(z.all, p.all)
  
}

#combine all lists into a single dataframe
singleDF <- rbind.fill(outData)

#add name of site
names.list <- names(a)
repeat.names <- rep(names.list, each = 2)

singleDF$Name <- repeat.names

#change name based on year
write.csv(singleDF,file=paste0("2015Indicies_zAndp_20190502v2.csv"))

##############################################################################

#CREATION OF NETWORK GRAPHS

#load libraries
library(bipartiteD3)
library(r2d3) #allows for saving D3 webs as png
library(webshot) #allows for saving D3 webs as png

QuadData <- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/AllFourWebs.csv")
names(QuadData) <- c("higher","lower","webID", "freq")
bipartite::frame2webs(QuadData)-> QuadWeb

#change which web you want by changing the # in the double brackets (1-4). Also rename the site.

ManualColoursa<- c(Chicory='lightgray', Everlasting.Pea='darkgrey', Narrow.Leaf.Plantain='gray', Other= 'dimgray', Red.Clover="black")
a<-bipartite_D3(data =QuadWeb[[1]], PrimaryLab = 'Flowers',
                SecondaryLab = 'Pollinators',colouroption = 'manual',
                NamedColourVector = ManualColoursa, ColourBy = 1,
                SiteNames = 'June Vacant Lot', filename = 'demo1')


ManualColoursb<- c(Chicory='lightgray', Other='darkgrey', Other.Planted='gray', Queen.Annes.Lace= 'dimgray', Red.Clover="black")
b<-bipartite_D3(data =QuadWeb[[2]], PrimaryLab = 'Flowers',
                SecondaryLab = 'Pollinators',colouroption = 'manual',
                NamedColourVector = ManualColoursb, ColourBy = 1,
                SiteNames = 'June Pocket Prairie', filename = 'demo1')

ManualColoursc<- c(Birdsfoot.Trefoil='lightgray', Chicory='darkgrey', Other='gray', Queen.Annes.Lace= 'dimgray', Red.Clover="black")
c<-bipartite_D3(data =QuadWeb[[3]], PrimaryLab = 'Flowers',
                SecondaryLab = 'Pollinators',colouroption = 'manual',
                NamedColourVector = ManualColoursc, ColourBy = 1,
                SiteNames = 'August Vacant Lot', filename = 'demo1')

ManualColoursd<- c(Chicory='lightgray', Other='darkgrey', Other.Planted='gray', Queen.Annes.Lace= 'dimgray', Red.Clover="black")
d<-bipartite_D3(data =QuadWeb[[4]], PrimaryLab = 'Flowers',
                SecondaryLab = 'Pollinators',colouroption = 'manual',
                NamedColourVector = ManualColoursd, ColourBy = 1,
                SiteNames = 'August Pocket Prairie', filename = 'demo1')
