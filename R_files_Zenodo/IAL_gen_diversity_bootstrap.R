rm(list=ls(all=T))
setwd('~/R')

require('reshape2')
source('IAL_popsim.R')

f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# bootstrap analysis including Ontario population

dat <- read.csv('IAL_genotypes.csv',header = T)
nregion <- 8
pops <- unique(dat$Site)
datpop <- list()
nj <- matrix(0,nregion,1)
for (j in 1:nregion) {
  datpop[[j]] <- dat[dat$Site==pops[j],]
  nj[j] <- nrow(datpop[[j]])
}

nrep <- 1000
nloci <- 10
regions <- unique(dat$Site)
#minsize <- min(tapply(dat$birdid,dat$Site,max))
popsize <- 13

probminshannon <- matrix(0,3,2)
probminehet <- matrix(0,3,2)
probminallele <- matrix(0,3,2)
nallelemn <- matrix(0,nregion,3)
nallelesd <- matrix(0,nregion,3)
rownames(probminshannon) <- c('9','10','11')
colnames(probminshannon) <- c('jacknife','bootstrap')
rownames(probminehet) <- c('9','10','11')
colnames(probminehet) <- c('jacknife','bootstrap')
rownames(probminallele) <- c('9','10','11')
colnames(probminallele) <- c('jacknife','bootstrap')
rownames(nallelemn) <- pops
colnames(nallelemn) <- c('9','10','11')
rownames(nallelesd) <- pops
colnames(nallelesd) <- c('9','10','11')



popsizes <- c(9,10,11)
boots <- c(F,T)
sites <- as.character(unique(dat$Site))

for (i1 in 1:3) {
  for (i2 in 1:2) {
    stats <- popsim(datpop,nj,popsizes[i1],boots[i2],nrep,nregion,nloci)
    shannon <- stats[[1]]
    ehet <- stats[[2]]
    nallele <- stats[[3]]
    ehet2 <- 1-1/nloci*ehet
    shmin <- apply(shannon,1,min)
    ismin <- (shannon[,1]==shmin)
    probminshannon[i1,i2] <- mean(ismin)
    hetmin <- apply(ehet2,1,min)
    ismin <- (ehet2[,1]==hetmin)
    probminehet[i1,i2] <- mean(ismin) 
    nallelemin <- apply(nallele,1,min)
    ismin <- (nallele[,1]==nallelemin)
    probminallele[i1,i2] <- mean(ismin)
    
    if (boots[i2]) {
      nallelemn[,i1] <- apply(nallele/nloci,2,mean)
      nallelesd[,i1] <- apply(nallele/nloci,2,sd)
    }
    
    shndat <- data.frame(shannon)
    names(shndat) <- sites
    shndat$id <- seq(nrep)
    shndat <- melt(shndat,id='id')
    ggplot(shndat,aes(y=value,x=variable)) + stat_summary(fun.data = f, geom="boxplot")
  }
}

probminshannon
probminehet
probminallele
nallelemn
nallelesd

# bootstrap analysis without Ontario

dat <- read.csv('IAL_genotypes.csv',header = T)
dat <- dat[dat$Site!='Ont',]
nregion <- length(unique(dat$Site))
pops <- unique(dat$Site)
datpop <- list()
nj <- matrix(0,nregion,1)
for (j in 1:nregion) {
  datpop[[j]] <- dat[dat$Site==pops[j],]
  nj[j] <- nrow(datpop[[j]])
}

nrep <- 1000
nloci <- 10
regions <- unique(dat$Site)

nallelemn <- matrix(0,nregion,3)
nallelesd <- matrix(0,nregion,3)
probminshannon <- matrix(0,3,2)
probminehet <- matrix(0,3,2)
probminallele <- matrix(0,3,2)
rownames(probminshannon) <- c('15','16','17')
colnames(probminshannon) <- c('jacknife','bootstrap')
rownames(probminehet) <- c('15','16','17')
colnames(probminehet) <- c('jacknife','bootstrap')
rownames(probminallele) <- c('15','16','17')
colnames(probminallele) <- c('jacknife','bootstrap')
rownames(nallelemn) <- pops
colnames(nallelemn) <- c('15','16','17')
rownames(nallelesd) <- pops
colnames(nallelesd) <- c('15','16','17')

popsizes <- c(15,16,17)
boots <- c(F,T)
sites <- as.character(unique(dat$Site))

for (i1 in 1:3) {
  for (i2 in 1:2) {
    stats <- popsim(datpop,nj,popsizes[i1],boots[i2],nrep,nregion,nloci)
    shannon <- stats[[1]]
    ehet <- stats[[2]]
    nallele <- stats[[3]]
    ehet2 <- 1-1/nloci*ehet
    shmin <- apply(shannon,1,min)
    ismin <- (shannon[,1]==shmin)
    probminshannon[i1,i2] <- mean(ismin)
    hetmin <- apply(ehet2,1,min)
    ismin <- (ehet2[,1]==hetmin)
    probminehet[i1,i2] <- mean(ismin)
    nallelemin <- apply(nallele,1,min)
    ismin <- (nallele[,1]==nallelemin)
    probminallele[i1,i2] <- mean(ismin)
    
    if (boots[i2]) {
      nallelemn[,i1] <- apply(nallele/nloci,2,mean)
      nallelesd[,i1] <- apply(nallele/nloci,2,sd)
    }
    
    shndat <- data.frame(shannon)
    names(shndat) <- sites
    shndat$id <- seq(nrep)
    shndat <- melt(shndat,id='id')
    ggplot(shndat,aes(y=value,x=variable)) + stat_summary(fun.data = f, geom="boxplot")
  }
}

probminshannon
probminehet
probminallele
nallelemn
nallelesd