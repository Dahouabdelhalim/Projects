rm(list=ls(all=T))
setwd('~/R')
require('BayesLogit')
require('reshape')
require('ggplot2')

birds <- read.csv('IAL_EPP_allpops.csv')

m1 <- glm(EPY~1,data=birds,family='binomial')
summary(m1)
m2 <- glm(EPY~Pop-1,data=birds,family='binomial')
summary(m2)

mles <- c(m2$coef,m1$coef)
areas <- unique(as.character(birds$Pop))
areas <- c(areas,"Continental mean")

#Random-effects model
nmc <- 11000
X <- model.matrix(EPY~Pop-1,data=birds)
p <- ncol(X)
bet <- matrix(0,p,1)
sigma2b <- 9
sigma2c <- 9
bet0b <- 0
ac <- 2
bc <- 2
y <- birds$EPY
kappa <- y-.5
bet0 <- 0
BETA <- matrix(0,nmc,p)
BET0 <- matrix(0,nmc,1)
sigma2 <- 1
n <- nrow(X)

for (t in 1:nmc) {
  if (t%%100 == 0) {
    print(t)
  }
  omega <- rpg.devroye(n,1,X%*%bet)
  Omega <- diag(omega)
  B <- c(sigma2b,rep(sigma2,p-1))
  V <- solve(t(X)%*%Omega%*%X+diag(1/B))
  b <- c(bet0b,rep(bet0,p-1))
  mv <- V%*%(t(X)%*%kappa+diag(1/B)%*%b)
  bet <- rnorm(p,mv,sqrt(diag(V)))
  
  vbet <- (1/sigma2c+(p-1)*1/sigma2)^(-1)
  mbet <- vbet*(p-1)*1/sigma2*mean(bet[2:p])
  bet0 <- rnorm(1,mbet,sqrt(vbet))
  
  ssq <- sum((bet[2:p]-bet0)^2)
  sigma2 <- rgamma(1,ac+(p-1)/2,rate=bc+ssq/2)
  
  BETA[t,] <- bet
  BET0[t] <- bet0
}

bdat <- data.frame(cbind(BETA[1000:11000,],BET0[1000:11000]))
names(bdat) <- areas
bdat$id <- seq(nmc-999)
bdat <- melt(bdat,id="id")
mles <- data.frame(mles)
names(mles)<-'z'
mles$variable<-areas

#print to screen
p<-ggplot(bdat,aes(x=value)) + geom_histogram() + facet_wrap(~variable) + geom_vline(aes(xintercept=z),mles)
p + labs(x="Value",y="Count") + theme(strip.text=element_text(size=14),axis.text = element_text(size = 14),axis.title = element_text(size = rel(1.5)))
