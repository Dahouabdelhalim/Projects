##BiSSE and HiSSE
library(diversitree)
library(hisse)
library(ape)
library(geiger)

##hisse analysis for the categorial characters
##HiSSE: The hidden state speciation and extinction model
##This model, proposed by Beaulieu and O'Meara, offers one possible solution to the problem 
##identified by Rabosky and Goldberg with the SSE methods. 
##The idea is to provide an additional model to test a state dependent model against. 
##In this model, the change in diversification is attributed to a second, unobserved character 
##that co-occurs with the observed character identified in a typical BiSSE analysis. 
##In this example, we are going to ask if transitions to coral reefs in haemulids change 
##the diversification rate.
##First lets get the data read.
library(hisse)
## Loading required package: GenSA
library(geiger)
##get the tree and data
tree <- read.tree("C_tree.nwk")
data <- read.csv("stigma.csv", row.names = 1)
head(data)
name.check(tree, data)

hd<-cbind(rownames(data), data[,"stigma"])


#1 BiSSE-like unequal turnover rate and extinction rate (BiSSE-like full model, BLFull)

BLFull<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=TransMatMaker.old())

BLFull.logL <- BLFull$loglik
BLFull.logL
BLFull.AIC <- BLFull$AIC
BLFull.AIC
BLFull.AICc <- BLFull$AICc
BLFull.AICc


#2 BiSSE-like equal turnover rate but unequal extinction rate (BiSSE-like equal.t model, BLEqual.t)

BLEqual.t<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,2,0,0), trans.rate=TransMatMaker.old())

BLEqual.t.logL <- BLEqual.t$loglik
BLEqual.t.logL
BLEqual.t.AIC <- BLEqual.t$AIC
BLEqual.t.AIC
BLEqual.t.AICc <- BLEqual.t$AICc
BLEqual.t.AICc


#3 BiSSE-like unequal turnover rate but equal extinction rate (BiSSE-like equal.e model, BLEqual.e)

BLEqual.e<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=TransMatMaker.old())

BLEqual.e.logL <- BLEqual.e$loglik
BLEqual.e.logL
BLEqual.e.AIC <- BLEqual.e$AIC
BLEqual.e.AIC
BLEqual.e.AICc <- BLEqual.e$AICc
BLEqual.e.AICc


#4 BiSSE-like equal transition (BiSSE-like equal.q model, BLEqual.q)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.allequal<-ParEqual(trans.rates, c(1,2))
trans.rates.allequal

BLEqual.q<-hisse.old(tree, hd, hidden.states=FALSE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.allequal)

BLEqual.q.logL <- BLEqual.q$loglik
BLEqual.q.logL
BLEqual.q.AIC <- BLEqual.q$AIC
BLEqual.q.AIC
BLEqual.q.AICc <- BLEqual.q$AICc
BLEqual.q.AICc


#5 BiSSE-like non transition from stateA to stateB (BiSSE-like qPD0 model, BL.qPD0)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BL.qPD0 <- hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.irrev)

BLEqual.qPD0.logL <- BLEqual.qPD0$loglik
BLEqual.qPD0.logL
BLEqual.qPD0.AIC <- BLEqual.qPD0$AIC
BLEqual.qPD0.AIC
BLEqual.qPD0.AICc <- BLEqual.qPD0$AICc
BLEqual.qPD0.AICc


#6 BiSSE-like equal turnover rate and extinction rate (BiSSE-like equal.te model, BLEqual.te)

BLEqual.te<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.allequal)

BLEqual.te.logL <- BLEqual.te$loglik
BLEqual.te.logL
BLEqual.te.AIC <- BLEqual.te$AIC
BLEqual.te.AIC
BLEqual.te.AICc <- BLEqual.te$AICc
BLEqual.te.AICc


#7 BiSSE-like equal turnover rate and transition rate (BiSSE-like equal.tq model, BLEqual.tq)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.allequal<-ParEqual(trans.rates, c(1,2))
trans.rates.allequal

BLEqual.tq<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.allequal)

BLEqual.tq.logL <- BLEqual.tq$loglik
BLEqual.tq.logL
BLEqual.tq.AIC <- BLEqual.tq$AIC
BLEqual.tq.AIC
BLEqual.tq.AICc <- BLEqual.tq$AICc
BLEqual.tq.AICc


#8 BiSSE-like equal extinction rate and transition rate (BiSSE-like equal.eq model, BLEqual.eq)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.allequal<-ParEqual(trans.rates, c(1,2))
trans.rates.allequal

BLEqual.eq<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.allequal)

BLEqual.eq.logL <- BLEqual.eq$loglik
BLEqual.eq.logL
BLEqual.eq.AIC <- BLEqual.eq$AIC
BLEqual.eq.AIC
BLEqual.eq.AICc <- BLEqual.eq$AICc
BLEqual.eq.AICc


#9 BiSSE-like equal turnover rate and non transition from polyploidy to diploidy (BiSSE-like equal.t, qPD0 model, BLEqual.t.qPD0)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BLEqual.t.qPD0<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.irrev)

BLEqual.t.qPD0.logL <- BLEqual.t.qPD0$loglik
BLEqual.t.qPD0.logL
BLEqual.t.qPD0.AIC <- BLEqual.t.qPD0$AIC
BLEqual.t.qPD0.AIC
BLEqual.t.qPD0.AICc <- BLEqual.t.qPD0$AICc
BLEqual.t.qPD0.AICc


#10 BiSSE-like equal extinction rate and non transition from polyploidy to diploidy (BiSSE-like equal.e, qPD0 model, BLEqual.e.qPD0)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BLEqual.e.qPD0<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.irrev)

BLEqual.e.qPD0.logL <- BLEqual.e.qPD0$loglik
BLEqual.e.qPD0.logL
BLEqual.e.qPD0.AIC <- BLEqual.e.qPD0$AIC
BLEqual.e.qPD0.AIC
BLEqual.e.qPD0.AICc <- BLEqual.e.qPD0$AICc
BLEqual.e.qPD0.AICc


#11 BiSSE-like equal turnover rate, extinction rate and transition rate (BiSSE-like equal.teq model, BLEqual.teq)

BLEqual.teq<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=TransMatMaker.old())

BLEqual.teq.logL <- BLEqual.teq$loglik
BLEqual.teq.logL
BLEqual.teq.AIC <- BLEqual.teq$AIC
BLEqual.teq.AIC
BLEqual.teq.AICc <- BLEqual.teq$AICc
BLEqual.teq.AICc


#12 BiSSE-like equal turnover rate and equal extinction rate, but non transition from polyploidy to diploidy (BiSSE-like equal.te, qPD0 model, BLEqual.te.qPD0)

trans.rates<-TransMatMaker.old()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BLEqual.te.qPD0<-hisse.old(tree, hd, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.irrev)

BLEqual.te.qPD0.logL <- BLEqual.te.qPD0$loglik
BLEqual.te.qPD0.logL
BLEqual.te.qPD0.AIC <- BLEqual.te.qPD0$AIC
BLEqual.te.qPD0.AIC
BLEqual.te.qPD0.AICc <- BLEqual.te.qPD0$AICc
BLEqual.te.qPD0.AICc


## HiSSE models (with hidden states)
#13 HiSSE full model (HFull)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates

HFull<-hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates)

HFull.logL <- HFull$loglik
HFull.logL
HFull.AIC <- HFull$AIC
HFull.AIC
HFull.AICc <- HFull$AICc
HFull.AICc


#14 HiSSE model without dual transitions between both the observed trait and the hidden trait (HiSSE nodual, H.nodual)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

H.nodual<-hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual)

H.nodual.logL <- H.nodual$loglik
H.nodual.logL
H.nodual.AIC <- H.nodual$AIC
H.nodual.AIC
H.nodual.AICc <- H.nodual$AICc
H.nodual.AICc


#15 HiSSE model with all equal transitions (HiSSE equal.q, H.Equal.q)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.allequal <- ParEqual(trans.rates, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,11,1,12))
trans.rates.allequal

H.Equal.q <- hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.allequal)

H.Equal.q.logL <- H.Equal.q$loglik
H.Equal.q.logL
H.Equal.q.AIC <- H.Equal.q$AIC
H.Equal.q.AIC
H.Equal.q.AICc <- H.Equal.q$AICc
H.Equal.q.AICc


#16 HiSSE model with equal transitions for nodual ones (HiSSE equal.q/nodual, H.nodual.Equal.q)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual.allequal <- ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

H.nodual.Equal.q <- hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.allequal)

H.nodual.Equal.q.logL <- H.nodual.Equal.q$loglik
H.nodual.Equal.q.logL
H.nodual.Equal.q.AIC <- H.nodual.Equal.q$AIC
H.nodual.Equal.q.AIC
H.nodual.Equal.q.AICc <- H.nodual.Equal.q$AICc
H.nodual.Equal.q.AICc


#17 HiSSE model with non transition from polyploidy to diploidy (HiSSE qPD0, H.qPD0)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.irrev2 <- ParDrop(trans.rates, c(1,3,8,9))
trans.rates.irrev2

H.qPD0 <- hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.irrev2)

H.qPD0.logL <- H.qPD0$loglik
H.qPD0.logL
H.qPD0.AIC <- H.qPD0$AIC
H.qPD0.AIC
H.qPD0.AICc <- H.qPD0$AICc
H.qPD0.AICc


#18 HiSSE model with non dual transition and non transition from polyploidy to diploidy (HiSSE qPD0/nodual, H.nodual.qPD0)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.irrev3 <- ParDrop(trans.rates.nodual, c(1,0,6))
trans.rates.irrev3

H.nodual.qPD0 <- hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.irrev3)

H.nodual.qPD0.logL <- H.nodual.qPD0$loglik
H.nodual.qPD0.logL
H.nodual.qPD0.AIC <- H.nodual.qPD0$AIC
H.nodual.qPD0.AIC
H.nodual.qPD0.AICc <- H.nodual.qPD0$AICc
H.nodual.qPD0.AICc


#19 HiSSE model with all equal transitions and non transition from polyploidy to diploidy (HiSSE euqal.q/qPD0, H.Equal.q.qPD0)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.irrev4 <- ParDrop(trans.rates, c(1,3,8,9))
trans.rates.irrev4
trans.rates.allequal2 <- ParEqual(trans.rates.irrev4, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.allequal2

H.Equal.q.qPD0 <- hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.allequal2)

H.Equal.q.qPD0.logL <- H.Equal.q.qPD0$loglik
H.Equal.q.qPD0.logL
H.Equal.q.qPD0.AIC <- H.Equal.q.qPD0$AIC
H.Equal.q.qPD0.AIC
H.Equal.q.qPD0.AICc <- H.Equal.q.qPD0$AICc
H.Equal.q.qPD0.AICc


#20 HiSSE model with equal transitions for nodual ones and non transition from polyploidy to diploidy (HiSSE equal.q/qPD0/nodual, H.nodual.Equal.q.qPD0)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual2 <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual2
trans.rates.irrev5 <- ParDrop(trans.rates.nodual2, c(1,0,6))
trans.rates.irrev5
trans.rates.nodual.allequal2 <- ParEqual(trans.rates.irrev5, c(1,2,1,3,1,4,1,5))
trans.rates.nodual.allequal2

H.nodual.Equal.q.qPD0 <- hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.allequal2)

H.nodual.Equal.q.qPD0.logL <- H.nodual.Equal.q.qPD0$loglik
H.nodual.Equal.q.qPD0.logL
H.nodual.Equal.q.qPD0.AIC <- H.nodual.Equal.q.qPD0$AIC
H.nodual.Equal.q.qPD0.AIC
H.nodual.Equal.q.qPD0.AICc <- H.nodual.Equal.q.qPD0$AICc
H.nodual.Equal.q.qPD0.AICc


#21 HiSSE model with equal turnover rate (HiSSE equal.t, H.Equal.t)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates

H.Equal.t<-hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,1,1,1), eps.anc=c(1,2,3,4), trans.rate=trans.rates)

H.Equal.t.logL <- H.Equal.t$loglik
H.Equal.t.logL
H.Equal.t.AIC <- H.Equal.t$AIC
H.Equal.t.AIC
H.Equal.t.AICc <- H.Equal.t$AICc
H.Equal.t.AICc


#22 HiSSE model with equal extinction rate (HiSSE equal.e, H.Equal.e)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates

H.Equal.e<-hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates)

H.Equal.e.logL <- H.Equal.e$loglik
H.Equal.e.logL
H.Equal.e.AIC <- H.Equal.e$AIC
H.Equal.e.AIC
H.Equal.e.AICc <- H.Equal.e$AICc
H.Equal.e.AICc


#23 HiSSE model with 3 rates (HiSSE 3 rates, H.3rates) 

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.nodual.threerates <- trans.rates.nodual

to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1

to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2

to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates

H.3rates<-hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.threerates)

H.3rates.logL <- H.3rates$loglik
H.3rates.logL
H.3rates.AIC <- H.3rates$AIC
H.3rates.AIC
H.3rates.AICc <- H.3rates$AICc
H.3rates.AICc


#24 HiSSE model with two rate irreversible (HiSSE 2 rate/irreversible, H.2rates.irrev)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.two.rates<- trans.rates.nodual
to.change <- cbind(c(1,3), c(2,4))
trans.rates.two.rates[to.change] = 1

to.change <- cbind(c(2,4), c(1,3))
trans.rates.two.rates[to.change] = 0

to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.two.rates[to.change] = 2
trans.rates.two.rates

H.2rates.irrev <-hisse.old(tree, hd, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.two.rates)

H.2rates.irrev.logL <- H.2rates.irrev$loglik
H.2rates.irrev.logL
H.2rates.irrev.AIC <- H.2rates.irrev$AIC
H.2rates.irrev.AIC
H.2rates.irrev.AICc <- H.2rates.irrev$AICc
H.2rates.irrev.AICc


# CID2
#25 HiSSE null-two model (HiSSE null-2,HNull2)

trans.rates.nodual[!is.na(trans.rates.nodual) & !trans.rates.nodual == 0] = 1

HNull2 <- hisse.old(tree, hd, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)

HNull2.logL <- HNull2$loglik
HNull2.logL
HNull2.AIC <- HNull2$AIC
HNull2.AIC
HNull2.AICc <- HNull2$AICc
HNull2.AICc


#26 HiSSE null-two model with non transition from polyploidy to diploidy (HiSSE null-2/qPD0, HNull2.qPD0)

HNull2.qPD0 <- hisse.old(tree, hd, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.irrev2)

HNull2.qPD0.logL <- HNull2.qPD0$loglik
HNull2.qPD0.logL
HNull2.qPD0.AIC <- HNull2.qPD0$AIC
HNull2.qPD0.AIC
HNull2.qPD0.AICc <- HNull2.qPD0$AICc
HNull2.qPD0.AICc


#27 HiSSE null-two model with non dual transition (HiSSE null-2/nodual, HNull2.nondual)

HNull2.nondual <- hisse.old(tree, hd, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)

HNull2.nondual.logL <- HNull2.nondual$loglik
HNull2.nondual.logL
HNull2.nondual.AIC <- HNull2.nondual$AIC
HNull2.nondual.AIC
HNull2.nondual.AICc <- HNull2.nondual$AICc
HNull2.nondual.AICc


# CID3
#28 HiSSE null-three model, HNull3

HNull3 <- hisse.null4(phy, states, f=sampfrac, trans.type = "three.rate")

HNull3.logL <- HNull3$loglik
HNull3.logL
HNull3.AIC <- HNull3$AIC
HNull3.AIC
HNull3.AICc <- HNull3$AICc
HNull3.AICc


# CID4 
#29 HiSSE null-four model, HNull4

HNull4 <- hisse.null4.old(phy, states, f=sampfrac, trans.type = "equal")

HNull4.logL <- HNull4$loglik
HNull4.logL
HNull4.AIC <- HNull4$AIC
HNull4.AIC
HNull4.AICc <- HNull4$AICc
HNull4.AICc


#30 HiSSE null-two model with three rates (HiSSE null-2/3 rates, HNull2.3rates)

HNull2.3rates<- hisse.null4.old(tree, hd, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual.threerates)

HNull2.3rates.logL <- HNull2.3rates$loglik
HNull2.3rates.logL
HNull2.3rates.AIC <- HNull2.3rates$AIC
HNull2.3rates.AIC
HNull2.3rates.AICc <- HNull2.3rates$AICc
HNull2.3rates.AICc

