
## read in master data frame
load('~/master.time.frame.final.rdata', verbose=T)
master<-master.time.frame# for easier use
head(master)


## ********************************************************************************
## Origination and extinction rates by period boundary-crosser method: Fig. 1
## ********************************************************************************

mtf <- master

#call plyr library for later use
library(plyr)

#identify extinction and origination by stage and then code as binary variables
orex <- ddply(mtf, c("taxon_genus"), function(df)c(max(df$stage), min(df$stage)))
names(orex) <- c("taxon_genus", "or.stage", "ex.stage")

mtf.1 <- merge(mtf, orex, by="taxon_genus")
mtf.1$or <- 0
mtf.1$or[mtf.1$stage==mtf.1$or.stage] <- 1
mtf.1$ex <- 0
mtf.1$ex[mtf.1$stage==mtf.1$ex.stage] <- 1

mtf.1$duration <- mtf.1$or.stage - mtf.1$ex.stage

mtf.1 <- subset(mtf.1, duration<150 & class!="" & (ex==0 | or==0)) #remove, long-ranging genera and singletons

#count number of modes by class for total dataset, Holocene, Permian, and Ordovician time slices #AMB: plus ends of other periods
mt <- ddply(mtf.1, c("class"), function(df)c(length(unique(na.omit(df$ecospace))), 
                                             length(unique(df$taxon_genus)),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==453]))), # Ordovician (before extinction)
                                             length(unique(df$taxon_genus[df$stage==453])),
                                             length(unique(df$taxon_genus[df$stage>=453])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==419.2]))), # Silurian, used to be 423
                                             length(unique(df$taxon_genus[df$stage==419.2])),
                                             length(unique(df$taxon_genus[df$stage>=419.2])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==358.9]))), # Devonian, used to be 372.2
                                             length(unique(df$taxon_genus[df$stage==358.9])),
                                             length(unique(df$taxon_genus[df$stage>=358.9])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==298.9]))), # Carboniferous, used to be 303.7
                                             length(unique(df$taxon_genus[df$stage==298.9])),
                                             length(unique(df$taxon_genus[df$stage>=298.9])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==252.2]))), # Permian, used to be 254.2
                                             length(unique(df$taxon_genus[df$stage==252.2])),
                                             length(unique(df$taxon_genus[df$stage>=252.2])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==201.3]))), # Triassic, used to be 208.5 
                                             length(unique(df$taxon_genus[df$stage==201.3])),
                                             length(unique(df$taxon_genus[df$stage>=201.3])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==145]))), # Jurassic, used to be 152.1
                                             length(unique(df$taxon_genus[df$stage==145])),
                                             length(unique(df$taxon_genus[df$stage>=145])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==66]))), # Cretaceous, used to be 72.1
                                             length(unique(df$taxon_genus[df$stage==66])),
                                             length(unique(df$taxon_genus[df$stage>=66])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==61.6]))), # Danian, used to be Paleocene (59.2)
                                             length(unique(df$taxon_genus[df$stage==61.6])),
                                             length(unique(df$taxon_genus[df$stage>=61.6])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==23.03]))), # Paleogene, used to be 28.1
                                             length(unique(df$taxon_genus[df$stage==23.03])),
                                             length(unique(df$taxon_genus[df$stage>=23.03])),
                                                                                        
                                             length(unique(na.omit(df$ecospace[df$stage==0.0117]))), #Holocene, used to be 0.0117
                                             length(unique(df$taxon_genus[df$stage==0.0117])),
                                             max(df$stage)))

names(mt) <- c("class", "modes.tot", "gen.tot", "modes.ord", "gen.ord", "gen.preord", "modes.sil", "gen.sil", "gen.presil", 
               "modes.dev", "gen.dev", "gen.predev", "modes.carb", "gen.carb", "gen.precarb", "modes.perm", "gen.perm", 
               "gen.preperm", "modes.tri", "gen.tri", "gen.pretri", "modes.jur", "gen.jur", "gen.prejur", 
               "modes.cret", "gen.cret", "gen.precret", "modes.dani", "gen.dani", "gen.predani", "modes.palg", "gen.palg", "gen.prepalg", "modes.mod", "gen.mod", "orig.class")

mtf.2 <- merge(mtf.1, mt, by="class")

min.gen.standing <- 2
min.gen.total <- 100


## MODERN ANALYSIS FOOTE RATES ##
mod <- ddply(mtf.2, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])), #extinction rate
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])), #origination rate
    length(df$ex[df$ex==0 & df$or==0]), #range through genera
    length(df$or[df$ex==0]), #surviving genera
    length(df$ex[df$or==0]))) #genera entering from previous stage
names(mod) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

mod.1 <- subset(mod, logex<100 & logor<100) #remove infinite rates and NA values
mod.2 <- ddply(mod.1, c("class"), function(df) #calculate average extinction and origination rates
  c(mean(df$logex), mean(df$logor)))
names(mod.2) <- c("class", "mean.ex", "mean.or")
mod.2$mean.or <- apply(mod.2[,2:3], MARGIN=1, FUN=max) 
mod.3 <- merge(mod.2, mt, by="class") #merge with other data by class (genus and mode counts)

mod.results <- subset(mod.3, gen.tot>=min.gen.total & gen.mod>=min.gen.standing & orig.class>400) #subset based on diversity criteria (standing and total diversity) 

## PALEOGENE ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=23.03)
palg <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(palg) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

palg.1 <- subset(palg, logex<100 & logor<100)
palg.2 <- ddply(palg.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(palg.2) <- c("class", "mean.ex", "mean.or")
palg.2$mean.or <- apply(palg.2[,2:3], MARGIN=1, FUN=max) 
palg.3 <- merge(palg.2, mt, by="class")

palg.results <- subset(palg.3, gen.tot>=min.gen.total & gen.palg>=min.gen.standing & orig.class>400)

## DANIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=61.6) 
dani <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(dani) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

dani.1 <- subset(dani, logex<100 & logor<100)
dani.2 <- ddply(dani.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(dani.2) <- c("class", "mean.ex", "mean.or")
dani.2$mean.or <- apply(dani.2[,2:3], MARGIN=1, FUN=max) 
dani.3 <- merge(dani.2, mt, by="class")

dani.results <- subset(dani.3, gen.tot>=min.gen.total & gen.dani>=min.gen.standing & orig.class>400)

## CRETACEOUS ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=66)
cret <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(cret) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

cret.1 <- subset(cret, logex<100 & logor<100)
cret.2 <- ddply(cret.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(cret.2) <- c("class", "mean.ex", "mean.or")
cret.2$mean.or <- apply(cret.2[,2:3], MARGIN=1, FUN=max) 
cret.3 <- merge(cret.2, mt, by="class")

cret.results <- subset(cret.3, gen.tot>=min.gen.total & gen.cret>=min.gen.standing & orig.class>400)

## JURASSIC ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=145)
jur <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(jur) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

jur.1 <- subset(jur, logex<100 & logor<100)
jur.2 <- ddply(jur.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(jur.2) <- c("class", "mean.ex", "mean.or")
jur.2$mean.or <- apply(jur.2[,2:3], MARGIN=1, FUN=max) 
jur.3 <- merge(jur.2, mt, by="class")

jur.results <- subset(jur.3, gen.tot>=min.gen.total & gen.jur>=min.gen.standing & orig.class>400)

## TRIASSIC ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=201.3)
tri <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(tri) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

tri.1 <- subset(tri, logex<100 & logor<100)
tri.2 <- ddply(tri.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(tri.2) <- c("class", "mean.ex", "mean.or")
tri.2$mean.or <- apply(tri.2[,2:3], MARGIN=1, FUN=max) 
tri.3 <- merge(tri.2, mt, by="class")

tri.results <- subset(tri.3, gen.tot>=min.gen.total & gen.tri>=min.gen.standing & orig.class>400)

## PERMIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=252.2) #changed from 252
perm <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(perm) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

perm.1 <- subset(perm, logex<100 & logor<100)
perm.2 <- ddply(perm.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(perm.2) <- c("class", "mean.ex", "mean.or")
perm.2$mean.or <- apply(perm.2[,2:3], MARGIN=1, FUN=max) 
perm.3 <- merge(perm.2, mt, by="class")

perm.results <- subset(perm.3, gen.tot>=min.gen.total & gen.perm>=min.gen.standing & orig.class>400)

## CARBONIFEROUS ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=298.9)
carb <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(carb) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

carb.1 <- subset(carb, logex<100 & logor<100)
carb.2 <- ddply(carb.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(carb.2) <- c("class", "mean.ex", "mean.or")
carb.2$mean.or <- apply(carb.2[,2:3], MARGIN=1, FUN=max) 
carb.3 <- merge(carb.2, mt, by="class")

carb.results <- subset(carb.3, gen.tot>=min.gen.total & gen.carb>=min.gen.standing & orig.class>400)

## DEVONIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=358.9)
dev <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(dev) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

dev.1 <- subset(dev, logex<100 & logor<100)
dev.2 <- ddply(dev.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(dev.2) <- c("class", "mean.ex", "mean.or")
dev.2$mean.or <- apply(dev.2[,2:3], MARGIN=1, FUN=max) 
dev.3 <- merge(dev.2, mt, by="class")

dev.results <- subset(dev.3, gen.tot>=min.gen.total & gen.dev>=min.gen.standing & orig.class>400)

## SILURIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=419.2)
sil <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(sil) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

sil.1 <- subset(sil, logex<100 & logor<100)
sil.2 <- ddply(sil.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(sil.2) <- c("class", "mean.ex", "mean.or")
sil.2$mean.or <- apply(sil.2[,2:3], MARGIN=1, FUN=max) 
sil.3 <- merge(sil.2, mt, by="class")

sil.results <- subset(sil.3, gen.tot>=min.gen.total & gen.sil>=min.gen.standing & orig.class>400)

## ORDOVICIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=453)
ord <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(ord) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

ord.1 <- subset(ord, logex<100 & logor<100)
ord.2 <- ddply(ord.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(ord.2) <- c("class", "mean.ex", "mean.or")
ord.2$mean.or <- apply(ord.2[,2:3], MARGIN=1, FUN=max) 
ord.3 <- merge(ord.2, mt, by="class")

ord.results <- subset(ord.3, gen.tot>=min.gen.total & gen.ord>=min.gen.standing & orig.class>400)




## FIGURE FOR GENUS DIVERSITY vs Rates##
pdf("Figure 20 panel dates.pdf", h=8, w=8)
par(mfrow=c(4,5), mar=c(4,4,1,1), mgp=c(1.75,.6,0), pin=c(.8,.8), tcl=-.3, lab=c(4,4,1), cex.axis=.8, cex.main=.8)

y.min <- 0
#y.max <- 1.2
y.max <- .6
x.min <- 0
x.max.1 <- log10(max(mod.results$gen.mod))+0.2
x.max.2 <- log10(max(mod.results$modes.mod))+0.2

if (F) { # turns off ordovician
plot(ord.results$mean.ex~log10(ord.results$gen.ord), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Ordovician",
     xlim=c(x.min,x.max.1), ylim=c(y.min, 1), type="n")
#points(ord.results$mean.or~log10(ord.results$gen.ord), pch=2, col="blue", cex=2)
clip(log10(min(ord.results$gen.ord)), log10(max(ord.results$gen.ord)), 0, .6)
abline(lm(ord.results$mean.ex~log10(ord.results$gen.ord)), col="red", lwd=2)
abline(lm(ord.results$mean.or~log10(ord.results$gen.ord)), col="blue", lwd=2)
}

plot(sil.results$mean.ex~log10(sil.results$gen.sil), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Silurian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(sil.results$mean.or~log10(sil.results$gen.sil), pch=2, col="blue", cex=2)
clip(log10(min(sil.results$gen.sil)), log10(max(sil.results$gen.sil)), 0, .6)
abline(lm(sil.results$mean.ex~log10(sil.results$gen.sil)), col="red", lwd=2)
abline(lm(sil.results$mean.or~log10(sil.results$gen.sil)), col="blue", lwd=2)

plot(dev.results$mean.ex~log10(dev.results$gen.dev), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Devonian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(dev.results$mean.or~log10(dev.results$gen.dev), pch=2, col="blue", cex=2)
clip(log10(min(dev.results$gen.dev)), log10(max(dev.results$gen.dev)), 0, .6)
abline(lm(dev.results$mean.ex~log10(dev.results$gen.dev)), col="red", lwd=2)
abline(lm(dev.results$mean.or~log10(dev.results$gen.dev)), col="blue", lwd=2)

plot(carb.results$mean.ex~log10(carb.results$gen.carb), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Carboniferous",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(carb.results$mean.or~log10(carb.results$gen.carb), pch=2, col="blue", cex=2)
clip(log10(min(carb.results$gen.carb)), log10(max(carb.results$gen.carb)), 0, .6)
abline(lm(carb.results$mean.ex~log10(carb.results$gen.carb)), col="red", lwd=2)
abline(lm(carb.results$mean.or~log10(carb.results$gen.carb)), col="blue", lwd=2)

plot(perm.results$mean.ex~log10(perm.results$gen.perm), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Permian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(perm.results$mean.or~log10(perm.results$gen.perm), pch=2, col="blue", cex=2)
clip(log10(min(perm.results$gen.perm)), log10(max(perm.results$gen.perm)), 0, .6)
abline(lm(perm.results$mean.ex~log10(perm.results$gen.perm)), col="red", lwd=2)
abline(lm(perm.results$mean.or~log10(perm.results$gen.perm)), col="blue", lwd=2)

plot(tri.results$mean.ex~log10(tri.results$gen.tri), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Triassic",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(tri.results$mean.or~log10(tri.results$gen.tri), pch=2, col="blue", cex=2)
clip(log10(min(tri.results$gen.tri)), log10(max(tri.results$gen.tri)), 0, .6)
abline(lm(tri.results$mean.ex~log10(tri.results$gen.tri)), col="red", lwd=2)
abline(lm(tri.results$mean.or~log10(tri.results$gen.tri)), col="blue", lwd=2)

plot(jur.results$mean.ex~log10(jur.results$gen.jur), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Jurassic",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(jur.results$mean.or~log10(jur.results$gen.jur), pch=2, col="blue", cex=2)
clip(log10(min(jur.results$gen.jur)), log10(max(jur.results$gen.jur)), 0, .6)
abline(lm(jur.results$mean.ex~log10(jur.results$gen.jur)), col="red", lwd=2)
abline(lm(jur.results$mean.or~log10(jur.results$gen.jur)), col="blue", lwd=2)

plot(cret.results$mean.ex~log10(cret.results$gen.cret), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Cretaceous",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(cret.results$mean.or~log10(cret.results$gen.cret), pch=2, col="blue", cex=2)
clip(log10(min(cret.results$gen.cret)), log10(max(cret.results$gen.cret)), 0, .6)
abline(lm(cret.results$mean.ex~log10(cret.results$gen.cret)), col="red", lwd=2)
abline(lm(cret.results$mean.or~log10(cret.results$gen.cret)), col="blue", lwd=2)

plot(dani.results$mean.ex~log10(dani.results$gen.dani), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Danian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(dani.results$mean.or~log10(dani.results$gen.dani), pch=2, col="blue", cex=2)
clip(log10(min(dani.results$gen.dani)), log10(max(dani.results$gen.dani)), 0, .6)
abline(lm(dani.results$mean.ex~log10(dani.results$gen.dani)), col="red", lwd=2)
abline(lm(dani.results$mean.or~log10(dani.results$gen.dani)), col="blue", lwd=2)

plot(palg.results$mean.ex~log10(palg.results$gen.palg), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Paleogene",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(palg.results$mean.or~log10(palg.results$gen.palg), pch=2, col="blue", cex=2)
clip(log10(min(palg.results$gen.palg)), log10(max(palg.results$gen.palg)), 0, .6)
abline(lm(palg.results$mean.ex~log10(palg.results$gen.palg)), col="red", lwd=2)
abline(lm(palg.results$mean.or~log10(palg.results$gen.palg)), col="blue", lwd=2)

plot(mod.results$mean.ex~log10(mod.results$gen.mod), pch=1, col="red", las=1, cex=2,
     xlab="log10 genera", ylab="mean rate (BC)", main="Modern (Holocene)",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="n")
#points(mod.results$mean.or~log10(mod.results$gen.mod), pch=2, col="blue", cex=2)
clip(log10(min(mod.results$gen.mod)), log10(max(mod.results$gen.mod)), 0, .6)
abline(lm(mod.results$mean.ex~log10(mod.results$gen.mod)), col="red", lwd=2)
abline(lm(mod.results$mean.or~log10(mod.results$gen.mod)), col="blue", lwd=2)


#dev.off() # disabled so genera and modes plot on same figure

## FIGURE FOR Modes vs Rates ##
#pdf("Figure Modes 10 panel.pdf", h=4, w=8) # disabled so genera and modes plot on same figure
#par(mfrow=c(2,5), mar=c(4,4,1,1), mgp=c(2.5,1,0), pin=c(1,1), tcl=-.5, lab=c(4,4,1))


y.min <- 0
y.max <- .6
x.min <- 0
x.max.1 <- log10(max(mod.results$gen.mod))+0.2
x.max.2 <- log10(max(mod.results$modes.mod))+0.2


plot(sil.results$mean.ex~log10(sil.results$modes.sil), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Silurian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(sil.results$mean.or~log10(sil.results$modes.sil), pch=2, col="blue", cex=2)
clip(log10(min(sil.results$modes.sil)), log10(max(sil.results$modes.sil)), 0, .6)
abline(lm(sil.results$mean.ex~log10(sil.results$modes.sil)), col="red", lwd=2)
abline(lm(sil.results$mean.or~log10(sil.results$modes.sil)), col="blue", lwd=2)

plot(dev.results$mean.ex~log10(dev.results$modes.dev), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Devonian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(dev.results$mean.or~log10(dev.results$modes.dev), pch=2, col="blue", cex=2)
clip(log10(min(dev.results$modes.dev)), log10(max(dev.results$modes.dev)), 0, .6)
abline(lm(dev.results$mean.ex~log10(dev.results$modes.dev)), col="red", lwd=2)
abline(lm(dev.results$mean.or~log10(dev.results$modes.dev)), col="blue", lwd=2)

plot(carb.results$mean.ex~log10(carb.results$modes.carb), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Carboniferous",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(carb.results$mean.or~log10(carb.results$modes.carb), pch=2, col="blue", cex=2)
clip(log10(min(carb.results$modes.carb)), log10(max(carb.results$modes.carb)), 0, .6)
abline(lm(carb.results$mean.ex~log10(carb.results$modes.carb)), col="red", lwd=2)
abline(lm(carb.results$mean.or~log10(carb.results$modes.carb)), col="blue", lwd=2)

plot(perm.results$mean.ex~log10(perm.results$modes.perm), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Permian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(perm.results$mean.or~log10(perm.results$modes.perm), pch=2, col="blue", cex=2)
clip(log10(min(perm.results$modes.perm)), log10(max(perm.results$modes.perm)), 0, .6)
abline(lm(perm.results$mean.ex~log10(perm.results$modes.perm)), col="red", lwd=2)
abline(lm(perm.results$mean.or~log10(perm.results$modes.perm)), col="blue", lwd=2)

plot(tri.results$mean.ex~log10(tri.results$modes.tri), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Triassic",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(tri.results$mean.or~log10(tri.results$modes.tri), pch=2, col="blue", cex=2)
clip(log10(min(tri.results$modes.tri)), log10(max(tri.results$modes.tri)), 0, .6)
abline(lm(tri.results$mean.ex~log10(tri.results$modes.tri)), col="red", lwd=2)
abline(lm(tri.results$mean.or~log10(tri.results$modes.tri)), col="blue", lwd=2)

plot(jur.results$mean.ex~log10(jur.results$modes.jur), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Jurassic",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(jur.results$mean.or~log10(jur.results$modes.jur), pch=2, col="blue", cex=2)
clip(log10(min(jur.results$modes.jur)), log10(max(jur.results$modes.jur)), 0, .6)
abline(lm(jur.results$mean.ex~log10(jur.results$modes.jur)), col="red", lwd=2)
abline(lm(jur.results$mean.or~log10(jur.results$modes.jur)), col="blue", lwd=2)

plot(cret.results$mean.ex~log10(cret.results$modes.cret), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Cretaceous",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(cret.results$mean.or~log10(cret.results$modes.cret), pch=2, col="blue", cex=2)
clip(log10(min(cret.results$modes.cret)), log10(max(cret.results$modes.cret)), 0, .6)
abline(lm(cret.results$mean.ex~log10(cret.results$modes.cret)), col="red", lwd=2)
abline(lm(cret.results$mean.or~log10(cret.results$modes.cret)), col="blue", lwd=2)

plot(dani.results$mean.ex~log10(dani.results$modes.dani), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Danian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(dani.results$mean.or~log10(dani.results$modes.dani), pch=2, col="blue", cex=2)
clip(log10(min(dani.results$modes.dani)), log10(max(dani.results$modes.dani)), 0, .6)
abline(lm(dani.results$mean.ex~log10(dani.results$modes.dani)), col="red", lwd=2)
abline(lm(dani.results$mean.or~log10(dani.results$modes.dani)), col="blue", lwd=2)

plot(palg.results$mean.ex~log10(palg.results$modes.palg), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Paleogene",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(palg.results$mean.or~log10(palg.results$modes.palg), pch=2, col="blue", cex=2)
clip(log10(min(palg.results$modes.palg)), log10(max(palg.results$modes.palg)), 0, .6)
abline(lm(palg.results$mean.ex~log10(palg.results$modes.palg)), col="red", lwd=2)
abline(lm(palg.results$mean.or~log10(palg.results$modes.palg)), col="blue", lwd=2)

plot(mod.results$mean.ex~log10(mod.results$modes.mod), pch=1, col="red", las=1, cex=2,
     xlab="log10 modes", ylab="mean rate (BC)", main="Modern (Holocene)",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="n")
#points(mod.results$mean.or~log10(mod.results$modes.mod), pch=2, col="blue", cex=2)
clip(log10(min(mod.results$modes.mod)), log10(max(mod.results$modes.mod)), 0, .6)
abline(lm(mod.results$mean.ex~log10(mod.results$modes.mod)), col="red", lwd=2)
abline(lm(mod.results$mean.or~log10(mod.results$modes.mod)), col="blue", lwd=2)



dev.off()

##*******************************************************************************
##Fig. S13, S14, Table S5
##*******************************************************************************


#identify extinction and origination by stage and then code as binary variables
orex <- ddply(mtf, c("taxon_genus"), function(df)c(max(df$stage), min(df$stage)))
names(orex) <- c("taxon_genus", "or.stage", "ex.stage")

mtf.1 <- merge(mtf, orex, by="taxon_genus")
mtf.1$or <- 0
mtf.1$or[mtf.1$stage==mtf.1$or.stage] <- 1
mtf.1$ex <- 0
mtf.1$ex[mtf.1$stage==mtf.1$ex.stage] <- 1

mtf.1$duration <- mtf.1$or.stage - mtf.1$ex.stage

mtf.1 <- subset(mtf.1, duration<150 & class!="" & (ex==0 | or==0)) #remove, long-ranging genera and singletons

#count number of modes by class for selected time slices
mt <- ddply(mtf.1, c("class"), function(df)c(length(unique(na.omit(df$ecospace))), 
                                             length(unique(df$taxon_genus)),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==453]))), # Ordovician
                                             length(unique(df$taxon_genus[df$stage==453])),
                                             length(unique(df$taxon_genus[df$stage>=453])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==419.2]))), # Silurian
                                             length(unique(df$taxon_genus[df$stage==419.2])),
                                             length(unique(df$taxon_genus[df$stage>=419.2])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==358.9]))), # Devonian
                                             length(unique(df$taxon_genus[df$stage==358.9])),
                                             length(unique(df$taxon_genus[df$stage>=358.9])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==298.9]))), # Carboniferous
                                             length(unique(df$taxon_genus[df$stage==298.9])),
                                             length(unique(df$taxon_genus[df$stage>=298.9])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==252.2]))), # Permian
                                             length(unique(df$taxon_genus[df$stage==252.2])),
                                             length(unique(df$taxon_genus[df$stage>=252.2])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==201.3]))), # Triassic 
                                             length(unique(df$taxon_genus[df$stage==201.3])),
                                             length(unique(df$taxon_genus[df$stage>=201.3])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==145]))), # Jurassic
                                             length(unique(df$taxon_genus[df$stage==145])),
                                             length(unique(df$taxon_genus[df$stage>=145])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==66]))), # Cretaceous
                                             length(unique(df$taxon_genus[df$stage==66])),
                                             length(unique(df$taxon_genus[df$stage>=66])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==61.6]))), # Danian
                                             length(unique(df$taxon_genus[df$stage==61.6])),
                                             length(unique(df$taxon_genus[df$stage>=61.6])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==23.03]))), # Paleogene
                                             length(unique(df$taxon_genus[df$stage==23.03])),
                                             length(unique(df$taxon_genus[df$stage>=23.03])),
                                                                                        
                                             length(unique(na.omit(df$ecospace[df$stage==0.0117]))), #Holocene
                                             length(unique(df$taxon_genus[df$stage==0.0117])),
                                             max(df$stage)))

names(mt) <- c("class", "modes.tot", "gen.tot", "modes.ord", "gen.ord", "gen.preord", "modes.sil", "gen.sil", "gen.presil", 
               "modes.dev", "gen.dev", "gen.pblueev", "modes.carb", "gen.carb", "gen.precarb", "modes.perm", "gen.perm", 
               "gen.preperm", "modes.tri", "gen.tri", "gen.pretri", "modes.jur", "gen.jur", "gen.prejur", 
               "modes.cret", "gen.cret", "gen.precret", "modes.dani", "gen.dani", "gen.pblueani", "modes.palg", "gen.palg", "gen.prepalg", "modes.mod", "gen.mod", "orig.class")

mtf.2 <- merge(mtf.1, mt, by="class")

min.gen.standing <- 2
min.gen.total <- 100


## MODERN ANALYSIS FOOTE RATES ##
mod <- ddply(mtf.2, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])), #extinction rate
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])), #origination rate
    length(df$ex[df$ex==0 & df$or==0]), #range through genera
    length(df$or[df$ex==0]), #surviving genera
    length(df$ex[df$or==0]))) #genera entering from previous stage
names(mod) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

mod.1 <- subset(mod, logex<100 & logor<100) #remove infinite rates and NA values
mod.2 <- ddply(mod.1, c("class"), function(df) #calculate average extinction and origination rates
  c(mean(df$logex), mean(df$logor)))
names(mod.2) <- c("class", "mean.ex", "mean.or")
mod.2$mean.or <- apply(mod.2[,2:3], MARGIN=1, FUN=max) 
mod.3 <- merge(mod.2, mt, by="class") #merge with other data by class (genus and mode counts)

mod.results <- subset(mod.3, gen.tot>=min.gen.total & gen.mod>=min.gen.standing & orig.class>400) #subset based on diversity criteria (standing and total diversity) 

## PALEOGENE ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=23.03)
palg <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(palg) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

palg.1 <- subset(palg, logex<100 & logor<100)
palg.2 <- ddply(palg.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(palg.2) <- c("class", "mean.ex", "mean.or")
palg.2$mean.or <- apply(palg.2[,2:3], MARGIN=1, FUN=max) 
palg.3 <- merge(palg.2, mt, by="class")

palg.results <- subset(palg.3, gen.tot>=min.gen.total & gen.palg>=min.gen.standing & orig.class>400)

## DANIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=61.6) 
dani <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(dani) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

dani.1 <- subset(dani, logex<100 & logor<100)
dani.2 <- ddply(dani.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(dani.2) <- c("class", "mean.ex", "mean.or")
dani.2$mean.or <- apply(dani.2[,2:3], MARGIN=1, FUN=max) 
dani.3 <- merge(dani.2, mt, by="class")

dani.results <- subset(dani.3, gen.tot>=min.gen.total & gen.dani>=min.gen.standing & orig.class>400)

## CRETACEOUS ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=66)
cret <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(cret) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

cret.1 <- subset(cret, logex<100 & logor<100)
cret.2 <- ddply(cret.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(cret.2) <- c("class", "mean.ex", "mean.or")
cret.2$mean.or <- apply(cret.2[,2:3], MARGIN=1, FUN=max) 
cret.3 <- merge(cret.2, mt, by="class")

cret.results <- subset(cret.3, gen.tot>=min.gen.total & gen.cret>=min.gen.standing & orig.class>400)

## JURASSIC ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=145)
jur <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(jur) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

jur.1 <- subset(jur, logex<100 & logor<100)
jur.2 <- ddply(jur.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(jur.2) <- c("class", "mean.ex", "mean.or")
jur.2$mean.or <- apply(jur.2[,2:3], MARGIN=1, FUN=max) 
jur.3 <- merge(jur.2, mt, by="class")

jur.results <- subset(jur.3, gen.tot>=min.gen.total & gen.jur>=min.gen.standing & orig.class>400)

## TRIASSIC ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=201.3)
tri <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(tri) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

tri.1 <- subset(tri, logex<100 & logor<100)
tri.2 <- ddply(tri.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(tri.2) <- c("class", "mean.ex", "mean.or")
tri.2$mean.or <- apply(tri.2[,2:3], MARGIN=1, FUN=max) 
tri.3 <- merge(tri.2, mt, by="class")

tri.results <- subset(tri.3, gen.tot>=min.gen.total & gen.tri>=min.gen.standing & orig.class>400)

## PERMIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=252.2) #changed from 252
perm <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(perm) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

perm.1 <- subset(perm, logex<100 & logor<100)
perm.2 <- ddply(perm.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(perm.2) <- c("class", "mean.ex", "mean.or")
perm.2$mean.or <- apply(perm.2[,2:3], MARGIN=1, FUN=max) 
perm.3 <- merge(perm.2, mt, by="class")

perm.results <- subset(perm.3, gen.tot>=min.gen.total & gen.perm>=min.gen.standing & orig.class>400)

## CARBONIFEROUS ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=298.9)
carb <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(carb) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

carb.1 <- subset(carb, logex<100 & logor<100)
carb.2 <- ddply(carb.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(carb.2) <- c("class", "mean.ex", "mean.or")
carb.2$mean.or <- apply(carb.2[,2:3], MARGIN=1, FUN=max) 
carb.3 <- merge(carb.2, mt, by="class")

carb.results <- subset(carb.3, gen.tot>=min.gen.total & gen.carb>=min.gen.standing & orig.class>400)

## DEVONIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=358.9)
dev <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(dev) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

dev.1 <- subset(dev, logex<100 & logor<100)
dev.2 <- ddply(dev.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(dev.2) <- c("class", "mean.ex", "mean.or")
dev.2$mean.or <- apply(dev.2[,2:3], MARGIN=1, FUN=max)
dev.3 <- merge(dev.2, mt, by="class")

dev.results <- subset(dev.3, gen.tot>=min.gen.total & gen.dev>=min.gen.standing & orig.class>400)

## SILURIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=419.2)
sil <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(sil) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

sil.1 <- subset(sil, logex<100 & logor<100)
sil.2 <- ddply(sil.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(sil.2) <- c("class", "mean.ex", "mean.or")
sil.2$mean.or <- apply(sil.2[,2:3], MARGIN=1, FUN=max) 
sil.3 <- merge(sil.2, mt, by="class")

sil.results <- subset(sil.3, gen.tot>=min.gen.total & gen.sil>=min.gen.standing & orig.class>400)

## ORDOVICIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=453)
ord <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(ord) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

ord.1 <- subset(ord, logex<100 & logor<100)
ord.2 <- ddply(ord.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(ord.2) <- c("class", "mean.ex", "mean.or")
ord.2$mean.or <- apply(ord.2[,2:3], MARGIN=1, FUN=max) 
ord.3 <- merge(ord.2, mt, by="class")

ord.results <- subset(ord.3, gen.tot>=min.gen.total & gen.ord>=min.gen.standing & orig.class>400)









## FIGURE FOR GENUS DIVERSITY vs Rates##
pdf("regr orig and diversification.pdf", h=8, w=8)
par(mfrow=c(4,5), mar=c(4,4,1,1), mgp=c(1.75,.6,0), pin=c(1,1), tcl=-.3, lab=c(4,4,1), cex.axis=.8, cex.main=.8)

y.min <- 0
#y.max <- 1.2
y.max <- .8
x.min <- 0
x.max.1 <- log10(max(mod.results$gen.mod))+0.2
x.max.2 <- log10(max(mod.results$modes.mod))+0.2

pval_orig_rate <- rep(NA,10)
pval_div_rate <- rep(NA,10)
counter <- 1

plot((sil.results$mean.or-sil.results$mean.ex)~log10(sil.results$gen.sil), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Silurian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points((sil.results$mean.or)~log10(sil.results$gen.sil), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(sil.results$gen.sil)), log10(max(sil.results$gen.sil)), 0, .6)
abline(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$gen.sil)), col="black", lwd=2)
abline(lm((sil.results$mean.or)~log10(sil.results$gen.sil)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$gen.sil)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(sil.results$mean.or~log10(sil.results$gen.sil)))$coefficients[2,4]
counter <- counter+1

plot(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$gen.dev), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Devonian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(dev.results$mean.or~log10(dev.results$gen.dev), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(dev.results$gen.dev)), log10(max(dev.results$gen.dev)), 0, .6)
abline(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$gen.dev)), col="black", lwd=2)
abline(lm(dev.results$mean.or~log10(dev.results$gen.dev)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$gen.dev)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(dev.results$mean.or~log10(dev.results$gen.dev)))$coefficients[2,4]
counter <- counter+1

plot(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$gen.carb), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Carboniferous",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(carb.results$mean.or~log10(carb.results$gen.carb), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(carb.results$gen.carb)), log10(max(carb.results$gen.carb)), 0, .6)
abline(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$gen.carb)), col="black", lwd=2)
abline(lm(carb.results$mean.or~log10(carb.results$gen.carb)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$gen.carb)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(carb.results$mean.or~log10(carb.results$gen.carb)))$coefficients[2,4]
counter <- counter+1

plot(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$gen.perm), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Permian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(perm.results$mean.or~log10(perm.results$gen.perm), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(perm.results$gen.perm)), log10(max(perm.results$gen.perm)), 0, .6)
abline(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$gen.perm)), col="black", lwd=2)
abline(lm(perm.results$mean.or~log10(perm.results$gen.perm)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$gen.perm)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(perm.results$mean.or~log10(perm.results$gen.perm)))$coefficients[2,4]
counter <- counter+1


plot(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$gen.tri), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Triassic",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(tri.results$mean.or~log10(tri.results$gen.tri), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(tri.results$gen.tri)), log10(max(tri.results$gen.tri)), 0, .6)
abline(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$gen.tri)), col="black", lwd=2)
abline(lm(tri.results$mean.or~log10(tri.results$gen.tri)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$gen.tri)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(tri.results$mean.or~log10(tri.results$gen.tri)))$coefficients[2,4]
counter <- counter+1

plot(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$gen.jur), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Jurassic",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(jur.results$mean.or~log10(jur.results$gen.jur), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(jur.results$gen.jur)), log10(max(jur.results$gen.jur)), 0, .6)
abline(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$gen.jur)), col="black", lwd=2)
abline(lm(jur.results$mean.or~log10(jur.results$gen.jur)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$gen.jur)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(jur.results$mean.or~log10(jur.results$gen.jur)))$coefficients[2,4]
counter <- counter+1

plot(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$gen.cret), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Cretaceous",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(cret.results$mean.or~log10(cret.results$gen.cret), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(cret.results$gen.cret)), log10(max(cret.results$gen.cret)), 0, .6)
abline(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$gen.cret)), col="black", lwd=2)
abline(lm(cret.results$mean.or~log10(cret.results$gen.cret)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$gen.cret)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(cret.results$mean.or~log10(cret.results$gen.cret)))$coefficients[2,4]
counter <- counter+1

plot(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$gen.dani), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Danian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(dani.results$mean.or~log10(dani.results$gen.dani), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(dani.results$gen.dani)), log10(max(dani.results$gen.dani)), 0, .6)
abline(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$gen.dani)), col="black", lwd=2)
abline(lm(dani.results$mean.or~log10(dani.results$gen.dani)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$gen.dani)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(dani.results$mean.or~log10(dani.results$gen.dani)))$coefficients[2,4]
counter <- counter+1

plot(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$gen.palg), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Paleogene",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(palg.results$mean.or~log10(palg.results$gen.palg), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(palg.results$gen.palg)), log10(max(palg.results$gen.palg)), 0, .6)
abline(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$gen.palg)), col="black", lwd=2)
abline(lm(palg.results$mean.or~log10(palg.results$gen.palg)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$gen.palg)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(palg.results$mean.or~log10(palg.results$gen.palg)))$coefficients[2,4]
counter <- counter+1

plot((mod.results$mean.or-mod.results$mean.ex)~log10(mod.results$gen.mod), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Pleistocene",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points((mod.results$mean.or)~log10(mod.results$gen.mod), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(mod.results$gen.mod)), log10(max(mod.results$gen.mod)), 0, .6)
abline(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$gen.mod)), col="black", lwd=2)
abline(lm((mod.results$mean.or)~log10(mod.results$gen.mod)), col="blue", lwd=2)
pval_div_rate[counter] <- summary(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$gen.mod)))$coefficients[2,4]
pval_orig_rate[counter] <- summary(lm(mod.results$mean.or~log10(mod.results$gen.mod)))$coefficients[2,4]
counter <- counter+1


## FIGURE FOR Modes vs Rates ##


pval_modes_orig_rate <- rep(NA,10)
pval_modes_div_rate <- rep(NA,10)
counter <- 1

y.min <- 0
x.min <- 0
x.max.1 <- log10(max(mod.results$gen.mod))+0.2
x.max.2 <- log10(max(mod.results$modes.mod))+0.2

plot((sil.results$mean.or-sil.results$mean.ex)~log10(sil.results$modes.sil), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Silurian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points((sil.results$mean.or)~log10(sil.results$modes.sil), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(sil.results$modes.sil)), log10(max(sil.results$modes.sil)), 0, .6)
abline(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$modes.sil)), col="black", lwd=2)
abline(lm((sil.results$mean.or)~log10(sil.results$modes.sil)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$modes.sil)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(sil.results$mean.or~log10(sil.results$modes.sil)))$coefficients[2,4]
counter <- counter+1

plot(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$modes.dev), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Devonian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(dev.results$mean.or~log10(dev.results$modes.dev), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(dev.results$modes.dev)), log10(max(dev.results$modes.dev)), 0, .6)
abline(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$modes.dev)), col="black", lwd=2)
abline(lm(dev.results$mean.or~log10(dev.results$modes.dev)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$modes.dev)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(dev.results$mean.or~log10(dev.results$modes.dev)))$coefficients[2,4]
counter <- counter+1

plot(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$modes.carb), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Carboniferous",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(carb.results$mean.or~log10(carb.results$modes.carb), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(carb.results$modes.carb)), log10(max(carb.results$modes.carb)), 0, .6)
abline(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$modes.carb)), col="black", lwd=2)
abline(lm(carb.results$mean.or~log10(carb.results$modes.carb)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$modes.carb)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(carb.results$mean.or~log10(carb.results$modes.carb)))$coefficients[2,4]
counter <- counter+1


plot(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$modes.perm), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Permian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(perm.results$mean.or~log10(perm.results$modes.perm), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(perm.results$modes.perm)), log10(max(perm.results$modes.perm)), 0, .6)
abline(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$modes.perm)), col="black", lwd=2)
abline(lm(perm.results$mean.or~log10(perm.results$modes.perm)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$modes.perm)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(perm.results$mean.or~log10(perm.results$modes.perm)))$coefficients[2,4]
counter <- counter+1

plot(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$modes.tri), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Triassic",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(tri.results$mean.or~log10(tri.results$modes.tri), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(tri.results$modes.tri)), log10(max(tri.results$modes.tri)), 0, .6)
abline(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$modes.tri)), col="black", lwd=2)
abline(lm(tri.results$mean.or~log10(tri.results$modes.tri)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$modes.tri)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(tri.results$mean.or~log10(tri.results$modes.tri)))$coefficients[2,4]
counter <- counter+1

plot(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$modes.jur), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Jurassic",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(jur.results$mean.or~log10(jur.results$modes.jur), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(jur.results$modes.jur)), log10(max(jur.results$modes.jur)), 0, .6)
abline(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$modes.jur)), col="black", lwd=2)
abline(lm(jur.results$mean.or~log10(jur.results$modes.jur)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$modes.jur)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(jur.results$mean.or~log10(jur.results$modes.jur)))$coefficients[2,4]
counter <- counter+1

plot(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$modes.cret), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Cretaceous",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(cret.results$mean.or~log10(cret.results$modes.cret), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(cret.results$modes.cret)), log10(max(cret.results$modes.cret)), 0, .6)
abline(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$modes.cret)), col="black", lwd=2)
abline(lm(cret.results$mean.or~log10(cret.results$modes.cret)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$modes.cret)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(cret.results$mean.or~log10(cret.results$modes.cret)))$coefficients[2,4]
counter <- counter+1

plot(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$modes.dani), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Danian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(dani.results$mean.or~log10(dani.results$modes.dani), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(dani.results$modes.dani)), log10(max(dani.results$modes.dani)), 0, .6)
abline(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$modes.dani)), col="black", lwd=2)
abline(lm(dani.results$mean.or~log10(dani.results$modes.dani)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$modes.dani)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(dani.results$mean.or~log10(dani.results$modes.dani)))$coefficients[2,4]
counter <- counter+1

plot(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$modes.palg), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Paleogene",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(palg.results$mean.or~log10(palg.results$modes.palg), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(palg.results$modes.palg)), log10(max(palg.results$modes.palg)), 0, .6)
abline(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$modes.palg)), col="black", lwd=2)
abline(lm(palg.results$mean.or~log10(palg.results$modes.palg)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$modes.palg)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(palg.results$mean.or~log10(palg.results$modes.palg)))$coefficients[2,4]
counter <- counter+1

plot((mod.results$mean.or-mod.results$mean.ex)~log10(mod.results$modes.mod), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Pleistocene",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points((mod.results$mean.or)~log10(mod.results$modes.mod), pch=23, col="black", bg="blue", cex=1)
clip(log10(min(mod.results$modes.mod)), log10(max(mod.results$modes.mod)), 0, .6)
abline(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$modes.mod)), col="black", lwd=2)
abline(lm((mod.results$mean.or)~log10(mod.results$modes.mod)), col="blue", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$modes.mod)))$coefficients[2,4]
pval_modes_orig_rate[counter] <- summary(lm(mod.results$mean.or~log10(mod.results$modes.mod)))$coefficients[2,4]
counter <- counter+1

dev.off()

p_vals <- cbind(pval_div_rate, pval_orig_rate, pval_modes_div_rate, pval_modes_orig_rate)
write.csv(p_vals, "p_values for regressions origination.csv")

#identify extinction and origination by stage and then code as binary variables
orex <- ddply(mtf, c("taxon_genus"), function(df)c(max(df$stage), min(df$stage)))
names(orex) <- c("taxon_genus", "or.stage", "ex.stage")

mtf.1 <- merge(mtf, orex, by="taxon_genus")
mtf.1$or <- 0
mtf.1$or[mtf.1$stage==mtf.1$or.stage] <- 1
mtf.1$ex <- 0
mtf.1$ex[mtf.1$stage==mtf.1$ex.stage] <- 1

mtf.1$duration <- mtf.1$or.stage - mtf.1$ex.stage

mtf.1 <- subset(mtf.1, duration<150 & class!="" & (ex==0 | or==0)) #remove, long-ranging genera and singletons

#count number of modes by class for selected time slices
mt <- ddply(mtf.1, c("class"), function(df)c(length(unique(na.omit(df$ecospace))), 
                                             length(unique(df$taxon_genus)),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==453]))), # Ordovician
                                             length(unique(df$taxon_genus[df$stage==453])),
                                             length(unique(df$taxon_genus[df$stage>=453])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==419.2]))), # Silurian
                                             length(unique(df$taxon_genus[df$stage==419.2])),
                                             length(unique(df$taxon_genus[df$stage>=419.2])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==358.9]))), # Devonian
                                             length(unique(df$taxon_genus[df$stage==358.9])),
                                             length(unique(df$taxon_genus[df$stage>=358.9])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==298.9]))), # Carboniferous
                                             length(unique(df$taxon_genus[df$stage==298.9])),
                                             length(unique(df$taxon_genus[df$stage>=298.9])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==252.2]))), # Permian
                                             length(unique(df$taxon_genus[df$stage==252.2])),
                                             length(unique(df$taxon_genus[df$stage>=252.2])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==201.3]))), # Triassic 
                                             length(unique(df$taxon_genus[df$stage==201.3])),
                                             length(unique(df$taxon_genus[df$stage>=201.3])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==145]))), # Jurassic
                                             length(unique(df$taxon_genus[df$stage==145])),
                                             length(unique(df$taxon_genus[df$stage>=145])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==66]))), # Cretaceous
                                             length(unique(df$taxon_genus[df$stage==66])),
                                             length(unique(df$taxon_genus[df$stage>=66])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==61.6]))), # Danian
                                             length(unique(df$taxon_genus[df$stage==61.6])),
                                             length(unique(df$taxon_genus[df$stage>=61.6])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==23.03]))), # Paleogene
                                             length(unique(df$taxon_genus[df$stage==23.03])),
                                             length(unique(df$taxon_genus[df$stage>=23.03])),
                                             
                                             length(unique(na.omit(df$ecospace[df$stage==0.0117]))), #Holocene
                                             length(unique(df$taxon_genus[df$stage==0.0117])),
                                             max(df$stage)))

names(mt) <- c("class", "modes.tot", "gen.tot", "modes.ord", "gen.ord", "gen.preord", "modes.sil", "gen.sil", "gen.presil", 
               "modes.dev", "gen.dev", "gen.predev", "modes.carb", "gen.carb", "gen.precarb", "modes.perm", "gen.perm", 
               "gen.preperm", "modes.tri", "gen.tri", "gen.pretri", "modes.jur", "gen.jur", "gen.prejur", 
               "modes.cret", "gen.cret", "gen.precret", "modes.dani", "gen.dani", "gen.predani", "modes.palg", "gen.palg", "gen.prepalg", "modes.mod", "gen.mod", "orig.class")

mtf.2 <- merge(mtf.1, mt, by="class")

min.gen.standing <- 2
min.gen.total <- 100


## MODERN ANALYSIS FOOTE RATES ##
mod <- ddply(mtf.2, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])), #extinction rate
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])), #origination rate
    length(df$ex[df$ex==0 & df$or==0]), #range through genera
    length(df$or[df$ex==0]), #surviving genera
    length(df$ex[df$or==0]))) #genera entering from previous stage
names(mod) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

mod.1 <- subset(mod, logex<100 & logor<100) #remove infinite rates and NA values
mod.2 <- ddply(mod.1, c("class"), function(df) #calculate average extinction and origination rates
  c(mean(df$logex), mean(df$logor)))
names(mod.2) <- c("class", "mean.ex", "mean.or")
mod.2$mean.or <- apply(mod.2[,2:3], MARGIN=1, FUN=max) 
mod.3 <- merge(mod.2, mt, by="class") #merge with other data by class (genus and mode counts)

mod.results <- subset(mod.3, gen.tot>=min.gen.total & gen.mod>=min.gen.standing & orig.class>400) #subset based on diversity criteria (standing and total diversity) 

## PALEOGENE ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=23.03)
palg <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(palg) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

palg.1 <- subset(palg, logex<100 & logor<100)
palg.2 <- ddply(palg.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(palg.2) <- c("class", "mean.ex", "mean.or")
palg.2$mean.or <- apply(palg.2[,2:3], MARGIN=1, FUN=max) 
palg.3 <- merge(palg.2, mt, by="class")

palg.results <- subset(palg.3, gen.tot>=min.gen.total & gen.palg>=min.gen.standing & orig.class>400)

## DANIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=61.6) 
dani <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(dani) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

dani.1 <- subset(dani, logex<100 & logor<100)
dani.2 <- ddply(dani.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(dani.2) <- c("class", "mean.ex", "mean.or")
dani.2$mean.or <- apply(dani.2[,2:3], MARGIN=1, FUN=max) 
dani.3 <- merge(dani.2, mt, by="class")

dani.results <- subset(dani.3, gen.tot>=min.gen.total & gen.dani>=min.gen.standing & orig.class>400)

## CRETACEOUS ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=66)
cret <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(cret) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

cret.1 <- subset(cret, logex<100 & logor<100)
cret.2 <- ddply(cret.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(cret.2) <- c("class", "mean.ex", "mean.or")
cret.2$mean.or <- apply(cret.2[,2:3], MARGIN=1, FUN=max) 
cret.3 <- merge(cret.2, mt, by="class")

cret.results <- subset(cret.3, gen.tot>=min.gen.total & gen.cret>=min.gen.standing & orig.class>400)

## JURASSIC ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=145)
jur <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(jur) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

jur.1 <- subset(jur, logex<100 & logor<100)
jur.2 <- ddply(jur.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(jur.2) <- c("class", "mean.ex", "mean.or")
jur.2$mean.or <- apply(jur.2[,2:3], MARGIN=1, FUN=max) 
jur.3 <- merge(jur.2, mt, by="class")

jur.results <- subset(jur.3, gen.tot>=min.gen.total & gen.jur>=min.gen.standing & orig.class>400)

## TRIASSIC ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=201.3)
tri <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(tri) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

tri.1 <- subset(tri, logex<100 & logor<100)
tri.2 <- ddply(tri.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(tri.2) <- c("class", "mean.ex", "mean.or")
tri.2$mean.or <- apply(tri.2[,2:3], MARGIN=1, FUN=max) 
tri.3 <- merge(tri.2, mt, by="class")

tri.results <- subset(tri.3, gen.tot>=min.gen.total & gen.tri>=min.gen.standing & orig.class>400)

## PERMIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=252.2) #changed from 252
perm <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(perm) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

perm.1 <- subset(perm, logex<100 & logor<100)
perm.2 <- ddply(perm.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(perm.2) <- c("class", "mean.ex", "mean.or")
perm.2$mean.or <- apply(perm.2[,2:3], MARGIN=1, FUN=max) 
perm.3 <- merge(perm.2, mt, by="class")

perm.results <- subset(perm.3, gen.tot>=min.gen.total & gen.perm>=min.gen.standing & orig.class>400)

## CARBONIFEROUS ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=298.9)
carb <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(carb) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

carb.1 <- subset(carb, logex<100 & logor<100)
carb.2 <- ddply(carb.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(carb.2) <- c("class", "mean.ex", "mean.or")
carb.2$mean.or <- apply(carb.2[,2:3], MARGIN=1, FUN=max) 
carb.3 <- merge(carb.2, mt, by="class")

carb.results <- subset(carb.3, gen.tot>=min.gen.total & gen.carb>=min.gen.standing & orig.class>400)

## DEVONIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=358.9)
dev <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(dev) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

dev.1 <- subset(dev, logex<100 & logor<100)
dev.2 <- ddply(dev.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(dev.2) <- c("class", "mean.ex", "mean.or")
dev.2$mean.or <- apply(dev.2[,2:3], MARGIN=1, FUN=max)
dev.3 <- merge(dev.2, mt, by="class")

dev.results <- subset(dev.3, gen.tot>=min.gen.total & gen.dev>=min.gen.standing & orig.class>400)

## SILURIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=419.2)
sil <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(sil) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

sil.1 <- subset(sil, logex<100 & logor<100)
sil.2 <- ddply(sil.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(sil.2) <- c("class", "mean.ex", "mean.or")
sil.2$mean.or <- apply(sil.2[,2:3], MARGIN=1, FUN=max) 
sil.3 <- merge(sil.2, mt, by="class")

sil.results <- subset(sil.3, gen.tot>=min.gen.total & gen.sil>=min.gen.standing & orig.class>400)

## ORDOVICIAN ANALYSIS FOOTE RATES ##
mtf.3 <- subset(mtf.2, stage>=453)
ord <- ddply(mtf.3, c("class", "stage"), function(df)
  c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])),
    -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])),
    length(df$ex[df$ex==0 & df$or==0]),
    length(df$or[df$ex==0]),
    length(df$ex[df$or==0])))
names(ord) <- c("class", "stage", "logex", "logor", "Nbt", "Nt", "Nb")

ord.1 <- subset(ord, logex<100 & logor<100)
ord.2 <- ddply(ord.1, c("class"), function(df)
  c(mean(df$logex), mean(df$logor)))
names(ord.2) <- c("class", "mean.ex", "mean.or")
ord.2$mean.or <- apply(ord.2[,2:3], MARGIN=1, FUN=max) 
ord.3 <- merge(ord.2, mt, by="class")

ord.results <- subset(ord.3, gen.tot>=min.gen.total & gen.ord>=min.gen.standing & orig.class>400)









## FIGURE FOR GENUS DIVERSITY vs Rates##
pdf("regr ext and diversification.pdf", h=8, w=8)
par(mfrow=c(4,5), mar=c(4,4,1,1), mgp=c(1.75,.6,0), pin=c(1,1), tcl=-.3, lab=c(4,4,1), cex.axis=.8, cex.main=.8)

y.min <- 0
y.max <- .8
x.min <- 0
x.max.1 <- log10(max(mod.results$gen.mod))+0.2
x.max.2 <- log10(max(mod.results$modes.mod))+0.2

pval_ext_rate <- rep(NA,10)
pval_div_rate <- rep(NA,10)
counter <- 1

plot((sil.results$mean.or-sil.results$mean.ex)~log10(sil.results$gen.sil), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Silurian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points((sil.results$mean.ex)~log10(sil.results$gen.sil), pch=23, col="black", bg="red", cex=1)
clip(log10(min(sil.results$gen.sil)), log10(max(sil.results$gen.sil)), 0, .6)
abline(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$gen.sil)), col="black", lwd=2)
abline(lm((sil.results$mean.ex)~log10(sil.results$gen.sil)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$gen.sil)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(sil.results$mean.ex~log10(sil.results$gen.sil)))$coefficients[2,4]
counter <- counter+1

plot(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$gen.dev), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Devonian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(dev.results$mean.ex~log10(dev.results$gen.dev), pch=23, col="black", bg="red", cex=1)
clip(log10(min(dev.results$gen.dev)), log10(max(dev.results$gen.dev)), 0, .6)
abline(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$gen.dev)), col="black", lwd=2)
abline(lm(dev.results$mean.ex~log10(dev.results$gen.dev)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$gen.dev)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(dev.results$mean.ex~log10(dev.results$gen.dev)))$coefficients[2,4]
counter <- counter+1

plot(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$gen.carb), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Carboniferous",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(carb.results$mean.ex~log10(carb.results$gen.carb), pch=23, col="black", bg="red", cex=1)
clip(log10(min(carb.results$gen.carb)), log10(max(carb.results$gen.carb)), 0, .6)
abline(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$gen.carb)), col="black", lwd=2)
abline(lm(carb.results$mean.ex~log10(carb.results$gen.carb)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$gen.carb)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(carb.results$mean.ex~log10(carb.results$gen.carb)))$coefficients[2,4]
counter <- counter+1

plot(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$gen.perm), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Permian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(perm.results$mean.ex~log10(perm.results$gen.perm), pch=23, col="black", bg="red", cex=1)
clip(log10(min(perm.results$gen.perm)), log10(max(perm.results$gen.perm)), 0, .6)
abline(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$gen.perm)), col="black", lwd=2)
abline(lm(perm.results$mean.ex~log10(perm.results$gen.perm)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$gen.perm)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(perm.results$mean.ex~log10(perm.results$gen.perm)))$coefficients[2,4]
counter <- counter+1


plot(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$gen.tri), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Triassic",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(tri.results$mean.ex~log10(tri.results$gen.tri), pch=23, col="black", bg="red", cex=1)
clip(log10(min(tri.results$gen.tri)), log10(max(tri.results$gen.tri)), 0, .6)
abline(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$gen.tri)), col="black", lwd=2)
abline(lm(tri.results$mean.ex~log10(tri.results$gen.tri)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$gen.tri)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(tri.results$mean.ex~log10(tri.results$gen.tri)))$coefficients[2,4]
counter <- counter+1

plot(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$gen.jur), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Jurassic",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(jur.results$mean.ex~log10(jur.results$gen.jur), pch=23, col="black", bg="red", cex=1)
clip(log10(min(jur.results$gen.jur)), log10(max(jur.results$gen.jur)), 0, .6)
abline(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$gen.jur)), col="black", lwd=2)
abline(lm(jur.results$mean.ex~log10(jur.results$gen.jur)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$gen.jur)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(jur.results$mean.ex~log10(jur.results$gen.jur)))$coefficients[2,4]
counter <- counter+1

plot(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$gen.cret), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Cretaceous",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(cret.results$mean.ex~log10(cret.results$gen.cret), pch=23, col="black", bg="red", cex=1)
clip(log10(min(cret.results$gen.cret)), log10(max(cret.results$gen.cret)), 0, .6)
abline(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$gen.cret)), col="black", lwd=2)
abline(lm(cret.results$mean.ex~log10(cret.results$gen.cret)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$gen.cret)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(cret.results$mean.ex~log10(cret.results$gen.cret)))$coefficients[2,4]
counter <- counter+1

plot(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$gen.dani), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Danian",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(dani.results$mean.ex~log10(dani.results$gen.dani), pch=23, col="black", bg="red", cex=1)
clip(log10(min(dani.results$gen.dani)), log10(max(dani.results$gen.dani)), 0, .6)
abline(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$gen.dani)), col="black", lwd=2)
abline(lm(dani.results$mean.ex~log10(dani.results$gen.dani)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$gen.dani)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(dani.results$mean.ex~log10(dani.results$gen.dani)))$coefficients[2,4]
counter <- counter+1

plot(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$gen.palg), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Paleogene",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points(palg.results$mean.ex~log10(palg.results$gen.palg), pch=23, col="black", bg="red", cex=1)
clip(log10(min(palg.results$gen.palg)), log10(max(palg.results$gen.palg)), 0, .6)
abline(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$gen.palg)), col="black", lwd=2)
abline(lm(palg.results$mean.ex~log10(palg.results$gen.palg)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$gen.palg)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(palg.results$mean.ex~log10(palg.results$gen.palg)))$coefficients[2,4]
counter <- counter+1

plot((mod.results$mean.or-mod.results$mean.ex)~log10(mod.results$gen.mod), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 genera", ylab="mean rate (BC)", main="Pleistocene",
     xlim=c(x.min,x.max.1), ylim=c(y.min, y.max), type="p")
points((mod.results$mean.ex)~log10(mod.results$gen.mod), pch=23, col="black", bg="red", cex=1)
clip(log10(min(mod.results$gen.mod)), log10(max(mod.results$gen.mod)), 0, .6)
abline(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$gen.mod)), col="black", lwd=2)
abline(lm((mod.results$mean.ex)~log10(mod.results$gen.mod)), col="red", lwd=2)
pval_div_rate[counter] <- summary(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$gen.mod)))$coefficients[2,4]
pval_ext_rate[counter] <- summary(lm(mod.results$mean.ex~log10(mod.results$gen.mod)))$coefficients[2,4]
counter <- counter+1


## FIGURE FOR Modes vs Rates ##


pval_modes_ext_rate <- rep(NA,10)
pval_modes_div_rate <- rep(NA,10)
counter <- 1

y.min <- 0
x.min <- 0
x.max.1 <- log10(max(mod.results$gen.mod))+0.2
x.max.2 <- log10(max(mod.results$modes.mod))+0.2

plot((sil.results$mean.or-sil.results$mean.ex)~log10(sil.results$modes.sil), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Silurian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points((sil.results$mean.ex)~log10(sil.results$modes.sil), pch=23, col="black", bg="red", cex=1)
clip(log10(min(sil.results$modes.sil)), log10(max(sil.results$modes.sil)), 0, .6)
abline(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$modes.sil)), col="black", lwd=2)
abline(lm((sil.results$mean.ex)~log10(sil.results$modes.sil)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(sil.results$mean.or-sil.results$mean.ex~log10(sil.results$modes.sil)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(sil.results$mean.ex~log10(sil.results$modes.sil)))$coefficients[2,4]
counter <- counter+1

plot(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$modes.dev), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Devonian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(dev.results$mean.ex~log10(dev.results$modes.dev), pch=23, col="black", bg="red", cex=1)
clip(log10(min(dev.results$modes.dev)), log10(max(dev.results$modes.dev)), 0, .6)
abline(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$modes.dev)), col="black", lwd=2)
abline(lm(dev.results$mean.ex~log10(dev.results$modes.dev)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(dev.results$mean.or-dev.results$mean.ex~log10(dev.results$modes.dev)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(dev.results$mean.ex~log10(dev.results$modes.dev)))$coefficients[2,4]
counter <- counter+1

plot(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$modes.carb), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Carboniferous",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(carb.results$mean.ex~log10(carb.results$modes.carb), pch=23, col="black", bg="red", cex=1)
clip(log10(min(carb.results$modes.carb)), log10(max(carb.results$modes.carb)), 0, .6)
abline(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$modes.carb)), col="black", lwd=2)
abline(lm(carb.results$mean.ex~log10(carb.results$modes.carb)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(carb.results$mean.or-carb.results$mean.ex~log10(carb.results$modes.carb)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(carb.results$mean.ex~log10(carb.results$modes.carb)))$coefficients[2,4]
counter <- counter+1


plot(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$modes.perm), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Permian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(perm.results$mean.ex~log10(perm.results$modes.perm), pch=23, col="black", bg="red", cex=1)
clip(log10(min(perm.results$modes.perm)), log10(max(perm.results$modes.perm)), 0, .6)
abline(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$modes.perm)), col="black", lwd=2)
abline(lm(perm.results$mean.ex~log10(perm.results$modes.perm)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(perm.results$mean.or-perm.results$mean.ex~log10(perm.results$modes.perm)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(perm.results$mean.ex~log10(perm.results$modes.perm)))$coefficients[2,4]
counter <- counter+1

plot(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$modes.tri), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Triassic",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(tri.results$mean.ex~log10(tri.results$modes.tri), pch=23, col="black", bg="red", cex=1)
clip(log10(min(tri.results$modes.tri)), log10(max(tri.results$modes.tri)), 0, .6)
abline(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$modes.tri)), col="black", lwd=2)
abline(lm(tri.results$mean.ex~log10(tri.results$modes.tri)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(tri.results$mean.or-tri.results$mean.ex~log10(tri.results$modes.tri)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(tri.results$mean.ex~log10(tri.results$modes.tri)))$coefficients[2,4]
counter <- counter+1

plot(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$modes.jur), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Jurassic",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(jur.results$mean.ex~log10(jur.results$modes.jur), pch=23, col="black", bg="red", cex=1)
clip(log10(min(jur.results$modes.jur)), log10(max(jur.results$modes.jur)), 0, .6)
abline(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$modes.jur)), col="black", lwd=2)
abline(lm(jur.results$mean.ex~log10(jur.results$modes.jur)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(jur.results$mean.or-jur.results$mean.ex~log10(jur.results$modes.jur)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(jur.results$mean.ex~log10(jur.results$modes.jur)))$coefficients[2,4]
counter <- counter+1

plot(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$modes.cret), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Cretaceous",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(cret.results$mean.ex~log10(cret.results$modes.cret), pch=23, col="black", bg="red", cex=1)
clip(log10(min(cret.results$modes.cret)), log10(max(cret.results$modes.cret)), 0, .6)
abline(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$modes.cret)), col="black", lwd=2)
abline(lm(cret.results$mean.ex~log10(cret.results$modes.cret)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(cret.results$mean.or-cret.results$mean.ex~log10(cret.results$modes.cret)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(cret.results$mean.ex~log10(cret.results$modes.cret)))$coefficients[2,4]
counter <- counter+1

plot(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$modes.dani), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Danian",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(dani.results$mean.ex~log10(dani.results$modes.dani), pch=23, col="black", bg="red", cex=1)
clip(log10(min(dani.results$modes.dani)), log10(max(dani.results$modes.dani)), 0, .6)
abline(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$modes.dani)), col="black", lwd=2)
abline(lm(dani.results$mean.ex~log10(dani.results$modes.dani)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(dani.results$mean.or-dani.results$mean.ex~log10(dani.results$modes.dani)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(dani.results$mean.ex~log10(dani.results$modes.dani)))$coefficients[2,4]
counter <- counter+1

plot(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$modes.palg), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Paleogene",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points(palg.results$mean.ex~log10(palg.results$modes.palg), pch=23, col="black", bg="red", cex=1)
clip(log10(min(palg.results$modes.palg)), log10(max(palg.results$modes.palg)), 0, .6)
abline(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$modes.palg)), col="black", lwd=2)
abline(lm(palg.results$mean.ex~log10(palg.results$modes.palg)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(palg.results$mean.or-palg.results$mean.ex~log10(palg.results$modes.palg)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(palg.results$mean.ex~log10(palg.results$modes.palg)))$coefficients[2,4]
counter <- counter+1

plot((mod.results$mean.or-mod.results$mean.ex)~log10(mod.results$modes.mod), pch=23, col="black", bg="yellow", las=1, cex=1,
     xlab="log10 modes", ylab="mean rate (BC)", main="Pleistocene",
     xlim=c(x.min,x.max.2), ylim=c(y.min, y.max), type="p")
points((mod.results$mean.ex)~log10(mod.results$modes.mod), pch=23, col="black", bg="red", cex=1)
clip(log10(min(mod.results$modes.mod)), log10(max(mod.results$modes.mod)), 0, .6)
abline(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$modes.mod)), col="black", lwd=2)
abline(lm((mod.results$mean.ex)~log10(mod.results$modes.mod)), col="red", lwd=2)
pval_modes_div_rate[counter] <- summary(lm(mod.results$mean.or-mod.results$mean.ex~log10(mod.results$modes.mod)))$coefficients[2,4]
pval_modes_ext_rate[counter] <- summary(lm(mod.results$mean.ex~log10(mod.results$modes.mod)))$coefficients[2,4]
counter <- counter+1

dev.off()

p_vals <- cbind(pval_div_rate, pval_ext_rate, pval_modes_div_rate, pval_modes_ext_rate)
write.csv(p_vals, "p_values for regressions extinction.csv")




## ********************************************************************************
## Assess the relationship between # ecomodes and # genera through time: Figs. 2, S1, S2, Tables S1,S2
## ********************************************************************************
#rm(list=ls())

library(reshape2)
library(viridis)
library(nlme)
library(plotrix)

## ************************************************************************

## Establish core phyla, and remove all others:
core.phyla<-c('Hemichordata','Porifera','Cnidaria','Bryozoa','Gastropoda','Cephalopoda',
				'Bivalvia','Mollusca','Echinodermata','Brachiopoda','Arthropoda',
				'Chordata')
master.core.phyla<-master[master$phylum %in% core.phyla,]

# Drop levels, and create a bin category for use with aggregate function to 
# collapse df across stages
master.core.phyla$phylum<-factor(master.core.phyla$phylum)
master.core.phyla$class<-factor(master.core.phyla$class)
master.core.phyla$bin<-1 # create a bin for counting durring data aggregation

## Create array of number of genera in each class X stage X ecomode combination
XX<-acast(master.core.phyla, class ~  stage ~ ecospace , sum)	

## Establish data arrays from XX
aa.paleo<-XX[-1,,] # Eliminate all genera with unknown classes (blank class designation)

## Convert aa.paleo to summary of number of modes & number of genera
mm.gen<-apply(aa.paleo, c(1,2), sum)
mm.modes<-apply(aa.paleo>0, c(1,2), sum)

# Establish counts and colors to be used below
nclass<-dim(aa.paleo)[1]
nstage<-dim(aa.paleo)[2]

# Convert matrixes to dataframes for analysis
df.modes<-melt(mm.modes)
df.gen<-melt(mm.gen)
df.paleo<-data.frame(df.modes, genera=df.gen$value)
names(df.paleo)<-c('class','stage','modes', 'genera')

## Generate new dataframe with relevant information for analysis
df.paleo$lmodes<-log10(df.paleo$modes)
df.paleo$lgenera<-log10(df.paleo$genera)
df.paleo$fstage<-factor(df.paleo$stage)
df.paleo$lgenera2<-log10(df.paleo$genera)^2

# Name the different time periods
df.paleo.narm<-df.paleo[-which(df.paleo$lmodes==-Inf),] # Remove -Infs
df.paleo.narm$period<-""
df.paleo.narm$period[df.paleo.narm$stage <444]<-"Sil-Dev"
df.paleo.narm$period[df.paleo.narm$stage <358]<-"Carb-Perm"
df.paleo.narm$period[df.paleo.narm$stage <252]<-"Triassic"
df.paleo.narm$period[df.paleo.narm$stage <201]<-"Jur-Cret"
df.paleo.narm$period[df.paleo.narm$stage <66]<-"Ceno"

# Remove stages before Sil-Dev
df.paleo.period<-df.paleo.narm[!df.paleo.narm$period=="",]
df.paleo.period$period<-factor(df.paleo.period$period)


# Assign relevant colors
cols<-rev(viridis(5))[c(4,1,2,5,3)]


plot.stage.lin.0<-function(stage=1){
	hold<-df.paleo.period[df.paleo.period$fstage==levels(df.paleo.period$fstage)[stage],]
	varf<-varExp(1, ~ lgenera)

 	mod<-gls(lmodes~0+lgenera, data=hold, weights=varf)
	
	colzz.hold<-cols[as.numeric(hold$period)[1]]
	preds<-predict(mod)
	ii<-order(hold$lgenera)
	points(hold$lgenera[ii], preds[ii], type='l', col=colzz.hold)
	
	c(0,coef(mod),confint(mod)[1,])
	
}


full.modern <- read.table(file="class.csv", sep=",", header=T)

varf<-varExp(1, ~ lgenera)
mod.modern.0<-gls(lmodes~0+lgenera, data= full.modern, weights=varf)

modern.coefs.0<-c(0, coef(mod.modern.0),confint(mod.modern.0)[1,])

period.names<-levels(df.paleo.period$period)[c(2,3,5,1,4)]
period.names	[1]<-"Paleo-Neo"
full.names<-paste(period.names, c("(66.0 - 0.0117 Ma)", "(201.3 - 66.0 Ma)", "(252.2 - 201.3 Ma)" ,"(358.9 - 252.2 Ma)", "(443.4 - 358.9 Ma)"))



## ********************************************************************************


pdf("figure-modern.pdf", width=6, height=6)

	par(mfrow=c(1,1))
	par(mar=c(5,6,1.5,1), oma=c(0,1,1,0))
		
	plot(full.modern $lgenera, full.modern $lmodes, pch=16, cex=1, ylim=c(0,1.6), main="", xlab="", ylab="", las=1, xlim=c(0,4))

	mtext(expression("Modern Ecological Diversity"), 2, line=4, cex=1.5) 
	mtext(expression((Log[10]*" "*modes)), 2, line=2.3, cex=1.5)
	mtext(expression("Modern Richness ("*Log[10]*" "*genera*")"),1, cex=1.5, line=3.5)
	
	lab.force<-sapply(c('sar','rhy'), function(i) which(full.modern $abbrev==i))
	force.vec<-c(0.04)
	lab.force2<-sapply(c('cep'), function(i) which(full.modern $abbrev==i))
	force.vec2<-c(0.04)
	lab.force3<-sapply(c('cho'), function(i) which(full.modern $abbrev==i))
	force.vec3<-c(-0.04)
	
	all.force<-c(lab.force, lab.force2, lab.force3)
	
	preds.modern<-predict(mod.modern.0)
	ii<-order(full.modern $lgenera)
	points(full.modern $lgenera[ii], preds.modern[ii], type='l', col='dodgerblue', lwd=3)
	
	
	thigmophobe.labels(full.modern $lgenera[-all.force], full.modern $lmodes[-all.force], full.modern $abbrev[-all.force], cex=1)
	text(full.modern $lgenera[lab.force], full.modern $lmodes[lab.force]+force.vec, full.modern $abbrev[lab.force], cex=1)
	text(full.modern $lgenera[lab.force2]+force.vec2, full.modern $lmodes[lab.force2]+force.vec2, full.modern $abbrev[lab.force2], cex=1)
	text(full.modern $lgenera[lab.force3]+4*force.vec3, full.modern $lmodes[lab.force3]+force.vec3, full.modern $abbrev[lab.force3], cex=1)

dev.off()



pdf("figure2-twopanel-no-varcovar.pdf", width=11, height=5.5)

## Panel A
	par(mfrow=c(1,2))
	par(mar=c(5,6,1.5,1), oma=c(0,1,1,0))


plot(c(0,4), c(0,1), col='white',
			xlab="", ylab="", cex.lab=1.5, las=1, ylim=c(0,1.6))
	
	mtext(expression("Ecological Diversity"), 2, line=4, cex=1.5) 
	mtext(expression((Log[10]*" "*modes)), 2, line=2.3, cex=1.5)
	mtext(expression("Richness ("*Log[10]*" "*genera*")"),1, cex=1.5, line=3.5)


coef.df<-data.frame(t(plot.stage.lin.0(1)))

for(i in 1:80) {coef.df[i,]<-plot.stage.lin.0(i)}

	legend('topleft', legend=c('Modern', full.names), 
			lwd=4, lty=c( rep(1, 6)), col= c('dodgerblue',rev(viridis(5))), seg.len = 3)

	preds.modern<-predict(mod.modern.0)
	ii<-order(full.modern $lgenera)
	points(full.modern $lgenera[ii], preds.modern[ii], type='l', col='dodgerblue', lwd=3)

	mtext("a)", side=3, at=-0.6, cex=2)

	
	## Panel B
	
	stage.col<-character()
	for (stage in 1:80){
		hold<-df.paleo.period[df.paleo.period$fstage==levels(df.paleo.period$fstage)[stage],]
		colzz.hold<-cols[as.numeric(hold$period)[1]]
		stage.col[stage]<-colzz.hold
	}
	
	plot(unique(df.paleo.period$stage), coef.df[,2], ylim=c(0, 0.45), ylab='', xlab="", col=stage.col, xlim=c(450,0), las=1)
	mtext(expression("Slope"), 2, line=3, cex=1.5)
	mtext(expression("Time (Ma)"),1, cex=1.7, line=3.5)
	
	points(0, modern.coefs.0[2], col='dodgerblue', lwd=3)
	
	arrows(x0=unique(df.paleo.period$stage), y0=coef.df[,3], y1=coef.df[,4], length=0, col=stage.col)
	arrows(x0=0, y0= modern.coefs.0[3], y1= modern.coefs.0[4], length=0, col='dodgerblue', lwd=3)
	abline(h=0, lty=1)
	
		mtext("b)", side=3, at=520, cex=2)
	
	abline(v=c(444, 358, 252, 201, 66), lty=2)

dev.off()

## ************************************************************************


## Secondary analyses to see what drives changes in linear effect
df.stages<-aggregate(modes ~ period + stage , FUN=mean, data=df.paleo.period)
df.stages$coef<-coef.df[,2] # Extract slope term

# Put period levels in order for easier interpretation
df.stages$period<-factor(df.stages$period, levels=c(levels(df.stages$period)[c(2,3,5,1,4)]))

mod.coef<-lm(coef~period + stage, data= df.stages)
summary(mod.coef) # Table S1

drop1(mod.coef, test='F') # Table S2



## ********************************************************************************
## Assess if early functional differentiation predicts later genus richness 
##and vice-versa: Figs. S11, S15
## ********************************************************************************

mtf <- master.time.frame

classes <- sort(unique(mtf$class)) #get a list of unique classes for later use
classes <- classes[classes!=""] #remove genera without class assignment
stages <- sort(unique(mtf$stage)) #get a list of unique stages for later use
stages <- stages[stages<486] #limit analysis to Ordovician through Recent

genera <- data.frame(matrix(nrow=length(classes), ncol=length(stages))) #define matrix of genera per class and stage
names(genera) <- stages #give column names to the matrix, which are the stages in Mya
rownames(genera) <- classes #give row names to the matrix, which are Linnaean classes

### loop through all classes and all stages to find out how many genera fall in each combination
for (i in 1:length(classes)) {
  data.i <- subset(mtf, class==classes[i]) #subset data to class of interest
  for (j in 1:length(stages)) {
    data.ij <- subset(data.i, stage==stages[j]) #subset further to the stage of interest
    genera[i,j] <- length(data.ij$taxon_genus) #record the number of genera for the class by stage combination
    print(i) #print i to see that code is running
    print(j) #print j to see that code is running
  }
}

ecomodes <- data.frame(matrix(nrow=length(classes), ncol=length(stages))) #define matrix of ecomodes per class and stage
names(ecomodes) <- stages #give column names to the matrix, which are stages in Mya
rownames(ecomodes) <- classes #give row names to the matrix, which are Linnaean classes

for (i in 1:length(classes)) {
  data.i <- subset(mtf, class==classes[i]) #subset data to class of interest
  for (j in 1:length(stages)) {
    data.ij <- subset(data.i, stage==stages[j]) #subset further to the stage of interest
    ecomodes[i,j] <- length(unique(data.ij$ecospace)) #record the number of unique ecomodes for the class by stage combination
    print(i) #print i to see that code is running
    print(j) #print j to see that code is running
  }
}

genera[genera<10] <- NA #turn values smaller than 10 into NAs so that they are not used in the correlation analysis
ecomodes[ecomodes==0] <- NA #turn zero values into NAs so that they are not used in the correlation analysis

gen <- log10(genera) #log-transform data for analysis
eco <- log10(ecomodes) #log-transform data for analysis

correlations <- data.frame(matrix(nrow=length(stages), ncol=4)) #define output matrix for correlation analysis
cor.total <- data.frame(matrix(nrow=1, ncol=4)) #define final output matrix

for (i in 1:length(stages)) {
  for (j in 1:length(stages)) {
    cor.test.i <- cor.test(gen[,i], eco[,j], method="pearson", use = "pairwise.complete.obs")
    correlations[j,1] <- cor.test.i$estimate #pearson r calculation
    correlations[j,2] <- stages[i] #stage from which genus diversity is used
    correlations[j,3] <- stages[j] #stage from which ecomode count is used
    correlations[j,4] <- cor.test.i$p.value #extract p.value
  }
  cor.total <- rbind(cor.total, correlations) #combine all output
}

names(cor.total) <- c("pearson", "genera.stage", "ecomode.stage", "p.value") #name output variables
cor.total$color <- "blue" #define color for positive correlations
cor.total$color[cor.total$pearson<0] <- "red" #define color for inverse correlations

##plot correlation coefficients as a scatterplot where point size reflects strength of correlation
par(xpd=FALSE)
par(mar=c(4,4,4,1))
par(mfrow=c(1,1))
pdf("pearson correlations.pdf", h=5, w=5)
plot(cor.total$ecomode.stage ~ cor.total$genera.stage, pch=16, cex=abs(cor.total$pearson)*1.5, col=cor.total$color,
     xlim=c(487,0), ylim=c(487,0), las=1, ylab="Interval for ecomode count (Mya)", xlab="Interval for genus count (Mya)")
segments(252,485.4,252,0, col="black", lwd=1)
segments(66,485.4,66,0, col="black", lwd=1)
segments(485.4,252,0,252, col="black", lwd=1)
segments(485.4,66,0,66, col="black", lwd=1)
segments(485.4,485.4,0,0, col="black")
par(xpd=NA)
legend(500, -70, pch=16, col="black", pt.cex=0.3, cex=0.8, legend=c("r=0.2"), bty="n")
legend(400, -70, pch=16, col="black", pt.cex=0.75, cex=0.8, legend=c("r=0.5"), bty="n")
legend(300, -70, pch=16, col="black", pt.cex=1.5, cex=0.8, legend=c("r=1.0"), bty="n")
legend(200, -70, pch=16, col="blue", pt.cex=0.75, cex=0.8, legend=c("r>0"), bty="n")
legend(100, -70, pch=16, col="red", pt.cex=0.75, cex=0.8, legend=c("r<0"), bty="n")

dev.off()


names(cor.total) <- c("pearson", "genera.stage", "ecomode.stage", "p.value") #name output variables
cor.total$color <- "blue" #define color for positive correlations
cor.total$color[cor.total$pearson<0] <- "red" #define color for inverse correlations
cor.total$color[cor.total$p.value>0.05] <- NA #do not plot points that are not significant at a=0.05
##plot correlation coefficients as a scatterplot where point size reflects strength of correlation
par(xpd=FALSE)
par(mar=c(4,4,4,1))
par(mfrow=c(1,1))
#pdf("pearson correlations significant.pdf", h=5, w=5)
plot(cor.total$ecomode.stage ~ cor.total$genera.stage, pch=16, cex=abs(cor.total$pearson)*1.5, col=cor.total$color,
     xlim=c(487,0), ylim=c(487,0), las=1, ylab="Interval for ecomode count (Mya)", xlab="Interval for genus count (Mya)")
segments(252,485.4,252,0, col="black", lwd=1)
segments(66,485.4,66,0, col="black", lwd=1)
segments(485.4,252,0,252, col="black", lwd=1)
segments(485.4,66,0,66, col="black", lwd=1)
segments(485.4,485.4,0,0, col="black")
par(xpd=NA)
legend(500, -70, pch=16, col="black", pt.cex=0.3, cex=0.8, legend=c("r=0.2"), bty="n")
legend(400, -70, pch=16, col="black", pt.cex=0.75, cex=0.8, legend=c("r=0.5"), bty="n")
legend(300, -70, pch=16, col="black", pt.cex=1.5, cex=0.8, legend=c("r=1.0"), bty="n")
legend(200, -70, pch=16, col="blue", pt.cex=0.75, cex=0.8, legend=c("r>0"), bty="n")
legend(100, -70, pch=16, col="red", pt.cex=0.75, cex=0.8, legend=c("r<0"), bty="n")

dev.off()


## ********************************************************************************
## Assess is pre-existing functional differentiation predicts survival during 
## mass extinction better than during background intervals: Figs. S11,S5,S6,S7
## ********************************************************************************

## Read stage data
stages<-read.csv("/Users/knope/Desktop/timescale_final.csv")
# get rid of Pleistocene and Cambrian 
stages.reduced <- subset(stages, period != 'Cambrian' & interval_name != 'Pleistocene')


## ************************************************************************

## Establish core phyla, and remove all others:
core.phyla<-c('Hemichordata','Porifera','Cnidaria','Bryozoa','Mollusca','Echinodermata','Brachiopoda','Arthropoda','Chordata')
master.core.phyla<-master[master$phylum %in% core.phyla & master$class != '',] # drop taxa without class assignments

# Drop levels, and create a bin category for use with aggregate function to 
# collapse df across stages
master.core.phyla$phylum<-factor(master.core.phyla$phylum)
master.core.phyla$class<-factor(master.core.phyla$class)
master.core.phyla$bin<-1

## add column that marks the first and last occurrence
master.core.phyla$first_occurrence <- 0 # set all to zero at first
genus_fad <- data.frame(fad=tapply(master.core.phyla$stage, master.core.phyla$taxon_genus, max)) # calculate the fad interval for each genus
genus_fad$taxon_genus <- rownames(genus_fad)
master.core.phyla$first_occurrence[with(master.core.phyla, paste(taxon_genus, stage)) %in% with(genus_fad, paste(taxon_genus, fad)) ] <- 1

master.core.phyla$last_occurrence <- 0
genus_lad <- data.frame(lad=tapply(master.core.phyla$stage, master.core.phyla$taxon_genus, min)) # calculate the lad interval for each genus
genus_lad$taxon_genus <- rownames(genus_lad)
master.core.phyla$last_occurrence[with(master.core.phyla, paste(taxon_genus, stage)) %in% with(genus_lad, paste(taxon_genus, lad)) ] <- 1

## add column to stages for counting bins.
stages$bin.number <- rev(1:nrow(stages))

## create data.frame with pertenant genus-level data
genera <- unique(master.core.phyla[,match(c('phylum','class','taxon_genus','ecospace'),colnames(master.core.phyla))])
## calcualte genus ranges & apply total number of ecomodes for each class
genera$fad_age=tapply(master.core.phyla$stage, master.core.phyla$taxon_genus, max)
genera$lad_age=tapply(master.core.phyla$stage, master.core.phyla$taxon_genus, min)
## adjust fad_age, to account for master 'stage' variable is the interval top
genera$fad_age = stages$age_bottom[match(genera$fad_age, stages$age_top)]

# number of unique ecomodes per phylum-class combination (combined with phylum b/c there are a few genera with no class assignments)
class.modes <- data.frame(modes=tapply(master.core.phyla$ecospace, paste(master.core.phyla$phylum, master.core.phyla$class), function(x){length(unique(x))}))
class.modes$phylum_class <- rownames(class.modes)
genera$class_modes <- class.modes$modes[match(paste(genera$phylum, genera$class), class.modes$phylum_class) ]


## ********************************************************************************
## LOGISTIC REGRESSIONS
## ********************************************************************************
# minimum number of events to run regression
minEvents <- 5

# extinction & origination
class.diff.q <- data.frame()
class.diff.p <- data.frame()
class.diff.q.sansceph <- data.frame()
class.diff.p.sansceph <- data.frame()
cor.s.mode <- data.frame()
class.modes.time <- data.frame(matrix(NA, nrow=nrow(stages.reduced), ncol=length(levels(master.core.phyla$class)), dimnames=list(stages.reduced$interval_name, levels(master.core.phyla$class))))
class.genera.time <- class.modes.time
class.extinction.time <- class.modes.time
class.origination.time <- class.modes.time
cor.ext.gen.mode <- data.frame()
cor.orig.gen.mode <- data.frame()
cor.ext.gen.mode.sansceph <- data.frame()
cor.orig.gen.mode.sansceph <- data.frame()
for(i in 1:nrow(stages.reduced)) {
	# get extant taxa
	temp <- subset(master.core.phyla, stage==stages.reduced$age_top[i])
	classCounts <- as.numeric(table(temp$class))
	modeCounts <- as.numeric(tapply(temp$ecospace, temp$class, function(x){return(length(unique(x)))}))
	tempExtinction <- tapply(temp$last_occurrence, temp$class, sum, na.rm=T)/table(temp$class)
	tempOrigination <- tapply(temp$first_occurrence, temp$class, sum, na.rm=T)/table(temp$class)
	class.genera.time[i,] <- as.numeric(classCounts)
	class.modes.time[i,] <- as.numeric(modeCounts)
	class.extinction.time[i,] <- as.numeric(tempExtinction)
	class.origination.time[i,] <- as.numeric(tempOrigination)
	
	temp$classCount <- scale(classCounts[match(temp$class, levels(temp$class))])[,1]
	temp$modeCount <- scale(modeCounts[match(temp$class, levels(temp$class))])[,1]

	nTaxa <- nrow(temp)
	
	classInclude <- classCounts > 0
	# extinction
	if(var(tempExtinction[classInclude]) > 0) {
		y <- tempExtinction[classInclude]
		xMode <- scale(modeCounts[classInclude])[,1]
		xGenera <- scale(classCounts[classInclude])[,1]
		myGlm <- lm(y ~ xMode + xGenera, na.action='na.omit')
		confInts <- confint(myGlm)
		result <- data.frame(nTaxa, 'modeCoef'=myGlm$coefficients[match('xMode',names(myGlm$coefficients))], 'modeMinus'=confInts[match('xMode', rownames(confInts)),1], 'modePlus'=confInts[match('xMode', rownames(confInts)),2], 'generaCoef'=myGlm$coefficients[match('xGenera',names(myGlm$coefficients))], 'generaMinus'=confInts[match('xGenera', rownames(confInts)),1], 'generaPlus'=confInts[match('xGenera', rownames(confInts)),2])	
	} else {
		result <- data.frame(nTaxa, 'modeCoef'=NA, 'modeMinus'=NA, 'modePlus'=NA, 'generaCoef'=NA, 'generaMinus'=NA, 'generaPlus'=NA)
	}
	cor.ext.gen.mode <- rbind(cor.ext.gen.mode, result)
	
	# origination
	if(var(tempOrigination[classInclude]) > 0) {
		y <- tempOrigination[classInclude]
		xMode <- scale(modeCounts[classInclude])[,1]
		xGenera <- scale(classCounts[classInclude])[,1]
		myGlm <- lm(y ~ xMode + xGenera, na.action='na.omit')
		confInts <- confint(myGlm)
		result <- data.frame(nTaxa, 'modeCoef'=myGlm$coefficients[match('xMode',names(myGlm$coefficients))], 'modeMinus'=confInts[match('xMode', rownames(confInts)),1], 'modePlus'=confInts[match('xMode', rownames(confInts)),2], 'generaCoef'=myGlm$coefficients[match('xGenera',names(myGlm$coefficients))], 'generaMinus'=confInts[match('xGenera', rownames(confInts)),1], 'generaPlus'=confInts[match('xGenera', rownames(confInts)),2])	
	} else {
		result <- data.frame(nTaxa, 'modeCoef'=NA, 'modeMinus'=NA, 'modePlus'=NA, 'generaCoef'=NA, 'generaMinus'=NA, 'generaPlus'=NA)
	}
	cor.orig.gen.mode <- rbind(cor.orig.gen.mode, result)
	
	### exclude cephalopods
	# extinction
	classInclude[match("Cephalopoda", rownames(tempExtinction))] <- FALSE
	if(var(tempExtinction[classInclude]) > 0) {
		y <- tempExtinction[classInclude]
		xMode <- scale(modeCounts[classInclude])[,1]
		xGenera <- scale(classCounts[classInclude])[,1]
		myGlm <- lm(y ~ xMode + xGenera, na.action='na.omit')
		confInts <- confint(myGlm)
		result <- data.frame(nTaxa, 'modeCoef'=myGlm$coefficients[match('xMode',names(myGlm$coefficients))], 'modeMinus'=confInts[match('xMode', rownames(confInts)),1], 'modePlus'=confInts[match('xMode', rownames(confInts)),2], 'generaCoef'=myGlm$coefficients[match('xGenera',names(myGlm$coefficients))], 'generaMinus'=confInts[match('xGenera', rownames(confInts)),1], 'generaPlus'=confInts[match('xGenera', rownames(confInts)),2])	
	} else {
		result <- data.frame(nTaxa, 'modeCoef'=NA, 'modeMinus'=NA, 'modePlus'=NA, 'generaCoef'=NA, 'generaMinus'=NA, 'generaPlus'=NA)
	}
	cor.ext.gen.mode.sansceph <- rbind(cor.ext.gen.mode.sansceph, result)
	
	# origination
	if(var(tempOrigination[classInclude]) > 0) {
		y <- tempOrigination[classInclude]
		xMode <- scale(modeCounts[classInclude])[,1]
		xGenera <- scale(classCounts[classInclude])[,1]
		myGlm <- lm(y ~ xMode + xGenera, na.action='na.omit')
		confInts <- confint(myGlm)
		result <- data.frame(nTaxa, 'modeCoef'=myGlm$coefficients[match('xMode',names(myGlm$coefficients))], 'modeMinus'=confInts[match('xMode', rownames(confInts)),1], 'modePlus'=confInts[match('xMode', rownames(confInts)),2], 'generaCoef'=myGlm$coefficients[match('xGenera',names(myGlm$coefficients))], 'generaMinus'=confInts[match('xGenera', rownames(confInts)),1], 'generaPlus'=confInts[match('xGenera', rownames(confInts)),2])	
	} else {
		result <- data.frame(nTaxa, 'modeCoef'=NA, 'modeMinus'=NA, 'modePlus'=NA, 'generaCoef'=NA, 'generaMinus'=NA, 'generaPlus'=NA)
	}
	cor.orig.gen.mode.sansceph <- rbind(cor.orig.gen.mode.sansceph, result)
	####
	
	# extinction
	logOddsS <- NA
	pValueS <- NA
	ciPlusS <- NA
	ciMinusS <- NA
	logOddsMode <- NA
	pValueMode <- NA
	ciPlusMode <- NA
	ciMinusMode <- NA
		
	nEvents <- sum(temp$last_occurrence)
	pct <- nEvents/nTaxa
	if(nEvents >= minEvents & nTaxa - nEvents >= minEvents) {
		myGlm <- glm(last_occurrence ~ classCount + modeCount, data=temp, family='binomial')
		confInts <- confint(myGlm)
		logOddsS <- summary(myGlm)$coefficients[2,1] 
		pValueS <- summary(myGlm)$coefficients[2,4]
		ciPlusS <- confInts[2,1]
		ciMinusS <- confInts[2,2]
		
		logOddsMode <- summary(myGlm)$coefficients[3,1] 
		pValueMode <- summary(myGlm)$coefficients[3,4]
		ciPlusMode <- confInts[3,2]
		ciMinusMode <- confInts[3,1]
	}
	result <- data.frame(nTaxa, nEvents, pct, logOddsS, pValueS, ciPlusS, ciMinusS, logOddsMode, pValueMode, ciPlusMode, ciMinusMode, row.names=NULL)
	class.diff.q <- rbind(class.diff.q, result)
	
	# origination
	logOddsS <- NA
	pValueS <- NA
	ciPlusS <- NA
	ciMinusS <- NA
	logOddsMode <- NA
	pValueMode <- NA
	ciPlusMode <- NA
	ciMinusMode <- NA
	
	nEvents <- sum(temp$first_occurrence)
	pct <- nEvents/nTaxa
	if(nEvents >= minEvents & nTaxa - nEvents >= minEvents) {
		myGlm <- glm(first_occurrence ~ classCount + modeCount, data=temp, family='binomial')
		confInts <- confint(myGlm)
		logOddsS <- summary(myGlm)$coefficients[2,1] 
		pValueS <- summary(myGlm)$coefficients[2,4]
		ciPlusS <- confInts[2,1]
		ciMinusS <- confInts[2,2]
		
		logOddsMode <- summary(myGlm)$coefficients[3,1] 
		pValueMode <- summary(myGlm)$coefficients[3,4]
		ciPlusMode <- confInts[3,2]
		ciMinusMode <- confInts[3,1]
	}
	result <- data.frame(nTaxa, nEvents, pct, logOddsS, pValueS, ciPlusS, ciMinusS, logOddsMode, pValueMode, ciPlusMode, ciMinusMode, row.names=NULL)
	class.diff.p <- rbind(class.diff.p, result)
	
	
	#### without cephalopods
	# extinction
	logOddsS <- NA
	pValueS <- NA
	ciPlusS <- NA
	ciMinusS <- NA
	logOddsMode <- NA
	pValueMode <- NA
	ciPlusMode <- NA
	ciMinusMode <- NA
		
	nEvents <- sum(temp$last_occurrence)
	pct <- nEvents/nTaxa
	if(nEvents >= minEvents & nTaxa - nEvents >= minEvents) {
		myGlm <- glm(last_occurrence ~ classCount + modeCount, data=temp[temp$class != 'Cephalopoda',], family='binomial')
		confInts <- confint(myGlm)
		logOddsS <- summary(myGlm)$coefficients[2,1] 
		pValueS <- summary(myGlm)$coefficients[2,4]
		ciPlusS <- confInts[2,1]
		ciMinusS <- confInts[2,2]
		
		logOddsMode <- summary(myGlm)$coefficients[3,1] 
		pValueMode <- summary(myGlm)$coefficients[3,4]
		ciPlusMode <- confInts[3,2]
		ciMinusMode <- confInts[3,1]
	}
	result <- data.frame(nTaxa, nEvents, pct, logOddsS, pValueS, ciPlusS, ciMinusS, logOddsMode, pValueMode, ciPlusMode, ciMinusMode, row.names=NULL)
	class.diff.q.sansceph <- rbind(class.diff.q.sansceph, result)
	
	# origination
	logOddsS <- NA
	pValueS <- NA
	ciPlusS <- NA
	ciMinusS <- NA
	logOddsMode <- NA
	pValueMode <- NA
	ciPlusMode <- NA
	ciMinusMode <- NA
	
	nEvents <- sum(temp$first_occurrence)
	pct <- nEvents/nTaxa
	if(nEvents >= minEvents & nTaxa - nEvents >= minEvents) {
		myGlm <- glm(first_occurrence ~ classCount + modeCount, data=temp[temp$class != 'Cephalopoda',], family='binomial')
		confInts <- confint(myGlm)
		logOddsS <- summary(myGlm)$coefficients[2,1] 
		pValueS <- summary(myGlm)$coefficients[2,4]
		ciPlusS <- confInts[2,1]
		ciMinusS <- confInts[2,2]
		
		logOddsMode <- summary(myGlm)$coefficients[3,1] 
		pValueMode <- summary(myGlm)$coefficients[3,4]
		ciPlusMode <- confInts[3,2]
		ciMinusMode <- confInts[3,1]
	}
	result <- data.frame(nTaxa, nEvents, pct, logOddsS, pValueS, ciPlusS, ciMinusS, logOddsMode, pValueMode, ciPlusMode, ciMinusMode, row.names=NULL)
	class.diff.p.sansceph <- rbind(class.diff.p.sansceph, result)
	####
	
	# class mode correlation
	temp.pear <- cor.test(temp$classCount, temp$modeCount)
	temp.spear <- cor.test(temp$classCount, temp$modeCount, type="spearman")
	
	cor.s.mode <- rbind(cor.s.mode, data.frame(nTaxa, 's'=temp.pear$estimate, 'sMinus'=temp.pear$conf.int[1], 'sPlus'=temp.pear$conf.int[2], 'rho'=temp.spear$estimate, 'rhoMinus'=temp.spear$conf.int[1], 'rhoPlus'=temp.spear$conf.int[2]))
}
row.names(class.diff.q) <- stages.reduced$interval_name
row.names(class.diff.p) <- stages.reduced$interval_name
row.names(class.diff.q.sansceph) <- stages.reduced$interval_name
row.names(class.diff.p.sansceph) <- stages.reduced$interval_name
row.names(cor.s.mode) <- stages.reduced$interval_name
row.names(cor.ext.gen.mode) <- stages.reduced$interval_name
row.names(cor.orig.gen.mode) <- stages.reduced$interval_name
row.names(cor.ext.gen.mode.sansceph) <- stages.reduced$interval_name
row.names(cor.orig.gen.mode.sansceph) <- stages.reduced$interval_name


## ********************************************************************************
## VOLATILITY (Gilinsky 1994), Fig. 4
## ********************************************************************************

classCounts <- table(master.core.phyla$class)
abundClasses <- names(classCounts[classCounts >= 100 & !is.element(names(classCounts), c("Irregulares","Regulares"))])
abundTaxa <- subset(master.core.phyla, is.element(master.core.phyla$class, abundClasses) & stage <= 477.7)
abundTaxa$class <- factor(as.character(abundTaxa$class), levels=abundClasses) # eliminate low diversity classes as factors
abundTaxa$turnover <- 0 # set value for genera that originate OR go extinct in intervals
abundTaxa$turnover[abundTaxa$first_occurrence == 1 | abundTaxa$last_occurrence == 1] <- 1
volatility <- data.frame(matrix(NA, nrow=nrow(stages.reduced), ncol=length(abundClasses), dimnames=list(stages.reduced$interval_name, abundClasses)))
volatilityGilinsky <- volatility
classDurations <- vector()
for(i in 1:length(abundClasses)) {
	taxa <- subset(abundTaxa, class == abundClasses[i])
	
	timeMax <- max(taxa$stage)
	timeMin <- min(taxa$stage)
	duration <- subset(stages.reduced, age_top >= timeMin & age_top <= timeMax) # total duration of class
	
	taxa$stage <- factor(taxa$stage, levels=duration$age_top) # convert stage number to factor

	
	tn <- duration$age_bottom - duration$age_top

	vntot <- as.numeric(table(taxa$stage))
	vn1 <- as.numeric(vntot - tapply(taxa$last_occurrence, taxa$stage, sum)) # 'number of taxa entering stage i + 1' or taxa that survive
	vn <- as.numeric(vntot - tapply(taxa$first_occurrence, taxa$stage, sum)) # 'number of taxa entering stage 1' or taxa that do not originate
	
	xFLxFtxbL <- as.numeric(tapply(taxa$turnover, taxa$stage, sum)) # Volatility = (XFL + XFt + XbL)/X, numerator is total number of taxa that do not pass through
	vol <- xFLxFtxbL/vntot
	volatility[match(duration$interval_name, rownames(volatility)), match(abundClasses[i], colnames(volatility))] <- vol

	vol <- abs(vn1 - vn) / (vntot * tn)
	volatilityGilinsky[match(duration$interval_name, rownames(volatility)), match(abundClasses[i], colnames(volatility))] <- vol
	
	classDurations[i] <- nrow(duration)
}
classVolatility <- data.frame('class'=levels(abundTaxa$class), 
	'modes'=tapply(abundTaxa$ecospace, abundTaxa$class, function(x){return(length(unique(x)))}),
	'genera'=as.numeric(table(abundTaxa$class)),
	'duration'=stages$age_bottom[match(tapply(abundTaxa$stage, abundTaxa$class, max), stages$age_top)] - tapply(abundTaxa$stage, abundTaxa$class, min),
	'volatility'=apply(volatility, 2, sum, na.rm=TRUE)/classDurations,
	'volatilityGilinsky'=apply(volatilityGilinsky, 2, sum, na.rm=TRUE)/classDurations
)

vol.cor <- data.frame()
vol.cor.gil <- vol.cor
for(i in 1:nrow(stages.reduced)) {
	# get extant taxa
	temp <- subset(abundTaxa, stage==stages.reduced$age_top[i])

	genusCounts <- as.numeric(table(temp$class))
	modeCounts <- as.numeric(tapply(temp$ecospace, temp$class, function(x){return(length(unique(x)))}))

	nTaxa <- nrow(temp)
	
	classInclude <- genusCounts > 0
	# extinction
	if(var(as.numeric(volatility[i,classInclude])) > 0) {
		y <- as.numeric(volatility[i,classInclude])
		xMode <- scale(modeCounts[classInclude])[,1]
		xGenera <- scale(log10(genusCounts[classInclude]))[,1]
		myGlm <- lm(y ~ xMode + xGenera, na.action='na.omit')
		confInts <- confint(myGlm)
		stdErr <- summary(myGlm)$coefficients[-1,]
		result <- data.frame(nTaxa, 'modeCoef'=myGlm$coefficients[match('xMode',names(myGlm$coefficients))], 'modeMinus'=confInts[match('xMode', rownames(confInts)),1], 'modePlus'=confInts[match('xMode', rownames(confInts)),2], 'generaCoef'=myGlm$coefficients[match('xGenera',names(myGlm$coefficients))], 'generaMinus'=confInts[match('xGenera', rownames(confInts)),1], 'generaPlus'=confInts[match('xGenera', rownames(confInts)),2])	

		y <- as.numeric(volatilityGilinsky[i,classInclude])
		myGlm <- lm(y ~ xMode + xGenera, na.action='na.omit')
		confInts <- confint(myGlm)
		result.gil <- data.frame(nTaxa, 'modeCoef'=myGlm$coefficients[match('xMode',names(myGlm$coefficients))], 'modeMinus'=confInts[match('xMode', rownames(confInts)),1], 'modePlus'=confInts[match('xMode', rownames(confInts)),2], 'generaCoef'=myGlm$coefficients[match('xGenera',names(myGlm$coefficients))], 'generaMinus'=confInts[match('xGenera', rownames(confInts)),1], 'generaPlus'=confInts[match('xGenera', rownames(confInts)),2])	

	} else {
		result <- data.frame(nTaxa, 'modeCoef'=NA, 'modeMinus'=NA, 'modePlus'=NA, 'generaCoef'=NA, 'generaMinus'=NA, 'generaPlus'=NA)
		result.gil <- data.frame(nTaxa, 'modeCoef'=NA, 'modeMinus'=NA, 'modePlus'=NA, 'generaCoef'=NA, 'generaMinus'=NA, 'generaPlus'=NA)
	}
	vol.cor <- rbind(vol.cor, result)
	vol.cor.gil <- rbind(vol.cor.gil, result.gil)
}
row.names(vol.cor) <- stages.reduced$interval_name
row.names(vol.cor.gil) <- stages.reduced$interval_name

vol.mean.cor <- data.frame(matrix(NA, nrow=2, ncol=3, dimnames=list(c('mode','genera'),c('mean','n','ci'))))
vol.mean.cor$mean[1] <- mean(vol.cor$modeCoef, na.rm=T)
vol.mean.cor$n[1] <- length(vol.cor$modeCoef[!is.na(vol.cor$modeCoef)])
vol.mean.cor$ci[1] <- mean(vol.cor$modeCoef, na.rm=T) * sd(vol.cor$modeCoef, na.rm=T) / sqrt(length(vol.cor$modeCoef[!is.na(vol.cor$modeCoef)]))

vol.mean.cor$mean[2] <- mean(vol.cor$generaCoef, na.rm=T)
vol.mean.cor$n[2] <- length(vol.cor$generaCoef[!is.na(vol.cor$generaCoef)])
vol.mean.cor$ci[2] <- mean(vol.cor$generaCoef, na.rm=T) * sd(vol.cor$generaCoef, na.rm=T) / sqrt(length(vol.cor$generaCoef[!is.na(vol.cor$generaCoef)]))


## ********************************************************************************
## FUNCTIONS
## ********************************************************************************

source("/Users/knope/Desktop/timescalePlotFunction.r");



## fix few missing CI
for(i in 1:nrow(class.diff.p)) {
	if(is.na(class.diff.p$ciPlusMode[i]) & !is.na(class.diff.p$ciMinusMode[i])) {
		class.diff.p$ciPlusMode[i] <- class.diff.p$logOddsMode[i] + (class.diff.p$ciMinusMode[i] - class.diff.p$logOddsMode[i]) * -1
	} else if(!is.na(class.diff.p$ciPlusMode[i]) & is.na(class.diff.p$ciMinusMode[i])) {
		class.diff.p$ciMinusMode[i] <- class.diff.p$logOddsMode[i] - (class.diff.p$ciPlusMode[i] - class.diff.p$logOddsMode[i])
	}
	
	if(is.na(class.diff.p$ciPlusS[i]) & !is.na(class.diff.p$ciMinusS[i])) {
		class.diff.p$ciPlusS[i] <- class.diff.p$logOddsS[i] + (class.diff.p$ciMinusS[i] - class.diff.p$logOddsS[i]) * -1
	} else if(!is.na(class.diff.p$ciPlusS[i]) & is.na(class.diff.p$ciMinusS[i])) {
		class.diff.p$ciMinusS[i] <- class.diff.p$logOddsS[i] - (class.diff.p$ciPlusS[i] - class.diff.p$logOddsS[i])
	}
}


## ********************************************************************************
## PLOT DATA
## ********************************************************************************

library(berryFunctions)

# extinction and origination intensity
time.plot.mult(nrow=2, ncol=1, pdf.name="FigS5_CorrelationOrigExt.pdf", plot.width=10, plot.height=8.5, bottom.mar=4.25, left.mar=4, top.mar=1.5, timescale='post-cambrian', mgpX=c(2.2, 0.75, 0), mgpY=c(3, 0.75, 0))
par(las=1, pch=16)
	
	#extinction
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="Linear regression coefficient", ylim=range(c(cor.ext.gen.mode$modeCoef, cor.ext.gen.mode$generaCoef),na.rm=TRUE)+c(-0.1,0.1))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, cor.ext.gen.mode$generaCoef, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, cor.ext.gen.mode$generaPlus, stages.reduced$age_mid, cor.ext.gen.mode$generaMinus, col="#ff7285")
	points(stages.reduced$age_mid, cor.ext.gen.mode$modeCoef, col="blue")
	segments(stages.reduced$age_mid, cor.ext.gen.mode$modePlus, stages.reduced$age_mid, cor.ext.gen.mode$modeMinus, col="blue")
	legend("topleft", legend=c("genus richness","ecological diversity"), pch=c(18,16), col=c("#ff7285","blue"), pt.cex=c(1.25,1), bty="n")
	axis(side=1, at=seq(0,400,100), labels=FALSE)
	text(x=485.4/2, y=max(c(cor.ext.gen.mode$modeCoef, cor.ext.gen.mode$generaCoef),na.rm=TRUE)+0.12, pos=1, labels="Extinction", cex=1.25, font=2)
	
	# origination 
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="Linear regression coefficient", ylim=range(c(cor.orig.gen.mode$modeCoef, cor.orig.gen.mode$generaCoef),na.rm=TRUE)+c(-0.1,0.1))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, cor.orig.gen.mode$generaCoef, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, cor.orig.gen.mode$generaPlus, stages.reduced$age_mid, cor.orig.gen.mode$generaMinus, col="#ff7285")
	points(stages.reduced$age_mid, cor.orig.gen.mode$modeCoef, col="blue")
	segments(stages.reduced$age_mid, cor.orig.gen.mode$modePlus, stages.reduced$age_mid, cor.orig.gen.mode$modeMinus, col="blue")
	text(x=485.4/2, y=max(c(cor.orig.gen.mode$modeCoef, cor.orig.gen.mode$generaCoef),na.rm=TRUE)+0.12, pos=1, labels="Origination", cex=1.25, font=2)
dev.off()


# extinction and origination intensity WITHOUT CEPHALOPODS
time.plot.mult(nrow=2, ncol=1, pdf.name="FigS6_CorrelationOrigExtSansCephalopods.pdf", plot.width=10, plot.height=8.5, bottom.mar=4.25, left.mar=4, top.mar=1.5, timescale='post-cambrian', mgpX=c(2.2, 0.75, 0), mgpY=c(3, 0.75, 0))
par(las=1, pch=16)
	
	#extinction
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="Linear regression coefficient", ylim=range(c(cor.ext.gen.mode.sansceph$modeCoef, cor.ext.gen.mode.sansceph$generaCoef),na.rm=TRUE)+c(-0.1,0.1))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, cor.ext.gen.mode.sansceph$generaCoef, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, cor.ext.gen.mode.sansceph$generaPlus, stages.reduced$age_mid, cor.ext.gen.mode.sansceph$generaMinus, col="#ff7285")
	points(stages.reduced$age_mid, cor.ext.gen.mode.sansceph$modeCoef, col="blue")
	segments(stages.reduced$age_mid, cor.ext.gen.mode.sansceph$modePlus, stages.reduced$age_mid, cor.ext.gen.mode.sansceph$modeMinus, col="blue")
	legend("topleft", legend=c("genus richness","ecological diversity"), pch=c(18,16), col=c("#ff7285","blue"), pt.cex=c(1.25,1), bty="n")
	axis(side=1, at=seq(0,400,100), labels=FALSE)
	text(x=485.4/2, y=max(c(cor.ext.gen.mode.sansceph$modeCoef, cor.ext.gen.mode.sansceph$generaCoef),na.rm=TRUE)+0.12, pos=1, labels="Extinction", cex=1.25, font=2)
	
	# origination 
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="Linear regression coefficient", ylim=range(c(cor.orig.gen.mode.sansceph$modeCoef, cor.orig.gen.mode.sansceph$generaCoef),na.rm=TRUE)+c(-0.1,0.1))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, cor.orig.gen.mode.sansceph$generaCoef, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, cor.orig.gen.mode.sansceph$generaPlus, stages.reduced$age_mid, cor.orig.gen.mode.sansceph$generaMinus, col="#ff7285")
	points(stages.reduced$age_mid, cor.orig.gen.mode.sansceph$modeCoef, col="blue")
	segments(stages.reduced$age_mid, cor.orig.gen.mode.sansceph$modePlus, stages.reduced$age_mid, cor.orig.gen.mode.sansceph$modeMinus, col="blue")
	text(x=485.4/2, y=max(c(cor.orig.gen.mode.sansceph$modeCoef, cor.orig.gen.mode.sansceph$generaCoef),na.rm=TRUE)+0.12, pos=1, labels="Origination", cex=1.25, font=2)
dev.off()

	
# extinction and origination selectivity
time.plot.mult(nrow=2, ncol=1, pdf.name="FigS7_LogOddsOrigExt.pdf", plot.width=10, plot.height=8.5, bottom.mar=4.25, left.mar=4, top.mar=1.5, timescale='post-cambrian', mgpX=c(2.2, 0.75, 0), mgpY=c(3, 0.75, 0))
par(las=1, pch=16)
	
	#extinction
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="log-odds of extinction", ylim=range(c(class.diff.q$logOddsS, class.diff.q$logOddsMode),na.rm=TRUE)+c(-0.1,0.1))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, class.diff.q$logOddsS, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, class.diff.q$ciPlusS, stages.reduced$age_mid, class.diff.q$ciMinusS, col="#ff7285")
	points(stages.reduced$age_mid, class.diff.q$logOddsMode, col="blue")
	segments(stages.reduced$age_mid, class.diff.q$ciPlusMode, stages.reduced$age_mid, class.diff.q$ciMinusMode, col="blue")
	legend("topleft", legend=c("genus richness","ecological diversity"), pch=c(18,16), col=c("#ff7285","blue"), pt.cex=c(1.25,1), bty="n")
	axis(side=1, at=seq(0,400,100), labels=FALSE)
	text(x=485.4/2, y=max(c(class.diff.q$logOddsS, class.diff.q$logOddsMode),na.rm=TRUE)+0.40, pos=1, labels="Extinction", cex=1.25, font=2)
	
	# origination 
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="log-odds of origination", ylim=c(min(range(class.diff.p$logOddsS, class.diff.p$logOddsMode,na.rm=TRUE)), 5)+c(-0.1,0.2))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, class.diff.p$logOddsS, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, class.diff.p$ciPlusS, stages.reduced$age_mid, class.diff.p$ciMinusS, col="#ff7285")
	points(stages.reduced$age_mid, class.diff.p$logOddsMode, col="blue")
	segments(stages.reduced$age_mid, class.diff.p$ciPlusMode, stages.reduced$age_mid, class.diff.p$ciMinusMode, col="blue")
	text(x=485.4/2, y=5.5, pos=1, labels="Origination", cex=1.25, font=2)
dev.off()


# volatility correlations with modes and genera
myBreaks <- seq(-0.2, 0.3, 0.01)
nInts <- length(vol.cor$generaCoef[!is.na(vol.cor$generaCoef)])

time.plot.mult(nrow=2, ncol=2, pdf.name="Fig4_CorrelationVolatilityModesGenera.pdf", plot.width=c(10,3.25), plot.height=8.75, bottom.mar=4.65, left.mar=4.25, top.mar=1.5, right.mar=0.45, 
  timescale='post-cambrian', time.height=1.35, all.cols=1, mgpX=c(2.2, 0.75, 0), mgpY=c(3, 0.75, 0))
  
  par(las=1, pch=16)
	
	# Genera Time Series
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="Linear regression coefficient", ylim=c(-0.3, 0.4))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, vol.cor$generaCoef, col="#ff7285", pch=18, cex=1.25)
	segments(stages.reduced$age_mid, vol.cor$generaPlus, stages.reduced$age_mid, vol.cor$generaMinus, col="#ff7285")
	axis(side=1, at=seq(0,400,100), labels=FALSE)
	text(x=485.4/2, y=0.425, pos=1, labels="Genus Richness", cex=1.25, font=2, xpd=TRUE)
	
	# Genera histogram
	par(mar=c(0, 0, 1.5, 0.65), xaxs="r", lwd=0.5)
	temp <- horizHist(vol.cor$generaCoef, breaks=myBreaks, xlim=c(0,12), ylim=c(-0.3, 0.4), xlab="", main="", col="#ff7285", labelat=NA, axes=FALSE, xaxs="r")
	tempStats <- c(mean(vol.cor$generaCoef, na.rm=TRUE), sd(vol.cor$generaCoef, na.rm=TRUE))

	par(new=TRUE, lwd=1)
	plot(1:10, type="n", axes=FALSE, xlab="",ylab="", xlim=c(0,12), ylim=c(-0.3, 0.4))
	axis(side=1, at=seq(0,12,3), labels=NA, tick=TRUE)
	abline(h=0, lty=2, lwd=1)
	points(11.5, tempStats[1], col="#ff7285", pch=16, cex=0.5)
	arrows(11.5, tempStats[1]-1.96*tempStats[2]/sqrt(nInts), 11.5, tempStats[1]+1.96*tempStats[2]/sqrt(nInts), code=3, angle=90, lwd=1, length=0.075, col="#ff7285")
  abline(v=0, lwd=1)
  
	# Modes Time Series
	par(mar=c(0, 4.25, 1.5, 0.45), xaxs="i")
	plot(1:10, type="n", xaxt="n", xlab="", xlim=c(485.4,0), ylab="Linear regression coefficient", ylim=c(-0.3, 0.4))
	abline(h=0, lty=2)
	points(stages.reduced$age_mid, vol.cor$modeCoef, col="blue")
	segments(stages.reduced$age_mid, vol.cor$modePlus, stages.reduced$age_mid, vol.cor$modeMinus, col="blue")
	text(x=485.4/2, y=0.425, pos=1, labels="Ecological Diversity", cex=1.23, font=2, xpd=TRUE)	

	# Modes histogram
	par(mar=c(0, 0, 1.5, 0.65), xaxs="r", lwd=0.5)
	temp <- horizHist(vol.cor$modeCoef, breaks=myBreaks, xlim=c(0,12), ylim=c(-0.3, 0.4), xlab="", main="", col="blue", labelat=NA, axes=FALSE, xaxs="r")
	tempStats <- c(mean(vol.cor$modeCoef, na.rm=TRUE), sd(vol.cor$modeCoef, na.rm=TRUE))
	
	par(new=TRUE, lwd=1)
	plot(1:10, type="n", axes=FALSE, xlab="",ylab="", xlim=c(0,12), ylim=c(-0.3, 0.4))
	axis(side=1, at=seq(0,12,3), line=0.97, cex.axis=1.2152, tcl=-0.475)
	mtext("Number of stages", side=1, cex=1.5, line=4)
	abline(h=0, lty=2, lwd=1)
	points(11.5, tempStats[1], col="blue", pch=16, cex=0.5)
	arrows(11.5, tempStats[1]-1.96*tempStats[2]/sqrt(nInts), 11.5, tempStats[1]+1.96*tempStats[2]/sqrt(nInts), code=3, angle=90, lwd=1, length=0.075, col="blue")
	par(xpd=NA)
	lines(c(0,0), c(-0.38,0.428))
dev.off()



##############################################################################
#  Assess clade age vs. taxonomic richness and ecological diversity: Fig. S3. 
##############################################################################
source("timescalePlotFunction.r")

library(vegan)
#install.packages("vegan")

master<-master.time.frame # make shorter variable name


fossil.only <- droplevels(subset(master, stage > 2.588 & class != "")) #exclude Pleistocene and Holocene genera and genera that have not been assigned to a Class
class.max <- tapply(fossil.only$stage, fossil.only$class, max, na.rm=T) # get maximum age for each class

class.counts <- table(fossil.only$class) # count genera in each class
# all genera, including Pleistocene & Holocene


modes.by.class <- as.data.frame.matrix(table(fossil.only$ecospace, fossil.only$class))
modes.by.class[modes.by.class >= 1] <- 1 # convert counts to presence absence (of mode)
class.modes <- apply(modes.by.class, 2, sum) # sum modes 

# plot
quartz(height=7, width=14)
par(mfrow=c(1,2), pch=16)
# genus richness
plot(as.numeric(class.max), log10(as.numeric(class.counts)), main="Genus Richness",
	xlab="Age of class origination", 
	ylab=expression(paste("log"[10],"(Total number of genera)")), xlim=c(541,0), las=1, col="dark green")
# number of modes
plot(as.numeric(class.max), log10(class.modes), main="Ecological Diversity", 
	xlab="Age of class origination", 
	ylab=expression(paste("log"[10],"(Total number of ecological modes)")), xlim=c(541,0), las=1, col="dark green" )
	
## correlation tests
cor.test(class.counts, class.max, method='spearman')
cor.test(class.modes, class.max, method='spearman')



##############################################################################
# Taxonomic evenness (number of genera) in ecological modes of life across time calculated as Shannonâ€™s Diversity Index and the Probability of Interspecific Encounter (inverse Simpson Diversity Index normalized for functional richness): Fig. S10
##############################################################################

library(vegan)
ecoSTAGE<-xtabs(~master$stage+ master$ecospace) # make time by ecospace

# plot shannon's diversity and PIE in two frame plot with spiffy Phanerozoic timescale as x-axis
time.plot.mult(nrow=1, ncol=2, time.height=1, left.mar = 4.2, mgpX=c(2.2, 0.75, 0), mgpY=c(3, 0.75, 0), timescale='post-cambrian')

div.ind<-diversity(ecoSTAGE,index="shannon",MARGIN=1) # calculate Shannon's Diversity

# plot Shannon
plot(as.numeric(names(div.ind)), div.ind, xlim=c(485.5,0), lwd=3, xaxt="n", xlab="", ylab="Shannon's Diversity Index", type="l", las=1)
polygon(x=c(0,as.numeric(names(div.ind)),444) , y=c(0,div.ind,0), col="blue" )

abline(v=66,col=1,lty=5)
abline(v=200,col=1,lty=3)
abline(v=252.5,col=1,lty=5)
abline(v=359,col=1,lty=3)
#abline(v=445,col=1,lty=3)
text(465,2.9,"(A)", font=2, cex=1.5)

#Probability of Interspecific Encounter (PIE): (S/(S-1))*(1-D), where D=Simpson's Diversity (vegan returns 1-D), S=richness
S <- specnumber(ecoSTAGE)
div.ind <- (S/(S-1)) * diversity(ecoSTAGE,index="simpson",MARGIN=1)

plot(as.numeric(names(div.ind)), div.ind, xlim=c(485.5,0), lwd=3, xaxt="n", xlab="", ylab="Probability of Interspecific Encounter (PIE)", type="l", las=1)
polygon(x=c(0,as.numeric(names(div.ind)),444) , y=c(0,div.ind,0), col="blue" )

abline(v=66,col=1,lty=5)
abline(v=200,col=1,lty=3)
abline(v=252.5,col=1,lty=5)
abline(v=359,col=1,lty=3)
#abline(v=445,col=1,lty=3)
text(465,0.934,"(B)", font=2, cex=1.5)


##############################################################################
General linear model for ecological covariates of extinction: Table S3 
##############################################################################

library(plyr)
library(lme4)

## Read stage data
stages<-read.csv("/timescale_final.csv")

## Pull in modern data
dd.modern.raw<-read.csv("/Users/knope/Desktop/genera_modes_modern_final.csv")

## Establish periods for some downstream applications
period<-character(length(stages$age_top))

## Establish core phyla, and remove all others:
core.phyla<-c('Hemichordata','Porifera','Cnidaria','Bryozoa','Gastropoda','Cephalopoda',
				'Bivalvia','Mollusca','Echinodermata','Brachiopoda','Arthropoda',
				'Chordata')
master.core.phyla<-master[master$phylum %in% core.phyla,]

# Drop levels, and create a bin category for use with aggregate function to 
# collapse df across stages
master.core.phyla$phylum<-factor(master.core.phyla$phylum)
master.core.phyla$class<-factor(master.core.phyla$class)

classes <- sort(unique(master.core.phyla$class)) #get a list of unique classes for later use
stages <- sort(unique(master.core.phyla$stage)) #get a list of unique stages for later use


## ************************************************************************

##find when genera go extinct and add this to the master file##
ex.data <- ddply(master.core.phyla, c("taxon_genus"), function(df)c(min(df$stage), max(df$stage), max(df$stage)-min(df$stage)))
names(ex.data) <- c("taxon_genus", "ex.stage", "or.stage", "surv.time")
genus.data <- merge(master.core.phyla, ex.data, by="taxon_genus")
genus.data$ex <- 0
genus.data$ex[genus.data$ex.stage==genus.data$stage] <- 1
genus.data$or <- 0
genus.data$or[genus.data$or.stage==genus.data$stage] <- 1

##get counts of genera and modes by class through the Phanerozoic##
class.data <- ddply(master.core.phyla, c("class", "stage"), function(df)c(length(df$taxon_genus), length(unique(df$ecospace))))
names(class.data) <- c("class", "stage", "genera", "modes")

ngenera<-length(unique(genus.data$taxon_genus))


## ************************************************************************
## Run an analysis on all stages at once in a GLMM framework
## ************************************************************************
genus.data$tiering<-factor(genus.data$tiering)
genus.data$motility_three<-factor(genus.data$motility_three)
genus.data$buffering_three<-factor(genus.data$buffering_three)

# Compile master datasheet with stage specific ecospace counts

tmp<-genus.data[genus.data$stage== unique(genus.data$stage)[order(unique(genus.data$stage))][2],]


# Calculate modes of class during a given time slice
tmp.classes<-unique(tmp$class)
tmp.ecos.class<-sapply(tmp.classes, function(x) length(unique(tmp[tmp$class==x,'ecospace'])))

names(tmp.ecos.class)<-tmp.classes	
tmp$n.ecos.class<-tmp.ecos.class[match(tmp$class, names(tmp.ecos.class))]
tmp$n.ecos.class.ls <-scale(log(tmp$n.ecos.class))
old<-tmp

nstage<-length(unique(genus.data$stage))

for (i in 3:nstage) {
		tmp<-genus.data[genus.data$stage== unique(genus.data$stage)[order(unique(genus.data$stage))][i],]
	
	
	# Calculate modes of class during a given time slice
	tmp.classes<-unique(tmp$class)
	tmp.ecos.class<-sapply(tmp.classes, function(x) 
		length(unique(tmp[tmp$class==x,'ecospace'])))
	
	names(tmp.ecos.class)<-tmp.classes	
	tmp$n.ecos.class<-tmp.ecos.class[match(tmp$class, names(tmp.ecos.class))]
	tmp$n.ecos.class.ls <-scale(log(tmp$n.ecos.class))

	old<-rbind(tmp, old)

}

# Eliminate unknown classes from dataset
old1<-old[old$class!="",]

# Note that a number of buffer values are NAs
old2<-old1[!is.na(old1$buffering_binary),]

## ********************************************************************************

## Final model set

a.full<-glmer(ex ~ n.ecos.class.ls + buffering_binary*motility_three*feeding_binary+tiering_binary + (1|stage) + (1|class) + (1|taxon_genus), data=old2, family='binomial', nAGQ=0)
summary(a.full)



## ********************************************************************************
##The relationship between origination and extinction rates across geologic time Fig S8 and table S4
## ********************************************************************************

#install.packages("viridis")
library(viridis)
cols<-viridis(5)
#modern is "dodgerblue"

mtf <- master.time.frame

#call plyr library for later use
library(plyr)

#identify extinction and origination by stage and then code as binary variables
orex <- ddply(mtf, c("taxon_genus"), function(df)c(max(df$stage), min(df$stage)))
names(orex) <- c("taxon_genus", "or.stage", "ex.stage")

mtf.1 <- merge(mtf, orex, by="taxon_genus")
mtf.1$or <- 0
mtf.1$or[mtf.1$stage==mtf.1$or.stage] <- 1
mtf.1$ex <- 0
mtf.1$ex[mtf.1$stage==mtf.1$ex.stage] <- 1

mtf.1$duration <- mtf.1$or.stage - mtf.1$ex.stage

mtf.1 <- subset(mtf.1, duration<150 & class!="" & (ex==0 | or==0)) #remove, long-ranging genera and singletons

#get stage-by-stage extinction, origination rates, counts of genera and modes
  mod.i <- ddply(mtf.1, c("class", "stage"), function(df)
    c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])), #extinction rate
      -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])), #origination rate
      length(unique(df$taxon_genus)),
      length(unique(df$ecospace))))
  names(mod.i) <- c("class", "stage", "logex", "logor", "genera", "modes")
  mod.i <- subset(mod.i, logex<100 & logor<100) #remove infinite rates and NA values
  
#calculate mean rate by class for stages >= to stage of interest
  out <- data.frame(matrix(nrow=1, ncol=4))
  names(out) <- c("class", "mean.ex", "mean.or", "stage")
  #get list of unique stages
  stages <- sort(unique(mod.i$stage))
for (i in 1: length(stages)) {
  data.i <- subset(mod.i, stage>=stages[i])
  mod.ii <- ddply(data.i, c("class"), function(df) #calculate average extinction and origination rates
    c(mean(df$logex), 
      mean(df$logor)))
    names(mod.ii) <- c("class", "mean.ex", "mean.or")
  mod.ii$stage <- stages[i]
  out <- rbind(out, mod.ii)
}
  
out.2 <- na.omit(merge(out, mod.i, by=c("stage", "class")))

### plot extinction versus origination ###
pdf("extinction vs origination.pdf", h=5, w=10)
par(mfrow=c(2,5), mar=c(4,4,1,1))
plot.stages <- c(419.2, 358.9, 298.9, 252.2, 201.3, 145, 66, 61.6, 23.03, 2.588)
names <- c("Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", 
           "Jurassic", "Cretaceous", "Danian", "Paleogene", "Pleistocene")
min.gen <- 2
out.5 <- data.frame(matrix(ncol=4, nrow=length(plot.stages)))
names(out.5) <- c("Period", "Coefficient", "R.squared", "p.value")
for (i in 1:length(plot.stages)) {
plot(out.2$mean.ex[out.2$stage==plot.stages[i] & out.2$genera>=min.gen]
     ~out.2$mean.or[out.2$stage==plot.stages[i] & out.2$genera>=min.gen],
     las=1.5,
     xlim=c(0,1.0),
     ylim=c(0,1.0),
     xlab="origination rate", 
     ylab="extinction rate")
segments(0,0,1.5,1.5)
mtext(side=3, names[i], line=-1, cex=0.75)
temp <- lm(out.2$mean.ex[out.2$stage==plot.stages[i] & out.2$genera>=min.gen]
           ~out.2$mean.or[out.2$stage==plot.stages[i] & out.2$genera>=min.gen])
t1 <- summary(temp)
out.5[i,1] <- names[i]
out.5[i,2] <- t1$coefficients[2]
out.5[i,3] <- t1$r.squared
out.5[i,4] <- t1$coefficients[2,4]
}

dev.off()

write.table(out.5, file="exor.corr.coefficients.txt", sep="\\t")

 
## ************************************************************************************************************
##The effect of number of ecological modes and number of genera on diversification rate shifts Fig. S8
## ************************************************************************************************************


library(scales)

# call plyr library for later use
library(plyr)

# identify extinction and origination by stage and then code as binary variables
orex <- ddply(mtf, c("taxon_genus"), function(df)c(max(df$stage), min(df$stage)))
names(orex) <- c("taxon_genus", "or.stage", "ex.stage")

mtf.1 <- merge(mtf, orex, by="taxon_genus")
mtf.1$or <- 0
mtf.1$or[mtf.1$stage==mtf.1$or.stage] <- 1
mtf.1$ex <- 0
mtf.1$ex[mtf.1$stage==mtf.1$ex.stage] <- 1

mtf.1$duration <- mtf.1$or.stage - mtf.1$ex.stage

mtf.1 <- subset(mtf.1, duration<150 & class!="" & (ex==0 | or==0)) # remove, long-ranging genera and singletons

# get stage-by-stage extinction, origination rates, counts of genera and modes
  mod.i <- ddply(mtf.1, c("class", "stage"), function(df)
    c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])), #extinction rate
      -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])), #origination rate
      length(unique(df$taxon_genus)),
      length(unique(df$ecospace))))
  names(mod.i) <- c("class", "stage", "logex", "logor", "genera", "modes")
  mod.i <- subset(mod.i, logex<100 & logor<100) #remove infinite rates and NA values
  stages <- sort(unique(mod.i$stage))



## Create div rates
mod.i$div<-mod.i$logor - mod.i$logex

# Adjust variables
mod.i$lmodes<-log(mod.i$modes)
mod.i$lgen<-log(mod.i$genera)
mod.i$fstage<-factor(mod.i$stage)

## Create temporary dataframe, to and run analysis on first stage of data
tmp<-mod.i[mod.i$fstage==levels(mod.i$fstage)[1],]

mod.full<-lm(div ~ lmodes + lgen, data=tmp)
sds<-1.96*summary(mod.full)[[4]][2:3,2]
ps<-summary(mod.full)[[4]][2:3,4]

mod.modes<-lm(div ~ lmodes, data=tmp)

mod.gen<-lm(div ~ lgen, data=tmp)

xx<-data.frame(lmod.full=coef(mod.full)["lmodes"], lgen.full=coef(mod.full)["lgen"],
	lgen.gen=coef(mod.gen)["lgen"], lmod.modes=coef(mod.modes)["lmodes"],
	se.lmodes=sds[1], se.lgen=sds[2],
	p.lmodes=ps[1], p.lgen=ps[2])
rownames(xx)<-levels(mod.i$fstage)[1]

# Loop through all stages
for (i in 2:length(levels(mod.i$fstage))) {
	tmp<-mod.i[mod.i$fstage==levels(mod.i$fstage)[i],]
	
	mod.full<-lm(div ~ lmodes + lgen, data=tmp)
	sds<-1.96*summary(mod.full)[[4]][2:3,2]
	ps<-summary(mod.full)[[4]][2:3,4]
	mod.modes<-lm(div ~ lmodes, data=tmp)
	
	mod.gen<-lm(div ~ lgen, data=tmp)
	
	xxx<-data.frame(lmod.full=coef(mod.full)["lmodes"], lgen.full=coef(mod.full)["lgen"],
		lgen.gen=coef(mod.gen)["lgen"], lmod.modes=coef(mod.modes)["lmodes"], 
		se.lmodes=sds[1], se.lgen=sds[2],
		p.lmodes=ps[1], p.lgen=ps[2])
	rownames(xxx)<-levels(mod.i$fstage)[i]
	
	xx<-rbind(xx,xxx)
}


# Mixed effect model for overall slope
mod.time<-lmer(div ~ stage*lgen*lmodes + (1|fstage), data=mod.i)
summary(mod.time)
drop1(mod.time, test='Chi')

mod.time1<-lmer(div ~ stage*lgen + stage*lmodes + lmodes*lgen + (1|fstage), data=mod.i)
drop1(mod.time1, test='Chi')

mod.time2<-lmer(div ~ stage*lgen + stage*lmodes + (1|fstage), data=mod.i)
drop1(mod.time2, test='Chi')

mod.time.final<-mod.time2
summary(mod.time.final)



pdf(file="diversification-eachstage-controlled-alltime-allgenera.pdf", width=8, height=4)

	par(mfrow=c(1,2), mar=c(5,6,2,1))
	
	plot(as.numeric(rownames(xx)), xx$lmod.full, type='l', xlab='Time (Ma)', ylab='Effect on\\nDiversification Rate', main='N. modes', cex.lab=1.5, xlim=c(530,0), las=1); abline(h=0, lty=2)
	points(0:500, fixef(mod.time.final)["lmodes"] + 0:500*fixef(mod.time.final)["stage:lmodes"], type='l', lwd=5)
	
	
	plot(as.numeric(rownames(xx)), xx$lgen.full, type='l', xlab='Time (Ma)', ylab='Effect on\\nDiversification Rate', main='N. genera', cex.lab=1.5, xlim=c(530,0), las=1); abline(h=0, lty=2)
	points(0:500, fixef(mod.time.final)["lgen"] + 0:500*fixef(mod.time.final)["stage:lgen"], col='black', type='l', lwd=5)

dev.off()


## ************************************************************************************************************
## Fig S12 and Table S4
## ************************************************************************************************************

#install.packages("viridis")
library(viridis)
cols<-viridis(5)
#modern is "dodgerblue"

mtf <- master.time.frame

#call plyr library for later use
library(plyr)

#identify extinction and origination by stage and then code as binary variables
orex <- ddply(mtf, c("taxon_genus"), function(df)c(max(df$stage), min(df$stage)))
names(orex) <- c("taxon_genus", "or.stage", "ex.stage")

mtf.1 <- merge(mtf, orex, by="taxon_genus")
mtf.1$or <- 0
mtf.1$or[mtf.1$stage==mtf.1$or.stage] <- 1
mtf.1$ex <- 0
mtf.1$ex[mtf.1$stage==mtf.1$ex.stage] <- 1

mtf.1$duration <- mtf.1$or.stage - mtf.1$ex.stage

mtf.1 <- subset(mtf.1, duration<150 & class!="" & (ex==0 | or==0)) #remove, long-ranging genera and singletons

#get stage-by-stage extinction, origination rates, counts of genera and modes
  mod.i <- ddply(mtf.1, c("class", "stage"), function(df)
    c(-log(length(df$ex[df$ex==0 & df$or==0])/length(df$or[df$or==0])), #extinction rate
      -log(length(df$or[df$ex==0 & df$or==0])/length(df$ex[df$ex==0])), #origination rate
      length(unique(df$taxon_genus)),
      length(unique(df$ecospace))))
  names(mod.i) <- c("class", "stage", "logex", "logor", "genera", "modes")
  mod.i <- subset(mod.i, logex<100 & logor<100) #remove infinite rates and NA values
  
#calculate mean rate by class for stages >= to stage of interest
  out <- data.frame(matrix(nrow=1, ncol=4))
  names(out) <- c("class", "mean.ex", "mean.or", "stage")
  #get list of unique stages
  stages <- sort(unique(mod.i$stage))
for (i in 1: length(stages)) {
  data.i <- subset(mod.i, stage>=stages[i])
  mod.ii <- ddply(data.i, c("class"), function(df) #calculate average extinction and origination rates
    c(mean(df$logex), 
      mean(df$logor)))
    names(mod.ii) <- c("class", "mean.ex", "mean.or")
  mod.ii$stage <- stages[i]
  out <- rbind(out, mod.ii)
}
  
out.2 <- na.omit(merge(out, mod.i, by=c("stage", "class")))

### plot extinction versus origination ###
pdf("extinction vs origination.pdf", h=5, w=10)
par(mfrow=c(2,5), mar=c(4,4,1,1))
plot.stages <- c(419.2, 358.9, 298.9, 252.2, 201.3, 145, 66, 61.6, 23.03, 2.588)
names <- c("Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", 
           "Jurassic", "Cretaceous", "Danian", "Paleogene", "Pleistocene")
min.gen <- 2
out.5 <- data.frame(matrix(ncol=4, nrow=length(plot.stages)))
names(out.5) <- c("Period", "Coefficient", "R.squared", "p.value")
for (i in 1:length(plot.stages)) {
plot(out.2$mean.ex[out.2$stage==plot.stages[i] & out.2$genera>=min.gen]
     ~out.2$mean.or[out.2$stage==plot.stages[i] & out.2$genera>=min.gen],
     las=1.5,
     xlim=c(0,1.0),
     ylim=c(0,1.0),
     xlab="origination rate", 
     ylab="extinction rate")
segments(0,0,1.5,1.5)
mtext(side=3, names[i], line=-1, cex=0.75)
temp <- lm(out.2$mean.ex[out.2$stage==plot.stages[i] & out.2$genera>=min.gen]
           ~out.2$mean.or[out.2$stage==plot.stages[i] & out.2$genera>=min.gen])
t1 <- summary(temp)
out.5[i,1] <- names[i]
out.5[i,2] <- t1$coefficients[2]
out.5[i,3] <- t1$r.squared
out.5[i,4] <- t1$coefficients[2,4]
}

dev.off()

write.table(out.5, file="exor.corr.coefficients.txt", sep="\\t")

## ************************************************************************************************************
## Fig S9 heatmap of genus diversity in ecological modes across time
## ************************************************************************************************************

master<-master.time.frame

nstages<-length(unique(master$stage))
necos<-length(levels(factor(master$ecospace)))

phylum<-matrix(nrow= nstages, ncol= necos )#creates empty matrix
class<-matrix(nrow= nstages, ncol= necos)
family<-matrix(nrow= nstages, ncol= necos)
genus<-matrix(nrow= nstages, ncol= necos)


levels(factor(master$phylum))[1:5]
levels(factor(master$class))
unique(master$taxon_genus)[1:6]


#to run at all taxonomic levels
for (i in 1:nstages) {
	for (j in 1:necos){
		
		hold<-master[which(master$ecospace==unique(master$ecospace)[j] & master$stage==unique(master$stage)[i]),]
		phylum[i,j]<-length(levels(factor(hold$phylum)))
		class[i,j]<-length(levels(factor(hold$class[!hold$class==""])))
		family[i,j]<-length(levels(factor(hold$family[!hold$family==""])))
		genus[i,j]<-length(levels(factor(hold$taxon_genus)))
	}
}



row.names(phylum)<-row.names(class)<-row.names(family)<-row.names(genus)<-unique(master$stage)
colnames(phylum)<-colnames(class)<-colnames(family)<-colnames(genus)<-unique(master$ecospace)



phylum<-t(phylum)### same here 
class<-t(class)### same here 
family<-t(family)### same here 
genus<-t(genus)### same here 

 
#order by total diversity in each ecospace across phanerozoic 
apply(genus, 1, sum)  #apply to array or matrix, object, margin = 1 = rows, but margin = 2 = columns etc, then function in this case sum

#organize by greatest total diversity at top
totaldiversity<-apply(genus, 1, sum)  #apply sum function to all rows


index2<-order(totaldiversity)   #get actual reverse of order
genus<-genus[index2,]  #reorder rows in ecos as defined by index2
apply(genus, 1, sum)  #to verify that it worked



#change numbers to words
mot<-character(length(row.names(genus)))
feed<-character(length(row.names(genus)))
hab<-character(length(row.names(genus)))


eco<-row.names(genus)
#assigning names to number codes for ecospaces

hab[grep("^1", eco)]<-"pelagic"
hab[grep("^2", eco)]<-"erect"
hab[grep("^3", eco)]<-"surficial"
hab[grep("^4", eco)]<-"semi-inf"
hab[grep("^5", eco)]<-"shallow"
hab[grep("^6", eco)]<-"deep"

feed[grep("1$", eco)]<-"suspension"
feed[grep("2$", eco)]<-"deposit"
feed[grep("3$", eco)]<-"mining"
feed[grep("4$", eco)]<-"grazing"
feed[grep("5$", eco)]<-"predatory"
feed[grep("6$", eco)]<-"other"

mot[grep("\\\\d1\\\\d", eco)]<-"fast"
mot[grep("\\\\d2\\\\d", eco)]<-"slow"
mot[grep("\\\\d3\\\\d", eco)]<-"fac. unattached"
mot[grep("\\\\d4\\\\d", eco)]<-"fac. attached"
mot[grep("\\\\d5\\\\d", eco)]<-"non-mot unattached"
mot[grep("\\\\d6\\\\d", eco)]<-"non-mot attached"

names<-paste(hab, mot, feed)

row.names(phylum)<-names
row.names(class)<-names
row.names(family)<-names
row.names(genus)<-names


## ********************************************************************************
## Make Final Figure
## ********************************************************************************
my.heat<-rev(heat.colors(100, alpha=1))

labels<-round(exp(seq(log(c(1, range(genus)[2]))[1], log(c(1, range(genus)[2]))[2],l=5)))
labels[1]<-0

genus.post.cam<-genus[,as.numeric(as.character(colnames(genus)))<486]

stg.names<-colnames(genus.post.cam)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
stg.names<-as.character(specify_decimal(as.numeric(stg.names), 1))
stg.names[seq(1,length(stg.names), by=2)]<-""
nstages<-length(stg.names)

pdf("FigS9_heatmap_no_cam.pdf", width=10, height=7) 
	layout(matrix(1:2,ncol=2), width = c(1,5), height = c(1,1))
	par(oma=c(4,0,2,10), mar=c(7,2,5,0))
	
	legend_image <- as.raster(matrix(rev(c("white",my.heat)), ncol=1))
	plot(c(0,2),c(0,1), type = 'n', axes = F,xlab = '', ylab = '')
	mtext("Number of Genera", cex=1.2, side=2)
	text(x=1.5, y = seq(0,1,l=5), labels = labels )
	polygon(c(0,0,1,1), y=c(0,1,1,0), col='black', lwd=3)
	rasterImage(legend_image, 0, 0, 1,1)


	par(mar=c(1,1,1,1))
	image(log(t(genus.post.cam)), xaxt='n', yaxt='n', col= my.heat)
	axis(4, seq(0,1, length.out=necos), labels=rownames(genus), las=1, cex.axis=0.6)
	# axis(1, seq(0,1, length.out=nstages), labels= stg.names, las=2, cex.axis=0.7)

	axis(1, seq(0,1, length.out=nstages)[stg.names!=""], labels= stg.names[stg.names!=""], las=2, cex.axis=0.75)
	axis(1, seq(0,1, length.out=nstages), labels= rep("", nstages), las=2, cex.axis=0.1, tck=-0.01)

	era.bounds<-which(colnames(genus.post.cam) %in% c(66.0,252.2))
	stg.pos<-seq(0,1, length.out=nstages)
	era.bounds.pos<-stg.pos[era.bounds]+diff(stg.pos)[1]/2

	abline(v=era.bounds.pos, lty=3)
	mtext("Geologic Time (Ma)", 1, line=3, cex=1.3)

	text(c(0.22, 0.62, 0.9), c(0.05, 0.05, 0.05), c("Paleozoic", "Mesozoic","Cenozoic"))

dev.off()




