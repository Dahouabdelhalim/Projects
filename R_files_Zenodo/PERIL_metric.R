## Calculating PERIL for modern-day Bivalvia
## Katie S. Collins et al
## 2018

## Clear up env ----
rm(list=ls())

## libraries ----
library(reshape)    # data management
library(data.table) # data management
library(GGally)     # plot extras
library(RcmdrMisc)  # helper functions
library(gtable)     # plot extras
library(grid)       # plot extras
library(gridExtra)  # plot extras
library(tidyverse)  # data management
library(forcats)    # data management of factors
library(Hmisc)      # pairwise correlations

# first make sure you've downloaded PERIL_rawdata.csv
# don't forget to set your working directory to the place you saved the csv

## read in data ----
PERIL <- read.csv("PERIL_rawdata.csv")

## data prep for PERIL calculation ----

## qhat = extinction rate 
## CHull = convex hulls of ranges
## CHullS = convex hulls of ranges clipped to the 200m shelf break (for comparison)
## range_SST = temperature range of realized thermal niche

## "r" = range-scaled
## "l" = logged 
## "i" = inverted 
## "rilCHull" or "rilCHullS" (using shelf range) is the final geog variable to use 
## "rQhatH" is the final extinction rate variable to use
## "riTemp" is the final thermal niche variable to use

PERIL <- ungroup(PERIL)%>%
  mutate(rQhatH=rescaler(qhatH, type="range"),
         lCHullS=log(CHullS),
         ilCHullS=lCHullS*-1,
         rilCHullS=rescaler(ilCHullS, type="range"), 
         lCHull=log(CHull),
         ilCHull=lCHull*-1,
         rilCHull=rescaler(ilCHull, type="range"),
         iTemp=range_SST*-1,
         riTemp=rescaler(iTemp, type="range"))

## PERIL ----
## add everything up
## "geogPERIL" is geography (full chulls) plus qhat <-- as used for Pliocene
## "shelfPERIL" is geography (shelf chulls) plus qhat
## "tempPERIL" is temperature range plus qhat
## "comboPERIL" is everything (full chulls) <-- PREFERRED METRIC
## "comboPERILshelf" is everything (shell chulls)
PERIL <- mutate(PERIL, 
                geogPERIL=rescaler(rQhatH+rilCHull, type="range"),
                shelfPERIL = rescaler(rQhatH+rilCHullS, type="range"),
                tempPERIL=rescaler(rQhatH+riTemp, type="range"), 
                comboPERILfull=rescaler(rQhatH+rilCHull+riTemp, type="range"),
                comboPERILshelf=rescaler(rQhatH+rilCHullS+riTemp, type="range"))

## correlations ----
correlates <- dplyr::select(PERIL, rQhatH, rilCHull, riTemp)%>%
  as.matrix()
PERILcorr <- rcorr(correlates, type="spearman")
PERILcorr$r
#write.csv(EVILcorr$r, file="PERILcomp_spearmancorr.csv", row.names=F)

## Percentiles and thresholds ----
## find the species in the top percentiles of interest
## 80th <-- preferred (codenamed "threatened") USED IN THE MAIN PAPER
## Comparison thresholds as used in the SOM:
## 50th (codenamed "vulnerable")
## 'td' is the 'Tridacna' threshold - the score of Tridacna derasa, 
## the lowest PERIL-scoring species officially assessed as Vulnerable by the IUCN
## (This is a different cutoff for each kind of PERIL as it has a slightly
## different score in each formulation)
PERIL <-  mutate(PERIL, 
                 cPERIL80=quantile(comboPERILfull, probs=0.80),
                 cPERIL50=quantile(comboPERILfull, probs=0.50),
                 cPERILtd=0.34592261,
                 csPERIL80=quantile(comboPERILshelf, probs=0.80),
                 csPERIL50=quantile(comboPERILshelf, probs=0.50),
                 csPERILtd=0.3528753,
                 sPERIL80=quantile(shelfPERIL, probs=0.80),
                 sPERIL50=quantile(shelfPERIL, probs=0.50),
                 sPERILtd=0.3572033,
                 gPERIL80=quantile(geogPERIL, probs=0.80),
                 gPERIL50=quantile(geogPERIL, probs=0.50),
                 gPERILtd=0.3384421,
                 tPERIL80=quantile(tempPERIL, probs=0.80),
                 tPERIL50=quantile(tempPERIL, probs=0.50),
                 tPERILtd=0.48338908)

## find which taxa are above or below the three thresholds
PERIL <-  mutate(PERIL,
                 threatenedC=comboPERILfull>cPERIL80,
                 vulnerableC=comboPERILfull>cPERIL50, tridacnaC=comboPERILfull>cPERILtd,
                 threatenedCS=comboPERILshelf>csPERIL80,
                 vulnerableCS=comboPERILshelf>csPERIL50, tridacnaCS=comboPERILshelf>csPERILtd,
                 threatenedG=geogPERIL>gPERIL80,
                 vulnerableG=geogPERIL>gPERIL50, tridacnaG=geogPERIL>gPERILtd,
                 threatenedS=shelfPERIL>sPERIL80,
                 vulnerableS=shelfPERIL>sPERIL50, tridacnaS=shelfPERIL>sPERILtd,
                 threatenedT=tempPERIL>tPERIL80,
                 vulnerableT=tempPERIL>tPERIL50, tridacnaT=tempPERIL>tPERILtd)

## numbers of threatened species

thrFamS <- filter(PERIL, threatenedS==T)%>%
  group_by(family)%>%
  mutate(total_thrS=n())%>%
  dplyr::select(family, total_thrS)%>%
  distinct()
PERIL <- left_join(PERIL, thrFamS, by="family")
rm(thrFamS)

thrFamC <- filter(PERIL, threatenedC==T)%>%
  group_by(family)%>%
  mutate(total_thrC=n())%>%
  dplyr::select(family, total_thrC)%>%
  distinct()
PERIL <- left_join(PERIL, thrFamC, by="family")
rm(thrFamC)

thrFamCS <- filter(PERIL, threatenedCS==T)%>%
  group_by(family)%>%
  mutate(total_thrCS=n())%>%
  dplyr::select(family, total_thrCS)%>%
  distinct()
PERIL <- left_join(PERIL, thrFamCS, by="family")
rm(thrFamCS)

thrFamG <- filter(PERIL, threatenedG==T)%>%
  group_by(family)%>%
  mutate(total_thrG=n())%>%
  dplyr::select(family, total_thrG)%>%
  distinct()
PERIL <- left_join(PERIL, thrFamG, by="family")
rm(thrFamG)

thrFamT <- filter(PERIL, threatenedT==T)%>%
  group_by(family)%>%
  mutate(total_thrT=n())%>%
  dplyr::select(family, total_thrT)%>%
  distinct()
PERIL <- left_join(PERIL, thrFamT, by="family")
rm(thrFamT)

## replace NAs with zero
PERIL[c("total_thrS", "total_thrC", "total_thrCS", "total_thrG", "total_thrT")][is.na(PERIL[c("total_thrS", "total_thrC", "total_thrCS", "total_thrG", "total_thrT")])] <- 0

## per family, proportionally, how many are threatened?
PERIL <- ungroup(PERIL)%>%
  mutate(threatPropS=total_thrS/totalExtantSpp, 
         threatPropC=total_thrC/totalExtantSpp, 
         threatPropCS=total_thrCS/totalExtantSpp, 
         threatPropG=total_thrG/totalExtantSpp, 
         threatPropT=total_thrT/totalExtantSpp)

## output the data
write.csv(PERIL, file="PERIL_metric.csv", row.names=F)

## FIGURES ----
## For Text Figures 2, 3 and SOM Figures 2, 3, contact KSC
## For Text Figure 1 and SOM Figure 4, see codefile 'paleoPERIL_metric.R'

## SOM Figure 1 ----
## Geogr ranges: convex hulls versus convex hulls clipped to shelf break
PERIL$clade <- fct_relevel(PERIL$clade, "Protobranchia", "Pteriomorpha", "Paleoheterodonta", "Heterodonta")

grid.arrange(
  ggplot(PERIL, aes(log(size), family))+
    geom_point(col="gray80")+
    geom_point(data=subset(PERIL, comboPERILfull > cPERIL80), col="black")+
    facet_wrap(~clade, scales = "free_y", nrow=1)+
    theme(axis.text.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()),
  ggplot(PERIL, aes(log(size), family))+
    geom_point(col="gray80")+
    geom_point(data=subset(PERIL, comboPERILshelf > sPERIL80, col="black"))+
    facet_wrap(~clade, scales="free_y", nrow=1)+
    theme(axis.text.y = element_blank()),
  nrow=2
)

## SOM Figure 5 ----
## K-S tests of size

pc80=quantile(PERIL$comboPERILfull, 0.80) 

above80=filter(PERIL, PERIL$comboPERILfull>=pc80)%>%
  mutate(threshold=paste("above"))
below80=filter(PERIL, PERIL$comboPERILfull < pc80)%>%
  mutate(threshold=paste("below"))

PERIL2 <- rbind(above80, below80)%>%
  mutate(clade=paste0("Bivalvia"))

ks80=ks.test(subset(PERIL2$size, PERIL2$threshold=="above"), subset(PERIL2$size, PERIL2$threshold=="below")) ## do a K-S test

fams <- dplyr::select(PERIL2, family, size, threshold)%>%
  group_by(family, threshold)%>%
  filter(n()>=2)%>%
  droplevels()%>%
  ungroup()

fams <- PERIL2 %>%
  group_by(family) %>%
  filter(any(threshold=="above"))%>%
  dplyr::select(clade, family, operational_genus, species, size, threshold)%>%
  mutate(nAbove=length(which(threshold=="above")), nBelow= length(which(threshold=="below")))%>%
  filter(!nBelow==0)

sp <- split(fams, fams$family, drop=TRUE)

KST <- lapply(sp, function(oo) {
  # oo <- sp[[14]
  oo <- data.frame(oo)
  kst <- ks.test(as.vector(oo[oo$threshold=="above", "size"]), as.vector(oo[oo$threshold=="below", "size"]))
  data.frame(family=unique(oo$family), kst$statistic, kst$p.value, nAbove=length(which(oo$threshold=="above")), nBelow= length(which(oo$threshold=="below")))
}) %>% bind_rows()

signif <- KST %>%
  ungroup()%>%
  mutate(significant=ifelse(kst.p.value>0.05, "no", "yes"), pvalue=kst.p.value)%>%
  dplyr::select(family, significant, pvalue)

fams <- left_join(fams, signif, by="family")%>%
  unique()

grid.arrange(
  ggplot(PERIL2, aes(threshold, log(size)))+
    geom_boxplot()+
    theme(axis.title.x=element_blank())+
    facet_wrap(~clade),
  ggplot(data=subset(fams, significant=="yes"), aes(threshold, log(size), group=threshold))+
    geom_jitter(col="grey50", size=0.5, alpha=0.5)+
    geom_boxplot(position = "dodge")+
    facet_wrap(~family, nrow=4),
  ncol=2, widths=1:2)
