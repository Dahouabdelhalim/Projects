# Code associated with the paper "Advancing our understanding of ecological stability". S. Kéfi, V. Dominguez-Garcia, I. Donohue, C. Fontaine, E. Thébault, V. Dakos. Ecology Letters 

#-------------------------------------------------------------------------------
## Fig 1B
#-------------------------------------------------------------------------------

dist.metrics=read.csv("fig_hist_metricsfreq.csv",sep=",",header=F) 
hist(dist.metrics$V2,breaks=20,xlim=c(0,140),ylim=c(0,30),xlab="# times metrics used across papers",ylab="frequency",main="")

#-------------------------------------------------------------------------------
## Fig. 3 Alluvial plot
#-------------------------------------------------------------------------------
library(alluvial)

dat=read.csv("pert.csv",sep=";",header=T) 
str(dat)
dat$freq=1
dat$decade=round(dat$date / 10) * 10

unique(dat$p_natu_cat)

dat$p_natu_cat[which(dat$p_natu_cat=="env_dist")]="abiotic"
dat$p_natu_cat[which(dat$p_natu_cat=="process")]="abiotic"

unique(dat$p_type)
dat$p_type[which(dat$p_type=="pulset")]="pulse"
dat$p_type[which(dat$p_type=="pulse_t")]="pulse"
dat$p_type[which(dat$p_type=="Pulse_t")]="pulse"
dat$p_type[which(dat$p_type=="pulse_freq")]="pulse"
dat$p_type[which(dat$p_type=="pulsefreq")]="pulse"
dat$p_type[which(dat$p_type=="Press_t")]="press"
dat$p_type[which(dat$p_type=="press_t")]="press"
dat$p_type[which(dat$p_type=="pressrate")]="press"
dat$p_type[which(dat$p_type=="treat")]="null"

dim(dat[dat$p_type=="null",])
dim(dat[dat$p_type=="press",])
dim(dat[dat$p_type=="noise",])
dim(dat[dat$p_type=="pulse",])

dupli=dat[grep("\\\\.{1}",dat$compo_scatego),]
extra=gsub(".+\\\\.","",dupli$compo_scatego)
extra=c(extra,gsub("\\\\..+","",dupli$compo_scatego))
dupli=rbind.data.frame(dupli,dupli)
dupli$compo_scatego=as.factor(extra)

tripli=dat[grep("\\\\..+\\\\.",dat$compo_scatego),]
extra=gsub(".+\\\\.","",tripli$compo_scatego)
extra1=gsub("\\\\..+","",tripli$compo_scatego)
extra2=gsub(".*\\\\.(.*)\\\\..*", "\\\\1",tripli$compo_scatego)
extra=c(extra,extra1,extra2)
tripli=rbind.data.frame(tripli,tripli,tripli)
tripli$compo_scatego=as.factor(extra)

singl=dat[-grep("\\\\.",dat$compo_scatego),]

datt=rbind.data.frame(tripli,dupli,singl)

dat8=aggregate(freq~p_type+p_natu_cat+freq,data=datt,FUN=sum)
alluvial(dat8[,c(1:3)],freq=dat8$freq,
         col = as.numeric(dat8$p_type),
         border = as.numeric(dat8$p_type),
         ordering = list(
           NULL,
           NULL,
           NULL))

#-------------------------------------------------------------------------------
## Fig. S2
#-------------------------------------------------------------------------------

dat2=read.csv("papers.csv",sep=";",header=T) 
str(dat2)
temp = as.numeric(as.character(dat2$Size))
dat2$SizeNb = temp

dim(dat2[dat2$category_type=="field",])
dim(dat2[dat2$category_type=="theo",])
dim(dat2[dat2$category_type=="mixed",])

dat3 = dat2[,c("SizeNb","category_type")]
dat3 <- na.omit(dat3)
dat_emp = dat3[dat3$category_type=="field",]
dat_th = dat3[dat3$category_type=="theo",]
dat_mixed = dat3[dat3$category_type=="mixed",]

x<-hist(dat_th$SizeNb,breaks=2000,xlim=c(0,1000),ylim=c(0,140),xlab="System size",ylab="frequency",main="")
#print(x$breaks)
#print(x$counts) 

y = hist(dat_emp$SizeNb,breaks=80,xlim=c(0,400),ylim=c(0,100),xlab="System size",ylab="frequency",main="")
#print(y$breaks)
#print(y$counts) 

z = hist(dat_mixed$SizeNb,breaks=10,xlim=c(0,21),ylim=c(0,5),xlab="System size",ylab="frequency",main="")
#print(z2bis$breaks)
#print(z2bis$counts) 

#-------------------------------------------------------------------------------
## Fig. S3
#-------------------------------------------------------------------------------
library('stringr')
library('plyr')
library('ggplot2')

dat2=read.csv("papers.csv",sep=";",header=T) 

str(dat2) 
dat2$decade=round(dat2$date/10) * 10

# Number of metrics per paper
smetrics = dat2$compo_smetric
nbmetrics = as.numeric(lapply(smetrics,function(x) str_count(x,"\\\\.")+1))
dat2$nbmetrics=nbmetrics

dfc <- ddply(dat2, c("dat2$date", "dat2$nbmetrics"), "nrow", .drop = FALSE)
colnames(dfc) <- c("year","nbmetrics","nrow")
n = which(dfc$nrow==0)
dfc2 = dfc[-n,]
colnames(dfc2) <- c("year","nbmetrics","nrow")

ggplot(data = dfc2, aes(x = dfc2$year, y = dfc2$nbmetrics)) +
  geom_point(aes(size = dfc2$nrow))+ 
  labs(y = "# metrics per paper", x = "year")+ 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())

# Nb of perturbations per paper
pert=read.csv("pert.csv",sep=";",header=T) 
#print(unique(pert$ID))
x = seq(from = -0.5, to = 458.5 , by = 1)
hist.nbpert = hist(pert$ID,breaks=x)
nbpert = hist.nbpert$counts
dat2$nbpert = nbpert

dfc <- ddply(dat2, c("dat2$date", "dat2$nbpert"), "nrow", .drop = FALSE)
colnames(dfc) <- c("year","nbpert","nrow")
n = which(dfc$nrow==0)
dfc2 = dfc[-n,]
colnames(dfc2) <- c("year","nbpert","nrow")
ggplot(data = dfc2, aes(x = dfc2$year, y = dfc2$nbpert)) +
  geom_point(aes(size = dfc2$nrow))+ 
  labs(y = "# pert", x = "decade")+ 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())

### Nb of scales measured per paper
sscale = dat$compo_sscale
nbsscale = as.numeric(lapply(sscale,function(x) str_count(x,"\\\\.")+1))
dat$nbsscale=nbsscale

dfc <- ddply(dat, c("dat$date", "dat$nbsscale"), "nrow", .drop = FALSE)
colnames(dfc) <- c("year","nbsscale","nrow")
n = which(dfc$nrow==0)
dfc2 = dfc[-n,]
colnames(dfc2) <- c("year","nbsscale","nrow")

ggplot(data = dfc2, aes(x = dfc2$year, y = dfc2$nbsscale)) +
  geom_point(aes(size = dfc2$nrow))+ 
  labs(y = "# scales measured", x = "decade")+ 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())

# Nb of scales perturbed per paper
pscale = dat$compo_pscale
nbpscale = as.numeric(lapply(pscale,function(x) str_count(x,"\\\\.")+1))
dat$nbpscale=nbpscale

dfc <- ddply(dat, c("dat$date", "dat$nbpscale"), "nrow", .drop = FALSE)
colnames(dfc) <- c("year","nbpscale","nrow")
n = which(dfc$nrow==0)
dfc2 = dfc[-n,]
colnames(dfc2) <- c("year","nbpscale","nrow")

ggplot(data = dfc2, aes(x = dfc2$year, y = dfc2$nbpscale)) +
  geom_point(aes(size = dfc2$nrow))+ 
  labs(y = "# scales perturbed", x = "decade")+ 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())




