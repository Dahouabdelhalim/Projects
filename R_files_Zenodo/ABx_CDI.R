# The following R code re-creates the analyses and output for the results section of (Brown KA, Khanafer N, Daneman N, Fisman DN. Antibiotics and the risk of community-associated Clostridium difficile infection (CDI): a meta-analysis. Antimicrob. Agents Chemother. 2013). The analyses are divided into sections A-F which correspond to the subheadings in the original article's results section.

# Load necessary packages
library(meta) # For metagen function
library(metafor) # Used for tau2 and meta-regression models

# Function: gets a table of effects with confidence intervals from an rma object
rma.eff<-function(x,expo=T,r=2){
  b<-x$b
  ci.lb<-x$ci.lb
  ci.ub<-x$ci.ub
  table<-cbind(b,ci.lb,ci.ub)
  ifelse(expo==T,out<-exp(table),out<-table)
  round(out,r)
}

# Load Datasets
primary<-read.csv("C:/Kevin/Dropbox/Dropbox/Kevin_2011/UT_RESEARCH/Systematic review/AAC/DRYAD.DATA/ABX_CDI_1.csv")
secondary<-read.csv("C:/Kevin/Dropbox/Dropbox/Kevin_2011/UT_RESEARCH/Systematic review/AAC/DRYAD.DATA/ABX_CDI_2.csv")

# Create factor variables
primary$f1 <- ifelse(primary$abx.type == 1, 1, 0)
primary$f2 <- ifelse(primary$abx.type == 2, 1, 0)
primary$f3 <- ifelse(primary$abx.type == 3, 1, 0)
primary$f4 <- ifelse(primary$abx.type == 4, 1, 0)
primary$f5 <- ifelse(primary$abx.type == 5, 1, 0)
primary$f6 <- ifelse(primary$abx.type == 6, 1, 0)
primary$f7 <- ifelse(primary$abx.type == 7, 1, 0)
primary$a0 <- ifelse(primary$author == "Delaney (2007)", 1, 0)
primary$a1 <- ifelse(primary$author == "Dial (2008)", 1, 0)
primary$a2 <- ifelse(primary$author == "Kuntz (2011)", 1, 0)
primary$a3 <- ifelse(primary$author == "Naggie (2011)", 1, 0)
primary$a4 <- ifelse(primary$author == "Wilcox (2008)", 1, 0)

# A. Pooled effects
a<-with(primary,metagen(lneff,selneff,studlab=author,comb.fixed=F,sm="OR")); 
summary(a)

# B. Antibiotic types
b<-with(primary,metagen(lneff,selneff,studlab=author,byvar=abx.type,comb.fixed=F,sm="OR")); 
summary(b)

# Calculate reduction in heterogeneity (tausq) from A to B.
m <- rma(lneff, sei=selneff, intercept=T, method="DL",  data = primary); m;
m <- rma(lneff, sei=selneff, mods = cbind(f1,f2,f3,f4,f5,f6,f7), intercept=F, method="DL",  data = primary); m;
# .62-->.27 = 0.56; 
(.62-.27)/.62

# C. Meta Regression 

# Accounting for study-level effects
res <- rma(lneff, sei=selneff, mods = cbind(f1,f2,f3,f4,f5,f6,f7,a1,a2,a3,a4), intercept=F, btt=c(8,9,10,11), method="DL",  data = primary); res; plot(res); rma.eff(res);
# tausq = 0.1

# Dial (2008) effect larger than other effects
res <- rma(lneff, sei=selneff, mods = cbind(f1,f2,f3,f4,f5,f6,f7,a1), intercept=F, data = primary, knha = TRUE); res; plot(res);
exp(.6579);exp(.2605);exp(1.0553)

# Removal of Dial (2008)
a<-with(subset(primary,a1==0),metagen(lneff,selneff,studlab=author,comb.fixed=F,sm="OR")); summary(a)
res <- rma(lneff, sei=selneff, mods = cbind(f1,f2,f3,f4,f5,f6,f7,a2,a3,a4), subset=a1==0, intercept=F,method="DL",btt=c(8,9,10), data = primary); res; plot(res); rma.eff(res);

# Subset of high study quality: Delaney (2007) and Kuntz (2011)
res <- rma(lneff, sei=selneff, mods = cbind(f1,f2,f3,f4,f5,f6,f7,a2), intercept=F, method="DL", btt=c(8), data = subset(primary,author %in% c("Delaney (2007)","Kuntz (2011)"))); res; rma.eff(res);  plot(res); 

a<-with(subset(primary,author %in% c("Delaney (2007)","Kuntz (2011)")),metagen(lneff,selneff,studlab=author,byvar=abx.type,comb.fixed=F,sm="OR")); summary(a)

# D. Publication Bias
res <- rma(lneff, sei=selneff, mods = cbind(f1,f2,f3,f4,f5,f6,f7), intercept=F, method="DL", data = primary); res; rma.eff(res);  plot(res); 
regtest(res, predictor = "vi")

# E. Risk Index
res <- rma(lneff, sei=selneff, mods = cbind(abx.index-1), intercept=F, method="DL", data = primary); res; rma.eff(res);  plot(res); 
exp(.88);exp(.76);exp(1.01)

# F. Secondary Analysis
a<-with(subset(secondary,abx.type %in% c(1,3,4,5,6)),metagen(lneff,selneff,studlab=author,byvar=abx.type,comb.fixed=F,comb.random=T,sm="OR")); summary(a); 
