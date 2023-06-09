#### Asreml analysis for each LG #########################

rm(list=ls())

library(asreml)
library(nadiv)


#read phenotypic data
data <- read.delim( "Pheno.txt", header=T)
 
#select trait
trait <- "shellLength"


# add column for second ID, is necessary for including two relationship matrices
data$ID2 <- data$ID

# Genomic relationships (with IDs as row/column-names); method:Yang 2009
Gmat <- read.table(file="GRM/GRM.txt",
                   header=T)
# inverse of relationship matrix  in Asreml-format
Gmat.inv <- read.table("GRM/GRM.inv.txt", header=T)
attr(Gmat.inv, "rowNames") <- as.character(colnames(Gmat))


###### Models according to Robinson et al. 2013 #####################################
# for details see http://doi.wiley.com/10.1111/mec.12375
# model 3: G with all markers
mod3 <- asreml(fixed=eval(as.name(trait))~1+batch+sex, random=~ped(ID), 
                                ginverse = list(ID = Gmat.inv ), data=data, na.method.X="include")
# variance components: additive genetic variance (ped(ID)!ped)
# and residual variance (R!variance)
summary(mod3)$varcomp



for(j in 1:17){
  #select GRM based on focal linkage group/ genomic region
  input <- paste("GRM/Ginv_LG", j, ".txt", sep="")
  G.reg.inv <- read.table(input, header=T)
  #add attributes (= sample IDs)
  attr(G.reg.inv , "rowNames") <- as.character(colnames(Gmat))
  
  ### G without markers of focal region
  input.No <- paste("GRM/Ginv_LG", j,"NOT", ".txt", sep="")
  G.regNOT.inv <- read.table(input.No, header=T)
  attr(G.regNOT.inv, "rowNames") <- as.character(colnames(Gmat))
  
  ### run different models ##########################
  name <- paste("LG",j, "_", trait, sep="")
  # model 1: G excluding markers on focal region
  mod1 <- asreml(fixed=eval(as.name(trait)) ~1+batch+sex, random=~ped(ID), data=data, ginverse=list(ID=G.regNOT.inv), na.method.X="include")
  assign(paste( "mod1_",name, sep=""),mod1)
  
  # model 2: G without markers from focal region plus G with markers of focal region
  mod2 <- asreml(fixed=eval(as.name(trait))~1+batch+sex, random=~ped(ID)+giv(ID2), 
                 ginverse = list(ID = G.regNOT.inv , ID2 = G.reg.inv), data=data, na.method.X="include")
  assign(paste( "mod2_",name, sep=""),mod2)
  
  # model 4: G with all markers plus G with only regional markers
  mod4 <- asreml(fixed=eval(as.name(trait))~1+batch+sex, random=~ped(ID)+giv(ID2), 
                 ginverse = list(ID = Gmat.inv, ID2 = G.reg.inv), data=data, na.method.X="include")
  assign(paste( "mod4_",name, sep=""),mod4)
  
  
}


####### regional heritabilities and standard errors  ##########################
h <- c()
h.se <- c()
for(i in 1:17){
  mod.output <- eval(as.name(paste("mod2_LG",i, "_", trait, sep="")))
  h[i] <- as.numeric(pin(mod.output, h~ V2/(V1+V2+V3))[1]) 
  h.se[i] <-as.numeric(pin(mod.output, h~ V2/(V1+V2+V3))[2]) 
}


### significance testing using likelihood ratio tests ####
# test whether a genomic region explains significant amounts of variance (p.var.explained)
# test whether a genomic region explains significantly more variance than expected given its size (p.var.greater)

p.var.explained <- c()
p.var.greater <- c()

for(i in 1: 17){
  #select model output for respective region
  mod2 <- eval(as.name(paste("mod2_LG",i, "_", trait, sep="")))
  mod1 <- eval(as.name(paste("mod1_LG",i, "_", trait, sep="")))
  mod4 <- eval(as.name(paste("mod4_LG",i, "_", trait, sep="")))
  
  #LRT test 1
  p.var.explained[i] <- 1-pchisq(2*(mod2$loglik-mod1$loglik),1)
  #LRT test 2
  p.var.greater[i] <-  1-pchisq(2*(mod4$loglik-mod3$loglik),1)
}

#one-tailed test: variance cannot be below zero-> divide pvalues by 2
p.var.explained <- p.var.explained/2
p.var.greater <- p.var.greater/2


# Renumber linkage groups to be consistent with the  previous map that was based on
# a cross between Crab individuals
crab.lg <- c(12,2,1,3,17,6,4,7,9,14,10,5,13,15,11,8,16)

#results output
LG <- paste("LG", 1:17, sep="")
crabLG <- paste("LG.crab", crab.lg, sep="")
output.results <- cbind(LG, crabLG, h, h.se, p.var.explained, p.var.greater)


##### Correlation between variance explained  and linkage group size ######################################
# size = sum of contig length that were assigned to a linkage group (in Mb)
size <- c(187.86459, 208.70784, 177.81884, 120.92790, 127.53224, 106.70604, 102.98803,
          92.88234,  73.90762,  75.03283, 84.65308,  86.67678,  63.71682,  55.99972,  51.06867,
          43.45757,  31.99672 )

cor.test(size, h)



