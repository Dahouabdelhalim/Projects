###########################################################################
########      STATISTICAL ANALYSES on niche inferences  (part 3)    #######
###########################################################################


########################
###  Load libraries  ###

library(gtools)
library(vegan) 
library(stringr)
library(adegenet)




###########################################
###  Set working directories and paths  ###

mainDir = "PATH/to/workingDirectory/"
dataDir = "PATH/to/data/"
outDir = "PATH/to/outputs/"
setwd(mainDir) 

# load R objects generated for niche modeling (refer to Rscript "EcoGen-1-Niche_modeling.RDat")
load("PATH/to/EcoGen-1-Niche_modeling.RData")





## Load individual genotypes ##

# Distance matrix based on Rousset's a copmuted in SPAGEDI (Hardy & Vekemans 2002)
gdataDir = "PATH/to/data/spagedi_mcdu_indDist_ARousset.txt"
agen = read.table(gdataDir,sep="\\t",stringsAsFactors=F, h=T, row.names=1)
agen[is.na(agen)] = 0

# List of raw genotypes for RDA tests
gdataDir = "PATH/to/data/spagedi_mcdu.txt"
rskip = 3 ; nb = grep("END",readLines(gdataDir,warn=F))[1]-(rskip+2)
agen2 = read.table(gdataDir, skip=rskip, nrows=nb, sep="\\t", colClasses="character", stringsAsFactors=F,h=T)






#################################################################
###  Format biogeographical and ecological distance matrices  ###
###  -----------  for Mantel correlation tests   -----------  ###
#################################################################

# Formatting of ecological, geographic and genetic data
# NB: All genotyped accessions('aeg') are listed before the occurrence accessions('clim') in 'aegc' (See Rscript 'ELE-1-Niche_modeling.R'). Therefore, #ROW of genotyped accessions in 'biodata' & 'pca.cal$li' EQUALS #ROW in 'aegc'.

# Formatting individual genotypes into allelic frequency matrix (from - PopGen analyses -)
# 30 loci subset
loc = c("3Piso","AcylcoA","Aldolase","BRCA1","CENP_E","Cyclin","DMC1","E2F","eif3k","elongTS","galac","Glusyn","GTP","HS","KinMot","lipase","mtPorin","PCNA_2L","PHY_C","RFC1","Ribprot","RING","Sec24","SMC4","sucP","TransFac","vesATPase","VHS_GAT","waxy","zinc")
colnames(agen2)[which(colnames(agen2)=="X3Piso")] = "3Piso"
agen2$Pop = str_extract(agen2$Ind, "[A-Za-z]+" ) ; agen2$Cat = agen2$Sp

# input data frame
# Species names were abbreviated in the analyses for simplicity as follows:
di = c("Ca","Co","Ta","Um") # diploids: Ae. caudata (Ca) / Ae. comosa (Co) / Ae. tauschii (Ta) / Ae. umbellulata (Um)
tetra = c("Cr","Cy","Ge","Tr") # allopolyploids: Ae. crassa (Cr) / Ae. cylindrica (Cy) / Ae. geniculata (Ge) / Ae. triuncialis (Tr)
mcdu = c(di,tetra)

df = agen2[,c("Ind","Cat","Pop",loc)]
df.2x = df[which(df$Pop %in% di),]
df.4x = df[grep("Ge|Tr|Cy|Cr193|Cr195|Cr228|Cr375|Cr392",df$Ind),]
df.6x = df[which(df$Ind %in% c("Cr90","Cr137","Cr319","Cr408","Cr422")),]
# df to genind
obj.2x = df2genind(df.2x[,loc],ploidy=2,sep="/",ncode=3,NA.char="-9",ind.names=df.2x$Ind,pop=factor(df.2x$Pop))
obj.4x = df2genind(df.4x[,loc],ploidy=4,sep="/",ncode=3,NA.char="-9",ind.names=df.4x$Ind,pop=factor(df.4x$Pop))
obj.6x = df2genind(df.6x[,loc],ploidy=6,sep="/",ncode=3,NA.char="-9",ind.names=df.6x$Ind,pop=factor(df.6x$Pop))
obj = repool(obj.2x,obj.4x,obj.6x)
# matrix of alleles based on genotypes
X = tab(obj,freq=T,NA.method="mean")




#############################################################
###           Extract accessions that were used           ###
###        in both genetical & ecological analyses        ###
#############################################################


# Filter out non-genotyped accessions
aegc$ID = as.character(aegc$ID)
sp$ID = as.character(sp$ID)

nind = as.character(aegc$ID[grep(paste(mcdu,collapse="|"), aegc$ID)])
nind = mixedsort(intersect(nind,rownames(X))) #for allelic matrix derived from genotypes using adegenet



#############################################################################
###   Prepare the ecological / genetic / geographic individual datasets   ###
###                         for correlation tests                         ###
#############################################################################


for(i in c(mcdu,"pCr","pCy","pGe","pTr")){
  
  z = get(paste0("z.",i))
  
  nid = nind[grep(i,nind)] #extract samples of species i
  
  if(i=="pCr"){ nid = nind[grep("Ta|Co",nind)] #extract samples of corresponding diploid species for expected polyploid niche
  }else if(i=="pCy"){ nid = nind[grep("Ta|Ca",nind)]
  }else if(i=="pGe"){ nid = nind[grep("Um|Co",nind)]
  }else if(i=="pTr"){nid = nind[grep("Um|Ca",nind)]}

  nrow = which(aegc$ID %in% nid) #identifying corresponding rows in ecovar table
  
  
  ###  Extract data to generate lists for RDA tests  ###
  
  # ecological data
  ecoList = pca.cal$li[nrow,] ; rownames(ecoList) = aegc$ID[nrow] ; colnames(ecoList) = c("X","Y")
  ecoList = ecoList[nid,]
  assign(paste0("ecoList.",i), ecoList)
  
  # geographical data
  geoList = sp[which(sp$ID %in% nid), c("ID","X","Y")] ; rownames(geoList) = geoList$ID
  geoList = geoList[nid,-c(1)]
  assign(paste0("geoList.",i), geoList)
  
  # genetic data
  genList = X[nid,] # ensure all lists are similarly ordered
  genList = genList[, colSums(genList != 0) > 0]
  assign(paste0("genList.",i), genList)
  
  #print(identical(rownames(ecoList), rownames(genList))) #to ensure all lists are similarly ordered
  #print(identical(rownames(ecoList), rownames(geoList)))
  
  
  ###  Compute Euclidean distance between pairwise individuals' PCA coordinates for Mantel tests  ###
  
  ecoDist = as.matrix(dist(ecoList, method="euclidean")) # ecological data
  assign(paste0("ecoDist.",i), ecoDist)
  
  geoDist = as.matrix(dist(geoList, method="euclidean")) # geographical data
  assign(paste0("geoDist.",i), geoDist)
  
  genDist = agen[nid,nid] # genetic data (Rousset's a distance matrix)
  assign(paste0("genDist.",i), genDist)
  
  
  #print(identical(rownames(ecoDist), rownames(genDist))) #to ensure all matrices are similarly ordered
  #print(identical(rownames(ecoDist), rownames(geoDist)))
  #print(identical(colnames(ecoDist), colnames(genDist)))
  #print(identical(colnames(ecoDist), colnames(geoDist)))
  #print(identical(rownames(ecoDist), colnames(geoDist)))
  
  
  ###  Write results for use in Mantel and PCNM tests  ###
  write.table(ecoList, paste0(alyDir,"eco_xyPCA_",i,".txt"), sep="\\t", quote=F)
  write.table(geoList, paste0(alyDir,"geo_xyCoor_",i,".txt"), sep="\\t", quote=F)
  write.table(genList, paste0(alyDir,"gen_allfreq_",i,".txt"), sep="\\t", quote=F)
  write.table(ecoDist, paste0(alyDir,"eco_dPCA_",i,".txt"), sep="\\t", quote=F)
  write.table(geoDist, paste0(alyDir,"geo_dCoor_",i,".txt"), sep="\\t", quote=F)
  write.table(genDist, paste0(alyDir,"gen_dRousset_",i,".txt"), sep="\\t", quote=F)
  
}








##################################################
###  ----------------------------------------  ###
###  --------  Compute MANTEL TESTS  --------  ###
###  ----------------------------------------  ###
##################################################

pm.res = matrix(nrow=8,ncol=3,dimnames=list(mcdu,c("geo","eco","eco*geo")))


for(x in mcdu){
  
  geo = get(paste0("geoDist.",x))
  eco = get(paste0("ecoDist.",x))
  gen = get(paste0("genDist.",x))
  
  
  #mantel test: gen vs geo
  m.ac = mantel(xdis=gen, ydis=geo, permutations=20000)
  pm.res[x,"geo"] = paste0(round(m.ac$statistic,3)," (",round(m.ac$signif,3),")")
  
  #mantel test: gen vs eco
  m.ab = mantel(xdis=gen, ydis=eco, permutations=20000)
  pm.res[x,"eco"] = paste0(round(m.ab$statistic,3)," (",round(m.ab$signif,3),")")
  
  #partial mantel test: gen vs eco | geo
  pm = mantel.partial(xdis=gen, ydis=eco, zdis=geo, permutations=20000)
  pm.res[x,"eco*geo"] = paste0(round(pm$statistic,3)," (",round(pm$signif,3),")")
  
}

pm.res
write.table(pm.res, file=paste0(alyDir,"Mantel_tests.txt"), sep="\\t", quote=F)





###################################################
###  -----------------------------------------  ###
###  --------    Compute RDA tests    --------  ###
###  -----------------------------------------  ###
###################################################

###  Perform a multidimensional scaling on genetic data  ###
###  -- for use as community matrix in RDA/CCA tests --  ###

rda.res = matrix(nrow=length(mcdu), ncol=3, 
                 dimnames=list(mcdu,c("geo","eco","eco*geo")))


for(x in mcdu){
  
  gen = as.matrix(get(paste0("genList.",x)))
  geo = as.matrix(get(paste0("geoList.",x)))
  eco = as.matrix(get(paste0("ecoList.",x)))
  
  
  ## GEN ~ GEO
  dr1 = rda(gen ~ geo)
  aocc1 = anova.cca(dr1,permutations=20000, model="full")
  rda.res[x,"geo"] = paste0(round(aocc1$`F`[1],3)," (",round(aocc1$`Pr(>F)`[1],3),")")
  
  ## GEN ~ ECO
  dr2 = rda(gen ~ eco)
  aocc2 = anova.cca(dr2,permutations=20000, model="full")
  rda.res[x,"eco"] = paste0(round(aocc2$`F`[1],3)," (",round(aocc2$`Pr(>F)`[1],3),")")
  
  ## GEN ~ ECO | GEO
  dr3 = rda(gen ~ eco + Condition(geo))
  aocc3 = anova.cca(dr3,permutations=20000, model="full")
  rda.res[x,"eco*geo"] = paste0(round(aocc3$`F`[1],3)," (",round(aocc3$`Pr(>F)`[1],3),")")
  
}

rda.res
write.table(rda.res, paste0(alyDir,"rda_tests.txt"), quote=F, sep="\\t")
















