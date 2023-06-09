# Processing raw data

set.seed(6202021)

#setwd("C:/Users/drhod/Desktop/ruminant_allometry")

family.colors <- c("#0080BD","#904010","#00BBE3","#BB662A","#005EA0","#60290C") # blue -> dark brown, with even darker tagged onto the back
names(family.colors) <- c("Cervidae", "Bovidae", "Tragulidae", "Moschidae", "DarkBlue", "DarkBrown")
col.pal.sunset <- colorRampPalette(GetColors(256, scheme = "sunset"))

# universal color blind palette, with a couple shades of gray at the end
CBP  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "Gray30", "Gray70")
col.pal <-  colorRampPalette(Blue2DarkRed18Steps)


# reading in metadata from Haber 2016
load("data/SpecimenInfo.Rdata")
load("data/comptreesL.Rdata")
load("data/ildefL.Rdata")
full.inventory1 <- read.csv(file = "data/full_inventory.csv",row.names = 1) # spreadsheet of covariate data (i.e., %grass, taxonomy, body mass estimations, etc.)

# reading in landmark data from Haber 2016
ruminant.lm.2d <- read.table(file = 'data/Haber2016_lm.tab',row.names = 1) # reading in the raw data
ruminant.lm <- arrayspecs(A = ruminant.lm.2d, p = dim(ruminant.lm.2d)[2]/3, k = 3) # making it an array


# for some reason the z axis is flipped in a lot of the specimens
# this part flips those specimens
r1 <- gpagen(ruminant.lm[8:47,,])
r1.pca <- gm.prcomp(r1$coords)
l <- which(sign(r1.pca$x[,1])==-1)
ruminant.lm[,3,l] <- -ruminant.lm[,3,l]



# the landmark data for Saiga tatarica and Megaloceros giganteus, rotated onto the first specimen
saiga <- read.pts(file = "data/saiga.pts")
saiga <- rotonto(x = ruminant.lm[1:25,,1], y = saiga, scale = T)$yrot

Megaloceros <- read.pts(file = "data/Megaloceros_noantlers.pts")
Megaloceros <- rotonto(x = ruminant.lm[1:25,,1], y = Megaloceros, scale = T)$yrot


# array of all the landmark data
ruminants <- bindArr(ruminant.lm[1:25,,],saiga,Megaloceros,along=3)

# all specimen names
si.FV <- c(as.character(specimenInfo$species.FV),"Saiga","Megaloceros") # index of the species shorthand for each specimen
si.FV.original <- si.FV
si.FV.ID <- c(as.character(names(specimenInfo$species.FV)),"Saiga","Megaloceros") # index of the species shorthand for each specimen

si.FV[319:447] <- "Odo_he_sp" # I'm binning the Odocoileus subspecies into a single species because the phylogeny I'm using (Upham 2019) doesn't have divergence time estimates for the subspecies
si.FV[448:572] <- "Odo_vi_sp"

dimnames(ruminants)[[3]] <- si.FV


# the phylogeny that Haber 2016 used (that I won't be using - just taking the IDs)
haber.tree <- comptreesL$FV

# putting the new Odocoileus species shorthand into the tree Annat Haber 2016 used
haber.tree2 <- add.tips(tree = haber.tree, tips = c("Odo_he_sp","Odo_vi_sp","Saiga","Megaloceros"), where = 240)
# the topology here is wrong but it doesn't matter because I'm not using this topology, I'm using Upham 2019, this part is just for sorting out names
haber.tree3 <- drop.tip(phy = haber.tree2, tip = haber.tree2$tip.label[28:33])
haber.tree <- haber.tree3
names(haber.tree$tip.label)[210] <- "Odocoileus hemionus"
names(haber.tree$tip.label)[211] <- "Odocoileus virginianus"
names(haber.tree$tip.label)[212] <- "Saiga tatarica"
names(haber.tree$tip.label)[213] <- "Megaloceros giganteus"


# making a key between the species' shorthand from Haber 2016 and the scientific names
name.id <- matrix(data = NA, nrow = length(haber.tree$tip.label),ncol=3)
name.id[,1] <- haber.tree$tip.label
name.id[,2] <- names(haber.tree$tip.label)
name.id[,3] <- str_replace(names(haber.tree$tip.label), pattern = " ", replacement = "_")


# making a vector of all the species present in the dataset
specs.FV <- unique(si.FV)
# putting names on this vector
for(i in 1:length(specs.FV)){
  name <- name.id[which(name.id[,1] == specs.FV[i]), 3]
  names(specs.FV)[i] <- str_replace(name, patter = " ", replacement = "_")
}


# processing Upham et al 2019's mammal phylogenies

# upham's phylogeny
# downloading 1000 trees from the posterior distribution of their estimate
# http://vertlife.org/phylosubsets/
upham.pd <- read.tree(file = "data/ruminant_distribution_UphamEtAl2019.phy", keep.multi = TRUE)

# taking the maximum clade credibility tree from the distribution
upham.tree <- mcc(upham.pd)

# Axis_calamianensis, Neotragus_moschatus, Rusa_marianna, Damaliscus_korrigum need to be corrected
## add Axis calamianensis (node 229), Damaliscus_korrigum (node 362), Rusa_marianna (222) into the phylogeny
names(specs.FV)[28] <- "Neotragus_moschatus"

upham.tree2 <- add.tips(tree = upham.tree, tips = c("Axis_calamianensis", "Damaliscus_korrigum", "Rusa_marianna"), where = c(229,362,222))
upham.tree2$tip.label[60] <- "Neotragus_moschatus"
# > identical(rownames(full.inventory1), upham.tree2$tip.label[match(rownames(full.inventory1),upham.tree2$tip.label)])
# [1] TRUE
upham.tree.full <- ladderize(upham.tree2)


#
#making vector of the full scientific names
haber.names <- names(specs.FV)

# tree with only the specimens in Annat Haber's dataset
R.tree <- keep.tip(phy = upham.tree2, tip = haber.names)

# making sure things are in the correct order for my own sanity
specs.FV <- specs.FV[match(R.tree$tip.label, haber.names)]


w <- table(si.FV)
# plot(sort(w),ylab="specimen count", xlab="species") ; grid()



# making a list of all the species with at least 27 specimens
Nmin <- 27
specs20 <- c()
for(i in 1:length(specs.FV)){
  if(length(which(si.FV == specs.FV[[i]])) > (Nmin-1)){specs20[[i]] <- specs.FV[[i]]}
specs20 <- specs20[lengths(specs20) != 0]
}
specs20 <- unlist(specs20)
specs20 <- specs20[-which(specs20 == "Cep_ni7")]
# there's an ontogenetic signal in Cep_ni7, it's not just conspecific adults


for(i in 1:length(specs20)){
  names(specs20)[i] <- name.id[which(name.id[,1]==specs20[i]),3]
}
# just to confirm that everything is in the correct order :-)

# the tree with only species represented by at least 20 specimens
tree.20 <- keep.tip(R.tree, tip = names(specs20))
specs20 <- specs20[match(tree.20$tip.label, names(specs20))]



# making a list of superimpositions of each of the species with >27 specimens
proc.list.20 <- c()
proc.list.20.ID <- c()
Ps.20 <- array(NA, dim = c(75,75,length(specs20)))
for(i in 1:length(specs20)){
proc.list.20[[i]] <- gpagen(A = ruminants[1:25,,which(si.FV == specs20[[i]])],print.progress = FALSE, verbose = T)
Ps.20[,,i] <- var(two.d.array(proc.list.20[[i]]$coords))
proc.list.20.ID[[i]] <- si.FV.ID[which(si.FV == specs20[[i]])]
}
names(proc.list.20) <- names(specs20)


# list of superimpositions of all species, regardless of whether they have over 20 specimens
proc.list.all <- c()
proc.list.ID.all <- c()
for(i in 1:length(specs.FV)){
# if there's only one specimen representing the species, just superimpose itself by duplicating it, so that all the objects in 'proc.list.all' are gpagen objects
if(length(which(si.FV==specs.FV[i]))==1){
  proc.list.all[[i]] <- gpagen(A = ruminants[1:25,,c(which(si.FV == specs.FV[i]),which(si.FV == specs.FV[i]))],print.progress = FALSE)
  }else{
proc.list.all[[i]] <- gpagen(A = ruminants[1:25,,which(si.FV == specs.FV[i])],print.progress = FALSE)}
proc.list.ID.all[[i]] <- si.FV.ID[which(si.FV == specs.FV[[i]])]
  }


# for all the species
R.all <- array(dim = c(25, 3, length(specs.FV)))
R.size <- c()
# the Saiga has a different scale becuase I downloaded it from morphosource
for(i in 1:length(specs.FV)){
  R.all[,,i] <- proc.list.all[[i]]$consensus
  R.size[[i]] <- mean(proc.list.all[[i]]$Csize)
}


R <- gpagen(R.all)$coords
# > identical(R.tree$tip.label,names(specs.FV))
# > TRUE
dimnames(R)[[3]] <- names(specs.FV)
names(R.size) <- names(specs.FV)
R.size <- unlist(R.size)
R <- R[,,R.tree$tip.label]
R.size <- R.size[R.tree$tip.label]


#full.inventory1[names(R.size),"max_BM"]
mBM <- full.inventory1[,"max_BM"]
names(mBM) <- rownames(full.inventory1)
#mBM <- mBM[-which(mBM == 0)]
# first, need to predict the saiga & irish elk centroid sizes from  the body size ~ centroid size regression of all the other species
BM.CS.names <- intersect(names(R.size)[-which(names(R.size) == "Saiga_tatarica" | names(R.size) == "Megaloceros_giganteus")], names(mBM))

skull.size.body.mass.reg <- lm(log(R.size)[BM.CS.names]~log(mBM[BM.CS.names]))

# y = mx + b
# body.size = coef * CS size + intercept
R.size <- log(R.size)

# estimating the centroid sizes of Saiga and the irish elk based off of their maximum body masses
# there's a linear relationship between centroid size of the skull and body mass
R.size["Saiga_tatarica"] <- log(mBM["Saiga_tatarica"])*skull.size.body.mass.reg$coefficients[2] + skull.size.body.mass.reg$coefficients[1]
R.size["Megaloceros_giganteus"] <- log(mBM["Megaloceros_giganteus"])*skull.size.body.mass.reg$coefficients[2] + skull.size.body.mass.reg$coefficients[1]



# principal components analysis of mean shape for each species
R.pca <- gm.prcomp(A=R,phy=R.tree)

# all the individual specimens in PC morphospace
as.proc <- gpagen(ruminants)
allspecs <- as.proc$coords
allspecs.size <- as.proc$Csize

# projecting the individual specimens into the interspecific morphospace
allspecs.scores <- matrix(nrow=dim(allspecs)[3],ncol=dim(R.pca$rotation)[2])
for(i in 1:dim(allspecs)[3]){
  x <- rotonto(mshape(R),allspecs[,,i],scale=F)$yrot - mshape(R)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  allspecs.scores[i,] <- a.preds[1,] %*% R.pca$rotation
}




###
# plotTree(R.tree,node.numbers=TRUE,fsize=0.6)

ancestors <- R.pca$ancestors[c('137','172','131'),]
ancestors <- arrayspecs(ancestors, p = 25, k = 3)
cervid.ancestor <- ancestors[,,1]
bovid.ancestor <- ancestors[,,2]
common.ancestor <- ancestors[,,3]
# shapes of the ancestors of Bovidae, Cervidae, and Ruminantia


# list of species and their family
Cervidae2 <- Descendants(R.tree, node = 137, type =  "tips")[[1]]
Bovidae2 <- Descendants(R.tree, node = 172, type =  "tips")[[1]]
Moschidae2 <- Descendants(R.tree, node = 170, type =  "tips")[[1]]
f.fam2 <- dimnames(R)[[3]]
f.fam2[Cervidae2] <- 'Cervidae'
f.fam2[Bovidae2] <- 'Bovidae'
f.fam2[Moschidae2] <- 'Moschidae'
f.fam2[-c(Cervidae2,Bovidae2,Moschidae2)] <- 'Tragulidae'
f.fam2.nl <- f.fam2
f.fam2 <- as.factor(f.fam2)
fam <- f.fam2
tm.drop <- which(fam=="Tragulidae"|fam=="Moschidae")
f.fam3 <- as.factor(f.fam2.nl[-tm.drop])
fam.drop <- as.factor(f.fam3)
names(fam) <- dimnames(R)[[3]]
names(fam.drop) <- dimnames(R[,,-tm.drop])[[3]]

fam <- fam[names(specs.FV)]

# subset of tree w 20 specs
#plotTree(tree.20,node.numbers=TRUE)
Cervidae <- Descendants(tree.20, node = 53, type =  "tips")[[1]]
Bovidae <- Descendants(tree.20, node = 65, type =  "tips")[[1]]
f.fam <- tree.20$tip.label
names(f.fam) <- tree.20$tip.label
f.fam[Cervidae] <- 'Cervidae'
f.fam[Bovidae] <- 'Bovidae'
f.fam[-c(Cervidae,Bovidae)] <- 'Tragulidae'
f.fam <- as.factor(f.fam)

f.fam <- f.fam[names(specs20)]


# pruning data to just cervids
cervid.tree <- keep.tip(R.tree, tip = R.tree$tip.label[which(fam=='Cervidae')])
cervid.lm <- R[,,which(fam=='Cervidae')]
cervid.lm <- cervid.lm[,,cervid.tree$tip.label] # making sure things are in the correct order
cervid.size <- R.size[which(fam=='Cervidae')]
cervid.size <- cervid.size[cervid.tree$tip.label] # making sure things are in the correct order
cervid.pca <- gm.prcomp(A = cervid.lm, phy = cervid.tree) # interspecific pca of cervid dataset

# bovids
bovid.tree <- keep.tip(R.tree, tip = R.tree$tip.label[which(fam=='Bovidae')])
bovid.lm <- R[,,which(fam=='Bovidae')]
bovid.lm <- bovid.lm[,,bovid.tree$tip.label]
bovid.size <- R.size[which(fam=='Bovidae')]
bovid.size <- bovid.size[bovid.tree$tip.label]
bovid.pca <- gm.prcomp(A = bovid.lm, phy = bovid.tree)

# dataframes


# R.df
# full 130 species with morphometric data
# centroid size, face length, family, subfamily, tribe, phylogeny, procrustes coordinates
# make sure they're all in the right order
R.tree <- ladderize(R.tree)
R.df.n <- R.tree$tip.label
FL <- vapply(X = names(R.size), FUN = function(x){face.length(R[,,x])}, FUN.VALUE = 1) # calculate face length for each species
NR <- vapply(X = names(R.size), FUN = function(x){nasal.retraction(R[,,x])}, FUN.VALUE = 1) # calculate face length for each species
R.df <- geomorph.data.frame(coords = R[,,R.df.n],
                            csize = R.size[R.df.n],
                            facelength = FL[R.df.n],
                            nasalretraction = NR[R.df.n],
                            family = as.factor(full.inventory1[R.df.n,"family"]),
                            subfamily = as.factor(full.inventory1[R.df.n,"subfamily"]),
                            tribe = as.factor(full.inventory1[R.df.n,"tribe"]),
                            tree = R.tree)


 M <- mshape(R.df$coords)
# just the taxonomy
ruminant.taxonomy <- as.data.frame(R.df[c("family", "subfamily", "tribe")]); rownames(ruminant.taxonomy) <- R.df.n
rt <- apply(ruminant.taxonomy,2,as.character)
rownames(rt) <- rownames(ruminant.taxonomy)
ruminant.taxonomy2 <- ruminant.taxonomy[names(specs20),]
rt2 <- ruminant.taxonomy2



#R.test <- procD.pgls(f1 = coords~csize, phy = tree, data = R.df); summary(R.test)
#R.test2 <- procD.lm(f1 = coords~csize*tribe, data = R.df); summary(R.test2)


# R2.df - ruminant dataframe #2
# interspecific datasets with only the species that have both landmark, body mass, % grass, and horn length measurements
ad <- rownames(full.inventory1[which(full.inventory1$max_BM > 0 & !is.na(full.inventory1$percent_grass) & !is.na(full.inventory1$HL_male_cm)),])
all.data.names <- intersect(dimnames(R)[[3]], ad)

R2.tree <- ladderize(keep.tip(R.tree, tip = all.data.names))
R2.df.n <- R2.tree$tip.label

# horn length residuals
# hl.log <- log(full.inventory1[R2.df.n,"HL_male_cm"])
# names(hl.log) <- R2.df.n
# size.hl <- lm(hl.log~R.size[R2.df.n])
# #plot(R.size[R2.df.n],log(full.inventory1[R2.df.n,"HL_male_cm"]), xlab = "skull size", ylab = "horn length") ; abline(size.hl)
#
# hl.resids <- size.hl$residuals
# names(hl.resids) <- R2.df.n

R2.df <- geomorph.data.frame(coords = R[,,R2.df.n],
                            csize = R.size[R2.df.n],
                            maxBM = log(full.inventory1[R2.df.n,"max_BM"]),
                            percentgrass = full.inventory1[R2.df.n,"percent_grass"],
                           # hlresids = hl.resids[R2.df.n],
                            family = as.factor(full.inventory1[R2.df.n,"family"]),
                            subfamily = as.factor(full.inventory1[R2.df.n,"subfamily"]),
                            tribe = as.factor(full.inventory1[R2.df.n,"tribe"]),
                            tree = R2.tree)

#R2.test <- procD.pgls(f1 = coords~csize*hlresids, phy = tree, data = R2.df); anova(R2.test) # relative horn length has a small but significant effect on cranial morphology

# full tree with body mass data
size.full <- full.inventory1[intersect(upham.tree.full$tip.label,rownames(full.inventory1)),"max_BM"] ; names(size.full) <- upham.tree.full$tip.label ; size.full <- log(size.full)




# first, interspecific allometry... to measure intraspecific PC1s against
# PGLS - axis of evolutionary allometry of cranial shape in ruminants
R.size.pgls <- procD.pgls(f1 = coords~csize, phy = tree, data = R.df)
summary(R.size.pgls)

R.size.ls <- procD.lm(f1 = coords~csize, data = R.df)
summary(R.size.ls)



pts.col <- c()
pts.pch <- c()
for(i in 1:dim(R)[3]){
  if(fam[i] == "Cervidae"){pts.col[[i]] <- addTrans(family.colors[1],255*1) ; pts.pch[[i]] <- 21}else{
  if(fam[i] == "Bovidae"){pts.col[[i]] <- addTrans(family.colors[2],255*1) ; pts.pch[[i]] <- 22}else{
  if(fam[i] == "Moschidae"){pts.col[[i]] <- addTrans(family.colors[4],255*1) ; pts.pch[[i]] <- 25}
  else{pts.col[[i]] <- addTrans(family.colors[3],255*1) ; pts.pch[[i]] <- 24}}}
}
pts.col <- unlist(pts.col) ; names(pts.col) <- dimnames(R)[[3]]
pts.pch <- unlist(pts.pch) ; names(pts.pch) <- dimnames(R)[[3]]





# INTRASPECIFIC METRICS
# eigenvalue dispersion, relative standard deviation of eigenvalues
rSDE <- lapply(proc.list.20, function(x){
  A <- x$coords
  eigen.var(A, sd = TRUE, rel = TRUE, sample = dim(A)[3])
}) %>% unlist


# effect size of Vrel, Connaway & Adams 2022
Vrel <- lapply(proc.list.20, function(x){
  integration.Vrel(x$coords)$ZR
}) %>% unlist

# allometry trajectory
allo.vector <- R.size.pgls$pgls.coefficients[2,]

# projected variance of CREA
Mr <- mshape(R.df$coords)
Vproj <- lapply(proc.list.20, function(x){
  A <- x$coords
  for(i in 1:dim(A)[3]){
    A[,,i] <- rotonto(x = Mr, y = A[,,i])$yrot
    }
  vcv <- var(two.d.array(A))
  projected.variance(vcv,allo.vector)
}) %>% unlist


# calculating family-specific allometric trajectories:
fams <- c("Bovidae", "Cervidae")
allo.fams <- matrix(nrow=length(fams), ncol = length(allo.vector))
Vproj.specific <- numeric(length(Vproj))
for(i in 1:nrow(allo.fams)){
  NAMES <- R.df.n[which(R.df$family == fams[i])]
  PHYfoo <- keep.tip(R.df$tree, tip = NAMES)
  COORDS <- R.df$coords[,,NAMES] ; SIZE <- R.df$csize[NAMES]
  pglsF <- procD.pgls(f1 = COORDS~SIZE, phy = PHYfoo)
  VEC <- pglsF$pgls.coefficients[2,]
  PART <- which(ruminant.taxonomy2[,1] == fams[i])
  for(j in PART){
      A1 <- proc.list.20[[j]]$coords
      for(k in 1:dim(A1)[3]){
        A1[,,k] <- rotonto(x = Mr, y = A1[,,k])$yrot
      }
      vcv <- var(two.d.array(A=A1))
      Vproj.specific[j] <- projected.variance(vcv,VEC)
  }
}
names(Vproj.specific) <- names(Vproj)



# which nodes are the ancestors to certain taxonomic groups
pts <- R.pca$x
anc <- R.pca$anc.x
ind <- match(R.tree$tip.label, rownames(pts))
ind.N <- rownames(pts)
pts <- pts[ind, ]
z <- rbind(pts,anc)
edges <- as.matrix(R.tree$edge)

fam.MRCA <- c(132,137,170,172)
names(fam.MRCA) <- unique(rt[,1])

subfam.MRCA <- c(Tragulinae = 132,Cervinae = 138,
                 Capreolinae = 150,Hydropotinae = 166,
                 Moschinae = 170,Neotraginae = 175,
                 Aepycerotinae = 174, Reduncinae = 177,
                 Antilopinae = 185, Cephalophinae = 211,
                 Pantholopinae = 224, Caprinae = 225,
                 Alcelaphinae = 242, Hippotraginae = 245, Bovinae = 249)

# calculating the divergence magnitudes
subfam.divergence <- numeric(length(nrow(pts)))
for(i in 1:nrow(pts)){
  NM <- rownames(pts)[i]
  sf <- rt[NM,2]
  N <- which(edges[,2] == i)
  subfam.divergence[i] <- euc.dist(z[edges[N,2],], z[subfam.MRCA[sf],])
}
names(subfam.divergence) <- rownames(pts)


# While accounting for clade age
# not presenting this though, the results are nearly identical
# this is == evolutionary rate

# subfam.divergence2 <- numeric(length(nrow(pts)))
# for(i in 1:nrow(pts)){
#   NM <- rownames(pts)[i]
#   sf <- rt[NM,2]
#   N <- which(edges[,2] == i)
#   R.d <- distRoot(R.df$tree)[1]
#   age1 <- distNodes(R.df$tree, node = c(131, subfam.MRCA[sf]))[1,2]
#   age <- R.d - age1
#   subfam.divergence2[i] <- euc.dist(z[edges[N,2],], z[subfam.MRCA[sf],]) / age
# }
# names(subfam.divergence2) <- rownames(pts)


# the angle between the direction of divergence and CREA (the allometric trajectory)
div.traj <- matrix(NA,nrow=nrow(R.pca$x),ncol=length(allo.vector))
rownames(div.traj) <- rownames(R.pca$x)
for(i in 1:nrow(div.traj)){
  NM <- rownames(pts)[i]
  sf <- rt[NM,2]
  desc.shape <- flatten(R.df$coords[,,NM])
  anc.shape <- R.pca$ancestors[as.factor(subfam.MRCA[sf]),]
  div.traj[i,] <-  desc.shape - anc.shape
}
div.crea.angle.sf <- apply(div.traj,1,FUN=function(x)vector.angle(x,allo.vector))


Vproj.div <- numeric(length(specs20))
names(Vproj.div) <- names(specs20)
for(i in 1:length(specs20)){
  A <- proc.list.20[[i]]$coords
  vcv <- var(two.d.array(A))
  Vproj.div[i] <- projected.variance(vcv,div.traj[names(specs20)[i],])
}



## Computing Projected Variance distributions within Ruminantia, Bovidae, and Cervidae:
#
R.Ppca <- gm.prcomp(R.df$coords, phy = R.df$tree, GLS= TRUE, transform = TRUE)
R.paca <- gm.prcomp(R.df$coords, phy = R.df$tree, align.to.phy = TRUE)

# all ruminants
##########

PC1 <- R.pca$rotation[,1] ; PC2 <- R.pca$rotation[,2]
pPC1 <- R.Ppca$rotation[,1] ; pPC2 <- R.Ppca$rotation[,2]
PaC1 <- R.paca$rotation[,1] ; PaC2 <- R.paca$rotation[,2]

nsim <- 499
Nmin <- 27 # minimum of 27 specimens
# originally I bootstrapped this, but now for each species I'm just taking the average of 499 random skewers through P
Vproj.random1 <- c()
for(i in 1:nsim){
  x <- numeric(length(proc.list.20))
    for(j in 1:length(x)){
      l <- dim(proc.list.20[[j]]$coords)[3]
        vcv <- var(two.d.array(A=proc.list.20[[j]]$coords))#[,,sample(1:l, Nmin)]))
        random.Z <- scalar1(sample(seq(-1,1,length.out=200),size=75))
        x[j] <- projected.variance(vcv,random.Z)
    }
  Vproj.random1[[i]] <- x
}

Vproj.random1 <- Vproj.random1 %>% do.call(rbind, .)
Vproj.random <- apply(Vproj.random1,2,mean)



  Vproj.1 <- setNames(numeric(length(proc.list.20)), nm = names(specs20))
  PC1D <- Vproj.1
  PC2D <- Vproj.1
  pPC1D <- Vproj.1
  pPC2D <- Vproj.1
  PAC1D <- Vproj.1
  PAC2D <- Vproj.1

    for(j in 1:length(proc.list.20)){
      A2d <- two.d.array(A=proc.list.20[[j]]$coords)
      vcv <- var(A2d)
      Vproj.1[j] <- projected.variance(vcv,allo.vector, scale=T)
      pPC1D[j] <- projected.variance(vcv,pPC1, scale=T)
      pPC2D[j] <- projected.variance(vcv,pPC2, scale=T)
      PAC1D[j] <- projected.variance(vcv,PaC1, scale=T)
      PAC2D[j] <- projected.variance(vcv,PaC2, scale=T)
      PC1D[j] <- projected.variance(vcv,PC1, scale=T)
      PC2D[j] <- projected.variance(vcv,PC2, scale=T)
  }


vproj.df <- data.frame(Random = Vproj.random, CREA = Vproj, Divergence = Vproj.div, PC1 = PC1D, PC2 = PC2D,
                       pPC1 = pPC1D, pPC2 = pPC2D, PAC1 = PAC1D, PAC2 = PAC2D)


##########


Ps <- array(NA,dim=c(75,75,length(specs20)))
for(i in 1:length(specs20)){
  Ps[,,i] <- var(two.d.array(A=proc.list.20[[i]]$coords))
}
dimnames(Ps)[[3]] <- names(specs20)



# Cervidae
##########

cervid.pPCA <- gm.prcomp(A=cervid.lm,phy=cervid.tree,
                         GLS=T,transform = T)
cervid.PACA <- gm.prcomp(A=cervid.lm,phy=cervid.tree,align.to.phy = T)

PC1 <- cervid.pca$rotation[,1] ; PC2 <- cervid.pca$rotation[,2]
pPC1 <- cervid.pPCA$rotation[,1] ; pPC2 <- cervid.pPCA$rotation[,2]
PaC1 <- cervid.PACA$rotation[,1] ; PaC2 <- cervid.PACA$rotation[,2]

Ps.cervidae <- Ps[,,which(rt2[,1]=="Cervidae")]

nsim <- 499
Nmin <- 27 # minimum of 27 specimens
Vproj.random1 <- c()
Cervid.L <- dim(Ps.cervidae)[3]
for(i in 1:nsim){
  x <- numeric(Cervid.L)
    for(j in 1:dim(Ps.cervidae)[3]){
      l <- dim(proc.list.20[[j]]$coords)[3]
        vcv <- Ps.cervidae[,,j]
        random.Z <- scalar1(sample(seq(-1,1,length.out=200),size=75))
        x[j] <- projected.variance(vcv,random.Z)
    }
  Vproj.random1[[i]] <- x
}

Vproj.random1 <- Vproj.random1 %>% do.call(rbind, .)
Vproj.random <- apply(Vproj.random1,2,mean)



  Vproj.1 <- setNames(numeric(Cervid.L), nm = dimnames(Ps.cervidae)[[3]])
  PC1D <- Vproj.1
  PC2D <- Vproj.1
  pPC1D <- Vproj.1
  pPC2D <- Vproj.1
  PAC1D <- Vproj.1
  PAC2D <- Vproj.1

    for(j in 1:dim(Ps.cervidae)[3]){
      vcv <- Ps.cervidae[,,j]
      Vproj.1[j] <- projected.variance(vcv,allo.vector, scale=T)
      pPC1D[j] <- projected.variance(vcv,pPC1, scale=T)
      pPC2D[j] <- projected.variance(vcv,pPC2, scale=T)
      PAC1D[j] <- projected.variance(vcv,PaC1, scale=T)
      PAC2D[j] <- projected.variance(vcv,PaC2, scale=T)
      PC1D[j] <- projected.variance(vcv,PC1, scale=T)
      PC2D[j] <- projected.variance(vcv,PC2, scale=T)
  }


  vproj.df.cervidae <- data.frame(Random = Vproj.random, CREA = Vproj.1, Divergence = Vproj.div[names(Vproj.1)], PC1 = PC1D, PC2 = PC2D,
                       pPC1 = pPC1D, pPC2 = pPC2D, PAC1 = PAC1D, PAC2 = PAC2D)


##########



# Bovidae
##########

bovid.pPCA <- gm.prcomp(A=bovid.lm,phy=bovid.tree,
                         GLS=T,transform = T)
bovid.PACA <- gm.prcomp(A=bovid.lm,phy=bovid.tree,align.to.phy = T)

PC1 <- bovid.pca$rotation[,1] ; PC2 <- bovid.pca$rotation[,2]
pPC1 <- bovid.pPCA$rotation[,1] ; pPC2 <- bovid.pPCA$rotation[,2]
PaC1 <- bovid.PACA$rotation[,1] ; PaC2 <- bovid.PACA$rotation[,2]

Ps.bovidae <- Ps[,,which(rt2[,1]=="Bovidae")]

nsim <- 499
Nmin <- 27 # minimum of 27 specimens
Vproj.random1 <- c()
bovid.L <- dim(Ps.bovidae)[3]
for(i in 1:nsim){
  x <- numeric(bovid.L)
    for(j in 1:dim(Ps.bovidae)[3]){
      l <- dim(proc.list.20[[j]]$coords)[3]
        vcv <- Ps.bovidae[,,j]
        random.Z <- scalar1(sample(seq(-1,1,length.out=200),size=75))
        x[j] <- projected.variance(vcv,random.Z)
    }
  Vproj.random1[[i]] <- x
}

Vproj.random1 <- Vproj.random1 %>% do.call(rbind, .)
Vproj.random <- apply(Vproj.random1,2,mean)


  Vproj.1 <- setNames(numeric(bovid.L), nm = dimnames(Ps.bovidae)[[3]])
  PC1D <- Vproj.1
  PC2D <- Vproj.1
  pPC1D <- Vproj.1
  pPC2D <- Vproj.1
  PAC1D <- Vproj.1
  PAC2D <- Vproj.1

    for(j in 1:dim(Ps.bovidae)[3]){
      vcv <- Ps.bovidae[,,j]
      Vproj.1[j] <- projected.variance(vcv,allo.vector, scale=T)
      pPC1D[j] <- projected.variance(vcv,pPC1, scale=T)
      pPC2D[j] <- projected.variance(vcv,pPC2, scale=T)
      PAC1D[j] <- projected.variance(vcv,PaC1, scale=T)
      PAC2D[j] <- projected.variance(vcv,PaC2, scale=T)
      PC1D[j] <- projected.variance(vcv,PC1, scale=T)
      PC2D[j] <- projected.variance(vcv,PC2, scale=T)
  }


  vproj.df.bovidae <- data.frame(Random = Vproj.random, CREA = Vproj.1, Divergence = Vproj.div[names(Vproj.1)], PC1 = PC1D, PC2 = PC2D,
                       pPC1 = pPC1D, pPC2 = pPC2D, PAC1 = PAC1D, PAC2 = PAC2D)


##########



getMeaningfulPCs(R.pca$d,n = dim(R.df$coords)[3]) # use the first 5 PCs

R.bm <- mvBM(tree = R.df$tree, data = R.pca$x[,1:5], model = "BM1")


save.image(file = "ruminant_data.RData") ; save.image(file = "data/ruminant_data.RData")
