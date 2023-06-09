# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# Vector correlation
# Andrew M. Scott

#### NOTE before running:
#### Script depends on objects generated from running FC_Plasticity_DE and FC_Artificial_Selection_DE Scripts!

library(scales)
set.seed(1234)

## Generate full list of plasticity hits based on both the main effect of treatment and the ET-IT contrast
edger.hits.ITET <- edger.hits[!(row.names(edger.hits) %in% row.names(edger.trt.hits)), ]
plast.hits.vc <- rbind(edger.trt.hits[,-6], edger.hits.ITET[,-6]) # 375 total hits

#### Significant artificial selection genes
#### FROM Nebula - NOTE, logFCs are in the same direction for AS and plasticity
####    Meaning: positive logFC values = higher FC group > low FC group (in counts)

AS.hits.vc <- neb.as.hits # 919 total hits


## Below portion of script from Zinna et al. (2018) (Mol. Ecol.), code written by Ian Dworkin
#__________________________________________________________________________________________________________________#
#The function below calculate Euclidean Distances (the L2 norm) so we can use unit vectors in the estimation of the vector correlation.
#These functions follow the logic of Kuruvilla et al 2002, and were adapted from Pitchers et al 2013
PD <- function(x) { 
  sqrt(t(x)%*%x)}
#this function gives me the vector correlation and angle, and vector magnitude ratio, alpha, between two vectors
#alpha of 1 means that vector 2 is the same length as vector 1. lower than one means that it is smaller, greater than 1 means that vector 2 is larger
ang.vec.alph <- function(vec1, vec2) {
  vec1 <- vec1 - mean(vec1)
  vec2 <- vec2 - mean(vec2)
  vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
  vec.angle <- acos(vec.cor)*(180/pi)
  vec.alpha <- PD(vec1)/PD(vec2) # ID > RZ The absolute value was unnecessary since these are both magnitudes.
  vec.ED <- PD(vec2-vec1) #Subtract vector one from vector two and then calculate PD for Euclidean Distance. 
  ## Do not use vec.ED unless examining mean gene expr vectors. Specifically do not use for fold change
  return(c(vector.cor=vec.cor, vec.angle=vec.angle, vec.alpha=vec.alpha, vector.ED=vec.ED))} 
#The function below will generate two vectors of the same length as the experimental vectors provided to it with the argument 'sig_vec'
RandomVectorSampler <- function(vec1, vec2, sig_vec) {
  sig_vector_length <- length(sig_vec)
  
  random_index <- sample(1:length(vec1), size = sig_vector_length)
  sub_vec1 <- na.omit(vec1[random_index])
  sub_vec2 <- na.omit(vec2[random_index])
  
  ang.vec.alph(sub_vec1, sub_vec2)
}
comment(RandomVectorSampler) <- c("vec1 = vector of complete genelist comparison 1",
                                  "vec2 = vector of complete genelist comparison 2",
                                  "sig_vec = vector of significant genes, to get length")
#__________________________________________________________________________________________________________________#

edger.trt$table$logFC # plasticity full list
neb.mod.1$summary$logFC_trt.nebU # AS full list

#######################################################################
######## VC - the overlap genes

overlapping.hits.expr <- get(load(file = "overlapping_hits_expr.Rdata"))

# observed vc, alpha

vc.overlap.result <- as.data.frame(rbind(
  ang.vec.alph(overlapping.hits.expr$logFC_Plast, overlapping.hits.expr$logFC_AS)), 
  rownames="estimates")

# Sampling to generate empirical distributions to compare our observed vector correlation, alphas to

# align full lists by FBgn to allow for correct sampling

full.logFC.list <- merge(edger.trt$table, neb.mod.1$summary, by = "row.names", all.x = F)
row.names(full.logFC.list) <- full.logFC.list$Row.names
full.logFC.list <- full.logFC.list[, -1]

# Sampling to generate empirical distributions to compare our observed vector correlation, alphas to
# Random sampling

random_vectors_overlap <- t(replicate(10000,
                                      RandomVectorSampler(full.logFC.list$logFC, 
                                                          full.logFC.list$logFC_trt.nebU, 
                                                          overlapping.hits.expr$logFC_AS)))

# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

obs.vs.empir.dist.overlap <- rbind(vc.overlap.result, 
                                   apply(random_vectors_overlap, MARGIN = 2, 
                                         quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))

########################################################################
########## VC - unshrunken AS estimates
##########
#### Vector correlation of significant plasticity genes, and corresponding genes in AS list

# corresponding genes in AS list

plast.hits.in.AS.list <- neb.mod.1$summary[(row.names(neb.mod.1$summary) %in% row.names(plast.hits.vc)), ]

# combine lists to line up genes

vc.sig.plast.genes <- merge(plast.hits.vc, plast.hits.in.AS.list, by = "row.names")
row.names(vc.sig.plast.genes) <- vc.sig.plast.genes$Row.names
vc.sig.plast.genes <- vc.sig.plast.genes[, -1]

# Vector correlation on "observed" logFC values between sig plast hits and corresponding genes in AS analysis
# NOTE - alphas larger than 1 means that vector 1 is bigger (I.e. first argument is bigger)

vc.sig.plast.result <- as.data.frame(rbind(
  ang.vec.alph(vc.sig.plast.genes$edger_logFC, vc.sig.plast.genes$logFC_trt.nebU)), 
  rownames="estimates")

# Random sampling

random_vectors_Plast_AS <- t(replicate(10000,
                                       RandomVectorSampler(full.logFC.list$logFC, 
                                                           full.logFC.list$logFC_trt.nebU, 
                                                           vc.sig.plast.genes$edger_logFC)))

plot(density(random_vectors_Plast_AS[,1]))


# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

Plast.AS.VC.res.vs.empir.dist <- rbind(vc.sig.plast.result, 
                                       apply(random_vectors_Plast_AS, MARGIN = 2, 
                                             quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))


###########
#### Vector correlation of significant AS genes, and corresponding genes in Plasticity list

AS.hits.in.plast.list <- edger.trt$table[(row.names(edger.trt$table) %in% row.names(AS.hits.vc)), ]

# 40 genes not present in plasticity list, so these will be removed during combining/aligning step:

vc.sig.AS.genes <- merge(AS.hits.vc, AS.hits.in.plast.list, by = "row.names", all.x = F)
row.names(vc.sig.AS.genes) <- vc.sig.AS.genes$Row.names
vc.sig.AS.genes <- vc.sig.AS.genes[, -1]

# Vector correlation on "observed" logFC values between sig AS hits and corresponding genes in plasticity analysis

vc.sig.AS.result <- as.data.frame(rbind(
  ang.vec.alph(vc.sig.AS.genes$logFC, vc.sig.AS.genes$logFC_trt.nebU)), 
  rownames="estimates")


# Random sampling

random_vectors_AS_hits <- t(replicate(10000,
                                      RandomVectorSampler(full.logFC.list$logFC, 
                                                          full.logFC.list$logFC_trt.nebU, 
                                                          vc.sig.AS.genes$logFC)))

plot(density(random_vectors_AS_hits[,1]))


# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

Plast.AS.VC.res.vs.empir.dist.2 <- rbind(vc.sig.AS.result, 
                                         apply(random_vectors_AS_hits, MARGIN = 2, 
                                               quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))

########################################################################
########## VC - shrunken AS estimates (APEGLM)
##########
#### Vector correlation of significant plasticity genes, and corresponding genes in AS (shrunken) list

edger.trt$table$logFC # plasticity full list
apeglm.AS.MD$trt.UvD.logFC # AS full list

# corresponding genes in AS list

plast.hits.in.AS_shrunken.list <- apeglm.AS.MD[(row.names(apeglm.AS.MD) %in% row.names(plast.hits.vc)), ]

# combine lists to line up genes

vc.sig.plast.genes_ASshrunk <- merge(plast.hits.vc, plast.hits.in.AS_shrunken.list, by = "row.names")
row.names(vc.sig.plast.genes_ASshrunk) <- vc.sig.plast.genes_ASshrunk$Row.names
vc.sig.plast.genes_ASshrunk <- vc.sig.plast.genes_ASshrunk[, -1]

# observed vc
vc.sig.plast.result_ASshrunk <- as.data.frame(rbind(
  ang.vec.alph(vc.sig.plast.genes_ASshrunk$edger_logFC, vc.sig.plast.genes_ASshrunk$trt.UvD.logFC)), 
  rownames="estimates")

# generate emprical dist

full.logFC.list_ASshrunk <- merge(edger.trt$table, apeglm.AS.MD, by = "row.names", all.x = F)
row.names(full.logFC.list_ASshrunk) <- full.logFC.list_ASshrunk$Row.names
full.logFC.list_ASshrunk <- full.logFC.list_ASshrunk[, -1]

# Random sampling

random_vectors_Plast_AS_shrunk <- t(replicate(10000,
                                              RandomVectorSampler(full.logFC.list_ASshrunk$logFC, 
                                                                  full.logFC.list_ASshrunk$trt.UvD.logFC, 
                                                                  vc.sig.plast.genes_ASshrunk$edger_logFC)))
plot(density(random_vectors_Plast_AS_shrunk[,1]))

# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

Plast.AS_shrunk.VC.res.vs.empir.dist <- rbind(vc.sig.plast.result_ASshrunk, 
                                              apply(random_vectors_Plast_AS_shrunk, MARGIN = 2, 
                                                    quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))


#### Vector correlation of significant AS (shrunken) genes, and corresponding genes in Plast list

# corresponding genes in Plast list

AS_shrunken.hits.in.plast.list <- edger.trt$table[(row.names(edger.trt$table) %in% row.names(apeglm.AS.hits)), ]

# combine lists to line up genes

vc.sig.ASshrunk.genes_wPlast <- merge(apeglm.AS.hits, AS_shrunken.hits.in.plast.list, by = "row.names")
row.names(vc.sig.ASshrunk.genes_wPlast) <- vc.sig.ASshrunk.genes_wPlast$Row.names
vc.sig.ASshrunk.genes_wPlast <- vc.sig.ASshrunk.genes_wPlast[, -1]

# observed vc
vc.sig.ASshrunk.result <- as.data.frame(rbind(
  ang.vec.alph(vc.sig.ASshrunk.genes_wPlast$logFC, vc.sig.ASshrunk.genes_wPlast$trt.UvD.logFC)), 
  rownames="estimates")

# Random sampling

random_vectors_AS_shrunk_Plast <- t(replicate(10000,
                                              RandomVectorSampler(full.logFC.list_ASshrunk$logFC, 
                                                                  full.logFC.list_ASshrunk$trt.UvD.logFC, 
                                                                  vc.sig.ASshrunk.genes_wPlast$trt.UvD.logFC)))

# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

AS_shrunk_plast.VC.res.vs.empir.dist <- rbind(vc.sig.ASshrunk.result, 
                                              apply(random_vectors_AS_shrunk_Plast, MARGIN = 2, 
                                                    quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))


########################################################################
########## VC - EdgeR of AS and Plasticity
##########
#### Vector correlation of significant plasticity genes, and corresponding genes in AS (shrunken) list

edger.trt$table$logFC
edgeR.AS.trt.p.0.25.cutoff$table$logFC

# full edgeR AS list
edgeR.AS.trt.fullList$table$logFC

# corresponding genes in AS list

plast.hits.in.AS_edgeR.list <- edgeR.AS.trt.fullList$table[(row.names(edgeR.AS.trt.fullList$table) %in% row.names(plast.hits.vc)), ]

# combine lists to line up genes

vc.sig.plast.genes_ASedgeR <- merge(plast.hits.vc, plast.hits.in.AS_edgeR.list, by = "row.names")
row.names(vc.sig.plast.genes_ASedgeR) <- vc.sig.plast.genes_ASedgeR$Row.names
vc.sig.plast.genes_ASedgeR <- vc.sig.plast.genes_ASedgeR[, -1]

# observed vc
vc.sig.plast.result_ASedgeR <- as.data.frame(rbind(
  ang.vec.alph(vc.sig.plast.genes_ASedgeR$edger_logFC, vc.sig.plast.genes_ASedgeR$logFC)), 
  rownames="estimates")

# generate emprical dist

full.logFC.list_ASedgeR <- merge(edger.trt$table, edgeR.AS.trt.fullList$table, by = "row.names", all.x = F)
row.names(full.logFC.list_ASedgeR) <- full.logFC.list_ASedgeR$Row.names
full.logFC.list_ASedgeR <- full.logFC.list_ASedgeR[, -1]

# Random sampling

random_vectors_Plast_AS_edgeR <- t(replicate(10000,
                                             RandomVectorSampler(full.logFC.list_ASedgeR$logFC.x, 
                                                                 full.logFC.list_ASedgeR$logFC.y, 
                                                                 vc.sig.plast.genes_ASedgeR$edger_logFC)))

# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

Plast.AS_edgeR.VC.res.vs.empir.dist <- rbind(vc.sig.plast.result_ASedgeR, 
                                             apply(random_vectors_Plast_AS_edgeR, MARGIN = 2, 
                                                   quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))


#### Vector correlation of significant AS (shrunken) genes, and corresponding genes in Plast list

# corresponding genes in Plast list

AS_edgeR.hits.in.plast.list <- edger.trt$table[(row.names(edger.trt$table) %in% row.names(edgeR.AS.trt.p.0.25.cutoff$table)), ]

# combine lists to line up genes

vc.sig.ASedgeR.genes_wPlast <- merge(edgeR.AS.trt.p.0.25.cutoff$table, AS_shrunken.hits.in.plast.list, by = "row.names")
row.names(vc.sig.ASedgeR.genes_wPlast) <- vc.sig.ASedgeR.genes_wPlast$Row.names
vc.sig.ASedgeR.genes_wPlast <- vc.sig.ASedgeR.genes_wPlast[, -1]

# observed vc
vc.sig.ASedgeR.result <- as.data.frame(rbind(
  ang.vec.alph(vc.sig.ASedgeR.genes_wPlast$logFC.y, vc.sig.ASedgeR.genes_wPlast$logFC.x)), 
  rownames="estimates")

# Random sampling

random_vectors_AS_edgeR_Plast <- t(replicate(10000,
                                             RandomVectorSampler(full.logFC.list_ASedgeR$logFC.y, 
                                                                 full.logFC.list_ASedgeR$logFC.x, 
                                                                 vc.sig.ASedgeR.genes_wPlast$logFC.y)))

# Comparison of observed VC, alpha, etc to empirical dist 95% conf int

AS_edgeR_plast.VC.res.vs.empir.dist <- rbind(vc.sig.ASedgeR.result, 
                                             apply(random_vectors_AS_edgeR_Plast, MARGIN = 2, 
                                                   quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))


