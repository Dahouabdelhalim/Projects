require(OutFLANK)
require(vcfR)
require(adegenet)
require(gdm)
require(gradientForest)
require(foreach)
require(doParallel)
require(pbapply)
require(gdata)
require(data.table)
require(PresenceAbsence)
require(ROCR)
require(modEvA)
require(dplyr)
require(grid)
require(gridExtra)
require(gtools)
require(stringr)
require(reshape2)
require(hierfstat)
require(clusterGeneration)


options(scipen = 999)

#cores <- 7
#cl <- makeCluster(cores)
#registerDoParallel(cl)

setwd("E:/Research_AJL/TTT_Offset_Vulnerability_GF_Sims/")

seed_table <- read.table("seeds_source_R.txt")

for(w in 1:nrow(seed_table)){

  seed_table <- read.table("seeds_source_R.txt")

  
seed=seed_table$V2[w]
#seed=8435562691650
#seed=commandArgs(trailingOnly=T)
fit<-read.table(paste("results/SLiM_output/Sim_sum/",seed,"_Freq_ML_WF.txt",sep=""), header=T)

N<-data.frame(fit$n[1]*100)
colnames(N)<-"N"

if(ncol(fit)==613){
  (specs<-data.frame(seed,fit[1,1:2],N,fit[1,3:12]))
  plotTitle <- paste(colnames(specs)[1],":",specs[[1]],", ", colnames(specs)[2],":",specs[2],", ",colnames(specs)[3],":",specs[3],", ",colnames(specs)[4],":",specs[4],",",colnames(specs)[5],":",specs[5],"\\n",colnames(specs)[6],":",specs[6],", ",colnames(specs)[7],":",specs[7],", ",colnames(specs)[8],":",specs[8],", ",colnames(specs)[9],":",specs[9],", ",colnames(specs)[10],":",specs[10], sep="")
  fitt<-data.frame(t(fit[,-1:-13]))
  } else {
  (specs<-data.frame(seed,fit[1,1:2],N,fit[1,3:13])) # for case4
  plotTitle <- paste(colnames(specs)[1],":",specs[[1]],", ", colnames(specs)[2],":",specs[2],", ",colnames(specs)[3],":",specs[3],", ",colnames(specs)[4],":",specs[4],",",colnames(specs)[5],":",specs[5],"\\n",colnames(specs)[6],":",specs[6],", ",colnames(specs)[7],":",specs[7],", ",colnames(specs)[8],":",specs[8],", ",colnames(specs)[9],":",specs[9],", ",colnames(specs)[10],":",specs[10],", ",colnames(specs)[11],":",specs[11],sep="") # for case4
  fitt<-data.frame(t(fit[,-1:-14])) # for case4
}

gen_nam <- paste("Gen",fit$Generation,sep="")

colnames(fitt)<-gen_nam

fitt$Location <- factor(rep(paste("P",1:100,sep=""),6), levels = unique(rep(paste("P",1:100,sep=""),6)))

fitt$Type <- as.factor(c(rep("Fit",100),rep("Freq",100),rep("Phen1",100), rep("Phen2",100),rep("Env1",100), rep("Env2",100)))

#VCF files are filtered with vcftools, as it is much faster than R. The filtering for MAF > 0.05 is accomplished with the following code:
vcf1_filt <- read.vcfR(paste("results/SLiM_output/VCF_files/T1_",seed,"_filtered_subset.recode.vcf",sep=""))
geno1_filt <- vcf1_filt@gt[,-1] # Remove 1st column, which is 'Format'
position1_filt <- as.numeric(getPOS(vcf1_filt)) # Positions in bp
chromosome1_filt <- as.numeric(getCHROM(vcf1_filt)) # Chromosome information

No_A<-unname(dim(vcf1_filt)[2]) #Get the number of filtered alleles

rm(vcf1_filt)
gc()

### Convert VCF to 012 format ####
# Character matrix containing the genotypes
# individuals in columns

#Create Genotype matrix
G1f <- matrix(NA, nrow = nrow(geno1_filt), ncol = ncol(geno1_filt))
G1f[geno1_filt %in% c("0/0", "0|0")] <- 0
G1f[geno1_filt %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G1f[geno1_filt %in% c("1/1", "1|1")] <- 2

#Occasionally have a "2|#" genotype show up because of overlapping mutations from SLiM. This procuded NA genotypes which breaks 
#the per locus FST calc. later, so this just gets rid of overlapped sites 
if(any(geno1_filt[!geno1_filt %in% c("0/0", "0|0", "0/1", "1/0", "1|0", "0|1", "1/1", "1|1")])){
  writeLines("Replaced unusual genotypes with '0|0'")
}
G1f[!geno1_filt %in% c("0/0", "0|0", "0/1", "1/0", "1|0", "0|1", "1/1", "1|1")] <- 0

#stopCluster(cl)

#Check number of duplicate positions
sum(duplicated(position1_filt))
position1_filt[duplicated(position1_filt)]
duplicated(position1_filt[duplicated(position1_filt)])
#Read in positions of m2 mutations
#pos_T1 <- scan(paste("results/SLiM_output/Sim_sum/",seed,"_T1_Pos_ML_WF.txt",sep=""))

muts <- scan(paste("results/SLiM_output/",seed,"_casual_muts2.txt",sep=""))
if(any(duplicated(muts))){
  print("Duplicates in positions!")
  break
}

#select_pos_num <-unique(pos_T1)
select_pos <- paste("M",muts,sep="")

duplic_pos <- position1_filt[duplicated(position1_filt)]
#select_pos <- paste("M",unique(pos_T1),sep="")

#select_pos[select_pos%in%duplic_pos]

for(i in 1:length(position1_filt)){
  if(duplicated(position1_filt)[i]){
    position1_filt[i]<-position1_filt[i]+0.5
    #print(position1_filt[i])
  }
}

Start<-seq(1,ncol(geno1_filt),10)
#Stop<-seq(100,ncol(geno1_filt),10)
Stop<-seq(10,ncol(geno1_filt),10)

Pop_afreq1<-NULL
for(i in 1:100){
  Pop_afreq1<-rbind(Pop_afreq1,rowSums(G1f[,Start[i]:Stop[i]])/(2*ncol(G1f[,Start[i]:Stop[i]])))
  #print(paste("Start:", Start[i],", Stop:", Stop[i]))
}
Pop_afreq1<-data.frame(Pop_afreq1)
colnames(Pop_afreq1)<-paste("M",position1_filt,sep="")

#######################################################
# Set up environmental vartiables for Gradient Forest
#######################################################
cg<-read.table(paste("results/SLiM_output/CG_files/",seed,"_fitnessmat_pop.txt",sep=""),header=F)

colnames(cg)<-rep(paste("H",seq(1:100),sep=""))
cg$Transplant<-paste("T",seq(1:100),sep="")

cg_df <- melt(cg,  id.vars ="Transplant",  variable.name = "Home",
              value.name = "Fitness")

cg_df$Transplant<-factor(cg_df$Transplant,levels=unique(cg_df$Transplant))

#Subset the environmental variables for the generation you're considering (make sure the M2 AF and environmental data are not being compare across generations)
envPop1<-data.frame(fitt[fitt$Type=="Env1",gen_nam[length(gen_nam)-3]]) #300 years prior to the end of the simulation is taken as the "before environmental shift" time
names(envPop1) <- "envPop1"

fakeEnv1 <- envPop1$envPop1 + rnorm(nrow(envPop1),0,1.3)

cor(envPop1$envPop1, fakeEnv1)

envPop2<-data.frame(fitt[fitt$Type=="Env2",gen_nam[length(gen_nam)-3]]) #300 years prior to the end of the simulation is taken as the "before environmental shift" time
names(envPop2) <- "envPop2"

fakeEnv2 <- envPop2$envPop2 + rnorm(nrow(envPop2),0,1.3)
# this standard deviation generally produces a correlation between 0.3 and 0.6
cor(envPop2$envPop2, fakeEnv2)

Pop <- levels(fitt$Location)
Popsenv <- data.frame(Pop, envPop1,envPop2,fakeEnv1,fakeEnv2)

nfake <- 10

Popsenv[,6:(5+nfake)] <- NA

cov1 <- genPositiveDefMat(nfake,covMethod="unifcorrmat" )
head(cov1)

a<- mvrnorm(nrow(Popsenv),mu=rep(0, nfake), Sigma=cov1$Sigma)

Popsenv[,6:(5+nfake)] <- a

# head(Popsenv)
# tail(Popsenv)

sel_env_cols <- 2:3
all_env_cols <- 2:ncol(Popsenv)
sel_env_cols_plus2 <- 2:5

cov_allEnv <- cov(Popsenv[,all_env_cols])
cov_selEnv <- cov(Popsenv[,sel_env_cols])
cov_selEnv_plus2 <- cov(Popsenv[,sel_env_cols_plus2])

start_time <- Sys.time()
for (i in 1:nrow(cg_df)){
  # get the row in PopsenvStnd for the common garden genotype
  row1 = which(Popsenv==gsub("T","P",as.character(cg_df$Transplant[i])))
  # get the row in PopsenvStnd for the genotype source
  row2 = which(Popsenv==gsub("H","P",as.character(cg_df$Home[i])))
  
  # Look up the envi (all var)
  (envpop1_all <- Popsenv[row1,all_env_cols])
  (envpop2_all <- Popsenv[row2,all_env_cols])
  
  # Look up the envi selected only
  (envpop1_sel <- Popsenv[row1,sel_env_cols])
  (envpop2_sel <- Popsenv[row2,sel_env_cols])
  
  # Look up the envi selected plus 2 correlated env
  (envpop1_sel_plus2 <- Popsenv[row1,sel_env_cols_plus2])
  (envpop2_sel_plus2 <- Popsenv[row2,sel_env_cols_plus2])  
  
  
  ### Calculate the environmental distance between the two rows
  
  # Euclidean distance for selective environments
  cg_df$EdSelEnv[i] <- dist(rbind(envpop1_sel,
                                  envpop2_sel))
  
  # Mahalanobis distance for selective environments
  cg_df$MdSelEnv[i] <- mahalanobis(as.numeric(envpop1_sel), 
                                   as.numeric(envpop2_sel), 
                                   cov_selEnv)
  
  # Euclidean distance for ALL environments
  cg_df$EdAllEnv[i] <- dist(rbind(envpop1_all,
                                  envpop2_all))
  
  # Mahalanobis distance for ALL environments
  cg_df$MdAllEnv[i] <- mahalanobis(as.numeric(envpop1_all), 
                                   as.numeric(envpop2_all), 
                                   cov_allEnv)
  
  # Euclidean distance for selective environments + 2 env
  cg_df$EdSelEnvPlus2[i] <- dist(rbind(envpop1_sel_plus2,
                                       envpop2_sel_plus2))
  
  # Mahalanobis distance for selective environments + 2 env
  cg_df$MdSelEnvPlus2[i] <- mahalanobis(as.numeric(envpop1_sel_plus2), 
                                        as.numeric(envpop2_sel_plus2), 
                                        cov_selEnv_plus2) 
}
end_time <- Sys.time()
writeLines("Environmental distance calculation:")
print(end_time - start_time)
# this take a few minutes

#write.csv(cg_df,"CGfit_Dist.csv",row.names = F)

# Env_sel1<-rep(envTab[,"envPop1"], 100)
# Env_sel2<-rep(envTab[,"envPop2"], 100)
# Env_nonsel1<-rep(envTab[,"fakeEnv1"], 100)
# Env_nonsel2<-rep(envTab[,"fakeEnv2"], 100)
# 
# cg_df$Env_sel1 <- Env_sel1
# cg_df$Env_sel2 <- Env_sel2
# cg_df$Env_nonsel1 <- Env_nonsel1
# cg_df$Env_nonsel2 <- Env_nonsel2


envTab <- Popsenv[,-1]
#envTab <- cbind(envPop1,envPop2,fakeEnv1,fakeEnv2)

#envTab <- cbind(envTab$envPop1/sd(envPop1$envPop1),envTab$envPop2/sd(envPop2$envPop2),envTab$fakeEnv1/sd(fakeEnv1),envTab$fakeEnv2/sd(fakeEnv2))
#colnames(envTab) <- c("envPop1","envPop2","fakeEnv1","fakeEnv2")
#envPop1.shift<-data.frame(fitt[fitt$Type=="Env1",gen_nam[length(gen_nam)]])
#names(envPop1.shift) <- "envSelect"

#Merge the population specific allele frequencies of all neutral (M1) alleles with the population specific frequency of the selected (M2) allele
#alFreq<-cbind(data.frame(Pop_afreq1))
alFreq<-Pop_afreq1

##############################################
# Chunk to fit GF models to minor allele frequencies at the level of
# populations
# GF is fit to each SNP individually to
# ease computational / memory burden

#cl <- makeCluster(cores)
#registerDoParallel(cl)

##### added by MCF, running all loci in on model ##########

vars <- colnames(envTab)


`getSplitImprove` <-function(fit, X) {
  #   return a data-frame: var name, rsq, var number, split value, improvement
  trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
  ok <- sapply(trees, class) != "try-error"
  ok2 <- apply(ok, 2, sum)==nrow(ok)
  tmp <- do.call("rbind", lapply((1:fit$ntree)[ok2], function(k) cbind(tree = k, trees[[k]])))
  tmp <- tmp[tmp[,"status"]==-3 & zapsmall(tmp[,"improve"]) > 0,c("split var","split point","improve")]
  colnames(tmp) <- c("var_n","split","improve")
  rownames(tmp)<-NULL     #S.J. Smith 11/05/2009
  res <- cbind(data.frame(var=names(X)[tmp[,"var_n"]],rsq=rep(fit$rsq[fit$ntree],nrow(tmp))),tmp)
  ok <- zapsmall(res[,"improve"]) > 0
  res[ok,] 
}

`getSplitImproveClassCompact` <- function(fit, bins, err0) {
  #   Return a data-frame: var name, rsq, split value, improvement
  #   Compact the splits into bins defined by bins matrix
  #   The i'th bin for predictor p is the interval (bin[i,p],bin[i+1,p])
  #   Every predictor is split into the same number of bins (nrow(bins)-1)
  
  #   extract all trees to a matrix and select for splits with some improvement
  trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
  ok <- sapply(trees, class) != "try-error"
  ok2 <- apply(ok, 2, sum)==nrow(ok)
  tmp <- do.call("rbind", lapply((1:fit$ntree)[ok2], function(k) cbind(tree = k, trees[[k]])))
  tmp <- tmp[tmp[,"status"]== 1 & zapsmall(tmp[,"improve"]) > 0,c("split var","split point","improve")]
  colnames(tmp) <- c("var_n","split","improve")
  rownames(tmp) <- NULL
  
  #   assign the split to the appropriate bin and aggregate importance in each bin
  Xnames <- colnames(bins)
  tmp <- data.frame(var=Xnames[tmp[,"var_n"]], tmp, bin=rep(0,nrow(tmp)))
  for(p in Xnames) {
    if(any(sub <- with(tmp,var==p)))
      tmp$bin[sub] <- as.numeric(cut(tmp$split[sub], bins[,p], include=TRUE, ordered=TRUE))
  }
  tmp <- with(tmp[tmp$bin>0,],agg.sum(improve,list(var,bin),sort.it=TRUE))
  names(tmp) <- c("var","bin","improve")
  
  #   Set the split value to the bin centre, but retain the bin number in case
  #   the bin centre is not appropriate value
  tmp <- cbind(tmp,split=rep(NA,nrow(tmp)),rsq=rep((err0-fit$err.rate[fit$ntree, "OOB"])/err0, nrow(tmp)))
  for(p in Xnames) {
    if(any(sub <- with(tmp,var==p)))
      tmp$split[sub] <- midpoints(bins[,p])[tmp$bin[sub]]
  }
  tmp[,c("var","rsq","split","improve","bin")]
}

`getSplitImproveClass` <- function(fit, X, err0)
  {
    #   return a data-frame: var name, rsq, var number, split value, improvement
    trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
    ok <- sapply(trees, class) != "try-error"
    ok2 <- apply(ok, 2, sum)==nrow(ok)
    tmp <- do.call("rbind", lapply((1:fit$ntree)[ok2], function(k) cbind(tree = k, trees[[k]])))
    tmp <- tmp[tmp[,"status"]==1,c("split var","split point","improve")]
    dimnames(tmp) <- list(NULL,c("var_n","split","improve"))
    res<-cbind(data.frame(var=names(X)[tmp[,"var_n"]],rsq=rep((err0-fit$err.rate[fit$ntree,"OOB"])/err0,nrow(tmp))),tmp)
    res
  }

`Impurity.based.measures` <-function(obj)
{
  #becomes an internal function not usually used by users
  #Modified 07/10/2009 by SJS re: NE changes for classification trees.
  dens <- lapply(names(obj$X), function(i) density(na.omit(obj$X[,i]),from=min(na.omit(obj$X[,i])),to=max(na.omit(obj$X[,i]))))
  dens <- lapply(dens,whiten,lambda=0.90) # hard-coded whitening
  names(dens) <- names(obj$X)
  res <- do.call("rbind", lapply(names(obj$result), function(spec) cbind(spec=spec,obj$result[[spec]]))) #added by Smith 13/05/2009
  res$spec <- as.factor(res$spec)
  res$var <- as.factor(res$var)
  res$improve <- pmax(0,res$improve)
  res$rsq <- pmax(0,res$rsq)   #added by Ellis 12/05/2009
  res$improve.tot <- tapply(res$improve,res$spec,sum)[res$spec]
  res$improve.tot.var <- tapply(res$improve,interaction(res$spec,res$var),sum)[interaction(res$spec,res$var)]
  res$improve.norm <- with(res,improve/improve.tot*rsq)
  nodup <- !duplicated(res[,1:2])
  res.u <- res[nodup, c("spec","var","rsq","improve.tot","improve.tot.var")]
  res.u$rsq <- with(res.u, ifelse(is.na(rsq), 0, rsq))
  res.u$rsq.var <- with(res.u,rsq*improve.tot.var/improve.tot)
  list(res=res,res.u=res.u,dens=dens)
}

`whiten` <-  function(dens, lambda=0.9)
{
  # add a small uniform value to the density to avoid zeroes when taking inverse
  dens$y <- lambda*dens$y + (1-lambda)/diff(range(dens$x))
  dens
}

gradientForest2<- function (data, predictor.vars, response.vars, ntree = 10, mtry = NULL, 
                            transform = NULL, maxLevel = 0, corr.threshold = 0.5, compact = FALSE, 
                            nbin = 101, trace = TRUE) 
{
  if (!inherits(data, "data.frame")) 
    stop("'data' must be a data.frame")
  X <- data[predictor.vars]
  Y <- data[response.vars]
  if (compact) {
    bins <- do.call("cbind", lapply(X, function(x) bin(x, 
                                                       nbin = nbin)))
  }
  if (!is.null(transform)) {
    Y <- apply(Y, 2, transform)
  }
  imp <- matrix(0, 0, 2, dimnames = list(NULL, c("%IncMSE", 
                                                 "IncNodePurity")))
  if (is.null(mtry)) 
    fitcmd <- quote(randomForest(Species ~ rhs, data = cbind(Y, 
                                                             X), maxLevel = maxLevel, keep.forest = TRUE, importance = TRUE, 
                                 ntree = ntree, keep.group = TRUE, keep.inbag = TRUE, 
                                 corr.threshold = corr.threshold, na.action = na.omit))
  else fitcmd <- quote(randomForest(Species ~ rhs, data = cbind(Y, 
                                                                X), maxLevel = maxLevel, keep.forest = TRUE, importance = TRUE, 
                                    ntree = ntree, mtry = mtry, keep.group = TRUE, keep.inbag = TRUE, 
                                    corr.threshold = corr.threshold, na.action = na.omit))
  result <- list()
  species.pos.rsq <- 0
  form.rhs <- as.formula(paste("Y ~ ", paste(predictor.vars, 
                                             collapse = "+")))
  if (trace) {
    spcount <- 0
    cat("Calculating forests for", length(response.vars), 
        "species\\n")
  }
  ##################################
  #Here is the problem
  ##################################
  for (spec in response.vars) {
    if (trace) 
      cat(if ((spcount <- spcount + 1)%%options("width")$width == 
              0) 
        "\\n."
        else ".")
    try({
      # spec <- response.vars[1]
      thisfitcmd <- do.call("substitute", list(fitcmd, 
                                               list(Species = as.name(spec), SpeciesName = spec, 
                                                    ntree = ntree, rhs = form.rhs[[3]])))
      fit <- eval(thisfitcmd)
      if (fit$type == "regression") {
        if (!is.na(fit$rsq[fit$ntree])) {
          if (fit$rsq[fit$ntree] > 0) {
            species.pos.rsq <- species.pos.rsq + 1
            if (compact) {
              result[[spec]] <- getSplitImproveCompact(fit, 
                                                       bins)
            }
            else {
              result[[spec]] <- getSplitImprove(fit, 
                                                X)
            }
            imp <- rbind(imp, fit$importance)
          }
        }
      }
      else if (fit$type == "classification") {
        if (!is.na(fit$err.rate[fit$ntree, "OOB"])) {
          p <- sum(Y[[spec]] == levels(Y[[spec]])[1])/length(Y[[spec]])
          err0 <- 2 * p * (1 - p)
          if (fit$err.rate[fit$ntree, "OOB"] < 
              2 * p * (1 - p)) {
            species.pos.rsq <- species.pos.rsq + 1
            if (compact) {
              result[[spec]] <- getSplitImproveClassCompact(fit, 
                                                            bins, err0)
            }
            else {
              result[[spec]] <- getSplitImproveClass(fit, 
                                                     X, err0)
            }
            nclass <- length(levels(Y[[spec]]))
            imp <- rbind(imp, fit$importance[, -(1:nclass)])
          }
        }
      }
      else stop(paste("unknown randomForest type:", 
                      fit$type))
    }, silent = FALSE)
  }
  if (!length(result)) {
    warning("No species models provided a positive R^2. \\nThe gradient forest is empty")
    return(NULL)
  }
  rsq <- sapply(result, function(x) x$rsq[1])
  imp.rsq <- matrix(imp[, 1], length(predictor.vars), dimnames = list(predictor.vars, 
                                                                      names(result)))
  imp.rsq[imp.rsq < 0] <- 0
  imp.rsq <- sweep(imp.rsq, 2, colSums(imp.rsq, na.rm = T), 
                   "/")
  imp.rsq <- sweep(imp.rsq, 2, rsq, "*")
  overall.imp <- tapply(imp[, 1], dimnames(imp)[[1]], mean, 
                        na.rm = T)
  overall.imp2 <- tapply(imp[, 2], dimnames(imp)[[1]], mean, 
                         na.rm = T)
  out1 <- list(X = X, Y = Y, result = result, overall.imp = overall.imp, 
               overall.imp2 = overall.imp2, ntree = ntree, imp.rsq = imp.rsq, 
               species.pos.rsq = species.pos.rsq, ranForest.type = fit$type)
  out2 <- Impurity.based.measures(out1)
  out1$result <- rsq
  out <- c(out1, out2, call = match.call())
  class(out) <- c("gradientForest", "list")
  out
}

# data(CoMLsimulation)
# preds <- colnames(Xsimulation)
# specs <- colnames(Ysimulation)
# start_time <- Sys.time()
# f1 <- gradientForest2(data.frame(Ysimulation,Xsimulation), preds, specs, ntree=10)
# end_time <- Sys.time()
# writeLines("GF genome calculation:")
# print(end_time - start_time)

start_time <- Sys.time()
gfMod_all <- gradientForest(data=data.frame(envTab[, vars], alFreq),
                        predictor.vars=vars,
                        response.vars=colnames(alFreq),
                        corr.threshold=0.5 ,
                        ntree=500,
                        trace=T)
end_time <- Sys.time()
writeLines("GF genome calculation:")
print(end_time - start_time)

GF_dist_gen <- NULL

for(i in 1:ncol(envTab)){
  tmp <- data.frame(envTab[,i])
  colnames(tmp) <- colnames(envTab)[i]
  
  writeLines(paste("Predicting ",colnames(envTab)[i]))
  gftmp <- predict(gfMod_all, tmp)
  
  GF_dist_gen <- data.frame(c(GF_dist_gen,gftmp))
}

# GF_dist_bef<-data.frame(cbind(gfTrans1e1,gfTrans1e2, gfTrans1e3, gfTrans1e4))
# 
# gfTrans1e1 <- predict(gfMod_all, envPop1)
# colnames(gfTrans1e1)<-"C.Imp_genome_envPop1"
# 
# gfTrans1e2 <- predict(gfMod_all, envPop2)
# colnames(gfTrans1e2)<-"C.Imp_genome_envPop2"
# 
# gfTrans1e3 <- predict(gfMod_all, data.frame(fakeEnv1))
# colnames(gfTrans1e3)<-"C.Imp_genome_fakeEnv1"
# 
# gfTrans1e4 <- predict(gfMod_all, data.frame(fakeEnv2))
# colnames(gfTrans1e4)<-"C.Imp_genome_fakeEnv2"
# 
# GF_dist_bef<-data.frame(cbind(gfTrans1e1,gfTrans1e2, gfTrans1e3, gfTrans1e4))
# colnames(GF_dist_bef)<-c("CI_env1", "CI_env2", "CI_fakeEnv1","CI_fakeEnv2")


#gfTrans1e1 <- predict(gfMod1, envPop1)
#colnames(gfTrans1e1)<-"C.Imp_genome_before"

# cImp1 <- cumimp(gfMod_all, "envPop1", type="Species")
# cImp1 <- data.frame(rbindlist(cImp1, idcol="allele"))
# 
# cImp2 <- cumimp(gfMod_all, "envPop2", type="Species")
# cImp2 <- data.frame(rbindlist(cImp2, idcol="allele"))
# 
# cImp3 <- cumimp(gfMod_all, "fakeEnv1", type="Species")
# cImp3 <- data.frame(rbindlist(cImp3, idcol="allele"))
# 
# cImp4 <- cumimp(gfMod_all, "fakeEnv2", type="Species")
# cImp4 <- data.frame(rbindlist(cImp4, idcol="allele"))


# ##### added by MCF, running all loci in on model ##########
# 
# gfMod2 <- gradientForest(data=data.frame(envPop2, alFreq),
#                          predictor.vars=colnames(envPop2),
#                          response.vars=colnames(alFreq),
#                          corr.threshold=0.5,
#                          ntree=500,
#                          trace=T)
# # Calculate genomic offset
# # note that I am doing this for the avearge across all alleles since
# # GF was fit to all alleles simultaneously
# # The more correct way is to calculate offset for adaptive alleles only,
# # either individually or for a model fit to just those alleles.

# offset needs to be considered using absolute values ()

alFreq_sel<-alFreq[colnames(alFreq)%in%select_pos]

# gfMod_sel1 <- gradientForest(data=data.frame(envPop1, alFreq_sel),
                        # predictor.vars=colnames(envPop1),
                        # response.vars=colnames(alFreq_sel),
                        # corr.threshold=0.5,
                        # ntree=500,
                        # trace=T)
start_time <- Sys.time()
gfMod_sel2 <- gradientForest(data=data.frame(envTab[, vars], alFreq_sel),
                        predictor.vars=vars,
                        response.vars=colnames(alFreq_sel),
                        corr.threshold=0.5 ,
                        ntree=500,
                        trace=T)
end_time <- Sys.time()
writeLines("GF causal calculation:")
print(end_time - start_time)

GF_dist_caus <- NULL

for(i in 1:ncol(envTab)){
  tmp <- data.frame(envTab[,i])
  colnames(tmp) <- colnames(envTab)[i]
  
  writeLines(paste("Predicting ",colnames(envTab)[i]))
  gftmp <- predict(gfMod_sel, tmp)
  
  GF_dist_caus <- data.frame(c(GF_dist_caus,gftmp))
}


# gfTrans1_sel1 <- predict(gfMod_sel, envPop1)
# colnames(gfTrans1_sel1)<-"C.Imp_causal_envPop1"
# 
# gfTrans1_sel2 <- predict(gfMod_sel, envPop2)
# colnames(gfTrans1_sel2)<-"C.Imp_causal_envPop2"
# 
# gfTrans1_sel3 <- predict(gfMod_sel, data.frame(fakeEnv1))
# colnames(gfTrans1_sel3)<-"C.Imp_causal_fakeEnv1"
# 
# gfTrans1_sel4 <- predict(gfMod_sel, data.frame(fakeEnv2))
# colnames(gfTrans1_sel4)<-"C.Imp_causal_fakeEnv2"

# cImp1_sel <- cumimp(gfMod_sel, "envPop1", type="Species")
# cImp1_sel <- data.frame(rbindlist(cImp1_sel, idcol="allele"))
# 
# cImp2_sel <- cumimp(gfMod_sel, "envPop2", type="Species")
# cImp2_sel <- data.frame(rbindlist(cImp2_sel, idcol="allele"))
# 
# cImp3_sel <- cumimp(gfMod_sel, "fakeEnv1", type="Species")
# cImp3_sel <- data.frame(rbindlist(cImp3_sel, idcol="allele"))
# 
# cImp4_sel <- cumimp(gfMod_sel, "fakeEnv2", type="Species")
# cImp4_sel <- data.frame(rbindlist(cImp4_sel, idcol="allele"))



##############################################################################
#Get  Weir & Cockerham F_ST values from the VCF files and use population data

#Create an object listing every population in the whole dataset
PopsALL <- NULL
for(j in rep(1:100)){
  for(i in rep(j,10)){
    PopsALL <- c(PopsALL,i)
  }
}

# #Create an object splitting a single population into a Pre ("T1") and Post ("T2") "population"
# PopsP <- c(rep("T1",10),rep("T2",10))
# 
# # cores<-3
# # cl <- makeCluster(cores)
# # registerDoParallel(cl)
# 
# #Create unfiltered Genotype matrix
# G1 <- matrix(NA, nrow = nrow(geno1), ncol = ncol(geno1))
# G1[geno1 %in% c("0/0", "0|0")] <- 0
# G1[geno1  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
# G1[geno1 %in% c("1/1", "1|1")] <- 2
# 
# 
# G2 <- matrix(NA, nrow = nrow(geno2), ncol = ncol(geno2))
# G2[geno2 %in% c("0/0", "0|0")] <- 0
# G2[geno2  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
# G2[geno2 %in% c("1/1", "1|1")] <- 2
# 
# #Calculate allele frequencies across the whole meta population
# a_freq1 <- rowSums(G1)/(2*ncol(G1))
# a_freq2 <- rowSums(G2)/(2*ncol(G2))
# 
# #Prepare the Pre (G1) and Post (G2) G matrices for FST calculation
# Gt1<-t(G1)
# rownames(Gt1)<-PopsALL
# #colnames(Gt1)<-MID_pre
# colnames(Gt1)<-paste("M",position1,sep="")
# 
# Gt2<-t(G2)
# rownames(Gt2)<-PopsALL
# #colnames(Gt2)<-MID_post
# colnames(Gt2)<-paste("M",position2,sep="")
# 
# #Add loop to iterate across each TRUE population, not just x-location
# #Filter by population (same pop. before and after env. shift)
# listGt<-list()
# for(i in 1:100){
#   Gt1_i<-data.frame(Gt1[rownames(Gt1)==i,])
#   Gt2_i<-data.frame(Gt2[rownames(Gt2)==i,])
#   Gt1_i[setdiff(names(Gt2_i), names(Gt1_i))] <- 0
#   Gt2_i[setdiff(names(Gt1_i), names(Gt2_i))] <- 0
#   listGt[[i]]<-rbind(Gt1_i,Gt2_i)
# }
# 
# #Filter the files based on Major and Minor AF filtration, then seed the files so that all variants are present at both time points to be compared
# a_freq<-list()
# a_freq_filt<-list()
# listGt_filt<-list()
# for(i in 1:100){
#   a_freq[[i]] <- colSums(listGt[[i]])/(2*nrow(listGt[[i]]))
#   a_freq_filt[[i]] <- a_freq[[i]][a_freq[[i]]>0.05 & a_freq[[i]]<0.95]
#   listGt_filt[[i]]<-listGt[[i]][colnames(listGt[[i]])%in%names(a_freq_filt[[i]])]
# }
# 
# #Calculate per locus FST values
# listPfst<-list()
# for(i in 1:100){
#   listPfst[[i]]<-MakeDiploidFSTMat(SNPmat = listGt_filt[[i]], locusNames = colnames(listGt_filt[[i]]), popNames = PopsP)
# }
# 
# #Filter out NA values
# listPfst_noNa<-list()
# for(i in 1:100){
#   listPfst_noNa<-lapply(listPfst,function(x) x[!is.na(x$FST),])
# }
# 
# #Calculate FST values for each populations
# FST_genome_pop<-NULL
# for(i in 1:100){
#   FST_genome_pop<-c(FST_genome_pop,mean(listPfst_noNa[[i]]$T1)/mean(listPfst_noNa[[i]]$T2))
# }
# 
getFSTs_diploids = function(popNameList, SNPDataColumn){
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]

  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0

  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }

  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)}
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }

  out = WC_FST_Diploids_2Alleles(HetCounts)
  return(out)
}

MakeDiploidFSTMat_2<-function(SNPmat,locusNames,popNames){
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  if(any(!(snplevs%in%c(0,1,2,9)))==TRUE) {
    print("Error: Your snp matrix has a character other than 0,1,2 or 9")
    break
  }
  if (dim(SNPmat)[1] != length(popname)) {
    print("Error: your population names do not match your SNP matrix")
    break
  }
  if (dim(SNPmat)[2] != length(locusname)) {
    print("Error:  your locus names do not match your SNP matrix")
    break
  }
  writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
  for (i in 1:nloci) {
    FSTmat[i, ] = unlist(getFSTs_diploids(popname, SNPmat[,i]))
    if (i%%10000 == 0) {
      print(paste(i, "done of", nloci))
    }
  }
  outTemp = as.data.frame(FSTmat)
  outTemp = cbind(locusname, outTemp)
  colnames(outTemp) = c("LocusName", "He", "FST", "T1", "T2",
                        "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return(outTemp)
}

########################################################
#Calculate per locus FST values pre environmental shift for Edge and cCore populations

#Whole genome
Gt1f<-t(G1f)
colnames(Gt1f)<-paste("M",position1_filt,sep="")

Pre_geno <- data.frame(PopsALL,Gt1f)
colnames(Pre_geno)[1]<-"Locality"

Edges <- c(1,5,6,10,12,19,41,50,51,60,82,89,91,95,96,100)

Pre_geno_Edge <- Pre_geno[Pre_geno$Locality%in%Edges,]

start_time <- Sys.time()
Pre_FST_Edge<-pairwise.WCfst(Pre_geno_Edge,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_Edge, na.rm=T)

Cores <- c(34, 35, 36, 37, 44, 45, 46, 47, 54, 55, 56, 57, 64, 65, 66, 67)

Pre_geno_Cores <- Pre_geno[Pre_geno$Locality%in%Cores,]

start_time <- Sys.time()
Pre_FST_Cores<-pairwise.WCfst(Pre_geno_Cores,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_Cores, na.rm=T)

#Genome under selection
Gt1f_sel<-Gt1f[,colnames(Gt1f)%in%select_pos]
Pre_geno_sel<-data.frame(PopsALL,Gt1f_sel)
colnames(Pre_geno_sel)[1]<-"Locality"

Pre_geno_sel_Edge <- Pre_geno_sel[Pre_geno_sel$Locality%in%Edges,]

start_time <- Sys.time()
Pre_FST_sel_Edge<-pairwise.WCfst(Pre_geno_sel_Edge,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_sel_Edge, na.rm=T)

Cores <- c(34, 35, 36, 37, 44, 45, 46, 47, 54, 55, 56, 57, 64, 65, 66, 67)

Pre_geno_sel_Cores <- Pre_geno_sel[Pre_geno_sel$Locality%in%Cores,]

start_time <- Sys.time()
Pre_FST_sel_Cores<-pairwise.WCfst(Pre_geno_sel_Cores,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_sel_Cores, na.rm=T)
# # Post_geno<-data.frame(PopsALL,Gt2f)
# # colnames(Post_geno)[1]<-"Locality"
# 
# 

# start_time <- Sys.time()
# Pre_FST<-pairwise.WCfst(Pre_geno,diploid = T)
# end_time <- Sys.time()
# print(end_time - start_time)

# 
# # start_time <- Sys.time()
# # Post_FST<-pairwise.WCfst(Post_geno,diploid = T)
# # end_time <- Sys.time()
# # print(end_time - start_time)
# #

#write.table(Pre_FST,file=paste("results/R_results/",seed,"_FST.csv",sep=""),sep=",",col.names = F,row.names = F,quote = F)

# # write.table(Post_FST,file="~/Desktop/Post_FST_2.txt",sep=",",col.names = F,row.names = F,quote = F)
# 
# Gt1f_sel<-Gt1f[,colnames(Gt1f)%in%select_pos]
# Pre_geno_sel<-data.frame(PopsALL,Gt1f_sel)
# colnames(Pre_geno_sel)[1]<-"Locality"
# 
# start_time <- Sys.time()
# Pre_FST_sel<-pairwise.WCfst(Pre_geno_sel,diploid = T)
# end_time <- Sys.time()
# print(end_time - start_time)

#write.table(Pre_FST,file=paste("results/R_results/",seed,"_FST_sel.csv",sep=""),sep=",",col.names = F,row.names = F,quote = F)
# ########################################################

#Get per population FST pre and post environmental shift
Gt1fd <- data.frame(Gt1f)
Pfst_pre_filt<-MakeDiploidFSTMat_2(SNPmat = Gt1fd, locusNames = colnames(Gt1fd), popNames = PopsALL)

#Filter out NA values
Pfst_pre_noNa<-Pfst_pre_filt[!is.na(Pfst_pre_filt$FST),]

#ink_bef<-data.frame(cor(Gt1_m2,Gt1f[,colnames(Gt1f)!=M2_MID],method="pearson"))

#Heterozgosity per allele before env. shift
Het_bef<-Pfst_pre_filt$He

#Calculate FST value for each allele pre environmental shift
F_ST_ll1<-Pfst_pre_noNa$T1/Pfst_pre_noNa$T2

#Calculate FST values averaged across each allele pre environmental shift
F_ST_l1<-mean(Pfst_pre_noNa$T1)/mean(Pfst_pre_noNa$T2)

# Gt2f<-t(G2f)
# rownames(Gt2f)<-PopsALL
# #colnames(Gt2f)<-MID_post_filt
# colnames(Gt2f)<-paste("M",position2_filt,sep="")
# #Calculate per locus FST values post environmental shift
# Pfst_post<-MakeDiploidFSTMat(SNPmat = Gt2f, locusNames = colnames(Gt2f), popNames = PopsALL)
# 
# #Filter out NA values
# Pfst_post_noNa<-Pfst_post[!is.na(Pfst_post$FST),]
# 
# #Heterozgosity per allele after env. shift
# Het_aft<-Pfst_post_noNa$He
# 
# #Calculate FST value for each allele pre environmental shift
# F_ST_ll2<-Pfst_post_noNa$T1/Pfst_post_noNa$T2
# 
# #Calculate FST values for each  allele pre environmental shift
# F_ST_l2<-mean(Pfst_post_noNa$T1)/mean(Pfst_post_noNa$T2)
# 
# #Pop_afreq2 is not filtered for MAF in order to properly compare all AF shifts from Pop_afreq1
# Pop_afreq2<-NULL
# for(i in 1:100){
#   Pop_afreq2<-rbind(Pop_afreq2,rowSums(G2[,Start[i]:Stop[i]])/(2*ncol(G2[,Start[i]:Stop[i]])))
# }
# 
# Pop_afreq2<-data.frame(Pop_afreq2)
# colnames(Pop_afreq2)<-paste("M",position2,sep="")

#stopCluster(cl)

##################################
#Population specific summary stats
##################################

#Location values
Loc <- NULL
for(j in 1:10){
  for(i in 1:10){
    Loc <- c(Loc,paste("A",i,sep=""))
  }
}
Loc<-factor(Loc,levels=Loc[1:10])

#Population values
Pop <- NULL
for(i in 1:100){
  Pop <- c(Pop,paste("P",i,sep=""))
}

X <- NULL
for(j in 1:10){
  for(i in 1:10){
    X <- c(X,i)
  }
}

Y <- NULL
for(j in 1:10){
  Y<-c(Y,rep(j,10))
}

# Env1<-envPop1$envPop1
# Env2<-envPop2$envPop2

#Env_after<-envPop.shift$envSelect

#Diff_env<-Env_after-Env_before

#Env_range<- envPop$envSelect%in%round(envPop.shift$envSelect,1)&round(envPop.shift$envSelect,1)%in%envPop$envSelect

# M2_AF_before<-data.frame(Pop_afreq1[,c(which(colnames(Pop_afreq1)==positionM2))])
# colnames(M2_AF_before)<-"M2_AF_before"
#
# M2_AF_after<-data.frame(Pop_afreq2[,c(which(colnames(Pop_afreq2)==positionM2))])
# colnames(M2_AF_after)<-"M2_AF_after"
#
# M2_AF_diff<-M2_AF_after-M2_AF_before
# colnames(M2_AF_diff)<-"M2_AF_diff"

# M1_AF_before_all<-Pop_afreq1[,-which(colnames(Pop_afreq1)==positionM2)]
# M1_AF_after_all<-Pop_afreq2[,-which(colnames(Pop_afreq2)==positionM2)]
#
# M1_AF_before_shared<-M1_AF_before_all[colnames(M1_AF_before_all)%in%colnames(M1_AF_after_all)]
# M1_AF_after_shared<-M1_AF_after_all[colnames(M1_AF_after_all)%in%colnames(M1_AF_before_all)]

AF_all<-data.frame(rowMeans(Pop_afreq1))
# colnames(M1_AF_before)<-"M1_AF_before"
# 
# M1_AF_after<-data.frame(rowMeans(M1_AF_after_shared))
# colnames(M1_AF_after)<-"M1_AF_after"
# 
# M1_AF_diff<-M1_AF_after-M1_AF_before
# colnames(M1_AF_diff)<-"M1_AF_diff"

# F_ST_genome_bef.aft.<-data.frame(FST_genome_pop)
# colnames(F_ST_genome_bef.aft.)<-"F_ST_genome_bef.aft."

# F_ST_M2_bef.aft.<-data.frame(FST_M2_pop)
# colnames(F_ST_M2_bef.aft.)<-"F_ST_M2_bef.aft."

Rel_Fit <- data.frame(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)-30]])
colnames(Rel_Fit)<-"Rel_Fit"

# Rel_Fit_after <- data.frame(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
# colnames(Rel_Fit_after)<-"Rel_Fit_after"
# 
# Rel_Fit_diff<-Rel_Fit_after-Rel_Fit_before
# colnames(Rel_Fit_diff)<-"Rel_Fit_diff"

Pop_size<-rep(specs$n,100)
# Summary_Pop<-cbind(Pop,X,Y,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1,gfM2Trans2,M2offset,M2_AF_before,M2_AF_after,M2_AF_diff,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)

#If neutral:
# Pop_size<-unlist(rep(unname(specs[9:18]),10))
# Summary_Pop<-cbind(Pop,X,Y,Pop_size,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1=0,gfM2Trans2=0,M2offset=0,M2_AF_before=0,M2_AF_after=0,M2_AF_diff=0,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.=0,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)

Summary_Pop<-cbind(Pop,X,Y,envTab, GF_dist_gen, AF_all, Rel_Fit)
write.csv(Summary_Pop,file=paste("results/R_results/",seed,"_summary_Pop.csv",sep=""),row.names=F)
#Summary_Pop<-read.csv("Summary_Pop_1576675870126.csv")

##################################
# Allele specific summary stats
##################################

R2_genome<-data.frame(gfMod_all$result) #Those allelese which had an R2 value > 0
colnames(R2_genome)<-"R2"

#R2_e2_genome<-data.frame(gfMod2$result) #Those allelese which had an R2 value > 0
#colnames(R2_e2_genome)<-"R2_e2_genome"

R2_sel<-data.frame(gfMod_sel$result) #Those allelese which had an R2 value > 0
colnames(R2_sel)<-"R2_adapt"

# R2_e2_sel<-data.frame(gfMod_sel2$result) #Those allelese which had an R2 value > 0
# colnames(R2_e2_sel)<-"R2_e2_genome"


R0<-colnames(alFreq[,!colnames(alFreq)%in%names(gfMod_all$result)]) #Get all alleles, regardless of R2 value
R0<-data.frame(rep(0,length(R0)),row.names = colnames(alFreq[,!colnames(alFreq)%in%names(gfMod_all$result)])) #Filter out those who we already have saved in R2
colnames(R0)<-"R2"
R2_all<-rbind(R2_genome,R0) #merge them so we have all alleles with accompanying R2 values

R2_all$Pos <- substring(row.names(R2_all),2)

Rho_Env1<-cor(as.matrix(alFreq),envPop1,method = "pearson")
Rho_Env2<-cor(as.matrix(alFreq),envPop2,method = "pearson")

temp<-merge(Rho_Env1, Rho_Env2, by=0)
colnames(temp)<-c("Locus","Rho_Env1","Rho_Env2")
row.names(temp)<-temp$Locus
temp<-merge(temp, R2_all, by=0)
temp <- temp[,-1]
colnames(temp)<-c("Locus","Rho_Env1","Rho_Env2","R2","Position")

Linkage<-NULL
position1_filt_scaled<-NULL
for(i in 1:length(position1_filt)){
  if(position1_filt[i]>=0 & position1_filt[i]<50001){
    Linkage<-c(Linkage,1)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i])
  }
  if(position1_filt[i]>50000 & position1_filt[i]<100001){
    Linkage<-c(Linkage,2)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-50000)
  }
  if(position1_filt[i]>100000 & position1_filt[i]<150001){
    Linkage<-c(Linkage,3)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-100000)
  }
  if(position1_filt[i]>150000 & position1_filt[i]<200001){
    Linkage<-c(Linkage,4)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-150000)
  }
  if(position1_filt[i]>200000 & position1_filt[i]<250001){
    Linkage<-c(Linkage,5)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-200000)
  }
  if(position1_filt[i]>250000 & position1_filt[i]<300001){
    Linkage<-c(Linkage,6)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-250000)
  }
  if(position1_filt[i]>300000 & position1_filt[i]<350001){
    Linkage<-c(Linkage,7)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-300000)
  }
  if(position1_filt[i]>350000 & position1_filt[i]<400001){
    Linkage<-c(Linkage,8)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-350000)
  }
  if(position1_filt[i]>400000 & position1_filt[i]<450001){
    Linkage<-c(Linkage,9)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-400000)
  }
  if(position1_filt[i]>450000 & position1_filt[i]<500001){
    Linkage<-c(Linkage,10)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-450000)
  }
}

#PreN$LG<-Linkage
#PreN[rownames(PreN)!=Pfst_pre_filt$LocusName,]
#PreN$FST<-Pfst_pre_filt$FST
#PreN$DistM2<-abs(PreN$PP-PreN$PP[PreN$LT=="M2"])
#temp$Position[300]=113353.51

Summary_Locus<-cbind(temp[order(as.numeric(temp$Position)),],Linkage, F_ST_ll1)

write.csv(Summary_Locus,file=paste("results/R_results/",seed,"_summary_Loc.csv",sep=""),row.names=F)

##################################
# Simulation specific summary stats
##################################

#No_A<-unname(dim(vcf1_filt)[1])
PR2<-gfMod_all$species.pos.rsq/No_A
#PR2_UL
F_ST_l1
#F_ST_l2
Rel.Fit_l1<-mean(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)-30]])
#Rel.Fit_l2<-mean(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
Rho_Env1R2<-cor(temp$Rho_Env1,temp$R2,method = "pearson", use="complete.obs")
Rho_Env2R2<-cor(temp$Rho_Env2,temp$R2,method = "pearson", use="complete.obs")

Summary_Sim<-cbind(seed,No_A,PR2,F_ST_l1,Rel.Fit_l1,Rho_Env1R2,Rho_Env2R2)

write.csv(Summary_Sim,file=paste("results/R_results/",seed,"_summary_Sim.csv",sep=""),row.names=F)

###############################################
#Visualize R^2>0 compared to Spearman correltation of alFreq to each env.
###############################################

#EnvCor<-cor(as.matrix(R2),envPop,method = "spearman")
#EnvCor<-data.frame(rownames(EnvCor),EnvCor)
#colnames(EnvCor)<-c("MID","rho")

#R2MID<-unique(data.frame(impDat$allele,impDat$r2))
#colnames(R2MID)<-c("MID","R2")

#Comp<-merge(EnvCor,R2MID, by="MID")

#Link<-NULL
#for(i in 1:length(Comp$MID)){
#  if(Comp$MID[i]%in%linked_MID){
#    Link<-c(Link,"Linked")
#    }
#  else{
#    Link<-c(Link,"Unlinked")
#  }
#}

#Comp$Link<-Link

#r2<-data.frame(gfMod$result)
#colnames(r2)<-"r2"

####################################################
#Get difference in Cumulative Importance (GF Offset)
####################################################

CI<-NULL
for(i in 1:nrow(GF_dist_gen)){
  for(j in 1:nrow(GF_dist_gen)){
    CI<-c(CI,dist(rbind(GF_dist_gen[i,],GF_dist_gen[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_allEnv<-CI

GFcaus_dist_gen<-GF_dist_gen[,c(1,2)]
CI_causal<-NULL
for(i in 1:nrow(GFcaus_dist_gen)){
  for(j in 1:nrow(GFcaus_dist_gen)){
    CI_causal<-c(CI_causal,dist(rbind(GFcaus_dist_gen[i,],GFcaus_dist_gen[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_causEnv<-CI_causal


CI_sel<-NULL
for(i in 1:nrow(GF_dist_caus)){
  for(j in 1:nrow(GF_dist_caus)){
    CI_sel<-c(CI_sel,dist(rbind(GF_dist_caus[i,],GF_dist_caus[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_caus_allEnv<-CI_sel

GFcaus_dist_caus<-GF_dist_caus[,c(1,2)]
CIcaus_causal<-NULL
for(i in 1:nrow(GFcaus_dist_caus)){
  for(j in 1:nrow(GFcaus_dist_caus)){
    CIcaus_causal<-c(CIcaus_causal,dist(rbind(GFcaus_dist_caus[i,],GFcaus_dist_caus[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_caus_causEnv<-CIcaus_causal


options(scipen = 0)

dat_sum<-cor.test(x=cg_df$D_CI_allEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$D_CI_allEnv, y=cg_df$Fitness)
paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$D_CI_causEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$D_CI_causEnv, y=cg_df$Fitness)
paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$D_CI_caus_allEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$D_CI_caus_allEnv, y=cg_df$Fitness)
paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$D_CI_caus_causEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$D_CI_caus_causEnv, y=cg_df$Fitness)
paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$EdSelEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$EdSelEnv, y=cg_df$Fitness)
paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$MdSelEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$MdSelEnv, y=cg_df$Fitness)
paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$EdAllEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$EdAllEnv, y=cg_df$Fitness)
paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$MdAllEnv, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$MdAllEnv, y=cg_df$Fitness)
paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$EdSelEnvPlus2, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$MdAllEnv, y=cg_df$Fitness)
paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

dat_sum<-cor.test(x=cg_df$MdSelEnvPlus2, y=cg_df$Fitness, method = "spearman")
plot(x=cg_df$MdSelEnvPlus2, y=cg_df$Fitness)
paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")

Pre_FST_Cores[is.na(Pre_FST_Cores)]<-0
Pre_FST_Edge[is.na(Pre_FST_Edge)]<-0
Pre_FST_sel_Cores[is.na(Pre_FST_sel_Cores)]<-0
Pre_FST_sel_Edge[is.na(Pre_FST_sel_Edge)]<-0

cg_df_Edges <- cg_df[as.numeric(substr(as.character(cg_df$Home),2,nchar(as.character(cg_df$Home))))%in%c(Edges,Cores),]
cg_df_Edges <- cg_df_Edges[as.numeric(substr(as.character(cg_df_Edges$Transplant),2,nchar(as.character(cg_df_Edges$Transplant))))%in%Edges,]

cg_df_Edges$Pre_FST <- as.vector(Pre_FST_Edge)
cg_df_Edges$Pre_FSt_sel <- as.vector(Pre_FST_sel_Edge)

cg_df_Cores <- cg_df[as.numeric(substr(as.character(cg_df$Home),2,nchar(as.character(cg_df$Home))))%in%c(Edges,Cores),]
cg_df_Cores <- cg_df_Cores[as.numeric(substr(as.character(cg_df_Cores$Transplant),2,nchar(as.character(cg_df_Cores$Transplant))))%in%Cores,]

cg_df_Cores$Pre_FST <- as.vector(Pre_FST_Cores)
cg_df_Cores$Pre_FSt_sel <- as.vector(Pre_FST_sel_Cores)

#compare edge pops
dat_sum<-cor.test(x=cg_df_Edges$D_CI_allEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
#summary(dat_sum)
GF_off_genome_allEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$D_CI_causEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_causEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$D_CI_caus_allEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_allEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$D_CI_caus_causEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_causEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$EdSelEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_caus_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$MdSelEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_caus_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$EdAllEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_all_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$MdAllEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_all_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$EdSelEnvPlus2, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_causP2_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$MdSelEnvPlus2, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_causP2_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$Pre_FST, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_genome_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$Pre_FSt_sel, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_causal_Edges <- round(dat_sum$estimate[[1]],3)

#compare core pops
dat_sum<-cor.test(x=cg_df_Cores$D_CI_allEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_allEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI_causEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_causEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI_caus_allEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_allEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI_caus_causEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_causEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$EdSelEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_caus_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$MdSelEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_caus_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$EdAllEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_all_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$MdAllEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_all_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$EdSelEnvPlus2, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_causP2_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$MdSelEnvPlus2, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_causP2_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$Pre_FST, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_genome_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$Pre_FSt_sel, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_causal_Cores <- round(dat_sum$estimate[[1]],3)

case <- data.frame(GF_off_genome_allEnv_Edges, 
                   GF_off_genome_causEnv_Edges, 
                   GF_off_causal_allEnv_Edges, 
                   GF_off_causal_causEnv_Edges, 
                   Env_ED_caus_env_Edges,
                   Env_MD_caus_env_Edges,
                   Env_ED_all_env_Edges,
                   Env_MD_all_env_Edges,
                   Env_ED_causP2_env_Edges,
                   Env_MD_causP2_env_Edges, 
                   FST_genome_Edges, 
                   FST_causal_Edges,
                   GF_off_genome_allEnv_Cores, 
                   GF_off_genome_causEnv_Cores, 
                   GF_off_causal_allEnv_Cores, 
                   GF_off_causal_causEnv_Cores, 
                   Env_ED_caus_env_Cores,
                   Env_MD_caus_env_Cores,
                   Env_ED_all_env_Cores,
                   Env_MD_all_env_Cores,
                   Env_ED_causP2_env_Cores,
                   Env_MD_causP2_env_Cores, 
                   FST_genome_Cores, 
                   FST_causal_Cores)

write.csv(case, paste("results/R_results/Case",seed_table$V5[w],"_",seed,"_corr_allEnv.csv",sep=""),row.names = F)

save.image(paste("results/R_results/",seed,".RData",sep=""))


rm(list=ls())
gc()
}
#stopCluster(cl)
