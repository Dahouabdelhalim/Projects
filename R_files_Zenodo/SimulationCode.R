###Code to accompany "The role of behavioural flexibility in primate diversification"
##Maria Creighton, email: maria.creighton@mail.mcgill.ca
rm(list=ls())

#Set your working directory
setwd("")



####1.0 Technical Innovation and Social Learning Lineage-Richness Sampling Bias####
LineageData <- read.csv("LineageData.csv") #Read in lineage data

library(dplyr)
library(ape)
tree <- read.tree("consensus_tree_10kTrees") #Read in tree
LineageData <- LineageData %>%
  filter(!Genus %in% c("Presbytis", "Galagoides", "Otolemur", "Galago", "Euoticus", "Semnopithecus")) #Non-monophyletic genera (see Methods)
tree<- drop.tip(tree,setdiff(tree$tip.label,LineageData$Trees_Name)) #Drop genera with unreliable node ages from lineage data

#Genus data
genus_data_p<- read.csv("GenusData.csv") #Read in genus data

library(ape)
genus<- distinct(LineageData,Genus,.keep_all=TRUE) #This makes each genus have one row
tips<- as.character(genus$Trees_Name) #extract the unique names that correspond to the tree
g_tree<- drop.tip(tree,setdiff(tree$tip.label,tips)) #Drop off the excess tips so you have one per genus
g_tree$tip.label<- gsub('\\\\_.*','',g_tree$tip.label) #Rename tip labels to the Genus name (this is purely aesthetic)

g_tree<- drop.tip(g_tree,setdiff(g_tree$tip.label,genus_data_p$Genus)) #Drop genera with unreliable node ages from tree
m <- match(genus_data_p$Genus,g_tree$tip.label) #Match gensus to tip labels
summary(m)

genus_data_p<- genus_data_p[m,] #Put data in order in appears on tree

#PGLS
library(nlme)
genus_data_p<- as.data.frame(genus_data_p)
rownames(genus_data_p)<- genus_data_p$Genus

sl<- subset(genus_data_p,is.na(Genus_Social_Learning)==F)
tree_sl<- drop.tip(g_tree,setdiff(g_tree$tip.label,sl$Genus))
l1.gr<-gls(G_ER~as.factor(Genus_Social_Learning),correlation=corPagel(1,form = ~Genus,phy=tree_sl,fixed=F),data=sl,method='ML')
l1.glr<- gls(G_LER~as.factor(Genus_Social_Learning),correlation=corPagel(0,form = ~Genus,phy=tree_sl,fixed=T),data=sl,method='ML')
l1.gr.p<-gls(G_ER~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=genus_data_p,method='ML')
l1.glr.p<- gls(G_LER~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=genus_data_p,method='ML')


ti<- subset(genus_data_p,is.na(Genus_Technical_Innovation)==F)
tree_ti<- drop.tip(g_tree,setdiff(g_tree$tip.label,ti$Genus))
l2.gr<-gls(G_ER~as.factor(Genus_Technical_Innovation),correlation=corPagel(1,form = ~Genus,phy=tree_ti,fixed=F),data=ti,method='ML')
l2.glr<- gls(G_LER~as.factor(Genus_Technical_Innovation),correlation=corPagel(0,form = ~Genus,phy=tree_ti,fixed=T),data=ti,method='ML')
l2.gr.p<-gls(G_ER~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=genus_data_p,method='ML')
l2.glr.p<- gls(G_LER~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=genus_data_p,method='ML')


tisl<- subset(genus_data_p,is.na(G_TISL)==F)
tree_tisl<- drop.tip(g_tree,setdiff(g_tree$tip.label,tisl$Genus))
l3.gr<-gls(G_ER~as.factor(G_TISL),correlation=corPagel(1,form = ~Genus,phy=tree_tisl,fixed=F),data=tisl,method='ML')
l3.glr<- gls(G_LER~as.factor(G_TISL),correlation=corPagel(0,form = ~Genus,phy=tree_tisl,fixed=T),data=tisl,method='ML')
l3.gr.p<-gls(G_ER~as.factor(Imputed_G_TISL),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=genus_data_p,method='ML')
l3.glr.p<- gls(G_LER~as.factor(Imputed_G_TISL),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=genus_data_p,method='ML')


library(geiger); library(nlme)
SL<- LineageData$Imputed_Social_Learning_Binary
names(SL)<- LineageData$Trees_Name
TI<- LineageData$Imputed_Technical_Innovation_Binary
names(TI)<- LineageData$Trees_Name
TISL<- LineageData$Imputed_TISL
names(TISL)<- LineageData$Trees_Name

SL<- SL[complete.cases(SL)]
TI<- TI[complete.cases(TI)]
TISL<- TISL[complete.cases(TISL)]

SL<- as.factor(SL)
TI<- as.factor(TI)
TISL<- as.factor(TISL)

tree_SL<- drop.tip(tree,setdiff(tree$tip.label,names(SL)))
tree_TI<- drop.tip(tree,setdiff(tree$tip.label,names(TI)))
tree_TISL<- drop.tip(tree,setdiff(tree$tip.label,names(TISL)))

tree_SL<- multi2di(tree_SL, random = TRUE) #Resolve polytomy
tree_TI<- multi2di(tree_TI, random = TRUE)
tree_TISL<- multi2di(tree_TISL, random = TRUE)
tree.p <- multi2di(tree, random = TRUE)

#Symmetrical models
M.sym.sl<- fitDiscrete(tree_SL,SL,model='SYM',transform=c('none')) #Fit symmetric Mk models
M.sym.ti<- fitDiscrete(tree_TI,TI,model='SYM',transform=c('none')) 
M.sym.tisl<- fitDiscrete(tree_TISL,TISL,model='SYM',transform=c('none'))

par_rates_sl<- matrix(ncol=2,nrow=2) #recreating q-matrix
par_rates_sl[1,2]<- M.sym.sl$opt$q12
par_rates_sl[2,1]<- M.sym.sl$opt$q21
par_rates_sl[1,1]<- -sum(par_rates_sl[1,2])
par_rates_sl[2,2]<- -sum(par_rates_sl[2,1])

par_rates_ti<- matrix(ncol=2,nrow=2)
par_rates_ti[1,2]<- M.sym.ti$opt$q12
par_rates_ti[2,1]<- M.sym.ti$opt$q21
par_rates_ti[1,1]<- -sum(par_rates_ti[1,2])
par_rates_ti[2,2]<- -sum(par_rates_ti[2,1])

par_rates_tisl<- matrix(ncol=2,nrow=2)
par_rates_tisl[1,2]<- M.sym.tisl$opt$q12
par_rates_tisl[2,1]<- M.sym.tisl$opt$q21
par_rates_tisl[1,1]<- -sum(par_rates_tisl[1,2])
par_rates_tisl[2,2]<- -sum(par_rates_tisl[2,1])

sim.sl<- sim.char(tree.p,par=par_rates_sl,model='discrete',root=1,nsim=1000)
sim.ti<- sim.char(tree.p,par=par_rates_ti,model='discrete',root=1,nsim=1000)
sim.tisl<- sim.char(tree.p,par=par_rates_tisl,model='discrete',root=1,nsim=1000)

b_sl.gr.p<- NA
b_sl.glr.p<- NA

pgl_sl.gr.p<- NA
pgl_sl.glr.p<- NA

b_ti.gr.p<- NA
b_ti.glr.p<- NA

pgl_ti.gr.p<- NA
pgl_ti.glr.p<- NA

b_tisl.gr.p<- NA
b_tisl.glr.p<- NA

b_pgls_sl.gr.p<- NA
b_pgls_sl.glr.p<- NA
b_pgls_ti.gr.p<- NA
b_pgls_ti.glr.p<- NA
b_pgls_tisl.gr.p<- NA
b_pgls_tisl.glr.p<- NA

n.SL.p<- NA
n.TI.p<- NA
n.TISL.p<- NA

for(i in 1:1000){
  m<- match(LineageData$Trees_Name,rownames(sim.sl))
  test<- data.frame(Genus=LineageData$Genus,Richness=LineageData$Richness,SL.p=sim.sl[m,,i],TI.p=sim.ti[m,,i],TISL.p=sim.tisl[m,,i])
  
  test_g<- group_by(test,Genus)
  test_g<- summarize(test_g,Genus_Richness=sum(Richness), Genus_Lineage_Richness=n(), Imputed_Genus_Technical_Innovation_Binary=max(TI.p), Imputed_Genus_Social_Learning_Binary=max(SL.p),Imputed_G_TISL=max(TISL.p))#Summarizing the data at Genus level
  test_g<- test_g %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), NA, x))
  test_g<- as.data.frame(test_g)
  rownames(test_g)<- test_g$Genus
  m<- match(test_g$Genus,genus_data_p$Genus)
  test_g[,7:8]<- genus_data_p[m,5:6]
  n.SL.p[i]<- unname(summary(as.factor(test$SL.p))[2])
  if(is.na(n.SL.p[i])==TRUE){next}
  n.TI.p[i]<-  unname(summary(as.factor(test$TI.p))[2])
  if(is.na(n.TI.p[i])==TRUE){next}
  n.TISL.p[i]<-  unname(summary(as.factor(test$TISL.p))[2])
  if(is.na(n.TISL.p[i])==TRUE){next}
  
  pgls_sl.glr.p<- gls(G_LER~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=test_g,method='ML')
  
  pgls_ti.glr.p<- gls(G_LER~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=test_g,method='ML')
  
  pgls_tisl.glr.p<- gls(G_LER~as.factor(Imputed_G_TISL),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=test_g,method='ML')
  
  skip_to_next <- FALSE
  tryCatch({pgls_tisl.gr.p<- gls(G_ER~as.factor(Imputed_G_TISL),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=test_g,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  
  tryCatch({pgls_ti.gr.p<- gls(G_ER~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=test_g,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  
  tryCatch({pgls_sl.gr.p<- gls(G_ER~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=test_g,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  
  b_pgls_sl.gr.p[i]<- pgls_sl.gr.p$coefficients[2]
  b_pgls_sl.glr.p[i]<- pgls_sl.glr.p$coefficients[2]
  b_pgls_ti.gr.p[i]<- pgls_ti.gr.p$coefficients[2]
  b_pgls_ti.glr.p[i]<- pgls_ti.glr.p$coefficients[2]
  b_pgls_tisl.gr.p[i]<- pgls_tisl.gr.p$coefficients[2]
  b_pgls_tisl.glr.p[i]<- pgls_tisl.glr.p$coefficients[2]
  
  n.SL.p[i]<- unname(summary(as.factor(test$SL.p))[2])
  n.TI.p[i]<-  unname(summary(as.factor(test$TI.p))[2])
  n.TISL.p[i]<-  unname(summary(as.factor(test$TISL.p))[2])
}

hist(n.SL.p)
abline(v=sum(LineageData$Imputed_Social_Learning_Binary))
hist(n.TI.p)
abline(v=sum(LineageData$Imputed_Technical_Innovation_Binary))
hist(n.TISL.p)
abline(v=sum(LineageData$Imputed_TISL))

#TI Taxa per Genus DR
hist(b_pgls_ti.gr.p,breaks=20,xlab = "Effect size")
abline(v=l2.gr.p$coefficients[2])
length(b_pgls_ti.gr.p[b_pgls_ti.gr.p>l2.gr.p$coefficients[2]])/1000 #PSEUDO P VALUE= 0.351
print(median(b_pgls_ti.gr.p)) #Median effect size= 0.009

#TI Lineage per Genus DR
hist(b_pgls_ti.glr.p,breaks=20,xlab = "Effect size")
abline(v=l2.glr.p$coefficients[2])
length(b_pgls_ti.glr.p[b_pgls_ti.glr.p>l2.glr.p$coefficients[2]])/1000 #PSEUDO P VALUE= 0.017
print(median(b_pgls_ti.glr.p)) #Median effect size= 0.020

#SL Taxa per Genus DR
hist(b_pgls_sl.gr.p,breaks=20,xlab = "Effect size")
abline(v=l1.gr.p$coefficients[2])
length(b_pgls_sl.gr.p[b_pgls_sl.gr.p>l1.gr.p$coefficients[2]])/1000 #PSEUDO P VALUE= 0.425
print(median(b_pgls_sl.gr.p)) #Median effect size= 0.025

#SL Lineage per Genus DR
hist(b_pgls_sl.glr.p,breaks=20,xlab = "Effect size")
abline(v=l1.glr.p$coefficients[2])
length(b_pgls_sl.glr.p[b_pgls_sl.glr.p>l1.glr.p$coefficients[2]])/1000 #PSEUDO P VALUE= 0.003
print(median(b_pgls_sl.glr.p)) #Median effect size= 0.031

#TISL Taxa per Genus DR
hist(b_pgls_tisl.gr.p,breaks=20,xlab = "Effect size")
abline(v=l3.gr.p$coefficients[2])
length(b_pgls_tisl.gr.p[b_pgls_tisl.gr.p>l3.gr.p$coefficients[2]])/1000 #PSEUDO P VALUE= 0.255
print(median(b_pgls_tisl.gr.p)) #Median effect size= 0.012

#TISL Lineage per Genus DR
hist(b_pgls_tisl.glr.p,breaks=20,xlab = "Effect size")
abline(v=l3.glr.p$coefficients[2])
length(b_pgls_tisl.glr.p[b_pgls_tisl.glr.p>l3.glr.p$coefficients[2]])/1000 #PSEUDO P VALUE= 0.003
print(median(b_pgls_tisl.glr.p)) #Median effect size= 0.023





####2.0 Research Effort Bias####
RE_dat<- subset(LineageData,is.na(Research_Effort)==F)

plot(Technical_Innovation_Binary~log(Research_Effort),data=RE_dat,bty='l')
abline(v=2.07,lty=5)

plot(Social_Learning_Binary~log(Research_Effort),data=RE_dat,bty='l')
abline(v=2.07,lty=5)

re_tree<- drop.tip(tree,setdiff(tree$tip.label,RE_dat$Trees_Name))
re_tree<- multi2di(re_tree, random = TRUE)
re<- log(RE_dat$Research_Effort)
names(re)<- RE_dat$Trees_Name

re.lambda<- geiger::fitContinuous(re_tree,re,model=c('lambda')) #Calculate lambda in log(Research Effort)


#Defining branch pruning and rescaling functions

pruningwise.branching.times <- function(phy) {   	
  ## calculates branching times = node ages, for ultrametric tree in pruningwise order
  ## !warning! no test that tree is in pruningwise order and ultrametric.
  ## branching time = age = time from the node to tips
  xx <- numeric(phy$Nnode) # times from root to nodes
  nt = length(phy$tip.label)
  interns <- which(phy$edge[, 2] > nt)
  for (i in rev(interns)) { # assumes preorder = 'intern'al nodes in reverse
    ## time of descendant = time of ancestor + branch length:
    xx[phy$edge[i,2] - nt] <- xx[phy$edge[i,1] - nt] + phy$edge.length[i]
  }
  depth <- xx[phy$edge[1, 1] - nt] + phy$edge.length[1] # total tree height
  xx <- depth - xx
  names(xx) <- if (is.null(phy$node.label)) (nt + 1):(nt + phy$Nnode) else phy$node.label
  return(xx)
}

pruningwise.distFromRoot <- function(phy) {
  ## distance from root to all nodes, for tree in pruningwise order
  ## !warning! no test that tree is in pruningwise order.
  nt = length(phy$tip.label)
  xx <- numeric(phy$Nnode + nt)
  for (i in length(phy$edge.length):1)
    xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
  names(xx) <- if (is.null(phy$node.label)) 1:(nt + phy$Nnode) else
    c(phy$tip.label, phy$node.label)
  return(xx)
}

transf.branch.lengths <-
  function(phy,
           model = c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"), 
           parameters = NULL, check.pruningwise = TRUE, check.ultrametric=TRUE, D=NULL, check.names = TRUE)	
  {
    if (!inherits(phy, "phylo")) stop("object \\"phy\\" is not of class \\"phylo\\".")
    model = match.arg(model)	
    if (model=="trend")
      if (is.ultrametric(phy))
        stop("the trend is unidentifiable for ultrametric trees.")	
    if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
    if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
    tol = 1e-10	
    if (check.pruningwise) phy = reorder(phy,"pruningwise")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    externalEdge = (des <= n)
    if (!is.null(phy$root.edge))
      if (phy$root.edge>0)
        stop("the tree is supposed to have no root edge (or of length 0).")
    
    ## Default parameters
    parameters.default = c(0,1,1,1,0,0)
    names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate", "sigma2_error")
    
    ## User defined parameters
    if (is.null(parameters)) {
      parameters = parameters.default
    } else { 
      if (class(parameters)!= "list") {
        stop("please specify parameters as a list().")
      } else {
        specified <- !c(is.null(parameters$alpha),
                        is.null(parameters$lambda),
                        is.null(parameters$kappa),
                        is.null(parameters$delta),
                        is.null(parameters$rate),
                        is.null(parameters$sigma2_error))
        parameters.user <- c(parameters$alpha,
                             parameters$lambda,
                             parameters$kappa,
                             parameters$delta,
                             parameters$rate,
                             parameters$sigma2_error)
        names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate", "sigma2_error")
        parameters <- parameters.default
        parameters[specified] <- parameters.user 
      }				
    }
    p = list(alpha = parameters[1],
             lambda = parameters[2],
             kappa = parameters[3],
             delta = parameters[4],
             rate = parameters[5],
             sigma2_error = parameters[6]) # note that sigma2_error = true_sigma2_error/sigma2
    
    root.edge = 0 # default, holds for most models. Assumes original tree has no root edge.
    diagWeight = rep(1,n)
    errEdge = rep(p$sigma2_error,n)
    
    ## BM model
    if (model %in% c("BM","trend")) {
      edge.length = phy$edge.length
    }	
    ## OU models
    OU = c("OUrandomRoot","OUfixedRoot")
    if (model %in% OU) {
      if (check.ultrametric){
        D = numeric(n) # adjustments to external branck lengths
        if (!is.ultrametric(phy)){
          flag = 1
          dis = pruningwise.distFromRoot(phy) # has all nodes
          D = max(dis[1:n]) - dis[1:n]
          D = D - mean(D)
          phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
        }
        ## phy is now ultrametric
      } else {
        if (is.null(D)) stop("Provide D if you choose check.ultrametric=F")
        if (length(D)!=n) stop("D should be a vector with one term for each tip in the tree")
        if (check.names) {
          if (is.null(names(D)))  stop("D is lacking names (tip labels)")
          ordr = match(phy$tip.label, names(D))
          if (sum(is.na(ordr))>0) stop("names of D do not match the tree tip labels.")
          D = D[ordr,drop=F]
        }
      }
      times <- pruningwise.branching.times(phy) # has internal nodes only
      Tmax <- max(times)
      alpha = p$alpha
      errEdge = errEdge*exp(-2*alpha*D[des[externalEdge]]) # adjust measurement errors for OU models
      ## OUrandomRoot model	
      if (model=="OUrandomRoot") {
        distFromRoot <-  exp(-2*alpha*times) # fixit: divide by 2 alpha??
        d1 = distFromRoot[anc-n] # distFromRoot has internal nodes only, not the n external nodes.
        d2 = numeric(N)
        d2[externalEdge]  = exp(-2*alpha*D[des[externalEdge]])
        d2[!externalEdge] = distFromRoot[des[!externalEdge]-n]
      }
      ## OUfixedRoot model
      if (model=="OUfixedRoot") {	
        distFromRoot <-  exp(-2*alpha*times)*(1 - exp(-2*alpha*(Tmax-times))) # fixit: divide by 2 alpha?
        d1 = distFromRoot[anc-n]
        d2 = numeric(N)
        d2[externalEdge] = exp(-2*alpha*D[des[externalEdge]]) * (1-exp(-2*alpha*(Tmax-D[des[externalEdge]])))
        d2[!externalEdge]= distFromRoot[des[!externalEdge]-n]
      }
      edge.length = d2 - d1
      root.edge = min(distFromRoot)
      diagWeight = exp(alpha*D)
    }
    ## lambda model
    if (model=="lambda") {
      lambda = p$lambda
      distFromRoot <- pruningwise.distFromRoot(phy)
      edge.length = phy$edge.length * lambda 
      edge.length[externalEdge] = edge.length[externalEdge] + (1-lambda)*distFromRoot[des[externalEdge]]
    }
    ## kappa model
    if (model=="kappa") {
      kappa = p$kappa
      edge.length = phy$edge.length^kappa
    }
    ## delta model
    if (model=="delta") {
      delta = p$delta
      distFromRoot <- pruningwise.distFromRoot(phy)
      depth = max(distFromRoot)
      edge.length = (distFromRoot[des]^delta - distFromRoot[anc]^delta)*depth^(1-delta)
    }
    ## early burst model
    if (model=="EB") {
      rate = p$rate
      if (rate==0) edge.length = phy$edge.length
      else {
        distFromRoot <- pruningwise.distFromRoot(phy)
        edge.length = (exp(rate*distFromRoot[des])-exp(rate*distFromRoot[anc]))/rate
      }			
    }
    
    edge.length[externalEdge] = edge.length[externalEdge] + errEdge # add measurement errors to the tree
    phy$edge.length = edge.length
    phy$root.edge = root.edge
    names(diagWeight) = phy$tip.label
    return(list(tree = phy, diagWeight = diagWeight))
  }

three.point.compute <-
  function(phy, P, Q = NULL, diagWeight = NULL, 
           check.pruningwise = TRUE, check.names = TRUE) 
  {
    ## For an extra tip Variance, like diagonal measurement error E:
    ## V = DVoD + E = D Vn D with Vn = Vo + D^{-1}ED^{-1} still from tree.
    ## D^{-1}ED^{-1} is diagonal --> add these terms to external branch length
    if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo".')
    if (check.pruningwise)	phy = reorder(phy,"pruningwise")
    if ((!check.names)&(check.pruningwise))
      stop("check.names has to be TRUE unless check.pruningwise=FALSE")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    root.edge = if (is.null(phy$root.edge)) 0 else phy$root.edge 
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    if (is.null(diagWeight)) {
      diagWeight = rep(1,n)
      names(diagWeight) = phy$tip.label
    } else {
      if (any(abs(diagWeight) <= .Machine$double.eps ^ 0.8))
        stop ("diagonal weights need to be non-zero.")
    }
    flag = 0
    if (is.null(Q)) {
      flag = 1
      Q = rep(1,n)
      names(Q) = phy$tip.label
    }
    P = as.matrix(P)
    Q = as.matrix(Q)
    if (check.names) {
      if (is.null(rownames(P))) stop("P needs to have row names.")
      ordr = match(phy$tip.label, rownames(P))
      if (sum(is.na(ordr))>0)
        stop("row names of P do not match the tree tip labels.")
      P = P[ordr,,drop=F]
      if (is.null(rownames(Q))) stop("Q needs to have row names.")
      ordr = match(phy$tip.label, rownames(Q))
      if (sum(is.na(ordr))>0)
        stop("row names of Q do not match the tree tip labels.")
      Q = Q[ordr,,drop=F]
      if (is.null(names(diagWeight))) stop("diagWeight needs to have names.")
      ordr = match(phy$tip.label, names(diagWeight))
      if (sum(is.na(ordr))>0)
        stop("names of diagWeight do not match the tree tip labels.")
      diagWeight = diagWeight[ordr,drop=F]
    }
    if (nrow(P)!=n)
      stop("the number of rows in P needs to be the same as the number of tips in the tree.")
    if (nrow(Q)!=n)
      stop("the number of rows in Q needs to be the same as the number of tips in the tree.")
    if (length(diagWeight)!=n)
      stop("the length of diagWeight needs to be the same as the number of tips in the tree.")
    ## now doing: Q'V^{-1}P = Q' (DVoD)^{-1} P = (D^{-1}Q)' Vo^{-1} (D^{-1}P)
    ##        and log|V|= log|Vo| + 2log|D|.
    P = cbind(rep(1,n),P)
    Q = cbind(rep(1,n),Q)
    P = P/diagWeight
    Q = Q/diagWeight
    colP = ncol(P)
    colQ = ncol(Q)
    nout = 2 + colP + colP^2 + colQ + colQ^2 + colP*colQ
    tmp=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
           as.integer(colP), as.integer(colQ), as.integer(ROOT), as.double(root.edge),
           as.double(phy$edge.length),as.integer(des), as.integer(anc),
           as.double(as.vector(P)), as.double(as.vector(Q)),
           result=double(nout))$result # P=y in threepoint.c, and Q=X
    logd=tmp[1] + 2*sum(log(diagWeight)) # vec11, P1 and Q1 not needed here
    PP=matrix(tmp[2+colP+ 1:(colP^2)], colP, colP)
    QQ=matrix(tmp[2+colP+colP^2+colQ+ 1:(colQ^2)], colQ, colQ)
    QP=matrix(tmp[2+colP+colP^2+colQ+colQ^2 + 1:(colP*colQ)], colQ, colP)
    
    if (flag==1) # Q absent
      return(list(vec11=PP[1,1], P1=PP[1,-1], PP=PP[-1,-1], logd=logd))
    return(list(vec11=PP[1,1], P1=PP[1,-1], PP=PP[-1,-1],
                Q1=QQ[1,-1], QQ=QQ[-1,-1], QP=QP[-1,-1], logd=logd)) 
    
  }


rescaled_tree<- transf.branch.lengths(re_tree,model=c('lambda'),parameters=list(lambda=re.lambda$opt$lambda))$tree #Rescale the phylogeny based on lambda for log(Research Effort)
bm_rescaled<- geiger::fitContinuous(rescaled_tree,re,model=c('BM')) #Calculate BM on rescaled phylogeny

tree2<-  multi2di(tree, random = TRUE)
rescaled_tree_full<- transf.branch.lengths(tree2,model=c('lambda'),parameters=list(lambda=re.lambda$opt$lambda))$tree #Rescale the phylogeny based on lambda for log(Research Effort)

##
b_sl.gr.p<- NA
b_sl.glr.p<- NA

pgl_sl.gr.p<- NA
pgl_sl.glr.p<- NA

b_ti.gr.p<- NA
b_ti.glr.p<- NA

pgl_ti.gr.p<- NA
pgl_ti.glr.p<- NA

b_tisl.gr.p<- NA
b_tisl.glr.p<- NA

b_pgls_sl.gr.p<- NA
b_pgls_sl.glr.p<- NA
b_pgls_ti.gr.p<- NA
b_pgls_ti.glr.p<- NA
b_pgls_tisl.gr.p<- NA
b_pgls_tisl.glr.p<- NA

n.SL.p<- NA
n.TI.p<- NA
n.TISL.p<- NA

sim.sl<- sim.char(tree.p,par=par_rates_sl,model='discrete',root=1,nsim=1000)
sim.ti<- sim.char(tree.p,par=par_rates_ti,model='discrete',root=1,nsim=1000)
sim.tisl<- sim.char(tree.p,par=par_rates_tisl,model='discrete',root=1,nsim=1000)
for(i in 1:1000){
  sim.Research_Effort<- diversitree::sim.character(rescaled_tree_full, par=c(bm_rescaled$opt$sigsq,0,6.7),x0=bm_rescaled$opt$z0,model='bbm') #Let research effort "evolve" across all lineages
  m<- match(LineageData$Trees_Name,rownames(sim.sl))
  m2<-  match(LineageData$Trees_Name,names(sim.Research_Effort))
  test<- data.frame(Genus=LineageData$Genus,Richness=LineageData$Richness,SL.p=sim.sl[m,,i],TI.p=sim.ti[m,,i],TISL.p=sim.tisl[m,,i],RE=sim.Research_Effort[m2])
  
  for(z in 1:nrow(test)){
    if(test$RE[z]<2.07){ #Cut-off for <8 studies (the minimum research effort where technical innovation has been observed) 
      test[z,3:5]<- 1 #If under-studied, turn into a zero
    }
  }
  
  test_g<- group_by(test,Genus)
  test_g<- summarize(test_g,Genus_Richness=sum(Richness),Genus_Lineage_Richness=n(), Genus_total_RE=sum(RE),Genus_max_RE=max(RE),Imputed_Genus_Technical_Innovation_Binary=max(TI.p), Imputed_Genus_Social_Learning_Binary=max(SL.p),Imputed_G_TISL=max(TISL.p))#Summarizing the data at Genus level
  test_g<- test_g %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), NA, x))
  test_g<- as.data.frame(test_g)
  rownames(test_g)<- test_g$Genus
  m<- match(test_g$Genus,genus_data_p$Genus)
  test_g[,9:10]<- genus_data_p[m,5:6]
  n.SL.p[i]<- unname(summary(as.factor(test$SL.p))[2])
  if(is.na(n.SL.p[i])==TRUE){next}
  n.TI.p[i]<-  unname(summary(as.factor(test$TI.p))[2])
  if(is.na(n.TI.p[i])==TRUE){next}
  n.TISL.p[i]<-  unname(summary(as.factor(test$TISL.p))[2])
  if(is.na(n.TISL.p[i])==TRUE){next}
  
  pgls_sl.glr.p<- gls(G_LER~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=test_g,method='ML')
  
  pgls_ti.glr.p<- gls(G_LER~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=test_g,method='ML')
  
  pgls_tisl.glr.p<- gls(G_LER~as.factor(Imputed_G_TISL),correlation=corPagel(0,form = ~Genus,phy=g_tree,fixed=T),data=test_g,method='ML')
  
  skip_to_next <- FALSE
  tryCatch({pgls_tisl.gr.p<- gls(G_ER~as.factor(Imputed_G_TISL),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=test_g,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  
  tryCatch({pgls_ti.gr.p<- gls(G_ER~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=test_g,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  
  tryCatch({pgls_sl.gr.p<- gls(G_ER~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(1,form = ~Genus,phy=g_tree,fixed=F),data=test_g,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  
  b_pgls_sl.gr.p[i]<- pgls_sl.gr.p$coefficients[2]
  b_pgls_sl.glr.p[i]<- pgls_sl.glr.p$coefficients[2]
  b_pgls_ti.gr.p[i]<- pgls_ti.gr.p$coefficients[2]
  b_pgls_ti.glr.p[i]<- pgls_ti.glr.p$coefficients[2]
  b_pgls_tisl.gr.p[i]<- pgls_tisl.gr.p$coefficients[2]
  b_pgls_tisl.glr.p[i]<- pgls_tisl.glr.p$coefficients[2]
  
}

hist(n.SL.p)
abline(v=sum(LineageData$Imputed_Social_Learning_Binary))
hist(n.TI.p)
abline(v=sum(LineageData$Imputed_Technical_Innovation_Binary))
hist(n.TISL.p)
abline(v=sum(LineageData$Imputed_TISL))

#TI Taxa per Genus DR
hist(b_pgls_ti.gr.p,breaks=20,xlab = "Effect size")
abline(v=l2.gr.p$coefficients[2])
length(b_pgls_ti.gr.p[b_pgls_ti.gr.p>l2.gr.p$coefficients[2]])/length(na.omit(b_pgls_ti.gr.p)) #PSEUDO P VALUE= 0.408
print(median(na.omit(b_pgls_ti.gr.p))) #Median effect size= 0.023

#TI Lineage per Genus DR
hist(b_pgls_ti.glr.p,breaks=20,xlab = "Effect size")
abline(v=l2.glr.p$coefficients[2])
length(b_pgls_ti.glr.p[b_pgls_ti.glr.p>l2.glr.p$coefficients[2]])/length(na.omit(b_pgls_ti.glr.p)) #PSEUDO P VALUE= 0.109
print(median(na.omit(b_pgls_ti.glr.p))) #Median effect size= 0.033

#SL Taxa per Genus DR
hist(b_pgls_sl.gr.p,breaks=20,xlab = "Effect size")
abline(v=l1.gr.p$coefficients[2])
length(b_pgls_sl.gr.p[b_pgls_sl.gr.p>l1.gr.p$coefficients[2]])/length(na.omit(b_pgls_sl.gr.p)) #PSEUDO P VALUE= 0.439
print(median(na.omit(b_pgls_sl.gr.p))) #Median effect size= 0.028

#SL Lineage per Genus DR
hist(b_pgls_sl.glr.p,breaks=20,xlab = "Effect size")
abline(v=l1.glr.p$coefficients[2])
length(b_pgls_sl.glr.p[b_pgls_sl.glr.p>l1.glr.p$coefficients[2]])/length(na.omit(b_pgls_sl.glr.p)) #PSEUDO P VALUE= 0.019
print(median(na.omit(b_pgls_sl.glr.p))) #Median effect size= 0.040

#TISL Taxa per Genus DR
hist(b_pgls_tisl.gr.p,breaks=20,xlab = "Effect size")
abline(v=l3.gr.p$coefficients[2])
length(b_pgls_tisl.gr.p[b_pgls_tisl.gr.p>l3.gr.p$coefficients[2]])/length(na.omit(b_pgls_tisl.gr.p)) #PSEUDO P VALUE= 0.281
print(median(na.omit(b_pgls_tisl.gr.p))) #Median effect size= 0.020

#TISL Lineage per Genus DR
hist(b_pgls_tisl.glr.p,breaks=20,xlab = "Effect size")
abline(v=l3.glr.p$coefficients[2])
length(b_pgls_tisl.glr.p[b_pgls_tisl.glr.p>l3.glr.p$coefficients[2]])/length(na.omit(b_pgls_tisl.glr.p)) #PSEUDO P VALUE= 0.027
print(median(na.omit(na.omit(b_pgls_tisl.glr.p)))) #Median effect size= 0.032
