## Termite body size evolution
## last update: 3/26/2021, N. Mizumoto
## format and comment for submission: 6/11/2021 N. Mizumoto
## script with R projects

# Environment -------------------------------------------------------------
# get path for r project
{
  RPROJ <- list(PROJHOME = normalizePath(getwd()))
  attach(RPROJ)
  rm(RPROJ)
}

# Setup + Read all data ---------------------------------------------------
{
  rm(list = ls())
  
  today <- Sys.Date()
  
  ## Packages
  {
    library(ape)
    library(geiger)
    library(ggplot2)
    library(mnormt)
    library(paleoTS)
    library(parallel)
    library(phangorn)
    library(phytools)
    library(qpcR)
    library(stringr)
    library(motmot)
    library(RRphylo)
    library(viridis)
    library( cowplot )
    library(qwraps2)
    library(sjPlot)
    library(tidyr)
    library(dplyr)
    library(nlme)
    library(car)
    library(RColorBrewer)
    library(ggExtra)
    library(caper)
    library(extrafont)
    font_import(pattern="PT")
    loadfonts()
  }
  
  ## functions
  source( file.path(PROJHOME, "scripts\\\\210222_Functions.R") )

  ## data
  {
    ## Head width
    {
      load(file.path(PROJHOME, "data\\\\HeadWidthData.rda"))
      
      ## Make species level dataset (as original datasets include several data points for each species. here calucurate mean value for each species)
      data <- d[,c("Family","Sub.Family", "Genus", "Species", "Subspecies", "caste",
                   "castenote", "Value", "WoodNester", "Forager", "AbeNesting", 
                   "TrueWorker", "Soldierless", "Fossil")]
      names <- c("Cryptocercidae", "Mastotermitidae", "GroupA", "Termopsidae",
                 "Archotermopsidae", "Hodotermitidae", "Stolotermitidae","GroupB",
                 "Kalotermitidae", "Archeorhinotermitidae", "Stylotermitidae",
                 "Serritermitidae", "Rhinotermitidae", "Termitidae","Sphaerotermitinae",
                 "Macrotermitinae", "Foraminitermitinae", "Apicotermitinae",
                 "Syntermitinae", "Nasutitermitinae", "Cubitermitinae", "Termitinae")
      data$Sub.Family <- factor(data$Sub.Family, levels=rev(names) )
      data$caste <- as.factor(data$caste)
      data$FullName <- paste(data$Genus, data$Species, data$Subspecies)
      Res <- NULL
      for(i in unique(data$FullName)){
        df <- data[data$FullName==i,]
        df <- df[!df$castenote=="small",]
        Size <- (tapply(df$Value, df$caste, mean))
        Size.Ratio <- c(Size[1]/Size[3], Size[2]/Size[3], Size[1]/Size[2])
        Res <- rbind(Res, data.frame(df[1,c("Family","Sub.Family", "Genus", "Species",
                                            "Subspecies", "WoodNester","Forager",
                                            "AbeNesting", "TrueWorker",
                                            "Soldierless", "Fossil")],
                                     Imago=Size[1], Soldier=Size[2], Worker=Size[3],
                                     Imago_Worker = Size.Ratio[1], 
                                     Soldier_Worker = Size.Ratio[2], 
                                     Imago_Soldier = Size.Ratio[3]))
      }
    }
    
    ## output data summary
    {
      Res.Living <- Res[Res$Fossil==0,]
      Res.Fossil <- Res[Res$Fossil==1,]
      
      ## output number of species
      ResPrint <- Res[Res$Genus!="Cryptocercus",]
      dim(ResPrint)
      print(paste("living species", dim(ResPrint[ResPrint$Fossil==0,])[1]))
      print(paste("fossil species", dim(ResPrint[ResPrint$Fossil==1,])[1]))
      print(paste("Imago", length(na.omit(ResPrint$Imago))))
      print(paste("Soldier", length(na.omit(ResPrint$Soldier))))
      print(paste("Worker", length(na.omit(ResPrint$Worker))))
      print(paste("Worker", length(na.omit(ResPrint$Worker))))
      
      print(paste("Imago-Living", length(na.omit(Res[Res$Fossil==0,]$Imago))))
      print(paste("Soldier-Living", length(na.omit(Res[Res$Fossil==0,]$Soldier))))
      print(paste("Worker-Living", length(na.omit(Res[Res$Fossil==0,]$Worker))))
      print(paste("TrueWorker-Living", length(na.omit(Res[Res$Fossil==0 & Res$TrueWorker==1,]$Worker))))
      print(paste("FalseWorker-Living", length(na.omit(Res[Res$Fossil==0 & Res$TrueWorker==0,]$Worker))))
      
      print(paste("Imago-Fossil", length(na.omit(Res[Res$Fossil==1,]$Imago))))
      print(paste("Soldier-Fossil", length(na.omit(Res[Res$Fossil==1,]$Soldier))))
      print(paste("Worker-Fossil", length(na.omit(Res[Res$Fossil==1,]$Worker))))
      
      print(paste("Number of Genus", length(unique(ResPrint$Genus))))
      
      ## Range of head width
      print(range(Res[Res$Family!="Cryptocercidae",]$Imago, na.rm=T))
      print(range(Res[Res$Family!="Cryptocercidae",]$Soldier, na.rm=T))
      print(range(Res[Res$Family!="Cryptocercidae",]$Worker, na.rm=T))
    }
  }
  
  ## Phylogeny
  {
    tree.file <- "path to DataS4-tree.tre"
    tree_no3rd <- read.nexus(file.path(PROJHOME, tree.file))
    
    # drop outgroup
    label <- tree_no3rd$tip.label
    label[label=="Coptotermes_priscus"] <-  "Coptotermes_priscus_fossil"
    tree_no3rd$tip.label <- label
    out.group <- is.na(str_locate(label, "Ter")[,1] | str_locate(label, "ter")[,1] |
                         str_locate(label, "cercus")[,1])
    termite.tree <- drop.tip(tree_no3rd, label[out.group])
    termite.tree <- ladderize(termite.tree)
    plot(termite.tree, cex=0.5)
    
    label<- termite.tree$tip.label
    fossil <- !is.na(str_locate(label, "fossil")[,1])
    
    ## label cleaning 
    ## to only include genus name for living species, and correct typo in tree file
    {
      n <- length(label)
      clean.label <- label
      for(i in 1:n){
        clean.label[i] <- str_replace(clean.label[i], "_sp", "")
        clean.label[i] <- str_replace(clean.label[i], "_A", "")
        clean.label[i] <- str_replace(clean.label[i], "_Cameroon", "")
        clean.label[i] <- str_replace(clean.label[i], "_Burundi", "")
        clean.label[i] <- str_replace(clean.label[i], "_C", "")
        clean.label[i] <- str_replace(clean.label[i], "_New_Guinea", "")
        clean.label[i] <- str_replace(clean.label[i], "_French_Guiana", "")
        clean.label[i] <- str_replace(clean.label[i], "_G", "")
        clean.label[i] <- str_replace(clean.label[i], "_EMPTY", "")
        clean.label[i] <- str_replace(clean.label[i], "_Japan", "")
        clean.label[i] <- str_replace(clean.label[i], "_Thailand", "")
        clean.label[i] <- str_replace(clean.label[i], "1", "")
        clean.label[i] <- str_replace(clean.label[i], "2", "")
        clean.label[i] <- str_replace(clean.label[i], "_n_nr", "")
        clean.label[i] <- str_replace(clean.label[i], "_nr_", "_")
        clean.label[i] <- str_replace(clean.label[i], "_cf", "")
        clean.label[i] <- str_replace(clean.label[i], "_living", "")
        if(!is.na(str_locate(clean.label[i], "_")[1])){
          if(str_locate(clean.label[i], "_")[1] == str_length(clean.label[i])){
            clean.label[i] <- str_replace(clean.label[i], "_", "")
          }
        }
        clean.label[i] <- str_replace(clean.label[i], "acanthotorax", "acanthothorax")
        clean.label[i] <- str_replace(clean.label[i], "singaporensis", "singaporiensis")
        clean.label[i] <- str_replace(clean.label[i], "wolfschwennigeri", 
                                      "wolfschwenningeri")
        clean.label[i] <- str_replace(clean.label[i], "electrodominicanus",
                                      "electrodominicus")
        clean.label[i] <- str_replace(clean.label[i], "Krishnatermes", 
                                      "Krishnatermes_yaddha")
        clean.label[i] <- str_replace(clean.label[i], "Coptotermes_priscus", 
                                      "Coptotermes_priscus_fossil")
      }
      clean.label[is.na(str_locate(clean.label, "_")[,1])] <- 
        paste0(clean.label[is.na(str_locate(clean.label, "_")[,1])], "_sp")
      termite.tree$tip.label <- clean.label
    }
  }
}

# Plot overall head width data --------------------------------------------
{
  # histogram for each caste
  {
    Res3 <- data.frame(Sub.Family = Res$Sub.Family, 
                       caste=rep(c("Imago","Soldier","Worker"), each= dim(Res)[1]), 
                       HW = c(Res$Imago, Res$Soldier, Res$Worker), Fossil=Res$Fossil,
                       Workertype = Res$TrueWorker)
    Res4 <- Res3[!is.na(Res3$HW),]
    Res4[Res4$Sub.Family=="Cryptocercidae",]$caste <- "xCockroaches"
    
    ggplot(Res4[Res4$caste=="Imago" | Res4$caste=="xCockroaches",], aes(x=HW, fill=interaction(caste,Fossil))) + 
      geom_histogram(alpha=0.6,  binwidth = 0.05, position = "stack")+
      scale_fill_viridis(discrete=TRUE) + theme_bw() + 
      theme(aspect.ratio = 3/10, legend.position = c(0.8,0.7)) +
      scale_x_continuous(limits = c(0,7))+
      scale_y_continuous(limits = c(0,80))+
      xlab("Head width (mm)") + ylab("Number of species")
    ggsave(paste0(today,"_ImagoDistribution.pdf"), width=4, height = 3)
    
    ggplot(Res4[Res4$caste=="Soldier",], aes(x=HW, fill=interaction(caste,Fossil))) + 
      geom_histogram(alpha=0.6,  binwidth = 0.05, position = "stack")+
      scale_fill_viridis(discrete=TRUE) + theme_bw() + 
      theme(aspect.ratio = 3/6, legend.position = c(0.8,0.7)) +
      scale_x_continuous(limits = c(0,7))+
      scale_y_continuous(limits = c(0,80))+
      xlab("Head width (mm)") + ylab("Number of species")
    ggsave(paste0(today,"_SoldierDistribution.pdf"), width=4, height = 3)
    
    Res4$Workertype[is.na(Res4$Workertype)] <- -1
    ggplot(Res4[Res4$caste=="Worker",], aes(x=HW, fill=interaction(-Workertype))) + 
      geom_histogram(alpha=0.6,  binwidth = 0.05, position = "stack")+
      scale_fill_viridis(discrete=TRUE) + theme_bw() + 
      theme(aspect.ratio = 3/6, legend.position = c(0.8,0.7)) +
      scale_x_continuous(limits = c(0,7))+
      scale_y_continuous(limits = c(0,80))+
      xlab("Head width (mm)") + ylab("Number of species")
    ggsave(paste0(today,"_WorkerDistribution.pdf"), width=4, height = 3)
  }
}

# Functions for Model fitting ---------------------------------------------
{
  ## Create data set for model-fitting
  # getmean = F (default): sample one species for living genus
  # getmean = T: get mean value for living genus
  data.create <- function(df=Res, AnalyzeCaste="Image", Logarithmï¼�T, tree = termite.tree, 
                          drop.Cryptocercus=T, drop.fossil.species=F, 
                          getmean = F){
    n = length(tree$tip.label)
    label = tree$tip.label
    fossil = !is.na(str_match(label, "fossil"))
    ## head width matching
    HW.value = Subfamily = Family = AbeNesting = TrueWorker = rep(0,n)
    if(getmean){
      HW.se = rep(0,n)
      HW.n = rep(1,n)  
    }
    new.label = label
    if(Logarithm){
      df[,c("Imago","Soldier","Worker")] <- log10(df[,c("Imago","Soldier","Worker")])
    }
    df.Living <- df[df$Fossil==0,]
    df.Fossil <- df[df$Fossil==1,]
    
    ## Fossil species: species level matching -> random sample for genus level id
    ## Living species: random sample for each genus
    for(i in 1:n){
      if(fossil[i]){
        new.label[i] <- str_replace(label[i], "_fossil", "")
        if(!is.na(str_locate(new.label[i], "_")[1])){ # have species name?
          Genus <- str_sub(new.label[i], 1, str_locate(new.label[i], "_")[1]-1)
          Species <- str_sub(new.label[i], 
                             str_locate(new.label[i], "_")[1]+1, str_length(new.label[i]))
          Family[i] = as.character(df.Fossil[df.Fossil$Genus==Genus,"Family"][1])
          Subfamily[i] = as.character(df.Fossil[df.Fossil$Genus==Genus,2][1])
          AbeNesting[i] = NA
          TrueWorker[i] = NA
          
          ## match HW data
          df.df <- df.Fossil[df.Fossil$Genus == Genus & df.Fossil$Species == Species,]
          if(dim(df.df)[1] > 0){
            ## with Imago data
            HW.value[i] <- df.df[,AnalyzeCaste]
          } else {
            HW.value[i] = NA
          }
        } else{ # without species name
          Genus <- new.label[i]
          Family[i] = as.character(df[df$Genus==Genus,"Family"][1])
          Subfamily[i] = as.character(df[df$Genus==Genus,2][1])
          AbeNesting[i] = NA
          TrueWorker[i] = NA
          
          HW.genus <- na.omit(df.Fossil[df.Fossil$Genus == Genus,AnalyzeCaste])
          if(length(HW.genus) > 0){
            if(getmean){
              HW.value[i] <- mean(HW.genus, na.rm=T)
              HW.se[i] <- se(HW.genus)
              HW.n[i] <- length(na.omit(HW.genus))
            } else {
              HW.value[i] <- sample(HW.genus, 1)
            }
          } else {
            HW.value[i] = NA
          }
          
        }
        new.label[i] <- paste0(new.label[i], "_fossil")
        
      } 
      else {
        Genus <- str_sub(new.label[i], 1, str_locate(new.label[i], "_")[1]-1)
        new.label[i] = Genus
        Family[i] = as.character(df.Living[df.Living$Genus==Genus,"Family"][1])
        Subfamily[i] = as.character(df.Living[df.Living$Genus==Genus,"Sub.Family"][1])
        AbeNesting[i] = as.character(df.Living[df.Living$Genus==Genus,"AbeNesting"][1])
        TrueWorker[i] = df.Living[df.Living$Genus==Genus,"TrueWorker"][1]
        
        HW.genus <- na.omit( df.Living[Res.Living$Genus == Genus,AnalyzeCaste] )
        if(length(HW.genus) > 0){
          if(getmean){
            HW.value[i] <- mean(HW.genus, na.rm=T)
            HW.se[i] <- se(HW.genus)
            HW.n[i] <- length(na.omit(HW.genus))
          } else {
            HW.value[i] <- sample(HW.genus, 1)
          }
        } else {
          HW.value[i] = NA
        }
      }
    }
    
    if(getmean){
      HWdata <- data.frame(Family, Subfamily, new.label, HW.value,
                           HW.se, HW.n, fossil,
                           AbeNesting, TrueWorker)
    } else {
      HWdata <- data.frame(Family, Subfamily, new.label, HW.value, fossil,
                           AbeNesting, TrueWorker)
    }
    termite.tree$tip.label <- new.label
    termite.tree2 <- drop.tip(termite.tree, new.label[is.na(HWdata$HW.value)])
    HWdata2 <- HWdata[!is.na(HWdata$HW.value),]
    row.names(HWdata2) <- 1:dim(HWdata2)[1]
    
    if(drop.Cryptocercus){
      termite.tree2 <- drop.tip(termite.tree2, "Cryptocercus")
      HWdata2 <- HWdata2[HWdata2$Family!="Cryptocercidae",]
    }
    
    if(drop.fossil.species){
      termite.tree2 <- drop.tip(termite.tree2, HWdata2[HWdata2$fossil,"new.label"])
      HWdata2 <- HWdata2[!HWdata2$fossil,]
    }
  
    return(list(termite.tree2, HWdata2))
  }

  ## Create data set for relationship (caste ratio)
  # getmean = F: sample one species for living genus
  # getmean = T (default): get mean value for living genus
  data.create2 <- function(AnalyzeCaste1, AnalyzeCaste2, tree = termite.tree, getmean=T){
    n = length(tree$tip.label)
    label = tree$tip.label
    fossil = !is.na(str_match(label, "fossil"))
    
    ## head width matching
    HW.ratio.sampled = HW.N = Subfamily = Family = AbeNesting = TrueWorker =
      HW.ratio.mean = HW.mean1 = HW.mean2 = rep(NA,n)
    new.label = clean.label
    Res.Living <- Res[Res$Fossil==0,]
    for(i in 1:n){
      if(!fossil[i]){
        # without species name
        Genus <- str_sub(clean.label[i], 1, str_locate(clean.label[i], "_")[1]-1)
        Family[i] = as.character(Res.Living[Res.Living$Genus==Genus,"Family"][1])
        Subfamily[i] = as.character(Res.Living[Res.Living$Genus==Genus,"Sub.Family"][1])
        AbeNesting[i] = as.character(Res.Living[Res.Living$Genus==Genus,"AbeNesting"][1])
        TrueWorker[i] = Res.Living[Res.Living$Genus==Genus,"TrueWorker"][1]
        new.label[i] = Genus
        
        if(length(na.omit(
          Res.Living[Res.Living$Genus == Genus, AnalyzeCaste1]/
          Res.Living[Res.Living$Genus == Genus,AnalyzeCaste2])) > 0){
          HW.value1 <- Res.Living[Res.Living$Genus == Genus,AnalyzeCaste1]
          HW.value2 <- Res.Living[Res.Living$Genus == Genus,AnalyzeCaste2]
          HW.ratio <- (HW.value1 / HW.value2)
          if(getmean){
            HW.ratio.mean[i] = mean(HW.ratio, na.rm=T)
            HW.mean1[i] = mean(HW.value1[!is.na(HW.value2)], na.rm=T)
            HW.mean2[i] = mean(HW.value2[!is.na(HW.value1)], na.rm=T)
          } else {
            if(length(na.omit(HW.ratio))==1){
              HW.ratio.sampled[i] = HW.ratio[!is.na(HW.ratio)]
              HW.mean1[i] = HW.value1[!is.na(HW.ratio)]
              HW.mean2[i] = HW.value2[!is.na(HW.ratio)]
            } else {
              x <- 1:length(HW.ratio[!is.na(HW.ratio)])
              x <- sample(x, 1)
              HW.ratio.sampled[i] <- HW.ratio[!is.na(HW.ratio)][x]
              HW.mean1[i] = HW.value1[!is.na(HW.ratio)][x]
              HW.mean2[i] = HW.value2[!is.na(HW.ratio)][x]
            }
          }
        } else {
          HW.ratio.sampled[i] = NA
        }
      }
    }
    
    if(getmean){
      HWdata <- data.frame(Family, Subfamily, new.label, 
                           AnalyzeCasteRatio = HW.ratio.mean,
                           HW.mean1 = HW.mean1, HW.mean2 = HW.mean2,
                           AbeNesting, TrueWorker)
    } else {
      HWdata <- data.frame(Family, Subfamily, new.label, 
                           AnalyzeCasteRatio = HW.ratio.sampled,
                           HW.mean1 = HW.mean1, HW.mean2 = HW.mean2,
                           AbeNesting, TrueWorker)
    }
    termite.tree$tip.label <- new.label
    termite.tree2 <- drop.tip(termite.tree, new.label[is.na(HWdata$AnalyzeCasteRatio)])
    HWdata2 <- HWdata[!is.na(HWdata$AnalyzeCasteRatio),]
    row.names(HWdata2) <- 1:dim(HWdata2)[1]
    
    return(list(termite.tree2, HWdata2))
  }
  
  ## Model fitting
  # bn = branch at which modes for evolution change (for trend-BM model)
  Model.fit <- function(tree, HWdata, HW.se, Fixed.ancester, x0=NULL,
                        EB.lower = -0.1, only.living=F, bn=1){
    
    Model = AncestralState = Variance1 = Variance2 = Mean1 = Mean2 = 
      Parameter = n.Param = AICc <- rep(NA, 8)
    
    # prep
    VV<- vcv(tree)  # the expected variances and covariances of a continuous trait
    nt<- Ntip(tree) # number of tips
    vs.init<- var(HWdata)/mean(diag(VV))  # very rough guess
    trm <- tree
    for(i in 1:length(bn)){
      trm <- paintSubTree(trm, node=bn, state=2, stem=1)
    }
    trm <- orderMappedEdge(trm, ordering="numerical")
    CC <- multiC(trm)
    mmi <- matrix(nrow=2, ncol=nt)
    mmi[1,]<- diag(CC[[1]])  # converts to a matrix
    mmi[2,]<- diag(CC[[2]])  # converts to a matrix
    n = length(HWdata)
    
    ### Unbiased Brownian Motion
    p0 <- c(log(vs.init), phyMean(HWdata, VV))
    names(p0)<- c("vs", "x0")
    w <- optim(p0, logL.uBM, method="L", control=cl<- list(fnscale= -1),
               hessian=F, x=HWdata, VV=VV, meserr=HW.se)
    k = 2
    Model[1] = "uBM"
    Variance1[1] = exp(w$par[1])
    AncestralState[1] = w$par[2]
    AICc[1] = -2*w$value+2*k*n/(n-k-1)
    n.Param[1] = k
    
    ### Ornstein-Uhlenbeck (Single Stationary Peak model)
    fit.ou <- fitContinuous(tree, HWdata, SE=HW.se, model="OU")
    Model[2] = "OU"
    Variance1[2] = fit.ou$opt$sigsq
    Parameter[2] = fit.ou$opt$alpha
    AncestralState[2] = fit.ou$opt$z0
    AICc[2] = fit.ou$opt$aicc
    n.Param[2] = 3
    
    ### Brownian Motion with Directional trend
    if(!only.living){
      p0 <- c(log(vs.init), lm(HWdata~diag(VV))$coefficients[2], phyMean(HWdata, VV))
      names(p0)<- c("vs", "ms", "x0")
      w <- optim(p0, logL.trend, method="L", control=cl<- list(fnscale= -1),
                 hessian=F, x=HWdata, VV=VV, meserr=HW.se)
      Model[3] = "Trend"
      Variance1[3] = exp(w$par[1])
      Mean1[3] = w$par[2]
      AncestralState[3] = w$par[3]
      k=3; n.Param[3] = k
      AICc[3] = -2*w$value+2*k*n/(n-k-1)
    }
    
    ### White noize
    p0 <- c(log(vs.init), phyMean(HWdata, VV))
    names(p0)<- c("vs", "x0")
    w <- optim(p0, lnl.noise, method="L", control=cl<- list(fnscale= -1),
               hessian=F, x=HWdata, VV=VV, meserr=HW.se)
    Model[4] = "White"
    Variance1[4] = exp(w$par[1])
    AncestralState[4] = w$par[2]
    k = 2; n.Param[4] = k
    AICc[4] =-2*w$value+2*k*n/(n-k-1)
    
    ### lambda
    p0 <- c(log(vs.init), 0, phyMean(HWdata, VV))
    names(p0)<- c("vs", "lambda", "x0")
    w <- optim(p0, lnl.lambda, method="L", control=cl<- list(fnscale= -1),
               hessian=F, x=HWdata, VV=VV, meserr=HW.se, upper = c(Inf,1,Inf),
               lower = c(-Inf,exp(-500),0))
    Model[5] = "lambda"
    Variance1[5] = exp(w$par[1])
    Parameter[5] = (w$par[2])
    AncestralState[5] = w$par[3]
    k = 3; n.Param[5] = k
    AICc[5] =-2*w$value+2*k*n/(n-k-1)
    
    ### kappa
    p0 <- c(log(vs.init), 0.01, phyMean(HWdata, VV))
    names(p0)<- c("vs", "kappa", "x0")
    w <- optim(p0, lnl.kappa, method="L", control=cl<- list(fnscale= -1),
               hessian=F, x=HWdata, VV=VV, meserr=HW.se, tree=tree)
    n = length(HWdata); k = 3
    Model[6] = "kappa"
    Variance1[6] = exp(w$par[1])
    AncestralState[6] = w$par[3]
    Parameter[6] = exp(w$par[2])
    AICc[6] = -2*w$value+2*k*n/(n-k-1)
    n.Param[6] = k
    
    ### delta
    p0 <- c(log(vs.init), 0.01, phyMean(HWdata, VV))
    names(p0)<- c("vs", "delta", "x0")
    w <- optim(p0, lnl.delta, method="L", control=cl<- list(fnscale= -1),
               hessian=F, x=HWdata, VV=VV, meserr=HW.se, tree=tree)
    n = length(HWdata); k = 3
    Model[7] = "delta"
    Variance1[7] = exp(w$par[1])
    AncestralState[7] = w$par[3]
    Parameter[7] = exp(w$par[2])
    AICc[7] = -2*w$value+2*k*n/(n-k-1)
    n.Param[7] = k
    
    if(length(bn)>0){
      
      ### Trend to unbiased BM (separate variance)
      p0<- c(0, rep(log(vs.init),2), phyMean(HWdata, VV))
      names(p0)<- c("ms1","vs1", "vs2", "x0")
      w <- optim(p0, logL.trendShift.vs, method="Nelder", control=list(fnscale= -1), 
                 hessian=F, bn=bn, mmi=mmi, x=HWdata, CC=CC, meserr=HW.se)	
      Model[8] = "Trend-uBM"
      Variance1[8] = exp(w$par[2])
      Variance2[8] = exp(w$par[3])
      Mean1[8] = w$par[1]
      AncestralState[8] = w$par[4]
      k = 4
      AICc[8] = -2*w$value+2*k*nt/(nt-k-1)
      
    }
    
    df <- data.frame(Model, AICc, AncestralState, Variance1, Variance2, Mean1, Mean2, Parameter)
    return(df)
  }
}

# Imago body size analysis ------------------------------------------------
{
  AnalyzeCaste <- "Imago"
  
  ## For plot the data
  {
    dl <- data.create(df=Res, AnalyzeCaste, tree=termite.tree, Logarithm=F, 
                      drop.Cryptocercus=T, drop.fossil.species=F,
                      getmean = T )
    HWdata <- dl[[2]]
    HW.tree <- dl[[1]]
    n <- length(HW.tree$tip.label)
    HW.tree <- multi2di(HW.tree, random=F)
    HWdata.mean <- HWdata$HW.value
    HWdata.se <- HWdata$HW.se
    HWdata.n <- HWdata$HW.n
    names(HWdata.mean) =names(HWdata.se) = names(HWdata.n) = HWdata$new.label
    
    ## phenogram (traitgram)
    {
      plot(HW.tree, cex=0.5)
      nodelabels(text=1:HW.tree$Nnode+Ntip(HW.tree),
                 node=1:HW.tree$Nnode+Ntip(HW.tree), cex=0.5)
      HW.tree<- paintSubTree(HW.tree, node=244, state=2, stem=1)
      HW.tree<- paintSubTree(HW.tree, node=151, state=4, stem=1)
      HW.tree<- paintSubTree(HW.tree, node=156, state=3, stem=1)
      plot(HW.tree, fsize=0.5)
      
      cols = c("#101010AA", "#1B9E77", "#D95F02", "#7570B3"); names(cols) <- 1:4
      cols = viridis(4); names(cols) <- 1:4
      
      #pdf(paste0(today,"Phenogram-Imago-All.pdf"))
      par(pin=c(4,3))
      phenogram(HW.tree, HWdata.mean, fsize=0.6, 
                ylab="Head width (mm)", xlab="Million years",
                log = "y", lty=1, lwd=2, 
                ylim=c(0.5,5), colors=cols, las=1,
                label=F, ftype="off", type="b", xlim=c(-9.4, 145))
      axis(1, at=c(-9.4, -9.4+50, -9.4+100, -9.4+150), labels=c(-150,-100,-50,0))
      #dev.off()  
    }
    
    ## Phylogeny
    ## create a split plot
    {
      HW.tree.plot <- HW.tree
      label.plot <- HW.tree.plot$tip.label
      label.plot <- str_replace(label.plot, "_", " ")
      label.plot <- str_replace(label.plot, "fossil", "")
      HW.tree.plot$tip.label <- label.plot
      HWdata.mean.plot <- HWdata.mean
      names(HWdata.mean.plot) <- label.plot
      
      is_tip <- HW.tree$edge[,2] <= length(HW.tree$tip.label)
      ordered_tips <- HW.tree$edge[is_tip, 2]
      
      #pdf(paste0(today,"-Phylogeny-Imago-All.pdf"), width=5, height=9, family = "PT Serif", paper = "a4")
      layout(matrix(c(1,2),1,2),c(0.2,0.2))
      plotTree(HW.tree.plot, mar=c(2, 2, 1, 0), fsize=0.4, col=cols, 
           xlim=c(-5,150), ftype="i")
      tiplabels(pch=20, cex=log((HWdata.n))/4, col="#00000099")
      axisPhylo(las=1)
      col2 <- viridis(2, option="inferno",alpha = 0.7)
      XonFig <- barplot(HWdata.mean.plot[ordered_tips],
              horiz=TRUE,width=0.5,space=1,
              ylim=c(1,length(HW.tree.plot$tip.label))-0.5, names="",
              xlim=c(0,4.5), col = col2[(HWdata$fossil*1+1)[ordered_tips]])
      arrows(HWdata.mean.plot[ordered_tips], XonFig,
             HWdata.mean.plot[ordered_tips]+ HWdata$HW.se[ordered_tips],XonFig, angle=90, length=0.025)
      
      abline(v=c(0,1,2,3,4))
      #dev.off()
    }

    ## Only living species
    {
    HW.tree <- dl[[1]]
    HW.Living <- HWdata[!HWdata$fossil,]
    HW.Living.Mean <- HW.Living$HW.value
    names(HW.Living.Mean) = HW.Living$new.label
    Living.tree <- drop.tip(HW.tree, HWdata$new.label[HWdata$fossil])
    
    ## phenogram
    dev.off()
    plot(Living.tree, cex=0.5)
    nodelabels(text=1:Living.tree$Nnode+Ntip(Living.tree),
               node=1:Living.tree$Nnode+Ntip(Living.tree), cex=0.5)
    Living.tree<- paintSubTree(Living.tree, node=203, state=2, stem=1)
    Living.tree<- paintSubTree(Living.tree, node=117, state=4, stem=1)
    Living.tree<- paintSubTree(Living.tree, node=122, state=3, stem=1)
    
    cols = c("#101010AA", "#1B9E77", "#D95F02", "#7570B3"); names(cols) <- 1:4
    cols = viridis(4); names(cols) <- 1:4    
    
    #pdf(paste0(today,"-Phenogram-Imago-Living.pdf"), family = "PT Sans", paper = "a4")
    par(pin=c(4,3))
    phenogram(Living.tree, HW.Living.Mean, fsize=0.6, 
              ylab="Head width (mm)", xlab="Million years",
              log = "y", lty=1, 
              ylim=c(0.5,5), colors=cols,
              label=T, ftype="off", type="b", xlim=c(-9.4, 145))
    axis(1, at=c(-9.4, -9.4+50, -9.4+100, -9.4+150), labels=c(-150,-100,-50,0))
    #dev.off()  
    }
  }
  
  ## For fitting
  {
    ## Analyze
    if(T){
      Res.Fitting.Sum <- NULL
      for(i in 1:100){
        dl <- data.create(df=Res, AnalyzeCaste, tree=termite.tree, Logarithm=T, drop.Cryptocercus=T)
        HWdata <- dl[[2]]
        HW.tree <- dl[[1]]
        HW.tree <- multi2di(HW.tree, random=F)
        n <- length(HW.tree$tip.label)
        HWdata.mean <- HWdata$HW.value
        HWdata.se <- rep(0, n)
        names(HWdata.mean) = names(HWdata.se) <- HWdata$new.label
        
        Res.all.Fitting <- Model.fit(tree = HW.tree, HWdata = HWdata.mean,
                                     HW.se = HWdata.se,
                                     Fixed.ancester = F, x0 = log10(4.8726667),
                                     EB.lower= 1e-10, bn=c(156))
        row.names(Res.all.Fitting) <- Res.all.Fitting$Model
        Res.all.Fitting <- Res.all.Fitting[c("uBM", "OU", "lambda", "delta", "kappa", "Trend", "Trend-uBM"),]
        plot(as.factor(Res.all.Fitting$Model), aicw(Res.all.Fitting$AICc)$w, main=paste("all", i), ylim=c(0,1))
      
        ## Only living species
        HW.Living <- HWdata[!HWdata$fossil,]
        HW.Living.Mean <- HW.Living$HW.value
        HW.Living.SE <- rep(0, length(HW.Living.Mean))
        names(HW.Living.Mean) = names(HW.Living.SE) = HW.Living$new.label
        Living.tree <- drop.tip(HW.tree, HWdata$new.label[HWdata$fossil])
        
        names(HWdata.mean) = names(HWdata.se) <- HW.Living$new.label
        
        Res.living.fitting <- Model.fit(tree = Living.tree, 
                                        HWdata = HW.Living.Mean, 
                         HW.se = HW.Living.SE,
                         Fixed.ancester = F, x0 = log10(4.8726667),
                         EB.lower= 1e-2, only.living=T, bn=c(122))
        
        Res.living.fitting <- Res.living.fitting[!is.na(Res.living.fitting$Model),]
        row.names(Res.living.fitting) <- Res.living.fitting$Model
        Res.living.fitting <- Res.living.fitting[c("uBM", "OU", "lambda", "delta", "kappa", "Trend-uBM"),]
        plot(as.factor(Res.living.fitting$Model), aicw(Res.living.fitting$AICc)$w, main=paste("living", i), ylim=c(0,1))
        
        df <- data.frame(
          Data = c(rep("Living+Fossil", 
                       dim(Res.all.Fitting)[1]), rep("Living", dim(Res.living.fitting)[1])),
          Model = c(Res.all.Fitting$Model, Res.living.fitting$Model),
          AICc = c(Res.all.Fitting$AICc, Res.living.fitting$AICc),
          wAIC = c(aicw(Res.all.Fitting$AICc)$w, aicw(Res.living.fitting$AICc)$w),
          AncestralState = c(Res.all.Fitting$AncestralState, Res.living.fitting$AncestralState),
          Variance1 = c(Res.all.Fitting$Variance1, Res.living.fitting$Variance1),
          Variance2 = c(Res.all.Fitting$Variance2, Res.living.fitting$Variance2),
          Mean1 = c(Res.all.Fitting$Mean1, Res.living.fitting$Mean1),
          Parameter = c(Res.all.Fitting$Parameter, Res.living.fitting$Parameter)
        )
        df$Model <- factor(df$Model, levels=(Res.all.Fitting$Model) )
        Res.Fitting.Sum <- rbind(Res.Fitting.Sum, df)
        print(i)
      }
      #save(Res.Fitting.Sum, file=paste0("data/", today, "-Res.Fitting.Sum.rda"))
    }
    
    ## plot
    {
      load("data/2021-03-24-Res.Fitting.Sum.rda")
      Res.Fitting.Sum.Sum <- data_summrize(Res.Fitting.Sum, "wAIC", c("Model", "Data"))
      Res.Fitting.Sum.Sum$Model <- rep(levels(Res.Fitting.Sum$Model), 2)
      ggplot(Res.Fitting.Sum, aes(x=Model, y=wAIC, col=Data)) + 
        geom_point(alpha = 0.2, position = position_jitter(.05), size = 1) +
        scale_colour_brewer(palette = "Dark2")+
        guides(fill=FALSE) +
        theme_bw() + 
        scale_y_continuous(expand = c(0, 0.05), limits = c(0,1)) +
        theme(legend.position = c(0.1, 0.9), aspect.ratio = 1,
              axis.text.x = element_text(size=6, angle = 30, hjust = 1)) +
        geom_errorbar( data=Res.Fitting.Sum.Sum, 
                       aes(x=Model, y=wAIC.mean, ymin=wAIC.mean,
                           ymax=wAIC.mean, col=Data),
                       width = 0.5)
      #
      ggsave(paste0(today, "-Fitting-Imago.pdf"), width=5, height = 3, family="PT Sans", paper="a4")
      
     
      ## table
      data <- Res.Fitting.Sum
      data$AncestralStateValue <- 10^data$AncestralState
      data$P <- data$Parameter
      data$P[is.na(data$P)] <- data$Mean1[is.na(data$P)]
      data$P[is.na(data$P)] <- data$Variance1[is.na(data$P)]
      
      sdtable <- data %>%
        dplyr::group_by(Data, Model) %>% 
        dplyr::summarise(
          "AICc" = mean_sd(AICc, denote_sd = "paren", 
                          markup = getOption("qwraps2_markup", "markdown")),
          "wAIC" = paste0(formatC(mean(wAIC), format = "g", digits = 2),
                          " (", formatC(sd(wAIC), format = "g", digits = 2),
                          ")"),
          "Root" = mean_sd(AncestralStateValue, denote_sd = "paren",
                          markup = getOption("qwraps2_markup", "markdown")),
          "Parameter" = paste0(formatC(mean(P), format = "g", digits = 2),
                          " (", formatC(sd(P), format = "g", digits = 2),
                          ")")
          )
      
      tab_df(sdtable)
      #tab_df(sdtable,  file = paste0(today, "_fitting-table.html"))
      #tab_df(sdtable,  file = paste0(today, "_fitting-table.doc"))
      #save(Res.Fitting.Sum, file="Res.Fitting.Sum.rda")
    }
  }
  
  ## Comparison of variance between nesting type or true worker
  {
    ## use all data
    {
      Res.Imago <- Res[!is.na(Res$Imago) & !Res$Fossil & Res$Family!="Cryptocercidae",]
      ggplot(Res.Imago, aes(x=as.factor(TrueWorker), y=Imago, fill=as.factor(TrueWorker))) + 
        geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth=0.04) +
        scale_fill_viridis(discrete = T, end=0.5, alpha=0.5)+
        theme_bw() + ylab("Imago head width (mm)") +
        theme(legend.position = "none", aspect.ratio = 1) +
        stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                     geom="pointrange", color="red") 
      ggsave(paste0(today, "_TrueWorkerImago.pdf"), height=3, width=3, family="PT Sans", paper="a4")
      
      bartlett.test(Imago~as.factor(TrueWorker), data=Res.Imago)
      tapply(Res.Imago$Imago, Res.Imago$TrueWorker, var)
      
      Res.Imago$AbeNesting<- factor(Res.Imago$AbeNesting, levels=c("OP","MP","CP"))
      ggplot(Res.Imago, aes(x=(AbeNesting), y=Imago, fill=as.factor(AbeNesting))) + 
        geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth=0.04) +
        scale_fill_viridis(discrete = T, end=1, alpha=0.5, direction=-1)+
        scale_x_discrete(limits = c("OP", "MP", "CP")) +
        theme_bw() + ylab("Imago head width (mm)") +
        theme(legend.position = "none", aspect.ratio = 1) +
        stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                     geom="pointrange", color="red") 
      ggsave(paste0(today, "_NestingImago.pdf"), height=3, width=3, family="PT Sans", paper="a4")
      
      bartlett.test(Imago~as.factor(AbeNesting), data=Res.Imago)
      tapply(Res.Imago$Imago, Res.Imago$AbeNesting, var)
      bartlett.test(Imago~as.factor(AbeNesting), data=Res.Imago[Res.Imago$AbeNesting!="OP",])
      bartlett.test(Imago~as.factor(AbeNesting), data=Res.Imago[Res.Imago$AbeNesting!="MP",])
      bartlett.test(Imago~as.factor(AbeNesting), data=Res.Imago[Res.Imago$AbeNesting!="CP",])
    }
  }
  
  ## Comparison of mean between nesting type or true worker
  {
    ## use mean
    dl <- data.create(df=Res, AnalyzeCaste, tree=termite.tree, Logarithm=F,
                      drop.Cryptocercus=T, getmean = T)
    HWdata <- dl[[2]]
    HW.Living <- HWdata[!HWdata$fossil,]
    HW.Living <- HW.Living[,c(1:4,8:9)]
    HW.tree <- dl[[1]]
    Living.tree <- drop.tip(HW.tree, HWdata$new.label[HWdata$fossil])
    
    res1=res2 =NULL
    
    HW.Living$TrueWorker <- as.factor(HW.Living$TrueWorker)
    comp.data<-comparative.data(Living.tree, HW.Living, names.col="new.label")
    
    # True worker
    model1<-pgls(HW.value~TrueWorker, data=comp.data)
    model2<-pgls(HW.value~TrueWorker, data=comp.data, lambda="ML")
    model3<-pgls(HW.value~TrueWorker, data=comp.data, kappa="ML")
    model4<-pgls(HW.value~TrueWorker, data=comp.data, delta="ML")
    model5<-pgls(HW.value~TrueWorker, data=comp.data, lambda="ML", kappa="ML", delta="ML")
    AIC(model1, model2, model3, model4, model5)
    summary(model5)
    anova(model5) 
    
    # AbeNesting
    model1<-pgls(HW.value~AbeNesting, data=comp.data)
    model2<-pgls(HW.value~AbeNesting, data=comp.data, lambda="ML")
    model3<-pgls(HW.value~AbeNesting, data=comp.data, kappa="ML")
    model4<-pgls(HW.value~AbeNesting, data=comp.data, delta="ML")
    model5<-pgls(HW.value~AbeNesting, data=comp.data, lambda="ML", kappa="ML", delta="ML")
    AIC(model1, model2, model3, model4, model5)
    summary(model5)
    anova(model5)
    
  }

  ## Correlation with colony size
  {
    load(file.path(PROJHOME, "data\\\\ColonySize.rda"))
    Res.Imago <- Res[!is.na(Res$Imago) & !Res$Fossil & Res$Family!="Cryptocercidae",]
    
    Res.colony.size <- NULL
    for(i in 1:dim(d.cs)[1]){
      s.match <- Res.Imago$Genus == d.cs$Genus[i] & Res.Imago$Species == d.cs$Species[i]
      if(sum(s.match)==1){
        df <- cbind(Res.Imago[s.match,], Colony.Size = d.cs[i,3])
        Res.colony.size <- rbind(Res.colony.size, df)
      }
    }
    Res.colony.size$AbeNesting = factor(Res.colony.size$AbeNesting, levels=c("OP","MP","CP"))
    
    ggplot(Res.colony.size, aes(x=Colony.Size, y=Imago, col=AbeNesting))+
      geom_point(size=2)+
      scale_color_manual(values= c("#78bbbad2", "#43005399", "#d95f02ff"))+
      scale_x_log10(limits=c(100,20000000)) + scale_y_log10(limits=c(0.7,4))+
      theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
    #ggsave(paste0(today, "ColonySizeImago.pdf"), width=3, height=3)
    
    Imago.mean.cs <- tapply(Res.colony.size$Imago, Res.colony.size$Genus, mean)
    CS.mean.cs <- tapply(Res.colony.size$Colony.Size, Res.colony.size$Genus, mean)
    d.label <- names(CS.mean.cs)
  
    dl <- data.create(df=Res, AnalyzeCaste, tree=termite.tree, Logarithm=T, drop.Cryptocercus=T)
    HW.tree <- dl[[1]]
    HW.tree <- multi2di(HW.tree, random=F)
    HWdata <-  dl[[2]]
    HW.Living <- HWdata[HWdata$fossil==0,]
    Living.tree <- drop.tip(HW.tree, HWdata$new.label[HWdata$fossil])
    
    t.label <- Living.tree$tip.label
    
    include1 <- rep(F, length(t.label))
    include2 <- rep(F, length(d.label))
    for(i in 1:length(d.label)){
      a <- !is.na(str_locate(t.label, d.label[i])[,1])
      if(sum(a)>0){
        include1[a] <- T
        include2[i] <- T
      } else {
        #print(i)
      }
    }
  
  
    df <- data.frame(label = t.label[include1], Imago  = Imago.mean.cs[include2],
                     Colony.Size = CS.mean.cs[include2])
    TrueWorker = HW.Living$TrueWorker
    AbeNesting = HW.Living$AbeNesting
    names(AbeNesting) = names(TrueWorker) = HW.Living$new.label
    df$TrueWorker <- TrueWorker[df$label]
    df$AbeNesting <- AbeNesting[df$label]
    df$Colony.Size <- log10(df$Colony.Size)
    cs.tree <- drop.tip(Living.tree, t.label[!include1])
    
    comp.data<-comparative.data(cs.tree, df, names.col="label", vcv=T, 
                                vcv.dim=3, warn.dropped=TRUE)
    model1<-pgls(Imago~(Colony.Size), data=comp.data)
    model2<-pgls(Imago~(Colony.Size), data=comp.data, lambda="ML", bounds=list(lambda=c(0.1,1)))
    model3<-pgls(Imago~(Colony.Size), data=comp.data, kappa="ML", bounds=list(kappa=c(0.1,1)))
    model4<-pgls(Imago~(Colony.Size), data=comp.data, delta="ML")
    model5<-pgls(Imago~(Colony.Size), data=comp.data, lambda="ML", kappa="ML", delta="ML",
                 bounds=list(kappa=c(0.1,1), lambda=c(0.1,1)))
    AIC(model1, model2, model3, model4, model5)
    summary(model5)
    anova(model5) 
  }
}


# Caste comparison --------------------------------------------------------
# data prep
{
  Res.Living.plot <- Res.Living
  Res.Living.plot[Res.Living.plot$Soldierless==1,"Soldier_Worker"] = 
    rnorm(length(Res.Living.plot[Res.Living.plot$Soldierless==1,"Soldier_Worker"]), 0, 0.05)
  Res.Living.plot <- Res.Living.plot[!is.na(Res.Living.plot$Imago_Worker)&!
                                       is.na(Res.Living.plot$Soldier_Worker),]
  Res.Living.plot$Sub.Family <- as.character(Res.Living.plot$Sub.Family)
  
  Res.Living.plot[Res.Living.plot$Sub.Family=="Stylotermitidae","Sub.Family"] <- "Rhinotermitidae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Stolotermitidae","Sub.Family"] <- "Archotermopsidae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Mastotermitidae","Sub.Family"] <- "Archotermopsidae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Hodotermitidae","Sub.Family"] <- "Archotermopsidae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Foraminitermitinae","Sub.Family"] <- "Macrotermitinae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Sphaerotermitinae","Sub.Family"] <- "Macrotermitinae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Sphaerotermitinae","Sub.Family"] <- "Macrotermitinae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Archotermopsidae","Sub.Family"] <- "Kalotermitidae"
  Res.Living.plot[Res.Living.plot$Sub.Family=="Kalotermitidae","Sub.Family"] <- "non.Neoisoptera"
  
  names <- c("non.Neoisoptera", "Rhinotermitidae", "Macrotermitinae",
             "Apicotermitinae", "Syntermitinae", "Nasutitermitinae",
             "Cubitermitinae", "Termitinae")
  Res.Living.plot$Sub.Family <- factor(Res.Living.plot$Sub.Family, levels=rev(names) )
  cols = c(rev(brewer.pal(n=8, name="Dark2")))
  cols <- adjustcolor(cols,0.7)
}

## overall plot
{
  #pdf(paste0(today, "-CasteProportion.pdf"),width=5,height = 5)
  g <- ggplot(Res.Living.plot, aes(x=Imago_Worker, y=Soldier_Worker, col=Sub.Family)) + 
    geom_point(size=2) + 
    scale_color_manual(values = rev(cols)) + 
    #scale_x_continuous(limits = c(0.45,2.5)) + scale_y_continuous(limits = c(0.45,2.5)) + 
    theme_classic() +
    theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank()) +
    geom_hline(yintercept = seq(0.5,2.5,0.5), lty=2) + 
    geom_vline(xintercept = seq(0.5,2.5,0.5), lty=2) +
    ggtitle("Rhinotermitidae+Stylotermitidae,
            Macrotermitinae+Sphaerotermitinae+Foraminitermitinae")
  g
  ggMarginal(g, Res.Living.plot[Res.Living.plot$Soldierless==0,], 
             x=Imago/Worker, y=Soldier/Worker,
             type = "density", margins = "both", size = 7,
             groupColour = F, groupFill = F)
  #dev.off()
  #ggMarginal(g, Res.Living.plot, x=Imago/Worker, y=Soldier/Worker,
  #           type = "boxplot", margins = "both", size = 4,
  #           groupColour = F, groupFill = T)
}
  
## Imago-Worker
{
  ## count number of species
  sum(Res.Living[Res.Living$TrueWorker==0,]$Imago <= Res.Living[Res.Living$TrueWorker==0,]$Worker, na.rm=T)
  sum(Res.Living[Res.Living$TrueWorker==0,]$Imago > Res.Living[Res.Living$TrueWorker==0,]$Worker, na.rm=T)
  11/21
  
  sum(Res.Living[Res.Living$TrueWorker==1,]$Imago <= Res.Living[Res.Living$TrueWorker==1,]$Worker, na.rm=T)
  sum(Res.Living[Res.Living$TrueWorker==1,]$Imago > Res.Living[Res.Living$TrueWorker==1,]$Worker, na.rm=T)
  400/467
  
  ## analysis
  dl <- data.create2(AnalyzeCaste1 = "Imago",
                     AnalyzeCaste2 = "Worker")
  HWdata <- dl[[2]]
  HW.tree <- dl[[1]]
  HW.tree <- multi2di(HW.tree, random=F)
  {
    HWdata$Subfamily <- as.character(HWdata$Subfamily)
    HWdata[HWdata$Subfamily=="Stylotermitidae","Subfamily"] <- "Rhinotermitidae"
    HWdata[HWdata$Subfamily=="Stolotermitidae","Subfamily"] <- "Archotermopsidae"
    HWdata[HWdata$Subfamily=="Mastotermitidae","Subfamily"] <- "Archotermopsidae"
    HWdata[HWdata$Subfamily=="Hodotermitidae","Subfamily"] <- "Archotermopsidae"
    HWdata[HWdata$Subfamily=="Foraminitermitinae","Subfamily"] <- "Macrotermitinae"
    HWdata[HWdata$Subfamily=="Sphaerotermitinae","Subfamily"] <- "Macrotermitinae"
    HWdata[HWdata$Subfamily=="Sphaerotermitinae","Subfamily"] <- "Macrotermitinae"
    HWdata[HWdata$Subfamily=="Archotermopsidae","Subfamily"] <- "Kalotermitidae"
    HWdata[HWdata$Subfamily=="Kalotermitidae","Subfamily"] <- "non.Neoisoptera"
    names <- c("non.Neoisoptera", "Rhinotermitidae", "Macrotermitinae",
               "Apicotermitinae", "Syntermitinae", "Nasutitermitinae",
               "Cubitermitinae", "Termitinae")
    HWdata$Subfamily <- factor(HWdata$Subfamily, levels=rev(names) )
  }
  
  ## correlation (with true-worker)
  {
    HWdata2 <- HWdata[HWdata$TrueWorker==1,]
    HW.tree2 <- drop.tip(HW.tree, HWdata[HWdata$TrueWorker==0,]$new.label)
    
    comp.data<-comparative.data(HW.tree2, HWdata2, names.col="new.label", vcv=T, 
                                vcv.dim=3, warn.dropped=TRUE)
    model1<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data)
    model2<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML")
    model3<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, kappa="ML")
    model4<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, delta="ML")
    model5<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML", kappa="ML", delta="ML")
    
    AIC(model1, model2, model3, model4, model5)
    
    summary(model5)
    anova(model5)
    
    ggplot(Res.Living.plot[Res.Living.plot$TrueWorker==1,], aes(x=Worker, y=Imago, col=Sub.Family)) +
      geom_point(size=0.5) + scale_color_manual(values = rev(cols)) +
      coord_fixed() +
      geom_point(data=HWdata2, aes(x=HW.mean2, y=HW.mean1, col=Subfamily), pch=3)+
      geom_abline(slope=1, intercept = 0) +
      geom_abline(slope=model4$model$coef[2], intercept = model4$model$coef[1], col=2)+
      scale_x_log10(limit=c(0.4,5)) + scale_y_log10(limit=c(0.4,5)) + 
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank())+
      annotation_logticks()
    ggsave(paste0(today, "-Imago-TrueWorker.pdf"), width = 3, height=4, family="PT Sans", paper="a4")
  }
  
  ## correlation (with pseudergate)
  {
    HWdata2 <- HWdata[HWdata$TrueWorker==0,]
    HW.tree2 <- drop.tip(HW.tree, HWdata[HWdata$TrueWorker==1,]$new.label)
    
    comp.data<-comparative.data(HW.tree2, HWdata2, names.col="new.label", vcv=T, 
                                vcv.dim=3, warn.dropped=TRUE)
    model1<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data)
    model2<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML")
    model3<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, kappa="ML")
    model4<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, delta="ML")
    model5<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML", kappa="ML", delta="ML")
    
    AIC(model1, model2, model3, model4, model5)
    
    summary(model5)
    anova(model5)
    
    ggplot(Res.Living.plot[Res.Living.plot$TrueWorker==0,], aes(x=Worker, y=Imago, col=Sub.Family)) +
      geom_point(size=0.5) + scale_color_manual(values = rev(cols)) +
      coord_fixed() +
      geom_point(data=HWdata2, aes(x=HW.mean2, y=HW.mean1, col=Subfamily), pch=3)+
      geom_abline(slope=1, intercept = 0) +
      geom_abline(slope=model4$model$coef[2], intercept = model4$model$coef[1], col=2)+
      scale_x_log10(limit=c(0.4,5)) + scale_y_log10(limit=c(0.4,5)) + 
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank())+
      annotation_logticks()
    ggsave(paste0(today, "-Imago-Pseudogate.pdf"), width = 3, height=4, family="PT Sans", paper="a4")
  }
  
  ## Comparison between worker types
  {
    HWdata$pdisp <- (HWdata$HW.mean1-HWdata$HW.mean2)/HWdata$HW.mean2
    HWdata$TrueWorker <- as.factor(HWdata$TrueWorker)
    comp.data<-comparative.data(HW.tree, HWdata, names.col="new.label")
    model1<-pgls(pdisp~TrueWorker, data=comp.data)
    model2<-pgls(pdisp~TrueWorker, data=comp.data, lambda="ML")
    model3<-pgls(pdisp~TrueWorker, data=comp.data, kappa="ML")
    model4<-pgls(pdisp~TrueWorker, data=comp.data, delta="ML")
    model5<-pgls(pdisp~TrueWorker, data=comp.data, lambda="ML", kappa="ML", delta="ML")
    AIC(model1, model2, model3, model4, model5)
    summary(model5)
    anova(model5)
    
    HWdata$TrueWorker <- factor(HWdata$TrueWorker, levels = c(1,0))
    comp.data<-comparative.data(HW.tree, HWdata, names.col="new.label", vcv=T, 
                                vcv.dim=3, warn.dropped=TRUE)
    model5<-pgls(pdisp~TrueWorker, data=comp.data, lambda="ML", kappa="ML", delta="ML")
    summary(model5)
    
    bartlett.test(pdisp~TrueWorker, data=HWdata)
    
    Res.Living.plot$pdisp <- 
      (Res.Living.plot$Imago-Res.Living.plot$Worker)/Res.Living.plot$Worker
    bartlett.test(pdisp~TrueWorker, data=Res.Living.plot)
    ggplot(Res.Living.plot, aes(x=TrueWorker, y=pdisp, fill=as.factor(TrueWorker))) +
      geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.02, dotsize = 1, alpha = .5)+
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                   geom="pointrange", color="red") +
      scale_fill_manual(values = cols[c(1,2)]) +
      theme_bw() + ylab("(Imago - Worker) / Worker") +
      theme(legend.position = "none", aspect.ratio = 1)
    ggsave(paste0(today, "_TrueWorker-Imago-Worker.pdf"),
           height=3, width=3, family="PT Sans", paper="a4")
    
    bartlett.test(pdisp~AbeNesting, data=Res.Living.plot)
    
    ## ancestral state reconstruction
    dl <- data.create2(AnalyzeCaste1 = "Imago",
                       AnalyzeCaste2 = "Worker")
    HWdata <- dl[[2]]
    HW.tree <- dl[[1]]
    HW.tree <- multi2di(HW.tree, random=F)
    
    TrueWorker <- HWdata$TrueWorker
    names(TrueWorker) <- HW.tree$tip.label
    TrueWorker <- factor(TrueWorker, levels = c("0", "1"))
    HWdata$pdisp <- (HWdata$HW.mean1 - HWdata$HW.mean2 ) / HWdata$HW.mean1
    
    mtrees <- make.simmap(HW.tree, TrueWorker, model="ER", nsim=100)
    pd<-summary(mtrees,plot=FALSE)
    
    
    is_tip <- HW.tree$edge[,2] <= length(HW.tree$tip.label)
    ordered_tips <- HW.tree$edge[is_tip, 2]
    
    cols = viridis(2, end=0.5)
    names(cols) <- c("0", "1")
    {
      pdf(paste0(today,"-Phylogeny-Trueworker.pdf"), width=5, height=7,
          family = "PT Serif", paper = "a4")
      layout(matrix(c(1,2),1,2),c(0.2,0.2))
      plot(mtrees[2], mar=c(5, 2, 3, 0), fsize=0.4, colors=cols, 
           xlim=c(-5,150), ftype="i")
      nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
      tiplabels(pie=to.matrix(TrueWorker,sort(unique(TrueWorker))),piecol=cols,cex=0.4)
      axisPhylo(las=1)
      col2 <- viridis(2, option="inferno",alpha = 0.7)
      XonFig <- barplot(HWdata$AnalyzeCasteRatio[ordered_tips],
                        horiz=TRUE,width=0.5,space=1,
                        ylim=c(1,length(HW.tree$tip.label))-0.5, names="",
                        xlim=c(0,3), col = cols[HWdata$TrueWorker[ordered_tips]+1])
      abline(v=c(0,1,2))
      dev.off()
    }
  }
}

## Soldier-Worker
{
  cols = c(rev(brewer.pal(n=8, name="Dark2")))
  cols <- adjustcolor(cols,0.7)
  
  dl <- data.create2(AnalyzeCaste1 = "Soldier",
                     AnalyzeCaste2 = "Worker")
  HWdata <- dl[[2]]
  HW.tree <- dl[[1]]
  HW.tree <- multi2di(HW.tree, random=F)
  {
    HWdata$Subfamily <- as.character(HWdata$Subfamily)
    HWdata[HWdata$Subfamily=="Stylotermitidae","Subfamily"] <- "Rhinotermitidae"
    HWdata[HWdata$Subfamily=="Stolotermitidae","Subfamily"] <- "Archotermopsidae"
    HWdata[HWdata$Subfamily=="Mastotermitidae","Subfamily"] <- "Archotermopsidae"
    HWdata[HWdata$Subfamily=="Hodotermitidae","Subfamily"] <- "Archotermopsidae"
    HWdata[HWdata$Subfamily=="Foraminitermitinae","Subfamily"] <- "Macrotermitinae"
    HWdata[HWdata$Subfamily=="Sphaerotermitinae","Subfamily"] <- "Macrotermitinae"
    HWdata[HWdata$Subfamily=="Sphaerotermitinae","Subfamily"] <- "Macrotermitinae"
    HWdata[HWdata$Subfamily=="Archotermopsidae","Subfamily"] <- "Kalotermitidae"
    HWdata[HWdata$Subfamily=="Kalotermitidae","Subfamily"] <- "non.Neoisoptera"
    names <- c("non.Neoisoptera", "Rhinotermitidae", "Macrotermitinae",
               "Apicotermitinae", "Syntermitinae", "Nasutitermitinae",
               "Cubitermitinae", "Termitinae")
    HWdata$Subfamily <- factor(HWdata$Subfamily, levels=rev(names) )
  }
  comp.data<-comparative.data(HW.tree, HWdata, names.col="new.label", vcv=T, 
                              vcv.dim=3, warn.dropped=TRUE)
  model1<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data)
  model2<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML")
  model3<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, kappa="ML")
  model4<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, delta="ML")
  model5<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML", kappa="ML", delta="ML")
  
  AIC(model1, model2, model3, model4, model5)
  
  summary(model5)
  anova(model5)
  
  ggplot(Res.Living.plot, aes(x=Worker, y=Soldier, col=Sub.Family)) +
    geom_point(size=0.5) + scale_color_manual(values = rev(cols)) +
    coord_fixed() +
    geom_point(data=HWdata, aes(x=HW.mean2, y=HW.mean1, col=Subfamily), pch=3)+
    geom_abline(slope=1, intercept = 0) +
    geom_abline(slope=model4$model$coef[2], intercept = model4$model$coef[1], col=2)+
    scale_x_log10(limit=c(0.4,7)) + scale_y_log10(limit=c(0.4,7)) + 
    theme_classic() +
    theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank())+
    annotation_logticks()
  ggsave(paste0(today, "-Soldier-Worker.pdf"), width = 3, height=4, family="PT Sans", paper="a4")
  
}

## Soldier-Imago
{
  ## count species
  {
    sum(Res.Living[Res.Living$TrueWorker==0,]$Imago <= Res.Living[Res.Living$TrueWorker==0,]$Soldier, na.rm=T)
    sum(Res.Living[Res.Living$TrueWorker==0,]$Imago > Res.Living[Res.Living$TrueWorker==0,]$Soldier, na.rm=T)
    119/129
    
    sum(Res.Living[Res.Living$TrueWorker==1,]$Imago <= Res.Living[Res.Living$TrueWorker==1,]$Soldier, na.rm=T)
    sum(Res.Living[Res.Living$TrueWorker==1,]$Imago > Res.Living[Res.Living$TrueWorker==1,]$Soldier, na.rm=T)
    280/607
  }
  
  ## correlation
  {
    dl <- data.create2(AnalyzeCaste1 = "Soldier",
                       AnalyzeCaste2 = "Imago")
    HWdata <- dl[[2]]
    HW.tree <- dl[[1]]
    HW.tree <- multi2di(HW.tree, random=F)
    {
      HWdata$Subfamily <- as.character(HWdata$Subfamily)
      HWdata[HWdata$Subfamily=="Stylotermitidae","Subfamily"] <- "Rhinotermitidae"
      HWdata[HWdata$Subfamily=="Stolotermitidae","Subfamily"] <- "Archotermopsidae"
      HWdata[HWdata$Subfamily=="Mastotermitidae","Subfamily"] <- "Archotermopsidae"
      HWdata[HWdata$Subfamily=="Hodotermitidae","Subfamily"] <- "Archotermopsidae"
      HWdata[HWdata$Subfamily=="Foraminitermitinae","Subfamily"] <- "Macrotermitinae"
      HWdata[HWdata$Subfamily=="Sphaerotermitinae","Subfamily"] <- "Macrotermitinae"
      HWdata[HWdata$Subfamily=="Sphaerotermitinae","Subfamily"] <- "Macrotermitinae"
      HWdata[HWdata$Subfamily=="Archotermopsidae","Subfamily"] <- "Kalotermitidae"
      HWdata[HWdata$Subfamily=="Kalotermitidae","Subfamily"] <- "non.Neoisoptera"
      HWdata[HWdata$Subfamily=="Serritermitidae","Subfamily"] <- "Rhinotermitidae"
      names <- c("non.Neoisoptera", "Rhinotermitidae", "Macrotermitinae",
                 "Apicotermitinae", "Syntermitinae", "Nasutitermitinae",
                 "Cubitermitinae", "Termitinae")
      HWdata$Subfamily <- factor(HWdata$Subfamily, levels=rev(names) )
    }
    comp.data<-comparative.data(HW.tree, HWdata, names.col="new.label")
    model1<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data)
    model2<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML")
    model3<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, kappa="ML")
    model4<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, delta="ML")
    model5<-pgls(log10(HW.mean1)~log10(HW.mean2), data=comp.data, lambda="ML", kappa="ML", delta="ML")
    
    AIC(model1, model2, model3, model4, model5)
    
    summary(model5)
    anova(model5)
    
    ggplot(Res.Living.plot, aes(x=Imago, y=Soldier, col=Sub.Family)) +
      geom_point(size=0.5) + scale_color_manual(values = rev(cols)) +
      coord_fixed() +
      geom_point(data=HWdata, aes(x=HW.mean2, y=HW.mean1, col=Subfamily), pch=3)+
      geom_abline(slope=1, intercept = 0) +
      geom_abline(slope=model4$model$coef[2], intercept = model4$model$coef[1], col=2)+
      scale_x_log10(limit=c(0.3,7)) + scale_y_log10(limit=c(0.3,7)) + 
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank())+
      annotation_logticks()
    ggsave(paste0(today, "-ImagoSoldier.pdf"), width = 3, height=4, family="PT Sans", paper="a4")
    
  }
  
  ## comparison working caste type
  {
    HWdata$pdisp <- (HWdata$HW.mean2-HWdata$HW.mean1)/HWdata$HW.mean1 
    HWdata$TrueWorker <- factor(HWdata$TrueWorker )
    comp.data<-comparative.data(HW.tree, HWdata, names.col="new.label", vcv=T, 
                                vcv.dim=3, warn.dropped=TRUE)
    model1<-pgls(pdisp~TrueWorker, data=comp.data)
    model2<-pgls(pdisp~TrueWorker, data=comp.data, lambda="ML")
    model3<-pgls(pdisp~TrueWorker, data=comp.data, kappa="ML")
    model4<-pgls(pdisp~TrueWorker, data=comp.data, delta="ML")
    model5<-pgls(pdisp~TrueWorker, data=comp.data, lambda="ML", kappa="ML", delta="ML")
    AIC(model1, model2, model3, model4, model5)
    summary(model5)
    anova(model5)
    
    HWdata$TrueWorker <- factor(HWdata$TrueWorker, levels = c("1", "0"))
    comp.data<-comparative.data(HW.tree, HWdata, names.col="new.label", vcv=T, 
                                vcv.dim=3, warn.dropped=TRUE)
    model5<-pgls(pdisp~TrueWorker, data=comp.data, lambda="ML", kappa="ML", delta="ML")
    summary(model5)
    
    
    ggplot(Res.Living.plot, aes(x=TrueWorker, y=pdisp, fill=as.factor(TrueWorker))) +
      geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.02, dotsize = 1.5, alpha = .5)+
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                   geom="pointrange", color="red") +
      scale_fill_manual(values = cols) +
      theme_bw() + ylab("(Imago - Soldier) / Soldier") +
      theme(legend.position = "none", aspect.ratio = 1)
    ggsave(paste0(today, "_TrueWorker-Soldier-Imago.pdf"),
           height=3, width=3, family="PT Sans", paper="a4")
    
    bartlett.test(pdisp~TrueWorker, data=Res.Living.plot)
  }  
  
  ## ancestral state reconstruction
  {
    dl <- data.create2(AnalyzeCaste1 = "Soldier",
                       AnalyzeCaste2 = "Imago")
    HWdata <- dl[[2]]
    HW.tree <- dl[[1]]
    HW.tree <- multi2di(HW.tree, random=F)
    
    TrueWorker <- HWdata$TrueWorker
    names(TrueWorker) <- HW.tree$tip.label
    TrueWorker <- factor(TrueWorker, levels = c(0,1))
    
    mtrees <- make.simmap(HW.tree, TrueWorker, nsim=100)
    pd<-summary(mtrees,plot=FALSE)
    
    
    is_tip <- HW.tree$edge[,2] <= length(HW.tree$tip.label)
    ordered_tips <- HW.tree$edge[is_tip, 2]
    
    cols = viridis(2, end=0.5)
    names(cols) <- c("0", "1")
    {
      pdf(paste0(today,"-Phylogeny-TrueWorker-sol-imago.pdf"), width=5, height=7,
          family = "PT Serif", paper = "a4")
      layout(matrix(c(1,2),1,2),c(0.2,0.2))
      plot(mtrees[1], mar=c(5, 2, 3, 0), fsize=0.4, colors=cols, 
           xlim=c(-5,150), ftype="i")
      nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
      tiplabels(pie=to.matrix(TrueWorker,sort(unique(TrueWorker))),piecol=cols[c(1,3,2)],cex=0.4)
      axisPhylo(las=1)
      col2 <- viridis(2, option="inferno",alpha = 0.7)
      XonFig <- barplot((HWdata$HW.mean2/HWdata$HW.mean1)[ordered_tips],
                        horiz=TRUE,width=0.5,space=1,
                        ylim=c(1,length(HW.tree$tip.label))-0.5, names="",
                        xlim=c(0,3), col = 
                          cols[as.numeric(as.factor(HWdata$TrueWorker))[ordered_tips]])
      abline(v=c(0,1,2))
      dev.off()
    }
  }
}
  