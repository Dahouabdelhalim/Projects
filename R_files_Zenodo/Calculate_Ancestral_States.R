
# CALCULATE ANCESTRAL STATES (CHRIS COONEY - 2017):

rm(list=ls())

library(ape)
library(BAMMtools)
library(parallel)
library(phytools)
library(geiger)
library(caper)
library(hexbin)
library(RColorBrewer)

setwd("~/Dropbox/Multivariate_disparity_through_time/")

# TREE:
raw.tree <- read.tree("Data/Raw_tree.tre")

# DATA:
scores <- read.csv("Data/Traits.csv", strings=F)
rownames(scores) <- scores$TipLabel

# Add branching times to raw tree:
raw.tree <- BAMMtools:::getStartStopTimes(raw.tree)


# ------------------------------------------------------------------------------------ #


traits <- c("trait1","trait2","trait3") # add names of trait variables here


# ------------------------------------------------------------------------------------ #

# Funcion to estimate ancestral states from rate scaled trees

estimate.ancestral.states <- function (i) {
  x <- focalscores
  x <- c(x, fastAnc(focaltrees[[i]], x)) # Insert RATETREE here - use phytools fastAnc to calculate ancestral states using rate scaled tree
  x[1:length(raw.tree$tip.label)] <- x[raw.tree$tip.label]
  names(x)[1:length(raw.tree$tip.label)] <- 1:length(raw.tree$tip.label)
  X <- matrix(x[as.character(raw.tree$edge)], nrow(raw.tree$edge), ncol(raw.tree$edge))
  return(X) # Produces a matrix equivalent to 'edge' in tree object with ancestral states of ancestor / descendent
}

# ------------------------------------------------------------------------------------ #


# Estimate Ancestral States - loop through traits axes:

ntrees <- 1 # Number of trees to estimate anc states for (1 = only mean rate scaled tree)

for (reps in 1:length(traits)) {
  cat("\\r", reps, "of", length(traits))
  scaledTree <- read.tree(paste("Data/Scaled_tree_", traits[reps], ".tre", sep="")) # Read in mean rate scaled tree
  scaledTrees <- c() # Read in posterior distribution of rate scaled trees (if applicable)
  allTrees <- as.list(rep(NA, length(scaledTrees) + 1))
  allTrees[[1]] <- scaledTree
  # for (j in 1:length(scaledTrees)) {
  #   allTrees[[j+1]] <- scaledTrees[[j]]
  # }
  focaltrees <- allTrees[1:ntrees]
  focalscores <- scores[,traits[reps]]; names(focalscores) <- scores$TipLabel
  anc.states <- mclapply(c(1:length(focaltrees)), estimate.ancestral.states, mc.cores = 1) # Run in parallel - adjust cores if necessary
  saveRDS(anc.states, paste("Outputs/", traits[reps], "_ancestral_states.rds", sep=""))
}


# ------------------------------------------------------------------------------------ #


# Calculate Disparity Through Time (indicating which are Meliphagidae lineages):

H <- nodeHeights(raw.tree)
start <- min(H)
end <- max(H)
Bt <- seq(start, floor(end), by = 0.5) # Set time slices through tree here

ntrees <- 1 # i.e. 1 = mean scaled tree

# Loop through each trait (and trees if applicable):
# Outputs a file with the estimated trait values of all the lineages alive at each chosen timeslice:
for (tree in 1:ntrees) {
  for (reps in 1:length(traits)) {
    cat("\\r", reps, "of", length(traits))
    anc.states <- readRDS(paste("Outputs/", traits[reps], "_ancestral_states.rds", sep=""))
    A <- anc.states[[tree]]
    cols <- c("t","traitval")
    write(cols, paste("Outputs/Res/", tree, "_", traits[reps], "_ancestral_states_through_time.csv", sep=""), ncolumns=3, sep=",")
    for (i in 1:length(Bt)) {
      fBt <- Bt[i]
      test1 <- H[,1] <= fBt # Branch begins before or at timeslice
      test2 <- H[,2] > fBt # Branch ends after timeslice
      fbranches <- which(test1 & test2)
      for (j in 1:length(fbranches)) {
        fbranch <- fbranches[j]
        if (H[fbranch,1] == fBt) {
          traitval <- A[fbranch,1]
          fres <-  c(fBt, traitval)
          write(fres, paste("Outputs/Res/", tree, "_", traits[reps], "_ancestral_states_through_time.csv", sep=""), ncolumns=3, sep=",", append=T)
        } else {
          timediff <- fBt - H[fbranch,1]
          timeprop <- timediff / diff(H[fbranch,])
          traitval <- A[fbranch,1] + (diff(A[fbranch,]) * timeprop)
          fres <-  c(fBt, traitval)
          write(fres, paste("Outputs/Res/", tree, "_", traits[reps], "_ancestral_states_through_time.csv", sep=""), ncolumns=3, sep=",", append=T)
        }
      }
    }
  }
}


# ------------------------------------------------------------------------------------ #


# Summarise Disparity Through Time:

ntrees <- 1
sumResList <- as.list(rep(NA, ntrees))

# Function to calculate disparity (variance or sum of variances) for each trait at each time slice:

summarise.disparity.through.time.all <- function(tree) {
  
  withinTreeRes <- as.list(rep(NA, length(traits)+2))
  names(withinTreeRes) <- c(traits, "MULTIVAR")
  for (reps in 1:length(traits)) {
    res <- read.csv(paste("Outputs/Res/", tree, "_", traits[reps], "_ancestral_states_through_time.csv", sep=""))
    
    # Univariate
    sumRes <- data.frame()
    for (i in 1:length(unique(res$t))) {
      fdata <- res[res$t == unique(res$t)[i], -1]
      sumvar <- var(fdata) # Variance for univariate traits
      nlin <- length(fdata)
      sumRes <- rbind(sumRes, c(unique(res$t)[i], sumvar, nlin))
    }
    colnames(sumRes) <- c("t", "sumvar", "nlin")
    sumRes <- rbind(sumRes, c(max(branching.times(raw.tree)),
                              var(scores[,traits[reps]]),
                              length(scores[,traits[reps]])
                              )
                    ) # Add present-day disparity
    withinTreeRes[[traits[reps]]] <- sumRes
    if (reps == 1) {
      allRes <- res
    } else {
      allRes <- cbind(allRes, res[,2])
    }
  }
  colnames(allRes) <- c("t", traits)
  
  # Multivariate
  sumRes <- data.frame()
  for (i in 1:length(unique(allRes$t))) {
    fdata <- allRes[allRes$t == unique(allRes$t)[i], -1]
    sumvar <- sum(diag(cov(fdata))) # Sum of variances for multivariate traits
    nlin <- nrow(fdata)
    sumRes <- rbind(sumRes, c(unique(allRes$t)[i], sumvar, nlin))
  }
  colnames(sumRes) <- c("t", "sumvar", "nlin")
  sumRes <- rbind(sumRes, c(max(branching.times(raw.tree)),
                            sum(diag(cov(scores[,traits]))),
                            nrow(scores[,traits])
                            )
                  ) # Add present-day disparity
  withinTreeRes[["MULTIVAR"]] <- sumRes
  
  return(withinTreeRes) # Produces a list of the disparity at each times slice for each trait (and number of alive lineages)
  
}

sumResList <- mclapply(c(1:ntrees), summarise.disparity.through.time.all, mc.cores = 1) # Run in parallel - adjust cores if necessary

saveRDS(sumResList, "Outputs/Empirical_DTT_data.rds")


# ------------------------------------------------------------------------------------ #


# Plot Ancestral States Through Time:

# Plot using heatmap:

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
for (reps in 1:length(traits)) {
  res <- read.csv(paste("Outputs/Res/", tree, "_", traits[reps], "_ancestral_states_through_time.csv", sep=""))
  if (reps == 1) {
    allRes <- res
  } else {
    allRes <- cbind(allRes, res[,2])
  }
}
colnames(allRes) <- c("t", traits)
h1 <- hexbinplot(trait1 ~ t, data=allRes, colramp=rf, trans=log, inv=exp, asp=1)
h2 <- hexbinplot(trait2 ~ t, data=allRes, colramp=rf, trans=log, inv=exp, asp=1)
h3 <- hexbinplot(trait3 ~ t, data=allRes, colramp=rf, trans=log, inv=exp, asp=1)
maxcount <- max(c(h1$panel.args.common$maxcnt,
                  h2$panel.args.common$maxcnt,
                  h3$panel.args.common$maxcnt))
pdf("Plots/Ancestral_states_heatmap_through_time.pdf", width = 5, height = 5)
    hexbinplot(trait1 ~ t, data=allRes, colramp=rf, trans=log, inv=exp, asp=1, cex.lab=0.75, cex.title=0.8) # maxcnt=maxcount
    hexbinplot(trait2 ~ t, data=allRes, colramp=rf, trans=log, inv=exp, asp=1, cex.lab=0.75, cex.title=0.8) # maxcnt=maxcount
    hexbinplot(trait3 ~ t, data=allRes, colramp=rf, trans=log, inv=exp, asp=1, cex.lab=0.75, cex.title=0.8) # maxcnt=maxcount
dev.off()


# ------------------------------------------------------------------------------------ #


# Plot Disparity Metrics States Through Time:


pdf("Plots/Ancestral_disparity_metrics_against_time.pdf", width = (length(traits)+1)*2, height = 2)

  par(mfcol=c(1,length(traits)+1))
  par(mar=c(2.5,4,2.5,1))
  par(cex.axis=0.9)
  par(cex.lab=1.25)
  
  plot(sumvar ~ t, sumResList[[1]]$MULTIVAR, type="l", lwd=2, log="", xlab="", ylab="Disparity"); box(lwd=2); title(main = "Multivariate")
  for (reps in 1:length(traits)) {
    plot(sumvar ~ t, sumResList[[1]][[traits[reps]]], type="l", lwd=2, log="", xlab="", ylab=""); box(lwd=2); title(main = traits[reps])
  }

dev.off()


# ------------------------------------------------------------------------------------ #


# SIMULATIONS


# Simulate BM Disparity Through Time:

nsims <- 50 # Set number of simulations

# Function to estimate ancestral states from simulated data using raw tree (i.e. BM model)
estimate.simulated.ancestral.states.bm <- function (i) {
  x <- simdat[,i]
  x <- c(x, fastAnc(raw.tree, x))
  x[1:length(raw.tree$tip.label)] <- x[raw.tree$tip.label]
  names(x)[1:length(raw.tree$tip.label)] <- 1:length(raw.tree$tip.label)
  X <- matrix(x[as.character(raw.tree$edge)], nrow(raw.tree$edge), ncol(raw.tree$edge))
  return(X)
}

# Function to estimate ancestral states from simulated data using rate scaled tree (i.e. VR model)
estimate.simulated.ancestral.states.vr <- function (i) {
  x <- simdat[,i]
  x <- c(x, fastAnc(scaledTree, x))
  x[1:length(raw.tree$tip.label)] <- x[raw.tree$tip.label]
  names(x)[1:length(raw.tree$tip.label)] <- 1:length(raw.tree$tip.label)
  X <- matrix(x[as.character(raw.tree$edge)], nrow(raw.tree$edge), ncol(raw.tree$edge))
  return(X)
}

# BM model - simulate traits under BM and estimate ancestral states:
for (reps in 1:length(traits)) {
  dat <- scores[,traits[reps]]; names(dat) <- scores$TipLabel
  fit <- fitContinuous(raw.tree, dat, model = "BM")
  simdat <- fastBM(raw.tree, a = fit$opt$z0, sig2 = fit$opt$sigsq, nsim = nsims) # Simulate data using parameters from BM model fitted to the raw data in each axis
  sim.anc.states <- mclapply(c(1:nsims), estimate.simulated.ancestral.states.bm, mc.cores = 20) # Run in parallel - adjust cores if necessary
  saveRDS(simdat, paste("Outputs/Data_sim/", traits[reps], "_sim_tip_states_BM.rds", sep=""))
  saveRDS(sim.anc.states, paste("Outputs/Data_sim/", traits[reps], "_sim_ancestral_states_BM.rds", sep=""))
}

# Variable rates model:
for (reps in 1:length(traits)) {
  dat <- scores[,traits[reps]]; names(dat) <- scores$TipLabel
  scaledTree <- read.tree(paste("Data/Scaled_tree_", traits[reps], ".tre", sep="")) # Read in rate scaled tree
  fit <- fitContinuous(scaledTree, dat, model = "BM") # Fit BM model on rate scaled tree
  simdat <- fastBM(scaledTree, a = fit$opt$z0, sig2 = fit$opt$sigsq, nsim = nsims) # Simulate on rate scaled tree using simulated parameters
  sim.anc.states <- mclapply(c(1:nsims), estimate.simulated.ancestral.states.vr, mc.cores = 20)
  saveRDS(simdat, paste("Outputs/Data_sim/", traits[reps], "_sim_tip_states_VR.rds", sep=""))
  saveRDS(sim.anc.states, paste("Outputs/Data_sim/", traits[reps], "_sim_ancestral_states_VR.rds", sep=""))
}


# Calculate simulated disparity through time (this part is quite slow):

H <- nodeHeights(raw.tree)
start <- min(H)
end <- max(H)
Bt <- seq(start, floor(end), by = 0.5) # Set timeslices (use same as empirical)

mod.types <- c("BM","VR")

# Loop through simulations 
for (sim in 1:nsims) {
  for (reps in 1:length(traits)) {
    for (type in 1:length(mod.types)) {
      cat("\\r", sim, "of", nsims)
      sim.anc.states <- readRDS(paste("Outputs/Data_sim/", traits[reps], "_sim_ancestral_states_", mod.types[type], ".rds", sep=""))
      A <- sim.anc.states[[sim]]
      cols <- c("t","traitval")
      write(cols, paste("Outputs/Res_sim/", sim, "_", traits[reps], "_sim_ancestral_states_through_time_", mod.types[type], ".csv", sep=""), ncolumns=2, sep=",")
      for (i in 1:length(Bt)) {
        #print(i)
        fBt <- Bt[i]
        test1 <- H[,1] <= fBt # Branch begins before or at timeslice
        test2 <- H[,2] > fBt # Branch ends after timeslice
        fbranches <- which(test1 & test2)
        for (j in 1:length(fbranches)) {
          fbranch <- fbranches[j]
          if (H[fbranch,1] == fBt) {
            traitval <- A[fbranch,1]
            fres <-  c(fBt, traitval)
            write(fres, paste("Outputs/Res_sim/", sim, "_", traits[reps], "_sim_ancestral_states_through_time_", mod.types[type], ".csv", sep=""), ncolumns=2, sep=",", append=T)
          } else {
            timediff <- fBt - H[fbranch,1]
            timeprop <- timediff / diff(H[fbranch,])
            traitval <- A[fbranch,1] + (diff(A[fbranch,]) * timeprop)
            fres <-  c(fBt, traitval)
            write(fres, paste("Outputs/Res_sim/", sim, "_", traits[reps], "_sim_ancestral_states_through_time_", mod.types[type], ".csv", sep=""), ncolumns=2, sep=",", append=T)
          }
        }
      }
    }
  }
}


# Summarise Disparity Through Time:


summarise.simulated.disparity.through.time <- function(sim, type = NA) {
  
  withinTreeRes <- as.list(rep(NA, length(traits)+2))
  names(withinTreeRes) <- c(traits, "MULTIVAR")
  
  # Collate info on tips and ancestral states
  for (reps in 1:length(traits)) {
    tips <- readRDS(paste("Outputs/Data_sim/", traits[reps], "_sim_tip_states_", type, ".rds", sep=""))
    res <- read.csv(paste("Outputs/Res_sim/", sim, "_", traits[reps], "_sim_ancestral_states_through_time_", type, ".csv", sep=""))
    if (reps == 1) {
      allTips <- tips[,sim,drop=F]
      allRes <- res
    } else {
      allTips <- cbind(allTips, tips[,sim])
      allRes <- cbind(allRes, res[,2])
    }
  }
  colnames(allRes) <- c("t", traits)
  
  # Univariate disparity:
  for (reps in 1:length(traits)) {
    tips <- readRDS(paste("Outputs/Data_sim/", traits[reps], "_sim_tip_states_", type, ".rds", sep=""))
    res <- read.csv(paste("Outputs/Res_sim/", sim, "_", traits[reps], "_sim_ancestral_states_through_time_", type, ".csv", sep=""))
    sumRes <- data.frame()
    for (i in 1:length(unique(res$t))) {
      fdata <- res[res$t == unique(res$t)[i], -1]
      sumvar <- var(fdata)
      nlin <- length(fdata)
      sumRes <- rbind(sumRes, c(unique(res$t)[i], sumvar, nlin))
    }
    colnames(sumRes) <- c("t", "sumvar", "nlin")
    sumRes <- rbind(sumRes, c(max(branching.times(raw.tree)),
                              var(tips[,sim]),
                              length(tips[,sim])
                              )
                    ) # Add present-day disparity
    withinTreeRes[[traits[reps]]] <- sumRes
  }
  
  # Multivariate disparity:
  sumRes <- data.frame()
  for (i in 1:length(unique(allRes$t))) {
    fdata <- allRes[allRes$t == unique(allRes$t)[i], -1]
    sumvar <- sum(diag(cov(fdata)))
    nlin <- nrow(fdata)
    sumRes <- rbind(sumRes, c(unique(allRes$t)[i], sumvar, nlin))
  }
  colnames(sumRes) <- c("t", "sumvar", "nlin")
  sumRes <- rbind(sumRes, c(max(branching.times(raw.tree)),
                            sum(diag(cov(allTips))),
                            nrow(allTips)
                            )
                  ) # Add present-day disparity
  withinTreeRes[["MULTIVAR"]] <- sumRes
  
  return(withinTreeRes)
  
}

sumResListSim_BM <- mclapply(c(1:nsims), summarise.simulated.disparity.through.time, type="BM", mc.cores = 1) # Run in parallel -adjust cores if necessary
sumResListSim_VR <- mclapply(c(1:nsims), summarise.simulated.disparity.through.time, type="VR", mc.cores = 1)

saveRDS(sumResListSim_BM, "Outputs/Simulated_DTT_data_BM.rds")
saveRDS(sumResListSim_VR, "Outputs/Simulated_DTT_data_VR.rds")


# ----------------- #


# Compare estimated and and simulated disparity through time:

library(MCMCglmm)

pdf("Plots/Ancestral_disparity_metrics_against_time_with_BM_null.pdf", width = (length(traits)+1)*2, height = 2)
  
  par(mfcol=c(1,length(traits)+1))
  par(mar=c(2.5,4,2.5,1))
  par(cex.axis=0.9)
  par(cex.lab=1.25)
  
  # Multivariate:
  ts <- sumResList[[1]]$MULTIVAR$t
  sumvar.vals <- c()
  for (i in 1:length(sumResListSim_BM)) {
    sumvar.vals <- rbind(sumvar.vals, sumResListSim_BM[[i]]$MULTIVAR$sumvar / max(sumResListSim_BM[[i]]$MULTIVAR$sumvar)) # Scale between 0 - 1
  }
  sumvar.lci <- HPDinterval(as.mcmc(sumvar.vals))[,1]; sumvar.uci <- HPDinterval(as.mcmc(sumvar.vals))[,2]; sumvar.sim.mean <- colMeans(sumvar.vals) # Calculate confidence intervals
  plot(1, type="n", xlab="t", ylab="Disparity", xlim=c(0,max(ts)), ylim=c(0,1)); box(lwd=2); title(main = "Multivariate")
  polygon(c(ts,rev(ts)), y = c(sumvar.lci,rev(sumvar.uci)), col = "grey90", border = "white")
  lines(I(sumvar / max(sumvar)) ~ t, sumResList[[1]]$MULTIVAR, lwd=2, col="black")
  
  # Unvariate disparity:
  for (reps in 1:length(traits)) {
    sumvar.vals <- c()
    for (i in 1:length(sumResListSim_BM)) {
      sumvar.vals <- rbind(sumvar.vals, sumResListSim_BM[[i]][[traits[reps]]]$sumvar / max(sumResListSim_BM[[i]][[traits[reps]]]$sumvar))
    }
    sumvar.lci <- HPDinterval(as.mcmc(sumvar.vals))[,1]; sumvar.uci <- HPDinterval(as.mcmc(sumvar.vals))[,2]; sumvar.sim.mean <- colMeans(sumvar.vals)
    plot(1, type="n", xlab="t", ylab="", xlim=c(0,max(ts)), ylim=c(0,1)); box(lwd=2); title(main = traits[reps])
    polygon(c(ts,rev(ts)), y = c(sumvar.lci,rev(sumvar.uci)), col = "grey90", border = "white")
    lines(I(sumvar / max(sumvar)) ~ t, sumResList[[1]][[traits[reps]]], lwd=2, col="black")
  }

dev.off()


pdf("Plots/Ancestral_disparity_metrics_against_time_with_VR_null.pdf", width = (length(traits)+1)*2, height = 2)

  par(mfcol=c(1,length(traits)+1))
  par(mar=c(2.5,4,2.5,1))
  par(cex.axis=0.9)
  par(cex.lab=1.25)
  
  # Multivariate
  ts <- sumResList[[1]]$MULTIVAR$t
  sumvar.vals <- c()
  for (i in 1:length(sumResListSim_VR)) {
    sumvar.vals <- rbind(sumvar.vals, sumResListSim_VR[[i]]$MULTIVAR$sumvar / max(sumResListSim_VR[[i]]$MULTIVAR$sumvar))
  }
  sumvar.lci <- HPDinterval(as.mcmc(sumvar.vals))[,1]; sumvar.uci <- HPDinterval(as.mcmc(sumvar.vals))[,2]; sumvar.sim.mean <- colMeans(sumvar.vals)
  plot(1, type="n", xlab="t", ylab="Disparity", xlim=c(0,max(ts)), ylim=c(0,1)); box(lwd=2); title(main = "Multivar")
  polygon(c(ts,rev(ts)), y = c(sumvar.lci,rev(sumvar.uci)), col = "grey90", border = "white")
  lines(I(sumvar / max(sumvar)) ~ t, sumResList[[1]]$MULTIVAR, lwd=2, col="black")
  
  # Unvariate disparity:
  for (reps in 1:length(traits)) {
    sumvar.vals <- c()
    for (i in 1:length(sumResListSim_VR)) {
      sumvar.vals <- rbind(sumvar.vals, sumResListSim_VR[[i]][[traits[reps]]]$sumvar / max(sumResListSim_VR[[i]][[traits[reps]]]$sumvar))
    }
    sumvar.lci <- HPDinterval(as.mcmc(sumvar.vals))[,1]; sumvar.uci <- HPDinterval(as.mcmc(sumvar.vals))[,2]; sumvar.sim.mean <- colMeans(sumvar.vals)
    plot(1, type="n", xlab="t", ylab="", xlim=c(0,max(ts)), ylim=c(0,1)); box(lwd=2); title(main = traits[reps])
    polygon(c(ts,rev(ts)), y = c(sumvar.lci,rev(sumvar.uci)), col = "grey90", border = "white")
    lines(I(sumvar / max(sumvar)) ~ t, sumResList[[1]][[traits[reps]]], lwd=2, col="black")
  }

dev.off()


# --------------------- #


# Replicate Fig 2a,b plot (with additional intermediate plot) for multivariate data (all traits):

pdf("Plots/Disparity_through_time_BM_VR_Multivar.pdf", height=2.4, width=7.2)

  par(mfrow=c(1,3))
  par(mar=c(3.5,3.5,1,1))
  par(cex.axis=0.9)
  par(cex.lab=1)
  
  pal <- brewer.pal(12, "Paired")
  
  # --------- #
  
  # Plot MULTIVAR:
  
  ts <- sumResList[[1]]$MULTIVAR$t
  sumvar.obs <- sumResList[[1]]$MULTIVAR$sumvar # re-scale observed data
  sumvar.obs <- (sumvar.obs - min(sumvar.obs)) / (max(sumvar.obs) - min(sumvar.obs))
  sumvar.vals.BM <- c()
  sumvar.vals.VR <- c()
  for (i in 1:length(sumResListSim_BM)) {
    sumvar.vals.BM <- rbind(sumvar.vals.BM, sumResListSim_BM[[i]]$MULTIVAR$sumvar)
    sumvar.vals.VR <- rbind(sumvar.vals.VR, sumResListSim_VR[[i]]$MULTIVAR$sumvar)
  }
  for (i in 1:length(sumResListSim_BM)) {
    sumvar.vals.BM[i,] <- (sumvar.vals.BM[i,] - min(sumvar.vals.BM[i,])) / (max(sumvar.vals.BM[i,]) - min(sumvar.vals.BM[i,])) # re-scale
    sumvar.vals.VR[i,] <- (sumvar.vals.VR[i,] - min(sumvar.vals.VR[i,])) / (max(sumvar.vals.VR[i,]) - min(sumvar.vals.VR[i,])) # re-scale
  }
  sumvar.sim.BM.mean <- colMeans(sumvar.vals.BM) # find average
  sumvar.sim.VR.mean <- colMeans(sumvar.vals.VR) # find average
  sumvar.obs.mod <- loess(sumvar.obs ~ ts) # Smooth empirical curves
  sumvar.sim.BM.mod <- loess(sumvar.sim.BM.mean ~ ts)
  sumvar.sim.VR.mod <- loess(sumvar.sim.VR.mean ~ ts)
  sumvar.obs.loess <- predict(sumvar.obs.mod, data.frame(ts = ts), se = TRUE)$fit
  sumvar.sim.BM.mean.loess <- predict(sumvar.sim.BM.mod, data.frame(ts = ts), se = TRUE)$fit
  sumvar.sim.VR.mean.loess <- predict(sumvar.sim.VR.mod, data.frame(ts = ts), se = TRUE)$fit
  st <- max(branching.times(raw.tree))
  xvals <- st - ts
  
  # Plot disparity accumulation curves:
  plot(sumvar.obs ~ xvals, xlim=c(max(branching.times(raw.tree)),0), ylim=c(0,1), type="n", lty=1, lwd=2, xlab="", ylab="", axes=F); box(lwd=1.5)
  lines(sumvar.obs.loess ~ xvals)
  lines(sumvar.sim.BM.mean.loess ~ xvals, col=pal[2])
  lines(sumvar.sim.VR.mean.loess ~ xvals, col=pal[6])
  lines(sumvar.obs ~ xvals, lty=1, lwd=2)
  lines(sumvar.sim.BM.mean ~ xvals, lty=2, lwd=2, col=pal[2])
  lines(sumvar.sim.VR.mean ~ xvals, lty=3, lwd=2, col=pal[6])
  axis(1, at=seq(0,100,5), labels=seq(0,100,5), padj=-1); title(xlab="Time (MYA)", line=2)
  axis(2, at=seq(0,1,0.2), labels=seq(0,1,0.2), padj=0.8); title(ylab="Disparity (re-scaled)", line=2)
  legend("bottomright", bty="n", legend=c("Observed", "Null (BM)", "Null (VR)"), lty=c(1,2,3), col=c("black",pal[2],pal[6]), cex=0.75, lwd=1.5)
  
  # --------- #
  
  # Calculate moving slope gradient:
  
  twindow <- 2 # Set the time window to calculate the slope of disparity accumulation
  
  sumvar.slopes.obs <- sumvar.slopes.sim.BM <- sumvar.slopes.sim.VR <- c()
  for (i in 1:ceiling(length(ts)/twindow)) {
    fslices <- (((i-1)*twindow)+1) : (((i-1)*twindow)+twindow)
    fvals.obs <- sumvar.obs.loess[fslices]
    fvals.sim.BM <- sumvar.sim.BM.mean.loess[fslices]
    fvals.sim.VR <- sumvar.sim.VR.mean.loess[fslices]
    fts <- ts[fslices]
    fslope.obs <- coef(lm(fvals.obs ~ ts[fslices]))[2]
    fslope.sim.BM <- coef(lm(fvals.sim.BM ~ ts[fslices]))[2]
    fslope.sim.VR <- coef(lm(fvals.sim.VR ~ ts[fslices]))[2]
    names(fslope.obs) <- mean(fts)
    names(fslope.sim.BM) <- mean(fts)
    names(fslope.sim.VR) <- mean(fts)
    sumvar.slopes.obs <- c(sumvar.slopes.obs, fslope.obs)
    sumvar.slopes.sim.BM <- c(sumvar.slopes.sim.BM, fslope.sim.BM)
    sumvar.slopes.sim.VR <- c(sumvar.slopes.sim.VR, fslope.sim.VR)
  }
  sumvar.slopes.sim.matrix.BM <- sumvar.slopes.sim.matrix.VR <- c()
  for (j in 1:length(sumResListSim_BM)) {
    f.sumvar.sim.BM <- sumvar.vals.BM[j,]
    f.sumvar.sim.VR <- sumvar.vals.VR[j,]
    f.sumvar.sim.mod.BM <- loess(f.sumvar.sim.BM ~ ts)
    f.sumvar.sim.mod.VR <- loess(f.sumvar.sim.VR ~ ts)
    f.sumvar.sim.vals.BM <- predict(f.sumvar.sim.mod.BM, data.frame(ts = ts), se = TRUE)$fit
    f.sumvar.sim.vals.VR <- predict(f.sumvar.sim.mod.VR, data.frame(ts = ts), se = TRUE)$fit
    f.sumvar.slopes.sim.BM <- f.sumvar.slopes.sim.VR <- c()
    for (i in 1:ceiling(length(ts)/twindow)) {
      fslices <- (((i-1)*twindow)+1) : (((i-1)*twindow)+twindow)
      fvals.sim.BM <- f.sumvar.sim.vals.BM[fslices]
      fvals.sim.VR <- f.sumvar.sim.vals.VR[fslices]
      fts <- ts[fslices]
      fslope.sim.BM <- coef(lm(fvals.sim.BM ~ ts[fslices]))[2]
      fslope.sim.VR <- coef(lm(fvals.sim.VR ~ ts[fslices]))[2]
      names(fslope.sim.BM) <- mean(fts)
      names(fslope.sim.VR) <- mean(fts)
      f.sumvar.slopes.sim.BM <- c(f.sumvar.slopes.sim.BM, fslope.sim.BM)
      f.sumvar.slopes.sim.VR <- c(f.sumvar.slopes.sim.VR, fslope.sim.VR)
    }
    sumvar.slopes.sim.matrix.BM <- rbind(sumvar.slopes.sim.matrix.BM, f.sumvar.slopes.sim.BM)
    sumvar.slopes.sim.matrix.VR <- rbind(sumvar.slopes.sim.matrix.VR, f.sumvar.slopes.sim.VR)
  }
  sim.l95.BM <- apply(sumvar.slopes.sim.matrix.BM, 2, quantile, 0.025)
  sim.u95.BM <- apply(sumvar.slopes.sim.matrix.BM, 2, quantile, 0.955)
  sim.l95.VR <- apply(sumvar.slopes.sim.matrix.VR, 2, quantile, 0.025)
  sim.u95.VR <- apply(sumvar.slopes.sim.matrix.VR, 2, quantile, 0.955)
  
  # Plot slope values:
  yrange <- range(c(sumvar.slopes.obs, sumvar.slopes.sim.BM, sim.l95.BM, sim.u95.BM, sumvar.slopes.sim.VR, sim.l95.VR, sim.u95.VR))
  xvals <- c(as.numeric(names(sumvar.slopes.obs)))
  xvals <- st - xvals
  xvals.poly <- c(xvals, rev(xvals))
  plot(sumvar.slopes.obs ~ xvals, type="n", lwd=2, xlim=c(max(branching.times(raw.tree)),0), ylim=yrange, ylab="", xlab="", axes=F); box(lwd=1.5)
  polygon(xvals.poly, c(sim.l95.VR,rev(sim.u95.VR)), col=paste(pal[5],"50",sep=""), border=NA)
  lines(sumvar.slopes.sim.VR ~ xvals, lty=3, lwd=2, col=pal[6])
  polygon(xvals.poly, c(sim.l95.BM,rev(sim.u95.BM)), col=paste(pal[1],"50",sep=""), border=NA)
  lines(sumvar.slopes.sim.BM ~ xvals, lty=2, lwd=2, col=pal[2])
  lines(sumvar.slopes.obs ~ xvals, lty=1, lwd=2)
  axis(1, at=seq(0,100,5), labels=seq(0,100,5), padj=-1); title(xlab="Time (MYA)", line=2)
  axis(2, at=seq(-0.4,0.4,0.02), labels=seq(-0.4,0.4,0.02), padj=0.8); title(ylab="Slope (disparity ~ time)", line=2)
  #legend("topleft", bty="n", legend=c("Observed", "Null (BM)"), lty=c(1,2), cex=0.6, lwd=1.5)
  
  # --------- #
  
  # Plot difference in slope values (observed - simulated)
  
  sumvar.slopes.diff.matrix.BM <- sumvar.slopes.diff.matrix.VR <- c()
  for (j in 1:length(sumResListSim_BM)) {
    fdiff.BM <- sumvar.slopes.obs - sumvar.slopes.sim.matrix.BM[j,]
    fdiff.VR <- sumvar.slopes.obs - sumvar.slopes.sim.matrix.VR[j,]
    sumvar.slopes.diff.matrix.BM <- rbind(sumvar.slopes.diff.matrix.BM, fdiff.BM)
    sumvar.slopes.diff.matrix.VR <- rbind(sumvar.slopes.diff.matrix.VR, fdiff.VR)
  }
  diff.l95.BM <- apply(sumvar.slopes.diff.matrix.BM, 2, quantile, 0.025)
  diff.u95.BM <- apply(sumvar.slopes.diff.matrix.BM, 2, quantile, 0.955)
  diff.l95.VR <- apply(sumvar.slopes.diff.matrix.VR, 2, quantile, 0.025)
  diff.u95.VR <- apply(sumvar.slopes.diff.matrix.VR, 2, quantile, 0.955)
  yrange <- range(c(sumvar.slopes.obs - sumvar.slopes.sim.BM, diff.l95.BM, diff.u95.BM, sumvar.slopes.obs - sumvar.slopes.sim.VR, diff.l95.VR, diff.u95.VR))
  xvals <- c(as.numeric(names(sumvar.slopes.obs)))
  xvals <- st - xvals
  xvals.poly <- c(xvals, rev(xvals))
  
  # Plot difference in slope values:
  plot(sumvar.slopes.obs - sumvar.slopes.sim.BM, type="n", xlim=c(max(branching.times(raw.tree)),0), ylim=yrange, lwd=2, ylab="", xlab="", axes=F); box(lwd=1.5)
  abline(0,0, col="grey80", lty=2)
  polygon(xvals.poly, c(diff.l95.VR,rev(diff.u95.VR)), col=paste(pal[5],"50",sep=""), border=NA)
  lines(sumvar.slopes.obs - sumvar.slopes.sim.VR ~ xvals, lwd=2, lty=3, col=pal[6])
  polygon(xvals.poly, c(diff.l95.BM,rev(diff.u95.BM)), col=paste(pal[1],"50",sep=""), border=NA)
  lines(sumvar.slopes.obs - sumvar.slopes.sim.BM ~ xvals, lwd=2, lty=2, col=pal[2])
  axis(1, at=seq(0,100,5), labels=seq(0,100,5), padj=-1); title(xlab="Time (MYA)", line=2)
  axis(2, at=seq(-0.03,0.03,0.01), labels=seq(-0.03,0.03,0.01), padj=0.8); title(ylab="Difference in slope (observed - null)", line=2)

dev.off()


# ============================================================================== #

