### Load packages and simple functions

# The following packages are needed to run the script

library(reshape)
library(ggplot2)
library(plyr)
library(dplyr)
library(Hmisc)
library(EloRating)

summarise_bootstraps <- function(bootstrap, report=TRUE){
  sb <- data.frame(driver = c("climate_change", "direct_exploitation", "ias", "land_sea_use", "pollution", "other"),
                   mean=rep(NA, 6), lower_ci=rep(NA, 6), upper_ci=rep(NA, 6))
  upper <- ifelse(exclude.other.drivers == "Yes", 5, 6)
  for (i in 1:upper){
    sb$mean[i] <- mean(bootstrap[[i]])
    sb$lower_ci[i] <- quantile(bootstrap[[i]], 0.025)
    sb$upper_ci[i] <- quantile(bootstrap[[i]], 0.975)
    if (report==TRUE){
      cat(paste(sb$driver[i]),": mean = ", round(sb$mean[i],3), "; 95% CI:", round(sb$lower_ci[i], 3), " - ", round(sb$upper_ci[i], 3), "\\n", sep="")
    }
  }
  if (report==TRUE) cat("\\n")
  sb <- sb[1:upper,]
  return(sb)
}

write.results.csv <- function(observed.david.score=observed.david.score, sb=sb, ps=ps, filestem){
  sb <- sb[1:5,]
  p.boot <- ps$uncorrected.p
  p.adj <- ps$adjusted.p
  
  results <- data.frame(Driver = sort(observed.david.score$david_scores$ID),
                        David_score = observed.david.score$david_scores$normDS[order(observed.david.score$david_scores$ID)],
                        DS0.025 = sb$lower_ci[1:5],
                        DS0.975 = sb$upper_ci[1:5],
                        p_diff_CC = c(NA, p.boot[1:4]),
                        p_diff_DE = c(p.boot[1], NA, p.boot[5:7]),
                        p_diff_IAS = c(p.boot[c(2, 5)], NA, c(p.boot[8:9])),
                        p_diff_LU = c(p.boot[c(3, 6, 8)], NA, p.boot[10]),
                        p_diff_PO = c(p.boot[c(4, 7, 9, 10)], NA),
                        padj_diff_CC = c(NA, p.adj[1:4]),
                        padj_diff_DE = c(p.adj[1], NA, p.adj[5:7]),
                        padj_diff_IAS = c(p.adj[c(2, 5)], NA, c(p.adj[8:9])),
                        padj_diff_LU = c(p.adj[c(3, 6, 8)], NA, p.adj[10]),
                        padj_diff_PO = c(p.adj[c(4, 7, 9, 10)], NA))
  filename <- paste("output/DSresults_", make.names(paste(indicator.reweights, "_", sub(" ", "", filestem),
                                                "_", Sys.time(), ".csv", sep="")))
  write.csv(results, file=filename)
  print(paste("Results written to", filename))
  return(results)
}


