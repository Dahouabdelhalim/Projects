# Function to test whether normalised David's scores from different subsets of the data differ significantly

compare_subsets <- function(pooled, nrands=10000, plot=TRUE, output.each = TRUE, pairwise="Conditional", what="among subsets"){
  # Test for significance of among-subset differences in normalised DS scores
  # Whether all pairs are then compared depends on value of pairwise: "Yes" -> yes; "No" -> no; 
  #   "Conditional" = only if p < 0.05 & n.subsets > 2
  
  # 1. Calculate the observed value of the test statistic
  pooled <- droplevels(pooled) #necessary given how loop is constructed
  
  for (j in 1:nlevels(pooled$source)){
    this <- subset(pooled, source == levels(pooled$source)[j])
    this.DS <- analyse.DS.mat(this, study.drivers, IR=indicator.reweights, summarise.evidence=FALSE, show.result=FALSE)
    this.DS$david_scores$source <- levels(pooled$source)[j]
    if (j == 1){
      this.obs <- this.DS$david_scores
    }else{
      this.obs <- rbind(this.obs, this.DS$david_scores)
    }
    if (output.each == TRUE){
      # Run the bootstrap
      # Still need to add the code
    }
  }
  m1 <- lm(normDS ~ ID, data = this.obs)
  observed.sigma <- summary(m1)$sigma
  
  # 2. Do randomisations
  null.sigma <- rep(NA, nrands)
  shuffled <- pooled
  for (i in 1:nrands){
    shuffled$source <- sample(shuffled$source)
    for (j in 1:nlevels(shuffled$source)){
      this <- subset(shuffled, source == levels(shuffled$source)[j])
      this.DS <- analyse.DS.mat(this, study.drivers, IR=indicator.reweights, summarise.evidence=FALSE, show.result=FALSE)
      if (j == 1){
        this.rand <- this.DS$david_scores
      }else{
        this.rand <- rbind(this.rand, this.DS$david_scores)
      }
    }
    m.rand <- lm(normDS ~ ID, data = this.rand)
    null.sigma[i] <- summary(m.rand)$sigma
  }
  
  # 3. Test
  p.value <- sum(null.sigma >= observed.sigma)/length(null.sigma)
  if (plot==TRUE){
    plot(density(null.sigma), main = paste("Test for sig diff ", what, ": p = ", p.value, sep=""), 
         xlab="Residual standard error",
         xlim=c(min(c(observed.sigma, null.sigma)), max(c(observed.sigma, null.sigma))))
    abline(v=observed.sigma, col="red")
  }
  
  # 4. Do pairwise comparisons among levels if appropriate
  if (pairwise == "Conditional" & p.value < 0.05 & nlevels(pooled$source) > 2) pairwise <- TRUE
  if (pairwise == TRUE){
    # Make all pairwise comparisons
    comparisons <- NULL
    p.values <- NULL
    k <- nlevels(pooled$source)
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        relevant <- levels(pooled$source)[c(i, j)]
        comparison <- paste(levels(pooled$source)[c(i, j)], collapse=" v ")
        print(comparison)
        compdata <- subset(pooled, source %in% relevant)
        pairdiff <- compare_subsets(compdata, what=comparison, nrands=1000) #Not 10,000 here - life's too short
        comparisons <- append(comparisons, comparison)
        p.values <- append(p.values, pairdiff$p.value)
      }
    }
    pairdiffs <- data.frame(comparison=comparisons, uncorrected.p=p.values)
    pairdiffs$corrected.p <- p.adjust(pairdiffs$uncorrected.p, method="BY")
    to.return <- list(matrices=this.obs, sigma=observed.sigma, p.value=p.value, null.sigma=null.sigma, pairwise=pairdiffs)
  }else{
    to.return <- list(matrices=this.obs, sigma=observed.sigma, p.value=p.value, null.sigma=null.sigma)
  }
  
  # 4. Return results
  return(to.return)
}
