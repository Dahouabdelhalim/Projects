################################################################################################################################################################
# This code is for derived analyses of the model
# from 'Hierarchical variation in phenotypic flexibility across timescales and associated survival selection shape the dynamics of partial seasonal migration' 
# by Paul Acker, Francis Daunt, Sarah Wanless, Sarah J. Burthe, Mark A. Newell, Michael P. Harris, Carrie Gunn, Robert Swann, Ana Payo-Payo, and Jane M. Reid.
################################################################################################################################################################


### Load the posterior samples, we will use them in derived calculations below
load('posterior_samples.Rdata')


### Define variables used throughout
N_sex <- 2
N_year <- 9
N_tactic <- 3
tactics <- c('Full-R','Mixed-RM','Full-M')
sex <- c('F','M')
sexes <- c("Females","Males") 
yrlabs <- years <- 2010:2018
fullyrlabs <- paste0(years,'-',11:19)
N_iter <- 5000


### Define useful functions
## Function for computing basic posterior summaries
summarize <- function(x) { 
  output <- round(c(mean(x),quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975)),mean(x>0)),3) 
  names(output) <- c('mean','2.5%','25%','50%','75%','97.5%','Pr(x>0)')
  return(output)
}
## Function for computing cross-year grand mean and standard deviation
crossyr  <- function(x) {
  SUM <- rowSums(x)
  GM <- SUM/dim(x)[2]
  SD <- sqrt((1/dim(x)[2])*rowSums((x-GM)^2))
  out <- list(GM=GM,SD=SD)
  return(out)
}
## Function for computing summaries of cross-year grand mean and standard deviation
crossyrsumm <- function(x) {
  SUM <- rowSums(x)
  GM <- SUM/dim(x)[2]
  SD <- sqrt((1/dim(x)[2])*rowSums((x-GM)^2))
  out <- rbind(summarize(GM),summarize(SD))
  rownames(out) <- c('GM','SD')
  return(out)
}


### Additional summaries of the relative tactic frequencies in the beginning of the year (F_1)
### Results subsection "Tactic frequencies and flexibility within years" and OSM S5.2

## Basic cross-year summaries per sex (provided in main text)
F1_crossyr <- array(dim=c(2,6,N_tactic*2))
dimnames(F1_crossyr) <- list(c('GM','SD'),
                             c('mean','2.5%','25%','50%','75%','97.5%'),
                             paste(rep(tactics,each=2),sexes,sep='::'))
cnt <- 0
for(k in c(1,3,2)) {
  for (s in 1:2) {
    cnt <- cnt+1
    F1_crossyr[,,cnt] <- round(crossyrsumm(F_1[,s,,k])[,-7],2)
  }
}
print(F1_crossyr)

## Detailed summary of between-tactic differences in F_1 (mentioned in main text, provided in OSM S5.2 Table S3)
Table_S3 <- array(dim=c((N_year+2)*2,N_tactic*2+1))
colnames(Table_S3) <- paste(c('sex','FRvsMRM-delta','FRvsMRM-delta-P','FRvsFM-delta','FRvsFM-P','MRMvsFM-delta','MRMvsFM-P'))
rownames(Table_S3) <- rep(c(years,'GM','SD'),each=2)
comp <- cbind(c(1,3),c(1,2),c(3,2))
for (c in 1:3) {
  ic <- c(1,3,5)[c]
  for (s in 1:2) {
    A <- F_1[,s,,comp[1,c]]
    B <- F_1[,s,,comp[2,c]]
    summ <- apply(A-B,2,summarize)
    summ <- cbind(summ,summarize(crossyr(A)[['GM']]-crossyr(B)[['GM']]))
    summ <- cbind(summ,summarize(crossyr(A)[['SD']]-crossyr(B)[['SD']]))
    summ <- round(summ,2)
    for (y in 1:11) {
      ir <- seq(1,22,2)[y]-1+s
      Table_S3[ir,1] <- sex[s]
      Table_S3[ir,ic+1] <- paste0(sprintf("%.2f",summ['mean',y]),' [',sprintf("%.2f",summ['2.5%',y]),',',sprintf("%.2f",summ['97.5%',y]),']')
      Table_S3[ir,ic+2] <- sprintf("%.2f",summ['Pr(x>0)',y])
    }
  }
}
print(Table_S3,quote=F)

## Detailed summary of sex differences in F_1 in each tactic (mentioned in main text, provided in OSM S5.2 Table S4)
Table_S4 <- array(dim=c(N_year+2,N_tactic*2))
colnames(Table_S4) <- paste(c('FullR-delta','FullR-Pr(delta>0)','MixedRM-delta','MixedRM-Pr(delta>0)','FullM-delta','FullM-Pr(delta>0)'))
rownames(Table_S4) <- rep(c(years,'GM','SD'))
for (c in 1:3) {
  k <- c(1,3,2)[c]
  ic <- c(1,3,5)[c]
  A <- F_1[,1,,k]
  B <- F_1[,2,,k]
  summ <- apply(A-B,2,summarize)
  summ <- cbind(summ,summarize(crossyr(A)[['GM']]-crossyr(B)[['GM']]))
  summ <- cbind(summ,summarize(crossyr(A)[['SD']]-crossyr(B)[['SD']]))
  summ <- round(summ,2)
  for (y in 1:11) {
    Table_S4[y,ic] <- paste0(sprintf("%.2f",summ['mean',y]),' [',sprintf("%.2f",summ['2.5%',y]),',',sprintf("%.2f",summ['97.5%',y]),']')
    Table_S4[y,ic+1] <- sprintf("%.2f",summ['Pr(x>0)',y])
  }
}
print(Table_S4,quote=F)

## Evidence for pairwise year difference in F_1 P(delta>0) (data underlying heatmaps of Fig. S20)
nyrs <- length(yrlabs)
difference <- array(dim=c(N_sex,N_tactic,N_year,N_year))
for (s in 1:2) {
  for (k in 1:3) {
    for (y in 1:9) {
      for (z in 1:9) {
        difference[s,k,y,z] <- mean((F_1[,s,y,k]-F_1[,s,z,k])>0)
      }
    }
  }
}
# for full-R first:'heatmapS20_fullR'
pmatrix_females <- round(difference[1,1,,],2)
pmatrix_males <- round(difference[2,1,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS20_fullR <- pmatrix_females
heatmapS20_fullR[upper.tri(heatmapS20_fullR)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS20_fullR[1:9,] <- heatmapS20_fullR[9:1,] # upper triangle: females, lower triangle: males
print(heatmapS20_fullR)
# then for mixed-RM:'heatmapS20_mixedRM'
pmatrix_females <- round(difference[1,3,,],2)
pmatrix_males <- round(difference[2,3,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS20_mixedRM <- pmatrix_females
heatmapS20_mixedRM[upper.tri(heatmapS20_mixedRM)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS20_mixedRM[1:9,] <- heatmapS20_mixedRM[9:1,] # upper triangle: females, lower triangle: males
print(heatmapS20_mixedRM)
# then for full-M:'heatmapS20_fullM'
pmatrix_females <- round(difference[1,2,,],2)
pmatrix_males <- round(difference[2,2,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS20_fullM <- pmatrix_females
heatmapS20_fullM[upper.tri(heatmapS20_fullM)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS20_fullM[1:9,] <- heatmapS20_fullM[9:1,] # upper triangle: females, lower triangle: males
print(heatmapS20_fullM)

### Additional summaries of between-year tactic switching probabilities
### Results subsection "Supraflexibility between years" and OSM S5.3

## Basic cross-year summaries per sex (provided in main text)
switching_summary <- array(dim=c(N_tactic,N_tactic,N_sex))
dimnames(switching_summary) <- list(paste('from',tactics),paste('to',tactics),sexes)
for (s in 1:2) {
  cnt1 <- 0
  for (k in c(1,3,2)) {
    cnt1 <- cnt1 + 1
    cnt2 <- 0
    for (l in c(1,3,2)) {
      cnt2 <- cnt2 + 1
      summ <- crossyrsumm(kappa[,s,k,2:9,l])
      switching_summary[cnt1,cnt2,s] <- paste0(round(summ[1,1],2),' [',round(summ[1,2],2),',',round(summ[1,6],2),']')
    }
  }
}
print(switching_summary)

## Detailed summary of between-tactic differences in the probability of repeating the same tactic (OSM S5.3 Table S5)
Table_S5 <- array(dim=c((N_year-1+2)*2,N_tactic*2+1))
colnames(Table_S5) <- paste(c('sex','FRvsMRM-delta','FRvsMRM-delta-P','FRvsFM-delta','FRvsFM-P','MRMvsFM-delta','MRMvsFM-P'))
rownames(Table_S5) <- rep(c(years[-1],'GM','SD'),each=2)
comp <- cbind(c(1,3),c(1,2),c(3,2))
for (c in 1:3) {
  ic <- c(1,3,5)[c]
  for (s in 1:2) {
    A <- kappa[,s,comp[1,c],2:9,comp[1,c]]
    B <- kappa[,s,comp[2,c],2:9,comp[2,c]]
    summ <- apply(A-B,2,summarize)
    summ <- cbind(summ,summarize(crossyr(A)[['GM']]-crossyr(B)[['GM']]))
    summ <- cbind(summ,summarize(crossyr(A)[['SD']]-crossyr(B)[['SD']]))
    summ <- round(summ,2)
    for (y in 1:10) {
      ir <- seq(1,20,2)[y]-1+s
      Table_S5[ir,1] <- sex[s]
      Table_S5[ir,ic+1] <- paste0(sprintf("%.2f",summ['mean',y]),' [',sprintf("%.2f",summ['2.5%',y]),',',sprintf("%.2f",summ['97.5%',y]),']')
      Table_S5[ir,ic+2] <- sprintf("%.2f",summ['Pr(x>0)',y])
    }
  }
}
print(Table_S5,quote=F)

## Detailed summary of within-tactic differences in the probability of switching to an alternative tactic (OSM S5.3 Table S6)
Table_S6 <- array(dim=c((N_year-1+2)*2,N_tactic*2+1))
sex <- c('F','M')
colnames(Table_S6) <- paste(c('sex',rep(tactics,each=2)))
rownames(Table_S6) <- rep(c(years[-1],'GM','SD'),each=2)
comp <- cbind(c(3,2),c(1,3),c(1,2))
for (c in 1:3) {
  ic <- c(1,5,3)[c]
  for (s in 1:2) {
    A <- kappa[,s,c,2:9,comp[1,c]]
    B <- kappa[,s,c,2:9,comp[2,c]]
    summ <- apply(A-B,2,summarize)
    summ <- cbind(summ,summarize(crossyr(A)[['GM']]-crossyr(B)[['GM']]))
    summ <- cbind(summ,summarize(crossyr(A)[['SD']]-crossyr(B)[['SD']]))
    summ <- round(summ,2)
    for (y in 1:10) {
      ir <- seq(1,20,2)[y]-1+s
      Table_S6[ir,1] <- sex[s]
      Table_S6[ir,ic+1] <- paste0(sprintf("%.2f",summ['mean',y]),' [',sprintf("%.2f",summ['2.5%',y]),',',sprintf("%.2f",summ['97.5%',y]),']')
      Table_S6[ir,ic+2] <- sprintf("%.2f",summ['Pr(x>0)',y])
    }
  }
}
print(Table_S6,quote=F)

## Detailed summary of between-sex differences in the probability of repeating the same tactic (OSM S5.3 Table S7)
Table_S7 <- array(dim=c((N_year-1+2)*3,N_tactic*2))
colnames(Table_S7) <- paste(c('FR-delta','FR-delta-P','MRM-delta','MRM-P','FM-delta','FM-P'))
rownames(Table_S7) <- rep(c(years[-1],'GM','SD'),3)
for (c in 1:3) {
  cnt <- 0
  for (c2 in 1:3) {
    k <- c(1,3,2)[c]
    k2 <- c(1,3,2)[c2]
    ic <- c(1,3,5)[c]
    A <- kappa[,1,k,2:9,k2]
    B <- kappa[,2,k,2:9,k2]
    summ <- apply(A-B,2,summarize)
    summ <- cbind(summ,summarize(crossyr(A)[['GM']]-crossyr(B)[['GM']]))
    summ <- cbind(summ,summarize(crossyr(A)[['SD']]-crossyr(B)[['SD']]))
    summ <- round(summ,2)
    for (y in 1:10) {
      cnt <- cnt+1
      Table_S7[cnt,ic] <- paste0(sprintf("%.2f",summ['mean',y]),' [',sprintf("%.2f",summ['2.5%',y]),',',sprintf("%.2f",summ['97.5%',y]),']')
      Table_S7[cnt,ic+1] <- sprintf("%.2f",summ['Pr(x>0)',y])
    }
  }
}
print(Table_S7,quote=F)

## Evidence for pairwise year difference in probabilities of repeating the same tactic (data underlying heatmaps of Fig. S21)
yrlabs <- years[-9]
nyrs <- length(yrlabs)
difference <- array(dim=c(N_sex,N_tactic,N_year-1,N_year-1))
for (s in 1:2) {
  for (k in 1:3) {
    for (y in 2:9) {
      for (z in 2:9) {
        difference[s,k,y-1,z-1] <- mean((kappa[,s,k,y,k]-kappa[,s,k,z,k])>0)
      }
    }
  }
}
# for full-R first:'heatmapS21_fullR'
pmatrix_females <- round(difference[1,1,,],2)
pmatrix_males <- round(difference[2,1,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS21_fullR <- pmatrix_females
heatmapS21_fullR[upper.tri(heatmapS21_fullR)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS21_fullR[1:8,] <- heatmapS21_fullR[8:1,] # upper triangle: females, lower triangle: males
print(heatmapS21_fullR)
# then for mixed-RM:'heatmapS21_mixedRM'
pmatrix_females <- round(difference[1,3,,],2)
pmatrix_males <- round(difference[2,3,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS21_mixedRM <- pmatrix_females
heatmapS21_mixedRM[upper.tri(heatmapS21_mixedRM)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS21_mixedRM[1:8,] <- heatmapS21_mixedRM[8:1,] # upper triangle: females, lower triangle: males
print(heatmapS21_mixedRM)
# then for full-M:'heatmapS21_fullM'
pmatrix_females <- round(difference[1,2,,],2)
pmatrix_males <- round(difference[2,2,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS21_fullM <- pmatrix_females
heatmapS21_fullM[upper.tri(heatmapS21_fullM)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS21_fullM[1:8,] <- heatmapS21_fullM[8:1,] # upper triangle: females, lower triangle: males
print(heatmapS21_fullM)

## Evidence for pairwise year difference in probabilities of switching to an alternative tactic (data underlying heatmaps of Fig. S22)
rho <- rhosigned <- array(dim=c(N_iter,2,3,10))
for (s in 1:2) {
  for (k in 1:3) {
    to <- ifelse(k==1,2,ifelse(k==2,3,2))
    p <- kappa[,s,k,,to]/(1-kappa[,s,k,,k]) # probability of switching to the rightest tactic on the R-M continuum
    rhosigned[,s,k,] <- p-0.5 # negative when preference to the left, positive when preference to the right
    rho[,s,k,] <- 2*abs(rhosigned[,s,k,])
  }
}
yrlabs <- years[-1]
nyrs <- length(yrlabs)
difference <- array(dim=c(N_sex,N_tactic,N_year-1,N_year-1))
for (s in 1:2) {
  for (k in 1:3) {
    for (y in 2:9) {
      for (z in 2:9) {
        difference[s,k,y-1,z-1] <- mean((rhosigned[,s,k,y]*rhosigned[,s,k,z])<0)
      }
    }
  }
}
# for full-R first:'heatmapS22_fullR'
pmatrix_females <- round(difference[1,1,,],2)
pmatrix_males <- round(difference[2,1,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS22_fullR <- pmatrix_females
heatmapS22_fullR[upper.tri(heatmapS22_fullR)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS22_fullR[1:8,] <- heatmapS22_fullR[8:1,] # upper triangle: females, lower triangle: males
print(heatmapS22_fullR)
# then for mixed-RM:'heatmapS22_mixedRM'
pmatrix_females <- round(difference[1,3,,],2)
pmatrix_males <- round(difference[2,3,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS22_mixedRM <- pmatrix_females
heatmapS22_mixedRM[upper.tri(heatmapS22_mixedRM)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS22_mixedRM[1:8,] <- heatmapS22_mixedRM[8:1,] # upper triangle: females, lower triangle: males
print(heatmapS22_mixedRM)
# then for full-M:'heatmapS22_fullM'
pmatrix_females <- round(difference[1,2,,],2)
pmatrix_males <- round(difference[2,2,,],2)
colnames(pmatrix_females) <- rownames(pmatrix_females) <- colnames(pmatrix_males) <- rownames(pmatrix_males) <- as.character(yrlabs)
pmatrix_females[upper.tri(pmatrix_females)] <- pmatrix_males[lower.tri(pmatrix_males)] <- NA
for (i in 1:nyrs) pmatrix_females[i,i] <- pmatrix_males[i,i] <- NA
heatmapS22_fullM <- pmatrix_females
heatmapS22_fullM[upper.tri(heatmapS22_fullM)] <- pmatrix_males[upper.tri(pmatrix_males)]
heatmapS22_fullM[1:8,] <- heatmapS22_fullM[8:1,] # upper triangle: females, lower triangle: males
print(heatmapS22_fullM)


### Additional summaries on survival selection on tactics
### Results subsection "Survival selection on tactics and flexibility"

## Between-tactic differences in survival
comp <- cbind(c(1,3),c(1,2),c(3,2))
survival_differences <- array(dim=c(N_tactic*2,7,4))
dimnames(survival_differences) <-   list(paste(rep(c('fullR-mixedRM','fullR-fullM','mixedRM-fullM'),each=2),sexes,sep="::"),
                                         c('mean','2.5%','25%','50%','75%','97.5%','Pr(x>0)'),
                                         c('Non-extreme years','2012-13','2013-14','2017-18'))
for (y in 1:4) {
  cnt <- 0
  for (c in 1:3) {
    for (s in 1:2) {
      cnt <- cnt+1
      survival_differences[cnt,,y] <- round(summarize(yphi[,s,comp[1,c],y]-yphi[,s,comp[2,c],y]),2)
    }
  }
}
print(survival_differences)

## Posterior probabilities of selection shapes
shape_proba <- array(dim=c(4,4,N_sex))
dimnames(shape_proba) <- list(c('Stabilising','Disruptive','Directional towards full-R','Directional towards full-M'),
                              c('Non-extreme years','2012-13','2013-14','2017-18'),
                              sexes)
for (s in 1:2) {
  for (y in 1:4) {
    shape_proba[1,y,s] <- round(mean(yphi[,s,1,y]<yphi[,s,3,y] & yphi[,s,3,y]>yphi[,s,2,y]),2)
    shape_proba[2,y,s] <- round(mean(yphi[,s,1,y]>yphi[,s,3,y] & yphi[,s,3,y]<yphi[,s,2,y]),2)
    shape_proba[3,y,s] <- round(mean(yphi[,s,1,y]>yphi[,s,3,y] & yphi[,s,3,y]>yphi[,s,2,y]),2)
    shape_proba[4,y,s] <- round(mean(yphi[,s,1,y]<yphi[,s,3,y] & yphi[,s,3,y]<yphi[,s,2,y]),2)
  }
}
print(shape_proba)


### Additional summaries on phenotypic dynamics
### Results subsection "Net effects on phenotypic dynamics" and OSM S5.4

## Detailed cross-year summary of sequential effects on year-to-year variation E_1, E_2 and E_3 (OSM S5.4 Table S8)

Table_S8 <- array(dim=c(N_year,5,N_sex))

summ <- array(dim=c(3,3,3,N_tactic,N_sex))
dimnames(summ) <- list(c('mean','2.5%','97.5%'),c('GM','SD','GM(abs)'),c('E_1','E_2','E_3'),tactics,sexes)

for (s in 1:2) {
  cnt <- 0
  for (k in c(1,3,2)) {
    cnt <- cnt+1
    summ[,1:2,1,cnt,s] <- t(crossyrsumm(E_1[,s,,k])[,c('mean','2.5%','97.5%')])
    summ[,1:2,2,cnt,s] <- t(crossyrsumm(E_2[,s,,k])[,c('mean','2.5%','97.5%')])
    summ[,1:2,3,cnt,s] <- t(crossyrsumm(E_3[,s,,k])[,c('mean','2.5%','97.5%')])
    summ[,3,1,cnt,s] <- t(crossyrsumm(abs(E_1[,s,,k]))[1,c('mean','2.5%','97.5%')])
    summ[,3,2,cnt,s] <- t(crossyrsumm(abs(E_2[,s,,k]))[1,c('mean','2.5%','97.5%')])
    summ[,3,3,cnt,s] <- t(crossyrsumm(abs(E_3[,s,,k]))[1,c('mean','2.5%','97.5%')])
    for (e in 1:3) {
      Table_S8[e,1:2,s] <- c('GM',paste0('E',e))
      Table_S8[e,2+cnt,s] <- paste0(sprintf("%.2f",summ['mean','GM',e,cnt,s]),' [',sprintf("%.2f",summ['2.5%','GM',e,cnt,s]),',',sprintf("%.2f",summ['97.5%','GM',e,cnt,s]),']')
      Table_S8[3+e,1:2,s] <- c('SD',paste0('E',e))
      Table_S8[3+e,2+cnt,s] <- paste0(sprintf("%.2f",summ['mean','SD',e,cnt,s]),' [',sprintf("%.2f",summ['2.5%','SD',e,cnt,s]),',',sprintf("%.2f",summ['97.5%','SD',e,cnt,s]),']')
      Table_S8[6+e,1:2,s] <- c('GM(abs)',paste0('|E',e,'|'))
      Table_S8[6+e,2+cnt,s] <- paste0(sprintf("%.2f",summ['mean','GM(abs)',e,cnt,s]),' [',sprintf("%.2f",summ['2.5%','GM(abs)',e,cnt,s]),',',sprintf("%.2f",summ['97.5%','GM(abs)',e,cnt,s]),']')
    }
  }
}
dimnames(Table_S8) <- list(NULL,c('Quantity','Effect',tactics),sexes)
print(Table_S8,quote=F)

## Detailed summary of sequential effects on year-to-year variation (E_1, E2, E_3) in females (OSM S5.4 Table S9)
Table_S9 <- array(dim=c(25,8))
summ <- array(dim=c(4,3,N_tactic,N_year))
dimnames(summ) <- list(c('mean','2.5%','97.5%','Pr(x>0)'),c('E_1','E_2','E_3'),tactics,fullyrlabs)
s <- 1 # females
cntrow <- 0
for (y in 1:9) {
  for (e in 1:3) {
    if (!((y==9 & e==3)|(y==9 & e==2))) {
      cntrow <- cntrow+1
      cnt <- 0
      for (k in c(1,3,2)) {
        cnt <- cnt+1
        summ[,e,cnt,y] <- summarize(get(paste0('E_',e))[,s,y,k])[c('mean','2.5%','97.5%','Pr(x>0)')]
        Table_S9[cntrow,(3:4)+ifelse(cnt==1,0,ifelse(cnt==2,2,4))] <- c(paste0(sprintf("%.2f",summ['mean',e,cnt,y]),' [',sprintf("%.2f",summ['2.5%',e,cnt,y]),',',sprintf("%.2f",summ['97.5%',e,cnt,y]),']'),sprintf("%.2f",summ['Pr(x>0)',e,cnt,y]))
      }
      Table_S9[cntrow,1:2] <- c(fullyrlabs[y],paste0('E',e))

    }
  }
}
colnames(Table_S9) <- c('Year','Effect',paste0(rep(tactics,each=2),':',c('Value','Pr(E>0)')))
print(Table_S9,quote=F)

## Detailed summary of sequential effects on year-to-year variation (E_1, E2, E_3) in males (OSM S5.4 Table S10)
Table_S10 <- array(dim=c(25,8))
summ <- array(dim=c(4,3,N_tactic,N_year))
dimnames(summ) <- list(c('mean','2.5%','97.5%','Pr(x>0)'),c('E_1','E_2','E_3'),tactics,fullyrlabs)
s <- 2 # males
cntrow <- 0
for (y in 1:9) {
  for (e in 1:3) {
    if (!((y==9 & e==3)|(y==9 & e==2))) {
      cntrow <- cntrow+1
      cnt <- 0
      for (k in c(1,3,2)) {
        cnt <- cnt+1
        summ[,e,cnt,y] <- summarize(get(paste0('E_',e))[,s,y,k])[c('mean','2.5%','97.5%','Pr(x>0)')]
        Table_S10[cntrow,(3:4)+ifelse(cnt==1,0,ifelse(cnt==2,2,4))] <- c(paste0(sprintf("%.2f",summ['mean',e,cnt,y]),' [',sprintf("%.2f",summ['2.5%',e,cnt,y]),',',sprintf("%.2f",summ['97.5%',e,cnt,y]),']'),sprintf("%.2f",summ['Pr(x>0)',e,cnt,y]))
      }
      Table_S10[cntrow,1:2] <- c(fullyrlabs[y],paste0('E',e))
      
    }
  }
}
colnames(Table_S10) <- c('Year','Effect',paste0(rep(tactics,each=2),':',c('Value','Pr(E>0)')))
print(Table_S10,quote=F)

## Detailed cross-year summary of sequential effects on year-to-year variation (OSM S5.4 Table S11)
Table_S11 <- array(dim=c(2*3,N_tactic*2,N_sex))
dimnames(Table_S11) <- list(c('GM-1vs2','GM-1vs3','GM-2vs3','SD-1vs2','SD-1vs3','SD-2vs3'),
                            paste(c('FullR-delta','FullR-P(delta>0)','MixedRM-delta','MixedRM-P(delta>0)','FullM-delta','FullM-P(delta>0)')),
                            sexes)
for (s in 1:2) {
  for (k in c(1,3,2)) {
    lk <- c(1,5,3)[k]
    ABdif <- crossyr(abs(E_1[,s,-9,k]))[['GM']] - crossyr(abs(E_2[,s,,k]))[['GM']]
    ACdif <- crossyr(abs(E_1[,s,-9,k]))[['GM']] - crossyr(abs(E_3[,s,,k]))[['GM']]
    BCdif <- crossyr(abs(E_2[,s,,k]))[['GM']] - crossyr(abs(E_3[,s,,k]))[['GM']]
    summ <- summarize(ABdif)
    summ <- cbind(summ,summarize(ACdif))
    summ <- cbind(summ,summarize(BCdif))
    ABdif <- crossyr(abs(E_1[,s,-9,k]))[['SD']] - crossyr(abs(E_2[,s,,k]))[['SD']]
    ACdif <- crossyr(abs(E_1[,s,-9,k]))[['SD']] - crossyr(abs(E_3[,s,,k]))[['SD']]
    BCdif <- crossyr(abs(E_2[,s,,k]))[['SD']] - crossyr(abs(E_3[,s,,k]))[['SD']]
    summ <- cbind(summ,summarize(ABdif))
    summ <- cbind(summ,summarize(ACdif))
    summ <- cbind(summ,summarize(BCdif))
    for (c in 1:6) {
      Table_S11[c,lk,s] <- paste0(sprintf("%.2f",summ['mean',c]),' [',sprintf("%.2f",summ['2.5%',c]),',',sprintf("%.2f",summ['97.5%',c]),']')
      Table_S11[c,lk+1,s] <- sprintf("%.2f",summ['Pr(x>0)',c])
    }
  }
}
print(Table_S11)

