library(ape)
mims <- read.tree("rooted_tree.tre")
plot(mims, cex = 0.25)
nodelabels(cex = 0.5)
#Use these node numbers to extract clades.
mimallonidae.tree  <- extract.clade(mims, 143)
plot(mimallonidae.tree, cex =0.35)

#PGLS
library(geiger)
library(nlme)
library(phytools)
data<-read.csv("FWL_frenulum_data.csv", row.names = 1)
name.check(mimallonidae.tree, data)
pglsModel <- gls(FWL ~ frenulum, correlation = corBrownian(phy=mimallonidae.tree),data = data, method = "ML")
> pglsModel <- gls(FWL ~ frenulum, correlation = corBrownian(phy=mimallonidae.tree),data = data, method = "ML")
> anova(pglsModel)
Denom. DF: 119 
            numDF   F-value p-value
(Intercept)     1 21.698370  <.0001
frenulum        1  7.983379  0.0055
> coef(pglsModel)
(Intercept)    frenulum 
  16.198502   -3.571894 

#caper
library(caper)
mims_continuous_discrete<-treedata(mimallonidae.tree,data)
mims_continuous_discrete_phy<-mims_continuous_discrete$phy
mims_continuous_discrete_dat<-mims_continuous_discrete$data
mims_continuous_discrete_dat_df<-as.data.frame(mims_continuous_discrete$data)
write.table(mims_continuous_discrete_dat_df, file="mims_continuous_discrete_dat_df.txt", sep = "\\t")
#added "taxon" column header in textwrangler
mims_continuous_discrete_dat_df<-read.table("mims_continuous_discrete_dat_df.txt",header = TRUE, sep = "\\t")
str(mims_continuous_discrete_dat_df)

mims_continuous_discrete_dat_df_as_factors<-mims_continuous_discrete_dat_df
mims_continuous_discrete_dat_df_as_factors$frenulum<-as.factor(mims_continuous_discrete_dat_df_as_factors$frenulum)
str(mims_continuous_discrete_dat_df_as_factors)
caper<-comparative.data(mimallonidae.tree,mims_continuous_discrete_dat_df_as_factors, names.col = taxon)
> caper
Comparative dataset of 121 taxa:
Phylogeny: mimallonidae.tree 
   121 tips, 120 internal nodes
   chr [1:121] "LEP24852_Mimallonidae_Cicinninae_Cicinnini_Gonogramma_submarcata" ...
Data: mims_continuous_discrete_dat_df_as_factors 
   $ frenulum: Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
   $ FWL     : num [1:121] 17.8 22.5 14.2 20 18.2 23.9 22.3 24.3 26 19.4 ...
> comparison<-brunch(frenulum~FWL,data=caper)
> summary(comparison)

Call:
lm(frenulum ~ FWL - 1, data = contrData)

Residuals:
    Min      1Q  Median      3Q     Max 
-0.5215 -0.1299  0.9862  2.8316  5.0702 

Coefficients:
    Estimate Std. Error t value Pr(>|t|)   
FWL -0.18708    0.04672  -4.004  0.00393 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 2.473 on 8 degrees of freedom
Multiple R-squared:  0.6672,	Adjusted R-squared:  0.6256 
F-statistic: 16.04 on 1 and 8 DF,  p-value: 0.003925

###threshBayes
library(coda)
data<-read.csv("FWL_frenulum_data.csv", row.names = 1)
sample <- 1000 # sample every 1000 steps
ngen <- 8000000 # chain length, > 2 million is suggested
burnin <- 0.2 * ngen # 20% of all data is discarded as burnin
head(data)
matrix_data <-as.matrix(data)
head(matrix_data)
thresh<-threshBayes(mimallonidae.tree, matrix_data, types=c("discrete","continuous"),ngen = ngen, control = list(sample = sample))
thresh$par
plot(thresh$par[,"gen"],thresh$par[,"logL"],type="l",
     xlab="generation",ylab="logL")
plot(density(thresh$par[(burnin/sample+1):nrow(thresh$par),"r"],bw=0.1),xlab="r",main="posterior density for r")
abline(v=mean(thresh$par[(burnin/sample + 1):nrow (thresh$par), "r"]),col="red")
r_thresh <- thresh$par[(burnin/sample + 1):nrow(thresh$par), "r"]
class(r_thresh) <- "mcmc"
> effectiveSize(r_thresh)
    var1 
440.5992
> HPDinterval(r_thresh)
         lower      upper
var1 -0.747149 -0.1752151
attr(,"Probability")
[1] 0.9500078
> thresh

Object of class "threshBayes" consisting of a matrix (L) of
sampled liabilities for the tips of the tree & a second matrix
(par) with the sample model parameters & correlation.

Mean correlation (r) from the posterior sample is: -0.49814.

Ordination of discrete traits:

	Trait 1: 0 <-> 1
