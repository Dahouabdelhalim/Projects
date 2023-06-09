##data preparation
#plink --vcf XXX.vcf --recode A --allow-extra-chr --out gm.326ind ###for stacks v1.4

#plink --vcf ps.all.203.obs0.8.maf0.05.19634.recode.vcf --make-bed --allow-extra-chr --out ps.203ind  ###for stacks 2.0
#plink --bfile ps.203ind --recode A --allow-extra-chr --out ps.203ind   ###for stacks 2.0
#install.packages(c("psych","vegan"), dependencies=TRUE)
#install.packages("adegenet")
library(psych)    
library(vegan)    

gm.snp<- read.table("ps.203.maf0.05.19634.raw",header = T,row.names = 2)[,-(1:5)] #read into R
dim(gm.snp)    #check
sum(is.na(gm.snp)) #NAs in the matrix 
write.table(gm.snp,"ps.snp.maf0.05.19634.txt",quote = F)
###replace NAs with the most common genotype at each SNP across all individuals
gm.snp.imp <- apply(gm.snp,2,function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gm.snp.imp))

gm.env <- read.table("ps.ind.env", header = T)
str(gm.env)
gm.env$IID <- as.character(gm.env$IID) # Make individual names characters (not factors)


coord <-read.table("sample.coord.txt", header=T, stringsAsFactors=F)
coord

pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
str(pcnm.keep)
write.table(pcnm.keep, "pcnm.keep", sep="\\t", quote=F, row.names=F) 

gm.env <- cbind(gm.env, pcnm.keep)
str(gm.env)
# Confirm that genotypes and environmental data are in the same order, it must be TURE
identical(rownames(gm.snp.imp), gm.env[,2]) 

pdf("pairs.panels.23.pdf", width = 16, height = 16)
pairs.panels(gm.env[,3:23],scale = T)
dev.off()

#remove strong correlated varables
pred <- gm.env[,3:ncol(gm.env)]
str(pred)
pred <- subset(pred, select = - c(bio11,bio6,bio13, bio16, bio18,bio8,bio2,bio7,bio19,bio1,bio12,bio17,bio4,PCNM1,bio10))
str(pred)

pdf("pairs.panels.6vars.pdf", width = 16, height = 16)
pairs.panels(pred,scale = T)
dev.off()


##full
gm.rda <- rda(gm.snp.imp ~., dat =pred,scale scale = T)
##partial
gm.rda <- rda(gm.snp.imp ~ bio3+bio5+bio9+bio15+bio14+Condition(PCNM2), data = pred, scale = T)
gm.rda  ###show results
###The proportion of the variance explained by the environmental predictors is 
###given under the “Proportion” column for “Constrained”; this is equivalent to 
###the R2 of a multiple regression. Just like in multiple regression, 
###this R2 will be biased and should be adjusted based on the number of predictors. 
###We can calculate the adjusted R2 using:

RsquareAdj(gm.rda)
####Our constrained ordination explains about 13.5% of the variation (OFM data)
###The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:
summary(eigenvals(gm.rda, model = "constrained"))
###We can visualize this information using a screeplot of the canonical eigenvalues by calling screeplot:
pdf("screeplot of the canonical eigenvalues_6vars.pdf",width = 16,height = 10)
screeplot(gm.rda)
dev.off()

###Here, we can see that the first two (OFM data) or three constrained axes explain most of the variance.
###The screeplot provides an informal (and quick) way to determine how many constrained axes to include 
###when we search for candidate SNPs (below). We could start by investigating RDA axes that explain 
###the most variance (excluding those after the “drop off” point in the screeplot.)

###Now let’s check our RDA model for significance using formal tests. 
###We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). 
###The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. 

signif<- anova.cca(gm.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif

###The full model is significant, but that doesn’t tell us much. We can check each constrained axis for significance 
###using the code below. For this test, each constrained axis is tested using all previous constrained axes as conditions. 
###The purpose here is to determine which constrained axes we should investigate for candidate loci.

signif.axis <- anova.cca(gm.rda, by="axis", parallel=getOption("mc.cores"))  ##taking up to a few hours
signif.axis

###vegan has a simple function for checking Variance Inflation Factors for the predictor variables used in the model:
vif.cca(gm.rda)

###0< VIF <=5 no multicollinearity, 5 < VIF <= 10 slight multicollinearity, 10 < VIF <= 100 moderate multicollinearity, VIF > 100 strong multicollinearity
####All values are below 10, and most are below 5, which indicates that multicollinearity among these predictors shouldn’t be a problem for the model. 
###We could remove one of the XXX variables (bio? or bio?) if we were concerned about these higher VIF values
pdf("simple.plots.gm.rda.6vars.pdf", width = 16, height = 16)
plot(gm.rda, scaling = 3)
plot(gm.rda, choices = c(1,3), scaling = 3)
dev.off()

###Here, the SNPs are in red (in the center of each plot), and the individuals are the black circles. 
###The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space 
###reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.

###color the individual points based on their pop, which we can find in the env data set.
levels(gm.env$pop) <- c("AHHF", "FJPT", "GDGZ", "GXNN", "HNCS", "HNYY", "JSWX", "JXGZ", "JXNC", "XJTL", "ZJHZ")
pop <- gm.env$pop
#library("RColorBrewer")
##bg <- brewer.pal(n = 11, name = 'RdBu') ###Creates nice looking color palettes especially for thematic map
bg <- rainbow(11)

pdf("ps.rda with snp pop and 6vars on axes12.pdf", width = 16, height = 16)
plot(gm.rda, type="n", scaling=3)
points(gm.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3) #the SNPs
points(gm.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[pop]) #the pop
text(gm.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

pdf("ps.rda with snp pop and 9vars on axes13.pdf", width = 16, height = 16)
plot(gm.rda, type="n", scaling=3, choices=c(1,3))
points(gm.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(gm.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[pop], choices=c(1,3))
text(gm.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

#####We’ll use the loadings of the SNPs in the ordination space to determine which SNPs are candidates 
###for local adaptation. The SNP loadings are stored as species in the RDA object. We’ll extract the 
###SNP loadings from the three significant constrained axes:
load.rda <- scores(gm.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

pdf("snp loadings on RDA_6vars.pdf",width = 8,height = 4)
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 
dev.off()
###SNPs loading at the center of the distribution are not showing a relationship with the environmental 
###predictors; those loading in the tails are, and are more likely to be under selection as a function 
###of those predictors (or some other predictor correlated with them).

### a function to identify SNPs that load in the tails of these distributions.
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

##apply it to each significant constrained axis:
cand1 <- outliers(load.rda[,1],3) # 32
cand2 <- outliers(load.rda[,2],3) # 5
cand3 <- outliers(load.rda[,3],3) # 75
####standard deviation cutoff can be 2.5, 3, or 3.5. higher sd identify loci under stronger selection.

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
###For OFM data, we have 32 candidates on axis 1, 5 on axis 2, and 75 on axis 3, for a total of 112 candidate SNPs

###organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor:
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

###add in the correlations of each candidate SNP with the eight environmental predictors:
ofm <- matrix(nrow=(ncand), ncol=6)  # 10 columns for 10 predictors
colnames(ofm) <- c("bio3", "bio5", "bio9", "bio14","bio15","PCNM2")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gm.snp.imp[,nam]
  ofm[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,ofm)  
head(cand)


###Now we have a data frame of 141 candidate SNPs and their correlation with our 8 environmental predictors.

###looking for duplicate detections. These are SNPs that are identified as candidates on more than one RDA axis.
length(cand$snp[duplicated(cand$snp)])

ofm <- cbind(cand$axis, duplicated(cand$snp)) 
table(ofm[ofm[,1]==1,2]) # duplicates on axis 1
table(ofm[ofm[,1]==2,2]) # duplicates on axis 2
table(ofm[ofm[,1]==3,2]) # duplicates on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections, if any
head(cand)

###which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,13] <- names(which.max(abs(bar[4:12]))) # gives the variable
  cand[i,14] <- max(abs(bar[4:12]))              # gives the correlation
}

colnames(cand)[13] <- "predictor"
colnames(cand)[14] <- "correlation"

table(cand$predictor) 
write.csv(cand,"snp.candiates.3.1.partial.csv")
###Note that, in some cases, correlations may be strong for multiple variables (depending on collinearity among predictors).
###It may be useful to consider how candidate SNPs are correlated with multiple predictors.

### focus in on the SNPs in the ordination space.
###We can color code the SNPs based on the predictor variable that they are most strongly correlated with. 

sel <- cand$snp
env <- cand$predictor
##display.brewer.pal(n = 3, name = 'Paired')  ##brewer.pal(n = 3, name = 'Paired')###Creates nice looking color

#"bio11", "bio12", "bio13", "bio15", "bio19", "bio2", "bio3", "bio4", "bio5", "bio8"

env[env=="bio3"] <- 'red'
env[env=="bio5"] <- 'green'
env[env=="bio9"] <- 'blue'
env[env=="bio14"] <- '#f1eef6'
env[env=="bio15"] <- '#f1eef6'
env[env=="PCNM2"] <- '#f1eef6'
#env[env=="PCNM4"] <- '#f1eef6'
#env[env=="PCNM5"] <- '#f1eef6'
#env[env=="PCNM6"] <- '#f1eef6'

# color by predictor:
col.pred <- rownames(gm.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  ofm <- match(sel[i],col.pred)
  col.pred[ofm] <- env[i]
}

col.pred[grep("X",col.pred)] <- '#f1eef6' # non-candidate SNPs, grep chr or X, depending on your data!
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c("red", "green", "blue", "#f1eef6", "#f1eef6", "#f1eef6", "#f1eef6", "#f1eef6", "#f1eef6")

# plot on axes 1 & 2
pdf("climate 9-6variables correlated3 SNPs axes12.pdf")
plot(gm.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(gm.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(gm.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(gm.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("bio5",  "bio6",  "bio15"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

# axes 1 & 3
pdf("climate 9-6variables correlated3 SNPs axes13.pdf")
plot(gm.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(gm.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(gm.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(gm.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("bio5",  "bio6",  "bio15"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()


