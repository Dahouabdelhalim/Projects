library(Hmisc)
library(vegan)
all <- t(Phenig$G)
rownames(all)<-Phenig$id
env<- Phenig$env

env<- subset(env, select = c(ind, bio01, bio03, bio05, bio07, bio12, bio09))
colnames(env)[-1]<- c("Ann MT","Isothermality", "MaxT war mon", "T ann rge", "Ann ppt", "MT dry qua")

identical(rownames(all), env$ind) 
#The data are imputed
gen_all.imp <- apply(all, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
colnames(gen_all.imp)<- Phenig$SNP

all.rda <- rda(gen_all.imp~ ., data=sam_env_ordered[,-1], scale=T)
all.rda

signif.all <- anova.cca(all.rda, parallel=getOption("mc.cores"), permutations = 99) 
signif.all

signif_all.axis <- anova.cca(all.rda, by="axis", parallel=getOption("mc.cores"),permutations = 99 )
signif_all.axis
# The five first axes are significant

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

load_all.rda <- scores(all.rda, choices=c(1:5), display="species")

cand1 <- outliers(load_all.rda[,1],3) # 38
cand2 <- outliers(load_all.rda[,2],3) # 69
cand3 <- outliers(load_all.rda[,3],3) # 34
cand4 <- outliers(load_all.rda[,4],3)
cand5 <- outliers(load_all.rda[,5],3)

ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5)  
ncand
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))

colnames(cand1)<-colnames(cand2)<- colnames(cand3) <- colnames(cand4) <- colnames(cand5) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3, cand4, cand5)
cand$snp <- as.character(cand$snp)
snp_all<- cand$snp
snp_all<- unique(cand$snp)
foo <- matrix(nrow=(length(snp_all)), ncol=6)  
colnames(foo) <- colnames(env[,c(2:7)])
for (i in 1:length(snp_all)) {
  nam <- snp_all[i]
  snp.gen <- gen_all.imp[,nam]
  foo[i,] <- apply(env[,c(2:7)],2,function(x) cor(x,snp.gen))
}
cand_all <- data.frame("snp"=snp_all,foo)
head(cand_all)
for (i in 1:length(cand_all$snp)) {
  bar <- cand_all[i,]
  cand_all[i,8] <- names(which.max(abs(bar[2:7]))) # gives the variable
  cand_all[i,9] <- max(abs(bar[2:7]))              # gives the correlation
}
colnames(cand_all)[8] <- "predictor"
colnames(cand_all)[9] <- "correlation"
table(cand_all$predictor)
View(cand_all)

