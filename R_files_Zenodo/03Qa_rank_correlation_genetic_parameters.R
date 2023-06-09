setwd("G:/data_analysis/GradientForest20210517/genetic_diversity_rank_correlation")
library(corrplot)

rawdata <- read.table("QaNS_intra_paraeter.txt6", head=T)
head(rawdata)

data <- subset(rawdata, s12>= 1000 & N12 >= 10 & fst12 != "NaN")#& CDS_len != "NaN"
dim(data)


mean(data$sum_rho1)
sd(data$sum_rho1)/42499^0.5

mean(data$sum_rho2)
sd(data$sum_rho2)/42499^0.5

data$Dxy <- data$dxy12/data$s12
data$pi1 <- data$tP1/data$s1
data$pi2 <- data$tP2/data$s2

head(data)


intradata <- data[c("Dxy","fst12","pi1","pi2","sum_rho1","sum_rho2")]

head(intradata)
corr <- cor(intradata)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(mtcars)


pdf("inter_intrapop_data_correlation_20211104.pdf")
col<- colorRampPalette(c("#0000FF", "white", "#FF0000"))(20)
corrplot(corr, method = "circle", type="lower",col=col, diag=F, p.mat = p.mat, sig.level = 0.01, addCoef.col = "black", title = "")
dev.off()











