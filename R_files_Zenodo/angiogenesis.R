#Figure 1C
y <- readRDS("smallBRCADataset.rds")
mat <- y$data
labels <- y$pheno$SeuratCluster

#expression normalization
mat <- mat[,which(labels == "Basal")]
mat <- t(t(mat)/apply(mat, 2, sum))*1000000
mat <- log(mat+1,2)
lmo2 <- mat["LMO2",]

genelists <- fread("/data1/ggulati/mberger/validationSets/FullTables_validationVersion2.0/GenesetAnalysis/msigdb_Tirosh_Table.txt", header = F)
genelists <- genelists[,-1]

library(GSVA)
#Hallmark angiogenesis
hallangio <- genelists[grep("HALLMARK_ANGIOGENESIS", genelists$V2),]
hallangio <- as.character(unlist(hallangio$V3))
hallangio <- gsva(data.matrix(mat), list(hallangio), method = "ssgsea", parallel.sz = 12)

t.test(hallangio[which(lmo2>0)], hallangio[-which(lmo2>0)], var.equal = T)


p = 0
for(i in 1:10000){
set.seed(i)
temp <- hallmark[,sample(which(lmo2 == 0), 7)]
temp2 <- apply(temp, 1, mean)
temp3 <- apply(hallmark[,which(lmo2>0)], 1, mean)
temp4 <- as.numeric(temp3>temp2)
ifelse(i ==1, temp5 <- temp4, temp5 <- temp5 + temp4) 
}
names(temp5) <- rownames(hallmark)
pval <- 1 - temp5/10000

hallmark["HALLMARK_ANGIOGENESIS",] -> angio

sup <- factor(lmo2>0)
levels(sup) <- c("LMO2neg", "LMO2pos")
sup <- relevel(sup, ref = "LMO2neg")


##
tog = -0.1
tog2 = 0.4
##BOXPLOT
pdf("LMO2_angiogenesis.pdf", width = 4, height = 7, useDingbats = FALSE)
par(oma = c(2,3,4,2))
par(mar = c(5,6.75,1,1), xpd = NA)
cols2 <- rep(c(adjustcolor("darkred", alpha.f = 0.6), adjustcolor("navyblue", alpha.f = 0.6)),  3)
cols3 <- rep(c("darkred", "navyblue"), 3)
boxplot(angio~sup, outline = F, xaxt = "n", yaxt = "n",
        staplelwd = 1,
        medlwd = 2.5,
        whisklty = 1,
        border = "black",
        col = adjustcolor(cols2, alpha.f = 0.6),
        ylab = "",
        xlab = "",
        cex.lab = 1.75,
        las = 1,
        frame.plot = F,
        ylim = c(tog,tog2),
        cex.axis = 2)
mtext("ssGSEA",
      cex = 2, side = 2, line = 5.5)
title(main = "",font.main = 1, cex.main = 2.25, line = 3)
axis(1, pos = tog, at = 1:length(unique(sup)), labels = FALSE)
axis(2, pos = 0, at = signif(as.numeric(formatC(seq(tog, tog2, 0.1), format = "f"), 2)), las = 1, cex.axis = 2)
segments(-0.15,tog,length(unique(sup))+0.5,tog)
text(seq_along(1:length(unique(sup))), par("usr")[3] - 0.0125, labels = levels(sup),
     srt = 45, adj = 1, xpd = TRUE, col = cols, cex = 2)
#set.seed(333)
#stripchart(angio~sup, vertical = TRUE,
#           method = "jitter", jitter = 0.3, add = TRUE, pch = 16,
#           col = adjustcolor(cols3,alpha.f = 0.5),lwd = 1, cex = 0.85)
#set.seed(333)
#stripchart(angio~sup, vertical = TRUE,
#           method = "jitter", jitter = 0.3, add = TRUE, pch = 1,
#           col = cols3,lwd = 1, cex = 0.85)
dev.off()

