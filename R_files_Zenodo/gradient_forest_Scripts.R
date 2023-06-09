require(raster)
require(rgdal)

snp <- read.table("snp.forR", header = T, row.names = 1)

sample.coord <-read.table("cs.sample.location.txt", header=T, stringsAsFactors=F)

library(vegan)
coord <- sample.coord[,c("longitude","latitude")]
pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
write.table(pcnm.keep, "pcnm.keep", sep="\\t", quote=F, row.names=F) 
pcnm.keep

clim.points <- read.table("cs.6climate", header = T, stringsAsFactors=F)

env.gf <- cbind(clim.points, pcnm.keep)

library(gradientForest)

maxLevel <- log2(0.368*nrow(env.gf)/2)
gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=1000, maxLevel=maxLevel, trace=T, corr.threshold=0.50)

save.image("/ipm1/ljcao/cs.rad/2_analysis_and_processing_of_data/5_outlier/gradientforest/reduced/bios6/cs.gf.data")

###plot bar graphs depicting the importance of each spatial and climate variable.
pdf("GF_VariableImportance.pdf")
plot(gf, plot.type = "O")
dev.off()

###plot the "turnover functions" showing how allelic composition changes along the spatial or environmental gradients. 
by.importance <- names(importance(gf))

pdf("GF_TurnoverFunctions.pdf")
plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
dev.off()

cu.sp.bio02 <- cumimp(gf, "bio01", "Overall")
cu.sp.bio03 <- cumimp(gf, "bio02", "Overall")
cu.sp.bio04 <- cumimp(gf, "bio03", "Overall")
cu.sp.bio05 <- cumimp(gf, "bio10", "Overall")
cu.sp.bio15 <- cumimp(gf, "bio13", "Overall")
cu.sp.bio18 <- cumimp(gf, "bio15", "Overall")
cu.sp.PCNM1 <- cumimp(gf, "PCNM1", "Overall")
cu.sp.PCNM2 <- cumimp(gf, "PCNM2", "Overall")
cu.sp.PCNM3 <- cumimp(gf, "PCNM3", "Overall")
cu.sp.PCNM4 <- cumimp(gf, "PCNM4", "Overall")

write.table(cu.sp.bio02,file="cs.bio01.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio03,file="cs.bio02.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio04,file="cs.bio03.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio05,file="cs.bio10.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio15,file="cs.bio13.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio18,file="cs.bio15.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.PCNM1,file="cs.PCNM1.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.PCNM2,file="cs.PCNM2.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.PCNM3,file="cs.PCNM3.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.PCNM4,file="cs.PCNM4.cumimp.PD.res",sep = " ",quote = F)

gm.imp <- as.data.frame(t(gf$imp.rsq))

write.table(gm.imp,"cs.gf.snp.imp",sep = " ",row.names=T,col.names = T)

for (i in 1:length(gm.imp[,"bio02"])) {
  bar <- gm.imp[i,]
  gm.imp[i,11] <- names(which.max(abs(bar[1:10]))) # gives the variable, j in gm.imp[i,j] is 3+biovars+PCNMvars+1, m and n in bar[m:n] is m=4, n= 4+vasrs(biovars + PCNMvars)
  gm.imp[i,12] <- max(abs(bar[1:10]))              # gives the correlation,j in gm.imp[i,j] is 3+biovars+PCNMvars+2, m and n in bar[m:n] is m=4, n= 4+vasrs(biovars + PCNMvars)
}

colnames(gm.imp)[11] <- "predictor"
colnames(gm.imp)[12] <- "correlation"

table(gm.imp$predictor)
gm.imp.bios <- gm.imp[which(gm.imp$predictor%in%c("bio01","bio02","bio03","bio10","bio13","bio15")),]
#write.table(gm.imp.bios,"gm.gf.snp.imp.7bios",sep = " ",row.names=T,col.names = T)

quants.cor <- quantile(gm.imp.bios[,"correlation"],probs = c(0.5,0.9,0.95,0.99,0.999),names = T)
quants.cor[1]
gm.imp.bios[gm.imp.bios[,"correlation"]<quants.cor[1],] <- NA

gm.imp.bios.high <- gm.imp.bios[apply(gm.imp.bios[,1:12],1,function(gm.imp.bios)!all(is.na(gm.imp.bios))),]
write.table(gm.imp.bios.high,"cs.gf.imp.snp.bios.cor.0.5.bios",quote = F,sep = " ",row.names = T,col.names = T)

