rm( list = ls() )
workdir <- "/Users/awwood/Desktop/r_dir/waynyi/pcangsd"
setwd(workdir)

colors <- scan( file = "population_color_list.txt", what = "character" )

matrix <- as.matrix( read.table ( file="whole_genome.cov", stringsAsFactors = F ) )
e <- eigen( matrix )

eigen_sum <- sum(e$values)
pc3 <- e$values[3] / eigen_sum
pc3 <- round(pc3, digits=3)
pc3 <- sprintf("%.1f%%", 100*pc3)
pc4 <- e$values[4] / eigen_sum
pc4 <- round(pc4, digits=3)
pc4 <- sprintf("%.1f%%", 100*pc4)
pc2 <- e$values[2] / eigen_sum
pc2 <- round(pc2, digits=3)
pc2 <- sprintf("%.1f%%", 100*pc2)

# Calculate the percent variance explained by each principal component.

all_pc <- c(rep(0, length(e$values)) )

for (i in 1 : length(all_pc) ){

pc <- e$values[i] / eigen_sum
pc <- round(pc, digits=3)
#pc <- 100*pc

all_pc[i] <- pc

}

names(all_pc) <- as.character(1:29)
barplot(all_pc, col = "Black", names = labels, xlab)

x <- paste( "PC2 (", pc2,")", sep = "")
y <- paste( "PC3 (", pc3,")", sep = "")

lgnd <- c("TN", "AR", "Uwharrie", "IN","NY", "NC [S. v. waynei]" )
lgnd_fill <- c( "dodgerblue", "cyan", "darkorchid", "cadetblue", "blue", "gold" )

pdf( file = "pc2_3_population_color_Whole_Genome.pdf" )
par(mar = c(6,6,2,2))
plot(e$vectors[,2:3],lwd=2,ylab=y,xlab=x,main="Setophaga virens waynei whole genome", bg=colors, pch=21, col="black", cex=3, cex.lab=1.0, cex.main=1.0 )
legend("bottomright", legend = lgnd , fill = lgnd_fill )
dev.off()
	
