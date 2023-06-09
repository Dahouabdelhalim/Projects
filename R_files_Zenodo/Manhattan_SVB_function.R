###
# Functions to adjust scaffold positions, written by Steven Van Belleghem
###

#calculate chromosome coordinates
chrom.coords <- function(scafL,chromNames, gap = 2000000) {
  chromosome = vector()
  chromLengths  = vector()
  chromStarts = vector()
  chromEnds = vector()
  chromMid = vector()
  chrom = 1
  endLast = 0
  scafCurrent <- subset(scafL, chromosome == chromNames[1])
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom]/2
  endLast = chromEnds[chrom]
  chrom = chrom + 1
  for (i in 2:length(chromNames)) {
    chromosome[chrom] <- chrom
    scafCurrent <- subset(scafL, chromosome == chromNames[i])
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom]/2
    endLast = chromEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(chromosome,chromLengths,chromStarts,chromEnds,chromMid))
  return(table)
}

#calculate scaffold coordinates
scaf.coords <- function(scafL,gap = 0) {
  scaffold = vector()
  scafStarts = vector()
  scafEnds = vector()
  chrom = 1
  endLast = 0
  for (e in 1:nrow(scafL)) {
    scaffold[chrom] <- levels(scafL$scaffold)[e]
    scafStarts[chrom] <- endLast + gap + 1
    scafEnds[chrom] <- endLast + gap + scafL$length[e]
    endLast = scafEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(scaffold,scafStarts,scafEnds))
  return(table)
}


###
scafL<-read.table("scaffold_lengths_chrom_melp2.5.txt",h=T)
chromNames <- c(1:21)
chrom_coords <- chrom.coords(scafL,chromNames)
scaf_coords <- scaf.coords(scafL)

scaf_coords2 <- merge(scafL,chrom_coords,by="chromosome", all.x=TRUE)
###

###
# Adjust your positions
###

fileX <- read.csv("Melpomene.geno.popPairDist5000.csv")

fileX_merged <- merge(fileX, scaf_coords2,by = "scaffold", all.x=TRUE)

fileX_merged_pos <- cbind(fileX_merged, chromPos = fileX_merged$mid + fileX_merged$scafStart + fileX_merged$chromStarts-2)


###
# plot
###

# I put my files in a list. If I have multiple, it can loop and plot one underneath the other (adjust mfrow!)
comp=list(fileX_merged_pos)
names=c("fileXX_name")
col=c("black")


par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 1
bot = 0

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[21]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$Fst_P1_P2, type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=1),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Fst")), cex=0.5,line = 0.8)
}

# plot chromosome labels
axis(1, at=chrom_coords[,5][1:21]/1000000, labels=(1:21),lwd=0, lwd.ticks=0)

# plot scale (you might have to adjust this to show up)
segments(c(390,390,400), c(250.07,250.0702,250.0702), c(400,390,400) ,c(250.07,250.0702,250.0702), lwd = 1)
text(395,300.086,labels = "10Mb", cex = 1)
