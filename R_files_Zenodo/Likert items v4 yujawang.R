### Wildlife gardening: an urban nexus of social and ecological relationships
### Laura Mumaw and Luis Mata
### https://ecoevorxiv.org/9rkhm/
### Code by Luis Mata | v4 yujawang | 16 February 2021

# Read in the data
data = read.csv("data.csv", header=TRUE, sep=",", na.strings=TRUE)

# Data wrangling
w = as.data.frame(data$well)
for (i in 1:10){
  w[,i+1] = data[,54+i]
}  
w

ww = matrix(NA, 5, 10)

for (i in 1:10){
  ww[1,i] = sum(w[,i+1]==-10)
  ww[2,i] = sum(w[,i+1]==-5)
  ww[3,i] = sum(w[,i+1]==-0)
  ww[4,i] = sum(w[,i+1]==5)
  ww[5,i] = sum(w[,i+1]==10)
}  
ww
l = dim(w)[1]
xx = (ww*100)/l 
xx

colnames(xx) = c("kcomm","joinpro","attgarden","attknox","attnat","sensepur",
                 "helpwild","pridegar","learnnew","posfeel")
rownames(xx) = c("-10","-5","0","5","10")

zz = xx[,c("helpwild","posfeel","attnat","sensepur","pridegar","attgarden","learnnew","attknox","joinpro","kcomm")]

zz = zz[c(5,4,3,2,1),c(10:1)]

# Plot
folder = "likertitems.jpg"
jpeg(filename=folder, width=1100, height=700)

  par(mai=c(1,2,1,2), cex=1)
  pal = c("royalblue","dodgerblue","grey77","grey55","grey33")  
  barplot(zz, col=pal , border="grey77", xlab="", axes=F, xlim=c(0,100), horiz=T)
  axis(1, at=c(0,20,40,60,80,100), labels=c(0,"20%","40%","60%","80%","100%"), 
       las=1, cex.axis=2.5, padj=.5)

dev.off()

# Legend
folder = "legendlikertitems.jpg"
jpeg(filename=folder, width=1100, height=700)

  par(mai=c(1,2,1,2), cex=1)
  barplot(zz, col=c("white","white","white","white","white") , border="white", xlab="", 
          ylab="", axes=F, xlim=c(0,100), horiz=T, axisnames=F)
  legend(x=5, y=12, legend=rownames(zz), fill=pal, bty="n", cex=3.5,
         horiz=T, text.col="white")

dev.off()


