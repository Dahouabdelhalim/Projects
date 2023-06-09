###############
library(RColorBrewer)
library(gplots)
library('corrplot')

make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
} 

wraplabels <- function(longnames, w=10) {
  retval=c()
  for (n in longnames) {
    sn = paste(strwrap(n, width = w), collapse="\\n")
    sn = substr(sn,0,30)
    retval=c(retval, sn)
  }
  return(retval)
}

trunclabs <- function(longnames, stopAt=30) {
  retval=c()
  for (n in longnames) {
    sn=n
    if (nchar(n)>stopAt) {
      sn = paste(substr(n, 0, stopAt-3), "...")  
    }
    if(sn=="")
      sn="N/A"
    retval=c(retval, sn)
  }
  return(retval)
}

###############

###  START  ###
print("Importing data")
data = read.csv("selected-results.csv", header = TRUE, sep=",", quote = "\\"", check.names=FALSE)

### categorical ordinal answers, answers replaced by integer code to enable correlation analysis
dataCORR = read.csv("catordQUCODANSFULL_process.csv", header = TRUE, sep = ",",
 stringsAsFactors = TRUE, quote = "\\"", check.names=FALSE)

paste("Discovered", nrow(data), "completed responses answering", ncol(data), "questions each")

### OUTPUT FOLDER  ###
print("Creating output folder")
make.dir("out")
setwd("out")

#https://stackoverflow.com/questions/23913276/texture-in-barplot-for-7-bars-in-r
cols = c("firebrick2", "blue3","forestgreen","purple2","darkorange2","darkturquoise","cornsilk4","navy")

# questions' strings manipulation
qnames = names(data)
fnames = make.names(qnames)

for (i in 1:length(data)) {
  # again strings manipulation
  fullqn = qnames[i]
  codes = unlist(strsplit(fullqn, "[.]"))[1]
  qtitle = substr(fullqn, nchar(codes)+3, nchar(fullqn))
  
  seccod = substr(codes,0,2)
  qcod = substr(codes,3,nchar(codes))

  qn = substr(paste(strwrap(qtitle, width = 60), collapse = "\\n"),0,160)
  qn = ""
  
  fn = fnames[i]
  fn = substr(gsub("\\\\.", "_", fn),0,25)
  
  # Section folder
  if(!file.exists(seccod))
    dir.create(seccod)
  
  cdir = paste(seccod,"/",qcod,sep="")
  make.dir(cdir) 
  setwd(cdir)
  
  toplot = sort(table(data[i]), decreasing = TRUE)
  rn = rownames(toplot)
  rn = trunclabs(rn,250)
  mylabels = wraplabels(rn,10)
  leglabs = trunclabs(rn,40)
  
  # Plot histograms with and without titles
  for (dev in 1:2) {
    
    if(dev==1) {
      pdf(paste(codes,"_h.pdf",sep=""))
    }
    if(dev==2) {
      qn = qnames[i]
      qn = substr(paste(strwrap(qn, width = 60), collapse ="\\n"),0,160)
      pdf(paste(codes,"_hWT.pdf",sep=""))
    }
    
    barplot(toplot/sum(toplot)*100, 
            las=2,
            main=qn,
            ylab = "Percentage [%]",
            ylim=c(0,100),
            col=cols,
            legend = leglabs,
            names.arg="",
            angle=c(45,0,90), density=seq(10,30,20)
    )
    dev.off()
  }
  setwd("../..")
}
  
# Preparing folders for extra plots
make.dir("extra")

setwd("./extra")
make.dir("corrplotMAT")
make.dir("corrplotUPPER")
make.dir("heatmap_rainbow")

### correlation diagrams
M <- cor(dataCORR)

for (i in 1:2) {
  
  if (i %% 2) {
    mytitle = "Pearson Correlation Index"
    myfname = "corrplotMAT/corrplotMAT_WT.pdf"
  } else {
    mytitle = ""
    myfname = "corrplotMAT/corrplotMAT.pdf"
  }
  pdf(myfname)
  par(xpd=TRUE)
  corrplot(M, method = "square", col = colorRampPalette(c("navyblue","white","red3"))(100),
           tl.col = "black",title=mytitle, mar = c(2, 2, 2, 2), tl.cex=0.7)
  dev.off()
  
  if (i %% 2) {
    mytitle = "Pearson Correlation Index"
    myfname = "corrplotUPPER/corrplotUPPER_WT.pdf"
  } else {
    mytitle = ""
    myfname = "corrplotUPPER/corrplotUPPER.pdf"
  }
  pdf(myfname)
  corrplot(M, type = "upper", order = "hclust", tl.col = "black",
           title=mytitle, tl.srt = 75, tl.cex=0.7, mar = c(2, 2, 2, 2))
  dev.off()
  
  if (i %% 2) {
    mytitle = "Heatmap"
    myfname = "heatmap_rainbow/heatmap_rainbow_WT.pdf"
  } else {
    mytitle = ""
    myfname = "heatmap_rainbow/heatmap_rainbow.pdf"
  }
  pdf(myfname)
  par(mar=c(7,4,4,2)+0.1) 
  heatmap.2(x=M, Rowv=NULL,Colv="Rowv",
            main=mytitle,
            cexRow=1,
            cexCol=1.2,
            col = rev(rainbow(20*10, start = 0/6, end = 4/6)), 
            scale="none",
            margins=c(8,2.5), # ("margin.Y", "margin.X")
            trace='none', 
            symkey=FALSE, 
            symbreaks=FALSE, 
            dendrogram='none',
            density.info='histogram', 
            denscol="black",
            keysize=1, 
            key.par=list(mar=c(3.5,0,3,0)),
            lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
  )
  dev.off()
  
}