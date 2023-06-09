library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)

#make sure an input is given
args<-commandArgs(trailingOnly = T)
if(length(args)==0){
  print("Input file missing")
  q()
}
fileName<-args[1]

#read in the populations.txt file
files<-read.csv(file=fileName, header = FALSE)
poplist<-dplyr::pull(files, V1)

##calculate mean R2 for increasing distances in increments of 20kb (script following https://www.biostars.org/p/300381/)
for (i in poplist){
  name<-paste (i,".maf_20.ld.summary",sep="")
  output<-paste(i,".txt",sep="")
  dfr <- read.delim(name, sep="",header=F)
  colnames(dfr) <- c("dist","rsq")
  dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=20000))
  dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq))
  dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\\\(\\\\)\\\\[\\\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\\\(\\\\)\\\\[\\\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  write.csv(dfr1, file=output,row.names = FALSE)
}