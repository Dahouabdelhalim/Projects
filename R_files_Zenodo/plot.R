library(ggplot2)
library(reshape2)
library(tidyverse)
library(hrbrthemes)

chrorder <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
             "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
             "chr18","chr19","chr20","chr21","chr22","chrX")

dfhm <- read.csv("9606.protein.links.v11.5_zscores_edges.txt",
                sep="\\t", skip=1, header=TRUE)
colnames(dfhm) <- c("var1", "var2", "zscore")
ggplot(dfhm, aes(x=var1, y=var2, fill=zscore))+
        geom_tile()+
        scale_fill_distiller(name="z-score", palette="Spectral")+
        scale_x_discrete(limits=chrorder)+
        scale_y_discrete(limits=chrorder)+
        theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1),
            legend.position="left",
            legend.key.height = unit(1.5, "cm"))


dfhg <- read.csv("9606.protein.links.v11.5_zscores_singletons.txt",
                sep="\\t", skip=1, header=TRUE)
colnames(dfhg) <- c("chr", "n", "m", "zscore")
ggplot(dfhg, aes(x=chr, y=zscore))+
        geom_bar(stat="identity", fill="#69b3a2", color="#e9ecef", alpha=0.9)+
        scale_x_discrete(limits=chrorder)+
        theme_ipsum()+
        theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=10, angle=45, hjust=1))+
        xlab("z-score")+
        coord_flip()
