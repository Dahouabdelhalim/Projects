#load the required library
library(ggplot2)
library(tidyverse)

###plot FST for independent population comparisons; the HAM x MARi population pair is given as example
HAM_MARi<-read.csv(file="HAM_MARi.windowed.weir.fst", header = TRUE, sep='\\t')
HAM_MARi_CHR7<-HAM_MARi[HAM_MARi$CHROM == '7',]

#get median values of genome wide FST (used to order the population pairs)
median(HAM_MARi$WEIGHTED_FST)

#manhattan plot genomewide using ggplot2 (code adapted from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/)
###HAM x MARi genomewide
nCHR <- length(unique(HAM_MARi$CHROM))
HAM_MARi$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(HAM_MARi$CHROM)){
  nbp[i] <- max(HAM_MARi[HAM_MARi$CHROM == i,]$BIN_START)
  HAM_MARi[HAM_MARi$CHROM == i,"BPcum"] <- HAM_MARi[HAM_MARi$CHROM == i,"BIN_START"] + s
  s <- s + nbp[i]
}
axis.set <- HAM_MARi %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)


ggplot(HAM_MARi, aes(x=BPcum, y=WEIGHTED_FST)) +
  geom_point(aes(color=as.factor(CHROM)), alpha = 0.75, size = 1.25, shape=21) +
  geom_smooth(span=0.05,se = FALSE, colour="red", size=0.8, method = "loess")+
  scale_color_manual(values = rep(c("black","grey70"), nCHR)) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  labs(x = NULL, y = "FST") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))

###HAM x MARi chromosome 7
#the outlier positions input file was obtained using the "FST_per_pop_pair_50kb.sh" script
HAM_MARi_outliers<-read.csv(file="HAM_MARi.outlier_positions.txt", header = FALSE, sep='\\t')

ggplot(data = HAM_MARi_CHR7,
       aes(x = BIN_START, 
           y = WEIGHTED_FST))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=0, ymax=1, alpha=0.1, fill="black", size=0)+
  geom_point(size=1.5, colour="black", shape=21)+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=0.8, method = "loess")+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  geom_point(data=HAM_MARi_outliers, aes(x=V1, y=1.2), shape=124, size=2, colour="red")+
  ylim(0,1.3)+
  labs(x = NULL, y = "FST") +
  theme_bw()+
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))


Cuba_average_FST_genomic<-read.csv(file="average_FST.Cuba.txt", header = TRUE, sep=' ')
Cuba_average_FST_CHR7<-Cuba_average_FST_genomic[Cuba_average_FST_genomic$CHROM == '7',]

Florida_average_FST_genomic<-read.csv(file="average_FST.Florida.txt", header = TRUE, sep=' ')
Florida_average_FST_CHR7<-Florida_average_FST_genomic[Florida_average_FST_genomic$CHROM == '7',]

###plot average FST for Cuba population comparisons, chromosomes 1-10 (Fig. 2A)
nCHR <- length(unique(Cuba_average_FST_genomic$CHROM))
Cuba_average_FST_genomic$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(Cuba_average_FST_genomic$CHROM)){
  nbp[i] <- max(Cuba_average_FST_genomic[Cuba_average_FST_genomic$CHROM == i,]$BP)
  Cuba_average_FST_genomic[Cuba_average_FST_genomic$CHROM == i,"BPcum"] <- Cuba_average_FST_genomic[Cuba_average_FST_genomic$CHROM == i,"BP"] + s
  s <- s + nbp[i]
}
axis.set <- Cuba_average_FST_genomic %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

ggplot(Cuba_average_FST_genomic, aes(x=BPcum, y=FST)) +
  geom_point(aes(color=as.factor(CHROM)), size = 1.25) +
  scale_color_manual(values = rep(c("#576cadff","#6ab7c1ff"), nCHR)) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.75)) +
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank())

###plot average FST for Florida population comparisons, chromosomes 1-10 (Fig. 2A)
nCHR <- length(unique(Florida_average_FST_genomic$CHROM))
Florida_average_FST_genomic$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(Florida_average_FST_genomic$CHROM)){
  nbp[i] <- max(Florida_average_FST_genomic[Florida_average_FST_genomic$CHROM == i,]$BP)
  Florida_average_FST_genomic[Florida_average_FST_genomic$CHROM == i,"BPcum"] <- Florida_average_FST_genomic[Florida_average_FST_genomic$CHROM == i,"BP"] + s
  s <- s + nbp[i]
}
axis.set <- Florida_average_FST_genomic %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

ggplot(Florida_average_FST_genomic, aes(x=BPcum, y=FST)) +
  geom_point(aes(color=as.factor(CHROM)), size = 1.25) +
  scale_color_manual(values = rep(c("#576cadff","#6ab7c1ff"), nCHR)) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.75)) +
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank())


####plot average FST for Cuba and Florida population comparisons, chromosome 7 (Fig. 2B)
ggplot(data = Cuba_average_FST_CHR7,
       aes(x = BP, 
           y = FST))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=0, ymax=0.75, alpha=0.1, fill="black", size=0)+
  geom_point(colour="#576cadff",shape=21, size=1)+
  geom_point(data=Florida_average_FST_CHR7, colour="#576cadff", size=1)+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=1, method = "loess")+
  geom_smooth(data=Florida_average_FST_CHR7, span=0.2,se = FALSE, colour="black", size=1, method = "loess")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  labs(x = NULL, y = NULL) +
  theme_bw()+
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank())


###plot average FST for Cuba population comparisons, chromosomes 1-10, pairs selected based on geography
Cuba_average_FST_genomic_geography<-read.csv(file="average_FST.Cuba.pairs_by_geography.txt", header = TRUE, sep=' ')
Cuba_average_FST_CHR7_geography<-Cuba_average_FST_genomic_geography[Cuba_average_FST_genomic_geography$CHROM == '7',]


nCHR <- length(unique(Cuba_average_FST_genomic_geography$CHROM))
Cuba_average_FST_genomic_geography$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(Cuba_average_FST_genomic_geography$CHROM)){
  nbp[i] <- max(Cuba_average_FST_genomic_geography[Cuba_average_FST_genomic_geography$CHROM == i,]$BP)
  Cuba_average_FST_genomic_geography[Cuba_average_FST_genomic_geography$CHROM == i,"BPcum"] <- Cuba_average_FST_genomic_geography[Cuba_average_FST_genomic_geography$CHROM == i,"BP"] + s
  s <- s + nbp[i]
}
axis.set <- Cuba_average_FST_genomic_geography %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

library(ggplot2)
ggplot(Cuba_average_FST_genomic_geography, aes(x=BPcum, y=FST)) +
  geom_point(aes(color=as.factor(CHROM)), alpha = 0.75, size = 1.25) +
  scale_color_manual(values = rep(c("#576cadff","#6ab7c1ff"), nCHR)) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.75)) +
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank())



####plot average FST for Cuba population comparisons, chromosome 7, pairs selected based on geography
ggplot(data = Cuba_average_FST_CHR7_geography,
       aes(x = BP, 
           y = FST))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=0, ymax=0.75, alpha=0.1, fill="black", size=0)+
  geom_point(colour="#2e4057",shape=21, size=1)+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=1, method = "loess")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  labs(x = NULL, y = NULL) +
  theme_bw()+
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank())


####barplots of per-chromosome outlier content
results_geography <- read.csv("results_outlier_proportions_geography.csv", sep = ',')
results_geography$Chromosome = factor(results_geography$Chromosome, levels=c("1", "2","3","4","5","6","7","8","9","10"))

library(ggplot2)
ggplot(data=results_geography, aes(Chromosome,Percent_outlier_windows_5percent))+
  geom_bar(stat="identity",fill="white",colour="black", width=0.8)+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")+
  labs(x = NULL, y = "Percent")

library(ggplot2)
ggplot(data=results_geography, aes(Chromosome,Percent_outlier_windows_1percent))+
  geom_bar(stat="identity",fill="white",colour="black", width=0.8)+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")+
  labs(x = NULL, y = "Percent")

