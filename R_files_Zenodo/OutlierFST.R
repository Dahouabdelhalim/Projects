#  title: "Process_FST"
#author: "Giorgia Auteri"
#date: "August 5, 2018"
#Script modified from Auteri, Giorgia G. and L. Lacey Knowles. 2020. Decimated little brown bats
# show potential for adaptive change. Scientific reports: 10, 3023. doi https://doi.org/10.1038/s41598-020-59797-4

########################################################################################
#Overvieww and processing of the fst data

#read in data
fdat <- read.table("populations.fst_1-6.tsv",
                   header = TRUE, 
                   sep = "\\t")
colnames(fdat) <- c('LocusID',	'Pop1',	'Pop2',	'CHROM',	'POS', 'Column',	'OverallPi',	'AMOVAFst',	'FishersP',	'OddsRatio',	'CILow',	'CIHigh.All',	'LOD',	'CorrectedAMOVAFst', 'SmoothedAMOVAFst',	'SmoothedAMOVAFst P-value',	'WindowSNPCount')

#Check
dim(fdat) #row number is number of loci/SNPs

length(unique(fdat$LocusID))

#Distribution of fst values
hist(fdat$AMOVAFst, breaks=1000); abline(v=0.0045, col = "red"); abline(v=0.05, col = "blue")

mFst <- mean(fdat$AMOVAFst);mFst  #mean fst
sdFst <- sd(fdat$AMOVAFst); sdFst #one standard deviation of Fst
cutoff9 <- mFst + (sdFst * 9); cutoff9 #the cutoff for three SDs from the mean

summary(fdat) #overall summary

############################################################################################
#Plot Fst values

library(ggplot2) 

scaffolds <- unique(fdat$CHROM) #list of unique scaffold names

chrom.cols <- c(rep(c("#F95045" , "#3CE0FF"), floor(length(scaffolds)/2))) #alternating colors for scaffolds

cutoff <- cutoff9 #Assign cutoff

#Making highlighter lines
selLoc <- (fdat$LocusID[fdat$AMOVAFst >= cutoff])/1000000
highlights <- data.frame(LocusID=selLoc, AMOVAFst=c(rep(0.7, length(selLoc))))

ggplot(fdat, aes(x=LocusID/1000000, y=(AMOVAFst), colour=CHROM)) +
  geom_point(shape = 19, size = 2)+
  scale_y_continuous(limits=c(0,0.6))+
  scale_color_manual(values=chrom.cols)+
  xlab("Alignment position (Mbp)")+
  ylab(bquote(paste(italic('F'['ST']*' ')))) +
  geom_hline(yintercept = c(cutoff), linetype = "dashed", size = 1, color = "gray35") +
  theme(legend.position="none", text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#################################################################################################
#Write output file with selected loci

selectedLoci <- na.omit(fdat[fdat$AMOVAFst >= cutoff,])

dim(selectedLoci)

write.csv(selectedLoci, "SelectedLoci.csv")
