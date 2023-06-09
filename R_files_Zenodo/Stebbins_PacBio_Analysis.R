#load packages
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(phyloseq); packageVersion("phyloseq")
library("microbiome")

#set local working directory, load phyloseq objects (no meta data)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
bact16S<-readRDS("16S_DADA2output.RDS")
fungITS<-readRDS("ITS_DADA2output.RDS")

#load meta data
BBSCC_meta<-read.csv("sampleSubmit_RECORD_2.csv")
SF_meta<-read.csv("Sunflower collections - Sheet1.csv")

#clean SCC and BB meta data
BBSCC_meta$sample_submit_ID
BBSCC_meta$sample_submit_ID<-gsub( "\\\\.","-", BBSCC_meta$sample_submit_ID)
BBSCC_meta$Site

#loop through each samplea and pull relevant data
out<-list()
for(i in unique(BBSCC_meta$sample_submit_ID)){
  
  sub<-BBSCC_meta[BBSCC_meta$sample_submit_ID==i , ]
  
  numFlwrSp<-nrow(sub)
  
  date<-as.POSIXct(sub$Date, format = "%m/%d/%y")
  
  site<-sub$Site
  
  out[[i]]<-data.frame(sample = i, numFlwrSp=numFlwrSp[1], date = date[1], site = site[1])
}
#put together data into single data.frame
summaryBBSCC<-do.call(rbind , out)

#sample_data for non sunflower
sample_data(bact16S)<-sample_data(summaryBBSCC[match(rownames(otu_table(bact16S)) , summaryBBSCC$sample),])
sample_data(fungITS)<-sample_data(summaryBBSCC[match(rownames(otu_table(fungITS)) , summaryBBSCC$sample),])


#get rid of plant DNA
fungITS<-subset_taxa(fungITS , Kingdom!="k__Viridiplantae")


# look at contamination
data<-psmelt(fungITS)
p<-ggplot(data=data, aes(x= sample, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack")+facet_wrap(site~1, scales="free")
contam <- colnames(otu_table(fungITS)[,which(otu_table(fungITS)["PCR-C",]>1) ])
fungITS_contam<-prune_taxa(contam, fungITS)
tax_table(fungITS_contam)

# removes Acremonium, Metschnikowia gruessii, Penicillium, Mycosphaerella, and Acremonium and Aspergillus. I wonder if this is spillover(?)


data<-psmelt(subset_taxa(bact16S , Phylum!="Cyanobacteria"))
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") 


#remove anything that was in PCR control (6 taxa)
Notcontam<-colnames(otu_table(fungITS)[,-which(otu_table(fungITS)["PCR-C",]>1) ])
fungITS_prune<-prune_taxa(Notcontam, fungITS)
fungITS_prune<-subset_taxa(fungITS_prune ,  Kingdom!="k__Viridiplantae")
fungITS_prune<-subset_taxa(fungITS_prune ,  Kingdom=="k__Fungi")

relfungITS <- transform_sample_counts(fungITS_prune, function(x) x / sum(x))

sum(otu_table(fungITS_prune))/sum(otu_table(fungITS))

setwd("..")
data<-psmelt(relfungITS)
data<-data[-which(data$site!="BB"&data$site!="SCC"),]
p <- ggplot(data=data, aes(x=date, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack")+facet_grid(site~1)
ggsave("figures/siteXdate_ITS.pdf", width = 18, height = 6)



p <- ggplot(data=data, aes(x=Kingdom, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack")
ggsave("figures/summary_ITS.pdf")

p <- ggplot(data=data, aes(x=Kingdom, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack")
ggsave("figures/summary16S.pdf")

# which are the top 5 genera? 
TopNOTUs = names(taxa_sums(relfungITS)[1:9])
ent10 = prune_taxa(TopNOTUs, relfungITS)
plot_bar(ent10, fill = "Genus")

##########
bact16S_prune<-subset_taxa(bact16S ,  Genus!="NA")

relbact16S <- transform_sample_counts(bact16S_prune, function(x) x / sum(x))
data<-psmelt(relbact16S)
p <- ggplot(data=data, aes(x=date, y=Abundance, fill=OTU))
p + geom_bar(aes(), stat="identity", position="stack")+#facet_wrap(~Genus~1)+
  theme(legend.position="none")
ggsave("figures/siteXdate_16S.pdf")

data<-psmelt(relbact16S)
p <- ggplot(data=data, aes(x=date, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack")+facet_grid(site~1)

p <- ggplot(data=data, aes(x=Kingdom, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack")
ggsave("figures/summary16S.pdf")

