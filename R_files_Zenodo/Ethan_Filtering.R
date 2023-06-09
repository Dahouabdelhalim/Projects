BiocManager::install("phyloseq")
BiocManager::install("IRanges")

library(phyloseq)
library(ggplot2)
library(vegan)

otu_mat<-read.csv("ethan_otu_table2.csv", row.names = 1)
#stay the same
tax<-read.csv("ethanhollow_edna_species7.csv", row.names = 1)
tax
tax_mat<-as.matrix(tax)
#done
sample_df<-read.csv("ethanhollow_samplelist.csv", row.names = 1)
#done


##Make phyloseq object
otu<-otu_table(otu_mat, taxa_are_rows = T)
tax<-tax_table(tax_mat)
samples<-sample_data(sample_df)

holl<-phyloseq(otu, tax, samples)
holl
sample_data(holl)
tax_table(holl)

#Check ASVs found in extraction control
extracts<-subset_samples(holl, Type =="Extract")
taxa_sums(extracts)<2
keeptaxa<-taxa_sums(extracts)<2
holl2<-prune_taxa(keeptaxa, holl)
holl2<-subset_samples(holl2, Type != "Extract")
holl2
#These ASVs are human and will be removed in the next step anyway

#remove low abundance samples
#Check rarefaction curve
#rarecurve(t(otu_table(holl)), step=50, cex=0.5)
#rarecurve(t(otu_table(holl2)), step=50, cex=0.5, xlim=c(0, 10000))
sums<-as.data.frame(sample_sums(holl))
holl2<-subset_samples(holl, sample_sums(holl)>642)

#rarecurve(t(otu_table(holl2)), step=50, cex=0.5, xlim=c(0, 10000))
#rarecurve(t(otu_table(holl2)), step=50, cex=0.5)


#Filtering by Relative abundance within samples, taxa making up less than 0.5% of the sample have counts reverted to 0
#You can change the level of filtering here, or just remove singletons later if you think this is too strict
holl2.df<-psotu2veg(holl2)

holl2.df<-t(holl2.df)
holl2.df=data.frame(holl2.df)
colnames(holl2.df)
holl2.df2=holl2.df
for (i in 1:ncol(holl2.df)){
  
  holl2.df[,i]
  temp=holl2.df[,i]
  which(temp<sum(holl2.df[,i])*0.005)
  temp[which(temp<sum(holl2.df[,i])*0.005)] =0
  holl2.df[,i]=temp
  
}


##Make phyloseq object again
otu2<-otu_table(holl2.df, taxa_are_rows = T)

holl2.fil<-phyloseq(otu2, tax_table(holl2), sample_data(holl2))
sums2<-as.data.frame(sample_sums(holl2.fil))
holl2
tax_table(holl2.fil)

#Select only Chordata ASVs
holl3<-subset_taxa(holl2.fil, phylum=="Chordata")
holl3

#Remove contamination ASVs
holl3<-subset_taxa(holl3, Contamination=="N")
holl3
sums3<-as.data.frame(sample_sums(holl3))
#Also remove A LOT of sequences
sum(sums3)#5,616,468
sum(sums)#15,849,075



#Remove the now extract samples
holl4<-subset_samples(holl3, Type != "Extract")


################################################################
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
####################################################################

#########Aggregate taxa
holl5<-tax_glom(holl4, taxrank="final_name")
holl5
tax_table(holl5)
taxa_names(holl5)<-tax_table(holl5)[,9]
taxa_names(holl5)
?taxa_names
#attempting to duplicate multiple taxa after transposing holl2
#At this point holl4 is the filtered frame with ASVs, and holl5 is with the taxa aggregated

###########################################
#Transformations
library(microbiome)
holl5.log<-transform(holl5, "log")
holl5.pa<-transform(holl5, "pa")
holl4.log<-transform(holl4, "log")
holl4.pa<-transform(holl4, "pa")
###########################################

#Load camera data
cam.otu<-read.csv("ethan_camera_table2.csv", row.names=1)
cam.otu #picumnus present
cam.env<-read.csv("ethan_camera_map.csv", row.names=1)
cam.env
cam.tax<-read.csv("cam_tax_filled3.csv", row.names = 1)
cam.tax
cam.tax_mat<-as.matrix(cam.tax)

##Make phyloseq object
cam.otu<-otu_table(cam.otu, taxa_are_rows = T)
cam.tax<-tax_table(cam.tax_mat)
cam.samples<-sample_data(cam.env)

cam<-phyloseq(cam.otu, cam.tax, cam.samples)
cam
tax_table(cam)
cam.pa<-transform(cam, "pa")
cam.tax
cam.otu
cam.samples
