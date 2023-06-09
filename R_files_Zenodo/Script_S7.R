# This is code to replicate the analysis from Alleman et al. "Tandem-Running and Scouting Behavior is Characterized by Up-Regulation of Learning and Memory formation genes within the Ant Brain"
# Homologue analysis between T. americanus and T. longispinosus
# Transcriptomes constructed using Trinity
# DEGs obtained by using DESeq2
# Orthogroups created using Orthofinder
# Script by M. Stoldt

#####################################################################################################################
## Loading libraries
#####################################################################################################################

library(lme4)
library(car)
library(MASS)
library(ggplot2)
library(multcomp)
library(plyr)
library(VennDiagram)

#####################################################################################################################
################################################### Part I ##########################################################
######################################## Orthogroups between transcriptomes #########################################
#####################################################################################################################

#####################################################################################################################
## Extraction of the single-copy orthogroups found by Orthofinder
#####################################################################################################################

# Specify working directory for project
setwd("/Directory/Orthofinder/Transcriptomes/")

# Read in file with the number of contigs per species per orthogroup
geneCount<-read.table("Orthogroups.GeneCount.csv", header=TRUE, sep=",", quote="")

# Take only the ones where both species has only 1 contig present
geneCountSC<-subset(geneCount, geneCount$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1 & geneCount$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1)
contigsSC<-as.data.frame(geneCountSC$X)

# Read in file with contigs assigned to each orthogroup
orthogroups<-read.table("Orthogroups.csv", sep=",",
                        quote="", header=TRUE)
# take only the orthogroups (and according contigs) that we extracted before
SC_Orthogroups<-orthogroups[match(contigsSC$`geneCountSC$X`, orthogroups$X),]
colnames(SC_Orthogroups)<-c("X", "Y","Z")

#####################################################################################################################
## Extraction of all the remaining orthogroups
#####################################################################################################################

# Take only the ones where at least one species has only 1 contig present
geneCount_Tam<-subset(geneCount,geneCount$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1&geneCount$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder>1)
geneCount_TL<-subset(geneCount,geneCount$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1&geneCount$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder>1)
contigs_rest_tam<-as.data.frame(geneCount_Tam$X)
contigs_rest_tl<-as.data.frame(geneCount_TL$X)

# Take only the orthogroups (and according contigs) that we extracted before
Orthogroups_rest_tam<-orthogroups[match(contigs_rest_tam$`geneCount_Tam$X`,orthogroups$X),]

Orthogroups_rest_tl<-orthogroups[match(contigs_rest_tl$`geneCount_TL$X`,orthogroups$X),]
write.table(Orthogroups_rest_tl,file="Orthogroups_OneToMany_TL.csv", sep="\\t",row.names = FALSE)

# Now edit Orthogroups_OneToMany_TL manually and move Tlongi-entry to the start
Orthogroups_rest_tl_edit<-read.table("Orthogroups_OneToMany_TL_edited.csv",
             sep=",", quote = "", header = TRUE)

# now concatenate Orthogroup-files again
Orthogroups_rest_new<-rbind(Orthogroups_rest_tam, Orthogroups_rest_tl_edit)

#####################################################################################################################
## Creation of new single-copy orthogroups from the orthogroups extracted before
#####################################################################################################################

# Read in information which contig is assigned to which ID 
# Before opening delete uninteresting information (we just need ID and contig name) and replace ": " by ":"
contigToID<-read.table("SequenceIDs.txt",header=FALSE, sep=" ", quote="")

# Read in BLAST files
blast_1<-read.table("Blast0_1.txt", sep="\\t", quote="")
blast_2<-read.table("/Blast1_0.txt",sep="\\t", quote="")
# Only take the necessary information: query, hit and score
blast_1_1<-cbind(blast_1[1], blast_1[2], blast_1[12])
colnames(blast_1_1)<-c("query", "hit", "score")
blast_2_1<-cbind(blast_2[1], blast_2[2], blast_2[12])
colnames(blast_2_1)<-c("query", "hit", "score")

# Create dataframe to write new single-copy orthogroups in
new_orthogroups<-as.data.frame(Orthogroups_rest_new$V2)
new_orthogroups$V2<-c(rep(1:length(new_orthogroups)))

# Now assign to each single copy species entry the entry from the other species with the highest blast score
for (i in 1:length(Orthogroups_rest_new_edit$V1) ){
  # Get ID of the contig of interest
  blast_x<-as.character(contigToID$V1[match(as.character(Orthogroups_rest_new_edit$V2[i]),as.character(contigToID$V2))])
  # Get all the entries in the BLAST list whose query matches with the single copy contig
  blast_list<-blast[which(!is.na(match(as.character(blast$query),blast_x))),]
  # Order entries by hit score from lowest to highest
  blast_x_ordered<-blast_list[order(blast_list$score),]
  if(length(blast_list$query)>0){
    for(j in 1:length(blast_x_ordered$query)){
      blast_hit<-
        contigToID$V2[match(as.character(blast_x_ordered$hit[j]),as.character(contigToID$V1))
                      ]
      if(as.character(blast_hit) %in% t(as.matrix(Orthogroups_rest_new[i,]))){
        best_ortho<-as.character(Orthogroups_rest_new_edit[i,match(blast_hit, t(as.matrix(Orthogroups_rest_new_edit[i,])))])new_orthogroups$V2[i]<-best_ortho
      }
    }
  }
}

colnames(new_orthogroups)<-c("X", "Y", "Z")

# Combine lists of orthogroups
all<-rbind(SC_Orthogroups, new_orthogroups)

#####################################################################################################################
################################################### Part II #########################################################
###################### Orthogroups between DEGs of T. americanus and T. longispinosus transcriptome #################
#####################################################################################################################

#####################################################################################################################
## Extraction of the single-copy orthogroups found by Orthofinder
#####################################################################################################################

# Specify working directory for project
setwd("/Directory/Orthofinder/DEG_Tam_vs_Tlongi_Transcriptome/")

# Read in file with the number of contigs per species per orthogroup
geneCount<-read.table("Orthogroups.GeneCount.csv", header=TRUE, sep="\\t", quote="")

# Take only the ones where both species has only 1 contig present
geneCount_sc<-subset(geneCount, geneCount$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1 & geneCount$degs_sequences_tam==1)

contigsSC<-as.data.frame(geneCount_sc$X)

# Read in file with contigs assigned to each orthogroup
orthogroups<-read.table("Orthogroups.csv", sep="\\t", quote="", header=TRUE)

# Take only the orthogroups (and according contigs) that we extracted before
SC_Orthogroups<-orthogroups[match(contigsSC$`geneCount_sc$X`, orthogroups$X),]
colnames(SC_Orthogroups)<-c("X", "Y", "Z")

#####################################################################################################################
## Extraction of all the remaining orthogroups
#####################################################################################################################

# Take only the ones where at least one species has only 1 contig present
geneCount_rest<-subset(geneCount, xor(geneCount$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1, geneCount$degs_sequences_tam==1))

contigs_rest<-as.data.frame(geneCount_rest$X)

# Now create a data structure with all single copy orthologues in first column
geneCount_sc_tam<-subset(geneCount, geneCount$degs_sequences_tam ==1&geneCount$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder>1)

contigs_sc_tam<-as.data.frame(geneCount_sc_tam$X)

# Take only the orthogroups (and according contigs) that we extracted before (have to cheat a little to combine tables ;))
Orthogroups_sc_tam<-orthogroups[match(contigs_sc_tam$`geneCount_sc_tam$X`, orthogroups$X),]
Orthogroups_sc_tam_invert<-cbind(as.character(Orthogroups_sc_tam$X), as.character(Orthogroups_sc_tam$degs_sequences_tam), as.character(Orthogroups_sc_tam$Tlongi_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder))

#####################################################################################################################
## Creation of new single-copy orthogroups from the orthogroups extracted before
#####################################################################################################################

# Read in information which contig is assigned to which ID (important for BLAST output later)
# Before opening delete uninteresting information (we just need ID and contig name) and replace ": " by ":"
contigToID<-read.table("SequenceIDs_edit.txt", header=FALSE, sep=":", quote="")

# Read in BLAST files
difftam_tl<-read.table("Blast1_0.txt", sep="\\t", quote="")
difftl_tam<-read.table("Blast0_1.txt", sep="\\t", quote="")

# Only take the necessary information: query, hit and score
difftam_tlx<-cbind(difftam_tl[1], difftam_tl[2], difftam_tl[12])
colnames(difftam_tlx)<-c("query", "hit", "score")
difftl_tamx<-cbind(difftl_tam[1], difftl_tam[2], difftl_tam[12])
colnames(difftl_tamx)<-c("query", "hit", "score")

# Combine both BLAST outputs
blast_2<-rbind(difftam_tlx, difftl_tamx)

# Create empty dataframe for new orthogroups
new_orthogroups<-Orthogroups_sc_tam_invert[,1:3]
new_orthogroups[,3]<-c(rep(1:length(new_orthogroups[,3])))

Orthogroups_sc_tam_invert<-as.data.frame(Orthogroups_sc_tam_invert)

# Now assign to each single copy species entry the entry from the other species with the highest blast score
for (i in 1:length(Orthogroups_sc_tam_invert[,2]) ){
  # Get ID of the contig of interest
  blast<-as.character(contigToID$V1[match(as.character(Orthogroups_sc_tam_invert$V2[i]),as.character(contigToID$V2))])
  
  # Get all the entries in the blast list whose query matches with the single copy contig
  blast_list<-blast_2[which(!is.na(match(as.character(blast_2$query),blast))),]
  
  # Order entries by hit score from highest to lowest
  blast_ordered<-blast_list[order(blast_list$score),]
  
  for(j in 1:length(blast_ordered$query)){
    
    blast_hit<-contigToID$V2[match(as.character(blast_ordered$hit[j]),as.character(contigToID$V1))]
    if(as.character(blast_hit) %in% strsplit(t(as.matrix(as.data.frame(Orthogroups_sc_tam_invert[i,]$V3))), split=", ")[[1]]){
      
      best_ortho<-strsplit(t(as.matrix(Orthogroups_sc_tam_invert[i,]$V3)),split=", ")[[1]][match(as.character(blast_hit), strsplit(t(as.matrix(Orthogroups_sc_tam_invert[i,]$V3)),split=", ")[[1]])]
      new_orthogroups[i,3]<-as.character(best_ortho)
    }
  }
}

colnames(new_orthogroups)<-c("X", "Y", "Z")

# Combine orthogroups
all<-rbind(SC_Orthogroups, new_orthogroups)

####################################################################################################################
## Check for matches with differentially expressed contigs in other species
####################################################################################################################

# Read in the files containing the all DEGs in T. longispinosus
deg_follower_tl<-read.table("DESeq2_All_Follower_Tlongi.csv", header = T)
deg_leader_tl<-read.table("DESeq2_All_Leader_Tlongi.csv", header = T)
deg_scout_tl<-read.table("DESeq2_All_Scout_Tlongi.csv", header=T)

# Read in the files containing the all DEGs of T. americanus
deg_follower_tam<-read.table("DESeq2_All_Follower_Tam.csv", header = T)
deg_leader_tam<-read.table("DESeq2_All_Leader_Tam.csv", header = T)
deg_scout_tam<-read.table("DESeq2_All_Scout_Tam.csv", header=T)

x<-which(match(deg_follower_tam$ID, all$Y)!="NA")
orthos_follower_tam<-all[x,]
y<-which(match(deg_leader_tam$ID, all$Y)!="NA")
orthos_leader_tam<-all[y,]
z<-which(match(deg_scout_tam$ID, all$Y)!="NA")
orthos_scout_tam<-all[z,]

# Plot Venn diagrams depicting overlap between DEGs in T. longispinosus and orthologues of DEGs in T. americanus
png(file="Follower_Tam_Orthologues_vs_DEGs_Tlongi.png", width = 820, height = 830)
venn.plot<-venn.diagram(list(deg_follower_tl$ID, orthos_follower_tam$V3), NULL, lwd=c(2,2), col=c("black","black"),category.names=c("DEGs TL","Orthos Tam"),fill=c("midnightblue", "mediumvioletred"), alpha=c(0.6, 0.6), cex=3, cat.fontface=4, cat.cex=c(1,1))
grid.draw(venn.plot)
dev.off()

png(file="Leader_Tam_Orthologues_vs_DEGs_Tlongi.png", width = 820, height = 830)
venn.plot<-venn.diagram(list(deg_leader_tl$ID, orthos_leader_tam$V3), NULL, lwd=c(2,2), col=c("black","black"),category.names=c("DEGs TL","Orthos Tam"),fill=c("midnightblue", "mediumvioletred"), alpha=c(0.6, 0.6), cex=3, cat.fontface=4, cat.cex=c(1,1))
grid.draw(venn.plot)
dev.off()

png(file="Scout_Tam_Orthologues_vs_DEGs_Tlongi.png", width = 820, height = 830)
venn.plot<-venn.diagram(list(deg_scout_tl$ID, orthos_scout_tam$V3), NULL, lwd=c(2,2), col=c("black","black"),category.names=c("DEGs TL","Orthos Tam"),fill=c("midnightblue", "mediumvioletred"), alpha=c(0.6, 0.6), cex=3, cat.fontface=4, cat.cex=c(1,1))
grid.draw(venn.plot)
dev.off()

####################################################################################################################
## Plot expression of the orthologues of the DEGs in the other species and test
####################################################################################################################

# Load the countsmatrix of the other (the transcriptome) species (=T. longispinosus)
load("countsMatrix_Tlongi.Rdata")

# Extract the rows of the countsmatrix that match with the orthologues of the DEGs of T. americanus in T. longispinosus
orthos_deg_tl_fo<-all[na.omit(match(deg_follower_tam$ID, all$Y)),]
orthos_deg_tl_le<-all[na.omit(match(deg_leader_tam$ID, all$Y)),]
orthos_deg_tl_sc<-all[na.omit(match(deg_scout_tam$ID, all$Y)),]
counts_deg_tl_fo<-countsMatrix[match(orthos_deg_tl_fo$V3, rownames(countsMatrix)),]
counts_deg_tl_le<-countsMatrix[match(orthos_deg_tl_le$V3, rownames(countsMatrix)),]
counts_deg_tl_sc<-countsMatrix[match(orthos_deg_tl_sc$V3, rownames(countsMatrix)),]

# Assign the right groups to the colums of this subset matrix
sample_info<-data.frame(group=c("Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout"))

# Orthologues of DEGs up-regulated in followers of T. americanus
# Remove rows with rowsum = 0, because they cause error and are not of interest
counts_deg_tl_fo<-counts_deg_tl_fo[rowSums(counts_deg_tl_fo)>0,]

# Create an empty data frame for the test results
results_deg_tl_fo<-data.frame(gene=c(rep(3,length(rownames(counts_deg_tl_fo)))), pvalue=c(rep(3,length(rownames(counts_deg_tl_fo)))), p_fo_le=c(rep(3,length(rownames(counts_deg_tl_fo)))), p_fo_sc=c(rep(3,length(rownames(counts_deg_tl_fo)))), p_le_sc=c(rep(3,length(rownames(counts_deg_tl_fo)))))

# Iterate over the contigs and test expression according to group
for (i in 1:length(rownames(counts_deg_tl_fo))){
  results_deg_tl_fo$gene[i]<-as.character(rownames(counts_deg_tl_fo)[i])
  x<-cbind(sample_info,as.numeric(t(counts_deg_tl_fo[i,])))
  colnames(x)<-c("group", "expression")
  y<-kruskal.test(expression~group, data=x)
  results_deg_tl_fo$pvalue[i]<-y$p.value
  
  PT = dunnTest(expression ~ group,data=x,method="bh")
  results_deg_tl_fo$p_fo_le[i]<-PT$res$P.adj[1]
  results_deg_tl_fo$p_fo_sc[i]<-PT$res$P.adj[2]
  results_deg_tl_fo$p_le_sc[i]<-PT$res$P.adj[3]
  
  z<-ggboxplot(x, x = "group", y = "expression", 
               fill = "group", palette = c("paleturquoise3","palevioletred2", "steelblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = "Number of raw reads", title=rownames(counts_deg_tl_le)[i], xlab = FALSE, font.tickslab=c(14,"plain", "black"))+stat_compare_means()+guides(fill=FALSE)+theme(plot.margin = unit(c(1.2,1.2,1.2,0.8), "lines"), plot.title = element_text(hjust = 0.5))
  
  ggexport(z,filename = paste("Boxplot_expression_Ortho_of_Tam_DEG_in_Tlongi_Follower",rownames(counts_deg_tl_fo)[i],".png", sep=""))
}

write.table(results_deg_tl_fo, file="Results_KruskalWallis_expression_Ortho_of_Tam_DEG_in_Tlongi_Follower.txt", sep="\\t", row.names = FALSE, quote = FALSE)

# Orthologues of DEGs up-regulated in leaders of T. americanus
# Remove rows with rowsum = 0, because they cause error and are not of interest
counts_deg_tl_le<-counts_deg_tl_le[rowSums(counts_deg_tl_le)>0,]

# Create an empty data frame for the test results
results_deg_tl_le<-data.frame(gene=c(rep(3,length(rownames(counts_deg_tl_le)))), pvalue=c(rep(3,length(rownames(counts_deg_tl_le)))), p_fo_le=c(rep(3,length(rownames(counts_deg_tl_le)))), p_fo_sc=c(rep(3,length(rownames(counts_deg_tl_le)))), p_le_sc=c(rep(3,length(rownames(counts_deg_tl_le)))))

# Iterate over the contigs and test expression according to group
for (i in 1:length(rownames(counts_deg_tl_le))){
  results_deg_tl_le$gene[i]<-as.character(rownames(counts_deg_tl_le)[i])
  x<-cbind(sample_info,as.numeric(t(counts_deg_tl_le[i,])))
  colnames(x)<-c("group", "expression")
  y<-kruskal.test(expression~group, data=x)
  results_deg_tl_le$pvalue[i]<-y$p.value
  
  PT = dunnTest(expression ~ group,data=x,method="bh")
  results_deg_tl_le$p_fo_le[i]<-PT$res$P.adj[1]
  results_deg_tl_le$p_fo_sc[i]<-PT$res$P.adj[2]
  results_deg_tl_le$p_le_sc[i]<-PT$res$P.adj[3]
  
  z<-ggboxplot(x, x = "group", y = "expression", 
               fill = "group", palette = c("paleturquoise3","palevioletred2", "steelblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = "Number of raw reads", title=rownames(counts_deg_tl_le)[i], xlab = FALSE, font.tickslab=c(14,"plain", "black"))+stat_compare_means()+guides(fill=FALSE)+theme(plot.margin = unit(c(1.2,1.2,1.2,0.8), "lines"), plot.title = element_text(hjust = 0.5))
  
  ggexport(z,filename = paste("Boxplot_expression_Ortho_of_Tam_DEG_in_Tlongi_Leader",rownames(counts_deg_tl_le)[i],".png", sep=""))
}

write.table(results_deg_tl_le, file="Results_KruskalWallis_expression_Ortho_of_Tam_DEG_in_Tlongi_Leader.txt", sep="\\t", row.names = FALSE, quote = FALSE)

# Orthologues of DEGs up-regulated in scouts of T. americanus
# Remove rows with rowsum = 0, because they cause error and are not of interest
counts_deg_tl_sc<-counts_deg_tl_sc[rowSums(counts_deg_tl_sc)>0,]

# Create an empty data frame for the test results
results_deg_tl_sc<-data.frame(gene=c(rep(3,length(rownames(counts_deg_tl_sc)))), pvalue=c(rep(3,length(rownames(counts_deg_tl_sc)))), p_fo_le=c(rep(3,length(rownames(counts_deg_tl_sc)))), p_fo_sc=c(rep(3,length(rownames(counts_deg_tl_sc)))), p_le_sc=c(rep(3,length(rownames(counts_deg_tl_sc)))))
# Iterate over the contigs and test expression according to group
for (i in 1:length(rownames(counts_deg_tl_sc))){
  results_deg_tl_sc$gene[i]<-as.character(rownames(counts_deg_tl_sc)[i])
  x<-cbind(sample_info,as.numeric(t(counts_deg_tl_sc[i,])))
  colnames(x)<-c("group", "expression")
  y<-kruskal.test(expression~group, data=x)
  results_deg_tl_sc$pvalue[i]<-y$p.value
  
  PT = dunnTest(expression ~ group,data=x,method="bh")
  results_deg_tl_sc$p_fo_le[i]<-PT$res$P.adj[1]
  results_deg_tl_sc$p_fo_sc[i]<-PT$res$P.adj[2]
  results_deg_tl_sc$p_le_sc[i]<-PT$res$P.adj[3]
  
  z<-ggboxplot(x, x = "group", y = "expression", 
               fill = "group", palette = c("paleturquoise3","palevioletred2", "steelblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = "Number of raw reads", title=rownames(counts_deg_tl_sc)[i], xlab = FALSE, font.tickslab=c(14,"plain", "black"))+stat_compare_means()+guides(fill=FALSE)+theme(plot.margin = unit(c(1.2,1.2,1.2,0.8), "lines"), plot.title = element_text(hjust = 0.5))
  
  ggexport(z,filename = paste("Boxplot_expression_Ortho_of_Tam_DEG_in_Tlongi_Scout",rownames(counts_deg_tl_sc)[i],".png", sep=""))
}

write.table(results_deg_tl_sc, file="Results_KruskalWallis_expression_Ortho_of_Tam_DEG_in_Tlongi_Scout.txt", sep="\\t", row.names = FALSE, quote = FALSE)


#####################################################################################################################
################################################### Part III ########################################################
###################### Orthogroups between DEGs of T. longispinosus and T. americanus transcriptome #################
#####################################################################################################################

#####################################################################################################################
## Extraction of the single-copy orthogroups found by Orthofinder
#####################################################################################################################

# Specify working directory for project
setwd("/Directory/Orthofinder/DEG_Tlongi_vs_Tam_Transcriptome/")

# Read in file with the number of contigs per species per orthogroup
geneCount<-read.table("Orthogroups.GeneCount.csv", header=TRUE, sep="\\t", quote="")

# Take only the ones where both species has only 1 contig present
geneCount_sc<-subset(geneCount, geneCount$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1 & geneCount$degs_sequences_tlongi==1)

contigsSC<-as.data.frame(geneCount_sc$X)

# Read in file with contigs assigned to each orthogroup
orthogroups<-read.table("Orthogroups.csv", sep="\\t", quote="", header=TRUE)

# Take only the orthogroups (and according contigs) that we extracted before
SC_Orthogroups<-orthogroups[match(contigsSC$`geneCount_sc$X`, orthogroups$X),]
colnames(SC_Orthogroups)<-c("X", "Y", "Z")

#####################################################################################################################
## Extraction of all the remaining orthogroups
#####################################################################################################################

# Take only the ones where at least one species has only 1 contig present
geneCount_rest<-subset(geneCount, xor(geneCount$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder==1, geneCount$degs_sequences_tlongi==1))

contigs_rest<-as.data.frame(geneCount_rest$X)

# Now create a data structure with all single copy orthologues in first column
geneCount_sc_tlongi<-subset(geneCount, geneCount$degs_sequences_tlongi ==1&geneCount$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder>1)

contigs_sc_tlongi<-as.data.frame(geneCount_sc_tlongi$X)

# Take only the orthogroups (and according contigs) that we extracted before (have to cheat a little to combine tables ;))
Orthogroups_sc_tlongi<-orthogroups[match(contigs_sc_tlongi$`geneCount_sc_tlongi$X`, orthogroups$X),]
Orthogroups_sc_tlongi_invert<-cbind(as.character(Orthogroups_sc_tlongi$X), as.character(Orthogroups_sc_tlongi$degs_sequences_tlongi), as.character(Orthogroups_sc_tlongi$Tamer_TrinityOnly.ShortContigNames_Transcriptome.fasta.transdecoder))

#####################################################################################################################
## Creation of new single-copy orthogroups from the orthogroups extracted before
#####################################################################################################################

# Read in information which contig is assigned to which ID (important for BLAST output later)
# Before opening delete uninteresting information (we just need ID and contig name) and replace ": " by ":"
contigToID<-read.table("SequenceIDs_edit.txt", header=FALSE, sep=":", quote="")

# Read in BLAST files
difftam_tl<-read.table("Blast0_1.txt", sep="\\t", quote="")
difftl_tam<-read.table("Blast1_0.txt", sep="\\t", quote="")

# Only take the necessary information: query, hit and score
difftam_tlx<-cbind(difftam_tl[1], difftam_tl[2], difftam_tl[12])
colnames(difftam_tlx)<-c("query", "hit", "score")
difftl_tamx<-cbind(difftl_tam[1], difftl_tam[2], difftl_tam[12])
colnames(difftl_tamx)<-c("query", "hit", "score")

# Combine both BLAST outputs
blast_2<-rbind(difftam_tlx, difftl_tamx)

# Create empty dataframe for new orthogroups
new_orthogroups<-Orthogroups_sc_tlongi_invert[,1:3]
new_orthogroups[,3]<-c(rep(1:length(new_orthogroups[,3])))

Orthogroups_sc_tlongi_invert<-as.data.frame(Orthogroups_sc_tlongi_invert)

# Now assign to each single copy species entry the entry from the other species with the highest blast score
for (i in 1:length(Orthogroups_sc_tlongi_invert[,2]) ){
  # Get ID of the contig of interest
  blast<-as.character(contigToID$V1[match(as.character(Orthogroups_sc_tlongi_invert$V2[i]),as.character(contigToID$V2))])
  
  # Get all the entries in the blast list whose query matches with the single copy contig
  blast_list<-blast_2[which(!is.na(match(as.character(blast_2$query),blast))),]
  
  # Order entries by hit score from highest to lowest
  blast_ordered<-blast_list[order(blast_list$score),]
  
  for(j in 1:length(blast_ordered$query)){
    
    blast_hit<-contigToID$V2[match(as.character(blast_ordered$hit[j]),as.character(contigToID$V1))]
    if(as.character(blast_hit) %in% strsplit(t(as.matrix(as.data.frame(Orthogroups_sc_tlongi_invert[i,]$V3))), split=", ")[[1]]){
      
      best_ortho<-strsplit(t(as.matrix(Orthogroups_sc_tlongi_invert[i,]$V3)),split=", ")[[1]][match(as.character(blast_hit), strsplit(t(as.matrix(Orthogroups_sc_tlongi_invert[i,]$V3)),split=", ")[[1]])]
      new_orthogroups[i,3]<-as.character(best_ortho)
    }
  }
}

colnames(new_orthogroups)<-c("X", "Y", "Z")

# Combine orthogroups
all<-rbind(SC_Orthogroups, new_orthogroups)

####################################################################################################################
## Check for matches with differentially expressed contigs in other species
####################################################################################################################

# Read in the files containing the all DEGs in T. longispinosus
deg_follower_tl<-read.table("DESeq2_All_Follower_Tlongi.csv", header = T)
deg_leader_tl<-read.table("DESeq2_All_Leader_Tlongi.csv", header = T)
deg_scout_tl<-read.table("DESeq2_All_Scout_Tlongi.csv", header=T)

# Read in the files containing the all DEGs of T. americanus
deg_follower_tam<-read.table("DESeq2_All_Follower_Tam.csv", header = T)
deg_leader_tam<-read.table("DESeq2_All_Leader_Tam.csv", header = T)
deg_scout_tam<-read.table("DESeq2_All_Scout_Tam.csv", header=T)

x<-which(match(deg_follower_tlongi$ID, all$Y)!="NA")
orthos_follower_tlongi<-all[x,]
y<-which(match(deg_leader_tlongi$ID, all$Y)!="NA")
orthos_leader_tlongi<-all[y,]
z<-which(match(deg_scout_tlongi$ID, all$Y)!="NA")
orthos_scout_tlongi<-all[z,]

# Plot Venn diagrams depicting overlap between DEGs in T. americanus and orthologues of DEGs in T. longispinosus
png(file="Follower_Tlongi_Orthologues_vs_DEGs_Tam.png", width = 820, height = 830)
venn.plot<-venn.diagram(list(deg_follower_tl$ID, orthos_follower_tam$V3), NULL, lwd=c(2,2), col=c("black","black"),category.names=c("DEGs TL","Orthos Tam"),fill=c("midnightblue", "mediumvioletred"), alpha=c(0.6, 0.6), cex=3, cat.fontface=4, cat.cex=c(1,1))
grid.draw(venn.plot)
dev.off()

png(file="Leader_Tlongi_Orthologues_vs_DEGs_Tam.png", width = 820, height = 830)
venn.plot<-venn.diagram(list(deg_leader_tl$ID, orthos_leader_tam$V3), NULL, lwd=c(2,2), col=c("black","black"),category.names=c("DEGs TL","Orthos Tam"),fill=c("midnightblue", "mediumvioletred"), alpha=c(0.6, 0.6), cex=3, cat.fontface=4, cat.cex=c(1,1))
grid.draw(venn.plot)
dev.off()

png(file="Scout_Tlongi_Orthologues_vs_DEGs_Tam.png", width = 820, height = 830)
venn.plot<-venn.diagram(list(deg_scout_tl$ID, orthos_scout_tam$V3), NULL, lwd=c(2,2), col=c("black","black"),category.names=c("DEGs TL","Orthos Tam"),fill=c("midnightblue", "mediumvioletred"), alpha=c(0.6, 0.6), cex=3, cat.fontface=4, cat.cex=c(1,1))
grid.draw(venn.plot)
dev.off()

####################################################################################################################
## Plot expression of the orthologues of the DEGs in the other species and test
####################################################################################################################

# Load the countsmatrix of the other (the transcriptome) species (=T. americanus)
load("countsMatrix_Tam.Rdata")

# Extract the rows of the countsmatrix that match with the orthologues of the DEGs of T. longispinosus in T. americanus
orthos_deg_tam_fo<-all[na.omit(match(deg_follower_tl$ID, all$Y)),]
orthos_deg_tam_le<-all[na.omit(match(deg_leader_tl$ID, all$Y)),]
orthos_deg_tam_sc<-all[na.omit(match(deg_scout_tl$ID, all$Y)),]
counts_deg_tam_fo<-countsMatrix[match(orthos_deg_tam_fo$V3, rownames(countsMatrix)),]
counts_deg_tam_le<-countsMatrix[match(orthos_deg_tam_le$V3, rownames(countsMatrix)),]
counts_deg_tam_sc<-countsMatrix[match(orthos_deg_tam_sc$V3, rownames(countsMatrix)),]

# Assign the right groups to the colums of this subset matrix
sample_info<-data.frame(group=c("Follower", "Leader", "Scout","Scout", "Leader", "Follower","Leader", "Scout", "Follower","Leader", "Scout", "Follower","Leader", "Follower", "Scout"))

# Orthologues of DEGs up-regulated in followers of T. americanus
# Remove rows with rowsum = 0, because they cause error and are not of interest
counts_deg_tam_fo<-counts_deg_tam_fo[rowSums(counts_deg_tam_fo)>0,]

# Create an empty data frame for the test results
results_deg_tam_fo<-data.frame(gene=c(rep(3,length(rownames(counts_deg_tam_fo)))), pvalue=c(rep(3,length(rownames(counts_deg_tam_fo)))), p_fo_le=c(rep(3,length(rownames(counts_deg_tam_fo)))), p_fo_sc=c(rep(3,length(rownames(counts_deg_tam_fo)))), p_le_sc=c(rep(3,length(rownames(counts_deg_tam_fo)))))

# Iterate over the contigs and test expression according to group
for (i in 1:length(rownames(counts_deg_tam_fo))){
  results_deg_tam_fo$gene[i]<-as.character(rownames(counts_deg_tam_fo)[i])
  x<-cbind(sample_info,as.numeric(t(counts_deg_tam_fo[i,])))
  colnames(x)<-c("group", "expression")
  y<-kruskal.test(expression~group, data=x)
  results_deg_tam_fo$pvalue[i]<-y$p.value
  
  PT = dunnTest(expression ~ group,data=x,method="bh")
  results_deg_tam_fo$p_fo_le[i]<-PT$res$P.adj[1]
  results_deg_tam_fo$p_fo_sc[i]<-PT$res$P.adj[2]
  results_deg_tam_fo$p_le_sc[i]<-PT$res$P.adj[3]
  
  z<-ggboxplot(x, x = "group", y = "expression", 
               fill = "group", palette = c("mediumturquoise","mediumvioletred", "midnightblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = "Number of raw reads", title=rownames(counts_deg_tam_le)[i], xlab = FALSE, font.tickslab=c(14,"plain", "black"))+stat_compare_means()+guides(fill=FALSE)+theme(plot.margin = unit(c(1.2,1.2,1.2,0.8), "lines"), plot.title = element_text(hjust = 0.5))
  
  ggexport(z,filename = paste("Boxplot_expression_Ortho_of_Tlongi_DEG_in_Tam_Follower",rownames(counts_deg_tl_fo)[i],".png", sep=""))
}

write.table(results_deg_tam_fo, file="Results_KruskalWallis_expression_Ortho_of_Tam_DEG_in_Tlongi_Follower.txt", sep="\\t", row.names = FALSE, quote = FALSE)

# Orthologues of DEGs up-regulated in leaders of T. longispinosus
# Remove rows with rowsum = 0, because they cause error and are not of interest
counts_deg_tam_le<-counts_deg_tam_le[rowSums(counts_deg_tam_le)>0,]

# Create an empty data frame for the test results
results_deg_tam_le<-data.frame(gene=c(rep(3,length(rownames(counts_deg_tam_le)))), pvalue=c(rep(3,length(rownames(counts_deg_tam_le)))), p_fo_le=c(rep(3,length(rownames(counts_deg_tam_le)))), p_fo_sc=c(rep(3,length(rownames(counts_deg_tam_le)))), p_le_sc=c(rep(3,length(rownames(counts_deg_tam_le)))))

# Iterate over the contigs and test expression according to group
for (i in 1:length(rownames(counts_deg_tam_le))){
  results_deg_tam_le$gene[i]<-as.character(rownames(counts_deg_tam_le)[i])
  x<-cbind(sample_info,as.numeric(t(counts_deg_tam_le[i,])))
  colnames(x)<-c("group", "expression")
  y<-kruskal.test(expression~group, data=x)
  results_deg_tam_le$pvalue[i]<-y$p.value
  
  PT = dunnTest(expression ~ group,data=x,method="bh")
  results_deg_tam_le$p_fo_le[i]<-PT$res$P.adj[1]
  results_deg_tam_le$p_fo_sc[i]<-PT$res$P.adj[2]
  results_deg_tam_le$p_le_sc[i]<-PT$res$P.adj[3]
  
  z<-ggboxplot(x, x = "group", y = "expression", 
               fill = "group", palette = c("mediumturquoise","mediumvioletred", "midnightblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = "Number of raw reads", title=rownames(counts_deg_tam_le)[i], xlab = FALSE, font.tickslab=c(14,"plain", "black"))+stat_compare_means()+guides(fill=FALSE)+theme(plot.margin = unit(c(1.2,1.2,1.2,0.8), "lines"), plot.title = element_text(hjust = 0.5))
  
  ggexport(z,filename = paste("Boxplot_expression_Ortho_of_Tlongi_DEG_in_Tam_Leader",rownames(counts_deg_tam_le)[i],".png", sep=""))
}

write.table(results_deg_tam_le, file="Results_KruskalWallis_expression_Ortho_of_Tlongi_DEG_in_Tam_Leader.txt", sep="\\t", row.names = FALSE, quote = FALSE)

# Orthologues of DEGs up-regulated in scouts of T. longispinosus
# Remove rows with rowsum = 0, because they cause error and are not of interest
counts_deg_tam_sc<-counts_deg_tam_sc[rowSums(counts_deg_tam_sc)>0,]

# Create an empty data frame for the test results
results_deg_tam_sc<-data.frame(gene=c(rep(3,length(rownames(counts_deg_tam_sc)))), pvalue=c(rep(3,length(rownames(counts_deg_tam_sc)))), p_fo_le=c(rep(3,length(rownames(counts_deg_tam_sc)))), p_fo_sc=c(rep(3,length(rownames(counts_deg_tam_sc)))), p_le_sc=c(rep(3,length(rownames(counts_deg_tam_sc)))))
# Iterate over the contigs and test expression according to group
for (i in 1:length(rownames(counts_deg_tam_sc))){
  results_deg_tam_sc$gene[i]<-as.character(rownames(counts_deg_tam_sc)[i])
  x<-cbind(sample_info,as.numeric(t(counts_deg_tam_sc[i,])))
  colnames(x)<-c("group", "expression")
  y<-kruskal.test(expression~group, data=x)
  results_deg_tam_sc$pvalue[i]<-y$p.value
  
  PT = dunnTest(expression ~ group,data=x,method="bh")
  results_deg_tam_sc$p_fo_le[i]<-PT$res$P.adj[1]
  results_deg_tam_sc$p_fo_sc[i]<-PT$res$P.adj[2]
  results_deg_tam_sc$p_le_sc[i]<-PT$res$P.adj[3]
  
  z<-ggboxplot(x, x = "group", y = "expression", 
               fill = "group", palette = c("mediumturquoise","mediumvioletred", "midnightblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = "Number of raw reads", title=rownames(counts_deg_tam_sc)[i], xlab = FALSE, font.tickslab=c(14,"plain", "black"))+stat_compare_means()+guides(fill=FALSE)+theme(plot.margin = unit(c(1.2,1.2,1.2,0.8), "lines"), plot.title = element_text(hjust = 0.5))
  
  ggexport(z,filename = paste("Boxplot_expression_Ortho_of_Tlongi_DEG_in_Tam_Scout",rownames(counts_deg_tam_sc)[i],".png", sep=""))
}

write.table(results_deg_tam_sc, file="Results_KruskalWallis_expression_Ortho_of_Tlongi_DEG_in_Tam_Scout.txt", sep="\\t", row.names = FALSE, quote = FALSE)
