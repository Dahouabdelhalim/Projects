#####Discriminant Analysis of Principal Components on single-read SNPs on Omy5
##This script assigns each individual a resident-migratory genotype using SNPs on Omy5
##This code was written and run on a PC
##Contact Suzanne Kelson, skelson@berkeley.edu, for questions

rm(list = ls())
library(adegenet)
library(reshape2)
library(dplyr)
R19198 <- "omy05.28579373"
put_last_dim_first <- function(df){df <- df[,c(dim(df)[2], 1:dim(df)[2]-1)]}

##SET WORKING DIRECTORY WHERE YOU HAVE SAVED FILES
setwd("")

#fish data
fish <- read.csv("platemaps_omy5genotypes.csv") ##Note that originally, we read this file in without the Omy5 genotypes
blanks <- subset(fish, reason_included == "blank")
fish <- subset(fish, reason_included != "blank") 

##Use IBS matrix (single read SNPs) with SNPs that are missing data at max 50% of individuals 
single_reads <- read.csv("single_read_snps_maxmis50p.csv")
#format data frame
single_reads_snps <- data.frame(colnames(single_reads[-1]))
colnames(single_reads_snps)[1]<- "ChrPos"
single_reads_snps$Chr_Pos<-single_reads_snps$ChrPos
single_reads_snps$ChrPos <- gsub("_", ".",single_reads_snps$ChrPos)
single_reads_snps$Chr <- unlist(strsplit(single_reads_snps$ChrPos,"[.]"))[seq(from=1,to=7427*2,by=2)]
single_reads_snps$Chr_num <- as.numeric(as.factor((single_reads_snps$Chr)))

##subset for Omy5
omy5_snps <- droplevels(single_reads_snps$Chr_Pos[single_reads_snps$Chr_num==5])
omy5cols <- colnames(single_reads)%in%omy5_snps
single_reads_omy5<- single_reads[,omy5cols==TRUE]
single_reads_omy5$sample_ID <-single_reads$sample_ID
single_reads_omy5<- put_last_dim_first(single_reads_omy5)
##remove blanks from singe reads data frame
single_reads_omy5 <-single_reads_omy5[!(single_reads_omy5$sample_ID %in% blanks$sample_ID),] ##remove blanks
single_reads_omy5<- subset(single_reads_omy5, sample_ID != "pA01_wA01.1") ##duplicated in genotype data frame)

##these are ghost individuals that are in the single read file
not_samples <- c("pB04_wE11", "pB04_wE09", "pB04_wC02", "pB04_wG03", "pB04_wE08", "pB04_wF05", "pB04_wC04","pB04_wH08",
                  "pB04_wE02","pB04_wG09" ,"pB04_wH07" ,"pB04_wF03" ,"pB04_wC11", "pB04_wF09", "pB04_wH06","pB04_wC09",
                  "pB04_wD02", "pB04_wE12","pB04_wE03" ,"pB04_wD10", "pB04_wE05", "pB04_wC10", "pB04_wC05" ,"pB04_wD04",
                  "pB04_wF08", "pB04_wD09", "pB04_wD11", "pB04_wC03","pB04_wC08" ,"pB04_wE06", "pB04_wE01", "pB04_wD05",
                  "pB04_wF01", "pB04_wD03", "pB04_wD12", "pB04_wF02" ,"pB04_wC12", "pB04_wD07","pB04_wC07", "pB04_wD08",
                  "pB04_wD06", "pB04_wC01","pB04_wC06", "pB04_wD01")
single_reads_omy5 <- single_reads_omy5[!(single_reads_omy5$sample_ID %in% not_samples),]



single_reads_omy5 <-  droplevels(merge(select(fish, year, FID, location, sample_ID, pool, FL, wt, PIT, age_class, recap), single_reads_omy5,
                            by = "sample_ID",all.y = T))


#create genind object for analyses in adegenet
genind <- df2genind(single_reads_omy5[,c(11:dim(single_reads_omy5)[2])], ploidy=2, sep="",
                    loc.names = colnames(single_reads_omy5)[11:dim(single_reads_omy5)[2]], 
                    ind.names = single_reads_omy5$sample_ID,
                    pop = single_reads_omy5$location,
                    NA.char = -1)

##include 300 PCS ##this next line is slow, wait for it to finish before moving on
group <- find.clusters(genind, max.n.clust = 10, n.pca = 300) ##choose 3 clusters
##this next line is slow, wait for it to finish before moving on
dapc <- dapc(genind, group$grp) ##used 300 PCs and 2 Discriminant functions

##make a data frame with groups and PCs
groups_df <- data.frame(group = group$grp)
groups_df$sample_ID <- as.factor(rownames(groups_df))
pc_tab <- data.frame(dapc$tab)
pc_tab <- pc_tab[,c(1:10)]
pc_tab$sample_ID <- as.factor(row.names(pc_tab))

groups_df <- merge(groups_df, pc_tab, by = "sample_ID")
rm(pc_tab)

plot(PCA.pc.2 ~ PCA.pc.1, col = group, data=groups_df)
groups_df <- merge(fish, groups_df, all.y = T, by = "sample_ID")

#scatter(dapc, scree.da = T)

##pca in adegenet
table(groups_df$group, groups_df$location)
#### IMPORTANT NOTE: Every time this code is run, the naming of each group (1,2,3) will change. 
#Code here represents the naming for the file that was used.
groups_df$genotype <- "Heterozygote"
groups_df$genotype[groups_df$group ==1]<- "Resident"
groups_df$genotype[groups_df$group ==3]<- "Migratory"
groups_df$genotype<- as.factor(groups_df$genotype)
#groups_df <- merge(groups_df, single_reads_colsnps_omy5[,c(1,41)], by = "sample_ID", all.x = T)
#3 is het, #1s migratory #2 is resident ##(This will change for every run)
##plot the DAPC
redcol = c( "orange", "orangered2","dark red")
scatter(dapc, scree.da = FALSE, solid = 0.5, pch = c(15,16,17), col = redcol, clab = 1, leg = T,
        txt.leg = c( "Resident","Heterozygote", "Migratory"))
loadingplot(dapc$var.contr)#

##create a missing data column by counting the number of -1s
single_reads_omy5$missing_data <- 0
for (i in 1:dim(single_reads_omy5)[1]){
  single_reads_omy5$missing_data[i] <- sum(single_reads_omy5[i,c(7:dim(single_reads_omy5)[2])] == -1)} 

groups_df <- merge(groups_df, single_reads_omy5[,c(1,422)], by = "sample_ID", all.x = T)

##Group probabilities
group_probabilities <- data.frame(dapc$posterior)
group_probabilities$sample_ID <- rownames(group_probabilities)
colnames(group_probabilities)[c(1:3)]<- c("post_prob_1", "post_prob_2","post_prob_3")

groups_df<- merge(groups_df, group_probabilities, by = "sample_ID")


##export file with DAPC groups, DAPC genotypes (same thing), number of missing SNPs, and R19198
#write.csv(groups_df, row.names = F)


