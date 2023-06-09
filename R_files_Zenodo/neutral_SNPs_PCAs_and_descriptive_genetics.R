###Exploration of neutral SNPs in South Fork Eel River and its tributaries
##Contact Suzanne Kelson, skelson@berkeley.edu, for questions
##This code was written and run on a PC

rm(list = ls())
library(adegenet)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(hierfstat)
library(plyr);library(dplyr)


setwd("")
fish <- read.csv("platemaps_omy5genotypes.csv")
fish <- subset(fish, reason_included != "blank") ##4517 samples TOTAL 

###Single read SNPs (IBS matrix) where SNPs that have a maximum of 20% missing per individual are included###
single_reads_colsnps <- read.csv("single_read_snps_maxmis20p.csv")

##remove the blanks from the dataframe
not_samples <- c("pB04_wE11", "pB04_wE09", "pB04_wC02", "pB04_wG03", "pB04_wE08", "pB04_wF05", "pB04_wC04",
                 "pB04_wH08", "pB04_wE02","pB04_wG09" ,"pB04_wH07" ,"pB04_wF03" ,"pB04_wC11", "pB04_wF09", 
                 "pB04_wH06","pB04_wC09","pB04_wD02", "pB04_wE12" ,"pB04_wE03" ,"pB04_wD10", "pB04_wE05", 
                 "pB04_wC10", "pB04_wC05" ,"pB04_wD04","pB04_wF08", "pB04_wD09", "pB04_wD11", "pB04_wC03",
                 "pB04_wC08" ,"pB04_wE06", "pB04_wE01", "pB04_wD05", "pB04_wF01", "pB04_wD03", "pB04_wD12",
                 "pB04_wF02" ,"pB04_wC12", "pB04_wD07", "pB04_wC07", "pB04_wD08", "pB04_wD06", "pB04_wC01",
                 "pB04_wC06", "pB04_wD01","NA_20E12", "NA_20F06", "NA_4H12", "NA_6F06","NA_8E03","NA_8G11",
                 "NA_8G12","NA_9H07","pA02_wE02","pA09_wF01", "pD02_wD09",   "pA01_wA01.1")
single_reads_colsnps <- single_reads_colsnps[!(single_reads_colsnps$sample_ID %in% not_samples),]
single_reads_colsnps <-  merge(select(fish,sample_ID,FID,year,date,location,reason_included,pool,FL,wt,age_class), single_reads_colsnps, by = "sample_ID",all.y = T)

########
##genind for analyses in adegenet WITHOUT omy5 SNPs
#####

snp_list <- data.frame(colnames(single_reads_colsnps)[11:dim(single_reads_colsnps)[2]])
colnames(snp_list)[1]<- "chr_pos"
snp_list$chrpos <- snp_list$chr_pos
snp_list$chrpos<- gsub("_",".",snp_list$chrpos)
snp_list$chr <- unlist(strsplit(snp_list$chrpos,split="[.]"))[seq(from=1,to=726*2,by=2)]
snp_list$chr_num <- as.numeric(as.factor(snp_list$chr))
omy5_snps <- droplevels(snp_list$chr_pos[snp_list$chr_num==5]) ##34 SNPs

##remove SNPs on Omy5
single_reads_colsnps_woomy5 <- single_reads_colsnps[,!(colnames(single_reads_colsnps)%in%omy5_snps)]

genind_wo5 <- df2genind(single_reads_colsnps_woomy5[,c(11:dim(single_reads_colsnps_woomy5)[2])], ploidy=2, sep="",
                        loc.names = colnames(single_reads_colsnps_woomy5)[11:dim(single_reads_colsnps_woomy5)[2]], 
                        ind.names = single_reads_colsnps_woomy5$sample_ID,
                        pop = single_reads_colsnps_woomy5$location,
                        NA.char = -1)
genind_wo5_scaled <- scaleGen(genind_wo5, NA.method="mean")

pca_wo5 <- dudi.pca(genind_wo5_scaled,scannf = FALSE, nf = 300)#select 300 axes

##merge dataframe to pull in covariates
pcatab <- data.frame(pca_wo5$li)
pcatab$sample_ID <- rownames(pcatab)
tab_woomy5 <- merge(pcatab[,c(1:10, 301)], fish, by = "sample_ID")

##relabel and re order location variable
colnames(tab_woomy5)[14]<- "Location"
levels(tab_woomy5$Location)[c(1,3)]<- c("Elder Above", "Elder Below")
tab_woomy5$Location<-factor(tab_woomy5$Location,levels=levels(tab_woomy5$Location)[c(4,3,1,8,9,10,2,5,6,7)])

##eigenvalue bar plot
eig_tab_wo5 <- data.frame(PC = seq(1:10), Eigenvalue = pca_wo5$eig[1:10],Variance = pca_wo5$eig[1:10]/sum(pca_wo5$eig)*100)
bp<- ggplot(data=eig_tab_wo5,aes(x=as.factor(PC),y=Variance))+geom_bar(stat="identity")+
  labs(x = "Principal Components", y = "% Variance")+theme_classic(10)
bp

ggsave("figS2_barplot.jpg", plot = bp,
       scale = 1, width = 2.5, height = 1.8, units = "in",
       dpi = 500 )

within_circle <- subset(tab_woomy5, Axis1>-4 & Axis2 <4)
3607/4505

##loading plots
loadings_wo5 <- data.frame(pca_wo5$c1)
loadings_wo5$ChrPos <- as.factor(rownames(loadings_wo5))
loadings_wo5$Chr <- as.factor(substr(loadings_wo5$ChrPos,1,5))
##lOading plot with colored bars, no threshold
par(mfrow = c(1,1), mar = c(2,4,3,2),oma = c(1,1,1,1))
loadingplot(loadings_wo5[,1:300]^2, col = loadings_wo5$Chr, main = "Loadings for Principal Components",lab=NA,
            xlab = "SNPs", ylab = "Loadings",threshold = F)


##scatter plot by location - Fig S2
scatterpl <- ggplot(data=subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"),aes(x=Axis1,y=Axis2,color=Location))+
  geom_point(alpha=3/5,size = 0.75)+
  scale_color_brewer(palette ="Set1")+labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")+
  theme_classic(10)
scatterpl

location_meanse <- ddply(subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"), 
                         .(Location),summarize,
                          mean_axis1 = mean(Axis1,na.rm = T),
                          sd_axis1 = sd(Axis1, na.rm = T),
                          mean_axis2 = mean(Axis2,na.rm = T),
                          sd_axis2 = sd(Axis2,na.rm = T))
#location_colors <- RColorBrewer::brewer.pal(6, "Set1")

scatterpl_location_wmeanse <- ggplot()+
  geom_point(data=subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"),
             aes(x = Axis1, y = Axis2, color = Location),alpha=3/5, size = 0.75)+
  geom_errorbar(data=location_meanse, width=0.4,
                aes(group = Location,x=mean_axis1,
                    ymin = (mean_axis2 - sd_axis2), ymax = (mean_axis2 + sd_axis2)))+
  geom_errorbarh(data=location_meanse,
                 aes(group = Location, xmin = (mean_axis1-sd_axis1), xmax = (mean_axis1+sd_axis1),
                     y=mean_axis2),height=0.4)+
  geom_point(data = location_meanse, aes(x = mean_axis1,y=mean_axis2, fill = Location),
             size=3,color="black", shape = 21)+
  theme_classic(10)+
  labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")+
  scale_color_brewer(palette ="Set1")+
  scale_fill_brewer(palette ="Set1")+
  theme(legend.position = "top",legend.key.size = unit(0.2, "cm"),legend.text=element_text(size=8))
scatterpl_location_wmeanse


 ##plotting by omy5 genotype
tab_woomy5$genotype[tab_woomy5$missing_data>250]<- NA  ##remove genotypes for fish that were missing data at over 250 SNPs
tab_woomy5$genotype<- factor(tab_woomy5$genotype, levels = levels(tab_woomy5$genotype)[c(2,1,3)])

scatterpl_genotypes <- ggplot(subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"), aes(x = Axis1, y = Axis2, color = genotype))+
  geom_point(alpha=3/5, size = 0.75)+
  theme_classic(10)+
  labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")+
  scale_color_manual(values = c("firebrick","purple","darkblue"), name = "Genotype",na.translate=F)
scatterpl_genotypes

##add means and standard deviation
genotypes_meanse <- ddply(subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"),
                          .(genotype),summarize,
                          mean_axis1 = mean(Axis1,na.rm = T),sd_axis1 = sd(Axis1, na.rm = T),
                          mean_axis2 = mean(Axis2,na.rm = T), sd_axis2 = sd(Axis2,na.rm = T))


scatterpl_genotypes_wmeanse <- ggplot()+
  geom_point(data=subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"), 
             aes(x = Axis1, y = Axis2, color = genotype),alpha=3/5, size = 0.75)+
  geom_errorbar(data=genotypes_meanse, width=0.4,
                aes(group = genotype,x=mean_axis1,
                    ymin = (mean_axis2 - sd_axis2), ymax = (mean_axis2 + sd_axis2)))+
  geom_errorbarh(data=genotypes_meanse,
                 aes(group = genotype, xmin = (mean_axis1-sd_axis1), xmax = (mean_axis1+sd_axis1),
                     y=mean_axis2),height=0.4)+
  geom_point(data = genotypes_meanse, aes(x = mean_axis1,y=mean_axis2, fill = genotype),
             size=3,color="black", shape = 21)+
  theme_classic(10)+
  labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")+
  scale_color_manual(values = c("firebrick","purple","darkblue"), name = "Genotype",na.translate=F)+
  scale_fill_manual(values = c("firebrick","purple","darkblue"), name = "Genotype",na.translate=F)+
  theme(legend.position = "top",legend.key.size = unit(0.2, "cm"),legend.text=element_text(size=8))
scatterpl_genotypes_wmeanse


##scatter plot by year
scatterpl_year <- ggplot(subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"),

                                                  aes(x = Axis1, y = Axis2, color = as.factor(year)))+
  geom_point(alpha=3/5,size=0.75)+
  theme_classic(10)+
  scale_color_manual(values = c("wheat3","darkred","skyblue3","darkblue"),na.translate = F,name = "Year")+
  labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")
scatterpl_year

year_meanse <- ddply(subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"),
                          .(year),summarize,
                          mean_axis1 = mean(Axis1,na.rm = T),sd_axis1 = sd(Axis1, na.rm = T),
                          mean_axis2 = mean(Axis2,na.rm = T),sd_axis2 = sd(Axis2,na.rm = T))

scatterpl_year_wmeanse <- ggplot()+
  geom_point(data=subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"), 
             aes(x = Axis1, y = Axis2, color = as.factor(year)),alpha=3/5, size = 0.75)+
  geom_errorbar(data=year_meanse, width=0.4,
                aes(group = year,x=mean_axis1,
                    ymin = (mean_axis2 - sd_axis2), ymax = (mean_axis2 + sd_axis2)))+
  geom_errorbarh(data=year_meanse,
                 aes(group = year, xmin = (mean_axis1-sd_axis1), xmax = (mean_axis1+sd_axis1),
                     y=mean_axis2),height=0.4)+
  geom_point(data = year_meanse, aes(x = mean_axis1,y=mean_axis2, fill = as.factor(year)),
             size=3,color="black", shape = 21)+
  theme_classic(10)+
  labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")+
  scale_color_manual(values = c("wheat3","darkred","skyblue3","darkblue"), name = "Year",na.translate=F)+
  scale_fill_manual(values = c("wheat3","darkred","skyblue3","darkblue"), name = "Year",na.translate=F)+
  theme(legend.position = "top",legend.key.size = unit(0.2, "cm"),legend.text=element_text(size=8))
scatterpl_year_wmeanse


###plot by missing data
missing_data_summary <- data.frame(sample_ID = single_reads_colsnps_woomy5$sample_ID, missing_data_allSNPs = NA)
for(i in 1:dim(single_reads_colsnps_woomy5)[1]){
  missing_data_summary$missing_data_allSNPs[i]<- sum(single_reads_colsnps_woomy5[i,-c(1:11)]==-1)
  single_reads_colsnps_woomy5}
##merge missing data frame to tab_woomy5 data frame
tab_woomy5<- merge(missing_data_summary, tab_woomy5, by = "sample_ID")

scatterpl_missingdata <- ggplot(subset(tab_woomy5, stream == "Elder"| stream == "Fox" | stream == "South Fork"), 
aes(x = Axis1, y = Axis2, color = missing_data_allSNPs))+
  geom_point(alpha=3/5,size=0.75)+
  theme_classic(10)+
  labs(x= "Axis 1 - 1.5%", y = "Axis 2 - 1.2%")+
  scale_color_continuous(name = "Num. SNPs\\nMissing Data")+
  theme(legend.position = "top",legend.text=element_text(size=8))
scatterpl_missingdata

###Group all scatter plots together (Fig S2)
allscatterplots <- ggarrange(scatterpl_location_wmeanse, scatterpl_genotypes_wmeanse,scatterpl_year_wmeanse, scatterpl_missingdata, nrow=2,ncol = 2)
allscatterplots
ggsave("figS2_scatterplots.pdf",
   plot = allscatterplots,scale = 1, width = 7, height =7.25, units = c("in"), dpi = 350 )

ggsave("figS2_scatterplots.jpg",
       plot = allscatterplots,scale = 1, width = 7, height =7.25, units = c("in"), dpi = 350 )


######################################################
##  Descriptive Genetics with called genotypes dataset: Fst Values and Heterozygosity ###
########################################################

put_last_dim_first <- function(df){df <- df[,c(dim(df)[2],1:(dim(df)[2]-1))]}
calledgenos <- read.csv("calledgenos_minInd50_95p.csv")
calledgenos$ChrPos <- paste(calledgenos$Chr, calledgenos$Pos, sep = ".")
calledgenos <- put_last_dim_first(calledgenos)

num_samples <- (dim(calledgenos)[2]-3)/3
##pulling out the AG column for adegenet
calledgenos <- calledgenos[c(1:3, seq(from=5,by=3,length.out=num_samples))]

##remove individuals that are missing data at all SNPs - calculate number of levels
nlevels <- data.frame(colnum = seq(from=4,to=4608, by = 1), nlevels = NA)
for(i in 1:4605){nlevels[i,2]<- nlevels(calledgenos[,(i+3)])}
levels_NN <- subset(nlevels, nlevels == 1)
##remove these columns
calledgenos <- calledgenos[,-c(levels_NN$colnum)];rm(levels_NN, nlevels)

snps <- calledgenos[,c(1:3)]
snps$ChrPos <- as.factor(snps$ChrPos)
snps$Chrnum <- as.numeric(snps$Chr)
snps$locnames <- paste(snps$Chrnum, snps$Pos, sep = "-")

##change to SNPs in column format
calledgenos_colsnps <- data.frame(t(calledgenos[,c(4:dim(calledgenos)[2])]))
colnames(calledgenos_colsnps)<- snps$locnames
calledgenos_colsnps$sample_ID <- colnames(calledgenos[4:dim(calledgenos)[2]])
calledgenos_colsnps$sample_ID<- as.factor(calledgenos_colsnps$sample_ID)
calledgenos_colsnps<- put_last_dim_first(calledgenos_colsnps)
for(i in 1:dim(calledgenos_colsnps)[2]){attr(calledgenos_colsnps[,i], "names")<- NULL} ##remove attributes

##add location and year data
calledgenos_colsnps <-  merge(select(fish,sample_ID,year,location,reason_included, pool, FL, wt), 
                                     calledgenos_colsnps, by = "sample_ID",all.y = T)
calledgenos_colsnps<- droplevels(calledgenos_colsnps[!(is.na(calledgenos_colsnps$location)),]) ##remove the blanks
##Should have 4441 individuals here 

##remove Mill Creek Barnwell and hatchery fish for analyses
calledgenos_colsnps<- calledgenos_colsnps[!(calledgenos_colsnps$location %in% c("Barnwell","hatchery","Mill Above","Mill Below")),]
calledgenos_colsnps<- droplevels(calledgenos_colsnps)

##SUBSTITUTE nn FOR na
calledgenos_colsnps <- data.frame(lapply(calledgenos_colsnps, function(x){gsub("NN", NA,x)}))
colnames(calledgenos_colsnps)[-c(1:7)]<- snps$locnames
###########
######## Genind object without Omy5 snps to calculate pariwise Fst between populations______####
##########
omy5snps_colnames <-  c("5-28579373","5-29641956","5-31476724","5-31476741","5-31675268","5-31675278","5-33918491","5-35157471", 
"5-45479654","5-45735028","5-47278098","5-49756324","5-49756334", "5-60143089","5-60143127","5-63158919","5-64034136,5-64867490","5-64867523")  
calledgenos_colsnps_woOmy5 <- calledgenos_colsnps[,!(colnames(calledgenos_colsnps)%in% omy5snps_colnames)]
#calledgenos_colsnps_woOmy5 <- calledgenos_colsnps[,c(1:78,98:477)]
#colnames(calledgenos_colsnps_woOmy5)
loci <- calledgenos_colsnps_woOmy5[,-c(1:7)]
loci <- data.frame(lapply(loci, function(x){gsub("NN", NA,x)}))
snps_woOmy5 <- droplevels(subset(snps, Chr != "omy05"))
colnames(loci)<- snps_woOmy5$locnames
sampleid <- calledgenos_colsnps_woOmy5$sample_ID
pop <- calledgenos_colsnps_woOmy5$location
locnames <- colnames(calledgenos_colsnps_woOmy5[,-c(1:7)])
##make the genind object
mygendata_woOmy05 <- df2genind(loci, ploidy=2, sep="",ind.names = sampleid, pop = pop,loc.names=locnames)

##calculate pairwise Fst using hierFst
pairwisefst <- pairwise.fst(mygendata_woOmy05) ##all pairwise Fst values are under 0.02 when rounded up
levels(mygendata_woOmy05$pop)

###########
######## Genind object with all snps to calculate observed vs expected heterozygosity______####
##########
loci_all <- calledgenos_colsnps[,-c(1:7)] ##remove first seven information columns
loci_all <- data.frame(lapply(loci_all, function(x){gsub("NN", NA,x)}))##substitute NN for NA
colnames(loci_all)<- snps$locnames
sampleid <- calledgenos_colsnps$sample_ID
pop <- calledgenos_colsnps$location
locnames <- colnames(calledgenos_colsnps[,-c(1:7)])
mygendata <- df2genind(loci_all, ploidy=2, sep="",ind.names = sampleid, pop = pop,loc.names=locnames)

##summarize data
summary <- summary(mygendata)
summary_df <- data.frame(Hobs = summary$Hobs, Hexp = summary$Hexp) ##create a dataframe with Hexp and Hobs
summary_df$locnames <- rownames(summary_df)
summary_df <- merge(summary_df, snps, by = "locnames") ##put chromosome number in the data frame
summary_df<- summary_df[order(summary_df$Chrnum),]
barplot(summary$Hobs-summary$Hexp, main="Heterozygosity: observed-expected",
        ylab="Hobs - Hexp")
t.test(summary$Hexp,summary$Hobs,pair=T,var.equal=TRUE, alt = "greater") ##Test if Expected Heterozygosity is greater than Observed Heterozygosity


cols1 <- rep(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00","#F0E442","#CC79A7"),4)
cols_bw <- c(rep(c("black","gray50"),15))

levels(summary_df$Chr)<-c("Omy01","Omy02","Omy03","Omy04","Omy05","Omy06","Omy07","Omy08","Omy09","Omy10","Omy11","Omy12",
                          "Omy13","Omy14","Omy15","Omy16","Omy17","Omy18","Omy19","Omy20","Omy21","Omy22","Omy23","Omy24",
                          "Omy25","Omy26","Omy27","Omy28","Omy29")
###with legend
Heterozygositypl <- ggplot(data = summary_df, aes(x = Chrnum, y = Hobs-Hexp))+geom_jitter(aes(color = Chr), size = 1)+
  theme_classic(10)+labs(x = "Chromosome")+
  geom_hline(aes(yintercept = 0), color = "gray50")+
  scale_color_manual(values = cols1, name = "Chromosome")+
  scale_x_continuous(breaks = c(1,5,10,15,20,25,29))
Heterozygositypl

##no legend
Heterozygositypl <- ggplot(data = summary_df, aes(x = Chrnum, y = Hobs-Hexp))+geom_jitter(aes(color = Chr), size = 1)+
  theme_classic(10)+labs(x = "Chromosome")+
  geom_hline(aes(yintercept = 0), color = "gray50")+
  scale_color_manual(values = cols1, name = "Chromosome")+
  scale_x_continuous(breaks = c(1:29))+guides(color = F)
Heterozygositypl

ggsave("figS1_wlegend.jpg", plot = Heterozygositypl,
   scale = 1, width = 6, height = 5, units = c("in"),
  dpi = 500 )


