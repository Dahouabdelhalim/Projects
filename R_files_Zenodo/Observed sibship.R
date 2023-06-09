### This code reproduces the results reported for 'Experiment 5' in the manuscript: 
### Burgess SC, Powell J, Bueno M. Dispersal, kin aggregation, and the fitness consequences of not spreading sibling larvae.
# Code finalized Jun 2022
# Any comments or error reporting, please contact Scott Burgess. sburgess@bio.fsu.edu

# R version 4.0.5 (2021-03-31) -- "Shake and Throw"

# Import data
# set.wd() # set working directory, or file path, before importing data
dat <- read.csv("Sibships.csv",colClasses=c(rep('factor',2),rep('numeric',2)))
sampleIDfile <- read.csv("sampleIDfile_msats.csv",colClasses='factor')

summary(dat)

# # with(dat[dat$Run==1,],hist(Probability,breaks=seq(0,1,0.05)))
# with(dat[dat$Run==2,],hist(Probability,breaks=seq(0,1,0.05)))
# with(dat[dat$Run==3,],hist(Probability,breaks=seq(0,1,0.05)))

# Get full or half-sibs with highest confidence (Probability >= threshold in all runs)
threshold <- 0.9 # CHANGE THIS TO EITHER 0.8, 0.9, or 0.95 to replicate results in paper
dat_sub <- dat[dat$Probability>=threshold,]

dat_sub$dyad <- factor(interaction(dat_sub$OffspringID1,dat_sub$OffspringID2))
dat_sub_table <- with(dat_sub, table(dyad))
dat_sub_table <- data.frame(dat_sub_table)

sibdyads <- dat_sub_table[dat_sub_table$Freq %in% max(dat_sub$Run,na.rm=T),1] 

sibships <- dat_sub[dat_sub$dyad %in% sibdyads,]

sibships <- sibships[!duplicated(sibships$dyad),]

# tmp <- unlist(list(sibships$OffspringID1,(sibships$OffspringID2)))

# Add locations for the offspring 
sib.groups <- cbind.data.frame(
	sampleIDfile[match(sibships$OffspringID1, sampleIDfile$Ind),],
	sampleIDfile[match(sibships$OffspringID2, sampleIDfile$Ind),]
	)
names(sib.groups) <- c("Ind1","Sample1","Group1","Ind2","Sample2","Group2")

sib.groups$SameDiff <- ifelse(sib.groups$Group1=="","Different",
							ifelse(sib.groups$Group1==sib.groups$Group2,"Same","Different"))

# Results
# threshold=0.8: 2 out of 14 dyads in the seagrass sample where on the same blade
sib.groups
       # Ind1 Sample1 Group1   Ind2 Sample2 Group2  SameDiff
# 2    ML_A02     A02      1 ML_T19     T19        Different
# 6    ML_A06     A06      1 ML_A77     A77     14 Different
# 7    ML_A07     A07      1 ML_A22     A22      4 Different
# 14   ML_A14     A14      3 ML_A69     A69     12 Different
# 19   ML_A19     A19      4 ML_A55     A55      9 Different
# 20   ML_A20     A20      4 ML_A70     A70     12 Different
# 22   ML_A22     A22      4 ML_A61     A61      9 Different
# 26   ML_A27     A27      4 ML_A56     A56      9 Different
# 26.1 ML_A27     A27      4 ML_A67     A67     11 Different
# 28   ML_A29     A29      4 ML_A64     A64     10 Different
# 28.1 ML_A29     A29      4 ML_T01     T01        Different
# 31   ML_A33     A33      5 ML_T42     T42        Different
# 34   ML_A36     A36     15 ML_T19     T19        Different
# 36   ML_A38     A38     15 ML_T25     T25        Different
# 37   ML_A39     A39      7 ML_T39     T39        Different
# 43   ML_A45     A45      7 ML_A71     A71     13 Different
# 46   ML_A48     A48      7 ML_A74     A74     14 Different
# 48   ML_A50     A50      8 ML_T20     T20        Different
# 49   ML_A51     A51      8 ML_A52     A52      8      Same
# 49.1 ML_A51     A51      8 ML_A59     A59      9 Different
# 52   ML_A54     A54      9 ML_A59     A59      9      Same
# 55   ML_A57     A57      9 ML_T41     T41        Different
# 61   ML_A64     A64     10 ML_T01     T01        Different
# 61.1 ML_A64     A64     10 ML_T54     T54        Different
# 65   ML_A68     A68     11 ML_T09     T09        Different
# 69   ML_A72     A72     13 ML_T55     T55        Different
# 77   ML_T03     T03        ML_T35     T35        Different
# 78   ML_T04     T04        ML_T33     T33        Different
# 104  ML_T31     T31        ML_T40     T40        Different
# 110  ML_T37     T37        ML_T55     T55        Different
# 119  ML_T47     T47        ML_T51     T51        Different


sib.groups[sib.groups$SameDiff=="Same",] # 2 dyads each on different blades (8,9)
# 103   ML_A51     A51      8 ML_A52     A52      8     Same
# 106   ML_A54     A54      9 ML_A59     A59      9     Same


# threshold=0.9: 1 out of 4 dyads in the seagrass sample where on the same blade
sib.groups
       # Ind1 Sample1 Group1   Ind2 Sample2 Group2  SameDiff
# 19   ML_A19     A19      4 ML_A55     A55      9 Different
# 20   ML_A20     A20      4 ML_A70     A70     12 Different
# 28   ML_A29     A29      4 ML_A64     A64     10 Different
# 28.1 ML_A29     A29      4 ML_T01     T01        Different
# 31   ML_A33     A33      5 ML_T42     T42        Different
# 36   ML_A38     A38     15 ML_T25     T25        Different
# 48   ML_A50     A50      8 ML_T20     T20        Different
# 52   ML_A54     A54      9 ML_A59     A59      9      Same
# 61   ML_A64     A64     10 ML_T01     T01        Different
# 61.1 ML_A64     A64     10 ML_T54     T54        Different
# 65   ML_A68     A68     11 ML_T09     T09        Different
# 77   ML_T03     T03        ML_T35     T35        Different
# 78   ML_T04     T04        ML_T33     T33        Different

sib.groups[sib.groups$SameDiff=="Same",] # 1 dyad on blade 9
# 106  ML_A54     A54      9 ML_A59     A59      9     Same


# threshold=0.95: 1 out of 3 dyads in the seagrass sample where on the same blade
sib.groups
      # Ind1 Sample1 Group1   Ind2 Sample2 Group2  SameDiff
# 19   ML_A19     A19      4 ML_A55     A55      9 Different
# 20   ML_A20     A20      4 ML_A70     A70     12 Different
# 28   ML_A29     A29      4 ML_A64     A64     10 Different
# 28.1 ML_A29     A29      4 ML_T01     T01        Different
# 31   ML_A33     A33      5 ML_T42     T42        Different
# 36   ML_A38     A38     15 ML_T25     T25        Different
# 48   ML_A50     A50      8 ML_T20     T20        Different
# 52   ML_A54     A54      9 ML_A59     A59      9      Same
# 61   ML_A64     A64     10 ML_T01     T01        Different
# 61.1 ML_A64     A64     10 ML_T54     T54        Different
# 65   ML_A68     A68     11 ML_T09     T09        Different
# 77   ML_T03     T03        ML_T35     T35        Different
# 78   ML_T04     T04        ML_T33     T33        Different

sib.groups[sib.groups$SameDiff=="Same",] # 1 dyad on blade 9
# 06 ML_A54     A54      9 ML_A59     A59      9     Same

