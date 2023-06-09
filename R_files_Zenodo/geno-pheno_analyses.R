library("ggplot2")
library("MuMIn")
library("lme4")
library("lmerTest")

##### load and tidy data 

all_data <- read.csv("pipit_sample_information.csv", header = T)
head(all_data)

all_data <- subset(all_data, Species == "berthelots") # just want to look at the Berthelots data set

subset(all_data, Head_length < 30) # correct this error
all_data$Head_length[all_data$Sample == 250] <- 33.65

# create "Newsample" column 
all_data$Newsample <- paste0("PIP_",all_data$Sample)

# just look at data where no individuals have NAs 
all_data2 <- all_data[complete.cases(all_data[,c("Weight","Wing_length","Bill_width","Head_length","Tarsus_length","Bill_height","Bill_length")]),]
all_data2$Newsample <- paste("PIP",all_data2$Sample,sep = "_")

morph <- all_data2[,c("Weight","Wing_length","Bill_width","Head_length","Tarsus_length","Bill_height","Bill_length")]

colnames(morph) <- c("WT", "WL", "BW", "HL", "TL", "BH", "BL")
  


### load in the allele frequencies for the ADAM12 SNP 

SNP_1585s94_ped <- read.table("Berthelots_SNP1585s94.ped")
str(SNP_1585s94_ped)

#Check sample IDs match with original data
all(all_data$Newsample == SNP_1585s94_ped$V2)

#If so subset based on morphometric NAs
SNP_1585s94_ped2 <- SNP_1585s94_ped[complete.cases(all_data[,c("Weight","Wing_length","Bill_width","Head_length","Tarsus_length","Bill_height","Bill_length")]),]

#Check it's worked
all(all_data2$Newsample == SNP_1585s94_ped2$V2)


#Other SNP

SNP_1585s112_ped <- read.table("Berthelots_SNP1585s112.ped")

SNP_1585s112_ped2 <- SNP_1585s112_ped[complete.cases(all_data[,c("Weight","Wing_length","Bill_width","Head_length","Tarsus_length","Bill_height","Bill_length")]),]

all(all_data2$Newsample == SNP_1585s112_ped2$V2)


# replace genotype numbers with genotype letters 
SNP_1585s94_ped2$V7 <- gsub("2", "C", SNP_1585s94_ped2$V7)
SNP_1585s94_ped2$V7 <- gsub("4", "T", SNP_1585s94_ped2$V7)

SNP_1585s94_ped2$V8 <- gsub("2", "C", SNP_1585s94_ped2$V8)
SNP_1585s94_ped2$V8 <- gsub("4", "T", SNP_1585s94_ped2$V8)


SNP_1585s112_ped2$V7 <- gsub("3", "G", SNP_1585s112_ped2$V7)
SNP_1585s112_ped2$V7 <- gsub("1", "A", SNP_1585s112_ped2$V7)

SNP_1585s112_ped2$V8 <- gsub("3", "G", SNP_1585s112_ped2$V8)
SNP_1585s112_ped2$V8 <- gsub("1", "A", SNP_1585s112_ped2$V8)


# combine the final two columns to give the overall genotype for each individual 

SNP_1585s94_ped2$Genotype <- paste0(SNP_1585s94_ped2$V7, SNP_1585s94_ped2$V8)
SNP_1585s112_ped2$Genotype <- paste0(SNP_1585s112_ped2$V7, SNP_1585s112_ped2$V8)

### SNP 158s94
y <- data.frame(mypca$scores)
y$Pop <- SNP_1585s94_ped2$V1
y$Genotype <- SNP_1585s94_ped2$Genotype

y$Geno_n <- as.numeric(factor(y$Genotype))-1

### remane and order the islands 
y$Pop <- gsub("C_01GRA", "GRA", y$Pop)
y$Pop <- gsub("C_02LZ", "LZ", y$Pop)
y$Pop <- gsub("C_03FV", "FV", y$Pop)
y$Pop <- gsub("C_04GC", "GC", y$Pop)
y$Pop <- gsub("C_05TF", "TF", y$Pop)
y$Pop <- gsub("C_06TEID", "TEID", y$Pop)
y$Pop <- gsub("C_07GOM", "GOM", y$Pop)
y$Pop <- gsub("C_08LP", "LP", y$Pop)
y$Pop <- gsub("C_09EH", "EH", y$Pop)
y$Pop <- gsub("M_11M", "M", y$Pop)
y$Pop <- gsub("M_12PS", "PS", y$Pop)
y$Pop <- gsub("M_13DG", "DG", y$Pop)
y$Pop <- gsub("S_10SG", "SG", y$Pop)

y$Pop <- factor(y$Pop, levels = c("EH","GOM","LP","TEID","TF","GC","FV","LZ","GRA","M","PS","DG", "SG"))
head(y)



# for the other SNP##

x <- data.frame(mypca$scores)
x$Pop <- SNP_1585s112_ped2$V1
x$Genotype <- SNP_1585s112_ped2$Genotype

x$Geno_n <- as.numeric(factor(x$Genotype))-1

x$Pop <- gsub("C_01GRA", "GRA", x$Pop)
x$Pop <- gsub("C_02LZ", "LZ", x$Pop)
x$Pop <- gsub("C_03FV", "FV", x$Pop)
x$Pop <- gsub("C_04GC", "GC", x$Pop)
x$Pop <- gsub("C_05TF", "TF", x$Pop)
x$Pop <- gsub("C_06TEID", "TEID", x$Pop)
x$Pop <- gsub("C_07GOM", "GOM", x$Pop)
x$Pop <- gsub("C_08LP", "LP", x$Pop)
x$Pop <- gsub("C_09EH", "EH", x$Pop)
x$Pop <- gsub("M_11M", "M", x$Pop)
x$Pop <- gsub("M_12PS", "PS", x$Pop)
x$Pop <- gsub("M_13DG", "DG", x$Pop)
x$Pop <- gsub("S_10SG", "SG", x$Pop)

x$Pop <- factor(x$Pop, levels = c("EH","GOM","LP","TEID","TF","GC","FV","LZ","GRA","M","PS","DG", "SG"))

head(x)

# add the ADAM12 genotypes into the dataframe as a column 
all_data2$Geno_n_x <- x$Geno_n
all_data2$Geno_n_y <- y$Geno_n



# GLLM --------------------------------------------------------------------

# Mixed effects model for SNP 1585s112 then replicate for the other SNP using "y", genotype as continous (number of copies of minor allele) ------------------------------

# mixed-effects model with gaussian error family for both SNPs. Model the genotype as continuous
head(all_data2)

#Tarsus length
model1_x <- lmer(Tarsus_length ~ Geno_n_x (1|Archipelago/Population) + Sex + Age,
                 REML = TRUE, data = all_data2, na.action = na.exclude) 

summary(model1_x)
r.squaredGLMM(model1_x)

#Head length
model1_x <- lmer(Head_length~ Geno_n_x + (1|Archipelago/Population) + Sex + Age,
                 REML = TRUE, data = all_data2, na.action = na.exclude) # significance

summary(model1_x)
r.squaredGLMM(model1_x)

# Check the assumptions of the models 
sresid <- resid(model1_x, type = "pearson") 
hist(sresid) # they look normal 
plot(model3_x)

# * FIGURE S11 *
ggplot(all_data2,aes(x = Genotype,y = Head_length))+
  theme_claudia()+
  ylab("Head length")+
  xlab("SNP 1585s112 Genotype")+
  geom_jitter(width = 0.1, cex = 0.7)+
  facet_wrap(~V1)

ggplot(all_data2,aes(y = head_length,x = Geno_n_x,col = Pop))+
  stat_smooth(method = "lm",se = F)

#Bill length
model1_x <- lmer(Bill_length ~ Geno_n_y + (1|Archipelago/Population) + Sex + Age,
                 REML = TRUE, data = all_data2, na.action = na.exclude) 

summary(model1_x)
r.squaredGLMM(model1_x)


# Check the correlation of the genotype variaibles
# co linear SNP variables. say that you checked for conlinearity of explanitary variables. >98%  
model_both_geno_1 <- lmer(Head_length ~ Geno_n_x + Geno_n_y + (1|Archipelago/Population) + Sex + Age,
                           REML = TRUE, data = all_data2, na.action = na.exclude) 

plot(all_data2$Head_length ~ all_data2$Sex)


