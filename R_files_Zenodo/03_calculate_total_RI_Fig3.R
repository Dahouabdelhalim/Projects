### 0. import data -------------------------------------------------------------
setwd("~/Documents/Post_Doc/01_Reproductive_Isolation_Review/4_Manuscript/files_for_Dryad_REVISED//")
list.files(pattern = ".csv")
dat <- read.csv(file = "RI_data_FINAL.csv", header = T, stringsAsFactors = F)


### 01. find average barrier strength for each taxa pair -----------------------

# add individual sterility columns 
  # combination of Pollen Sterility and Ovule Fertility 
dat$F1Sterility1 <- apply(dat[, c("F1PollenSterility1", "F1OvuleFertility1")], MARGIN = 1, mean, na.rm = T)
dat$F1Sterility2 <- apply(dat[, c("F1PollenSterility2", "F1OvuleFertility2")], MARGIN = 1, mean, na.rm = T)

# add individual mating system columns (for later asymmetry analyses)
  # combination of Mating System and Differential Pollen 
dat$MatingSystem1 <- apply(dat[, c("MatingSystem1", "DifferentialPollen1")], MARGIN = 1, mean, na.rm = T)
dat$MatingSystem2 <- apply(dat[, c("MatingSystem2", "DifferentialPollen2")], MARGIN = 1, mean, na.rm = T)


### 01. build function to calculate total RI -----------------------------------
# following Equation RI4E from Sobel and Chen 2014
# modeled after "total isolation" spreadsheet 
# provided in supplementary Excel file

# this function is designed to be applied across rows,
# in which rows represent taxa pairs,
# and columns represent individual strengths of different isolating barriers


# write function to calculate total RI

Calc.Total.RI <- function(RI_values){
    # here we assess 12 potential isolating barriers
    
    ### prezygotic barriers that affect co-occurrence ###
        # Barrier 1 - Ecogeographic Isolation
        # Barrier 2 - Immigrant Inviability
        # Barrier 3 - Phenology
    
    # Barrier 1
    RI1 = RI_values[1]
    shared1 = 1 - RI1
    unshared1 = RI1
    hetero1 = 0.5 * shared1
    con1 = 1 - hetero1
    cum1 = 1 - 2*(hetero1/(hetero1 + con1))
    
    # Barrier 2
    RI2 = RI_values[2]
    shared2 = 1 - RI2
    unshared2 = RI2
    hetero2 = 0.5 * shared1 * shared2
    con2 = 1 - hetero2
    cum2 = 1 - 2*(hetero2/(hetero2 + con2))
    
    # Barrier 3
    RI3 = RI_values[3]
    shared3 = 1 - RI3
    unshared3 = RI3
    hetero3 = 0.5 * shared1 * shared2 * shared3
    con3 = 1 - hetero3
    cum3 = 1 - 2*(hetero3/(hetero3 + con3))
    
    
    ### prezyogitc barriers that do not affect co-occurrence ###
        # Barrier 4 - Mating System
        # Barrier 5 - Floral Isolation
        # Barrier 6 - Pollen Pistil Interactions   
         
    # Barrier 4
    RI4 = RI_values[4]
    H4 = (1 - RI4) / 2
    C4 = 1 - H4
    hetero4 = (shared1*shared2*shared3 * H4) / (H4 + C4)  
    con4 = 1 - hetero4
    cum4 = 1 - 2*(hetero4 / (hetero4 + con4))
        
    # Barrier 5
    RI5 = RI_values[5]
    H5 = (1 - RI5) / 2
    C5 = 1 - H5
    hetero5 = (shared1*shared2*shared3 * H4*H5) / (H4*H5 + C4*C5)  
    con5 = 1 - hetero5
    cum5 = 1 - 2*(hetero5 / (hetero5 + con5))
    
    # Barrier 6
    RI6 = RI_values[6]
    H6 = (1 - RI6) / 2
    C6 = 1 - H6
    hetero6 = (shared1*shared2*shared3 * H4*H5*H6) / (H4*H5*H6 + C4*C5*C6)  
    con6 = 1 - hetero6
    cum6 = 1 - 2*(hetero6 / (hetero6 + con6))
    

    ### postzygotic barriers ###
        # Barrier 7 - Fruit Production 
        # Barrier 8 - Seed Production  
        # Barrier 9 - F1 Germination
        # Barrier 10 - F1 Viability 
        # Barrier 11 - F1 Sterility
        # Barrier 12 - Extrinsic Post

    # Barrier 7
    RI7 = RI_values[7]
    H7 = (1 - RI7) / 2
    C7 = 1 - H7
    hetero7 = hetero6*H7 
    con7 = con6 * C7
    cum7 = 1 - 2*(hetero7 / (hetero7 + con7)) 
    
    # Barrier 8
    RI8 = RI_values[8]
    H8 = (1 - RI8) / 2
    C8 = 1 - H8
    hetero8 = hetero6*H7*H8
    con8 = con6*C7*C8
    cum8 = 1 - 2*(hetero8 / (hetero8 + con8)) 
    
    # Barrier 9
    RI9 = RI_values[9]
    H9 = (1 - RI9) / 2
    C9 = 1 - H9
    hetero9 = hetero6*H7*H8*H9
    con9 = con6*C7*C8*C9
    cum9 = 1 - 2*(hetero9 / (hetero9 + con9)) 
    
    # Barrier 10
    RI10 = RI_values[10]
    H10 = (1 - RI10) / 2
    C10 = 1 - H10
    hetero10 = hetero6*H7*H8*H9*H10
    con10 = con6*C7*C8*C9*C10
    cum10 = 1 - 2*(hetero10 / (hetero10 + con10)) 
    
    # Barrier 11
    RI11 = RI_values[11]
    H11 = (1 - RI11) / 2
    C11 = 1 - H11
    hetero11 = hetero6*H7*H8*H9*H10*H11
    con11 = con6*C7*C8*C9*C10*C11
    cum11 = 1 - 2*(hetero11 / (hetero11 + con11))
    
    # Barrier 12
    RI12 = RI_values[12]
    H12 = (1 - RI12) / 2
    C12 = 1 - H12
    hetero12 = hetero6*H7*H8*H9*H10*H11*H12
    con12 = con6*C7*C8*C9*C10*C11*C12
    cum12 = 1 - 2*(hetero12 / (hetero12 + con12))
    
    # return total RI
    TOTAL_RI = cum12
    return(TOTAL_RI)
    
    }


### 02. write function to calculate total prezygotic RI ------------------------

Calc.Total.PRE <- function(RI_values){
  # here we assess 6 potential isolating barriers
  
  ### prezygotic barriers that affect co-occurrence ###
  # Barrier 1 - Ecogeographic Isolation
  # Barrier 2 - Immigrant Inviability
  # Barrier 3 - Phenology
  
  # Barrier 1
  RI1 = RI_values[1]
  shared1 = 1 - RI1
  unshared1 = RI1
  hetero1 = 0.5 * shared1
  con1 = 1 - hetero1
  cum1 = 1 - 2*(hetero1/(hetero1 + con1))
  
  # Barrier 2
  RI2 = RI_values[2]
  shared2 = 1 - RI2
  unshared2 = RI2
  hetero2 = 0.5 * shared1 * shared2
  con2 = 1 - hetero2
  cum2 = 1 - 2*(hetero2/(hetero2 + con2))
  
  # Barrier 3
  RI3 = RI_values[3]
  shared3 = 1 - RI3
  unshared3 = RI3
  hetero3 = 0.5 * shared1 * shared2 * shared3
  con3 = 1 - hetero3
  cum3 = 1 - 2*(hetero3/(hetero3 + con3))
  
  
  ### prezyogitc barriers that do not affect co-occurrence ###
  # Barrier 4 - Mating System
  # Barrier 5 - Floral Isolation
  # Barrier 6 - Pollen Pistil Interactions   
  
  # Barrier 4
  RI4 = RI_values[4]
  H4 = (1 - RI4) / 2
  C4 = 1 - H4
  hetero4 = (shared1*shared2*shared3 * H4) / (H4 + C4)  
  con4 = 1 - hetero4
  cum4 = 1 - 2*(hetero4 / (hetero4 + con4))
  
  # Barrier 5
  RI5 = RI_values[5]
  H5 = (1 - RI5) / 2
  C5 = 1 - H5
  hetero5 = (shared1*shared2*shared3 * H4*H5) / (H4*H5 + C4*C5)  
  con5 = 1 - hetero5
  cum5 = 1 - 2*(hetero5 / (hetero5 + con5))
  
  # Barrier 6
  RI6 = RI_values[6]
  H6 = (1 - RI6) / 2
  C6 = 1 - H6
  hetero6 = (shared1*shared2*shared3 * H4*H5*H6) / (H4*H5*H6 + C4*C5*C6)  
  con6 = 1 - hetero6
  cum6 = 1 - 2*(hetero6 / (hetero6 + con6))
  
  # return total RI
  TOTAL_PRE_RI = cum6
  return(TOTAL_PRE_RI)
  
}


### 03. write function to calculate total postzygotic RI -----------------------

Calc.Total.POST <- function(RI_values){
  
  ### to calculate total POST-zygotic RI
  ### assume that prezygotic barriers are absent 
  hetero6 = 0.5
  con6 = 0.5 
  
  ### postzygotic barriers ###
  # Barrier 7 - Fruit Production 
  # Barrier 8 - Seed Production  
  # Barrier 9 - F1 Germination
  # Barrier 10 - F1 Viability 
  # Barrier 11 - F1 Sterility
  # Barrier 12 - Extrinsic Post
  
  # Barrier 7
  RI7 = RI_values[7]
  H7 = (1 - RI7) / 2
  C7 = 1 - H7
  hetero7 = hetero6*H7 
  con7 = con6*C7
  cum7 = 1 - 2*(hetero7 / (hetero7 + con7)) 
  
  # Barrier 8
  RI8 = RI_values[8]
  H8 = (1 - RI8) / 2
  C8 = 1 - H8
  hetero8 = hetero6*H7*H8
  con8 = con6*C7*C8
  cum8 = 1 - 2*(hetero8 / (hetero8 + con8)) 
  
  # Barrier 9
  RI9 = RI_values[9]
  H9 = (1 - RI9) / 2
  C9 = 1 - H9
  hetero9 = hetero6*H7*H8*H9
  con9 = con6*C7*C8*C9
  cum9 = 1 - 2*(hetero9 / (hetero9 + con9)) 
  
  # Barrier 10
  RI10 = RI_values[10]
  H10 = (1 - RI10) / 2
  C10 = 1 - H10
  hetero10 = hetero6*H7*H8*H9*H10
  con10 = con6*C7*C8*C9*C10
  cum10 = 1 - 2*(hetero10 / (hetero10 + con10)) 
  
  # Barrier 11
  RI11 = RI_values[11]
  H11 = (1 - RI11) / 2
  C11 = 1 - H11
  hetero11 = hetero6*H7*H8*H9*H10*H11
  con11 = con6*C7*C8*C9*C10*C11
  cum11 = 1 - 2*(hetero11 / (hetero11 + con11))
  
  # Barrier 12
  RI12 = RI_values[12]
  H12 = (1 - RI12) / 2
  C12 = 1 - H12
  hetero12 = hetero6*H7*H8*H9*H10*H11*H12
  con12 = con6*C7*C8*C9*C10*C11*C12
  cum12 = 1 - 2*(hetero12 / (hetero12 + con12))
  
  # return total RI
  TOTAL_POST_RI = cum12
  return(TOTAL_POST_RI)
  
}


### 04. prepare data.frame for analysis ----------------------------------------

# isolate first taxon in pair
RI1_dat <- dat[, 
                c("Ecogeo1", "ImmigrantInviability1", "Pheno1",
                  "MatingSystem1", "FloralIsolation1", "PollenPistil1",
                  "FruitProduction1", "SeedProduction1", "F1Germination1",
                  "F1Viability1", "F1Sterility1", "ExtrinsicPost1")]

# isolate second taxon in pair
RI2_dat <- dat[, 
               c("Ecogeo2", "ImmigrantInviability2", "Pheno2",
                 "MatingSystem2", "FloralIsolation2", "PollenPistil2",
                 "FruitProduction2", "SeedProduction2", "F1Germination2",
                 "F1Viability2", "F1Sterility2", "ExtrinsicPost2")]

# convert NA values to zero
RI1_dat[is.na(RI1_dat)] <- 0
RI2_dat[is.na(RI2_dat)] <- 0


### 05.  test function ---------------------------------------------------------
Calc.Total.RI(c(0.57, 0.3, 0.1, 0, 0.5, 0.1, 0.2, 0.4, 0, 0.2, 0.4, 0.5))
Calc.Total.RI(c(as.numeric(RI1_dat[1,])))
apply(RI1_dat, MARGIN = 1, FUN = Calc.Total.RI)[1]


### 06. summarize pre-, post-, and total RI for all species pairs --------------
RI_df <- data.frame(
    total1 = apply(RI1_dat, MARGIN = 1, FUN = Calc.Total.RI),
    total2 = apply(RI2_dat, MARGIN = 1, FUN = Calc.Total.RI),
    pre1 = apply(RI1_dat, MARGIN = 1, FUN = Calc.Total.PRE),
    pre2 = apply(RI2_dat, MARGIN = 1, FUN = Calc.Total.PRE),
    post1 = apply(RI1_dat, MARGIN = 1, FUN = Calc.Total.POST),
    post2 = apply(RI2_dat, MARGIN = 1, FUN = Calc.Total.POST))

# find average values
RI_df$avg_total_RI <- 
    apply(RI_df[,c("total1", "total2")], MARGIN = 1, FUN = mean, na.rm = T)

RI_df$avg_pre_RI <- 
    apply(RI_df[,c("pre1", "pre2")], MARGIN = 1, FUN = mean, na.rm = T)

RI_df$avg_post_RI <- 
    apply(RI_df[,c("post1", "post2")], MARGIN = 1, FUN = mean, na.rm = T)

# add geography and taxa type
RI_df$geography <- dat$geography
RI_df$taxa_type <- dat$Taxa_type2

# RI_df$taxa_type <- factor(RI_df$taxa_type,
 #                          levels = c("ecotypes", "subspecies", "cytotypes", "species"))



### 07. calculate summary statistics -------------------------------------------

# total RI
summary(RI_df$avg_total_RI)
quantile(RI_df$avg_total_RI, probs = c(0.25, 0.75)) # 0.7994242 0.9996365 

length(RI_df$avg_total_RI[RI_df$avg_total_RI == 1]) 
# 24 of 89 taxa pairs showed complete RI
24/89

length(RI_df$avg_total_RI[RI_df$avg_total_RI > 0.95]) 
# 47 of 89 taxa pairs showed near complete RI
47/89


# do the average strengths of pre- vs. post differ?
wilcox.test(RI_df$avg_pre_RI, RI_df$avg_post_RI) # p < 0.001


# summary of prezygotic RI
summary(RI_df$avg_pre_RI)
quantile(RI_df$avg_pre_RI, probs = c(0.25, 0.75)) # 0.5967985 0.9840483

length(RI_df$avg_pre_RI[RI_df$avg_pre_RI == 1]) 
# 17 of 89 taxa pairs showed complete RI
17/89

length(RI_df$avg_pre_RI[RI_df$avg_pre_RI > 0.95]) 
# 32 of 89 taxa pairs showed near complete RI
32/89

length(RI_df$avg_pre_RI[RI_df$avg_pre_RI > 0.75]) 
# 54 of 89 taxa pairs showed substantial RI
54/89


# summary of postzygotic RI
summary(RI_df$avg_post_RI)
quantile(RI_df$avg_post_RI, probs = c(0.25, 0.75)) # 0.02553242 0.74302804  

length(RI_df$avg_post_RI[RI_df$avg_post_RI == 1]) 
# 7 of 89 taxa pairs showed complete RI
7/89

length(RI_df$avg_post_RI[RI_df$avg_post_RI > 0.95]) 
# 10 of 89 taxa pairs showed near complete RI
10/89

length(RI_df$avg_post_RI[RI_df$avg_post_RI > 0.75]) 
# 21 of 89 taxa pairs showed substantial RI
21/89

mean(RI_df$avg_post_RI) # 0.395858
median(RI_df$avg_post_RI) # 0.4133782


# how many pairs have prezygotic or postzygotic RI > 0.75?
length(RI_df$avg_pre_RI[RI_df$avg_pre_RI > 0.75]) # 53 prezygotic
length(RI_df$avg_post_RI[RI_df$avg_post_RI > 0.75]) # 22 prezygotic

# how many pairs have prezygotic and postzygotic RI > 0.75?
nrow(RI_df[RI_df$avg_pre_RI > 0.75 & RI_df$avg_post_RI > 0.75 , ]) # 15
15/89

nrow(RI_df[RI_df$avg_pre_RI < 0.75 & RI_df$avg_post_RI < 0.75 , ]) # 29
29/89


# are pre- and postzygotic RI correlated?
cor.test(RI_df$avg_pre_RI, RI_df$avg_post_RI, method = "pearson") # p-value = 0.27, r = -0.12


### 08. Plot pre- vs. post for individual taxa pairs (FIGURE 3) ----------------
library(ggplot2)

dev.new()

ggplot(RI_df, aes(x = avg_pre_RI, y = avg_post_RI)) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "solid") +
  annotate("rect", xmin = 0.75, xmax = 1, ymin = 0.75, ymax = 1, alpha = 0.5, fill = "black") +
  annotate("rect", xmin = 0.75, xmax = 1, ymin = -1, ymax = 0.75, alpha = 0.25, fill = "black") +
  annotate("rect", xmin = -0.5, xmax = 0.75, ymin = 0.75, ymax = 1, alpha = 0.25, fill = "black") +
  annotate("rect", xmin = -0.5, xmax = 0.75, ymin = -1, ymax = 0.75, alpha = 0.1, fill = "black") +
  geom_point(aes(fill = taxa_type, shape = geography), size = 2.5, alpha = 0.85) +
  scale_shape_manual(values = c(24,22,21), na.translate = F) +
  guides(fill = guide_legend(override.aes=list(shape=23))) +
  theme_minimal() +
  xlim(-0.5,1) +
  ylim(-1,1)+
  labs(x = "Total prezygotic RI", y = "Total postzygotic RI") +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.title.x = element_text(face = "bold", size = 16)) +
  theme(legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size = 12)) + 
  labs(shape = "Geography", fill = "Taxa type")

# quartz.save(file = "Figure3_FINAL.jpg", type = "jpg", dpi = 600)



### 09. assess correlation between 
### number of barriers assessed and overall RI ---------------------------------

# find average barrier strength for both taxa
dat$Ecogeo_mean <- apply(dat[, c("Ecogeo1", "Ecogeo2")], MARGIN = 1, mean, na.rm = T)
dat$ImmigrantInviability_mean <- apply(dat[, c("ImmigrantInviability1", "ImmigrantInviability2")], MARGIN = 1, mean, na.rm = T)
dat$Pheno_mean <- apply(dat[, c("Pheno1", "Pheno2")], MARGIN = 1, mean, na.rm = T)

# combine "Mating System" and "Differential Pollen (Production)" columns 
dat$MatingSystem_mean <- apply(dat[, c("MatingSystem1", "MatingSystem2", "DifferentialPollen1", "DifferentialPollen2")], MARGIN = 1, mean, na.rm = T)

dat$FloralIsolation_mean <- apply(dat[, c("FloralIsolation1", "FloralIsolation2")], MARGIN = 1, mean, na.rm = T)
dat$PollenPistil_mean <- apply(dat[, c("PollenPistil1", "PollenPistil2")], MARGIN = 1, mean, na.rm = T)
dat$FruitProduction_mean <- apply(dat[, c("FruitProduction1", "FruitProduction2")], MARGIN = 1, mean, na.rm = T)
dat$SeedProduction_mean <- apply(dat[, c("SeedProduction1", "SeedProduction2")], MARGIN = 1, mean, na.rm = T)
dat$F1Germination_mean <- apply(dat[, c("F1Germination1", "F1Germination2")], MARGIN = 1, mean, na.rm = T)
dat$F1Viability_mean <- apply(dat[, c("F1Viability1", "F1Viability2")], MARGIN = 1, mean, na.rm = T)

# combine "Pollen Sterility" and "Ovule Fertility" columns 
dat$F1Sterility_mean <- apply(dat[, c("F1PollenSterility1", "F1PollenSterility2", "F1OvuleFertility1", "F1OvuleFertility2")], MARGIN = 1, mean, na.rm = T)

dat$ExtrinsicPost_mean <- apply(dat[, c("ExtrinsicPost1", "ExtrinsicPost2")], MARGIN = 1, mean, na.rm = T)

# create count number of barriers each study measured
names(dat)
d2 <- dat[, 62:73]

n_total <- apply(d2[, c(1:12)], MARGIN = 1, function(x) length(x[!is.na(x)]))
n_pre <- apply(d2[, c(1:6)], MARGIN = 1, function(x) length(x[!is.na(x)]))
n_post = apply(d2[, c(7:12)], MARGIN = 1, function(x) length(x[!is.na(x)]))


# assess correlation
cor.test(RI_df$avg_pre_RI, n_pre) # p-value = 0.01336, cor = 0.2613714  
cor.test(RI_df$avg_post_RI, n_post) # p-value = 0.04094, cor = 0.2171569 
cor.test(RI_df$avg_total_RI, n_total) # p-value = 0.8077, cor = 0.02617111 


### 10. Compare barrier strength estimates to Lowry et al. 2008 ----------------

# subset taxa pairs in Lowry et al. 2008
RI_in_Lowry <- RI_df[grepl(pattern = "YES", x = dat$in_Lowry_2008),]

mean(RI_in_Lowry$avg_total_RI) # 0.797565
mean(RI_df$avg_total_RI) # 0.8424227

mean(RI_in_Lowry$avg_pre_RI) # 0.6953776
mean(RI_df$avg_pre_RI) # 0.733024

mean(RI_in_Lowry$avg_post_RI) # 0.426434
mean(RI_df$avg_post_RI) # 0.4061396

# find number of barriers measured
n_barriers_measured <- data.frame(
  in_Lowry_2008 = dat$in_Lowry_2008,
  n_measured = n_total
)

temp1 <- subset(n_barriers_measured, in_Lowry_2008 == "YES")
mean(temp1$n_measured) # 4.4
sd(temp1$n_measured) # 2.07

temp2 <- subset(n_barriers_measured, in_Lowry_2008 == "NO")
mean(temp2$n_measured) # 4.51
sd(temp2$n_measured) # 1.77




