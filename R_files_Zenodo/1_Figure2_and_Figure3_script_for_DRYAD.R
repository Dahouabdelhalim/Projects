### Analysis script for:

# "Asymmetrical reproductive barriers in sympatric Jewelflowers: 
# are floral isolation, genetic incompatibilities, 
# and floral trait displacement connected?â€� (Manuscript ID: BJLS-6658)

# K. Christie, J.P. Doan, W.C. McBride, S.Y. Strauss 2021

# Biological Journal of the Linnean Society; (Manuscript ID: BJLS-6658)



### Figure 2 and Figure 3

# import data
setwd("~/Desktop/Dryad_data/")
dat <- read.csv(file = "1_preference_and_constancy_data_for_DRYAD.csv", header = T, stringsAsFactors = F)


### 01. summarize data ---------------------------------------------------------
nrow(dat) # 778 pollinator bouts

# total flowers visited
sum(dat$total_flws_visited) # 6122 total floral visits
sum(dat$B_flws_visited) # 3931 S.breweri flowers visited
sum(dat$H_flws_visited) # 2175 S.hesperidis flowers visited

# total plants visited
sum(dat$total_plants) # 1680 total floral visits
sum(dat$B_plants_visited) # 1169 S.breweri flowers visited
sum(dat$H_plants_visited) # 511 S.hesperidis flowers visited

# calculate proportion of plants visited during each pollinator bout
dat$pro_B_plants <- dat$B_plants_visited / dat$total_plants
dat$pro_H_plants <- dat$H_plants_visited / dat$total_plants

mean(dat$pro_B_plants) # 66.6%
mean(dat$pro_H_plants) # 33.4%


### 02. test for preference ----------------------------------------------------

# conduct chi.sq test across all all habitats
all_habitat <-  matrix(c(sum(dat$B_plants_visited),sum(dat$H_plants_visited)), nrow = 1, byrow = T)
all <- chisq.test(all_habitat)


# conduct individual chi.sq tests for each habitat


# breweri sites
b_dat <- subset(dat, outcrop == "breweri")
nrow(b_dat) # 359 foraging bouts

sum(b_dat$B_plants_visited) # 427
sum(b_dat$H_plants_visited) # 207

b_habitat <- matrix(c(sum(b_dat$B_plants_visited), sum(b_dat$H_plants_visited)), nrow = 1, byrow = T)
b <- chisq.test(b_habitat) ### X2 = 76.34 (correct on Table 1)
p.adjust(p = b$p.value, method = "bonferroni", n = 3) ### p = 7.16E-18 (correct on Table 1)


# hesperidis sites
h_dat <- subset(dat, outcrop == "hesperidis")
nrow(h_dat) # 250 foraging bouts

sum(h_dat$B_plants_visited) # 475
sum(h_dat$H_plants_visited) # 206

h_habitat <- matrix(c(sum(h_dat$B_plants_visited), sum(h_dat$H_plants_visited)), nrow = 1, byrow = T)
h <- chisq.test(h_habitat)
p.adjust(p = h$p.value, method = "bonferroni", n = 3)


# intermediate sites
in_dat <- subset(dat, outcrop == "between")
nrow(in_dat) # 169 

sum(in_dat$B_plants_visited) # 267
sum(in_dat$H_plants_visited) # 97

in_habitat <- matrix(c(sum(in_dat$B_plants_visited), sum(in_dat$H_plants_visited)), nrow = 1, byrow = T)
in_hab <- chisq.test(in_habitat)
p.adjust(p = in_hab$p.value, method = "bonferroni", n = 3)


### 03. test for difference in preference based on biotic environment ----------
library(DescTools)

# show number of plants visited across all observations in each habitat
sum(b_dat$B_plants_visited) # 427
sum(b_dat$H_plants_visited) # 207
sum(h_dat$B_plants_visited) # 475
sum(h_dat$H_plants_visited) # 206
sum(in_dat$B_plants_visited) # 267
sum(in_dat$H_plants_visited) # 97



table1 <- as.table(rbind(c(427, 475, 267), c(207,206,97)))
dimnames(table1) <- 
    list(
        species = c("breweri", "hesperidis"),
        site = c("breweri", "hesperidis", "intermediate"))

GTest(table1) # G = 3.9859, X-squared df = 2, p-value = 0.1363


# 04. plot Figure 2 ------------------------------------------------------------

dev.new()

library(scales)
par(font.lab = 2, cex.lab = 1.5, mar = c(5,5,5,5))

mp <- barplot(
    c(
        mean(b_dat$pro_B_plants), mean(b_dat$pro_H_plants),
        mean(h_dat$pro_B_plants), mean(h_dat$pro_H_plants),
        mean(in_dat$pro_B_plants), mean(in_dat$pro_H_plants)),
    ylim = c(0,1), col = alpha(c("blue","darkgreen"), alpha = 0.65),
    names = c("S. breweri", "S. hesperidis", "S. breweri", "S. hesperidis", "S. breweri", "S. hesperidis"),
    ylab = "Proportion of visits")


# add error bars
m1 <- mean(b_dat$pro_B)
se1 <- sd(b_dat$pro_B) / sqrt(nrow(b_dat)) # standard error

arrows(x0 = mp[1], y0 = (m1 - 2*se1),
       x1 = mp[1], y1 = (m1 + 2*se1),
       lwd = 2, angle = 90, code = 3)

m2 <- mean(b_dat$pro_H)
se2 <- sd(b_dat$pro_H) / sqrt(nrow(b_dat)) # standard error

arrows(x0 = mp[2], y0 = (m2 - 2*se2),
       x1 = mp[2], y1 = (m2 + 2*se2),
       lwd = 2, angle = 90, code = 3)

m3 <- mean(h_dat$pro_B)
se3 <- sd(h_dat$pro_B) / sqrt(nrow(h_dat)) # standard error

arrows(x0 = mp[3], y0 = (m3 - 2*se3),
       x1 = mp[3], y1 = (m3 + 2*se3),
       lwd = 2, angle = 90, code = 3)

m4 <- mean(h_dat$pro_H)
se4 <- sd(h_dat$pro_H) / sqrt(nrow(h_dat)) # standard error

arrows(x0 = mp[4], y0 = (m4 - 2*se4),
       x1 = mp[4], y1 = (m4 + 2*se4),
       lwd = 2, angle = 90, code = 3)

m5 <- mean(in_dat$pro_B)
se5 <- sd(in_dat$pro_B) / sqrt(nrow(in_dat)) # standard error

arrows(x0 = mp[5], y0 = (m5 - 2*se4),
       x1 = mp[5], y1 = (m5 + 2*se4),
       lwd = 2, angle = 90, code = 3)

m6 <- mean(in_dat$pro_H)
se6 <- sd(in_dat$pro_H) / sqrt(nrow(in_dat)) # standard error

arrows(x0 = mp[6], y0 = (m6 - 2*se4),
       x1 = mp[6], y1 = (m6 + 2*se4),
       lwd = 2, angle = 90, code = 3)

# add text
text(x = (mp[1] + mp[2])/2, y = 0.9, "S. breweri locations", font = 2, cex = 1.25)
text(x = (mp[3] + mp[4])/2, y = 0.9, "S. hesperidis locations", font = 2, cex = 1.25)
text(x = (mp[5] + mp[6])/2, y = 0.9, "Unoccupied locations", font = 2, cex = 1.25)

# add divider lines
segments(x0 = 2.5, x1 = 2.5, y0 = 0, y1 = 1, lwd = 2, col = "gray")
segments(x0 = 4.9, x1 = 4.9, y0 = 0, y1 = 1, lwd = 2, col = "gray")

# add p-values
text(x = (mp[1] + mp[2])/2, y = 0.8, "***", font = 2, cex = 1.25)
text(x = (mp[3] + mp[4])/2, y = 0.8, "***", font = 2, cex = 1.25)
text(x = (mp[5] + mp[6])/2, y = 0.8, "***", font = 2, cex = 1.25)



### 05. summarize transition frequencies ---------------------------------------

# isolate foraging bouts in which there was at least a single inter-plant transition
trans_dat <- subset(dat, n_trans > 0)

nrow(trans_dat) # 278 foraging bouts with at least one transition

# total transitions
sum(colSums(trans_dat[,15:18])) # 905 inter-plant transitions


# conduct chi.square test for transition frequencies

# breweri transitions
sum(trans_dat$B.B) # 524 breweri-breweri transitions
sum(trans_dat$B.H) # 132 breweri-hesperidis transitions

chisq.test(matrix(c(sum(trans_dat$B.B),sum(trans_dat$B.H)), nrow = 1, byrow = T))

# hesperidis transitions
sum(trans_dat$H.H) # 121 hesperidis-hesperidis transitions
sum(trans_dat$H.H) # 128 hesperidis-breweri transitions

chisq.test(matrix(c(sum(trans_dat$H.H),sum(trans_dat$H.H)), nrow = 1, byrow = T))


### is there a difference in transition probability based on biotic enviro? ----

# calculate average number of each transitions in each environment
b_dat_trans <- subset(trans_dat, outcrop == "breweri")
h_dat_trans <- subset(trans_dat, outcrop == "hesperidis")
in_dat_trans <- subset(trans_dat, outcrop == "between")

# average transitions
colMeans(b_dat_trans[,15:18])
colMeans(h_dat_trans[,15:18])
colMeans(in_dat_trans[,15:18])

table2 <- as.table(
    rbind(c(1.5377358, 0.3584906, 0.4716981, 0.2641509), 
          c(2.1929825, 0.5263158, 0.5350877, 0.5263158),
          c(1.9137931, 0.5862069, 0.1724138, 0.6896552)))

dimnames(table2) <- 
    list(
        species = c("breweri-sites", "hesperidis-sites", "intermediate-sites"),
        transition_type = c("B-B", "B-H", "H-H", "H-B"))

# no difference in transition rate based on background environment
GTest(table2) # G = 0.38199, X-squared df = 6, p-value = 0.999



### 06. Figure 3A (Transition frequencies) -------------------------------------
# plot average transitions per bout
dev.new()
par(font.lab = 2, cex.lab = 1.5, mar = c(5,5,5,1))
par(mfrow = c(1,2))

mp2 <- barplot(colMeans(trans_dat[,15:18]), ylim = c(0,2.5),
    names = c("B-B", "B-H", "H-H", "H-B"),
    ylab = "Transitions per foraging bout",
    col = alpha(c("blue","darkgray", "darkgreen", "darkgray"), alpha = 0.65)) # average number of transitions per foraging bout

# add error bars
m1 <- mean(trans_dat$B.B)
se1 <- sd(trans_dat$B.B) / sqrt(nrow(trans_dat)) # standard error

arrows(x0 = mp2[1], y0 = (m1 - 2*se1),
       x1 = mp2[1], y1 = (m1 + 2*se1),
       lwd = 2, angle = 90, code = 3)

m2 <- mean(trans_dat$B.H)
se2 <- sd(trans_dat$B.H) / sqrt(nrow(trans_dat)) # standard error

arrows(x0 = mp2[2], y0 = (m2 - 2*se2),
       x1 = mp2[2], y1 = (m2 + 2*se2),
       lwd = 2, angle = 90, code = 3)

m3 <- mean(trans_dat$H.H)
se3 <- sd(trans_dat$H.H) / sqrt(nrow(trans_dat)) # standard error

arrows(x0 = mp2[3], y0 = (m3 - 2*se3),
       x1 = mp2[3], y1 = (m3 + 2*se3),
       lwd = 2, angle = 90, code = 3)

m4 <- mean(trans_dat$H.B)
se4 <- sd(trans_dat$H.B) / sqrt(nrow(trans_dat)) # standard error

arrows(x0 = mp2[4], y0 = (m4 - 2*se4),
       x1 = mp2[4], y1 = (m4 + 2*se4),
       lwd = 2, angle = 90, code = 3)

# add p-values for chi-square test
text(x = 1.3, y = 2.4, "***", font = 2, cex = 1.75)
text(x = 3.7, y = 0.8, "NS", font = 2, cex = 1.25)



### 07. calculate RI based on average transition frequency ---------------------
colMeans(trans_dat[,15:18])

Calc.RI <- function(H,C){
    RI = 1 - 2 * (H / (H+C))
    
    return(RI)
}

### RI values in arrays, not necessarily in nature

# breweri (ovules) receive breweri pollen (C)
# breweri (ovules) receive hesperidis pollen (H)
Calc.RI(H = 0.4604317, C = 1.884) # 0.6072125

# hesperidis (ovules) receive hesperidis pollen (C)  
# hesperidis (ovules) receive hesperidis pollen (H)
Calc.RI(H = 0.4748201, C = 0.4352518) # -0.04347821

# breweri pollen goes to breweri ovule (C)
# breweri pollen goes to hesperidis ovule (H)
Calc.RI(H = 0.4748201, C = 1.884) # 0.5974088

# hesperidis pollen goes to hesperidis ovule (C)
# hesperidis pollen goes to breweri ovule (H)
Calc.RI(H = 0.4604317 , C = 0.4352518) # -0.0281125



### 07. calculate overall floral constancy -------------------------------------
# among all bouts with inter-plant transitions, 
# find plant species-specific constancy

# for a single plant species
# let c = observed proportion of conspecific transitions
# let e = expected proprtion of conspecific transitions
# based on overall frequency with which floral visitors visit focal species (preference)


# function to calculate Gegear and Thompson's Constancy Index (CI) -------------
Calc.CI <- function(c,e){
    # c = observed proportion of conspecific movements
    # e = expected proportion of conspecific movements
    # based on overall frequency of visitation
    CI = (c-e)/(c+e-2*c*e)
    return(CI)
}

# find "c" values (observed conspecific transition rate) per foraging bout
trans_dat$HH_c <- trans_dat$H.H / (trans_dat$H.H + trans_dat$H.B)
trans_dat$BB_c <- trans_dat$B.B / (trans_dat$B.B + trans_dat$B.H)


# find "e" values (expected based on preference; or proportion of plants visited)
BB_e <- sum(trans_dat$B_plants_visited) / (sum(trans_dat$B_plants_visited) + sum(trans_dat$H_plants_visited))
    # 0.7254237
HH_e <- sum(trans_dat$H_plants_visited) / (sum(trans_dat$B_plants_visited) + sum(trans_dat$H_plants_visited))
    # 0.2745763


# calculate species-specific CI values
trans_dat$CI_breweri  <- Calc.CI(c = c(trans_dat$BB_c), e = 0.7373102)
trans_dat$CI_hesperidis  <- Calc.CI(c = c(trans_dat$HH_c), e = 0.2626898)

mean(trans_dat$CI_breweri, na.rm = T) # 0.56
mean(trans_dat$CI_hesperidis, na.rm = T) # 0.11

# calculate standard error
mean_CI_B <- mean(trans_dat$CI_breweri, na.rm = T)
se_CI_B <- sd(trans_dat$CI_breweri, na.rm = T) / sqrt(nrow(trans_dat)) # standard error
lower_95_B <- mean_CI_B - 2*se_CI_B
upper_95_B <- mean_CI_B + 2*se_CI_B

mean_CI_H <- mean(trans_dat$CI_hesperidis, na.rm = T)
se_CI_H <- sd(trans_dat$CI_hesperidis, na.rm = T) / sqrt(nrow(trans_dat)) # standard error
lower_95_H<- mean_CI_H - 2*se_CI_H
upper_95_H <- mean_CI_H + 2*se_CI_H


### 08. Figure 3B (Constancy) --------------------------------------------------
# barplot
# dev.new()
# par(font.lab = 2, cex.lab = 1.5, mar = c(5,5,5,5))

bp1 <- barplot(c(mean_CI_B, mean_CI_H), ylim = c(0, 0.7),
    col = alpha(c("blue","darkgreen"), alpha = 0.65),
    names = c("S. breweri", "S. hesperidis"),
    ylab = "Constancy Index (CI)")


arrows(x0 = 0.7, y0 = lower_95_B,
       x1 = 0.7, y1 = upper_95_B,
       lwd = 2, angle = 90, code = 3)

arrows(x0 = 1.9, y0 = lower_95_H,
       x1 = 1.9, y1 = upper_95_H,
       lwd = 2, angle = 90, code = 3)






















