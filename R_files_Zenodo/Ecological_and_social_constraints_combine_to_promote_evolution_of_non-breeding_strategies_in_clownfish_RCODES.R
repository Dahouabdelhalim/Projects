setwd("~/Desktop/PNG")
library(readxl)

data <- read_excel("ratio_R2-R3.xlsx")
View(data)
 
hist(data$ratio,main= "R2-R3 ratio", xlab="frequency") 

hist(data$ratio,xlim=c(0.4,0.9),ylim=c(0, 30), main= "R2-R3 ratio", xlab="R2-R3 RATIO") 


# ECOLOGICAL CONSTRAINTS FIRST EXPERIMENT: 
#TEST 1 (MALE VS EMPTY) -EFFECT OF VARIATION IN THE ALTERNATIVE OPTION   - DO NON BREEDRS/R3 DISPERSE TO AN OUTSIDE OPTION PLACED AT 0.5 M? 
TEST1 <- matrix(c(16, 15, 0,1), ncol = 2)
colnames(TEST1) <- c("NOT moved TEST 1","moved TEST 1") # <- gives column names
rownames(TEST1) <- c("EMPTY", "MALE") # <- gives row names
fisher.test(TEST1, conf.int = TRUE, conf.level = 0.99)
barplot(TEST1, beside=TRUE)


# ECOLOGICAL CONSTRAINTS SECOND EXPERIMENT: 
#TEST 2 (MALE VS EMPTY) -EFFECT OF VARIATION IN THE ALTERNATIVE OPTION   - DO NON BREEDRS/R3 RETURN TO THEIR HOME ANEMONE WHEN DISPLACED INTO OUTSIDE OPTIONS PLACED AT 0.5 M?
TEST2 <- matrix(c(9, 13, 7,3), ncol = 2)
colnames(TEST2) <- c("moved TEST 2","NOT moved TEST 2") # <- gives column names
rownames(TEST2) <- c("MALE", "EMPTY") # <- gives row names
fisher.test(TEST2, conf.int = TRUE, conf.level = 0.99)
barplot(TEST2, beside=TRUE, legend=TRUE)


# TEST 3 (MALE VS EMPTY) -EFFECT OF VARIATION IN THE ALTERNATIVE OPTION   - DO NON BREEDRS/R3 RETURN TO THEIR HOME ANEMONE WHEN DISPLACED INTO OUTSIDE OPTIONS PLACED AT 5 M?
TEST3<- matrix(c(0, 0, 16,16), ncol = 2)
colnames(TEST3) <- c("moved TEST 2","NOT moved TEST 2") # <- gives column names
rownames(TEST3) <- c("MALE", "EMPTY") # <- gives row names
fisher.test(TEST3, conf.int = TRUE, conf.level = 0.99)
barplot(TEST3, beside=TRUE, legend=TRUE)


# TEST 4 EMPTY OPTION ONLY: EFFECT OF DISTANCE 0.5 M VS 5 M DISTANCE -DO NON BREEDRS/R3 RETURN TO THEIR HOME ANEMONE WHEN DISPLACED INTO EMPTY OUTSIDE OPTIONS PLACED AT 5 M VS O.5 M?
TEST4<- matrix(c(13, 0, 3,16), ncol = 2)
colnames(TEST4) <- c("EMPTY moved"," EMPTY NOT moved") # <- gives column names
rownames(TEST4) <- c("TEST2", "TEST3") # <- gives row names
fisher.test(TEST4, conf.int = TRUE, conf.level = 0.99)
barplot(TEST4, beside=TRUE, legend=TRUE)

# TEST 5 BREEDING MALE OPTION ONLY: EFFECT OF DISTANCE 0.5 M VS 5 M DISTANCE -DO NON BREEDRS/R3 RETURN TO THEIR HOME ANEMONE WHEN DISPLACED INTO OUTSIDE OPTIONS WITH A BREEDING MALE PLACED AT 5 M VS O.5 M?
TEST5<- matrix(c(9, 0, 7,16), ncol = 2)
colnames(TEST5) <- c("MALE moved"," MALE NOT moved") # <- gives column names
rownames(TEST5) <- c("TEST2", "TEST3") # <- gives row names
fisher.test(TEST5, conf.int = TRUE, conf.level = 0.99)
barplot(TEST5, beside=TRUE, legend=TRUE)



# ECOLOGICAL CONSTRAINTS: EXPERIMENT 1 (MALE+EMPTY) VS EXPERIMENT 2   -MOVEMENTS IN THE 1ST EXP VS MOVEMENTS IN THE SECOND EXP A5 O.5. M   
TESTE<- matrix(c(1, 22, 31,10), ncol = 2)
colnames(TESTE) <- c(" moved"," NOT moved") # <- gives column names
rownames(TESTE) <- c("TEST1", "TEST2") # <- gives row names
fisher.test(TESTE, conf.int = TRUE, conf.level = 0.99)
barplot(TESTE, beside=TRUE, legend=TRUE)


# SOCIAL CONTRAINTS: EVICTIONS/SMALL VS BIG
x <- matrix(c(12, 3, 4,13), ncol = 2)
colnames(x) <- c("EVICTED","TOLLERATED") # <- gives column names
rownames(x) <- c("Bigger R3", "Smaller R3") # <- gives row names
fisher.test(x, conf.int = TRUE, conf.level = 0.99)
barplot(x, legend=TRUE,beside=TRUE, main= "Eviction status by size", ylab = "Frequency of outocmes", ylim=c(0,16))
                                                                                         