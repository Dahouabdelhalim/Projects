#############################################################################################################################################-
##
## Title: Floral scents of a deceptive plant are hyperdiverse and under population-specific phenotypic selection
## Authors: GFRERER Eva*, LAINA Danae, GIBERNAU Marc, FUCHS Roman, HAPP Martin**, TOLASCH Till, TRUTSCHNIG Wolfgang, HÖRGER Anja C., COMES Hans Peter, and DÖTTERL Stefan
## Date: 2021-05-20
## Description: File to simulate impact of non-detects on estimate in elastic net
## 
## Questions regarding the study contact *eva.gfrerer[at]stud.sbg.ac.at or **stefan.doetterl[at]plus.ac.at
## Questions regarding the code contact *eva.gfrerer[at]stud.sbg.ac.at or **martin.happ[at]plus.ac.at
#############################################################################################################################################-


library(pacman)
pacman::p_load(readxl, MASS, glmnet)


#---------------------------------------------------------------------------#
## 0. Read in Data, Define Helper Functions and Simulation Settings----
#---------------------------------------------------------------------------#

A <- read_xlsx(path="data/Data_Arum_scent_selection_Gfrerer_etal.xlsx",
               sheet = 1, skip = 2) #transform data or choose relevant sheet (2 = JOS, 3 = DAO)
N <- subset(N, N$popul == "JOS")
S <- subset(S, S$popul == "DAO")



f <- function(x, n){ 
  z <- x[x!=0]
  sample(x = z, size = n, replace = TRUE)
}

g <- function(x){ 
  z <- x[x!=0]
  return(length(z))
}

# Simulation Settings
alpha <- 0
prob <- 0.5
R <- 10^4


#---------------------------------------------------------------------------#
## 1. Simulation for northern JOS population ----
#---------------------------------------------------------------------------#

B <- N[, -c(1:8)]
B <- as.data.frame(B)


 select <- apply(B, 2, g) >= 5  
 #sum(select)
 
 B <- B[, select]
 d <- dim(B)[2]
 n <- dim(B)[1]
 
 



result <- as.data.frame(matrix(0, nrow = R, ncol = d))
result_tilde <- as.data.frame(matrix(0, nrow = R, ncol = d))






nn <- c(0.3,0.4,0.5,0.6)
error_north <- data.frame(zero_percentage = nn, completeSens = 0, zerosSens = 0, 
                          completeSpec = 0, zerosSpec = 0, completeAcc = 0, zerosAcc = 0)

set.seed(0)
for(k in 1:length(nn)){
  result <- rep(0,R)
  result0 <- rep(0,R)
  result_coef <- rep(0,R)
  result0_coef <- rep(0,R)
  result_acc <- rep(0,R)
  result0_acc <- rep(0,R)
  
  for(i in 1:R) {


    a <- runif(dim(B)[2])*10
    b <- 3
    
    x <- apply(B, 2, f, n)
    e <- rnorm(n, 0, sqrt(2))
    
    a0 <- rbinom(dim(B)[2], 1, prob = prob)
    a2 <- a*a0
    i0 <- (1:dim(B)[2]*(1-a0))
    i0 <- i0[i0>0]
    i1 <- (1:dim(B)[2]*a0)
    i1 <- i1[i1>0]
    
    y <- a2%*%t(x)+b+e
    l <- cv.glmnet(x, t(y), family = "gaussian", alpha = alpha)$lambda.min
    m <- glmnet(x, t(y), family = "gaussian", alpha = alpha, lambda = l)
    cm <- coef(m)[,1]
    c1 <- ifelse(cm > 0, 1, 0)
    
    
    df <- data.frame(a2=a2, c1=c1[-1])
    df$both <- df$a2*df$c1
    df$a22 <- ifelse(a2 > 0, 1, 0)
    df$only_c1 <- df$c1*(1-df$a22)
    df$TN <- (1-df$a22)*(1-df$c1)
    
    result[i] <- sum(df$both>0)/sum(df$a2>0)    
    cm <- cm[-1]
    result_coef[i] <- sum(df$only_c1)/sum(df$a2 == 0)
    result_acc <- (sum(df$both > 0) + sum(df$TN > 0))/dim(df)[1]
    
    
    x <- ifelse(x <= quantile(x, nn[k]), 0, x)
    
    l <- cv.glmnet(x, t(y), family = "gaussian", alpha = alpha)$lambda.min
    m <- glmnet(x, t(y), family = "gaussian", alpha = alpha, lambda = l)
    cm <- coef(m)[,1]
    c1 <- ifelse(cm > 0, 1, 0)
    
    
    df <- data.frame(a2=a2, c1=c1[-1])
    df$both <- df$a2*df$c1
    df$a22 <- ifelse(a2 > 0, 1, 0)
    df$only_c1 <- df$c1*(1-df$a22)
    df$TN <- (1-df$a22)*(1-df$c1)
    
    result0[i] <- sum(df$both>0)/sum(df$a2>0)
    cm <- cm[-1]
    result0_coef[i] <- sum(df$only_c1)/sum(df$a2 == 0)
    result0_acc <- (sum(df$both > 0) + sum(df$TN > 0))/dim(df)[1]
    
  }
  error_north[k, 2] <- mean(result)
  error_north[k, 3] <- mean(result0)
  error_north[k, 4] <- mean(result_coef)
  error_north[k, 5] <- mean(result0_coef)
}

# relative errors for the complete data and for the data with non-detects
error_north


#---------------------------------------------------------------------------#
## 2. Simulation for southern DAO population ---- 
#---------------------------------------------------------------------------#

B <- S[, -c(1:8)]
B <- as.data.frame(B)


select <- apply(B, 2, g) >= 40 #hier nicht 5?

B <- B[, select]
d <- dim(B)[2]
n <- dim(B)[1]


result <- as.data.frame(matrix(0, nrow = R, ncol = d))
result_tilde <- as.data.frame(matrix(0, nrow = R, ncol = d))


a <- runif(dim(B)[2])*10
b <- 3

nn <- c(0.3,0.4,0.5,0.6)
error_south <- data.frame(zero_percentage = nn, completeSens = 0, zerosSens = 0, completeSpec = 0, zerosSpec = 0, completeAcc = 0, zerosAcc = 0)


for(k in 1:length(nn)){
  set.seed(k)
  result <- rep(0,R)
  result0 <- rep(0,R)
  result_coef <- rep(0,R)
  result0_coef <- rep(0,R)
  result_acc <- rep(0,R)
  result0_acc <- rep(0,R)
  for(i in 1:R) {
    
    
    a <- runif(dim(B)[2])*10
    b <- 3
    
    x <- apply(B, 2, f, n)
    e <- rnorm(n, 0, sqrt(2))
    
    a0 <- rbinom(dim(B)[2], 1, prob = prob)
    a2 <- a*a0
    i0 <- (1:dim(B)[2]*(1-a0))
    i0 <- i0[i0>0]
    i1 <- (1:dim(B)[2]*a0)
    i1 <- i1[i1>0]
    
    y <- a2%*%t(x)+b+e
    l <- cv.glmnet(x, t(y), family = "gaussian", alpha = alpha)$lambda.min
    m <- glmnet(x, t(y), family = "gaussian", alpha = alpha, lambda = l)
    cm <- coef(m)[,1]
    c1 <- ifelse(cm > 0, 1, 0)
    
    
    df <- data.frame(a2=a2, c1=c1[-1])
    df$both <- df$a2*df$c1
    df$a22 <- ifelse(a2 > 0, 1, 0)
    df$only_c1 <- df$c1*(1-df$a22)
    df$TN <- (1-df$a22)*(1-df$c1)
    
    result[i] <- sum(df$both>0)/sum(df$a2>0)    #PPV, sensitivity
    cm <- cm[-1]
    result_coef[i] <- sum(df$only_c1)/sum(df$a2 == 0) #specificity
    result_acc <- (sum(df$both > 0) + sum(df$TN > 0))/dim(df)[1]
    
    
    
    x <- ifelse(x <= quantile(x, nn[k]), 0, x)
    
    l <- cv.glmnet(x, t(y), family = "gaussian", alpha = alpha)$lambda.min
    m <- glmnet(x, t(y), family = "gaussian", alpha = alpha, lambda = l)
    cm <- coef(m)[,1]
    c1 <- ifelse(cm > 0, 1, 0)
    
    
    df <- data.frame(a2=a2, c1=c1[-1])
    df$both <- df$a2*df$c1    
    df$a22 <- ifelse(a2 > 0, 1, 0)
    df$only_c1 <- df$c1*(1-df$a22)
    df$TN <- (1-df$a22)*(1-df$c1)
    
    result0[i] <- sum(df$both>0)/sum(df$a2>0)
    cm <- cm[-1]
    result0_coef[i] <- sum(df$only_c1)/sum(df$a2 == 0)
    result0_acc <- (sum(df$both > 0) + sum(df$TN > 0))/dim(df)[1]
    
  }
  error_south[k, 2] <- mean(result)
  error_south[k, 3] <- mean(result0)
  error_south[k, 4] <- mean(result_coef)
  error_south[k, 5] <- mean(result0_coef)
  error_south[k, 6] <- mean(result_acc)
  error_south[k, 7] <- mean(result0_acc)
}

# relative errors for the complete data and for the data with non-detects
error_south