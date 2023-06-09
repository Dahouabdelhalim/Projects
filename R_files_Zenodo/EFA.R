if (!require("psych")){
  install.packages("psych", dep = T)
}
require("psych")

if (!require("GPArotation")){
  install.packages("GPArotation", dep = T)
}
require("GPArotation")

source("import_anon_survey2020.r")
source("filtering.r")

ds<-filter_finished(ds)

ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_02", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_09","OP01_10","OP01_11","OP01_12", "OP01_13", "OP01_14", "OP01_15"))]

#Remove rows with missing values (N=80)
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)

#descriptive stats
describe(ds.beforeEFA)

#remove OP01_10 and OP01_11 since they are skewed
#repeat the entire loading process, since there might now be more complete rows
ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_02", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_09","OP01_12", "OP01_13", "OP01_14", "OP01_15"))]
#Remove rows with missing values
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)

#check again - we have five more complete answers (N=85)
describe(ds.beforeEFA)

#run KMO
KMO(ds.beforeEFA)

#Low overall MSA 0.68, let's remove all individuals with MSA <0.4 (OP01_12)
#again, we reload the whole thing to get more complete answers
ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_02", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_09", "OP01_13", "OP01_14", "OP01_15"))]
#Remove rows with missing values
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)

#check again - we have four more complete answers (N=89)
describe(ds.beforeEFA)

#redo KMO - 0.73, acceptable starting value!
KMO(ds.beforeEFA)

#Shapiro-Wilk to check for normality - not normal -> we'll run PAF
apply(ds.beforeEFA, 2, shapiro.test)

#Bartlet sphericity test - significant
cortest.bartlett(ds.beforeEFA)

#Determining no of factors
#Scree, Kaiser's rule, and parallel analysis
scree(ds.beforeEFA)
parallel = fa.parallel(ds.beforeEFA, fm = "pa", fa = "fa")
#Kaiser's eigenvalue rule says 2 factors, scree plot as well
#Parallel says 3
#--> Let's start at two factors and make our way up.

#PAF FA with two factors and promax rotation (we don't expect orthogonal factors)
fa =fa(ds.beforeEFA, nfactors = 2, fm = "pa", rotate = "promax")
print(fa, cut = 0.3, digits = 3)

#We exclude OP01_09 based on factor loadings (and repeat our procedures from before: loading data, removing missing)
ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_02", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_13", "OP01_14", "OP01_15"))]
#Remove rows with missing values
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)
#check again - we have one more complete answer (N=90)
describe(ds.beforeEFA)

#PAF FA with two factors and promax rotation (we don't expect orthogonal factors)
fa =fa(ds.beforeEFA, nfactors = 2, fm = "pa", rotate = "promax")
print(fa, cut = 0.3, digits = 3)

#We remove now variables with low h2 (<.25), one at a time (starting with lowest, OP01_02)
ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_13", "OP01_14", "OP01_15"))]
#Remove rows with missing values
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)
#check again - we have eight more complete answers (N=98)
describe(ds.beforeEFA)
#PAF FA with two factors and promax rotation (we don't expect orthogonal factors)
fa =fa(ds.beforeEFA, nfactors = 2, fm = "pa", rotate = "promax")
print(fa, cut = 0.3, digits = 3)

#next, OP01_13
ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_14", "OP01_15"))]
#Remove rows with missing values
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)
#check again - we have eight more complete answers (N=106)
describe(ds.beforeEFA)
fa =fa(ds.beforeEFA, nfactors = 2, fm = "pa", rotate = "promax")
print(fa, cut = 0.3, digits = 3)

#next, OP01_14
ds.beforeEFA <- ds[, which(names(ds) %in% c("OP01_01", "OP01_03", "OP01_04", "OP01_05", "OP01_06", "OP01_07", "OP01_08", "OP01_15"))]
#Remove rows with missing values
ds.beforeEFA[ds.beforeEFA == -1] <- NA
ds.beforeEFA <- na.omit(ds.beforeEFA)
#check again - we have one more complete answer (N=107)
describe(ds.beforeEFA)

#We do this in retrospective to report the KMO and sphericity numbers, since we know it's the last round.
#redo KMO - 0.83, Good
KMO(ds.beforeEFA)
#Shapiro-Wilk to check for normality - nothing has changed -> we'll run PAF
apply(ds.beforeEFA, 2, shapiro.test)
#Bartlet sphericity test - significant
cortest.bartlett(ds.beforeEFA)

#Final FA with 2 factors
fa =fa(ds.beforeEFA, nfactors = 2, fm = "pa", rotate = "promax")
print(fa, cut = 0.3, digits = 3)

#all h2 are now >.25, factor loadings are large enough >.32, no cross-loadings. NICE!
#cumulative variance explained: 0.534
#correlation PA1 and PA2 is 0.495
#RMSR=0.033, Chi Square prob <0.271 (non sign), TLI 0.9811, RMSEA 0.0493
#           PA1    PA2    h2    u2  com
#OP01_01  0.539        0.525 0.475 1.52
#OP01_03        -0.566 0.496 0.504 1.30
#OP01_04  0.840        0.754 0.246 1.01
#OP01_05  0.797        0.636 0.364 1.00
#OP01_06  0.624        0.429 0.571 1.02
#OP01_07         0.797 0.670 0.330 1.01
#OP01_08  0.721        0.411 0.589 1.18
#OP01_15         0.657 0.348 0.652 1.14

#We stop here!

#Cronbach alphas
PA1<-c("OP01_01","OP01_04","OP01_05","OP01_06","OP01_08")
PA2<-c("OP01_03","OP01_07","OP01_15")
alpha.pa1 <- psych::alpha(ds.beforeEFA[PA1])
alpha.pa2 <- psych::alpha(ds.beforeEFA[PA2],check.keys=TRUE)
#pa1: raw alpha 0.84, cannot be increased by removing, r.cor>0.5
#pa2: raw alpha 0.72, cannot be increased by removing, r.cor>0.5
#Nice!

#Final dataframes
PA1.df <- ds.beforeEFA[PA1]
PA2.df <- ds.beforeEFA[PA2]
print("PA1")
c(c(1, NCOL(PA1.df) * 5), c(min(rowSums(PA1.df)), max(rowSums(PA1.df))))
print("PA2")
c(c(1, NCOL(PA2.df) * 5), c(min(rowSums(PA2.df)), max(rowSums(PA2.df))))

after.efa.df <- data.frame(PA1 = rowSums(PA1.df), PA2 = rowSums(PA2.df))
describe(after.efa.df)

error.bars(after.efa.df, ylab = "Score", xlab = "Factor", main = "Confidence intervals \\n for PA1, PA2", eyes = F, sd = F)