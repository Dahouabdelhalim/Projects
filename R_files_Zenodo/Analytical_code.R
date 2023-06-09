# If using this code please cite both the archived dataset and the published article.

# According to your available memory, you may need to remove models as you progress through the code.

# Long-term email address: s.r.evans.02@cantab.net


##########################################
# Using data that includes male identity #
##########################################

rm(list = ls())
library(asreml)
library(nadiv)

# Attach phenotypic data
phenotypic.data <- read.csv("~/YOURDIRECTORY/great tit breeding data.csv", 
                            header = T, 
                            na.strings = "NA")
# Exclude experimental nests
phenotypic.data <- subset(phenotypic.data, experimental_manipulation == 0)

phenotypic.data <- phenotypic.data[!with(phenotypic.data, is.na(phenotypic.data$laydate)&is.na(phenotypic.data$clutch_size)),]  # Omit records where both laydate and clutch size values are missing
dim(phenotypic.data)

range(phenotypic.data$laydate, na.rm = T)  # Range of laydate values
range(phenotypic.data$clutch_size, na.rm = T)  # Range of clutch size values

# REMOVE EXTREME VALUES
quantile(phenotypic.data$laydate, na.rm = T, prob=c(0.005, 0.995))  # Range of laydate values
quantile(phenotypic.data$clutch_size, na.rm = T, prob=c(0.005, 0.995))  # Range of clutch size values

phenotypic.data <- phenotypic.data[phenotypic.data$clutch_size >=4 &
                                     phenotypic.data$clutch_size <=13 &
                                     phenotypic.data$laydate >=5 &
                                     phenotypic.data$laydate <=55,]
range(phenotypic.data$laydate, na.rm = T)  # Revised range of laydate values
range(phenotypic.data$clutch_size, na.rm = T)  # Revised range of clutch size values

dim(phenotypic.data)

phenotypic.data$female.animal      <- as.factor(phenotypic.data$female_ID)
phenotypic.data$female.ID          <- as.factor(phenotypic.data$female_ID)
phenotypic.data$YEAR               <- as.factor(phenotypic.data$year)
phenotypic.data$LAYDATE            <- as.numeric(phenotypic.data$laydate)
phenotypic.data$CLUTCHSIZE         <- as.numeric(phenotypic.data$clutch_size)
phenotypic.data$male.animal        <- as.factor(phenotypic.data$male_ID)
phenotypic.data$male.ID            <- as.factor(phenotypic.data$male_ID)
phenotypic.data$female.natal.brood <- as.factor(phenotypic.data$female_natal_brood)
phenotypic.data$male.natal.brood   <- as.factor(phenotypic.data$male_natal_brood)
phenotypic.data$female.maternal.ID <- as.factor(phenotypic.data$female_maternal_ID)
phenotypic.data$male.maternal.ID   <- as.factor(phenotypic.data$male_maternal_ID)

phenotypic.data$F.AGE.CLASS <- phenotypic.data$female.age
phenotypic.data$F.AGE.CLASS[phenotypic.data$F.AGE.CLASS > 1] <- "adult"
phenotypic.data$F.AGE.CLASS[phenotypic.data$F.AGE.CLASS == 1] <- "juvenile"
phenotypic.data$F.AGE.CLASS <- as.factor(phenotypic.data$F.AGE.CLASS)

phenotypic.data$M.AGE.CLASS <- phenotypic.data$male.age
phenotypic.data$M.AGE.CLASS[phenotypic.data$M.AGE.CLASS > 1] <- "adult"
phenotypic.data$M.AGE.CLASS[phenotypic.data$M.AGE.CLASS == 1] <- "juvenile"
phenotypic.data$M.AGE.CLASS <- as.factor(phenotypic.data$M.AGE.CLASS)

# Make a vector detailing all brood identities that appear in the dataset
brood.IDs <- union(levels(phenotypic.data$female.natal.brood), levels(phenotypic.data$male.natal.brood))
length(brood.IDs) # Number of identified broods

# 'Relevel' both of the brood indices
phenotypic.data$female.natal.brood <- factor(phenotypic.data$female.natal.brood, levels = brood.IDs)
phenotypic.data$male.natal.brood <- factor(phenotypic.data$male.natal.brood, levels = brood.IDs)

# Make a vector detailing all maternal identities that appear in the dataset
maternal.IDs <- union(levels(phenotypic.data$female.maternal.ID), levels(phenotypic.data$male.maternal.ID))
length(maternal.IDs)  # Number of identified mothers

# 'Relevel' both of the maternal identity indices
phenotypic.data$female.maternal.ID <- factor(phenotypic.data$female.maternal.ID, levels = maternal.IDs)
phenotypic.data$male.maternal.ID <- factor(phenotypic.data$male.maternal.ID, levels = maternal.IDs)

# Create variable detailing nestbox ID and allocate nestboxes missing from the location data the identity of a neighbouring nestbox
phenotypic.data$nestbox <- substring(phenotypic.data$nest_ID, 6, 20)
phenotypic.data$nestbox[phenotypic.data$nestbox=="C142"] <- "C42"
phenotypic.data$nestbox[phenotypic.data$nestbox=="C160"] <- "C160A"
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX0"] <- NA
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX3A"] <- "EX3"
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX23"] <- "EX22"
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX24"] <- "EX25"
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX32A"] <- "EX32"
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX64A"] <- "EX64"
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX65A"] <- "EX65" 
phenotypic.data$nestbox[phenotypic.data$nestbox=="EX87"] <- "EX86"
phenotypic.data$nestbox[phenotypic.data$nestbox=="O89"] <- "EX89"
phenotypic.data$nestbox[phenotypic.data$nestbox=="O90"] <- "EX90"
phenotypic.data$nestbox[phenotypic.data$nestbox=="O92"] <- "EX92"
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW12A"] <- "SW12"
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW93"] <- NA
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW94"] <- NA
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW96"] <- NA
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW98"] <- NA
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW99"] <- NA
phenotypic.data$nestbox[phenotypic.data$nestbox=="SW104B"] <- "SW104"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W10A"] <- "W10"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W10B"] <- "W10"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W20B"] <- "W20"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W20A"] <- "W20"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W22A"] <- "W22"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W42A"] <- "W42"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W43A"] <- "W43"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W48"] <- "W47"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W61A"] <- "W61"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W70"] <- "W69"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W81A"] <- "W81"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W90"] <- "W89"
phenotypic.data$nestbox[phenotypic.data$nestbox=="W98A"] <- "W98"

# Create variable detailing nestbox round
phenotypic.data$round <- substr(phenotypic.data$nestbox, 1, nchar(phenotypic.data$nestbox)-1)  # Remove final alphanumeric, to deal with 'A', 'B', etc variants of nestboxes
phenotypic.data$round <- gsub('[[:digit:]]+', '', phenotypic.data$round)
table(phenotypic.data$round)

phenotypic.data$round[phenotypic.data$round == "B"] <- "Marley"
phenotypic.data$round[phenotypic.data$round == "C"] <- "Bean"
phenotypic.data$round[phenotypic.data$round == "CP"] <- "Common Piece"
phenotypic.data$round[phenotypic.data$round == "EX"] <- "Extra"
phenotypic.data$round[phenotypic.data$round == "MP"] <- "Marley Plantation"
phenotypic.data$round[phenotypic.data$round == "O"] <- "Broad Oak"
phenotypic.data$round[phenotypic.data$round == "P"] <- "Pasticks"
phenotypic.data$round[phenotypic.data$round == "SW"] <- "Singing Way"
phenotypic.data$round[phenotypic.data$round == "W"] <- "Great Wood"
table(phenotypic.data$round)

# Add information on nestbox locations
nestbox.locations <- read.csv("~/YOUR DIRECTORY/nestbox coordinates.csv", header = T)
phenotypic.data <- merge(x = phenotypic.data,
                         y = nestbox.locations,
                         by.x = c("nestbox", "round"),  # Add round so that it is not duplicated
                         by.y = c("box1", "round"),
                         all.x = T,
                         all.y = F)
sum(is.na(phenotypic.data$x))  # Number of records lacking location data for the nestbox
phenotypic.data<- phenotypic.data[!is.na(phenotypic.data$x),]  # Remove records from nestboxes of unknown location
phenotypic.data$round <- as.factor(phenotypic.data$round)
dim(phenotypic.data)

phenotypic.data <- subset(phenotypic.data, year < 2017)  # Remove records from 2017 since they lack data on male IDs
phenotypic.data <- phenotypic.data[!is.na(phenotypic.data$male.ID),]  # Omit data where male identity is unavailable
# phenotypic.data <- phenotypic.data[!is.na(phenotypic.data$LAYDATE),]  # Omit nests where LAYDATE information is missing
dim(phenotypic.data)

length(unique(phenotypic.data$YEAR))  # Number of years represented
min(phenotypic.data$year)  # Earliest year represented
max(phenotypic.data$year)  # Most recent year represented
mean(table(phenotypic.data[, "YEAR"])[1:(max(phenotypic.data$year-min(phenotypic.data$year)+1))])  # Mean annual sample size
range(table(phenotypic.data[, "YEAR"])[1:(max(phenotypic.data$year-min(phenotypic.data$year)+1))])  # Minimum and maximum annual sample sizes
sum(!is.na(phenotypic.data$LAYDATE))  # Number of laydate observations
sum(!is.na(phenotypic.data$CLUTCHSIZE))  # Number of clutch size observations
length(unique(phenotypic.data$female.ID))  # Number of female great tits sampled
dim(phenotypic.data)[1]/length(unique(phenotypic.data$female.ID)) # Average number of samples per female
max(table(phenotypic.data[, "female.ID"]))  # Highest number of samples from an individual female
length(unique(phenotypic.data$male.ID))  # Number of male great tits represented
dim(phenotypic.data)[1]/length(unique(phenotypic.data$male.ID)) # Average number of samples per male
max(table(phenotypic.data[, "male.ID"]))  # Highest number of samples from an individual male
known.female.maternal <- phenotypic.data[!is.na(phenotypic.data$female.maternal.ID),]
length(unique(known.female.maternal$female.ID))  # Number of females with known maternal identity
length(unique(known.female.maternal$female.ID))/length(unique(phenotypic.data$female.ID))  # Proportion of females with known maternal identity
known.male.maternal <- phenotypic.data[!is.na(phenotypic.data$male.maternal.ID),]
length(unique(known.male.maternal$male.ID))  # Number of males with known maternal identity
length(unique(known.male.maternal$male.ID))/length(unique(phenotypic.data$male.ID))  # Proportion of males with known maternal identity
dim(known.female.maternal[!is.na(known.female.maternal$male.maternal.ID),])[1]  # Number of data for which both male and female maternal identity is known
dim(known.female.maternal[!is.na(known.female.maternal$male.maternal.ID),])[1]/dim(phenotypic.data)[1]  # Proportion of nests where maternal identity is known for both parents

pedigree <- read.csv("~/IGEs in GTs/Evolution/great tit pedigree (to 2016).csv",
                     header = T,
                     na.strings = NA)

ainv <- ainverse(pedigree)



###############################
###############################
## CLASSICAL BIVARIATE MODEL ##
###############################
###############################

classical.bivariate.model <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                    + at(trait, 1):F.AGE.CLASS
                                    + at(trait, 2):F.AGE.CLASS,
                                    random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                    + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                    + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                    + corgh(trait):vm(female.animal, ainv),
                                    residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                    data = phenotypic.data,
                                    workspace = 32e6,  # Boost memory
                                    maxit = 20)

summary(classical.bivariate.model)$varcomp[,c("component", "std.error", "bound")]
summary(classical.bivariate.model,  coef=T)$coef.fixed
wald.asreml(classical.bivariate.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

classic.LD.Vp <- vpredict(classical.bivariate.model, LD.Vp ~ (V2 + V8 + V11 + V15))  # Excluding year
classic.LD.plot2 <- vpredict(classical.bivariate.model, LD.plot2 ~ V2 / (V2+V8+V11+V15))
classic.LD.y2 <- vpredict(classical.bivariate.model, LD.y2 ~ V5 / (V2+V8+V11+V15))
classic.LD.pe2 <- vpredict(classical.bivariate.model,LD.pe2 ~ V8 / (V2+V8+V11+V15))
classic.LD.h2 <- vpredict(classical.bivariate.model, LD.h2 ~ V11 / (V2+V8+V11+V15))
classic.LD.r2 <- vpredict(classical.bivariate.model, LD.r2 ~ V15 / (V2+V8+V11+V15))
classic.LD.i <- vpredict(classical.bivariate.model, LD.i ~ V8+V11)
classic.LD.i2 <- vpredict(classical.bivariate.model, LD.i2 ~ (V8+V11)/(V2+V8+V11+V15))

classic.CS.Vp <- vpredict(classical.bivariate.model, CS.Vp ~ V3 + V9 + V12 + V16)
classic.CS.plot2 <- vpredict(classical.bivariate.model, CS.plot2 ~ V3 / (V3+V9+V12+V16))
classic.CS.y2 <- vpredict(classical.bivariate.model, CS.y2 ~ V6 / (V3+V9+V12+V16))
classic.CS.pe2 <- vpredict(classical.bivariate.model,CS.pe2 ~ V9 / (V3+V9+V12+V16))
classic.CS.h2 <- vpredict(classical.bivariate.model, CS.h2 ~ V12 / (V3+V9+V12+V16))
classic.CS.r2 <- vpredict(classical.bivariate.model, CS.r2 ~ V16 / (V3+V9+V12+V16))
classic.CS.i <- vpredict(classical.bivariate.model, CS.i ~ V9+V12)
classic.CS.i2 <- vpredict(classical.bivariate.model, CS.i2 ~ (V9+V12)/(V3+V9+V12+V16))

classic.COVplot <- vpredict(classical.bivariate.model, classic.COVplot ~ V1*sqrt(V2*V3))
classic.COVy <- vpredict(classical.bivariate.model, classic.COVy ~ V4*sqrt(V5*V6))
classic.COVfpe <- vpredict(classical.bivariate.model, classic.COVfpe ~ V7*sqrt(V8*V9))
classic.COVa <- vpredict(classical.bivariate.model, classic.COVa ~ V10*sqrt(V11*V12))
classic.COVr <- vpredict(classical.bivariate.model, classic.COVr ~ V14*sqrt(V15*V16))

classic.LD.Vp
classic.LD.plot2
classic.LD.y2
classic.LD.pe2
classic.LD.h2
classic.LD.r2
classic.LD.i
classic.LD.i2

classic.CS.Vp
classic.CS.plot2
classic.CS.y2
classic.CS.pe2
classic.CS.h2
classic.CS.r2
classic.CS.i
classic.CS.i2

classic.COVplot
classic.COVy
classic.COVfpe
classic.COVa
classic.COVr


##########################
# TESTING VARIANCE TERMS #
##########################

# FUNCTION FOR CALCULATING P-VALUE FROM A MIXTURE OF TWO CHI-SQUARED DISTRIBUTIONS
p.chimixer <- function(chi.sq, df1, df2){
  P.1 <- 1-pchisq(chi.sq, df1)
  P.2 <- 1-pchisq(chi.sq, df2)
  P_value <- (P.1 + P.2)/2
  return(P_value)
}

# PLOT
# Plot covariance
classical.bivariate.model_no.COVplot <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                               + at(trait, 1):F.AGE.CLASS
                                               + at(trait, 2):F.AGE.CLASS,
                                               random = ~ idh(trait, init = c(2, 0.05)):round
                                               + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                               + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                               + corgh(trait):vm(female.animal, ainv),
                                               residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                               data = phenotypic.data,
                                               workspace = 32e6,  # Boost memory
                                               maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.COVplot$loglik - classical.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Plot effect on laydate
fix.LD.plot <- c(1e-8, 0.05)
names(fix.LD.plot) <- c("F", "P")
classical.bivariate.model_no.LD.Vplot <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS,
                                                random = ~ idh(trait, init = fix.LD.plot):round
                                                + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                                + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                                + corgh(trait):vm(female.animal, ainv),
                                                residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                                data = phenotypic.data,
                                                workspace = 32e6,  # Boost memory
                                                maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.LD.Vplot$loglik - classical.bivariate.model_no.COVplot$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Plot effect on clutch size
fix.CS.plot <- c(2, 1e-8)
names(fix.CS.plot) <- c("P", "F")
classical.bivariate.model_no.CS.Vplot <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS,
                                                random = ~ idh(trait, init = fix.CS.plot):round
                                                + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                                + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                                + corgh(trait):vm(female.animal, ainv),
                                                residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                                data = phenotypic.data,
                                                workspace = 32e6,  # Boost memory
                                                maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.CS.Vplot$loglik - classical.bivariate.model_no.COVplot$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# YEAR
# covariance
classical.bivariate.model_no.COVy <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                            + at(trait, 1):F.AGE.CLASS
                                            + at(trait, 2):F.AGE.CLASS,
                                            random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                            + idh(trait, init = c(40, 0.5)):YEAR
                                            + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                            + corgh(trait):vm(female.animal, ainv),
                                            residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                            data = phenotypic.data,
                                            workspace = 32e6,  # Boost memory
                                            maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.COVy$loglik - classical.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Year effect on laydate
fix.LD.year <- c(1e-8, 0.5)
names(fix.LD.year) <- c("F", "P")
classical.bivariate.model_no.LD.Vy <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                             + idh(trait, init = fix.LD.year):YEAR
                                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                             + corgh(trait):vm(female.animal, ainv),
                                             residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.LD.Vy$loglik - classical.bivariate.model_no.COVy$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Year effect on clutch size
fix.CS.year <- c(40, 1e-8)
names(fix.CS.year) <- c("P", "F")
classical.bivariate.model_no.CS.Vy <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                             + idh(trait, init = fix.CS.year):YEAR
                                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                             + corgh(trait):vm(female.animal, ainv),
                                             residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.CS.Vy$loglik - classical.bivariate.model_no.COVy$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value


# PERMANENT ENVIRONMENT
# covariance
classical.bivariate.model_no.COVpe <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                             + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                             + idh(trait, init = c(4, 0.5)):female.ID 
                                             + corgh(trait):vm(female.animal, ainv),
                                             residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.COVpe$loglik - classical.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Permanent environment effect on laydate
fix.LD.PE <- c(1e-8, 0.5)
names(fix.LD.PE) <- c("F", "P")
classical.bivariate.model_no.LD.Vpe <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS,
                                              random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                              + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                              + idh(trait, init = fix.LD.PE):female.ID 
                                              + corgh(trait):vm(female.animal, ainv),
                                              residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                              data = phenotypic.data,
                                              workspace = 32e6,  # Boost memory
                                              maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.LD.Vpe$loglik - classical.bivariate.model_no.COVpe$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Permanent environment effect on clutch size
fix.CS.PE <- c(4, 1e-8)
names(fix.CS.PE) <- c("P", "F")
classical.bivariate.model_no.CS.Vpe <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS,
                                              random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                              + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                              + idh(trait, init = fix.CS.PE):female.ID 
                                              + corgh(trait):vm(female.animal, ainv),
                                              residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                              data = phenotypic.data,
                                              workspace = 32e6,  # Boost memory
                                              maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.CS.Vpe$loglik - classical.bivariate.model_no.COVpe$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# ADDITIVE GENETIC
classical.bivariate.model_no.COVa <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                            + at(trait, 1):F.AGE.CLASS
                                            + at(trait, 2):F.AGE.CLASS,
                                            random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                            + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                            + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                            + idh(trait):vm(female.animal, ainv),
                                            residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                            data = phenotypic.data,
                                            workspace = 32e6,  # Boost memory
                                            maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.COVa$loglik - classical.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Additive genetic effect on laydate
fix.LD.a <- c(1e-8, 0.6)
names(fix.LD.a) <- c("F", "P")
classical.bivariate.model_no.LD.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                             + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                             + idh(trait, init = fix.LD.a):vm(female.animal, ainv),
                                             residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.LD.Va$loglik - classical.bivariate.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Additive genetic effect on clutch size
fix.CS.a <- c(8, 1e-8)
names(fix.CS.a) <- c("P", "F")
classical.bivariate.model_no.CS.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 2, 0.05)):round
                                             + corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID 
                                             + idh(trait, init = fix.CS.a):vm(female.animal, ainv),
                                             residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20)
chi.sq <- -2*(classical.bivariate.model_no.CS.Va$loglik - classical.bivariate.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value



##############################
##############################
## EXTENDED BIVARIATE MODEL ##
##############################
##############################

extended.bivariate.model <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                   + at(trait, 1):F.AGE.CLASS
                                   + at(trait, 2):F.AGE.CLASS
                                   + at(trait, 1):M.AGE.CLASS
                                   + at(trait, 2):M.AGE.CLASS,
                                   random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                   + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                   + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                   + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                   + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                         ~corgh(4):vm(female.animal, ainv)),
                                   residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                   data = phenotypic.data, 
                                   workspace = 6e8,  # Boost memory
                                   maxit = 20)
extended.bivariate.model <- update(extended.bivariate.model)
summary(extended.bivariate.model)$varcomp[,c("component", "std.error", "bound")]
summary(extended.bivariate.model, coef=T)$coef.fixed
wald.asreml(extended.bivariate.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

ext.LD.Vp       <- vpredict(extended.bivariate.model, LD.Vphenotypic ~ V5 + V8 + V11 + V19 + V21 + V25)  # Year excluded from phenotypic variance
ext.LD.y2       <- vpredict(extended.bivariate.model, LD.y2 ~ V2 / (V5+V8+V11+V19+V21+V25))
ext.LD.plot2    <- vpredict(extended.bivariate.model, LD.plot2 ~ V5 / (V5+V8+V11+V19+V21+V25))
ext.LD.f.pe2    <- vpredict(extended.bivariate.model, LD.f.pe2 ~ V8 / (V5+V8+V11+V19+V21+V25))
ext.LD.m.pe2    <- vpredict(extended.bivariate.model, LD.m.pe2 ~ V11/ (V5+V8+V11+V19+V21+V25))
ext.LD.f.h2     <- vpredict(extended.bivariate.model, LD.f.h2 ~ V19 / (V5+V8+V11+V19+V21+V25))
ext.LD.m.h2     <- vpredict(extended.bivariate.model, LD.m.h2 ~ V21 / (V5+V8+V11+V19+V21+V25))
ext.LD.cov.a    <- vpredict(extended.bivariate.model, LD.cov.a ~ V14 * sqrt(V19*V21))
ext.LD.tbv      <- vpredict(extended.bivariate.model, LD.tbv ~ V19+2*(V14*sqrt(V19*V21))+V21)
ext.LD.total.h2 <- vpredict(extended.bivariate.model, LD.total.h2 ~ (V19+2*(V14*sqrt(V19*V21))+V21) / (V5+V8+V11+V19+V21+V25))
ext.LD.r2       <- vpredict(extended.bivariate.model, LD.r2 ~ V25 / (V5+V8+V11+V19+V21+V25))
ext.LD.f.i      <- vpredict(extended.bivariate.model, LD.i ~ V8+V19)
ext.LD.m.i      <- vpredict(extended.bivariate.model, LD.i ~ V11+V21)
ext.LD.f.i2     <- vpredict(extended.bivariate.model, LD.f.i2 ~ (V8+V19) / (V5+V8+V11+V19+V21+V25))
ext.LD.m.i2     <- vpredict(extended.bivariate.model, LD.m.i2 ~ (V11+V21) / (V5+V8+V11+V19+V21+V25))

ext.CS.Vp       <- vpredict(extended.bivariate.model, CS.Vphenotypic ~ V6+V9+V12+V20+V22+V26)
ext.CS.y2       <- vpredict(extended.bivariate.model, CS.y2 ~ V3 / (V6+V9+V12+V20+V22+V26))
ext.CS.plot2    <- vpredict(extended.bivariate.model, CS.plot2 ~ V6 / (V6+V9+V12+V20+V22+V26))
ext.CS.f.pe2    <- vpredict(extended.bivariate.model, CS.f.pe2 ~ V9 / (V6+V9+V12+V20+V22+V26))
ext.CS.m.pe2    <- vpredict(extended.bivariate.model, CS.m.pe2 ~ V12 / (V6+V9+V12+V20+V22+V26))
ext.CS.f.h2     <- vpredict(extended.bivariate.model, CS.f.h2 ~ V20 / (V6+V9+V12+V20+V22+V26))
ext.CS.m.h2     <- vpredict(extended.bivariate.model, CS.m.h2 ~ V22 / (V6+V9+V12+V20+V22+V26))
ext.CS.cov.a    <- vpredict(extended.bivariate.model, CS.cov.a ~ V17 * sqrt(V20*V22))
ext.CS.tbv      <- vpredict(extended.bivariate.model, CS.tbv ~ V20+2*(V17*sqrt(V20*V22))+V22)
ext.CS.total.h2 <- vpredict(extended.bivariate.model, CS.h2 ~ (V20+2*(V17*sqrt(V20*V22))+V22) / (V6+V9+V12+V20+V22+V26))
ext.CS.r2       <- vpredict(extended.bivariate.model, CS.r2 ~ V26 / (V6+V9+V12+V15+V22+V26))
ext.CS.f.i      <- vpredict(extended.bivariate.model, CS.f.i ~ V9+V20)
ext.CS.m.i      <- vpredict(extended.bivariate.model, CS.m.i ~ V12+V22)
ext.CS.f.i2     <- vpredict(extended.bivariate.model, CS.f.i2 ~ (V9+V20)/(V6+V9+V12+V20+V22+V26))
ext.CS.m.i2     <- vpredict(extended.bivariate.model, CS.m.i2 ~ (V12+V22)/(V6+V9+V12+V20+V22+V26))

ext.COV.plot  <- vpredict(extended.bivariate.model, ext.COV.plot ~ V4*sqrt(V5*V6))
ext.COV.y     <- vpredict(extended.bivariate.model, ext.COV.y ~ V1*sqrt(V2*V3))
ext.COV.fpe   <- vpredict(extended.bivariate.model, ext.COV.fpe ~ V7*sqrt(V8*V9))
ext.COV.mpe   <- vpredict(extended.bivariate.model, ext.COV.mpe ~ V10*sqrt(V11*V12))
ext.COV.a.fge <- vpredict(extended.bivariate.model, ext.COV.a.fge ~ V13*sqrt(V19*V20))  # Additive genetic covariance of female effects
ext.COV.a.mge <- vpredict(extended.bivariate.model, ext.COV.a.mge ~ V18*sqrt(V21*V22))  # Additive genetic covariance of male effects
ext.COV.a.fld_mcs <- vpredict(extended.bivariate.model, ext.COV.a.fld_mcs ~ V16*sqrt(V19*V22))  # Cross-sex cross-trait genetic covariance
ext.COV.a.fcs_mld <- vpredict(extended.bivariate.model, ext.COV.a.fcs_mld ~ V15*sqrt(V20*V21))
ext.COV.r     <- vpredict(extended.bivariate.model, ext.COV.r ~ V24*sqrt(V25*V26))

LD.m.h2_including.year <- vpredict(extended.bivariate.model, LD.m.h2_including.year ~ V21 / (V2+V5+V8+V11+V19+V21+V25))

ext.LD.Vp
ext.LD.y2
ext.LD.plot2
ext.LD.f.pe2
ext.LD.m.pe2
ext.LD.f.h2
ext.LD.m.h2
ext.LD.cov.a
ext.LD.tbv
ext.LD.total.h2
ext.LD.r2
ext.LD.f.i
ext.LD.m.i
ext.LD.f.i2
ext.LD.m.i2

ext.CS.Vp
ext.CS.y2
ext.CS.plot2
ext.CS.f.pe2
ext.CS.m.pe2
ext.CS.f.h2
ext.CS.m.h2
ext.CS.cov.a
ext.CS.tbv
ext.CS.total.h2
ext.CS.r2
ext.CS.f.i
ext.CS.m.i
ext.CS.f.i2
ext.CS.m.i2

ext.COV.plot
ext.COV.y
ext.COV.fpe
ext.COV.mpe
ext.COV.a.fge
ext.COV.a.mge
ext.COV.a.fld_mcs
ext.COV.a.fcs_mld
ext.COV.r

LD.m.h2_including.year

###############################
# TESTING VARIANCE COMPONENTS #
###############################

# PLOT
# covariance
extended.bivariate.model_no.plot.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                               + at(trait, 1):F.AGE.CLASS
                                               + at(trait, 2):F.AGE.CLASS
                                               + at(trait, 1):M.AGE.CLASS
                                               + at(trait, 2):M.AGE.CLASS,
                                               random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                               + idh(trait, init = c(2, 0.05)):round 
                                               + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                               + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                               + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                     ~corgh(4):vm(female.animal, ainv)),
                                               residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                               data = phenotypic.data, 
                                               workspace = 6e8,  # Boost memory
                                               maxit = 20)
extended.bivariate.model_no.plot.COV <- update(extended.bivariate.model_no.plot.COV)
chi.sq <- -2*(extended.bivariate.model_no.plot.COV$loglik - extended.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Plot effect on laydate
fix.LD.plot <- c(1e-8, 0.05)
names(fix.LD.plot) <- c("F", "P")
extended.bivariate.model_no.LD.plot <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                              + idh(trait, init = fix.LD.plot):round 
                                              + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                              + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 20)
extended.bivariate.model_no.LD.plot <- update(extended.bivariate.model_no.LD.plot)
chi.sq <- -2*(extended.bivariate.model_no.LD.plot$loglik - extended.bivariate.model_no.plot.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Plot effect on clutch size
fix.CS.plot <- c(2, 1e-8)
names(fix.CS.plot) <- c("P", "F")
extended.bivariate.model_no.CS.plot <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                              + idh(trait, init = fix.CS.plot):round 
                                              + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                              + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 20)
extended.bivariate.model_no.CS.plot <- update(extended.bivariate.model_no.CS.plot)
chi.sq <- -2*(extended.bivariate.model_no.CS.plot$loglik - extended.bivariate.model_no.plot.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# YEAR
# Year effect on laydate
# Not fitting year covariance
extended.bivariate.model_no.year.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                               + at(trait, 1):F.AGE.CLASS
                                               + at(trait, 2):F.AGE.CLASS
                                               + at(trait, 1):M.AGE.CLASS
                                               + at(trait, 2):M.AGE.CLASS,
                                               random = ~ idh(trait, init = c(40, 0.5)):YEAR 
                                               + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                               + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                               + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                               + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                     ~corgh(4):vm(female.animal, ainv)),
                                               residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                               data = phenotypic.data, 
                                               workspace = 6e8,  # Boost memory
                                               maxit = 30)
extended.bivariate.model_no.year.COV <- update(extended.bivariate.model_no.year.COV)
chi.sq <- -2*(extended.bivariate.model_no.year.COV$loglik - extended.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Laydate
fix.LD.year <- c(1e-8, 0.5)
names(fix.LD.year) <- c("F", "P")
extended.bivariate.model_no.LD.year <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ idh(trait, init = fix.LD.year):YEAR 
                                              + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                              + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                              + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 50)
extended.bivariate.model_no.LD.year <- update(extended.bivariate.model_no.LD.year)
chi.sq <- -2*(extended.bivariate.model_no.LD.year$loglik - extended.bivariate.model_no.year.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# YEAR EFFECT ON CLUTCH SIZE
fix.CS.year <- c(40, 1e-8)
names(fix.CS.year) <- c("P", "F")
extended.bivariate.model_no.CS.year <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ idh(trait, init = fix.CS.year):YEAR 
                                              + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                              + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                              + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 50)
extended.bivariate.model_no.CS.year <- update(extended.bivariate.model_no.CS.year)
chi.sq <- -2*(extended.bivariate.model_no.CS.year$loglik - extended.bivariate.model_no.year.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# FEMALE PERMANENT ENVIRONMENT EFFECT
# Omitting female permanent environmental covariance
extended.bivariate.model_no.f.PE.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                               + at(trait, 1):F.AGE.CLASS
                                               + at(trait, 2):F.AGE.CLASS
                                               + at(trait, 1):M.AGE.CLASS
                                               + at(trait, 2):M.AGE.CLASS,
                                               random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                               + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                               + idh(trait, init = c(4, 0.5)):female.ID
                                               + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                               + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                     ~corgh(4):vm(female.animal, ainv)),
                                               residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                               data = phenotypic.data, 
                                               workspace = 6e8,  # Boost memory
                                               maxit = 20)
extended.bivariate.model_no.f.PE.COV <- update(extended.bivariate.model_no.f.PE.COV)
chi.sq <- -2*(extended.bivariate.model_no.f.PE.COV$loglik - extended.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Female permanent environment effect on laydate
fix.LD.fPE <- c(1e-8, 0.5)
names(fix.LD.fPE) <- c("F", "U")
extended.bivariate.model_no.LD.fPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                             + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                             + idh(trait, init = fix.LD.fPE):female.ID
                                             + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv)),
                                             residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                             data = phenotypic.data, 
                                             workspace = 6e8,  # Boost memory
                                             maxit = 20)
extended.bivariate.model_no.LD.fPE <- update(extended.bivariate.model_no.LD.fPE)
chi.sq <- -2*(extended.bivariate.model_no.LD.fPE$loglik - extended.bivariate.model_no.f.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment effect on clutch size
fix.CS.fPE <- c(2, 1e-8)
names(fix.CS.fPE) <- c("P", "F")
extended.bivariate.model_no.CS.fPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                             + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                             + idh(trait, init = fix.CS.fPE):female.ID
                                             + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv)),
                                             residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                             data = phenotypic.data, 
                                             workspace = 6e8,  # Boost memory
                                             maxit = 20)
extended.bivariate.model_no.CS.fPE <- update(extended.bivariate.model_no.CS.fPE)
chi.sq <- -2*(extended.bivariate.model_no.CS.fPE$loglik - extended.bivariate.model_no.f.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# MALE PERMANENT ENVIRONMENT EFFECT
# Not fitting male permanent environmental covariance
extended.bivariate.model_no.m.PE.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                               + at(trait, 1):F.AGE.CLASS
                                               + at(trait, 2):F.AGE.CLASS
                                               + at(trait, 1):M.AGE.CLASS
                                               + at(trait, 2):M.AGE.CLASS,
                                               random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                               + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                               + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                               + idh(trait, init = c(2, 0.05)):male.ID
                                               + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                     ~corgh(4):vm(female.animal, ainv)),
                                               residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                               data = phenotypic.data, 
                                               workspace = 6e8,  # Boost memory
                                               maxit = 20)
extended.bivariate.model_no.m.PE.COV <- update(extended.bivariate.model_no.m.PE.COV)
extended.bivariate.model_no.m.PE.COV <- update(extended.bivariate.model_no.m.PE.COV)
chi.sq <- -2*(extended.bivariate.model_no.m.PE.COV$loglik - extended.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Male permanent environment effect on laydate
fix.LD.mPE <- c(1e-8, 0.04)
names(fix.LD.mPE) <- c("F", "P")
extended.bivariate.model_no.LD.mPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                             + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                             + idh(trait, init = fix.LD.mPE):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv)),
                                             residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                             data = phenotypic.data, 
                                             workspace = 6e8,  # Boost memory
                                             maxit = 20)
extended.bivariate.model_no.LD.mPE <- update(extended.bivariate.model_no.LD.mPE)
extended.bivariate.model_no.LD.mPE <- update(extended.bivariate.model_no.LD.mPE)
chi.sq <- -2*(extended.bivariate.model_no.LD.mPE$loglik - extended.bivariate.model_no.m.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male permanent environment effect on clutch size
fix.CS.mPE <- c(2, 1e-8)
names(fix.CS.mPE) <- c("P", "F")
extended.bivariate.model_no.CS.mPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                             + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                             + idh(trait, init = fix.CS.mPE):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv)),
                                             residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                             data = phenotypic.data, 
                                             workspace = 6e8,  # Boost memory
                                             maxit = 20)
extended.bivariate.model_no.CS.mPE <- update(extended.bivariate.model_no.CS.mPE)
chi.sq <- -2*(extended.bivariate.model_no.CS.mPE$loglik - extended.bivariate.model_no.m.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value


# CONSTRANING COMPONENTS OF THE G-MATRIX
# Extract starting parameters
gamma.extraction.model <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                 + at(trait, 1):F.AGE.CLASS
                                 + at(trait, 2):F.AGE.CLASS
                                 + at(trait, 1):M.AGE.CLASS
                                 + at(trait, 2):M.AGE.CLASS,
                                 random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                 + corgh(trait, init = c(-0.1, 2, 0.05)):round 
                                 + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                 + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                 + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                       ~corgh(4):vm(female.animal, ainv)),
                                 residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                 data = phenotypic.data, 
                                 workspace = 6e8,
                                 start.values = "gammas.txt") 

# Female additive genetic effect on laydate
# Fixing covariances with female additive genetic variance in laydate to zero
gam.LD.f.COVs <- gamma.extraction.model$vparameters.table 
gam.LD.f.COVs[c(13,14,16),]$Value <- 0
gam.LD.f.COVs[c(13,14,16),]$Constraint <- "F"
extended.bivariate.model_no.LD.f.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS
                                                + at(trait, 1):M.AGE.CLASS
                                                + at(trait, 2):M.AGE.CLASS,
                                                random = ~ corgh(trait):YEAR 
                                                + corgh(trait):round 
                                                + corgh(trait):female.ID
                                                + corgh(trait):male.ID
                                                + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                      ~corgh(4):vm(female.animal, ainv)),
                                                residual = ~units:corgh(trait),
                                                data = phenotypic.data, 
                                                workspace = 6e8,  # Boost memory
                                                maxit = 30,
                                                G.param = gam.LD.f.COVs)
extended.bivariate.model_no.LD.f.COVs <- update(extended.bivariate.model_no.LD.f.COVs)
chi.sq <- -2*(extended.bivariate.model_no.LD.f.Va$loglik - extended.bivariate.model_no.LD.f.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Additionally fixing female additive genetic variance in laydate close to zero
gam.LD.f.Va <- gam.LD.f.COVs 
gam.LD.f.Va[19,]$Value <- 1e-8
gam.LD.f.Va[19,]$Constraint <- "F"
extended.bivariate.model_no.LD.f.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait):YEAR 
                                              + corgh(trait):round 
                                              + corgh(trait):female.ID
                                              + corgh(trait):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 30,
                                              G.param = gam.LD.f.Va)
extended.bivariate.model_no.LD.f.Va <- update(extended.bivariate.model_no.LD.f.Va)
chi.sq <- -2*(extended.bivariate.model_no.LD.f.Va$loglik - extended.bivariate.model_no.LD.f.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female additive genetic effect on clutch size
# Fixing covariances with female additive genetic variance in clutch size to zero
gam.CS.f.COVs <- gamma.extraction.model$vparameters.table 
gam.CS.f.COVs[c(13,15,17),]$Value <- 0
gam.CS.f.COVs[c(13,15,17),]$Constraint <- "F"
extended.bivariate.model_no.CS.f.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS
                                                + at(trait, 1):M.AGE.CLASS
                                                + at(trait, 2):M.AGE.CLASS,
                                                random = ~ corgh(trait):YEAR 
                                                + corgh(trait):round 
                                                + corgh(trait):female.ID
                                                + corgh(trait):male.ID
                                                + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                      ~corgh(4):vm(female.animal, ainv)),
                                                residual = ~units:corgh(trait),
                                                data = phenotypic.data, 
                                                workspace = 6e8,  # Boost memory
                                                maxit = 30,
                                                G.param = gam.CS.f.COVs)
extended.bivariate.model_no.CS.f.COVs <- update(extended.bivariate.model_no.CS.f.COVs)

# Additionally fixing female additive genetic variance of clutch size close to zero
gam.CS.f.Va <- gam.CS.f.COVs 
gam.CS.f.Va[20,]$Value <- 1e-8
gam.CS.f.Va[20,]$Constraint <- "F"
extended.bivariate.model_no.CS.f.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait):YEAR 
                                              + corgh(trait):round 
                                              + corgh(trait):female.ID
                                              + corgh(trait):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 30,
                                              G.param = gam.CS.f.Va)
extended.bivariate.model_no.CS.f.Va <- update(extended.bivariate.model_no.CS.f.Va)
chi.sq <- -2*(extended.bivariate.model_no.CS.f.Va$loglik - extended.bivariate.model_no.CS.f.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic effect on laydate
# Fixing covariances with male additive genetic variance in laydate to zero
gam.LD.m.COVs <- gamma.extraction.model$vparameters.table 
gam.LD.m.COVs[c(14,15,18),]$Value <- 0
gam.LD.m.COVs[c(14,15,18),]$Constraint <- "F"
extended.bivariate.model_no.LD.m.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS
                                                + at(trait, 1):M.AGE.CLASS
                                                + at(trait, 2):M.AGE.CLASS,
                                                random = ~ corgh(trait):YEAR 
                                                + corgh(trait):round 
                                                + corgh(trait):female.ID
                                                + corgh(trait):male.ID
                                                + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                      ~corgh(4):vm(female.animal, ainv)),
                                                residual = ~units:corgh(trait),
                                                data = phenotypic.data, 
                                                workspace = 6e8,  # Boost memory
                                                maxit = 50,
                                                G.param = gam.LD.m.COVs)
extended.bivariate.model_no.LD.m.COVs <- update(extended.bivariate.model_no.LD.m.COVs)

# Additionally fixing male additive genetic variance in laydate close to zero
gam.LD.m.Va <- gam.LD.m.COVs 
gam.LD.m.Va[21,]$Value <- 1e-8
gam.LD.m.Va[21,]$Constraint <- "F"
extended.bivariate.model_no.LD.m.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait):YEAR 
                                              + corgh(trait):round 
                                              + corgh(trait):female.ID
                                              + corgh(trait):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 50,
                                              G.param = gam.LD.m.Va)
extended.bivariate.model_no.LD.m.Va <- update(extended.bivariate.model_no.LD.m.Va)
chi.sq <- -2*(extended.bivariate.model_no.LD.m.Va$loglik - extended.bivariate.model_no.LD.m.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic effect on clutch size
# Fixing covariances with male additive genetic variance in clutch size to zero
gam.CS.m.COVs <- gamma.extraction.model$vparameters.table 
gam.CS.m.COVs[c(16,17,18),]$Value <- 0
gam.CS.m.COVs[c(16,17,18),]$Constraint <- "F"
extended.bivariate.model_no.CS.m.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS
                                                + at(trait, 1):M.AGE.CLASS
                                                + at(trait, 2):M.AGE.CLASS,
                                                random = ~ corgh(trait):YEAR 
                                                + corgh(trait):round 
                                                + corgh(trait):female.ID
                                                + corgh(trait):male.ID
                                                + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                      ~corgh(4):vm(female.animal, ainv)),
                                                residual = ~units:corgh(trait),
                                                data = phenotypic.data, 
                                                workspace = 6e8,  # Boost memory
                                                maxit = 50,
                                                G.param = gam.CS.m.COVs)
extended.bivariate.model_no.CS.m.COVs <- update(extended.bivariate.model_no.CS.m.COVs)

# Additionally fixing male additive genetic variance in laydate close to zero
gam.CS.m.Va <- gam.CS.m.COVs 
gam.CS.m.Va[22,]$Value <- 1e-8
gam.CS.m.Va[22,]$Constraint <- "F"
extended.bivariate.model_no.CS.m.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait):YEAR 
                                              + corgh(trait):round 
                                              + corgh(trait):female.ID
                                              + corgh(trait):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 50,
                                              G.param = gam.CS.m.Va)
extended.bivariate.model_no.CS.m.Va <- update(extended.bivariate.model_no.CS.m.Va)
chi.sq <- -2*(extended.bivariate.model_no.CS.m.Va$loglik - extended.bivariate.model_no.CS.m.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Covariance between female additive genetic effects on laydate and clutch size
gam.f.COVa <- gamma.extraction.model$vparameters.table 
gam.f.COVa[13,]$Value <- 0
gam.f.COVa[13,]$Constraint <- "F"
extended.bivariate.model_no.f.COVa <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait):YEAR 
                                             + corgh(trait):round 
                                             + corgh(trait):female.ID
                                             + corgh(trait):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv)),
                                             residual = ~units:corgh(trait),
                                             data = phenotypic.data, 
                                             workspace = 6e8,  # Boost memory
                                             maxit = 50,
                                             G.param = gam.f.COVa)
extended.bivariate.model_no.f.COVa <- update(extended.bivariate.model_no.f.COVa)
chi.sq <- -2*(extended.bivariate.model_no.f.COVa$loglik - extended.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Covariance between male additive genetic effects on laydate and clutch size
gam.m.COVa <- gamma.extraction.model$vparameters.table 
gam.m.COVa[18,]$Value <- 0
gam.m.COVa[18,]$Constraint <- "F"
extended.bivariate.model_no.COVa.mge <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                               + at(trait, 1):F.AGE.CLASS
                                               + at(trait, 2):F.AGE.CLASS
                                               + at(trait, 1):M.AGE.CLASS
                                               + at(trait, 2):M.AGE.CLASS,
                                               random = ~ corgh(trait):YEAR 
                                               + corgh(trait):round 
                                               + corgh(trait):female.ID
                                               + corgh(trait):male.ID
                                               + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                     ~corgh(4):vm(female.animal, ainv)),
                                               residual = ~units:corgh(trait),
                                               data = phenotypic.data, 
                                               workspace = 6e8,  # Boost memory
                                               maxit = 50,
                                               G.param = gam.m.COVa)
extended.bivariate.model_no.COVa.mge <- update(extended.bivariate.model_no.COVa.mge)
chi.sq <- -2*(extended.bivariate.model_no.COVa.mge$loglik - extended.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Covariance between female and male additive genetic effects on laydate
gam.COVa.LD <- gamma.extraction.model$vparameters.table 
gam.COVa.LD[14,]$Value <- 0
gam.COVa.LD[14,]$Constraint <- "F"
extended.bivariate.model_no.COVa.LD <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait):YEAR 
                                              + corgh(trait):round 
                                              + corgh(trait):female.ID
                                              + corgh(trait):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 50,
                                              G.param = gam.COVa.LD)
extended.bivariate.model_no.COVa.LD <- update(extended.bivariate.model_no.COVa.LD)
chi.sq <- -2*(extended.bivariate.model_no.COVa.LD$loglik - extended.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Covariance between female and male additive genetic effects on clutch size
gam.COVa.CS <- gamma.extraction.model$vparameters.table 
gam.COVa.CS[17,]$Value <- 0
gam.COVa.CS[17,]$Constraint <- "F"
extended.bivariate.model_no.COVa.CS <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                              + at(trait, 1):F.AGE.CLASS
                                              + at(trait, 2):F.AGE.CLASS
                                              + at(trait, 1):M.AGE.CLASS
                                              + at(trait, 2):M.AGE.CLASS,
                                              random = ~ corgh(trait):YEAR 
                                              + corgh(trait):round 
                                              + corgh(trait):female.ID
                                              + corgh(trait):male.ID
                                              + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                    ~corgh(4):vm(female.animal, ainv)),
                                              residual = ~units:corgh(trait),
                                              data = phenotypic.data, 
                                              workspace = 6e8,  # Boost memory
                                              maxit = 50,
                                              G.param = gam.COVa.CS)
extended.bivariate.model_no.COVa.CS <- update(extended.bivariate.model_no.COVa.CS)
chi.sq <- -2*(extended.bivariate.model_no.COVa.CS$loglik - extended.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Additive genetic covariance between female laydate and male clutch size
gam.COVa.fLD_mCS <- gamma.extraction.model$vparameters.table 
gam.COVa.fLD_mCS[16,]$Value <- 0
gam.COVa.fLD_mCS[16,]$Constraint <- "F"
extended.bivariate.model_no.COVa.fLD_mCS <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                   + at(trait, 1):F.AGE.CLASS
                                                   + at(trait, 2):F.AGE.CLASS
                                                   + at(trait, 1):M.AGE.CLASS
                                                   + at(trait, 2):M.AGE.CLASS,
                                                   random = ~ corgh(trait):YEAR 
                                                   + corgh(trait):round 
                                                   + corgh(trait):female.ID
                                                   + corgh(trait):male.ID
                                                   + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                         ~corgh(4):vm(female.animal, ainv)),
                                                   residual = ~units:corgh(trait),
                                                   data = phenotypic.data, 
                                                   workspace = 6e8,  # Boost memory
                                                   maxit = 50,
                                                   G.param = gam.COVa.fLD_mCS)
extended.bivariate.model_no.COVa.fLD_mCS <- update(extended.bivariate.model_no.COVa.fLD_mCS)
chi.sq <- -2*(extended.bivariate.model_no.COVa.fLD_mCS$loglik - extended.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Additive genetic covariance between female clutch size and male laydate
gam.COVa.fCS_mLD <- gamma.extraction.model$vparameters.table 
gam.COVa.fCS_mLD[15,]$Value <- 0
gam.COVa.fCS_mLD[15,]$Constraint <- "F"
extended.bivariate.model_no.COVa.fCS_mLD <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                   + at(trait, 1):F.AGE.CLASS
                                                   + at(trait, 2):F.AGE.CLASS
                                                   + at(trait, 1):M.AGE.CLASS
                                                   + at(trait, 2):M.AGE.CLASS,
                                                   random = ~ corgh(trait):YEAR 
                                                   + corgh(trait):round 
                                                   + corgh(trait):female.ID
                                                   + corgh(trait):male.ID
                                                   + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                         ~corgh(4):vm(female.animal, ainv)),
                                                   residual = ~units:corgh(trait),
                                                   data = phenotypic.data, 
                                                   workspace = 6e8,  # Boost memory
                                                   maxit = 50,
                                                   G.param = gam.COVa.fCS_mLD)
extended.bivariate.model_no.COVa.fCS_mLD <- update(extended.bivariate.model_no.COVa.fCS_mLD)
chi.sq <- -2*(extended.bivariate.model_no.COVa.fCS_mLD$loglik - extended.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# TESTING ENTIRE LAYDATE G-MATRIX
# Fixing all laydate-clutch size covariances to zero
gam.LDCS.COVs <- gamma.extraction.model$vparameters.table
gam.LDCS.COVs[c(13,15,16,18),]$Value <- 0  # Fix relevant covariances to zero
gam.LDCS.COVs[c(13,15,16,18),]$Constraint <- "F"

extended.bivariate.model_no.LDCS.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                                + at(trait, 1):F.AGE.CLASS
                                                + at(trait, 2):F.AGE.CLASS
                                                + at(trait, 1):M.AGE.CLASS
                                                + at(trait, 2):M.AGE.CLASS,
                                                random = ~ corgh(trait):YEAR 
                                                + corgh(trait):round 
                                                + corgh(trait):female.ID
                                                + corgh(trait):male.ID
                                                + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                      ~corgh(4):vm(female.animal, ainv)),
                                                residual = ~units:corgh(trait),
                                                data = phenotypic.data, 
                                                workspace = 6e8,  # Boost memory
                                                maxit = 50,
                                                G.param = gam.LDCS.COVs)
extended.bivariate.model_no.LDCS.COVs <- update(extended.bivariate.model_no.LDCS.COVs)

# Additionally fixing laydate (co)variances to zero
gam.LD.G <- gam.LDCS.COVs 
gam.LD.G[c(19,21),]$Value <- 1e-8  # Fix relevant variances to very small value
gam.LD.G[c(19,21),]$Constraint <- "F"
gam.LD.G[14,]$Value <- 0  # Fix male-female additive genetic covariance for laydate to zero
gam.LD.G[14,]$Constraint <- "F"
extended.bivariate.model_no.LD.G <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                           + at(trait, 1):F.AGE.CLASS
                                           + at(trait, 2):F.AGE.CLASS
                                           + at(trait, 1):M.AGE.CLASS
                                           + at(trait, 2):M.AGE.CLASS,
                                           random = ~ corgh(trait):YEAR 
                                           + corgh(trait):round 
                                           + corgh(trait):female.ID
                                           + corgh(trait):male.ID
                                           + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                 ~corgh(4):vm(female.animal, ainv)),
                                           residual = ~units:corgh(trait),
                                           data = phenotypic.data, 
                                           workspace = 6e8,  # Boost memory
                                           maxit = 50,
                                           G.param = gam.LD.G)
extended.bivariate.model_no.LD.G <- update(extended.bivariate.model_no.LD.G)
chi.sq <- -2*(extended.bivariate.model_no.LD.G$loglik - extended.bivariate.model_no.LDCS.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value

# Testing entire clutch size G-matrix
gam.CS.G <- gam.LDCS.COVs 
gam.CS.G[c(20,20),]$Value <- 1e-8  # Fix relevant variances to very small value
gam.CS.G[c(20,22),]$Constraint <- "F"
gam.CS.G[17,]$Value <- 0  # Fix intersexual genetic covariance for clutch size to zero
gam.CS.G[17,]$Constraint <- "F"
extended.bivariate.model_no.CS.G <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                           + at(trait, 1):F.AGE.CLASS
                                           + at(trait, 2):F.AGE.CLASS
                                           + at(trait, 1):M.AGE.CLASS
                                           + at(trait, 2):M.AGE.CLASS,
                                           random = ~ corgh(trait):YEAR 
                                           + corgh(trait):round 
                                           + corgh(trait):female.ID
                                           + corgh(trait):male.ID
                                           + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                 ~corgh(4):vm(female.animal, ainv)),
                                           residual = ~units:corgh(trait),
                                           data = phenotypic.data, 
                                           workspace = 6e8,  # Boost memory
                                           maxit = 50,
                                           G.param = gam.CS.G)
extended.bivariate.model_no.CS.G <- update(extended.bivariate.model_no.CS.G)
chi.sq <- -2*(extended.bivariate.model_no.CS.G$loglik - extended.bivariate.model_no.LDCS.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value



###################################
###################################
## SPATIAL AUTOCORRELATION MODEL ##
###################################
###################################

phenotypic.data$col <- as.factor(floor(phenotypic.data$x/50))
phenotypic.data$row <- as.factor(floor(phenotypic.data$y/50))
length(table(phenotypic.data$row))
length(table(phenotypic.data$col))

SA.bivariate.model <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                             + at(trait, 1):F.AGE.CLASS
                             + at(trait, 2):F.AGE.CLASS
                             + at(trait, 1):M.AGE.CLASS
                             + at(trait, 2):M.AGE.CLASS,
                             random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                             + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                             + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                             + str(~corgh(trait):vm(female.animal, ainv) 
                                   + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                             + trait:ar1v(col):ar1(row),
                             residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                             data = phenotypic.data,
                             workspace = 8e8,  # Boost memory
                             maxit = 100)
SA.bivariate.model <- update(SA.bivariate.model)

summary(SA.bivariate.model)$varcomp[c("component", "std.error", "bound")]
summary(SA.bivariate.model,  coef=T)$coef.fixed
wald.asreml(SA.bivariate.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

SA.LD.Vp       <- vpredict(SA.bivariate.model, SA.LD.Vphenotypic ~ V5 + V8 + V16 + V18 + V21 + V25)
SA.LD.y2       <- vpredict(SA.bivariate.model, SA.LD.y2 ~ V2 / (V5+V8+V16+V18+V21+V25))
SA.LD.f.pe2    <- vpredict(SA.bivariate.model, SA.LD.f.pe2 ~ V5 / (V5+V8+V16+V18+V21+V25))
SA.LD.m.pe2    <- vpredict(SA.bivariate.model, SA.LD.m.pe2 ~ V8/ (V5+V8+V16+V18+V21+V25))
SA.LD.f.h2     <- vpredict(SA.bivariate.model, SA.LD.f.h2 ~ V16 / (V5+V8+V16+V18+V21+V25))
SA.LD.m.h2     <- vpredict(SA.bivariate.model, SA.LD.m.h2 ~ V18 / (V5+V8+V16+V18+V21+V25))
SA.LD.cov.a    <- vpredict(SA.bivariate.model, SA.LD.cov.a ~ V11*sqrt(V16*V18))
SA.LD.tbv      <- vpredict(SA.bivariate.model, SA.LD.tbv ~ V16+2*(V11*sqrt(V16*V18))+V18)
SA.LD.total.h2 <- vpredict(SA.bivariate.model, SA.LD.total.h2 ~ (V16+2*(V11*sqrt(V16*V18))+V18) / (V5+V8+V16+V18+V21+V25))
SA.LD.sa2       <- vpredict(SA.bivariate.model, SA.LD.SA ~ V21 / (V5+V8+V16+V18+V21+V25))
SA.LD.r2       <- vpredict(SA.bivariate.model, SA.LD.r2 ~ V25 / (V5+V8+V16+V18+V21+V25))
SA.LD.f.i      <- vpredict(SA.bivariate.model, SA.LD.i ~ V5+V16)
SA.LD.m.i      <- vpredict(SA.bivariate.model, SA.LD.i ~ V8+V18)
SA.LD.f.i2     <- vpredict(SA.bivariate.model, SA.LD.f.i2 ~ (V5+V16) / (V5+V8+V16+V18+V21+V25))
SA.LD.m.i2     <- vpredict(SA.bivariate.model, SA.LD.m.i2 ~ (V8+V18) / (V5+V8+V16+V18+V21+V25))

SA.CS.Vp       <- vpredict(SA.bivariate.model, SA.CS.Vphenotypic ~ V6+V9+V17+V19+V21+V26)
SA.CS.y2       <- vpredict(SA.bivariate.model, SA.CS.y2 ~ V3 / (V6+V9+V17+V19+V21+V26))
SA.CS.f.pe2    <- vpredict(SA.bivariate.model, SA.CS.f.pe2 ~ V6 / (V6+V9+V17+V19+V21+V26))
SA.CS.m.pe2    <- vpredict(SA.bivariate.model, SA.CS.m.pe2 ~ V9 / (V6+V9+V17+V19+V21+V26))
SA.CS.f.h2     <- vpredict(SA.bivariate.model, SA.CS.f.h2 ~ V17 / (V6+V9+V17+V19+V21+V26))
SA.CS.m.h2     <- vpredict(SA.bivariate.model, SA.CS.m.h2 ~ V19 / (V6+V9+V17+V19+V21+V26))
SA.CS.cov.a    <- vpredict(SA.bivariate.model, SA.CS.cov.a ~ V14*sqrt(V17*V19))
SA.CS.tbv      <- vpredict(SA.bivariate.model, SA.CS.tbv ~ V17+2*(V14*sqrt(V17*V19))+V19)
SA.CS.total.h2 <- vpredict(SA.bivariate.model, SA.CS.h2 ~ (V17+2*(V14*sqrt(V17*V19))+V19) / (V6+V9+V17+V19+V21+V26))
SA.CS.sa2       <- vpredict(SA.bivariate.model, SA.CS.SA ~ V21 / (V6+V9+V17+V19+V21+V26))
SA.CS.r2       <- vpredict(SA.bivariate.model, SA.CS.r2 ~ V26 / (V6+V9+V17+V19+V21+V26))
SA.CS.f.i      <- vpredict(SA.bivariate.model, SA.CS.f.i ~ V6+V17)
SA.CS.m.i      <- vpredict(SA.bivariate.model, SA.CS.m.i ~ V9+V19)
SA.CS.f.i2     <- vpredict(SA.bivariate.model, SA.CS.f.i2 ~ (V6+V17)/(V6+V9+V17+V19+V21+V26))
SA.CS.m.i2     <- vpredict(SA.bivariate.model, SA.CS.m.i2 ~ (V9+V19)/(V6+V9+V17+V19+V21+V26))

SA.cov.y     <- vpredict(SA.bivariate.model, SA.cov.y ~ V1*sqrt(V2*V3))
SA.cov.fpe   <- vpredict(SA.bivariate.model, SA.cov.fpe ~ V4*sqrt(V5*V6))
SA.cov.mpe   <- vpredict(SA.bivariate.model, SA.cov.mpe ~ V7*sqrt(V8*V9))
SA.cov.A.fge <- vpredict(SA.bivariate.model, SA.cov.A.fge ~ V10*sqrt(V16*V17))  # Additive genetic covariance of female effects
SA.cov.A.mge <- vpredict(SA.bivariate.model, SA.cov.A.mge ~ V15*sqrt(V18*V19))  # Additive genetic covariance of male effects
SA.COV.a.fld_mcs <- vpredict(SA.bivariate.model, SA.COV.a.fld_mcs ~ V13*sqrt(V16*V19))  # Cross-sex cross-trait genetic covariance
SA.COV.a.fcs_mld <- vpredict(SA.bivariate.model, SA.COV.a.fcs_mld ~ V12*sqrt(V17*V18))
SA.COV.r <- vpredict(SA.bivariate.model, SA.COV.r ~ V24*sqrt(V25*V26))

SA.LD.Vp
SA.LD.y2
SA.LD.f.pe2
SA.LD.m.pe2
SA.LD.f.h2
SA.LD.m.h2
SA.LD.cov.a
SA.LD.tbv
SA.LD.total.h2
SA.LD.sa2
SA.LD.r2
SA.LD.f.i
SA.LD.m.i
SA.LD.f.i2
SA.LD.m.i2

SA.CS.Vp
SA.CS.y2
SA.CS.f.pe2
SA.CS.m.pe2
SA.CS.f.h2
SA.CS.m.h2
SA.CS.cov.a
SA.CS.tbv
SA.CS.total.h2
SA.CS.sa2
SA.CS.r2
SA.CS.f.i
SA.CS.m.i
SA.CS.f.i2
SA.CS.m.i2

SA.cov.y
SA.cov.fpe
SA.cov.mpe
SA.cov.A.fge
SA.cov.A.mge
SA.COV.a.fld_mcs
SA.COV.a.fcs_mld
SA.COV.r


###############################
# TESTING VARIANCE COMPONENTS #
###############################

# YEAR
# Not fitting year covariance
SA.bivariate.model_no.year.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                         + at(trait, 1):F.AGE.CLASS
                                         + at(trait, 2):F.AGE.CLASS
                                         + at(trait, 1):M.AGE.CLASS
                                         + at(trait, 2):M.AGE.CLASS,
                                         random = ~ idh(trait, init = c(40, 0.5)):YEAR 
                                         + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                         + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                         + str(~corgh(trait):vm(female.animal, ainv) 
                                               + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                         + trait:ar1v(col):ar1(row),
                                         residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                         data = phenotypic.data,
                                         workspace = 8e8,  # Boost memory
                                         maxit = 100)
SA.bivariate.model_no.year.COV <- update(SA.bivariate.model_no.year.COV)
chi.sq <- -2*(SA.bivariate.model_no.year.COV$loglik - SA.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Year effect on laydate
fix.LD.year <- c(1e-8, 0.05)
names(fix.LD.year) <- c("F", "P")
SA.bivariate.model_no.LD.year <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ idh(trait, init = fix.LD.year):YEAR 
                                        + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                        + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) 
                                              + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                        data = phenotypic.data,
                                        workspace = 8e8,  # Boost memory
                                        maxit = 100)
SA.bivariate.model_no.LD.year <- update(SA.bivariate.model_no.LD.year)
chi.sq <- -2*(SA.bivariate.model_no.LD.year$loglik - SA.bivariate.model_no.year.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Year effect on clutch size
fix.CS.year <- c(40, 1e-8)
names(fix.CS.year) <- c("P", "F")
SA.bivariate.model_no.CS.year <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ idh(trait, init = fix.CS.year):YEAR 
                                        + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                        + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) 
                                              + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                        data = phenotypic.data,
                                        workspace = 8e8,  # Boost memory
                                        maxit = 100)
SA.bivariate.model_no.CS.year <- update(SA.bivariate.model_no.CS.year)
chi.sq <- -2*(SA.bivariate.model_no.CS.year$loglik - SA.bivariate.model_no.year.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# FEMALE PERMANENT ENVIRONMENT EFFECT
# Omitting female permanent environmental covariance
SA.bivariate.model_no.f.PE.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                         + at(trait, 1):F.AGE.CLASS
                                         + at(trait, 2):F.AGE.CLASS
                                         + at(trait, 1):M.AGE.CLASS
                                         + at(trait, 2):M.AGE.CLASS,
                                         random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                         + idh(trait, init = c(4, 0.5)):female.ID
                                         + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                         + str(~corgh(trait):vm(female.animal, ainv) 
                                               + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                         + trait:ar1v(col):ar1(row),
                                         residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                         data = phenotypic.data,
                                         workspace = 8e8,  # Boost memory
                                         maxit = 30)
SA.bivariate.model_no.f.PE.COV <- update(SA.bivariate.model_no.f.PE.COV)
chi.sq <- -2*(SA.bivariate.model_no.f.PE.COV$loglik - SA.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Female permanent environment effect on laydate
fix.LD.fPE <- c(1e-8, 0.5)
names(fix.LD.fPE) <- c("F", "P")
SA.bivariate.model_no.LD.fPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                       + at(trait, 1):F.AGE.CLASS
                                       + at(trait, 2):F.AGE.CLASS
                                       + at(trait, 1):M.AGE.CLASS
                                       + at(trait, 2):M.AGE.CLASS,
                                       random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                       + idh(trait, init = fix.LD.fPE):female.ID
                                       + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                       + str(~corgh(trait):vm(female.animal, ainv) 
                                             + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                       + trait:ar1v(col):ar1(row),
                                       residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                       data = phenotypic.data,
                                       workspace = 8e8,  # Boost memory
                                       maxit = 30)
SA.bivariate.model_no.LD.fPE <- update(SA.bivariate.model_no.LD.fPE)
chi.sq <- -2*(SA.bivariate.model_no.LD.fPE$loglik - SA.bivariate.model_no.f.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment effect on clutch size
fix.CS.fPE <- c(2, 1e-8)
names(fix.CS.fPE) <- c("P", "F")
SA.bivariate.model_no.CS.fPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                       + at(trait, 1):F.AGE.CLASS
                                       + at(trait, 2):F.AGE.CLASS
                                       + at(trait, 1):M.AGE.CLASS
                                       + at(trait, 2):M.AGE.CLASS,
                                       random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                       + idh(trait, init = fix.CS.fPE):female.ID
                                       + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                       + str(~corgh(trait):vm(female.animal, ainv) 
                                             + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                       + trait:ar1v(col):ar1(row),
                                       residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                       data = phenotypic.data,
                                       workspace = 8e8,  # Boost memory
                                       maxit = 30)
SA.bivariate.model_no.CS.fPE <- update(SA.bivariate.model_no.CS.fPE)
chi.sq <- -2*(SA.bivariate.model_no.CS.fPE$loglik - SA.bivariate.model_no.f.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# MALE PERMANENT ENVIRONMENT EFFECT
# Not fitting male permanent environmental covariance
SA.bivariate.model_no.m.PE.COV <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                         + at(trait, 1):F.AGE.CLASS
                                         + at(trait, 2):F.AGE.CLASS
                                         + at(trait, 1):M.AGE.CLASS
                                         + at(trait, 2):M.AGE.CLASS,
                                         random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                         + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                         + idh(trait, init = c(2, 0.05)):male.ID
                                         + str(~corgh(trait):vm(female.animal, ainv) 
                                               + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                         + trait:ar1v(col):ar1(row),
                                         residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                         data = phenotypic.data,
                                         workspace = 8e8,  # Boost memory
                                         maxit = 30)
SA.bivariate.model_no.m.PE.COV <- update(SA.bivariate.model_no.m.PE.COV)
chi.sq <- -2*(SA.bivariate.model_no.m.PE.COV$loglik - SA.bivariate.model$loglik)
chi.sq
1-pchisq(chi.sq, 1)

# Male permanent environment effect on laydate
fix.LD.mPE <- c(1e-8, 0.04)
names(fix.LD.mPE) <- c("F", "P")
SA.bivariate.model_no.LD.mPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                       + at(trait, 1):F.AGE.CLASS
                                       + at(trait, 2):F.AGE.CLASS
                                       + at(trait, 1):M.AGE.CLASS
                                       + at(trait, 2):M.AGE.CLASS,
                                       random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                       + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                       + idh(trait, init = fix.LD.mPE):male.ID
                                       + str(~corgh(trait):vm(female.animal, ainv) 
                                             + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                       + trait:ar1v(col):ar1(row),
                                       residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                       data = phenotypic.data,
                                       workspace = 8e8,  # Boost memory
                                       maxit = 30)
chi.sq <- -2*(SA.bivariate.model_no.LD.mPE$loglik - SA.bivariate.model_no.m.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male permanent environment effect on clutch size
fix.CS.mPE <- c(2, 1e-8)
names(fix.CS.mPE) <- c("P", "F")
SA.bivariate.model_no.CS.mPE <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                       + at(trait, 1):F.AGE.CLASS
                                       + at(trait, 2):F.AGE.CLASS
                                       + at(trait, 1):M.AGE.CLASS
                                       + at(trait, 2):M.AGE.CLASS,
                                       random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                       + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                       + idh(trait, init = fix.CS.mPE):male.ID
                                       + str(~corgh(trait):vm(female.animal, ainv) 
                                             + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv))
                                       + trait:ar1v(col):ar1(row),
                                       residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                       data = phenotypic.data,
                                       workspace = 8e8,  # Boost memory
                                       maxit = 30)
chi.sq <- -2*(SA.bivariate.model_no.CS.mPE$loglik - SA.bivariate.model_no.m.PE.COV$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value


# CONSTRANING COMPONENTS OF THE G-MATRIX
# Extract starting parameters
SA.gamma.extraction.model <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                    + at(trait, 1):F.AGE.CLASS
                                    + at(trait, 2):F.AGE.CLASS
                                    + at(trait, 1):M.AGE.CLASS
                                    + at(trait, 2):M.AGE.CLASS,
                                    random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR
                                    + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                    + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                    + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                          ~corgh(4):vm(female.animal, ainv))
                                    + trait:ar1v(col):ar1(row),
                                    residual = ~units:corgh(trait, init=c(-0.1, 25, 1)),
                                    data = phenotypic.data, 
                                    workspace = 6e8,
                                    start.values = "gammas.txt") 

# Female additive genetic effect on laydate
# Fixing covariances with female additive genetic variance in laydate to zero
gam.LD.f.COVs <- SA.gamma.extraction.model$vparameters.table  
gam.LD.f.COVs[c(10,11,13),]$Value <- 0
gam.LD.f.COVs[c(10,11,13),]$Constraint <- "F"
SA.bivariate.model_no.LD.f.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                          + at(trait, 1):F.AGE.CLASS
                                          + at(trait, 2):F.AGE.CLASS
                                          + at(trait, 1):M.AGE.CLASS
                                          + at(trait, 2):M.AGE.CLASS,
                                          random = ~ corgh(trait):YEAR
                                          + corgh(trait):female.ID
                                          + corgh(trait):male.ID
                                          + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                ~corgh(4):vm(female.animal, ainv))
                                          + trait:ar1v(col):ar1(row),
                                          residual = ~units:corgh(trait),
                                          data = phenotypic.data, 
                                          workspace = 6e8,
                                          maxiter = 30,
                                          G.param = gam.LD.f.COVs) 
SA.bivariate.model_no.LD.f.COVs <- update(SA.bivariate.model_no.LD.f.COVs)

# Additionally fixing female additive genetic variance in laydate close to zero
gam.LD.f.Va <-  gam.LD.f.COVs
gam.LD.f.Va[16,]$Value <- 1e-8
gam.LD.f.Va[16,]$Constraint <- "F"
SA.bivariate.model_no.LD.f.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ corgh(trait):YEAR
                                        + corgh(trait):female.ID
                                        + corgh(trait):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                              ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait),
                                        data = phenotypic.data, 
                                        workspace = 6e8,
                                        maxiter = 30,
                                        G.param = gam.LD.f.Va)
SA.bivariate.model_no.LD.f.Va <- update(SA.bivariate.model_no.LD.f.Va)
chi.sq <- -2*(SA.bivariate.model_no.LD.f.Va$loglik - SA.bivariate.model_no.LD.f.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female additive genetic effect on clutch size
# Fixing covariances with female additive genetic variance in clutch size to zero
gam.CS.f.COVs <- SA.gamma.extraction.model$vparameters.table
gam.CS.f.COVs[c(10,12,14),]$Value <- 0
gam.CS.f.COVs[c(10,12,14),]$Constraint <- "F"
SA.bivariate.model_no.CS.f.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                          + at(trait, 1):F.AGE.CLASS
                                          + at(trait, 2):F.AGE.CLASS
                                          + at(trait, 1):M.AGE.CLASS
                                          + at(trait, 2):M.AGE.CLASS,
                                          random = ~ corgh(trait):YEAR
                                          + corgh(trait):female.ID
                                          + corgh(trait):male.ID
                                          + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                ~corgh(4):vm(female.animal, ainv))
                                          + trait:ar1v(col):ar1(row),
                                          residual = ~units:corgh(trait),
                                          data = phenotypic.data, 
                                          workspace = 6e8,
                                          maxiter = 30,
                                          G.param = gam.CS.f.COVs)
SA.bivariate.model_no.CS.f.COVs <- update(SA.bivariate.model_no.CS.f.COVs)

# Additionally fixing female additive genetic variance in clutch size close to zero
gam.CS.f.Va <-  gam.CS.f.COVs 
gam.CS.f.Va[17,]$Value <- 1e-8
gam.CS.f.Va[17,]$Constraint <- "F"
SA.bivariate.model_no.CS.f.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ corgh(trait):YEAR
                                        + corgh(trait):female.ID
                                        + corgh(trait):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                              ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait),
                                        data = phenotypic.data, 
                                        workspace = 6e8,
                                        maxiter = 30,
                                        G.param = gam.CS.f.Va)
chi.sq <- -2*(SA.bivariate.model_no.CS.f.Va$loglik - SA.bivariate.model_no.CS.f.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic effect on laydate
# Fixing covariances with male additive genetic variance in laydate to zero
gam.LD.m.COVs <- SA.gamma.extraction.model$vparameters.table
gam.LD.m.COVs[c(11,12,15),]$Value <- 0
gam.LD.m.COVs[c(11,12,15),]$Constraint <- "F"
SA.bivariate.model_no.LD.m.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                          + at(trait, 1):F.AGE.CLASS
                                          + at(trait, 2):F.AGE.CLASS
                                          + at(trait, 1):M.AGE.CLASS
                                          + at(trait, 2):M.AGE.CLASS,
                                          random = ~ corgh(trait):YEAR
                                          + corgh(trait):female.ID
                                          + corgh(trait):male.ID
                                          + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                ~corgh(4):vm(female.animal, ainv))
                                          + trait:ar1v(col):ar1(row),
                                          residual = ~units:corgh(trait),
                                          data = phenotypic.data, 
                                          workspace = 6e8,
                                          maxiter = 30,
                                          G.param = gam.LD.m.COVs)
SA.bivariate.model_no.LD.m.COVs <- update(SA.bivariate.model_no.LD.m.COVs)

# Additionally fixing male additive genetic variance in laydate close to zero
gam.LD.m.Va <- gam.LD.m.COVs 
gam.LD.m.Va[18,]$Value <- 1e-8
gam.LD.m.Va[18,]$Constraint <- "F"
SA.bivariate.model_no.LD.m.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ corgh(trait):YEAR
                                        + corgh(trait):female.ID
                                        + corgh(trait):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                              ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait),
                                        data = phenotypic.data, 
                                        workspace = 6e8,
                                        maxiter = 30,
                                        G.param = gam.LD.m.Va)
SA.bivariate.model_no.LD.m.Va <- update(SA.bivariate.model_no.LD.m.Va)
chi.sq <- -2*(SA.bivariate.model_no.LD.m.Va$loglik - SA.bivariate.model_no.LD.m.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic effect on clutch size
# Fixing covariances with male additive genetic variance in clutch size to zero
gam.CS.m.COVs <- SA.gamma.extraction.model$vparameters.table
gam.CS.m.COVs[c(13,14,15),]$Value <- 0
gam.CS.m.COVs[c(13,14,15),]$Constraint <- "F"
SA.bivariate.model_no.CS.m.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                          + at(trait, 1):F.AGE.CLASS
                                          + at(trait, 2):F.AGE.CLASS
                                          + at(trait, 1):M.AGE.CLASS
                                          + at(trait, 2):M.AGE.CLASS,
                                          random = ~ corgh(trait):YEAR
                                          + corgh(trait):female.ID
                                          + corgh(trait):male.ID
                                          + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                ~corgh(4):vm(female.animal, ainv))
                                          + trait:ar1v(col):ar1(row),
                                          residual = ~units:corgh(trait),
                                          data = phenotypic.data, 
                                          workspace = 6e8,
                                          maxiter = 30,
                                          G.param = gam.CS.m.COVs)
SA.bivariate.model_no.CS.m.COVs <- update(SA.bivariate.model_no.CS.m.COVs)

# Additionally fixing male additive genetic variance in clutch size close to zero
gam.CS.m.Va <- gam.CS.m.COVs
gam.CS.m.Va[19,]$Value <- 1e-8
gam.CS.m.Va[19,]$Constraint <- "F"
SA.bivariate.model_no.CS.m.Va <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ corgh(trait):YEAR
                                        + corgh(trait):female.ID
                                        + corgh(trait):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                              ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait),
                                        data = phenotypic.data, 
                                        workspace = 6e8,
                                        maxiter = 30,
                                        G.param = gam.CS.m.Va)
SA.bivariate.model_no.CS.m.Va <- update(SA.bivariate.model_no.CS.m.Va)
chi.sq <- -2*(SA.bivariate.model_no.CS.m.Va$loglik - SA.bivariate.model_no.CS.m.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Covariance between female additive genetic effects on laydate and clutch size
gam.f.COVa <- SA.gamma.extraction.model$vparameters.table
gam.f.COVa[10,]$Value <- 0
gam.f.COVa[10,]$Constraint <- "F"
SA.bivariate.model_no.f.COVa <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                       + at(trait, 1):F.AGE.CLASS
                                       + at(trait, 2):F.AGE.CLASS
                                       + at(trait, 1):M.AGE.CLASS
                                       + at(trait, 2):M.AGE.CLASS,
                                       random = ~ corgh(trait):YEAR
                                       + corgh(trait):female.ID
                                       + corgh(trait):male.ID
                                       + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                             ~corgh(4):vm(female.animal, ainv))
                                       + trait:ar1v(col):ar1(row),
                                       residual = ~units:corgh(trait),
                                       data = phenotypic.data, 
                                       workspace = 6e8,
                                       maxiter = 30,
                                       G.param = gam.f.COVa)
SA.bivariate.model_no.f.COVa <- update(SA.bivariate.model_no.f.COVa)
chi.sq <- -2*(SA.bivariate.model_no.f.COVa$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Covariance between male additive genetic effects on laydate and clutch size
gam.m.COVa <- SA.gamma.extraction.model$vparameters.table
gam.m.COVa[15,]$Value <- 0
gam.m.COVa[15,]$Constraint <- "F"
SA.bivariate.model_no.COVa.mge <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                         + at(trait, 1):F.AGE.CLASS
                                         + at(trait, 2):F.AGE.CLASS
                                         + at(trait, 1):M.AGE.CLASS
                                         + at(trait, 2):M.AGE.CLASS,
                                         random = ~ corgh(trait):YEAR
                                         + corgh(trait):female.ID
                                         + corgh(trait):male.ID
                                         + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                               ~corgh(4):vm(female.animal, ainv))
                                         + trait:ar1v(col):ar1(row),
                                         residual = ~units:corgh(trait),
                                         data = phenotypic.data, 
                                         workspace = 6e8,
                                         maxiter = 30,
                                         G.param = gam.m.COVa)
SA.bivariate.model_no.COVa.mge <- update(SA.bivariate.model_no.COVa.mge)
chi.sq <- -2*(SA.bivariate.model_no.COVa.mge$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Covariance between female and male additive genetic effects on laydate
gam.COVa.LD <- SA.gamma.extraction.model$vparameters.table
gam.COVa.LD[11,]$Value <- 0
gam.COVa.LD[11,]$Constraint <- "F"
SA.bivariate.model_no.COVa.LD <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ corgh(trait):YEAR
                                        + corgh(trait):female.ID
                                        + corgh(trait):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                              ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait),
                                        data = phenotypic.data, 
                                        workspace = 6e8,
                                        maxiter = 30,
                                        G.param = gam.COVa.LD)
SA.bivariate.model_no.COVa.LD <- update(SA.bivariate.model_no.COVa.LD)
chi.sq <- -2*(SA.bivariate.model_no.COVa.LD$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Covariance between female and male additive genetic effects on clutch size
gam.COVa.CS <- SA.gamma.extraction.model$vparameters.table
gam.COVa.CS[14,]$Value <- 0
gam.COVa.CS[14,]$Constraint <- "F"
SA.bivariate.model_no.COVa.CS <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                        + at(trait, 1):F.AGE.CLASS
                                        + at(trait, 2):F.AGE.CLASS
                                        + at(trait, 1):M.AGE.CLASS
                                        + at(trait, 2):M.AGE.CLASS,
                                        random = ~ corgh(trait):YEAR
                                        + corgh(trait):female.ID
                                        + corgh(trait):male.ID
                                        + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                              ~corgh(4):vm(female.animal, ainv))
                                        + trait:ar1v(col):ar1(row),
                                        residual = ~units:corgh(trait),
                                        data = phenotypic.data, 
                                        workspace = 6e8,
                                        maxiter = 30,
                                        G.param = gam.COVa.CS)
SA.bivariate.model_no.COVa.CS <- update(SA.bivariate.model_no.COVa.CS)
chi.sq <- -2*(SA.bivariate.model_no.COVa.CS$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Additive genetic covariance between female laydate and male clutch size
gam.COVa.fLD_mCS <- SA.gamma.extraction.model$vparameters.table
gam.COVa.fLD_mCS[13,]$Value <- 0
gam.COVa.fLD_mCS[13,]$Constraint <- "F"
SA.bivariate.model_no.COVa.fLD_mCS <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait):YEAR
                                             + corgh(trait):female.ID
                                             + corgh(trait):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv))
                                             + trait:ar1v(col):ar1(row),
                                             residual = ~units:corgh(trait),
                                             data = phenotypic.data, 
                                             workspace = 6e8,
                                             maxiter = 30,
                                             G.param = gam.COVa.fLD_mCS)
SA.bivariate.model_no.COVa.fLD_mCS <- update(SA.bivariate.model_no.COVa.fLD_mCS)
chi.sq <- -2*(SA.bivariate.model_no.COVa.fLD_mCS$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# Additive genetic covariance between female clutch size and male laydate
gam.COVa.mLD_fCS <- SA.gamma.extraction.model$vparameters.table
gam.COVa.mLD_fCS[12,]$Value <- 0
gam.COVa.mLD_fCS[12,]$Constraint <- "F"
SA.bivariate.model_no.COVa.mLD_fCS <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                             + at(trait, 1):F.AGE.CLASS
                                             + at(trait, 2):F.AGE.CLASS
                                             + at(trait, 1):M.AGE.CLASS
                                             + at(trait, 2):M.AGE.CLASS,
                                             random = ~ corgh(trait):YEAR
                                             + corgh(trait):female.ID
                                             + corgh(trait):male.ID
                                             + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                   ~corgh(4):vm(female.animal, ainv))
                                             + trait:ar1v(col):ar1(row),
                                             residual = ~units:corgh(trait),
                                             data = phenotypic.data, 
                                             workspace = 6e8,
                                             maxiter = 30,
                                             G.param = gam.COVa.mLD_fCS)
SA.bivariate.model_no.COVa.mLD_fCS <- update(SA.bivariate.model_no.COVa.mLD_fCS)
chi.sq <- -2*(SA.bivariate.model_no.COVa.mLD_fCS$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)

# TESTING ENTIRE LAYDATE G-MATRIX
# Fixing all laydate-clutch size covariances to zero
gam.LDCS.COVs <- SA.gamma.extraction.model$vparameters.table
gam.LDCS.COVs[c(10,12,13,15),]$Value <- 0  # Fix relevant covariances to zero
gam.LDCS.COVs[c(10,12,13,15),]$Constraint <- "F"
SA.bivariate.model_no.LDCS.COVs <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                          + at(trait, 1):F.AGE.CLASS
                                          + at(trait, 2):F.AGE.CLASS
                                          + at(trait, 1):M.AGE.CLASS
                                          + at(trait, 2):M.AGE.CLASS,
                                          random = ~ corgh(trait):YEAR
                                          + corgh(trait):female.ID
                                          + corgh(trait):male.ID
                                          + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                                ~corgh(4):vm(female.animal, ainv))
                                          + trait:ar1v(col):ar1(row),
                                          residual = ~units:corgh(trait),
                                          data = phenotypic.data, 
                                          workspace = 6e8,
                                          maxiter = 30,
                                          G.param = gam.LDCS.COVs)
SA.bivariate.model_no.LDCS.COVs <- update(SA.bivariate.model_no.LDCS.COVs)

# Additionally fixing laydate (co)variances to zero
gam.LD.G <- gam.LDCS.COVs 
gam.LD.G[c(16,18),]$Value <- 1e-8  # Fix relevant variances to very small value
gam.LD.G[c(16,18),]$Constraint <- "F"
gam.LD.G[11,]$Value <- 0  # Fix intersexual genetic correlation for laydate to zero
gam.LD.G[11,]$Constraint <- "F"
SA.bivariate.model_no.LD.G <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                     + at(trait, 1):F.AGE.CLASS
                                     + at(trait, 2):F.AGE.CLASS
                                     + at(trait, 1):M.AGE.CLASS
                                     + at(trait, 2):M.AGE.CLASS,
                                     random = ~ corgh(trait):YEAR
                                     + corgh(trait):female.ID
                                     + corgh(trait):male.ID
                                     + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                           ~corgh(4):vm(female.animal, ainv))
                                     + trait:ar1v(col):ar1(row),
                                     residual = ~units:corgh(trait),
                                     data = phenotypic.data, 
                                     workspace = 6e8,
                                     maxiter = 30,
                                     G.param = gam.LD.G)
chi.sq <- -2*(SA.bivariate.model_no.LD.G$loglik - SA.bivariate.model_no.LDCS.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value

# Testing entire clutch size G-matrix
gam.CS.G <- gam.LDCS.COVs 
gam.CS.G[c(17,19),]$Value <- 1e-8  # Fix relevant variances to very small value
gam.CS.G[c(17,19),]$Constraint <- "F"
gam.CS.G[14,]$Value <- 0  # Fix intersexual genetic correlation for clutch size to zero
gam.CS.G[14,]$Constraint <- "F"
SA.bivariate.model_no.CS.G <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                     + at(trait, 1):F.AGE.CLASS
                                     + at(trait, 2):F.AGE.CLASS
                                     + at(trait, 1):M.AGE.CLASS
                                     + at(trait, 2):M.AGE.CLASS,
                                     random = ~ corgh(trait):YEAR
                                     + corgh(trait):female.ID
                                     + corgh(trait):male.ID
                                     + str(~corgh(trait):vm(female.animal, ainv) + corgh(trait):vm(male.animal, ainv), 
                                           ~corgh(4):vm(female.animal, ainv))
                                     + trait:ar1v(col):ar1(row),
                                     residual = ~units:corgh(trait),
                                     data = phenotypic.data, 
                                     workspace = 6e8,
                                     maxiter = 30,
                                     G.param = gam.CS.G)
SA.bivariate.model_no.CS.G <- update(SA.bivariate.model_no.CS.G)
chi.sq <- -2*(SA.bivariate.model_no.CS.G$loglik - SA.bivariate.model_no.LDCS.COVs$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value

# Testing variance attributable to spatial autocorrelation
SA.bivariate.model_no.SA <- asreml(fixed = cbind(LAYDATE, CLUTCHSIZE) ~ trait 
                                   + at(trait, 1):F.AGE.CLASS
                                   + at(trait, 2):F.AGE.CLASS
                                   + at(trait, 1):M.AGE.CLASS
                                   + at(trait, 2):M.AGE.CLASS,
                                   random = ~ corgh(trait, init = c(-0.1, 40, 0.5)):YEAR 
                                   + corgh(trait, init = c(-0.1, 4, 0.5)):female.ID
                                   + corgh(trait, init = c(-0.1, 2, 0.05)):male.ID
                                   + str(~corgh(trait):vm(female.animal, ainv) 
                                         + corgh(trait):vm(male.animal, ainv), ~corgh(4):vm(female.animal, ainv)),
                                   residual = ~units:corgh(trait, init = c(-0.1, 25, 1)),
                                   data = phenotypic.data,
                                   workspace = 8e8,  # Boost memory
                                   maxit = 30)
SA.bivariate.model_no.SA <- update(SA.bivariate.model_no.SA)
SA.bivariate.model_no.SA <- update(SA.bivariate.model_no.SA)
chi.sq <- -2*(SA.bivariate.model_no.SA$loglik - SA.bivariate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 2, 3)



#############
#############
##         ##
## FIGURES ##
##         ##
#############
#############

######################################
# FIGURE 1: PHENOTYPIC DISTRIBUTIONS #
######################################

par(mfrow = c(1, 2), mar = c(4, 4.5, 1, 1))
hist(phenotypic.data$LAYDATE, las = 1, breaks = 60, main = "", xlab = "", ylab = "", axes = "F", right = F)
mtext("Laydate (Aprilday)", side = 1, line = 2.5)
mtext("Frequency", side = 2, line = 3.3)
axis(side = 2, labels=c("",""), at = c(0, 420), lwd.ticks=0)
axis(2, at=seq(0, 400, by =100), las = 1)
axis(side = 1, at = seq(10, max(phenotypic.data$LAYDATE), by = 10))

par(mar = c(4, 1, 1, 4.5))
hist(phenotypic.data$CLUTCHSIZE, las = 1, ylim = c(0, 2000), breaks = 10, main = "", xlab = "", ylab = "", axes = "F", right = F)
mtext("Clutch size (eggs)", side = 1, line = 2.5)
mtext("Frequency", side = 4, line = 3.3)
axis(side = 4, at = seq(0, 2000, by = 500), las = 1)
axis(side = 1, at = seq(min(phenotypic.data$CLUTCHSIZE), max(phenotypic.data$CLUTCHSIZE), by = 2))


#############################
# FIGURE 2: BLOCKING DESIGN #
#############################

phenotypic.data$row.jitt <- jitter(as.numeric(phenotypic.data$row), 1.7)
phenotypic.data$col.jitt <- jitter(as.numeric(phenotypic.data$col), 1.7)
par(mfrow = c(1,1), mar = c(4.5,4.5,0.5,0.5))
plot(x = phenotypic.data$col.jitt, y = phenotypic.data$row.jitt, 
     bty = "n",
     pch = ".", 
     las = 1, 
     ylab = "Row", 
     xlab = "Column",
     cex.axis = 1.2,
     cex.lab = 1.2,
     ylim = c(0, 61), 
     xlim = c(1, 72))
axis(side = 1, at = c(0, 72), labels = c("", ""), lwd.ticks = 0)
axis(side = 2, at = c(0, 61), labels = c("", ""), lwd.ticks = 0)


#################################################################
# FIGURE 3: ADDITIVE GENETIC VARIANCE ESTIMATES FROM EACH MODEL #
#################################################################

# USING VARIANCE VALUES
classical.LD <- c(summary(classical.bivariate.model)$varcomp[11,"component"], 0, 0)

classical.CS <- c(summary(classical.bivariate.model)$varcomp[12,"component"], 0, 0)

extended.LD <- c(summary(extended.bivariate.model)$varcomp[19,"component"],
                 summary(extended.bivariate.model)$varcomp[21,"component"],
                 2*summary(extended.bivariate.model)$varcomp[14,"component"]*sqrt(summary(extended.bivariate.model)$varcomp[19,"component"]*summary(extended.bivariate.model)$varcomp[21,"component"]))

extended.CS <- c(summary(extended.bivariate.model)$varcomp[20,"component"],
                 summary(extended.bivariate.model)$varcomp[22,"component"],
                 2*summary(extended.bivariate.model)$varcomp[17,"component"]*sqrt(summary(extended.bivariate.model)$varcomp[20,"component"]*summary(extended.bivariate.model)$varcomp[22,"component"]))

spatial.LD <- c(summary(SA.bivariate.model)$varcomp[16,"component"],
                summary(SA.bivariate.model)$varcomp[18,"component"],
                2*summary(SA.bivariate.model)$varcomp[11,"component"]*sqrt(summary(SA.bivariate.model)$varcomp[16,"component"]*summary(SA.bivariate.model)$varcomp[18,"component"]))

spatial.CS <- c(summary(SA.bivariate.model)$varcomp[17,"component"],
                summary(SA.bivariate.model)$varcomp[19,"component"],
                2*summary(SA.bivariate.model)$varcomp[14,"component"]*sqrt(summary(SA.bivariate.model)$varcomp[17,"component"]*summary(SA.bivariate.model)$varcomp[19,"component"]))


LD.data <- data.frame(classical.LD, extended.LD, spatial.LD)
rownames(LD.data) <- c("female.Va", "male.Va", "mf.COVa")
LD.data <- t(LD.data)
LD.data <- t(LD.data)

CS.data <- data.frame(classical.CS, extended.CS, spatial.CS)
rownames(CS.data) <- c("female.Va", "male.Va", "mf.COVa")
CS.data <- t(CS.data)
CS.data <- t(CS.data)

par(mfrow = c(2,1), mar = c(0.5, 5, 4.5, 0.2))
barplot(LD.data, 
        space = c(0, 0.4, 0.4),
        las = 1,
        names.arg = c("", "", ""),
        legend = c("female additive genetic variance", "male additive genetic variance", "male-female additive genetic covariance"),
        args.legend = list(bty = "n", cex = 0.85, x = 3.8, y = 14.5),
        ylim = range(pretty(c(0, (max(sum(LD.data[,1]), sum(LD.data[,2]), sum(LD.data[,3])))))))
mtext("Additive genetic variance", side = 2, line = 4, cex = 1)
mtext(expression(paste("of laydate (days" ^ "2",")")), side = 2, line = 2.7, cex = 1)

par(mar = c(2, 5, 2.5, 0.2))
barplot(CS.data,
        space = c(0, 0.4,0.4),
        las = 1, 
        names.arg = c("classic", "extended", "SA"),
        ylim = c(0, 0.6)) # range(pretty(c(0,(max(sum(CS.data[,1]), sum(CS.data[,2]), sum(CS.data[,3])))))))
axis(side = 2, at = c(0, 0.61), labels = F, tick = F)
mtext("Additive genetic variance", side = 2, line = 4, cex = 1)
mtext(expression(paste("of clutch size (eggs" ^ "2",")")), side = 2, line = 2.7, cex = 1)



##############################
##############################
##                          ##
##  SUPPLEMENTARY ANALYSES  ##
##                          ##
##############################
##############################

##############################
# CLASSICAL MODEL OF LAYDATE #
##############################

classical.laydate.model <- asreml(fixed = LAYDATE ~ F.AGE.CLASS,
                                  random = ~round
                                  + YEAR
                                  + female.ID
                                  + female.maternal.ID
                                  + vm(female.animal, ainv),
                                  residual = ~idv(units),
                                  data = phenotypic.data,
                                  workspace = 32e6,  # Boost memory
                                  maxit = 20,
                                  na.action = na.method(x = "include", y = "omit"))

summary(classical.laydate.model)$varcomp[,c("component", "std.error", "bound")]
summary(classical.laydate.model,  coef=T)$coef.fixed
wald.asreml(classical.laydate.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

uni.LD.Vp       <- vpredict(classical.laydate.model, LD.Vphenotypic ~ V1+V3+V4+V5+V6)  # Year excluded from phenotypic variance
uni.LD.plot2    <- vpredict(classical.laydate.model, LD.plot2 ~ V1 / (V1+V3+V4+V5+V6))
uni.LD.y2       <- vpredict(classical.laydate.model, LD.y2 ~ V2 / (V1+V3+V4+V5+V6))
uni.LD.f.pe2    <- vpredict(classical.laydate.model, LD.f.pe2 ~ V4 / (V1+V3+V4+V5+V6))
uni.LD.f.mat2    <- vpredict(classical.laydate.model, LD.m.pe2 ~ V3/ (V1+V3+V4+V5+V6))
uni.LD.f.h2     <- vpredict(classical.laydate.model, LD.f.h2 ~ V5 / (V1+V3+V4+V5+V6))
uni.LD.r2       <- vpredict(classical.laydate.model, LD.r2 ~ V6 / (V1+V3+V4+V5+V6))

uni.LD.Vp
uni.LD.plot2
uni.LD.y2
uni.LD.f.pe2
uni.LD.f.mat2
uni.LD.f.h2
uni.LD.r2

# TESTING VARIANCE COMPONENTS
# Plot
classical.laydate.model_no.plot <- asreml(fixed = LAYDATE ~ F.AGE.CLASS,
                                          random = ~YEAR
                                          + female.ID
                                          + female.maternal.ID
                                          + vm(female.animal, ainv),
                                          residual = ~idv(units),
                                          data = phenotypic.data,
                                          workspace = 32e6,  # Boost memory
                                          maxit = 20,
                                          na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.laydate.model_no.plot$loglik - classical.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Year
classical.laydate.model_no.year <- asreml(fixed = LAYDATE ~ F.AGE.CLASS,
                                          random = ~round
                                          + female.ID
                                          + female.maternal.ID
                                          + vm(female.animal, ainv),
                                          residual = ~idv(units),
                                          data = phenotypic.data,
                                          workspace = 32e6,  # Boost memory
                                          maxit = 20,
                                          na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.laydate.model_no.year$loglik - classical.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment
classical.laydate.model_no.fPE <- asreml(fixed = LAYDATE ~ F.AGE.CLASS,
                                         random = ~round
                                         + YEAR
                                         + female.maternal.ID
                                         + vm(female.animal, ainv),
                                         residual = ~idv(units),
                                         data = phenotypic.data,
                                         workspace = 32e6,  # Boost memory
                                         maxit = 20,
                                         na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.laydate.model_no.fPE$loglik - classical.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female maternal environment
classical.laydate.model_no.fmat <- asreml(fixed = LAYDATE ~ F.AGE.CLASS,
                                          random = ~round
                                          + YEAR
                                          + female.ID
                                          + vm(female.animal, ainv),
                                          residual = ~idv(units),
                                          data = phenotypic.data,
                                          workspace = 32e6,  # Boost memory
                                          maxit = 20,
                                          na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.laydate.model_no.fmat$loglik - classical.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female additive genetic
classical.laydate.model_no.fVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS,
                                         random = ~round
                                         + YEAR
                                         + female.ID
                                         + female.maternal.ID,
                                         residual = ~idv(units),
                                         data = phenotypic.data,
                                         workspace = 32e6,  # Boost memory
                                         maxit = 20,
                                         na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.laydate.model_no.fVa$loglik - classical.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value


#############################
# EXTENDED MODEL OF LAYDATE #
#############################

extended.laydate.model <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                 + M.AGE.CLASS,
                                 random = ~ round
                                 + YEAR 
                                 + female.ID
                                 + male.ID
                                 + female.maternal.ID
                                 + male.maternal.ID
                                 + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                       ~corgh(2):vm(female.animal, ainv)),
                                 residual = ~idv(units),
                                 data = phenotypic.data, 
                                 workspace = 6e8,  # Boost memory
                                 maxit = 20,
                                 na.action = na.method(x = "include", y = "omit"))
summary(extended.laydate.model)$varcomp[,c("component", "std.error", "bound")]
summary(extended.laydate.model, coef=T)$coef.fixed
wald.asreml(extended.laydate.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

ext.LD.Vp       <- vpredict(extended.laydate.model, ext.LD.Vphenotypic ~ V1+V3+V4+V5+V6+V7+V8+V9+V10)  # Year excluded from phenotypic variance
ext.LD.plot2    <- vpredict(extended.laydate.model, ext.LD.plot2 ~ V1 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.y2       <- vpredict(extended.laydate.model, ext.LD.y2 ~ V2 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.f.pe2    <- vpredict(extended.laydate.model, ext.LD.f.pe2 ~ V3 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.m.pe2    <- vpredict(extended.laydate.model, ext.LD.m.pe2 ~ V4 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.f.mat2   <- vpredict(extended.laydate.model, ext.LD.f.mat2 ~ V5 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.m.mat2   <- vpredict(extended.laydate.model, ext.LD.m.mat2 ~ V6 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.f.h2     <- vpredict(extended.laydate.model, ext.LD.f.h2 ~ V8 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.m.h2     <- vpredict(extended.laydate.model, ext.LD.m.h2 ~ V9 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.cov.a    <- vpredict(extended.laydate.model, ext.LD.cov.a ~ V7 * sqrt(V8*V9))
ext.LD.tbv      <- vpredict(extended.laydate.model, ext.LD.tbv ~ V8 + 2*(V7*sqrt(V8*V9)) + V9)
ext.LD.total.h2 <- vpredict(extended.laydate.model, ext.LD.total.h2 ~  (V8+2*(V7*sqrt(V8*V9))+V9)/ (V1+V3+V4+V5+V6+V7+V8+V9+V10))
ext.LD.r2       <- vpredict(extended.laydate.model, ext.LD.r2 ~ V10 / (V1+V3+V4+V5+V6+V7+V8+V9+V10))

ext.LD.Vp
ext.LD.plot2
ext.LD.y2
ext.LD.f.pe2
ext.LD.m.pe2
ext.LD.f.mat2
ext.LD.m.mat2
ext.LD.f.h2
ext.LD.m.h2
ext.LD.cov.a
ext.LD.tbv
ext.LD.total.h2
ext.LD.r2

# TESTING VARIANCE COMPONENTS
# Plot
extended.laydate.model_no.plot <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                         + M.AGE.CLASS,
                                         random = ~ YEAR 
                                         + female.ID
                                         + male.ID
                                         + female.maternal.ID
                                         + male.maternal.ID
                                         + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                               ~corgh(2):vm(female.animal, ainv)),
                                         residual = ~idv(units),
                                         data = phenotypic.data, 
                                         workspace = 6e8,  # Boost memory
                                         maxit = 20,
                                         na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.plot$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Year
extended.laydate.model_no.year <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                         + M.AGE.CLASS,
                                         random = ~ round
                                         + female.ID
                                         + male.ID
                                         + female.maternal.ID
                                         + male.maternal.ID
                                         + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                               ~corgh(2):vm(female.animal, ainv)),
                                         residual = ~idv(units),
                                         data = phenotypic.data, 
                                         workspace = 6e8,  # Boost memory
                                         maxit = 20,
                                         na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.year$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment
extended.laydate.model_no.fPE <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                        + M.AGE.CLASS,
                                        random = ~ round
                                        + YEAR 
                                        + male.ID
                                        + female.maternal.ID
                                        + male.maternal.ID
                                        + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                              ~corgh(2):vm(female.animal, ainv)),
                                        residual = ~idv(units),
                                        data = phenotypic.data, 
                                        workspace = 6e8,  # Boost memory
                                        maxit = 20,
                                        na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.fPE$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male permanent environment
extended.laydate.model_no.mPE <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                        + M.AGE.CLASS,
                                        random = ~ round
                                        + YEAR 
                                        + female.ID
                                        + female.maternal.ID
                                        + male.maternal.ID
                                        + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                              ~corgh(2):vm(female.animal, ainv)),
                                        residual = ~idv(units),
                                        data = phenotypic.data, 
                                        workspace = 6e8,  # Boost memory
                                        maxit = 20,
                                        na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.mPE$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female maternal environment
extended.laydate.model_no.fmat <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                         + M.AGE.CLASS,
                                         random = ~ round
                                         + YEAR 
                                         + female.ID
                                         + male.ID
                                         + male.maternal.ID
                                         + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                               ~corgh(2):vm(female.animal, ainv)),
                                         residual = ~idv(units),
                                         data = phenotypic.data, 
                                         workspace = 6e8,  # Boost memory
                                         maxit = 20,
                                         na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.fmat$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male maternal environment - bound at zero

# Intersexaul additive genetic covariance
extended.laydate.model_no.COVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                         + M.AGE.CLASS,
                                         random = ~ round
                                         + YEAR 
                                         + female.ID
                                         + male.ID
                                         + female.maternal.ID
                                         + male.maternal.ID
                                         + vm(female.animal, ainv) 
                                         + vm(male.animal, ainv),
                                         residual = ~idv(units),
                                         data = phenotypic.data, 
                                         workspace = 6e8,  # Boost memory
                                         maxit = 20,
                                         na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.COVa$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)  # P-value

# Female additive genetic
extended.laydate.model_no.fVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                        + M.AGE.CLASS,
                                        random = ~ round
                                        + YEAR 
                                        + female.ID
                                        + male.ID
                                        + female.maternal.ID
                                        + male.maternal.ID
                                        + vm(male.animal, ainv),
                                        residual = ~idv(units),
                                        data = phenotypic.data, 
                                        workspace = 6e8,  # Boost memory
                                        maxit = 20,
                                        na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.fVa$loglik - extended.laydate.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic
extended.laydate.model_no.mVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                        + M.AGE.CLASS,
                                        random = ~ round
                                        + YEAR 
                                        + female.ID
                                        + male.ID
                                        + female.maternal.ID
                                        + male.maternal.ID
                                        + vm(female.animal, ainv),
                                        residual = ~idv(units),
                                        data = phenotypic.data, 
                                        workspace = 6e8,  # Boost memory
                                        maxit = 20,
                                        na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.mVa$loglik - extended.laydate.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# G-matrix
extended.laydate.model_no.G <-asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                     + M.AGE.CLASS,
                                     random = ~ round
                                     + YEAR 
                                     + female.ID
                                     + male.ID
                                     + female.maternal.ID
                                     + male.maternal.ID,
                                     residual = ~idv(units),
                                     data = phenotypic.data, 
                                     workspace = 6e8,  # Boost memory
                                     maxit = 20,
                                     na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.laydate.model_no.G$loglik - extended.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value

##################################
# CLASSICAL MODEL OF CLUTCH SIZE #
##################################

classical.clutchsize.model <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS,
                                     random = ~round
                                     + YEAR
                                     + female.ID
                                     + female.maternal.ID
                                     + vm(female.animal, ainv),
                                     residual = ~idv(units),
                                     data = phenotypic.data,
                                     workspace = 32e6,  # Boost memory
                                     maxit = 20,
                                     na.action = na.method(x = "include", y = "omit"))
classical.clutchsize.model <- update(classical.clutchsize.model)
summary(classical.clutchsize.model)$varcomp[,c("component", "std.error", "bound")]
summary(classical.clutchsize.model,  coef=T)$coef.fixed
wald.asreml(classical.clutchsize.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

uni.CS.Vp       <- vpredict(classical.clutchsize.model, CS.Vphenotypic ~ V1+V3+V4+V5+V6)  # Year excluded from phenotypic variance
uni.CS.plot2    <- vpredict(classical.clutchsize.model, CS.plot2 ~ V1 / (V1+V3+V4+V5+V6))
uni.CS.y2       <- vpredict(classical.clutchsize.model, CS.y2 ~ V2 / (V1+V3+V4+V5+V6))
uni.CS.f.pe2    <- vpredict(classical.clutchsize.model, CS.f.pe2 ~ V4 / (V1+V3+V4+V5+V6))
uni.CS.f.mat2    <- vpredict(classical.clutchsize.model, CS.m.pe2 ~ V3/ (V1+V3+V4+V5+V6))
uni.CS.f.h2     <- vpredict(classical.clutchsize.model, CS.f.h2 ~ V5 / (V1+V3+V4+V5+V6))
uni.CS.r2       <- vpredict(classical.clutchsize.model, CS.r2 ~ V6 / (V1+V3+V4+V5+V6))

uni.CS.Vp
uni.CS.plot2
uni.CS.y2
uni.CS.f.pe2
uni.CS.f.mat2
uni.CS.f.h2
uni.CS.r2

# TESTING VARIANCE COMPONENTS
# Plot
classical.clutchsize.model_no.plot <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS,
                                             random = ~YEAR
                                             + female.ID
                                             + female.maternal.ID
                                             + vm(female.animal, ainv),
                                             residual = ~idv(units),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20,
                                             na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.clutchsize.model_no.plot$loglik - classical.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Year
classical.clutchsize.model_no.year <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS,
                                             random = ~round
                                             + female.ID
                                             + female.maternal.ID
                                             + vm(female.animal, ainv),
                                             residual = ~idv(units),
                                             data = phenotypic.data,
                                             workspace = 32e6,  # Boost memory
                                             maxit = 20,
                                             na.action = na.method(x = "include", y = "omit"))
classical.clutchsize.model_no.year <- update(classical.clutchsize.model_no.year)
chi.sq <- -2*(classical.clutchsize.model_no.year$loglik - classical.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment
classical.clutchsize.model_no.fPE <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS,
                                            random = ~round
                                            + YEAR
                                            + female.maternal.ID
                                            + vm(female.animal, ainv),
                                            residual = ~idv(units),
                                            data = phenotypic.data,
                                            workspace = 32e6,  # Boost memory
                                            maxit = 20,
                                            na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.clutchsize.model_no.fPE$loglik - classical.clutchsize.model$loglik)
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female maternal environment - bound at zero

# Female additive genetic
classical.clutchsize.model_no.fVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS,
                                            random = ~round
                                            + YEAR
                                            + female.ID
                                            + female.maternal.ID,
                                            residual = ~idv(units),
                                            data = phenotypic.data,
                                            workspace = 32e6,  # Boost memory
                                            maxit = 20,
                                            na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(classical.clutchsize.model_no.fVa$loglik - classical.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

################################
# EXTENDED MODEL OF CLUTCHSIZE #
################################

extended.clutchsize.model <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                    + M.AGE.CLASS,
                                    random = ~ round
                                    + YEAR 
                                    + female.ID
                                    + male.ID
                                    + female.maternal.ID
                                    + male.maternal.ID
                                    + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                          ~corgh(2):vm(female.animal, ainv)),
                                    residual = ~idv(units),
                                    data = phenotypic.data, 
                                    workspace = 6e8,  # Boost memory
                                    maxit = 20,
                                    na.action = na.method(x = "include", y = "omit"))
summary(extended.clutchsize.model)$varcomp[,c("component", "std.error", "bound")]
summary(extended.clutchsize.model, coef=T)$coef.fixed
wald.asreml(extended.clutchsize.model, ssType = "conditional", denDF = "numeric")  # P-value for fixed effects

ext.CS.Vp       <- vpredict(extended.clutchsize.model, ext.CS.Vphenotypic ~ V1+V3+V4+V5+V6+V8+V9+V10)  # Year excluded from phenotypic variance
ext.CS.plot2    <- vpredict(extended.clutchsize.model, ext.CS.plot2 ~ V1 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.y2       <- vpredict(extended.clutchsize.model, ext.CS.y2 ~ V2 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.f.pe2    <- vpredict(extended.clutchsize.model, ext.CS.f.pe2 ~ V3 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.m.pe2    <- vpredict(extended.clutchsize.model, ext.CS.m.pe2 ~ V4 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.f.mat2   <- vpredict(extended.clutchsize.model, ext.CS.f.mat2 ~ V5 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.m.mat2   <- vpredict(extended.clutchsize.model, ext.CS.m.mat2 ~ V6 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.f.h2     <- vpredict(extended.clutchsize.model, ext.CS.f.h2 ~ V8 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.m.h2     <- vpredict(extended.clutchsize.model, ext.CS.m.h2 ~ V9 / (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.cov.a    <- vpredict(extended.clutchsize.model, ext.CS.cov.a ~ V7 * sqrt(V8*V9))
ext.CS.tbv      <- vpredict(extended.clutchsize.model, ext.CS.tbv ~ V8 + 2*(V7*sqrt(V8*V9)) + V9)
ext.CS.total.h2 <- vpredict(extended.clutchsize.model, ext.CS.total.h2 ~  (V8+2*(V7*sqrt(V8*V9))+V9)/ (V1+V3+V4+V5+V6+V8+V9+V10))
ext.CS.r2       <- vpredict(extended.clutchsize.model, ext.CS.r2 ~ V10 / (V1+V3+V4+V5+V6+V8+V9+V10))

ext.CS.Vp
ext.CS.plot2
ext.CS.y2
ext.CS.f.pe2
ext.CS.m.pe2
ext.CS.f.mat2
ext.CS.m.mat2
ext.CS.f.h2
ext.CS.m.h2
ext.CS.cov.a
ext.CS.tbv
ext.CS.total.h2
ext.CS.r2

# TESTING VARIANCE COMPONENTS
# Plot
extended.clutchsize.model_no.plot <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                            + M.AGE.CLASS,
                                            random = ~ YEAR 
                                            + female.ID
                                            + male.ID
                                            + female.maternal.ID
                                            + male.maternal.ID
                                            + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                                  ~corgh(2):vm(female.animal, ainv)),
                                            residual = ~idv(units),
                                            data = phenotypic.data, 
                                            workspace = 6e8,  # Boost memory
                                            maxit = 20,
                                            na.action = na.method(x = "include", y = "omit"))
extended.clutchsize.model_no.plot <- update(extended.clutchsize.model_no.plot)
chi.sq <- -2*(extended.clutchsize.model_no.plot$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Year
extended.clutchsize.model_no.year <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                            + M.AGE.CLASS,
                                            random = ~ round
                                            + female.ID
                                            + male.ID
                                            + female.maternal.ID
                                            + male.maternal.ID
                                            + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                                  ~corgh(2):vm(female.animal, ainv)),
                                            residual = ~idv(units),
                                            data = phenotypic.data, 
                                            workspace = 6e8,  # Boost memory
                                            maxit = 20,
                                            na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.clutchsize.model_no.year$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment
extended.clutchsize.model_no.fPE <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                           + M.AGE.CLASS,
                                           random = ~ round
                                           + YEAR 
                                           + male.ID
                                           + female.maternal.ID
                                           + male.maternal.ID
                                           + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                                 ~corgh(2):vm(female.animal, ainv)),
                                           residual = ~idv(units),
                                           data = phenotypic.data, 
                                           workspace = 6e8,  # Boost memory
                                           maxit = 20,
                                           na.action = na.method(x = "include", y = "omit"))
extended.clutchsize.model_no.fPE <- update(extended.clutchsize.model_no.fPE)
chi.sq <- -2*(extended.clutchsize.model_no.fPE$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male permanent environment
extended.clutchsize.model_no.mPE <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                           + M.AGE.CLASS,
                                           random = ~ round
                                           + YEAR 
                                           + female.ID
                                           + female.maternal.ID
                                           + male.maternal.ID
                                           + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                                 ~corgh(2):vm(female.animal, ainv)),
                                           residual = ~idv(units),
                                           data = phenotypic.data, 
                                           workspace = 6e8,  # Boost memory
                                           maxit = 20,
                                           na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.clutchsize.model_no.mPE$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female maternal environment - bound at zero

# Male maternal environment
extended.clutchsize.model_no.mmat <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                            + M.AGE.CLASS,
                                            random = ~ round
                                            + YEAR 
                                            + female.ID
                                            + male.ID
                                            + female.maternal.ID
                                            + str(~vm(female.animal, ainv) + vm(male.animal, ainv), 
                                                  ~corgh(2):vm(female.animal, ainv)),
                                            residual = ~idv(units),
                                            data = phenotypic.data, 
                                            workspace = 6e8,  # Boost memory
                                            maxit = 20,
                                            na.action = na.method(x = "include", y = "omit"))
extended.clutchsize.model_no.mmat <- update(extended.clutchsize.model_no.mmat)
chi.sq <- -2*(extended.clutchsize.model_no.mmat$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Intersexaul additive genetic covariance
extended.clutchsize.model_no.COVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                            + M.AGE.CLASS,
                                            random = ~ round
                                            + YEAR 
                                            + female.ID
                                            + male.ID
                                            + female.maternal.ID
                                            + male.maternal.ID
                                            + vm(female.animal, ainv) 
                                            + vm(male.animal, ainv),
                                            residual = ~idv(units),
                                            data = phenotypic.data, 
                                            workspace = 6e8,  # Boost memory
                                            maxit = 20,
                                            na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.clutchsize.model_no.COVa$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)  # P-value

# Female additive genetic
extended.clutchsize.model_no.fVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                           + M.AGE.CLASS,
                                           random = ~ round
                                           + YEAR 
                                           + female.ID
                                           + male.ID
                                           + female.maternal.ID
                                           + male.maternal.ID
                                           + vm(male.animal, ainv),
                                           residual = ~idv(units),
                                           data = phenotypic.data, 
                                           workspace = 6e8,  # Boost memory
                                           maxit = 20,
                                           na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.clutchsize.model_no.fVa$loglik - extended.clutchsize.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic
extended.clutchsize.model_no.mVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                           + M.AGE.CLASS,
                                           random = ~ round
                                           + YEAR 
                                           + female.ID
                                           + male.ID
                                           + female.maternal.ID
                                           + male.maternal.ID
                                           + vm(female.animal, ainv),
                                           residual = ~idv(units),
                                           data = phenotypic.data, 
                                           workspace = 6e8,  # Boost memory
                                           maxit = 20,
                                           na.action = na.method(x = "include", y = "omit"))
extended.clutchsize.model_no.mVa <- update(extended.clutchsize.model_no.mVa)
chi.sq <- -2*(extended.clutchsize.model_no.mVa$loglik - extended.clutchsize.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# G-matrix
extended.clutchsize.model_no.G <-asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                        + M.AGE.CLASS,
                                        random = ~ round
                                        + YEAR 
                                        + female.ID
                                        + male.ID
                                        + female.maternal.ID
                                        + male.maternal.ID,
                                        residual = ~idv(units),
                                        data = phenotypic.data, 
                                        workspace = 6e8,  # Boost memory
                                        maxit = 20,
                                        na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(extended.clutchsize.model_no.G$loglik - extended.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value


############################################
# SPATIAL AUTOCORRELATION MODEL OF LAYDATE #
############################################

SA.laydate.model <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                           + M.AGE.CLASS,
                           random = ~ YEAR 
                           + female.ID
                           + male.ID
                           + str(~vm(female.animal, ainv) 
                                 + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                           + ar1v(col):ar1(row),
                           residual = ~idv(units),
                           data = phenotypic.data,
                           workspace = 8e8,  # Boost memory
                           maxit = 30,
                           na.action = na.method(x = "include", y = "omit"))
summary(SA.laydate.model)$varcomp[,c("component", "std.error", "bound")]
summary(SA.laydate.model, coef=T)$coef.fixed
wald.asreml(SA.laydate.model, ssType = "conditional", denDF = "numeric")

SA.LD.Vp       <- vpredict(SA.laydate.model, SA.LD.Vphenotypic ~ V2+V3+V5+V6+V8+V10)  # Year excluded from phenotypic variance
SA.LD.y2       <- vpredict(SA.laydate.model, SA.LD.y2 ~ V1 / (V2+V3+V5+V6+V8+V10))
SA.LD.f.pe2    <- vpredict(SA.laydate.model, SA.LD.f.pe2 ~ V2 / (V2+V3+V5+V6+V8+V10))
SA.LD.m.pe2    <- vpredict(SA.laydate.model, SA.LD.m.pe2 ~ V3 / (V2+V3+V5+V6+V8+V10))
SA.LD.f.h2     <- vpredict(SA.laydate.model, SA.LD.f.h2 ~ V5 / (V2+V3+V5+V6+V8+V10))
SA.LD.m.h2     <- vpredict(SA.laydate.model, SA.LD.m.h2 ~ V6 / (V2+V3+V5+V6+V8+V10))
SA.LD.cov.a    <- vpredict(SA.laydate.model, SA.LD.cov.a ~ V4 * sqrt(V5*V6))
SA.LD.tbv      <- vpredict(SA.laydate.model, SA.LD.tbv ~ V5 + 2*(V4*sqrt(V5*V6)) + V6)
SA.LD.total.h2 <- vpredict(SA.laydate.model, SA.LD.total.h2 ~ (V5+2*(V4*sqrt(V5*V6))+V6)/ (V2+V3+V5+V6+V8+V10))
SA.LD.sa2      <- vpredict(SA.laydate.model, SA.LD.sa2 ~ V8/(V2+V3+V5+V6+V8+V10))
SA.LD.r2       <- vpredict(SA.laydate.model, SA.LD.r2 ~ V10 / (V2+V3+V5+V6+V8+V10))

SA.LD.Vp
SA.LD.y2
SA.LD.f.pe2
SA.LD.m.pe2
SA.LD.f.h2
SA.LD.m.h2
SA.LD.cov.a
SA.LD.tbv
SA.LD.total.h2
SA.LD.sa2
SA.LD.r2

# TESTING VARIANCE COMPONENTS
# Year
SA.laydate.model_no.year <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                   + M.AGE.CLASS,
                                   random = ~ female.ID
                                   + male.ID
                                   + str(~vm(female.animal, ainv) 
                                         + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                                   + ar1v(col):ar1(row),
                                   residual = ~idv(units),
                                   data = phenotypic.data,
                                   workspace = 8e8,  # Boost memory
                                   maxit = 30,
                                   na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.year$loglik - SA.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment
SA.laydate.model_no.fPE <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                  + M.AGE.CLASS,
                                  random = ~ YEAR 
                                  + male.ID
                                  + str(~vm(female.animal, ainv) 
                                        + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                                  + ar1v(col):ar1(row),
                                  residual = ~idv(units),
                                  data = phenotypic.data,
                                  workspace = 8e8,  # Boost memory
                                  maxit = 30,
                                  na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.fPE$loglik - SA.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male permanent environment
SA.laydate.model_no.mPE <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                  + M.AGE.CLASS,
                                  random = ~ YEAR 
                                  + female.ID
                                  + str(~vm(female.animal, ainv) 
                                        + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                                  + ar1v(col):ar1(row),
                                  residual = ~idv(units),
                                  data = phenotypic.data,
                                  workspace = 8e8,  # Boost memory
                                  maxit = 30,
                                  na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.mPE$loglik - SA.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Intersexaul additive genetic covariance
SA.laydate.model_no.COVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                   + M.AGE.CLASS,
                                   random = ~ YEAR 
                                   + female.ID
                                   + male.ID
                                   + vm(female.animal, ainv)
                                   + vm(male.animal, ainv)
                                   + ar1v(col):ar1(row),
                                   residual = ~idv(units),
                                   data = phenotypic.data,
                                   workspace = 8e8,  # Boost memory
                                   maxit = 30,
                                   na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.COVa$loglik - SA.laydate.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)  # P-value

# Female additive genetic
SA.laydate.model_no.fVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                  + M.AGE.CLASS,
                                  random = ~ YEAR 
                                  + female.ID
                                  + male.ID
                                  + vm(male.animal, ainv)
                                  + ar1v(col):ar1(row),
                                  residual = ~idv(units),
                                  data = phenotypic.data,
                                  workspace = 8e8,  # Boost memory
                                  maxit = 30,
                                  na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.fVa$loglik - SA.laydate.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic
SA.laydate.model_no.mVa <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                  + M.AGE.CLASS,
                                  random = ~ YEAR 
                                  + female.ID
                                  + male.ID
                                  + vm(female.animal, ainv)
                                  + ar1v(col):ar1(row),
                                  residual = ~idv(units),
                                  data = phenotypic.data,
                                  workspace = 8e8,  # Boost memory
                                  maxit = 30,
                                  na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.mVa$loglik - SA.laydate.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# G-matrix
SA.laydate.model_no.G <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                + M.AGE.CLASS,
                                random = ~ YEAR 
                                + female.ID
                                + male.ID
                                + ar1v(col):ar1(row),
                                residual = ~idv(units),
                                data = phenotypic.data,
                                workspace = 8e8,  # Boost memory
                                maxit = 30,
                                na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.G$loglik - SA.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value

# Spatial autocorrelation - row
SA.laydate.model_no.SA <- asreml(fixed = LAYDATE ~ F.AGE.CLASS
                                 + M.AGE.CLASS,
                                 random = ~ YEAR 
                                 + female.ID
                                 + male.ID
                                 + str(~vm(female.animal, ainv) 
                                       + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv)),
                                 residual = ~idv(units),
                                 data = phenotypic.data,
                                 workspace = 8e8,  # Boost memory
                                 maxit = 30,
                                 na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.laydate.model_no.SA$loglik - SA.laydate.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 2, 3)  # P-value


################################################
# SPATIAL AUTOCORRELATION MODEL OF CLUTCH SIZE #
################################################

SA.clutchsize.model <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                              + M.AGE.CLASS,
                              random = ~ YEAR 
                              + female.ID
                              + male.ID
                              + str(~vm(female.animal, ainv) 
                                    + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                              + ar1v(col):ar1(row),
                              residual = ~idv(units),
                              data = phenotypic.data,
                              workspace = 8e8,  # Boost memory
                              maxit = 100)
summary(SA.clutchsize.model)$varcomp[,c("component", "std.error", "bound")]
summary(SA.clutchsize.model, coef=T)$coef.fixed
wald.asreml(SA.clutchsize.model, ssType = "conditional", denDF = "numeric")

SA.CS.Vp       <- vpredict(SA.clutchsize.model, SA.CS.Vphenotypic ~ V2+V3+V5+V6+V8+V10)  # Year excluded from phenotypic variance
SA.CS.y2       <- vpredict(SA.clutchsize.model, SA.CS.y2 ~ V1 / (V2+V3+V5+V6+V8+V10))
SA.CS.f.pe2    <- vpredict(SA.clutchsize.model, SA.CS.f.pe2 ~ V2 / (V2+V3+V5+V6+V8+V10))
SA.CS.m.pe2    <- vpredict(SA.clutchsize.model, SA.CS.m.pe2 ~ V3 / (V2+V3+V5+V6+V8+V10))
SA.CS.f.h2     <- vpredict(SA.clutchsize.model, SA.CS.f.h2 ~ V5 / (V2+V3+V5+V6+V8+V10))
SA.CS.m.h2     <- vpredict(SA.clutchsize.model, SA.CS.m.h2 ~ V6 / (V2+V3+V5+V6+V8+V10))
SA.CS.cov.a    <- vpredict(SA.clutchsize.model, SA.CS.cov.a ~ V4 * sqrt(V5*V6))
SA.CS.tbv      <- vpredict(SA.clutchsize.model, SA.CS.tbv ~ V5 + 2*(V4*sqrt(V5*V6)) + V6)
SA.CS.total.h2 <- vpredict(SA.clutchsize.model, SA.CS.total.h2 ~ (V5+2*(V4*sqrt(V5*V6))+V6)/ (V2+V3+V5+V6+V8+V10))
SA.CS.sa2      <- vpredict(SA.clutchsize.model, SA.CS.sa2 ~ V8/(V2+V3+V5+V6+V8+V10))
SA.CS.r2       <- vpredict(SA.clutchsize.model, SA.CS.r2 ~ V10 / (V2+V3+V5+V6+V8+V10))

SA.CS.Vp
SA.CS.y2
SA.CS.f.pe2
SA.CS.m.pe2
SA.CS.f.h2
SA.CS.m.h2
SA.CS.cov.a
SA.CS.tbv
SA.CS.total.h2
SA.CS.sa2
SA.CS.r2

# TESTING VARIANCE COMPONENTS
# Year
SA.clutchsize.model_no.year <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                      + M.AGE.CLASS,
                                      random = ~ female.ID
                                      + male.ID
                                      + str(~vm(female.animal, ainv) 
                                            + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                                      + ar1v(col):ar1(row),
                                      residual = ~idv(units),
                                      data = phenotypic.data,
                                      workspace = 8e8,  # Boost memory
                                      maxit = 30,
                                      na.action = na.method(x = "include", y = "omit"))
SA.clutchsize.model_no.year <- update(SA.clutchsize.model_no.year)
chi.sq <- -2*(SA.clutchsize.model_no.year$loglik - SA.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Female permanent environment
SA.clutchsize.model_no.fPE <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                     + M.AGE.CLASS,
                                     random = ~ YEAR 
                                     + male.ID
                                     + str(~vm(female.animal, ainv) 
                                           + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                                     + ar1v(col):ar1(row),
                                     residual = ~idv(units),
                                     data = phenotypic.data,
                                     workspace = 8e8,  # Boost memory
                                     maxit = 30,
                                     na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.fPE$loglik - SA.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male permanent environment
SA.clutchsize.model_no.mPE <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                     + M.AGE.CLASS,
                                     random = ~ YEAR 
                                     + female.ID
                                     + str(~vm(female.animal, ainv) 
                                           + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv))
                                     + ar1v(col):ar1(row),
                                     residual = ~idv(units),
                                     data = phenotypic.data,
                                     workspace = 8e8,  # Boost memory
                                     maxit = 30,
                                     na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.mPE$loglik - SA.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Intersexaul additive genetic covariance
SA.clutchsize.model_no.COVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                      + M.AGE.CLASS,
                                      random = ~ YEAR 
                                      + female.ID
                                      + male.ID
                                      + vm(female.animal, ainv)
                                      + vm(male.animal, ainv)
                                      + ar1v(col):ar1(row),
                                      residual = ~idv(units),
                                      data = phenotypic.data,
                                      workspace = 8e8,  # Boost memory
                                      maxit = 30,
                                      na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.COVa$loglik - SA.clutchsize.model$loglik)  
chi.sq  # Test statistic
1-pchisq(chi.sq, 1)  # P-value

# Female additive genetic
SA.clutchsize.model_no.fVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                     + M.AGE.CLASS,
                                     random = ~ YEAR 
                                     + female.ID
                                     + male.ID
                                     + vm(male.animal, ainv)
                                     + ar1v(col):ar1(row),
                                     residual = ~idv(units),
                                     data = phenotypic.data,
                                     workspace = 8e8,  # Boost memory
                                     maxit = 30,
                                     na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.fVa$loglik - SA.clutchsize.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# Male additive genetic
SA.clutchsize.model_no.mVa <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                     + M.AGE.CLASS,
                                     random = ~ YEAR 
                                     + female.ID
                                     + male.ID
                                     + vm(female.animal, ainv)
                                     + ar1v(col):ar1(row),
                                     residual = ~idv(units),
                                     data = phenotypic.data,
                                     workspace = 8e8,  # Boost memory
                                     maxit = 30,
                                     na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.mVa$loglik - SA.clutchsize.model_no.COVa$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 0, 1)  # P-value

# G-matrix
SA.clutchsize.model_no.G <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                   + M.AGE.CLASS,
                                   random = ~ YEAR 
                                   + female.ID
                                   + male.ID
                                   + ar1v(col):ar1(row),
                                   residual = ~idv(units),
                                   data = phenotypic.data,
                                   workspace = 8e8,  # Boost memory
                                   maxit = 30,
                                   na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.G$loglik - SA.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 1, 3)  # P-value

# Spatial autocorrelation
SA.clutchsize.model_no.SA <- asreml(fixed = CLUTCHSIZE ~ F.AGE.CLASS
                                    + M.AGE.CLASS,
                                    random = ~ YEAR 
                                    + female.ID
                                    + male.ID
                                    + str(~vm(female.animal, ainv) 
                                          + vm(male.animal, ainv), ~corgh(2):vm(female.animal, ainv)),
                                    residual = ~idv(units),
                                    data = phenotypic.data,
                                    workspace = 8e8,  # Boost memory
                                    maxit = 30,
                                    na.action = na.method(x = "include", y = "omit"))
chi.sq <- -2*(SA.clutchsize.model_no.SA$loglik - SA.clutchsize.model$loglik)  
chi.sq  # Test statistic
p.chimixer(chi.sq, 2, 3)  # P-value