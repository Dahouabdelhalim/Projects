# Required libraries
library(xlsx)
library(car)
library(AER)
library(MASS)
library(ggplot2)

# Load Excel file directly or use file.choose()
myDataFC <- read.xlsx("StageAppearanceData.xlsx", sheetName="FirstCauliflower")

# Define HostAge as a fixed factor
myDataFC$HostAge <- factor(myDataFC$HostAge)
str(myDataFC)

# Fit generalized linear model (Poisson distribution)
fitFC <- glm(FirstAppearance ~ HostAge, family="poisson", data=myDataFC)
summary(fitFC)
Anova(fitFC)

# Confirm lack of overdispersion
dispersiontest(fitFC)

# Display box plot
boxPlotFC <- ggplot(myDataFC, aes(x=HostAge, y=FirstAppearance))
boxPlotFC + geom_boxplot()

# Load Excel file directly or use file.choose()
myDataGS <- read.xlsx("StageAppearanceData.xlsx", sheetName="GrapeSeed")

# Define HostAge as a fixed factor
myDataGS$HostAge <- factor(myDataGS$HostAge)
str(myDataGS)

# Fit generalized linear model (Poisson distribution)
fitGS <- glm(FirstAppearance ~ HostAge, family="poisson", data=myDataGS)
summary(fitGS)
Anova(fitGS)

# Confirm lack of overdispersion
dispersiontest(fitGS)

# Fit generalized linear model (Negative Binomial distribution) due to overdispersion
fitGS <- glm.nb(FirstAppearance ~ HostAge, data=myDataGS)
summary(fitGS)
Anova(fitGS)

# Display box plot
boxPlotGS <- ggplot(myDataGS, aes(x=HostAge, y=FirstAppearance))
boxPlotGS + geom_boxplot()

# Load Excel file directly or use file.choose()
myDataMS <- read.xlsx("StageAppearanceData.xlsx", sheetName="MatureSpores")

# Define HostAge as a fixed factor
myDataMS$HostAge <- factor(myDataMS$HostAge)
str(myDataMS)

# Fit generalized linear model (Poisson distribution)
fitMS <- glm(FirstAppearance ~ HostAge, family="poisson", data=myDataMS)
summary(fitMS)
Anova(fitMS)

# Confirm lack of overdispersion
dispersiontest(fitMS)

# Display box plot
boxPlotMS <- ggplot(myDataMS, aes(x=HostAge, y=FirstAppearance))
boxPlotMS + geom_boxplot()

# Load Excel file directly or use file.choose()
myDataSC <- read.xlsx("StageAppearanceData.xlsx", sheetName="SecondCauliflower")

# Define HostAge as a fixed factor
myDataSC$HostAge <- factor(myDataSC$HostAge)
str(myDataSC)

# Fit generalized linear model (Poisson distribution)
fitSC <- glm(FirstAppearance ~ HostAge, family="poisson", data=myDataSC)
summary(fitSC)
Anova(fitSC)

# Confirm lack of overdispersion
dispersiontest(fitSC)

# Display box plot
boxPlotSC <- ggplot(myDataSC, aes(x=HostAge, y=FirstAppearance))
boxPlotSC + geom_boxplot()
