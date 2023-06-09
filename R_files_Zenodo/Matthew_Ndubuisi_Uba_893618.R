# This is the code for the analysis of the data set "Kebbi.csv" for the 
# Environmental Data Analysis project (ENVDA -MOD 1) titled Deforestation:
# An Implication of High Rate of Unemployed Household in Rural Communities 
# in Kebbi State, Nigeria.

# I first set the working directory: Click on session tab, select working directory,
# choose working directory and click on the folder where I save the data set 
# "Kebbi.csv".

# Reading the data set "Kebbi.csv"
Kebbitable<-read.csv("Kebbi.csv")
head(Kebbitable)

# Check is there are missing value in the data set "Kebbi.csv".
summary(Kebbitable)

# Omitting missing values from the data set "Kebbi.csv".
Kebbitable.1<-na.omit(Kebbitable)
summary(Kebbitable.1)
View(Kebbitable.1)

attach(Kebbitable.1)


par(mfcol=c(1,2))

# Get the frequency of employed members of the households in Kebbi
hist(Employed_members,freq = TRUE,col = "red",main = paste("Employed Members"))
abline(v=mEmploy<-mean(Employed_members),lwd=2,lty="dashed")

# Get the frequency of the annual income of households in Kebbi
hist(Annual_income,freq = TRUE,col = "red",main = paste("Annual Income"))
abline(v=mIncome<-mean(Annual_income),lwd=2,lty="dashed")

par(mfcol=c(1,2))

# Get the frequency of the different sources of lighting fuel in households in Kebbi
hist(Lighting_fuel_source,freq = TRUE,col = "red",main = paste("Sources of Lighting Fuel"))
abline(v=mLight<-mean(Lighting_fuel_source),lwd=2,lty="dashed")


# Get the frequency of the different sources of cooking fuel in households in Kebbi
hist(Cooking_fuel_source,freq = TRUE,col = "red",main = paste("Sources of Cooking Fuel"))
abline(v=mCook<-mean(Cooking_fuel_source),lwd=2,lty="dashed")

par(mfcol=c(1,2))

# Calculating the mean of the employed members of the households in Kebbi.
boxplot(Employed_members, main = "Mean of Employed Members", horizontal = FALSE,col = "red")
mEmploy

# Calculating the mean of annual income of the households in Kebbi.
boxplot(Annual_income, main = "Mean of Annual Income", horizontal = FALSE,col = "red")
mIncome

par(mfcol=c(1,2))

# Calculating the mean of sources of lighting fuel of the households in Kebbi.
boxplot(Lighting_fuel_source, main = "Mean of Sources of Lighting Fuel", horizontal = FALSE,col = "red")
mLight

# Calculating the mean of sources of cooking fuel of the households in Kebbi.
boxplot(Cooking_fuel_source, main = "Mean of Sources of Cooking Fuel", horizontal = FALSE,col = "red")
mCook


# Hypothesis Testing
levels(Lighting_fuel_source)<-c("collected firewood","purchased firewood","grass","kerosene","phcn electricity","generator","gas","battery","candle")
levels(Cooking_fuel_source)<-c("collected firewood","purchased firewood","coal","grass","kerosene","phcn electricity","generator","gas","others")
levels(Annual_income)<-c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15","N16","N17","N18","N19","20")


## Testing for Null Hypothesis H0 (1a) The member of household that has worked in a business or self-employment in the past 12 months has no effect on the source of lighting fuel 
boxplot(Employed_members~Lighting_fuel_source,main="emplyed members vs lighting fuel",cex.lab=1.5,col="lightblue")
summary(aov(Employed_members~Lighting_fuel_source))

## Testing for Null Hypothesis H0 (1b) The member of household that has worked in a business or self-employment in the past 12 months has no effect on the source of cooking fuel 
boxplot(Employed_members~Cooking_fuel_source,main="emplyed members vs cooking fuel",cex.lab=1.5,col="lightblue")
summary(aov(Employed_members~Cooking_fuel_source))

# Testing for Null Hypothesis Ho (2) The member of household that has worked in a business or self-employment in the past 12 months has no effect on the number of naira for annual income of the member of the household 
boxplot(Employed_members~Annual_income,main="employed members vs annual income",cex.lab=1.5,col="lightblue")
summary(aov(Employed_members~Annual_income))

## Testing for Null Hypothesis H0 (3a) The number of naira for annual income of the member of the household has no effect on the source of lightning and cooking fuel 
boxplot(Annual_income~Lighting_fuel_source,main="annual income vs lighting fuel",cex.lab=1.5,col="lightblue")
summary(aov(Annual_income~Lighting_fuel_source))

## Testing for Null Hypothesis H0 (3b) The number of naira for annual income of the member of the household has no effect on the source of lightning and cooking fuel 
boxplot(Annual_income~Cooking_fuel_source,main="annual income vs cooking fuel",cex.lab=1.5,col="lightblue")
summary(aov(Annual_income~Cooking_fuel_source))
