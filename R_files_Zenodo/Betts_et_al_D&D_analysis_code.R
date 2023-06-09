# Analysis code for:
# --- 'Old-growth forests buffer climate-sensitive bird populations from warming' ---
# by Betts MG, Phalan B, Frey SJK, Rousseau JS, Yang Z
# Published in Diversity & Distributions December 2017

setwd("This will depend on your computer and where you put the files")

# load libraries
library(stringi)
library(MuMIn)

# import and subset data

# BBS bird trend data
BBS_birdtrends <- read.csv("BBS_birdtrends.csv")

# subset to routes that start before 2006
BBS_birdtrends06 <- subset(BBS_birdtrends, Start<2006)

# import predictors (climate, old-growth forest, elevation, lat/long, forest cover)
predictors <- read.csv("predictors.csv")
head(predictors)


# -- Run models and create AIC table for choosing the month of temperature to use for each species --

# species list (common name)
species <- c("Brown Creeper", "Chestnut-backed Chickade", "Golden-crowned Kinglet", "Hammond's Flycatcher", "Hermit Thrush", "Pacific Wren", "Cordilleran & Pacific-sl", "Red-breasted Sapsucker", "Townsend's Warbler", "Varied Thrush", "Vaux's Swift", "Wilson's Warbler", "Hermit Warbler")
length(species)

# Create temporary matrix to fill with AIC values to pick month of temperature for each species
AIC_temp <- matrix(NA, 13, 14)

colnames(AIC_temp) <- c("JuneAIC", "JulyAIC", "JuneEst", "JuneSE",  "JulyEst", "JulySE", "JuneP", "JulyP", "JuneCI", "JulyCI", "monthAIC", "AICdif", "monthP", "monthCI")
rownames(AIC_temp) <- species
AIC_temp

AIC_temp <- as.data.frame(AIC_temp)

# For loop for running the models for each fpes 
for (s in species) {
	
	sp <- subset (BBS_birdtrends06, SpeciesCommon == s) #Select species by code
	pred_sp <- merge(predictors, sp, by.x = "Route") #merge with climate data
	bbs_pred <- subset(pred_sp, EEQ_Variance!= "NaN") #remove NAs

# subset to EEQ_NumYrs > 6
	trend_pred <- subset(bbs_pred, EEQ_NumYrs > 6)

	trend <- trend_pred$EEQ_TrendEst
	
# pick June or July
	tempJun <- trend_pred$JunTMax
	tempJul <- trend_pred$JulTMax
	
# old growth
	OG <- trend_pred$ogsi

# Compute the linear regression (scale variables)
	fitJun <- glm(trend~scale(tempJun)*scale(OG), weights=1/trend_pred$EEQ_Variance^2)
	fitJul <- glm(trend~scale(tempJul)*scale(OG), weights=1/trend_pred$EEQ_Variance^2)

	#fitJul <- glm(trend~tempJul*OG, weights=1/trend_pred$EEQ_Variance^2)

# save the models with species names
	m <- stri_replace_all_fixed(s, " ", "")
	m <- stri_replace_all_fixed(m, "-", "")
	m <- stri_replace_all_fixed(m, "'", "")
	m <- stri_replace_all_fixed(m, "&", "")

	assign(paste(m, "fitJun", sep=""), fitJun)
	assign(paste(m, "fitJul", sep=""), fitJul)

# fill in table
	AIC_temp[s, "JuneAIC"] <- as.numeric(fitJun$aic)
	AIC_temp[s, "JuneEst"] <- coef(fitJun)[2]
	AIC_temp[s, "JuneSE"] <- coef(summary(fitJun))[2,2]
	AIC_temp[s, "JuneP"] <- coef(summary(fitJun))[2,4]

# calculate CIs
	AIC_temp[s, "JuneCI"] <- abs(confint(fitJun,2)[1])+abs(confint(fitJun,2)[2])


# fill in table
	AIC_temp[s, "JulyAIC"] <- as.numeric(fitJul$aic)
	AIC_temp[s, "JulyEst"] <- coef(fitJul)[2]
	AIC_temp[s, "JulySE"] <- coef(summary(fitJul))[2,2]
	AIC_temp[s, "JulyP"] <- coef(summary(fitJul))[2,4]
	
	AIC_temp[s, "JulyCI"] <- abs(confint(fitJul,2)[1])+abs(confint(fitJul,2)[2])

	AIC_temp[s, "AICdif"] <- abs(as.numeric(fitJun$aic) - as.numeric(fitJul$aic))
	
# pick month with lowest AIC
	if (fitJun$aic > fitJul$aic)
	{AIC_temp[s, "monthAIC"] <- "July"
		} else {AIC_temp[s, "monthAIC"] <- "June"}
		
# pick month with smallest P-value
	if (coef(summary(fitJun))[2,4] > coef(summary(fitJul))[2,4])
	{AIC_temp[s, "monthP"] <- "July"
		} else {AIC_temp[s, "monthP"] <- "June"}
	
	CIjune <- abs(confint(fitJun,2)[1])+abs(confint(fitJun,2)[2])
	CIjuly <- abs(confint(fitJul,2)[1])+abs(confint(fitJul,2)[2])
	
# pick month with smallest CI
	if (CIjune > CIjuly)
	{AIC_temp[s, "monthCI"] <- "July"
		} else {AIC_temp[s, "monthCI"] <- "June"}

}
	
AIC_temp
str(AIC_temp)

# export table
getwd()
write.csv(AIC_temp, file="AIC_temp.csv", row.names=TRUE)


# --- ASSIGN MONTH FOR EACH SPECIES

# Brown creeper *JULY*
# June
summary(BrownCreeperfitJun)
confint(BrownCreeperfitJun)
# July
summary(BrownCreeperfitJul)
confint(BrownCreeperfitJul)

BRCRJuly <- BrownCreeperfitJul


# Chestnut-backed Chickadee *JUNE*
# June
summary(ChestnutbackedChickadefitJun)
# July
summary(ChestnutbackedChickadefitJul)

CBCHJune <- ChestnutbackedChickadefitJun

# Golden-crowned Kinglet *JULY*
# June
summary(GoldencrownedKingletfitJun)
confint(GoldencrownedKingletfitJun)
# July
summary(GoldencrownedKingletfitJul)
confint(GoldencrownedKingletfitJul)

GCKIJuly <- GoldencrownedKingletfitJul

# Hammond's Flycatcher *JUNE*
# June
summary(HammondsFlycatcherfitJun)
# July
summary(HammondsFlycatcherfitJul)

HAFLJune <- HammondsFlycatcherfitJun

# Hermit Thrush *JUNE*
# June
summary(HermitThrushfitJun)
# July
summary(HermitThrushfitJul)

HETHJune <- HermitThrushfitJun

# Pacific Wren *JULY*
# June
summary(PacificWrenfitJun)
confint(PacificWrenfitJun)
# July
summary(PacificWrenfitJul)
confint(PacificWrenfitJul)

PAWRJuly <- PacificWrenfitJul

# Pacific-slope Flyccatcher *JUNE*
# June
summary(CordilleranPacificslfitJun)
confint(CordilleranPacificslfitJun)
# July
summary(CordilleranPacificslfitJul)
confint(CordilleranPacificslfitJul)

PSFLJune <- CordilleranPacificslfitJun

# Red-breasted Sapsucker *JUNE*
# June
summary(RedbreastedSapsuckerfitJun)
# July
summary(RedbreastedSapsuckerfitJul)

RBSAJune <- RedbreastedSapsuckerfitJun

# Townsend's Warbler *JUNE*
# June
summary(TownsendsWarblerfitJun)
confint(TownsendsWarblerfitJun)
# July
summary(TownsendsWarblerfitJul)
confint(TownsendsWarblerfitJul)

TOWAJune <- TownsendsWarblerfitJun

# Varied Thrush *JULY*
# June
summary(VariedThrushfitJun)
# July
summary(VariedThrushfitJul)

VATHJuly <- VariedThrushfitJul

# Vaux's Swift *JULY*
# June
summary(VauxsSwiftfitJun)
# July
summary(VauxsSwiftfitJul)

VASWJuly <- VauxsSwiftfitJul


# Wilson's Warbler *JUNE*
# June
summary(WilsonsWarblerfitJun)
confint(WilsonsWarblerfitJun)
# July
summary(WilsonsWarblerfitJul)
confint(WilsonsWarblerfitJul)

WIWAJune <- WilsonsWarblerfitJun


# Hermit Warbler *JULY*
# June
summary(HermitWarblerfitJun)
# July
summary(HermitWarblerfitJul)

HEWAJuly <- HermitWarblerfitJul

# Table 1 for temp*OG models (using the best month for each species)

# List of species and months
sp_month <- c("BRCRJuly", "CBCHJune", "GCKIJuly", "HAFLJune", "HETHJune", "HEWAJuly", "PAWRJuly", "PSFLJune", "RBSAJune", "TOWAJune", "VATHJuly", "VASWJuly", "WIWAJune")

# Create temporary matrix to fill with parameter estimates and SE for each species' best model
temp_month <- matrix(NA, 13, 18)

colnames(temp_month) <- c("Species", "SpCode", "Month", "Temp", "TempSE", "TempLCL", "TempUCL", "TempP",  "OG", "OGSE", "OGLCL", "OGUCL", "OGP", "TempOG", "TempOGSE", "TempOGLCL", "TempOGUCL", "TempOGP")
rownames(temp_month) <- sp_month
temp_month

temp_month <- as.data.frame(temp_month)

temp_month$Species <- c("Brown Creeper", "Chestnut-backed Chickadee", "Golden-crowned Kinglet", "Hammond's Flycatcher", "Hermit Thrush", "Hermit Warbler", "Pacific Wren", "Pacific-slope Flycatcher", "Red-breasted Sapsucker", "Townsend's Warbler", "Varied Thrush", "Vaux's Swift", "Wilson's Warbler")

temp_month$SpCode <- c("BRCR", "CBCH", "GCKI", "HAFL", "HETH", "HEWA", "PAWR", "PSFL", "RBSA", "TOWA", "VATH", "VASW", "WIWA")

temp_month$Month <- c("July", "June", "July", "June", "June", "July", "June", "June", "June", "July", "July", "June", "July")

# look at dataframe 
temp_month

# For loop 

for (m in sp_month){

model <- get(m)

# Temp
temp_month[m, "Temp"]  <- coef(summary(model))[2,1]
temp_month[m, "TempSE"]  <- coef(summary(model))[2,2]
temp_month[m, "TempLCL"]  <- confint(model)[2,1]
temp_month[m, "TempUCL"]  <- confint(model)[2,2]
temp_month[m, "TempP"]  <- coef(summary(model))[2,4]

# OG
temp_month[m, "OG"]  <- coef(summary(model))[3,1]
temp_month[m, "OGSE"]  <- coef(summary(model))[3,2]
temp_month[m, "OGLCL"]  <- confint(model)[3,1]
temp_month[m, "OGUCL"]  <- confint(model)[3,2]
temp_month[m, "OGP"]  <- coef(summary(model))[3,4]

# OGxTemp
temp_month[m, "TempOG"]  <- coef(summary(model))[4,1]
temp_month[m, "TempOGSE"]  <- coef(summary(model))[4,2]
temp_month[m, "TempOGLCL"]  <- confint(model)[4,1]
temp_month[m, "TempOGUCL"]  <- confint(model)[4,2]
temp_month[m, "TempOGP"]  <- coef(summary(model))[4,4]	
	
}

temp_month

# export table
getwd()
write.csv(temp_month, file="TempOG_models.csv", row.names=TRUE)




# --- Global models for species with significant Temp*OG interactions --- (for supplementary tables S1 & S2)


# HEWA
# global model
spHEWA <- subset(BBS_birdtrends06, SpeciesCommon == "Hermit Warbler") #Select species by code
bbs_spHEWA <- merge(predictors, spHEWA, by.x = "Route") #merge with climate data
bbs_spHEWA <- subset(bbs_spHEWA, EEQ_Variance!= "NaN") #remove NAs
# subset to EEQ_NumYrs > 6
bbs_spHEWA <- subset(bbs_spHEWA, EEQ_NumYrs > 6)
head(bbs_spHEWA)

# *MODEL*
HEWAglobal <- glm(EEQ_TrendEst ~ scale(JulTMax)*scale(lat) + scale(JulTMax)*scale(mean_elev) + scale(JulTMax)*scale(stdDev) + scale(JulTMax)*scale(Av_for_cover) + scale(JulTMax)*scale(ogsi), weights = 1/EEQ_Variance^2, na.action = "na.fail", data=bbs_spHEWA)
summary(HEWAglobal)

# model selection
aHEWA <- dredge(HEWAglobal)
model.avg(aHEWA)
importance(aHEWA)
summary(model.avg(aHEWA, subset = delta < 4))
HEWAma <- summary(model.avg(aHEWA, subset = delta < 4))

# export results
str(HEWAma)
HEWAma4 <- as.data.frame(HEWAma$coefmat.subset)
HEWAma4
write.csv(HEWAma4, "HEWAglobal_ma_results.csv", row.names=TRUE)



# WIWA
# global model
spWIWA <- subset(BBS_birdtrends06, SpeciesCommon == "Wilson's Warbler") #Select species by code
bbs_spWIWA <- merge(predictors, spWIWA, by.x = "Route") #merge with climate data
bbs_spWIWA <- subset(bbs_spWIWA, EEQ_Variance!= "NaN") #remove NAs
# subset to EEQ_NumYrs > 6
bbs_spWIWA <- subset(bbs_spWIWA, EEQ_NumYrs > 6)
head(bbs_spWIWA)


# *MODEL*
WIWAglobal <- glm(EEQ_TrendEst ~ scale(JunTMax)*scale(lat) + scale(JunTMax)*scale(mean_elev) + scale(JunTMax)*scale(stdDev) + scale(JunTMax)*scale(Av_for_cover) + scale(JunTMax)*scale(ogsi), weights = 1/EEQ_Variance^2, na.action = "na.fail", data=bbs_spWIWA)
summary(WIWAglobal)

# model selection
aWIWA <- dredge(WIWAglobal)
model.avg(aWIWA)
importance(aWIWA)
summary(model.avg(aWIWA, subset = delta < 4))
WIWAma <- summary(model.avg(aWIWA, subset = delta < 4))

# export results
str(WIWAma)
WIWAma4 <- as.data.frame(WIWAma$coefmat.subset)
WIWAma4
write.csv(WIWAma4, "WIWAglobal_ma_results.csv", row.names=TRUE)

# % change per year with 1degC

# HEWA 
trendHEWA <- bbs_spHEWA$EEQ_TrendEst
tempJulHEWA <- bbs_spHEWA$JulTMax
ogHEWA <- bbs_spHEWA$ogsi

JunHEWA <- lm(trendHEWA~tempJulHEWA*ogHEWA, weights=1/bbs_spHEWA$EEQ_Variance^2)
summary(JunHEWA)
confint(JunHEWA)

# WIWA 
trendWIWA <- bbs_spWIWA$EEQ_TrendEst
tempJunWIWA <- bbs_spWIWA$JunTMax
ogWIWA <- bbs_spWIWA$ogsi

JunWIWA <- lm(trendWIWA~tempJunWIWA*ogWIWA, weights=1/bbs_spWIWA$EEQ_Variance^2)
summary(JunWIWA)
confint(JunWIWA)



