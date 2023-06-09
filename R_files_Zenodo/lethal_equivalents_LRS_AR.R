#Script for estimating lethal equivalents for AR and LRS in Helgeland house sparrows 
#Alina Niskanen & Stefanie Muff
#alina.niskanen@gmail.com
#April 2020


library(INLA)
library(MCMCpack)
library(MasterBayes)
library(gridGraphics)
library(ggplot2)
library(nadiv)
library(MCMCglmm)
library(pedigreemm)
library(pryr)
library(cowplot)

#Store island names and codes into a data frame to use in plotting
Island_names <- as.data.frame(c("Nesøy","Myken","Træna","Selvær","Gjerøy","Hestmannøy","Indre Kvarøy","Aldra"))
Island_names$code <- c("20","22","23","24","26","27","28","38")
colnames(Island_names)[1] <- "Island"

#Import the phenotypic data and inbreeding estimates
data_LRS <- read.csv("Data.txt", header = T, stringsAsFactors = F, sep = "\\t")

###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_LRS[which(data_LRS$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)


###
###Lifetime recruit production
###

#Only those individuals that have LRS_OK estimate and have been adults on one of the study islands
data_ok <- data_adult_temp[!is.na(data_adult_temp$LRS_OK),] 


#Remove temp files
rm(list=ls(pattern="temp"))

###
###Prepare the variables for the dataset that has correct LRS data for adult islands
###

#Scale FROH
data_ok$z_FROH <- scale(data_ok$FROH, center = T, scale = T)

# make island factor covariate
data_ok$f_laflok <- as.factor(data_ok$laflok)

#Add combined variable for hatch year and adult island
data_ok$IAY <- as.factor(paste(data_ok$all_hatchyears,data_ok$laflok,sep="_"))

#Import the pedigree
ped <- read.table("pedigree.txt", header = T, stringsAsFactors = F)

# need an ordered Pedigree for the inverseA() function:
ped <- orderPed(ped)

# introduce the "ID", a new identity for individuals, enumerated from 1 to number of ind.:
ped$ID <- 1:(nrow(ped))

# fathers and mothers are replaced by these new IDs
d.map <- ped[,c("id","ID")]
ped$dam.id <- d.map[match(ped$dam, d.map$id),"ID"]
ped$sire.id <- d.map[match(ped$sire, d.map$id),"ID"]


# compute A inverse, using the new IDs:
Cmatrix <- inverseA(ped[,c("ID","dam.id","sire.id")])$Ainv

# Update ID column of the data file:
data_ok$ID <- d.map[match(data_ok$id, d.map$id), "ID"]

#Make factors out of necessary variables
data_ok$f_laflok2 <- as.factor(data_ok$f_laflok2)
data_ok$IAY2 <- as.factor(data_ok$IAY2) 

#order factor f_laflok2 to get the islands in logical order from INLA sampling output
data_ok$f_laflok2 = factor(data_ok$f_laflok2,levels=c("20","22","23","24","26","27","28","38"))

##########
#Lifetime reproductive success
##########

###Estimate lethal equivalents using standardised FROH, i.e. the model used in the main inbreeding depression analysis

#Including FROH:laflok and FROH:year_island interactions in the random term, FROH scaled
formula_7.1_FROH = LRS_OK ~ z_FROH + gen_sex +  
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2,z_FROH,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FROH, model = "iid", hyper=list( prec=list(initial=log(36), prior="pc.prec",param=c(0.1,0.01)))) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  ) 

#Run INLA
model.zip_7.1_FROH = inla(formula=formula_7.1_FROH, family="zeroinflatedpoisson1",
                          data=data_ok,
                          control.family=list(link='log'),
                          control.compute=list(dic=T, config=T), 
                          control.predictor=list(link=1, compute=TRUE)
)
summary(model.zip_7.1_FROH)


#Estimate lethal equivalents per island

###Extracting island specific inbreeding effects
### Sample from posterior distribution
nsamples <- 1000
model <- model.zip_7.1_FROH
sample <- inla.posterior.sample(n=nsamples,model, add.names = TRUE, use.improved.mean=TRUE, intern=FALSE)
rm(model)

#Take the posterior values for FROH effect and each island effect separately and restore them in a data frame
island_samp <- as.data.frame(rep(NA,nsamples))

for (i in 1:nsamples){
  island_samp[i,1] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:1",1]
  island_samp[i,2] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:2",1]
  island_samp[i,3] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:3",1]
  island_samp[i,4] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:4",1]
  island_samp[i,5] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:5",1]
  island_samp[i,6] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:6",1]
  island_samp[i,7] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:7",1]
  island_samp[i,8] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:8",1]
}

###Effect sizes
mean_effect_LRS_FROH<-apply(island_samp,2,mean)
lc_LRS_FROH<-apply(island_samp,2,quantile, 0.025)
uc_LRS_FROH<-apply(island_samp,2,quantile, 0.975)

par(mai=c(1.02,1.15,0.82,0.42))
plot(c(1:8)~mean_effect_LRS_FROH, xlim=c(-1,0.5), col=rainbow(8), lwd=3, yaxt="n", ylab="", xlab="Effect size of standardised FROH", main="Effect of inbreeding on lifetime reproductive success")
abline(v=0, lty=3)
segments(uc_LRS_FROH, 1:8, lc_LRS_FROH, 1:8, col=rainbow(8), lwd=3)
axis(2, c(1:8), Island_names$Island, las=2)

#Estimate lethal equivalents from the standardised FROH results
#Estimate lethal equivalents and their upper and lower confidence limits using the mean of the posterior sample for each island
LE_LRS <- as.data.frame(apply(island_samp,2,mean))
LE_LRS[,2:7] <- NA 
std_Froh_LRS <- sd(data_ok$FROH)
colnames(LE_LRS) <- c("mean_effect_size", "lc_LRS_FROH","uc_LRS_FROH","scaled_LE", "lethal_equivalents","lc_LRS_lethal_equivalents","uc_LRS_lethal_equivalents")
rownames(LE_LRS) <- c("Nesøy","Myken","Træna","Selvær","Gjerøy","Hestmannøy","Indre Kvarøy","Aldra")
LE_LRS$scaled_LE <- -LE_LRS$mean_effect_size*2
LE_LRS$lethal_equivalents <- -LE_LRS$mean_effect_size*2/std_Froh_LRS
LE_LRS$lc_LRS_FROH<- apply(island_samp,2,quantile, 0.025)
LE_LRS$lc_LRS_lethal_equivalents <- (LE_LRS$scaled_LE-((LE_LRS$mean_effect_size-LE_LRS$lc_LRS_FROH)*2))/std_Froh_LRS
LE_LRS$uc_LRS_FROH <- apply(island_samp,2,quantile, 0.975)
LE_LRS$uc_LRS_lethal_equivalents <- (((LE_LRS$uc_LRS_FROH-LE_LRS$mean_effect_size)*2)+LE_LRS$scaled_LE)/std_Froh_LRS

#Write table of lethal equivalents in LRS per island
write.table(LE_LRS,"lethal_equivalents_LRS_z_FROH.txt", col.names = T, row.names = T, quote = F)


rm(sample)
rm(island_samp)



##########
#Yearly reproductive success
##########


#Import the phenotypic data and inbreeding estimates
data_LRS <- read.csv("Data.txt", header = T, stringsAsFactors = F, sep = "\\t")

###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_LRS[which(data_LRS$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)

###
###Per year recruit production
###

#Start choosing the correct data from "data_adult_temp" that includes all individuals that were observed as adults (laflok) on one of the 8 study islands
data_no_YR <- data_adult_temp[(is.na(data_adult_temp$X1998) & is.na(data_adult_temp$X1999) & is.na(data_adult_temp$X2000) & is.na(data_adult_temp$X2001) & is.na(data_adult_temp$X2002) & is.na(data_adult_temp$X2003) & is.na(data_adult_temp$X2004) & is.na(data_adult_temp$X2005) & is.na(data_adult_temp$X2006) & is.na(data_adult_temp$X2007) & is.na(data_adult_temp$X2008) & is.na(data_adult_temp$X2009) & is.na(data_adult_temp$X2010) & is.na(data_adult_temp$X2011) & is.na(data_adult_temp$X2012)),]

#Make the final yearly estimate dataset, includes only individuals with laflok on one of the 8 study islands and YR estimate available
data_YR_temp <- data_adult_temp[!data_adult_temp$id %in% data_no_YR$id,]

#Scale FROH
data_YR_temp$z_FROH <- scale(data_YR_temp$FROH, center = T, scale = T)

#Transpose the data to include each individual yearly estimate on its own row
data_YR <- gather(data_YR_temp, "Year", "recruits", 22:36)
data_YR$Year <- gsub("X","", data_YR$Year)

#Remove all rows with missing YR to reduce the file size
data_YR <- data_YR[!is.na(data_YR$recruits),]
table(data_YR$recruits)
hist(data_YR$recruits)
#Add age of the individual when having the offspring
data_YR$Year <- as.integer(data_YR$Year)
data_YR$age <- data_YR$Year-data_YR$all_hatchyears

#Remove temp files
rm(list=ls(pattern="temp"))


###
###Prepare the variables for the dataset that has correct YR data for adult islands
###

# make island factor
data_YR$f_laflok <- as.factor(data_YR$laflok)

#Combine age classes 6 and above
data_YR$age_comb <- data_YR$age
data_YR$age_comb[which(data_YR$age>=6)] <- 6
#Center combined age
data_YR$c_age_comb <- scale(data_YR$age_comb, center = T, scale = F)
#add squared age for quadratic effect
data_YR$c_age_comb2 <- data_YR$c_age_comb^2

#Add combined variable for hatch year and adult island
data_YR$IAY <- as.factor(paste(data_YR$all_hatchyears,data_YR$laflok,sep="_"))
#Add combined variable for hatch year and adult island to be used in
data_YR$IAY2 <- data_YR$IAY

# also need to add ID column to data file:
data_YR$ID <- d.map[match(data_YR$id, d.map$id), "ID"]
data_YR$ID2 <- data_YR$ID


###Estimate lethal equivalents for AR using standardised FROH, i.e. the model used in the main inbreeding depression analysis

#Model including centered age as quadratic term, FROH:laflok and FROH:hatchyear_laflok random interaction, FROH scaled
formula_13_YR_FROH = recruits ~ z_FROH + gen_sex + c_age_comb + c_age_comb2 +
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(f_laflok2,z_FROH,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FROH, model = "iid",param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )) +
  f(ID2,model="iid",param = c(0.1,0.01),
    constr=TRUE)


#Run INLA
model.poi_YR_13_FROH = inla(formula=formula_13_YR_FROH, family="poisson",
                            data=data_YR,
                            control.compute=list(dic=T, config=T))
summary(model.poi_YR_13_FROH)


#Estimate lethal equivalents per island

###Extracting island specific inbreeding effects
### Sample from posterior distribution
nsamples <- 1000
model <- model.poi_YR_13_FROH
sample <- inla.posterior.sample(n=nsamples,model, add.names = TRUE, use.improved.mean=TRUE, intern=FALSE)
rm(model)

#Take the posterior values for FROH effect and each island effect separately and restore them in a data frame
island_samp <- as.data.frame(rep(NA,nsamples))
for (i in 1:nsamples){
  island_samp[i,1] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:1",1]
  island_samp[i,2] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:2",1]
  island_samp[i,3] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:3",1]
  island_samp[i,4] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:4",1]
  island_samp[i,5] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:5",1]
  island_samp[i,6] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:6",1]
  island_samp[i,7] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:7",1]
  island_samp[i,8] <- sample[[i]]$latent["z_FROH:1",1] + sample[[i]]$latent["f_laflok2:8",1]
}

rm(sample)

###Effect sizes
mean_effect_YR_FROH<-apply(island_samp,2,mean)
lc_YR_FROH<-apply(island_samp,2,quantile, 0.025)
uc_YR_FROH<-apply(island_samp,2,quantile, 0.975)

par(mai=c(1.02,1.15,0.82,0.42))
plot(c(1:8)~mean_effect_YR_FROH, xlim=c(-0.7,0.1), col=rainbow(8), lwd=3, yaxt="n", ylab="", xlab = expression(paste("Effect size of standardized ", italic("F"["ROH"]))), main="Effect of inbreeding on annual reproductive success")
abline(v=0, lty=3)
segments(uc_YR_FROH, 1:8, lc_YR_FROH, 1:8, col=rainbow(8), lwd=3)
axis(2, c(1:8), Island_names$Island, las=2)


#Estimate lethal equivalents from the standardised FROH results
#Estimate lethal equivalents and their upper and lower confidence limits using the mean of the posterior sample for each island
std_Froh_AR <- 0.03333521 #Estimated from sd(data_YR_temp$FROH), before transposing the data of individuals into long form
LE_AR <- as.data.frame(apply(island_samp,2,mean))
LE_AR[,2:7] <- NA
colnames(LE_AR) <- c("mean_effect_size", "lc_AR_FROH","uc_AR_FROH","scaled_LE", "lethal_equivalents","lc_AR_lethal_equivalents","uc_AR_lethal_equivalents")
rownames(LE_AR) <- c("Nesøy","Myken","Træna","Selvær","Gjerøy","Hestmannøy","Indre Kvarøy","Aldra")
LE_AR$scaled_LE <- -LE_AR$mean_effect_size*2
LE_AR$lethal_equivalents <- -LE_AR$mean_effect_size*2/std_Froh_AR
LE_AR$lc_AR_FROH<- apply(island_samp,2,quantile, 0.025)
LE_AR$lc_AR_lethal_equivalents <- (LE_AR$scaled_LE-((LE_AR$mean_effect_size-LE_AR$lc_AR_FROH)*2))/std_Froh_AR
LE_AR$uc_AR_FROH <- apply(island_samp,2,quantile, 0.975)
LE_AR$uc_AR_lethal_equivalents <- (((LE_AR$uc_AR_FROH-LE_AR$mean_effect_size)*2)+LE_AR$scaled_LE)/std_Froh_AR


#Write table of lethal equivalents in AR per island
write.table(LE_AR,"lethal_equivalents_AR_z_FROH.txt", col.names = T, row.names = T, quote = F)

rm(island_samp)

