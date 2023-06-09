#Variance partitioning for inbreeding in Helgeland house sparrows
#######################################################################
#Alina Niskanen
#alina.niskanen@gmail.com
#April 2020

library(INLA)
library(MCMCglmm)


#Import the phenotypic data and inbreeding estimates
data_LRS <- read.csv("Data.txt", header = T, stringsAsFactors = F, sep = "\\t")

#Import population size estimate per island for all study years
Pop_size_w_SNPs <- read.csv("Pop_size_1998_2013.csv", header = T, stringsAsFactors = F, sep = "\\t") 

###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_LRS[which(data_LRS$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)

#Add presence of an island in 2013, for which no offspring data is available based on SNPs
data_adult_temp$X2013 <- NA
data_adult_temp$X2013[which(data_adult_temp$minadyear<=2013 & data_adult_temp$maxadyear >=2013)] <- 1

#Start choosing the correct individuals that were observed as adults (laflok) on one of the 8 study islands
data_absence <- data_adult_temp[(is.na(data_adult_temp$X1998) & is.na(data_adult_temp$X1999) & is.na(data_adult_temp$X2000) & is.na(data_adult_temp$X2001) & is.na(data_adult_temp$X2002) & is.na(data_adult_temp$X2003) & is.na(data_adult_temp$X2004) & is.na(data_adult_temp$X2005) & is.na(data_adult_temp$X2006) & is.na(data_adult_temp$X2007) & is.na(data_adult_temp$X2008) & is.na(data_adult_temp$X2009) & is.na(data_adult_temp$X2010) & is.na(data_adult_temp$X2011) & is.na(data_adult_temp$X2012) & is.na(data_adult_temp$X2013)),]

#Make the final dataset, includes only individuals with laflok on one of the 8 study islands 
data_inbreeding <- data_adult_temp[!data_adult_temp$id %in% data_absence$id,]

#Remove temp files
rm(list=ls(pattern="temp"))

# make island factor
data_inbreeding$f_laflok <- as.factor(data_inbreeding$laflok)

#Add combined variable for hatch year and adult island
data_inbreeding$IAY <- as.factor(paste(data_inbreeding$all_hatchyears,data_inbreeding$laflok,sep="_"))
#Add combined variable for hatch year and adult island
data_inbreeding$IAY2 <- data_inbreeding$IAY


############################################################
### Animal model to account for relatedness among individuals
##############################################################

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

# also need to add ID column to data files:
data_inbreeding$ID <- d.map[match(data_inbreeding$id, d.map$id), "ID"]
data_inbreeding$ID2 <- data_inbreeding$ID

#Scale FGRM
data_inbreeding$z_FGRM <- scale(data_inbreeding$FGRM, center = T, scale = T)


#Inbreeding variance including sex as a fixed effect, and adult island and hatchyear nested within island as random effects, FGRM standardised
formula_0.1_inbr = z_FGRM ~ gen_sex +  
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01))


#Inbreeding variance including sex and island status as fixed effects, and adult island and hatchyear nested within island as random effects, FGRM standardised
formula_3.1_inbr = z_FGRM ~ gen_sex + laflok_status + 
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01))


#Run INLA

#Inbreeding variance including sex as a fixed effect, and adult island and hatchyear nested within island as random effects, FGRM standardised
model.inbr_0.1 = inla(formula=formula_0.1_inbr, family="gaussian",
                    data=data_inbreeding,
                    control.compute=list(dic=T, config=T)
)
summary(model.inbr_0.1)
inlaPosteriors_inbr_0(model.inbr_0.1)


#Inbreeding variance including sex and island status as fixed effects, and adult island and hatchyear nested within island as random effects, FGRM standardised
model.inbr_3.1 = inla(formula=formula_3.1_inbr, family="gaussian",
                    data=data_inbreeding,
                    control.compute=list(dic=T, config=T)
)
summary(model.inbr_3.1)
inlaPosteriors_inbr_0(model.inbr_3.1)

#Plot
plot(model.inbr_0.1$marginals.hyperpar$'Precision for f_laflok')
plot(model.inbr_0.1$marginals.hyperpar$'Precision for IAY')


#Formula to transform precision to variance
inlaPosteriors_inbr_0 <- function(model){
  sigma.island = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for f_laflok")
  m.island=inla.mmarginal(sigma.island)
  e.island=inla.emarginal(function(x) x, sigma.island)
  
  sigma.IAY = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for IAY")
  m.IAY=inla.mmarginal(sigma.IAY)
  e.IAY=inla.emarginal(function(x) x, sigma.IAY)
  
  sigma.e = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for the Gaussian observations")
  m.e =inla.mmarginal(sigma.e)
  e.e =inla.emarginal(function(x) x, sigma.e)
  
  results_tab <- cbind(rbind(VarE = e.e, island = e.island, Island_year=e.IAY),
                       rbind(VarE = m.e, island = m.island, Island_year=m.IAY),
                       1/model$summary.hyperpar[,c(5,3)])
  
  names(results_tab) <- c("Mean","Mode", "2.5%", "97.5 %")
  return(results_tab)
}


#Make output file
sink("Inbreeding_variance.txt")
cat("###Model summaries inbreeding variance partitioning.###", sep = "\\n")
cat("\\n")

cat("####Inbreeding variance for sex, adult island and island_year, FGRM standardised.###")
cat("\\n")
summary(model.inbr_0.1)
cat("\\n")
inlaPosteriors_inbr_0(model.inbr_0.1)
cat("\\n")
cat("\\n")
cat("####Inbreeding variance for sex, habitat type, adult island and island_year, FGRM standardised.###")
cat("\\n")
summary(model.inbr_3.1)
cat("\\n")
inlaPosteriors_inbr_0(model.inbr_3.1)
cat("\\n")

sink()

