#this script uses MCMCglmm and QGglmm to run a model and extract the quant gen parameters 

require(MCMCglmm)
require(QGglmm)
require(ggplot2)
require(tidyverse)

#load the data

map_pop<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/IndivPlast_MappingPops/TomalesBay_capsuledata.csv", header=TRUE, sep=",")
SvN<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/SoCALvsNorCAL/SoCalvNorCal_Data.csv", header=TRUE, sep=",")

#get it all in one dataframe

p_map_pop<-subset(map_pop, gen %in% "p")

p_map_pop$date<-as.Date(p_map_pop$date, "%m/%d/%Y")

SvN$date<-as.Date(SvN$date, "%m/%d/%Y")

#rbind SvN and maping pops to make a master high versus low (hvl) data set

HvL<-rbind(SvN[,c("Capsule", "original", "date", "Family", "ID")],p_map_pop[, c("Capsule", "original", "date", "Family", "ID")])

HvL$date<-as.Date(HvL$date, "%m/%d/%Y")

HvL$type_150 <- ifelse(HvL$Capsule>0.150,1,0)


#get data for capsule size as individual means

ID_mean<-aggregate(HvL$Capsule, list(HvL$ID) , mean)
colnames(ID_mean)<-c("ID", "mean_cap")

#now merge them together with the original data so we can get the Family etc info
ID_mean_tog<-merge(ID_mean, HvL, by=c("ID"))
ID_mean_tog<-subset(ID_mean_tog, select = -c(Capsule, date))
ID_mean_tog_uni<-unique.data.frame(ID_mean_tog)


#get type, but make it a binary variable where 1 = lecithotrophy and 0 = planktotrophy

ID_mean_tog_uni$type_150 <- ifelse(ID_mean_tog_uni$mean_cap>0.150,1,0)

#now we have our data

#look at how the responses are distributed
ggplot(data = HvL) + geom_density(aes(x = Capsule), size = 2)


ggplot(data = HvL) + geom_histogram(aes(x = type_150))

#models run with MCMCglmm and QGglmm


#make a dataframe for my data with the ID renamed as animal
colnames(HvL)<-c("Capsule", "original", "date", "Family", "animal", "eggmass")

#get region for each family
HvL$region<-ifelse(grepl("SC*", HvL$Family), "Long Beach", ifelse(grepl("W128|13[0-4]", HvL$Family), "Tomales", "Mill Valley"))

#get region for family for the mean per egg mass dataset
ID_mean_tog_uni$region<-ifelse(grepl("SC*", ID_mean_tog_uni$Family), "Long Beach", ifelse(grepl("W128|13[0-4]", ID_mean_tog_uni$Family), "Tomales", "Mill Valley"))


#make a pedigree with the animal (ID) and the Family info
pedigree<-HvL[!duplicated(HvL$animal),c("animal", "Family")]
colnames(pedigree)<-c("animal", "dam")
pedigree$sire<-NA

#add in the parental families and rbind them to the data

parents<-as.data.frame(cbind(as.vector(unique(HvL$Family)), NA, NA))
colnames(parents)<-c("animal", "dam", "sire")

pedigree_all<-rbind(parents, pedigree)

#the data for capsule size is not normally distributed, but the residules are (mostly) when we account for salinity (original)

#there is no latent scale here. This is a gaussian trait. 

#make a prior
prior <- list(R = list(V = 1, nu = 1),
              G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                       G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                       G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                       G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))



#make a new prior for looking at type_150

prior_bin <- list(R = list(V = 1, fix = 1),
                  G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                  G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                  G3 = list (V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                  G4 = list (V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

#to run a nested model in MCMCglmm you need to made the nested factors uniquely identifiable 

HvL$Familyregion<-paste0(HvL$Family, HvL$region)


# separate the animal column by letters A-D to find the dish name and add dish to the models
#first change the SCW label to something without a C
ID<-sub("SCW", "SXW", HvL$animal)

dish<-strsplit(ID, c("A", "B", "C", "D"), 1)

HvL$dish_f<-sapply(dish, "[", 1)


model_3_animal <- MCMCglmm(Capsule ~ 1 + original, random = ~animal + Family + dish_f + Familyregion, 
                    data = HvL, prior=prior, pedigree=pedigree_all, family = "gaussian", 
                    nit = 50000, burnin = 10000, thin = 10) #sample size greater than 1000!

model3.3<- MCMCglmm(eggmass ~ 1 + original, random = ~animal + Family + dish_f + Familyregion, 
                     data = HvL, prior=prior_bin, pedigree = pedigree_all, family = "threshold", 
                     nit = 603000, burnin = 10000, thin = 10)#sample size is greater than 1000!

#run for categorical trait
parental<-sel_all_mlp_tog[sel_all_mlp_tog$gen == "p",]
colnames(parental)<-c("date", "animal", "eggmass", "original", "x", "is_l", 
                      "is_p", "is_m", "eggmass_type", "Family", "current", "gen", "region")


parental$Familyregion<-paste0(parental$Family, parental$region)


model_mlp<-MCMCglmm(eggmass_type ~ 1+original,  random = ~animal + Family + dish_f + Familyregion,
                    data = parental, pedigree= pedigree_all, prior=prior_bin, family = "ordinal", 
                    nit = 603000, burnin = 10000, thin = 10)

save(model_3_animal, file = "model3_capsule_AmNat.Rdata")
load("model3_binary.Rdata")
summary(model3.3)

#plots of the postior destributions

library(reshape2)
ref <- data.frame(model3.3$VCV)
rdf <- melt(ref, value.name="original")
ggplot(rdf, aes(x=sqrt(original),color=variable))+geom_density()

plot(log(model_3_animal$VCV))

#pull out genetic variance components

mu_gaus <- model_3_animal[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_gaus <- model_3_animal[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_gaus <-model_3_animal[["VCV"]][ , "units"]
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_gaus<-va_gaus/(va_gaus+vr_gaus)

posterior.mode(h2_gaus)
median(h2_gaus)		
HPDinterval(h2_gaus)	#0.5044 0.5796

#a few parameters
median(mu_gaus)
HPDinterval(mu_gaus)

median(va_gaus)
HPDinterval(va_gaus)

median(vr_gaus)
HPDinterval(vr_gaus)



#for the binary model
mu_bin_lat <- model3.3[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_bin_lat <- model3.3[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_bin_lat <-(model3.3[["VCV"]][ , "animal"] + model3.3[["VCV"]][ , "units"] - 1)
var.p <- (model3.3[["VCV"]][ , "animal"] + model3.3[["VCV"]][ , "units"] - 1)
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_bin_lat <- va_bin_lat / (va_bin_lat + var.p)

posterior.mode(h2_bin_lat)
median(h2_bin_lat)		
HPDinterval(h2_bin_lat)	#0.445 0.4688

## Computing the parameters on the liability of a threshold model
#We just need to add the "link variance" of 1 to the above calculation
h2_bin_liab <- va_bin_lat / (va_bin_lat + vr_bin_lat + 1)
posterior.mode(h2_bin_liab)	#Posterior Mode: 0.456
HPDinterval(h2_bin_liab)	#Credible Interval: [0.445 - 0.468]

## Computing the parameters on the data scale
#Using QGglmm
out_bin <- QGparams(mu = posterior.mode(mu_bin_lat),
                var.p = posterior.mode(va_bin_lat + vr_bin_lat),
                var.a = posterior.mode(va_bin_lat),
                model = "binom1.probit")
h2_bin <- out[["h2.obs"]]	#0.2281



#model stuff for the mlp model

mu_cat_lat <- model_mlp[["Sol"]][,"(Intercept)"]
#Cut-point
cp_cat <- model_mlp[["CP"]][ , 1]
#Additive genetic variance on the latent scale
va_cat_lat <- model_mlp[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_cat_lat <- model_mlp[["VCV"]][ , "units"]
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_cat_lat <- va_cat_lat / (va_cat_lat + vr_cat_lat)
posterior.mode(h2_cat_lat)    #Posterior Mode: 0.93
median(h2_cat_lat)            #Posterior Median: 0.93
HPDinterval(h2_cat_lat)       #Credible Interval: [0.91 - 0.95]

## Computing the parameters on a threshold model liability scale
#Since we used MCMCglmm's "ordinal" family, we can relate it to a 
#  multiple threshold model by adding the "link variance" (i.e. 1)
h2_perso_liab <- va_perso_lat / (va_perso_lat + vr_perso_lat + 1)
posterior.mode(h2_perso_liab)   #0.23

## Computing the parameters on the data scale
#Using QGglmm
out_mlp <- QGparams(mu = posterior.mode(mu_perso_lat),
                var.p = posterior.mode(va_perso_lat + vr_perso_lat),
                var.a = posterior.mode(va_perso_lat),
                model = "ordinal", 
                cut.points = c(-Inf, 0, posterior.mode(cp_perso), Inf))
h2_perso <- out[["h2.obs"]] #0.12 / 9.66e-05 / 0.12





####re-run models with pedigrees that contain full sibs or half sibs


#make self pedigree
pedigree_self<-pedigree_all[,c(1,2)]
pedigree_self$sire<-pedigree_all$dam
head(pedigree_self)

model_3s <- MCMCglmm(Capsule ~ 1+original+region, random = ~Family +animal + Familyregion, 
                    data = HvL, prior=prior, pedigree=pedigree_self, family = "gaussian", 
                    nit = 500000, burnin = 10000, thin = 100) #sample size all greater than 1000!

model3.3s <- MCMCglmm(eggmass ~ 1 + original + region, random = ~Family +animal + Familyregion, 
                     data = HvL, prior=prior_bin, pedigree = pedigree_self, family = "threshold", 
                     nit = 500000, burnin = 10000, thin = 10)#sample size is greater than 1000!

mu_s_gaus <- model_3s[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_s_gaus <- model_3s[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_s_gaus <-model_3s[["VCV"]][ , "units"]
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_s_gaus<-va_s_gaus/(va_s_gaus+vr_s_gaus)

posterior.mode(h2_s_gaus)
median(h2_s_gaus)		
HPDinterval(h2_s_gaus)	#0.5044 0.5796

#a few parameters
median(mu_s_gaus)#0.129
HPDinterval(mu_s_gaus)

median(va_s_gaus)
HPDinterval(va_s_gaus)

median(vr_s_gaus)
HPDinterval(vr_s_gaus)

#for the threshold model
#for the binary model
mu_s_bin_lat <- model3.3s[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_s_bin_lat <- model3.3s[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_s_bin_lat <-(model3.3s[["VCV"]][ , "animal"] + model3.3s[["VCV"]][ , "units"] - 1)
var.p_s <- (model3.3s[["VCV"]][ , "animal"] + model3.3s[["VCV"]][ , "units"] - 1)
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_bin_lat_s <- va_s_bin_lat / (va_s_bin_lat + var.p_s)

posterior.mode(h2_bin_lat_s)
median(h2_bin_lat_s)		
HPDinterval(h2_bin_lat_s)	#0.445 0.4688

## Computing the parameters on the liability of a threshold model
#We just need to add the "link variance" of 1 to the above calculation
h2_bin_liab_s <- va_s_bin_lat / (va_s_bin_lat + vr_s_bin_lat + 1)
posterior.mode(h2_bin_liab_s)	#Posterior Mode: 0.456
HPDinterval(h2_bin_liab_s)	#Credible Interval: [0.445 - 0.468]

## Computing the parameters on the data scale
#Using QGglmm
out_bin_s <- QGparams(mu = posterior.mode(mu_s_bin_lat),
                    var.p = posterior.mode(va_s_bin_lat + vr_s_bin_lat),
                    var.a = posterior.mode(va_s_bin_lat),
                    model = "binom1.probit")
h2_bin_s <- out[["h2.obs"]]	#0.2281



#make full-sib pedigree. Generate some new sires for each maternal line

pedigree_full<-pedigree_self[,c(1,2)]

#get unique list of moms
moms<-unique(pedigree_full$dam)

#make a vector for sires
pedigree_full$sire[pedigree_self$dam == moms[3]]<-"S3"

#run a for loop to populate each dam with a unique sire

for(i in 1:17){
  pedigree_full$sire[pedigree_self$dam == moms[i]]<-paste("s", i, sep="")
}


#add in the new sires
parents_full<-as.data.frame(cbind(as.vector(unique(pedigree_full$sire)), NA, NA))
parents_full
colnames(parents_full)<-c("animal", "dam", "sire")

pedigree_fall<-rbind(parents_full[c(-1),], pedigree_full)

#then rerun the model

model_3full <- MCMCglmm(Capsule ~ 1+original, random = ~Family +animal + Familyregion, 
                     data = HvL, prior=prior, pedigree=pedigree_fall, family = "gaussian", 
                     nit = 5000, burnin = 3000, thin = 10) #this looks better. Sample size is still low

model3.3f <- MCMCglmm(eggmass ~ 1 + original + region, random = ~Family +animal + Familyregion,
                      data = HvL, prior=prior_bin, pedigree = pedigree_fall, family = "threshold", 
                      nit = 5000, burnin = 3000, thin = 10)#sample size is greater than 1000!


#for the gausian model
mu_f_gaus <- model_3full[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_f_gaus <- model_3full[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_f_gaus <-model_3full[["VCV"]][ , "units"]
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_f_gaus<-va_f_gaus/(va_f_gaus + vr_f_gaus)

posterior.mode(h2_f_gaus)
median(h2_f_gaus)		
HPDinterval(h2_f_gaus)	#0.5044 0.5796

#a few parameters
median(mu_f_gaus)
HPDinterval(mu_f_gaus)

median(va_f_gaus)
HPDinterval(va_f_gaus)

median(vr_f_gaus)
HPDinterval(vr_f_gaus)


##For the threshold model
mu_bin_lat_f <- model3.3f[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_bin_lat_f <- model3.3f[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_bin_lat_f <-(model3.3f[["VCV"]][ , "animal"] + model3.3f[["VCV"]][ , "units"] - 1)
var.p_f <- (model3.3f[["VCV"]][ , "animal"] + model3.3f[["VCV"]][ , "units"] - 1)
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_bin_lat_f <- va_bin_lat_f / (va_bin_lat_f + var.p_f)

posterior.mode(h2_bin_lat_f)
median(h2_bin_lat_f)		
HPDinterval(h2_bin_lat_f)	#0.445 0.4688

## Computing the parameters on the liability of a threshold model
#We just need to add the "link variance" of 1 to the above calculation
h2_bin_liab_f <- va_bin_lat_f / (va_bin_lat_f + vr_bin_lat_f + 1)
posterior.mode(h2_bin_liab_f)	#Posterior Mode: 0.456
HPDinterval(h2_bin_liab_f)	#Credible Interval: [0.445 - 0.468]

## Computing the parameters on the data scale
#Using QGglmm
out_f <- QGparams(mu = posterior.mode(mu_bin_lat_f),
                var.p = posterior.mode(va_bin_lat_f + vr_bin_lat_f),
                var.a = posterior.mode(va_bin_lat_f),
                model = "binom1.probit")
h2_bin_f <- out_f[["h2.obs"]]	#0.2281
h2_bin_f

#make a pedigree that is half full-sib and half-half sib. 

#sample 
#use pedigree_fall and pedigree_s
#make a mixed pedigree
pedigree_mix<-pedigree_full[order(pedigree_full$dam),c(1,2)]
head(pedigree_mix)

#remove the dams for now

pedigree_mix2<-pedigree_mix[1:433,]

#this makes a list of sires combined from full and self
test<-c(rbind(pedigree_full$sire, as.character(pedigree_self$sire)))
test<-test[!is.na(test)]

#make for the dams. We will use this to loop through
mat<-c(rbind(as.character(pedigree_full$dam), as.character(pedigree_full$dam)))
mat<-mat[!is.na(mat)]

#goal is to sample the list of sires by the number of individuals in each maternal family

pedigree_mix$sire
#this is a stupid way to do this, but the indexing is complicated
one<-sample(test[1:42], table(pedigree_fall$sire)[1])
two<-sample(test[43:81], table(pedigree_fall$sire)[2])
three<-sample(test[82:122], table(pedigree_fall$sire)[3])
four<-sample(test[123:154], table(pedigree_fall$sire)[4])
five<-sample(test[155:180], table(pedigree_fall$sire)[5])
six<-sample(test[182:214], table(pedigree_fall$sire)[6], replace=T)
seven<-sample(test[215:264], table(pedigree_fall$sire)[7])
eight<-sample(test[265:308], table(pedigree_fall$sire)[8])
nine<-sample(test[309:346], table(pedigree_fall$sire)[9])
ten<-sample(test[347:386], table(pedigree_fall$sire)[10], replace=T)
eleven<-sample(c("W128", "s12"), table(pedigree_fall$sire)[11], replace=T)
twelve<-sample(c("W129", "s13"), table(pedigree_fall$sire)[12], replace=T)
thirteen<-sample(c("W130", "s14"), table(pedigree_fall$sire)[13], replace=T)
fourteen<-sample(c("W131", "s15"), table(pedigree_fall$sire)[14], replace=T)
fifteen<-sample(c("W132", "s16"), table(pedigree_fall$sire)[15], replace=T)
sixteen<-sample(c("W133", "s17"), table(pedigree_fall$sire)[16]+5, replace=T)

#make a string of na's for the top
n<-rep(NA, 33)


#put it all together

pedigree_mix2$sire<-as.factor(c(one, two, three, four, five, six, seven, eight, nine,ten, eleven,sixteen, twelve, thirteen, fourteen, fifteen))

#add in the dam and sire ids
pedigree_mix<-rbind(pedigree_fall[1:33,], pedigree_mix2)
#need to make the pedigree in terms of the dams. Need to have an easy way to specify the ratio of full versus self or half



#run the models

model_3mix <- MCMCglmm(Capsule ~ 1+original, random = ~Family +animal + Familyregion, 
                        data = HvL, prior=prior, pedigree=pedigree_mix, family = "gaussian", 
                        nit = 5000, burnin = 3000, thin = 10) #this looks better. Sample size is still low

model3.3mix <- MCMCglmm(eggmass ~ 1 + original + region, random = ~Family +animal + Familyregion,
                      data = HvL, prior=prior_bin, pedigree = pedigree_mix, family = "threshold", 
                      nit = 5000, burnin = 3000, thin = 10)#sample size is greater than 1000!

#for the gausian model
mu_m_gaus <- model_3mix[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_m_gaus <- model_3mix[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_m_gaus <-model_3mix[["VCV"]][ , "units"]
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_m_gaus<-va_m_gaus/(va_m_gaus + vr_m_gaus)

posterior.mode(h2_m_gaus)
median(h2_m_gaus)		
HPDinterval(h2_m_gaus)	#0.5044 0.5796

#a few parameters
median(mu_m_gaus)
HPDinterval(mu_m_gaus)

median(va_m_gaus)
HPDinterval(va_m_gaus)

median(vr_m_gaus)
HPDinterval(vr_m_gaus)


##For the threshold model
mu_bin_lat_m <- model3.3mix[["Sol"]][,"(Intercept)"]
#Additive genetic variance on the latent scale
va_bin_lat_m <- model3.3mix[["VCV"]][ , "animal"]
#Residual variance on the latent scale
vr_bin_lat_m <-(model3.3mix[["VCV"]][ , "animal"] + model3.3mix[["VCV"]][ , "units"] - 1)
var.p_m <- (model3.3mix[["VCV"]][ , "animal"] + model3.3mix[["VCV"]][ , "units"] - 1)
#Heritability on the latent scale (h2 = Va/ (Va + Vr))
h2_bin_lat_m <- va_bin_lat_m / (va_bin_lat_m + var.p_m)

posterior.mode(h2_bin_lat_m)
median(h2_bin_lat_m)		
HPDinterval(h2_bin_lat_m)	#0.445 0.4688

## Computing the parameters on the liability of a threshold model
#We just need to add the "link variance" of 1 to the above calculation
h2_bin_liab_m <- va_bin_lat_m / (va_bin_lat_m + vr_bin_lat_m + 1)
posterior.mode(h2_bin_liab_m)	#Posterior Mode: 0.456
HPDinterval(h2_bin_liab_m)	#Credible Interval: [0.445 - 0.468]

## Computing the parameters on the data scale
#Using QGglmm
out_m <- QGparams(mu = posterior.mode(mu_bin_lat_m),
                  var.p = posterior.mode(va_bin_lat_m + vr_bin_lat_m),
                  var.a = posterior.mode(va_bin_lat_m),
                  model = "binom1.probit")
h2_bin_m <- out_m[["h2.obs"]]	#0.2281
h2_bin_m
