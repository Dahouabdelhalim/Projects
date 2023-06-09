#############################################################
#PETSCH, DK; BERTONCIN, A.P.S; ORTEGA, J.C.G; THOMAZ, S.M.
#Exotic species drive biotic homogenization, but it# 
#depends on the realm, dimension and study design. Oikos#####
#############################################################
#R script last checked in October, 20th, 2021, by Authors
#Import 'metafor' package and data:
library(metafor)
meta_homog_data <- read.csv("meta_homog.csv", h = T)
#"meta_homog.csv" is a '.csv' dataplan with "Studies Data" tab from Supplementary Material

#Compute effect size estimates (log response ratio):
ef_si_ori <- escalc(measure="ROM", #"ROM" is the log response ratio code in 'escalc' function
                    m1i= Mean_Invaded, sd1i= SD_Invaded, n1i= n_Invaded,#mean, sd and sample size for treatment (samples with non-native species) 
                    m2i= Mean_Control, sd2i= SD_Control, n2i= n_Control,#same for control (samples without non-native species)
                    data = meta_homog_data)
names(ef_si_ori) #Checking column names
dim(ef_si_ori) #Checking matrix dimension
head(ef_si_ori) #Checking first lines

#Some studies presented similarities instead of dissimilarities. We changed the sign of these effect sizes:
yi <- ifelse (ef_si_ori$Similarity == "Similarity", ef_si_ori$yi * -1, ef_si_ori$yi)#convertendo similaridade em dissimilaridade

#Then we replaced the column with effect sizes (yi) with the effect sizes corrected for their sign:
ef_si <- cbind(ef_si_ori[ ,-20], yi)

#95 confidence interval (95 CI) estimates for each effect size:
low_b <- ef_si$yi - (1.96*sqrt(ef_si$vi)) #lower bound
upp_b <- ef_si$yi + (1.96*sqrt(ef_si$vi)) #upper bound

#Computing how many positive, negative and non-significant effects:
sum(low_b > 0)
sum(upp_b < 0)
sum(low_b < 0 & upp_b > 0)

########################
#Meta-analytical models#
########################

#Code the study-level ("Study_number") and between-study heterogeneity ("Bet_stu").
#These objects are needed for "rma.mv" below to correctly estimate the multi-level estimates:
Study_number <- ef_si$Study_number #study-level heterogeneity
Bet_stu <- seq(from = 1, to = nrow(ef_si)) #between-study heterogeneity

#Mean weighted effect size (a model with only the intercept):
res.redu <- rma.mv(data= ef_si,  mods= ~ 1,
                   yi= yi, V= vi, 
                   random=list( ~1|Bet_stu, ~1|Study_number), 
                   method="ML")
res.redu
#95 CI for the weighthed effect size:
(res.redu$ci.ub - res.redu$ci.lb)/2

#Meta-regression model assessing the effect of only spatial extent:
plot(ef_si$yi ~ ef_si$Spatialextent_km)
plot(ef_si$yi ~ log(ef_si$Spatialextent_km)) #Transforming by log improve the relationship reducing the difference in scale of spatial extent

res.spa <- rma.mv(data= ef_si,  mods= ~ log(Spatialextent_km),
                  yi= yi, V= vi, 
                  random=list( ~1|Bet_stu, ~1|Study_number), 
                  method="ML")
res.spa

#We ran a model with only the intercept but excluding the "Native_range" level from subgroup analysis due to the low number of studies.
#This intercept-only model is necessary to compute Pseudo-R^2 below:
res.redu.1 <- rma.mv(data= ef_si,  mods= ~ 1,
                     yi= yi, V= vi,
                     subset= ef_si$Type_of_control!="Native_range",
                     random=list( ~1|Bet_stu, ~1|Study_number), 
                     method="ML") #just intercept

#Full subgroup analysis with dimension (taxonomic, functional or phylogenetic),
#type of control used in the primary studies, and realm as moderator variables:
res.full <- rma.mv(data= ef_si,  
                   mods= ~ Dimension  + Type_of_control + Realm,
                   yi=yi,V=vi,
                   subset= ef_si$Type_of_control!="Native_range",
                   random = list( ~1|Bet_stu, ~1|Study_number), 
                   method="ML")
res.full

#Pseudo-R^2:
round(1-(res.full$sigma2/res.redu.1$sigma2),3)
#The first value correspond to between-study heterogeneity and the second to the study-level heterogeneity

#I^2 statistic for "rma.mv" function:
s2_fun<-function(vi, data){
  #"vi" is a vector with effect size variances
  #"data" is the object with effect size variances
  wi<-1/vi
  wi.2<-wi^2
  k<-nrow(data)
  first.term<-(k-1)*sum(wi)
  second.term<-(sum(wi)^2)-(sum(wi.2))
  s2_obs<-first.term/second.term
  return(s2_obs)
}

error_meta <- s2_fun(ef_si$vi, data = ef_si) #Sampling error estimate
vartotal <- error_meta + sum(res.redu$sigma2) #Total variance (heterogeneity)
i2 <- res.redu$sigma2/vartotal #I^2 statistic
i2 *100 #The first value correspond to between-study heterogeneity and the second to the study-level heterogeneity

############################
#End meta-analitycal models#
############################

##################
#Publication bias#
##################
#We computed mean effect size and variance by study to avoid multiplicity:
mean_yi <- as.numeric(unlist(by(ef_si$yi, INDICES = ef_si$Study_number, FUN = mean)))
mean_vi <- as.numeric(unlist(by(ef_si$vi, INDICES = ef_si$Study_number, FUN = mean)))

#Orwin fail-safe number:
mean(mean_yi) #Reference value for publication bias assessment

fsn(mean_yi, mean_vi, type ="Orwin", target = -0.1401*.25) #Number of studies required to reduce to 25% of the observed mean effect
fsn(mean_yi, mean_vi, type ="Orwin", target = -0.1401*.5) #Number of studies required to reduce to 50% of the observed mean effect
fsn(mean_yi, mean_vi, type ="Orwin", target = -0.1401*.75) #Number of studies required to reduce to 75% of the observed mean effect

#Trim-and-fill procedure:
funnel_resu <- rma.uni(mean_yi, mean_vi,method="ML")

trim <- trimfill(funnel_resu)
trim

#Just to see where are the missing studies:
funnel(trim)
####################End