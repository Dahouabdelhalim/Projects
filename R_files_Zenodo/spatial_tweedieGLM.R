#compare lightning frequency among ecosystem types
#Helene suggested including continent as a fixed effect so that we can estimate different intercepts and slopes for each continent

#open packages I need
library(dplyr)
library(geoR)
library(lme4)
library(Hmisc)


#clear workspace
rm(list=ls())

#import data 
agg_data#import the ENTLN data and save it as "agg_data"
str(agg_data)

######################################################
#create a new df that only includes rows representing more than 80% of their cells
#this makes sure that cells primarily represent individual ecosystem types
agg_data_ov80<-agg_data[agg_data$prop.cover.type>.8,]

##########################################
#create different spatial variables to account for spatial autocorrelation

#first create a variable at the 10 degree scale
breaks_10deg_long_cut<-seq(-180,180,10)
drop_cells_under10$longitude_10degbins<-cut(drop_cells_under10$Longitude,breaks_10deg_long_cut)
summary(as.factor(drop_cells_under10$longitude_10degbins))

#do the same with latitude
breaks_10deg_cut<-c(-23.5,seq(-16.5,23.5,10))
drop_cells_under10$latitude_10degbins<-cut(drop_cells_under10$Latitude,breaks_10deg_cut)
summary(as.factor(drop_cells_under10$latitude_10degbins))

#create bin aggregating all of these values in a grid
drop_cells_under10$bins_10deg<-paste(as.character(drop_cells_under10$longitude_10degbins),
                                     as.character(drop_cells_under10$latitude_10degbins),sep="_")
table(drop_cells_under10$bins_10deg)
length(levels(as.factor(drop_cells_under10$bins_10deg)))#107 10x10deg bins

#also create 0.5 degree bins
breaks_lat_halfdeg_cut<-seq(-23.5,23.5,0.5)
drop_cells_under10$latitude_halfdegbins<-cut(drop_cells_under10$Latitude,breaks_lat_halfdeg_cut)
breaks_long_halfdeg_cut<-seq(-180,180,0.5)
drop_cells_under10$longitude_halfdegbins<-cut(drop_cells_under10$Longitude,breaks_long_halfdeg_cut)
drop_cells_under10$bins_halfdeg<-paste(drop_cells_under10$longitude_halfdegbins,drop_cells_under10$latitude_halfdegbins,sep=".")
length(levels(as.factor(drop_cells_under10$bins_halfdeg)))#15594 0.5x0.5deg bins

#convert these bin groupings into factors
drop_cells_under10$bins_10deg<-as.factor(drop_cells_under10$bins_10deg)
drop_cells_under10$bins_halfdeg<-as.factor(drop_cells_under10$bins_halfdeg)

#filter out 0.5 degree bins with fewer than 10 values within each 
over10cells_per_halfdeg_bins<-drop_cells_under10%>%
  group_by(bins_halfdeg)%>%
  filter(n() >=10)

#how many cells were dropped in total?
#how many cells were dropped?
nrow(agg_data_ov80)-nrow(over10cells_per_halfdeg_bins)#5086 total cells removed
(nrow(agg_data_ov80)-nrow(over10cells_per_halfdeg_bins))/nrow(agg_data_ov80)*100#this is only 0.54% of all cells

#look at how data are distributed among ecosystem types and continents
table(interaction(over10cells_per_halfdeg_bins$Continent,over10cells_per_halfdeg_bins$aggregated_ecosystems))


##########################################################################
#confirm that there is variation within groups
View(over10cells_per_halfdeg_bins%>%
       group_by(bins_10deg,bins_halfdeg)%>%
       summarise(sd_CG_freq = sd(CG_FRD),
                 sd_CG_freq = sd(CG_FRD),
                 N_CG_freq = length(CG_count)))
####all but 1 cell has variation, so we should drop this cell to avoid issues with singularity

#drop the bin with no variation at all among its cells - "(-77.5,-77].(-12,-11.5]"
#this is Lima, which an all-Urban bin on the west coast of South America
over10cells_per_halfdeg_bins<-over10cells_per_halfdeg_bins[over10cells_per_halfdeg_bins$bins_halfdeg!="(-77.5,-77].(-12,-11.5]",]


#run a model using a Tweedie distribution (compound Poisson-gamma)
#the cpglm function works great for our purposes
#try analyzing these data with a Tweedie model
library(tweedie)
library(statmod)
library(cplm)

############################################
#run compound poisson-gamma models with 0.5 degree bins w/i 10 degree bins

#link = 0 means that we are using a log-link for the "link.power" function (this is from the "tweedie" function)
cplm_halfdeg_fit1<-cpglmm(CG_FRD~aggregated_ecosystems + (1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins)
summary(cplm_halfdeg_fit1)
cplm_halfdeg_fit2<-cpglmm(CG_FRD~ 1+(1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins)
anova(cplm_halfdeg_fit1,cplm_halfdeg_fit2)
AIC(cplm_halfdeg_fit1);AIC(cplm_halfdeg_fit2)

#save results
#creat relevant vector for aggregated ecosystems
agg_ecosystems<-levels(over10cells_per_halfdeg_bins$aggregated_ecosystems)
AllContinents_cpglmm_results<-data.frame(agg_ecosystems,"Continent" = rep("Pantropical",6),as.data.frame(summary(cplm_halfdeg_fit1)[5])[,1:2])

##########################################################
# quantile residual plot
u <- tweedie::ptweedie(cplm_halfdeg_fit1$y, cplm_halfdeg_fit1$p, fitted(cplm_halfdeg_fit1), cplm_halfdeg_fit1$phi)
u[cplm_halfdeg_fit1$y == 0] <- runif(sum(cplm_halfdeg_fit1$y == 0), 0, u[cplm_halfdeg_fit1$y == 0])
r2 <- qnorm(u)
qqnorm(r2, cex = 0.5)#the fit looks ok

#######################################################
#look at model residuals with variogram
subset_data<-data.frame(over10cells_per_halfdeg_bins[,c(6,5)],residuals = resid(cplm_halfdeg_fit1))
subset_variogram<-sample_n(subset_data,20000,replace=FALSE)
subset_geoR_dat<-as.geodata(subset_variogram,coords.col=c(1,2),data.col=3)

#create vector for separating values by degrees
deg_vec = seq(0,25,1)

#create variogram
subset_geoR_variog<-variog(subset_geoR_dat,max.dist = 20,uvec = deg_vec)
plot(subset_geoR_variog)

#now create a small-scale semivariogram

#create vector for separating values by degrees
deg_vec2 = seq(0,1,.05)

#create variogram
subset_geoR_variog<-variog(subset_geoR_dat,max.dist = 1,uvec = deg_vec2)
plot(subset_geoR_variog)


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#run the same model with just the Americas
over10cells_per_halfdeg_bins_Americas<-over10cells_per_halfdeg_bins[over10cells_per_halfdeg_bins$Continent=="Americas",]

cplm_halfdeg_fit1_Americas<-cpglmm(CG_FRD~aggregated_ecosystems + (1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Americas)
summary(cplm_halfdeg_fit1_Americas)
cplm_halfdeg_fit2_Americas<-cpglmm(CG_FRD~ 1+(1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Americas)
anova(cplm_halfdeg_fit1_Americas,cplm_halfdeg_fit2_Americas)
AIC(cplm_halfdeg_fit1_Americas);AIC(cplm_halfdeg_fit2_Americas)

#save results
Americas_cpglmm_results<-data.frame("agg_ecosystems"=agg_ecosystems,"Continent" = rep("Americas",6),as.data.frame(summary(cplm_halfdeg_fit1_Americas)[5])[,1:2])

##########################################################
# quantile residual plot
u <- tweedie::ptweedie(cplm_halfdeg_fit1_Americas$y, cplm_halfdeg_fit1_Americas$p, fitted(cplm_halfdeg_fit1_Americas), cplm_halfdeg_fit1_Americas$phi)
u[cplm_halfdeg_fit1_Americas$y == 0] <- runif(sum(cplm_halfdeg_fit1_Americas$y == 0), 0, u[cplm_halfdeg_fit1_Americas$y == 0])
r2 <- qnorm(u)
qqnorm(r2, cex = 0.5)#the fit looks ok, but not great
#######################################################
#look at model residuals with variogram
subset_data_Americas<-data.frame(over10cells_per_halfdeg_bins_Americas[,c(6,5)],residuals = resid(cplm_halfdeg_fit1_Americas))
subset_variogram_Americas<-sample_n(subset_data_Americas,20000,replace=FALSE)
subset_geoR_dat_Americas<-as.geodata(subset_variogram_Americas,coords.col=c(1,2),data.col=3)

#create vector for separating values by degrees
deg_vec = seq(0,25,1)

#create variogram
subset_geoR_variog_Americas<-variog(subset_geoR_dat_Americas,max.dist = 20,uvec = deg_vec)
plot(subset_geoR_variog_Americas)
####################################################
#now create a small-scale semivariogram

#create vector for separating values by degrees
deg_vec2 = seq(0,1,.05)

#create variogram
subset_geoR_variog<-variog(subset_geoR_dat,max.dist = 1,uvec = deg_vec2)
plot(subset_geoR_variog)


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#run the same model with just the Africa
over10cells_per_halfdeg_bins_Africa<-over10cells_per_halfdeg_bins[over10cells_per_halfdeg_bins$Continent=="Africa",]

cplm_halfdeg_fit1_Africa<-cpglmm(CG_FRD~aggregated_ecosystems + (1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Africa)
summary(cplm_halfdeg_fit1_Africa)
cplm_halfdeg_fit2_Africa<-cpglmm(CG_FRD~ 1+(1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Africa)
anova(cplm_halfdeg_fit1_Africa,cplm_halfdeg_fit2_Africa)
AIC(cplm_halfdeg_fit1_Africa);AIC(cplm_halfdeg_fit2_Africa)

#save results
Africa_cpglmm_results<-data.frame("agg_ecosystems"=agg_ecosystems,"Continent" = rep("Africa",6),as.data.frame(summary(cplm_halfdeg_fit1_Africa)[5])[,1:2])

##########################################################
# quantile residual plot
u <- tweedie::ptweedie(cplm_halfdeg_fit1_Africa$y, cplm_halfdeg_fit1_Africa$p, fitted(cplm_halfdeg_fit1_Africa), cplm_halfdeg_fit1_Africa$phi)
u[cplm_halfdeg_fit1_Africa$y == 0] <- runif(sum(cplm_halfdeg_fit1_Africa$y == 0), 0, u[cplm_halfdeg_fit1_Africa$y == 0])
r2 <- qnorm(u)
qqnorm(r2, cex = 0.5)#the fit looks quite good!
#######################################################
#look at model residuals with variogram
subset_data_Africa<-data.frame(over10cells_per_halfdeg_bins_Africa[,c(6,5)],residuals = resid(cplm_halfdeg_fit1_Africa))
subset_variogram_Africa<-sample_n(subset_data_Africa,20000,replace=FALSE)
subset_geoR_dat_Africa<-as.geodata(subset_variogram_Africa,coords.col=c(1,2),data.col=3)

#create vector for separating values by degrees
deg_vec = seq(0,25,1)

#create variogram
subset_geoR_variog_Africa<-variog(subset_geoR_dat_Africa,max.dist = 20,uvec = deg_vec)
plot(subset_geoR_variog_Africa)

####################################################
#now create a small-scale semivariogram

#create vector for separating values by degrees
deg_vec2 = seq(0,1,.05)

#create variogram
subset_geoR_variog<-variog(subset_geoR_dat,max.dist = 1,uvec = deg_vec2)
plot(subset_geoR_variog)


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#run the same model with just the Asia
over10cells_per_halfdeg_bins_Asia<-over10cells_per_halfdeg_bins[over10cells_per_halfdeg_bins$Continent=="Asia",]

cplm_halfdeg_fit1_Asia<-cpglmm(CG_FRD~aggregated_ecosystems + (1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Asia)
summary(cplm_halfdeg_fit1_Asia)
cplm_halfdeg_fit2_Asia<-cpglmm(CG_FRD~ 1+(1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Asia)
anova(cplm_halfdeg_fit1_Asia,cplm_halfdeg_fit2_Asia)
AIC(cplm_halfdeg_fit1_Asia);AIC(cplm_halfdeg_fit2_Asia)

#save results
Asia_cpglmm_results<-data.frame("agg_ecosystems"=agg_ecosystems,"Continent" = rep("Asia",6),as.data.frame(summary(cplm_halfdeg_fit1_Asia)[5])[,1:2])

##########################################################
# quantile residual plot
u <- tweedie::ptweedie(cplm_halfdeg_fit1_Asia$y, cplm_halfdeg_fit1_Asia$p, fitted(cplm_halfdeg_fit1_Asia), cplm_halfdeg_fit1_Asia$phi)
u[cplm_halfdeg_fit1_Asia$y == 0] <- runif(sum(cplm_halfdeg_fit1_Asia$y == 0), 0, u[cplm_halfdeg_fit1_Asia$y == 0])
r2 <- qnorm(u)
qqnorm(r2, cex = 0.5)#the fit is great

#######################################################
#look at model residuals with variogram
subset_data_Asia<-data.frame(over10cells_per_halfdeg_bins_Asia[,c(6,5)],residuals = resid(cplm_halfdeg_fit1_Asia))
subset_variogram_Asia<-sample_n(subset_data_Asia,20000,replace=FALSE)
subset_geoR_dat_Asia<-as.geodata(subset_variogram_Asia,coords.col=c(1,2),data.col=3)

#create vector for separating values by degrees
deg_vec = seq(0,25,1)

#create variogram
subset_geoR_variog_Asia<-variog(subset_geoR_dat_Asia,max.dist = 20,uvec = deg_vec)
plot(subset_geoR_variog_Asia)

####################################################
#now create a small-scale semivariogram

#create vector for separating values by degrees
deg_vec2 = seq(0,1,.05)

#create variogram
subset_geoR_variog<-variog(subset_geoR_dat,max.dist = 1,uvec = deg_vec2)
plot(subset_geoR_variog)

####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#run the same model with just the Australia
over10cells_per_halfdeg_bins_Australia<-over10cells_per_halfdeg_bins[over10cells_per_halfdeg_bins$Continent=="Australia",]

cplm_halfdeg_fit1_Australia<-cpglmm(CG_FRD~aggregated_ecosystems + (1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Australia)
summary(cplm_halfdeg_fit1_Australia)
cplm_halfdeg_fit2_Australia<-cpglmm(CG_FRD~ 1+(1|bins_10deg) + (1|bins_halfdeg), link = 0,data = over10cells_per_halfdeg_bins_Australia)
anova(cplm_halfdeg_fit1_Australia,cplm_halfdeg_fit2_Australia)
AIC(cplm_halfdeg_fit1_Australia);AIC(cplm_halfdeg_fit2_Australia)

#save results
Australia_cpglmm_results<-data.frame("agg_ecosystems"=agg_ecosystems[1:5],"Continent" = rep("Australia",5),as.data.frame(summary(cplm_halfdeg_fit1_Australia)[5])[,1:2])

##########################################################
# quantile residual plot
u <- tweedie::ptweedie(cplm_halfdeg_fit1_Australia$y, cplm_halfdeg_fit1_Australia$p, fitted(cplm_halfdeg_fit1_Australia), cplm_halfdeg_fit1_Australia$phi)
u[cplm_halfdeg_fit1_Australia$y == 0] <- runif(sum(cplm_halfdeg_fit1_Australia$y == 0), 0, u[cplm_halfdeg_fit1_Australia$y == 0])
r2 <- qnorm(u)
qqnorm(r2, cex = 0.5)#the fit is pretty great

#######################################################
#look at model residuals with variogram
subset_data_Australia<-data.frame(over10cells_per_halfdeg_bins_Australia[,c(6,5)],residuals = resid(cplm_halfdeg_fit1_Australia))
subset_variogram_Australia<-sample_n(subset_data_Australia,20000,replace=FALSE)
subset_geoR_dat_Australia<-as.geodata(subset_variogram_Australia,coords.col=c(1,2),data.col=3)

#create vector for separating values by degrees
deg_vec = seq(0,25,1)

#create variogram
subset_geoR_variog_Australia<-variog(subset_geoR_dat_Australia,max.dist = 20,uvec = deg_vec)
plot(subset_geoR_variog_Australia)
####################################################
#now create a small-scale semivariogram

#create vector for separating values by degrees
deg_vec2 = seq(0,1,.05)

#create variogram
subset_geoR_variog<-variog(subset_geoR_dat,max.dist = 1,uvec = deg_vec2)
plot(subset_geoR_variog)

####################################################

#save model results
cpglmm_model_results_halfdegree<-rbind.data.frame(AllContinents_cpglmm_results,Americas_cpglmm_results,Africa_cpglmm_results,Asia_cpglmm_results,Australia_cpglmm_results)
write.csv(cpglmm_model_results_halfdegree,"cpglmm_model_results_halfdegree.csv");getwd()
####################################################

