# Predicting activity patterns of Gila monster using daily microclimatic data and niche modeling
# Toro-Cardona et al. 2022 
# Journal of Animal Ecology
# The following code includes all process made in our study.
#
#############################################################

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(NicheMapR)


setwd("D:/Proyectos/En proceso/MAESTRIA/Manuscrito/JOAE_review_coments/Analisis/Modelos_0m/")

# Get nc data
get.global.climate(folder="Microclimate_input")

nc_data <- nc_open('global_climate.nc')

# Read input data to biuld microclimatic models

data<-read.table("Data_Input_NMR.txt", header = T, sep="\\t")
names(data)
head(data)
summary(data)
# Creo vector con los dias de cada mes para un año
diasmes<-rep(1:12, times =c(31,28,31,30,31,30,31,31,30,31,30,31))
diasmes2<-rep(diasmes, times=10)
year<-rep(1:10, each=365)
length(diasmes2)
# creo tabla
filas <-  12*24
names(data)
clima <- data.frame(longitud= rep(NA,filas), latitud=rep(NA,filas) ,Mes=rep(NA,filas), Hora=rep(NA,filas),Tmin=rep(NA,filas),Tmax=rep(NA,filas), Tprom=rep(NA,filas), Hrmin=rep(NA,filas), Hrmax=rep(NA,filas), Hrprom=rep(NA,filas))

ind <- seq(0,dim(clima)[1], 12)

# Function (loop). We defined some metrics for temp and humid.

micromod <- function(x) {
    temp <- micro_global(loc= x[2:3], timeinterval= 365, nyears= 10, soiltype= x[9], elev= x[5], slope= x[6], aspect= x[4], minshade= x[7], maxshade=x[8], Usrhyt= 0.5)$metout
    
    horas <- unique(temp[,2])
    meses<-rep(c(1,2,3,4,5,6,7,8,9,10,11,12), times=24, length.out=288)
    
    for(i in 1:24) {  
        
        tiempox <- subset(data.frame(temp), TIME==horas[i])
        tmean <- apply(tapply(tiempox$TALOC, INDEX=list(year,diasmes2), FUN=mean),2,FUN=mean)
        tmin <- apply(tapply(tiempox$TALOC, INDEX=list(year,diasmes2), FUN=min), 2,FUN=mean)
        tmax <- apply(tapply(tiempox$TALOC, INDEX=list(year,diasmes2), FUN=max), 2,FUN=mean)
        hrmean <- apply(tapply(tiempox$RHLOC, INDEX=list(year,diasmes2), FUN=mean),2,FUN=mean)
        hrmin <- apply(tapply(tiempox$RHLOC, INDEX=list(year,diasmes2), FUN=min),2,FUN=mean)
        hrmax <- apply(tapply(tiempox$RHLOC, INDEX=list(year,diasmes2), FUN=max),2,FUN=mean)
        
        
        clima[(ind[i]+1):ind[i+1],1] <- x[2]
        clima[(ind[i]+1):ind[i+1],2] <- x[3] 
        clima[(ind[i]+1):ind[i+1],3] <- 1:12
        clima[(ind[i]+1):ind[i+1],4] <- horas[i]
        clima[(ind[i]+1):ind[i+1],5] <- tmin
        clima[(ind[i]+1):ind[i+1],6] <- tmax
        clima[(ind[i]+1):ind[i+1],7] <- tmean
        clima[(ind[i]+1):ind[i+1],8] <- hrmin
        clima[(ind[i]+1):ind[i+1],9] <- hrmax
        clima[(ind[i]+1):ind[i+1],10] <- hrmean

    }
    return(clima)
}

# Apply the funtion to the data

models_nmr <- apply(data, 1, micromod)


# check dimensions
dim(models_nmr[[1]])


# Getting the final data
df_nmr <- do.call('rbind', models_nmr)
dim(df_nmr)

# Saving output 

write.table (tabla, file="Results/result_mod_0m_10y.txt", sep= "\\t")
getwd()
dim(tabla)
summary(tabla)

########################################################################
########################################################################

library(dplyr)

data<-read.csv("Data_Input_NMR.csv", header=T)

microclim<-read.table("result_mod_0m_10y.txt", header=T, sep="\\t")
head(microclim)
coll<- data %>% filter(Data=="Preserved specimen")
inat<- data %>% filter(Data=="Observation")

# Loop to get the data according each point ID for inta data

dat1 <- vector("list",dim(inat)[1])
for(i in 1: dim(inat)[1]) {
    dat1[[i]] <- subset(microclim, longitud==inat[i,"Long"] & latitud == inat[i,"Lat"] & Mes==inat[i,"Month"] & Hora==60*(inat[i,"Hour"]), select=c(4,5,6,8))
}

dat1_tab <- do.call("rbind",dat1) %>% as.data.frame() %>% cbind(inat)



# Loop to get the data according each point ID for inta collected data


dat2 <- vector("list",dim(coll)[1])

for(i in 1: dim(coll)[1]) {
    dat2[[i]] <- subset(microclim, longitud==coll[i,"Long"] & latitud == coll[i,"Lat"] & Mes==coll[i,"Month"], select=c(4,5,6,8))
}


dat3<-lapply(dat2, FUN=function(x) x[1:24,])

dat3[[2]]

dat3_tab <- do.call("rbind",dat3) %>% as.data.frame() %>% cbind(inat)

# loop 2

lista2<-list()
for (i in 1:dim(colectas)[1]){
    lista2[[i]]<-data.frame(matrix(coll[i,], nrow=24, ncol=dim(coll)[2], byrow=T))
}

# Final dataframe
df<-do.call("rbind",lista2) %>% as.data.frame()
df2<-cbind(dat3_tab, df)

head(df2)
summary(df2)

df3<-as.matrix(df2) %>% rbind(dat1_tab)
df4<-as.matrix(df2)

##########################################################################
##########################################################################

# Data Analysis

library(ellipsenm)
library(dplyr)
library(raster)
library(ggplot2)
library(rasterVis)
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(cowplot)
library(rgl)


# Generating niche ellipsoids ---------------------------------------------
# Database with observed and modeled records including microclimatic data

both_data_t<-read.csv("D:/Proyectos/En proceso/MAESTRIA/Manuscrito/JOAE_review_coments/Analisis/Modelos_0m/Microclim_values_full_data.csv", header=T)

head(both_data_t)

# Subset for modeling and observed data

mod_t<-both_data_t %>% filter(Data=="Preserved specimen")
inat_t<-both_data_t %>% filter(Data=="Observation")

# Subset for traditional modeling approach

tradmod_t<-mod_t %>% filter(Month==5) %>% filter(Hour%in% c (8, 9, 10, 11)) %>% as.data.frame()  # Datos con res variable


# Niche estimation using activity range from May to August at activity hours following Beck(2005). First we filtered the hours from summer and joint the result with the traditional approach

actmod_t<-mod_t %>% filter(Month%in% c(6,7,8)) %>% filter(Hour%in% c (8,9,10, 11, 16,17,18,19,20,21)) %>% as.data.frame() %>% rbind(tradmod_t)

# For niches ellipsoids we use only 3 variables we exclude some variables that ellipsenm dont use

trad_ellipsoid_t<-tradmod_t %>% dplyr::select(3,4,5,9,10)
acti_ellipsoid_t<-actmod_t %>% dplyr::select(3,4,5,9,10)


# Create overlap ellipsoids to plot and compare their volumes

elip_trad_t<-overlap_object(data=trad_ellipsoid_t, species="Heloderma suspectum",longitude="Long", latitude = "Lat", method = "mve1", level = 95)

elip_acti_t<-overlap_object(data=acti_ellipsoid_t, species="Heloderma suspectum",longitude="Long", latitude = "Lat", method = "mve1", level = 95)


# Ploting the overlap

# Overlaping both niche ellipsoids at 50 a 0 cm
overlap_t<-ellipsoid_overlap(elip_trad_t, elip_acti_t)

plot_overlap(object=overlap_t, niches = c(1,2), data = T,  data_col = c("black","red2"),  background = F, change_labels = FALSE, xlab = "Hrmin", ylab = "Tmax", zlab = "Tmin", legend = F, niche_col = c("black","red2"))

# loop 100 iterations to perform niche overlap 
vec_t<-c(rep(NA,100))
for (i in 1:100){
    vec_t[i]<-ellipsoid_overlap(elip_trad_t, elip_acti_t, overlap_type="all")@full_overlap$overlap 
}

mean(vec_t)

# Estimating the proportion of records out of both ellipsoids using iNaturalist data. Values can vary slightly 

test.trad_ellip<-ellipsoid_fit(trad_ellipsoid_t, longitude="Long", latitude="Lat", method = "mve1", level = 95, raster_layers = NULL)

# Volume= 164.2124
#Slot "centroid":
#    Tmin     Tmax    Hrmin 
#32.02034 39.01495 11.78782 

test.acti_ellip<-ellipsoid_fit(acti_ellipsoid_t, longitude="Long", latitude="Lat", method = "mve1", level = 95, raster_layers = NULL)
# Volume= 9433.393
#Slot "centroid":
#    Tmin     Tmax    Hrmin 
#32.42523 35.86412 20.60118  


# Predicting suitability for iNaturalist occurences. Low Suitability values (=0) represents occurrences out of the ellipsoid

inat_pred<-inat_t %>% dplyr::select(3,4,5)

suit_inat_trad<-predict(test.trad_ellip, projection_variables=inat_pred, prediction = "suitability", truncate = TRUE, return_numeric=T, tolerance = 1e-60, name = NULL, format, overwrite = FALSE, force_return = FALSE)

suit_inat_acti<-predict(test.acti_ellip, projection_variables=inat_pred, prediction = "suitability", truncate = TRUE, return_numeric=T, tolerance = 1e-60, name = NULL, format, overwrite = FALSE, force_return = FALSE)

# Ordering data into one dataframe
inat_suit_results<-cbind(inat_t, suit_inat_trad@suitability, suit_inat_acti@suitability)
head(inat_suit_results)

names(inat_suit_results)
colnames(inat_suit_results)[14]<-"Suitability_trad"
colnames(inat_suit_results)[15]<-"Suitability_act_range"

# Calculate the proportion of the records out of both ellipsoids. Occurrences with suitability==0 are out of ellipsoid 

out_elip_trad<- inat_suit_results %>% filter(Suitability_trad==0)
out_elip_acti<- inat_suit_results %>% filter(Suitability_act_range==0)

prop_oe_trad<-(length(out_elip_trad$Suitability_trad)/(length(inat_t$Data)))
# 86 - 87.04% de los registros quedan por fuera con la aproximación tadicional
prop_oe_acti<-(length(out_elip_acti$Suitability_act_range)/(length(inat_t$Data)))
# 5.09 - 6% de los resgistros quedan por fuera con el elipsoide horario


# Generating response curves ----------------------------------------------
# To achive this, we take the range o the variables and test each variable when other two are always equal to the mean

# finding minimum, maximum and mean values of the variables using the occurences for modeling according selected elipsoid. 

# minimum
niche<-actmod_t

min_tmax<-min(niche$Tmax)
min_tmin<-min(niche$Tmin)
min_hrmin<-min(niche$Hrmin)

# maximum
max_tmax<-max(niche$Tmax)
max_tmin<-max(niche$Tmin)
max_hrmin<-max(niche$Hrmin)

# Setting the range of all variables to the same  number of rows 
tmin<-c(seq(10, 55, 1.5)) %>% as.data.frame()
tmax<-c(seq(10, 55, 1.5)) %>% as.data.frame()
hrmin<-c(seq(4, 68, 2.1)) %>% as.data.frame()
64/31
# All variables have 29 rows
length(tmin$.)
length(tmax$.)
length(hrmin$.)

# Calculating the mean of all variables, then create a dataframe with the same number of rows
med_tmax<-rep(mean(niche$Tmax), 31) %>% as.data.frame()
med_tmin<-rep(mean(niche$Tmin), 31) %>% as.data.frame()
med_hrmin<-rep(mean(niche$Hrmin), 31) %>% as.data.frame()

# Creating a dataframe to test each variable. Each variable will have that range o values, and the other variables will be equal to their mean

tmin_data<-cbind(tmin, med_tmax, med_hrmin)
head(tmin_data)
colnames(tmin_data)[1]="Tmin"
colnames(tmin_data)[2]="Tmax"
colnames(tmin_data)[3]="Hrmin"

# dataframe para probar tmax
tmax_data<-cbind(med_tmin, tmax, med_hrmin)
head(tmax_data)
colnames(tmax_data)[1]="Tmin"
colnames(tmax_data)[2]="Tmax"
colnames(tmax_data)[3]="Hrmin"

# dataframe para probar hrmin
hrmin_data<-cbind(med_tmin, med_tmax, hrmin)
head(hrmin_data)
colnames(hrmin_data)[1]="Tmin"
colnames(hrmin_data)[2]="Tmax"
colnames(hrmin_data)[3]="Hrmin"

# Projecting the niche ellipsoid selected into each variable dataframe

suit_tmin<-predict(test.acti_ellip, projection_variables=tmin_data, prediction = "both", truncate = TRUE, return_numeric=T, tolerance = 1e-60, name = NULL, format, overwrite = FALSE, force_return = FALSE) 

suit_tmax<-predict(test.acti_ellip, projection_variables=tmax_data, prediction = "both", truncate = TRUE, return_numeric=T, tolerance = 1e-60, name = NULL, format, overwrite = FALSE, force_return = FALSE)

suit_hrmin<-predict(test.acti_ellip, projection_variables=hrmin_data, prediction = "both", truncate = TRUE, return_numeric=T, tolerance = 1e-60, name = NULL, format, overwrite = FALSE, force_return = FALSE)

# Results 
s_tmin<-as.data.frame(suit_tmin@suitability)
dis_tmin<-as.data.frame(suit_tmin@mahalanobis)
rc_tmin<-cbind(tmin, s_tmin, dis_tmin)

s_tmax<-as.data.frame(suit_tmax@suitability)
dis_tmax<-as.data.frame(suit_tmax@mahalanobis)
rc_tmax<-cbind(tmax, s_tmax, dis_tmax)

s_hrmin<-as.data.frame(suit_hrmin@suitability)
dis_hrmin<-as.data.frame(suit_hrmin@mahalanobis)
rc_hrmin<-cbind(hrmin, s_hrmin, dis_hrmin)

# Ploting
tmin_rc<-ggplot()+
    geom_line(data=rc_tmin, aes(x=., y=suit_tmin@suitability), size=1.3, col="black")+
    theme_bw()+
    geom_vline(xintercept=17.4, color= "gray61", size=1, linetype= "dashed")+
    geom_vline(xintercept=36.8, color= "gray61", size=1, linetype= "dashed")+
    theme(axis.title=element_text(size=9), axis.text.x  = element_text(size = 9), axis.text.y  = element_text(size = 9))+
    labs(x="Minimum temperature (°C)", y="Suitability", size=10)
tmin_rc

tmax_rc<-ggplot()+
    geom_line(data=rc_tmax, aes(x=., y=suit_tmax@suitability), size=1.3, col="black")+
    theme_bw()+
    geom_vline(xintercept=17.4, color= "gray61", size=1, linetype= "dashed")+
    geom_vline(xintercept=36.8, color= "gray61", size=1, linetype= "dashed")+
    theme(axis.title=element_text(size=9), axis.text.x  = element_text(size = 9), axis.text.y  = element_text(size = 9))+
    labs(x="Maximum temperature (°C)", y="Suitability", size=10)
tmax_rc

hrmin_rc<-ggplot()+
    geom_line(data=rc_hrmin, aes(x=., y=suit_hrmin@suitability), size=1.3, col="black")+
    theme_bw()+
    theme(axis.title=element_text(size=9), axis.text.x  = element_text(size = 9), axis.text.y  = element_text(size = 9))+
    labs(x="Minimum relative humidity (%)", y="Suitability", size=10)

hrmin_rc

res_curve_plot<-plot_grid(hrmin_rc, tmin_rc, tmax_rc, nrow = 1)

# Predicting suitability for model occurrences using the selected ellipsoid -----------
suit_datos<- both_data_t %>% dplyr::select(3,4,5) %>% as.data.frame()
head(suit_datos)

suit_occ<-predict(test.acti_ellip, projection_variables=suit_datos, prediction = "both", truncate = TRUE, return_numeric=T, tolerance = 1e-60, name = NULL, format, overwrite = FALSE, force_return = FALSE)

suit_occ<-as.data.frame(suit_occ@suitability)
result_occ<-cbind(both_data_t, suit_occ)
names(result_occ)
colnames(result_occ)[14]<-"suitability"
names(result_occ)

write.csv(result_occ, "D:/Proyectos/En proceso/MAESTRIA/Manuscrito/JOAE_review_coments/Final_products/Suitability_values_act_model_0m_FINAL.csv")


# Validating the model using the observed data and the seasonal hypothesis of activity following Beck and physiological limits --------------------------------

# Suitability data for occurences by season (not include occurrences outside the ellipsoid = 5%)
summer_m<-result_occ %>% filter(Month%in% c(6,7,8)) %>% filter(Hour%in% c(8:23)) %>% filter(suitability>0) 
spring_m<- result_occ %>% filter(Month%in% c(3,4,5)) %>% filter(Hour%in% c(8:19)) %>% filter(suitability>0)
fallwin_m<-result_occ %>% filter(Month%in% c(1,12,11,10,9,2)) %>% filter(Hour%in% c(9:16)) %>% filter(suitability>0)
head(fallwin_m)
# load Becks data
# preparing Becks data

beck<-read.csv("D:/Proyectos/En proceso/MAESTRIA/Manuscrito/Data/Becks_seasonal_observation.csv", header=T)

beck_sum<- beck %>% filter(Season=="Summer")
beck_spr<- beck %>% filter(Season=="Spring")
beck_fw<- beck %>% filter(Season=="Fall_winter")


# loading and subseting Data from iNaturalist by season (proportion of observation by hours within each season)

inat_prop<-read.csv("D:/Proyectos/En proceso/MAESTRIA/Manuscrito/JOAE_review_coments/Final_products/Prop_mod_seasons.csv", header=T)
summer_inat<-inat_prop %>% filter(Season=="Summer") 
spring_inat<- inat_prop %>% filter(Season=="Spring") 
fallwin_inat<-inat_prop %>% filter(Season=="Fall and Winter") 


# Plot becks patterns according the evaluated hour range for each season

becks_inv<-ggplot()+
    geom_smooth(data= beck_fw, aes(x= Hour, y= Becks_obs),col="black", size=1.5,se=F)+
    theme_bw()+
    labs(title = "Fall/Winter", y="Beck's Activity\\nObservations", x=" ", size=10)+
    theme(axis.title=element_text(size=9), axis.text.x  = element_blank(), axis.text.y  = element_text(size = 9), plot.title = element_text(hjust = 0.5))+
    coord_cartesian(xlim=c(9, 16), ylim=c(0, 35))

becks_sum<-ggplot()+
    geom_smooth(data= beck_sum, aes(x= Hour, y= Becks_obs),col="black", size=1.5,se=F)+
    theme_bw()+
    labs(title="Summer", y="Beck's Activity\\nObservations", x=" ", size=10)+
    theme(axis.title=element_text(size=9), axis.text.x  = element_blank(), axis.text.y  = element_text(size = 9), plot.title = element_text(hjust = 0.5))+
    coord_cartesian(ylim=c(0, 35))

becks_spr<-ggplot()+
    geom_smooth(data= beck_spr, aes(x= Hour, y= Becks_obs),col="black", size=1.5,se=F)+
    theme_bw()+
    labs(title="Spring", y="Beck's Activity\\nObservations", x=" ", size=10)+
    theme(axis.title=element_text(size=9), axis.text.x  = element_blank(), axis.text.y  = element_text(size = 9), plot.title = element_text(hjust = 0.5))+
    coord_cartesian(ylim=c(0, 35))

# Plot the suitability trend of the modeled data and proportion of inaturalist records.

Sui_inta_wint<-ggplot()+
    geom_bar(data=fallwin_inat, aes(x=Hour, y=Freq), col="black", fill= "gray25", stat="identity", width = 0.8)+
    geom_smooth(data=fallwin_m, aes(x=Hour, y=suitability), size=1.5, col="red3")+
    theme_bw()+
    labs(x="", y="Suitability(red line) \\n Proportion of occurences (gray bars)", size=10)+
    theme(axis.title=element_text(size=9), axis.text.x  = element_blank(), axis.text.y  = element_text(size = 9))+
    #scale_y_continuous(sec.axis= sec_axis(~., name= "Obsevations proportion"))+
    coord_cartesian(xlim=c(8.7, 16), ylim = c(0,0.4))

Sui_inta_sum<-ggplot()+
    geom_bar(data=summer_inat, aes(x=Hour, y=Freq), col="black", fill= "gray25", stat="identity", width = 0.8)+
    geom_smooth(data=summer_m, aes(x=Hour, y=suitability), size=1.5, col="red3")+
    theme_bw()+
    labs(x="", y="Suitability (red line) \\n Proportion of occurences (gray bars)", size=10)+
    theme(axis.title=element_text(size=9), axis.text.x  = element_blank(), axis.text.y  = element_text(size = 9))+
    #scale_y_continuous(sec.axis= sec_axis(~., name= "Obsevations proportion"))+
    coord_cartesian(ylim = c(0,0.55))

Sui_inta_spr<-ggplot()+
    geom_bar(data=spring_inat, aes(x=Hour, y=Freq), col="black", fill= "gray25", stat="identity", width = 0.8)+
    geom_smooth(data=spring_m, aes(x=Hour, y=suitability), size=1.5, col="red3")+
    theme_bw()+
    labs(x="", y="Suitability (red line) \\n Proportion of occurences (gray bars)", size=10)+
    theme(axis.title=element_text(size=9), axis.text.x  = element_blank(), axis.text.y  = element_text(size = 9))+
    #scale_y_continuous(sec.axis= sec_axis(~., name= "Obsevations proportion"))+
    coord_cartesian(ylim = c(0,0.45))

# Boxplot maximum temperature of both datasets with temperature activity range as a reference

fw_box<-ggplot()+
    geom_hline(yintercept=17.4, color= "gray61", size=1.1, linetype= "dashed")+
    geom_hline(yintercept=36.8, color= "gray61", size=1.1, linetype= "dashed")+
    geom_boxplot(data=fallwin_m, aes(x=as.factor(Hour), y=Tmax, fill=Data), col="black", position="dodge2", outlier.size =1, outlier.alpha=0.5)+
    theme_bw()+
    scale_fill_manual(values = c("gray25", "red3"))+
    labs(x="Hour", y="Maximum\\ntemperature (°C)", size=10)+
    coord_cartesian(ylim=c(0, 50))+
    theme(axis.title=element_text(size=9), axis.text.x  = element_text(size = 9), axis.text.y  = element_text(size = 9), legend.position = "none")

sum_box<-ggplot()+
    geom_hline(yintercept=17.4, color= "gray61", size=1.1, linetype= "dashed")+
    geom_hline(yintercept=36.8, color= "gray61", size=1.1, linetype= "dashed")+
    geom_boxplot(data=summer_m, aes(x=as.factor(Hour), y=Tmax, fill=Data), col="black", position="dodge2", outlier.size =1, outlier.alpha=0.5)+
    theme_bw()+
    scale_fill_manual(values = c("gray25", "red3"))+
    labs(x="Hour", y="Maximum\\ntemperature (°C)", size=10)+
    coord_cartesian(ylim=c(0, 50))+
    theme(axis.title=element_text(size=9), axis.text.x  = element_text(size = 9), axis.text.y  = element_text(size = 9), legend.position = "none")

spr_box<-ggplot()+
    geom_hline(yintercept=17.4, color= "gray61", size=1.1, linetype= "dashed")+
    geom_hline(yintercept=36.8, color= "gray61", size=1.1, linetype= "dashed")+
    geom_boxplot(data=spring_m, aes(x=as.factor(Hour), y=Tmax, fill=Data), col="black", position="dodge2",  outlier.size =1, outlier.alpha=0.5)+
    theme_bw()+
    scale_fill_manual(values = c("gray25", "red3"))+
    labs(x="Hour", y="Maximum\\ntemperature (°C)", size=10)+
    coord_cartesian(ylim=c(0, 50))+
    theme(axis.title=element_text(size=9), axis.text.x  = element_text(size = 9), axis.text.y  = element_text(size = 9), legend.position = "none")


# Correlation analisys

# vector inat occs by hour
inat_sum <- as.vector(summer_inat$Count)
inat_spr <- as.vector(spring_inat$Count)
inat_fw  <- as.vector(fallwin_inat$Count)


# vector mean suitability by hours
summer_suit<-result_occ %>% filter(Data=="Preserved specimen") %>%  filter(Month%in% c(6,7,8)) %>% filter(Hour%in% c(8:23)) %>% group_by(Hour) %>% summarise(Mean_suit=mean(suitability)) %>% as.data.frame() %>% dplyr::select(2)
summ_suit<-as.vector(summer_suit$Mean_suit)

spring_suit<- result_occ %>% filter(Data=="Preserved specimen") %>% filter(Month%in% c(3,4,5)) %>% filter(Hour%in% c(8:19)) %>% filter(suitability>0) %>% group_by(Hour) %>% summarise(Mean_suit=mean(suitability)) %>% as.data.frame %>% dplyr::select(2) %>% as.vector()
spri_suit<-as.vector(spring_suit$Mean_suit)

fallwin_suit<-result_occ %>% filter(Data=="Preserved specimen") %>%filter(Month%in% c(1,12,11,10,9,2)) %>%  filter(Hour%in% c(9:16)) %>% group_by(Hour) %>% summarise(Mean_suit=mean(suitability)) %>% as.data.frame() %>% dplyr::select(2) %>% as.vector()
fawin_suit<-as.vector(fallwin_suit$Mean_suit)

# vector beck observations
beck_sum<- beck %>% filter(Season=="Summer") %>% dplyr::select(3) 
bsum<-as.vector(beck_sum[[1]])

beck_spr<- beck %>% filter(Season=="Spring") %>% dplyr::select(3)
bspr<-as.vector(beck_spr[[1]])

beck_fw<- beck %>% filter(Season=="Fall_winter") %>% dplyr::select(3)
bfw<-as.vector(beck_fw[[1]])

# Correaltion test i naturalist and suitability
i_s_sum<-cor.test(inat_sum, summ_suit, method = "spearman")
i_s_spr<-cor.test(inat_spr, spri_suit, method = "spearman")
i_s_fw <-cor.test(inat_fw, fawin_suit, method = "spearman")

# correlation test suitability and becks data
b_s_sum<-cor.test(bsum, summ_suit, method = "spearman")
b_s_spr<-cor.test(bspr, spri_suit, method = "spearman")
b_s_fw <-cor.test(bfw, fawin_suit, method = "spearman")


#cross correlation using spearman
# Correaltion test i naturalist and suitability

i_s_sum_1<-cor.test(inat_sum[2:16], summ_suit[1:15], method = "spearman")
i_s_sum_2<-cor.test(inat_sum[3:16], summ_suit[1:14], method = "spearman")
i_s_sum_3<-cor.test(inat_sum[4:16], summ_suit[1:13], method = "spearman")
i_s_sum_m1<-cor.test(inat_sum[1:15], summ_suit[2:16], method = "spearman")
i_s_sum_m2<-cor.test(inat_sum[1:14], summ_suit[3:16], method = "spearman")
i_s_sum_m3<-cor.test(inat_sum[1:13], summ_suit[4:16], method = "spearman")

i_s_spr_1<-cor.test(inat_spr[2:12], spri_suit[1:11], method = "spearman")
i_s_spr_2<-cor.test(inat_spr[3:12], spri_suit[1:10], method = "spearman")
i_s_spr_3<-cor.test(inat_spr[4:12], spri_suit[1:9], method = "spearman")
i_s_spr_m1<-cor.test(inat_spr[1:11], spri_suit[2:12], method = "spearman")
i_s_spr_m2<-cor.test(inat_spr[1:10], spri_suit[3:12], method = "spearman")
i_s_spr_m3<-cor.test(inat_spr[1:9], spri_suit[4:12], method = "spearman")

i_s_fw_1<-cor.test(inat_fw[2:8], fawin_suit[1:7], method = "spearman")
i_s_fw_2<-cor.test(inat_fw[3:8], fawin_suit[1:6], method = "spearman")
i_s_fw_3<-cor.test(inat_fw[4:8], fawin_suit[1:5], method = "spearman")
i_s_fw_m1<-cor.test(inat_fw[1:7], fawin_suit[2:8], method = "spearman")
i_s_fw_m2<-cor.test(inat_fw[1:6], fawin_suit[3:8], method = "spearman")
i_s_fw_m3<-cor.test(inat_fw[1:5], fawin_suit[4:8], method = "spearman")


# Correaltion test Becks data and suitability

b_s_sum_1<-cor.test(bsum[2:16], summ_suit[1:15], method = "spearman")
b_s_sum_2<-cor.test(bsum[3:16], summ_suit[1:14], method = "spearman")
b_s_sum_3<-cor.test(bsum[4:16], summ_suit[1:13], method = "spearman")
b_s_sum_m1<-cor.test(bsum[1:15], summ_suit[2:16], method = "spearman")
b_s_sum_m2<-cor.test(bsum[1:14], summ_suit[3:16], method = "spearman")
b_s_sum_m3<-cor.test(bsum[1:13], summ_suit[4:16], method = "spearman")

b_s_spr_1<-cor.test(bspr[2:12], spri_suit[1:11], method = "spearman")
b_s_spr_2<-cor.test(bspr[3:12], spri_suit[1:10], method = "spearman")
b_s_spr_3<-cor.test(bspr[4:12], spri_suit[1:9], method = "spearman")
b_s_spr_m1<-cor.test(bspr[1:11], spri_suit[2:12], method = "spearman")
b_s_spr_m2<-cor.test(bspr[1:10], spri_suit[3:12], method = "spearman")
b_s_spr_m3<-cor.test(bspr[1:9], spri_suit[4:12], method = "spearman")

b_s_fw_1<-cor.test(bfw[2:8], fawin_suit[1:7], method = "spearman")
b_s_fw_2<-cor.test(bfw[3:8], fawin_suit[1:6], method = "spearman")
b_s_fw_3<-cor.test(bfw[4:8], fawin_suit[1:5], method = "spearman")
b_s_fw_m1<-cor.test(bfw[1:7], fawin_suit[2:8], method = "spearman")
b_s_fw_m2<-cor.test(bfw[1:6], fawin_suit[3:8], method = "spearman")
b_s_fw_m3<-cor.test(bfw[1:5], fawin_suit[4:8], method = "spearman")




