#ClarkWolf_et_al_2022_code.R
#
# This script will reproduce the supplemental figures presented in 
# Wolf et al. (2022). Performed in R version 4.1.1.
#
# Wolf, K.D., Higuera, P.E., Davis, K.T. 2022. Conifer seedling demography reveals
# mechanisms of initial forest resilience to wildfires in the northern Rocky Mountains.
# Forest ecology and management.
#
# Data File Requirements: 
# (1) ClarkWolf_et_al_2022_PlotData.csv - plot geospatial data, field measurements, with all predictor variables used in statistical modeling (averaged to plot level)
# (2) ClarkWolf_et_al_2022_Microclimate.csv - microclimate data (temperature and VPD) aggregated to daily timestep, with fitted values based on field measurements at a subset of sites
# (3) ClarkWolf_et_al_2022_SoilData.csv - soil data, including ammonium (NH4) and nitrate (NO3) in mineral (M) and organic (O) soil layers by year
# (4) ClarkWolf_et_al_2022_OverstoryData.csv - plot-level tree density and basal area measurements by species
# (5) ClarkWolf_et_al_2022_SubplotData.csv - ground cover and canopy cover measurements in each subplot by year
# (6) ClarkWolf_et_al_2022_RecruitmentData.csv - annual counts of germination-year seedlings in each plot by species
# (7) ClarkWolf_et_al_2022_RegenMonitoring.csv - individual seedling height, age, and living status by year and species for subplots within plots 
# (8) ClarkWolf_et_al_2022_gridMET.csv - daily temperature and vpd data extracted from gridMET (Abatzoglou et al. 2013, climatologylab.org/gridmet)
#
#
# Created by: Kyra D. Wolf
# Created on: December 2021
# University of Montana, PaleoEcology and Fire Ecology Lab
#
# kyra.clark-wolf@umontana.edu
# philip.higuera@umontana.edu


# Set working directory ---------------------------------------------------
#setwd("~/Desktop") # The directory should be changed to reflect your data file location
#setwd("L:/3_labMembers/Kyra/Manuscripts/Seedlings_GRIN/Data_archive")

# Set up libraries, data, and functions --------------------------------------------------------------

# Load required packages
x = c("ggplot2","dplyr","tidyr","mgcv","lme4","lmerTest","effects",
      "patchwork","ade4","FactoMineR","factoextra","gridExtra", "zoo", 
      "sp","rgdal","rgeos","maptools","raster","rasterVis","dismo",
      "jtools","gstat","ggstance","performance", "MuMIn", "MASS", "scales",
      "pROC","DHARMa", "car","stringr", "officer","rvg","rstatix")

# If packages have not been previously installed, uncomment the following code,
# or install these packages manually:
# lapply(x,install.packages,character.only=T)

lapply(x,library,character.only = T)

############Read in data

#plot data
plot <- read.csv("ClarkWolf_et_al_2022_PlotData.csv")

#seedling monitoring
sd <- read.csv("ClarkWolf_et_al_2022_RegenMonitoring.csv")

#subplot data
sp <- read.csv("ClarkWolf_et_al_2022_SubplotData.csv")

#microclimate and gridded climate data
mc <- read.csv("ClarkWolf_et_al_2022_Microclimate.csv")
gm <- read.csv("ClarkWolf_et_al_2022_gridMET.csv")

#recruitment data
rc <- read.csv("ClarkWolf_et_al_2022_RecruitmentData.csv")

#soil data
sl <- read.csv("ClarkWolf_et_al_2022_SoilData.csv")
  
#overstory data
ov <- read.csv("ClarkWolf_et_al_2022_OverstoryData.csv")

theme_set(theme_classic())
p = c("#D55E00",  "#F0E442", "#009E73")

####Functions

#Crossvalidation function, for glmms
CV_plot = function(model,response,data){
  P = levels(factor(data$Plot))
  acc.metrics<-data.frame(Plot = P, RMSE = NA)
  for(i in P){
    data.train<-data[data$Plot != i,]
    data.valid<-data[data$Plot == i,]
    temp.model<-update(model,data=data.train)
    pred<-predict(temp.model,newdata=data.valid, re.form=NA) 
    obs<-subset(data.valid,select=response)
    rmse<-sqrt(sum((pred-obs)^2)/length(data.valid[,1]))
    acc.metrics[acc.metrics$Plot == i, 2]<-rmse
  }
  accuracy<<-as.numeric(acc.metrics$RMSE)
  mean.acc<<-mean(accuracy)
  sd.acc<<-sd(accuracy)
  print(mean.acc)
  print(sd.acc)
  hist(accuracy)
}

####Figure A.2 ##################################################################################
###Modeling daily microclimate conditions based on field measurements from a subset of sites

######################################Set up data

###Select and rescale predictor variables
#canopy cover
sb$Canopy.tot_sqrt <- sqrt(sb$Canopy.tot) #reduce skewness
sb1 <- pivot_wider(sb[,c("Plot","Postfire_year","Canopy.tot_sqrt")],
                   names_from = Postfire_year, values_from = Canopy.tot_sqrt,names_prefix = "Canopy.tot_")
#plot data
preds <- merge(plot[,c("Plot","Axis1","Axis2","DEF","HLI","Elev")],
               sb1[,c("Plot","Canopy.tot_1","Canopy.tot_2","Canopy.tot_3")], by = "Plot")
#gridMET data
preds <- merge(preds,gm,by = "Plot",all=T)
#rescale and center predictors
preds[,c("Canopy.tot_1","Canopy.tot_2","Canopy.tot_3","Axis1","Axis2","DEF","HLI","Elev","Tmax.gm")] = 
  scale(preds[,c("Canopy.tot_1","Canopy.tot_2","Canopy.tot_3","Axis1","Axis2","DEF","HLI","Elev","Tmax.gm")],center =TRUE)
preds[which(is.na(preds$Canopy.tot_1)),"Canopy.tot_1"] <- #fill in a few missing values with averages
  apply(preds[which(is.na(preds$Canopy.tot_1)),c("Canopy.tot_2","Canopy.tot_3")],1,mean)

###Merge microclimate measurements with predictors
preds$Date <- as.Date(preds$Date)
mc$Date <- as.Date(mc$Date)
mc <- separate(mc,"Plot",into =c("Fire","Rep","Severity"),sep = "_",remove=F)
d=merge(mc[-which(is.na(mc$Meas_Tmax)),],preds,by=c("Plot","Date")) 
d$Canopy.tot <- ifelse(d$Year==2018,d$Canopy.tot_1,ifelse(d$Year==2019, d$Canopy.tot_2, d$Canopy.tot_3))

train <- d[d$Year == 2019,] #select training data
val <- d[d$Year == 2018,] #select validation data

##########################################Temporal autocorrelation & data subsetting

#vector of plot names
v = unique(train$Plot)

#empty dataframe of the maximum lag with significant autocorrelation
lags = data.frame(Plot = v,lag.Tmax = NA)

#loop through to calculate significant lag for each plot (max Temperature)
for (i in v){
  d1 = train[train$Plot == i,]
  a = acf(d1$Meas_Tmax, main = i)
  temp = setNames(data.frame(matrix(NA,nrow = dim(a$lag),ncol = 4)), c("lag","cor","cutoff","sig"))
  temp[,2] = a$acf[,,1]
  temp[,1] = a$lag[,,1]
  temp$cutoff = 1.96/sqrt(dim(d1)[1]-temp$lag)
  temp$sig = ifelse(temp$cor>temp$cutoff,1,0)
  lags[lags$Plot==i,"lag.Tmax"] = sum(temp$sig) -1
}
mean(lags$lag.Tmax)

#######Subset data to retain one day out of every 6

#how many dates are there
dates = seq(from = min(as.Date(train$Date)),to = max(as.Date(train$Date)),by=1)
l = length(dates)

#create a vector to select 1 out of every 6 days)
keep = c(rep(c(0,1,0,0,0,0),l/6),0,0)
ds3 = merge(train,data.frame(Date=dates,keep=keep), by = "Date") %>% filter(keep == 1)

######Fit random structure
m1 = lmer(Meas_Tmax~(DEF+HLI+poly(Tmax.gm,2,raw=T)+Canopy.tot+Elev+Axis1+Axis2)^2 + (1|Plot),data = ds1)
m2 = update(m1,.~.-(1|Plot) + (1|Fire/Plot))
m3 = update(m1,.~.-(1|Plot) + (1|Fire)) 

#Compare different random structures
BIC(m1,m2,m3) #m1 has lowest BIC
AIC(m1,m2,m3) #m1 has lowest AIC

######Fit fixed structure
m.full = update(m1,REML = F,data=ds3)
s = get_model(step(m.full))
summary(s)

#Further reduce fixed structure to minimize RMSE
m1 = lmer(Meas_Tmax~Canopy.tot+Axis1+Axis2+DEF+HLI+Elev+poly(Tmax.gm,2,raw=T)+
            Axis1:Elev+
            Axis2:HLI+
            DEF:Elev+
            HLI:poly(Tmax.gm,2,raw=T)+
            HLI:Canopy.tot+
            HLI:Elev+
            poly(Tmax.gm,2,raw=T):Canopy.tot+
            Canopy.tot:Elev+
            (1|Plot),data=ds3,REML=F)
summary(m1)
CV_plot(m1,"Meas_Tmax",ds3) #3.10

m2 = lmer(Meas_Tmax~Canopy.tot+Axis1+Axis2+HLI+Elev+poly(Tmax.gm,2,raw=T)+
            Axis2:HLI+
            HLI:Tmax.gm+
            HLI:Elev+
            Canopy.tot:Elev+
            (1|Plot),data=ds3,REML=F)
anova(m1,m2)
CV_plot(m2,"Meas_Tmax",ds3) #3.04

m3 = lmer(Meas_Tmax~Canopy.tot+HLI+Elev+poly(Tmax.gm,2,raw=T)+Axis1+
            HLI:Tmax.gm+
            HLI:Elev+
            (1|Plot),data=ds3,REML=F)
anova(m2,m3)
CV_plot(m3,"Meas_Tmax",ds3) #3.00

#Validate using 2018 microclimate measurements
sqrt(mean((val$Meas_Tmax-predict(m3,newdata=val,re.form=NA))^2)) #2.60

#Check residuals
plot(m3, type = c("p", "smooth"))
hist(resid(m3))

#Predict across all plots and years
fitted <- predict(m3,newdata = preds,allow.new.levels =TRUE)

#################Plotting

####Fig. A.2
gm$Date <- as.Date(gm$Date)
dd = merge(mc,gm,by=c("Plot","Date"))
dd = dd %>% pivot_longer(cols = c("Meas_Tmax","fitted_Tmax","Tmax.gm"),names_to = "type",values_to = "Tmax")
dd$type = as.factor(dd$type)
levels(dd$type) = c("fitted","measured","gridMET")

p1 <- ggplot(dd[dd$Year==2019,], aes(x = Date, y =  Tmax, color = type))+facet_wrap(~Severity)+
  geom_smooth(method = "gam",lwd = 1)+
  labs(color = NULL, title = "2019 calibration data (n = 40 plots)")

p2 <- ggplot(dd[dd$Year==2018,], aes(x = Date, y =  Tmax, color = type))+facet_wrap(~Severity)+
  geom_smooth(method = "gam",lwd = 1)+
  labs(color = NULL, title = "2018 validation data (n = 15 plots)")

p1/p2
ggsave("Figure_A-2.png", height = 10, width = 8, units = "in")

#
#
#
####Figure A.3 ####################################################################
####Principal components analysis to characterize field-based fire severity

###pre-treat data to reduce skewness
#sqrt transform ground cover data & canopy cover & basal area
d <- pivot_longer(plot,cols = starts_with("cover."),names_prefix = "cover.", #convert to long format to make this faster
                  names_to = c("CoverType","type"),names_sep = "_",values_to = "cover")
d$cover = sqrt(d$cover)
d = pivot_wider(d,values_from = cover, names_from = CoverType)

d$DeadCanopy = sqrt(d$Canopy.dead_Avg)
d$LiveCanopy = sqrt(d$Canopy.green_Avg)
d$LiveBA = sqrt(d$LiveBA)
d$DeadBA = sqrt(d$DeadBA)
d$CanopyTotal = sqrt(d$Canopy.tot_Avg)

###run the PCA
PCA<-dudi.pca(da[,c("BG","Wood","Litter","Moss","Forb","Grass","Shrub","DeadCanopy","LiveCanopy",
                    "LiveBA","DeadBA","SH","mort", "PP")], center = TRUE, scale=TRUE, scannf=FALSE, nf = 2)

###Appendix S1: Fig. S2
b<-fviz_pca_biplot(PCA,
                   geom.ind ="point",
                   pointsize=1,
                   labelsize=4,
                   palette = c(),
                   col.var="black",
                   col.ind="grey",
                   repel=TRUE,
                   mean.point=FALSE,
                   arrowsize=.5)

p1 = b+labs(title = "Averaged")+theme(axis.text=element_text(size=11),axis.title=element_text(size=13))#,plot.title = element_blank())
p1

b<-fviz_pca_biplot(PCA,
                   geom.ind ="point",
                   pointsize=2,
                   labelsize=0,
                   palette =c("#D55E00", "#E69F00", "#009E73"),
                   col.var=NA,
                   repel=TRUE,
                   mean.point=FALSE,
                   arrowsize=0,
                   habillage=da$Sev,
                   addEllipses = TRUE)
p2 = b+ labs(title = "")+theme(axis.text=element_text(size=12),axis.title=element_text(size=13), legend.text = element_text(size = 12),
                               legend.title = element_blank(),legend.position = c(.15,.1), 
                               legend.background = element_rect(color = "black"), legend.margin = margin(c(0,0.2,0.1,0.1),unit = "cm"))
p2

p1+p2

ggsave("AppendixS1_FigS2.png", dpi = 600, width = 10, height = 5, units = "in" )

#Axis 1 and 2 values (PCA$li) are used as field-based fire-severity metrics in further analyses

#
#
#
####Figure B.1 - B.2 ##################################################################
####Analysis of seedling density in burned and unburned plots

#calculate total and species density of regenerating seedlings after fire
plot$Regen_StudySpp_ct <- (plot$Regen_Abies_ct+plot$Regen_LAOC_ct+plot$Regen_PICO_ct+plot$Regen_PIEN_ct+plot$Regen_PIPO_ct+plot$Regen_PSME_ct)
regen_dens <- plot %>%  dplyr::select(c("Plot","Severity","TransectSize","Regen_StudySpp_ct","Regen_AllSpp_ct","Regen_Abies_ct","Regen_LAOC_ct","Regen_PICO_ct","Regen_PIEN_ct",
                                        "Regen_PIPO_ct","Regen_PSME_ct"))%>%
  mutate(Regen_StudySpp = Regen_StudySpp_ct/TransectSize*10000,
         Regen_AllSpp = Regen_AllSpp_ct/TransectSize*10000,
         Regen_Abies = Regen_Abies_ct/TransectSize*10000,
         Regen_LAOC = Regen_LAOC_ct/TransectSize*10000,
         Regen_PICO = Regen_PICO_ct/TransectSize*10000,
         Regen_PIEN = Regen_PIEN_ct/TransectSize*10000,
         Regen_PIPO = Regen_PIPO_ct/TransectSize*10000,
         Regen_PSME = Regen_PSME_ct/TransectSize*10000) %>%
  pivot_longer(cols = c("Regen_AllSpp","Regen_Abies","Regen_LAOC","Regen_PICO","Regen_PIEN","Regen_PIPO","Regen_PSME"),
               values_to = "Regen",names_to = "Species",names_prefix = "Regen_")

regen_dens$Status = ifelse(regen_dens$Severity == "Unburned","Unburned","Burned")
regen_dens$Species <- factor(regen_dens$Species,levels = c("AllSpp","Abies","LAOC","PICO","PIEN","PIPO","PSME"))

######Calculate how many plots have seedling density equal to or exceeding reconstructed prefire tree density
nrow(plot[plot$Regen_AllSpp_ct/plot$TransectSize*10000 >= plot$Mature_AllSpp,]) #42/69 or 60% all plots at replacement
b <- plot[plot$Severity != "Unburned",]
nrow(b[which(b$Regen_AllSpp_ct/b$TransectSize*10000 >= b$Mature_AllSpp),]) #32/47 or 68% burned plots at replacement
#25/47 or 53% burned plots at 2x replacement

###################Figure B.1
sim <- data.frame(x = seq(100,10000,100),y = seq(100,10000,100))
g1 <- ggplot(b,aes(x=Mature_AllSpp,y=Regen_AllSpp_ct/TransectSize*10000))+
  geom_point()+
  labs(y =expression("Post-fire seedling density (# ha"^-1*")"),x=expression("Pre-fire tree density (# ha"^-1*")"))+
  scale_y_log10(labels = comma, limits = c(10,550000))+
  scale_x_log10(labels = comma,limits = c(100,10000))+
  geom_line(data = sim, aes(x=x, y=y),linetype = "dashed")+
  annotate("text",label = "1:1",x=7500, y = 12000)
g1
#ggsave("Figure_B-1.png",dpi = 600, width = 3.5, height = 3, units = "in")


#######Compare regeneration density with pre-fire tree density by species

#calculate total and species density of mature trees (live and dead) to estimate pre-fire overstory density
ov$Mature <- ov$dens_Mature_Dead+ov$dens_Mature_Live
m <- ov %>% pivot_wider(id_cols = Plot,names_from = Species,values_from = Mature,names_prefix = "Mature_") %>%
  mutate(Mature_Abies = Mature_ABGR+Mature_ABLA,
         Mature_StudySpp = Mature_Abies+Mature_LAOC+Mature_PICO+Mature_PIEN+Mature_PIPO+Mature_PSME) %>%
  pivot_longer(cols = c("Mature_Abies","Mature_LAOC","Mature_PICO","Mature_PIEN","Mature_PIPO","Mature_PSME"),
               values_to = "Mature",names_to = "Species",names_prefix = "Mature_")

#calculate species composition of trees & regenerating seedlings
r <- regen_dens[-which(regen_dens$Species == "AllSpp"),]
r$Regen_perc <- r$Regen/r$Regen_StudySpp
m$Mature_perc <- m$Mature/m$Mature_StudySpp
dd <- merge(r,m,by=c("Plot","Species"))

dd <- dd %>% pivot_longer(cols = c("Regen_perc","Mature_perc"),values_to = "Percent",names_to = "Class")

########################FigureB.2
ggplot(dd,aes(x=Class,y = Percent*100,fill = Species))+
  geom_boxplot()+
  scale_fill_manual(values = c("#a63d40","#fb8500","#09814a","#6754A0","#C1A5A9","#25AED0"))+
  labs(y = "Species composition (%)",x=NULL)+
  scale_x_discrete(labels = c("Pre-fire mature","Post-fire regeneration"))+
  scale_y_continuous(expand = c(0.01,0))+
  theme(axis.text.x = element_text(size =10.4, color = "black"))
#ggsave("Figure_B-2.png",dpi = 600, width = 5.5, height = 3.5, units = "in")



#
#
#
####See 'ClarkWolf_et_al_2022_code.R' for Figures B.3, B.4, and B.5 ##################################################
####Figure B.6 #########################################################################
####Soil nitrogen data

sl <- merge(sl,plot[c("Plot","Severity")])%>%
  mutate(Status = ifelse(Severity == "Unburned","Unburned","Burned"))

s1 <- sl %>% dplyr::select(Plot,Severity,Status,ends_with("day")) %>%
  pivot_longer(cols = starts_with("soil"),values_to = "Resin_n",names_to = "Type",names_prefix = "soil_") %>%
  mutate(Type = factor(Type, levels = c("ug.N.day","NH4.N.day","NO3.N.day")))
levels(s1$Type)

g1 <-  ggplot(s1, aes(x = Type, y = Resin_n, fill = Status))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey50","grey80"))+
  scale_x_discrete(labels = c("Total","NH4","NO3"))+
  coord_cartesian(ylim = c(0,2))+
  labs(x= NULL, y = expression("Resin nitrogen ("~mu*"g N"~day^-1~")"),fill=NULL)+
  ggtitle("Aug 2018 - Jun 2019")+
  annotate("text", x = 1, y = 2, label = "***")+
  annotate("text", x = 2, y = 2, label = "***")+
  annotate("text", x = 3, y = 2, label = "***")+
  #coord_equal(ratio = 3, ylim = c(0,2))+
  theme(axis.text.x = element_text(size = 10, color = "black"))
g1                  

s2 <- sl %>% dplyr::select(Plot,Severity,Status,starts_with("soil_M_20")) %>%
  mutate(soil_M_2018_Total = soil_M_2018_NH4+soil_M_2018_NO3, soil_M_2019_Total = soil_M_2019_NH4+soil_M_2019_NO3,) %>%
  pivot_longer(cols = starts_with("soil"),values_to = "Soil_N",names_to = c("Year","Type"),names_prefix = "soil_M_",names_sep = "_") %>%
  mutate(Type = factor(Type, levels = c("Total","NH4","NO3")))

g2 <-s2 %>% filter(Year == "2018") %>%
  ggplot(aes(x = Type, y = Soil_N, fill = Status))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey50","grey80"))+
  scale_x_discrete(labels = c("Total","NH4","NO3"))+
  scale_y_continuous(limits = c(0,31))+
  #coord_cartesian(ylim = c(0,2))+
  labs(x= NULL, y = expression("Mineral soil inorganic nitrogen ("~mu*"g N"~g^-1~"soil )"),fill=NULL)+
  ggtitle("2018")+
  annotate("text", x = 1, y = 31, label = "***")+
  annotate("text", x = 2, y = 31, label = "***")+
  annotate("text", x = 3, y = 31, label = "***")+
  theme(axis.text.x = element_text(size = 10, color = "black"))
g2                

g3 <- s2 %>% filter(Year == "2019") %>%
  ggplot(aes(x = Type, y = Soil_N, fill = Status))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey50","grey80"))+
  scale_x_discrete(labels = c("Total","NH4","NO3"))+
  scale_y_continuous(limits = c(0,31))+
  labs(x= NULL, y = NULL,fill=NULL)+
  ggtitle("2019")+
  annotate("text", x = 1, y = 31, label = "*")+
  annotate("text", x = 2, y = 31, label = "*")+
  theme(axis.text.x = element_text(size = 10, color = "black"))
g3                

g4 <- ggplot(sl, aes(x=Status, y = soil_M_perc_C,fill=Status))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey50","grey80"))+
  labs(x= "", y = "Mineral soil carbon (%)",fill=NULL)+
  theme(axis.text.x = element_text(size = 10, color = "black"))
g4

g5 <- ggplot(sl, aes(x=Status, y = soil_M_perc_N,fill=Status))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey50","grey80"))+
  labs(x= "", y = "Mineral soil total nitrogen (%)",fill=NULL)+
  annotate("text", x = 1.5, y = 0.3, label = "*")+
  theme(axis.text.x = element_text(size = 10, color = "black"))
g5

(g1|(g4|g5)/(g2|g3))+plot_layout(guides = "collect",width = c(40,60)) & theme(legend.text = element_text(size = 11))
ggsave("Figure_B-6.jpg",dpi = 300, width = 8, height = 6.5, units = "in")






