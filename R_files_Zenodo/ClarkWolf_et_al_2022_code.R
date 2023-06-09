#ClarkWolf_et_al_2022_code.R
#
# This script will reproduce the results figures presented in 
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
      "jtools","gstat","ggstance","performance", "MuMIn", "MASS",
      "pROC","DHARMa", "car","stringr", "officer","rvg")

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

#microclimate data
mc <- read.csv("ClarkWolf_et_al_2022_Microclimate.csv")

#recruitment data
rc <- read.csv("ClarkWolf_et_al_2022_RecruitmentData.csv")

####Functions

#Crossvalidation function, for glmms
CV = function(model,response,data,k,z) #k is the number of folds and z the number of iterations
{
  acc.metrics<-data.frame(id = seq(1,z), RMSE = NA, r = NA)
  for(i in 1:z){
    tryCatch({
      temp = data
      temp$K = sample(seq(1,k),nrow(temp),replace=T)
      data.train<-temp[temp$K != 1,]
      data.valid<-temp[temp$K == 1,]
      temp.model<-update(model,data=data.train,control = glmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
      pred<-predict(temp.model,type = "response",newdata=data.valid, re.form=NA) 
      obs<-as.vector(subset(data.valid,select=response)[,1])
      rmse<-sqrt(sum((log(obs)-log(pred))^2)/length(data.valid[,1]))
      acc.metrics[acc.metrics$id == i, 2]<-rmse
      acc.metrics[acc.metrics$id == i, 3]<-cor.test(obs,pred,method = "spearman")$est[[1]]
    },error = function(e){cat("ERROR :",conditionMessage(e), "\\n")})
  }
  accuracy<<-as.numeric(acc.metrics$RMSE)
  correlation <<- as.numeric(acc.metrics$r)
  mean.acc<<-mean(accuracy,na.rm = T)
  mean.cor<<- mean(correlation,na.rm = T)
  median.cor<<- median(correlation,na.rm = T)
  print(mean.acc)
  print(mean.cor)
  print(median.cor)
  hist(correlation)
}

#Crossvalidation function, for glms
CV.glm = function(model,response,data,k,z) #k is the number of folds and z the number of iterations
{
  acc.metrics<-data.frame(id = seq(1,z), RMSE = NA, r = NA)
  for(i in 1:z){
    tryCatch({
      temp = data
      temp$K = sample(seq(1,k),nrow(temp),replace=T)
      data.train<-temp[temp$K != 1,]
      data.valid<-temp[temp$K == 1,]
      temp.model<-update(model,data=data.train)
      pred<-predict(temp.model,type = "response",newdata=data.valid) 
      obs<-as.vector(subset(data.valid,select=response)[,1])
      rmse<-sqrt(sum((obs-pred)^2)/length(data.valid[,1]))
      acc.metrics[acc.metrics$id == i, 2]<-rmse
      acc.metrics[acc.metrics$id == i, 3]<-cor.test(obs,pred,method = "spearman")$est[[1]]
    },error = function(e){cat("ERROR :",conditionMessage(e), "\\n")})
  }
  accuracy<<-as.numeric(acc.metrics$RMSE)
  correlation <<- as.numeric(acc.metrics$r)
  mean.acc<<-mean(accuracy,na.rm = T)
  mean.cor<<- mean(correlation,na.rm = T)
  median.cor<<- median(correlation,na.rm = T)
  print(mean.acc)
  print(mean.cor)
  print(median.cor)
  hist(correlation)
}

#Overdispersion function
overdisp_fun = function(model) {
  sum( residuals(model, type = "pearson")^2)/df.residual(model)
}

#Logit to probability function
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#
#
#
#
####Initial data setup ################################################

#############################Set up data

d <- plot

####Transform highly skewed predictors
transform = c("PP_dWt","DSS","ba_PICO", "BA_PICO_L","BA_LAOC_L","BA_PSME_L",
              "cover.BG_Avg","cover.Moss_Avg","cover.Shrub_Avg","cover.Grass_Avg","cover.Wood_Avg",
              "Canopy.green_Avg","NO3","NH4","soilN")

#visualize distributions
op <- par(mfrow = c(4,4),mar = c(4,4,1,2))
for (i in 1:length(transform)){
  predictor = transform[i]
  hist(d[,paste0(predictor)], main = paste0(predictor))
  hist(sqrt(d[,paste0(predictor)]), main = paste0(predictor))
}
par(op)

#Square root transform
d1 = d[,transform] %>%
  mutate_all(sqrt)
names(d1) = paste0(names(d1),"_sqrt")
d = cbind(d,d1)
d$logArea = log(d$TransectSize)

#Rename predictor variables for ease in typing later on
preds = c("PP_dWt_sqrt","DSS_sqrt","ba_PICO_sqrt","BA_PICO_L_sqrt","BA_PSME_L_sqrt","BA_LAOC_L_sqrt",
          "HLI","ppt_ann","ppt_JJAS","ppt_pf","DEF","Tmax_abs","Tmax_avg","Axis1","Axis2","dnbr",
          "Canopy.tot_Avg","Canopy.green_Avg_sqrt", "cover.BG_Avg_sqrt","cover.Forb_Avg","cover.Moss_Avg_sqrt",
          "cover.Shrub_Avg_sqrt","cover.Grass_Avg_sqrt","cover.Wood_Avg_sqrt","NO3_sqrt","NH4_sqrt","soilN_sqrt")
preds.names = c("PP","DSS","ba_PICO","Live_ba_PICO","Live_ba_PSME","Live_ba_LAOC",
                "HLI","ppt_ann","ppt_JJAS","ppt_pf","DEF","Tmax","Tmax_avg","Axis1","Axis2","dnbr",
                "Canopy","CanopyLive", "BG","Forb","Moss","Shrub","Grass","Wood","NO3","NH4","soilN")
name =  paste0("z.",preds.names)
d1 = d[,preds]; names(d1) = name
d = cbind(d,d1)

#Separate out data by severity classes for separate models
u = d[which(d$Severity == "Unburned"),]
b = d[which(d$Severity != "Unburned"),]
a = d

#Scale predictor variables (z-score)
a[,c(name)] = scale(a[,c(name)],center = T)
u[,c(name)] = scale(u[,c(name)],center = T)
b[,c(name)] = scale(b[,c(name)],center = T)

####Figure 3#########################################################
theme_set(theme_classic())

m <- plot %>% dplyr::select(Plot,Fire,Severity,Date_1,n_live_V1_c1:germ_2019_tot_c2)

#Calculate survivorship at each sampling interval
m$Surv_v1v3_c1 <- m$n_live_V3_c1/m$n_live_V1_c1
m[which(m$Surv_v1v3_c1 == Inf),"Surv_v1v3_c1"] <- NA
m$Surv_v2v3_c2 <- m$n_live_V3_c2/m$germ_2019_tot_c2
m$Surv_v3v4_c1 <- m$n_live_V4_c1/m$n_live_V3_c1
m$Surv_v3v4_c2 <- m$n_live_V4_c2/m$n_live_V3_c2

#Calculate survivorship by seedling age
A <- pivot_longer(m,cols = starts_with("Surv"),names_to = c("time","cohort"),names_sep = "_",
                  names_prefix = "Surv_",values_to = "Surv")

A$Age <- NA
A[which(A$time=="v1v3" & A$cohort=="c1"),"Age"] <- "1-2"
A[which(A$time=="v2v3" & A$cohort=="c2"),"Age"] <- "<1"
A[which(A$time=="v3v4" & A$cohort=="c2"),"Age"] <- "1-2"
A[which(A$time=="v3v4" & A$cohort=="c1"),"Age"] <- "2-3"

A$Year <- NA
A$Year <- factor(A$time,levels = c("v2v3","v1v3","v3v4"),labels = c("summer 2019","2018-2019","2019-2020"))

A <- A[-which(is.na(A$Surv)),]
A$Date_1 <- as.Date(A$Date_1,"%m/%d/%Y")
A <- A[-which(A$Date_1< "2018-07-15" & A$cohort == "c1" & A$Year == "2018-2019"),]

########Plot survivorship by seedling age transition
p1 <- ggplot(A,aes(x=as.factor(Age),y=Surv*100))+#,fill = Year))+
  geom_boxplot()+
  labs(x="Seedling age transition (yr)",y = "Survivorship (%)")+
  theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 13, color = "black"),
        axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 14))
p1

#######Pool across cohorts - total survivorship in each plot
m$n_tot <- m$n_live_V1_c1+m$n_live_V3_c2
m$n_live <- m$n_live_V4_c1+m$n_live_V4_c2
m$Surv <- m$n_live/m$n_tot

m$Burned <- ifelse(m$Severity=="Unburned","Unburned","Burned")

p2 <- ggplot(m,aes(x=Burned,y=Surv*100))+
  geom_boxplot()+
  labs(x="",y = "Survivorship (%)")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14, color = "black"),# angle = -45, hjust = 0, vjust = 1),
        axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 14))
p2

#######Plot regeneration by burned/unburned
d <- plot
d$Status = ifelse(d$Sev == "Unburned","Unburned","Burned")
d1 <- d %>% mutate(across(Regen_AllSpp_ct:Regen_Abies_ct,~(.x/TransectSize)*10000)) %>%#calculate seedling densities
  rename_with(~str_remove(., '_ct')) %>%
  dplyr::select(c("Status", starts_with("Regen_"),-ends_with("UNKN")))%>%
  pivot_longer(cols = starts_with("Regen_"),names_to = "Spp",values_to = "Regen",names_prefix = "Regen_")

d1$Spp <- factor(d1$Spp, levels = c("AllSpp","Abies","LAOC","PICO","PIEN","PIPO","PSME"),
                 labels = c("All","Abies","LAOC","PICO","PIEN","PIPO","PSME"))
p3 <- ggplot(d1,aes(x=Spp,y = Regen+1, fill = Status))+
  geom_boxplot()+
  labs(x = "Species",y = expression("Total regeneration  ( # "~ha^-1~")"))+
  scale_y_log10(expand = c(0,0.1),breaks = c(1,10,100,1000,10000,100000),
                labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  scale_fill_manual(values = c("grey50","grey80"))+
  annotate("text",x=4, y = 6*10^5, label = "***")+
  annotate("text",x=3, y = 6*10^5, label = "**")+
  annotate("text",x=6, y = 6*10^5, label = "**")+
  theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 14))
p3

#######Plot annual recruitment
r1 <- rc[rc$Spp=="AllSpp",] #select total recruitment (all species)
r1$Germ <- (r1$Germ/r1$TransectSize)*10000

#remove year 1 data from plots sampled early in the season 
r1 <- merge(r1,plot[,c("Plot","Severity","Date_1")],by="Plot"); r1$Date_1 <- as.Date(r1$Date_1,"%m/%d/%Y")
r1 <- r1[-which(r1$Postfire_year==1 & r1$Date_1<"2018-07-15"),]
r1$Postfire_year <- as.factor(r1$Postfire_year)
r1$Status <- ifelse(r1$Severity == "Unburned","Unburned","Burned")

p4 <- ggplot(r1, aes(y=Germ+1, fill = Postfire_year, x = Status))+
  geom_boxplot(alpha = 0.85)+
  scale_fill_manual(values = c( "#34a0a4","#1a759f","#184e77"))+
  labs(y = expression("Annual recruitment  ( # "~ha^-1~")"), x = NULL, fill = "Year")+
  scale_y_log10(expand = c(0,0.1),breaks = c(1,10,100,1000,10000,100000),
                labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  annotate("text",x=1, y = 6*10^5, label = "***",size = 4)+  
  annotate("text",x=2, y = 6*10^5, label = "*",size = 4)+
  theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 14))
p4

######Height
sd[which(sd$Spp %in% c("ABLA", "ABGR", "ABIES")), "Spp"] <- "Abies"
sd <- sd[-which(sd$Spp %in% c("UNKN","POTR","TSME")),]
sd$Spp <- factor(sd$Spp, levels = c("Abies","LAOC","PICO","PIEN","PIPO","PSME"))

p5 <- ggplot(sd,aes(x=as.factor(Age),y=Ht,fill = Spp))+
  geom_boxplot()+
  labs(y = "Height (cm)",fill = "Species",x = "Seedling age (yr)")+
  scale_fill_manual(values = c("#a63d40","#fb8500","#09814a","#6754A0","#C1A5A9","#25AED0"))+
  theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 13, color = "black"),
        axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 14))
p5

zz <- p3+p4+plot_layout(widths = c(13,7))
zy <- p2+p1+p5+plot_layout(widths = c(5,5,9))
myplot <- zz/zy+plot_layout(guides = "collect")+plot_annotation(tag_levels = "A") &
  theme(legend.box.background = element_rect(color = "black", fill = NULL),
        legend.box.margin = margin(2,1,2,2,unit = "pt"), legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))
myplot
#ggsave("Figure_3.png",dpi = 300, width =10, height = 7.5, units = "in")

#Can fix the formatting by exporting an editable figure into powerpoint
editable_graph <- dml(ggobj = myplot)
doc <- read_pptx("newplot.pptx")
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_fullsize() )
print(doc, target = "newplot.pptx")

####Total regeneration statistical model (Figures B.3 - B.4) ##########################################################################################

###########################################Total regeneration in burned plots 
theme_set(theme_classic(base_size = 8))

m <-glmer.nb(Regen_AllSpp_ct~offset(logArea)+z.DSS+z.ba_PICO+z.Tmax+
               z.Axis1+z.Axis2+z.soilN+ (1|Fire),data=b)

#Model diagnostics
summary(m)
vif(m)
CV(m,"Regen_AllSpp_ct",b,4,1000) #rmse (log) = 1.338, rho = 0.632, median 0.665
cor.test(predict(m),log(b$Regen_AllSpp_ct), method = "spearman")
r.squaredGLMM(m)

#Check residuals
check_model(m)
hist(resid(m))
op <- par(mfrow=c(1,2),mar=c(4,4,1,2))
scatter.smooth(predict(m),resid(m))
abline(a=0,b=0,col="red")
plot(x=predict(m),y=log(b$Regen_AllSpp_ct))
l = lm(log(b$Regen_AllSpp_ct)~predict(m))
abline(a=l$coefficients[[1]], b= l$coefficients[[2]])
abline(a=0,b=1,col="red") 

#Check spatial correlation
points = plot[which(plot$Plot %in% b$Plot),]
coordinates(points) <- c("LAT_WGS84" ,"LON_WGS84")
crs(points) <- "+proj=longlat +datum=WGS84 +no_defs"

E = resid(m,type = "pearson")
newdata = data.frame(E, points$LON_WGS84, points$LAT_WGS84)
coordinates(newdata) <- ~points.LON_WGS84 + points.LAT_WGS84
bubble(newdata, "E", col = c("red", "blue"), main = "Residuals", xlab = "longitude", ylab = "latitude",alpha=0.5)

Vario1 = variogram(E ~ 1, newdata)  # Creates a semivariogram to examine the autocorrelation of the residuals
plot(Vario1)

vario2 = variogram(E ~ 1, newdata, alpha = c(0,45,90,135))  # Splits the semivariogram to check for autocorrelation in all directions
plot(vario2)

#
#
#
###############################################Plotting: Figure B.3
off <- mean(b$TransectSize) #Average transect size to use as offset

eff <- effect("z.ba_PICO",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g1 <- ggplot(b,aes(x  = z.ba_PICO, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3)+
  geom_ribbon(data = effect, aes(x = z.ba_PICO, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.ba_PICO,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 0.5)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(limits = c(-1,3))+
  labs(x = "Basal area PICO (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g1

eff <- effect("z.DSS",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g2 <- ggplot(b,aes(x  = z.DSS, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3)+
  geom_ribbon(data = effect, aes(x = z.DSS, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.DSS,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 0.5)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  labs(x = "Distance to seed source (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g2

eff <- effect("z.Tmax",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g3 <- ggplot(b,aes(x  = z.Tmax, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3)+
  geom_ribbon(data = effect, aes(x = z.Tmax, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.Tmax,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 0.5)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  labs(x = expression("Microclimate ("*italic(T[max])*" z-score)"),y = expression("Predicted regeneration ( # "~ha^-1~")"))
g3

eff <- effect("z.Axis1",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g4 <- ggplot(b,aes(x  = z.Axis1, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3)+
  geom_ribbon(data = effect, aes(x = z.Axis1, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.Axis1,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 0.5)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  labs(x =expression("Overstory fire severity ("*italic("Axis1")*" z-score)"),y = expression("Predicted regeneration ( # "~ha^-1~")"))
g4

eff <- effect("z.Axis2",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g5 <- ggplot(b,aes(x  = z.Axis2, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3)+
  geom_ribbon(data = effect, aes(x = z.Axis2, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.Axis2,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 0.5)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_reverse()+
  labs(x = expression("Moss + forb cover ("*italic("Axis2")*" z-score)"),
       y = expression("Predicted regeneration ( # "~ha^-1~")"))
g5

eff <- effect("z.soilN",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g6 <- ggplot(b,aes(x  = z.soilN, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3)+
  geom_ribbon(data = effect, aes(x = z.soilN, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.soilN,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 0.5)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  labs(x = "Soil N availability (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g6

p1 = ggplot(data.frame( x = 1, y = 1)) +
  annotate("text", label = expression("Expected regeneration ( # "~ha^-1~")"), angle = 90, y = 1, x = 1,size = 3.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

g1$labels$y <-g2$labels$y<-g3$labels$y<-g5$labels$y<-g6$labels$y<-g4$labels$y <- NULL

p1+(g1+g2+g3)/(g4+g5+g6)+ plot_layout(widths = c(1,25))+ plot_annotation(tag_levels = list(c("","A","B","C","D","E","F"))) & theme(
  plot.tag.position = c(0,1)
)
#ggsave("Figure_B-3.png",dpi = 600, width = 6.5, height = 4, units = "in")

#
#
#
#################################################Total regeneration in unburned plots

m =glm.nb(Regen_AllSpp_ct~offset(logArea) 
          +z.DSS #negative - fewer seedlings in sparse, low productivity forests
          +z.CanopyLive #negative - light limitation
          +z.Moss #positive - more seedlings in sites with better seedbeds
          +z.soilN #positive - more seedlings in sites with more nutrients (or closer to burned areas)
          ,data=u)
summary(m)

off <- mean(u$TransectSize) #Average transect size to use as offset

eff <- effect("z.DSS",m,se=T,offset = log(off),xlevels = 100)
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g1 <- ggplot(u,aes(x  = z.DSS, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.DSS, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.DSS,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 1)+
  scale_y_log10(breaks = c(10,100,1000,10000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4)))+
  coord_cartesian(ylim = c(100,30000))+
  labs(x = "DSS (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g1

eff <- effect("z.CanopyLive",m,se=T,offset = log(off),xlevels = 100)
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g2 <- ggplot(u,aes(x  = z.CanopyLive, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.CanopyLive, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.CanopyLive,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 1)+
  scale_y_log10(breaks = c(10,100,1000,10000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4)))+
  coord_cartesian(ylim = c(100,30000))+
  labs(x = "Live canopy cover (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g2

eff <- effect("z.Moss",m,se=T,offset = log(off),xlevels = 100)
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g3 <- ggplot(u,aes(x  = z.Moss, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.Moss, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.Moss,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 1)+
  scale_y_log10(breaks = c(10,100,1000,10000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4)))+
  coord_cartesian(ylim = c(100,30000))+
  labs(x = "Moss cover (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g3

eff <- effect("z.soilN",m,se=T,offset = log(off),xlevels = 100)
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g4 <- ggplot(u,aes(x  = z.soilN, y = (Regen_AllSpp_ct/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.soilN, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = z.soilN,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 1)+
  scale_y_log10(breaks = c(10,100,1000,10000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4)))+
  coord_cartesian(ylim = c(100,30000))+
  labs(x = "Soil N availability (z-score)",y = expression("Predicted regeneration ( # "~ha^-1~")"))
g4

p1 = ggplot(data.frame( x = 1, y = 1)) +
  annotate("text", label = expression("Expected regeneration ( # "~ha^-1~")"), angle = 90, y = 1, x = 1,size = 4) + 
  theme_void() +
  coord_cartesian(clip = "off")

g1$labels$y <-g2$labels$y<-g3$labels$y<-g5$labels$y<-g6$labels$y<-g4$labels$y <- " "

p1+(g1+g2)/(g3+g4)+ plot_layout(widths = c(1,25))
#ggsave("Figure_B-4.jpg",dpi = 600, width = 5.5, height = 4.5, units = "in")

#
#
#
####Hypothesis testing #################################################################################
#########Test effects of seed availability, climate, & fire severity variables on seedling regeneration in burned plots

#final (reduced) model for total seedling regeneration 3 yrs postfire
m = glmer.nb(Regen_AllSpp_ct ~ offset(logArea)+z.ba_PICO+z.DSS +z.Tmax_2019+z.Axis1 +z.Axis2  +z.soilN+(1|Fire),data = b)

#null model
null = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) +(1|Fire),data = b, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#seed availability
seed = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) +z.ba_PICO +z.DSS +(1|Fire),data = b, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(null,seed)
AIC(seed);BIC(seed)
cor.test(predict(seed),log(b$Regen_AllSpp_ct),method = "spearman")#0.419
CV(seed,"Regen_AllSpp_ct",b,4,500) #2567, 0.35, 0.37

#seed + macroclimate
macro = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) +z.ba_PICO +z.DSS+z.HLI +z.ppt_JJAS+z.DEF+(1|Fire),data = b)
anova(seed,macro)
AIC(macro);BIC(macro)
cor.test(predict(macro),log(b$Regen_AllSpp_ct),method = "spearman")#0.615
CV(macro,"Regen_AllSpp_ct",b,4,500) #1470, 0.491, 0.516

#seed+microclimate
micro = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) +z.ba_PICO +z.DSS +z.HLI+z.ppt_JJAS+z.Tmax +(1|Fire),data = b)
anova(seed,micro)
AIC(micro);BIC(micro)
cor.test(predict(micro),log(b$Regen_AllSpp_ct),method = "spearman")#0.627
CV(micro,"Regen_AllSpp_ct",b,4,500) #814, 0.497, 0.518

#seed+macroclimate+fire
fire = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) +z.ba_PICO  +z.DSS +z.HLI+z.ppt_JJAS  +z.DEF+z.Axis1 +z.Axis2+z.soilN
                +(1|Fire),data = b, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(macro,fire)
AIC(fire);BIC(fire)
cor.test(predict(fire),log(b$Regen_AllSpp_ct),method = "spearman")#0.837
CV(fire,"Regen_AllSpp_ct",b,4,1000) #614, 0.6048, 0.646

#seed+microclimate+fire
mcfire = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) +z.ba_PICO +z.DSS +z.HLI+z.ppt_JJAS+z.Tmax+z.Axis1+z.Axis2
                  +z.soilN+(1|Fire),data = b, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(mcfire);isSingular(mcfire)
anova(micro,mcfire)
AIC(mcfire);BIC(mcfire)
cor.test(predict(mcfire),log(b$Regen_AllSpp_ct),method = "spearman")#0.8618
CV(mcfire,"Regen_AllSpp_ct",b,4,500) #513, 0.6097, 0.648

#########summarize cross-validation results for all models (to use for plotting)
All <- data.frame(Species = rep("All",6000),Model = c(rep("Seed",1000),rep("Climate",1000),rep("MicroClim",1000),
                                                      rep("Fire",1000),rep("MC_fire",1000),rep("Reduced",1000)),num = rep(seq(1,1000),6),rho = NA)
models = list(seed,macro,micro,fire,mcfire,m)

#iterate crossvalidation and store for each subset model - this will hake a while!
for (i in 1:6){
  model = models[[i]]
  modelName = unique(All$Model)[i]
  response = "Regen_AllSpp_ct"
  for(j in 1:1000){
    tryCatch({
      temp = b
      temp$K = sample(seq(1,4),nrow(temp),replace=T)
      data.train<-temp[temp$K != 1,]
      data.valid<-temp[temp$K == 1,]
      temp.model<-update(model,data=data.train,control = glmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
      pred<-predict(temp.model,type = "response",newdata=data.valid, re.form=NA) 
      obs<-as.vector(subset(data.valid,select=response)[,1])
      All[which(All$Model == modelName & All$num ==j), "rho"]<-cor.test(obs,pred,method = "spearman")$est[[1]]
    },error = function(e){cat("ERROR :",conditionMessage(e), "\\n")})
  }
}

All_rho = data.frame(Models =  c("Seed","Climate","MicroClim","Fire","MC_fire","Reduced"), rho= NA)

for (i in 1:6){
  model = models[[i]]
  All_rho[i,"rho"] <- cor.test(predict(model,type = "response"),d1$Regen_AllSpp_ct,method = "spearman")$est[[1]]
}

###########################################Repeat for individual species models
##################################PSME
b$PSME_present = as.factor(b$PSME_present)

m = glmer.nb(Regen_PSME_ct ~ offset(logArea) 
             +PSME_present 
             +z.Live_ba_PSME
             +z.Tmax
             +z.dnbr
             +poly(z.Axis2,2,raw=T)
             +(1|Fire)
             ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

seed = glmer.nb(Regen_PSME_ct ~ offset(logArea) 
                +PSME_present 
                +z.Live_ba_PSME
                +(1|Fire)
                ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

macro = glmer.nb(Regen_PSME_ct ~ offset(logArea) 
                 +PSME_present 
                 +z.Live_ba_PSME 
                 +z.HLI
                 +z.ppt_pf
                 +z.DEF
                 +(1|Fire)
                 ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

micro = glmer.nb(Regen_PSME_ct ~ offset(logArea) 
                 +PSME_present 
                 +z.Live_ba_PSME 
                 +z.HLI
                 +z.ppt_pf
                 +z.Tmax
                 +(1|Fire)
                 ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

fire = glmer.nb(Regen_PSME_ct ~ offset(logArea) 
                +PSME_present 
                +z.Live_ba_PSME 
                +z.HLI
                +z.ppt_pf
                +z.DEF
                +z.dnbr
                +z.Axis2
                +z.soilN
                +(1|Fire)
                ,data = b, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

mcfire = glmer.nb(Regen_PSME_ct ~ offset(logArea) 
                  +PSME_present 
                  +z.Live_ba_PSME 
                  +z.HLI
                  +z.ppt_pf
                  +z.Tmax
                  +z.dnbr
                  +z.Axis2
                  +z.soilN
                  +(1|Fire)
                  ,data = b, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

PSME <- data.frame(Species = rep("PSME",6000),Model = c(rep("Seed",1000),rep("Climate",1000),rep("MicroClim",1000),rep("Fire",1000),
                                                        rep("MC_fire",1000),rep("Reduced",1000)),num = rep(seq(1,1000),6),rho = NA)
models = list(seed,macro,micro,fire,mcfire,m)

for (i in 1:6){
  model = models[[i]]
  modelName = unique(PSME$Model)[i]
  response = "Regen_PSME_ct"
  for(j in 1:1000){
    tryCatch({
      temp = b
      temp$K = sample(seq(1,4),nrow(temp),replace=T)
      data.train<-temp[temp$K != 1,]
      data.valid<-temp[temp$K == 1,]
      temp.model<-update(model,data=data.train,control = glmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
      pred<-predict(temp.model,type = "response",newdata=data.valid, re.form=NA) 
      obs<-as.vector(subset(data.valid,select=response)[,1])
      PSME[which(PSME$Model == modelName & PSME$num ==j), "rho"]<-cor.test(obs,pred,method = "spearman")$est[[1]]
    },error = function(e){cat("ERROR :",conditionMessage(e), "\\n")})
  }
}

PSME_rho = data.frame(Models =  c("Seed","Climate","MicroClim","Fire","MC_fire","Reduced"), rho= NA)

for (i in 1:6){
  model = models[[i]]
  PSME_rho[i,"rho"] <- cor.test(predict(model,type = "response"),b$Regen_PSME_ct,method = "spearman")$est[[1]]
}


####################################PICO
b$PICO_present=as.factor(b$PICO_present)

m = glm.nb(Regen_PICO_ct ~ offset(logArea) 
           +PICO_present
           +z.PP
           +z.ba_PICO
           +z.HLI
           +z.DEF
           +poly(z.Canopy,2,raw=T)
           +z.Wood
           +z.soilN
           +z.Wood:z.DEF
           ,data = b,control = glm.control(maxit =100, trace = F))

seed = glm.nb(Regen_PICO_ct ~ offset(logArea) 
              +PICO_present
              +z.PP
              +z.ba_PICO
              ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

macro = glm.nb(Regen_PICO_ct ~ offset(logArea) 
               +PICO_present
               +z.ba_PICO
               +z.PP
               +z.HLI
               +z.ppt_JJAS
               +z.DEF
               ,data = b,control = glm.control(maxit =100, trace = F)) 

micro = glm.nb(Regen_PICO_ct ~ offset(logArea) 
               +PICO_present
               +z.ba_PICO
               +z.PP
               +z.HLI
               +z.ppt_JJAS
               +z.Tmax
               ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

fire = glm.nb(Regen_PICO_ct ~ offset(logArea) 
              +PICO_present
              +z.PP
              +z.ba_PICO
              +z.HLI
              +z.ppt_JJAS
              +z.DEF
              +z.Axis1
              +z.Wood
              +z.soilN
              ,data = b, control =glm.control(maxit = 100,trace=F))

mcfire = glm.nb(Regen_PICO_ct ~ offset(logArea) 
                +PICO_present 
                +z.PP
                +z.ba_PICO
                +z.HLI
                +z.ppt_JJAS
                +z.Tmax
                +z.Axis1
                +z.Wood
                +z.soilN
                ,data = b, control =glm.control(maxit = 100,trace=F))

PICO <- data.frame(Species = rep("PICO",6000),Model = c(rep("Seed",1000),rep("Climate",1000),rep("MicroClim",1000),rep("Fire",1000),
                                                        rep("MC_fire",1000),rep("Reduced",1000)),num = rep(seq(1,1000),6),rho = NA)
models = list(seed,macro,micro,fire,mcfire,m)

for (i in 1:6){
  model = models[[i]]
  modelName = unique(PICO$Model)[i]
  response = "Regen_PICO_ct"
  for(j in 1:1000){
    tryCatch({
      temp = b
      temp$K = sample(seq(1,4),nrow(temp),replace=T)
      data.train<-temp[temp$K != 1,]
      data.valid<-temp[temp$K == 1,]
      temp.model<-update(model,data=data.train)
      pred<-predict(temp.model,type = "response",newdata=data.valid) 
      obs<-as.vector(subset(data.valid,select=response)[,1])
      PICO[which(PICO$Model == modelName & PICO$num ==j), "rho"]<-cor.test(obs,pred,method = "spearman")$est[[1]]
    },error = function(e){cat("ERROR :",conditionMessage(e), "\\n")})
  }
}

PICO_rho = data.frame(Models =  c("Seed","Climate","MicroClim","Fire","MC_fire","Reduced"), rho= NA)

for (i in 1:6){
  model = models[[i]]
  PICO_rho[i,"rho"] <- cor.test(predict(model,type = "response"),b$Regen_PICO_ct,method = "spearman")$est[[1]]
}

#####################################LAOC
b$LAOC_present=as.factor(b$LAOC_present)

m = glm.nb(Regen_LAOC_ct ~ offset(logArea) 
           +LAOC_present
           +z.PP
           +z.Tmax
           +z.Axis2 
           +z.NO3
           ,data = b,control = glm.control(maxit =100, trace = F))

seed = glm.nb(Regen_LAOC_ct ~ offset(logArea) 
              +LAOC_present
              +z.PP
              ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

macro = glm.nb(Regen_LAOC_ct ~ offset(logArea) 
               +LAOC_present
               +z.PP 
               +z.HLI
               +z.ppt_ann
               +z.DEF
               ,data = b,control = glm.control(maxit =200, trace = FALSE)) #didn't fully converge

micro = glm.nb(Regen_LAOC_ct ~ offset(logArea) 
               +LAOC_present 
               +z.PP 
               +z.HLI
               +z.ppt_ann
               +z.Tmax
               ,data = b)#,control = glm.control(maxit =100, trace = FALSE)) 

fire = glm.nb(Regen_LAOC_ct ~ offset(logArea) 
              +LAOC_present
              +z.PP
              +z.HLI
              +z.ppt_ann
              +z.DEF
              +z.Axis1
              +z.Axis2
              +z.NO3
              ,data = b)#, control =glm.control(maxit = 200,trace=FALSE))

mcfire = glm.nb(Regen_LAOC_ct ~ offset(logArea) 
                +LAOC_present 
                +z.PP
                +z.HLI
                +z.ppt_ann
                +z.Tmax
                +z.Axis1
                +z.Axis2
                +z.NO3
                ,data = b)#, control =glm.control(maxit = 200,trace=F))

LAOC <- data.frame(Species = rep("LAOC",6000),Model = c(rep("Seed",1000),rep("Climate",1000),rep("MicroClim",1000),rep("Fire",1000),
                                                        rep("MC_fire",1000),rep("Reduced",1000)),num = rep(seq(1,1000),6),rho = NA)
models = list(seed,macro,micro,fire,mcfire,m)

for (i in 1:6){
  model = models[[i]]
  modelName = unique(LAOC$Model)[i]
  response = "Regen_LAOC_ct"
  for(j in 1:1000){
    tryCatch({
      temp = b
      temp$K = sample(seq(1,4),nrow(temp),replace=T)
      data.train<-temp[temp$K != 1,]
      data.valid<-temp[temp$K == 1,]
      temp.model<-update(model,data=data.train,evaluate = TRUE)
      pred<-predict(temp.model,type = "response",newdata=data.valid) 
      obs<-as.vector(subset(data.valid,select=response)[,1])
      LAOC[which(LAOC$Model == modelName & LAOC$num == j), "rho"]<-cor.test(obs,pred,method = "spearman")$est[[1]]
    },error = function(e){cat("ERROR :",conditionMessage(e), "\\n")})
  }
}

LAOC_rho = data.frame(Models =  c("Seed","Climate","MicroClim","Fire","MC_fire","Reduced"), rho= NA)

for (i in 1:6){
  model = models[[i]]
  LAOC_rho[i,"rho"] <- cor.test(predict(model,type = "response"),b$Regen_LAOC_ct,method = "spearman")$est[[1]]
}


##########################Combine all of the crossvalidation results

#Combine all of the crossvalidation results
CV_comb <- rbind(All,PSME,PICO,LAOC)

#Calculate standard errors for the average crossvalidated correlation coefficient for each subset model
cc <- CV_comb %>% group_by(Species, Model)%>%
  summarize_at("rho",c(mean,median,sd),na.rm=T)
names(cc) <- c("Species","Model","rho","median","sd")

cc$Species <- factor(cc$Species,levels = c("All","PSME","LAOC","PICO"))
cc$Model <- factor(cc$Model,levels = c("Seed","Climate","MicroClim","Fire","MC_fire","Reduced"),labels = c("Seed","Seed+Clim","Seed+MC","Seed+Clim+Fire","Seed+MC+Fire","Reduced"))
cc$se2 <- 2*(cc$sd/sqrt(1000))
cc$low <- cc$rho - cc$se2
cc$hi <- cc$rho + cc$se2

ggplot(cc,aes(x=Model,y=rho,color = Species))+
  geom_point(position = position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = low, ymax = hi),position = position_dodge(width=0.5))+
  labs(y = expression("Crossvalidated spearman's "~rho~" (predicted vs. observed)"))

#Combine all of the non-crossvalidated correlation coefficients for nonsubset models 
All_rho$Species <- "All"
PSME_rho$Species <- "PSME"
LAOC_rho$Species <- "LAOC"
PICO_rho$Species <- "PICO"
rho_comb <- rbind(All_rho,PSME_rho,PICO_rho,LAOC_rho)
names(cc)[3] <- "rho_CV"
rho_comb$Model <- factor(rho_comb$Model,levels = c("Seed","Climate","MicroClim","Fire","MC_fire","Reduced"),labels = c("Seed","Seed+Clim","Seed+MC","Seed+Clim+Fire","Seed+MC+Fire","Reduced"))

#Merge crossvalidated and non-crosvalidated results
rho <- merge(cc,rho_comb,by=c("Species","Model"))
#rho <- read.csv("L:/1_projectsData/NoRockies_GRIN_seedlings/Analysis/Field_data/Models_rho_CV.csv")

rho$Species <- factor(rho$Species,levels = c("All","PSME","LAOC","PICO"))
rho$Model <- factor(rho$Model,levels = c("Seed","Seed+Clim","Seed+MC","Seed+Clim+Fire","Seed+MC+Fire","Reduced"),
                    labels = c("Seed","Seed+Clim","Seed+MC","Seed+Clim+Fire","Seed+MC+Fire","Best-fit"))


#
#
#
####Figure 4 ############################################################################################
###############################Plots of coefficient values and crosvalidation results
theme_set(theme_classic())

#Need data frame "b" with scaled predictor variables in burned plots (see Initial Data Setup section above)

#######Fit reduced (best-fit) models for total regeneration density and individual species models
AllSpp = glmer.nb(Regen_AllSpp_ct ~ offset(logArea) 
                  +z.ba_PICO 
                  +z.DSS 
                  +z.Tmax
                  +z.Axis1 
                  +z.Axis2 
                  +z.soilN
                  +(1|Fire)
                  ,data = b)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

PSME = glmer.nb(Regen_PSME_ct~offset(logArea)
                +PSME_present
                +z.Live_ba_PSME
                +z.Tmax
                +z.dnbr
                #+z.Axis2
                +poly(z.Axis2,2,raw=TRUE)
                +(1|Fire), data=b)#, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

PICO = glm.nb(Regen_PICO_ct~offset(logArea)
              +PICO_present
              +z.HLI
              +z.ba_PICO
              +z.PP
              +z.DEF
              #+z.Canopy
              +poly(z.Canopy,2,raw=TRUE)
              +z.Wood
              +z.soilN
              +z.DEF:z.Wood
              ,data=b, control = glm.control(maxit=100,trace = F))

LAOC = glm.nb(Regen_LAOC_ct~offset(logArea)
              +LAOC_present
              +z.PP
              +z.Tmax
              +z.Axis1
              +z.Axis2
              +z.NO3
              , data=b, control = glm.control(maxit=100,trace = F))

######Extract coefficients and confidence intervals from each model
s <- summ(AllSpp,confint=TRUE)
s1 <- summ(PSME,confint=TRUE)
s2 <- summ(PICO,confint=TRUE)
s3 <- summ(LAOC,confint=TRUE)

coefs <- data.frame(Model = "All",Variable = rownames(s$coeftable),Est = s$coeftable[,1],Low = s$coeftable[,2],High = s$coeftable[,3])
coefs$Variable = c("Intercept","BA_PICO","DSS","Tmax","Axis1","Axis2","soilN")
coefs1 <- data.frame(Model = "PSME",Variable = rownames(s1$coeftable),Est = s1$coeftable[,1],Low = s1$coeftable[,2],High = s1$coeftable[,3])
coefs1$Variable = c("Intercept","Spp_present","BA_PSME","Tmax","dnbr","Axis2","Axis2^2")
coefs2 <- data.frame(Model = "PICO",Variable = rownames(s2$coeftable),Est = s2$coeftable[,1],Low = s2$coeftable[,2],High = s2$coeftable[,3])
coefs2$Variable = c("Intercept","Spp_present","HLI","BA_PICO","PP","DEF","Canopy","Canopy^2","Wood","soilN","DEF:Wood")
coefs3 <- data.frame(Model = "LAOC",Variable = rownames(s3$coeftable),Est = s3$coeftable[,1],Low = s3$coeftable[,2],High = s3$coeftable[,3])
coefs3$Variable = c("Intercept","Spp_present","PP","Tmax","Axis1","Axis2","NO3")

#####Reorder and name coefficients for plotting
coef <- rbind(coefs,coefs1,coefs2,coefs3)
coef$Variable <- factor(coef$Variable, levels = c("NO3","soilN","DEF:Wood","Wood","Axis2^2","Axis2","Canopy^2","Canopy","dnbr","Axis1","HLI","DEF","Tmax","PP","DSS","BA_PSME","BA_PICO","Spp_present","Intercept"))
levels(coef$Variable)
coef$VariableType <- ifelse(coef$Variable %in% c("NO3","soilN"),"Soil",ifelse(coef$Variable %in% c("DEF:Wood","Wood","Axis2^2","Axis2"),"Understory",
                                                                              ifelse(coef$Variable %in% c("Canopy^2","Canopy","dnbr","Axis1"),"Fire severity",
                                                                                     ifelse(coef$Variable %in% c("HLI","DEF","Tmax"),"Climate",
                                                                                            ifelse(coef$Variable %in% c("PP","DSS","BA_PSME","BA_PICO"),"Seed avail.","Species")))))
coef$VariableType <- factor(coef$VariableType, levels = c("Species","Seed avail.","Climate","Fire severity","Understory","Soil"))
coef$Model <- factor(coef$Model,levels = c("PSME","LAOC","PICO","All"))

#####################Plotting: Figure 4A
c1 <- coef[-which(coef$Variable %in% c("Intercept","Spp_present","Canopy^2","Axis2^2","DEF:Wood")),]
g1 <- ggplot(c1,aes(y=Variable, x = Est, color = Model))+
  geom_point(aes(shape = Model),size = 2,position=position_dodge(width = 0.6))+
  geom_errorbarh(aes(xmin = Low, xmax = High),position=position_dodge(width = 0.6),height = 0)+
  geom_vline(xintercept = 0,linetype = 2)+
  scale_color_manual(values = c("#219ebc","#fb8500","#09814a","grey10"))+#,"#a63d40"))+
  scale_shape_manual(values = c(17,15,18,19))+
  labs(color = "Species",shape = "Species",x = "Estimate",y = NULL)+
  scale_x_continuous(limits = c(-2,2),breaks = c(-2,-1,0,1,2,3,4),minor_breaks = c(-2,-1,0,1,2,3,4),expand = c(0,0))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 8.5),legend.key.height = unit(0.4,"cm"),
        legend.spacing.y = unit(0.05,"cm"), panel.grid = element_blank())#,legend.box.background = element_rect(color = "grey20"))#,legend.box.margin = margin(1,0,0,1,unit = "pt"))
g1 
g1 <- g1+theme(legend.position = "none")


####################Plotting:Figure 4B - 4C (need results of hypothesis testing crossvalidation, see above section)
theme_set(theme_bw())
g2 <-  ggplot(rho, aes(x=Model,y=rho,color = Species, shape = Species))+
  geom_point(size = 3)+
  ggtitle(label = expression(rho))+
  scale_color_manual(values = c("grey10","#219ebc","#fb8500","#09814a"))+#,"#a63d40"))+
  scale_shape_manual(values = c(19,17,15,18))+
  scale_y_continuous(limits = c(0.3, 1),expand = c(0,0))+
  labs(y=expression("Spearman's "~rho~" (pred. vs. obs. density)"), x = "Model")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.text = element_text(size = 11),legend.title = element_text(size=12))
g2

g3 <-  ggplot(rho, aes(x=Model,y=rho_CV,color = Species,shape = Species))+
  geom_point(size = 3)+
  ggtitle(label = expression(rho~"-CV"))+
  scale_color_manual(values = c("grey10","#219ebc","#fb8500","#09814a"))+#,"#a63d40"))+
  scale_shape_manual(values = c(19,17,15,18))+
  scale_y_continuous(limits = c(0.3, 1),expand = c(0,0))+
  labs(y=expression("Spearman's "~rho~" (pred. vs. obs. density)"), x = "Model")+
  theme(axis.text.x = element_text(angle = -45,hjust=.1,size =11),axis.title.x=element_text(size=12),legend.text = element_text(size = 11),legend.title = element_text(size=12))
g3

g1|g2/g3+plot_layout(guides="collect")
#ggsave("Figure4.png",dpi = 600, height = 7.5, width = 7.5, units = "in")

#
#
#
####Figure 5 ################################################
#Model of annual seedling recruitment


###################################################Set up data
#need scaled plot data from burned plots (see "Initial Data Setup" above)

#Recruitment data
r <- rc[rc$Spp=="AllSpp",] #select total recruitment (all species)

#Merge recruitment data with scaled plot-level data
preds <- c("z.PP","z.DSS","z.ba_PICO","z.Live_ba_PICO","z.Live_ba_PSME","z.Live_ba_LAOC",
           "z.HLI","z.ppt_ann","z.ppt_JJAS","z.ppt_pf","z.DEF","z.Axis1","z.Axis2","z.dnbr",
           "z.NO3","z.NH4","z.soilN")
r <- merge(r,b[,c("Plot","Severity",paste(preds))])

#####Format annual cover data
sp.avg <- sp %>% dplyr::select(-Subplot) %>% group_by(Plot, Postfire_year)%>% #average across subplots within each plot
  summarize_all(function(x){mean(x,na.rm=T)} )
sp.avg <- sp.avg[,c(1:9,13:15)]
sp1 <- merge(sp.avg,plot[,c("Plot","Severity","cover.BG_Avg","cover.Wood_Avg","cover.Litter_Avg",
                               "cover.Moss_Avg","cover.Grass_Avg","cover.Forb_Avg","cover.Shrub_Avg",
                               "Canopy.tot_Avg","Canopy.green_Avg","Canopy.dead_Avg")])

#Use averages to fill in missing annual values
sp1[which(is.na(sp1$cover.BG)),"cover.BG"] = sp1[which(is.na(sp1$cover.BG)),"cover.BG_Avg"] 
sp1[which(is.na(sp1$cover.Wood)),"cover.Wood"] = sp1[which(is.na(sp1$cover.Wood)),"cover.Wood_Avg"] 
sp1[which(is.na(sp1$cover.Litter)),"cover.Litter"] = sp1[which(is.na(sp1$cover.Litter)),"cover.Litter_Avg"] 
sp1[which(is.na(sp1$cover.Grass)),"cover.Grass"] = sp1[which(is.na(sp1$cover.Grass)),"cover.Grass_Avg"] 
sp1[which(is.na(sp1$cover.Forb)),"cover.Forb"] = sp1[which(is.na(sp1$cover.Forb)),"cover.Forb_Avg"] 
sp1[which(is.na(sp1$cover.Shrub)),"cover.Shrub"] = sp1[which(is.na(sp1$cover.Shrub)),"cover.Shrub_Avg"] 
sp1[which(is.na(sp1$cover.Moss)),"cover.Moss"] = sp1[which(is.na(sp1$cover.Moss)),"cover.Moss_Avg"] 
sp1[which(is.na(sp1$Canopy.tot)),"Canopy.tot"] = sp1[which(is.na(sp1$Canopy.tot)),"Canopy.tot_Avg"] 
sp1[which(is.na(sp1$Canopy.green)),"Canopy.green"] = sp1[which(is.na(sp1$Canopy.green)),"Canopy.green_Avg"] 
sp1[which(is.na(sp1$Canopy.dead)),"Canopy.dead"] = sp1[which(is.na(sp1$Canopy.dead)),"Canopy.dead_Avg"] 

#Transform and scale
sp1 <- sp1[-which(sp1$Severity=="Unburned"),c(1:12)]
sp.z <- sp1 %>% pivot_wider(names_from = "Postfire_year",values_from = c("cover.BG","cover.Wood","cover.Litter",
                                                                          "cover.Moss","cover.Grass","cover.Forb","cover.Shrub",
                                                                          "Canopy.tot","Canopy.green","Canopy.dead") ) %>%
  mutate(across(c(2:31),~as.numeric(scale(sqrt(.x)))))%>%
  pivot_longer(cols = c(2:31),values_to = "Cover",names_to = c("CoverType","Postfire_year"),
               names_sep = "_") %>%
  pivot_wider(names_from = "CoverType",values_from = "Cover",names_prefix = "z.")

#Combine with recruitment data
r <- merge(r,as.data.frame(sp.z),by=c("Plot","Postfire_year"))

#####Format microclimate data
mc$Postfire_year <- mc$Year-2017
mc1 <- mc %>% dplyr::select(Plot,Postfire_year,fitted_Tmax) %>%
  group_by(Plot,Postfire_year)%>%
  summarize(Tmax.abs = max(fitted_Tmax),
            Tmax.ann = mean(fitted_Tmax))%>%
  pivot_wider(id_cols="Plot",names_from = "Postfire_year",values_from = c("Tmax.abs","Tmax.ann"))
mc1 <- merge(mc1,plot[,c("Plot","Severity")])
mc1 <- mc1[-which(mc1$Severity == "Unburned"),c(1:7)]
m.z <-mc1 %>% mutate(across(!Plot,~as.numeric(scale(.x)))) %>%
  pivot_longer(cols = starts_with("Tmax"),values_to = "Tmax",names_to = c("MCType","Postfire_year"), 
               names_prefix = "Tmax_",names_sep = "_") %>%
  pivot_wider(names_from = "MCType",values_from = "Tmax",names_prefix = "z.")
mc1 <- mc1 %>%
  pivot_longer(cols = starts_with("Tmax"),values_to = "Tmax",names_to = c("MCType","Postfire_year"), 
               names_prefix = "Tmax_",names_sep = "_") %>%
  pivot_wider(names_from = "MCType",values_from = "Tmax")

#Combine with recruitment data
r <- merge(r, as.data.frame(m.z), by=c("Plot","Postfire_year"))
  
#remove year 1 data from plots sampled early in the season 
r <- merge(r,plot[,c("Plot","Fire","Date_1")],by="Plot")
r$Date_1 <- as.Date(r$Date_1,"%m/%d/%Y")
r <- r[-which(r$Postfire_year==1 & r$Date_1<"2018-07-15"),]
r$logArea <- log(r$TransectSize)
r$Postfire_year <- as.factor(r$Postfire_year)


#
#
#
########################################################Final model 

theme_set(theme_classic())

m = glmer.nb(Germ~offset(logArea)
             +z.ba_PICO
             +z.DSS
             +z.HLI
             +z.ppt_ann
             +z.Tmax.ann
             +z.Axis1
             +z.cover.Forb
             +z.soilN
             +Postfire_year
             +z.ppt_ann:z.HLI
             +z.ba_PICO:Postfire_year
             +z.HLI:z.cover.Forb
             +(1|Fire/Plot)
             ,data = r)#, control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

#Model diagnostics
summary(m)
CV(m,"Germ",r,4,1000) #rho-CV = 0.60
#without ppt or HLI: mean cor. 0.58, median 0.59
cor.test(predict(m),log(r$Germ), method = "spearman") #rho = 0.89
#without ppt or HLI: 0.89
r.squaredGLMM(m) #marginal r2 = 0.72
#without ppt or HLI: marginal R2 = 0.67

#Check residuals
hist(resid(m))
op <- par(mfrow=c(1,2),mar=c(4,4,1,2))
scatter.smooth(predict(m),resid(m))
abline(a=0,b=0,col="red")
plot(x=predict(m),y=log(r$Germ+1))
l = lm(log(r$Germ+1)~predict(m))
abline(a=l$coefficients[[1]], b= l$coefficients[[2]])
abline(a=0,b=1,col="red") 

#
#
#
##############################################################Plotting

off <- mean(r$TransectSize) #Average transect size to use as offset

eff <- effect("z.ba_PICO:Postfire_year",m,se=T,offset = log(off))
effect <- data.frame(x = eff$x$z.ba_PICO, Postfire_year = eff$x$Postfire_year, fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(sqrt(b$ba_PICO*0.229568)))
mylabels <- c(0,10,50)
mybreaks <- scale(sqrt(mylabels),att$'scaled:center',att$'scaled:scale')[,1]
g1 <- ggplot(r,aes(x  = z.ba_PICO, y = (Germ/TransectSize)*10000))+#facet_wrap(~Yr, labeller = label_both,ncol=3)+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = x, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000,fill = Postfire_year),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = (exp(fit)/off)*10000,color = Postfire_year),method = "lm",se = F,size = 2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(expand = c(0,0),breaks = mybreaks, labels = mylabels)+#,
    #inverse = function(x){(sqrt(x/0.229568) - mean(b$ba_PICO_sqrt))*sd(b$ba_PICO_sqrt)})+
  scale_fill_manual(values = c( "#34a0a4","#1a759f","#184e77"))+
  scale_color_manual(values = c("#34a0a4","#1a759f","#184e77"))+
  #coord_equal(ratio = 1,ylim = c(10,500000))+
  labs(x = expression("BA PICO ("*m^2~ha^-2~")"),#Basal area PICO (z-score)",
       y = expression("Expected regeneration ( # "~ha^-1~")"),fill = "Year",color = "Year")+
  theme(legend.position = "top")#,legend.background = element_rect(color = "black"))
g1

eff <- predictorEffect("z.ppt_ann",m,se=T,offset = log(off),xlevels = list(z.HLI = c(-1.87,1.83)))
effect <- data.frame(x = eff$x$z.ppt_ann, HLI = as.factor(eff$x$z.HLI), fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(b$ppt_ann))
mylabels <- c(600,1200,1800)
mybreaks <- scale(mylabels,att$'scaled:center',att$'scaled:scale')[,1]
g2 <- ggplot(data = effect)+#facet_grid(~HLI,labeller = label_both)+
  geom_rug(r,mapping = aes(x  = z.ppt_ann, y = (Germ/TransectSize)*10000),size = 0.3,position = "jitter",inherit.aes = FALSE)+
  geom_ribbon(data = effect, aes(x = x, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000,fill=HLI ),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = (exp(fit)/off)*10000,color=HLI),method = "lm",se = F,size = 2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  scale_fill_manual(values = c("#fb8500","#a63d40"),labels = c("0.4","1"))+ #"#FFC914","#E4572E"
  scale_color_manual(values = c("#fb8500","#a63d40"),labels = c("0.4","1"))+
  labs(x = "Annual precip. (mm)",y = expression("Predicted regeneration ( # "~ha^-1~")"),color = "HLI",fill="HLI")+
  theme(legend.position = "top")#,legend.background = element_rect(color = "black"))
g2

eff <- predictorEffect("z.cover.Forb",m,se=T,offset = log(off),xlevels = list(z.HLI = c(-1.87,1.83)))
effect <- data.frame(x = eff$x$z.cover.Forb, HLI = as.factor(eff$x$z.HLI), fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(sqrt(sp$cover.Forb)))
mylabels <- c(1,10,30,50)
mybreaks <- scale(sqrt(mylabels),att$'scaled:center',att$'scaled:scale')[,1]
g3 <- ggplot(data = effect)+#facet_grid(~HLI,labeller = label_both)+
  geom_rug(r,mapping = aes(x  = z.cover.Forb, y = (Germ/TransectSize)*10000),size = 0.3,position = "jitter",inherit.aes = FALSE)+
  geom_ribbon(data = effect, aes(x = x, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000,fill = HLI),
              alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = (exp(fit)/off)*10000, color = HLI),method = "lm",se = F,size = 2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  scale_fill_manual(values = c("#fb8500","#a63d40"),labels = c("0.4","1"))+ #"#FFC914","#E4572E"
  scale_color_manual(values = c("#fb8500","#a63d40"),labels = c("0.4","1"))+
  labs(x = "Forb cover (%)",y = expression("Predicted regeneration ( # "~ha^-1~")"),color = "HLI",fill="HLI")+
  theme(legend.position = "none")#,legend.background = element_rect(color = "black"))
g3

eff <- effect("z.DSS",m,se=T,offset = log(off),xlevels = list(z.DSS = c(min(r$z.DSS),-2,-1,0,1,2,max(r$z.DSS))))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(sqrt(b$DSS)))
mylabels <- c(10,50,100)
mybreaks <- scale(sqrt(mylabels),att$'scaled:center',att$'scaled:scale')[,1]
g4 <- ggplot(r,aes(x  = z.DSS, y = (Germ/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.DSS, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3, color = "grey90")+
  geom_smooth(data = effect, aes(x = z.DSS,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(expand = c(0,0),limits = c(min(r$z.DSS-.4),max(r$z.DSS+.1)),labels = mylabels, breaks = mybreaks)+
  labs(x = "DSS (m)",
       y = expression("Predicted regeneration ( # "~ha^-1~")"))
g4

eff <- effect("z.Axis1",m,se=T,offset = log(off),xlevels = list(z.Axis1 = c(min(r$z.Axis1),-2,-1,0,1,2,max(r$z.Axis1))))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(b$Axis1))
mylabels <- seq(-2,4,2)
mybreaks <- scale(mylabels,att$'scaled:center',att$'scaled:scale')[,1]
g5 <- ggplot(r,aes(x  = z.Axis1, y = (Germ/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.Axis1, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3, color = "grey90")+
  geom_smooth(data = effect, aes(x = z.Axis1,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size =2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  labs(x = expression("Fire severity ("*italic("Axis1")*")"),
       y = expression("Predicted regeneration ( # "~ha^-1~")"))
g5

eff <- effect("z.Tmax.ann",m,se=T,offset = log(off),xlevels = list(z.Tmax.ann = c(min(r$z.Tmax.ann),-2,-1,0,1,2,max(r$z.Tmax.ann))))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(mc1$Tmax.ann))
mylabels <- seq(0,40,5)
mybreaks <- scale(mylabels,att$'scaled:center',att$'scaled:scale')[,1]
g6 <- ggplot(r,aes(x  = z.Tmax.ann, y = (Germ/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.Tmax.ann, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3, color = "grey90")+
  geom_smooth(data = effect, aes(x = z.Tmax.ann,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000))+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  labs(x =expression("Microclim. ("*italic(T[avg])*","~degree*"C)"),y = expression("Predicted regeneration ( # "~ha^-1~")"))
g6

eff <- effect("z.soilN",m,se=T,offset = log(off),xlevels = list(z.soilN = c(min(r$z.soilN),-1,0,1,2,2.5)))
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
att <- attributes(scale(sqrt(b$soilN)))
mylabels <- c(0.1,0.5,1,2)
mybreaks <- scale(sqrt(mylabels),att$'scaled:center',att$'scaled:scale')[,1]
g7 <- ggplot(r,aes(x  = z.soilN, y = (Germ/TransectSize)*10000))+
  geom_rug(size = 0.3,position = "jitter")+
  geom_ribbon(data = effect, aes(x = z.soilN, y = (exp(fit)/off)*10000, 
                                 ymin = (exp(lower)/off)*10000,ymax = (exp(upper)/off)*10000),
              alpha = 0.3, color = "grey90")+
  geom_smooth(data = effect, aes(x = z.soilN,y = (exp(fit)/off)*10000),method = "lm",se = F,color = "black",size = 2)+
  scale_y_log10(breaks = c(10,100,1000,10000,100000),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  coord_cartesian(ylim = c(10,500000),xlim = c(-1.4,2.5))+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  labs(x = expression("Soil N avail. ("*mu*"g N"~day^-1~")"),y = expression("Predicted regeneration ( # "~ha^-1~")"))
g7

g2$labels$y<-g3$labels$y<-g4$labels$y<-g5$labels$y<-g6$labels$y<-g7$labels$y <- NULL

p1 <- (g2/g3)|(g4/g6)|(g5/g7)
g1+p1+plot_layout(guides = "collect",design = layout,heights = c(1,1))+plot_annotation(tag_levels= list(c("A","B","E","C","F","D","G"))) & theme(
  legend.position = "right",plot.tag.position = c(0,1.02),
  axis.title.y = element_text(size = 13)) 
#ggsave("Figure5.png",dpi = 600, width = 10, height = 5, units = "in")

#
#
#



####Figure 6 ############################################################################################
#Statistical model of seedling survivorship (all species) in burned plots from year 2-3 post-fire
theme_set(theme_classic())

#####################################################Set up data

sd[,c("Subplot","Seedling_ID")] = sapply(sd[,c("Subplot","Seedling_ID")],as.factor)

#Remove first year seedlings (to focus on interannual mortality/survivorship)
d <- sd[-which(sd$Age==1),]

###########Merge with scaled plot data - need data frame "b" from Initial Data Setup section above
d <- merge(sd,b[,c("Plot","Severity","z.PP","z.DSS","z.ba_PICO","z.Live_ba_PICO","z.Live_ba_PSME","z.Live_ba_LAOC",
                   "z.HLI","z.ppt_ann","z.ppt_JJAS","z.ppt_pf","z.DEF","z.Tmax","z.Tmax_avg",
                   "z.Axis1","z.Axis2","z.dnbr","z.NO3","z.NH4","z.soilN")],by="Plot")

###########Set up subplot data (retain individual subplot data)
  #Calculate average cover values across all three measurement years
sp.avg <- sp[,-3] %>% group_by(Plot,Subplot)%>%
  summarize_all(function(x){mean(x,na.rm=T)} )
names(sp.avg)[3:15] <- paste0(names(sp.avg)[3:15],"_Avg")
sp <- merge(sp,sp.avg,by=c("Plot","Subplot"))

  #Use averages to fill in missing annual values
sp[which(is.na(sp$cover.BG)),"cover.BG"] = sp[which(is.na(sp$cover.BG)),"cover.BG_Avg"] 
sp[which(is.na(sp$cover.Wood)),"cover.Wood"] = sp[which(is.na(sp$cover.Wood)),"cover.Wood_Avg"] 
sp[which(is.na(sp$cover.Litter)),"cover.Litter"] = sp[which(is.na(sp$cover.Litter)),"cover.Litter_Avg"] 
sp[which(is.na(sp$cover.Grass)),"cover.Grass"] = sp[which(is.na(sp$cover.Grass)),"cover.Grass_Avg"] 
sp[which(is.na(sp$cover.Forb)),"cover.Forb"] = sp[which(is.na(sp$cover.Forb)),"cover.Forb_Avg"] 
sp[which(is.na(sp$cover.Shrub)),"cover.Shrub"] = sp[which(is.na(sp$cover.Shrub)),"cover.Shrub_Avg"] 
sp[which(is.na(sp$cover.Moss)),"cover.Moss"] = sp[which(is.na(sp$cover.Moss)),"cover.Moss_Avg"] 
sp[which(is.na(sp$cover.Seedling)),"cover.Seedling"] = sp[which(is.na(sp$cover.Seedling)),"cover.Seedling_Avg"]
sp[which(is.na(sp$cover.Sapling)),"cover.Sapling"] = sp[which(is.na(sp$cover.Sapling)),"cover.Sapling_Avg"] 
sp[which(is.na(sp$cover.Veg)),"cover.Veg"] = sp[which(is.na(sp$cover.Veg)),"cover.Veg_Avg"] 
sp[which(is.na(sp$Canopy.tot)),"Canopy.tot"] = sp[which(is.na(sp$Canopy.tot)),"Canopy.tot_Avg"] 
sp[which(is.na(sp$Canopy.green)),"Canopy.green"] = sp[which(is.na(sp$Canopy.green)),"Canopy.green_Avg"] 
sp[which(is.na(sp$Canopy.dead)),"Canopy.dead"] = sp[which(is.na(sp$Canopy.dead)),"Canopy.dead_Avg"] 

  #Calculate average Moss & Forb cover - analogous to Axis2 values, but annual and subplot-specific (more fine-grained data)
sp$MossForb <- (sp$cover.Moss+sp$cover.Forb)/2

  #Transform as necessary to reduce skewness
sq = c("MossForb","cover.BG","cover.Wood", "cover.Grass","cover.Forb",
       "cover.Shrub","cover.Moss")
sp[,sq] = sapply(sp[,sq],sqrt)

  #Restrict subplots to those with seedling monitoring data
sp$PlotSubplot <- paste0(sp$Plot,sp$Subplot,sp$Postfire_year)
sp <- sp[which(sp$PlotSubplot %in% unique(paste0(d$Plot,d$Subplot,d$Postfire_year))),]
sp.avg$PlotSubplot <- paste0(sp.avg$Plot,sp.avg$Subplot)
sp.avg <- sp.avg[which(sp.avg$PlotSubplot %in% unique(paste0(d$Plot,d$Subplot))),]

  #Scale subplot data (z scores)
spz <- cbind(sp[,1:3],scale(sp[,c(4:16,30)],center = T))
spz.avg <- cbind(as.data.frame(sp.avg[,1:2]),scale(sp.avg[,3:15],center = T))
spz <- merge(spz,spz.avg,by=c("Plot","Subplot"))
names(spz)[4:30]<- paste0("z.",names(spz)[4:30])

  #merge with regeneration data
d <- merge(d,spz, by = c("Plot","Subplot","Postfire_year") )

########Set up microclimate data
m <- mc %>% group_by(Plot, Year)%>%
  summarize(Tmax.yr = max(fitted_Tmax),
            Tmax.avg.yr = mean(fitted_Tmax))
m$Postfire_year <- m$Year - 2017
m$PlotYr <- paste0(m$Plot,m$Postfire_year)
m <- m[which(m$PlotYr %in% unique(paste0(d$Plot,d$Postfire_year))),]
m[,3:4] <- scale(m[,3:4],center=T)
names(m)[3:4] <- paste0("z.",names(m)[3:4])

  #merge with regeneration data
d <- merge(d,m[,c(1,3,4,5)],by=c("Plot","Postfire_year"))

########Calculate number of seedlings dead and alive in each year
d1 = d[which(d$Postfire_year == 3),]
d2 = d1 %>% group_by(Plot,Subplot,Age) %>%
  summarize_at("Alive",sum)
d3 = d1 %>% group_by(Plot,Subplot,Age) %>%
  count()
d2 = merge(d2,d3,by=c("Plot","Subplot","Age"))
d2$n_alive = d2$Alive
d2$n_dead = d2$n - d2$n_alive

d3 <- merge(d2,unique(d1[,c(1:3,16:64)]),by=c("Plot","Subplot"),all=F)


#
#
#
#########################################################Final model
m = glmer(cbind(n_alive,n_dead)~
            +z.DEF
          +z.MossForb #Subplot level (year specific)
          +z.Canopy.tot #Subplot level (year specific)
          +z.NO3
          +Age
          +Age:z.DEF
          +(1|Plot)
          ,data=d3, family = binomial(link="logit"),control =glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))

#Model diagnostics
summary(m)
check_convergence(m)
r= roc(d1$Alive,predict(m,newdata=d1,type="response",re.form=NA))
plot(r);r
dat.sim = simulateResiduals(m1)
plot(dat.sim) #these look good
1 - pchisq(sum(resid(m, type = "pearson")^2), summary(m)$AICtab[[5]]) #GOF test with pearson residuals - no significant patterns remaining
testUniformity(dat.sim) #residuals are uniform
testDispersion(simulateResiduals(m1))#, refit = T)) #no evidence of overdispersion
r.squaredGLMM(m) #marginal r2 = 0.53 (theoretical), 0.34 (delta - based on observation-level variance)

#
#
#
###################################################Plotting
eff <- predictorEffect("z.DEF",m,se=T,xlevels = 2)
effect <- data.frame(x = eff$x$z.DEF, Age = factor(eff$x$Age,labels = c("1-2","2-3")), fit = logit2prob(eff$fit)*100, lower = logit2prob(eff$lower)*100, upper = logit2prob(eff$upper)*100)
att <- attributes(scale(b$DEF))
mylabels <- seq(200,500,50)
mybreaks <- scale(mylabels,att$'scaled:center',att$'scaled:scale')[,1]
g1 <- ggplot()+#facet_wrap(~Yr, labeller = label_both,ncol=3)+
  geom_rug(d3,mapping = aes(x  = z.DEF, y = (n_alive/(n_alive+n_dead))*100,size= n_alive+n_dead),inherit.aes = FALSE)+
  geom_ribbon(data = effect, aes(x = x, y = fit , ymin = lower,ymax = upper,
                                 fill = Age),alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = fit,color = Age),se = F,size = 1)+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  scale_fill_manual(values = c( "#34a0a4","#184e77"), labels = c("1-2","2-3"))+
  scale_color_manual(values = c("#34a0a4","#184e77"), labels = c("1-2","2-3"))+
  labs(#x = expression("Cool-wet "%<->%" Warm-dry (DEF z-score)"),
       x = "DEF (mm)",y = "Predicted survivorship (%)",size = "n", color = "Age transition",fill = "Age transition")+
  theme(legend.position = c(0.25,0.35),legend.box = "horizontal", legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(1,1,1,1,unit= "pt") ,axis.title.x = element_text(size =14),axis.text = element_text(size = 11),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
g1

eff <- predictorEffect("z.MossForb",m,se=T,xlevels = list(z.MossForb = seq(min(d3$z.MossForb),max(d3$z.MossForb),.1)))
effect <- data.frame(x = eff$x$z.MossForb, fit = logit2prob(eff$fit)*100, lower = logit2prob(eff$lower)*100, upper = logit2prob(eff$upper)*100)
#att <- attributes(scale(sqrt((sp$MossForb^2)/2)))
att <- attributes(scale(sp$MossForb))
mylabels <- c(1,10,30,60)
mybreaks <- scale(sqrt(mylabels),att$'scaled:center',att$'scaled:scale')[,1]
g2 <- ggplot()+#facet_wrap(~Yr, labeller = label_both,ncol=3)+
  geom_rug(data =d3,mapping = aes(x  = z.MossForb, y = (n_alive/(n_alive+n_dead))*100,size= n_alive+n_dead),inherit.aes = FALSE,position="jitter")+
  geom_ribbon(data = effect, aes(x = x, y = fit , ymin = lower,ymax = upper),alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = fit),se = F,size = 1, color = "black")+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  labs(x = "Moss & forb cover (%)",y = "Predicted survivorship (%)",size = "n")+
  theme(legend.position = "none",axis.title.x = element_text(size =14),axis.text = element_text(size = 11))
g2

eff <- predictorEffect("z.Canopy.tot",m,se=T,xlevels = list())
effect <- data.frame(x = eff$x$z.Canopy.tot,fit = logit2prob(eff$fit)*100, lower = logit2prob(eff$lower)*100, upper = logit2prob(eff$upper)*100)
att <- attributes(scale(sp$Canopy.tot))
mylabels <- seq(0,80,20)
mybreaks <- scale(mylabels,att$'scaled:center',att$'scaled:scale')[,1]
g3 <- ggplot()+#facet_wrap(~Yr, labeller = label_both,ncol=3)+
  geom_rug(d3,mapping = aes(x  = z.Canopy.tot, y = (n_alive/(n_alive+n_dead))*100,size= n_alive+n_dead),inherit.aes = FALSE)+
  geom_ribbon(data = effect, aes(x = x, y = fit , ymin = lower,ymax = upper),alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = fit),se = F,size = 1, color = "black")+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  labs(x = "Canopy cover (%)",y = "Predicted survivorship (%)",size = "n")+
  theme(legend.position = "none",axis.title.x = element_text(size =14),axis.text = element_text(size = 11))
g3

eff <- predictorEffect("z.NO3",m,se=T,xlevels = list())
effect <- data.frame(x = eff$x$z.NO3, fit = logit2prob(eff$fit)*100, lower = logit2prob(eff$lower)*100, upper = logit2prob(eff$upper)*100)
att <- attributes(scale(sqrt(b$NO3)))
mylabels <- c(0,0.1,0.5,1)
mybreaks <- scale(sqrt(mylabels),att$'scaled:center',att$'scaled:scale')[,1]
g4 <- ggplot()+#facet_wrap(~Yr, labeller = label_both,ncol=3)+
  geom_rug(d3,mapping = aes(x  = z.NO3,y = (n_alive/(n_alive+n_dead))*100,size= n_alive+n_dead),inherit.aes = FALSE)+
  geom_ribbon(data = effect, aes(x = x, y = fit , ymin = lower,ymax = upper),alpha = 0.3)+
  geom_smooth(data = effect, aes(x = x,y = fit),se = F,size = 1, color = "black")+
  scale_x_continuous(expand = c(0,0),labels = mylabels, breaks = mybreaks)+
  labs(x = expression(NO[3]~"("*mu*"g N"~day^-1~")"),y = "Predicted survivorship (%)",size = "n")+
  theme(legend.position = "none",axis.title.x = element_text(size =14),axis.text = element_text(size = 11))
g4

p1 = ggplot(data.frame( x = 1, y = 1)) +
  annotate("text", label = "Expected survivorship (%)" , angle = 90, y = 1, x = 1,size = 6) + 
  theme_void() +
  coord_cartesian(clip = "off")

g1$labels$y <-g2$labels$y<-g3$labels$y<-g4$labels$y <- " "

p1+((g1|g2)/(g3|g4))+ plot_layout(widths = c(1,24),guides = "collect")+ plot_annotation(tag_levels = list(c("","A","B","C","D"))) & theme(
  plot.tag.position = c(0.1,1.05),plot.margin = margin(t=10, unit = "pt"))

#ggsave("Figure_6.png",dpi = 600, width = 8, height = 6.8, units = "in")


#
#
#
#
####Seedling height statistical model (Figure B.5) ######################################################################################
#############################################Model of seedling height 

##################Set up data - need data frame 'a' from section 'Initial data setup' above
h <- merge(sd,a[,c("Plot","Severity","z.PP","z.DSS","z.ba_PICO","z.Live_ba_PICO","z.Live_ba_PSME","z.Live_ba_LAOC",
                   "z.HLI","z.ppt_ann","z.ppt_JJAS","z.ppt_pf","z.DEF","z.Tmax","z.Tmax_avg",
                   "z.Axis1","z.Axis2","z.dnbr","z.NO3","z.NH4","z.soilN")],by="Plot")

h <-h %>% filter(Alive == 1) %>% filter(!is.na(Ht))
h[which(h$Spp %in% c("ABLA","ABGR")),"Spp"] <- "ABIES"
h <- h[which(h$Spp %in% c("ABIES","LAOC","PICO","PIEN","PIPO","PSME")),]
h <- h %>% filter(Postfire_year == 3)

m <- glmer(Ht ~ Spp 
           + Age 
           + z.Tmax_avg 
           + z.Axis2
           +z.Tmax_avg:Age
           + (1|Plot/Subplot)
           , data = h, family = Gamma(link = "log"),control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun =2e6)))

eff <- predictorEffect("Spp",m,se=T)
effect <- data.frame(Spp = as.factor(eff$x$Spp), fit =eff$fit, lower = eff$lower, upper = eff$upper)
g1 <- ggplot()+
  geom_errorbar(data = effect, aes(x = Spp, y = exp(fit), ymin = exp(lower),ymax = exp(upper)))+
  geom_point(data = effect, aes(x = Spp,y = exp(fit)),size = 1.5)+
  scale_y_continuous(expand = c(0,0),breaks = seq(2,8,2),limits = c(2,9))+
  scale_x_discrete(labels = c("Abies","LAOC","PICO","PIEN","PIPO","PSME"))+
  #scale_fill_manual(values = c("#09814a","#219ebc","#023047","#C1A5A9","#fb8500","#a63d40"))+
  #scale_color_manual(values = c("#09814a","#219ebc","#023047","#C1A5A9","#fb8500","#a63d40"))+
  labs(x = NULL,y = "Predicted height (cm)")+theme(axis.text.x = element_text(colour = "black",size = 10))
g1

eff <- effect("z.Axis2",m,se=T,xlevels=20)
effect <- data.frame(x = eff$x, fit =eff$fit, lower = eff$lower, upper = eff$upper)
g3 <- ggplot()+
  geom_ribbon(data = effect, aes(x = z.Axis2, y = exp(fit), ymin = exp(lower),ymax = exp(upper)),alpha = 0.2)+
  geom_smooth(data = effect, aes(x = z.Axis2,y = exp(fit)),se = F,size = 1.5,color="black")+
  # scale_y_log10(labels = comma)+
  scale_x_reverse(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),breaks = seq(2,8,2),limits = c(2,9))+
  # coord_cartesian(ylim = c(10,500000))+
  labs(x = expression("Moss + forb cover ("*italic(Axis2)*" z-score)"), y = "Predicted height (cm)")
g3

eff <- predictorEffect("Age",m,se=T, xlevels = list(z.Tmax_avg = c(-2,2)))
effect <- data.frame(x = eff$x$Age, Tmax = as.factor(eff$x$z.Tmax_avg), fit =eff$fit, lower = eff$lower, upper = eff$upper)
g4 <- ggplot()+
  geom_ribbon(data = effect, aes(x = x, y = exp(fit), ymin = exp(lower),ymax = exp(upper),fill = Tmax),alpha = 0.2)+
  geom_smooth(data = effect, aes(x = x,y = exp(fit),color = Tmax),se = F,size = 1.5)+
  scale_x_continuous(expand = c(0,0),breaks = c(1,2,3))+
  scale_y_continuous(expand = c(0,0),breaks = seq(2,14,2))+
  scale_fill_manual(values = c("#FFC914","#E4572E"))+
  scale_color_manual(values = c("#FFC914","#E4572E"))+
  labs(fill= expression(italic(T[avg])*" z-score"),color= expression(italic(T[avg])*" z-score"),
       x = "Age (year)", y = "Predicted height (cm)")+
  theme(legend.position = c(0.3,0.8),legend.background = element_rect(color = "black"))
g4

p1 = ggplot(data.frame( x = 1, y = 1)) +
  annotate("text", label = "Expected height (cm)", angle = 90, y = 1, x = 1,size = 5) + 
  theme_void() +
  coord_cartesian(clip = "off")

g1$labels$y <-  g3$labels$y <-g4$labels$y <- NULL

p1+ g1/(g3|g4) + plot_layout(widths = c(1,22)) + plot_annotation(tag_levels = list(c("","A","B","C")))
#ggsave("Figure_B-5.jpg",dpi = 600, width = 6, height = 6, units = "in")

