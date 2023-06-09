#This R script is part of the manuscript entitled:
#'Functional diversity of the avian assemblages better than 
#taxonomic diversity in capturing the effect of protection'

##Load R packages and input files
{
rm(list = ls()) 
set.seed (123)

library(tidyverse)    #For data handling and processing
library(vegan)        #For estimating diversity indices
library(codyn)        #For estimating differences in species abundance
library(ggbiplot)     #For plotting figures
library(gridExtra)    #For plotting figures
library(ggplot2)      #For plotting figures
library(ggpubr)       #For plotting p-values on boxplots
library(FD)           #For estimating Functional diversity
library(mgcv)         #For GAM analysis
library(tidymv)       #For plotting GAM results
library(indicspecies) #For Indicator species analysis
library(permute)      #For Indicator species analysis

#resolve the function conflict
filter <- dplyr::filter
select <- dplyr::select
summarise<-dplyr::summarise
mutate<-dplyr::mutate
rename<-dplyr::rename

#Read the KBA dataset and clean the data
#All "Spuhs" and "Slashes" are to be eliminated except few.
#Little spiderhunter is the only spiderhunter species in Kerala, hence spiderhunter sp. changed to 'Little spiderhunter'
#Green/Greenish warbler and Booted/Sykes warbler merged into single taxa, respectively. 
#Saunder's Tern and Square-tailed Drongo-Cuckoo are not found in Kerala
#eliminate all nocturnal species, ambiguous species (sps. and slashes)

file1<-readRDS("kba_data.rds")              #contains Kerala Bird Atlas details
file2<-readRDS("kba_names.rds")             #contains information on ambiguous taxa to be eliminated
}

##Kerala Bird Atlas Data Filtering
{
if(is.null(file2$Action) == F) {
  file1$Common.Name<- file2$Assigned.Name.for.Atlas[match(file1$Common.Name,file2$Common.Name)]
  file1$Common.Name[file1$Common.Name == ''] <- 'Removed'
}

##File for analysis
KBA <- file1 %>% filter(!Common.Name == 'Removed')%>%unique()%>%
  mutate(Common.Name=replace(Common.Name,Common.Name == "Green/Greenish Warbler","Greenish Warbler"),
         Common.Name=replace(Common.Name,Common.Name == "Booted/Sykes's Warbler","Booted Warbler"))

remove(file1,file2)
}

##Select Protected Area (PA) and Reserved forests (RF) subcells
#load subcell file
#exclude areas of elevation below 200m and above 2000 m
#randomly select 300 subcells within PAs and RFs

{
  PA<- readRDS("subcells.rds")%>%filter(Treatment=="PA")%>%select(Sub.cell,DEM)%>%
    filter(DEM<=2000 & DEM>=200)%>%slice_sample(n=300, replace=F)%>%select(Sub.cell)
  
  RF<- readRDS("subcells.rds")%>%filter(Treatment=="RF")%>%select(Sub.cell,DEM)%>%
    filter(DEM<=2000 & DEM>=200)%>%slice_sample(n=300, replace=F)%>%select(Sub.cell)
  
 PA.kba<-KBA%>%dplyr::filter(Season=="Dry")%>%
    group_by(Common.Name)%>%mutate(count=n())%>%filter(!count==1)%>%right_join(PA)%>%unique()%>%
    select(List.ID, Common.Name,Sub.cell)%>%unique()%>%mutate(pres_abs=1)
  
  RF.kba<-KBA%>%dplyr::filter(Season=="Dry")%>%
    group_by(Common.Name)%>%mutate(count=n())%>%filter(!count==1)%>%right_join(RF)%>%unique()%>%
    select(List.ID, Common.Name,Sub.cell)%>%unique()%>%mutate(pres_abs=1)
}

##Prepare file for various analysis 
#(Species [row] by sites matrix [column])
{
  PA_1<-PA.kba%>%group_by(Sub.cell)%>%mutate(list=n_distinct(List.ID))%>%
    group_by(Sub.cell,Common.Name)%>%summarise(freq=n()/list)%>%unique()%>%
    spread(key = Common.Name, value = freq)%>%replace(is.na(.), 0)%>%
    column_to_rownames(., var = "Sub.cell")
  
  RF_1<-RF.kba%>%group_by(Sub.cell)%>%mutate(list=n_distinct(List.ID))%>%
    group_by(Sub.cell,Common.Name)%>%summarise(freq=n()/list)%>%unique()%>%
    spread(key = Common.Name, value = freq)%>%replace(is.na(.), 0)%>%
    column_to_rownames(., var = "Sub.cell")
  
  file<-bind_rows(RF_1,PA_1)%>%replace(is.na(.), 0)
  groups=c(rep("RF",300),rep("PS",300)) #treatment
}

##Figure S1
# species accumulation curve across samples
{
  par(mfrow = c(1, 1), mar=c(3.5,3.5,1,1),mgp = c(3, 0.5, 0))
  plot(specaccum(RF_1), xlab = "", ylab = "",col = "Grey50",xlim=c(0,300),ylim=c(0,260))
  par(new=TRUE)
  plot(specaccum(PA_1), xlab = "", ylab = "", col = "Grey10",xlim=c(0,300),ylim=c(0,260),yaxt="n",xaxt="n",axes=FALSE)
  title(xlab = "sites count", line = 2.2)
  title(ylab = "Species count", line = 2.2)
  
  legend(x = "bottomright", c("Protected Areas", " Reserve Forests"), 
         fill=c("Grey10","Grey50")) 
}

##Figure 2
#Comparison across PA and RF sites
#Plot elevation, latitude and longitude of PA/RF
#LULC layers and Bioclimatic layers across PA and RF

{
  par(mfrow = c(3, 3), mar=c(2,2,2,1),mgp = c(3, 0.5, 0))
  a<-readRDS("subcells.rds")%>%right_join(PA)
  b<-readRDS("subcells.rds")%>%right_join(RF)
  
  atr<-rbind(a,b) 
  remove(a,b)
  
  plot1<-ggplot(atr,aes(x=Treatment, y=DEM))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Elevation (in metres)")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test")
  
  plot2<-ggplot(atr,aes(x=Treatment, y=lon))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Longitude")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+ 
    stat_compare_means(method = "wilcox.test")
  
  plot3<-ggplot(atr,aes(x=Treatment, y=lat))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Latitude")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+ 
    stat_compare_means(method = "wilcox.test")
  
  plot4<-ggplot(atr,aes(x=Treatment, y=Wet.forest))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Wet Forest")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+ 
    stat_compare_means(method = "wilcox.test")
  
  plot5<-ggplot(atr,aes(x=Treatment, y=Dry.forest))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Dry forest")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test")
  
  plot6<-ggplot(atr,aes(x=Treatment, y=Open.area))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Open area")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test")
  
  plot7<-ggplot(atr,aes(x=Treatment, y=Water.body))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Waterbody")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+ 
    stat_compare_means(method = "wilcox.test")
  
  plot8<-ggplot(atr,aes(x=Treatment, y=Plantations))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Plantations")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test")
  
  plot9<-ggplot(atr,aes(x=Treatment, y=Agriculture))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Agriculture")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+ 
    stat_compare_means(method = "wilcox.test")
  
  Chelsa<-atr%>%filter(!is.na(bio2))
  
  #Perform PCA
  
  dataset = data.frame(Chelsa[,c(13:30)])
  cat = factor(Chelsa[,1])
  pca = prcomp(dataset,scale=T)
  summary(pca)
  plot10<-ggbiplot(pca,choices=c(1,2),obs.scale=0,var.scale=1,groups=cat,ellipse=T,
                   circle=T,var.axes= F)+theme_test()+
    scale_colour_manual(labels = c("Protected areas (PA)", "Reserve forests (RF)"), 
                        values = c("Grey10", "Grey50"))+
    theme(legend.text = element_text(colour="black", size=11))+
    theme(axis.text = element_text(size = 10, face = "plain", colour = "black"))+
    theme(axis.title = element_text(size = 12, face = "plain"))+
    ggtitle("PCA on Bioclimatic variables")+ 
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.title = element_blank(),
          legend.box.margin = margin(0.2,0.2,0.2,0.2,"cm"),
          legend.spacing.y = unit(0, "mm"),
          panel.border = element_rect(colour = "black", fill=NA),
          aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
          legend.position = c(0.8, 0.067),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))

  #Default ggplot colors are Green #00BA38 Red #F8766D  Blue #619CFF
  
  plota<-grid.arrange(plot1, plot2, plot3, ncol=3)
  plotb<-grid.arrange(plot4, plot5, plot6, ncol=3)
  plotc<-grid.arrange(plot7, plot8, plot9, ncol=3)
  plotd<-grid.arrange(plota, plotb, plotc, nrow=3)
  plote<-grid.arrange(plotd, plot10, ncol=2)
  
}
remove(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plota,plotb,plotc,plotd)

##Figure 3
#Plot alpha and beta diversity indices
#Perform ANOSIM analsis to test if the two communities are statistically different

{
  
  Alpha1<-specnumber(file)%>%as.data.frame()%>%rename('MSC'='.')%>%
    rownames_to_column("Sub.cell")%>%left_join(readRDS("subcells.rds"))%>%
    select(Sub.cell,MSC,Treatment)%>%unique()
  
  Alpha2<-diversity(file, "invsimpson")%>%as.data.frame()%>%rename('IVS'='.')%>%
    rownames_to_column("Sub.cell")%>%left_join(Alpha1)
  
  Alpha3<-diversity(file, "shannon")%>%as.data.frame()%>%rename('SNN'='.')%>%
    rownames_to_column("Sub.cell")%>%left_join(Alpha2)
  
  Alpha5<-diversity(file)/log(specnumber(file))%>%as.data.frame()%>%rename('EVN'='.')
  Alpha5<-Alpha5%>%rownames_to_column("Sub.cell")%>%left_join(Alpha3)
  
  #Shapiro test of normality to check if the values are normally distributed
  #P value<0.05 suggests that data is not normal
  
  AlphaRF<-Alpha5%>%filter(Treatment=="RF")
  AlphaPA<-Alpha5%>%filter(Treatment=="PA") 
  
  shapiro.test(AlphaRF$MSC)
  shapiro.test(AlphaRF$EVN)
  shapiro.test(AlphaRF$SNN)
  shapiro.test(AlphaRF$IVS)
  
  shapiro.test(AlphaPA$MSC)
  shapiro.test(AlphaPA$EVN)
  shapiro.test(AlphaPA$SNN)
  shapiro.test(AlphaPA$IVS)
  
  remove(Alpha1,Alpha2,Alpha3,AlphaPA,AlphaRF)
  
  #As the data is not normal, use non-parametric test to compare means     
  
  #Estimate betadiversity
  #Run NMDS 
  kba_NMDS=metaMDS(file,dist = "bray", k=3,trymax=1000)
  
  
  #Plot a combined plot for alpha and betadiversity
  
  plot1<-ggplot(Alpha5,aes(x=Treatment, y=MSC))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Mean species count")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test",label.x = 1,label.y = 70)
  
  plot2<-ggplot(Alpha5,aes(x=Treatment, y=IVS))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Invsimpson index")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test",label.x = 1,label.y = 50)
  
  plot3<-ggplot(Alpha5,aes(x=Treatment, y=SNN))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("Shannon index")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test",label.x = 1,label.y = 3.5)
  
  plot4<-ggplot(Alpha5,aes(x=Treatment, y=EVN))+geom_boxplot()+theme_bw()+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    ggtitle("pielou's evenness")+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text =element_text(color = "black"))+
    stat_compare_means(method = "wilcox.test",label.x = 1,label.y = 0.98)
  
  plota<-grid.arrange(plot1, plot2,nrow=2)
  plotb<-grid.arrange(plot3, plot4, nrow=2)
  grid.arrange(plota, plotb, ncol=2)
  
  par(mfrow=c(2,1),mai = c(0.4, 0.5,0.1,0.1))
  
  plot5<-ordiplot(kba_NMDS, type = "none", xlim = c(-1,1),ylim = c(-1,1))
  points(plot5, "sites", pch = 19, col = "Grey10", select = groups == "PS")
  points(plot5, "sites", pch = 19, col = "Grey50", select = groups  == "RF")
  
  # add confidence ellipses around habitat types
  ordiellipse(kba_NMDS, groups, conf = 0.95, label = TRUE)
  legend(x = "topleft", c("Protected areas", " Reserve Forests"),fill=c("Grey10","Grey50")) 
  
  stressplot(kba_NMDS)
  
  #Test if there is statistical difference between the two communities
  
  group_n<-readRDS("subcells.rds")%>%select(Treatment,Sub.cell)%>%unique()
  
  file_m<-bind_rows(PA_1,RF_1)%>%replace(is.na(.), 0)%>%
    rownames_to_column("Sub.cell")%>%left_join(group_n)%>%relocate(Treatment)
  
  com = as.matrix(file_m[,3:ncol(file_m)])
  anosim(com, file_m$Treatment, distance = "bray", permutations = 9999)
  
}
remove(plot1,plot2,plot3,plot4,plot5,plota,plotb)

##Figure 4
#Estimate Functional diversity
#Run GAM model

{
  #Get DEM values
  a<-readRDS("subcells.rds")%>%right_join(PA)
  b<-readRDS("subcells.rds")%>%right_join(RF)
  
  atr<-rbind(a,b)%>%select(Treatment,Sub.cell,lon,lat,DEM)%>%unique() 
  remove(a,b)
  
  #Calculate functional diversity
  #for RF
  sp.RF<-RF.kba%>%select(Common.Name)%>%unique()
  trait.RF <- readRDS("kba_species.details.rds")%>%select(-c(2:11))%>%
    filter(duplicated(Common.Name) == FALSE)%>%right_join(sp.RF)
  trait.RF$Common.Name<-sort(trait.RF$Common.Name, decreasing=TRUE)
  trait.RF<-trait.RF%>%column_to_rownames(., var = "Common.Name")
  abund.RF<-RF_1%>%select(order(desc(colnames(.))))
  fd.RF<- dbFD(trait.RF, abund.RF,  w.abun = TRUE, corr = "cailliez",m="min",stand.x = TRUE,print.pco=TRUE)
  
  #for PA
  sp.PA<-PA.kba%>%select(Common.Name)%>%unique()
  trait.PA <- readRDS("kba_species.details.rds")%>%select(-c(2:11))%>%
    filter(duplicated(Common.Name) == FALSE)%>%right_join(sp.PA)
  trait.PA$Common.Name<-sort(trait.PA$Common.Name, decreasing=TRUE)
  trait.PA<-trait.PA%>%column_to_rownames(., var = "Common.Name")
  abund.PA<-PA_1%>%select(order(desc(colnames(.))))
  fd.PA<- dbFD(trait.PA, abund.PA,  w.abun = TRUE, corr = "cailliez",m="min",stand.x = TRUE,print.pco=TRUE)
  
  DF <- data.frame(SpRich = fd.PA$nbsp)
  DF <-rbind(DF, data.frame(SpRich = fd.RF$nbsp))%>%rownames_to_column("Sub.cell")
  
  
  DF1 <- data.frame(FRic = fd.PA$FRic)
  DF1 <-rbind(DF1, data.frame(FRic = fd.RF$FRic))%>%rownames_to_column("Sub.cell")%>%full_join(DF)
  
  DF2 <- data.frame(FEve = fd.PA$FEve)
  DF2 <-rbind(DF2, data.frame(FEve = fd.RF$FEve))%>%rownames_to_column("Sub.cell")%>%full_join(DF1)
  
  
  DF3 <- data.frame(FDiv = fd.PA$FDiv)
  DF3<-rbind(DF3, data.frame(FDiv = fd.RF$FDiv))%>%rownames_to_column("Sub.cell")%>%full_join(DF2)
  
  
  DF4 <- data.frame(FDis = fd.PA$FDis)
  DF4 <-rbind(DF4, data.frame(FDis = fd.RF$FDis))%>%rownames_to_column("Sub.cell")%>%full_join(DF3)
  
  #combine all the datasets
  file_atr<-rbind(atr)%>%full_join(DF4)
  
  file_atr$Treatment<-as.factor(file_atr$Treatment)
  
  #GAM for Species richness
  mod1<-gam(SpRich~Treatment+s(DEM)+s(lon)+s(lat),data = file_atr, method = "REML")
  
  summary(mod1)
  
  plot1<-plot_smooths(model = mod1,series = DEM,comparison = Treatment)+
    theme_test()+theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=1000, y=30, label= "Elevation",size = 4.5)
  plot2<-plot_smooths(model = mod1,series = lon,comparison = Treatment)+theme_test()+
    theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=76.5, y=35, label= "Longitude",size = 4.5)
  plot3<-plot_smooths(model = mod1,series = lat,comparison = Treatment)+theme_test()+ 
    theme(legend.position = c(0.2, 0.8))+theme(legend.title = element_blank())+
    xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=10.5, y=55, label= "Latitude",size = 4.5)
  
  grid.arrange(plot1, plot2,plot3, nrow=3)
  
  #GAM for Functional richness
  mod1<-gam(FRic~Treatment+s(DEM)+s(lon)+s(lat),data = file_atr, method = "REML")
  
  summary(mod1)
  
  plot1<-plot_smooths(model = mod1,series = DEM,comparison = Treatment)+
    theme_test()+theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=1000, y=0.30, label= "Elevation",size = 4.5)
  plot2<-plot_smooths(model = mod1,series = lon,comparison = Treatment)+theme_test()+
    theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=76.5, y=0.35, label= "Longitude",size = 4.5)
  plot3<-plot_smooths(model = mod1,series = lat,comparison = Treatment)+theme_test()+ 
    theme(legend.position = c(0.2, 0.8))+theme(legend.title = element_blank())+
    xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=10.5, y=0.3, label= "Latitude",size = 4.5)
  
  grid.arrange(plot1, plot2,plot3, nrow=3)
  
  #GAM for Functional eveness
  mod1<-gam(FEve~Treatment+s(DEM)+s(lon)+s(lat),data = file_atr, method = "REML")
  
  summary(mod1)
  
  plot1<-plot_smooths(model = mod1,series = DEM,comparison = Treatment)+
    theme_test()+theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=1000, y=0.865, label= "Elevation",size = 4.5)
  plot2<-plot_smooths(model = mod1,series = lon,comparison = Treatment)+theme_test()+
    theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=76.5, y=0.87, label= "Longitude",size = 4.5)
  plot3<-plot_smooths(model = mod1,series = lat,comparison = Treatment)+theme_test()+ 
    theme(legend.position = c(0.2, 0.8))+theme(legend.title = element_blank())+
    xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=10.5, y=0.88, label= "Latitude",size = 4.5)
  
  grid.arrange(plot1, plot2,plot3, nrow=3)
  
  
  #GAM for Functional Divergence
  mod1<-gam(FDiv~Treatment+s(DEM)+s(lon)+s(lat),data = file_atr, method = "REML")
  
  summary(mod1)
  
  plot1<-plot_smooths(model = mod1,series = DEM,comparison = Treatment)+
    theme_test()+theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=1000, y=0.77, label= "Elevation",size = 4.5)
  plot2<-plot_smooths(model = mod1,series = lon,comparison = Treatment)+theme_test()+
    theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=76.5, y=0.79, label= "Longitude",size = 4.5)
  plot3<-plot_smooths(model = mod1,series = lat,comparison = Treatment)+theme_test()+ 
    theme(legend.position = c(0.2, 0.8))+theme(legend.title = element_blank())+
    xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=10.5, y=0.82, label= "Latitude",size = 4.5)
  
  grid.arrange(plot1, plot2,plot3, nrow=3)
  
  #GAM for Functional Dispersion
  mod1<-gam(FDis~Treatment+s(DEM)+s(lon)+s(lat),data = file_atr, method = "REML")
  
  summary(mod1)
  
  plot1<-plot_smooths(model = mod1,series = DEM,comparison = Treatment)+
    theme_test()+theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=1000, y=0.114, label= "Elevation",size = 4.5)
  plot2<-plot_smooths(model = mod1,series = lon,comparison = Treatment)+theme_test()+
    theme(legend.position = 'null')+xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=76.5, y=0.11, label= "Longitude",size = 4.5)
  plot3<-plot_smooths(model = mod1,series = lat,comparison = Treatment)+theme_test()+ 
    theme(legend.position = c(0.2, 0.8))+theme(legend.title = element_blank())+
    xlab("")+
    theme(plot.margin=unit(c(0.01,0.01,0,0), "cm"))+
    theme(axis.text = element_text(color = "black"))+
    annotate("text", x=10.5, y=0.119, label= "Latitude",size = 4.5)
  
  grid.arrange(plot1, plot2,plot3, nrow=3)
  
}


##Figure 5
#Estimate differences in species abundance in PAs and RFs

{
#Prepare the input file
RF_2<-RF.kba%>%group_by(Sub.cell)%>%mutate(list=n_distinct(List.ID))%>%
  group_by(Sub.cell,Common.Name)%>%summarise(freq=n()/list)%>%unique()%>%
  mutate(treatment="RF")%>%rename(replicate=Sub.cell)
PA_2<-PA.kba%>%group_by(Sub.cell)%>%mutate(list=n_distinct(List.ID))%>%
  group_by(Sub.cell,Common.Name)%>%summarise(freq=n()/list)%>%unique()%>%
  mutate(treatment="PA")%>%rename(replicate=Sub.cell)

file_2<-bind_rows(RF_2,PA_2)

Abund.diff<-abundance_difference(df = file_2,species.var = "Common.Name",pool=T,
                                 abundance.var = "freq",treatment.var = "treatment",
                                 replicate.var="replicate")

#get the trophic guild, family, conservation status and habitat data
foragingstrate<-readRDS("kba_species.details.rds")%>%select(Common.Name, ForagingStrata)%>%unique()
feedingguild<-readRDS("kba_species.details.rds")%>%select(Common.Name, TrophicNiche)%>%unique()

Traits<-readRDS("kba_species.details.rds")%>%select(Common.Name,IUCN.Redlist.Status, SoIB.status, Distribution.Status)%>%
  unique()%>% mutate(Common.Name=replace(Common.Name,Common.Name == "Green Warbler","Greenish Warbler"),
                     IUCN.Redlist.Status=replace(IUCN.Redlist.Status,IUCN.Redlist.Status == "Critically Endangered","Threatened"),
                     IUCN.Redlist.Status=replace(IUCN.Redlist.Status,IUCN.Redlist.Status == "Endangered","Threatened"),
                     IUCN.Redlist.Status=replace(IUCN.Redlist.Status,IUCN.Redlist.Status == "Vulnerable","Threatened"),
                     Distribution.Status=replace(Distribution.Status,Distribution.Status == "Distribution Very Restricted","Distribution Restricted"))

T1<-Traits%>%select(-c(IUCN.Redlist.Status,Distribution.Status))%>%rename(trait=SoIB.status)
T2<-Traits%>%select(-c(SoIB.status,Distribution.Status))%>%rename(trait=IUCN.Redlist.Status)
T3<-Traits%>%select(-c(SoIB.status,IUCN.Redlist.Status))%>%rename(trait=Distribution.Status)
Traits<-rbind(T1,T2)%>%rbind(T3)
remove(T1,T2,T3)

file<-Abund.diff%>%group_by(Common.Name)%>%summarise(diff=mean(difference))%>%
  left_join(feedingguild)%>%left_join(Traits)%>%left_join(foragingstrate)%>%unique()
  
Diet<-file%>%group_by(TrophicNiche)%>%summarise(sum=mean(diff))%>%filter(!is.na(TrophicNiche))
IUCN<-file%>%group_by(trait)%>%summarise(sum=mean(diff))%>%filter(!is.na(trait))
forage<-file%>%group_by(ForagingStrata)%>%summarise(sum=mean(diff))%>%filter(!is.na(ForagingStrata))

colour <- c("gray50", "gray10")
  scale_fill_discrete <- function(...) {scale_fill_manual(..., values = colour)}
  
  plot1<-ggplot(Diet) +  geom_hline(yintercept = 0,lwd = 0.2,col = "gray") +
    geom_col(data = Diet,aes(x = TrophicNiche, y = sum,fill  = sum< 0), width = 0.2) +
    theme_test()+ coord_flip()+theme(legend.position = "None") +
    labs(x = "", y = "")+ geom_text(aes(x = TrophicNiche, y = 0, label = TrophicNiche),vjust = -1)+
    theme(axis.ticks.y = element_blank(), axis.text.x = element_text(color = "black"))+
    theme(axis.text.y = element_blank())+
    ggtitle("Feeding guild")+  theme(plot.title = element_text(hjust = 0.5))+
    ylim(-0.021,0.021)+theme(plot.margin=unit(c(0.1,0.03,0.1,0.03), "cm")) 
  
  IUCN$trait <- factor(IUCN$trait, levels = c("Threatened", "Near Threatened", "Least Concern",
                                              "SoIB High", "SoIB Moderate", "SoIB Low",
                                              "Distribution Very Large", "Distribution Large", "Distribution Moderate","Distribution Restricted"))
  
  plot2<-ggplot(IUCN) +  geom_hline(yintercept = 0,lwd = 0.2,col = "grey") +
    geom_col(data = IUCN,aes(x = factor(trait), y = sum,fill  = sum< 0), width = 0.2) +
    theme_test() + coord_flip()+theme(legend.position = "None") +
    labs(x = "", y = "")+ geom_text(aes(x = trait, y = 0, label = trait),vjust = -1)+
    theme(axis.ticks.y = element_blank(), axis.text.x = element_text(color = "black"))+
    theme(axis.text.y = element_blank())+
    ggtitle("Concern categories")+  theme(plot.title = element_text(hjust = 0.5))+
    geom_vline(xintercept=6.55, linetype="dashed", color = "black")+
    geom_vline(xintercept=3.6, linetype="dashed", color = "black")+
    ylim(-0.015,0.015)+theme(plot.margin=unit(c(0.1,0.03,0.1,0.03), "cm"))
  
  plot3<-ggplot(forage) +  geom_hline(yintercept = 0,lwd = 0.2,col = "grey") +
    geom_col(data = forage,aes(x = factor(ForagingStrata), y = sum,fill  = sum< 0), width = 0.2) +
    theme_test() + coord_flip()+theme(legend.position = "None") +
    labs(x = "", y = "")+ geom_text(aes(x = ForagingStrata, y = 0, label = ForagingStrata),vjust = -1)+
    theme(axis.ticks.y = element_blank(), axis.text.x = element_text(color = "black"))+
    theme(axis.text.y = element_blank())+
    ggtitle("Foraging strata")+  theme(plot.title = element_text(hjust = 0.5))+
    ylim(-0.006,0.006)+theme(plot.margin=unit(c(0.1,0.03,0.1,0.03), "cm"))
  
  grid.arrange(plot1,plot2, plot3, ncol=3)

}

#Figure S2
##For joint plotting family-wise differences between PA and RF

{
  fam<-readRDS("kba_species.details.rds")%>%select(Common.Name,Family)%>%unique()
  family<-Abund.diff%>%left_join(fam)%>%filter(!is.na(Family))%>%
    group_by(Family)%>%summarise(sum=round(mean(difference),digits = 3))
  
  A<-family[grepl("A|B|C|D|E",family$Family),] #22
  B<-family[grepl(c("F|G|H|I|J|K|L|O|P"),family$Family),]  #22
  C<-family[grepl(c("M|N|Q|R|S|T|U|V|W|X|Y|Z"),family$Family),] #22
  
  colour <- c("gray50", "gray10")
  scale_fill_discrete <- function(...) {scale_fill_manual(..., values = colour)}
  
  A1<-ggplot(A) +  geom_hline(yintercept = 0,lwd = 0.2,col = "grey") +
    geom_col(data = A,aes(x = reorder(Family, desc(sum)), y = sum,fill  = sum< 0), width = 0.2) +
    theme_test() + coord_flip()+theme(legend.position = "None") +
    labs(x = "", y = "")+theme(plot.margin=unit(c(0.1,0.3,0.1,0.3), "cm")) 
  B1<-ggplot(B) +  geom_hline(yintercept = 0,lwd = 0.2,col = "grey") +
    geom_col(data = B,aes(x = reorder(Family, desc(sum)), y = sum,fill  = sum< 0), width = 0.2) +
    theme_test() + coord_flip()+theme(legend.position = "None") +
    labs(x = "", y = "")+theme(plot.margin=unit(c(0.1,0.3,0.1,0.3), "cm")) 
  C1<-ggplot(C) +  geom_hline(yintercept = 0,lwd = 0.2,col = "grey") +
    geom_col(data = C,aes(x = reorder(Family, desc(sum)), y = sum,fill  = sum< 0), width = 0.2) +
    theme_test() + coord_flip()+theme(legend.position = "None") +
    labs(x = "", y = "")+theme(plot.margin=unit(c(0.1,0.3,0.1,0.3), "cm")) 
  grid.arrange(A1,B1,C1, ncol=3)
  
}

#Create a dataframe and fill it with details of PA and RF

{
  
  Compare = data.frame(matrix(vector(), ncol =  3,dimnames=list(c(),c("Attribute","PA","RF"))),stringsAsFactors=F)
  
  PA.kba<-KBA%>%dplyr::filter(Season=="Dry")%>%
    group_by(Common.Name)%>%mutate(count=n())%>%filter(!count==1)%>%right_join(PA)%>%unique()%>%
    mutate(Common.Name=replace(Common.Name,Common.Name == "Green Warbler","Greenish Warbler"))%>%
    unique()%>%left_join(readRDS("kba_species.details.rds"))
  
  RF.kba<-KBA%>%dplyr::filter(Season=="Dry")%>%
    group_by(Common.Name)%>%mutate(count=n())%>%filter(!count==1)%>%right_join(RF)%>%unique()%>%
    mutate(Common.Name=replace(Common.Name,Common.Name == "Green Warbler","Greenish Warbler"))%>%
    unique()%>%left_join(readRDS("kba_species.details.rds"))
  
  
  Compare[1,'Attribute']="Checklists"
  Compare[1,'PA']=n_distinct(PA.kba$List.ID)
  Compare[1,'RF']=n_distinct(RF.kba$List.ID)
  
  
  Compare[2,'Attribute']="Total Species count"
  Compare[2,'PA']=n_distinct(PA.kba$Common.Name)
  Compare[2,'RF']=n_distinct(RF.kba$Common.Name)
  
  Compare[3,'Attribute']="Total Family count"
  Compare[3,'PA']=n_distinct(PA.kba$Family)
  Compare[3,'RF']=n_distinct(RF.kba$Family)
  
  Compare[4,'Attribute']="Diversity 'inv-simpson'"
  Compare[4,'PA']=round(mean(diversity(PA_1, "invsimpson")),digits=2)
  Compare[4,'RF']=round(mean(diversity(RF_1, "invsimpson")),digits=2)
  
  Compare[5,'Attribute']="Diversity 'shannon'"
  Compare[5,'PA']=round(mean(diversity(PA_1, "shannon")),digits=2)
  Compare[5,'RF']=round(mean(diversity(RF_1, "shannon")),digits=2)
  
  Compare[6,'Attribute']="Mean species count"
  Compare[6,'PA']=round(mean(specnumber(PA_1)),digits=2)
  Compare[6,'RF']=round(mean(specnumber(RF_1)),digits=2)
  
  Compare[7,'Attribute']="Pielou's evenness"
  Compare[7,'PA']=round(mean(mean(diversity(PA_1))/mean(log(specnumber(PA_1)))),digits=2)
  Compare[7,'RF']=round(mean(mean(diversity(RF_1))/mean(log(specnumber(RF_1)))),digits=2)
  
  Compare[8,'Attribute']="IUCN 'Least Concern'"
  Compare[8,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$IUCN.Redlist.Status=="Least Concern"])
  Compare[8,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$IUCN.Redlist.Status=="Least Concern"])
  
  Compare[9,'Attribute']="IUCN 'Threatened/Near Threatened'"
  Compare[9,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$IUCN.Redlist.Status=="Critically Endangered"|
                                                   PA.kba$IUCN.Redlist.Status=="Endangered"|
                                                   PA.kba$IUCN.Redlist.Status=="Vulnerable"|
                                                   PA.kba$IUCN.Redlist.Status=="Near Threatened"])
  Compare[9,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$IUCN.Redlist.Status=="Critically Endangered"|
                                                   RF.kba$IUCN.Redlist.Status=="Endangered"|
                                                   RF.kba$IUCN.Redlist.Status=="Vulnerable"|
                                                   RF.kba$IUCN.Redlist.Status=="Near Threatened"])
  
  Compare[10,'Attribute']="Western Ghats Endemic"
  Compare[10,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$Endemicity=="Western Ghats"])
  Compare[10,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$Endemicity=="Western Ghats"])
  
  Compare[11,'Attribute']="SoIB status 'High'"
  Compare[11,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$SoIB.status=="SoIB High"])
  Compare[11,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$SoIB.status=="SoIB High"])
  
  Compare[12,'Attribute']="SoIB status 'Moderate'"
  Compare[12,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$SoIB.status=="SoIB Moderate"])
  Compare[12,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$SoIB.status=="SoIB Moderate"])
  
  Compare[13,'Attribute']="SoIB status 'Low'"
  Compare[13,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$SoIB.status=="SoIB Low"])
  Compare[13,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$SoIB.status=="SoIB Low"])
  
  Compare[14,'Attribute']="Distribution 'Large/Very Large'"
  Compare[14,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$Distribution.Status=="Distribution Very Large"|PA.kba$Distribution.Status=="Distribution Large"])
  Compare[14,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$Distribution.Status=="Distribution Very Large"|RF.kba$Distribution.Status=="Distribution Large"])
  
  Compare[15,'Attribute']="Distribution 'Moderate'"
  Compare[15,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$Distribution.Status=="Distribution Moderate"])
  Compare[15,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$Distribution.Status=="Distribution Moderate"])
  
  Compare[16,'Attribute']="Distribution 'Restricted/Very Restricted'"
  Compare[16,'PA']=n_distinct(PA.kba$Common.Name[PA.kba$Distribution.Status=="Distribution Very Restricted"|PA.kba$Distribution.Status=="Distribution Restricted"])
  Compare[16,'RF']=n_distinct(RF.kba$Common.Name[RF.kba$Distribution.Status=="Distribution Very Restricted"|RF.kba$Distribution.Status=="Distribution Restricted"])
  
  Compare$RF<- sub(".000", "", Compare$RF)
  Compare$PA<- sub(".000", "", Compare$PA)
}
write.csv(Compare,"Compare.csv",row.names = F)

##Identify Indicator species
{
  file<-RF_1%>%dplyr::bind_rows(PA_1)%>%replace(is.na(.), 0)
  indval = multipatt(file, groups,func = "IndVal.g",duleg=TRUE,control = how(nperm=999))
  summary(indval, indvalcomp=TRUE)
}  
  



