#Sara Hamimlton
#Feb 8, 2021
# Remnant/Persistant Population Analysis Figures
# Pycnopodia Project

library(tidyverse)
library(ggplot2)
########################################################################################################
# Measuring persistance of Pycnos in grid cells with annual survey data
########################################################################################################
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataandScripts/Data_ForGit/")
RemPop_02182021 <- read_csv("MasterPycno_02182021_SpatialJoin.csv")
post <- RemPop_02182021[RemPop_02182021$event_year>=2017,]
pre <- RemPop_02182021[RemPop_02182021$event_year>=2009 & RemPop_02182021$event_year<=2012,]
#NOTE: all the data is marked as shallow (no deep). We are only using shallow data for the persistant pop analysis.

# 1) the total number of Grid cells in each region
# 2) the total number that have data in all four years
# 3) the % of surveys in each 3-4 year grid cell that saw Pycnos 

#1) the total number of Grid cells in each region
Region =  c("WestGulfofAlaska","EastGulfofAlaska","NorthernCalifornia","BritishColumbia","Baja",
            "SalishSea","SoutheastAlaska","CentralCalifornia","Aleutians","Oregon", "SouthernCalifornia", "Washington")
Phase = rep("Pre",12)
data_pre = cbind(Region,Phase )
Phase = rep("Post",12)
data_post = cbind(Region, Phase)

#capture all cells in pre and post time frame that have data (without filtering for those with 3+ years of data)
#This is then exported so that it can be plotted on a map as background to give the viewer an idea of where data is availble. 
cells_toexport <- pre %>%
  distinct(GRID_ID, .keep_all = TRUE)%>%
  dplyr::select(region,GRID_ID)
write.csv(cells_toexport, "RemnantPopulations_PreOutbreak_TotalCells.csv")
cells_toexport <- post %>%
  distinct(GRID_ID, .keep_all = TRUE)%>%
  dplyr::select(region,GRID_ID)
write.csv(cells_toexport, "RemnantPopulations_PostOutbreak_TotalCells.csv")

#2) the total number that have data in all four years
pre_34years <- pre%>%
  group_by(GRID_ID)%>% 
  distinct(event_year, .keep_all = TRUE)%>%
  add_count()%>% 
  rename(numyears= n) %>%
  filter(numyears %in% c(3,4))
pre_34numcells <- pre_34years %>%
  group_by(region) %>%
  distinct(GRID_ID) %>%
  add_count() %>%
  rename(UniqueGrids = n) %>%
  distinct(region, .keep_all = TRUE)
data_pre <- merge(data_pre, pre_34numcells, by.x = "Region", by.y = "region", all.x = TRUE)
data_pre <- data_pre[,-3]
colnames(data_pre) <- c("Region","Phase","Cellswith3YearsData")

post_34years <- post%>%
  group_by(GRID_ID)%>% 
  distinct(event_year, .keep_all = TRUE) %>%
  add_count()%>% 
  rename(numyears= n) %>%
  filter(numyears%in% c(3,4))
post_34numcells <- post_34years %>%
  group_by(region) %>%
  distinct(GRID_ID) %>%
  add_count() %>%
  rename(UniqueGrids = n) %>%
  distinct(region, .keep_all = TRUE)
data_post<- merge(data_post,post_34numcells, by = "Region",  by.y = "region",all.x =TRUE)
data_post <- data_post[,-3]
colnames(data_post) <- c("Region","Phase","Cellswith3YearsData")

#3)the % of surveys in each 3-4 year grid cell that saw Pycnos 
pre_grids = unique(pre_34years$GRID_ID)
pre_percents <-pre %>%
  filter(GRID_ID %in% pre_grids) %>%
  group_by(GRID_ID) %>%
  mutate(npres = length(GRID_ID[pres_abs==1]),
         nabs = length(GRID_ID[pres_abs==0]))%>%
  mutate(percent_present = npres/(nabs+npres))
#Export this to be used in making maps 
pre_toexport <- pre_percents %>%
  distinct(GRID_ID, .keep_all = TRUE)%>%
  select(region, GRID_ID, npres, nabs, percent_present)
#Now going to transform the % of surveys that saw Pycnos into categorical information so that it can be more easily 
#plotted as a bar chart. Divisions are 0-0.05 = Absent, 0.5-0.5 = Patchy, 0.5 - 0.9 = Common, 0.9 - 1 
categories <- vector()
for(i in 1:nrow(pre_toexport)){
  if (pre_toexport[i,5] <= 0){
    categories[i] <- "Absent" 
  } else if (pre_toexport[i,5]>0 & pre_toexport[i,5]<0.3){
    categories[i] <- "Rare" 
  } else if (pre_toexport[i,5] >= 0.3 & pre_toexport[i,5]<0.9){
    categories[i] <- "Common" 
  }else if (pre_toexport[i,5] >= 0.9){
    categories[i] <- "Very Common" 
  }
}
pre_tograph <- data.frame(pre_toexport[,c(1,2,5)],categories)
colnames(pre_tograph)<- c("Region","Grid_ID","Percent_present","Categories")
pre_tograph$Region <- factor(pre_tograph$Region,levels = rev(c("EastGulfofAlaska","SoutheastAlaska", "BritishColumbia", "SalishSea","Washington",         
                                                                     "Oregon","NorthernCalifornia","CentralCalifornia","SouthernCalifornia","Baja")))


post_grids = unique(post_34years$GRID_ID)
post_percents <-post%>%
  filter(GRID_ID %in% post_grids) %>%
  group_by(GRID_ID) %>%
  mutate(npres = length(GRID_ID[pres_abs==1]),
         nabs = length(GRID_ID[pres_abs==0]))%>%
  mutate(percent_present = npres/(nabs+npres))
post_toexport <- post_percents %>%
  distinct(GRID_ID, .keep_all = TRUE)%>%
  select(region, GRID_ID, npres, nabs, percent_present)
categories <- vector()
for(i in 1:nrow(post_toexport)){
  if (post_toexport[i,5] <= 0){
    categories[i] <- "Absent" 
  } else if (post_toexport[i,5]>0 & post_toexport[i,5]<0.3){
    categories[i] <- "Rare" 
  } else if (post_toexport[i,5] >= 0.3 & post_toexport[i,5]<0.9){
    categories[i] <- "Common" 
  }else if (post_toexport[i,5] >= 0.9){
    categories[i] <- "Very Common" 
  }
}
post_tograph <- data.frame(post_toexport[,c(1,2,5)],categories)
colnames(post_tograph)<- c("Region","Grid_ID","Percent_present","Categories")
#ordering levels of the factor columns so that they appear correctly on graph
post_tograph$Region <- factor(post_tograph$Region,levels = rev(c("EastGulfofAlaska","SoutheastAlaska", "BritishColumbia", "SalishSea","Washington",         
                                               "Oregon","NorthernCalifornia","CentralCalifornia","SouthernCalifornia","Baja")))
post_tograph$Categories <- factor(post_tograph$Categories, levels = c("Very Common","Common","Rare","Absent"))


#########################################################################################################
#GRAPHING
#Decided to only use post-outbreak data to plot and include in paper 
#########################################################################################################
#annotations representing the total number of cells with 3+ years of data in each region. Drawn from data_post 
annotations = rev(c("N = 2","N = 5", "N = 29", "N = 124", "N = 6", "N = 2","N = 4","N = 12","N = 43","N = 11"))
ggplot(post_tograph, aes(x = Region, fill = Categories))+
  geom_bar(position = "fill", color = "gray", stat = "count")+
  scale_fill_manual(values = c("#f9d926","#eb7556", "#9524a1", "#3d1d98"), name = "")+
  scale_x_discrete(labels = rev(c("East Gulf of Alaska","Southeast Alaska", "British Columbia*", "Salish Sea","Washington outer coast*",         
                                  "Oregon","Northern California","Central California","Southern California","Baja")))+
  scale_y_continuous(name = "Frequency", limits = c(0,1))+
  coord_flip()+
  theme(axis.title = element_text(size = 24), axis.text = element_text(size=20),axis.title.y = element_blank(), 
        legend.text = element_text(size =20), legend.position = "top")+
  guides(fill = guide_legend(reverse=T))+
  annotate("text", x = 1:10, y = 0.1, label = annotations, size = 6, color = "white")


