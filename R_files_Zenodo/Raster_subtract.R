#Load required packages
library(raster)
library(dplyr)
library(dismo)
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(wallace)
library(ggplot2)
library(remotes)
library(wallace)
library(rgdal)
library(blockCV)
library(sf)
library(dplyr)
library(crs)
library(ggplot2)
library(usdm)
library(reshape2)
library(car)
library(BSDA)
library(plyr)
library(hrbrthemes)
library(grid)
library(svglite)
library(DataScienceR)
library(Cairo)
library(egg)

#Using the permutation csv files from your analyses, create a bargraph of permutation importance
#Read in your dynamic model permutation results 
permute_dynamic<-read.csv("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/permute_list_dynamic.csv", header = T, sep = ",")

#Subset datasheet to the variables involved, permutation importance, and model complexity (L and LQH)
permute_dynamic_subset<-as.data.frame(cbind(permute_dynamic$L.fc.L_rm.1.variable, permute_dynamic$L.fc.L_rm.1.permutation.importance, permute_dynamic$LQH.fc.LQH_rm.1.permutation.importance))
colnames(permute_dynamic_subset) <- c("enviromental.variable", "Simple Model (L)", " Complex Model (LQH)")

#Melt dataframe together for ggplot and make permutation importance a numeric value
permute_dynamic_subset<-melt(permute_dynamic_subset, value.name = c("permutation.importance"), id = 'enviromental.variable')
permute_dynamic_subset$permutation.importance<-as.numeric(permute_dynamic_subset$permutation.importance)

#Generate ggplot of permutation importance for the dynamic model
g <-ggplot(permute_dynamic_subset, aes(x=reorder(enviromental.variable, -permutation.importance), permutation.importance,  col = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  scale_fill_viridis_d() + 
  coord_flip() +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.title.y = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.text.y  = element_text(size = 20, color = "black", family = "Times New Roman"),
        axis.text.x  = element_text(size = 20, color = "black", family = "Times New Roman"),
        legend.text = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) +
  labs(fill = "Predicted Suitability") + theme (legend.title = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) + 
  ylab("Permutation Importance") + 
  xlab("Enviromental Variable") +
  facet_grid(.~variable) + 
  theme(
    strip.text.x = element_text(
      face="bold", size=20, color = "black", family = "Times New Roman"
    ),
    strip.text.y = element_text(
      face="bold", size=20, color = "black", family = "Times New Roman"
    )
  )

setwd("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Figures")
ggsave(file="Figure_5.svg", plot=g, width=11, height=6, device = grDevices::svg)

#Using the permutation csv files from your analyses, create a bargraph of permutation importance for your static models
permute_linear<-read.csv("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/permute_list.csv", header = T, sep = ",")
permute_linear_subset<-as.data.frame(cbind(permute_linear$fc.L_rm.1.variable, permute_linear$fc.L_rm.1.permutation.importance, permute_linear$fc.LQH_rm.1.permutation.importance, permute_linear$i))
colnames(permute_linear_subset) <- c("enviromental.variable", "permutation.importance_L", "permutation.importance_LQH", "year")
permute_linear_subset$permutation.importance_L<-as.numeric(permute_linear_subset$permutation.importance_L)
permute_linear_subset$permutation.importance_LQH<-as.numeric(permute_linear_subset$permutation.importance_LQH)

g<-ggplot(permute_linear_subset, aes(x=reorder(enviromental.variable, -permutation.importance_L), permutation.importance_L,  col = year, fill = year)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  scale_fill_viridis_d() + 
  coord_flip() +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.title.y = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.text.y  = element_text(size = 8, color = "black", family = "Times New Roman"),
        axis.text.x  = element_text(size = 8, color = "black", family = "Times New Roman"),
        legend.text = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) +
  labs(fill = "Predicted Suitability") + theme (legend.title = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) + 
  ylab("Permutation Importance") + 
  xlab("Enviromental Variable")  +
  facet_grid(.~year) + 
  theme(
    strip.text.x = element_text(
      face="bold", size = 20, color = "black", family = "Times New Roman"
    ),
    strip.text.y = element_text(
      face="bold", size = 20, color = "black", family = "Times New Roman"
    )
  )


setwd("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Figures")
ggsave(file="Figure_S4.svg", plot=g, width=11, height=6, device = grDevices::svg)

g<-ggplot(permute_linear_subset, aes(x=reorder(enviromental.variable, -permutation.importance_LQH), permutation.importance_LQH,  col = year, fill = year)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  scale_fill_viridis_d() + 
  coord_flip() +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.title.y = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.text.y  = element_text(size = 8, color = "black", family = "Times New Roman"),
        axis.text.x  = element_text(size = 8, color = "black", family = "Times New Roman"),
        legend.text = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) +
  labs(fill = "Predicted Suitability") + theme (legend.title = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) + 
  ylab("Permutation Importance") + 
  xlab("Enviromental Variable")  +
  facet_grid(.~year) + 
  theme(
    strip.text.x = element_text(
      face="bold", size = 20, color = "black", family = "Times New Roman"
    ),
    strip.text.y = element_text(
      face="bold", size = 20, color = "black", family = "Times New Roman"
    )
  )

setwd("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Figures")
ggsave(file="Figure_S5.svg", plot=g, width=11, height=6, device = grDevices::svg)

yearly<-stack(list.files("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Yearly_binary/L", pattern = ".tif", full.names = TRUE, recursive = TRUE))
dynamic<-stack(list.files("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Dynamic_binary/L", pattern = ".tif", full.names = TRUE, recursive = TRUE))
dynamic_subset<-stack(dynamic$X2008, dynamic$X2009, dynamic$X2010, dynamic$X2011, dynamic$X2012, dynamic$X2013, dynamic$X2014, dynamic$X2015, dynamic$X2017, dynamic$X2018, dynamic$X2019, dynamic$X2020)
subtract<-yearly-dynamic_subset
writeRaster(subtract, file.path("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Subtractions/L", names(yearly)), format = "GTiff", bylayer = T, overwrite = T)

#Stack your  raster predictions for your simple and complex dynamic and static models, subtract them from eachother, then save them
yearly<-stack(list.files("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Yearly_binary/LQH", pattern = ".tif", full.names = TRUE, recursive = TRUE))
dynamic<-stack(list.files("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Dynamic_binary/LQH", pattern = ".tif", full.names = TRUE, recursive = TRUE))
dynamic_subset<-stack(dynamic$X2008, dynamic$X2009, dynamic$X2010, dynamic$X2011, dynamic$X2012, dynamic$X2013, dynamic$X2014, dynamic$X2015, dynamic$X2017, dynamic$X2018, dynamic$X2019, dynamic$X2020)
subtract<-yearly-dynamic_subset
plot(subtract)
writeRaster(subtract, file.path("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Subtractions/LQH", names(yearly)), format = "GTiff", bylayer = T, overwrite = T)

#pull out the extracted background point environmental values for your top three environmental variables of highest permutation importance for 2009, 2019 and your dynamic model
#In the case of our analysis, it was mean percent tree cover, mean evapotranspiration, and mean percent non-vegetated cover
#here is an example of just mean evapotranspiration
dynamic_evapotranspiration<-bg.all.mean$evapotranspiration
bg_2009<-filter(bg.match, year == 2009)
bg_2019<-filter(bg.match, year == 2019)
bg_2009<-bg_2009$evapotranspiration
bg_2019<-bg_2019$evapotranspiration
evapotranspiration<-cbind(dynamic_evapotranspiration, bg_2009, bg_2019)
evapotranspiration<-as.data.frame(evapotranspiration)
evapotranspiration <- evapotranspiration[complete.cases(evapotranspiration), ]

#perform a kruskal-wallis test on your three distributions
ks.test(evapotranspiration$dynamic_evapotranspiration, evapotranspiration$bg_2019)

#Plot your distributions of all three of your variables, Remember to swap out a, b, and c for each variable used
evapotranspiration<-as.data.frame(melt(evapotranspiration))
mu <- ddply(evapotranspiration, "variable", summarise, grp.mean=mean(value))

a <-ggplot(evapotranspiration, aes(x=value, fill=variable)) +
  geom_density(alpha=.5) + 
  scale_fill_viridis_d() + 
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed") + 
  theme(axis.title.x = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.title.y = element_text(face="bold", size = 20, color = "black", family = "Times New Roman"), 
        axis.text.y  = element_text(size = 20, color = "black", family = "Times New Roman"),
        axis.text.x  = element_text(size = 20, color = "black", family = "Times New Roman"),
        legend.text = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) +
  labs(fill = "Background Distributions") + theme (legend.title = element_text(face="bold",size = 20, color = "black", family = "Times New Roman")) + 
  ylab("Density") + 
  #xlim(0,100) + 
  facet_wrap(.~variable) + 
  xlab("Percent Tree Cover") +
  theme(
    strip.text.x = element_text(
      face="bold", size = 20, color = "black", family = "Times New Roman"
    ),
    strip.text.y = element_text(
      face="bold", size = 20, color = "black", family = "Times New Roman"
    )
  )


figure <- ggarrange(a, b, c, labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)

setwd("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Figures")
ggsave(file="Figure_S6.svg", width=10, height=10, plot=figure)
 