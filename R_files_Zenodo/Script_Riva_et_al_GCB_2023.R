### Riva et al. 2023 - GCB 

#################################  
# OPEN R PACKAGES
#################################  

library("ecospat")
library("ENMTools")
library("rJava")
library("virtualspecies")
library("dismo")
library("spatialEco")
library("dplyr")
library("viridis")
library("maps")
library("raster")
library("sp")
library("envirem")
library("rgdal")
library("spatstat")
library("proj4")
library("SSDM")
library("RColorBrewer")
library("rgeos")
library("maptools")
library("ENMeval")
library("rangeModelMetadata")
library("enmSdm")
library('tidyr')
library('ggplot2')
library('ggthemes')
library("brms")
library("sjPlot")
library("ggspatial")
library("maps")
library("ggpubr")

#################################  
### LOAD DATA AND DATA PREPARATION
#################################  

load('PATH_TO_YOUR_DATA\\\\Data_Riva_et_al_GCB_2023.RData') # update the path to open the file .RDATA included in the repository

# vector with species names
names_ck_map1 <- sapply(ckmap1, "[[", 3) 
names_ck_map2 <- sapply(ckmap2, "[[", 3) 
names_ck_map3 <- sapply(ckmap3, "[[", 3) 
names_ck_map4 <- sapply(ckmap4, "[[", 3) 
names_ck_map5 <- sapply(ckmap5, "[[", 3) 

list_pred_t1 <- sapply(ckmap1,"[[",8)
list_pred_t2 <- sapply(ckmap2,"[[",8)
list_pred_t3 <- sapply(ckmap3,"[[",8)
list_pred_t4 <- sapply(ckmap4,"[[",8)
list_pred_t5 <- sapply(ckmap5,"[[",8)

#################################  
### ASSESSMENT OF ENVIRONMENTAL NICHE MODELS
#################################  

par(mfrow=c(1,5))

### FIGURES AND INFO IN SUPPLEMENTARY INFORMATION
# proportion of urban land cover in cells
hist(covariates_final1$urban, xlim = c(0,1))
hist(covariates_final2$urban, xlim = c(0,1))
hist(covariates_final3$urban, xlim = c(0,1))
hist(covariates_final4$urban, xlim = c(0,1))
hist(covariates_final5$urban, xlim = c(0,1))
# proportion of agricultural land cover in cells
hist(covariates_final1$agric, xlim = c(0,1))
hist(covariates_final2$agric, xlim = c(0,1))
hist(covariates_final3$agric, xlim = c(0,1))
hist(covariates_final4$agric, xlim = c(0,1))
hist(covariates_final5$agric, xlim = c(0,1))

# null models; # plot p values
p_ckmap <- function(list_ckmap){
  pval <- sapply(list_ckmap, "[[", 12)
  hist((pval), main = NULL, breaks=seq(0,1,0.05))
  return(list(mean((pval), na.rm = TRUE),sd((pval), na.rm = TRUE) ))
}
p_ckmap(ckmap1)
p_ckmap(ckmap2)
p_ckmap(ckmap3)
p_ckmap(ckmap4)
p_ckmap(ckmap5)

# count how many species per grain have p < 0.05
pval_ckmap <- function(list_ckmap){
  pval <- sapply(list_ckmap, "[[", 12)
  pval <- pval[!is.na(pval)]
  length(pval[pval<0.05])/length(pval)
}
pval_ckmap(ckmap1)
pval_ckmap(ckmap2)
pval_ckmap(ckmap3)
pval_ckmap(ckmap4)
pval_ckmap(ckmap5)

# count how many observations per species in the models
counts_ckmap <- function(list_ckmap){
  counts <- sapply(list_ckmap, "[[", 1)
  hist(log10(counts), main = NULL)
  return(list(mean((counts), na.rm = TRUE),sd((counts), na.rm = TRUE) ))
}

counts_ckmap(ckmap1)
counts_ckmap(ckmap2)
counts_ckmap(ckmap3)
counts_ckmap(ckmap4)
counts_ckmap(ckmap5)

# count how many species have at least 10 observationa
counts_ckmap_ten <- function(list_ckmap){
  counts_ten <- sapply(list_ckmap, "[[", 1)
  na_ten <- sapply(list_ckmap, "[[", 11)
  
  ten <- as.data.frame(t(rbind(counts_ten, na_ten)))
  ten <- ten[complete.cases(ten), ] 
  min(ten[,1])
  sum(ten$counts_ten > 10)/length(ten$counts_ten)
}

counts_ckmap_ten(ckmap1)
counts_ckmap_ten(ckmap2)
counts_ckmap_ten(ckmap3)
counts_ckmap_ten(ckmap4)
counts_ckmap_ten(ckmap5)

# vector with AUC
AUC_ckmap <- function(list_ckmap){
  AUC <- sapply(list_ckmap, "[[", 5)
  hist(AUC, main = NULL)
  abline(v=0.7,col="red",lwd=0.5)
  return(list(mean(AUC, na.rm = TRUE), sd(AUC, na.rm = TRUE)))
}

AUC_ckmap(ckmap1)
AUC_ckmap(ckmap2)
AUC_ckmap(ckmap3)
AUC_ckmap(ckmap4)
AUC_ckmap(ckmap5)


### RESULTS ON SPECIES RICHNESS 
names(list_pred_t1) <- NULL
names(list_pred_t2) <- NULL
names(list_pred_t3) <- NULL
names(list_pred_t4) <- NULL
names(list_pred_t5) <- NULL

list_pred_t1$fun <- sum
list_pred_t2$fun <- sum
list_pred_t3$fun <- sum
list_pred_t4$fun <- sum
list_pred_t5$fun <- sum

list_pred_t1$na.rm <- TRUE
list_pred_t2$na.rm <- TRUE
list_pred_t3$na.rm <- TRUE
list_pred_t4$na.rm <- TRUE
list_pred_t5$na.rm <- TRUE

# calculate species richness at different grains
rich1 <- do.call(mosaic, list_pred_t1)
rich2 <- do.call(mosaic, list_pred_t2)
rich3 <- do.call(mosaic, list_pred_t3)
rich4 <- do.call(mosaic, list_pred_t4)
rich5 <- do.call(mosaic, list_pred_t5)

# prepare dataframe for models
rich1_model <- as.data.frame(rasterToPoints(rich1))
rich2_model <- as.data.frame(rasterToPoints(rich2))
rich3_model <- as.data.frame(rasterToPoints(rich3))
rich4_model <- as.data.frame(rasterToPoints(rich4))
rich5_model <- as.data.frame(rasterToPoints(rich5))

cov1_model <- as.data.frame(rasterToPoints(covariates_final1))
cov2_model <- as.data.frame(rasterToPoints(covariates_final2))
cov3_model <- as.data.frame(rasterToPoints(covariates_final3))
cov4_model <- as.data.frame(rasterToPoints(covariates_final4))
cov5_model <- as.data.frame(rasterToPoints(covariates_final5))

m1 <- merge(rich1_model, cov1_model, by=c("x", "y"), all.x=T)
m1 <- m1[complete.cases(m1), ]

m2 <- merge(rich2_model, cov2_model, by=c("x", "y"), all.x=T)
m2 <- m2[complete.cases(m2), ]

m3 <- merge(rich3_model, cov3_model, by=c("x", "y"), all.x=T)
m3 <- m3[complete.cases(m3), ]

m4 <- merge(rich4_model, cov4_model, by=c("x", "y"), all.x=T)
m4 <- m4[complete.cases(m4), ]

m5 <- merge(rich5_model, cov5_model, by=c("x", "y"), all.x=T)
m5 <- m5[complete.cases(m5), ]

#################################  
### RESULT SECTION NUMBER ONE: ANALYSIS OF SPECIES RICHNESS PATTERNS
#################################  

model_richness_1_clim <- brms::brm(layer ~ mean_t * precip, data = m1,  iter = 8000, cores = 11)
model_richness_1_climlu <- brms::brm(layer ~ mean_t * precip + agric + urban, data = m1,  iter = 8000, cores = 11)

tab_model(model_richness_1_clim)
tab_model(model_richness_1_climlu)

model_richness_2_clim <- brms::brm(layer ~ mean_t * precip, data = m2,  iter = 8000, cores = 11)
model_richness_2_climlu <- brms::brm(layer ~ mean_t * precip + agric + urban, data = m2,  iter = 8000, cores = 11)

tab_model(model_richness_2_clim)
tab_model(model_richness_2_climlu)

model_richness_3_clim <- brms::brm(layer ~ mean_t * precip, data = m3,  iter = 8000, cores = 11)
model_richness_3_climlu <- brms::brm(layer ~ mean_t * precip + agric + urban, data = m3,  iter = 8000, cores = 11)

tab_model(model_richness_3_clim)
tab_model(model_richness_3_climlu)


model_richness_4_clim <- brms::brm(layer ~ mean_t * precip, data = m4,  iter = 8000, cores = 11)
model_richness_4_climlu <- brms::brm(layer ~ mean_t * precip + agric + urban, data = m4,  iter = 8000, cores = 11)

tab_model(model_richness_4_clim)
tab_model(model_richness_4_climlu)

model_richness_5_clim <- brms::brm(layer ~ mean_t * precip, data = m5, iter = 8000, cores = 11)
model_richness_5_climlu <- brms::brm(layer ~ mean_t * precip + agric + urban, data = m5, iter = 8000, cores = 11)

tab_model(model_richness_5_clim)
tab_model(model_richness_5_climlu)

## FIGURE 3
# maps
richness_gg1 <- data.frame(rasterToPoints(rich1))
colnames(richness_gg1) <- c("x","y", "Species_richness")
richness_gg1[richness_gg1 == 0] <- NA # removing the zeroes at high latitudes

richness_gg2 <- data.frame(rasterToPoints(rich2))
colnames(richness_gg2) <- c("x","y", "Species_richness")
richness_gg2[richness_gg2 == 0] <- NA # removing the zeroes at high latitudes

richness_gg5 <- data.frame(rasterToPoints(rich5))
colnames(richness_gg5) <- c("x","y", "Species_richness")
richness_gg5[richness_gg5 == 0] <- NA # removing the zeroes at high latitudes

richness_gg4 <- data.frame(rasterToPoints(rich4))
colnames(richness_gg4) <- c("x","y", "Species_richness")
richness_gg4[richness_gg4 == 0] <- NA # removing the zeroes at high latitudes

richness_gg3 <- data.frame(rasterToPoints(rich3))
colnames(richness_gg3) <- c("x","y", "Species_richness")
richness_gg3[richness_gg3 == 0] <- NA # removing the zeroes at high latitudes

ggrich1<- ggplot(data=richness_gg1, aes(y=y, x=x)) +
  geom_tile(aes(fill=Species_richness)) + 
  annotation_scale(plot_unit = "m")+
  labs(fill="Species richness") + 
  coord_equal()+
  scale_fill_gradientn(limits=c(0, 170), na.value = 'white', colours=c("khaki1", "gold", "orange", "red", "dark red")) +
  theme_few()+
  labs(x = "", y = "")+
  theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank())

ggrich2<- ggplot(data=richness_gg2, aes(y=y, x=x)) +
  geom_tile(aes(fill=Species_richness)) + 
  annotation_scale(plot_unit = "m")+
  labs(fill="Species richness") + 
  coord_equal()+
  scale_fill_gradientn(limits=c(0, 170), na.value = 'white', colours=c("khaki1", "gold", "orange", "red", "dark red")) +
  theme_few()+
  labs(x = "", y = "")+
  theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank())

ggrich3<-ggplot(data=richness_gg3, aes(y=y, x=x)) +
  geom_tile(aes(fill=Species_richness)) + 
  annotation_scale(plot_unit = "m")+
  labs(fill="Species richness") + 
  coord_equal()+
  scale_fill_gradientn(limits=c(0, 170), na.value = 'white', colours=c("khaki1", "gold", "orange", "red", "dark red")) +
  theme_few()+
  labs(x = "", y = "")+
  theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank())

ggrich4<-ggplot(data=richness_gg4, aes(y=y, x=x)) +
  geom_tile(aes(fill=Species_richness)) + 
  annotation_scale(plot_unit = "m")+
  labs(fill="Species richness") + 
  coord_equal()+
  scale_fill_gradientn(limits=c(0, 170), na.value = 'white', colours=c("khaki1", "gold", "orange", "red", "dark red")) +
  theme_few()+
  labs(x = "", y = "")+
  theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank())

ggrich5<-ggplot(data=richness_gg5, aes(y=y, x=x)) +
  geom_tile(aes(fill=Species_richness)) + 
  annotation_scale(plot_unit = "m")+
  labs(fill="Species richness") + 
  coord_equal()+
  scale_fill_gradientn(limits=c(0, 170), na.value = 'white', colours=c("khaki1", "gold", "orange", "red", "dark red")) +
  theme_few()+
  labs(x = "", y = "")+
  theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank())

legend_ggrich<-ggplot(data=richness_gg5, aes(y=y, x=x)) +
  geom_tile(aes(fill=Species_richness)) + 
  annotation_scale(plot_unit = "m")+
  labs(fill="Species richness") + 
  coord_equal()+
  scale_fill_gradientn(limits=c(0, 170), na.value = 'white', colours=c("khaki1", "gold", "orange", "red", "dark red")) +
  theme_few()+
  labs(x = "", y = "")

ggarrange(ggrich1, ggrich2, ggrich3, ggrich4, ggrich5, get_legend(legend_ggrich),
          labels = c("3-km grain", "6-km grain", "12-km grain", "24-km grain", "48-km grain", ""),
          ncol = 6, nrow = 1, widths = c(1,1,1,1,1,0.5), hjust = -1) + bgcolor("white")
# ggsave("fig3b.jpg", path = "C:\\\\Users\\\\feder\\\\OneDrive\\\\Desktop\\\\ckmap 2021-2022\\\\figures", 
#        width = 5000, height = 1300, units = "px",  device='jpeg', dpi=300)

## R2 barplot
grain <- (c("3-km grain", "3-km grain",
            "6-km grain", "6-km grain",
            "12-km grain", "12-km grain",
            "24-km grain", "24-km grain",
            "48-km grain", "48-km grain"))
condition <- (c("Climate", "Climate and land use",
                "Climate", "Climate and land use",
                "Climate", "Climate and land use",
                "Climate", "Climate and land use",
                "Climate", "Climate and land use"))

# R2 values based on models
value <- c(0.546, 0.668,
           0.585, 0.706,
           0.642, 0.727,
           0.664, 0.695,
           0.642, 0.669)
data_r2 <- data.frame(grain,condition,value)

data_r2$grain <- factor(data_r2$grain, levels=c("3-km grain",
                                                "6-km grain",
                                                "12-km grain",
                                                "24-km grain",
                                                "48-km grain"))

ggplot(data_r2, aes(fill=condition, y=value, x=grain)) + 
  geom_bar(position="dodge", stat="identity")+ 
  theme_few()+
  theme(legend.position="top",
        legend.title=element_blank(),
        axis.title.x = element_blank())+
  labs(y = bquote(~R^2))+
  scale_fill_manual(values=c("gold","red3"))

# ggsave("fig3b.jpg", path = "C:\\\\Users\\\\feder\\\\OneDrive\\\\Desktop\\\\ckmap 2021-2022\\\\figures", 
#        width = 2500, height = 1000, units = "px",  device='jpeg', dpi=300)

#################################  
### RESULT SECTION NUMBER TWO: ANALYSIS OF SINGLE-SPECIES MODELS
#################################  

# extract covariates that affected the distribution of species 
effects_ckmap <- function(list_ckmap){
  
  betas_ck_map <- sapply(list_ckmap, "[[", 4) # select the 4th element of all lists
  betas_ck_map <- purrr::map_df(betas_ck_map, ~as.data.frame(t(.))) # transform into a table
  betas_ck_map <- betas_ck_map[ , -which(names(betas_ck_map) %in% c("V1"))] # remove the column V1
  rownames(betas_ck_map) <- sapply(ckmap5, "[[", 3) # apply the name of the species to each row
  betas_ck_map <- betas_ck_map[,order(colnames(betas_ck_map))] #order alphabetically so that different scales can be assessed consistently
  
  colSums(!is.na(betas_ck_map)) # how many times each covariate is important?
  not_modelled <- sum(rowSums(is.na(betas_ck_map)) == 22) # how many species had no maxent model fitted? (all NAs)
  
  # 22 is the number of maximun NA (11 covariates and 11 quadratic effects); 47 species didn't fit the model; sample of 241 species 
  prop_species_affected <- colSums(!is.na(betas_ck_map))/ (288 - sum(rowSums(is.na(betas_ck_map)) == 22)) # how often  is each covariate important
  
  effects <- data.frame('agri NA' = sum(is.na(betas_ck_map$agric) & is.na(betas_ck_map$`I(agric^2)`), na.rm = TRUE),
                        'agri plus' = sum(betas_ck_map$agric > 0 & is.na(betas_ck_map$`I(agric^2)`), na.rm = TRUE),
                        'agri minus' = sum(betas_ck_map$agric < 0 & is.na(betas_ck_map$`I(agric^2)`), na.rm = TRUE),
                        'agri quadratic' = sum(!is.na(betas_ck_map$`I(agric^2)`), na.rm = TRUE), 
                        
                        'bare NA' = sum(is.na(betas_ck_map$bare) & is.na(betas_ck_map$`I(bare^2)`), na.rm = TRUE),
                        'bare plus' = sum(betas_ck_map$bare > 0 & is.na(betas_ck_map$`I(bare^2)`), na.rm = TRUE),
                        'bare minus' = sum(betas_ck_map$bare < 0 & is.na(betas_ck_map$`I(bare^2)`), na.rm = TRUE),
                        'bare quadratic' = sum(!is.na(betas_ck_map$`I(bare^2)`), na.rm = TRUE), 
                        
                        'forest NA' = sum(is.na(betas_ck_map$forest) & is.na(betas_ck_map$`I(forest^2)`), na.rm = TRUE),
                        'forest plus' = sum(betas_ck_map$forest > 0 & is.na(betas_ck_map$`I(forest^2)`), na.rm = TRUE),
                        'forest minus' = sum(betas_ck_map$forest < 0 & is.na(betas_ck_map$`I(forest^2)`), na.rm = TRUE),
                        'forest quadratic' = sum(!is.na(betas_ck_map$`I(forest^2)`), na.rm = TRUE), 
                        
                        'herb NA' = sum(is.na(betas_ck_map$herb) & is.na(betas_ck_map$`I(herb^2)`), na.rm = TRUE),
                        'herb plus' = sum(betas_ck_map$herb > 0 & is.na(betas_ck_map$`I(herb^2)`), na.rm = TRUE),
                        'herb minus' = sum(betas_ck_map$herb < 0 & is.na(betas_ck_map$`I(herb^2)`), na.rm = TRUE),
                        'herb quadratic' = sum(!is.na(betas_ck_map$`I(herb^2)`), na.rm = TRUE), 
                        
                        'urban NA' = sum(is.na(betas_ck_map$urban) & is.na(betas_ck_map$`I(urban^2)`), na.rm = TRUE),
                        'urban plus' = sum(betas_ck_map$urban > 0 & is.na(betas_ck_map$`I(urban^2)`), na.rm = TRUE),
                        'urban minus' = sum(betas_ck_map$urban < 0 & is.na(betas_ck_map$`I(urban^2)`), na.rm = TRUE),
                        'urban quadratic' = sum(!is.na(betas_ck_map$`I(urban^2)`), na.rm = TRUE), 
                        
                        'wet NA' = sum(is.na(betas_ck_map$wet) & is.na(betas_ck_map$`I(wet^2)`), na.rm = TRUE),
                        'wet plus' = sum(betas_ck_map$wet > 0 & is.na(betas_ck_map$`I(wet^2)`), na.rm = TRUE),
                        'wet minus' = sum(betas_ck_map$wet < 0 & is.na(betas_ck_map$`I(wet^2)`), na.rm = TRUE),
                        'wet quadratic' = sum(!is.na(betas_ck_map$`I(wet^2)`), na.rm = TRUE), 
                        
                        'mean_t NA' = sum(is.na(betas_ck_map$mean_t) & is.na(betas_ck_map$`I(mean_t^2)`), na.rm = TRUE),
                        'mean_t plus' = sum(betas_ck_map$mean_t > 0 & is.na(betas_ck_map$`I(mean_t^2)`), na.rm = TRUE),
                        'mean_t minus' = sum(betas_ck_map$mean_t < 0 & is.na(betas_ck_map$`I(mean_t^2)`), na.rm = TRUE),
                        'mean_t quadratic' = sum(!is.na(betas_ck_map$`I(mean_t^2)`), na.rm = TRUE), 
                        
                        'seasonality NA' = sum(is.na(betas_ck_map$seasonality) & is.na(betas_ck_map$`I(seasonality^2)`), na.rm = TRUE),
                        'seasonality plus' = sum(betas_ck_map$seasonality > 0 & is.na(betas_ck_map$`I(seasonality^2)`), na.rm = TRUE),
                        'seasonality minus' = sum(betas_ck_map$seasonality < 0 & is.na(betas_ck_map$`I(seasonality^2)`), na.rm = TRUE),
                        'seasonality quadratic' = sum(!is.na(betas_ck_map$`I(seasonality^2)`), na.rm = TRUE), 
                        
                        'precip NA' = sum(is.na(betas_ck_map$precip) & is.na(betas_ck_map$`I(precip^2)`), na.rm = TRUE),
                        'precip plus' = sum(betas_ck_map$precip > 0 & is.na(betas_ck_map$`I(precip^2)`), na.rm = TRUE),
                        'precip minus' = sum(betas_ck_map$precip < 0 & is.na(betas_ck_map$`I(precip^2)`), na.rm = TRUE),
                        'precip quadratic' = sum(!is.na(betas_ck_map$`I(precip^2)`), na.rm = TRUE), 
                        
                        'precip_coldest NA' = sum(is.na(betas_ck_map$precip_coldest) & is.na(betas_ck_map$`I(precip_coldest^2)`), na.rm = TRUE),
                        'precip_coldest plus' = sum(betas_ck_map$precip_coldest > 0 & is.na(betas_ck_map$`I(precip_coldest^2)`), na.rm = TRUE),
                        'precip_coldest minus' = sum(betas_ck_map$precip_coldest < 0 & is.na(betas_ck_map$`I(precip_coldest^2)`), na.rm = TRUE),
                        'precip_coldest quadratic' = sum(!is.na(betas_ck_map$`I(precip_coldest^2)`), na.rm = TRUE), 
                        
                        'PET_wettest_quarter NA' = sum(is.na(betas_ck_map$PET_wettest_quarter) & is.na(betas_ck_map$`I(PET_wettest_quarter^2)`), na.rm = TRUE),
                        'PET_wettest_quarter plus' = sum(betas_ck_map$PET_wettest_quarter > 0 & is.na(betas_ck_map$`I(PET_wettest_quarter^2)`), na.rm = TRUE),
                        'PET_wettest_quarter minus' = sum(betas_ck_map$PET_wettest_quarter < 0 & is.na(betas_ck_map$`I(PET_wettest_quarter^2)`), na.rm = TRUE),
                        'PET_wettest_quarter quadratic' = sum(!is.na(betas_ck_map$`I(PET_wettest_quarter^2)`), na.rm = TRUE) 
  )
  
  
  
}

# correct for a few fitting errors at the coarsest scale (NAs in raster cells)
ckmap5[[23]]$model_parameters <- NA 
ckmap5[[106]]$model_parameters <- NA 
ckmap5[[149]]$model_parameters <- NA 
ckmap5[[237]]$model_parameters <- NA 


effect1 <- effects_ckmap(ckmap1)
effect2 <- effects_ckmap(ckmap2)
effect3 <- effects_ckmap(ckmap3)
effect4 <- effects_ckmap(ckmap4)
effect5 <- effects_ckmap(ckmap5)

effects_all <- rbind(effect1, effect2, effect3, effect4, effect5)

effects_list <- list()
for (i in 1 : 11){
  effects_list[[i]] <- effects_all[, ((c(1,2,3,4) + 4*i) - 4) ]
  
}

effects_list <- purrr::map_df(effects_list, ~as.data.frame(t(.)))
colnames(effects_list) <- c("3 km", "6 km", "12 km", "24 km", "48 km")


modeled_ckmap <- function(list_ckmap){
  betas_ck_map <- sapply(list_ckmap, "[[", 4) # select the 4th element of all lists
  betas_ck_map <- purrr::map_df(betas_ck_map, ~as.data.frame(t(.))) # transform into a table
  betas_ck_map <- betas_ck_map[ , -which(names(betas_ck_map) %in% c("V1"))] # remove the column V1
  rownames(betas_ck_map) <- sapply(ckmap5, "[[", 3) # apply the name of the species to each row
  betas_ck_map <- betas_ck_map[,order(colnames(betas_ck_map))] #order alphabetically so that different scales can be assessed consistently
  modelled <- data.frame('modeled' = sum(rowSums(is.na(betas_ck_map)) == 22)) # how many species had no maxent model fitted? (all NAs)
}

modelled1 <- modeled_ckmap(ckmap1)
modelled2 <- modeled_ckmap(ckmap2)
modelled3 <- modeled_ckmap(ckmap3)
modelled4 <- modeled_ckmap(ckmap4)
modelled5 <- modeled_ckmap(ckmap5)
modelled_all <- c(modelled1, modelled2, modelled3, modelled4, modelled5)
modelled_all <- unlist(modelled_all)

effects_list_prop <- effects_list
for (i in 1 : 11){
  effects_list_prop[(1 + 4*i -4),] <- effects_list_prop[(1 + 4*i -4),] - modelled_all
  
}


effects_list_prop$`3 km` <- effects_list_prop$`3 km`/(288 - as.numeric(modelled1))
effects_list_prop$`6 km` <- effects_list_prop$`6 km`/(288 - as.numeric(modelled2))
effects_list_prop$`12 km` <- effects_list_prop$`12 km`/(288 - as.numeric(modelled3))
effects_list_prop$`24 km` <- effects_list_prop$`24 km`/(288 - as.numeric(modelled4))
effects_list_prop$`48 km` <- effects_list_prop$`48 km`/(288 - as.numeric(modelled5))

modelled_all_prop <- 1 - (modelled_all/288)
effects_list_prop$`3 km` <- effects_list_prop$`3 km` * modelled_all_prop[[1]]
effects_list_prop$`6 km` <- effects_list_prop$`6 km`* modelled_all_prop[[2]]
effects_list_prop$`12 km` <- effects_list_prop$`12 km`* modelled_all_prop[[3]]
effects_list_prop$`24 km` <- effects_list_prop$`24 km`* modelled_all_prop[[4]]
effects_list_prop$`48 km` <- effects_list_prop$`48 km`* modelled_all_prop[[5]]

effects_list_prop$covariate <- c(rep('Agriculture', 4),
                                 rep('Bare', 4),
                                 rep('Forest', 4),
                                 rep('Herbaceous', 4),
                                 rep('Urban', 4),
                                 rep('Water', 4),
                                 rep('Temperature', 4),
                                 rep('Seasonality', 4),
                                 rep('Precipitation', 4),
                                 rep('Precipitation coldest', 4),
                                 rep('PET wettest quarter', 4))

effects_list_prop$covariate_feature <- c(rep(c('No effect', 'Positive effect', 'Negative effect', 'Quadratic effect'), 11))




effects_list_plot <- gather(effects_list_prop, scale, measurement, `3 km`:`48 km`, factor_key=TRUE)
effects_list_plot$covariate <- factor(effects_list_plot$covariate, levels = c("Agriculture", "Urban", "Forest", "Herbaceous", "Bare", "Water", "Temperature", "Seasonality", "Precipitation", "Precipitation coldest", "PET wettest quarter"))

effects_list_plot$covariate_feature <- factor(effects_list_plot$covariate_feature, levels = c("Positive effect","Quadratic effect", "Negative effect","No effect" ))

effects_list_plot <- effects_list_plot[effects_list_plot$covariate != "Water", ] 

# figure 4
ggplot(effects_list_plot, aes(fill=covariate_feature, y=measurement, x=scale)) + 
  geom_bar(position="stack", stat="identity") +
  ylim(0,1) +
  #scale_fill_viridis(discrete = T) +
  #ggtitle("Studying 4 species..") +
  facet_wrap(~ covariate, ncol =5)+
  xlab("")+
  ylab("Proportion of species responding to each covariate")+
  theme_base()+
  #scale_fill_manual(values= c("chartreuse","deeppink","cyan", "black"))+
  scale_fill_manual(values= c("gold","orange","red3","black"))+ 
  theme(legend.title = element_blank())+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size=rel(1.5)))

# ggsave("fig4.jpg", path = "C:\\\\Users\\\\feder\\\\OneDrive\\\\Desktop\\\\ckmap 2021-2022\\\\figures", 
#        width = 5000, height = 2500, units = "px",  device='jpeg', dpi=300)

# model proportion of "no effect" across scales
no_eff <- effects_list_plot
no_eff <- no_eff[no_eff$covariate_feature == "No effect", ]
no_eff$type <- rep(c("LULC", "LULC","LULC","LULC","LULC",
                     "Climate", "Climate","Climate","Climate","Climate"), 5)
no_eff$scale_cont <- c(3,3,3,3,3,3,3,3,3,3,
                       6,6,6,6,6,6,6,6,6,6,
                       12,12,12,12,12,12,12,12,12,12,
                       24,24,24,24,24,24,24,24,24,24,
                       48,48,48,48,48, 48,48,48,48,48)

no_eff_model <- brm(measurement ~ type * scale + (1|covariate), data = no_eff)
summary(no_eff_model)
tab_model(no_eff_model)

#################################  
### RESULT SECTION NUMBER THREE: ANALYSIS OF FUNCTIONAL TRAITS 
#################################  

traits$Taxa.name <- (sub("_", " ", traits$Taxa.name))
traits_match <- traits[traits$Taxa.name %in% names_ck_map1, ]

dev.off()
##
ggplot(traits, aes(x=Wsp_female_avg)) + 
  geom_histogram(binwidth = 5) + 
  theme_bw()+
  labs(y= "Number of species", x = "Wingspan (mm)")

ggplot(traits, aes(x=total_flight_mo)) + 
  geom_histogram(binwidth = 1) + 
  theme_bw() + 
  xlim(0,13) +
  labs(y= "Number of species", x = "Length of the flight window (months)")


plot(traits$Wsp_female_avg ~ traits$total_flight_mo)

# effects of land use on species
sp_coeff <- function(ckmap){
  sp_coeff <- sapply(ckmap, "[[", 4) # select the 4th element of all lists
  sp_coeff <- purrr::map_df(sp_coeff, ~as.data.frame(t(.))) # transform into a table
  sp_coeff <- sp_coeff[ , -which(names(sp_coeff) %in% c("V1"))] # remove the column V1
  rownames(sp_coeff) <- sapply(ckmap, "[[", 3) # apply the name of the species to each row
  return(sp_coeff)
}

sp_coeff_1 <- sp_coeff(ckmap1)
sp_coeff_2 <- sp_coeff(ckmap2)
sp_coeff_3 <- sp_coeff(ckmap3)
sp_coeff_4 <- sp_coeff(ckmap4)
sp_coeff_5 <- sp_coeff(ckmap5)

sp_coeff_1 <- sp_coeff_1[rowSums(is.na(sp_coeff_1)) != ncol(sp_coeff_1), ]
sp_coeff_2 <- sp_coeff_2[rowSums(is.na(sp_coeff_2)) != ncol(sp_coeff_2), ]
sp_coeff_3 <- sp_coeff_3[rowSums(is.na(sp_coeff_3)) != ncol(sp_coeff_3), ]
sp_coeff_4 <- sp_coeff_4[rowSums(is.na(sp_coeff_4)) != ncol(sp_coeff_4), ]
sp_coeff_5 <- sp_coeff_5[rowSums(is.na(sp_coeff_5)) != ncol(sp_coeff_5), ]

#n species
nrow(sp_coeff_1)/length(names_ck_map1)
nrow(sp_coeff_2)/length(names_ck_map1)
nrow(sp_coeff_3)/length(names_ck_map1)
nrow(sp_coeff_4)/length(names_ck_map1)
nrow(sp_coeff_5)/length(names_ck_map1)

sp_coeff_1$species <- rownames(sp_coeff_1)
sp_coeff_1$scale <- rep(3, nrow(sp_coeff_1))

sp_coeff_2$species <- rownames(sp_coeff_2)
sp_coeff_2$scale <- rep(6, nrow(sp_coeff_2))

sp_coeff_3$species <- rownames(sp_coeff_3)
sp_coeff_3$scale <- rep(12, nrow(sp_coeff_3))

sp_coeff_4$species <- rownames(sp_coeff_4)
sp_coeff_4$scale <- rep(24, nrow(sp_coeff_4))

sp_coeff_5$species <- rownames(sp_coeff_5)
sp_coeff_5$scale <- rep(48, nrow(sp_coeff_5))

sp_coeff <- rbind(sp_coeff_1, sp_coeff_2, sp_coeff_3, sp_coeff_4, sp_coeff_5)

# effects of land use
# agriculture
agricultural_effect <- list()
for(i in 1:nrow(sp_coeff)){ 
  agricultural_effect[[i]] <- if (is.na(sp_coeff[i,]$`I(agric^2)`) & is.na(sp_coeff[i,]$agric)) {
    "No effect"
  } else {
    if (!is.na(sp_coeff[i,]$`I(agric^2)`)) {
      "Quadratic effect"
    } else {
      if (sp_coeff[i,]$agric > 0) {
        "Positive effect"
      } else {
        "Negative effect"
      }
    }
  }
}

agricultural_effect <- do.call(rbind, agricultural_effect)
sp_coeff$agricultural_effect <- agricultural_effect

# urban
urban_effect <- list()
for(i in 1:nrow(sp_coeff)){ 
  urban_effect[[i]] <- if (is.na(sp_coeff[i,]$`I(urban^2)` & is.na(sp_coeff[i,]$urban))) {
    "No effect"
  } else {
    if (!is.na(sp_coeff[i,]$`I(urban^2)`)) {
      "Quadratic effect"
    } else {
      if (sp_coeff[i,]$urban > 0) {
        "Positive effect"
      } else {
        "Negative effect"
      }
    }
  }
}

urban_effect <- do.call(rbind, urban_effect)
sp_coeff$urban_effect <- urban_effect

sp_coeff_model <- sp_coeff[,23:26]
names(traits)[[1]] <- "species"
names(number_obs)[[1]] <- "species"
sp_coeff_model <- merge(sp_coeff_model, traits, all.x = TRUE, by = "species")
sp_coeff_model <- merge(sp_coeff_model, number_obs, all.x = TRUE, by = "species")
sp_coeff_model$scale2 <- as.factor(sp_coeff_model$scale)

# models of response as a function of trait
fit_agri_mo <- brm(agricultural_effect ~ total_flight_mo * scale2 + (1|Family/species), data = sp_coeff_model, 
                   family = categorical(link = "logit"), #multinomial(), 
                   inits = 0,
                   seed   = 123,
                   iter = 14000,
                   chains = 4, cores = 11)

fit_urban_mo <- brm(urban_effect ~ total_flight_mo * scale2 + (1|Family/species), data = sp_coeff_model, 
                    family = categorical(link = "logit"), #multinomial(), 
                    inits = 0,
                    seed   = 123,
                    iter = 14000,
                    chains = 4, cores = 11)


fit_agri_wi <- brm(agricultural_effect ~ Wsp_female_avg * scale2 + (1|Family/species), data = sp_coeff_model, 
                   family = categorical(link = "logit"), #multinomial(), 
                   inits = 0,
                   seed   = 123,
                   iter = 14000,
                   chains = 4, cores = 11)

fit_urban_wi <- brm(urban_effect ~ Wsp_female_avg * scale2 + (1|Family/species), data = sp_coeff_model, 
                    family = categorical(link = "logit"), #multinomial(), 
                    inits = 0,
                    seed   = 123,
                    iter = 14000,
                    chains = 4, cores = 11)


tab_model(fit_agri_wi, transform = NULL)
tab_model(fit_agri_mo, transform = NULL)
tab_model(fit_urban_wi, transform = NULL)
tab_model(fit_urban_mo, transform = NULL)

plot_models(fit_agri_wi, transform = NULL)

agri_wi <- conditional_effects(fit_agri_wi, categorical = TRUE)
agri_mo <- conditional_effects(fit_agri_mo, categorical = TRUE)
urba_wi <- conditional_effects(fit_urban_wi, categorical = TRUE)
urba_mo <- conditional_effects(fit_urban_mo, categorical = TRUE)

#plot(conditional_effects(fit_urban_mo, categorical = TRUE), prob = 0.89, ask = FALSE)

conditions <- make_conditions(fit_urban_mo, vars = "scale2")

urban_wi_plot <- conditional_effects(
  fit_urban_wi, 
  conditions = conditions,
  prob = 0.89,
  categorical = TRUE,
  spaghetti = FALSE,
  ncol = 5)

urban_mo_plot <- conditional_effects(
  fit_urban_mo, 
  conditions = conditions,
  prob = 0.89,
  categorical = TRUE,
  spaghetti = FALSE,
  ncol = 5)

agri_wi_plot <- conditional_effects(
  fit_agri_wi, 
  conditions = conditions,
  prob = 0.89,
  categorical = TRUE,
  spaghetti = FALSE,
  ncol = 5)

agri_mo_plot <- conditional_effects(
  fit_agri_mo, 
  conditions = conditions,
  prob = 0.89,
  categorical = TRUE,
  spaghetti = FALSE,
  ncol = 5)
##

levels(urban_wi_plot[["Wsp_female_avg:cats__"]]$scale2) <- c("3-km grain", "6-km grain", "12-km grain", "24-km grain", "48-km grain")
levels(urban_mo_plot[["total_flight_mo:cats__"]]$scale2) <- c("3-km grain", "6-km grain", "12-km grain", "24-km grain", "48-km grain")
levels(agri_wi_plot[["Wsp_female_avg:cats__"]]$scale2) <- c("3-km grain", "6-km grain", "12-km grain", "24-km grain", "48-km grain")
levels(agri_mo_plot[["total_flight_mo:cats__"]]$scale2) <- c("3-km grain", "6-km grain", "12-km grain", "24-km grain", "48-km grain")

p1 <- plot(urban_wi_plot, plot = FALSE)[[1]] + 
  facet_wrap("scale2", ncol = 5) +
  scale_color_manual(values=c("red3", "black", "gold","orange"))+ 
  scale_fill_manual(values=c("red3", "black", "gold","orange"))+ 
  theme_few()+
  labs(x = "Wingspan", y = "Probability of effect")+
  theme(legend.position="none", legend.title=element_blank())+
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  scale_x_continuous(name="Wingspan", limits=c(20, 80), breaks=c(25, 50, 75))


p2 <- plot(urban_mo_plot, plot = FALSE)[[1]] + 
  facet_wrap("scale2", ncol = 5) +
  scale_color_manual(values=c("red3", "black", "gold","orange"))+ 
  scale_fill_manual(values=c("red3", "black", "gold","orange"))+ 
  theme_few()+
  labs(x = "Flight window", y = "Probability of effect")+
  theme(legend.position="none", legend.title=element_blank())+
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  scale_x_continuous(name="Flight window", limits=c(1, 12), breaks=c(4,8,12))

p3 <- plot(agri_wi_plot, plot = FALSE)[[1]] + 
  facet_wrap("scale2", ncol = 5) +
  scale_color_manual(values=c("red3", "black", "gold","orange"))+ 
  scale_fill_manual(values=c("red3", "black", "gold","orange"))+ 
  theme_few()+
  labs(x = "Wingspan", y = "Probability of effect")+
  theme(legend.position="none", legend.title=element_blank())+
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  scale_x_continuous(name="Wingspan", limits=c(20, 80), breaks=c(25, 50, 75))

p4 <- plot(agri_mo_plot, plot = FALSE)[[1]] + 
  facet_wrap("scale2", ncol = 5) +
  scale_color_manual(values=c("red3", "black", "gold","orange"))+ 
  scale_fill_manual(values=c("red3", "black", "gold","orange"))+ 
  theme_few()+
  labs(x = "Flight window", y = "Probability of effect")+
  theme(legend.position="none", legend.title=element_blank())+
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  scale_x_continuous(name="Flight window", limits=c(1, 12), breaks=c(4,8,12))


ggarrange(p1, p2, p3, p4, 
          labels = c("Urban", "Urban", "Agricultural", "Agricultural"),
          common.legend = TRUE,  legend = "bottom",
          font.label = list(size = 24, color = "black", face = "bold", family = NULL)) + bgcolor("white")
# ggsave("fig5.jpg", path = "C:\\\\Users\\\\feder\\\\OneDrive\\\\Desktop\\\\ckmap 2021-2022\\\\figures", 
#        width = 5000, height = 2000, units = "px",  device='jpeg', dpi=300)

#################################  
### SUPPLEMENTARY FUNCTIONS
#################################  

# ### FIT MAXENT MODELS
# eB_sdm <- function(species_name, coord_data, raster_cov, ecorg, bias_window, denominator_pseudoabsence, dispersal_distance){
#   
#   set.seed(123)
#   # remove double-observations in each raster cell based on resolution of "raster_cov"
#   occu <- base::subset(coord_data, coord_data$nome_sp == species_name)
#   occu <- occu[,2:3]
#   occu <- raster::rasterize(occu, raster_cov[[1]], fun='count')
#   occu[occu > 0] <- 1
#   occu <- raster::rasterToPoints(occu)
#   # as matrix necessary for species with only one obs
#   # when only one unique cell is represented in the data, need an extra line of code
#   if(nrow(occu) > 1){
#     occu <- occu[, 1:2]
#   } else {
#     occu <- occu[, 1:2]
#     occu <- t(as.data.frame(occu))
#   }
#   
#   tryCatch({
#     
#     # reduce study area
#     names_raster_cov <- names(raster_cov)
#     extract_regions <- raster::extract(ecorg, occu)
#     extract_regions <- colSums(extract_regions, na.rm = TRUE)
#     extract_regions[extract_regions > 0] <- 1
#     extract_regions <- t(as.data.frame(extract_regions))
#     
#     ecorg_drop <- raster::dropLayer(ecorg, (colnames(extract_regions)[colSums(pmax(extract_regions == 0)) > 0]) )
#     ecorg_drop <- sum(ecorg_drop)
#     ecorg_drop[ecorg_drop > 0] <- 1
#     ecorg_drop[ecorg_drop == 0] <- NA
#     
#     raster_cov <- raster_cov * ecorg_drop[[1]]
#     crop_extent <- extent(min(occu[,1])- dispersal_distance, max(occu[,1]) + dispersal_distance, min(occu[,2]) - dispersal_distance, max(occu[,2]) + dispersal_distance)
#     raster_cov <- crop(raster_cov, crop_extent)
#     names(raster_cov) <- names_raster_cov
#     raster_cov <- scale(raster_cov)
#     
#     # create pseudoabsences based on density of occurrence data
#     bias <- raster::rasterize(coord_data[, 2:3], raster_cov[[1]], fun='count') #
#     bias <- raster::focal(bias, w=matrix(1, bias_window, bias_window)  , fun=mean, pad=TRUE,  na.rm = TRUE )
#     bias <- raster::crop(bias, raster_cov[[1]])
#     bias <- raster::mask(bias, raster_cov[[1]])
#     bias <- log2(bias + 1)
#     pseudo_absences <- raster::xyFromCell(bias, sample(which(!is.na(values(bias))),
#                                                        (sum((!is.na(values(ecorg_drop))), na.rm=TRUE))/denominator_pseudoabsence, # denominator_pseudoabsence determines proportion of cells available in the raster that are not occupied
#                                                        prob=values(bias)[!is.na(values(bias))]))
# 
#     #fit maxent model
#     maxent_model <- ENMeval::ENMevaluate(occu,
#                                          raster_cov,
#                                          pseudo_absences,
#                                          tune.args = list(fc = c("L","LQ"), rm = 1:5),
#                                          partitions = "block",
#                                          other.settings = list(abs.auc.diff = FALSE, pred.type = "logistic", validation.bg = "partition"),
#                                          algorithm = "maxnet",
#                                          taxon.name = species_name,
#                                          overlap = TRUE,
#                                          parallel = FALSE) 
#     
#     if(all(is.na(maxent_model@results$AICc))){
#       best_model <- (maxent_model@models[[which.max(maxent_model@results$auc.val.avg)]])
#       pred_best_model <- predict(raster_cov, maxent_model@models[[which.max(maxent_model@results$auc.val.avg)]], type = 'cloglog')
#     } else {
#       best_model <- (maxent_model@models[[which.min(maxent_model@results$AICc)]])
#       pred_best_model <- predict(raster_cov, maxent_model@models[[which.min(maxent_model@results$AICc)]], type = 'cloglog')
#     }
#     
#     #find threshold of sensitivity equalling specificity
#     p <- raster::extract(pred_best_model, occu, method = "simple")
#     a <- raster::extract(pred_best_model, pseudo_absences, method = "simple")
#     e <- dismo::evaluate(p=p, a=a)
#     thresholds <- dismo::threshold(e)
#     
#     final_distr <- pred_best_model
#     final_distr[final_distr < thresholds$equal_sens_spec] <- 0
#     final_distr[final_distr >= thresholds$equal_sens_spec] <- 1
#     
#     mod.null <- ENMnulls(maxent_model, 
#                          mod.settings = list(fc = c("LQ"), rm = 3), 
#                          removeMxTemp = TRUE,
#                          userStats.signs = NULL,
#                          no.iter = 50)
# 
#     list_outputs <- list(n_obs = nrow(occu),
#                          occu = occu,
#                          name = species_name,
#                          model_parameters = best_model$betas,
#                          AUC = e@auc,
#                          model_pred = pred_best_model,
#                          threshold_equal_sens_spec = thresholds$equal_sens_spec,
#                          model_pred_thresh = final_distr,
#                          emp_mean = null.emp.results(mod.null)[1,2],
#                          null_mean = null.emp.results(mod.null)[3,2],
#                          znull = null.emp.results(mod.null)[5,2],
#                          pnull = null.emp.results(mod.null)[6,2])
#     return(list_outputs)
#     
#     
#   }, error = function(e) {
#     final_distr <- raster::rasterize(occu, raster_cov[[1]], fun='count')
#     final_distr[final_distr > 0] <- 1
#     plot(final_distr)
#     list_outputs <- list(n_obs = nrow(occu),
#                          occu = occu,
#                          name = species_name,
#                          model_parameters = NA,
#                          AUC = NA,
#                          model_pred = NA,
#                          threshold_equal_sens_spec = NA,
#                          model_pred_thresh = final_distr,
#                          emp_mean = NA,
#                          null_mean = NA,
#                          znull = NA,
#                          pnull = NA)
#     return(list_outputs)
#     
#   })
# }
# 


# ### VARIABLE IMPORTANCE ESTIMATES USING DISMO PACKAGE
# eB_var_imp <- function(species_name, coord_data, raster_cov, ecorg, bias_window, denominator_pseudoabsence, dispersal_distance){
#   
#   
#   set.seed(123)
#   occu <- base::subset(coord_data, coord_data$nome_sp == species_name)
#   occu <- occu[,2:3]
#   occu <- raster::rasterize(occu, raster_cov[[1]], fun='count')
#   occu[occu > 0] <- 1
#   occu <- raster::rasterToPoints(occu)
#   if(nrow(occu) > 1){
#     occu <- occu[, 1:2]
#   } else {
#     occu <- occu[, 1:2]
#     occu <- t(as.data.frame(occu))
#   }
#   
#   tryCatch({
#     
#     names_raster_cov <- names(raster_cov)
#     extract_regions <- raster::extract(ecorg, occu)
#     extract_regions <- colSums(extract_regions, na.rm = TRUE)
#     extract_regions[extract_regions > 0] <- 1
#     extract_regions <- t(as.data.frame(extract_regions))
#     
#     ecorg_drop <- raster::dropLayer(ecorg, (colnames(extract_regions)[colSums(pmax(extract_regions == 0)) > 0]) )
#     ecorg_drop <- sum(ecorg_drop)
#     ecorg_drop[ecorg_drop > 0] <- 1
#     ecorg_drop[ecorg_drop == 0] <- NA
#     
#     raster_cov <- raster_cov * ecorg_drop[[1]]
#     crop_extent <- extent(min(occu[,1])- dispersal_distance, max(occu[,1]) + dispersal_distance, min(occu[,2]) - dispersal_distance, max(occu[,2]) + dispersal_distance)
#     raster_cov <- crop(raster_cov, crop_extent)
#     names(raster_cov) <- names_raster_cov
#     raster_cov <- scale(raster_cov)
#     
#     bias <- raster::rasterize(coord_data[, 2:3], raster_cov[[1]], fun='count') #
#     bias <- raster::focal(bias, w=matrix(1, bias_window, bias_window)  , fun=mean, pad=TRUE,  na.rm = TRUE )
#     bias <- raster::crop(bias, raster_cov[[1]])
#     bias <- raster::mask(bias, raster_cov[[1]])
#     bias <- log2(bias + 1)
#     pseudo_absences <- raster::xyFromCell(bias, sample(which(!is.na(values(bias))),
#                                                        (sum((!is.na(values(ecorg_drop))), na.rm=TRUE))/denominator_pseudoabsence, 
#                                                        prob=values(bias)[!is.na(values(bias))]))
#    
#     
#     object <- dismo::maxent(raster_cov, p = occu, a = pseudo_absences,
#                             args=c("linear=true","quadratic=true","product=false","hinge=false", "threshold=false"))
#     
#     importance <- maxent_get_results(object, 'importance') 
#     
#     list_outputs <- list(name = species_name,
#                          var_importance = importance)
#     return(list_outputs)
#     
#     
#   }, error = function(e) {
#     list_outputs <- list(name = species_name,
#                          var_importance = rep(NA,11))
#     return(list_outputs)
#     
#   })
# }

# MORAN'S I TO EVALUATE SPATIAL AUTOCORRELATION IN SPECIES RICHNESS MODELS
# library(spdep)
# library(spatialreg)
# 
# #replace m1 with m2, m3, m4, m5
# cor_r <- pgirmess::correlog(coords=m1[,c("x", "y")],
#                             z= (lm(m1$layer ~ m1$mean_t * m1$precip))$residuals,
#                             method="Moran", nbclass=10)
# 
# cor_r
# 
# correlograms <- as.data.frame(cor_r)
# correlograms$variable <- "residuals_glm" 
# 
# ggplot(subset(correlograms, variable=="residuals_glm"), aes(dist.class, coef)) + 
#   geom_hline(yintercept = 0, col="grey") +
#   geom_line(col="steelblue") + 
#   geom_point(col="steelblue") +
#   xlab("distance") + 
#   ylab("Moran's coefficient")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))


#################################  
### INFORMATION ON THE R AND PACKAGES VERSIONS
################################# 

# #sessionInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                    LC_TIME=English_Canada.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.5.0             ggspatial_1.1.7          sjPlot_2.8.12            brms_2.18.0              Rcpp_1.0.9               ggthemes_4.2.4          
# [7] ggplot2_3.4.0            tidyr_1.2.1              enmSdm_0.5.3.7           rangeModelMetadata_0.1.4 ENMeval_2.0.3            magrittr_2.0.3          
# [13] maptools_1.1-5           rgeos_0.5-9              RColorBrewer_1.1-3       SSDM_0.2.8               proj4_1.0-12             spatstat_3.0-2          
# [19] spatstat.linnet_3.0-3    spatstat.model_3.0-2     rpart_4.1.19             spatstat.explore_3.0-5   nlme_3.1-152             spatstat.random_3.0-1   
# [25] spatstat.geom_3.0-3      spatstat.data_3.0-0      rgdal_1.6-2              envirem_2.3              palinsol_0.93            gsl_2.1-7.1             
# [31] maps_3.4.1               viridis_0.6.2            viridisLite_0.4.1        dplyr_1.0.10             spatialEco_2.0-0         virtualspecies_1.5.1    
# [37] rJava_1.0-6              ENMTools_1.0.7           dismo_1.3-9              raster_3.6-11            sp_1.5-1                 ecospat_3.4             
# 
# loaded via a namespace (and not attached):
#   [1] estimability_1.4.1     spacetime_1.2-8        coda_0.19-4            intervals_0.15.2       earth_5.3.1            nabor_0.5.0            knitr_1.41            
# [8] dygraphs_1.1.1.6       multcomp_1.4-20        data.table_1.14.6      inline_0.3.19          generics_0.1.3         callr_3.7.3            terra_1.6-41          
# [15] TH.data_1.1-1          proxy_0.4-27           httpuv_1.6.6           StanHeaders_2.21.0-7   assertthat_0.2.1       xfun_0.35              bayesplot_1.10.0      
# [22] promises_1.2.0.1       fansi_1.0.3            igraph_1.3.5           DBI_1.1.3              htmlwidgets_1.5.4      reshape_0.8.9          tensorA_0.36.2        
# [29] stats4_4.1.0           purrr_0.3.5            ellipsis_0.3.2         crosstalk_1.2.0        ks_1.14.0              backports_1.4.1        insight_0.18.8        
# [36] permute_0.9-7          markdown_1.4           RcppParallel_5.1.5     deldir_1.0-6           vctrs_0.5.1            sjlabelled_1.2.0       abind_1.4-5           
# [43] withr_2.5.0            checkmate_2.1.0        emmeans_1.8.2          vegan_2.6-4            xts_0.12.2             prettyunits_1.1.1      mclust_6.0.0          
# [50] goftest_1.2-3          cluster_2.1.2          gstat_2.1-0            dotCall64_1.0-2        ape_5.6-2              crayon_1.5.2           units_0.8-0           
# [57] pkgconfig_2.0.3        nnet_7.3-16            rlang_1.0.6            spThin_0.2.0           lifecycle_1.0.3        miniUI_0.1.1.1         colourpicker_1.2.0    
# [64] sandwich_3.0-2         PresenceAbsence_1.1.10 biomod2_4.1-2          modelr_0.1.10          distributional_0.3.1   randomForest_4.7-1.1   polyclip_1.10-4       
# [71] matrixStats_0.63.0     datawizard_0.6.4       Matrix_1.5-3           loo_2.5.1              carData_3.0-5          boot_1.3-28            zoo_1.8-11            
# [78] base64enc_0.1-3        gamm4_0.2-6            processx_3.8.0         png_0.1-8              shinydashboard_0.7.2   KernSmooth_2.23-20     spam_2.9-1            
# [85] pROC_1.18.0            classInt_0.4-8         stringr_1.4.1          maxnet_0.1.4           jpeg_0.1-10            shinystan_2.6.0        rstatix_0.7.1         
# [92] ggeffects_1.1.4        ggsignif_0.6.4         scales_1.2.1           plyr_1.8.8             threejs_0.3.3          compiler_4.1.0         rstantools_2.2.0      
# [99] plotrix_3.8-2          lme4_1.1-31            cli_3.4.1              ade4_1.7-20            ps_1.7.2               Brobdingnag_1.2-9      htmlTable_2.4.1       
# [106] Formula_1.2-4          MASS_7.3-54            mgcv_1.8-35            tidyselect_1.2.0       stringi_1.7.8          forcats_0.5.2          projpred_2.2.2        
# [113] latticeExtra_0.6-30    bridgesampling_1.1-2   grid_4.1.0             shapefiles_0.7.2       tools_4.1.0            parallel_4.1.0         rstudioapi_0.14       
# [120] foreach_1.5.2          foreign_0.8-81         gridExtra_2.3          posterior_1.3.1        farver_2.1.1           digest_0.6.30          FNN_1.1.3.1           
# [127] shiny_1.7.3            pracma_2.4.2           TeachingDemos_2.12     car_3.1-1              broom_1.0.1            performance_0.10.1     later_1.3.0           
# [134] mda_0.5-3              sf_1.0-9               sjstats_0.18.2         colorspace_2.0-3       tensor_1.5             splines_4.1.0          fields_14.1           
# [141] spatstat.utils_3.0-1   shinythemes_1.2.0      RSAGA_1.3.0            xtable_1.8-4           nloptr_2.0.3           rstan_2.21.7           poibin_1.5            
# [148] plotmo_3.6.2           R6_2.5.1               Hmisc_4.7-2            lhs_1.1.5              pillar_1.8.1           htmltools_0.5.3        mime_0.12             
# [155] glue_1.6.2             fastmap_1.1.0          minqa_1.2.5            DT_0.26                class_7.3-19           codetools_0.2-18       pkgbuild_1.4.0        
# [162] mvtnorm_1.1-3          utf8_1.2.2             lattice_0.20-44        spatstat.sparse_3.0-0  tibble_3.1.8           gbm_2.1.8.1            gtools_3.9.4          
# [169] shinyjs_2.1.0          interp_1.1-3           survival_3.2-11        munsell_0.5.0          e1071_1.7-12           iterators_1.0.14       sjmisc_2.8.9          
# [176] reshape2_1.4.4         gtable_0.3.1           bayestestR_0.13.0   