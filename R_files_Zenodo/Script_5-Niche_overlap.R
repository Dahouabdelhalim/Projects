# load image 
load('Script_4_image.RData')

# NMI and NICHE OVERLAP
# load functions from Broenniman et al
source("https://raw.githubusercontent.com/ecospat/NMI/master/functions/NMI_function.R")
source("https://raw.githubusercontent.com/ecospat/NMI/master/functions/mve_function.R")
### choose the parameters of the analysis
grain=10 # resolution of the climatic data. In Broennimann et al. 2020 NCOM: 10,30,60
envelope="kde" # choise of the envelop methode, "kde" or "mve"
level=99 # level of inclusion of rare climatic conditions for kde and mve. In Broennimann et al. 2020 NCOM: 99,95,90

current_climate <- raster::crop(current_climate, aus)
current_climate <- raster::mask(current_climate, aus)
clim.df<-na.exclude(getValues(current_climate))

pca<- dudi.pca(clim.df,scannf = FALSE, nf = 2)

# get current outline of Bd
bd_current <- raster('Result_rasters/Bd_current_projection_binary.tif')
bd_full <- raster('Result_rasters/Bd_current_projection.tif')
bd_current <- rasterToPolygons(bd_current, dissolve = TRUE)

sp.scores<-na.omit(suprow(pca,raster::extract(current_climate,bd_current))$li) 

info <- read.csv('Species_Information.csv')

bd_species <- gsub(' ', '_', subset(info, Bd_decline_Scheele_2019 == 'Yes')$Species)

results <- list()

for (i in 1:length(bd_species)){
  tryCatch({ 
    # select native niche of focal species and extract scores in PCA space
    print(i)
    print(Sys.time())
    # create scores of introductions in PCA space
    sp_current <- raster(paste('Result_rasters/',
                               bd_species[i],
                               '_current_projection_binary.tif',
                               sep = ''))
    sp <- xyFromCell(sp_current, cell = 1:ncell(sp_current))
    sp <- cbind(sp,raster::extract(sp_current, sp))
    sp <- sp[complete.cases(sp),]
    
    intros.scores<-suprow(pca,raster::extract(current_climate,sp[,1:2]))$li
    NArows<-is.na(intros.scores$Axis1)
    intros<-sp[!NArows,]
    intros.scores<-intros.scores[!NArows,]
    
    ### delineate niche margins
    # kde
    if(envelope=="kde"){
      fhat<-kde(sp.scores,compute.cont=TRUE)
      c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=fhat$cont[level])
      l99<-list()        
      for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
      sp.pol=SpatialPolygons(l99)
      sp.pol<-aggregate(sp.pol)
    }
    # mve
    if(envelope=="mve"){sp.pol <- mve(sp.scores, thresh=level/100, robust=F)}
    
    ### calculate NMI for introductions
    coordinates(intros.scores) <- cbind(intros.scores$Axis1 , intros.scores$Axis2) # DataFrame to SpatialPointsDataFrame
    intro.NMI<-NMI(foc.pop = intros.scores,niche=sp.pol)
    nmi<-intro.NMI$NMI
    nmi_val <- mean(nmi)
    
    # load full raster for niche overlap
    sp_full <- raster(paste('Result_rasters/',
                            bd_species[i],
                            '_current_projection.tif',
                            sep = ''))
    overlap_mets <- raster.overlap(bd_full, sp_full)
    results[[i]] <- data.frame(Species = bd_species[i],
                               Climate_Model = NA,
                               Emissions = NA,
                               Year = 'Current',
                               NMI = nmi_val,
                               D = overlap_mets$D,
                               I = overlap_mets$I,
                               Spearman = overlap_mets$rank.cor)
  }, error=function(e){cat("ERROR :", conditionMessage(e),"\\n")})
}
library(tidyverse)
current_results <- bind_rows(results)

# make a vector of climate models we will investigate for future time periods
climate_models <- list.files(path = 'Climate/Future',
                             pattern = 'ssp126_2021-2040.tif')
climate_models <- gsub('wc2.1_2.5m_bioc_',
                       '',
                       climate_models)
climate_models <- gsub('_ssp126_2021-2040.tif',
                       '',
                       climate_models)
# make a vector of time scenarios
time <- c('2021-2040',
          '2041-2060',
          '2061-2080',
          '2081-2100')
# and emissions scenarios
emissions <- c('ssp126', 'ssp370')
# read in elevation for later
elevation <- raster('Climate/wc2.1_2.5m_elev.tif')

# nor for future time periods
num = 1
period_results <- list()
for(i in 1:length(climate_models)) {
  for(j in 1:length(time)) {
    for(k in 1:length(emissions)) {
      current_climate <- stack(raster(paste('Climate/Future/wc2.1_2.5m_bioc_',climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''), band = 1),
                               raster(paste('Climate/Future/wc2.1_2.5m_bioc_',climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''), band = 2),
                               raster(paste('Climate/Future/wc2.1_2.5m_bioc_',climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''), band = 12),
                               raster(paste('Climate/Future/wc2.1_2.5m_bioc_',climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''), band = 16),
                               raster(paste('Climate/Future/wc2.1_2.5m_bioc_',climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''), band = 17),
                               raster('Climate/wc2.1_2.5m_elev.tif'))     
      current_climate <- raster::crop(current_climate, aus)
      current_climate <- raster::mask(current_climate, aus)
      clim.df<-na.exclude(getValues(current_climate))
      pca<- dudi.pca(clim.df,scannf = FALSE, nf = 2)
      # get current outline of Bd
      bd <- raster(paste('Result_rasters/Bd_', climate_models[i],'_',emissions[k],'_',time[j],'._binary.tif', sep = ''))
      bd_full <- raster(paste('Result_rasters/Bd_', climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''))
      bd <- rasterToPolygons(bd, dissolve = TRUE)
      sp.scores<-na.omit(suprow(pca,raster::extract(current_climate,bd))$li) 
      for (m in 1:nrow(current_results)){
        tryCatch({ 
          # select native niche of focal species and extract scores in PCA space
          print(m)
          print(Sys.time())
          # create scores of introductions in PCA space
          sp_current <- raster(paste('Result_rasters/',
                                     current_results$Species[m],
                                     '_', climate_models[i],'_',emissions[k],'_',time[j],'_binary.tif', sep = ''))
          sp <- xyFromCell(sp_current, cell = 1:ncell(sp_current))
          sp <- cbind(sp,raster::extract(sp_current, sp))
          sp <- sp[complete.cases(sp),]
          
          intros.scores<-suprow(pca,raster::extract(current_climate,sp[,1:2]))$li
          NArows<-is.na(intros.scores$Axis1)
          intros<-sp[!NArows,]
          intros.scores<-intros.scores[!NArows,]
          
          ### delineate niche margins
          # kde
          if(envelope=="kde"){
            fhat<-kde(sp.scores,compute.cont=TRUE)
            c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=fhat$cont[level])
            l99<-list()        
            for(p in 1:length(c99))l99[[p]]=Polygons(list(Polygon(c99[[p]][-1])),ID=p)
            sp.pol=SpatialPolygons(l99)
            sp.pol<-aggregate(sp.pol)
          }
          # mve
          if(envelope=="mve"){sp.pol <- mve(sp.scores, thresh=level/100, robust=F)}
          
          ### calculate NMI for introductions
          coordinates(intros.scores) <- cbind(intros.scores$Axis1 , intros.scores$Axis2) # DataFrame to SpatialPointsDataFrame
          intro.NMI<-NMI(foc.pop = intros.scores,niche=sp.pol)
          nmi<-intro.NMI$NMI
          nmi_val <- mean(nmi)
          
          sp_full <- raster(paste('Result_rasters/', current_results$Species[m],'_', climate_models[i],'_',emissions[k],'_',time[j],'.tif', sep = ''))
          overlap_mets <- raster.overlap(bd_full, sp_full)
          period_results[[num]] <- data.frame(Species = current_results$Species[m],
                                              Climate_Model = climate_models[i],
                                              Emissions = emissions[k],
                                              Year = time[j],
                                              NMI = nmi_val,
                                              D = overlap_mets$D,
                                              I = overlap_mets$I,
                                              Spearman = overlap_mets$rank.cor)
          num <- num+1
        }, error=function(e){cat("ERROR :", conditionMessage(e),"\\n")})
      }
    }
  }
} 

write.csv(current_results,
          'CURRENTRESULTS.csv')
write.csv(bind_rows(period_results),
          'PERIODRESULTS.csv')

save.image('Script_5_image.RData')

