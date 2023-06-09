#define a function to refine a range by land cover-based habitat suitability

refine.habitats <- function(range,type){
  #develop a binomial from the file name of the range
  binom <- gsub(".*./|\\\\.tif","",range)
  u_binom <- binom
  binom <- gsub("_"," ",binom)
  #use the IUCN Redlist API to collect habitat suitabilities for the species based on its name
  suit <- rl_habitats(suitability = c("Suitable","Marginal"),name = binom, key = "bcdf6849b4e6df1e0cecaf68490770dff7406cca9ae953a75174947117dd8d79")
  
  #take out the habitat codes
  suit <- suit$result
  habs <- suit$code
  
  #simplify them to first-level to match the crosswalk used in this analysis
  habs <- unique(substr(habs,1,1))
  
  #load the range as a raster
  pre.range <- raster(range)
  pre.range <- projectRaster(pre.range,lc,method = "ngb")
  if(min(values(pre.range),na.rm=T)==1){
    pre.range[is.na(pre.range)==T] <- 0
  }
  pre.range[is.na(pre.range)&!is.na(lc)] <- 0
  
  #relate suitable habitats to valid land covers using the crosswalk
  land.covers <- cw[cw$IUCN_habitat %in% habs,]
  land.covers <- unique(land.covers$ESA_CCI)
  
  #refine the range by suitable land covers as described above
  fullLC <- lc
  fullLC[!(fullLC %in% land.covers)] <- 0
  fullLC[(fullLC %in% land.covers)] <- 1
  fullLC[!(pre.range==1)] <- 0
  if(type=="mod"){
  writeRaster(fullLC, paste0("C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/mod_areas/",u_binom,".tif"),overwrite = TRUE)
  }
  if(type=="nat"){
  writeRaster(fullLC, paste0("C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/nat_areas/",u_binom,".tif"),overwrite = TRUE)
  }
}

lapply(modFiles[1:327],FUN = function(x1) refine.habitats(x1,"mod"))
lapply(natFiles[1:327],FUN = function(x1) refine.habitats(x1,"nat"))
