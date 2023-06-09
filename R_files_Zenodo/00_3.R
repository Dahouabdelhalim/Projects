mod.out <- "C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/mod_areas/"
post.modFiles <- paste0(mod.out,list.files(mod.out))

nat.out <- "C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/nat_areas/"
post.natFiles <- paste0(nat.out,list.files(nat.out))

getLoss <- function(modern,natural){
  binom <- gsub(".*./|\\\\.tif","",modern)
  mod.range <- raster(modern)
  nat.range <- raster(natural)
  loss <- nat.range
  values(loss) <- 0
  loss[is.na(lc)] <- NA
  loss[mod.range==0 & nat.range==1] <- 1
  writeRaster(loss,paste0("C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/loss_areas/",binom,".tif"),overwrite = TRUE)
  return <- loss
}

for(i in 1 :length(post.modFiles)){
  thisLoss <- getLoss(post.modFiles[i],post.natFiles[i])
}

#allMod <- lc
#values(allMod) <- 0
#allMod[is.na(lc)] <- NA
#for(i in 1:length(post.modFiles)){
  #mod.range <- raster(post.modFiles[i])
  #allMod <- allMod + mod.range
#}
#writeRaster(allMod,"C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/masking/allMod.tif")


allNat <- lc
values(allNat) <- 0
allNat[is.na(lc)] <- NA
for(i in 1:length(post.natFiles)){
  nat.range <- raster(post.natFiles[i])
  allNat <- allNat + nat.range
}
writeRaster(allNat,"C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/masking/allNat.tif")
