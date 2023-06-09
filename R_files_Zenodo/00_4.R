loss.out <- "C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/loss_areas/"
lossFiles <- paste0(loss.out,list.files(mod.out,pattern=".tif"))


allLoss <- lc
values(allLoss) <- 0
allLoss[is.na(lc)] <- NA

for(i in 1:length(lossFiles)){
  loss.range <- raster(lossFiles[i])
  allLoss <- allLoss + loss.range
}
writeRaster(allLoss,"C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/intact_areas/allLoss.tif")

wilmaOut <- allLoss
#loads the raster containing all modern species ranges
allMod <- raster("C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/masking/allMod.tif")

#loads the raster containing all natural species ranges
allNat <- raster("C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/masking/allNat.tif")


#negates any areas where no WILMA species currently live
wilmaOut[allMod == 0] <- NA

#negates any area where all "natural" species have been lost
wilmaOut[wilmaOut >= allNat] <- NA

plot(wilmaOut)


pureWilma <- wilmaOut

#filters areas by the paper's 1-3 species focus
pureWilma[wilmaOut > 3] <- NA


plot(pureWilma, col = c("#30BF4F","#8AEB6A","#D9FF5C","#F0E67F"),main = "Missing large mammal species")
writeRaster(pureWilma,"C:/Users/JoeG/OneDrive - WCMC/WILMA2/outputs/outputs_pres_only/intact_areas/draft_wilma.tif")
  