# ______________________________________________________________________________________ 
# ______________________________________________________________________________________ 
# ______________________________________________________________________________________ 
# |                                                                                    |
# |          SCRIPTS REALIZED BY Fabrice VINATIER fabrice.vinatier@inrae.fr            | 
# |                          and Denis FEURER denis.feurer@ird.fr                      | 
# |                              ----------------                                      | 
# |                              LICENCE CC-BY-SA                                      | 
# |                              ----------------                                      |
# | This license lets others remix, adapt, and build upon your work even for           |
# | commercial purposes, as long as they credit you and license their new creations    |
# | under the identical terms.                                                         |
# |                                                                                    |
# | The proposed code has a purely academic purpose, is valid under the conditions     |
# | of use of the scientific project for which it was funded and at the date of        |
# | acceptance of the article presenting the code. As with any research work, the      |
# | code is not free of possible errors, approximations, sub-optimisations or          |
# | defects in monitoring dependencies between libraries of the programme.             |
# |                                                                                    |
# ______________________________________________________________________________________ 
# |                                                                                    |
# | Cette licence permet à d'autres personnes de remixer, d'adapter et de              |
# | développer ce travail, même à des fins commerciales, à condition qu'elles          |
# | créditent l'auteur et accordent une licence pour leurs nouvelles créations aux     |
# | mêmes conditions.                                                                  |
# |                                                                                    |
# | Le code proposé a une visée purement académique, est valable dans les conditions   |
# | d'utilisation du projet scientifique pour lequel il a été financé et à la date de  |
# | d'acceptation de l'article de présentation du code.                                |
# | Comme tout travail de recherche, le code n'est pas exempt d'éventuelles erreurs,   |
# | approximations, sous-optimisations ou défauts de suivi des dépendances entre       |
# | sous-éléments du programme.                                                        |
# ______________________________________________________________________________________ 
# ______________________________________________________________________________________ 
# ______________________________________________________________________________________ 
# |                                                                                    |
# |                    PREPARATION OF THE WORK ENVIRONMENT                             |
# ______________________________________________________________________________________ ---------------------------------------------------------------------

rm(list=ls(all=TRUE))
# Libraries                                                                            |----
library(raster)
library(plyr)
library(rgeos)
library(scales)
library(prettymapr)
# ______________________________________________________________________________________ ---------------------------------------------------------------------
# |                                                                                    | ------------------------------------------------------------------------
# |                   MAIN FUNCTIONS                                                   | ------------------------------------------------------
# ______________________________________________________________________________________ ---------------------------------------------------------------------

# | Hole filling on a raster                     fill.na                             | ----------------------------------------------
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( mean(x, na.rm=TRUE) )
  } else {
    return( x[i] )
  }
}

# | Find zone of prospection    (CLOUDCOMPARE)   ProspectZone      Fabrice Vinatier  | ----------------------------------------------
ProspectZone=function(pathRef,pathCompared,pathCLOUDS){
  if(Sys.info()["sysname"]=="Windows")pathToProgram="C:/PROGRA~1/CloudCompare/CloudCompare"
  if(Sys.info()["sysname"]=="Linux")  pathToProgram="cloudcompare.CloudCompare"
  
  # Rasterization of each cloud
  system(paste(pathToProgram,"-SILENT",
               "-O",
               "-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),pathRef,
               "-SOR 6 1",
               #"AUTO_SAVE ON",
               "-RASTERIZE -GRID_STEP 0.001 -EMPTY_FILL INTERP -OUTPUT_RASTER_Z",collapse=" "),ignore.stderr = TRUE,ignore.stdout = TRUE)
  rREF=readAll(raster(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_BEFORE_SOR_RASTER_Z_"),sep=""),band=1))
  system(paste(pathToProgram,"-SILENT",
               "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),pathCompared,
               "-SOR 6 1",
               "-RASTERIZE -GRID_STEP 0.001 -EMPTY_FILL INTERP -OUTPUT_RASTER_Z",collapse=" "),ignore.stderr = TRUE,ignore.stdout = TRUE)
  rCOMP=readAll(raster(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_AFTER_SOR_RASTER_Z_"),sep=""),band=1))
  
  # Get new filenames
  listFilesRef=dir(pathCLOUDS,pattern=gsub(".ply","_SOR_",tail(strsplit(pathRef,"/")[[1]],1)))
  listFilesRef=listFilesRef[grep(".bin",listFilesRef)]
  pathRefNew=paste(pathCLOUDS,listFilesRef,sep="")
  listFilesCompared=dir(pathCLOUDS,pattern=gsub(".ply","_SOR_",tail(strsplit(pathCompared,"/")[[1]],1)))
  listFilesCompared=listFilesCompared[grep(".bin",listFilesCompared)]
  pathComparedNew=paste(pathCLOUDS,listFilesCompared,sep="")
  
  # Cloud to cloud difference using cloud compare, then rasterization of the cloud difference
  system(paste(pathToProgram,"-SILENT",
               "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),pathRef, 
               "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),pathCompared,
               "-c2c_dist -MAX_DIST 0.1",
               "-AUTO_SAVE OFF",
               "-RASTERIZE -GRID_STEP 0.001 -EMPTY_FILL INTERP -OUTPUT_RASTER_Z",collapse=" "),ignore.stderr = TRUE,ignore.stdout = TRUE)
  # Selection of the type of survey using the character string pathRef
  typeReleve=strsplit(tail(strsplit(pathRef,"/")[[1]],1),"_")[[1]][3]
  method=strsplit(pathRef,"/")[[1]][11]
  zoneStudy=strsplit(pathRef,"/")[[1]][1]
  
  # Estimation of the hole contour (in 2D)
  tmp=raster(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_C2C_DIST_MAX_DIST_0.1_RASTER_Z_"),sep=""),band=2) # downloading the raster of cloud difference
  tmp=aggregate(tmp,fact=10) # aggregation of the raster 
  tmp=focal(tmp,w=matrix(1,3,3),fun=fill.na,na.rm=F,pad=T) # Filling the little holes in the raster
  tmp=focal(tmp,w=matrix(1/25,nrow=5,ncol=5),fun=mean,na.rm=T) # Smoothing the values
  
  tmp
}
  
# | Delimitate border of excavated holes         delimitateHoles   Fabrice Vinatier  | ----------------------------------------------
delimitateHoles=function(rasterZone=rast_sel,tol_simplify=0.01,width_buffer=0.02,area_hole=0.008,thres_hole=0.00021){
  mask_zone=(rasterZone>thres_hole)      # Thresholding the smoothed raster to identify the hole
  rCLUMP=clump(mask_zone)         # detecting patches of connected cells with a unique ID
  dCLUMP=na.omit(data.frame(freq(rCLUMP))) # Estimation of areas of patches (in number of cells)
  maxV=dCLUMP[dCLUMP$count*mean(res(rCLUMP))^2>area_hole,"value"] # Identification of hole contour that reached a minimal area (multiple holes)
  rCLUMP[is.na(rCLUMP)]=0      # Selection of hole contour (1/3)
  rCLUMP[!rCLUMP[]%in%maxV]=NA # Selection of hole contour (2/3)
  rCLUMP[rCLUMP[]%in%maxV]=1   # Selection of hole contour (3/3)
  vCLUMP=rasterToPolygons(rCLUMP,dissolve=T) # Polygonization of the hole contour
  vCLUMP=gSimplify(vCLUMP,tol=tol_simplify)  # Reduction of stairs in the polygon
  vCLUMP=disaggregate(vCLUMP)                # Splitting multiple polygons
  vCLUMP=gBuffer(vCLUMP,width=width_buffer,byid=T) # Enlargement of polygon area
  vCLUMP=vCLUMP[order(coordinates(vCLUMP)[,1]),]   # Reordering the multiple polygons from left to right
  file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_C2C_"),sep="")) # Removing temporary CloudCompare files
  vCLUMP
}  

# | Volume calculation    (CLOUDCOMPARE)         calcVol           Fabrice Vinatier  | ----------------------------------------------
calcVol=function(pathRef,pathCompared,pathCLOUDS,holes=holes_sel){
  if(Sys.info()["sysname"]=="Windows")pathToProgram="C:/PROGRA~1/CloudCompare/CloudCompare"
  if(Sys.info()["sysname"]=="Linux")  pathToProgram="cloudcompare.CloudCompare"
  
  # Get new filenames
  listFilesRef=dir(pathCLOUDS,pattern=gsub(".ply","_SOR_",tail(strsplit(pathRef,"/")[[1]],1)))
  listFilesRef=listFilesRef[grep(".bin",listFilesRef)]
  pathRefNew=paste(pathCLOUDS,listFilesRef,sep="")
  listFilesCompared=dir(pathCLOUDS,pattern=gsub(".ply","_SOR_",tail(strsplit(pathCompared,"/")[[1]],1)))
  listFilesCompared=listFilesCompared[grep(".bin",listFilesCompared)]
  pathComparedNew=paste(pathCLOUDS,listFilesCompared,sep="")
  
  # Estimation of volume behind hole contours
  rCROPref=list() ; rCROPcomp=list()
  vols=sapply(1:length(holes),function(winSel){ # calculation on every polygons
    dWIN=holes[winSel,]@polygons[[1]]@Polygons[[1]]@coords # Selection of the hole
    dWIN=as.data.frame(t(t(dWIN))) # Recuperation of the numeric coordinates of the hole vertexes
    # Cropping the cloud around hole borders, meshing, then resampling at high density of points (Reference)
    system(paste(pathToProgram,"-SILENT",
                 "-C_EXPORT_FMT ASC",
                 "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),pathRefNew,
                 "-CROP2D Z",
                 dim(dWIN)[1],
                 paste(c(t(dWIN),recursive=T),collapse=" "),
                 "-DELAUNAY -BEST_FIT",
                 "-SAMPLE_MESH DENSITY 10000000",
                 "-RASTERIZE -GRID_STEP 0.001 -EMPTY_FILL INTERP -OUTPUT_RASTER_Z",
                 collapse=" "),ignore.stderr = TRUE,ignore.stdout = TRUE)
    
    # Cropping the cloud around hole borders, meshing, then resampling at high density of points (Comparison)
    system(paste(pathToProgram,"-SILENT",
                 "-C_EXPORT_FMT ASC",
                 "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),pathComparedNew,
                 "-CROP2D Z",
                 dim(dWIN)[1],
                 paste(c(t(dWIN),recursive=T),collapse=" "),
                 "-DELAUNAY -BEST_FIT",
                 "-SAMPLE_MESH DENSITY 10000000",
                 "-RASTERIZE -GRID_STEP 0.001 -EMPTY_FILL INTERP -OUTPUT_RASTER_Z",
                 collapse=" "),ignore.stderr = TRUE,ignore.stdout = TRUE)
    # Get the rasterized files
    rasterCROPED=dir(pathCLOUDS,pattern="_CROPPED_SAMPLED_POINTS_RASTER_Z_")
    rCROPref[[winSel]] <<-readAll(raster(paste(pathCLOUDS,rasterCROPED[grep("_BEFORE_",rasterCROPED)],sep=""),band=1)) # downloading the raster of cloud difference
    rCROPcomp[[winSel]]<<-readAll(raster(paste(pathCLOUDS,rasterCROPED[grep("_AFTER_",rasterCROPED)],sep=""),band=1))# downloading the raster of cloud difference
    file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_RASTER_"),sep=""))
    # Recuperation of the transformed clouds
    cloudCROPED=dir(pathCLOUDS,pattern="_CROPPED_SAMPLED_POINTS_")
    
    cropRef=paste(pathCLOUDS,cloudCROPED[grep("_BEFORE_",cloudCROPED)],sep="")
    cropCompared=paste(pathCLOUDS,cloudCROPED[grep("_AFTER_",cloudCROPED)],sep="")
    # Volume calculation between the comparison and the reference clouds
    system(paste("cloudcompare.CloudCompare","-SILENT",
                 "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),cropCompared, 
                 "-O","-GLOBAL_SHIFT",paste(c(0,0,0),collapse=" "),cropRef,
                 "-VOLUME -GRID_STEP 0.001",
                 collapse=" "),ignore.stderr = TRUE,ignore.stdout = TRUE)
    # Recuperation of the result
    vol=readLines(paste(pathCLOUDS,dir(pathCLOUDS,pattern="VolumeCalculationReport_"),sep=""),4)[4]
    vol=strsplit(vol," ")[[1]][3]
    vol=as.numeric(substr(vol,4,20))
    # Removing temporary CloudCompare files
    file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_CROPPED_"),sep=""))
    file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="VolumeCalculationReport_"),sep=""))
    vol})
  # Attribution of a ID for every volume calculated from multiple or single holes
  file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_SOR_"),sep=""))
  repD=paste("D",1:length(vols),sep="")
  # Results
  data.frame(Replicate=repD,Volume=vols)
}

# | Filtering multiple ply                       filterRRF         Fabrice Vinatier  | ----------------------------------------------
filterRRF=function(pathIN=pathRef,pathOUT="IN/2018-10_Mauguio/Flash_Lidar/Fosse1_H3_V_BEFORE.ply"){
  # Selection of all ply files
  filesPLY=sapply(dir(pathIN,pattern=".ply$"),function(file_sel){vcgPlyRead(paste(pathIN,"/",file_sel,sep=""))$vb},simplify="array")
  # Application of a median filter on depth
  medPLY=apply(filesPLY,1:2,median)
  tmp=t(medPLY)
  tmp=as.data.frame(tmp)
  tmp=tmp[tmp$V1!=0 & tmp$V2!=0 & tmp$V3!=0,]
  medPLY=t(tmp)
  # Formatting the resulting data.frame in a ply file for export
  medPLY=list(vb=medPLY,material=list())
  attr(medPLY,"class")="mesh3d"
  vcgPlyWrite(medPLY,filename = pathOUT)
}


# ______________________________________________________________________________________ ---------------------------------------------------------------------
# |                                                                                    | ------------------------------------------------------------------------
# |                     VOLUME CALCULATION                                             | ------------------------------------------------------
# ______________________________________________________________________________________ ---------------------------------------------------------------------

# Selection of the clouds (Ref: before excavation and Compared: after excavation)
pathRef=file.choose()
pathCompared=file.choose()

# Suppression of the residual cloudCompare files
endFile=strsplit(pathRef,"/")
pathCLOUDS=sub(endFile[[1]][length(endFile[[1]])],"",pathRef)
file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_C2C_"),sep=""),showWarnings=F)
file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_SOR"),sep=""),showWarnings=F)
file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="VolumeCalculationReport_"),sep=""),showWarnings=F)
file.remove(paste(pathCLOUDS,dir(pathCLOUDS,pattern="_RASTER_"),sep=""),showWarnings=F)

# Volume calculation
zone_tot=ProspectZone(pathRef,pathCompared,pathCLOUDS); plot(zone_tot) # Calculation of the prospecting zone
zone_sel=drawExtent()  # selection of the prospecting zone
rast_sel=crop(zone_tot,zone_sel)
holes_sel=delimitateHoles(rasterZone=rast_sel,tol_simplify=0.01,width_buffer=0.02,area_hole=0.008,thres_hole=0.00021)
plot(rast_sel)
plot(holes_sel,add=T)
calcVol(pathRef,pathCompared,pathCLOUDS,holes=holes_sel)
