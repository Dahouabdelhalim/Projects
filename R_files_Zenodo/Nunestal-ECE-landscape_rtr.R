################ rtr function used in:
#Nunes et al - Identifying mismatches between conservation area networks and vulnerable populations using spatial randomization

##landscape === polygon to do the RTR (e.g Focal Landscape)
##bios == raster of the background area to sample (e.g. State of Wisconsin with hostile habitat masked)
## tolerance == proportion of the landscape that is inside the raster (e.g.,not falling on open water or hostile habitat)

rtr<-function(landscape,bios,tolerance,rotation=TRUE,translation=TRUE){
  library('raster')
  library('maptools')
  library('dismo')
  library('rgeos') #for convexhull and centroid

  poly<-SpatialPoints(landscape) 
  crs(poly)<-crs(landscape)
  poly<-gConvexHull(poly) #make minimum convex polygon
  cent<-gCentroid(poly) #find centroid of polygon
  
  x_centre<-as.vector(extent(cent))[1]  #find x- centre of polygon
  y_centre<-as.vector(extent(cent))[3]  #find y- centre of polygon
  
  repeat{
    ###### rotation
    if(rotation==TRUE){ ##option to rotate the points, no translation
      deg<-runif(1,0,360) #random angle in degrees
      
      P.new<-maptools::elide(landscape, rotate=deg,center=c(x_centre,y_centre))
     
    }
    if(translation==TRUE){ #option to translate the points, no rotation
      #### translation 
      xmin=extent(bios)[1];xmax=extent(bios)[2];ymin=extent(bios)[3];ymax=extent(bios)[4]
      
      bb<-bbox(P.new)
      bb.x<-bb[1,]
      bb.y<-bb[2,]
      
      
      x_trans<-runif(1,-(max(xmin,bb.x[1])-min(xmin,bb.x[1])),max(xmax,bb.x[2])-min(xmax,bb.x[2])) #random long trans
      y_trans<-runif(1,-(max(ymin,bb.y[1])-min(ymin,bb.y[1])),max(ymax,bb.y[2])-min(ymax,bb.y[2])) #random long trans
      
      P.new<-maptools::elide(P.new, shift=c(x_trans,y_trans),center=c(x_centre,y_centre))
    }

    landscape.details = extract(bios, P.new) #if NA, then a point is in the ocean/outside study region - no variables
    if(isTRUE((length(which(is.na(landscape.details[[1]])==FALSE))/length(landscape.details[[1]]))>=tolerance)){  #if TRUE no NAs then break loop, keep RTR polygon
      break}}
  

gridcells<-extract(bios, P.new,cellnumbers=TRUE)[[1]][,1] #store gridcells to check for duplicates

sp.abundance<-sum(na.omit(landscape.details[[1]]))

prop.abundance<-sp.abundance/sum(na.omit(values(bios[[1]])))

    return(list(sp.abundance,prop.abundance,P.new,as.vector(gridcells)))} #return niche overlap of RTR and the gridcells occupied by the simulated points
