library(plotrix)

XArea = 64
YArea = 64
xyArea = 0
d = sqrt(XArea/pi)+sqrt(YArea/pi)

AxyD = function(XArea,YArea,d){

if(XArea>YArea) stop("XArea must be smaller or equal to YArea")
   
rX = sqrt(XArea/pi)  
rY = sqrt(YArea/pi)  

#if ((rX^2+d^2)<=rY^2){
  
thetaX = ((rX^2)+(d^2)-(rY^2))/(2*d*rX)
thetaX=  2 * acos(thetaX)

thetaY = ((rY^2)+(d^2)-(rX^2))/(2*d*rY)
thetaY=  2 * acos(thetaY)

AxArch = rX^2 * thetaX/2
segX = 2 * sin(thetaX/2) * rX
hightX = cos(thetaX/2) * rX

AxTriangle = -(segX*hightX/2)
  
AyArch = rY^2 * thetaY/2

segY = 2 * sin(thetaY/2) * rY
hightY = cos(thetaY/2) * rY

AyTriangle = -(segY*hightY/2)

d_calc = hightX + hightY

#}




Axy = AxArch + AxTriangle + AyArch + AyTriangle
  

results = c(Axy, d, rX, rY,thetaX, thetaY,segX,segY,hightX,hightY,d_calc, AxArch,AxTriangle, AyArch, AyTriangle )
names(results) = c("Axy", "d", "rX", "rY", "thetaX","thetaY","segX","segY","hightX","hightY","d_calc", "AxArch","AxTriangle","AyArch","AyTriangle")

return(results)
 
}



XY_Circles = function(XArea,YArea, xyArea){

AreaMin = min(XArea,YArea)
AreaMax = max(XArea,YArea)
Axy = xyArea



if(Axy>AreaMin) stop("Areaxy, larger than one of the total Areas")

opt_d = function(d){
  
calAxy = AxyD(XArea = AreaMin,YArea = AreaMax, d=d)["Axy"]
return(abs(calAxy-Axy)) 
  
}

rX = sqrt(XArea/pi)  
rY = sqrt(YArea/pi)  

opt_d_est = optimise (opt_d, lower = abs(rX-rY), upper = rX+rY)

return(AxyD(XArea = AreaMin,YArea = AreaMax, d=opt_d_est$minimum))

}



results = XY_Circles(XArea = 90,YArea = 45, xyArea = 45)



ABCVenn = function (AreaA, AreaB, AreaC, AreaAB, AreaAC, AreaBC){


AB_dist = XY_Circles(XArea = AreaA,YArea = AreaB, xyArea = AreaAB)
AC_dist = XY_Circles(XArea = AreaA,YArea = AreaC, xyArea = AreaAC)
BC_dist = XY_Circles(XArea = AreaB,YArea = AreaC, xyArea = AreaBC)

A_coord = c(0,0)
B_coord = c(AB_dist["d"],0)

C_coodX = ((AB_dist["d"]^2)+(AC_dist["d"]^2)-(BC_dist["d"]^2))/(2*AB_dist["d"])
C_coody = sqrt((AC_dist["d"]^2)-(C_coodX^2))

C_coord = c(C_coodX, C_coody)

out = data.frame(Circle=c("A","B","C"), X = c(A_coord[1],B_coord[1],C_coord[1]), Y = c(A_coord[2],B_coord[2],C_coord[2]), radius = c(sqrt(AreaA/pi),sqrt(AreaB/pi),sqrt(AreaC/pi)))
return(out)
}

biotic = bioticRaster
abiotic = abioticRaster
movement = movementRaster

BAMplot = function(biotic, abiotic, movement){

warning("rasters should be thresholded models with 0 and 1 values an all have the same cell with NA")

BA = biotic * abiotic
AM = abiotic * movement
BM = biotic * movement  
BAM = abiotic * biotic * movement
 
areaRaster = area(movement, na.rm = T)
 
totalArea = sum(areaRaster[], na.rm = T)

B_area = sum((areaRaster*biotic)[], na.rm = T)
A_area = sum((areaRaster*abiotic)[], na.rm = T)
M_area = sum((areaRaster*movement)[], na.rm = T)

BA_area = sum((areaRaster*BA)[], na.rm = T)
AM_area = sum((areaRaster*AM)[], na.rm = T)
BM_area = sum((areaRaster*BM)[], na.rm = T)

BAM_Areas = sum((areaRaster*BAM)[], na.rm = T)

BAM_ComplexArea = B_area + A_area + M_area - BA_area - AM_area - BM_area + BAM_Areas

BAM_To_Continent = c(totalArea,BAM_ComplexArea,BAM_ComplexArea/totalArea)
names(BAM_To_Continent) = c("Continent_Area","BAM_Area", "BAM_Proportion")

out = ABCVenn(AreaA = A_area,
              AreaB = B_area,
              AreaC = M_area,
              AreaAB = BA_area,
              AreaAC = AM_area,
              AreaBC = BM_area)

maxX = max(out[1,"X"]+out[1,"radius"],out[2,"X"]+out[2,"radius"],out[3,"X"]+out[3,"radius"])
minX = min(out[1,"X"]-out[1,"radius"],out[2,"X"]-out[2,"radius"],out[3,"X"]-out[3,"radius"])

maxY = max(out[1,"Y"]+out[1,"radius"],out[2,"Y"]+out[2,"radius"],out[3,"Y"]+out[3,"radius"])
minY = min(out[1,"Y"]-out[1,"radius"],out[2,"Y"]-out[2,"radius"],out[3,"Y"]-out[3,"radius"])

coordsRange = range(maxX,minX,maxY,minY )

return(list(out,coordsRange, BAM_To_Continent ))

}




