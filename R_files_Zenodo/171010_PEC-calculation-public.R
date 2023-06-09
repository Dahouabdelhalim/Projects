##  ##   ##  ##  ## ##  ## ##   ##   ##	   Concentrations  ##   ##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##
##  ##   ##  ##  ## ##  ## ##   ##   ##   ##   START  ##   ##   ##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

CF		<- matrix(0,6,2) ###Matrix for defined compartment volumes###
CF[,1]	<- c(1:6)

#Preparation of the compartment volumes


######## AIR
### EDIT HW: no AIR concentration calculated due to the assumed full dry/wet deposition of SiO2, but description has been left for explanatory uses.

##*10^12, transformation tons in µg 
##/365*10), residence time of 10 days was assumed corresponding to ultrafine particles (5).
##*1, the relevant air volume was calculated by
##multiplying the area of the particular region by the recommended value for atmospheric
##from http://europa.eu/about-eu/facts-figures/living/index_en.htm (4463000km2)
##*10^9, transformation km3 in m3 

#CF[1,2]	<- (10^12/365*10)/(4463000*1*10^9)


######## SEDIMENTS
##*10^9 transformation tons in mg
##4326337*(1-0.97), surface  of sediements (equivalent to surface water) in km2, 
##*10^6, transformation km2 in m2
##*0.03, sediment depth, the depth of the sediment of 0.03 m (1)
##0.2*1300, the density of dry sediment was calculated by subtracting the
##water content (80%) from the standardized value of 1300 kg m-3 recommended by the European
##Commission (1) resulting in 260 kg m-3.
# worst-case scenario 100% sedimentation

CF[1,2]	<- 10^9/(4326337*(1-0.97)*10^6*0.03*260) #See AIR for the explanation


######## SOIL (NATURAL AND URBAN)
##*10^12, transformation tons in µg 
##4326337 surface area of EU 27 in km2 ////  4,463,000 km2, according to http://europa.eu/about-eu/facts-figures/living/index_en.htm
##*0.97, proportion relevant soil not covered with concrete, 3% covered with concrete
##*10^6, transformation km2 in m2
##*0.2*0.27, soil depth of agricultural soil  times the proportion of such soil, see also sludge treated soil
##0.05*0.73, soil depth urban soil times the proportion of such soil, The soil
##volume was calculated by multiplying the soil depth depending on the mixing depth of
##different soil types (natural and urban soil: 0.05 m, agricultural soil: 0.2 m (1))
##(4326337*0.97*10^6*(0.2*0.27*0.99+0.05*0.73) Soil volume=346.6km3 (not 34 as indicated in support information) 
##*0.6*2500=1500 kg/m3 density of dry soil (for such volume), The density of dry
##soil was calculated by subtracting the water content from the standardized value of 1700 kg
##m-3 recommended by the European Commission (1) resulting in 1500 kg m-3; we used as equivalent calcualtion 0.6
##volume fraction solids in soil of the density of the solid phase resulting in 1500 kg m-3

CF[2,2]	<- 10^12/(4326337*0.97*10^6*(0.2*0.47+0.05*0.53)*1500) # According to the site below, 43% is agricultural soil, not 27%

#Denominator is equal to 7.585259e+14



######## SOIL (SLUDGE TREATED AGRICULTURAL SOIL)
## In 2009, 43% of the total EU suface was agricultural land. http://ec.europa.eu/eurostat/statistics-explained/index.php/Land_cover,_land_use_and_landscape
## According to Sun et al. 2014, 1% of agricultural soil is treated with sludge.
##0.2 m is the agricultural mixing depth

CF[3,2] <- (10^12/((9000000*0.55/20)*10^4*0.2*1500)) #Denominator equal to 7.425e+11


######## Surface water
##*10^12,transformation tons in µg 
##*(40/365), residence time in the system EU, according to the technical guidance document (1) a residence time of 40 days was applied
##4326337*(1-0.97), surface  of surface water in km2 (3%), the technical guidance document (1) 
##*10^6,transformation km2 in m2
##*3, water depth, water mixing depth of 3 m, the technical guidance document (1) 
##*1000, transformation of m3 in dm3=l

CF[4,2]	<-	10^12*(40/365)/(4326337*(1-0.97)*10^6*3*1000)


######## WIP waste
##*10^9 transformation tons in mg
##190'560'000 t Annual Municipal Solid Waste amount in EU27
## 0.2-0.3 municipal solid waste incineration rate in Europe
## *1000 transformation tons in kg

#CF[6,2]	<-	10^9/(190560000*0.3*1000)


######## WWTP effluent
##*10^12,transformation tons in µg 
##0.8*200*365 a daily water consumption per inhabitant of 200 L for Switzerland and the EU was assumed;the percentage of sources
##connected to sewage treatment facilities is 80% for the EU (17) 
##509000000 (http://data.worldbank.org/region/EUU)  // Again, I dont like that number
##According to http://ec.europa.eu/eurostat/tgm/table.do?tab=table&init=1&language=en&pcode=tps00001&plugin=1
##
CF[5,2]	<-	10^12/(0.8*200*365*509000000)


######## WWTP sludge
##*10^9 transformation kg in mg
##The considered annual dry sewage treatment sludge productions for the studied regions were:
## 9,000,000 t (4) for the EU.
## *1000transformation tons in kg

CF[6,2]	<- 10^9/(9000000*1000)


##  ##   ##  ##  ## ##  ## ##   ##   ##	   Concentrations  ##   ##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##
##  ##   ##  ##  ## ##  ## ##   ##   ##   ##   END  ##   ##   ##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##




####################
###Concentrations in technical and environmental compartments
####################

PECs	<- matrix(NA,6,SIM)

PECs[1,]	<- CF[1,2]*Mass["Sediments",]	#Sediments
PECs[2,]	<- CF[2,2]*Mass["Soil",]	#Soil (Natural & Urban)
PECs[3,]	<- CF[3,2]*Mass["Soil_sludge_treat",]+PECs[2,]	#Soil (Sludge Treated) + Soil concentration
PECs[4,]	<- CF[4,2]*Mass["Surfacewater",]	#Surface water
PECs[5,]	<- CF[5,2]*Mass["WWTP_effl",]	#WWTP effluent
PECs[6,]	<- CF[6,2]*Mass["Sewage_sludge_treat",]	#WWTP sludge

#########################################################################

TITLE	<- c("nano-SiO2 in sediments","nano-SiO2 in N&U soil","nano-SiO2 in ST soil","nano-SiO2 in surface water","nano-SiO2 in WWTP effluent","nano-SiO2 in WWTP sludge")
XLABEL	<- c("Concentration (µg/kg)","Concentration (µg/kg)","Concentration (µg/kg)","Concentration (µg/l)","Concentration (µg/l)","Concentration (mg/kg)")

for(i in 1:6)
{
  plot(density(PECs[i,], cut=10), main=TITLE[i],font.main=1,lwd=2,col="sandybrown", xlab=XLABEL[i],ylab="Density")
}


##  ##   ##  ##  ## ##   Export central measurements (q15,mode...) data to CSV file      	 ##  ##  ##  ##  ##  ## 
##  ##   ##  ##  ## ##  ## ##   ##   ##   ##   START  ##   ##   ##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  

#Mode calculation
Mode_Y	<- function(x)
{
  dens	<- density(x)
  ind		<- which(dens$y==max(dens$y))
  dens$x[ind]
}

####################
###Flows (preparing the data)
####################



flows <- matrix(NA,74,5)

for(i in 1:74)
{
  xx			<- Mass[i,]
  
  flows[i,1]	<- i
  flows[i,2]	<- quantile(xx,0.15)
  flows[i,3]	<- Mode_Y(xx)[1]
  flows[i,4]	<- mean(xx)
  flows[i,5]	<- quantile(xx,0.85)	
}

####################
###Concentrations (preparing the data)
####################

conc <- matrix(NA,6,5)

for(i in 1:6)
{
  xx	<- PECs[i,]
  
  conc[i,1]	<- i
  conc[i,2]	<- quantile(xx,0.15)
  conc[i,3]	<- Mode_Y(xx)[1]
  conc[i,4]	<- mean(xx)
  conc[i,5]	<- quantile(xx,0.85)	
}

####################
###Exporting. Provide the desired name and folder
####################

measures			<- flows
colnames(measures)	<- c("flow","Q15","Mode","Mean","Q85")
rownames(measures) <- names(inp)
measures_table		<- as.table(measures)
write.table(measures_table, "INSERTFILENAME")


measures			<- conc
colnames(measures)	<- c("Row","Q15","Mode","Mean","Q85")
measures_table		<- as.table(measures)
write.table(measures_table, "INSERTFILENAME")


pdf(file = "Output_of_Trial_PEC.pdf",
    height = 7.5,
    width = 7.5,
    pointsize = 10)

par(mfrow = c(3,3), mar = c(3,3,3,1), mgp = c(1.5,0.5,0), xpd = F)
color <- adjustcolor(rainbow(dim(PECs)[1]), alpha.f = 0.6)
for(co in 1:dim(PECs)[1]){
  density(PECs[co,], freq = F, xlab = paste(dimnames(PECs)[[1]][co], "(tonnes)"),
       ylab = "Probability density", main = dimnames(PECs)[[1]][co],
       col = color[co])
  legend("topright", c(paste("Mean:", round(mean(PECs[co,]), digits = 4)),
                       paste("SD:", round(sd(PECs[co,]), digits = 4))),
         bty = "n", cex = 0.85)
  box()
}



##  ##   ##  ##  ## ##   Export central measurements (q15,mode...) data to CSV file      	 ##  ##  ##  ##  ##  ## 
##  ##   ##  ##  ## ##  ## ##   ##   ##   ##   END  ##   ##   ##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 