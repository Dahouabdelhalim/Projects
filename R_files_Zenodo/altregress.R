#*** This analysis uses the data in Figure 1 and incorporates uncertainty into the correlation between Black Rail and VIRA turnover estimates ***#   

#Load packages
library(lme4)
require(deming) #Deming regression

#setwd()

rails <- read.csv("Both_rails_02_19_dynamics.csv")  #read in data from figure 1

str(rails)

#Ordinary correlations between turnover rates of black and Virginia rails with plots
plot(Lambda_BLRA~Lambda_VIRA,data=rails) #plot of Lambda estimates
with(rails,cor.test(Lambda_BLRA,Lambda_VIRA)) #simple correlation test

plot(Extinct_BLRA~Extinct_VIRA,data=rails) #plot of extinction rates
with(rails,cor.test(Extinct_BLRA,Extinct_VIRA)) #simple correlation test

plot(Colon_BLRA~Colon_VIRA,data=rails) #plot of colonization rates
with(rails,cor.test(Colon_BLRA,Colon_VIRA)) #simple correlation test

plot(Occupied_BLRA~Occupied_VIRA,data=rails) #plot of occupancy estimates
with(rails,cor.test(Occupied_BLRA,Occupied_VIRA)) #simple correlation test


#Deming regressions
d1=deming(Extinct_BLRA~Extinct_VIRA,data=rails,xstd=Se_Extin_VIRA, ystd=Se_Extin_BLRA);d1
d2=deming(Colon_BLRA~Colon_VIRA,data=rails,xstd=Se_Colon_VIRA, ystd=Se_Colon_BLRA);d2
d3=deming(Lambda_BLRA~Lambda_VIRA,data=rails,xstd=Se_Lambda_VIRA, ystd=Se_Lambda_BLRA);d3
d4=deming(Occupied_BLRA~Occupied_VIRA,data=rails,xstd=Occup_Se_VIRA, ystd=Occup_Se_BLRA);d4



