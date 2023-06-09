#This script splits the original database into the n main traits included in the BeechCOSTe52

#set working directory 
setwd("/WorkingDirectory/")

setwd("/home/marta/Dropbox/Database-check/")
Fsyl <- read.csv("finals for repository/Fsylvatica.csv")


#Read databases
Fsyl <- read.csv("Fsylvatica.csv")


#List available traits: 
#Tree height: 
H95cm,H96cm,H97cm,H98cm,H99cm,H00cm,H01cm,H02cm,H03cm,H04cm,H05cm,H06cm,H07cm,H08cm
#DBH
DBH02mm,DBH03mm,DBH04mm,DBH05mm,DBH06mm,DBH07mm,DBH08mm
#Basal Diameter
BasalD96mm,BasalD01mm,BasalD02mm,BasalD04mm,BasalD05mm,BasalD06mm,BasalD07mm
#Mortality
Mort95,Mort97,Mort99,Mort00,Mort01,Mort02,Mort07
#Spring phenology
#Autumn phenology

#Tree height: 
#H95cm,H96cm,H97cm,H98cm,H99cm,H00cm,H01cm,H02cm,H03cm,H04cm,H05cm,H06cm,H07cm,H08cm

height <- subset(Fsyl,(H95cm != "NA" | H96cm != "NA" | H97cm != "NA" |H98cm != "NA" |H99cm != "NA" |H00cm != "NA" |H01cm != "NA" |H02cm != "NA" |H03cm != "NA" |H04cm != "NA" |H05cm != "NA" |H06cm != "NA" |H07cm != "NA" |H08cm != "NA" ))
height[,(11:24)]
TH <- height[,c(1:10,11:24)]

dim(TH)# 148182 individual records of tree height
write.csv(TH, file = paste("traits/Fsylvatica_height.csv",sep="" ), sep="\\t", col.names=TRUE)

#DBH
#DBH02mm,DBH03mm,DBH04mm,DBH05mm,DBH06mm,DBH07mm,DBH08mm

DBH <- subset(Fsyl,(DBH02mm != "NA" | DBH03mm != "NA" | DBH04mm != "NA" |DBH05mm != "NA" |DBH06mm != "NA" |DBH07mm != "NA" |DBH08mm != "NA"  ))
DBH[,(25:31)]
diam <- DBH[,c(1:10,25:31)]

dim(diam)# 54051 individual records of DBH
write.csv(diam, file = paste("traits/Fsylvatica_DBH.csv",sep="" ), sep="\\t", col.names=TRUE)

#Basal Diameter
BasalD96mm,BasalD01mm,BasalD02mm,BasalD04mm,BasalD05mm,BasalD06mm,BasalD07mm

BasalD <- subset(Fsyl,(BasalD96mm != "NA" | BasalD01mm != "NA" | BasalD02mm != "NA" |BasalD04mm != "NA" |BasalD05mm != "NA" |BasalD06mm != "NA" |BasalD07mm != "NA"  ))
BasalD[,(32:38)]
BD <- BasalD[,c(1:10,32:38)]

dim(BD)# 20271 individual records of DBH
write.csv(BD, file = paste("traits/Fsylvatica_DBH.csv",sep="" ), sep="\\t", col.names=TRUE)

#Mortality
#Mort95,Mort97,Mort99,Mort00,Mort01,Mort02,Mort07

mortality <- subset(Fsyl,(Mort95 != "NA" | Mort97 != "NA" | Mort99 != "NA" |Mort00 != "NA" |Mort01 != "NA" |Mort02 != "NA" |Mort07 != "NA" ))


mortality[,(39:45)]
mort <- mortality[,c(1:10,39:45)]

dim(mort)# 45328 individual records of tree mortality
write.csv(mort, file = paste("traits/Fsylvatica_mortality.csv",sep="" ), sep="\\t", col.names=TRUE)

#Spring Phenology
flush <- subset(Fsyl, (SprPhen00C1 != "NA" | SprPhen00C1Date != "NA" |  SprPhen01C1 != "NA" |  SprPhen01C1Date != "NA" |  SprPhen02C1 != "NA" | SprPhen02C1Date != "NA" | SprPhen02C2 != "NA" | SprPhen02C2Date != "NA" | SprPhen04C1 != "NA" | SprPhen04C1Date != "NA" | SprPhen04C2 != "NA" | SprPhen04C2Date != "NA" |      
SprPhen04C3 != "NA" | SprPhen04C3Date != "NA" | SprPhen04C4 != "NA" | SprPhen04C4Date != "NA" | SprPhen05C06 != "NA" | SprPhen05C06Date != "NA" |      SprPhen05C07 != "NA" | SprPhen05C07Date != "NA" | SprPhen05C1 != "NA" | SprPhen05C1Date != "NA" | SprPhen05C2 != "NA" | SprPhen05C2Date != "NA" |  SprPhen05C3 != "NA" | SprPhen05C3Date != "NA" | SprPhen05C4 != "NA" | SprPhen05C4Date != "NA" | SprPhen05C5 != "NA" | SprPhen05C5Date != "NA" |  SprPhen06C1 != "NA" | SprPhen06C1Date != "NA" | SprPhen06C2 != "NA" | SprPhen06C2Date != "NA" | SprPhen06C3 != "NA" | SprPhen06C3Date != "NA" |      
SprPhen06C4 != "NA" | SprPhen06C4Date != "NA" | SprPhen06C5 != "NA" | SprPhen06C5Date != "NA" | SprPhen07C1 != "NA" | SprPhen07C10 != "NA" | SprPhen07C10Date != "NA" |      SprPhen07C11 != "NA" | SprPhen07C11Date != "NA" | SprPhen07C12 != "NA" | SprPhen07C12Date != "NA" | SprPhen07C13 != "NA" | SprPhen07C13Date != "NA" |   SprPhen07C1Date != "NA" | SprPhen07C2 != "NA" |  SprPhen07C2Date != "NA" | SprPhen07C3 != "NA" | SprPhen07C3Date != "NA" |  SprPhen07C4 != "NA" | SprPhen07C4Date != "NA" | SprPhen07C5 != "NA" |  SprPhen07C5Date != "NA" | SprPhen07C6 != "NA" |  SprPhen07C6Date != "NA" |      SprPhen07C7 != "NA" | SprPhen07C7Date != "NA" | SprPhen07C8 != "NA" |          SprPhen07C8Date != "NA" | SprPhen07C9 != "NA" | SprPhen07C9Date != "NA" |      SprPhen08C1 != "NA" | SprPhen08C1Date != "NA" | SprPhen08C2 != "NA" | SprPhen08C2Date != "NA" | SprPhen08C3 != "NA" |  SprPhen08C3Date != "NA" | SprPhen08C4 != "NA" | SprPhen08C4Date != "NA" | SprPhen08C5 != "NA" |       
SprPhen08C5Date != "NA" | SprPhen08C6 != "NA" | SprPhen08C6Date != "NA" | SprPhen08C7 != "NA" | SprPhen08C7Date != "NA" | SprPhen08C8 != "NA" |          
SprPhen08C8Date != "NA" | SprPhen09C1 != "NA" | SprPhen09C1Date != "NA" | SprPhen97C1 != "NA" | SprPhen97C1Date != "NA" | SprPhen97C2 != "NA" |      
SprPhen97C2Date != "NA" | SprPhen99C1 != "NA" | SprPhen99C1Date != "NA" | SprPhen99C2 != "NA" | SprPhen99C2Date != "NA" | SprPhen99C3 != "NA" |     SprPhen99C3Date != "NA" | SprPhen99C4 != "NA" | SprPhen99C4Date!= "NA" ))

flush[,(46:141)]
BB <- flush[,c(1:10,46:141)]

dim(BB)# 70849 individual records of spring phenology
write.csv(BB, file = paste("traits/Fsylvatica_SpringPheno.csv",sep="" ), sep="\\t", col.names=TRUE)


#Autumn phenology

senesc <- subset(Fsyl, (Autumn01C1 != "NA" | Autumn01C1Date != "NA" | Autumn04C1 != "NA" | Autumn04C1Date != "NA" | Autumn05C1 != "NA" | Autumn05C1Date != "NA" | Autumn05C2 != "NA" | Autumn05C2Date != "NA" | Autumn05C3 != "NA" | Autumn05C3Date != "NA" | Autumn05C4 != "NA" | Autumn05C4Date != "NA" | Autumn05C5 != "NA" | Autumn05C5Date != "NA" | Autumn05C6 != "NA" | Autumn05C6Date != "NA" | Autumn06C1 != "NA" | Autumn06C1Date != "NA" | Autumn06C2 != "NA" | Autumn06C2Date != "NA" |   Autumn06C3 != "NA" | Autumn06C3Date != "NA" | Autumn06C4 != "NA" | Autumn06C4Date != "NA" | Autumn06C5 != "NA" | Autumn06C5Date != "NA" | Autumn06C6 != "NA" | Autumn06C6Date != "NA" | Autumn07C1 != "NA" |  Autumn07C1Date != "NA" | Autumn07C2 != "NA" | Autumn07C2Date != "NA" | Autumn07C3 != "NA" | Autumn07C3Date != "NA" | Autumn07C4 != "NA" |          Autumn07C4Date != "NA" | Autumn07C5 != "NA" | Autumn07C5Date != "NA" | Autumn07C6 != "NA" | Autumn07C6Date != "NA" | Autumn08C1 != "NA" | Autumn08C10 != "NA" | Autumn08C10Date != "NA" | Autumn08C1Date != "NA" | Autumn08C2 != "NA" | Autumn08C2Date != "NA" | Autumn08C3 != "NA" | Autumn08C3Date != "NA" | Autumn08C4 != "NA" | Autumn08C4Date != "NA" | Autumn08C5 != "NA" | Autumn08C5Date != "NA" | Autumn08C6 != "NA" | Autumn08C6Date != "NA" | Autumn08C7 != "NA" | Autumn08C7Date != "NA" | Autumn08C8 != "NA" | Autumn08C8Date != "NA" | Autumn08C9 != "NA" | Autumn08C9Date != "NA"))


senesc[,(161:220)]
AP <- senesc[,c(1:10,161:220)]

dim(AP)# 1758 individual records of autumn phenology
write.csv(AP, file = paste("traits/Fsylvatica_AutumnPheno.csv",sep="" ), sep="\\t", col.names=TRUE)



