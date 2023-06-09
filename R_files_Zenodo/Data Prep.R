source("Otter_tg_Functions.R")
data <- read.csv("otter_tg_dataset.csv",header=T)

data$AC2 <- data$AGECLASS

data$AC2[which(data$AC2=="P")] <- "J"
data$AC2[which(data$AC2=="IM")] <- "J"
data$AC2[which(data$AC2=="A-AA")] <- "AA"
data$AC2 <- as.factor(data$AC2)
data$AC2[which(data$AC2=="")] <- NA
data$AC2 <- relevel(factor(data$AC2),ref="J")

data$AC3 <- data$AC2
data$AC3[which(data$AC3=="AA")] <- "A"
data$AC3 <- relevel(factor(data$AC3),ref="J")

data$AC4 <- data$AC3
data$AC4[which(data$AC4=="AA")] <- "A"
data$AC4 <- factor(data$AC4)
data$AC4 <- relevel(data$AC4,ref="J")

data$AC5 <- data$AC4
data$AC5[which(data$AC5=="J")]<-"SA"
data$AC5 <- relevel(factor(data$AC5),ref="SA")

data$PD10 <- data$PD/10
data$HD10 <- data$HD/10
data$RD10 <- data$RD/10

data$HDbin <- 1
data$HDbin[data$HD<5] <- 0

data$PD5cat <- ntiles(data$PD,5)
data$HD5cat <- ntiles(data$HD,5)
data$RD5cat <- ntiles(data$RD,5)

data$PD4cat <- ntiles(data$PD,4)
data$HD4cat <- ntiles(data$HD,4)
data$RD4cat <- ntiles(data$RD,4)

data$Devall <- with(data,rowSums(cbind(NLCDV21,NLCDV22,NLCDV23,NLCDV24)))
data$Devlowplus <- with(data,rowSums(cbind(NLCDV22,NLCDV23,NLCDV24)))
data$Devmedplus <- with(data,rowSums(cbind(NLCDV23,NLCDV24)))
data$Devhigh <- data$NLCDV24

data$imp10 <- data$imp/10
data$imp4 <- ntiles(data$imp)

data$impHI <- 0
data$impHI[data$imp>=1] <- 1

d <- data[,18:57]
d[is.na(d)] <- 0
data[,18:57] <- d

data$Dev10 <- 10*(data$Devall+data$LCIV17)
data$Frm10 <- 10*(data$NLCDV82+data$LCIV15)
data$Grz1 <- 10*data$NLCDV81
data$For10 <- 10*(data$NLCDV41+data$NLCDV42+data$LCIV1+data$LCIV5+data$LCIV6)
data$Wet1 <- 100*(data$NLCDV90+data$NLCDV95+data$LCIV14)
data$Scrub10 <- 10*(data$NLCDV51+data$NLCDV52+data$LCIV8)
data$Grass10 <- 10*(data$NLCDV71+data$LCIV10)
data$Othnat10 <- data$Grass10+data$Scrub10+10*(data$LCIV12+data$LCIV13)

data$length <- as.numeric(data$LENGTH)
data$length[data$length==0] <- NA
data$lengthper10cm <- data$length/10

data$AREA2 <- as.character(data$SITE)
calif <- subset(data,!is.na(data$ATOS))
calif$area <- ws_area(calif$SITE,calif$ATOS,calif$YEAR)
data$AREA2[!is.na(data$ATOS)] <- calif$area
data$AREA2 <- factor(data$AREA2)
data$toxo <- factor(data$TGPOS)

datastore <- data
