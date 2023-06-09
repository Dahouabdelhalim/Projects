### Analysis Tutorial_1
### Project Artenvielfalt erleben ###

# set working directory
setwd("/Users/juliamoritz/Desktop/Citizen_Science_PJ/Tutorialstudie/Zenodo/")
rm(list=ls())

library(dplyr)
library(lme4)
library(car)
library(emmeans)
library(afex)
source("https://sebastiansauer.github.io/Rcode/logit2prob.R")

# read data
dat <- read.csv("ORIGINALDATA_JAMZ20200527_Artenvielfalt_Tutorialstudie.csv", header=T)
dat <- dat[-c(1:2),]
labels <- read.csv("ORIGINALDATA_JAMZ20200527_Artenvielfalt_Tutorialstudie.csv", header=T)
labels <- labels[c(1:2),]
# add subject number
dat$VPNr <- 1:length(dat$ResponseId)

#### delete all columns with irrelevant information
dat$IPAddress <- NULL
dat$RecipientFirstName <- NULL
dat$RecipientLastName <- NULL
dat$RecipientEmail <- NULL
dat$ExternalReference <- NULL
dat$LocationLatitude <- NULL
dat$LocationLongitude <- NULL
dat$IPAddress <- NULL
dat$UserLanguage <- NULL 

# total duration of experiment -> StartDate; EndDate; Duration; RecordedDate 
dat$RecordedDate <- NULL # allways 2 sec after "EndDate"
dat$StartDate <- NULL 
dat$EndDate <- NULL 

# Fortschritt -> Progress in %; Finished in 1 vs. 0
#dat$Progress <- NULL
#dat$Finished <- NULL

# distribution channel
dat$DistributionChannel <- NULL # wie wurde es verteilt -> alle haben den selben Link; evtl herausfiltern von Testdurchlaeufen - falls noch welche gemacht wurden
dat$Status <- NULL # Testdurchlaeufe koennen auch hier herausgefilter werden

# consent and data retraction
Rueck <- subset(dat, dat$R端ckzug2 == 1 | dat$Einwilligung != 1, select = c("ResponseId", "VPNr")) # List with VP who wanted to retract their data or did not confirm the consent
dat <- subset(dat, dat$Einwilligung == 1) # Einwilligung mit 1
dat <- subset(dat, dat$R端ckzug2 != 1) # Datenrueckzug mit Antwort 1 der zweiten Frage 
dat$Einwilligung <- NULL
dat$R端ckzug1 <- NULL
dat$R端ckzug2 <- NULL


########## Die beiden Bedingungen (mit & ohne Tutorial) - kodieren und Variablen untereinander setzten ##########
# 1. Zwei Datensaetze (dat3, dat5) jeweils fuer eine Kondition erstellen, die jeweils nur die Variablen enthalten, die sich unterscheiden
# 2. jeweils in den beiden Datensaetzen die leeren Felder herausloeschen
# 3. Variable fuer die Bedingung erstellen
# 4. die beiden Datensaetze (dat3, dat5) zusammenfuegen (dat35)
# 5. Datensatz (dat0) erstellen, der nur Variablen beinhaltet, die fuer alle Versuchspersonen gleich sind
# 6. aus dat0 und dat35 einen finalen Datensatz (newData) zusammenfuegen

dat$cond <- "control"
dat$cond[is.na(as.numeric(as.character(dat$Video_3Legende_t_Page.Submit)))==F] <- "tutorial"
dat$cond <- as.factor(dat$cond)
dat_c <- dat[dat$cond=="control",-c(19:36)]
dat_t <- dat[dat$cond=="tutorial",-c(5:18)]

t_names <- names(dat_t)
t_names <- gsub("Video3", "X3", t_names)
t_names <- gsub("Video_3", "X3", t_names)
t_names <- gsub("Verstandnis", "Verstaendnis", t_names)

names(dat_t) <- t_names
dat_c$Video_t_First.Click <- NA        
dat_c$Video_t_Last.Click <- NA
dat_c$Video_t_Page.Submit <- NA
dat_c$Video_t_Click.Count <- NA

setdiff(names(dat_t), names(dat_c))
setdiff(names(dat_c), names(dat_t))

data <- rbind(dat_t, dat_c)
names(data) <- gsub("X3", "", names(data))

########## Klickfelder reduzieren auf jeweils 6 Variablen #############

names(data)[names(data) == 'Reg1_12'] <- '1,2,8' # 1,2,7 und 2,2,7 sind doppelt (2,8 fehlt jeweils)Reg1/2_12 eigentlich 2,8
names(data)[names(data) == 'Reg2_12'] <- '2,2,8'

RegtList1 <- names(select(data, starts_with("Reg1_t")))
RegtList2 <- names(select(data, starts_with("Reg2_t")))
RegList1 <- names(select(select(data, starts_with("Reg1")), -RegtList1)) # Liste, die alle Namen der Klickfelder beinhalten
RegList2 <- names(select(select(data, starts_with("Reg2")), -RegtList2))

# Variablen nach ihrer Matrix-Koordinate benennen
for (i in 1:length(names(data))) { # Iteration durch aller Variablen in dat35
  oldname <- names(data)[i] # aktueller Variablenname
  if(oldname %in% RegList1){ # nur Variablen waehlen, die ein Klickfeld beschreiben (hier fuer Region1)
    Erklaerung <- levels(data[,i])[3] # Auswahl der Beschreibung/Erklaerung der Variable
    MatrixPos <- substr(Erklaerung, nchar(Erklaerung)-2, nchar(Erklaerung)) # die letzten 3 Characters beinhalten die Matrixposition (Achtung die Reihe 10 wird zu 0)
    names(data)[names(data) == oldname] <- paste('1', MatrixPos, sep = ',') # Variable neu benennen nach ihrer Matrix-Koordinate, mit einer 1 fuer die Region-Karte vorangestellt
  }
  else if(oldname %in% RegList2){ # gleiches Vorgehen, nur fuer Region2
    Erklaerung <- levels(data[,i])[3]
    MatrixPos <- substr(Erklaerung, nchar(Erklaerung)-2, nchar(Erklaerung))
    names(data)[names(data) == oldname] <- paste('2', MatrixPos, sep = ',')
  }
}


data$FieldList <- list(character()) # neue Variable mit einer Liste fuer jede VP anlegen
MatrixFields <- names(select(data, c(starts_with("1,"), starts_with("2,")))) # Liste aller Variablennamen,die ein Klickfeld benennen

# Erstellen der Listen fuer jede VP
for (i in 1:length(names(data))){ # Iteration durch jede Variable in dat35
  name <- names(data)[i] # Name der aktuellen Variable
  if(name %in% MatrixFields){ # Nur Variablen waehlen, bei denen es sich um Klickfelder handelt
    for (j in 1:length(data$VPNr)) { # Iteration durch die VP
      if(data[j,i]=="On"){ # Auswahl der Klickfelder, die ausgewaehlt wurden
        data$FieldList[[j]] <- c(data$FieldList[[j]], name) # Anhaengen der ausgewaehlten Felder an FieldList der VP
      }
    }
  }
}

# Ausschliessen der Clickfelder aus dem Datensatz
data <- select(data, -MatrixFields)

# neue Variablen mit den x- & y-Matrix-Koordinaten von jedem gewaehlten Klickfeld, und fuer Region-Karte 1 & 2
for (i in 1:length(data$VPNr)){ # Iteration fuer jede VP
  Reg1x1 <- substring(data$FieldList[[i]][1], 3, 3) # FieldList enthaelt strings ("1,2,3" [region,x-Koordinate, y-Koordinate]); 
  Reg1y1 <- substring(data$FieldList[[i]][1], 5, 5) # der 5. string bezeichnet die y-Koordinate
  Reg1x2 <- substring(data$FieldList[[i]][2], 3, 3) # der 3. string bezeichnet die x-Koordinate
  Reg1y2 <- substring(data$FieldList[[i]][2], 5, 5) # insgesamt enthaelt jede Liste aus FieldList 6 strings mit jeweils einer Matrix-Koordinate
  Reg1x3 <- substring(data$FieldList[[i]][3], 3, 3)
  Reg1y3 <- substring(data$FieldList[[i]][3], 5, 5)
  data$Reg1[i] <- list(c(Reg1x1, Reg1y1, Reg1x2, Reg1y2, Reg1x3, Reg1y3))
  
  Reg2x1 <- substring(data$FieldList[[i]][4], 3, 3)
  Reg2y1 <- substring(data$FieldList[[i]][4], 5, 5)
  Reg2x2 <- substring(data$FieldList[[i]][5], 3, 3)
  Reg2y2 <- substring(data$FieldList[[i]][5], 5, 5)
  Reg2x3 <- substring(data$FieldList[[i]][6], 3, 3)
  Reg2y3 <- substring(data$FieldList[[i]][6], 5, 5)
  data$Reg2[i] <- list(c(Reg2x1, Reg2y1, Reg2x2, Reg2y2, Reg2x3, Reg2y3))
}

#Ausschliessen der Variable mit den Listen der 6 Matrixkoordinaten
data$FieldList <- NULL

newData <- data
########### Richtige Antworten codieren: ##########
##### fuer die Klickkarten
# 1. liegen die Daten nebeneinander?
# Funktion gibt TRUE zurueck, wenn zwei Felder nebeneinander liegen
neben <- function(x1, y1, x2, y2){
  if((x1 == x2+1 || x1 == x2-1) && y1 == y2){ # die Felder liegen nebeneinander, wenn x2 eins groesser oder kleiner als x1 ist
    return(TRUE);                             # und y1 gleich gross wie y2 ist -> in der selben y-Reihe
  }
  else if((y1 == y2+1 || y1 == y2-1) && x1 == x2){
    return(TRUE);
  }
  else{
    return(FALSE);
  }
}

# Funktion gibt TRUE zurueck, wenn die drei Felder nebeneinander liegen
verbunden <- function(l){ # Eingabe einer Liste l(x1,y1 , x2,y2 , x3,y3)
  if(is.na(l[1])){
    return(NA)
  }
  else{
    if(neben (l[1], l[2], l[3], l[4])){  # wenn Feld1(x1,y1) neben Feld2(x2,y2) liegt 
    return(neben(l[1], l[2], l[5], l[6]) || neben(l[3], l[4], l[5], l[6]))   # dann muss Feld3(x3,y3) neben Feld1 oder Feld2 liegen
    }
    else if(neben(l[1], l[2], l[5], l[6])){  # sonst wenn Feld1 neben Feld3 liegt 
      return(neben(l[3], l[4], l[5], l[6]))  # muss Feld2 neben Feld3 liegen (Feld1 kann in diesem Fall nicht neben Feld2 liegen!)
    }
    else{ # sonst liegen die Felder nicht nebeneinander
      return(FALSE)
    }
  }
}

# neue Variable, die enthaelt ob die Felder korrekter Weise nebeneinander liegen
data$Reg1 <- lapply(data$Reg1, as.numeric)
data$Reg1 <- lapply(data$Reg1, function(x) ifelse(x==0, 10, x)) # die 0 bedeutet Reihe 10
data$verb1 <- lapply(data$Reg1, verbunden)

data$Reg2 <- lapply(data$Reg2, as.numeric)
data$Reg2 <- lapply(data$Reg2, function(x) ifelse(x==0, 10, x))
data$verb2 <- lapply(data$Reg2, verbunden)

# 2. Wurden sie inhaltlich korrekt ausgewaehlt?
# Matrix mit Sichtungen: 0 - kein roter Punkt; 1 - ein roter Punkt
sichtungsM <- rbind(c(1,1,1,1,1,0,1,0),
                    c(1,0,0,0,0,0,0,0),
                    c(1,0,1,0,0,0,0,1),
                    c(1,0,0,1,1,1,1,1),
                    c(0,0,0,0,1,0,0,1),
                    c(0,0,0,0,0,0,0,0),
                    c(0,0,0,0,0,0,0,0),
                    c(1,1,1,0,0,0,0,1),
                    c(1,1,1,1,1,0,0,0),
                    c(0,0,1,1,1,1,1,0))

# Matrix mit Beobachtungsintensitaet: 0 - keine Beob.; 1 - ; 2 - ; 3 - ; 4 - ;
# alle >2 sind in der dreistufigen Bedingung 2 
intensitaetsM <- rbind(c(4,3,2,1,0,1,1,1),
                       c(3,1,1,0,1,1,1,1),
                       c(2,1,1,1,1,1,1,3),
                       c(3,1,1,3,3,4,3,1),
                       c(1,1,1,2,4,3,2,1),
                       c(1,1,1,1,1,1,3,1),
                       c(1,1,1,1,2,2,4,4),
                       c(4,4,2,1,3,1,2,2),
                       c(4,3,4,4,1,2,2,1),
                       c(3,2,3,3,3,2,2,2))

# Funktion gibt an, ob alle drei Felder keine Sichtung enthalten und die Intensitaet ueber 1 liegt
# Richtige Antwort fuer Reg1: Alle drei Felder duerfen keine Sichtung enthalten, die Intensitaet muss groesser 1 sein
CorrIntens <- function(l){
  return((intensitaetsM[l[1],l[2]] > 1 && intensitaetsM[l[3],l[4]] > 1 && intensitaetsM[l[5],l[6]] > 1) &&
           (sichtungsM[l[1],l[2]] == 0 && sichtungsM[l[3],l[4]] == 0 && sichtungsM[l[5],l[6]] == 0))
}

# neue Variable, die angibt, ob fuer Region 1 inhaltlich richtige Felder gewaehlt wurden
data$CorrAns1 <- lapply(data$Reg1, CorrIntens)

# neue Variable Reg1_CorrAns, die angibt, ob sowohl die Felder nebeneinander liegen, als auch die richtigen Felder gewaehlt wurden
for (i in 1:length(data$VPNr)){
  data$Reg1_CorrAns[i] <- data$verb1[[i]] && data$CorrAns1[[i]]
}

# Funktion gibt an, ob alle drei Felder eine Sichtung enthalten,
# Richtige Antwort fuer Reg2: Alle drei Felder muessen eine Sichtung enthalten
CorrSichtung <- function(l){
  return(sichtungsM[l[1],l[2]] == 1 && sichtungsM[l[3],l[4]] == 1 && sichtungsM[l[5],l[6]] == 1) 
}
# neue Variable, die angibt, ob fuer Region 2 inhaltlich richtige Felder gewaehlt wurden
data$CorrAns2 <- lapply(data$Reg2, CorrSichtung)

# neue Variable Reg1_CorrAns, die angibt, ob sowohl die Felder nebeneinander liegen, als auch die richtigen Felder gewaehlt wurden
for(i in 1:length(data$VPNr)){
  data$Reg2_CorrAns[i] <- data$verb2[[i]] && data$CorrAns2[[i]]
}


##### fuer die Karten

# Liste mit allen richtigen Antworten:
corrAnsList <- c(1,3,3,3,3,3,2,3,3,1,3,1,3,3,1)

# Variable fuer jede richtige Antwort mit Mapx.y_CorrA
l <- 1
for(i in 1:3){
  for(j in 1:5){
      corrA <- corrAnsList[l]
      var <- paste0("Map", i, ".", j, "_corrA")
      data$new <- corrA
      colnames(data)[colnames(data) == 'new'] <- var
      l <- l+1
  }
}

# Umbennenen der Variablennamen ohne ".Click", ".Submit", ".Count"
clickList <- names(select(data, ends_with(".Click")))
submitList <- names(select(data, ends_with(".Submit")))
countList <- names(select(data, ends_with(".Count")))
nameList <- names(data)

for (i in 1:length(nameList)) { # Iteration durch alle Variablen in data
  oldname <- nameList[i] # aktueller Variablenname
  if(oldname %in% clickList || oldname %in% countList){
    names(data)[names(data) == oldname] <- substr(oldname, 1, nchar(oldname)-6)
  }
  else if(oldname %in% submitList){
    names(data)[names(data) == oldname] <- substr(oldname, 1, nchar(oldname)-7)
  }
}

# Matrix, die alle Variablen enthaelt, die von reshape betroffen sein sollen
MapVars <- rbind(c("Map1.1", "Map1.1_t_First", "Map1.1_t_Last", "Map1.1_t_Page", "Map1.1_t_Click", "Map1.1_corrA"), 
                c("Map1.2", "Map1.2_t_First", "Map1.2_t_Last", "Map1.2_t_Page", "Map1.2_t_Click", "Map1.2_corrA"), 
                c("Map1.3", "Map1.3_t_First", "Map1.3_t_Last", "Map1.3_t_Page", "Map1.3_t_Click", "Map1.3_corrA"),
                c("Map1.4", "Map1.4_t_First", "Map1.4_t_Last", "Map1.4_t_Page", "Map1.4_t_Click", "Map1.4_corrA"),
                c("Map1.5", "Map1.5_t_First", "Map1.5_t_Last", "Map1.5_t_Page", "Map1.5_t_Click", "Map1.5_corrA"),
                c("Map2.1", "Map2.1_t_First", "Map2.1_t_Last", "Map2.1_t_Page", "Map2.1_t_Click", "Map2.1_corrA"),  
                c("Map2.2", "Map2.2_t_First", "Map2.2_t_Last", "Map2.2_t_Page", "Map2.2_t_Click", "Map2.2_corrA"),
                c("Map2.3", "Map2.3_t_First", "Map2.3_t_Last", "Map2.3_t_Page", "Map2.3_t_Click", "Map2.3_corrA"),  
                c("Map2.4", "Map2.4_t_First", "Map2.4_t_Last", "Map2.4_t_Page", "Map2.4_t_Click", "Map2.4_corrA"),   
                c("Map2.5", "Map2.5_t_First", "Map2.5_t_Last", "Map2.5_t_Page", "Map2.5_t_Click", "Map2.5_corrA"),
                c("Map3.1", "Map3.1_t_First", "Map3.1_t_Last", "Map3.1_t_Page", "Map3.1_t_Click", "Map3.1_corrA"),
                c("Map3.2", "Map3.2_t_First", "Map3.2_t_Last", "Map3.2_t_Page", "Map3.2_t_Click", "Map3.2_corrA"),
                c("Map3.3", "Map3.3_t_First", "Map3.3_t_Last", "Map3.3_t_Page", "Map3.3_t_Click", "Map3.3_corrA"),
                c("Map3.4", "Map3.4_t_First", "Map3.4_t_Last", "Map3.4_t_Page", "Map3.4_t_Click", "Map3.4_corrA"),
                c("Map3.5", "Map3.5_t_First", "Map3.5_t_Last", "Map3.5_t_Page", "Map3.5_t_Click", "Map3.5_corrA"),
                c("Verstaendnis1", "Verstaendnis1_t_First", "Verstaendnis1_t_Last", "Verstaendnis1_t_Page", "Verstaendnis1_t_Click", "Verstaendnis1"),
                c("Verstaendnis2", "Verstaendnis2_t_First", "Verstaendnis2_t_Last", "Verstaendnis2_t_Page", "Verstaendnis2_t_Click", "Verstaendnis2"))

# Reshapen der Daten in 'long' Format:
# dabei sollen aus den Variablen aus MapVars 6 neue Variablen mit jeweils 17 Zeilen entstehen
longData <- reshape(data, varying = list(MapVars[,1], MapVars[,2], MapVars[,3], MapVars[,4], MapVars[,5], MapVars[,6]),
                    direction="long", idvar="VPNr", timevar = "Map")

# Sortieren nach Versuchsperson
longData <- longData[order(longData$VPNr),]

# Variablen umbenennen: 
# Map1.1(Antworten); Map1.1_t_First; Map1.1_t_Last; Map1.1_t_Click; Map1.1_t_Page; Map1.1_corrA 
# in => Ans, FirstClick, LastClick, ClickCount, PageSubmit, CorrAns 
names(longData)[names(longData) == "Map1.1"] <- "Ans"
names(longData)[names(longData) == "Map1.1_t_First"] <- "FirstClick"
names(longData)[names(longData) == "Map1.1_t_Last"] <- "LastClick"
names(longData)[names(longData) == "Map1.1_t_Click"] <- "ClickCount"
names(longData)[names(longData) == "Map1.1_t_Page"] <- "PageSubmit"
names(longData)[names(longData) == "Map1.1_corrA"] <- "CorrAns"

# Map Variable beinhaltet Zahlen 1:17 fuer jede Karte einen Wert
# => umcodieren in die vorherige Kartenbezeichnung: 1.1, 1.2, 1.3 ..... 3.5, V1, V2 (fuer die beiden Verstaendnisfragen)
mapList <- c("1.1", "1.2", "1.3", "1.4", "1.5", "2.1", "2.2", "2.3", "2.4", "2.5", "3.1", "3.2", "3.3", "3.4", "3.5", "V1", "V2")
longData$Map <- lapply(longData$Map, function(x){mapList[x]})

# Neue Variable: CorrAnsBool, die ueberprueft, ob die richtige Antwort gegeben wurde:
for(i in 1:length(longData$VPNr)){
  longData$CorrAnsBool[i] <- longData$Ans[i] == longData$CorrAns[i]
}

### new Variables

# correct
longData$correct <- NA
longData$correct[longData$CorrAnsBool==FALSE] <- 0
longData$correct[longData$CorrAnsBool==TRUE] <- 1

# map
longData$map <- NA
longData$map <- substr(longData$Map,1,1)

# question
longData$question <- NA
longData$question <- substr(longData$Map,3,3)

# task (combination of map and question)
longData$task <- NA
longData$task <- substr(longData$Map,1,3)

longData$VPNr <- as.factor(longData$VPNr)
longData$Progress <- as.numeric(as.character(longData$Progress))
longData$Duration..in.seconds. <- as.numeric(as.character(longData$Duration..in.seconds.))
longData$TLX_mental_1 <- as.numeric(as.character(longData$TLX_mental_1))
longData$TLX_time_1 <- as.numeric(as.character(longData$TLX_time_1))
longData$TLX_performance_1 <- as.numeric(as.character(longData$TLX_performance_1))
longData$TLX_effort_1 <- as.numeric(as.character(longData$TLX_effort_1))
longData$TLX_frustration_1 <- as.numeric(as.character(longData$TLX_frustration_1))
longData$Map <- as.factor(unlist(longData$Map))
longData$Ans <- as.numeric(as.character(longData$Ans))
longData$FirstClick <- as.numeric(as.character(longData$FirstClick))
longData$LastClick <- as.numeric(as.character(longData$LastClick))
longData$PageSubmit <- as.numeric(as.character(longData$PageSubmit))
longData$ClickCount <- as.numeric(as.character(longData$ClickCount))
longData$CorrAns <- as.numeric(longData$CorrAns)
longData$map <- as.factor(longData$map)
longData$question <- as.factor(longData$question)
longData$task <- as.factor(longData$task)
longData$verb1 <- as.factor(unlist(longData$verb1))
longData$verb2 <- as.factor(unlist(longData$verb2))
longData$CorrAns1 <- as.factor(unlist(longData$CorrAns1))
longData$CorrAns2 <- as.factor(unlist(longData$CorrAns2))
data$TLX_mental_1 <- as.numeric(as.character(data$TLX_mental_1))
data$TLX_time_1 <- as.numeric(as.character(data$TLX_time_1))
data$TLX_performance_1 <- as.numeric(as.character(data$TLX_performance_1))
data$TLX_effort_1 <- as.numeric(as.character(data$TLX_effort_1))
data$TLX_frustration_1 <- as.numeric(as.character(data$TLX_frustration_1))

# data without first two questions on understanding of notation and only complete data
longData <- longData[longData$map!="V",]
longData <- longData[longData$Progress==100,]
data <- data[data$Progress==100,]

# for analyzing click maps

# click map 1
data$cm1_correct <- NA
data$cm1_correct[data$Reg1_CorrAns==TRUE] <- 1
data$cm1_correct[data$Reg1_CorrAns==FALSE] <- 0

# click map 2
data$cm2_correct <- NA
data$cm2_correct[data$Reg2_CorrAns==TRUE] <- 1
data$cm2_correct[data$Reg2_CorrAns==FALSE] <- 0

# export data for data management and storage
# dput(longData, file="RAWDATA_JAMZ20200527_Artenvielfalt_Tutorialstudie.txt")

# # read data in  
# data2 <- source("RAWDATA_JAMZ20200527_Artenvielfalt_Tutorialstudie.txt")
# View(data2$value)
# data2 <- data2$value
#########################
### get to know your data

# number of participants per condition

View(dat[dat$cond=="tutorial" & dat$Progress==100,]) # 41 (59)
View(dat[dat$cond=="control" & dat$Progress==100,]) # 48 (64)

summary(data[data$Progress==100,])
summary(data$Demo_ges[data$Progress==100])
summary(as.numeric(as.character(data$Demo_alt[data$Progress==100])))
summary(data$Demo_schule[data$Progress==100])
summary(data$ornitho_seit[data$Progress==100])
summary(data$ornitho_Experte[data$Progress==100])
summary(data$ornitho_oft[data$Progress==100])
summary(data$ornitho_Eintrag[data$Progress==100])
summary(data$ornitho_DB[data$Progress==100])
summary(data$ornitho_Fkt[data$Progress==100])
summary(data$KarteDigital1_1[data$Progress==100])
summary(data$KarteDigital1_2[data$Progress==100])
summary(data$KarteDigital2_1[data$Progress==100])
summary(data$KarteDigital2_2[data$Progress==100])
summary(data$KarteAnalog_1[data$Progress==100])
summary(data$KarteAnalog_2[data$Progress==100])
summary(data$KarteNavi_1[data$Progress==100])
summary(data$KartenKompetenz[data$Progress==100])
# plot for visualization, count (plyr) for listing

tapply(longData$correct, longData$cond, mean)
tapply(longData$correct, longData$cond, sd)

tapply(longData$correct, longData$map, mean)
tapply(longData$correct, longData$map, sd)

tapply(longData$correct, longData$task, mean)
tapply(longData$correct, longData$task, sd)

tapply(longData$correct[longData$cond=="tutorial"], longData$task[longData$cond=="tutorial"], mean)
tapply(longData$correct[longData$cond=="control"], longData$task[longData$cond=="control"], mean)

tapply(data$TLX_mental_1, data$cond, mean, na.rm=T)
tapply(data$TLX_time_1, data$cond, mean, na.rm=T)
tapply(data$TLX_performance_1, data$cond, mean, na.rm=T)
tapply(data$TLX_effort_1, data$cond, mean, na.rm=T)
tapply(data$TLX_frustration_1, data$cond, mean, na.rm=T)

t.test(data$TLX_mental_1 ~ data$cond)
t.test(data$TLX_time_1 ~ data$cond)
t.test(data$TLX_performance_1 ~ data$cond)
t.test(data$TLX_effort_1 ~ data$cond)
t.test(data$TLX_frustration_1 ~ data$cond)


###################################################
#<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>#
#                    ANALYSIS                     #
#<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>#
###################################################

########################################
#                                      #
##### 3 maps with 5 questions each #####
#                                      #
########################################

######################
# PROPORTION CORRECT #
######################

glmer_PC <- glmer(correct ~ cond  + (1|VPNr) + (1|task), family=binomial, longData)
summary(glmer_PC)
mixed_PC <- mixed(correct ~ cond + (1|VPNr) + (1|task), family=binomial, method="LRT", longData)
anova(mixed_PC)

emm_PC_1 <- emmeans(glmer_PC, pairwise ~ cond)
logit2prob(0.3176762) # control
logit2prob(0.1782899) # tutorial

logit2prob(0.4905411)
logit2prob(0.4957598)

### pool answers "information not given on map" and "is not true"
longData$CorrAns_pooled <- NA
longData$CorrAns_pooled[longData$Ans ==1 & longData$CorrAns==1] <- 1
longData$CorrAns_pooled[longData$Ans ==2 & longData$CorrAns==2] <- 1
longData$CorrAns_pooled[longData$Ans ==3 & longData$CorrAns==3] <- 1
longData$CorrAns_pooled[longData$Ans ==1 & longData$CorrAns==2] <- 0
longData$CorrAns_pooled[longData$Ans ==1 & longData$CorrAns==3] <- 0
longData$CorrAns_pooled[longData$Ans ==2 & longData$CorrAns==1] <- 0
longData$CorrAns_pooled[longData$Ans ==2 & longData$CorrAns==3] <- 1
longData$CorrAns_pooled[longData$Ans ==3 & longData$CorrAns==1] <- 0
longData$CorrAns_pooled[longData$Ans ==3 & longData$CorrAns==2] <- 1
longData$CorrAns_pooled[longData$Ans ==4 ] <- 0

glmer_PC_pooled <- glmer(CorrAns_pooled ~ cond + (1|VPNr) + (1|task), family=binomial, longData)
summary(glmer_PC_pooled)
Anova(glmer_PC_pooled, type=3)
mixed_PC_pooled <- mixed(CorrAns_pooled ~ cond + (1|VPNr) + (1|task), family=binomial, method="LRT", longData)
anova(mixed_PC_pooled)

# with pooled answers yes/no and info not given
longData$yn_ni[longData$CorrAns=="1"|longData$CorrAns=="2"] <- "1"
longData$yn_ni[longData$CorrAns=="3"] <- "0"

dat_aov_2 <- aggregate(list(correct=longData$correct), list(yn_ni=longData$yn_ni,VPNr=longData$VPNr),mean, na.rm=T)
tapply(dat_aov_2$correct, dat_aov_2$yn_ni, mean, na.rm=T)
tapply(dat_aov_2$correct, dat_aov_2$yn_ni, sd, na.rm=T)

mixed_PC_ni <- mixed(correct ~ cond * yn_ni + (1|VPNr) + (1|task), family=binomial, method="LRT", longData)
anova(mixed_PC_ni)

emm_PC_ni1 <- emmeans(mixed_PC_ni, pairwise ~ cond | yn_ni)
emm_PC_ni2 <- emmeans(mixed_PC_ni, pairwise ~ yn_ni)

logit2prob(-0.4855849)
logit2prob(-0.7644233)
logit2prob(1.8848500)
logit2prob(2.1823383)

logit2prob(0.4434292)
logit2prob(0.4511975)
logit2prob(0.6219359)
logit2prob(0.6358148)

logit2prob(-0.6250041)
logit2prob(2.0335942)

logit2prob(0.4261079)
logit2prob(0.6021431)


### ANOVA ###
dat_aov <- aggregate(list(correct=longData$correct), list(cond=longData$cond,VPNr=longData$VPNr),mean, na.rm=T)
aov <- aov(correct ~ cond , dat_aov)
summary(aov)

########################
#                      #
##### 2 click maps #####
#                      #
########################

######################
# PROPORTION CORRECT #
######################

# Clickmap 1 # 

glm_CM_1 <- glm(cm1_correct ~ cond , family=binomial, data[data$Progress==100,])
summary(glm_CM_1)
Anova(glm_CM_1, type=3)

emm_CM_1 <- emmeans(glm_CM_1, pairwise ~ cond)
logit2prob(-0.1670541) # control
logit2prob(0.8823892) # tutorial

logit2prob(0.2896827)
logit2prob(0.3432433)


# Clickmap 2 #

glm_CM_2 <- glm(cm2_correct ~ cond , family=binomial, data)
summary(glm_CM_2)
Anova(glm_CM_2, type=3)

emm_CM_2 <- emmeans(glm_CM_2, pairwise ~ cond)
logit2prob(2.397895) # control
logit2prob(2.970414) # tutorial

logit2prob(0.5222329)
logit2prob(0.7249655)

####################################
# SUBJECTIVE EVALUATION (NASA TLX) #
####################################

newData$TLX_mental_1 <- as.numeric(as.character(newData$TLX_mental_1))
newData$TLX_time_1 <- as.numeric(as.character(newData$TLX_time_1))
newData$TLX_performance_1 <- as.numeric(as.character(newData$TLX_performance_1))
newData$TLX_effort_1 <- as.numeric(as.character(newData$TLX_effort_1))
newData$TLX_frustration_1 <- as.numeric(as.character(newData$TLX_frustration_1))

# mental demand
lm_TLX_1 <- lm(TLX_mental_1 ~ cond, newData[newData$Progress==100,])
summary(lm_TLX_1)
emm_TLX_1_1 <- emmeans(lm_TLX_1, pairwise ~ cond)
emm_TLX_1_1
# time demand
lm_TLX_2 <- lm(TLX_time_1 ~ cond, newData[newData$Progress==100,])
summary(lm_TLX_2)
emm_TLX_2_1 <- emmeans(lm_TLX_2, pairwise ~ cond)
emm_TLX_2_1
# performance
lm_TLX_3 <- lm(TLX_performance_1 ~ cond, newData[newData$Progress==100,])
summary(lm_TLX_3)
emm_TLX_3_1 <- emmeans(lm_TLX_3, pairwise ~ cond)
emm_TLX_3_1
# effort
lm_TLX_4 <- lm(TLX_effort_1 ~ cond, newData[newData$Progress==100,])
summary(lm_TLX_4)
emm_TLX_4_1 <- emmeans(lm_TLX_4, pairwise ~ cond)
emm_TLX_4_1
# frustration
lm_TLX_5 <- lm(TLX_frustration_1 ~ cond, newData[newData$Progress==100,])
summary(lm_TLX_5)
emm_TLX_5_1 <- emmeans(lm_TLX_5, pairwise ~ cond)
emm_TLX_5_1
