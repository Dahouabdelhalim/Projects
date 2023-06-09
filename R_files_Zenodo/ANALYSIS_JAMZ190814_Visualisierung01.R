### Analysis Visualisierung_1
### Project Artenvielfalt erleben ###

# set working directory
setwd("/Users/.../data")
rm(list=ls())

library(dplyr)
library(lme4)
library(car)
library(emmeans)
library(afex)
library(tidyr)
source("https://sebastiansauer.github.io/Rcode/logit2prob.R")

# read data
dat_B <- read.csv("ORIGINALDATA_JAMZ190814_Artenvielfalt_Visualisierung01/ornitho-Visualization1_B_August 14, 2019_15.14.csv", header=T)
dat_E <- read.csv("ORIGINALDATA_JAMZ190814_Artenvielfalt_Visualisierung01/ornitho-Visualization1_E_August 14, 2019_15.13.csv", header=T)

# code expertise level
dat_B$expertise <- "beginner" 
dat_E$expertise <- "expert" 

# bind data (expert & beginner)
# dat <- rbind(dat_B[-c(1:2),], dat_E[-c(1:2),]) # we still need the information of line 3 for further processing
dat <- rbind(dat_B, dat_E[-c(1:2),])

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
Rueck <- subset(dat, dat$Rückzug2 == 1 | dat$Einwilligung != 1, select = c("ResponseId", "VPNr")) # List with VP who wanted to retract their data or did not confirm the consent
dat <- subset(dat, dat$Einwilligung == 1) # Einwilligung mit 1
dat <- subset(dat, dat$Rückzug2 != 1) # Datenrueckzug mit Antwort 1 der zweiten Frage 
dat$Einwilligung <- NULL
dat$Rückzug1 <- NULL
dat$Rückzug2 <- NULL


########## Die beiden Konditionen (dreistufig & fuenfstufig) - codieren und Variablen untereinander setzten ##########
# 1. Zwei Datensaetze (dat3, dat5) jeweils fuer eine Kondition erstellen, die jeweils nur die Variablen enthalten, die sich unterscheiden
# 2. jeweils in den beiden Datensaetzen die leeren Felder herausloeschen
# 3. Variable fuer die Kondition erstellen
# 4. die beiden Datensaetze (dat3, dat5) zusammenfuegen (dat35)
# 5. Datensatz (dat0) erstellen, der nur Variablen beinhaltet, die fuer alle Versuchspersonen gleich sind
# 6. aus dat0 und dat35 einen finalen Datensatz (newData) zusammenfuegen

dat3 <- select(dat, c(starts_with("X3"), "ResponseId", "VPNr")) # Spalten auswaehlen deren Name mit X3 beginnt, sowie ResponseId und VPNr, um die Daten spaeter identifizieren zu koennen.
dat3 <- subset(dat3, !is.na(as.numeric(levels(dat3$X3Legende_t_Page.Submit))[dat3$X3Legende_t_Page.Submit])) # Zeilen mit Kondition X5 loeschen, NA warnings ok
dat3$X3cond <- 3 # Neue Variable mit der Kondition einfuehren
dat3$X3cond <- as.factor(dat3$X3cond)
dat5 <- select(dat, c(starts_with("X5"), "ResponseId", "VPNr"))
dat5 <- subset(dat5, !is.na(as.numeric(levels(dat5$X5Legende_t_Page.Submit))[dat5$X5Legende_t_Page.Submit])) # NA warnings ok
dat5$X5cond <- 5
dat5$X5cond <- as.factor(dat5$X5cond)

### Datensatz ohne Karten-Fragen
dat0 <- select(dat, -c(starts_with("X3"), starts_with("X5")))

names(dat3)[names(dat3) != "ResponseId" & names(dat3) != "VPNr"] <- substring(names(dat3)[names(dat3) != "ResponseId" & names(dat3) != "VPNr"],3,10000) # Entfernen der ersten beiden Characters
names(dat5)[names(dat5) != "ResponseId" & names(dat5) != "VPNr"] <- substring(names(dat5)[names(dat5) != "ResponseId" & names(dat5) != "VPNr"],3,10000) # dadurch ergeben sich gleiche Namen fuer beide Daten
dat35 <- rbind(dat3, dat5) # Daten von beiden Bedingungen zusammenfuehren

########## Klickfelder reduzieren auf jeweils 6 Variablen #############

names(dat35)[names(dat35) == 'Reg1_12'] <- '1,2,8' # 1,2,7 und 2,2,7 sind doppelt (2,8 fehlt jeweils)Reg1/2_12 eigentlich 2,8
names(dat35)[names(dat35) == 'Reg2_12'] <- '2,2,8'

RegtList1 <- names(select(dat35, starts_with("Reg1_t")))
RegtList2 <- names(select(dat35, starts_with("Reg2_t")))
RegList1 <- names(select(select(dat35, starts_with("Reg1")), -RegtList1)) # Liste, die alle Namen der Klickfelder beinhalten
RegList2 <- names(select(select(dat35, starts_with("Reg2")), -RegtList2))

# Variablen nach ihrer Matrix-Koordinate benennen
for (i in 1:length(names(dat35))) { # Iteration durch aller Variablen in dat35
  oldname <- names(dat35)[i] # aktueller Variablenname
  if(oldname %in% RegList1){ # nur Variablen waehlen, die ein Klickfeld beschreiben (hier fuer Region1)
    Erklaerung <- levels(dat35[,i])[3] # Auswahl der Beschreibung/Erklaerung der Variable
    MatrixPos <- substr(Erklaerung, nchar(Erklaerung)-2, nchar(Erklaerung)) # die letzten 3 Characters beinhalten die Matrixposition (Achtung die Reihe 10 wird zu 0)
    names(dat35)[names(dat35) == oldname] <- paste('1', MatrixPos, sep = ',') # Variable neu benennen nach ihrer Matrix-Koordinate, mit einer 1 fuer die Region-Karte vorangestellt
  }
  else if(oldname %in% RegList2){ # gleiches Vorgehen, nur fuer Region2
    Erklaerung <- levels(dat35[,i])[3]
    MatrixPos <- substr(Erklaerung, nchar(Erklaerung)-2, nchar(Erklaerung))
    names(dat35)[names(dat35) == oldname] <- paste('2', MatrixPos, sep = ',')
  }
}


dat35$FieldList <- list(character()) # neue Variable mit einer Liste fuer jede VP anlegen
MatrixFields <- names(select(dat35, c(starts_with("1,"), starts_with("2,")))) # Liste aller Variablennamen,die ein Klickfeld benennen

# Erstellen der Listen fuer jede VP
for (i in 1:length(names(dat35))){ # Iteration durch jede Variable in dat35
  name <- names(dat35)[i] # Name der aktuellen Variable
  if(name %in% MatrixFields){ # Nur Variablen waehlen, bei denen es sich um Klickfelder handelt
    for (j in 1:length(dat35$VPNr)) { # Iteration durch die VP
      if(dat35[j,i]=="On"){ # Auswahl der Klickfelder, die ausgewaehlt wurden
        dat35$FieldList[[j]] <- c(dat35$FieldList[[j]], name) # Anhaengen der ausgewaehlten Felder an FieldList der VP
      }
    }
  }
}

# Ausschliessen der Clickfelder aus dem Datensatz
dat35 <- select(dat35, -MatrixFields)

# neue Variablen mit den x- & y-Matrix-Koordinaten von jedem gewaehlte Klickfeld, und fuer Region-Karte 1 & 2
for (i in 1:length(dat35$VPNr)){ # Iteration fuer jede VP
  Reg1x1 <- substring(dat35$FieldList[[i]][1], 3, 3) # FieldList enthaelt strings ("1,2,3" [region,x-Koordinate, y-Koordinate]); 
  Reg1y1 <- substring(dat35$FieldList[[i]][1], 5, 5) # der 5. string bezeichnet die y-Koordinate
  Reg1x2 <- substring(dat35$FieldList[[i]][2], 3, 3) # der 3. string bezeichnet die x-Koordinate
  Reg1y2 <- substring(dat35$FieldList[[i]][2], 5, 5) # insgesamt enthaelt jede Liste aus FieldList 6 strings mit jeweils einer Matrix-Koordinate
  Reg1x3 <- substring(dat35$FieldList[[i]][3], 3, 3)
  Reg1y3 <- substring(dat35$FieldList[[i]][3], 5, 5)
  dat35$Reg1[i] <- list(c(Reg1x1, Reg1y1, Reg1x2, Reg1y2, Reg1x3, Reg1y3))
  
  Reg2x1 <- substring(dat35$FieldList[[i]][4], 3, 3)
  Reg2y1 <- substring(dat35$FieldList[[i]][4], 5, 5)
  Reg2x2 <- substring(dat35$FieldList[[i]][5], 3, 3)
  Reg2y2 <- substring(dat35$FieldList[[i]][5], 5, 5)
  Reg2x3 <- substring(dat35$FieldList[[i]][6], 3, 3)
  Reg2y3 <- substring(dat35$FieldList[[i]][6], 5, 5)
  dat35$Reg2[i] <- list(c(Reg2x1, Reg2y1, Reg2x2, Reg2y2, Reg2x3, Reg2y3))
}

#Ausschliessen der Variable mit den Listen der 6 Matrixkoordinaten
dat35$FieldList <- NULL

# Zusammenfuehren dat35, dat0
newData <- merge(dat35, dat0, by=c("VPNr", "ResponseId")) # newData enthaelt 132 Variablen

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
newData$Reg1 <- lapply(newData$Reg1, as.numeric)
newData$Reg1 <- lapply(newData$Reg1, function(x) ifelse(x==0, 10, x)) # die 0 bedeutet Reihe 10
newData$verb1 <- lapply(newData$Reg1, verbunden)

newData$Reg2 <- lapply(newData$Reg2, as.numeric)
newData$Reg2 <- lapply(newData$Reg2, function(x) ifelse(x==0, 10, x))
newData$verb2 <- lapply(newData$Reg2, verbunden)

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

# Funktion gibt an, ob alle drei Felder keine Sichtung enthalten und die Insitaet ueber 1 liegt
# Richtige Antwort fuer Reg1: Alle drei Felder duerfen keine Sichtung enthalten, die Intensitaet muss groesser 1 sein
CorrIntens <- function(l){
  return((intensitaetsM[l[1],l[2]] > 1 && intensitaetsM[l[3],l[4]] > 1 && intensitaetsM[l[5],l[6]] > 1) &&
           (sichtungsM[l[1],l[2]] == 0 && sichtungsM[l[3],l[4]] == 0 && sichtungsM[l[5],l[6]] == 0))
}

# neue Variable, die angibt, fuer Region 1 inhaltlich richtige Felder gewaehlt wurden
newData$CorrAns1 <- lapply(newData$Reg1, CorrIntens)

# neue Variable Reg1_CorrAns, die angibt, ob sowohl die Felder nebeneinander liegen, als auch die richtigen Felder gewaehlt wurden
for (i in 1:length(newData$VPNr)){
  newData$Reg1_CorrAns[i] <- newData$verb1[[i]] && newData$CorrAns1[[i]]
}

# Funktion gibt an, ob alle drei Felder eine Sichtung enthalten,
# Richtige Antwort fuer Reg2: Alle drei Felder muessen eine Sichtung enthalten
CorrSichtung <- function(l){
  return(sichtungsM[l[1],l[2]] == 1 && sichtungsM[l[3],l[4]] == 1 && sichtungsM[l[5],l[6]] == 1) 
}
# neue Variable, die angibt, ob fuer Region 2 inhaltlich richtige Felder gewaehlt wurden
newData$CorrAns2 <- lapply(newData$Reg2, CorrSichtung)

# neue Variable Reg1_CorrAns, die angibt, ob sowohl die Felder nebeneinander liegen, als auch die richtigen Felder gewaehlt wurden
for(i in 1:length(newData$VPNr)){
  newData$Reg2_CorrAns[i] <- newData$verb2[[i]] && newData$CorrAns2[[i]]
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
      newData$new <- corrA
      colnames(newData)[colnames(newData) == 'new'] <- var
      l <- l+1
  }
}

# Umbennenen der Variablennamen ohne ".Click", ".Submit", ".Count"
clickList <- names(select(newData, ends_with(".Click")))
submitList <- names(select(newData, ends_with(".Submit")))
countList <- names(select(newData, ends_with(".Count")))
nameList <- names(newData)

for (i in 1:length(nameList)) { # Iteration durch alle Variablen in newData
  oldname <- nameList[i] # aktueller Variablenname
  if(oldname %in% clickList || oldname %in% countList){
    names(newData)[names(newData) == oldname] <- substr(oldname, 1, nchar(oldname)-6)
  }
  else if(oldname %in% submitList){
    names(newData)[names(newData) == oldname] <- substr(oldname, 1, nchar(oldname)-7)
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
                c("Verständnis1", "Verständnis1_t_First", "Verständnis1_t_Last", "Verständnis1_t_Page", "Verständnis1_t_Click", "Verständnis1"),
                c("Verständnis2", "Verständnis2_t_First", "Verständnis2_t_Last", "Verständnis2_t_Page", "Verständnis2_t_Click", "Verständnis2"))

# Reshapen der Daten in 'long' Format:
# dabei sollen aus den Variablen aus MapVars 6 neue Variablen mit jeweils 17 Zeilen entstehen
longData <- reshape(newData, varying = list(MapVars[,1], MapVars[,2], MapVars[,3], MapVars[,4], MapVars[,5], MapVars[,6]),
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
longData$TLX1_1 <- as.numeric(as.character(longData$TLX1_1))
longData$TLX2_1 <- as.numeric(as.character(longData$TLX2_1))
longData$TLX3_1 <- as.numeric(as.character(longData$TLX3_1))
longData$TLX4_1 <- as.numeric(as.character(longData$TLX4_1))
longData$TLX5_1 <- as.numeric(as.character(longData$TLX5_1))
longData$expertise <- as.factor(longData$expertise)
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

# data without first two questions on understanding of notation and only complete data
data <- longData[longData$map!="V",]
data <- data[data$Progress==100,]

# for analyzing click maps

# click map 1
newData$cm1_correct <- NA
newData$cm1_correct[newData$Reg1_CorrAns==TRUE] <- 1
newData$cm1_correct[newData$Reg1_CorrAns==FALSE] <- 0

# click map 2
newData$cm2_correct <- NA
newData$cm2_correct[newData$Reg2_CorrAns==TRUE] <- 1
newData$cm2_correct[newData$Reg2_CorrAns==FALSE] <- 0

# export data fpr data management and storage
# write.table(data[,-c(16,17)], "RAWDATA_JAMZ190814_Artenvielfalt_Visualisierung01.csv", sep="\\t", quote=F, row.names=F)

#########################
### get to know your data

# number of participants per condition
View(newData[newData$expertise=="beginner" & newData$cond=="3" & newData$Progress==100,]) # 43 finished (59 total)
View(newData[newData$expertise=="beginner" & newData$cond=="5" & newData$Progress==100,]) # 41 (60)
View(newData[newData$expertise=="expert" & newData$cond=="3" & newData$Progress==100,]) # 46 (58)
View(newData[newData$expertise=="expert" & newData$cond=="5" & newData$Progress==100,]) # 43 (58)

summary(newData[newData$Progress==100,])
summary(newData$Demo_ges[newData$Progress==100])
summary(as.numeric(as.character(newData$Demo_alt[newData$Progress==100])))
summary(newData$Demo_schule[newData$Progress==100])
summary(newData$ornitho_seit[newData$Progress==100])
summary(newData$ornitho_Experte[newData$Progress==100])
summary(newData$ornitho_oft[newData$Progress==100])
summary(newData$ornitho_Eintrag[newData$Progress==100])
summary(newData$ornitho_DB[newData$Progress==100])
summary(newData$ornitho_Fkt[newData$Progress==100])
summary(newData$KarteDigital1_1[newData$Progress==100])
summary(newData$KarteDigital1_2[newData$Progress==100])
summary(newData$KarteDigital2_1[newData$Progress==100])
summary(newData$KarteDigital2_2[newData$Progress==100])
summary(newData$KarteAnalog_1[newData$Progress==100])
summary(newData$KarteAnalog_2[newData$Progress==100])
summary(newData$KarteNavi_1[newData$Progress==100])
summary(newData$KartenKompetenz[newData$Progress==100])
# plot for visualization, count (plyr) for listing

tapply(data$correct, data$expertise, mean)

tapply(data$correct, data$cond, mean)
tapply(data$correct[data$expertise=="beginner"], data$cond[data$expertise=="beginner"], mean)
tapply(data$correct[data$expertise=="beginner"], data$cond[data$expertise=="beginner"], sd)
tapply(data$correct[data$expertise=="expert"], data$cond[data$expertise=="expert"], mean)
tapply(data$correct[data$expertise=="expert"], data$cond[data$expertise=="expert"], sd)

tapply(data$correct, data$map, mean)
tapply(data$correct[data$expertise=="beginner"], data$map[data$expertise=="beginner"], mean)
tapply(data$correct[data$expertise=="expert"], data$map[data$expertise=="expert"], mean)

tapply(data$correct, data$task, mean)
tapply(data$correct[data$expertise=="beginner"], data$task[data$expertise=="beginner"], mean)
tapply(data$correct[data$expertise=="expert"], data$task[data$expertise=="expert"], mean)

tapply(data$correct[data$cond=="3"], data$task[data$cond=="3"], mean)
tapply(data$correct[data$cond=="5"], data$task[data$cond=="5"], mean)

tapply(data$correct[data$expertise=="beginner" & data$cond=="3"], data$task[data$expertise=="beginner" & data$cond=="3"], mean)
tapply(data$correct[data$expertise=="beginner" & data$cond=="5"], data$task[data$expertise=="beginner" & data$cond=="5"], mean)
tapply(data$correct[data$expertise=="expert" & data$cond=="3"], data$task[data$expertise=="expert" & data$cond=="3"], mean)
tapply(data$correct[data$expertise=="expert" & data$cond=="5"], data$task[data$expertise=="expert" & data$cond=="5"], mean)

tapply(data$TLX1_1, data$expertise, mean, na.rm=T)
tapply(data$TLX2_1, data$expertise, mean, na.rm=T)
tapply(data$TLX3_1, data$expertise, mean, na.rm=T)
tapply(data$TLX4_1, data$expertise, mean, na.rm=T)
tapply(data$TLX5_1, data$expertise, mean, na.rm=T)

tapply(data$TLX1_1, data$cond, mean, na.rm=T)
tapply(data$TLX2_1, data$cond, mean, na.rm=T)
tapply(data$TLX3_1, data$cond, mean, na.rm=T)
tapply(data$TLX4_1, data$cond, mean, na.rm=T)
tapply(data$TLX5_1, data$cond, mean, na.rm=T)

plot(tapply(data$correct, data$VPNr, mean, na.rm=T))
hist(tapply(data$correct, data$VPNr, mean, na.rm=T),breaks=20, xlim=c(0,1))

# overview tasks
tapply(data$Ans, data$task, hist, breaks=c(0,1,2,3,4), xlim=c(0,4), ylim=c(0,180))

# click maps #
# CM 1
summary(newData$cm1_correct[newData$Progress==100])
tapply(newData$cm1_correct[newData$Progress==100], newData$cond[newData$Progress==100], mean)
tapply(newData$cm1_correct[newData$Progress==100], newData$expertise[newData$Progress==100], mean)

tapply(newData$cm1_correct[newData$Progress==100 & newData$cond=="3"], newData$expertise[newData$Progress==100 & newData$cond=="3"], mean)
tapply(newData$cm1_correct[newData$Progress==100 & newData$cond=="5"], newData$expertise[newData$Progress==100 & newData$cond=="5"], mean)

# CM 2
summary(newData$cm2_correct[newData$Progress==100])
tapply(newData$cm2_correct[newData$Progress==100], newData$cond[newData$Progress==100], mean)
tapply(newData$cm2_correct[newData$Progress==100], newData$expertise[newData$Progress==100], mean)

tapply(newData$cm2_correct[newData$Progress==100 & newData$cond=="3"], newData$expertise[newData$Progress==100 & newData$cond=="3"], mean)
tapply(newData$cm2_correct[newData$Progress==100 & newData$cond=="5"], newData$expertise[newData$Progress==100 & newData$cond=="5"], mean)

# gender

tapply(data$correct, data$Demo_ges, mean, na.rm=T)


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

glmer_PC <- glmer(correct ~ cond * expertise + (1|VPNr) + (1|task), family=binomial, data)
mixed_PC <- mixed(correct ~ cond * expertise + (1|VPNr) + (1|task), family=binomial, method="LRT", data)
anova(mixed_PC)

emm_PC_1 <- emmeans(glmer_PC, pairwise ~ expertise | cond)
logit2prob(0.3323165)
logit2prob(0.1734866)
logit2prob(0.0252538)
logit2prob(0.1196049)

logit2prob(0.4710691)
logit2prob(0.4686440)
logit2prob(0.4725326)
logit2prob(0.4709054)

emm_PC_2 <- emmeans(glmer_PC, pairwise ~ cond | expertise)
emm_PC_3 <- emmeans(glmer_PC, pairwise ~ cond)
emm_PC_4 <- emmeans(glmer_PC, pairwise ~ expertise)

### pool answers "information not given on map" and "is not true"
data$CorrAns_pooled <- NA
data$CorrAns_pooled[data$Ans ==1 & data$CorrAns==1] <- 1
data$CorrAns_pooled[data$Ans ==2 & data$CorrAns==2] <- 1
data$CorrAns_pooled[data$Ans ==3 & data$CorrAns==3] <- 1
data$CorrAns_pooled[data$Ans ==1 & data$CorrAns==2] <- 0
data$CorrAns_pooled[data$Ans ==1 & data$CorrAns==3] <- 0
data$CorrAns_pooled[data$Ans ==2 & data$CorrAns==1] <- 0
data$CorrAns_pooled[data$Ans ==2 & data$CorrAns==3] <- 1
data$CorrAns_pooled[data$Ans ==3 & data$CorrAns==1] <- 0
data$CorrAns_pooled[data$Ans ==3 & data$CorrAns==2] <- 1
data$CorrAns_pooled[data$Ans ==4 ] <- 0

glmer_PC_pooled <- glmer(CorrAns_pooled ~ cond * expertise + (1|VPNr) + (1|task), family=binomial, data)
summary(glmer_PC_pooled)
Anova(glmer_PC_pooled, type=3)
mixed_PC_pooled <- mixed(CorrAns_pooled ~ cond * expertise + (1|VPNr) + (1|task), family=binomial, method="LRT", data)
anova(mixed_PC_pooled)

# with pooled answers yes/no and info not given
data$yn_ni[data$CorrAns=="1"|data$CorrAns=="2"] <- "1"
data$yn_ni[data$CorrAns=="3"] <- "0"

mixed_PC_ni <- mixed(correct ~ cond * expertise + yn_ni + (1|VPNr) + (1|task), family=binomial, method="LRT", data)
anova(mixed_PC_ni)

emm_PC_ni1 <- emmeans(mixed_PC_ni, pairwise ~ expertise | yn_ni)
emm_PC_ni2 <- emmeans(mixed_PC_ni, pairwise ~ yn_ni)

logit2prob(-0.5788241)
logit2prob(-0.7488230)
logit2prob(1.6425601)
logit2prob(2.0087967)

logit2prob(0.4058920)
logit2prob(0.4050451)
logit2prob(0.5689254)
logit2prob(0.5706098)

logit2prob(-0.6638236)
logit2prob(1.8256784)

logit2prob(0.3935908)
logit2prob(0.5555788)

### ANOVA ###
dat_aov <- aggregate(list(correct=data$correct), list(cond=data$cond,expertise=data$expertise, VPNr=data$VPNr),mean, na.rm=T)
aov <- aov(correct ~ cond * expertise, dat_aov)
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

glm_CM_1 <- glm(cm1_correct ~ cond * expertise, family=binomial, newData[newData$Progress==100,])
summary(glm_CM_1)

emm_CM_1 <- emmeans(glm_CM_1, pairwise ~ expertise | cond)
logit2prob(-0.0465200)
logit2prob(0.0870114)
logit2prob(-0.3448405)
logit2prob(0.0465200)

logit2prob(0.3050796)
logit2prob(0.2951630)
logit2prob(0.3169889)
logit2prob(0.3050796)

tapply(newData$cm1_correct[newData$Progress==100], newData$cond[newData$Progress==100], mean)
tapply(newData$cm1_correct[newData$Progress==100], newData$expertise[newData$Progress==100], mean)

# Clickmap 2 #

glm_CM_2 <- glm(cm2_correct ~ cond * expertise, family=binomial, newData[newData$Progress==100,])
summary(glm_CM_2)

emm_CM_2 <- emmeans(glm_CM_2, pairwise ~ expertise | cond)

logit2prob(1.475907)
logit2prob(2.104134)
logit2prob(1.974081)
logit2prob(19.566069)

logit2prob(0.3919)
logit2prob(0.4737)
logit2prob(0.4773)
logit2prob(1639.9716)


####################################
# SUBJECTIVE EVALUATION (NASA TLX) #
####################################

newData$TLX1_1 <- as.numeric(as.character(newData$TLX1_1))
newData$TLX2_1 <- as.numeric(as.character(newData$TLX2_1))
newData$TLX3_1 <- as.numeric(as.character(newData$TLX3_1))
newData$TLX4_1 <- as.numeric(as.character(newData$TLX4_1))
newData$TLX5_1 <- as.numeric(as.character(newData$TLX5_1))

# mental demand
lm_TLX_1 <- lm(TLX1_1 ~ cond * expertise, newData[newData$Progress==100,])
summary(lm_TLX_1)
emm_TLX_1_1 <- emmeans(lm_TLX_1, pairwise ~ cond | expertise)
emm_TLX_1_1
# time demand
lm_TLX_2 <- lm(TLX2_1 ~ cond * expertise, newData[newData$Progress==100,])
summary(lm_TLX_2)
emm_TLX_2_1 <- emmeans(lm_TLX_2, pairwise ~ cond | expertise)
emm_TLX_2_1
# performance
lm_TLX_3 <- lm(TLX3_1 ~ cond * expertise, newData[newData$Progress==100,])
summary(lm_TLX_3)
emm_TLX_3_1 <- emmeans(lm_TLX_3, pairwise ~ cond | expertise)
emm_TLX_3_1
# effort
lm_TLX_4 <- lm(TLX4_1 ~ cond * expertise, newData[newData$Progress==100,])
summary(lm_TLX_4)
emm_TLX_4_1 <- emmeans(lm_TLX_4, pairwise ~ cond | expertise)
emm_TLX_4_1
# frustration
lm_TLX_5 <- lm(TLX5_1 ~ cond * expertise, newData[newData$Progress==100,])
summary(lm_TLX_5)
emm_TLX_5_1 <- emmeans(lm_TLX_5, pairwise ~ cond | expertise)
emm_TLX_5_1
