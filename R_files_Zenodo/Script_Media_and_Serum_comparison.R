#initialise les librairies -----------------------------------

rm(list=ls()) # Supprime les données présentes dans la ram et la rom de RStudio

packages<-c("here","sp") # Chemin du repetoire avec les donnees 


for(i in packages){if(!require(i,character.only = T)){install.packages(i,dep=T)}}
lapply(packages, library, character.only = T) #chargement de la liste des packages
library(rlang) # Charge la librairie rlang
library(ggplot2) #Charge la librairie ggplot2 (pour tracer les graphiques)
library(scales) #Charge la librairie scales
library(dplyr) #Charge la librairie dplyr
library(concatenate) #Charge la librairie concatenate
library(writexl) #Charge la librairie writexl
library(tidyr) #Charge la librairie tidyr
library(na.tools) #Charge la librairie na.tools
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggpubr)
library (psych)
library(multcomp)

# fonction  -----------------------------------
opentxt<-function(x) # Fonction utilisée pour ouvrir les fichiers textes
{
  data<-read.table(here("data",x),sep='\\t', header=T) # Ouvre et lit les fichiers txt pour les enregistrer dans une dataframe
}
txt<-dir("data",pattern = ".txt") #Enregistre le nombre de fichier texte et leur nom dans un tableau 

N<-length(txt) # Enregistre dans la variable N le nombre total de fichier texte
t_jours<-NA 
data_real<-structure(list(character()), class = "data.frame") #Création de la dataframe data
Jmax_table<-NA # création de la variable Jmax_table avec comme valeur par défaut NA

for (i in 1:N) # Ouvre tous les fichiers textes et l'enregistre dans un les dataframes data et data_real
{
  data_real<-rbind(data_real,assign(paste("data", i, sep = ""),opentxt(txt[i]))) # ouvre les fichiers txt et les enregistres dans data
  x<-assign(paste("data", i, sep = ""),opentxt(txt[i])) #pendant la boucle, enregistre le dernier fichier txt lu dans la dataframe x
  Jmax_table[i]<-max(x$j) # Récupère la durée maximale du dernier fichier txt afin de connaitre la durée maximale de toutes les expériences 
}
Jmax<-min(Jmax_table) # Enregistre dans cette variable de durée minimale de toutes les expériences

data<-data_real # Enregistre dans la dataframe data_real dans data

# fonction  -----------------------------------

N_of_variable<-function(Variable) # Fonction utilisée pour enregistrer dans un tableau l’ensemble des différentes valeurs existantes dans la colonne choisie de la dataframe en entrée de fonction (ici, Variable) (exemple : nom des milieux de culture, les jours pour chaque contrôles, …)
{
  compteur<-1 #initialise le compteur à la valeur 1 
  table<-NA # Création du tableau 
  table[compteur]<--1 # donne la valeur -1 a la première colonne du tableau
  Num<-NA # Création de la variable Num avec comme valeur NA
  exist<-0 # Création de la variable exist à 0
  for (i in 1:length(Variable)) # Pour chaque ligne contenue dans Variable  
  {
    Num<-Variable[i] # Enregistre dans Num la valeur contenue dans la ligne i
    for (j in 1:length(table)) # de la première colonne a la dernière colonne du tableau table
    {
      if(Num==table[j]) # vérifie si la valeur testée à la ligne j est déjà présente dans la table
      {
        j<-length(table) # si oui, arrêt la boucle for
        exist<-1 #indique que la valeur est déjà contenue dans la table
      }
    }
    if(exist==0) # si la valeur n'existe pas dans la table
    {
      table[compteur]<-Num #enregistre la valeur dans la colonne « compteur » de la table 
      compteur<-compteur+1 #change de colonne « compteur » de la table 
    }
    else(exist<-0) # Sinon, réinitialise la valeur exist a 0 pour la lecture de la ligne suivante
  }
  NVariableNumber<-(sort(table)) # range les valeurs contenues dans le tableau dans l'ordre croissant
  return(NVariableNumber) # renvoi les données de la table en sortie de la fonction.
}


# code  -----------------------------------

data$NDC<-data$NST-data$NDC # Pour chaque point de mesure, calcul le nombre de couples
data$NDC<-data$NDC*100/10 # Pour chaque point de mesure, calcul le pourcentage de couples
data$E2<-data$E2*100/data$NST # Pour chaque point de mesure, Calcul le pourcentage de vers à l'état E2
data$E1<-data$E1*100/data$NST # Pour chaque point de mesure, Calcul le pourcentage de vers à l'état E1
data$E0<-data$E0*100/data$NST # Pour chaque point de mesure, Calcul le pourcentage de vers à l'état E0
data$NST<-data$NST*100/20 # Pour chaque point de mesure, calcul le pourcentage de vers contenue dans le puit


for(i in 1:length(data$DDP)) # Pour chaque ligne de la dataframe data
{
  if(data$j[i]>Jmax) # donne la valeur NA à « j » si « j » et supérieur à Jmax (ici 29 jours)
  {
    data$j[i]<-NA
  }
}
data<-na.rm(data) # Supprime les lignes dans lesquelles j est égal à NA. (Ici, supprime les lignes de données après Jmax (29 jours))


data$line<-NA # ajoute une colonne ligne dans la dataframe data
for(i in 1:length(data$DDP))
{
  data$line[i]<-i # Donne un numéro de ligne a la colonne ligne dans data
}

data$DDP<-as.numeric(data$DDP) # informe que les valeurs contenues dans la colonne DDP sont numériques

Day_real<-rev(N_of_variable(data$j)) # utilise la fonction N_of_variable pour récupérer les différents jours de mesures
Media_real<-N_of_variable(data$Milieu) # utilise la fonction N_of_variable pour récupérer les différents numéros des milieux utilisés
Well_real<-N_of_variable(data$Puit) # utilise la fonction N_of_variable pour récupérer les différents numéros des puits 
Perfusion_day_real<-N_of_variable(data$DDP) # utilise la fonction N_of_variable pour récupérer les différentes dates de perfusion

Media<-N_of_variable(data$Milieu) # utilise la fonction N_of_variable pour récupérer les différents numéros des milieux 
for (Medias_number in 1:max(Media)) # Pour chaque milieu
{
  data_Media<-data[data$Milieu %in% c(Media[Medias_number]),] #Eclate les données par milieu de culture
  Well<-N_of_variable(data_Media$Puit) # Pour les données du milieu sélectionné, utilise la fonction N_of_variable pour récupérer les différents numéros de puits
  for(Wells_number in 1:max(Well)) #Pour chaque puit
  {
    data_Media_Well<-data_Media[data_Media$Puit %in% c(Well[Wells_number]),] #Eclate les données par numéro de puit
    Perfusion_day<-N_of_variable(data_Media_Well$DDP) #Pour chaque puit, utilise la fonction N_of_variable pour récupérer les différentes dates de perfusion
    for(Perfusions_number in 1:length(Perfusion_day)) # Pour chaque date de perfusion
    {
      data_Media_Well_Perfusions<-data_Media_Well[data_Media_Well$DDP %in% c(Perfusion_day[Perfusions_number]),] #Eclate les données par date de perfusion
      
      for(i in 1:length(data_Media_Well_Perfusions$DDP)) # Pour chaque ligne i du jeux de données éclaté data_Media_Well_Perfusions (exemple ligne 1 de data M199, puit N°1, date de perfusion N°1)
      {
        for(k in 1:length(data_Media_Well_Perfusions$DDP)) # Pour chaque ligne k du jeux de données éclaté data_Media_Well_Perfusions (exemple ligne 2 de data M199, puit N°1, date de perfusion N°1)
        {
          if(i!=k && data_Media_Well_Perfusions$j[k]==data_Media_Well_Perfusions$j[i]) #si la ligne i et la ligne K sont différentes et que les données contenues dans les lignes i et k sont identiques
          {
            print(paste("problem Two same days at line",data_Media_Well_Perfusions$line[k] ,"and line", data_Media_Well_Perfusions$line[i], "day",data_Media_Well_Perfusions$j[k])) # affiche qu'il y a un problème de double mesure (à corriger manuellement dans le jeu de donnée txt) 
          }
        }
      }
    }
  }
}

x<-NA
for (Medias_number in 1:max(Media)) # même chose que précédemment : Pour chaque milieu 
{
  data_Media<-data[data$Milieu %in% c(Media[Medias_number]),] # Eclate les données par milieu
  Well<-N_of_variable(data_Media$Puit)
  for(Wells_number in 1:max(Well)) # Pour chaque puit
  {
    data_Media_Well<-data_Media[data_Media$Puit %in% c(Well[Wells_number]),] # Eclate les données par numéro de puits
    Perfusion_day<-N_of_variable(data_Media_Well$DDP)
    for(Perfusions_number in 1:length(Perfusion_day)) # pour chaque date de perfusion
    {
      data_Media_Well_Perfusions<-data_Media_Well[data_Media_Well$DDP %in% c(Perfusion_day[Perfusions_number]),] # Eclate les données par date de perfusion (exemple data M199, puit N°1, date de perfusion N°1)
      for(i in (length(data_Media_Well_Perfusions$DDP)):1) # de la dernière ligne à la première ligne du jeu de données éclaté
      {
        if(data_Media_Well_Perfusions$E0[i]<data_Media_Well_Perfusions$E0[i-1]&&i!=1) # si des faux parasites morts ont été enregistrés dans le jeu de donnée
        {
          #si oui, affiche deux lignes concernées et leur contenu
          print("before correction") #if yes, print were it appears
          print(paste("Media",data_Media_Well_Perfusions$Milieu[i-1],
                      "DDP",data_Media_Well_Perfusions$DDP[i-1],
                      "/Day ",data_Media_Well_Perfusions$j[i-1],
                      "/E2 ",data_Media_Well_Perfusions$E2[i-1],
                      "/E1 ",data_Media_Well_Perfusions$E1[i-1],
                      "/E0 ",data_Media_Well_Perfusions$E0[i-1]))
          print(paste("Media",data_Media_Well_Perfusions$Milieu[i],
                      "DDP",data_Media_Well_Perfusions$DDP[i],
                      "/Day ",data_Media_Well_Perfusions$j[i],
                      "/E2 ",data_Media_Well_Perfusions$E2[i],
                      "/E1 ",data_Media_Well_Perfusions$E1[i],
                      "/E0 ",data_Media_Well_Perfusions$E0[i]))
          
          #Applique les corrections suivantes
          data_Media_Well_Perfusions$E1[i-1]<-data_Media_Well_Perfusions$E1[i-1]+(data_Media_Well_Perfusions$E0[i-1]-data_Media_Well_Perfusions$E0[i])
          data_Media_Well_Perfusions$E0[i-1]<-data_Media_Well_Perfusions$E0[i] # applique la correction suivante
          
          #affiche les deux lignes concernées et leur contenu après correction
          print("after correction")
          print(paste("Media",data_Media_Well_Perfusions$Milieu[i-1], # affiche la correction
                      "DDP",data_Media_Well_Perfusions$DDP[i-1],
                      "/Day ",data_Media_Well_Perfusions$j[i-1],
                      "/E2 ",data_Media_Well_Perfusions$E2[i-1],
                      "/E1 ",data_Media_Well_Perfusions$E1[i-1],
                      "/E0 ",data_Media_Well_Perfusions$E0[i-1]))
          print(paste("Media",data_Media_Well_Perfusions$Milieu[i],
                      "DDP",data_Media_Well_Perfusions$DDP[i],
                      "/Day ",data_Media_Well_Perfusions$j[i],
                      "/E2 ",data_Media_Well_Perfusions$E2[i],
                      "/E1",data_Media_Well_Perfusions$E1[i],
                      "/E0 ",data_Media_Well_Perfusions$E0[i]))
        }
      }
      x<-rbind(x,data_Media_Well_Perfusions) # rassemble les données éclatées dans la dataframe x après chaque traitement de la boucle
    }
  }
}
data<-na.rm(x) # supprimes les lignes NA ajoutées dans la dataframe x suite à la fonction Rbind  puis enregistre les données dans la dataframe data


#---------------Stat

#test J29 sans les "sans sérums" et sans les "humans"
Stat<-subset(data, select = -c(2,3,6,7,8,9,10,12))
Stat$Milieu <- factor(Stat$Milieu)

Stat$Day<-Stat$j
Stat$Media<-NA
Stat$Sera<-NA

for(i in 1:length(Stat$DDP))
{
  if(Stat$Milieu[i]==1)
  {
    Stat$Media[i]<-"M199"
    Stat$Sera[i]<-"Calf"
  }
  if(Stat$Milieu[i]==2)
  {
    Stat$Media[i]<-"DMEM"
    Stat$Sera[i]<-"Calf"
  }
  if(Stat$Milieu[i]==3)
  {
    Stat$Media[i]<-"RPMI"
    Stat$Sera[i]<-"Calf"
  }
  if(Stat$Milieu[i]==4)
  {
    Stat$Media[i]<-"M199"
    Stat$Sera[i]<-"Horse"
  }
  if(Stat$Milieu[i]==5)
  {
    Stat$Media[i]<-"DMEM"
    Stat$Sera[i]<-"Horse"
  }
  if(Stat$Milieu[i]==6)
  {
    Stat$Media[i]<-"RPMI"
    Stat$Sera[i]<-"Horse"
  }
  if(Stat$Milieu[i]==7)
  {
    Stat$Media[i]<-"M199"
    Stat$Sera[i]<-"Human"
  }
  if(Stat$Milieu[i]==8)
  {
    Stat$Media[i]<-"DMEM"
    Stat$Sera[i]<-"Human"
  }
  if(Stat$Milieu[i]==9)
  {
    Stat$Media[i]<-"RPMI"
    Stat$Sera[i]<-"Human"
  }
  if(Stat$Milieu[i]==10)
  {
    Stat$Media[i]<-"M199"
    Stat$Sera[i]<-"No serum"
  }
  if(Stat$Milieu[i]==11)
  {
    Stat$Media[i]<-"DMEM"
    Stat$Sera[i]<-"No serum"
  }
  if(Stat$Milieu[i]==10)
  {
    Stat$Media[i]<-"RPMI"
    Stat$Sera[i]<-"No serum"
  }
}



#selection du milieu
Table_numero_Milieux<-c(1,2,3,4,5,6,7,8,9,10,11,12)
Table_Nom_Milieux<-c("M199 Calf","DMEM Calf","RPMI Calf","M199 Horse","DMEM Horse","RPMI Horse","M199 Human","DMEM Human","RPMI Human","M199 No Serum","DMEM No Serum","RPMI No Serum")

Table_selection_numero<-NA
Table_selection_Nom<-NA

for (Numero_Selection in 1:12) #1 à 6 pour garder le sérum de veau et le sérum de cheval. 1 à 12 pour tous les prendre
{
  Table_selection_numero[Numero_Selection]<-Table_numero_Milieux[Numero_Selection]
  Table_selection_Nom[Numero_Selection]<-Table_Nom_Milieux[Numero_Selection]
}
Stat<-Stat[Stat$Milieu %in% Table_selection_numero,]
Stat$Milieu <- factor(Stat$Milieu, labels = Table_selection_Nom)


Stat$Mort<-Stat$E0
Stat$Jours<-Stat$j
Stat$Alive<-100-Stat$E0

#selection des milieux de culture

Stat$DDP <- factor(Stat$DDP)
Stat$Sera<- factor(Stat$Sera)
Time<-N_of_variable(Stat$Jours)

J29<-subset(Stat,Day==29)
#J29<-subset(J29,Sera!="No serum")
#J29<-subset(J29,Sera!="Human")
table(Stat$Milieu,Stat$DDP,useNA="ifany")
stripchart(Alive~Milieu,data=J29,method="jitter",cex.axis=0.8,vertical=TRUE)
AnovaJ29<-aov(Alive~Milieu+DDP,data=J29)
qqnorm(residuals(AnovaJ29))
qqline(residuals(AnovaJ29))
print(fitted(AnovaJ29))
plot(rstudent(AnovaJ29)~fitted(AnovaJ29))
abline(h=0,lty="dotted")
print(summary(AnovaJ29))

# Avec la fonction glht
comparaisons <- glht(AnovaJ29,linfct=mcp(Milieu="Tukey"))
# Afficher les comparaisons
summary(comparaisons)
cld(comparaisons)
tapply(J29$Alive,J29$Milieu,mean)
R2_for_J29 <- summary(AnovaJ29)[[1]][["Sum Sq"]]
R2_for_J29
sum(R2_for_J29[1:(length(R2_for_J29)-1)])/sum(R2_for_J29)

# Stat comparatif
Savestat<-Stat


Stat<-subset(Stat,Media=="M199")
for (t in length(Time):1)
{
  print(paste("Data_of_J",Time[t],sep=""))
  print(paste("Media : M199"))
  print(paste(""))
  assign(paste("j",Time[t],sep=""),Stat[Stat$Jours %in% c(Time[t]),])
  #print(table(get(paste("j",Time[t],sep=""))$Milieu,get(paste("j",Time[t],sep=""))$DDP,useNA="ifany"))
  stripchart(Mort~Milieu,data=get(paste("j",Time[t],sep="")),method="jitter",cex.axis=0.8,vertical=TRUE)
  
  
  assign(paste("AnovaJ",Time[t],sep=""),aov(Mort~Milieu+DDP,data=get(paste("j",Time[t],sep=""))))
  Anova<-get(paste("AnovaJ",Time[t],sep=""))
  #qqnorm(residuals(get(paste("AnovaJ",Time[t],sep=""))))
  #qqline(residuals(get(paste("AnovaJ",Time[t],sep=""))))
  qqnorm(residuals(Anova))
  qqline(residuals(Anova))
  #print(fitted(get(paste("AnovaJ",Time[t],sep=""))))
  
  #plot(rstudent(get(paste("AnovaJ",Time[t],sep="")))~fitted(get(paste("AnovaJ",Time[t],sep=""))))
  plot(rstudent(Anova)~fitted(Anova))
  abline(h=0,lty="dotted")
  
  print(summary(get(paste("AnovaJ",Time[t],sep=""))))
  
  # Avec la fonction glht
  assign(paste("comparaisons",Time[t],sep=""),glht(get(paste("AnovaJ",Time[t],sep="")),linfct=mcp(Milieu="Tukey")))
  # Afficher les comparaisons
  #print(summary(get(paste("comparaisons",Time[t],sep=""))))
  #print(cld(get(paste("comparaisons",Time[t],sep=""))))
  
  #print(tapply(get(paste("j",Time[t],sep=""))$Mort,get(paste("j",Time[t],sep=""))$Milieu,mean))
  
  
  assign(paste("R2_for_J",Time[t],sep=""),summary(get(paste("AnovaJ",Time[t],sep="")))[[1]][["Sum Sq"]])
  #print(get(paste("R2_for_J",Time[t],sep="")))
  #print(sum(get(paste("R2_for_J",Time[t],sep=""))[1:(length(get(paste("R2_for_J",Time[t],sep="")))-1)])/sum(get(paste("R2_for_J",Time[t],sep=""))))
}
Stat<-Savestat

Stat<-subset(Stat,Media=="DMEM")
for (t in length(Time):1)
{
  print(paste("Data_of_J",Time[t],sep=""))
  print(paste("Media : DMEM"))
  print(paste(""))
  assign(paste("j",Time[t],sep=""),Stat[Stat$Jours %in% c(Time[t]),])
  #print(table(get(paste("j",Time[t],sep=""))$Milieu,get(paste("j",Time[t],sep=""))$DDP,useNA="ifany"))
  stripchart(Mort~Milieu,data=get(paste("j",Time[t],sep="")),method="jitter",cex.axis=0.8,vertical=TRUE)
  
  
  assign(paste("AnovaJ",Time[t],sep=""),aov(Mort~Milieu+DDP,data=get(paste("j",Time[t],sep=""))))
  Anova<-get(paste("AnovaJ",Time[t],sep=""))
  #qqnorm(residuals(get(paste("AnovaJ",Time[t],sep=""))))
  #qqline(residuals(get(paste("AnovaJ",Time[t],sep=""))))
  qqnorm(residuals(Anova))
  qqline(residuals(Anova))
  #print(fitted(get(paste("AnovaJ",Time[t],sep=""))))
  
  #plot(rstudent(get(paste("AnovaJ",Time[t],sep="")))~fitted(get(paste("AnovaJ",Time[t],sep=""))))
  plot(rstudent(Anova)~fitted(Anova))
  abline(h=0,lty="dotted")
  
  print(summary(get(paste("AnovaJ",Time[t],sep=""))))
  
  # Avec la fonction glht
  assign(paste("comparaisons",Time[t],sep=""),glht(get(paste("AnovaJ",Time[t],sep="")),linfct=mcp(Milieu="Tukey")))
  # Afficher les comparaisons
  #print(summary(get(paste("comparaisons",Time[t],sep=""))))
  #print(cld(get(paste("comparaisons",Time[t],sep=""))))
  
  #print(tapply(get(paste("j",Time[t],sep=""))$Mort,get(paste("j",Time[t],sep=""))$Milieu,mean))
  
  
  assign(paste("R2_for_J",Time[t],sep=""),summary(get(paste("AnovaJ",Time[t],sep="")))[[1]][["Sum Sq"]])
  #print(get(paste("R2_for_J",Time[t],sep="")))
  #print(sum(get(paste("R2_for_J",Time[t],sep=""))[1:(length(get(paste("R2_for_J",Time[t],sep="")))-1)])/sum(get(paste("R2_for_J",Time[t],sep=""))))
}
Stat<-Savestat

Stat<-subset(Stat,Media=="RPMI")
for (t in length(Time):1)
{
  print(paste("Data_of_J",Time[t],sep=""))
  print(paste("Media : RPMI"))
  print(paste(""))
  assign(paste("j",Time[t],sep=""),Stat[Stat$Jours %in% c(Time[t]),])
  #print(table(get(paste("j",Time[t],sep=""))$Milieu,get(paste("j",Time[t],sep=""))$DDP,useNA="ifany"))
  stripchart(Mort~Milieu,data=get(paste("j",Time[t],sep="")),method="jitter",cex.axis=0.8,vertical=TRUE)
  
  
  assign(paste("AnovaJ",Time[t],sep=""),aov(Mort~Milieu+DDP,data=get(paste("j",Time[t],sep=""))))
  Anova<-get(paste("AnovaJ",Time[t],sep=""))
  #qqnorm(residuals(get(paste("AnovaJ",Time[t],sep=""))))
  #qqline(residuals(get(paste("AnovaJ",Time[t],sep=""))))
  qqnorm(residuals(Anova))
  qqline(residuals(Anova))
  #print(fitted(get(paste("AnovaJ",Time[t],sep=""))))
  
  #plot(rstudent(get(paste("AnovaJ",Time[t],sep="")))~fitted(get(paste("AnovaJ",Time[t],sep=""))))
  plot(rstudent(Anova)~fitted(Anova))
  abline(h=0,lty="dotted")
  
  print(summary(get(paste("AnovaJ",Time[t],sep=""))))
  
  # Avec la fonction glht
  assign(paste("comparaisons",Time[t],sep=""),glht(get(paste("AnovaJ",Time[t],sep="")),linfct=mcp(Milieu="Tukey")))
  # Afficher les comparaisons
  #print(summary(get(paste("comparaisons",Time[t],sep=""))))
  #print(cld(get(paste("comparaisons",Time[t],sep=""))))
  
  #print(tapply(get(paste("j",Time[t],sep=""))$Mort,get(paste("j",Time[t],sep=""))$Milieu,mean))
  
  
  assign(paste("R2_for_J",Time[t],sep=""),summary(get(paste("AnovaJ",Time[t],sep="")))[[1]][["Sum Sq"]])
  #print(get(paste("R2_for_J",Time[t],sep="")))
  #print(sum(get(paste("R2_for_J",Time[t],sep=""))[1:(length(get(paste("R2_for_J",Time[t],sep="")))-1)])/sum(get(paste("R2_for_J",Time[t],sep=""))))
}
Stat<-Savestat


#--------------Fin stat

data_save<-data # enregistre les données de data dans data_save

Tableau<-NA
data_save<-data # sauvegarde de la dataframe data
data$sum<-NA
data$sum<-data$E2+data$E1+data$E0 # vérifie que la somme de E2+E1+E0 = 100%

#création de nouvelles colonnes dans la dataframe data pour les prochains calculs (moyennes et des écarts types)

data$Day<-NA
data$Medium<-NA

data$Viable<-NA
data$Degraded<-NA
data$Mortality<-NA
data$Viability<-NA
data$Couple<-NA

data$deltaE2<-NA
data$deltaE1<-NA
data$deltaE0<-NA
data$deltaCouple<-NA

Tableau<-data
Tableau<-subset(data, select = -c(1,2,3,4,5,6,7))

print("0%")
ligne<-1
Media<-N_of_variable(data$Milieu) #utilise la fonction N_of_variable pour récupérer les différents numéros des milieux utilisés
for (Medias_number in 1:max(Media))
{
  data_Media<-data[data$Milieu %in% c(Media[Medias_number]),] #Eclate les données par milieu
  Day_to_verify<-rev(N_of_variable(data_Media$j)) #utilise la fonction N_of_variable pour récupérer les différents jours de mesure (exemple, jours 29)
  for(Day_number in 1:length(Day_to_verify))
  {
    data_Media_Day<-data_Media[data_Media$j %in% c(Day_real[Day_number]),] #Eclate les données par jours de mesure (exemple : M199 Horse au jours 29)

    Tableau$Viable[ligne]<-mean(data_Media_Day$E2) # calcul et enregistre la moyenne des vers viables dans la dataframe tableau
    Tableau$Degraded[ligne]<-mean(data_Media_Day$E1) # calcul et enregistre la moyenne des vers dégradés dans la dataframe tableau
    Tableau$Mortality[ligne]<-mean(data_Media_Day$E0) # calcul et enregistre la moyenne des vers morts dans la dataframe tableau
    Tableau$Viability[ligne]<-100-mean(data_Media_Day$E0) # calcul et enregistre la moyenne de la viabilité des vers dans la dataframe tableau

    Tableau$deltaE2[ligne]<-sd(data_Media_Day$E2) # calcul et enregistre  l'écart type des vers viables dans la dataframe tableau
    Tableau$deltaE1[ligne]<-sd(data_Media_Day$E1) # calcul et enregistre  l'écart type des vers dégradés dans la dataframe tableau
    Tableau$deltaE0[ligne]<-sd(data_Media_Day$E0) # calcul et enregistre  l'écart type des vers morts dans la dataframe tableau

    Tableau$Day[ligne]<-data_Media_Day$j # enregistre le jour de mesure dans la data frame tableau
    Tableau$Medium[ligne]<-data_Media_Day$Milieu # enregistre le type de milieu dans la data frame tableau

    Tableau$Couple[ligne]<-mean(data_Media_Day$NDC) #calcul et enregistre la moyenne des vers en couple dans la dataframe tableau
    Tableau$deltaCouple[ligne]<-sd(data_Media_Day$NDC) #calcul et enregistre l'a moyenne l'écart type  des vers en couple dans la dataframe tableau

    ligne<-ligne+1 # change la ligne à vérifier dans la dataframe tableau pour le prochain calcul de la boucle
  }
}
Tableau<-na.rm(Tableau) # supprime toutes les lignes de la dataframe tableau contenant un NA (ici, conserve uniquement les données après le calcul des moyennes et écarts type)

data<-Tableau # enregistre les données contenue dans tableau dans la dataframe data

data$Medium<-as.factor(data$Medium) # précise que la colonne des milieux est de type factor (de type caractère)
#précise le nom du milieu pour chaque numéro de milieu
levels(data$Medium)<-c("M199 Calf","DMEM Calf","RPMI Calf","M199 Horse","DMEM Horse","RPMI Horse","M199 Human","DMEM Human","RPMI Human","M199 No Serum","DMEM No Serum","RPMI No Serum")

data$Day<-as.numeric(data$Day) # précise que la colonne Day et de type numérique

#Tous les milieux sur un seul graphique ---------------------------

#Impose un ordre d'affichage de chaque milieu dans le futur graphique
data$Medium <- factor(data$Medium, levels = c("RPMI Human","DMEM Human","M199 Human","RPMI Horse","DMEM Horse","M199 Horse","RPMI Calf","DMEM Calf","M199 Calf","RPMI No Serum","DMEM No Serum","M199 No Serum"))

data_save<-data

#trace le graphique de viabilité des vers en fonction du temps et du type de milieu
a<-ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 102),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.1, 0.5),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+

  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))
a

ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 102),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),
        #legend.position = c(0.25, 0.35),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
  
  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))

ggsave(filename = file.path("figs",paste("Figure 4 supplementary","dose.TIFF")), width = 11.8, height = 8.3, units = "cm",dpi = 600)

ggsave(filename = file.path("figs",paste("Figure 4 supplementary","dose.pdf")), width = 11.8, height = 8.3, units = "cm",dpi = 600)


# Tous les milieux sur un seul graphique terminé ---------------------------

data<-data_save

x<-0.71 #c(0.71, 0.1981)
y<-0.28

# Tous les milieux RPMI sur un seul graphique ---------------------------

#trace le graphique de viabilité des vers en fonction du temps et des milieux RPMI
data<-data[data$Medium %in% c("RPMI Human","RPMI Horse","RPMI Calf","RPMI No Serum"),]
data$Medium <- factor(data$Medium, levels = c("RPMI Human","RPMI Horse","RPMI Calf","RPMI No Serum"))



b<-ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 127),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),

        legend.position = c(x, y-0.00005),
        #legend.position = c(0.78, 0.35),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
        
  geom_segment(aes(x=15, xend=29, y=105, yend=105), colour="black", size = 0.5, show.legend=F)+
        geom_segment(aes(x=15, xend=15, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
        geom_segment(aes(x=29, xend=29, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
        annotate("text", x = 22, y = 109+0.2, label = "p-value < 0.05 : Horse vs Calf", colour="black", size = 2)+
        
  
  geom_segment(aes(x=17, xend=29, y=105+8, yend=105+8), colour="black", size = 0.5, show.legend=F)+
        geom_segment(aes(x=17, xend=17, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
        geom_segment(aes(x=29, xend=29, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
        annotate("text", x = 21.8, y = 109+8.5, label = "p-value < 0.05 : Human vs Calf", colour="black", size = 2)+
        
  geom_segment(aes(x=26, xend=29, y=105+8+8.5, yend=105+8+8.5), colour="black", size = 0.5, show.legend=F)+
        geom_segment(aes(x=26, xend=26, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
        geom_segment(aes(x=29, xend=29, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
        annotate("text", x = 21.6, y = 109+8+8.7, label = "p-value < 0.05 : Human vs Horse", colour="black", size = 2)+

  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))


ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 127),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),
        
        legend.position = c(x, y-0.00005),
        #legend.position = c(0.78, 0.35),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
  
  geom_segment(aes(x=15, xend=29, y=105, yend=105), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=15, xend=15, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 22, y = 109+0.2, label = "p-value < 0.05 : Horse vs Calf", colour="black", size = 2)+
  
  
  geom_segment(aes(x=17, xend=29, y=105+8, yend=105+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=17, xend=17, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.8, y = 109+8.5, label = "p-value < 0.05 : Human vs Calf", colour="black", size = 2)+
  
  geom_segment(aes(x=26, xend=29, y=105+8+8.5, yend=105+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=26, xend=26, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.6, y = 109+8+8.7, label = "p-value < 0.05 : Human vs Horse", colour="black", size = 2)+
  
  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))

ggsave(filename = file.path("figs",paste("Figure 4 RPMI","sup.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)


ggsave(filename = file.path("figs",paste("Figure 4 RPMI","sup.pdf")), width = 8.3, height = 8.3, units = "cm",dpi = 600)



# Tous les milieux RPMI sur un seul graphique terminé ---------------------------


data<-data_save

# Tous les milieux DMEM sur un seul graphique ---------------------------

#trace le graphique de viabilité des vers en fonction du temps et des milieux DMEM

data<-data[data$Medium %in% c("DMEM Human","DMEM Horse","DMEM Calf","DMEM No Serum"),]
data$Medium <- factor(data$Medium, levels = c("DMEM Human","DMEM Horse","DMEM Calf","DMEM No Serum"))



c<-ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 127),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),

        legend.position = c(x+0.01,y),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
  geom_segment(aes(x=8, xend=29, y=105, yend=105), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=8, xend=8, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 22, y = 109+0.2, label = "p-value < 0.05 : Horse vs Calf", colour="black", size = 2)+
  
  
  geom_segment(aes(x=5, xend=29, y=105+8, yend=105+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=5, xend=5, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.8, y = 109+8.5, label = "p-value < 0.05 : Human vs Calf", colour="black", size = 2)+
  
  geom_segment(aes(x=15, xend=29, y=105+8+8.5, yend=105+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=15, xend=15, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.6, y = 109+8+8.7, label = "p-value < 0.05 : Human vs Horse", colour="black", size = 2)+

  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))



ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 127),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),
        
        legend.position = c(x+0.01,y),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
  geom_segment(aes(x=8, xend=29, y=105, yend=105), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=8, xend=8, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 22, y = 109+0.2, label = "p-value < 0.05 : Horse vs Calf", colour="black", size = 2)+
  
  
  geom_segment(aes(x=5, xend=29, y=105+8, yend=105+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=5, xend=5, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.8, y = 109+8.5, label = "p-value < 0.05 : Human vs Calf", colour="black", size = 2)+
  
  geom_segment(aes(x=15, xend=29, y=105+8+8.5, yend=105+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=15, xend=15, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.6, y = 109+8+8.7, label = "p-value < 0.05 : Human vs Horse", colour="black", size = 2)+
  
  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))

ggsave(filename = file.path("figs",paste("Figure 4 DMEM","sup.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

ggsave(filename = file.path("figs",paste("Figure 4 DMEM","sup.pdf")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Tous les milieux DMEM sur un seul graphique terminé ---------------------------

data<-data_save

# Tous les milieux M199 sur un seul graphique ---------------------------

#trace le graphique de viabilité des vers en fonction du temps et des milieux M199
data<-data[data$Medium %in% c("M199 Human","M199 Horse","M199 Calf","M199 No Serum"),]
data$Medium <- factor(data$Medium, levels = c("M199 Human","M199 Horse","M199 Calf","M199 No Serum"))

d<-ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 127),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),

        legend.position = c(x, y),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
  geom_segment(aes(x=12, xend=29, y=105, yend=105), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=12, xend=12, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 22, y = 109+0.2, label = "p-value < 0.05 : Horse vs Calf", colour="black", size = 2)+
  
  
  geom_segment(aes(x=10, xend=29, y=105+8, yend=105+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=10, xend=10, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.8, y = 109+8.5, label = "p-value < 0.05 : Human vs Calf", colour="black", size = 2)+
  
  geom_segment(aes(x=15, xend=29, y=105+8+8.5, yend=105+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=15, xend=15, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.6, y = 109+8+8.7, label = "p-value < 0.05 : Human vs Horse", colour="black", size = 2)+

  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))


ggplot(data,aes(x=Day, y=Viability, group = Medium,color = Medium, na.rm=TRUE))+
  geom_point(aes(shape=Medium),size=1.5)+
  geom_step()+
  #scale_color_grey(start=0.7, end=0)+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15,16,17,18))+
  geom_errorbar(aes(ymin=Viability-deltaE0, ymax=Viability+deltaE0),width=1)+
  scale_y_continuous( limits=c(0, 127),breaks=seq(0, 102, by = 20))+
  scale_x_continuous( limits=c(0, Jmax+1),breaks=seq(0, Jmax+1, by = 5))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.background = element_rect(fill = "transparent"),
        
        legend.position = c(x, y),
        #legend.position ="bottom",
        legend.direction="vertical",
        #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.key.size = unit(0.5, "lines"))+
  geom_segment(aes(x=12, xend=29, y=105, yend=105), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=12, xend=12, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103, yend=107), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 22, y = 109+0.2, label = "p-value < 0.05 : Horse vs Calf", colour="black", size = 2)+
  
  
  geom_segment(aes(x=10, xend=29, y=105+8, yend=105+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=10, xend=10, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8, yend=107+8), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.8, y = 109+8.5, label = "p-value < 0.05 : Human vs Calf", colour="black", size = 2)+
  
  geom_segment(aes(x=15, xend=29, y=105+8+8.5, yend=105+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=15, xend=15, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=29, xend=29, y=103+8+8.5, yend=107+8+8.5), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 21.6, y = 109+8+8.7, label = "p-value < 0.05 : Human vs Horse", colour="black", size = 2)+
  
  #ggtitle(paste("Worms survival depending on culture medium"))+
  labs(y= "Living worms (%)", x = paste("Time (days)"),colour = paste("Medium"))

ggsave(filename = file.path("figs",paste("Figure 4 M199","sup.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

ggsave(filename = file.path("figs",paste("Figure 4 M199","sup.pdf")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# # Tous les milieux DMEM sur un seul graphique terminé ---------------------------



figure <- ggarrange(b, d, c, a, labels = c("A", "B", "C","D"),ncol = 2, nrow = 2)

figure
ggsave(filename = file.path("figs",paste("Figure 4 All","sup.TIFF")), width = 8.3*2, height = 8.3*2, units = "cm",dpi = 600)
ggsave(filename = file.path("figs",paste("Figure 4 All","sup.pdf")), width = 8.3*2, height = 8.3*2, units = "cm",dpi = 600)

