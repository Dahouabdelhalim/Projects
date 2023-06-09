# Initialise les librairies-----------------------------------
rm(list=ls()) # Supprime les données présentes dans la ram et la rom de RStudio

packages<-c("here","sp") # Chemin du répertoire avec les données 

for(i in packages){if(!require(i,character.only = T)){install.packages(i,dep=T)}}
lapply(packages, library, character.only = T) # Chargement de la liste des packages
library(ggplot2)
library(scales)
library(dplyr)
library(na.tools)
library(concatenate)
library(writexl)
library(tidyr)

# fonction  -----------------------------------

opentxt<-function(x) # Fonction utilisée pour ouvrir les fichiers textes
  
{
  data<-read.table(here("data",x),sep='\\t', header=T) # Ouvre et lit les fichiers txt pour les enregistrer dans une dataframe
}
txt<-dir("data",pattern = ".txt") # Enregistre le nombre de fichier texte et leur nom dans un tableau 
N<-length(txt) # Enregistre dans la variable N le nombre total de fichier texte
data_real<-structure(list(character()), class = "data.frame") # Création de la dataframe data

for (i in 1:N) # Ouvre tous les fichiers textes et l'enregistre dans la dataframe data_real
{
  x<-opentxt(txt[i]) # Ouvre les fichiers txt et les enregistres dans x
  x$Number<-i # Ajoute une colonne Number dans la dataframe x
  data_real<-rbind(data_real,x) # Enregistre dans la dataframe data_real les données
  
}
data_save<-na.rm(data_real) # Supprimes tous les lignes de data_real contenants des NA 


# fonctions -----------------------------------
N_of_variable<-function(Variable) # Fonction utilisée pour enregistrer dans un tableau l’ensemble des différentes valeurs existantes dans la colonne choisi de la dataframe en entrée de fonction (ici, Variable) (exemple : nom des milieux de culture, les jours pour chaque contrôles, …)
{
  compteur<-1 # Initialise le compteur à la valeur 1
  table<-NA # Création du tableau 
  table[compteur]<--1 # Donne la valeur -1 à la première colonne du tableau
  Num<-NA # Création de la variable Num avec comme valeur NA
  exist<-0 # Création de la variable exist a 0
  for (i in 1:length(Variable)) # Pour chaque ligne contenue dans Variable 
  {
    Num<-Variable[i] # Enregistre dans Num la valeur contenue dans la ligne i
    for (j in 1:length(table)) # De la première colonne à la dernière colonne du tableau table
    {
      if(Num==table[j]) # Vérifie si la valeur est déjà présente dans la table
      {
        j<-length(table) # Si oui, arrêt la boucle for
        exist<-1 # Indique que la valeur est déjà contenue dans la table
      }
    }
    if(exist==0) # Si la valeur n'existe pas dans la table
    {
      table[compteur]<-Num # Enregistre la valeur dans la colonne « compteur » de la table 
      compteur<-compteur+1 # Change de colonne « compteur » de la table 
      print(Num) # Affiche la valeur contenue dans Num 
    }
    else(exist<-0) # Sinon, réinitialise la valeur exist a 0 pour la lecture de la ligne suivante
  }
  
  Ndose<-rev(sort(table)) # Range les valeurs contenues dans le tableau dans décroissant
  return(Ndose) # Renvoi les données de la table en sortie de la fonction
}


# Comme les vers ne sortent pas des puces pour les mêmes temps, (exemple : puce 1, 1 mâle de sortie à 5min 17 seconde, Puce 2, 1 mâle de sortie à 6 min et 25 secondes)
# J’ai créé les 3 fonction suivantes pour aligner les temps sur l’instant où la drogue a été injectée dans chaque puce

Vecteur_time<-function(data,Temps,TableDose) # Création d'une dataframe vide ou le temps est noté de 0 à 60 minutes avec pas de 1 secondes (0.01)
{
  t<-seq(0,Temps, by=0.01) # Création d'un vecteur temporel "t" de 0 à Temps (60 minutes) avec un pas de 0.01
  Time<-NA # Création d'une variable Time de valeur NA
  for(i in 1:(length(TableDose)*max(data$Number))) #  (Nombre de dose multiplie par nombre max expériences) afin qu'elle soit suficament grandes pour changer toutes les présent dans la dataframe data_real
  {
    Time<-c(Time,t) # Ajoute dans Time autant de vecteur temporel qu'il y a de doses et d'expériences (exemple : dose 400 nM => 4 expériences et dose 0 nM => 3 expériences, donc ajouté 16 vecteur temporelles)
  }
  Time<-na.rm(Time) # Supprime les NA apparues dans Time suite à la fonction C() (la valeur initiale de Time est NA à la ligne 1)
  Time<-as.numeric(Time) # Définie Time comme numérique
  df<-data.frame(Time) # Définie Time comme une dataframe et l'enregistre dans df 
  df$Flow<-NA # Création d'une colonne Flow
  df$Dose<-NA # Création d'une colonne Dose
  df$Male_out<-NA # Création d'une colonne Male_out
  df$Female_out<-NA # Création d'une colonne Female_out
  df$Couple_out<-NA # Création d'une colonne Couple_out
  df$Video<-NA # Création d'une colonne video
  df$Drug<-NA # Création d'une colonne Drug
  df$Delta<-NA # Création d'une colonne Delta
  df$Unit<-NA # Création d'une colonne unit
  df$Number<-NA # Création d'une colonne number
  df$Drug_Used<-NA # Création d'une colonne Drug_Used
  return(df) # renvoi la dataframe df vide en sortie de la fonction  
}
Vecteur_Dose_Number<-function(df,data,Ndose) # Après la fonction Vecteur Times, ajoute une valeur de dose et un numéro d'expérience pour chaque vecteur temporel (Exemple : Vecteur N°1 (0 à 60 minutes), Dose = 400 nM, Numéro = 1 ; Vecteur N°2, Dose = 400 nM, Numéro = 2 ; ...)
{
  print("traitement 1 start")
  Number<-1 # Création d'un variable Number a 1
  compteur_Number<-0 # Création d'une variable compteur a 0
  colonne_dose<-1 # Création d'une variable colonne dose a 1
  compteur_dose<-0 # Création d'une variable compteur dose a 0
  
  for(i in 1:length(df$Time)) # De la première ligne à la dernière ligne de la dataframe
  {
    compteur_Number<-compteur_Number+1 # Ajoute +1 à compteur_Number
    compteur_dose<-compteur_dose+1 # Ajoute +1 à compteur_dose
    df$Number[i]<-Number # Donne un numéro d'expérience dans la ligne i de la colonne Number
    df$Dose[i]<-Ndose[colonne_dose] # Donne une valeur de dose dans la ligne i de la colonne dose
    
    if(compteur_Number>=(length(df$Time)/max(data$Number))) # Quand la variable compteur_number est supérieur à la taille de la dataframe / par le nombre de d'expérience (exemple, 16 vecteur/ 4 expérience), change le numéro de d'expérience
    {
      Number<-Number+1 # Change le numéro d'expérience à noter dans la dataframe
      compteur_Number<-0 # Réinitialise la variable compteur Number à 0
    }
    if(compteur_dose>=(length(df$Time)/(max(data$Number)*length(Ndose)))) # Quand la variable compteur_Dose est supérieur à la taille de la dataframe / par le nombre de d'expérience et le Nombre de dose (exemple, 16 vecteur/ [4 expérience x 2 Dose]), change le numéro la valeur de la dose
    {
      compteur_dose<-0 # Change le numéro de dose à noter dans la dataframe 
      colonne_dose<-colonne_dose+1 # Réinitialise la valeur de colonne_dose 
      #print(paste("traitement 1",(i/length(df$Time)*100),"%"))
      if(colonne_dose>length(Ndose)) # Si la colonne_dose est supérieur au nombre de dose 
      {
        colonne_dose<-1 # Réinitialise colonne dose à 1 (redémarre à la première dose)
      }
    }
  }
  print("traitement 1 finished")
  return(df) # renvoi la dataframe df en sortie de la fonction 
}
Placement_Data_In_Vecteur<-function(df,data) # Après les fonctions Vecteur_time et Vecteur_Dose_Number, place les données contenues dans data dans la dataframe temporel 
{
  print("traitement 2 Start")
  df$Time<-round(df$Time,2) # Arrondie les valeurs de la colonne time à 2 chiffre après la virgules pour éviter les erreurs
  df$Time<-as.numeric(df$Time) # Définie la colonne type de df comme numérique
  df$Number<-as.numeric(df$Number) # Définie la number type de df comme numérique
  df$Dose<-as.numeric(df$Dose) # Définie la colonne Dose de df comme numérique
  data$Real_Time<-as.numeric(data$Real_Time) # Définie la colonne Real_time  de data comme numérique
  data$Number<-as.numeric(data$Number) # Définie la colonne Number de data comme numérique
  data$Dose<-as.numeric(data$Dose) # Définie la colonne Dose de data comme numérique
  
  for(i in 1:length(data$Real_Time)) # Pour toutes les lignes de data
  {
    for(j in 1:length(df$Time)) # Pour toutes les ligne de df 
    {
      if((df$Time[j]==data$Real_Time[i])&&(df$Number[j]==data$Number[i])) # Si ligne Time de df est identique à la ligne Real_time de data et si ligne numéro de puce est identique au numéro de puce présent dans data
      {
        
        if(df$Dose[j]==data$Dose[i]) # Si la ligne dose de df est identique à la ligne dose de data
        {
          # Copie la ligne contenue dans la dataframe data dans la ligne de la dataframe df 
          df$Flow[j]<-data$Flow[i]
          df$Dose[j]<-data$Dose[i]
          df$Male_out[j]<-data$Male_out[i]
          df$Female_out[j]<-data$Female_out[i]
          df$Couple_out[j]<-data$Couple_out[i]
          df$Video[j]<-data$Video[i]
          df$Drug[j]<-data$Drug[i]
          df$Delta[j]<-data$Delta[i]
          df$Unit[j]<-data$Unit[i]
          df$Number[j]<-data$Number[i]
          df$Drug_Used[j]<-data$Drug_Used[i]
          #print(paste("traitement 2",(i/length(data$Real_Time)*100),"%"))
        }
      }
    }
  }
  print("traitement 2 Done")
  return(df) # Renvoi la dataframe df en sortie de la fonction 
}

Normalisation_sortie_vers<-function(df1) # Après les fonctions Vecteur_time, Vecteur_Dose_Number et Placement_Data_In_Vecteur, calcul le nombre de vers présent dans les puces au cours du temps 
{
  print("traitement 3 Start %")
  Table1<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA) # Création d'un tableau pour enregistrer les valeurs suivantes au cours du traitement : Flow; Male_out; Female_out; Couple_out; Video; Drug; Delta; Unit; Drug_Used
  
  df1$Male<-NA # Création d'une colonne Male
  df1$Female<-NA # Création d'une colonne Female
  df1$Couple<-NA # Création d'une colonne couple
  
  Table2<-c(NA,NA,NA) # Création d'une table pour enregistrer les valeurs suivantes au cours du traitement : Male, Female, Couple
  
  for (i in 1:length(df1$Time)) # Pour chaque ligne de la dataframe df1
  {
    if((df1$Time[i]-floor(df1$Time[i]))<=0.59) # Prend en compte que les lignes où les secondes sont comprise entre 0 et 59 secondes (et ne prend pas en compte les lignes où les secondes dépasse 59 secondes (59 à 99)) 
    {
      if((FALSE==is.na(df1$Male_out[i]))||df1$Time[i]==0) # Si il s'agit de la ligne avec le nombre de vers initiaux dans la puce
      {
        Table1<-c(df1$Flow[i],df1$Male_out[i],df1$Female_out[i],df1$Couple_out[i], df1$Video[i], df1$Drug[i], df1$Delta[i],df1$Unit[i],df1$Drug_Used[i]) # Enregistre les données initiales du nombre de vers dans la puce dans la table 1
      }
      if(df1$Time[i]==0) # Si la ligne correspond à T0
      {
        df1$Male[i]<-df1$Male_out[i] # Enregistre le nombre initiale de mâle dans mâle 
        Table2[1]<-df1$Male_out[i] # Enregistre le nombre initiale de mâle dans Table2
        df1$Female[i]<-df1$Female_out[i] # Enregistre le nombre initiale de Femelle dans Female
        Table2[2]<-df1$Female_out[i] # Enregistre le nombre initiale de Femelle dans Table2
        df1$Couple[i]<-df1$Couple_out[i] # Enregistre le nombre initiale de Couple dans Couple
        Table2[3]<-df1$Couple_out[i] # Enregistre le nombre initiale de Couple dans Table2
      }
      
      
      if((FALSE==is.na(df1$Male_out[i]))&&df1$Time[i]!=0) # Si le temps est différent de zero et qu'il existe une valeur contenue dans la ligne Male_out
      {
        df1$Male[i]<-Table2[1]-df1$Male_out[i] # Le nombre de mâle est égale au nombre de mâle (Table2, colonne 1) - le nombre de mâle sorties de la puce
        Table2[1]<-Table2[1]-df1$Male_out[i] # Enregistre le nombre de mâle restant dans la table 2 (Table2, colonne 1)
        
        df1$Female[i]<-Table2[2]-df1$Female_out[i] # Le nombre de Femelle est égale au nombre de Femelle (Table2, colonne 2) - le nombre de Femelle sorties de la puce
        Table2[2]<-Table2[2]-df1$Female_out[i] # Enregistre le nombre de Femelle restante dans la table 2 (Table2, colonne 2)
        df1$Couple[i]<-Table2[3]-df1$Couple_out[i] # Le nombre de Couple est égale au nombre de Couple (Table2, colonne 3) - le nombre de couple sorties de la puce
        Table2[3]<-Table2[3]-df1$Couple_out[i] # Enregistre le nombre de couple restant dans la table 2 (Table2, colonne 3)
      }
      
      if (TRUE==is.na(df1$Male_out[i])) # Si aucune données n'est présente dans la ligne testée de la colonne Male_out
      {
        df1$Flow[i]<-Table1[1] # Enregistre le débit contenue dans la table1 dans la colonne Flow de la dataframe
        df1$Male_out[i]<-0 # Enregistre aucune sortie de mâle dans la colonne Male_out de la dataframe
        df1$Female_out[i]<-0 # Enregistre aucune sortie de femelle dans la colonne Female_out de la dataframe
        df1$Couple_out[i]<-0 # Enregistre aucune sortie de couple dans la colonne Couple_out de la dataframe
        df1$Video[i]<-Table1[5] # Enregistre le numéro de vidéo contenue dans la table1 dans la colonne video de la dataframe
        df1$Drug[i]<-Table1[6] # Enregistre le statut (présent ou absent) contenue dans la table1 dans la colonne drug de la dataframe
        df1$Delta[i]<-Table1[7] # Enregistre la valeur de Delta (temps entre les 5 minutes de lavage et l'injection de la drogue) contenue dans la table1 dans la colonne Drug de la dataframe
        df1$Unit[i]<-Table1[8] # Enregistre l'unité métrique de la drogue utilisée contenue dans la table1 dans la colonne Unit de la dataframe
        df1$Drug_Used[i]<-Table1[9] # Enregistre le nom de la drogue contenue dans la table1 dans la colonne Drug_Used de la dataframe
        
        
        df1$Male[i]<-Table2[1] # Enregistre le nombre de mâle restant dans la puce  dans la colonne mâle de la dataframe
        df1$Female[i]<-Table2[2] # Enregistre le nombre de Femelle restante dans la puce dans la colonne Female de la dataframe
        df1$Couple[i]<-Table2[3] # Enregistre le nombre de Couple restant dans la puce dans la colonne Couple de la dataframe
      }
      else
      {
        #print(paste("traitement 3",(i/length(df1$Time)*100),"%"))
      }
    }
  }
  
  print("traitement 3 done %")
  dfx<-na.rm(df1) 
  return(dfx) # Renvoi le contenue de la dataframe en sortie de la fonction.
}

Time_and_delta_Time<-function(data) # Fonction qui permets de créer un second axe temporelle ou T0 commences au moment de l'injection des drogues
{
  data$Delta<-as.numeric(data$Delta) # Définie la colonne delta comme numérique
  print("traitement delta time start")
  data$Real_Minutes<-NA # Création de la colonne Real_minutes 
  data$Real_Second<-NA # Création de la colonne Real_Second
  data$Delta_Minutes<-NA # Création de la colonne Delta_minutes
  data$Delta_Second<-NA # Création de la colonne Delta_Second
  
  
  data$Real_Minutes<-floor(data$Real_Time) # Enregistre dans la colonne Real_minutes les minutes contenues dans Real_Time (4min et 16 seconde donc 4 minutes)
  data$Real_Second<-data$Real_Time-data$Real_Minutes # Enregistre dans la colonne Real_Seconde les secondes contenues dans Real_Time (4min et 16 secondes donc 16 secondes)
  
  data$Delta_Minutes<-floor(data$Delta) # Enregistre dans la colonne Delta_minutes les minutes contenues dans Delta (4min et 16 seconde donc 4 minutes)
  data$Delta_Second<-data$Delta-data$Delta_Minutes # Enregistre dans la colonne Delta_Seconde les secondes contenues dans Delta (4min et 16 secondes donc 16 secondes)
  
  data$Delta_Second<-round(data$Delta_Second*100/60,2) # Passe les secondes de base 60 en base 100  (pour éviter les erreurs dans le calcul suivant)
  data$Real_Second<-round(data$Real_Second*100/60,2) # Passe les secondes de base 60 en base 100 ( pour éviter les erreurs dans le calcul suivant)
  
  data$Real_Time<-data$Real_Minutes+data$Real_Second # Enregistre dans Real_time la somme de Real_Minutes et Real_Secondes 
  data$Delta<-data$Delta_Minutes+data$Delta_Second # Enregistre dans Delta_time la somme de Delta_Minutes et Delta_Secondes 
  
  data$Time<-data$Real_Time-(5+data$Delta) # Enregistre dans la colonne Time le calcul suivant (Real_Time - (5+Delta_Time)) afin d'obtenir un axe où T0 commence au moment de l'injection de la drogue  
  
  print("traitement delta time done")
  return(data)
}

Pourcentage_Vers_PT0<-function(data_save) # Calcul de pourcentage de vers restant dans les puces microfluidiques depuis le démarrage de la pompe péristaltique (aligner sur l'axe temporel Real_Time)
{
  print("traitement PT0 Start")
  data$Vers<-NA # Création d'une colonne Vers pour enregistrer le nombre de parasites dans les puces au cours du temps
  data$Vers<-data$Male+data$Female+data$Couple*2 # Calcule le pourcentage de vers restants dans les puces microfluidiques au cours du temps
  data$VersP100<-NA # Création d'une colonne pour enregistrer le pourcentage de vers au cours du temps (aligner sur l'axe temporel Real_Time)  
  data_save$PT0_Vers<-NA # Création d'une colonne pour enregistrer le pourcentage de vers au cours du temps
  
  data<-data_save # Enregistre les données présentes dans data_save dans data
  
  male<-NA # Création d'une variable male
  female<-NA # Création d'une variable female
  PT0_couple<-NA # Création d'une variable PT0_couple 
  vers<-NA # Création d'une variable vers
  PT0max<-NA # Création d'une variable PT0max (enregistre le nombre de vers max à partir de T0)
  PT5Dmax<-NA # Création d'une variable PTDmax (enregistre le nombre de vers max à partir de T0+5minutes+delta)
  
  data$Number<-as.numeric(data$Number) # Défini la colonne Number comme numérique
  
  for (i in 1:length(data$Real_Time)) # Pour chaque ligne de data
  {
    if((data$Real_Time[i]==0)&&(data$Flow[i]==0.05)) # Si la ligne correspond à T0 et que le flux est à 0.05 mL/min (débit au démarrage de la pompe péristaltique)
    {
      data$Male[i]<-data$Male_out[i] # Enregistre le contenue de Male_out dans la ligne i de la colonne Male
      data$Female[i]<-data$Female_out[i] # Enregistre le contenue de Female_out dans la ligne i de la colonne Female
      data$Couple[i]<-data$Couple_out[i] # Enregistre le contenue de Couple_out dans la ligne i de la colonne Couple
      
      male<-data$Male[i] # Enregistre la valeur initiale de parasites mâle (non en couple) dans la variable mâle 
      female<-data$Female[i] # Enregistre la valeur initiale de parasites femelle (non en couple) dans la variable Female
      couple<-data$Couple[i] # Enregistre la valeur initiale de parasites en couple dans la variable couple
      data$Vers[i]<-data$Male[i]+data$Female[i]+data$Couple[i]*2 # Calcul le nombre de vers dans la puce microfluidique
      data$PT0_Vers[i]<-(data$Vers[i]/data$Vers[i])*100 # Calcul le pourcentage de vers dans la puce microfluidique
      
      for (j in 1:length(data$Real_Time)) # Pour chaque ligne de data
      {
        if(j!=i) # Si la ligne j de data est différente de la ligne i de data
        {
          if(data$Number[j]==data$Number[i]) # Si le numéro de puce est le même que pour les lignes i et j
          {
            
            if((data$Dose[i]==data$Dose[j]))#&&(data$Drug[i]==data$Drug[j])) # Si la dose est la même pour les lignes i et j
            {
              data$Male[j]<-male-data$Male_out[j] # Enregistre dans la ligne j de la colonne mâle la différence entre le nombre initiale de mâle et le nombre de mâle Sortie (Male_Out)
              data$Female[j]<-female-data$Female_out[j] # Enregistre dans la ligne j de la colonne Female la différence entre le nombre initiale de femelle et le nombre de Femelle Sortie (Female_Out)
              data$Couple[j]<-couple-data$Couple_out[j] # Enregistre dans la ligne j de la colonne Couple la différence entre le nombre initiale de couple et le nombre de Couple Sortie (Couple_Out)
              data$Vers[j]<-data$Male[j]+data$Female[j]+data$Couple[j]*2 # Enregistre dans la ligne j de la colonne mâle le nombre de vers presents dans la puce
              male<-data$Male[j] # Enregistre dans la variable mâle le nombre de mâle restants dans la puce
              female<-data$Female[j] # Enregistre dans la variable female le nombre de femmelle restantes dans la puce
              couple<-data$Couple[j] # Enregistre dans la variable couple le nombre de couple restantes dans la puce 
              data$PT0_Vers[j]<-(data$Vers[j]/data$Vers[i])*100 # Calcul le pourcentage de vers restant dans la puce 
            }
          }
        }
      }
    }
    #print(paste("traitement PT0",(i/length(data$Time)*100),"%"))
  }
  print("traitement PT0 Done")
  return(data) # renvoi le contenue de data en sortie de la fonction
}

Pourcentage_Vers_PT5D<-function(data) # Calcul de pourcentage de vers restant dans les puces microfluidiques après le lavage et l'injection de la drogue (aligner sur l'axe temporel Time)
{
  print("traitement PT5D Start")
  data$PT5D_Vers<-NA
  for (i in 1:length(data$Time)) # Pour chaque ligne de data
  {
    if(data$Time[i]<0) # Si le temps contenue dans la ligne i de Time est inférieurs à zéro (avant l'introduction de la drogue) 
    {
      data$PT5D_Vers[i]<-100 # Le pourcentage de vers dans la ligne i de PT5D (Pourcentage de vers à Temps après 5 minutes + Delta) est égale a 100%
    }
    if((data$Time[i]==0)&&(data$Flow[i]==1)) # Si la ligne correspond à T0 (après lavage de 5 minutes et introduction de la drogue) et que le flux est à 1 mL/min (débit au moment de l'introduction de la drogue) 
    {
      data$PT5D_Vers[i]<-(data$Vers[i]/data$Vers[i])*100 # Enregistre dans PT5D le pourcentage de vers présent dans la puce
      
      for (j in 1:length(data$Time)) # Pour chaque ligne de data
      {
        if(j!=i&&data$Time[j]>=0) # Si la ligne j de data est différente de la ligne i de data et que le temps à la ligne j de Time est supérieure ou égale à 0
        {
          if(data$Number[j]==data$Number[i]) # Si le numéro de puce est le même que pour les lignes i et j
          {
            if((data$Dose[i]==data$Dose[j]))#&&(data$Drug[i]==data$Drug[j])) # Si la dose est la même pour les lignes i et j
            {
              data$PT5D_Vers[j]<-(data$Vers[j]/data$Vers[i])*100 # Calcul le pourcentage de vers restant dans la puce
            }
          }
        }
      }
    }
    #print(paste("traitement PT5D",(i/length(data$Time)*100),"%"))
  }
  print("traitement PT5D Done")
  return(data) # renvoi de contenue de data en sortie de la fonction
}

Ecart_type<-function(dfx,NdataNumber,Ndose) #Calcul les écarts types pour expériences
{
  
  for (i in 1:max(NdataNumber)) # Pour chaque numéro d'expériences 
  {
    assign(paste("d", Ndose[i], sep = ""),dfx[dfx$Dose %in% c(Ndose[i]),]) # Eclate les données pour chaque dose en plusieur dataframe en les nommant "d"+valeur de la dose (exemple : d900, d700, d500 ...)
    for (j in 1:length(Ndose)) # Pour chaque dose
    {
      assign(paste("Split",Ndose[j]+i/1000,sep=""),get(paste("d", Ndose[j], sep = ""))[get(paste("d", Ndose[j], sep = ""))$Number %in% c(i),]) # Eclate les données de chaque dose (exemple : d900) pour chaque numéro d'expériences en les nommant ("Split"+Valeur de la dose + valeur de l'expérience /1000) (exemple : Split900.0001, Split900.0002, Split900.0003 ....)
      assign(paste("d",Ndose[j]+i/1000,sep=""),get(paste("d", Ndose[j], sep = ""))[get(paste("d", Ndose[j], sep = ""))$Number %in% c(i),]) # Eclate les données pour chaque dose (exemple : d900) pour chaque numéro d'expériences (exemple : d900.0001, d900.0002, d900.0003 ....)
    }
  }
  
  d<-get(paste("Split", Ndose[1]+1/1000, sep = "")) # Enregistre dans la dataframe d les données contenues dans Split Dose N°1 Expériences N°1
  if(length(d$Time)==0) # Si la dataframe d ne contient rien car il n'existe pas de test à une dose de 900 dans l'expérience N°1 (pas de données pour d900.0001 par exemple), recherche une autre dataframe où des données sont contenues
  {
    A<-1 # La variable A est égale à 1 (pour les doses)
    B<-1 # La variable B est égale à 1 (pour les Numéro d'expériences)
    while(length(d$Time)==0) # Temps que la dataframe d ne contient rien 
    {
      while(is.na(Ndose[B])) # Temps que la valeur contenue dans la case B du tableau Ndose est égale à NA
      {
        B<-B+1 # B = B+1
        if(B>max(Ndose)) # Si B est plus grand que la longueur du tableau NDose
        {
          B<-1 # On retourne à la valeur initiale du tableau Ndose (B=1)
          A<-A+1 # On change de numéro d'expérience (A=A+1)
        }
      }
      d<-get(paste("Split", Ndose[B]+A/1000, sep = "")) # Enregistre le contenue de la dataframe Split Dose B numéro A/1000 (exemple : Split900.0005 si Ndose = 900 et expérience=5)
      B<-B+1 # B=B+1
      if(B>max(NdataNumber)) # Si B est plus grand que la longueur du tableau NDose
      {
        B<-1 # On retourne à la valeur initiale du tableau Ndose (B=1)
        A<-A+1 # On change de numéro d'expérience (A=A+1)
      }
    }
  } # Une fois que d contient des données
  
  d$Number<-NA # Ajoute une colonne Number dans la dataframe d
  d$Dose<-NA # Ajoute une colonne Number dans la dataframe d
  d$delta_Vers_PT5D<-NA # Ajoute une colonne delta_Vers_PT5D dans la dataframe d
  d$delta_Vers_PT0<-NA # Ajoute une colonne delta_Vers_PT0 dans la dataframe d
  d$Mean_Vers_PT5D<-NA # Ajoute une colonne Mean_Vers_PT5D dans la dataframe d
  d$Mean_Vers_PT0<-NA # Ajoute une colonne Mean_Vers_PT0 dans la dataframe d
  
  for (i in 1:length(Ndose)) # Pour chaque dose 
  {
    assign(paste("Dose",Ndose[i],sep=""),get(paste("d"))) # Crée autant de dataframe Dose qu'il y a de dose dans le tableau Ndose en y copiant les valeur contenues dans d
  }
  
  TablePT5D<-NA # Crée un tableau TablePT5D
  TablePT0<-NA # Crée un tableau TablePT0
  for (i in 1:max(NdataNumber)) # Pour chaque numéro d'expériences
  {
    TablePT5D[i]<-NA # Défini la longueur du TablePT5D comme celle de NdataNumber puis donne les valeurs NA à chaque colonne
    TablePT0[i]<-NA # Défini la longueur du TablePT5D comme celle de NdataNumber puis donne les valeurs NA à chaque colonne
  }
  
  for(i in 1:length(Ndose)) # Pour chaque dose
  {
    for(j in 1:length(d$Time)) # Pour chaque Temps
    {
      for(k in 1:max(NdataNumber)) # Pour chaque numéro d'expérience
      {
        TablePT5D[k]<-get(paste("Split", Ndose[i]+k/1000, sep = ""))$PT5D_Vers[j] # Récupère la dataframe Split + Dose i + Numéro  (exemple Split900.0001) puis enregistre la valeur contenue dans colonne PT5D_Vers dans la colonne K de TablePT5D (exemple : enregistre dans TablePT5D colonne 1, la valeur contenue dans la ligne 1 de la colonne PT5D_Vers de Split900.0001)
        TablePT0[k]<-get(paste("Split", Ndose[i]+k/1000, sep = ""))$PT0_Vers[j] # Récupère la dataframe Split + Dose i + Numéro  (exemple Split900.0001) puis enregistre la valeur contenue dans PT5D_Vers dans la colonne K de TablePT5D (exemple : enregistre dans TablePT5D colonne 1 la valeur contenue dans la ligne 1 de la colonne PT5D_Versde Split900.0001)
      }
      d$Mean_Vers_PT5D[j]<-mean(TablePT5D,na.rm = TRUE) # Enregistre dans la ligne j de la colonne Mean_Vers_PT5D la valeur moyenne de TablePT5D
      d$delta_Vers_PT5D[j]<-sd(TablePT5D,na.rm = TRUE) # Enregistre dans la ligne j de la colonne delta_Vers_PT5D la valeur de l'écart type de TablePT5D
      
      d$Mean_Vers_PT0[j]<-mean(TablePT0,na.rm = TRUE) # Enregistre dans la ligne j la valeur moyenne de TablePT0
      d$delta_Vers_PT0[j]<-sd(TablePT0,na.rm = TRUE) # Enregistre dans la ligne j la valeur l'écart type de TablePT0
      
      d$Dose<-Ndose[i] # Enregistre la valeur de la dose contrôlée dans la colonne Dose de la dataframe d
      assign(paste("Dose",Ndose[i],sep=""),get(paste("d"))) # Enregistre le contenue de d dans une dataframe nommées Dose + Numéro de Dose
    }
    print(paste("traitement ecart type",(i/length(Ndose)*100),"%"))
  }
  print("traitement ecart type 100 %")
  
  #d<-get(paste("Dose", Ndose[i], sep = "")) # Après le traitement, 
  d<-NA # Vide le contenue de d
  
  for(i in 1: length(Ndose)) # Pour chaque dose
  {
    d<-rbind(d,get(paste("Dose", Ndose[i], sep = ""))) # Rassemble dans d toutes les données dataframe Dose+i
  }
  d<-subset(d, select = -c(Number,Male_out,Female_out,Couple_out,Video,Drug,Male,Female,Couple,Real_Minutes,Real_Second,Delta_Minutes,Delta_Second,PT0_Vers,Vers,PT5D_Vers)) # Supprime le colonnes inutiles de la dataframe d
  return(d) # Renvoi la dataframe d en sortie de la fonction
}

Suppression_dose<-function(Suppression,testvalDose,data) # Pour la dataframe testée, concerve qu'une seule dose si la fonction suppression est activée
{
  if(Suppression==1) # Si la fonction Suppression est activé
  {
    for(i in 1:length(data$Real_Time)) # Pour chaque temps
    {
      if(data$Dose[i]!=testvalDose) # Si la dose présentes dans la ligne i est différente de la dose à conserver
      {
        data$Dose[i]<-NA # La valeur de la dose contenue dans la ligne i de la colonne dose devient NA
      }
    }
    data<-na.rm(data) # Supprimes tous les lignes contenant un NA
  }
  return(data) # Renvoi le contenue de data en sortie de la fonction
}


Reduction_Du_Pas<-function(d,N) # Reduit le pas entre 2 points des échantillons pour éviter un graphique surcharger (car sinon, on a un point tous les 0.01 secondes)
{
  print("traitement reduction du pas Start")
  d2<-na.rm(d) # Supprimes les lignes contenants des NA
  A<-1 # Initialise une variable A à 1
  
  for(i in 1:length(d2$Time)) # Pour chaque temps 
  {
    if(A!=1) # Si A est différent de 1
    {
      d2$Time[i]<-NA # La ligne i de la colonne Time de d2 devient NA
    }
    A<-A+1 # Ajoute +1 à A
    
    if(A==N) # Si A est égale à N (N représente l'écart entre deux points) 
    {
      A<-1 # A est réinitialisé a 1
      #print(paste("traitement reduction du pas ",(i/length(d$Time)*100),"%"))
    }
  }
  print("traitement reduction du pas finished")
  d2<-na.rm(d2) # Supprimes les lignes de d2 contenant des NA
  return(d2) # Renvoi le contenue de d2 en sortie de la fonction
}


p100<-function(data) # Supprime les colonne 1 est 4 de la dataframe
{
  data<-subset(data,select = -c(1,4))
  return(data)
}

Print_graph_Time_PT5D<-function(data,Drug_used,Unit) # Trace l'évolution du pourcentage de vers PT5D en fonction du temps
{
  data$Unit<-Unit
  data$Drug_Used<-Drug_used
  A<-ggplot(data,aes(x=Time,y=Mean_Vers_PT5D,group=Dose,color=Dose,na.rm = TRUE))+
    geom_point(aes(shape=Dose),size=1)+
    scale_color_grey(start=0.7, end=0,name = paste("[",data$Drug_Used[1],"]",data$Unit[1]))+
    scale_shape_manual(values=c(11,12,13,14,15,16,17,18,19,20,21),name = paste("[",data$Drug_Used[1],"]",data$Unit[1]))+
    geom_errorbar(aes(ymin=Mean_Vers_PT5D-SDmin5D, ymax=Mean_Vers_PT5D+SDmax5D))+
    geom_step(aes(group=Dose))+
    #scale_y_continuous( limits=c(0, 100))+
    scale_x_continuous( limits=c(0, 45),breaks=seq(0, 45, by = 5))+
    theme_bw() + 
    theme(plot.title = element_text(size=7,face = "bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=6,face = "bold",colour = "black"),
          axis.text.y = element_text(size=6,face = "bold",colour = "black"),
          axis.title.x = element_text(size=6,face = "bold"),
          axis.title.y = element_text(size=6,face = "bold"),
          legend.title = element_text(size=4,face = "bold"),
          legend.text = element_text(face = "bold",size = 4),
          legend.position ="bottom",
          legend.direction="horizontal",
          legend.key.size = unit(0.5, "lines"))+
    guides(col = guide_legend(ncol = 13))+
    ggtitle(paste("Worms adhesion depending on",data$Drug_Used[1],"dose"))+
    labs(y= "Number of worms (%)", x = "Time (min)",colour = paste("Dose(",data$Unit[1],")"))
  
  ggsave(filename = file.path("figs",paste("Microfluidic Worms adhesion depending on",data$Drug_Used[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)
  return(A)
}

Print_graph_Time_PT0<-function(data,Drug_used,Unit) # Trace l'évolution du pourcentage de vers PT0 en fonction du temps
{
  data$Unit<-Unit
  data$Drug_Used<-Drug_used
  A<-ggplot(data,aes(x=Real_Time,y=Mean_Vers_PT0,group=Dose,color=Dose,na.rm = TRUE))+
    geom_point(aes(shape=Dose),size=2)+
    scale_color_grey(start=0.7, end=0,name = paste("[",data$Drug_Used[1],"]",data$Unit[1]))+
    scale_shape_manual(values=c(11,12,13,14,15,16,17,18,19,20,21),name = paste("[",data$Drug_Used[1],"]",data$Unit[1]))+
    geom_errorbar(aes(ymin=Mean_Vers_PT0-delta_Vers_PT0, ymax=Mean_Vers_PT0+delta_Vers_PT0))+
    geom_step(aes(group=Dose))+
    #scale_y_continuous( limits=c(0, 100))+
    scale_x_continuous( limits=c(0, 52.5),breaks=seq(0, 52.5, by = 5))+
    theme_bw() + 
    theme(plot.title = element_text(size=7,face = "bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=6,face = "bold",colour = "black"),
          axis.text.y = element_text(size=6,face = "bold",colour = "black"),
          axis.title.x = element_text(size=6,face = "bold"),
          axis.title.y = element_text(size=6,face = "bold"),
          legend.title = element_text(size=4,face = "bold"),
          legend.text = element_text(face = "bold",size = 4),
          legend.position ="bottom",
          legend.direction="horizontal",
          legend.key.size = unit(0.5, "lines"))+
    ggtitle(paste("Worms adhesion (PT0) depending on",data$Drug_Used[1],"dose"))+
    labs(y= "Number of worms (%)", x = "Time (min)",colour = paste("Dose(",data$Unit[1],")"))
  
  ggsave(filename = file.path("figs",paste("Worms adhesion (PT0) depending on",data$Drug_Used[1],"dose.TIFF")), width = 15, height = 10, units = "cm",dpi = 600)
  return(A)
}

Time_selection_IC50<-function(data,data_save2,start_time,gap) 
{
  Ndose<-N_of_variable(data_save2$Dose)
  compteur<-start_time 
  ligne<-1
  tableau <- data.frame(Dose=c(1:length(data$Time)), p_5=length(data$Time),Time=length(data$Time))
  tableau$Dose<-NA
  tableau$p_5<-NA
  tableau$Time<-NA
  while(compteur<=52.5)
  {
    for(i in 1:length(Ndose))
    {
      for(j in 1:length(data$Time))
      {
        if(data$Time[j]==compteur && data$Dose[j]==Ndose[i])
        {
          tableau$Dose[ligne]<-data$Dose[j]
          tableau$p_5[ligne]<-data$Mean_Vers_PT5D[j]
          tableau$Time[ligne]<-data$Time[j]
          tableau$delta[ligne]<-data$delta_Vers_PT5D[j]
          ligne<-ligne+1
        }
      }
    }
    compteur<-compteur+gap
  }
  
  tableau<-na.rm(tableau)
  tableau$Dose<-as.numeric(tableau$Dose)
  tableau$p_5<-as.numeric(tableau$p_5)
  tableau$Time<-as.factor(tableau$Time)
  return(tableau)
}

Print_graph_log_IC50<-function(tableau,Drug_used,Unit) # Trace l'évolution du pourcentage de vers PT0 en fonction de la dose est du temps en échelle logarithmique
{
  data$Unit<-Unit
  data$Drug_Used<-Drug_used
  A<-ggplot(tableau,aes(x=Dose,y=p_5,group=Time,color=Time,na.rm = TRUE))+
    geom_point(aes(shape=Time),size=2)+
    scale_color_grey(start=0.7, end=0,name = "Time (min)")+
    scale_shape_manual(values=c(11,12,13,14,15,16,17,18,19,20,21,22,23),name = "Time (min)")+
    geom_errorbar(aes(ymin=p_5-delta, ymax=p_5+delta))+
    geom_line()+
    #scale_x_log10(limits=c(1, 2.65), breaks=seq(1, 2.65 , by = 0.2))+
    scale_x_log10(limits=c(0, 1000))+
    annotation_logticks(sides = "b")+
    theme_bw() + 
    annotation_logticks(sides = "b")+scale_x_log10(limits=c(1, 1000))+
    theme(plot.title = element_text(size=7,face = "bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=6,face = "bold",colour = "black"),
          axis.text.y = element_text(size=6,face = "bold",colour = "black"),
          axis.title.x = element_text(size=6,face = "bold"),
          axis.title.y = element_text(size=6,face = "bold"),
          legend.title = element_text(size=4.2,face = "bold"),
          legend.text = element_text(face = "bold",size = 4),
          legend.position = c(0.87, 0.5),
          legend.direction="vertical",
          legend.key.size = unit(0.5, "lines"),
          legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
    ggtitle(paste("Microfluidic IC50 Worms adhesion depending on",data$Drug_Used[1],"dose"))+labs(y= "Number of worms (%)", x = paste("dose(",data$Unit[1],")"),colour = paste("log10(dose(",data$Unit[1],"))"))
  A
  ggsave(filename = file.path("figs",paste("Microfluidic IC50 Worms adhesion depending on",data$Drug_Used[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)
  return(A)
}

Print_graph_nolog_IC50<-function(tableau,Drug_used,Unit) # Trace l'évolution du pourcentage de vers PT0 en fonction de la dose est du temps
{
  data$Unit<-Unit
  data$Drug_Used<-Drug_used
  A<-ggplot(tableau,aes(x=Dose,y=p_5,group=Time,color=Time,na.rm = TRUE))+
    geom_point(aes(shape=Time),size=2)+
    scale_color_grey(start=0.7, end=0,name = "Time (min)")+
    scale_shape_manual(values=c(11,12,13,14,15,16,17,18,19,20,21,22,23),name = "Time (min)")+
    geom_line()+
    geom_errorbar(aes(ymin=p_5-delta, ymax=p_5+delta))+
    scale_x_continuous( limits=c(min(tableau$Dose)-20, max(tableau$Dose))+10,breaks=seq(min(tableau$Dose), max(tableau$Dose), by = 50))+
    theme_bw() + 
    theme(plot.title = element_text(size=7,face = "bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=6,face = "bold",colour = "black"),
          axis.text.y = element_text(size=6,face = "bold",colour = "black"),
          axis.title.x = element_text(size=6,face = "bold"),
          axis.title.y = element_text(size=6,face = "bold"),
          legend.title = element_text(size=4.2,face = "bold"),
          legend.text = element_text(face = "bold",size = 4),
          legend.position = c(0.87, 0.5),
          legend.direction="vertical",
          legend.key.size = unit(0.5, "lines"),
          legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
    ggtitle(paste("Microfluidic IC50 Worms adhesion depending on",data$Drug_Used[1],"dose"))+labs(y= "Number of worms (%)", x = paste("dose(",data$Unit[1],")"),colour = paste("log10(dose(",data$Unit[1],"))"))
  A
  ggsave(filename = file.path("figs",paste("Microfludic IC50 Worms adhesion depending on",data$Drug_Used[1],"dose (no log).TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)
  return(A)
}



#code -----------------------------------


data_save$Real_Time<-data_save$Real_Time+(60*(data_save$Video-1)) # Comme 1 vidéo est égale à 60 minutes, calcule de temps réel de chaque expérience en fonction du nombre d'enregistrement
data<-data_save # Sauvegarde les données dans data

Suppression<-1 # Active la variable suppression pour un temps de traitement plus rapide

Ndose_to_test<-N_of_variable(data$Dose) #Récupère le nombre de dose présent de le jeu de donnée
Ndose_real<-Ndose_to_test # Sauvegarde ces valeurs dans la variable Ndose_real
NdataNumber<-N_of_variable(data$Number) #Récupère le nombre de dose présent dans le jeu de donnée

Drug_used<-data_real$Drug_Used[1] # Récupère le nom de la drogue utilisée dans le jeu de donnée
Unit<-data_real$Unit[1] # Récupère l'unité métrique de la drogue utilisée durant les expériences

data$Delta<-as.numeric(data$Delta) # Définie la colonne Delta comme numérique

data_save2<-data # Sauvegarde le jeu de donné présent dans data dans data_save2
data2<-NA # Définie data2 comme NA
tableau_normalise<-NA # Définie tableau_normalise comme NA

for(i in 1:length(Ndose_real)) # Pour chaque dose
{
  data<-data_save2 # Recharge les données complète dans data
  
  if(Suppression==0) # Si la variable suppression est inactive
  {
    i<-length(Ndose_to_test) # Réalise le calcul de la boucle qu'une fois mais allonge considérablement le temps d'exécution des calcules
  }
  data<-Suppression_dose(Suppression,Ndose_real[i],data) # Si la variable suppression est active, conserve uniquement les données pour une valeur de dose défini (ex : conserve que les dose à 0 nM de PZQ)
  Ndose_to_test<-N_of_variable(data$Dose) # Affiche la valeur de la dose a tester après la suppression et l'enregistre dans Ndose_to_test (Exemple : reste que les doses à 0 nM)
  
  
  if(length(N_of_variable(data$Number))>=2) # Si cette dose a été testé au minimum en duplicate, exécute la boucle. Sinon, oublie le jeu de données  
  {
    NdataNumber<-N_of_variable(data$Number) # NdataNumber enregistre dans un tableau le nombre de fois ou la dose a été testé
    tableau_temps<-Vecteur_time(data,60,Ndose_to_test) # Création d'un tableau allant de 0 à 60 minutes avec un pas de 0.01 seconde  
    tableau_temps<-Vecteur_Dose_Number(tableau_temps,data,Ndose_to_test) # Concatène le tableau avec lui-même autant de fois qu'il y a de dose (fonction que si il n'y a pas eu de suppression) 
    df<-Placement_Data_In_Vecteur(tableau_temps,data) # Place les données de data dans le tableau temps  
    df$Real_Time<-df$Time
    data<-Normalisation_sortie_vers(df) #De T0 à T60 minutes, complète chaque ligne vide du tableau par les valeurs près remplit par la fonction préssédente  
    data<-Time_and_delta_Time(data) # Aligne les axes temporels de chaque expérience au moment de l'ajout de la drogue 
    data<-Pourcentage_Vers_PT0(data) # Calcul le pourcentage de vers restant dans les puces à partir du moment où la pompe est mise en route
    data<-Pourcentage_Vers_PT5D(data) # Calcul de pourcentage de vers restant dans les puces à partir du moment où la drogue est ajoutée dans les puces microfluidiques
    tableau_normalise<-rbind(tableau_normalise,data) # Ajoute dans tableau_normalise les données présente dans data  
    d<-Ecart_type(data,NdataNumber,Ndose_to_test) # Calcul les écarts types pour chaque expérience et l'enregistre dans d
    data2<-rbind(data2,d) # Ajoute dans data2 les valeur contenue dans d
  }
  
  print(paste("traitement total",(i/length(Ndose_real)*100),"%"))
}
data<-na.rm(data2) # Supprimes les lignes contenants des valeurs NA dans data2 et l'enregistre dans data
data_save3<-data
data<-data_save3
data$Dose<-as.factor(data$Dose)
data<-Reduction_Du_Pas(data,31) # Réduit le pas entre 2 points pour tous les échantillons
data_reduit<-data

data$Ymax5D<-NA
data$Ymin5D<-NA
data$SDmin5D<-NA
data$SDmax5D<-NA

for(i in 1:length(data$Mean_Vers_PT5D))
{
  
  data$Ymax5D[i]<-data$Mean_Vers_PT5D[i]+data$delta_Vers_PT5D[i]
  data$Ymin5D[i]<-data$Mean_Vers_PT5D[i]-data$delta_Vers_PT5D[i]
  if(data$Ymax5D[i]>100)
  {
    data$SDmax5D[i]<-100-data$Mean_Vers_PT5D[i]
  }
  else
  {
    data$SDmax5D[i]<-data$delta_Vers_PT5D[i]
  }
  if(data$Ymin5D[i]<0)
  {
    data$SDmin5D[i]<-sqrt((0-data$Mean_Vers_PT5D[i])^2)
  }
  else
  {
    data$SDmin5D[i]<-data$delta_Vers_PT5D[i]
  }
}
for(i in 1:length(data$Time))
{
  if(data$Time[i]>50)
  {
    data$Time[i]<-NA
  }
}
data<-na.rm(data)

Xlab<-expression (bold(Time~"(min)"))
Ylab<-expression (bold(Attached~worms~"(%)"))

ggplot(data,aes(x=Time,y=Mean_Vers_PT5D,group=Dose,color=Dose,na.rm = TRUE))+
  geom_step()+
  geom_errorbar(aes(ymin=Mean_Vers_PT5D-SDmin5D, ymax=Mean_Vers_PT5D+SDmax5D),width=0.085)+
  geom_point(aes(shape=Dose),size=1.5)+
  scale_shape_manual(values=c(11,12,13,14,15,16,17,18,19,20,21),name = paste("[",Drug_used,"]",Unit))+
  scale_y_continuous( limits=c(0, 100),breaks=seq(0, 100, by = 20))+
  scale_x_continuous(limits=c(0, 55),breaks=seq(0, 50, by = 5))+
  theme_bw() + 
  theme(panel.border = element_rect(fill = "transparent", color = "black", size = 0.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.position = c(0.94, 0.5),
        legend.background = element_rect(fill = "transparent"),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"))+
  labs(y= Ylab, x = Xlab,color=paste("[",Drug_used,"]",Unit))

ggsave(dpi=600,filename = file.path("figs",paste("Figure Microfluidic Worms adhesion depending on",Drug_used,"dose.TIFF")), width = 16.6, height = 8.3, units = "cm")

ggplot(data,aes(x=Time,y=Mean_Vers_PT5D,group=Dose,color=Dose,na.rm = TRUE))+
  geom_step()+
  geom_errorbar(aes(ymin=Mean_Vers_PT5D-SDmin5D, ymax=Mean_Vers_PT5D+SDmax5D),width=0.085)+
  geom_point(aes(shape=Dose),size=1.5)+
  scale_shape_manual(values=c(11,12,13,14,15,16,17,18,19,20,21),name = paste("[",Drug_used,"]",Unit))+
  scale_y_continuous( limits=c(0, 100),breaks=seq(0, 100, by = 20))+
  scale_x_continuous(limits=c(0, 55),breaks=seq(0, 50, by = 5))+
  theme_bw() + 
  theme(panel.border = element_rect(fill = "transparent", color = "black", size = 0.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black",size=0.2),
        axis.text.x = element_text(size=7.8,face = "bold",colour = "black"),
        axis.text.y = element_text(size=7.8,face = "bold",colour = "black"),
        axis.title.x = element_text(size=9,face = "bold"),
        axis.title.y = element_text(size=9,face = "bold"),
        legend.title = element_text(size=8.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 7.2),
        legend.position = c(0.94, 0.5),
        legend.background = element_rect(fill = "transparent"),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"))+
  labs(y= Ylab, x = Xlab,color=paste("[",Drug_used,"]",Unit))

ggsave(dpi=600,filename = file.path("figs",paste("Figure Microfluidic Worms adhesion depending on",Drug_used,"dose.pdf")), width = 16.6, height = 8.3, units = "cm")




data<-na.rm(data2) # Supprimes les lignes contenants des valeurs NA dans data2 et l'enregistre dans data
data_save3<-data
data<-data_save3
data$Dose<-as.factor(data$Dose)
data<-Reduction_Du_Pas(data,31) # Réduit le pas entre 2 points pour tous les échantillons
data_reduit<-data


Print_graph_Time_PT0(data,Drug_used,Unit) # Trace le graphique à partir du nombre de vers retenu dans les puce après le démarrage de la pompe 
Print_graph_Time_PT5D(data,Drug_used,Unit) # Trace le graphique à partir du nombre de vers retenu après l'ajout de la drogue


data<-data_save3

tableau<-Time_selection_IC50(data,data_save2,0,5)

Print_graph_nolog_IC50(tableau,Drug_used,Unit)
for(i in 1:length(tableau$Dose))
{
  if(tableau$Dose[i]==0)
  {
    tableau$Dose[i]<-1
  }
}
Print_graph_log_IC50(tableau,Drug_used,Unit)

data<-data_save3

data<-na.rm(tableau_normalise)

if(Drug_used=="ART") # Si la drogue est de l'ART, enregistre les données contenues dans tableau_normalise dans des fichiers Excel
{
  data0<-data[data$Time %in% c(0),] #enregistre dans data0 que les données à temps 0
  data5<-data[data$Time %in% c(5),] #enregistre dans data5 que les données à temps 5
  data10<-data[data$Time %in% c(10),] #enregistre dans data10 que les données à temps 10
  data15<-data[data$Time %in% c(15),] #enregistre dans data15 que les données à temps 15
  data20<-data[data$Time %in% c(20),] #enregistre dans data20 que les données à temps 20
  data25<-data[data$Time %in% c(25),] #enregistre dans data25 que les données à temps 25
  data30<-data[data$Time %in% c(30),] #enregistre dans data30 que les données à temps 30
  data35<-data[data$Time %in% c(35),] #enregistre dans data35 que les données à temps 35
  data40<-data[data$Time %in% c(40),] #enregistre dans data40 que les données à temps 40
  data45<-data[data$Time %in% c(45),] #enregistre dans data45 que les données à temps 45
  data50<-data[data$Time %in% c(50),] #enregistre dans data50 que les données à temps 50
  
  data0<-subset(data0,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22)) # Conserve que les colonnes Dose et le pourcentage de vers restant (PT5D_Vers) pour chaque expérience
  data5<-subset(data5,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data10<-subset(data10,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data15<-subset(data15,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data20<-subset(data20,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data25<-subset(data25,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data30<-subset(data30,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data35<-subset(data35,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data40<-subset(data40,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data45<-subset(data45,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data50<-subset(data50,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  
  data0<-data0 %>% spread(key= Number, value= PT5D_Vers)
  data5<-data5 %>% spread(key= Number, value= PT5D_Vers)
  data10<-data10 %>% spread(key= Number, value= PT5D_Vers)
  data15<-data15 %>% spread(key= Number, value= PT5D_Vers)
  data20<-data20 %>% spread(key= Number, value= PT5D_Vers)
  data25<-data25 %>% spread(key= Number, value= PT5D_Vers)
  data30<-data30 %>% spread(key= Number, value= PT5D_Vers)
  data35<-data35 %>% spread(key= Number, value= PT5D_Vers)
  data40<-data40 %>% spread(key= Number, value= PT5D_Vers)
  data45<-data45 %>% spread(key= Number, value= PT5D_Vers)
  data50<-data50 %>% spread(key= Number, value= PT5D_Vers)
  
  write_xlsx(data0,path = "ART_µM_data_0_minutes.xlsx",col_names = TRUE)
  write_xlsx(data5,path = "ART_µM_data_5_minutes.xlsx",col_names = TRUE)
  write_xlsx(data10,path = "ART_µM_data_10_minutes.xlsx",col_names = TRUE)
  write_xlsx(data15,path = "ART_µM_data_15_minutes.xlsx",col_names = TRUE)
  write_xlsx(data20,path = "ART_µM_data_20_minutes.xlsx",col_names = TRUE)
  write_xlsx(data25,path = "ART_µM_data_25_minutes.xlsx",col_names = TRUE)
  write_xlsx(data30,path = "ART_µM_data_30_minutes.xlsx",col_names = TRUE)
  write_xlsx(data35,path = "ART_µM_data_35_minutes.xlsx",col_names = TRUE)
  write_xlsx(data40,path = "ART_µM_data_40_minutes.xlsx",col_names = TRUE)
  write_xlsx(data45,path = "ART_µM_data_45_minutes.xlsx",col_names = TRUE)
  write_xlsx(data50,path = "ART_µM_data_50_minutes.xlsx",col_names = TRUE)
}

if(Drug_used=="PZQ") # Si la drogue est du PZQ, enregistre les données contenues dans tableau_normalise dans des fichiers Excel
{
  data0<-data[data$Time %in% c(0),]
  data5<-data[data$Time %in% c(5),]
  data10<-data[data$Time %in% c(10),]
  data15<-data[data$Time %in% c(15),]
  data20<-data[data$Time %in% c(20),]
  data25<-data[data$Time %in% c(25),]
  data30<-data[data$Time %in% c(30),]
  data35<-data[data$Time %in% c(35),]
  data40<-data[data$Time %in% c(40),]
  data45<-data[data$Time %in% c(45),]
  data50<-data[data$Time %in% c(50),]
  
  data0<-subset(data0,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data5<-subset(data5,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data10<-subset(data10,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data15<-subset(data15,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data20<-subset(data20,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data25<-subset(data25,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data30<-subset(data30,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data35<-subset(data35,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data40<-subset(data40,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data45<-subset(data45,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  data50<-subset(data50,select = -c(1,2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22))
  
  data0<-data0 %>% spread(key= Number, value= PT5D_Vers)
  data5<-data5 %>% spread(key= Number, value= PT5D_Vers)
  data10<-data10 %>% spread(key= Number, value= PT5D_Vers)
  data15<-data15 %>% spread(key= Number, value= PT5D_Vers)
  data20<-data20 %>% spread(key= Number, value= PT5D_Vers)
  data25<-data25 %>% spread(key= Number, value= PT5D_Vers)
  data30<-data30 %>% spread(key= Number, value= PT5D_Vers)
  data35<-data35 %>% spread(key= Number, value= PT5D_Vers)
  data40<-data40 %>% spread(key= Number, value= PT5D_Vers)
  data45<-data45 %>% spread(key= Number, value= PT5D_Vers)
  data50<-data50 %>% spread(key= Number, value= PT5D_Vers)
  
  
  write_xlsx(data0,path = "PZQ_nM_data_0_minutes.xlsx",col_names = TRUE)
  write_xlsx(data5,path = "PZQ_nM_data_5_minutes.xlsx",col_names = TRUE)
  write_xlsx(data10,path = "PZQ_nM_data_10_minutes.xlsx",col_names = TRUE)
  write_xlsx(data15,path = "PZQ_nM_data_15_minutes.xlsx",col_names = TRUE)
  write_xlsx(data20,path = "PZQ_nM_data_20_minutes.xlsx",col_names = TRUE)
  write_xlsx(data25,path = "PZQ_nM_data_25_minutes.xlsx",col_names = TRUE)
  write_xlsx(data30,path = "PZQ_nM_data_30_minutes.xlsx",col_names = TRUE)
  write_xlsx(data35,path = "PZQ_nM_data_35_minutes.xlsx",col_names = TRUE)
  write_xlsx(data40,path = "PZQ_nM_data_40_minutes.xlsx",col_names = TRUE)
  write_xlsx(data45,path = "PZQ_nM_data_45_minutes.xlsx",col_names = TRUE)
  write_xlsx(data50,path = "PZQ_nM_data_50_minutes.xlsx",col_names = TRUE)
}
