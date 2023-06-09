#Initialise les librairies -----------------------------------

rm(list=ls()) # Supprime les données présentes dans la ram et la rom de RStudio

packages<-c("here","sp") # Chemin du répertoire avec les données 

for(i in packages){if(!require(i,character.only = T)){install.packages(i,dep=T)}}
lapply(packages, library, character.only = T) # Chargement de la liste des packages
library(ggplot2) # Charge la librairie ggplot2 (pour tracer les graphiques)
library(scales) # Charge la librairie scales
library(dplyr) # Charge la librairie dplyr
library(na.tools) # Charge la librairie na.tools
library(concatenate) #Charge la librairie concatenate
library(writexl) # Charge la librairie writexl
library(tidyr) # Charge la librairie tidyr

# Fonction  -----------------------------------

opentxt<-function(x) # Fonction utilisée pour ouvrir les fichiers textes
  
{
  data<-read.table(here("data",x),sep='\\t', header=T) # Ouvre et lit les fichiers txt pour les enregistrer dans une dataframe
}
txt<-dir("data",pattern = ".txt") # Enregistre le nombre de fichier texte et leur nom dans un tableau  
N<-length(txt) # Enregistre dans la variable N le nombre total de fichier texte
data_real<-structure(list(character()), class = "data.frame") # Création d’une variable t_jours de valeur NA 

for (i in 1:N) # Ouvre tous les fichiers textes et les enregistrent les uns après les autres dans la dataframes excel
{
  excel<-opentxt(txt[i]) # Ouvre les fichiers txt et les enregistres dans excel
  excel$Number<-i # Donne un numéro pour chaque fichier excel ouvert dans la colonne Number
  data_real<-rbind(data_real,excel) # Après chaque ouverture des fichier txt, les enregistres les données dans data_real
  
}
data_real<-na.rm(data_real) # Supprime les lignes de données contenants un NA

data_real$H<-as.numeric(data_real$H) # Défini la colonne H comme numérique (H pour Heures)


N_of_variable<-function(variable) # Fonction utilisée pour enregistrer dans un tableau l’ensemble des différentes valeurs existantes dans la colonne choisi de la dataframe en entrée de fonction (ici, Variable) (exemple : nom des milieux de culture, les jours pour chaque contrôles, …)
{
  compteur<-1 # Initialise le compteur à la valeur 1
  table<-NA # Création du tableau 
  table[compteur]<--1 # Donne la valeur -1 à la première colonne du tableau
  Num<-NA # Création de la variable Num avec comme valeur NA
  exist<-0 # Création de la variable exist a 0
  for (i in 1:length(variable)) # Pour chaque ligne contenue dans Variable 
  {
    Num<-variable[i] # Enregistre dans Num la valeur contenue dans la ligne i
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
  
  #Ndose<-rev(sort(table))
  Ndose<-sort(table) # Range les valeurs contenues dans le tableau dans décroissant
  return(Ndose) # Renvoi les données de la table en sortie de la fonction
}

error_detection<-function(data) # Fonction qui permets la détection d’erreurs lors de l'analyse des schistosomes
{
  data$Couple<-NA # Crée une colonne Couple 
  data$Couple_Female<-NA # Crée une colonne Couple_Female (nombre de femelles en couple) 
  data$Couple_Male<-NA # Crée une colonne Couple_Male (nombre de Males en couple) 
  data$E4<-NA # Crée une colonne E4 (Etat de viabilité 4) 
  data$E3<-NA # Crée une colonne E3 (Etat de viabilité 3) 
  data$E2<-NA # Crée une colonne E2 (Etat de viabilité 2) 
  data$E1<-NA # Crée une colonne E1 (Etat de viabilité 1) 
  data$E0<-NA # Crée une colonne E0 (Etat de viabilité 0) 
  data$Error<-0 # Crée une colonne erreur (Permets d'enregistrer si il y a une erreur de notation pour chaque ligne de données) 
  
  data$Couple_Male<-data$CME4+data$CME3+data$CME2+data$CME1+data$CME0 # Calcul le nombre de mâle en couple
  data$Couple_Female<-data$CFE4+data$CFE3+data$CFE2+data$CFE1+data$CFE0 # Calcul de nombre de femelle en couple
  data$Couple<-data$Couple_Male+data$Couple_Female # Calcul le nombre de vers en couple
  
  data$E4<-data$ME4+data$FE4 # Calcul le nombre de schistosomes à état E4 
  data$E3<-data$ME3+data$FE3 # Calcul le nombre de schistosomes à état E3 
  data$E2<-data$ME2+data$FE2 # Calcul le nombre de schistosomes à état E2 
  data$E1<-data$ME1+data$FE1 # Calcul le nombre de schistosomes à état E1 
  data$E0<-data$ME0+data$FE0 # Calcul le nombre de schistosomes à état E0 
  Control<-0 # Crée une variable Control de valeur zéro
  for(i in 1: length(data$DDP)) # Pour chaque ligne du jeux de donnée
  {
    if(data$Couple_Female[i]!=data$Couple_Male[i]) # Si le nombre de Femelle en couple est différent du nombre de mâle en couple
    {
      print(paste("problem couple ligne",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$E4[i]+data$E3[i]+data$E2[i]+data$E1[i]+data$E0[i]!=data$NST[i]) # Si la somme de schistosomes évalués (E4+E3+E2+E1+E0) est différent du nombre de schistosomes totaux (NST)
    {
      print(paste("problem status ligne",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$FE4[i]+data$FE3[i]+data$FE2[i]+data$FE1[i]+data$FE0[i]!=data$NSTF[i]) # Si la somme des Femelles schistosomes evaluées (E4+E3+E2+E1+E0) est différent du nombre de schistosomes Femelles totale (NSTF)
    {
      print(paste("problem status Female ligne",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$ME4[i]+data$ME3[i]+data$ME2[i]+data$ME1[i]+data$ME0[i]!=data$NSTM[i]) # Si la somme des Males schistosomes evalués (E4+E3+E2+E1+E0) est différent du nombre de schistosomes mâle totaux (NSTM)
    {
      print(paste("problem status Male",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CME4[i]>data$ME4[i]) # Si le nombre de mâle en couple à l'état E4 est supérieur au nombre de mâle à l'état E4
    {
      print(paste("problem status 4 and couple 4 Male",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CME3[i]>data$ME3[i]) # Si le nombre de mâle en couple à l'état E3 est supérieur au nombre de mâle à l'état E3
    {
      print(paste("problem status 3 and couple 3 Male",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CME2[i]>data$ME2[i]) # Si le nombre de mâle en couple à l'état E2 est supérieur au nombre de mâle à l'état E2
    {
      print(paste("problem status 2 and couple 2 Male",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CME1[i]>data$ME1[i]) # Si le nombre de mâle en couple à l'état E1 est supérieur au nombre de mâle à l'état E1
    {
      print(paste("problem status 1 and couple 1 Male",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CME0[i]>data$ME0[i]) # Si le nombre de mâle en couple à l'état E0 est supérieur au nombre de mâle à l'état E0
    {
      print(paste("problem status 0 and couple 0 Male",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    
    if(data$CFE4[i]>data$FE4[i]) # Si le nombre de femelle en couple à l'état E4 est supérieur au nombre de femelle à l'état E4
    {
      print(paste("problem status 4 and couple 4 Female",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CFE3[i]>data$FE3[i]) # Si le nombre de femelle en couple à l'état E3 est supérieur au nombre de femelle à l'état E3
    {
      print(paste("problem status 3 and couple 3 Female",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CFE2[i]>data$FE2[i]) # Si le nombre de femelle en couple à l'état E2 est supérieur au nombre de femelle à l'état E2
    {
      print(paste("problem status 2 and couple 2 Female",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CFE1[i]>data$FE1[i]) # Si le nombre de femelle en couple à l'état E1 est supérieur au nombre de femelle à l'état E1
    {
      print(paste("problem status 1 and couple 1 Female",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
    if(data$CFE0[i]>data$FE0[i]) # Si le nombre de femelle en couple à l'état E0 est supérieur au nombre de femelle à l'état E0
    {
      print(paste("problem status 0 and couple 0 Female",i+1,"Hour :",data$H[i],"Media :",data$Milieu[i],"Drug :",data$Drug[i],"Dose :",data$Dose[i],"Well :",data$Puit[i])) # Affiche la ligne en défaut ainsi que leur contenue
      Control<-1 # Indique que la ligne testée possède une erreur
    }
  }
  data$Error<-Control # Enregistre qu'il y a eu au moins une erreur de notation dans la colonne Error
  print(paste("error detect :",Control)) # Affiche qu'il y a eu au moins une erreur de notation des schistosomes
  return(data) # Renvoi le contenu de data en sortie de la fonction
}

# programme main -----------------------------------


Ndose<-N_of_variable(data_real$Dose) # Récupère le nombre de dose présent de le jeu de donnée
Nmedia<-N_of_variable(data_real$Milieu) # Récupère le nombre de Milieu présent de le jeu de donnée
Nwell<-N_of_variable(data_real$Puit) # Récupère le nombre de Puit présent de le jeu de donnée
Ntime<-N_of_variable(data_real$H) # Récupère le nombre de d'Heures présent de le jeu de donnée

Drug<-levels(data_real$Drug) # Definie la colonne Drug comme de type levels (caractère)
Unit<-levels(data_real$Unit) # Definie la colonne Unit comme de type levels (caractère)
data<-data_real # Enregistre le contenu de data_real dans data
data_save_1<-data # Enregistre le contenu de data dans data_save_1

data<-error_detection(data) # Exécute la fonction error_detection pour vérifier si il y a eu des erreurs lors de la prise des données 

compteur<-0 # Initialise une variable compteur à zéro

data_error_correct_E1_E0<-NA # Crée une variable data_error_correct_E1_E0 de valeur NA

for (i in 1:length(Nmedia)) # Pour chaque milieu de culture
{
  Type_media<-Nmedia[i] # Enregistre dans Type_media la valeur du milieu en cours de vérification
  assign(paste("media",sep=""),data[data$Milieu %in% c(Nmedia[i]),]) # Récupérer dans "media" toutes les lignes de données avec le milieu en cours de vérification (exemple : Toutes les données pour le M199)
  
  Xdose<-N_of_variable(media$Dose) # Récupère le nombre de dose présent de le jeu de donnée
  
  for(j in 1:length(Xdose)) # Pour chaque dose
  {
    assign(paste("media_dose",sep=""),media[media$Dose %in% c(Xdose[j]),]) # Récupérer dans "media_dose" toutes les lignes de données avec la dose en cours de vérification (exemple : Toutes les données pour le M199 dose 900 nM de PZQ)
    Xwell<-N_of_variable(media$Puit) # Récupère le nombre de puit présent de le jeu de donnée
    for(k in 1:length(Xwell)) #Pour chaque puit
    {
      assign(paste("media_dose_well",sep=""),media_dose[media_dose$Puit %in% c(Xwell[k]),]) # Récupérer dans "media_dose_well" toutes les lignes de données avec le puit en cours de vérification (exemple : Toutes les données pour le M199 dose 900 nM de PZQ puit N°1)
      compteur<-compteur+1 # Ajoute plus 1 au compteur
      Xtime<-N_of_variable(media$H) # Récupère le nombre de Milieu d'heure présent de le jeu de donnée
      Time_to_controle<-rev(Xtime) # Range dans l'ordre décroissant les heures contenues dans le jeu de données media_dose_well 
      for(t in 1:length(Time_to_controle)) # Pour chaque heure de prise de données
      {
        for(ligne in 1: length(media_dose_well$DDP)) # Pour chaque ligne dans le jeu de données
        {
          if(media_dose_well$H[ligne]==Time_to_controle[t]) # Si l'heure présente dans la ligne contrôlée est égale au temps (Time_to-controle) à vérifier
          {
            E0_to_correct<-media_dose_well$E0[ligne] # Enregistre la valeur contenue dans la ligne E0 dans E0_to_correct
            E1_to_correct<-media_dose_well$E1[ligne] # Enregistre la valeur contenue dans la ligne E1 dans E1_to_correct
            for(ligne_to_controle in 1:length(media_dose_well$DDP)) # Pour chaque ligne du jeu de donnée 
            {
              if(t<length(Time_to_controle)) # Si t est inférieur à Time_to_controle
              {
                if(media_dose_well$H[ligne_to_controle]<media_dose_well$H[ligne]&&Time_to_controle[t+1]==media_dose_well$H[ligne_to_controle]) # Si l'heure à la ligne "ligne_to_controle" est inférieur à l'heure contenue dans la ligne "ligne" et que l'heure contenue dans "Time_to_controle[t+1]" est égale à l'heure présente dans la ligne "ligne_to_controle" (si media_dose_well$H[ligne_to_controle]= J5 ; media_dose_well$H[ligne] = J7 et Time_to_controle[t+1]= J7 ; media_dose_well$H[ligne_to_controle] = J7)  
                {
                  if(media_dose_well$E0[ligne_to_controle]>media_dose_well$E0[ligne]) # Si le nombre de mort à J5 est supérieur au nombre de mort à J7 
                  {
                    print("") # Affiche qu'il y a eu une erreur dans le nombre de mort lors de la prise de donnée et affiche les lignes conservées
                    print(paste("Time_before",media_dose_well$H[ligne_to_controle],"Time_after",media_dose_well$H[ligne]))
                    print(paste("E0 Value_before",media_dose_well$E0[ligne_to_controle],"E0 Value_after",media_dose_well$E0[ligne]))
                    print(paste("E1 Value_before",media_dose_well$E1[ligne_to_controle],"E1 Value_after",media_dose_well$E1[ligne]))
                    
                    # Applique la correction suivante :
                    media_dose_well$E1[ligne_to_controle]<-media_dose_well$E1[ligne_to_controle]+media_dose_well$E0[ligne_to_controle]-E0_to_correct  
                    media_dose_well$E0[ligne_to_controle]<-E0_to_correct
                    
                    # Affiche la correction appliquée 
                    print(paste("Time_before",media_dose_well$H[ligne_to_controle],"Time_after",media_dose_well$H[ligne]))
                    print(paste("E0 Value_before",media_dose_well$E0[ligne_to_controle],"E0 Value_after",media_dose_well$E0[ligne]))
                    print(paste("E1 Value_before",media_dose_well$E1[ligne_to_controle],"E1 Value_after",media_dose_well$E1[ligne]))
                    
                  }
                  
                }
              }
            }
          }
        }
      }
      
      data_error_correct_E1_E0<-na.rm(rbind(data_error_correct_E1_E0,media_dose_well)) # Enregistre les corrections dans data_error_correct_E1_E0 en supprimant les lignes avec des NA 
      
    }
  }
}

data<-data_error_correct_E1_E0 # Enregistre dans data_error_correct_E1_E0 dans data

data$Survie<-NA # Ajoute une colonne Survie 
data$Dettached<-NA # Ajoute une colonne Dettached
data$Inactive_Drug<-NA # Ajoute une colonne Inactive_Drug
data$Moderately_Active_Drug<-NA # Ajoute une colonne Moderately_Active_Drug 
data$Active_Drug<-NA # Ajoute une colonne Active_Drug

data$dSurvie<-NA # Ajoute une colonne dSurvie pour le calcul des écarts type
data$dDettached<-NA # Ajoute une colonne dDettached pour le calcul des écarts type
data$dActive_Drug<-NA # Ajoute une colonne dActive_Drug pour le calcul des écarts type
data$dModerately_Active_Drug<-NA # Ajoute une colonne dModerately_Active_Drug pour le calcul des écarts type
data$dActive_Drug<-NA # Ajoute une colonne dActive_Drug pour le calcul des écarts type

data$Survie<-(data$E4+data$E3+data$E2)*100/data$NST # Calcul le pourcentage de vers vivants dans chaque puit (E4+E3+E2)  
data$Dettached<-(data$E3+data$E2+data$E1+data$E0)*100/data$NST # Calcul le pourcentage de vers détachés de la surface de chaque puit (E1+E0) 

data$Inactive_Drug<-(data$E4+data$E3)*100/data$NST # Calcul le pourcentage le pourcentage de vers dont la drogue était inactive (E4+E3)
data$Moderately_Active_Drug<-(data$E3)*100/data$NST # Calcul le pourcentage le pourcentage de vers dont la drogue était modérément active (E4+E3)
data$Active_Drug<-(data$E1+data$E0)*100/data$NST # Calcul le pourcentage le pourcentage de vers dont la drogue était active (E4+E3)

data_save_1<-data # Enregistre le contenu de data dans data_save_1 

data_treat<-data # Enregistre le contenu de data dans data_treat
data_treat$NST<-NA # Remplace les valeurs contenues dans la colonne NST par des NA
data_treat$dNST<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dNSTM<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dNSTM<-NA # Création une colonne avec comme valeurs par défaut NA

data_treat$dCME4<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCME3<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCME2<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCME1<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCME0<-NA # Création une colonne avec comme valeurs par défaut NA

data_treat$dME4<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dME3<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dME2<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dME1<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dME0<-NA # Création une colonne avec comme valeurs par défaut NA

data_treat$dCFE4<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCFE3<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCFE2<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCFE1<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCFE0<-NA # Création une colonne avec comme valeurs par défaut NA

data_treat$dFE4<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dFE3<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dFE2<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dFE1<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dFE0<-NA # Création une colonne avec comme valeurs par défaut NA

data_treat$dCouple<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCouple_Female<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dCouple_Male<-NA # Création une colonne avec comme valeurs par défaut NA

data_treat$dE4<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dE3<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dE2<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dE1<-NA # Création une colonne avec comme valeurs par défaut NA
data_treat$dE0<-NA # Création une colonne avec comme valeurs par défaut NA

data_save_2<-data_treat # Enregistre le contenue de data_treat dans data_save_2

# exportation sur Excel des données survie, dettached, active drug, moderately active et inactive drug  -----------------------------------

data_treat$Time<-data_treat$H # Crée la colonne Time puis enregistre dans celle-ci le contenu de la colonne H 
data_RPMI_treat<-data_treat[data_treat$Milieu %in% c(6),] # récupère uniquement les données qui comprennent le milieu RPMI Horse
data_M199_treat<-data_treat[data_treat$Milieu %in% c(1),] # récupère uniquement les données qui comprennent le milieu M199 Calf

Drug<-N_of_variable(data_save_1$Drug) # récupère le noms des drogues présents dans le jeux de données


if(Drug[1]=="ART") # Si le jeu de données a été fait avec de l'Artemisinine, exporte les données traitées dans des fichiers Excel 
{
  data_RPMI_0<-data_RPMI_treat[data_RPMI_treat$Time %in% c(0),] # Récupère uniquement les données à T0
  data_RPMI_1<-data_RPMI_treat[data_RPMI_treat$Time %in% c(1),] # Récupère uniquement les données pour le temps 1 heure
  data_RPMI_12<-data_RPMI_treat[data_RPMI_treat$Time %in% c(12),] # Récupère uniquement les données pour le temps 12 heures
  data_RPMI_24<-data_RPMI_treat[data_RPMI_treat$Time %in% c(24),] # Récupère uniquement les données pour le temps 24 heures
  data_RPMI_48<-data_RPMI_treat[data_RPMI_treat$Time %in% c(48),] # Récupère uniquement les données pour le temps 48 heures
  data_RPMI_72<-data_RPMI_treat[data_RPMI_treat$Time %in% c(72),] # Récupère uniquement les données pour le temps 72 heures
  data_RPMI_96<-data_RPMI_treat[data_RPMI_treat$Time %in% c(96),] # Récupère uniquement les données pour le temps 96 heures
  data_RPMI_120<-data_RPMI_treat[data_RPMI_treat$Time %in% c(120),] # Récupère uniquement les données pour le temps 120 heures
  
  data_RPMI_0_Survie<-subset(data_RPMI_0,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants à T0
  data_RPMI_1_Survie<-subset(data_RPMI_1,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 1 heure
  data_RPMI_12_Survie<-subset(data_RPMI_12,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 12 heures
  data_RPMI_24_Survie<-subset(data_RPMI_24,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 24 heures
  data_RPMI_48_Survie<-subset(data_RPMI_48,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 48 heures
  data_RPMI_72_Survie<-subset(data_RPMI_72,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 72 heures
  data_RPMI_96_Survie<-subset(data_RPMI_96,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 96 heures
  data_RPMI_120_Survie<-subset(data_RPMI_120,select = c(7,9,43)) # Conserve que les colonnes doses, puits et le pourcentage de vers vivants pour le temps 120 heures
  
  data_RPMI_0_Dettached<-subset(data_RPMI_0,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit à T0
  data_RPMI_1_Dettached<-subset(data_RPMI_1,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 1 heure
  data_RPMI_12_Dettached<-subset(data_RPMI_12,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 12 heures
  data_RPMI_24_Dettached<-subset(data_RPMI_24,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 24 heures
  data_RPMI_48_Dettached<-subset(data_RPMI_48,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 48 heures
  data_RPMI_72_Dettached<-subset(data_RPMI_72,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 72 heures
  data_RPMI_96_Dettached<-subset(data_RPMI_96,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 96 heures
  data_RPMI_120_Dettached<-subset(data_RPMI_120,select = c(7,9,44)) # Conserve que les colonnes doses, puits et le pourcentage de vers détachés de la surface du puit pour le temps 120 heures
  
  data_RPMI_0_Inactive_Drug<-subset(data_RPMI_0,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) à T0
  data_RPMI_1_Inactive_Drug<-subset(data_RPMI_1,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 1 heure
  data_RPMI_12_Inactive_Drug<-subset(data_RPMI_12,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 12 heures
  data_RPMI_24_Inactive_Drug<-subset(data_RPMI_24,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 24 heures
  data_RPMI_48_Inactive_Drug<-subset(data_RPMI_48,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 48 heures
  data_RPMI_72_Inactive_Drug<-subset(data_RPMI_72,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 72 heures
  data_RPMI_96_Inactive_Drug<-subset(data_RPMI_96,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 96 heures
  data_RPMI_120_Inactive_Drug<-subset(data_RPMI_120,select = c(7,9,45)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été inactive (E4+E3) pour le temps 120 heures
  
  data_RPMI_0_Moderately_Active_Drug<-subset(data_RPMI_0,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) à T0
  data_RPMI_1_Moderately_Active_Drug<-subset(data_RPMI_1,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 1 heure
  data_RPMI_12_Moderately_Active_Drug<-subset(data_RPMI_12,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 12 heures
  data_RPMI_24_Moderately_Active_Drug<-subset(data_RPMI_24,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 24 heures
  data_RPMI_48_Moderately_Active_Drug<-subset(data_RPMI_48,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 48 heures
  data_RPMI_72_Moderately_Active_Drug<-subset(data_RPMI_72,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 72 heures
  data_RPMI_96_Moderately_Active_Drug<-subset(data_RPMI_96,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 96 heures
  data_RPMI_120_Moderately_Active_Drug<-subset(data_RPMI_120,select = c(7,9,46)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été modérément active (E2) pour le temps 120 heures
  
  data_RPMI_0_Active_Drug<-subset(data_RPMI_0,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) à T0
  data_RPMI_1_Active_Drug<-subset(data_RPMI_1,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 1 heure
  data_RPMI_12_Active_Drug<-subset(data_RPMI_12,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 12 heures
  data_RPMI_24_Active_Drug<-subset(data_RPMI_24,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 24 heures
  data_RPMI_48_Active_Drug<-subset(data_RPMI_48,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 48 heures
  data_RPMI_72_Active_Drug<-subset(data_RPMI_72,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 72 heures
  data_RPMI_96_Active_Drug<-subset(data_RPMI_96,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 96 heures
  data_RPMI_120_Active_Drug<-subset(data_RPMI_120,select = c(7,9,47)) # Conserve que les colonnes doses, puits et le pourcentage de vers où la drogue a été active (E4+E3) pour le temps 120 heures
  
  #Rearrangement des colonnes
  
  data_RPMI_0_Survie<-data_RPMI_0_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_1_Survie<-data_RPMI_1_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_12_Survie<-data_RPMI_12_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_24_Survie<-data_RPMI_24_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_48_Survie<-data_RPMI_48_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_72_Survie<-data_RPMI_72_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_96_Survie<-data_RPMI_96_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_120_Survie<-data_RPMI_120_Survie%>% spread(key= Puit, value= Survie)
  
  
  data_RPMI_0_Dettached<-data_RPMI_0_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_1_Dettached<-data_RPMI_1_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_12_Dettached<-data_RPMI_12_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_24_Dettached<-data_RPMI_24_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_48_Dettached<-data_RPMI_48_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_72_Dettached<-data_RPMI_72_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_96_Dettached<-data_RPMI_96_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_120_Dettached<-data_RPMI_120_Dettached%>% spread(key= Puit, value= Dettached)
  
  data_RPMI_0_Inactive_Drug<-data_RPMI_0_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_1_Inactive_Drug<-data_RPMI_1_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_12_Inactive_Drug<-data_RPMI_12_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_24_Inactive_Drug<-data_RPMI_24_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_48_Inactive_Drug<-data_RPMI_48_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_72_Inactive_Drug<-data_RPMI_72_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_96_Inactive_Drug<-data_RPMI_96_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_120_Inactive_Drug<-data_RPMI_120_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  
  data_RPMI_0_Moderately_Active_Drug<-data_RPMI_0_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_1_Moderately_Active_Drug<-data_RPMI_1_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_12_Moderately_Active_Drug<-data_RPMI_12_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_24_Moderately_Active_Drug<-data_RPMI_24_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_48_Moderately_Active_Drug<-data_RPMI_48_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_72_Moderately_Active_Drug<-data_RPMI_72_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_96_Moderately_Active_Drug<-data_RPMI_96_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_120_Moderately_Active_Drug<-data_RPMI_120_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  
  data_RPMI_0_Active_Drug<-data_RPMI_0_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_1_Active_Drug<-data_RPMI_1_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_12_Active_Drug<-data_RPMI_12_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_24_Active_Drug<-data_RPMI_24_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_48_Active_Drug<-data_RPMI_48_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_72_Active_Drug<-data_RPMI_72_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_96_Active_Drug<-data_RPMI_96_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_120_Active_Drug<-data_RPMI_120_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  
  # Sauvegarde dans ces fichiers Excel
  
  # Survie
  
  write_xlsx(data_RPMI_0_Survie,path = "ART_µM_data_RPMI_0_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Survie,path = "ART_µM_data_RPMI_1_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Survie,path = "ART_µM_data_RPMI_12_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Survie,path = "ART_µM_data_RPMI_24_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Survie,path = "ART_µM_data_RPMI_48_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Survie,path = "ART_µM_data_RPMI_72_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Survie,path = "ART_µM_data_RPMI_96_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Survie,path = "ART_µM_data_RPMI_120_Survie_heures.xlsx",col_names = TRUE)
  
  #Dettached
  
  write_xlsx(data_RPMI_0_Dettached,path = ,"ART_µM_data_RPMI_0_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Dettached,path = "ART_µM_data_RPMI_1_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Dettached,path = "ART_µM_data_RPMI_12_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Dettached,path = "ART_µM_data_RPMI_24_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Dettached,path = "ART_µM_data_RPMI_48_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Dettached,path = "ART_µM_data_RPMI_72_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Dettached,path = "ART_µM_data_RPMI_96_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Dettached,path = "ART_µM_data_RPMI_120_Dettached_heures.xlsx",col_names = TRUE)
  
  # Inactive_Drug
  
  write_xlsx(data_RPMI_0_Inactive_Drug,path = "ART_µM_data_RPMI_0_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Inactive_Drug,path = "ART_µM_data_RPMI_1_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Inactive_Drug,path = "ART_µM_data_RPMI_12_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Inactive_Drug,path = "ART_µM_data_RPMI_24_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Inactive_Drug,path = "ART_µM_data_RPMI_48_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Inactive_Drug,path = "ART_µM_data_RPMI_72_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Inactive_Drug,path = "ART_µM_data_RPMI_96_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Inactive_Drug,path = "ART_µM_data_RPMI_120_Inactive_Drug_heures.xlsx",col_names = TRUE)
  
  #Moderately_Active_Drug
  
  write_xlsx(data_RPMI_0_Moderately_Active_Drug,path = "ART_µM_data_RPMI_0_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Moderately_Active_Drug,path = "ART_µM_data_RPMI_1_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Moderately_Active_Drug,path = "ART_µM_data_RPMI_12_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Moderately_Active_Drug,path = "ART_µM_data_RPMI_24_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Moderately_Active_Drug,path = "ART_µM_data_RPMI_48_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Moderately_Active_Drug,path = "ART_µM_data_RPMI_72_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Moderately_Active_Drug,path = "ART_µM_data_RPMI_96_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Moderately_Active_Drug,path = "ART_µM_data_RPMI_120_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  
  #Active_Drug
  
  write_xlsx(data_RPMI_0_Active_Drug,path = "ART_µM_data_RPMI_0_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Active_Drug,path = "ART_µM_data_RPMI_1_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Active_Drug,path = "ART_µM_data_RPMI_12_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Active_Drug,path = "ART_µM_data_RPMI_24_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Active_Drug,path = "ART_µM_data_RPMI_48_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Active_Drug,path = "ART_µM_data_RPMI_72_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Active_Drug,path = "ART_µM_data_RPMI_96_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Active_Drug,path = "ART_µM_data_RPMI_120_Active_Drug_heures.xlsx",col_names = TRUE)
  
  
  #milieu M199 --------------------------------------------------
  # Même chose pour le milieu M199 Calf
  
  data_M199_0<-data_M199_treat[data_M199_treat$Time %in% c(0),]
  data_M199_1<-data_M199_treat[data_M199_treat$Time %in% c(1),]
  data_M199_12<-data_M199_treat[data_M199_treat$Time %in% c(12),]
  data_M199_24<-data_M199_treat[data_M199_treat$Time %in% c(24),]
  data_M199_48<-data_M199_treat[data_M199_treat$Time %in% c(48),]
  data_M199_72<-data_M199_treat[data_M199_treat$Time %in% c(72),]
  data_M199_96<-data_M199_treat[data_M199_treat$Time %in% c(96),]
  data_M199_120<-data_M199_treat[data_M199_treat$Time %in% c(120),]
  
  data_M199_0_Survie<-subset(data_M199_0,select = c(7,9,43))
  data_M199_1_Survie<-subset(data_M199_1,select = c(7,9,43))
  data_M199_12_Survie<-subset(data_M199_12,select = c(7,9,43))
  data_M199_24_Survie<-subset(data_M199_24,select = c(7,9,43))
  data_M199_48_Survie<-subset(data_M199_48,select = c(7,9,43))
  data_M199_72_Survie<-subset(data_M199_72,select = c(7,9,43))
  data_M199_96_Survie<-subset(data_M199_96,select = c(7,9,43))
  data_M199_120_Survie<-subset(data_M199_120,select = c(7,9,43))
  
  data_M199_0_Dettached<-subset(data_M199_0,select = c(7,9,44))
  data_M199_1_Dettached<-subset(data_M199_1,select = c(7,9,44))
  data_M199_12_Dettached<-subset(data_M199_12,select = c(7,9,44))
  data_M199_24_Dettached<-subset(data_M199_24,select = c(7,9,44))
  data_M199_48_Dettached<-subset(data_M199_48,select = c(7,9,44))
  data_M199_72_Dettached<-subset(data_M199_72,select = c(7,9,44))
  data_M199_96_Dettached<-subset(data_M199_96,select = c(7,9,44))
  data_M199_120_Dettached<-subset(data_M199_120,select = c(7,9,44))
  
  data_M199_0_Inactive_Drug<-subset(data_M199_0,select = c(7,9,45))
  data_M199_1_Inactive_Drug<-subset(data_M199_1,select = c(7,9,45))
  data_M199_12_Inactive_Drug<-subset(data_M199_12,select = c(7,9,45))
  data_M199_24_Inactive_Drug<-subset(data_M199_24,select = c(7,9,45))
  data_M199_48_Inactive_Drug<-subset(data_M199_48,select = c(7,9,45))
  data_M199_72_Inactive_Drug<-subset(data_M199_72,select = c(7,9,45))
  data_M199_96_Inactive_Drug<-subset(data_M199_96,select = c(7,9,45))
  data_M199_120_Inactive_Drug<-subset(data_M199_120,select = c(7,9,45))
  
  data_M199_0_Moderately_Active_Drug<-subset(data_M199_0,select = c(7,9,46))
  data_M199_1_Moderately_Active_Drug<-subset(data_M199_1,select = c(7,9,46))
  data_M199_12_Moderately_Active_Drug<-subset(data_M199_12,select = c(7,9,46))
  data_M199_24_Moderately_Active_Drug<-subset(data_M199_24,select = c(7,9,46))
  data_M199_48_Moderately_Active_Drug<-subset(data_M199_48,select = c(7,9,46))
  data_M199_72_Moderately_Active_Drug<-subset(data_M199_72,select = c(7,9,46))
  data_M199_96_Moderately_Active_Drug<-subset(data_M199_96,select = c(7,9,46))
  data_M199_120_Moderately_Active_Drug<-subset(data_M199_120,select = c(7,9,46))
  
  data_M199_0_Active_Drug<-subset(data_M199_0,select = c(7,9,47))
  data_M199_1_Active_Drug<-subset(data_M199_1,select = c(7,9,47))
  data_M199_12_Active_Drug<-subset(data_M199_12,select = c(7,9,47))
  data_M199_24_Active_Drug<-subset(data_M199_24,select = c(7,9,47))
  data_M199_48_Active_Drug<-subset(data_M199_48,select = c(7,9,47))
  data_M199_72_Active_Drug<-subset(data_M199_72,select = c(7,9,47))
  data_M199_96_Active_Drug<-subset(data_M199_96,select = c(7,9,47))
  data_M199_120_Active_Drug<-subset(data_M199_120,select = c(7,9,47))
  
  # Rearrangement
  
  data_M199_0_Survie<-data_M199_0_Survie%>% spread(key= Puit, value= Survie)
  data_M199_1_Survie<-data_M199_1_Survie%>% spread(key= Puit, value= Survie)
  data_M199_12_Survie<-data_M199_12_Survie%>% spread(key= Puit, value= Survie)
  data_M199_24_Survie<-data_M199_24_Survie%>% spread(key= Puit, value= Survie)
  data_M199_48_Survie<-data_M199_48_Survie%>% spread(key= Puit, value= Survie)
  data_M199_72_Survie<-data_M199_72_Survie%>% spread(key= Puit, value= Survie)
  data_M199_96_Survie<-data_M199_96_Survie%>% spread(key= Puit, value= Survie)
  data_M199_120_Survie<-data_M199_120_Survie%>% spread(key= Puit, value= Survie)
  
  
  data_M199_0_Dettached<-data_M199_0_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_1_Dettached<-data_M199_1_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_12_Dettached<-data_M199_12_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_24_Dettached<-data_M199_24_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_48_Dettached<-data_M199_48_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_72_Dettached<-data_M199_72_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_96_Dettached<-data_M199_96_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_120_Dettached<-data_M199_120_Dettached%>% spread(key= Puit, value= Dettached)
  
  data_M199_0_Inactive_Drug<-data_M199_0_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_1_Inactive_Drug<-data_M199_1_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_12_Inactive_Drug<-data_M199_12_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_24_Inactive_Drug<-data_M199_24_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_48_Inactive_Drug<-data_M199_48_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_72_Inactive_Drug<-data_M199_72_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_96_Inactive_Drug<-data_M199_96_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_120_Inactive_Drug<-data_M199_120_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  
  data_M199_0_Moderately_Active_Drug<-data_M199_0_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_1_Moderately_Active_Drug<-data_M199_1_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_12_Moderately_Active_Drug<-data_M199_12_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_24_Moderately_Active_Drug<-data_M199_24_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_48_Moderately_Active_Drug<-data_M199_48_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_72_Moderately_Active_Drug<-data_M199_72_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_96_Moderately_Active_Drug<-data_M199_96_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_120_Moderately_Active_Drug<-data_M199_120_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  
  data_M199_0_Active_Drug<-data_M199_0_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_1_Active_Drug<-data_M199_1_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_12_Active_Drug<-data_M199_12_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_24_Active_Drug<-data_M199_24_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_48_Active_Drug<-data_M199_48_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_72_Active_Drug<-data_M199_72_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_96_Active_Drug<-data_M199_96_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_120_Active_Drug<-data_M199_120_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  
  # Sauvegarde en excel
  
  # Survie
  
  write_xlsx(data_M199_0_Survie,path = "ART_µM_data_M199_0_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Survie,path = "ART_µM_data_M199_1_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Survie,path = "ART_µM_data_M199_12_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Survie,path = "ART_µM_data_M199_24_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Survie,path = "ART_µM_data_M199_48_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Survie,path = "ART_µM_data_M199_72_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Survie,path = "ART_µM_data_M199_96_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Survie,path = "ART_µM_data_M199_120_Survie_heures.xlsx",col_names = TRUE)
  
  # Dettached
  
  write_xlsx(data_M199_0_Dettached,path = ,"ART_µM_data_M199_0_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Dettached,path = "ART_µM_data_M199_1_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Dettached,path = "ART_µM_data_M199_12_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Dettached,path = "ART_µM_data_M199_24_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Dettached,path = "ART_µM_data_M199_48_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Dettached,path = "ART_µM_data_M199_72_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Dettached,path = "ART_µM_data_M199_96_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Dettached,path = "ART_µM_data_M199_120_Dettached_heures.xlsx",col_names = TRUE)
  
  # Inactive_Drug
  
  write_xlsx(data_M199_0_Inactive_Drug,path = "ART_µM_data_M199_0_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Inactive_Drug,path = "ART_µM_data_M199_1_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Inactive_Drug,path = "ART_µM_data_M199_12_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Inactive_Drug,path = "ART_µM_data_M199_24_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Inactive_Drug,path = "ART_µM_data_M199_48_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Inactive_Drug,path = "ART_µM_data_M199_72_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Inactive_Drug,path = "ART_µM_data_M199_96_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Inactive_Drug,path = "ART_µM_data_M199_120_Inactive_Drug_heures.xlsx",col_names = TRUE)
  
  # Moderately_Active_Drug
  
  write_xlsx(data_M199_0_Moderately_Active_Drug,path = "ART_µM_data_M199_0_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Moderately_Active_Drug,path = "ART_µM_data_M199_1_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Moderately_Active_Drug,path = "ART_µM_data_M199_12_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Moderately_Active_Drug,path = "ART_µM_data_M199_24_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Moderately_Active_Drug,path = "ART_µM_data_M199_48_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Moderately_Active_Drug,path = "ART_µM_data_M199_72_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Moderately_Active_Drug,path = "ART_µM_data_M199_96_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Moderately_Active_Drug,path = "ART_µM_data_M199_120_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  
  # Active_Drug
  
  write_xlsx(data_M199_0_Active_Drug,path = "ART_µM_data_M199_0_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Active_Drug,path = "ART_µM_data_M199_1_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Active_Drug,path = "ART_µM_data_M199_12_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Active_Drug,path = "ART_µM_data_M199_24_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Active_Drug,path = "ART_µM_data_M199_48_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Active_Drug,path = "ART_µM_data_M199_72_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Active_Drug,path = "ART_µM_data_M199_96_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Active_Drug,path = "ART_µM_data_M199_120_Active_Drug_heures.xlsx",col_names = TRUE)
  
}

if(Drug[1]=="PZQ") # Si le jeu de données a été fait avec de l'Artemicinine, exporte les données traitées dans des fichiers excels de la même manière que pour l'ART 
{
  data_RPMI_0<-data_RPMI_treat[data_RPMI_treat$Time %in% c(0),]
  data_RPMI_1<-data_RPMI_treat[data_RPMI_treat$Time %in% c(1),]
  data_RPMI_12<-data_RPMI_treat[data_RPMI_treat$Time %in% c(12),]
  data_RPMI_24<-data_RPMI_treat[data_RPMI_treat$Time %in% c(24),]
  data_RPMI_48<-data_RPMI_treat[data_RPMI_treat$Time %in% c(48),]
  data_RPMI_72<-data_RPMI_treat[data_RPMI_treat$Time %in% c(72),]
  data_RPMI_96<-data_RPMI_treat[data_RPMI_treat$Time %in% c(96),]
  data_RPMI_120<-data_RPMI_treat[data_RPMI_treat$Time %in% c(120),]
  
  data_RPMI_0_Survie<-subset(data_RPMI_0,select = c(7,9,43))
  data_RPMI_1_Survie<-subset(data_RPMI_1,select = c(7,9,43))
  data_RPMI_12_Survie<-subset(data_RPMI_12,select = c(7,9,43))
  data_RPMI_24_Survie<-subset(data_RPMI_24,select = c(7,9,43))
  data_RPMI_48_Survie<-subset(data_RPMI_48,select = c(7,9,43))
  data_RPMI_72_Survie<-subset(data_RPMI_72,select = c(7,9,43))
  data_RPMI_96_Survie<-subset(data_RPMI_96,select = c(7,9,43))
  data_RPMI_120_Survie<-subset(data_RPMI_120,select = c(7,9,43))
  
  data_RPMI_0_Dettached<-subset(data_RPMI_0,select = c(7,9,44))
  data_RPMI_1_Dettached<-subset(data_RPMI_1,select = c(7,9,44))
  data_RPMI_12_Dettached<-subset(data_RPMI_12,select = c(7,9,44))
  data_RPMI_24_Dettached<-subset(data_RPMI_24,select = c(7,9,44))
  data_RPMI_48_Dettached<-subset(data_RPMI_48,select = c(7,9,44))
  data_RPMI_72_Dettached<-subset(data_RPMI_72,select = c(7,9,44))
  data_RPMI_96_Dettached<-subset(data_RPMI_96,select = c(7,9,44))
  data_RPMI_120_Dettached<-subset(data_RPMI_120,select = c(7,9,44))
  
  data_RPMI_0_Inactive_Drug<-subset(data_RPMI_0,select = c(7,9,45))
  data_RPMI_1_Inactive_Drug<-subset(data_RPMI_1,select = c(7,9,45))
  data_RPMI_12_Inactive_Drug<-subset(data_RPMI_12,select = c(7,9,45))
  data_RPMI_24_Inactive_Drug<-subset(data_RPMI_24,select = c(7,9,45))
  data_RPMI_48_Inactive_Drug<-subset(data_RPMI_48,select = c(7,9,45))
  data_RPMI_72_Inactive_Drug<-subset(data_RPMI_72,select = c(7,9,45))
  data_RPMI_96_Inactive_Drug<-subset(data_RPMI_96,select = c(7,9,45))
  data_RPMI_120_Inactive_Drug<-subset(data_RPMI_120,select = c(7,9,45))
  
  data_RPMI_0_Moderately_Active_Drug<-subset(data_RPMI_0,select = c(7,9,46))
  data_RPMI_1_Moderately_Active_Drug<-subset(data_RPMI_1,select = c(7,9,46))
  data_RPMI_12_Moderately_Active_Drug<-subset(data_RPMI_12,select = c(7,9,46))
  data_RPMI_24_Moderately_Active_Drug<-subset(data_RPMI_24,select = c(7,9,46))
  data_RPMI_48_Moderately_Active_Drug<-subset(data_RPMI_48,select = c(7,9,46))
  data_RPMI_72_Moderately_Active_Drug<-subset(data_RPMI_72,select = c(7,9,46))
  data_RPMI_96_Moderately_Active_Drug<-subset(data_RPMI_96,select = c(7,9,46))
  data_RPMI_120_Moderately_Active_Drug<-subset(data_RPMI_120,select = c(7,9,46))
  
  data_RPMI_0_Active_Drug<-subset(data_RPMI_0,select = c(7,9,47))
  data_RPMI_1_Active_Drug<-subset(data_RPMI_1,select = c(7,9,47))
  data_RPMI_12_Active_Drug<-subset(data_RPMI_12,select = c(7,9,47))
  data_RPMI_24_Active_Drug<-subset(data_RPMI_24,select = c(7,9,47))
  data_RPMI_48_Active_Drug<-subset(data_RPMI_48,select = c(7,9,47))
  data_RPMI_72_Active_Drug<-subset(data_RPMI_72,select = c(7,9,47))
  data_RPMI_96_Active_Drug<-subset(data_RPMI_96,select = c(7,9,47))
  data_RPMI_120_Active_Drug<-subset(data_RPMI_120,select = c(7,9,47))
  
  # Rearrangement
  
  data_RPMI_0_Survie<-data_RPMI_0_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_1_Survie<-data_RPMI_1_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_12_Survie<-data_RPMI_12_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_24_Survie<-data_RPMI_24_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_48_Survie<-data_RPMI_48_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_72_Survie<-data_RPMI_72_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_96_Survie<-data_RPMI_96_Survie%>% spread(key= Puit, value= Survie)
  data_RPMI_120_Survie<-data_RPMI_120_Survie%>% spread(key= Puit, value= Survie)
  
  
  data_RPMI_0_Dettached<-data_RPMI_0_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_1_Dettached<-data_RPMI_1_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_12_Dettached<-data_RPMI_12_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_24_Dettached<-data_RPMI_24_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_48_Dettached<-data_RPMI_48_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_72_Dettached<-data_RPMI_72_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_96_Dettached<-data_RPMI_96_Dettached%>% spread(key= Puit, value= Dettached)
  data_RPMI_120_Dettached<-data_RPMI_120_Dettached%>% spread(key= Puit, value= Dettached)
  
  data_RPMI_0_Inactive_Drug<-data_RPMI_0_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_1_Inactive_Drug<-data_RPMI_1_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_12_Inactive_Drug<-data_RPMI_12_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_24_Inactive_Drug<-data_RPMI_24_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_48_Inactive_Drug<-data_RPMI_48_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_72_Inactive_Drug<-data_RPMI_72_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_96_Inactive_Drug<-data_RPMI_96_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_RPMI_120_Inactive_Drug<-data_RPMI_120_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  
  data_RPMI_0_Moderately_Active_Drug<-data_RPMI_0_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_1_Moderately_Active_Drug<-data_RPMI_1_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_12_Moderately_Active_Drug<-data_RPMI_12_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_24_Moderately_Active_Drug<-data_RPMI_24_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_48_Moderately_Active_Drug<-data_RPMI_48_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_72_Moderately_Active_Drug<-data_RPMI_72_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_96_Moderately_Active_Drug<-data_RPMI_96_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_RPMI_120_Moderately_Active_Drug<-data_RPMI_120_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  
  data_RPMI_0_Active_Drug<-data_RPMI_0_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_1_Active_Drug<-data_RPMI_1_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_12_Active_Drug<-data_RPMI_12_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_24_Active_Drug<-data_RPMI_24_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_48_Active_Drug<-data_RPMI_48_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_72_Active_Drug<-data_RPMI_72_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_96_Active_Drug<-data_RPMI_96_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_RPMI_120_Active_Drug<-data_RPMI_120_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  
  # Sauvegarde en excel
  
  # Survie
  
  write_xlsx(data_RPMI_0_Survie,path = "PZQ_µM_data_RPMI_0_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Survie,path = "PZQ_µM_data_RPMI_1_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Survie,path = "PZQ_µM_data_RPMI_12_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Survie,path = "PZQ_µM_data_RPMI_24_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Survie,path = "PZQ_µM_data_RPMI_48_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Survie,path = "PZQ_µM_data_RPMI_72_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Survie,path = "PZQ_µM_data_RPMI_96_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Survie,path = "PZQ_µM_data_RPMI_120_Survie_heures.xlsx",col_names = TRUE)
  
  # Dettached
  
  write_xlsx(data_RPMI_0_Dettached,path = ,"PZQ_µM_data_RPMI_0_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Dettached,path = "PZQ_µM_data_RPMI_1_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Dettached,path = "PZQ_µM_data_RPMI_12_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Dettached,path = "PZQ_µM_data_RPMI_24_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Dettached,path = "PZQ_µM_data_RPMI_48_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Dettached,path = "PZQ_µM_data_RPMI_72_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Dettached,path = "PZQ_µM_data_RPMI_96_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Dettached,path = "PZQ_µM_data_RPMI_120_Dettached_heures.xlsx",col_names = TRUE)
  
  # Inactive_Drug
  
  write_xlsx(data_RPMI_0_Inactive_Drug,path = "PZQ_µM_data_RPMI_0_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Inactive_Drug,path = "PZQ_µM_data_RPMI_1_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Inactive_Drug,path = "PZQ_µM_data_RPMI_12_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Inactive_Drug,path = "PZQ_µM_data_RPMI_24_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Inactive_Drug,path = "PZQ_µM_data_RPMI_48_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Inactive_Drug,path = "PZQ_µM_data_RPMI_72_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Inactive_Drug,path = "PZQ_µM_data_RPMI_96_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Inactive_Drug,path = "PZQ_µM_data_RPMI_120_Inactive_Drug_heures.xlsx",col_names = TRUE)
  
  # Moderately_Active_Drug
  
  write_xlsx(data_RPMI_0_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_0_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_1_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_12_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_24_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_48_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_72_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_96_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Moderately_Active_Drug,path = "PZQ_µM_data_RPMI_120_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  
  # Active_Drug
  
  write_xlsx(data_RPMI_0_Active_Drug,path = "PZQ_µM_data_RPMI_0_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_1_Active_Drug,path = "PZQ_µM_data_RPMI_1_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_12_Active_Drug,path = "PZQ_µM_data_RPMI_12_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_24_Active_Drug,path = "PZQ_µM_data_RPMI_24_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_48_Active_Drug,path = "PZQ_µM_data_RPMI_48_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_72_Active_Drug,path = "PZQ_µM_data_RPMI_72_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_96_Active_Drug,path = "PZQ_µM_data_RPMI_96_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_RPMI_120_Active_Drug,path = "PZQ_µM_data_RPMI_120_Active_Drug_heures.xlsx",col_names = TRUE)
  
  
  # Milieu M199 --------------------------------------------------
  
  data_M199_0<-data_M199_treat[data_M199_treat$Time %in% c(0),]
  data_M199_1<-data_M199_treat[data_M199_treat$Time %in% c(1),]
  data_M199_12<-data_M199_treat[data_M199_treat$Time %in% c(12),]
  data_M199_24<-data_M199_treat[data_M199_treat$Time %in% c(24),]
  data_M199_48<-data_M199_treat[data_M199_treat$Time %in% c(48),]
  data_M199_72<-data_M199_treat[data_M199_treat$Time %in% c(72),]
  data_M199_96<-data_M199_treat[data_M199_treat$Time %in% c(96),]
  data_M199_120<-data_M199_treat[data_M199_treat$Time %in% c(120),]
  
  data_M199_0_Survie<-subset(data_M199_0,select = c(7,9,43))
  data_M199_1_Survie<-subset(data_M199_1,select = c(7,9,43))
  data_M199_12_Survie<-subset(data_M199_12,select = c(7,9,43))
  data_M199_24_Survie<-subset(data_M199_24,select = c(7,9,43))
  data_M199_48_Survie<-subset(data_M199_48,select = c(7,9,43))
  data_M199_72_Survie<-subset(data_M199_72,select = c(7,9,43))
  data_M199_96_Survie<-subset(data_M199_96,select = c(7,9,43))
  data_M199_120_Survie<-subset(data_M199_120,select = c(7,9,43))
  
  data_M199_0_Dettached<-subset(data_M199_0,select = c(7,9,44))
  data_M199_1_Dettached<-subset(data_M199_1,select = c(7,9,44))
  data_M199_12_Dettached<-subset(data_M199_12,select = c(7,9,44))
  data_M199_24_Dettached<-subset(data_M199_24,select = c(7,9,44))
  data_M199_48_Dettached<-subset(data_M199_48,select = c(7,9,44))
  data_M199_72_Dettached<-subset(data_M199_72,select = c(7,9,44))
  data_M199_96_Dettached<-subset(data_M199_96,select = c(7,9,44))
  data_M199_120_Dettached<-subset(data_M199_120,select = c(7,9,44))
  
  data_M199_0_Inactive_Drug<-subset(data_M199_0,select = c(7,9,45))
  data_M199_1_Inactive_Drug<-subset(data_M199_1,select = c(7,9,45))
  data_M199_12_Inactive_Drug<-subset(data_M199_12,select = c(7,9,45))
  data_M199_24_Inactive_Drug<-subset(data_M199_24,select = c(7,9,45))
  data_M199_48_Inactive_Drug<-subset(data_M199_48,select = c(7,9,45))
  data_M199_72_Inactive_Drug<-subset(data_M199_72,select = c(7,9,45))
  data_M199_96_Inactive_Drug<-subset(data_M199_96,select = c(7,9,45))
  data_M199_120_Inactive_Drug<-subset(data_M199_120,select = c(7,9,45))
  
  data_M199_0_Moderately_Active_Drug<-subset(data_M199_0,select = c(7,9,46))
  data_M199_1_Moderately_Active_Drug<-subset(data_M199_1,select = c(7,9,46))
  data_M199_12_Moderately_Active_Drug<-subset(data_M199_12,select = c(7,9,46))
  data_M199_24_Moderately_Active_Drug<-subset(data_M199_24,select = c(7,9,46))
  data_M199_48_Moderately_Active_Drug<-subset(data_M199_48,select = c(7,9,46))
  data_M199_72_Moderately_Active_Drug<-subset(data_M199_72,select = c(7,9,46))
  data_M199_96_Moderately_Active_Drug<-subset(data_M199_96,select = c(7,9,46))
  data_M199_120_Moderately_Active_Drug<-subset(data_M199_120,select = c(7,9,46))
  
  data_M199_0_Active_Drug<-subset(data_M199_0,select = c(7,9,47))
  data_M199_1_Active_Drug<-subset(data_M199_1,select = c(7,9,47))
  data_M199_12_Active_Drug<-subset(data_M199_12,select = c(7,9,47))
  data_M199_24_Active_Drug<-subset(data_M199_24,select = c(7,9,47))
  data_M199_48_Active_Drug<-subset(data_M199_48,select = c(7,9,47))
  data_M199_72_Active_Drug<-subset(data_M199_72,select = c(7,9,47))
  data_M199_96_Active_Drug<-subset(data_M199_96,select = c(7,9,47))
  data_M199_120_Active_Drug<-subset(data_M199_120,select = c(7,9,47))
  
  # Rearrangement
  
  data_M199_0_Survie<-data_M199_0_Survie%>% spread(key= Puit, value= Survie)
  data_M199_1_Survie<-data_M199_1_Survie%>% spread(key= Puit, value= Survie)
  data_M199_12_Survie<-data_M199_12_Survie%>% spread(key= Puit, value= Survie)
  data_M199_24_Survie<-data_M199_24_Survie%>% spread(key= Puit, value= Survie)
  data_M199_48_Survie<-data_M199_48_Survie%>% spread(key= Puit, value= Survie)
  data_M199_72_Survie<-data_M199_72_Survie%>% spread(key= Puit, value= Survie)
  data_M199_96_Survie<-data_M199_96_Survie%>% spread(key= Puit, value= Survie)
  data_M199_120_Survie<-data_M199_120_Survie%>% spread(key= Puit, value= Survie)
  
  
  data_M199_0_Dettached<-data_M199_0_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_1_Dettached<-data_M199_1_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_12_Dettached<-data_M199_12_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_24_Dettached<-data_M199_24_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_48_Dettached<-data_M199_48_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_72_Dettached<-data_M199_72_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_96_Dettached<-data_M199_96_Dettached%>% spread(key= Puit, value= Dettached)
  data_M199_120_Dettached<-data_M199_120_Dettached%>% spread(key= Puit, value= Dettached)
  
  data_M199_0_Inactive_Drug<-data_M199_0_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_1_Inactive_Drug<-data_M199_1_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_12_Inactive_Drug<-data_M199_12_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_24_Inactive_Drug<-data_M199_24_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_48_Inactive_Drug<-data_M199_48_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_72_Inactive_Drug<-data_M199_72_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_96_Inactive_Drug<-data_M199_96_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  data_M199_120_Inactive_Drug<-data_M199_120_Inactive_Drug%>% spread(key= Puit, value= Inactive_Drug)
  
  data_M199_0_Moderately_Active_Drug<-data_M199_0_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_1_Moderately_Active_Drug<-data_M199_1_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_12_Moderately_Active_Drug<-data_M199_12_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_24_Moderately_Active_Drug<-data_M199_24_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_48_Moderately_Active_Drug<-data_M199_48_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_72_Moderately_Active_Drug<-data_M199_72_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_96_Moderately_Active_Drug<-data_M199_96_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  data_M199_120_Moderately_Active_Drug<-data_M199_120_Moderately_Active_Drug%>% spread(key= Puit, value= Moderately_Active_Drug)
  
  data_M199_0_Active_Drug<-data_M199_0_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_1_Active_Drug<-data_M199_1_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_12_Active_Drug<-data_M199_12_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_24_Active_Drug<-data_M199_24_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_48_Active_Drug<-data_M199_48_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_72_Active_Drug<-data_M199_72_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_96_Active_Drug<-data_M199_96_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  data_M199_120_Active_Drug<-data_M199_120_Active_Drug%>% spread(key= Puit, value= Active_Drug)
  
  # Sauvegarde en excel
  
  # Survie
  
  write_xlsx(data_M199_0_Survie,path = "PZQ_µM_data_M199_0_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Survie,path = "PZQ_µM_data_M199_1_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Survie,path = "PZQ_µM_data_M199_12_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Survie,path = "PZQ_µM_data_M199_24_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Survie,path = "PZQ_µM_data_M199_48_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Survie,path = "PZQ_µM_data_M199_72_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Survie,path = "PZQ_µM_data_M199_96_Survie_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Survie,path = "PZQ_µM_data_M199_120_Survie_heures.xlsx",col_names = TRUE)
  
  # Dettached
  
  write_xlsx(data_M199_0_Dettached,path = ,"PZQ_µM_data_M199_0_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Dettached,path = "PZQ_µM_data_M199_1_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Dettached,path = "PZQ_µM_data_M199_12_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Dettached,path = "PZQ_µM_data_M199_24_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Dettached,path = "PZQ_µM_data_M199_48_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Dettached,path = "PZQ_µM_data_M199_72_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Dettached,path = "PZQ_µM_data_M199_96_Dettached_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Dettached,path = "PZQ_µM_data_M199_120_Dettached_heures.xlsx",col_names = TRUE)
  
  # Inactive_Drug
  
  write_xlsx(data_M199_0_Inactive_Drug,path = "PZQ_µM_data_M199_0_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Inactive_Drug,path = "PZQ_µM_data_M199_1_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Inactive_Drug,path = "PZQ_µM_data_M199_12_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Inactive_Drug,path = "PZQ_µM_data_M199_24_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Inactive_Drug,path = "PZQ_µM_data_M199_48_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Inactive_Drug,path = "PZQ_µM_data_M199_72_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Inactive_Drug,path = "PZQ_µM_data_M199_96_Inactive_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Inactive_Drug,path = "PZQ_µM_data_M199_120_Inactive_Drug_heures.xlsx",col_names = TRUE)
  
  # Moderately_Active_Drug
  
  write_xlsx(data_M199_0_Moderately_Active_Drug,path = "PZQ_µM_data_M199_0_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Moderately_Active_Drug,path = "PZQ_µM_data_M199_1_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Moderately_Active_Drug,path = "PZQ_µM_data_M199_12_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Moderately_Active_Drug,path = "PZQ_µM_data_M199_24_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Moderately_Active_Drug,path = "PZQ_µM_data_M199_48_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Moderately_Active_Drug,path = "PZQ_µM_data_M199_72_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Moderately_Active_Drug,path = "PZQ_µM_data_M199_96_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Moderately_Active_Drug,path = "PZQ_µM_data_M199_120_Moderately_Active_Drug_heures.xlsx",col_names = TRUE)
  
  # Active_Drug
  
  write_xlsx(data_M199_0_Active_Drug,path = "PZQ_µM_data_M199_0_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_1_Active_Drug,path = "PZQ_µM_data_M199_1_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_12_Active_Drug,path = "PZQ_µM_data_M199_12_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_24_Active_Drug,path = "PZQ_µM_data_M199_24_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_48_Active_Drug,path = "PZQ_µM_data_M199_48_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_72_Active_Drug,path = "PZQ_µM_data_M199_72_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_96_Active_Drug,path = "PZQ_µM_data_M199_96_Active_Drug_heures.xlsx",col_names = TRUE)
  write_xlsx(data_M199_120_Active_Drug,path = "PZQ_µM_data_M199_120_Active_Drug_heures.xlsx",col_names = TRUE)
  
}


# Calcul des écarts type  -----------------------------------
data_treat<-data_save_2 # Recherche les données présente dans data_save_2 dans data_treat
data_treat$DDP<-NA

compteur<-0 # Initialise un compteur à zéro

for (i in 1:length(Nmedia)) # Pour chaque milieu
{
  assign(paste("media",sep=""),data[data$Milieu %in% c(Nmedia[i]),]) # Enregistre dans media les données avec le milieux de culture en cours de traitement (Exemple : Enregistres dans media les jeux de données pour le milieu M199 Calf)
  Xtime<-N_of_variable(media$H) # Récupère l'ensemble des temps d'évaluation
  
  for(j in 1:length(Xtime)) # Pour chaque temps d'évaluation
  {
    Xdose<-N_of_variable(media$Dose) # Récupère l'ensemble des dose présente dans le jeux de données
    assign(paste("media_time",sep=""),media[media$H %in% c(Xtime[j]),]) # Enregistre dans media_time les données avec l'heure en cours de traitement (Exemple : Enregistres les données pour le temps 72 heures)
    
    for(k in 1:length(Xdose)) # Pour chaque dose
    {
      compteur<-compteur+1 # Ajoute +1 au compteur
      assign(paste("media_time_dose",sep=""),media_time[media_time$Dose %in% c(Xdose[k]),]) # Enregistre dans media_time_dose les données avec la dose en cours de traitement (Exemple : Enregistres Milieu M199 Calf, 72h et dose de 50 nM)
      
      data_treat$DDP[compteur]<-media_time_dose$DDP[1] # Enregistre dans data_treat la date de perfusion 
      data_treat$DDE[compteur]<-media_time_dose$DDE[1] # Enregistre dans data_treat la date de d'expérience (de mesure)
      data_treat$H[compteur]<-media_time_dose$H[1] # Enregistre dans data_treat l'a date de perfusion l'heure de mesure (0, 1, 12, 24, 48, ..., 120h)
      data_treat$V[compteur]<-media_time_dose$V[1] # Enregistre dans data_treat le volume de milieu dans le puit
      data_treat$Milieu[compteur]<-media_time_dose$Milieu[1] # Enregistre dans data_treat le milieu 
      data_treat$Drug[compteur]<-media_time_dose$Drug[1] # Enregistre dans data_treat le nom de la drogue introduite
      data_treat$Dose[compteur]<-media_time_dose$Dose[1] # Enregistre dans data_treat la dose de drogue introduite
      data_treat$Unit[compteur]<-media_time_dose$Unit[1] # Enregistre dans data_treat l'unité métrique de la dose introduite
      data_treat$Puit[compteur]<-media_time_dose$Puit[1] # Enregistre dans la colonne puit la valeur 1
      
      data_treat$NST[compteur]<-mean(media_time_dose$NST) # Enregistre la moyenne du nombre de schistosome totaux
      data_treat$NSTM[compteur]<-mean(media_time_dose$NSTM) # Enregistre la moyenne du nombre de mâle schistosome totale
      data_treat$NSTM[compteur]<-mean(media_time_dose$NSTF) # Enregistre la moyenne du nombre de Femelle schistosome totale
      
      data_treat$CME4[compteur]<-mean(media_time_dose$CME4) # Enregistre la moyenne du nombre de mâle en couple à l'état E4
      data_treat$CME3[compteur]<-mean(media_time_dose$CME3) # Enregistre la moyenne du nombre de mâle en couple à l'état E3
      data_treat$CME2[compteur]<-mean(media_time_dose$CME2) # Enregistre la moyenne du nombre de mâle en couple à l'état E2
      data_treat$CME1[compteur]<-mean(media_time_dose$CME1) # Enregistre la moyenne du nombre de mâle en couple à l'état E1
      data_treat$CME0[compteur]<-mean(media_time_dose$CME0) # Enregistre la moyenne du nombre de mâle en couple à l'état E0
      
      data_treat$ME4[compteur]<-mean(media_time_dose$ME4) # Enregistre la moyenne du nombre de mâle à l'état E4
      data_treat$ME3[compteur]<-mean(media_time_dose$ME3) # Enregistre la moyenne du nombre de mâle à l'état E3
      data_treat$ME2[compteur]<-mean(media_time_dose$ME2) # Enregistre la moyenne du nombre de mâle à l'état E2
      data_treat$ME1[compteur]<-mean(media_time_dose$ME1) # Enregistre la moyenne du nombre de mâle à l'état E1
      data_treat$ME0[compteur]<-mean(media_time_dose$ME0) # Enregistre la moyenne du nombre de mâle à l'état E0
      
      data_treat$CFE4[compteur]<-mean(media_time_dose$CFE4) # Enregistre la moyenne du nombre de Femelle en couple à l'état E4
      data_treat$CFE3[compteur]<-mean(media_time_dose$CFE3) # Enregistre la moyenne du nombre de Femelle en couple à l'état E3
      data_treat$CFE2[compteur]<-mean(media_time_dose$CFE2) # Enregistre la moyenne du nombre de Femelle en couple à l'état E2
      data_treat$CFE1[compteur]<-mean(media_time_dose$CFE1) # Enregistre la moyenne du nombre de Femelle en couple à l'état E1
      data_treat$CFE0[compteur]<-mean(media_time_dose$CFE0) # Enregistre la moyenne du nombre de Femelle en couple à l'état E0
      
      data_treat$FE4[compteur]<-mean(media_time_dose$FE4) # Enregistre la moyenne du nombre de Femelle à l'état E4
      data_treat$FE3[compteur]<-mean(media_time_dose$FE3) # Enregistre la moyenne du nombre de Femelle à l'état E3
      data_treat$FE2[compteur]<-mean(media_time_dose$FE2) # Enregistre la moyenne du nombre de Femelle à l'état E2
      data_treat$FE1[compteur]<-mean(media_time_dose$FE1) # Enregistre la moyenne du nombre de Femelle à l'état E1
      data_treat$FE0[compteur]<-mean(media_time_dose$FE0) # Enregistre la moyenne du nombre de Femelle à l'état E0
      
      data_treat$Couple[compteur]<-mean(media_time_dose$Couple) # Enregistre la moyenne du nombre de vers en couple
      data_treat$Couple_Female[compteur]<-mean(media_time_dose$Couple_Female) # Enregistre la moyenne du nombre de Femelle en couple
      data_treat$Couple_Male[compteur]<-mean(media_time_dose$Couple_Male) # Enregistre la moyenne du nombre de mâle en couple
      
      data_treat$E4[compteur]<-mean(media_time_dose$E4) # Enregistre la moyenne du nombre de vers à l'état E4
      data_treat$E3[compteur]<-mean(media_time_dose$E3) # Enregistre la moyenne du nombre de vers à l'état E3
      data_treat$E2[compteur]<-mean(media_time_dose$E2) # Enregistre la moyenne du nombre de vers à l'état E2
      data_treat$E1[compteur]<-mean(media_time_dose$E1) # Enregistre la moyenne du nombre de vers à l'état E1
      data_treat$E0[compteur]<-mean(media_time_dose$E0) # Enregistre la moyenne du nombre de vers à l'état E0
      
      data_treat$Survie[compteur]<-mean(media_time_dose$Survie) # Enregistre la moyenne de survie des vers
      data_treat$Dettached[compteur]<-mean(media_time_dose$Dettached) # Enregistre le taux de détachement moyen des vers dans les puits
      
      data_treat$Inactive_Drug[compteur]<-mean(media_time_dose$Inactive_Drug) # Enregistre le taux d'inactivité moyen de la drogue
      data_treat$Moderately_Active_Drug[compteur]<-mean(media_time_dose$Moderately_Active) # Enregistre le taux d'activité modéré moyen de la drogue
      data_treat$Active_Drug[compteur]<-mean(media_time_dose$Active_Drug) # Enregistre le taux d'activité moyen de la drogue
      
      
      # Calcul les écarts types
      data_treat$dNST[compteur]<-sd(media_time_dose$NST)
      data_treat$dNSTM[compteur]<-sd(media_time_dose$NSTM)
      data_treat$dNSTM[compteur]<-sd(media_time_dose$NSTF)
      
      data_treat$dCME4[compteur]<-sd(media_time_dose$CME4)
      data_treat$dCME3[compteur]<-sd(media_time_dose$CME3)
      data_treat$dCME2[compteur]<-sd(media_time_dose$CME2)
      data_treat$dCME1[compteur]<-sd(media_time_dose$CME1)
      data_treat$dCME0[compteur]<-sd(media_time_dose$CME0)
      
      data_treat$dME4[compteur]<-sd(media_time_dose$ME4)
      data_treat$dME3[compteur]<-sd(media_time_dose$ME3)
      data_treat$dME2[compteur]<-sd(media_time_dose$ME2)
      data_treat$dME1[compteur]<-sd(media_time_dose$ME1)
      data_treat$dME0[compteur]<-sd(media_time_dose$ME0)
      
      data_treat$dCFE4[compteur]<-sd(media_time_dose$CFE4)
      data_treat$dCFE3[compteur]<-sd(media_time_dose$CFE3)
      data_treat$dCFE2[compteur]<-sd(media_time_dose$CFE2)
      data_treat$dCFE1[compteur]<-sd(media_time_dose$CFE1)
      data_treat$dCFE0[compteur]<-sd(media_time_dose$CFE0)
      
      data_treat$dFE4[compteur]<-sd(media_time_dose$FE4)
      data_treat$dFE3[compteur]<-sd(media_time_dose$FE3)
      data_treat$dFE2[compteur]<-sd(media_time_dose$FE2)
      data_treat$dFE1[compteur]<-sd(media_time_dose$FE1)
      data_treat$dFE0[compteur]<-sd(media_time_dose$FE0)
      
      data_treat$dCouple[compteur]<-sd(media_time_dose$Couple)
      data_treat$dCouple_Female[compteur]<-sd(media_time_dose$Couple_Female)
      data_treat$dCouple_Male[compteur]<-sd(media_time_dose$Couple_Male)
      
      data_treat$dE4[compteur]<-sd(media_time_dose$E4)
      data_treat$dE3[compteur]<-sd(media_time_dose$E3)
      data_treat$dE2[compteur]<-sd(media_time_dose$E2)
      data_treat$dE1[compteur]<-sd(media_time_dose$E1)
      data_treat$dE0[compteur]<-sd(media_time_dose$E0)
      
      data_treat$dSurvie[compteur]<-sd(media_time_dose$Survie)
      data_treat$dDettached[compteur]<-sd(media_time_dose$Dettached)
      
      data_treat$Inactive_Drug[compteur]<-sd(media_time_dose$Inactive_Drug)
      data_treat$dModerately_Active_Drug[compteur]<-sd(media_time_dose$Moderately_Active)
      data_treat$dActive_Drug[compteur]<-sd(media_time_dose$Active_Drug)
    }
  }
}

data_treat$Time<-data_treat$H # Enregistre le contenu de la colonne H dans la colonne Time
data_treat<-na.rm(data_treat) # Supprime les lignes contenants des NA
data_save_3<-data_treat # Sauvegarde le contenu de data_treat dans data_save_3

data_treat$Time<-as.factor(data_treat$Time) # Définis la colonne Time comme de type factor (caractère)

assign(paste("data_RPMI_treat",sep=""),data_treat[data_treat$Milieu %in% c(6),]) # Récupère uniquement les données qui comprennent le milieu RPMI Horse
assign(paste("data_M199_treat",sep=""),data_treat[data_treat$Milieu %in% c(1),]) # Récupère uniquement les données qui comprennent le milieu M199 Calf

# Graphique de la Survie des vers dans du RPMI Horse en fonction de la Dose et du temps -------------------

ggplot(data_RPMI_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0.7, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(11,12,13,14,15,16,17,18),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=40)+
  scale_x_continuous( limits=c(0-10, 900+25),breaks=seq(0, 900, by = 100))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
  ggtitle(paste("RPMI Horse IC50 Activity of drug depending on",data_RPMI_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_RPMI_treat$Unit[1],")"),colour = paste("log10(dose(",data_RPMI_treat$Unit[1],"))"))

ggsave(filename = file.path("figs",paste("RPMI Horse IC50 Activity of drug depending on",data_RPMI_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Graphique de la Survie des vers dans du M199 Calf en fonction de la Dose  et du temps -------------------

ggplot(data_M199_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0.7, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(11,12,13,14,15,16,17,18),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=40)+
  scale_x_continuous( limits=c(0-10, 900+25),breaks=seq(0, 900, by = 100))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
  
  ggtitle(paste("M199 Calf IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_M199_treat$Unit[1],")"),colour = paste("log10(dose(",data_M199_treat$Unit[1],"))"))

ggsave(filename = file.path("figs",paste("M199 Calf IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Pour travers les graphiques en échelle logarithmique, je modifie la valeur de dose 0 en valeur de dose 1

for(i in 1:length(data_M199_treat$Dose))
{
  if(data_M199_treat$Dose[i]==0)
  {
    data_M199_treat$Dose[i]<-1
  }
}

for(j in 1:length(data_RPMI_treat$Dose))
{
  if(data_RPMI_treat$Dose[j]==0)
  {
    data_RPMI_treat$Dose[j]<-1
  }
}

# Graphique de la Survie des vers dans du M199 Calf en fonction de la Dose et du temps (Echelle logarithmique) -------------------

ggplot(data_M199_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0.7, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(11,12,13,14,15,16,17,18),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=0.01)+
  scale_x_log10(limits=c(1, 1000))+
  annotation_logticks(sides = "b")+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
  
  
  ggtitle(paste("M199 Calf IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_M199_treat$Unit[1],")"),colour = paste("log10(dose(",data_M199_treat$Unit[1],"))"))

ggsave(filename = file.path("figs",paste("M199 Calf IC50 log10 Activity of drug depending on",data_M199_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Graphique de la Survie des vers dans du RPMI Horse en fonction de la Dose  et du temps (Echelle logarithmique) -------------------

ggplot(data_RPMI_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0.7, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(11,12,13,14,15,16,17,18),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=0.01)+
  annotation_logticks(sides = "b")+scale_x_log10(limits=c(1, 1000))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
  
  ggtitle(paste("RPMI Horse IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_M199_treat$Unit[1],")"),colour = paste("log10(dose(",data_M199_treat$Unit[1],"))"))
ggsave(filename = file.path("figs",paste("RPMI Horse IC50 log10 Activity of drug depending on",data_M199_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)



assign(paste("data_RPMI_treat",sep=""),data_treat[data_treat$Milieu %in% c(6),]) # Récupère uniquement les données qui comprennent le milieu RPMI Horse
assign(paste("data_M199_treat",sep=""),data_treat[data_treat$Milieu %in% c(1),]) # Récupère uniquement les données qui comprennent le milieu M199 Calf

data_RPMI_treat<-data_RPMI_treat[data_RPMI_treat$Time %in% c(120),] # Récupère uniquement les données qui correspondent au temps 120 heures avec le RPMI Horse
data_M199_treat<-data_M199_treat[data_M199_treat$Time %in% c(120),] # Récupère uniquement les données qui correspondent au temps 120 heures avec le M199 Calf

# Graphique de la Survie des vers dans du RPMI Horse au temps 120 heures et en fonction de la Dose  -------------------

ggplot(data_RPMI_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(17),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=40)+
  scale_x_continuous( limits=c(0-10, 900+25),breaks=seq(0, 900, by = 100))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
  ggtitle(paste("RPMI Horse 120h IC50 Activity of drug depending on",data_RPMI_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_RPMI_treat$Unit[1],")"),colour = paste("log10(dose(",data_RPMI_treat$Unit[1],"))"))

ggsave(filename = file.path("figs",paste("RPMI Horse 120h IC50 Activity of drug depending on",data_RPMI_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Graphique de la Survie des vers dans du M199 Calf au temps 120 heures et en fonction de la Dose  -------------------

ggplot(data_M199_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(17),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=40)+
  scale_x_continuous( limits=c(0-10, 900+25),breaks=seq(0, 900, by = 100))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+
  
  ggtitle(paste("M199 Calf 120h IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_M199_treat$Unit[1],")"),colour = paste("log10(dose(",data_M199_treat$Unit[1],"))"))

ggsave(filename = file.path("figs",paste("M199 Calf 120h IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Pour travers les graphiques en échelle logarithmique, je modifie la valeur de dose 0 en valeur de dose 1

for(i in 1:length(data_M199_treat$Dose))
{
  if(data_M199_treat$Dose[i]==0)
  {
    data_M199_treat$Dose[i]<-1
  }
}

for(j in 1:length(data_RPMI_treat$Dose))
{
  if(data_RPMI_treat$Dose[j]==0)
  {
    data_RPMI_treat$Dose[j]<-1
  }
}

# Graphique de la Survie des vers dans du M199 Calf au temps 120 heures et en fonction de la Dose (Echelle logarithmique) -------------------

ggplot(data_M199_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(17),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=0.01)+
  scale_x_log10(limits=c(1, 1000))+
  annotation_logticks(sides = "b")+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+  
  ggtitle(paste("M199 Calf 120h IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_M199_treat$Unit[1],")"),colour = paste("log10(dose(",data_M199_treat$Unit[1],"))"))

ggsave(filename = file.path("figs",paste("M199 Calf 120h IC50 log10 Activity of drug depending on",data_M199_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

# Graphique de la Survie des vers dans du RPMI Horse au temps 120 heures et en fonction de la Dose (Echelle logarithmique) -------------------

ggplot(data_RPMI_treat,aes(x=Dose,y=Survie,group=Time,color=Time,na.rm = TRUE))+
  geom_point(aes(shape=Time),size=2)+
  scale_color_grey(start=0, end=0,name = "Time (h)")+
  scale_shape_manual(values=c(17),name = "Time (h)")+
  geom_line(aes(group=Time))+
  geom_errorbar(aes(ymin=Survie-dSurvie, ymax=Survie+dSurvie),width=0.01)+
  annotation_logticks(sides = "b")+scale_x_log10(limits=c(1, 1000))+
  theme_bw() +
  theme(plot.title = element_text(size=7,face = "bold"),
        panel.border = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6,face = "bold",colour = "black"),
        axis.text.y = element_text(size=6,face = "bold",colour = "black"),
        axis.title.x = element_text(size=6,face = "bold"),
        axis.title.y = element_text(size=6,face = "bold"),
        legend.title = element_text(size=4.2,face = "bold"),
        legend.text = element_text(face = "bold",size = 4),
        legend.position = c(0.14, 0.25),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))+  
  ggtitle(paste("RPMI Horse 120h IC50 Activity of drug depending on",data_M199_treat$Drug[1],"dose"))+labs(y= "Number of worms (%)", x = paste("Dose(",data_M199_treat$Unit[1],")"),colour = paste("log10(dose(",data_M199_treat$Unit[1],"))"))
ggsave(filename = file.path("figs",paste("RPMI Horse 120h IC50 log10 Activity of drug depending on",data_M199_treat$Drug[1],"dose.TIFF")), width = 8.3, height = 8.3, units = "cm",dpi = 600)

