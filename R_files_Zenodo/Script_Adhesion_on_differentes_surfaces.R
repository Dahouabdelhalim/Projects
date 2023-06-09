# Initialise les librairies  -----------------------------------

rm(list=ls()) # Supprime les données présentes dans la ram et la rom de RStudio
packages<-c("here","sp") # Chemin du répertoire avec les données 

# Installation des packages (si besoin)
for(i in packages){if(!require(i,character.only = T)){install.packages(i,dep=T)}}
lapply(packages, library, character.only = T) # Chargement de la liste des packages

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

opentxt<-function(x)# Fonction utilisée pour ouvrir les fichiers textes
{
  data<-read.table(here("data",x),sep='\\t', header=T) # Ouvre et lit les fichiers txt pour les enregistrer dans une dataframe
}
txt<-dir("data",pattern = ".txt") #Enregistre le nombre de fichier texte et leur nom dans un tableau
N<-length(txt) # Enregistre dans la variable N le nombre total de fichier texte
t_jours<-NA
data_real<-structure(list(character()), class = "data.frame") # Création de la dataframe data_real

for (i in 1:N) # Ouvre tous les fichiers textes et les enregistre dans la dataframe data_real
{
  data_real<-rbind(data_real,assign(paste("data", i, sep = ""),opentxt(txt[i]))) # Ouvre les fichiers txt et les enregistres dans data_real
  x<-assign(paste("data", i, sep = ""),opentxt(txt[i])) # Pendant la boucle, enregistre le dernier fichier txt lu dans la dataframe x
}

N_Dose<-function(Variable) # Récupère le nombre dose dans un tableau
{
  compteur<-1 # Prépare un compteur à 1
  table<-NA # Création d'un tableau
  table[compteur]<--1 # Donne une valeur de -1 à la colonne compteur du tableau 
  dose<-NA #Donne une valeur de dose à NA
  exist<-0 # Création de la variable exist à 0
  for (i in 1:length(Variable)) # Pour chaque ligne contenue dans Variable  
  {
    dose<-Variable[i] # Enregistre dans dose la valeur dose présent dans variable à la ligne i 
    for (j in 1:length(table)) # De la première colonne à la dernière colonne du tableau table
    {
      if(dose==table[j]) # Vérifie si la dose est déjà présente dans la table
      {
        j<-length(table) # Si oui, arrêt la boucle for
        exist<-1 # Indique que la valeur est déjà contenue dans la table
      }
    }
    if(exist==0) # Si la valeur n'existe pas dans la table
    {
      table[compteur]<-dose # Enregistre la valeur dans la colonne « compteur » de la table 
      compteur<-compteur+1 # Change de colonne « compteur » de la table 
      print(dose) # Affiche la valeur ajoutée dans la table
    }
    else(exist<-0) # Sinon, réinitialise la valeur exist a 0 pour la lecture de la ligne suivante
  }
  
  Ndose<-(sort(table)) # Range les valeurs de dose contenues dans l'ordre croissant 
  return(Ndose) # Renvoi les données de la table de dose en sortie de la fonction.
}

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

#Code---------------

data<-na.rm(data_real) # Sauvegarde les données de data_real dans data en supprimant les lignes contenants des NA
data$vers<-data$Male.Total+data$Female.Total+data$Couple.totaux*2 # Calcul le nombre de vers total dans les puces (pour toutes les lignes de data)
data$Day<-as.numeric(data$Day) # Indique que la colonne Day de data est numérique
data$Couple.totaux<-data$Couple.totaux*2 # Multiplie par 2 pour récupérer le nombre exacte de vers en couple

Jmax<-max(data$Day) # Permet de récupérer le nombre d'expérience
Table<-c("Collagen","Glass","MDCK","Collagen + MDCK")# Table enregistre le nom de chaque type de surfaçage
data$Day<-as.numeric(data$Day)
data<-na.rm(data) # Supprimes les lignes avec des NA dans la dataframe
x<-NA

for(i in 1:max(data$Day)) # Pour chaque jour
{
  day<-data[data$Day %in% i,] # Eclates les données par jours
  for (j in 1:length(Table)) # Pour chaque type de surfaçage
  {
    Type<-day[day$Type %in% Table[j],] # Eclaté les données par type de surfaçage
    Type$vers<-(Type$vers/max(Type$vers))*100 # Calcul le pourcentage de vers présent dans chaque puce au cours du temps
    x<-rbind(x,Type) # Enregistre les données éclatées dans la dataframe x à la fin de chaque exécution de la boucle for
  }
}

data<-na.rm(x) # Supprimes les lignes avec des NA suites à la fonction Rbind puis enregistre les données dans data
save_for_sigmaplot<-data # Sauvegarde les données dans la dataframe save_for_sigmaplot

#---------------Stat

Stat<-data
Stat$Flow<-Stat$Flow.rate..mL.min.
Stat<-subset(Stat, select = -c(3,4,5,6,7,8,9,10,11,12))
Stat$Type <- factor(Stat$Type)
Stat$Day <- factor(Stat$Day)
Stat$Flow<- as.numeric(Stat$Flow)


NFlow<-N_of_variable(Stat$Flow)

#flow

Flow<-subset(Stat,Flow==1)#0.7 à 2.4 2.8
table(Flow$Type,Flow$Day,useNA="ifany")
stripchart(vers~Type,data=Flow,method="jitter",cex.axis=0.8,vertical=TRUE)
AnovaFlow<-aov(vers~Type+Day,data=Flow)#AnovaFlow<-aov(vers~Type,data=Flow)
qqnorm(residuals(AnovaFlow))
qqline(residuals(AnovaFlow))
print(fitted(AnovaFlow))
plot(rstudent(AnovaFlow)~fitted(AnovaFlow))
abline(h=0,lty="dotted")
print(summary(AnovaFlow))

# Avec la fonction glht
comparaisons <- glht(AnovaFlow,linfct=mcp(Type="Tukey"))
# Afficher les comparaisons
summary(comparaisons)
cld(comparaisons)
tapply(Flow$vers,Flow$Type,mean)
R2 <- summary(AnovaFlow)[[1]][["Sum Sq"]]
R2
sum(R2[1:(length(R2)-1)])/sum(R2)


#---------------Fin Stat





data<-save_for_sigmaplot
Nflow<-N_Dose(data$Flow.rate..mL.min.) # Récupéré chaque valeur de débit dans la colonne Nflow de tableau

Tableau<-data # Enregistre data dans Tableau
Tableau$Day<-NA # Change les données contenue dans la colonne Day par des NA
Tableau$SD<-NA # Création de la colonne des écarts type SD (NA comme valeur par défaut)
ligne<-1 # Initialise la variable ligne à 1

for (i in 1:length(Table)) # Pour chaque type de surfaçage
{
  Type<-data[data$Type %in% Table[i],] # Eclate les données par type de surfaçage
  for(j in 1:length(Nflow)) # Pour chaque valeur de débit
  {
    Flow<-Type[Type$Flow.rate..mL.min. %in% Nflow[j],] # Eclaté les données pour chaque valeur de débit
    Tableau$vers[ligne]<-mean(Flow$vers) # Calcul la moyenne du nombre de vers pour un débit et un type de surfaçage précis (ex: calcul de la moyenne pour toutes les puces au Collagène à un débit de 1mL/min)
    Tableau$SD[ligne]<-sd(Flow$vers) # Calcul l'écart type du nombre de vers pour un débit et un type de surfaçage donnée
    Tableau$Flow[ligne]<-Nflow[j] # Enregistre la valeur du débit testé
    Tableau$Day[ligne]<-0 # Enregistre dans la colonne Day une valeur à 0
    ligne<-ligne+1 # Change de ligne dans la dataframe Tableau pour le prochain calcul
  }
}
Tableau<-na.rm(Tableau) # Supprime tous les lignes contenants des NA 
Tableau<-subset(Tableau,select = c(1,13,14,15)) # Conserve que les colonnes Type, vers (%), Flow (mL/min) et les écart type (SD) 

data<-Tableau # Ecrase data par les données contenues dans Tableau
levels(data$Type)<-c("Collagen","Glass","MDCK","Collagen + MDCK") # Indique le nom de chaque type de surface en fonction du numéro contenue dans la colonne type
data$Type <- factor(data$Type, levels = c("Collagen + MDCK","Collagen","MDCK","Glass")) # Impose un ordre d'affichage des type de surfaçage dans la légende du graphique

data$Ymax<-NA
data$Ymin<-NA
data$SDmin<-NA
data$SDmax<-NA

for(i in 1:length(data$vers))
{
  data$Ymax[i]<-data$vers[i]+data$SD[i]
  data$Ymin[i]<-data$vers[i]-data$SD[i]
  if(data$Ymax[i]>100)
  {
    data$SDmax[i]<-data$SD[i]-(data$Ymax[i]-100)
  }
  else
  {
    data$SDmax[i]<-data$SD[i]
  }
  if(data$Ymin[i]<0)
  {
    data$SDmin[i]<-sqrt((0-data$vers[i])^2)
  }
  else
  {
    data$SDmin[i]<-data$SD[i]
  }
}

# Trace le graphique du nombre de vers dans les puces en fonction du débit et du type de surfaçage 
Xlab<-expression (bold(Flow~rate~(mL.min^{-1})))
Ylab<-expression (bold(Attached~worms~"(%)"))


ggplot(data,aes(x=Flow,y=vers, group=Type,color=Type,na.rm = TRUE))+
  geom_step()+
  geom_errorbar(aes(ymin=vers-SDmin, ymax=vers+SDmax),width=0.085)+
  geom_point(aes(shape=Type),size=1.5)+
  scale_shape_manual(values=c(8,15,19,17),name = "Type of coating")+
  scale_y_continuous( limits=c(-5, 100),breaks=seq(0, 100, by = 20))+
  scale_x_continuous(breaks=seq(0, 3.5, by = 0.5))+
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
        legend.position = c(0.78, 0.87),
        legend.background = element_rect(fill = "transparent"),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "lines"))+
  labs(y= Ylab, x = Xlab,color="Type of coating")+
  geom_segment(aes(x=0.7, xend=2.8, y=-4, yend=-4), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=0.7, xend=0.7, y=-5.0, yend=-3), colour="black", size = 0.5, show.legend=F)+
  geom_segment(aes(x=2.8, xend=2.8, y=-5.0, yend=-3), colour="black", size = 0.5, show.legend=F)+
  annotate("text", x = 1.75, y = -1.4, label = "p-value < 0.05 : Collagen vs Without Collagen", colour="black", size = 2)

ggsave(dpi=600,filename = file.path("figs","Figure 5.TIFF"), width = 8.3, height = 8.3, units = "cm")
ggsave(dpi=600,filename = file.path("figs","Figure 5.pdf"), width = 8.3, height = 8.3, units = "cm")

data<-save_for_sigmaplot # Charge les données sauvegardées dans save_for_sigmaplot dans la dataframe data
data$Flow<-data$Flow.rate..mL.min.


data<-subset(data,select = c(1,2,13,14)) # Conserve que les colonnes Type, vers (%), Flow (mL/min) et day
data1<-data[data$Day %in% c(1),] # Extrait que les données contenue au jour 1 et l'enregistre dans data1
data2<-data[data$Day %in% c(2),] #Extrait que les données contenue au jour 2 et l'enregistre dans data2

data1<-subset(data1, select = c(1,3,4)) # Supprime la colonne day
data2<-subset(data2, select = c(1,3,4)) # Supprime la colonne day

data1<-data1 %>% spread(key= Type, value= vers) # Oriente la disposition du tableau data1
data2<-data2 %>% spread(key= Type, value= vers) # Oriente la disposition du tableau data2

write_xlsx(data1,path = "feuille_exp_1_adhesion.xlsx",col_names = TRUE) # Exporte data1 dans un fichier Excel
write_xlsx(data2,path = "feuille_exp_2_adhesion.xlsx",col_names = TRUE) # Exporte data2 dans un fichier Excel

