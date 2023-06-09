Chemin=readline(prompt = "Please indicate the path of your working directory: ")
setwd(Chemin)

source(paste(Chemin,"/fonction_Animal.R", sep =""))


listparam = c("choix_plan", "choix_seq")

mydata <- read.table(paste(Chemin,"/fichier_de_parametres_defaut.csv", sep =""), sep =";", quote = "\\"",
                     header = TRUE, stringsAsFactors = FALSE)

vecteur = mydata$Nom%in%listparam


# CREER LA MATRICE DES COMBINAISONS D'ENTREES / SCENARIOS --> expand.grid

mat.param = expand.grid(c("Ad libitum", "Restriction a 2.5 kg/j"), c("Bi-phase AB","Multiphase en 10"))
colnames(mat.param)= listparam

id=1
system.time(list.result1 <- FUNCANIMAL(id, mat.param, vecteur, mydata, listparam, Chemin))

id=2
system.time(list.result2 <- FUNCANIMAL(id, mat.param, vecteur, mydata, listparam, Chemin))

id=3
system.time(list.result3 <- FUNCANIMAL(id, mat.param, vecteur, mydata, listparam, Chemin))

id=4
system.time(list.result4 <- FUNCANIMAL(id, mat.param, vecteur, mydata, listparam, Chemin))

mat.result = rbind(list.result1, list.result2, list.result3, list.result4)

colnames(mat.result)=c("Simu", "CC", "CC_F", "CC_H", "CC_M",
                       "AC", "AC_F", "AC_H", "AC_M",
                       "EU", "EU_F", "EU_H", "EU_M",
                       "CED", "CED_F", "CED_H", "CED_M",
                       "Land", "Land_F", "Land_H", "Land_M", 
                       "CC_Eng", "CC_F_Eng", "CC_H_Eng", "CC_M_Eng",
                       "AC_Eng", "AC_F_Eng", "AC_H_Eng", "AC_M_Eng",
                       "EU_Eng", "EU_F_Eng", "EU_H_Eng", "EU_M_Eng",
                       "CED_Eng", "CED_F_Eng", "CED_H_Eng", "CED_M_Eng",
                       "Land_Eng", "Land_F_Eng", "Land_H_Eng", "Land_M_Eng",
                       "prixvente", "produit",
                       "marge_Eng_reelle", "marge_Eng_estimee",
                       "CoutAlimTot", "CoutAlim_Eng", 
                       "IC_Eng", "RejeteN_Eng", "RejeteP_Eng",
                       "MOexc_Eng", "Resdig_Eng",
                       "PVFin", "AgeFin", "Gamme",
                       "Coeur_gamme", "Lourds", "Legers", "Perte")



mat.tot = cbind(mat.param,mat.result)
save('mat.tot', file = 'mat.tot_Animal.RData')

write.table(mat.tot,"mat_tot_Animal.csv",sep=";", col.names = TRUE,
            row.names = FALSE)

