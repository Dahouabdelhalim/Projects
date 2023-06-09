library(tidyverse)
library(stringr) 
library(naniar)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(readxl)
library(caret)
library(randomForest)
library(stats)
library(performance) 
library(ggeffects) 
library(dotwhisker) 
library(ggpubr)

# La encuesta

load("~/MEGA/Manu/UPO/Curso 5º-2/Soc/WVS_Cross-National_Wave_7_R_v2_0.rdata")
View(WVS_Cross_National_Wave_7_v2_0)

# Vamos a añadir dos variables nacionales (segundo nivel) que aparecen en la literatura científica y pueden ser muy relevantes:
# Vráblíková, K. (2014). How context matters? Mobilization, political opportunity structures, and nonelectoral political participation in old and new democracies. Comparative Political Studies, 47(2), 203-229.

# Indice de descentralizacion fiscal, de 0 (centralizado) a 1 (descentralizado)

FD <- read_excel("~/MEGA/Manu/UPO/Curso 5º-2/Soc/Fiscal_Decentralization.xlsx", sheet=2)
colnames(FD) <- c("B_COUNTRY_ALPHA", "FD")
FD

# Indice de autonomia politica regional - De 0 (centralizado) a 30 (descentralizado)

RAI <- read_excel("~/MEGA/Manu/UPO/Curso 5º-2/Soc/RAI_country-2021.xlsx")
RAI <- RAI %>% filter(year == 2017) %>% select(c("abbr_country", n_RAI))
colnames(RAI) <- c("B_COUNTRY_ALPHA", "RAI")
RAI

# Nuevo data frame

df <- merge(RAI, FD, by="B_COUNTRY_ALPHA", all.y=TRUE)
wvs0 <- merge(df, WVS_Cross_National_Wave_7_v2_0, by = "B_COUNTRY_ALPHA", all.y=TRUE)
View(wvs0)

# Variable peso
wvs0 <- mutate(wvs0, peso= W_WEIGHT*pwght)

# Variable Identificación partidista
wvs0$Q223[wvs0$Q223 < 10] <- 0
wvs0$Q223[wvs0$Q223 > 10] <- 1

# Variable denominación religiosa
wvs0$Q289[wvs0$Q289>0] <- 1

# variable etnia
wvs0$ethnic <- str_sub(wvs0$Q290, start=-1)
wvs0$ethnic <- as.numeric(wvs0$ethnic)
wvs0$ethnic[wvs0$ethnic == 1] <- 0
wvs0$ethnic[wvs0$ethnic < 0] <- NA
wvs0$ethnic[wvs0$ethnic > 0] <- 1

# Variable grupo_ocupacional

wvs0$Q281[wvs0$Q281 == 0] <- "No ha trabajado"
wvs0$Q281[wvs0$Q281 == 1 | wvs0$Q281 == 2] <- "Clase de servicio"
wvs0$Q281[wvs0$Q281 > 2 & wvs0$Q281 < 6] <- "Trabajadores no manuales"
wvs0$Q281[wvs0$Q281 > 5 & wvs0$Q281 < 9] <- "Clase obrera"
wvs0$Q281[wvs0$Q281 == 9 | wvs0$Q281 == 10] <- "Trabajadores del campo"
wvs0$Q281[wvs0$Q281 == 11 | wvs0$Q281 < 0] <- NA
wvs0$Q281 <- as.factor(wvs0$Q281)
summary(wvs0$Q281)

# Variable Tipo de empleo

wvs0$Q279[wvs0$Q279 == 1] <- "Trabajador a tiempo completo"
wvs0$Q279[wvs0$Q279 == 2] <- "Trabajador a tiempo parcial"
wvs0$Q279[wvs0$Q279 == 3] <- "Autónomo"
wvs0$Q279[wvs0$Q279 > 3 & wvs0$Q279 < 7] <- "Inactivo"
wvs0$Q279[wvs0$Q279 == 7] <- "Desempleado"
wvs0$Q279[wvs0$Q279 == 8 | wvs0$Q279 < 0] <- NA
wvs0$Q279 <- as.factor(wvs0$Q279)
summary(wvs0$Q279)

# Variable Nivel educativo

wvs0$Q275R[wvs0$Q275R == 1] <- "Primaria"
wvs0$Q275R[wvs0$Q275R == 2] <- "Secundaria"
wvs0$Q275R[wvs0$Q275R == 3] <- "Pos-secundaria"
wvs0$Q275R[wvs0$Q275R == 4] <- "Terciaria"
wvs0$Q275R[wvs0$Q275R < 0] <- NA
wvs0$Q275R <- as.factor(wvs0$Q275R)
summary(wvs0$Q275R)


# Marcamos los valores sin respuesta y no procede por NA

wvs0[,1:362][wvs0[,1:362] <0] <- NA
wvs0[,363:401][wvs0[,363:401]== -99] <- NA
wvs0[,402:515][wvs0[,402:515]== -9999] <- NA
wvs0[,516:550][wvs0[,516:550]== -999] <- NA


# Realizamos algunos cambios para que algunas variables sean positivas conforme el valor sea mayor

wvs0 <- mutate(wvs0, Q4 = (Q4-4)*-1) 
wvs0 <- mutate(wvs0, Q199 = (Q199-4)*-1)
wvs0 <- mutate(wvs0, Q200 = (Q200-3)*-1)
wvs0 <- mutate(wvs0, Q202 = (Q202-5)*-1)
wvs0 <- mutate(wvs0, Q206 = (Q206-5)*-1)
wvs0 <- mutate(wvs0, Q207 = (Q207-5)*-1)
wvs0 <- mutate(wvs0, Q224 = (Q224-4)*-1)
wvs0 <- mutate(wvs0, Q234 = (Q234-4)*-1)
wvs0 <- mutate(wvs0, Q235 = (Q235-4)*-1)


# Regiones de Norte y Sur de América

table(wvs0$B_COUNTRY_ALPHA, wvs0$regionWB)

wvs0 <- filter(wvs0, regionWB == 3 | regionWB == 5)


# Ahora hacemos la variable dependiente.
# Anduiza, Eva y Bosch, Antoni (2019) Comportamiento político y electoral, Barcelona, Editorial Ariel.
# Van Deth, Jan W. (2014). "A conceptual map of political participation", en Acta politica, vol. 43, no. 3, pp. 349-367.

unique(wvs0$Q98)
wvs0$partido[wvs0$Q98 == 0] <- 0
wvs0$partido[wvs0$Q98 > 0] <- 1
wvs0$partido <- as.numeric(wvs0$partido)
wvs0$partido # Activo o inactivo frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q209)
wvs0$Q209rec[wvs0$Q209 == 1] <- 1
wvs0$Q209rec[wvs0$Q209 > 1] <- 0
wvs0$Q209rec <- as.numeric(wvs0$Q209rec)
wvs0$Q209rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q210)
wvs0$Q210rec[wvs0$Q210 == 1] <- 1
wvs0$Q210rec[wvs0$Q210 > 1] <- 0
wvs0$Q210rec <- as.numeric(wvs0$Q210rec)
wvs0$Q210rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q211)
wvs0$Q211rec[wvs0$Q211 == 1] <- 1
wvs0$Q211rec[wvs0$Q211 > 1] <- 0
wvs0$Q211rec <- as.numeric(wvs0$Q211rec)
wvs0$Q211rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q212)
wvs0$Q212rec[wvs0$Q212 == 1] <- 1
wvs0$Q212rec[wvs0$Q212 > 1] <- 0
wvs0$Q212rec <- as.numeric(wvs0$Q212rec)
wvs0$Q212rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q213)
wvs0$Q213rec[wvs0$Q213 == 1] <- 1
wvs0$Q213rec[wvs0$Q213 > 1] <- 0
wvs0$Q213rec <- as.numeric(wvs0$Q213rec)
wvs0$Q213rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q214)
wvs0$Q214rec[wvs0$Q214 == 1] <- 1
wvs0$Q214rec[wvs0$Q214 > 1] <- 0
wvs0$Q214rec <- as.numeric(wvs0$Q214rec)
wvs0$Q214rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q215)
wvs0$Q215rec[wvs0$Q215 == 1] <- 1
wvs0$Q215rec[wvs0$Q215 > 1] <- 0
wvs0$Q215rec <- as.numeric(wvs0$Q215rec)
wvs0$Q215rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q216)
wvs0$activismo[wvs0$Q216 == 1] <- 1
wvs0$activismo[wvs0$Q216 > 1] <- 0
wvs0$activismo <- as.numeric(wvs0$activismo)
wvs0$activismo # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q217)
wvs0$Q217rec[wvs0$Q217 == 1] <- 1
wvs0$Q217rec[wvs0$Q217 > 1] <- 0
wvs0$Q217rec <- as.numeric(wvs0$Q217rec)
wvs0$Q217rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q218)
wvs0$Q218rec[wvs0$Q218 == 1] <- 1
wvs0$Q218rec[wvs0$Q218 > 1] <- 0
wvs0$Q218rec <- as.numeric(wvs0$Q218rec)
wvs0$Q218rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q219)
wvs0$Q219rec[wvs0$Q219 == 1] <- 1
wvs0$Q219rec[wvs0$Q219 > 1] <- 0
wvs0$Q219rec <- as.numeric(wvs0$Q219rec)
wvs0$Q219rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q220)
wvs0$Q220rec[wvs0$Q220 == 1] <- 1
wvs0$Q220rec[wvs0$Q220 > 1] <- 0
wvs0$Q220rec <- as.numeric(wvs0$Q220rec)
wvs0$Q220rec # Lo ha hecho frente a no lo ha hecho o podría hacerlo

unique(wvs0$Q221)
wvs0$voto_local[wvs0$Q221 == 4] <- NA
wvs0$voto_local[wvs0$Q221 == 3] <- 0
wvs0$voto_local[wvs0$Q221 == 1 | wvs0$Q221 == 2] <- 1
wvs0$voto_local <- as.numeric(wvs0$voto_local)
wvs0$voto_local # Siempre o frecuentemente vota frente a voto_Nacional

unique(wvs0$Q222)
wvs0$voto_nacional[wvs0$Q222 == 4] <- NA
wvs0$voto_nacional[wvs0$Q222 == 3] <- 0
wvs0$voto_nacional[wvs0$Q222 == 1 | wvs0$Q222 == 2] <- 1
wvs0$voto_nacional <- as.numeric(wvs0$voto_nacional)
wvs0$voto_nacional # Siempre o frecuentemente vota frente a nunca 

# Hacemos un ACP para ver como clasificar las variables

df1 <- PCA(wvs0[,c(553:567)])
df1$eig
df1$var

# La única separación razonable es voto vs resto de actividades

wvs0 <- mutate(wvs0, voto=((wvs0$voto_nacional + wvs0$voto_local)/2))
summary(wvs0$voto)

wvs0 <- mutate(wvs0, zPPOL = ((wvs0$partido + wvs0$Q213rec + wvs0$Q214rec + wvs0$Q209rec + wvs0$Q217rec 
                                + wvs0$Q218rec +  wvs0$Q219rec  + wvs0$Q220rec +wvs0$Q210rec + wvs0$Q211rec + 
                                  wvs0$Q212rec + wvs0$activismo)/12))
summary(wvs0$zPPOL)


PCA(wvs0[,c(568:569)])

# Crepaz, M. M., Jazayeri, K. B., & Polk, J. (2017). What's trust got to do with it? The effects of in-group and out-group trust on conventional and unconventional political participation. Social Science Quarterly, 98(1), 261-281.
# Ahora vamos a crear otras dos variables importantes en la literatura: la confianza endogena y exogena. Fuentes:

PCA(wvs0[,99:104])

# Todos valores de 1 (plenamente) a 4 (nada) pasados como 0 (nada) a 1 (plenamente)
wvs0 <- mutate(wvs0, endtrust = ((wvs0$Q58-4) + (wvs0$Q59-4) +(wvs0$Q60-4))/-9)
wvs0$endtrust

wvs0 <- mutate(wvs0, exotrust = ((wvs0$Q61-4) + (wvs0$Q62-4) +(wvs0$Q63-4))/-9)
wvs0$exotrust

# Hacemos ACP de nuevo

PCA(wvs0[,571:572])


# Participación en organizaciones no directamente políticas

unique(wvs0$Q94)
wvs0$Q94[wvs0$Q94>0] <- 1
wvs0$Q95[wvs0$Q95>0] <- 1
wvs0$Q96[wvs0$Q96>0] <- 1
wvs0$Q97[wvs0$Q97>0] <- 1
wvs0$Q99[wvs0$Q99>0] <- 1
wvs0$Q100[wvs0$Q100>0] <- 1
wvs0$Q101[wvs0$Q101>0] <- 1
wvs0$Q102[wvs0$Q102>0] <- 1
wvs0$Q103[wvs0$Q103>0] <- 1
wvs0$Q104[wvs0$Q104>0] <- 1
wvs0$Q105[wvs0$Q105>0] <- 1

wvs0 <- wvs0 %>% mutate(no_polpart = c(Q94 + Q95 + Q96 + Q97 + Q99 + Q100
                                       + Q101 + Q102 + Q103 + Q104 + Q105)) 
wvs0$no_polpart

# Confianza Política

unique(wvs0$Q71)
wvs0 <- wvs0 %>% mutate(poltrust = ((c(Q71 + Q72 + Q73 + Q76)-16)*-1))
wvs0$poltrust

# Eliminamos columnas
wvs0 <- wvs0[-c(4:23, 27:38, 44, 46, 49, 55, 65, 67:68, 71:72, 74, 76, 86, 99:104, 106, 110:114, 117, 124:139, 151:162,
                209:214, 228, 230, 235, 237:239, 241:242, 266:279, 281:282, 314, 321, 323:324, 328:330, 334,
                337:338, 340:341, 343:344, 346:347, 358, 361:362, 364:365, 368:402, 417, 431:432, 435:437, 444:445, 450, 459,
                465:467, 469, 516:521, 553:567)]

# Eliminamos variables con un alto porcentaje de NA
wvs0 <- wvs0[colSums(is.na(wvs0))/nrow(wvs0) < .3]


# Guardamos el archivo

save(wvs0, file="wvs0.rdata")
load(file = "~/MEGA/Manu/UPO/Curso 5º-2/Soc/wvs0.rdata")

# Creamos los descriptivos de la variable dependiente

df <- wvs0 %>% group_by(B_COUNTRY_ALPHA) %>% summarise(voto = mean(voto, na.rm=TRUE))
ggplot(df, aes(x=voto, y=reorder(B_COUNTRY_ALPHA, voto), fill=B_COUNTRY_ALPHA)) + 
  geom_col(show.legend = FALSE) + labs(y=NULL) + theme_minimal()

df <- wvs0 %>% group_by(B_COUNTRY_ALPHA) %>% summarise(zPPOL = mean(zPPOL, na.rm=TRUE))
ggplot(df, aes(x=zPPOL, y=reorder(B_COUNTRY_ALPHA, zPPOL), fill=B_COUNTRY_ALPHA)) + 
  geom_col(show.legend = FALSE) + labs(y=NULL) + theme_minimal()

# Recodificamos el voto como dicotómica
wvs0$voto[wvs0$voto > 0] <- 1
wvs0$voto <- as.factor(wvs0$voto) 


# Separamos el fichero

wvsNA <- filter(wvs0, regionWB==3)
wvsLA <- filter(wvs0, regionWB==5)

unique(wvsNA$B_COUNTRY_ALPHA)
unique(wvsLA$B_COUNTRY_ALPHA)

wvsNA <- wvsNA[,-237]
wvsLA <- wvsLA[,-237]


# Bosque Aleatorio

wvsNA <- wvsNA[,-1]
wvsLA <- wvsLA[,-1]

wvsNA <- na.roughfix(wvsNA)
wvsLA <- na.roughfix(wvsLA)

a <- wvsNA$peso
b <- wvsLA$peso

wvsNA <- wvsNA[,-309]
wvsLA <- wvsLA[,-309]


help("randomForest")
# MODELO 1 NA - voto

set.seed(25)
model1 <- randomForest(voto ~ ., data=wvsNA[,-311], weights=a, proximity=TRUE)
model1
i_scores <- varImp(model1, conditional=TRUE, scale=TRUE)
i_scores

model2 <- randomForest(voto ~ ., data=wvsLA[,-311], weights=b, proximity=TRUE)
model2
i_scores <- varImp(model2, conditional=TRUE, scale=TRUE)
i_scores

model3 <- randomForest(zPPOL ~ ., data=wvsNA[,-310], weights=a)
model3
i_scores <- varImp(model3, conditional=TRUE, scale=TRUE)
i_scores

model4 <- randomForest(zPPOL ~ ., data=wvsLA[,-310], weights=b)
model4
i_scores <- varImp(model4, conditional=TRUE, scale=TRUE)
i_scores


# Recodificación de principales variables



wvsNA <- wvsNA %>% mutate(Importancia_Politica = (Q4 - mean(Q4))/sd(Q4), 
                          Interés_Política = (Q199 - mean(Q199))/sd(Q199), 
                          Elecciones_honestas_afectan = (Q234 - mean(Q234))/sd(Q234), 
                          Identificación_partidista = Q223, 
                          Votos_justos = (Q224 - mean(Q224))/sd(Q224), 
                          Grupo_ocupacional = Q281, 
                          Edad = (Q262 - mean(Q262))/sd(Q262), 
                          Importancia_democracia = (Q250 - mean(Q250))/sd(Q250),
                          Valores_emancipatorios = (resemaval - mean(resemaval))/sd(resemaval),  
                          Valores_seculares = (sacsecval - mean(sacsecval))/sd(sacsecval), 
                          Líderes_elecciones_libres = (Q243 - mean(Q243))/sd(Q243), 
                          Afiliación_org_cívicas = (no_polpart - mean(no_polpart))/sd(no_polpart),
                          Tipo_empleo = Q279,
                          Información_TV = (Q202 - mean(Q202))/sd(Q202), 
                          N_hijos = (Q274 - mean(Q274))/sd(Q274), 
                          N_personas_hogar = (Q270 - mean(Q270))/sd(Q270), 
                          Derechos_civiles = (Q246 - mean(Q246))/sd(Q246),
                          Discusion_pol_amigos = (Q200 - mean(Q200))/sd(Q200),
                          Ideología = (Q240 - mean(Q240))/sd(Q240),
                          Líder_fuerte_autocrático = (Q235 - mean(Q235))/sd(Q235),
                          Información_RRSS = (Q207 - mean(Q207))/sd(Q207),
                          Información_Internet = (Q206 - mean(Q206))/sd(Q206),
                          Nivel_Educativo = Q275R)

wvsLA <- wvsLA %>% mutate(Importancia_Politica = (Q4 - mean(Q4))/sd(Q4), 
                          Interés_Política = (Q199 - mean(Q199))/sd(Q199), 
                          Elecciones_honestas_afectan = (Q234 - mean(Q234))/sd(Q234), 
                          Identificación_partidista = Q223, 
                          Votos_justos = (Q224 - mean(Q224))/sd(Q224), 
                          Grupo_ocupacional = Q281, 
                          Edad = (Q262 - mean(Q262))/sd(Q262), 
                          Importancia_democracia = (Q250 - mean(Q250))/sd(Q250),
                          Valores_emancipatorios = (resemaval - mean(resemaval))/sd(resemaval),  
                          Valores_seculares = (sacsecval - mean(sacsecval))/sd(sacsecval), 
                          Líderes_elecciones_libres = (Q243 - mean(Q243))/sd(Q243), 
                          Afiliación_org_cívicas = (no_polpart - mean(no_polpart))/sd(no_polpart),
                          Tipo_empleo = Q279,
                          Información_TV = (Q202 - mean(Q202))/sd(Q202), 
                          N_hijos = (Q274 - mean(Q274))/sd(Q274), 
                          N_personas_hogar = (Q270 - mean(Q270))/sd(Q270),
                          Derechos_civiles = (Q246 - mean(Q246))/sd(Q246),
                          Discusion_pol_amigos = (Q200 - mean(Q200))/sd(Q200),
                          Ideología = (Q240 - mean(Q240))/sd(Q240),
                          Líder_fuerte_autocrático = (Q235 - mean(Q235))/sd(Q235),
                          Información_RRSS = (Q207 - mean(Q207))/sd(Q207),
                          Información_Internet = (Q206 - mean(Q206))/sd(Q206),
                          Nivel_Educativo = Q275R)



# MODELO 1 NA - voto

model1 <- randomForest(voto ~ Interés_Política +  Elecciones_honestas_afectan +  Identificación_partidista + Edad  +
                         Importancia_democracia + Votos_justos + Importancia_Politica + Líderes_elecciones_libres +
                         Valores_seculares + Valores_emancipatorios, data=wvsNA[,-311], weights=a, proximity=TRUE)
model1

i_scores <- varImp(model1, conditional=TRUE, scale=TRUE)
i_scores <- i_scores %>% tibble::rownames_to_column("var") 
i_scores <- i_scores[order(-i_scores$Overall),]
i_scores$var<- i_scores$var %>% as.factor()
i_scores <- i_scores[1:10,]
i_scores %>% ggplot(aes(x=Overall, y=reorder(var, Overall), fill=var)) + 
  geom_col(show.legend = FALSE) + theme_minimal() + labs(x = NULL, y = NULL)



# MODELO 2 LA - voto


model2 <- randomForest(voto ~ Edad + Valores_seculares + Valores_emancipatorios + Información_TV + N_hijos + 
                         Grupo_ocupacional + Tipo_empleo + N_personas_hogar + endtrust + Elecciones_honestas_afectan, data=wvsLA[,-311], weights=b, proximity=TRUE)
model2

i_scores <- varImp(model2, conditional=TRUE, scale=TRUE)
i_scores <- i_scores %>% tibble::rownames_to_column("var") 
i_scores <- i_scores[order(-i_scores$Overall),]
i_scores$var<- i_scores$var %>% as.factor()
i_scores <- i_scores[1:10,]
i_scores %>% ggplot(aes(x=Overall, y=reorder(var, Overall), fill=var)) + 
  geom_col(show.legend = FALSE) + theme_minimal() + labs(x = NULL, y = NULL)


# MODELO 3 NA - PPOL

model3 <- randomForest(zPPOL ~ Valores_seculares + Valores_emancipatorios + Información_RRSS + Interés_Política + Derechos_civiles +
                         Ideología + Importancia_Politica + Líder_fuerte_autocrático + Afiliación_org_cívicas + Discusion_pol_amigos, data=wvsNA[,-310], weights=a)
model3


i_scores <- varImp(model3, conditional=TRUE, scale=TRUE)
i_scores <- i_scores %>% tibble::rownames_to_column("var") 
i_scores <- i_scores[order(-i_scores$Overall),]
i_scores$var<- i_scores$var %>% as.factor()
i_scores <- i_scores[1:10,]
i_scores %>% ggplot(aes(x=Overall, y=reorder(var, Overall), fill=var)) + 
  geom_col(show.legend = FALSE) + theme_minimal() + labs(x = NULL, y = NULL)

# MODELO 4 LA - PPOL


model4 <- randomForest(zPPOL ~ Edad + Valores_seculares + Valores_emancipatorios + Información_Internet + Interés_Política + 
                         Grupo_ocupacional + Importancia_Politica + Nivel_Educativo + no_polpart + Discusion_pol_amigos, data=wvsLA[,-310], weights=b)
model4

i_scores <- varImp(model4, conditional=TRUE, scale=TRUE)
i_scores <- i_scores %>% tibble::rownames_to_column("var") 
i_scores <- i_scores[order(-i_scores$Overall),]
i_scores$var<- i_scores$var %>% as.factor()
i_scores <- i_scores[1:10,]
i_scores %>% ggplot(aes(x=Overall, y=reorder(var, Overall), fill=var)) + 
  geom_col(show.legend = FALSE) + theme_minimal() + labs(x = NULL, y = NULL)



# modelo logit 

# MODELO 1 NA - voto

model10 <- glm(voto ~ Interés_Política +  Elecciones_honestas_afectan +  Identificación_partidista + Edad  +
                Importancia_democracia + Votos_justos + Importancia_Politica + Líderes_elecciones_libres +
                Valores_seculares + Valores_emancipatorios,
              family= quasibinomial, data=wvsNA, weights=a)
summary(model10)
r2_nagelkerke(model10)

model11 <- glm(voto ~ Interés_Política +  Elecciones_honestas_afectan +  Identificación_partidista + Edad  +
                Importancia_democracia + Votos_justos + Importancia_Politica + Líderes_elecciones_libres +
              Valores_seculares + Valores_emancipatorios + Edad*Elecciones_honestas_afectan + Edad*Valores_seculares,
              family= quasibinomial, data=wvsNA, weights=a)
summary(model11)
r2_nagelkerke(model11)

dwplot(list(model10, model11), ci_method="wald",
       vline = geom_vline(
         xintercept = 0,
         colour = "grey60",
         linetype = 2))


# MODELO 2 LA - voto



model20 <- glm(voto ~ Edad + Valores_seculares + Valores_emancipatorios + Información_TV + N_hijos +
                 Grupo_ocupacional + Tipo_empleo + N_personas_hogar + endtrust + Elecciones_honestas_afectan,
               family= quasibinomial, data=wvsLA, weights=b)
summary(model20)
r2_nagelkerke(model20)

model21 <- glm(voto ~ Edad + Valores_seculares + Valores_emancipatorios + Información_TV + N_hijos +
                Grupo_ocupacional + Tipo_empleo + N_personas_hogar + endtrust + Elecciones_honestas_afectan +
                Grupo_ocupacional*endtrust + Grupo_ocupacional*Elecciones_honestas_afectan + Edad*N_personas_hogar + Edad*Elecciones_honestas_afectan + Edad*Valores_seculares,
              family= quasibinomial, data=wvsLA, weights=b)
summary(model21)
r2_nagelkerke(model21)

dwplot(list(model20, model21), ci_method="wald",
       vline = geom_vline(
         xintercept = 0,
         colour = "grey60",
         linetype = 2))


# modelos lineales


wvsNA <- mutate(wvsNA, Ideología_sqr = Ideología^2)
model30 <- lm(zPPOL ~ Valores_seculares + Valores_emancipatorios + Información_RRSS + Interés_Política + Derechos_civiles +
                Ideología + Ideología_sqr + Importancia_Politica + Líder_fuerte_autocrático + Afiliación_org_cívicas + Discusion_pol_amigos, data=wvsNA[,-310], weights=a)
summary(model30)
r2(model30)

dwplot(model30, ci_method="wald",
       vline = geom_vline(
         xintercept = 0,
         colour = "grey60",
         linetype = 2))

par(mfrow = c(2, 2))
plot(model30)

# MODELO 4 LA - PPOL


model40 <- lm(zPPOL ~ Edad + Valores_seculares + Valores_emancipatorios + Información_Internet + Interés_Política + 
               Grupo_ocupacional + Importancia_Politica + Nivel_Educativo + no_polpart + Discusion_pol_amigos, data=wvsLA[,-310], weights=b, proximity=TRUE)
summary(model40)
r2(model40)

par(mfrow = c(2, 2))
plot(model40)


model41 <- lm(zPPOL ~ Edad + Valores_seculares + Valores_emancipatorios + Información_Internet + Interés_Política + 
                Grupo_ocupacional + Importancia_Politica + Nivel_Educativo + no_polpart + Discusion_pol_amigos +
                Edad*Valores_emancipatorios + Edad*Importancia_Politica + Edad*no_polpart +
                Grupo_ocupacional*Valores_emancipatorios + Grupo_ocupacional*Importancia_Politica + Grupo_ocupacional*no_polpart +
                Nivel_Educativo*Valores_emancipatorios + Nivel_Educativo*Importancia_Politica + Nivel_Educativo*no_polpart, data=wvsLA[,-310], weights=b, proximity=TRUE)
summary(model41)
r2(model41)


dwplot(list(model40, model41), ci_method="wald",
       vline = geom_vline(
         xintercept = 0,
         colour = "grey60",
         linetype = 2))
