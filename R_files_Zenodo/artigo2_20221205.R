#rodada 12: teste
library(epiDisplay)
library(rgdal)
library(spdep)
library(tmap)
library(tmaptools)
library(OpenStreetMap)
library(rgdal)
library(rgeos)
library(sf)
library(sjPlot)

#load("E:/IESC/Doutorado/analises/objetivo3/analises3_20191231.RData")
#load("E:/IESC/Doutorado/analises/objetivo2/analises2_10122020.RData")

load("E:/IESC/Doutorado/analises/objetivo3/analises3_20191231.RData")

dfdengue <- readOGR("E:/IESC/Eny/gwr/Mrj6.shp", stringsAsFactors=FALSE)

zktotal <- merge(zktotal, dfdengue[,c("CD_GEOCODI","TDS","TAS","TCS")], by="CD_GEOCODI", all.x=T)
summary(zktotal$TDS)

zktotal$lntxdeng <- ifelse(zktotal$TDS > 0, log(zktotal$TDS+1), 0)
zktotal$lntxdeng <- ifelse(is.na(zktotal$lntxdeng), 0, zktotal$lntxdeng)
summary(zktotal$lntxdeng)

#Tabela 1
summary(zktotal@data[,c("p_0a1sm","p_parpre","pr_analf","p_agua","p_esgoto","p_lixo","eq_esf15")])
zktotal@data %>% select(p_0a1sm,p_parpre,pr_analf,p_agua,p_esgoto,p_lixo,eq_esf15) %>% descr()

setwd("E:/IESC/Doutorado/analises/objetivo2/gwr")

tmap_options(check.and.fix = TRUE)

#análise do modelo GWR do Python
dfres_gwr12b <- st_read("zkgwr_setor12b.gpkg") #12b foi o modelo escolhido
names(dfres_gwr12b)
summary(dfres_gwr12b)

#tabela 1 do artigo
ax <- as.data.frame(dfres_gwr12b) %>% 
  dplyr::select(-geom) %>% 
  dplyr::select(lntxzk,p_0a1sm,p_parpr,pr_anlf,p_agua,p_esgot,p_lixo,eq_sf15,lntxdng) 

tab_corr(ax, corr.method="spearman", triangle = "upper")

dfres_gwr12b <- st_make_valid(dfres_gwr12b)
st_crs(dfres_gwr12b) <- 4674

boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$v0a1sm2)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$vagua2)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$vesgot2)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$vlixo2)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$eq_sf15)

#vegetal <- readOGR("E:/IESC/Doutorado/bases/shapes/cobertura_vegetal/sirgas2000_Cobertura_Vegetal.shp", 
#                   "sirgas2000_Cobertura_Vegetal", encoding="UTF-8")
vegetal <- st_read("E:/IESC/Doutorado/bases/shapes/cobertura_vegetal/sirgas2000_Cobertura_Vegetal.shp")

#classificação: 0=outros; 1=Mata; 2=água
tab1(vegetal@data$CLASSE, graph=F)
vegetal$novaclasse <- ifelse(substr(vegetal$CLASSE, 1, 5)=="Flore" | 
                                    substr(vegetal$CLASSE, 1, 5)=="Reflo", 1,
                                  ifelse(substr(vegetal$CLASSE, 1, 5)=="Corpo", 2, 0))
vegetal$novaclasse <- factor(vegetal$novaclasse, levels=2:0, 
                                  labels=c("Lagoon","Forest","Not information"))
tab1(vegetal$novaclasse, graph=F)

#para o not information não aparecer no mapa
vegetal$novaclasse2 <- ifelse(substr(vegetal$CLASSE, 1, 5)=="Flore" | 
                                    substr(vegetal$CLASSE, 1, 5)=="Reflo", 1,
                                  ifelse(substr(vegetal$CLASSE, 1, 5)=="Corpo", 2, NA))
vegetal$novaclasse2 <- factor(vegetal$novaclasse2, levels=2:1, 
                                  labels=c("Lagoon","Forest"))
tab1(vegetal$novaclasse2, graph=F)

st_crs(vegetal) <- 4674

#vegetal <- st_make_valid(vegetal)

dfmun <- st_read("E:/IESC/Doutorado/bases/shapes/munic_rj/novomun.shp")
#rj_osm <- read_osm(bb(dfmun, ext=1.5, projection ="longlat"), type="stamen-toner")
rj_osm <- read_osm(bb(dfmun, ext=1.5)) #projection ="longlat"))
st_crs(dfmun) <- 4674

nbsetor <- read.gal("E:/IESC/Doutorado/bases/shapes/setor_cens/zktotal_setor.gal", region.id=dfres_gwr12b$CD_GEOCODI)

mod_ols_mgwr <- lm(lntxzk ~ p_0a1sm+vagua2+vesgot2+eq_sf15+lntxdng, data=dfres_gwr12b) #sem padronizar
summary(mod_ols_mgwr)
moran.mc(residuals(mod_ols_mgwr), nb2listw(nbsetor, style="W"), nsim=999)

mod_sar_mgwr <- lagsarlm(lntxzk ~ v0a1sm2+vagua2+vesgot2+eq_sf15+lntxdng, 
                     data=dfres_gwr12b, nb2listw(nbsetor, style="W"), method="Matrix", quiet=T)
summary(mod_sar_mgwr, Nagelkerke=T)

mod_car_mgwr <- errorsarlm(lntxzk ~ v0a1sm2+vagua2+vesgot2+eq_sf15+lntxdng, 
                         data=dfres_gwr12b, nb2listw(nbsetor, style="W"), method="Matrix", quiet=T)
summary(mod_car_mgwr, Nagelkerke=T)

boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$vagua)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$vesgoto)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$vlixo)
boxplot(dfres_gwr12b$lntxzk ~ dfres_gwr12b$eq_sf15)

summary(dfres_gwr12b$mu)
summary(dfres_gwr12b$lntxzk)
dfres_gwr12b$novoresid <- dfres_gwr12b$lntxzk - dfres_gwr12b$mu

#Moran dos resíduos do modelo GWR
moran.mc(dfres_gwr12b$novoresid, nb2listw(nbsetor, style="W"), nsim=999) #0.09

moran.plot(dfres_gwr12b$novoresid, nb2listw(nbsetor, style="W"), labels=F, xlab="Resíduos",ylab="MGWR", 
           cex=0.8, spChk=F, ylim=c(-2,2), cex.axis=0.8, pch=19)

map_nz <- tm_shape(vegetal) + 
  tm_polygons("novaclasse2", border.alpha=0, title="Classes", palette=c("#00BFFF","#5F8c04","#BFBFBF"), 
              showNA=F, colorNA="#B8B8B8", legend.show=T) +
  tm_shape(dfres_gwr12b) +
  tm_polygons(col="novoresid",border.alpha=0, style="quantile",n=3, title = "MGWR Residuals", 
              palette=c("#FFF0A8","#FEAA38","#B74202"),showNA=F, colorNA="gray80") + 
  tm_shape(dfmun) +
  tm_borders(lwd=1) + 
  tm_legend(position = c("left", "top")) + 
  tm_scale_bar(position = c("right", "bottom"), text.size=0.7) + 
  tm_compass(position = c("right", "bottom"), size=0.95, text.size=0.7)
tmap_save(map_nz, "E:/IESC/Doutorado/projeto/Artigo2/nova_revista/Figure2.tiff", dpi=600) 

#filter t_vals do MGWR considera #default behavior using corrected alpha (valores = 0 são não significativos)
summary(dfres_gwr12b$t_v0a1sm2)
length(dfres_gwr12b$t_p0a1sm)
summary(abs(dfres_gwr12b$t_v0a1sm2[dfres_gwr12b$t_v0a1sm2 != 0]))
length(dfres_gwr12b$t_v0a1sm2[dfres_gwr12b$t_v0a1sm2 != 0])   #7233
length(dfres_gwr12b$t_vagua2[dfres_gwr12b$t_vagua2 != 0])     #834
length(dfres_gwr12b$t_vesgoto2[dfres_gwr12b$t_vesgoto2 != 0]) #733
length(dfres_gwr12b$t_eqesf15[dfres_gwr12b$t_eqesf15 != 0]) #1485

length(dfres_gwr12b$tcomum_p0a[dfres_gwr12b$tcomum_p0a != 0])   #2840
length(dfres_gwr12b$tcomum_pag[dfres_gwr12b$tcomum_pag != 0])     #2319
length(dfres_gwr12b$tcomum_pes[dfres_gwr12b$tcomum_pes != 0]) #3560
length(dfres_gwr12b$tcomum_eqe[dfres_gwr12b$tcomum_eqe != 0]) #6717

head(dfres_gwr12b[dfres_gwr12b$t_p0a1sm!=0, c("gwr_p0a1sm","t_p0a1sm")])

#avaliar se o coeficiente é significativo ou não de acordo com o t_
#filter t_vals do MGWR considera #default behavior using corrected alpha (valores = 0 são não significativos)

head(dfres_gwr12b[dfres_gwr12b$t_p0a1sm!=0, c("gwr_p0a1sm","t_p0a1sm")])

pal <- c("#FECC5C","#BD0026")

dfres_gwr12b$beta_p0a1sm <- ifelse(dfres_gwr12b$tcomum_p0a1sm != 0, dfres_gwr12b$gwr_p0a1sm, NA)
summary(dfres_gwr12b$beta_p0a1sm)
#summary(dfres_gwr12b$gwr_v0a1sm)
head(dfres_gwr12b[dfres_gwr12b$t_p0a1sm==0, c("beta_p0a1sm","gwr_p0a1sm","t_p0a1sm")])
tab1(cut(dfres_gwr12b$beta_p0a1sm, breaks=c(-5,0,5)), graph=F) #3 categorias

dfres_gwr12b$beta_pagua <- ifelse(dfres_gwr12b$tcomum_vagua2 != 0, dfres_gwr12b$gwr_vagua2, NA)
summary(dfres_gwr12b$beta_pagua)

dfres_gwr12b$beta_pesgoto <- ifelse(dfres_gwr12b$tcomum_vesgoto2 != 0, dfres_gwr12b$gwr_vesgot, NA)
summary(dfres_gwr12b$beta_pesgoto)

dfres_gwr12b$beta_eqesf15 <- ifelse(dfres_gwr12b$tcomum_eqesf15 != 0, dfres_gwr12b$gwr_eqesf15, NA)
summary(dfres_gwr12b$beta_eqesf15)

dfres_gwr12b$beta_lntxdng <- ifelse(dfres_gwr12b$tcomum_txdengue != 0, dfres_gwr12b$gwr_txdengue, NA)
summary(dfres_gwr12b$beta_lntxdng)

#setores com beta negativo para dengue no artigo
length(dfres_gwr12b$beta_lntxdng[dfres_gwr12b$beta_lntxdng < 0 & !is.na(dfres_gwr12b$beta_lntxdng)])

teste_nz1 <- tm_shape(vegetal) + 
  tm_polygons("novaclasse2", border.alpha=0, title="Classes", palette=c("#00BFFF","#5F8c04","#BFBFBF"), 
              showNA=F, colorNA="#B8B8B8", legend.show=T) +
  tm_shape(dfres_gwr12b) +
  tm_polygons(col="beta_p0a1sm",border.alpha = 0, breaks=c(-2,0,2),title="INCOME<1MW",
              palette=pal, showNA=TRUE, colorNA="#B8B8B8", textNA="not significant") + 
  tm_shape(dfmun) +
  tm_borders(lwd=1) + 
  tm_legend(position = c("left", "bottom")) + 
  tm_compass(position = c("right", "bottom"), size=0.6, text.size=0.5) + 
  tm_scale_bar(position = c("right", "bottom"), text.size=0.6) + 
  tm_layout(title="A", legend.title.size = 0.6, legend.text.size = 0.5, inner.margins=0.09)

teste_nz2 <- tm_shape(vegetal) +
  tm_polygons("novaclasse2", border.alpha=0, title="Classes", palette=c("#00BFFF","#5F8c04","#BFBFBF"), 
              showNA=F, colorNA="#B8B8B8", legend.show=T) +
  tm_shape(dfres_gwr12b) +
  tm_polygons(col="beta_pagua",border.alpha = 0, breaks=c(-30,0,8),title="RUNNING_WATER",palette=pal, showNA=TRUE,
              colorNA="#B8B8B8", textNA="not significant") + 
  tm_shape(dfmun) +
  tm_borders(lwd=1) + 
  tm_legend(position = c("left", "bottom")) + 
  tm_compass(position = c("right", "bottom"), size=0.6, text.size=0.5) + 
  tm_scale_bar(position = c("right", "bottom"), text.size=0.6) + 
  tm_layout(title="B", legend.title.size = 0.6, legend.text.size = 0.5, inner.margins=0.09)

teste_nz3 <- tm_shape(vegetal) +
  tm_polygons("novaclasse2", border.alpha=0, title="Classes", palette=c("#00BFFF","#5F8c04","#BFBFBF"), 
              showNA=F, colorNA="#B8B8B8", legend.show=T) +
  tm_shape(dfres_gwr12b) +
  tm_polygons(col="beta_pesgoto",border.alpha = 0, breaks=c(-4,0,8),title="SEWER_NETWORK",palette=pal, showNA=TRUE,
              colorNA="#B8B8B8", textNA="not significant") + 
  tm_shape(dfmun) +
  tm_borders(lwd=1) + 
  tm_legend(position = c("left", "bottom")) + 
  tm_compass(position = c("right", "bottom"), size=0.6, text.size=0.5) + 
  tm_scale_bar(position = c("right", "bottom"), text.size=0.6) + 
  tm_layout(title="C", legend.title.size = 0.6, legend.text.size = 0.5, inner.margins=0.09)

teste_nz4 <- tm_shape(vegetal) +
  tm_polygons("novaclasse2", border.alpha=0, title="Classes", palette=c("#00BFFF","#5F8c04","#BFBFBF"), 
              showNA=F, colorNA="#B8B8B8", legend.show=T) +
  tm_shape(dfres_gwr12b) +
  tm_polygons(col="beta_eqesf15",border.alpha = 0, breaks=c(-5,0,5),title="TEAM_FHS",palette=pal, showNA=TRUE,
              colorNA="#B8B8B8", textNA="not significant") + 
  tm_shape(dfmun) +
  tm_borders(lwd=1) + 
  tm_legend(position = c("left", "bottom")) + 
  tm_compass(position = c("right", "bottom"), size=0.6, text.size=0.5) + 
  tm_scale_bar(position = c("right", "bottom"), text.size=0.6) + 
  tm_layout(title="D", legend.title.size = 0.6, legend.text.size = 0.5, inner.margins=0.09)

teste_nz5 <- tm_shape(vegetal) +
  tm_polygons("novaclasse2", border.alpha=0, title="Classes", palette=c("#00BFFF","#5F8c04","#BFBFBF"), 
              showNA=F, colorNA="#B8B8B8", legend.show=T) +
  tm_shape(dfres_gwr12b) +
  tm_polygons(col="beta_lntxdng",border.alpha = 0, breaks=c(-1,0,2),title="DENV_RATE",palette=pal, showNA=TRUE,
              colorNA="#B8B8B8", textNA="not significant") + 
  tm_shape(dfmun) +
  tm_borders(lwd=1) + 
  tm_legend(position = c("left", "bottom")) + 
  tm_compass(position = c("right", "bottom"), size=0.6, text.size=0.5) + 
  tm_scale_bar(position = c("right", "bottom"), text.size=0.6) + 
  tm_layout(title="E", legend.title.size = 0.6, legend.text.size = 0.5, inner.margins=0.09)

tm1 <- tmap_arrange(teste_nz1, teste_nz2, teste_nz3, teste_nz4, teste_nz5)
tmap_save(tm1, "E:/IESC/Doutorado/projeto/Artigo2/plos_ntd/Figure3_all2.tiff", dpi=350) 


