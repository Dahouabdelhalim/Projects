## This is the code for the manuscript "Habitat loss shapes the structure,
##stability and the role of species in tropical mutualistic networks"
##By Fernando César Gonçalves Bonfim, Pavel Dodonov and Eliana Cazetta##
##E-mail: fernandouesb@gmail.com##

##Load your library

setwd("C:\\\\Users....")

library(spdep)
library(MuMIn)
library(MASS)
library(segmented)
library(mgcv)
library(ape)
library(nlme)
library(nls2)
library(bbmle)
library(PerformanceAnalytics)


dados<-read.table("analise.csv", header=T, sep=",")
str(dados)
dados$coords<-as.matrix(cbind(dados$lon,dados$lat))


##Test if sampling effor affect response variables: GAM models
#Robustness
robgam<-gam(robustness ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(robgam)##Not significant
plot(resid(robgam)) 

#Number of birds (NB)
nbgam<-gam(NB ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=poisson, na.action="na.fail")
summary(nbgam)##significant
plot(resid(nbgam)) 

#Number of plants(NP)
npgam<-gam(NP ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=poisson, na.action="na.fail")
summary(npgam)##significant
plot(resid(npgam)) 

#Number of interactions (NI)
nigam<-gam(NI ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=poisson, na.action="na.fail")
summary(nigam)##significant
plot(resid(nigam)) 

#Diversity of interactions (shanon)
divgam<-gam(shanon ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(divgam)##significant
plot(resid(divgam)) 

#Conectance z score (zconectancia)
congam<-gam(z_conectancia ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(congam)## significant
plot(resid(congam)) 

##Number of links by species- zscore (z_link)
linkgam<-gam(z_link ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(linkgam)## significant
plot(resid(linkgam)) 

##Nestedness- zscore (zNODF)
nodgam<-gam(zNODF ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(nodgam)##significant
plot(resid(nodgam)) 

#Body mass( mean)
bodygam<-gam(body_mass ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(bodygam)##Não significativo
plot(resid(bodygam)) 

#Bill width
billgam<-gam(bill_width ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(billgam)##Not significant
plot(resid(billgam)) 

##Seed diameter- diameter
seedgam<-gam(diameter ~ s(sam_eff, k= 3, fx = FALSE), data=dados , family=gaussian, na.action="na.fail")
summary(seedgam)## significant
plot(resid(seedgam)) 

##Test the correlation between forest cover and sampling effort
dados.pred <- dados[,5:21]
str(dados.pred)
chart.Correlation(dados.pred$sam_eff~dados.pred,  method = c("pearson"))

##################################
#First, we creat an object with all coverage.
dados.Fc <- dados[,5:20]
str(dados.Fc)

##Robustness
dados$resid.rob<-resid(robgam)

##Linear model
### Create a list to put the results of the model selection and shapiro test.
result.rob <- list()
shapiro.rob <- list()

### Start with a loop

for(i in 1:ncol(dados.Fc)) { ### Will repeat what the command is defining i=1, i=2, ..., i=15
  y.temp <- dados.Fc[,i] ### Catch the 1a, 2a..., 15a column
  name.temp <- names(dados.Fc)[i] ### Cathc the name of the same column
  mod.temp <- lm(dados$resid.rob ~ y.temp, na.action="na.fail") ### adjust the model, for the column selected as y.temp
  result.rob[[i]] <- mod.temp ### Put the result in a list
  names(result.rob)[i] <- name.temp ### Give the appropriate name
  shapiro.rob[[i]] <- shapiro.test(resid(mod.temp)) ### Put the Shapiro's result
  names(shapiro.rob)[i] <- name.temp ### Give the appropriate name
  print(i) # Put on the screem how many iteraction was done
}

model.sel(result.rob)


##Quadratic model

result.robq <- list()
shapiro.robq <- list()

for(i in 1:ncol(dados.Fc)) { 
  y.temp <- dados.Fc[,i] 
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.rob ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.robq[[i]] <- mod.temp 
  names(result.robq)[i] <- name.temp 
  shapiro.robq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.robq)[i] <- name.temp
  print(i)
}
model.sel(result.robq)


## Power law##
##First we create a function to with the power function
power<-function(a,x,b){
  a*x^b
}

##Than, we create a grid to search for initial values
grid.pow<-expand.grid(list(a = seq(0,100, by = 0.5),
                              b = seq(-1,1, by = 0.1)))

result.robp <- list()
shapiro.robp <- list()

##Some models do not converged
for(i in 3:ncol(dados.Fc)) { 
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  inicial.pow= nls2(resid.rob ~power(y.temp, a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.rob ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail") ### ajusta o modelo, pra coluna selecionada como y.temp
  result.robp[[i]] <- mod.temp 
  names(result.robp)[i] <- name.temp 
  shapiro.robp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.robp)[i] <- name.temp
  print(i) 
}

AICctab(result.robp[["Fc_700"]],result.robp[["Fc_800"]], result.robp[["Fc_900"]],
                result.robp[["Fc_1000"]], result.robp[["Fc_1100"]], result.robp[["Fc_1200"]],
                result.robp[["Fc_1300"]], result.robp[["Fc_1400"]], result.robp[["Fc_1500"]], 
                result.robp[["Fc_1600"]], result.robp[["Fc_1700"]], result.robp[["Fc_1800"]],
                result.robp[["Fc_1900"]], result.robp[["Fc_2000"]], nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.robpi <- list()
shapiro.robpi <- list()

resid.rob<-dados$resid.rob

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.rob ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.robpi[[i]] <- mod.temp.2
  names(result.robpi)[i] <- name.temp
  shapiro.robpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.robpi)[i] <- name.temp
  print(i) 
}

model.sel(result.robpi)

##Null model

null<- lm(resid.rob~1, data=dados)

##Test if the models have auto spatial correlation.
rob_500<-result.rob[["Fc_500"]]
robq_600<-result.robq[["Fc_600"]]
robl_600<-result.robl[["Fc_600"]]
robp_1600<-result.robp[["Fc_1600"]]
robpi_1900<-result.robpi[["Fc_1900"]]
robg_500<-result.robg[["Fc_500"]]

##To calculate auto spatial correlation
dists<-as.matrix(dist(cbind(dados$lon, dados$lat)))
dists.inv<-1/dists
diag(dists.inv)<-0
dists.inv[is.infinite(dists.inv)] <- 0

Moran.I(residuals(rob_500),dists.inv)
Moran.I(residuals(robq_600),dists.inv)
Moran.I(residuals(robp_1600),dists.inv)
Moran.I(residuals(robg_500),dists.inv)
Moran.I(residuals(robpi_1900),dists.inv)

AICctab(rob_500,robq_600, robp_1600, robpi_1900, null, nobs=25, base=T, weights=T, logLik=T)
summary(rob_500)

#####################################
##Number of birds
dados$resid.nb<- resid(nbgam)
##Linear model

result.nb <- list()
shapiro.nb <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.nb ~ y.temp, na.action="na.fail")
  result.nb[[i]] <- mod.temp
  names(result.nb)[i] <- name.temp
  shapiro.nb[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.nb)[i] <- name.temp
  print(i)
}

model.sel(result.nb)


##Quadratic model

result.nbq <- list()
shapiro.nbq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.nb ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.nbq[[i]] <- mod.temp
  names(result.nbq)[i] <- name.temp
  shapiro.nbq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.nbq)[i] <- name.temp
  print(i)
}

model.sel(result.nbq)



## Power law##
grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-1,1, by = 0.1)))
result.nbp <- list()
shapiro.nbp <- list()

##Models do not converged

for(i in 1:ncol(dados.Fc)) { ### Vai repetir o que os comandos definindo i=1, i=2, ..., i=15
  y.temp <- dados.Fc[,i] ### Pega a 1a, 2a..., 15a coluna
  name.temp <- names(dados.Fc)[i] ### Pega o nome da mesma coluna
  initial.pow= nls2(resid.np ~power(dados.Fc[,i], a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.np ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail") ### ajusta o modelo, pra coluna selecionada como y.temp
  result.nbp[[i]] <- mod.temp ### coloca o resultado na lista
  names(result.nbp)[i] <- name.temp ### dÃ¡ o nome apropriado
  shapiro.nbp[[i]] <- shapiro.test(resid(mod.temp)) ### coloca o resultado do shapiro
  names(shapiro.nbp)[i] <- name.temp ### dÃ¡ o nome apropriado
  print(i) # coloca na tela quantas iteraÃ§Ãµes jÃ¡ foram (talvez nÃ£o funcione em RStudio)
}

AICctab(result.robp, nobs=25, base=T, weights=T, logLik=T)


##Piecewise model

result.nbpi <- list()
shapiro.nbpi <- list()

resid.nb<-dados$resid.nb

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.nb ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.nbpi[[i]] <- mod.temp.2
  names(result.nbpi)[i] <- name.temp
  shapiro.nbpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.nbpi)[i] <- name.temp
  print(i)
}

model.sel(result.nbpi)


##Null model

null<- lm(resid.nb~1, data=dados)

##Test if the models have auto spatial correlation.
nb_500<-result.nb[["Fc_500"]]
nbq_500<-result.nbq[["Fc_500"]]
nbl_500<-result.nbl[["Fc_500"]]
##Power law models do not converged
nbpi_500<-result.nbpi[["Fc_500"]]
nbg_500<-result.nbg[["Fc_500"]]


dists<-as.matrix(dist(cbind(dados$lon, dados$lat)))
dists.inv<-1/dists
diag(dists.inv)<-0
dists.inv[is.infinite(dists.inv)] <- 0

Moran.I(residuals(nb_500),dists.inv)
Moran.I(residuals(nbq_500),dists.inv)
Moran.I(residuals(nbg_500),dists.inv)
Moran.I(residuals(nbpi_500),dists.inv)
plot(nbpi_500)

AICctab(nb_500,nbq_500, nbpi_500, null, nobs=25, base=T, weights=T, logLik=T)
summary(nbpi_500)


###########################################
##Number of plants
dados$resid.np<- resid(npgam)
##Linear model

### Create a list to put the results of the model selection and shapiro test.
result.np <- list()
shapiro.np <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.np ~ y.temp, na.action="na.fail")
  result.np[[i]] <- mod.temp
  names(result.np)[i] <- name.temp
  shapiro.np[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.np)[i] <- name.temp
  print(i)
}

model.sel(result.np)


##Quadratic model

result.npq <- list()
shapiro.npq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.np ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.npq[[i]] <- mod.temp
  names(result.npq)[i] <- name.temp
  shapiro.npq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.npq)[i] <- name.temp
  print(i)
}

model.sel(result.npq)



## Power law##
grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-1,1, by = 0.1)))
result.npp <- list()
shapiro.npp <- list()

##Models do not converged

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.np ~power(dados.Fc[,i], a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.np ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")
  result.npp[[i]] <- mod.temp
  names(result.npp)[i] <- name.temp
  shapiro.npp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.npp)[i] <- name.temp
  print(i)
}

AICctab(result.npp, nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.nppi <- list()
shapiro.nppi <- list()

resid.np<-dados$resid.np

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.np ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.nppi[[i]] <- mod.temp.2
  names(result.nppi)[i] <- name.temp
  shapiro.nppi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.nppi)[i] <- name.temp
  print(i)
}

model.sel(result.nppi)


##Null model

null<- lm(resid.np~1, data=dados)

##Test if the models have auto spatial correlation.
np_1700<-result.np[["Fc_1700"]]
npq_1500<-result.npq[["Fc_1500"]]
npl_1000<-result.npl[["Fc_1000"]]
##Power law models do not converged
nppi_1500<-result.nppi[["Fc_1500"]]
npg_1700<-result.npg[["Fc_1700"]]

##To calculate auto spatial correlation
Moran.I(residuals(np_1700),dists.inv)
Moran.I(residuals(npq_1500),dists.inv)
Moran.I(residuals(npg_1700),dists.inv)
Moran.I(residuals(nppi_1500),dists.inv)
plot(np_1700)

AICctab(np_1700,npq_1500, nppi_1500, null, nobs=25, base=T, weights=T, logLik=T)
summary(npq_1500)

###########################################
##Number of interactions
dados$resid.ni<- resid(nigam)
##Linear model

result.ni <- list()
shapiro.ni <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.ni ~ y.temp, na.action="na.fail")
  result.ni[[i]] <- mod.temp
  names(result.ni)[i] <- name.temp
  shapiro.ni[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.ni)[i] <- name.temp
  print(i)
}

model.sel(result.ni)


##Quadratic model

result.niq <- list()
shapiro.niq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.ni ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.niq[[i]] <- mod.temp
  names(result.niq)[i] <- name.temp
  shapiro.niq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.niq)[i] <- name.temp
  print(i)
}

model.sel(result.niq)


## Power law##
grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-1,1, by = 0.1)))

result.nip <- list()
shapiro.nip <- list()

##Some models do not converged

for(i in 3:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.ni ~power(dados.Fc[,i], a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.ni ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")
  result.nip[[i]] <- mod.temp
  names(result.nip)[i] <- name.temp
  shapiro.nip[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.nip)[i] <- name.temp
  print(i)
}

AICctab(result.nip[["Fc_700"]],result.nip[["Fc_800"]], result.nip[["Fc_900"]],
        result.nip[["Fc_1000"]], result.nip[["Fc_1100"]], result.nip[["Fc_1200"]],
        result.nip[["Fc_1300"]], result.nip[["Fc_1400"]], nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.nipi <- list()
shapiro.nipi <- list()

resid.ni<-dados$resid.ni

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.ni ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.nipi[[i]] <- mod.temp.2
  names(result.nipi)[i] <- name.temp
  shapiro.nipi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.nipi)[i] <- name.temp
  print(i)
}

model.sel(result.nipi)


##Null model

null<- lm(resid.ni~1, data=dados)

##Test if the models have auto spatial correlation.
ni_1600<-result.ni[["Fc_1600"]]
niq_600<-result.niq[["Fc_600"]]
nil_500<-result.nil[["Fc_500"]]
nip_500<-result.nip[["Fc_700"]]
nipi_1500<-result.nipi[["Fc_1500"]]
nig_600<-result.nig[["Fc_600"]]

##To calculate auto spatial correlation
Moran.I(residuals(ni_1600),dists.inv)
Moran.I(residuals(niq_600),dists.inv)
Moran.I(residuals(nip_700),dists.inv)
Moran.I(residuals(nig_600),dists.inv)
Moran.I(residuals(nipi_1500),dists.inv)
plot(ni_1600)

AICctab(ni_1600,niq_600, nip_700, nipi_1500, null, nobs=25, base=T, weights=T, logLik=T)
summary(niq_600)

###########################################
##Conectance
dados$resid.con<- resid(congam)

##Linear model

result.con <- list()
shapiro.con <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.con ~ y.temp, na.action="na.fail")
  result.con[[i]] <- mod.temp
  names(result.con)[i] <- name.temp
  shapiro.con[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.con)[i] <- name.temp
  print(i)
}

model.sel(result.con)


##Quadratic model

result.conq <- list()
shapiro.conq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.con ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.conq[[i]] <- mod.temp
  names(result.conq)[i] <- name.temp
  shapiro.conq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.conq)[i] <- name.temp
  print(i)
}

model.sel(result.conq)


## Power law##

grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-2,2, by = 0.2)))

result.conp <- list()
shapiro.conp <- list()

##Some models do not converged

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.con ~power(y.temp, a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.con ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")  
  result.conp[[i]] <- mod.temp
  names(result.conp)[i] <- name.temp
  shapiro.conp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.conp)[i] <- name.temp
  print(i)
}

AICctab(result.conp[["Fc_700"]], result.conp[["Fc_800"]], result.conp[["Fc_900"]],
        result.conp[["Fc_1000"]], result.conp[["Fc_1100"]], result.conp[["Fc_1200"]],
        result.conp[["Fc_1300"]], result.conp[["Fc_1400"]], result.conp[["Fc_1500"]],
        nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.conpi <- list()
shapiro.conpi <- list()

resid.con<-dados$resid.con

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.con ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.conpi[[i]] <- mod.temp.2
  names(result.conpi)[i] <- name.temp
  shapiro.conpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.conpi)[i] <- name.temp
  print(i)
}

model.sel(result.conpi)


##Null model

null<- lm(resid.con~1, data=dados)

##Test if the models have auto spatial correlation.
con_500<-result.con[["Fc_500"]]
conq_500<-result.conq[["Fc_500"]]
conl_1000<-result.conl[["Fc_1000"]]
conp_800<-result.conp[["Fc_800"]]
cong_500<-result.cong[["Fc_500"]]
conpi_1600<-result.conpi[["Fc_1600"]]

##To calculate auto spatial correlation
Moran.I(residuals(con_500),dists.inv)
Moran.I(residuals(conq_500),dists.inv)
Moran.I(residuals(conp_800),dists.inv)
Moran.I(residuals(cong_500),dists.inv)
Moran.I(residuals(conpi_1600),dists.inv)
plot(conpi_1300)

AICctab(con_500,conq_500, conp_800, conpi_1600, null, nobs=25, base=T, weights=T, logLik=T)
summary(con_500)

###########################################
##Number of links per species
dados$resid.link<- resid(linkgam)

##Linear model

result.link <- list()
shapiro.link <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.link ~ y.temp, na.action="na.fail")
  result.link[[i]] <- mod.temp
  names(result.link)[i] <- name.temp
  shapiro.link[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.link)[i] <- name.temp
  print(i)
}

model.sel(result.link)


##Quadratic model

result.linkq <- list()
shapiro.linkq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.link ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.linkq[[i]] <- mod.temp
  names(result.linkq)[i] <- name.temp
  shapiro.linkq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.linkq)[i] <- name.temp
  print(i)
}

model.sel(result.linkq)

## Power law##

grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-2,2, by = 0.2)))

result.linkp <- list()
shapiro.linkp <- list()

##Some models do not converged

for(i in 3:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.link ~power(y.temp, a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.link ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")  
  result.linkp[[i]] <- mod.temp
  names(result.linkp)[i] <- name.temp
  shapiro.linkp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.linkp)[i] <- name.temp
  print(i)
}

AICctab(result.linkp[["Fc_700"]], result.linkp[["Fc_800"]], result.linkp[["Fc_900"]],
        result.linkp[["Fc_1000"]], result.linkp[["Fc_1100"]], result.linkp[["Fc_1200"]],
        result.linkp[["Fc_1300"]], result.linkp[["Fc_1400"]], result.linkp[["Fc_1500"]],
        result.linkp[["Fc_1600"]], result.linkp[["Fc_1700"]], nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.linkpi <- list()
shapiro.linkpi <- list()

resid.link<-dados$resid.link

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.link ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.linkpi[[i]] <- mod.temp.2
  names(result.linkpi)[i] <- name.temp
  shapiro.linkpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.linkpi)[i] <- name.temp
  print(i)
}
model.sel(result.linkpi)

##Null model

null<- lm(resid.link~1, data=dados)

##Test if the models have auto spatial correlation.
link_500<-result.link[["Fc_500"]]
linkq_500<-result.linkq[["Fc_500"]]
linkl_1000<-result.linkl[["Fc_1000"]]
linkp_700<-result.linkp[["Fc_700"]]
linkg_500<-result.linkg[["Fc_500"]]
linkpi_1600<-result.linkpi[["Fc_1600"]]

##To calculate auto spatial correlation
Moran.I(residuals(link_500),dists.inv)
Moran.I(residuals(linkq_500),dists.inv)
Moran.I(residuals(linkp_700),dists.inv)
Moran.I(residuals(linkg_500),dists.inv)
Moran.I(residuals(linkpi_1600),dists.inv)
plot(linkpi_1600)

AICctab(link_500,linkq_500, linkp_700, linkpi_1600, null, nobs=25, base=T, weights=T, logLik=T)
summary(link_500)

###########################################
##Nestedness
dados$resid.nodf<- resid(nodgam)

##Linear model

result.nodf <- list()
shapiro.nodf <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.nodf ~ y.temp, na.action="na.fail")
  result.nodf[[i]] <- mod.temp
  names(result.nodf)[i] <- name.temp
  shapiro.nodf[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.nodf)[i] <- name.temp
  print(i)
}
model.sel(result.nodf)

##Quadratic model

result.nodfq <- list()
shapiro.nodfq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.nodf ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.nodfq[[i]] <- mod.temp
  names(result.nodfq)[i] <- name.temp
  shapiro.nodfq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.nodfq)[i] <- name.temp
  print(i)
}
model.sel(result.nodfq)

## Power law##

grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-3,2, by = 0.1)))

result.nodfp <- list()
shapiro.nodfp <- list()

##Models do not converged

for(i in 3:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.nodf ~power(y.temp, a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.nodf ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")  
  result.nodfp[[i]] <- mod.temp
  names(result.nodfp)[i] <- name.temp
  shapiro.nodfp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.nodfp)[i] <- name.temp
  print(i)
}
model.sel(result.nodfp)

##Piecewise model

result.nodfpi <- list()
shapiro.nodfpi <- list()

resid.nodf<-dados$resid.nodf

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.nodf ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.nodfpi[[i]] <- mod.temp.2
  names(result.nodfpi)[i] <- name.temp
  shapiro.nodfpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.nodfpi)[i] <- name.temp
  print(i)
}

model.sel(result.nodfpi)

##Null model

null<- lm(resid.nodf~1, data=dados)

##Test if the models have auto spatial correlation.
nodf_500<-result.nodf[["Fc_500"]]
nodfq_600<-result.nodfq[["Fc_600"]]
nodfl_1500<-result.nodfl[["Fc_1500"]]
#Power law models do not converged
nodfg_500<-result.nodfg[["Fc_500"]]
nodfpi_500<-result.nodfpi[["Fc_500"]]

##To calculate auto spatial correlation
Moran.I(residuals(nodf_500),dists.inv)
Moran.I(residuals(nodfq_600),dists.inv)
Moran.I(residuals(nodfg_500),dists.inv)
Moran.I(residuals(nodfpi_500),dists.inv)
plot(nodfpi_500)

AICctab(nodf_500,nodfq_600, nodfpi_500, null, nobs=25, base=T, weights=T, logLik=T)

summary(nodf_500)

###########################################
##Bill width
dados$resid.bill<- resid(billgam)

##Linear model

result.bill <- list()
shapiro.bill <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.bill ~ y.temp, na.action="na.fail")
  result.bill[[i]] <- mod.temp
  names(result.bill)[i] <- name.temp
  shapiro.bill[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.bill)[i] <- name.temp
  print(i)
}

model.sel(result.bill)

##Quadratic model

result.billq <- list()
shapiro.billq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.bill ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.billq[[i]] <- mod.temp
  names(result.billq)[i] <- name.temp
  shapiro.billq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.billq)[i] <- name.temp
  print(i)
}

model.sel(result.billq)

## Power law##

grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-5,5, by = 0.5)))

result.billp <- list()
shapiro.billp <- list()

##Some models do not converged

for(i in 3:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.bill ~power(y.temp, a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.bill ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")  
  result.billp[[i]] <- mod.temp
  names(result.billp)[i] <- name.temp
  shapiro.billp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.billp)[i] <- name.temp
  print(i)
}

bill_pow<-AICctab(result.billp[["Fc_700"]], result.billp[["Fc_800"]], result.billp[["Fc_900"]],
        result.billp[["Fc_1000"]], result.billp[["Fc_1100"]], result.billp[["Fc_1200"]],
        result.billp[["Fc_1300"]], result.billp[["Fc_1400"]], result.billp[["Fc_1500"]],
        result.billp[["Fc_1600"]], result.billp[["Fc_1700"]], result.billp[["Fc_1800"]],
        result.billp[["Fc_1900"]], result.billp[["Fc_2000"]], nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.billpi <- list()
shapiro.billpi <- list()

resid.bill<-dados$resid.bill

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.bill ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.billpi[[i]] <- mod.temp.2
  names(result.billpi)[i] <- name.temp
  shapiro.billpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.billpi)[i] <- name.temp
  print(i) # coloca na tela quantas iteraÃ§Ãµes jÃ¡ foram (talvez nÃ£o funcione em RStudio)
}

model.sel(result.billpi)


##Null model

null<- lm(resid.bill~1, data=dados)

##Test if the models have auto spatial correlation.
bill_2000<-result.bill[["Fc_2000"]]
billq_800<-result.billq[["Fc_800"]]
billl_1000<-result.billl[["Fc_1000"]]
billp_800<-result.billp[["Fc_800"]]
billg_800<-result.billg[["Fc_800"]]
billpi_800<-result.billpi[["Fc_800"]]

##To calculate auto spatial correlation
Moran.I(residuals(bill_2000),dists.inv)
Moran.I(residuals(billq_800),dists.inv)
Moran.I(residuals(billp_800),dists.inv)
Moran.I(residuals(billg_800),dists.inv)
Moran.I(residuals(billpi_800),dists.inv)
plot(body_2000)

AICctab(bill_2000,billq_800, billp_800, billpi_800, null, nobs=25, base=T, weights=T, logLik=T)
summary(billq_800)

###########################################
##Seed diameter
dados$resid.seed<- resid(seedgam)

##Linear model

result.seed <- list()
shapiro.seed <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.seed ~ y.temp, na.action="na.fail")
  result.seed[[i]] <- mod.temp
  names(result.seed)[i] <- name.temp
  shapiro.seed[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.seed)[i] <- name.temp
  print(i)
}

model.sel(result.seed)

##Quadratic model

result.seedq <- list()
shapiro.seedq <- list()

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp <- lm(dados$resid.seed ~ y.temp+I(y.temp^2), na.action="na.fail")
  result.seedq[[i]] <- mod.temp
  names(result.seedq)[i] <- name.temp
  shapiro.seedq[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.seedq)[i] <- name.temp
  print(i)
}

model.sel(result.seedq)


## Power law##

grid.pow<-expand.grid(list(a = seq(0,100, by = 1),
                           b = seq(-2,3, by = 0.5)))

result.seedp <- list()
shapiro.seedp <- list()

##Some models do not converged

for(i in 3:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  initial.pow= nls2(resid.seed ~power(y.temp, a, b), data=dados,
                    start=grid.pow, algorithm = "brute-force")  
  mod.temp <- gnls(resid.seed ~ a*(y.temp^b), control=gnlsControl(nlsTol=100), start=coef(initial.pow),data=dados, na.action = "na.fail")  
  result.seedp[[i]] <- mod.temp
  names(result.seedp)[i] <- name.temp
  shapiro.seedp[[i]] <- shapiro.test(resid(mod.temp))
  names(shapiro.seedp)[i] <- name.temp
  print(i)
}

AICctab(result.seedp[["Fc_700"]], result.seedp[["Fc_800"]], result.seedp[["Fc_900"]],
                result.seedp[["Fc_1000"]], result.seedp[["Fc_1100"]], result.seedp[["Fc_1200"]],
                result.seedp[["Fc_1300"]], result.seedp[["Fc_1400"]], result.seedp[["Fc_1500"]],
                result.seedp[["Fc_1600"]], result.seedp[["Fc_1700"]], result.seedp[["Fc_1800"]],
                result.seedp[["Fc_1900"]], result.seedp[["Fc_2000"]], nobs=25, base=T, weights=T, logLik=T)

##Piecewise model

result.seedpi <- list()
shapiro.seedpi <- list()

resid.seed<-dados$resid.seed

for(i in 1:ncol(dados.Fc)) {
  y.temp <- dados.Fc[,i]
  name.temp <- names(dados.Fc)[i]
  mod.temp.1<-lm(resid.seed ~ y.temp, na.action="na.fail")
  mod.temp.2<- segmented(mod.temp.1, seg.Z = ~ y.temp)
  result.seedpi[[i]] <- mod.temp.2
  names(result.seedpi)[i] <- name.temp
  shapiro.seedpi[[i]] <- shapiro.test(resid(mod.temp.2))
  names(shapiro.seedpi)[i] <- name.temp
  print(i)
}

model.sel(result.seedpi)

##Null model

null<- lm(resid.seed~1, data=dados)

##Test if the models have auto spatial correlation.
seed_500<-result.seed[["Fc_500"]]
seedq_500<-result.seedq[["Fc_500"]]
seedl_600<-result.seedl[["Fc_600"]]
seedp_1100<-result.seedp[["Fc_1100"]]
seedg_500<-result.seedg[["Fc_500"]]
seedpi_600<-result.seedpi[["Fc_600"]]

##To calculate auto spatial correlation
Moran.I(residuals(seed_500),dists.inv)
Moran.I(residuals(seedq_500),dists.inv)
Moran.I(residuals(seedp_1100),dists.inv)
Moran.I(residuals(seedg_500),dists.inv)
Moran.I(residuals(seedpi_600),dists.inv)
plot(seed_500)

AICctab(seed_500,seedq_500, seedp_1100, seedpi_600, null, nobs=25, base=T, weights=T, logLik=T)
summary(seedpi_600)
