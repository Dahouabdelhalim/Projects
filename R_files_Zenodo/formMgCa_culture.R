MgCa.dat.culture<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_culture.csv")

ruber<-MgCa.dat.culture[which(MgCa.dat.culture$Species=="G. ruber"),]
sac<-MgCa.dat.culture[which(MgCa.dat.culture$Species=="G. sacculifer"),]
bul<-MgCa.dat.culture[which(MgCa.dat.culture$Species=="G. bulloides"),]
uni<-MgCa.dat.culture[which(MgCa.dat.culture$Species=="O. universa"),]
incompta<-MgCa.dat.culture[which(MgCa.dat.culture$Species=="N. incompta"),]

#T,S, Omega
reg.ruber<-lm(log(Mg.Ca)~Temp+Salinity+log(Ω.Calcite),data=ruber)
reg.sac<-lm(log(Mg.Ca)~Temp+Salinity+log(Ω.Calcite),data=sac)
reg.bul<-lm(log(Mg.Ca)~Temp+Salinity+log(Ω.Calcite),data=bul)
reg.uni<-lm(log(Mg.Ca)~Temp+Salinity+log(Ω.Calcite),data=uni)
#reg.incompta<-lm(log(Mg.Ca)~Temp+Salinity+log(Ω.Calcite),data=incompta)

#T,S,pH
f=function(temp,salt,pH,a1,b1,b2,b3) {a1*exp(b1*temp+b2*salt)+b3*pH}
reg.ruber.pH<-nls(Mg.Ca~f(Temp,Salinity,pH,a1,b1,b2,b3),data=ruber,start=c(a1=1,b1=0.1,b2=0.1,b3=-1))
#reg.sac.pH<-nls(Mg.Ca~f(Temp,Salinity,pH,a1,b1,b2,b3),data=sac,start=c(a1=10,b1=0.01,b2=0.01,b3=-1))
reg.bul.pH<-nls(Mg.Ca~f(Temp,Salinity,pH,a1,b1,b2,b3),data=bul,start=c(a1=10,b1=0.01,b2=0.01,b3=-1))
reg.uni.pH<-nls(Mg.Ca~f(Temp,Salinity,pH,a1,b1,b2,b3),data=uni,start=c(a1=10,b1=0.01,b2=0.01,b3=-1))
#reg.incompta.pH<-nls(Mg.Ca~f(Temp,Salinity,pH,a1,b1,b2,b3),data=uni,start=c(a1=10,b1=0.01,b2=0.01,b3=-1))

ruber.T26S35<-ruber$Mg.Ca*exp(reg.ruber$coeff[2]*(26-ruber$Temp)+reg.ruber$coeff[3]*(35-ruber$Salinity))
sac.T26S35<-sac$Mg.Ca*exp(reg.sac$coeff[2]*(26-sac$Temp)+reg.sac$coeff[3]*(35-sac$Salinity))
bul.T26S35<-bul$Mg.Ca*exp(reg.bul$coeff[2]*(26-bul$Temp)+reg.bul$coeff[3]*(35-bul$Salinity))
uni.T26S35<-uni$Mg.Ca*exp(reg.uni$coeff[2]*(26-uni$Temp)+reg.uni$coeff[3]*(35-uni$Salinity))
#incompta.T26S35<-incompta$Mg.Ca*exp(reg.incompta$coeff[2]*(26-incompta$Temp)+reg.incompta$coeff[3]*(35-incompta$Salinity))

ruber.T26O5<-ruber$Mg.Ca*exp(reg.ruber$coeff[2]*(26-ruber$Temp))*(5/ruber$Ω.Calcite)^reg.ruber$coeff[4]
sac.T26O5<-sac$Mg.Ca*exp(reg.sac$coeff[2]*(26-sac$Temp))*(5/sac$Ω.Calcite)^reg.sac$coeff[4]
bul.T26O5<-bul$Mg.Ca*exp(reg.bul$coeff[2]*(26-bul$Temp))*(5/bul$Ω.Calcite)^reg.bul$coeff[4]
uni.T26O5<-uni$Mg.Ca*exp(reg.uni$coeff[2]*(26-uni$Temp))*(5/uni$Ω.Calcite)^reg.uni$coeff[4]
#incompta.T26O5<-incompta$Mg.Ca*exp(reg.incompta$coeff[2]*(26-incompta$Temp))*(5/incompta$Ω.Calcite)^reg.incompta$coeff[4]

ruber.S35O5<-ruber$Mg.Ca*exp(reg.ruber$coeff[3]*(35-ruber$Salinity))*(5/ruber$Ω.Calcite)^reg.ruber$coeff[4]
sac.S35O5<-sac$Mg.Ca*exp(reg.sac$coeff[3]*(35-sac$Salinity))*(5/sac$Ω.Calcite)^reg.sac$coeff[4]
bul.S35O5<-bul$Mg.Ca*exp(reg.bul$coeff[3]*(35-bul$Salinity))*(5/bul$Ω.Calcite)^reg.bul$coeff[4]
uni.S35O5<-uni$Mg.Ca*exp(reg.uni$coeff[3]*(35-uni$Salinity))*(5/uni$Ω.Calcite)^reg.uni$coeff[4]
#incompta.S35O5<-incompta$Mg.Ca*exp(reg.incompta$coeff[3]*(35-incompta$Salinity))*(5/incompta$Ω.Calcite)^reg.incompta$coeff[4]

ruber.T26S35.pH<-(ruber$Mg.Ca-ruber$pH*summary(reg.ruber.pH)$coeff[4])*exp(summary(reg.ruber.pH)$coeff[2]*(26-ruber$Temp)+summary(reg.ruber.pH)$coeff[3]*(35-ruber$Salinity))+ruber$pH*summary(reg.ruber.pH)$coeff[4]
bul.T26S35.pH<-(bul$Mg.Ca-bul$pH*summary(reg.bul.pH)$coeff[4])*exp(summary(reg.bul.pH)$coeff[2]*(26-bul$Temp)+summary(reg.bul.pH)$coeff[3]*(35-bul$Salinity))+bul$pH*summary(reg.bul.pH)$coeff[4]
uni.T26S35.pH<-(uni$Mg.Ca-uni$pH*summary(reg.uni.pH)$coeff[4])*exp(summary(reg.uni.pH)$coeff[2]*(26-uni$Temp)+summary(reg.uni.pH)$coeff[3]*(35-uni$Salinity))+uni$pH*summary(reg.uni.pH)$coeff[4]
incompta.T26S35.pH<-(incompta$Mg.Ca-incompta$pH*summary(reg.incompta.pH)$coeff[4])*exp(summary(reg.incompta.pH)$coeff[2]*(26-incompta$Temp)+summary(reg.incompta.pH)$coeff[3]*(35-incompta$Salinity))+incompta$pH*summary(reg.incompta.pH)$coeff[4]


dev.new(width=8, height=6);
par(mfrow=c(2,2),mar=c(4,4,1,4))

plot(ruber$Temp,ruber.S35O5,bg="gold3",las=1,xlab="temperature (ºC)",ylab="Mg/Ca @ 35 psu; Ω =5, (mmol/mol)",xlim=c(15,32),ylim=c(0,16),pch=21,col="black")
lines(seq(15,32),exp(reg.ruber$coeff[1])*exp(reg.ruber$coeff[2]*seq(15,32)+reg.ruber$coeff[3]*35)*5^reg.ruber$coeff[4],col="gold3")
points(sac$Temp,sac.S35O5,bg="palegreen3",pch=21,col="black")
lines(seq(22,32),exp(reg.sac$coeff[1])*exp(reg.sac$coeff[2]*seq(22,32)+reg.sac$coeff[3]*35)*5^reg.sac$coeff[4],col="palegreen3")
points(bul$Temp,bul.S35O5,bg="cornflowerblue",pch=21,col="black")
lines(seq(15,25),exp(reg.bul$coeff[1])*exp(reg.bul$coeff[2]*seq(15,25)+reg.bul$coeff[3]*35)*5^reg.bul$coeff[4],col="cornflowerblue")
points(uni$Temp,uni.S35O5,bg="firebrick",pch=21,col="black")
lines(seq(15,28),exp(reg.uni$coeff[1])*exp(reg.uni$coeff[2]*seq(15,28)+reg.uni$coeff[3]*35)*5^reg.uni$coeff[4],col="firebrick")
legend("topleft",c("O. universa","G. bulloides","G. ruber","G. sacculifer"),pch=c(21,21,21,21), col=c("black","black","black","black"), pt.bg=c("firebrick","cornflowerblue","gold3","palegreen3"),bty="n")

plot(ruber$Salinity,ruber.T26O5,bg="gold3",las=1,xlab="salinity (psu)",ylab="Mg/Ca @ 26ºC; Ω =5, (mmol/mol)",xlim=c(28,40),ylim=c(0,16),pch=21,col="black")
lines(seq(28,40,by=0.5),exp(reg.ruber$coeff[1])*exp(reg.ruber$coeff[2]*26+reg.ruber$coeff[3]*seq(28,40,by=0.5))*5^reg.ruber$coeff[4],col="gold3")
points(sac$Salinity,sac.T26O5,bg="palegreen3",pch=21,col="black")
lines(seq(28,40,by=0.5),exp(reg.sac$coeff[1])*exp(reg.sac$coeff[2]*26+reg.sac$coeff[3]*seq(28,40,by=0.5))*5^reg.sac$coeff[4],col="palegreen3")
points(bul$Salinity,bul.T26O5,bg="cornflowerblue",pch=21,col="black")
points(uni$Salinity,uni.T26O5,bg="firebrick",pch=21,col="black")
lines(seq(28,37,by=0.5),exp(reg.uni$coeff[1])*exp(reg.uni$coeff[2]*26+reg.uni$coeff[3]*seq(28,37,by=0.5))*5^reg.uni$coeff[4],col="firebrick")

plot(ruber$Ω.Calcite,ruber.T26S35,bg="gold3",las=1,xlab="Ω.Calcite",ylab="Mg/Ca @ 26ºC; 35 psu, (mmol/mol)",xlim=c(0,12),ylim=c(0,16),pch=21,col="black")
lines(seq(1,12,by=0.5),exp(reg.ruber$coeff[1])*exp(reg.ruber$coeff[2]*26+reg.ruber$coeff[3]*35)*seq(1,12,by=0.5)^reg.ruber$coeff[4],col="gold3")
points(sac$Ω.Calcite,sac.T26S35,bg="palegreen3",pch=21,col="black")
lines(seq(1,12,by=0.5),exp(reg.sac$coeff[1])*exp(reg.sac$coeff[2]*26+reg.sac$coeff[3]*35)*seq(1,12,by=0.5)^reg.sac$coeff[4],col="palegreen3")
points(bul$Ω.Calcite,bul.T26S35,bg="cornflowerblue",pch=21,col="black")
lines(seq(1,12,by=0.5),exp(reg.bul$coeff[1])*exp(reg.bul$coeff[2]*26+reg.bul$coeff[3]*35)*seq(1,12,by=0.5)^reg.bul$coeff[4],col="cornflowerblue")
points(uni$Ω.Calcite,uni.T26S35,bg="firebrick",pch=21,col="black")
lines(seq(1,12,by=0.5),exp(reg.uni$coeff[1])*exp(reg.uni$coeff[2]*26+reg.uni$coeff[3]*35)*seq(1,12,by=0.5)^reg.uni$coeff[4],col="firebrick")

plot(ruber$pH,ruber.T26S35.pH,bg="gold3",las=1,xlab="pH",ylab="Mg/Ca @ 26ºC; 35 psu, (mmol/mol)",xlim=c(7,9),ylim=c(0,16),pch=21,col="black")
lines(seq(7,9,by=0.2),summary(reg.ruber.pH)$coeff[1]*exp(summary(reg.ruber.pH)$coeff[2]*26+summary(reg.ruber.pH)$coeff[3]*35)+seq(7,9,by=0.2)*summary(reg.ruber.pH)$coeff[4],col="gold3")
points(bul$pH,bul.T26S35.pH,bg="cornflowerblue",pch=21,col="black")
lines(seq(7,9,by=0.2),summary(reg.bul.pH)$coeff[1]*exp(summary(reg.bul.pH)$coeff[2]*26+summary(reg.bul.pH)$coeff[3]*35)+seq(7,9,by=0.2)*summary(reg.bul.pH)$coeff[4],col="cornflowerblue")
points(uni$pH,uni.T26S35.pH,bg="firebrick",pch=21,col="black")
lines(seq(7,9,by=0.2),summary(reg.uni.pH)$coeff[1]*exp(summary(reg.uni.pH)$coeff[2]*26+summary(reg.uni.pH)$coeff[3]*35)+seq(7,9,by=0.2)*summary(reg.uni.pH)$coeff[4],col="firebrick")
