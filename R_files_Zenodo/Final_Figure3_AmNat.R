rm(list=ls())
#Load needed libraries. Use the following if you don't have a library installed, eg. lme4.
#install.packages("lme4")
{
  require(lme4)
  require(pracma)
  require(stats)
  require(boot)
  require(betareg)
  require(DHARMa)
  require(car)
}


#Transmission rate estimates from "Success, failure and ambiguity of the dilution effect among
#competitors"
betaB=7.216188e-7
betaA=2.475841e-6

#Set working directory wherever this data file was downloaded
#setwd("C:/Users/Jason/OneDrive - University of Pittsburgh/Documents/Thesis Chapter 2")
data_focal<-readRDS("data_focal.RDS")
#The same data is also provided as .csv files, for those who prefer that.

#Data is formatted such that Disease column indicates these tanks received spores
#Phosph column indicates that nutrient supply was 5 micrograms P/L or 50
#Genotype column indicates genotype treatment was 1=A43 (less resistant), 2=Bristol 10 (more resistant), or 4=A43+Bristol 10
#Tank column give tank ID
#Columns Juveniles through Males.1 are numbers of individual hosts in the 1 L sample
#Juveniles and InfJuveniles are the number of Juveniles and Infected Juveniles, respectively
#Swith0eggs is susceptible adults with 0 eggs in the brood chamber
#X1-X15 are susceptible adults with 1-15 eggs in the brood chamber
#Males and Males.1 are susceptible adult males and infected adult males, respectively
#Iwith0eggs is infected adults with zero eggs
#X1.1-X15.1 are infected adults with 1-15 eggs.
#Chlorophyll is micrograms of chlorophyll a per liter
#Each layer (third dimension) of data_focal represents a different sample date falling, in order, on the time vector below.
#E.g., data_focal[,,1] gives measurements on day 14.

#Daphnia added on day 1, 13 days before first sample data 
#Spores added July 21, 2015 (day = 28). 4th sample date is 
#the first one after spore addition
time=c(14,20,27,30,34,36,41,44,48,51,55,58,62,64,69,72,76,79,83,86)

#Only consider times two weeks or more after infection to the end of the experiment
timetrunc=time[9:17]
#truncate data in time to incorporate days 48 to 76 (inclusive)
data1=data_focal[,,9:17]

#Dropping problematic tubs, 9 for extinction and 36 for extremely and uncharacteristically low density
data2=data1[-which(data1[,4,1]%in%c(9,36)),,]

#A is A43, B is Bristol 10. AB is A43 and Bristol 10. The other
#genotype IDs are not used for this manuscript.
gen_ID=function(gen_num){
  c("A","B","S","AB","AS","BS")[gen_num]
}

#A function for finding tank ID #s for a given genotype (Gen), nutrient (Nut), and disease (Dis) treatment
replicates=function(Gen,Nut,Dis){
  data2[which(data2[,3,1]==Gen&data2[,2,1]==Nut&data2[,1,1]==Dis),4,1]
}

#Read in the genotype data. First column is sample date. Second is treatment. Third is tank ID. The fourth is the genotype ID of an individual animal
#There are up to 10 samples for a given tank and sample date
Genotype_data<-read.csv("Genotype_data.csv")

#Convert genotype identity to 1 (A43) or 0 (Bristol 10) for binomial model
#and exclude NA values
binary_data<-Genotype_data[!is.na(Genotype_data$genotype),1:4]
binary_data$A43_ID<-as.numeric(binary_data$genotype=="A")
#Convert treatments into 5s or 50s to represent nutrients
binary_data$Nutrient<-5+45*(binary_data$Treatment=="D_AB_50ug")

#Fig. 3A stats
m<-glmer(A43_ID~Nutrient+(1|Tank)+(1|Day),data=binary_data,family=binomial)
#Focal p-value for mixed effects binomial regression
summary(m)
plot(simulateResiduals(m))


#Can remove day as a random effect since there was extremely little
#variance associated with day. Doing so has very little effect on the
#model results, see m2 here. This is the version reported in the manuscript,
#but the results are extremely similar either way.
m2<-glmer(A43_ID~Nutrient+(1|Tank),data=binary_data,family=binomial)
summary(m2)
plot(simulateResiduals(m2))

#Raw percentages A43 of IDed individuals at each nutrient
mean(binary_data$A43_ID[binary_data$Nutrient==5])
mean(binary_data$A43_ID[binary_data$Nutrient==50])

#Now want to complement these with ecological plots
prev_tank=function(Tank_ID){
  row_num=which(data2[,4,1]==Tank_ID)
  trapz(timetrunc,ifelse(apply(data2[row_num,5:40,],2,sum)>0,apply(data2[row_num,23:40,],2,sum)/apply(data2[row_num,5:40,],2,sum),0))/diff(range(timetrunc))
}

H_tank=function(Tank_ID){
  row_num=which(data2[,4,1]==Tank_ID)
  trapz(timetrunc,apply(data2[row_num,5:40,],2,sum))/diff(range(timetrunc))
}
I_tank=function(Tank_ID){
  row_num=which(data2[,4,1]==Tank_ID)
  trapz(timetrunc,apply(data2[row_num,23:40,],2,sum))/diff(range(timetrunc))
}
Chl_tank=function(Tank_ID){
  row_num=which(data2[,4,1]==Tank_ID)
  trapz(timetrunc,data2[row_num,41,])/diff(range(timetrunc))
}

#Calculate average transmssion rate, beta, for a tank given the average
#frequency of A43 among genotyped individuals
beta_from_freq<-function(freqA){
  betaB+freqA*(betaA-betaB)
}

High_nut_betas=c(beta_from_freq(mean(binary_data[which(binary_data[,3]==replicates(4,50,1)[1]),5],na.rm=T)),beta_from_freq(mean(binary_data[which(binary_data[,3]==replicates(4,50,1)[2]),5],na.rm=T)),beta_from_freq(mean(binary_data[which(binary_data[,3]==replicates(4,50,1)[3]),5],na.rm=T)),betaA,betaA,betaA,betaB,betaB,betaB)
High_nut_prevs=c(prev_tank(replicates(4,50,1)[1]),prev_tank(replicates(4,50,1)[2]),prev_tank(replicates(4,50,1)[3]),prev_tank(replicates(1,50,1)[1]),prev_tank(replicates(1,50,1)[2]),prev_tank(replicates(1,50,1)[3]),prev_tank(replicates(2,50,1)[1]),NA,prev_tank(replicates(2,50,1)[2]))
High_nut_Hs=c(H_tank(replicates(4,50,1)[1]),H_tank(replicates(4,50,1)[2]),H_tank(replicates(4,50,1)[3]),H_tank(replicates(1,50,1)[1]),H_tank(replicates(1,50,1)[2]),H_tank(replicates(1,50,1)[3]),H_tank(replicates(2,50,1)[1]),NA,H_tank(replicates(2,50,1)[2]))
High_nut_Is=c(I_tank(replicates(4,50,1)[1]),I_tank(replicates(4,50,1)[2]),I_tank(replicates(4,50,1)[3]),I_tank(replicates(1,50,1)[1]),I_tank(replicates(1,50,1)[2]),I_tank(replicates(1,50,1)[3]),I_tank(replicates(2,50,1)[1]),NA,I_tank(replicates(2,50,1)[2]))

Low_nut_betas=c(beta_from_freq(mean(binary_data[which(binary_data[,3]==replicates(4,5,1)[1]),5],na.rm=T)),beta_from_freq(mean(binary_data[which(binary_data[,3]==replicates(4,5,1)[2]),5],na.rm=T)),beta_from_freq(mean(binary_data[which(binary_data[,3]==replicates(4,5,1)[3]),5],na.rm=T)),betaA,betaA,betaA,betaB,betaB,betaB)
Low_nut_prevs=c(prev_tank(replicates(4,5,1)[1]),prev_tank(replicates(4,5,1)[2]),prev_tank(replicates(4,5,1)[3]),NA,prev_tank(replicates(1,5,1)[1]),prev_tank(replicates(1,5,1)[2]),prev_tank(replicates(2,5,1)[1]),prev_tank(replicates(2,5,1)[2]),prev_tank(replicates(2,5,1)[3]))
Low_nut_Hs=c(H_tank(replicates(4,5,1)[1]),H_tank(replicates(4,5,1)[2]),H_tank(replicates(4,5,1)[3]),NA,H_tank(replicates(1,5,1)[1]),H_tank(replicates(1,5,1)[2]),H_tank(replicates(2,5,1)[1]),H_tank(replicates(2,5,1)[2]),H_tank(replicates(2,5,1)[3]))
Low_nut_Is=c(I_tank(replicates(4,5,1)[1]),I_tank(replicates(4,5,1)[2]),I_tank(replicates(4,5,1)[3]),NA,I_tank(replicates(1,5,1)[1]),I_tank(replicates(1,5,1)[2]),I_tank(replicates(2,5,1)[1]),I_tank(replicates(2,5,1)[2]),I_tank(replicates(2,5,1)[3]))

plot_predictor=High_nut_betas*1e6
plot_predictor_low=Low_nut_betas*1e6

#High nut I stats
lm_for_plot_I=lm(High_nut_Is~plot_predictor)
summary(lm_for_plot_I)
#plot(simulateResiduals(lm_for_plot_I))
#Low nut I stats
lm_for_plot_I_low=lm(Low_nut_Is~plot_predictor_low)
summary(lm_for_plot_I_low)
#plot(simulateResiduals(lm_for_plot_I_low))

#Fig. 3C stats
summary(mod4C<-betareg(Low_nut_prevs~plot_predictor_low))
#plot(mod4C)

#Fig. 3D stats
summary(mod4D<-lm(Low_nut_Hs~plot_predictor_low))
#plot(mod4D)

#Fig. 3E stats
mod_high_prev=betareg(High_nut_prevs~plot_predictor)
summary(mod_high_prev)
#plot(mod_high_prev)
new_beta=data.frame(plot_predictor=c(1.024547,1.69))
predict(mod_high_prev,new_beta)

#Fig. 3F stats
mod_high_H=lm(High_nut_Hs~plot_predictor)
summary(mod_high_H)
new_beta=data.frame(plot_predictor=c(1.02,1.69))
predict(mod_high_H,new_beta)
#plot(mod_high_H)

beta_example<-seq(0,3,.001)
new_beta_example<-data.frame(plot_predictor=beta_example)
predicted_prevs<-predict(mod_high_prev,new_beta_example)

#Make a simple function for calculating standard error
standard_error<-function(data_vec){
  sd(data_vec,na.rm = T)/sqrt(length(which(!is.na(data_vec))))  
}

#Calculate standard errors for putting bars on Figure 3.
mean_p_lowA=mean(Low_nut_prevs[4:6],na.rm=T)
mean_p_highA=mean(High_nut_prevs[4:6],na.rm=T)
mean_p_lowB=mean(Low_nut_prevs[7:9],na.rm=T)
mean_p_highB=mean(High_nut_prevs[7:9],na.rm=T)
mean_p_lowAB=mean(Low_nut_prevs[1:3],na.rm=T)
mean_p_highAB=mean(High_nut_prevs[1:3],na.rm=T)
se_p_lowA=standard_error(Low_nut_prevs[4:6])
se_p_highA=standard_error(High_nut_prevs[4:6])
se_p_lowB=standard_error(Low_nut_prevs[7:9])
se_p_highB=standard_error(High_nut_prevs[7:9])
se_p_lowAB=standard_error(Low_nut_prevs[1:3])
se_p_highAB=standard_error(High_nut_prevs[1:3])

#Calculate means for plotting on Figure 3.
mean_H_lowA=mean(Low_nut_Hs[4:6],na.rm=T)
mean_H_highA=mean(High_nut_Hs[4:6],na.rm=T)
mean_H_lowB=mean(Low_nut_Hs[7:9],na.rm=T)
mean_H_highB=mean(High_nut_Hs[7:9],na.rm=T)
mean_H_lowAB=mean(Low_nut_Hs[1:3],na.rm=T)
mean_H_highAB=mean(High_nut_Hs[1:3],na.rm=T)
se_H_lowA=standard_error(Low_nut_Hs[4:6])
se_H_highA=standard_error(High_nut_Hs[4:6])
se_H_lowB=standard_error(Low_nut_Hs[7:9])
se_H_highB=standard_error(High_nut_Hs[7:9])
se_H_lowAB=standard_error(Low_nut_Hs[1:3])
se_H_highAB=standard_error(High_nut_Hs[1:3])

shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}

# And here is an example of use:
# pdf(file="test2.pdf", width=2, height=2); par(mar=c(0,0,0,0)+.1)
plot(c(0,1), c(0,1), type="n", lwd=20, axes=FALSE, xlab="", ylab="")

rect(xleft = 0.5, xright = 1, ybottom = 0, ytop = 1, col=1)
text(1/6, 1/6, "Test 1")
shadowtext(2/6, 2/6, "Test 2", col='red', bg="blue")
shadowtext(3/6, 3/6, "Test 3", cex=2)

# `r` controls the width of the border
shadowtext(5/6, 5/6, "Test 4", col="black", bg="white", cex=4, r=0.2)

#Ch2_Mesocosm.png. Makes manuscript figure 3.
pdf("AmNat_Fig3.pdf",width=14,height=9,family="ArialMT",useDingbats=FALSE)
{
  #Since beta_from_freq is a linear function of frequency, we can take beta_from_freq
  #after finding mean frequency and standard error of the mean frequency
  mean_low=mean(binary_data$A43_ID[binary_data$Nutrient==5])
  se_low=sqrt((mean_low*(1-mean_low))/length(binary_data$A43_ID[binary_data$Nutrient==5]))
  mean_high=mean(binary_data$A43_ID[binary_data$Nutrient==50])
  se_high=sqrt((mean_high*(1-mean_high))/length(binary_data$A43_ID[binary_data$Nutrient==50]))

  cex_smallest_text=.5*1.2
  cex_minor_text=1*1.2
  cex_major_text=1.5*1.2
  cex_big_points=3
  lwd_minor=1.5*1.2
  lwd_major=2*1.2
  dens_choice=10
  par(mar=c(5,8.5,1,.75),cex.lab=cex_major_text,cex.axis=cex_major_text,cex=1,oma=c(0,0,0,0),font.lab=2,xaxs="i")
  
  m<-rbind(c(1,1,2,2,3,3),c(6,4,4,5,5,7))
  layout(m)

  plot(c(0,1),c(mean_low,mean_high),ylab="",xlab="",xlim=c(-.5,1.5),ylim=c(0,1),pch=21,bg="darkgrey",cex=cex_big_points,xaxt="n")
  polygon(c(-.1,1.1,1.1,-.1),c(.125,.125,.55,.55),col=rgb(0,0,1,0.25,maxColorValue=1),border=rgb(0,0,1,0.5,maxColorValue=1),fillOddEven="non-zero")
  segments(0,mean_low-se_low,0,mean_low+se_low,lwd=lwd_major)
  segments(1,mean_high-se_high,1,mean_high+se_high,lwd=lwd_major)
  points(c(0,1),c(mean_low,mean_high),ylab="",xlab="",xlim=c(-.5,1.5),ylim=c(0,1),pch=21,bg="darkgrey",cex=cex_big_points,xaxt="n")
  mtext(side=2,expression(atop("Freq. of low","resistance genotype")),cex=cex_major_text,line=2.6)
  text(0.48,.3,expression(atop("Resistance","is futile")),cex=cex_major_text,col=rgb(0,0,1,0.5,maxColorValue=1),adj=c(0.5,0.5))
  points(c(0,1),c(mean_low,mean_high),ylab=expression(atop("Evolved transmission rate:",italic(beta)~" (L parasite"^"-1"~" day"^"-1"~" x 10"^"-6"~")")),xlab="",xlim=c(-.5,1.5),ylim=1e6*c(betaB-1e-7,betaA+1e-7),pch=21,bg="darkgrey",cex=cex_big_points,xaxt="n")
  axis(side=1,at=c(0,1),labels=c("Low nut.","High nut."),font=1)
  points(c(0,1),1e6*c(mean_low,mean_high),ylab=expression(atop("Evolved transmission rate:",italic(beta)~" (L parasite"^"-1"~" day"^"-1"~" x 10"^"-6"~")")),xlab="",xlim=c(-.5,1.5),ylim=1e6*c(mean_low-se_low,mean_high+se_high),pch=21,bg="darkgrey",cex=cex_big_points,xaxt="n")
  points(c(-.05,.95),c(0,0),pch=24,cex=cex_big_points,bg="white",lwd=lwd_major)
  points(c(.05,1.05),c(1,1),cex=cex_big_points,pch=15,col="black")
  text((par("usr")[2]-par("usr")[1])*.95+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_major_text)
  par(xpd=NA)
  legend(-1.2,-0.1,c("Low resistance","Evolving","High resistance"),pch=c(15,21,24),pt.bg=c("black","darkgrey","white"),cex=cex_big_points,bty="n")
  par(xpd=F)
      
  plot(c(.05,1.05),c(mean_p_lowA,mean_p_highA),ylab="",ylim=c(0,mean_p_highA+se_p_highA),xaxt="n",cex=cex_big_points,pch=15,col="black",lwd=lwd_major,xlab="",xlim=c(-.5,1.5))
  axis(side=1,at=c(0,1),labels=c("Low nut.","High nut."),font=1)
  mtext(side=2,expression(atop("Prevalence:",italic(p)~" (unitless)")),line=2.6,cex=cex_major_text)
  mtext(side=1,expression("Nutrient treatment"),line=3.5,cex=cex_major_text)
  segments(-.05,mean_p_lowB-se_p_lowB,-.05,mean_p_lowB+se_p_lowB,lwd=lwd_major)
  segments(1-.05,mean_p_highB-se_p_highB,1-.05,mean_p_highB+se_p_highB,lwd=lwd_major)
  segments(.05,mean_p_lowA-se_p_lowA,.05,mean_p_lowA+se_p_lowA,lwd=lwd_major)
  segments(1+.05,mean_p_highA-se_p_highA,1+.05,mean_p_highA+se_p_highA,lwd=lwd_major)
  segments(0,mean_p_lowAB-se_p_lowAB,0,mean_p_lowAB+se_p_lowAB,lwd=lwd_major)
  segments(1,mean_p_highAB-se_p_highAB,1,mean_p_highAB+se_p_highAB,lwd=lwd_major)
  points(c(-.05,.95),c(mean_p_lowB,mean_p_highB),pch=24,cex=cex_big_points,bg="white",lwd=lwd_major)
  points(c(0,1),c(mean_p_lowAB,mean_p_highAB),pch=21,bg="darkgrey",lwd=lwd_major,cex=cex_big_points)
  text((par("usr")[2]-par("usr")[1])*.95+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_major_text)

  
  plot(c(.05,1.05),c(mean_H_lowA,mean_H_highA),ylab="",ylim=c(0,mean_H_highB+se_H_highB),xaxt="n",cex=cex_big_points,pch=15,col="black",lwd=lwd_major,xlab="",xlim=c(-.5,1.5))
  segments(-.05,mean_H_lowB-se_H_lowB,-.05,mean_H_lowB+se_H_lowB,lwd=lwd_major)
  segments(1-.05,mean_H_highB-se_H_highB,1-.05,mean_H_highB+se_H_highB,lwd=lwd_major)
  segments(.05,mean_H_lowA-se_H_lowA,.05,mean_H_lowA+se_H_lowA,lwd=lwd_major)
  segments(1+.05,mean_H_highA-se_H_highA,1+.05,mean_H_highA+se_H_highA,lwd=lwd_major)
  segments(0,mean_H_lowAB-se_H_lowAB,0,mean_H_lowAB+se_H_lowAB,lwd=lwd_major)
  segments(1,mean_H_highAB-se_H_highAB,1,mean_H_highAB+se_H_highAB,lwd=lwd_major)
  mtext(side=2,expression(atop("Host density:",italic(H)~"(Hosts L"^"-1"~")")),line=2.6,cex=cex_major_text)
  axis(side=1,at=c(0,1),labels=c("Low nut.","High nut."),font=1)
  points(c(-.05,.95),c(mean_H_lowB,mean_H_highB),pch=24,cex=cex_big_points,bg="white",lwd=lwd_major)
  points(c(0,1),c(mean_H_lowAB,mean_H_highAB),pch=21,bg="darkgrey",lwd=lwd_major,cex=cex_big_points)
  text((par("usr")[2]-par("usr")[1])*.95+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_major_text)
  
  plot(High_nut_betas*1e6,High_nut_prevs,xlim=c(betaB-1e-7,betaA+1e-7)*1e6,ylim=c(0,.3),xlab="",ylab=""
       ,pch=c(21,21,21,15,15,15,24,24,24),bg=c("darkgrey","darkgrey","darkgrey","black","black","black","white","white","white"),
       cex=cex_major_text,lwd=lwd_minor,yaxt="n")
  predict_point<-(mean_low*(betaA-betaB)+betaB)*1e6
  new_point<-data.frame(plot_predictor=predict_point)
  prev_point<-predict(mod_high_prev,new_point)
  points((mean_high*(betaA-betaB)+betaB)*1e6,mean(High_nut_prevs[1:3]),cex=cex_major_text*2,pch=21,bg="darkgrey")
  arrows(x0=(mean_low*(betaA-betaB)+betaB)*1e6,y0=prev_point,x1=(mean_high*(betaA-betaB)+betaB)*1e6,y1=mean(High_nut_prevs[1:3]),col=rgb(0,0,1,0.5,maxColorValue=1),lwd=lwd_major*2)
  axis(side=2,cex=cex_major_text,at=c(0,.1,.2,.3),labels=c(0,0.1,0.2,0.3))
  points(beta_example,predicted_prevs,lwd=lwd_minor,type="l")
  text((par("usr")[2]-par("usr")[1])*.95+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_major_text)
  text(1.05,0.3/20,cex=cex_major_text,col="gray40",expression(atop("High","resistance")))
  text(2.25,0.3/20,cex=cex_major_text,col="gray40",expression(atop("Low","resistance")))
  mtext(side=2,expression(atop("Prevalence:",italic(p)~" (unitless)")),line=2.6,cex=cex_major_text)
  text(1,0.17,"Consequences",cex=cex_major_text,col=rgb(0,0,1,0.5,maxColorValue=1))
  points((mean_low*(betaA-betaB)+betaB)*1e6,prev_point,pch=21,bg=rgb(0,0,1,0.5,maxColorValue=1),cex=cex_big_points*1.2)
  
  plot(High_nut_betas*1e6,High_nut_Hs,xlim=c(betaB-1e-7,betaA+1e-7)*1e6,ylim=c(0,182),xlab="",ylab="",type="l",col="white",cex=cex_major_text)
  H_point<-predict(mod_high_H,new_point)
  points((mean_high*(betaA-betaB)+betaB)*1e6,mean(High_nut_Hs[1:3]),pch=21,bg="darkgrey",cex=cex_major_text*2)
  arrows(x0=(mean_low*(betaA-betaB)+betaB)*1e6,y0=H_point,x1=(mean_high*(betaA-betaB)+betaB)*1e6,y1=mean(High_nut_Hs[1:3]),col=rgb(0,0,1,0.5,maxColorValue=1),lwd=lwd_major*2)
  abline(mod_high_H,lwd=lwd_major)
  points(High_nut_betas*1e6,High_nut_Hs,xlim=c(betaB,betaA)*1e6,ylim=c(0,177),xlab="",ylab=""
         ,pch=c(21,21,21,15,15,15,24,24,24),bg=c("darkgrey","darkgrey","darkgrey","black","black","black","white","white","white"),
         cex=cex_major_text,lwd=lwd_minor)
  mtext(side=1,expression("Transmission rate: "~italic(beta)~" (L parasite"^"-1"~" day"^"-1"~" x 10"^"-6"~")"),line=3.5,adj=1.35,cex=cex_major_text)
  text((par("usr")[2]-par("usr")[1])*.95+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"E",cex=cex_major_text)
  text(1,177/20,cex=cex_major_text,col="gray40",expression(atop("High","resistance")))
  text(2.25,177/20,cex=cex_major_text,col="gray40",expression(atop("Low","resistance")))
  mtext(side=2,expression(atop("Host density:",italic(H)~"(Hosts L"^"-1"~")")),line=2.6,cex=cex_major_text)
  text(1,105,"Consequences",cex=cex_major_text,col=rgb(0,0,1,0.5,maxColorValue=1))
  points((mean_low*(betaA-betaB)+betaB)*1e6,H_point,pch=21,bg=rgb(0,0,1,0.5,maxColorValue=1),cex=cex_big_points*1.2)
}
dev.off()

#Test the effect of nutrients and genotype on prevalence, stats for Fig. 3B
Prevalences<-c(Low_nut_prevs[1:3],High_nut_prevs[1:3],Low_nut_prevs[5:6],High_nut_prevs[4:6],Low_nut_prevs[7:9],High_nut_prevs[c(7,9)])
Hosts<-c(Low_nut_Hs[1:3],High_nut_Hs[1:3],Low_nut_Hs[5:6],High_nut_Hs[4:6],Low_nut_Hs[7:9],High_nut_Hs[c(7,9)])
Gens=c(rep("AB",6),rep("A",5),rep("B",5))
#Ensure nutrients are centered
Nuts=c(-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,-1,1,1)
levels(Nuts)=c(-1,1)
summary(Prev_model<-betareg(Prevalences~Nuts))
summary(H_model<-lm(Hosts~Nuts))
#plot(H_model)


