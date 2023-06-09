{
#Clear working environment
rm(list=ls())
#Load the needed libraries. Will need to install these libraries if they
#have not already been installed.
library(deSolve)
library(pracma)
}


#Define the differential equation
SIZR2 <- function(t,y,params){
  S = y[1]; I = y[2];S2 = y[3]; I2 = y[4]; Z = y[5]; R = y[6];
  with(
    as.list(params),
    {
      dS =e*f*R*(S+th*I)-B*S*Z-d*S
      dI =B*S*Z-(d+v)*I
      dS2 =e*f2*R*(S2+th*I2)-B2*S2*Z-d*S2
      dI2 =B2*S2*Z-(d+v)*I2
      dZ =s*(d+v)*(I+I2)-m*Z
      dR =w*R*(1-R/K)-R*(f*(S+I)+f2*(S2+I2))
      res = c(dS,dI,dS2,dI2,dZ,dR)
      list(res)
    }
  )
}


#Define a list of time values
times1=linspace(0,800,n=1000)
times2=linspace(0,62,n=100)

#Define parameter values that we will use.
#Note that th in the code is theta in the manuscript
#beta is written out as a variable name here
#s in the code is sigma in the manuscript
#w,K,e,th,d,fmax,fmin,fscale,v,s,m
w_ex=0.9
K_ex1=25
K_ex2=30
K_ex3=39.9
e_ex=.18
th_ex=.65
betaA=2.475841e-06
betaB=7.216188e-07
d_ex=0.05
v_ex=0.05
s_ex=1e+05
m_ex=1.5
fA=0.01231144
fB=0.009

#Run a simulation that finds frequencies of genotypes with the given traits,
#depending on the duration of the simulation (duration) and the carrying capacity of
#the resource (K).
freq_finder=function(duration,K){
    timestemp=linspace(0,duration,duration)
    temp=rk(c(S=5,I=0,S2=5,I2=0,Z=3600,R=K),
            timestemp,SIZR2,c(w=w_ex,K=K,e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
            ,method="rk45dp7")
    freqs=sum(temp[length(timestemp),2:3])/sum(temp[length(timestemp),2:5])
    freqs_av=mean(rowSums(temp[round(length(timestemp)*2/3):length(timestemp),2:3])/rowSums(temp[round(length(timestemp)*2/3):length(timestemp),2:5]))
    Hs_av=mean(rowSums(temp[round(length(timestemp)*2/3):length(timestemp),2:5]))
    prevs_av=mean(rowSums(temp[round(length(timestemp)*2/3):length(timestemp),c(3,5)])/rowSums(temp[round(length(timestemp)*2/3):length(timestemp),2:5]))
    prevs_A=mean(temp[round(length(timestemp)*2/3):length(timestemp),3]/rowSums(temp[round(length(timestemp)*2/3):length(timestemp),2:3]))
    prevs_B=mean(temp[round(length(timestemp)*2/3):length(timestemp),5]/rowSums(temp[round(length(timestemp)*2/3):length(timestemp),4:5]))
    R_av=mean(temp[round(length(timestemp)*2/3):length(timestemp),7])
    
  #Returns the frequency of genotype A at the end, and averages over the last
  #third of the time series. These averages, in order, are frequency of
  #genotype A, total host density, prevalence of infection, the prevalence
  #of infection experienced by genotype A, the prevalence of infection 
  #experienced by genotype B, and the density of resources.
  return(c(freqs,freqs_av,Hs_av,prevs_av,prevs_A,prevs_B,R_av))
}

#Similar to freq_finder but runs a simulation for genotype A alone
#and for genotype B alone.
solo_finder=function(duration,K){
  timestemp=linspace(0,duration,duration)
  #Run a simulation with just genotype A
  tempA=rk(c(S=10,I=0,S2=0,I2=0,Z=3600,R=K),
          timestemp,SIZR2,c(w=w_ex,K=K,e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
          ,method="rk45dp7")
  
  #Run a simulation with just genotype B
  tempB=rk(c(S=0,I=0,S2=10,I2=0,Z=3600,R=K),
           timestemp,SIZR2,c(w=w_ex,K=K,e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
           ,method="rk45dp7")
  
    Hs_avA=mean(rowSums(tempA[round(length(timestemp)*2/3):length(timestemp),2:5]))
    prevs_avA=mean(rowSums(tempA[round(length(timestemp)*2/3):length(timestemp),c(3,5)])/rowSums(tempA[round(length(timestemp)*2/3):length(timestemp),2:5]))
    Hs_avB=mean(rowSums(tempB[round(length(timestemp)*2/3):length(timestemp),2:5]))
    prevs_avB=mean(rowSums(tempB[round(length(timestemp)*2/3):length(timestemp),c(3,5)])/rowSums(tempB[round(length(timestemp)*2/3):length(timestemp),2:5]))

    #Return averages over the last third of the time series.
    #Respectively, these are the average density of genotype A,
    #The average infection prevalence experienced by genotype A,
    #and the same for genotype B.
  return(c(Hs_avA,prevs_avA,Hs_avB,prevs_avB))
}

#A similar function but for just one genotype and allows
#the user to input the traits of that genotype.
solo_finderfB=function(duration,K,f_in,B_in,th_in){
  timestemp=linspace(0,duration,duration)
    tempA=rk(c(S=10,I=0,S2=0,I2=0,Z=3600,R=K),
             timestemp,SIZR2,c(w=w_ex,K=K,e=e_ex,th=th_in,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=f_in,B=B_in,f2=fB,B2=betaB)
             ,method="rk45dp7")
    Hs_avA=mean(rowSums(tempA[round(length(timestemp)*2/3):length(timestemp),2:5]))
    prevs_avA=mean(rowSums(tempA[round(length(timestemp)*2/3):length(timestemp),c(3,5)])/rowSums(tempA[round(length(timestemp)*2/3):length(timestemp),2:5]))
    Ss_avA=mean(tempA[round(length(timestemp)*2/3):length(timestemp),2])
  return(c(Hs_avA,prevs_avA,Ss_avA))
}

#Some example K values for simulating, for plotting purposes.
K_ex1=40
K_ex2=63
K_ex3=150

K_range=linspace(25,150,200)

freq_av=array(NA,dim=length(K_range))
prev_av=array(NA,dim=length(K_range))
prev_A=array(NA,dim=length(K_range))
prev_B=array(NA,dim=length(K_range))
H_av=array(NA,dim=length(K_range))
R_av=array(NA,dim=length(K_range))
A_adv=array(NA,dim=length(K_range))
B_adv=array(NA,dim=length(K_range))

#Normally, duration is 150
for(i in 1:length(K_range)){
  temp=NA
  temp=freq_finder(150,K_range[i])
  freq_av[i]=temp[2]
  H_av[i]=temp[3]
  prev_av[i]=temp[4]
  prev_A[i]=temp[5]
  prev_B[i]=temp[6]
  R_av[i]=temp[7]
  print(i)
}

#Show some example simulations
initvar=0
K_example1=K_ex2
times_example1=linspace(0,5000,5000)
example1_simA=rk(c(S=10,I=0,S2=0,I2=0,Z=3600,R=K_example1),
                times_example1,SIZR2,c(w=w_ex,K=K_example1,e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
                ,method="rk45dp7")
example1_simB=rk(c(S=0,I=0,S2=10,I2=0,Z=3600,R=K_example1),
                times_example1,SIZR2,c(w=w_ex,K=K_example1,e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
                ,method="rk45dp7")
example1_simAB=rk(c(S=5,I=0,S2=5,I2=0,Z=3600,R=K_example1),
                 times_example1,SIZR2,c(w=w_ex,K=K_example1,e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
                 ,method="rk45dp7")


#Creates Figure A3 of the manuscript.
#png("Ch2_app_Simulations.png",width=10,height=6,units="in",res=300)
{
  cex_smallest_text=.5*1.1
  cex_minor_text=1*1.3
  cex_major_text=1.5*1.2
  cex_big_points=3
  lwd_minor=1.5*1.2
  lwd_major=2*1.2
  dens_choice=10
  par(cex.lab=cex_major_text,cex.axis=cex_major_text,cex=1,oma=c(3,0,2,0),font.lab=2,xaxs="i")
  par(mfcol=c(2,3),mar=c(2,5,1,1))
plot(log10(times_example1),example1_simA[,2],type="l",xaxt="n",xlab="",ylab=expression("Low-resistance"~italic(S)),ylim=range(c(example1_simA[,2:5],example1_simB[,2:5],example1_simAB[,2:5])),col="gray40",lty=2,lwd=lwd_major)
points(log10(times_example1),example1_simB[,2],type="l",col="gray40",lty=3,lwd=lwd_major)
points(log10(times_example1),example1_simAB[,2],type="l",lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_major_text)
text(log10(7),30,"Low resistance",cex=cex_minor_text,col="gray40")
text(log10(7),2,"High resistance",cex=cex_minor_text,col="gray40")
text(log10(3),10,"Evolving",cex=cex_minor_text,srt=25)

plot(log10(times_example1),example1_simA[,3],type="l",xaxt="n",xlab="",ylab=expression("Low-resistance"~italic(I)),ylim=range(c(example1_simA[,2:5],example1_simB[,2:5],example1_simAB[,2:5])),col="gray40",lty=2,lwd=lwd_major)
points(log10(times_example1),example1_simB[,3],type="l",col="gray40",lty=3,lwd=lwd_major)
points(log10(times_example1),example1_simAB[,3],type="l",lwd=lwd_major)
axis(side=1,at=log10(c(1,10,100,1000)),labels = c("1","10","100","1000"))
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_major_text)

plot(log10(times_example1),example1_simA[,4],type="l",xaxt="n",xlab="",ylab=expression("High-resistance"~italic(S)),ylim=range(c(example1_simA[,2:5],example1_simB[,2:5],example1_simAB[,2:5])),col="gray40",lty=2,lwd=lwd_major)
points(log10(times_example1),example1_simB[,4],type="l",col="gray40",lty=3,lwd=lwd_major)
points(log10(times_example1),example1_simAB[,4],type="l",lwd=lwd_major)
mtext(side=3,"Averages taken from days 100-150 for Fig. 2",line=1,cex=1.3)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_major_text)

plot(log10(times_example1),example1_simA[,5],type="l",xaxt="n",xlab="",ylab=expression("High-resistance"~italic(I)),ylim=range(c(example1_simA[,2:5],example1_simB[,2:5],example1_simAB[,2:5])),col="gray40",lty=2,lwd=lwd_major)
points(log10(times_example1),example1_simB[,5],type="l",col="gray40",lty=3,lwd=lwd_major)
points(log10(times_example1),example1_simAB[,5],type="l",lwd=lwd_major)
axis(side=1,at=log10(c(1,10,100,1000)),labels = c("1","10","100","1000"))
mtext(side=1,expression("Time on log scale:"~italic(t)~"(day)"),line=4,cex=1.3)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_major_text)

plot(log10(times_example1),example1_simA[,6],type="l",xaxt="n",xlab="",ylab=expression("Parasite propagules:"~italic(Z)),ylim=range(c(example1_simA[,6],example1_simB[,6],example1_simAB[,6])),col="gray40",lty=2,yaxt="n",lwd=lwd_major)
points(log10(times_example1),example1_simB[,6],type="l",col="gray40",lty=3,lwd=lwd_major)
points(log10(times_example1),example1_simAB[,6],type="l",lwd=lwd_major)
axis(side=2,at=c(0,50000,100000),labels=c(0,expression("5 x 10"^"4"),expression("1 x 10"^"5")))
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"E",cex=cex_major_text)

plot(log10(times_example1),example1_simA[,7],type="l",xaxt="n",xlab="",ylab=expression("Resources:"~italic(R)),ylim=range(c(example1_simA[,7],example1_simB[,7],example1_simAB[,7])),col="gray40",lty=2,lwd=lwd_major)
points(log10(times_example1),example1_simB[,7],type="l",col="gray40",lty=3,lwd=lwd_major)
points(log10(times_example1),example1_simAB[,7],type="l",lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"F",cex=cex_major_text)
axis(side=1,at=log10(c(1,10,100,1000)),labels = c("1","10","100","1000"))
}
#dev.off()


#Makes Figure A4 of the manuscript.
#png("Ch2_app_Resources_Prevalence.png",width=6,height=9,units="in",res=300)
{
  cex_smallest_text=.5*1.2
  cex_minor_text=1*1.2
  cex_major_text=1.5*1.5
  cex_big_points=3
  lwd_minor=1.5*1.2
  lwd_major=2*1.2
  dens_choice=10
  par(cex.lab=cex_major_text,cex.axis=cex_major_text,cex=1,oma=c(4,0,0,0),font.lab=2,xaxs="i")
  par(mfrow=c(3,1),mar=c(1,8,1,1))
plot(K_range,R_av,type="l",xlab="",ylab=expression(atop("Resource density:",italic(R)~"("~mu~"g chl"~italic(a)~"L"^"-1"~")")),lwd=lwd_major,xaxt="n")
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_major_text)

plot(K_range,prev_A,xlab="",ylab=expression(atop("Prevalence:",italic(p)~"(day"^"-1"~")")),type="l",col="gray40",lty=2,lwd=lwd_major,xaxt="n")
points(K_range,prev_B,type="l",col="gray40",lty=3,lwd=lwd_major)
legend(80,.35,c(expression(italic(p)["Low resistance"]),expression(italic(p)["High resistance"])),lty=c(2,3),lwd=lwd_major,col="gray40",cex=2)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_major_text)

plot(K_range,prev_A-prev_B,type="l",xlab="",ylab=expression(atop("Prevalence:",italic(p)["Low"]~"-"~italic(p)["High"]~"(day"^"-1"~")")),lwd=lwd_major)
mtext(side=1,expression("Carrying capacity:"~italic(K)~"("~mu~"g chl"~italic(a)~"L"^"-1"~")"),line=4,cex=1.5)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_major_text)
}
#dev.off()

#Check which K value maximizes the prevalence difference
K_range[which.max(prev_A-prev_B)]

freq_av2=array(NA,dim=length(K_range))
prev_av2=array(NA,dim=length(K_range))
H_av2=array(NA,dim=length(K_range))


#Normally, duration is 150. Run the simulations for 10 tims longer
#for comparison.
for(i in 1:length(K_range)){
  temp=NA
  temp=freq_finder(1500,K_range[i])
  freq_av2[i]=temp[2]
  H_av2[i]=temp[3]
  prev_av2[i]=temp[4]
  print(i)
}

#Running for 10 times longer, we see the simulation averages are qualitatively
#similar except that 
plot(K_range,freq_av2,type="l")
plot(K_range,prev_av2,type="l")
plot(K_range,H_av2,type="l")

#Can save these objects as .RDS objects for plotting on Figure 2.
#saveRDS(freq_av,"freq_av.RDS")
#saveRDS(H_av,"H_av.RDS")
#saveRDS(prev_av,"prev_av.RDS")

#It's also helpful to get some simulation outputs for each genotype
#alone.
prev_av_solo=array(NA,dim=c(length(K_range),2))
H_av_solo=array(NA,dim=c(length(K_range),2))
for(i in 1:length(K_range)){
  temp=NA
  temp=solo_finder(150,K_range[i])
  H_av_solo[i,1]=temp[1]
  H_av_solo[i,2]=temp[3]
  prev_av_solo[i,1]=temp[2]
  prev_av_solo[i,2]=temp[4]
  print(i)
}

#Save these for the purposes of plotting Figure 2.
#saveRDS(prev_av_solo,"prev_av_solo.RDS")
#saveRDS(H_av_solo,"H_av_solo.RDS")


#Run single genotype simulations at the focal values of carrying capacity.
K1_solos=solo_finder(150,K_ex1)
K2_solos=solo_finder(150,K_ex2)
K3_solos=solo_finder(150,K_ex3)

#Construct regressions of host density and prevalence and save outputs for
#plotting in Figure 3 of the manuscript.
H_line0=c(mean(K2_solos[3]),H_av[which.min(abs(K_range-K_ex2))],mean(K2_solos[1]))
prev_line0=c(mean(K2_solos[4]),prev_av[which.min(abs(K_range-K_ex2))],mean(K2_solos[2]))

H_line=c(mean(K3_solos[3]),H_av[which.min(abs(K_range-K_ex3))],mean(K3_solos[1]))
prev_line=c(mean(K3_solos[4]),prev_av[which.min(abs(K_range-K_ex3))],mean(K3_solos[2]))
beta_line=c(betaB,freq_av[which.min(abs(K_range-K_ex3))]*betaA+(1-freq_av[which.min(abs(K_range-K_ex3))])*betaB,betaA)

H_lm=lm(H_line~beta_line)

library(betareg)
prev_curve=betareg(prev_line~beta_line)

new_data=data.frame(beta_line=freq_av[which.min(abs(K_range-K_ex2))]*betaA+(1-freq_av[which.min(abs(K_range-K_ex2))])*betaB)
new_data2=data.frame(beta_line=linspace(min(beta_line)/2,max(beta_line)*1.5,1e4))
predict_H=predict(H_lm,new_data)
predict_prev=predict(prev_curve,new_data)
predict_prev_pretty=predict(prev_curve,new_data2)

summary(H_lm)
int=6.216e1
slope=-7.505e6

#Can save these for plotting Figure 3.
#saveRDS(H_line0,"H_line0.RDS")
#saveRDS(prev_line0,"prev_line0.RDS")
#saveRDS(H_line,"H_line.RDS")
#saveRDS(prev_line,"prev_line.RDS")
#saveRDS(beta_line,"beta_line.RDS")
#saveRDS(as.numeric(predict_prev_pretty),"prev_curve.RDS")
#saveRDS(as.numeric(new_data2$beta_line),"beta_curve.RDS")
#saveRDS(c(new_data$beta_line,predict_prev),"predict_prev.RDS")
#saveRDS(c(new_data$beta_line,predict_H),"predict_H.RDS")


K_length=200
K_series2=linspace(04,150,K_length)
example_middle=0.01513478
example_top=0.02064133
example_middleB=5.216614e-06
example_topB=4.091753e-05

#Here are the K values for which there can be oscillations for the middle
#value of beta and f as well as the top value of beta and f. Simulate
#these for plotting in Figure 2.
middle_Ks=seq(107,200,1)
top_Ks=seq(4,200,1)

#First column is prev, second is H, third is S
Fig1_outputs_middle=array(NA,dim=c(K_length,3))
Fig1_outputs_top=array(NA,dim=c(K_length,3))
for(i in 1:K_length){
  temp1=NA    
  if(i %in% middle_Ks){
    temp1=solo_finderfB(100000,K_series2[i],example_middle,example_middleB,0)
  }

  temp2=NA
  if(i %in% top_Ks){
    temp2=solo_finderfB(100000,K_series2[i],example_top,example_topB,0)    
  }
  
  Fig1_outputs_middle[i,1]=temp1[2]
  Fig1_outputs_middle[i,2]=temp1[1]
  Fig1_outputs_middle[i,3]=temp1[3]
  
  Fig1_outputs_top[i,1]=temp2[2]
  Fig1_outputs_top[i,2]=temp2[1]
  Fig1_outputs_top[i,3]=temp2[3]
  print(i)
}

#Set rows with tiny host density to have zero prevalence.
Fig1_outputs_middle[which(Fig1_outputs_middle[,2]<1e-16),1]=NA
Fig1_outputs_top[which(Fig1_outputs_top[,2]<1e-16),1]=NA

#Save these for plotting two of the gray curves on the left
#panels of Figure 1.
#saveRDS(Fig1_outputs_middle,"Fig1_outputs_middle.RDS")
#saveRDS(Fig1_outputs_top,"Fig1_outputs_top.RDS")


#Check change in last one-fifth. Note that this will take a LONG time.
timestemp=linspace(0,500000,500000)

#Each row corresponds to a K value. The first column
#is the maximum change that occurs in the last one-fifth
#of the time series for genotype A alone in a simulation.
#The second column is for genotype B, the third is for both
#genotypes together.
change_array=array(NA,dim=c(length(K_range),3))

for(k in 1:length(K_range)){
  tic()
  A_alone=NA
  B_alone=NA
  AB_together=NA
  #A_alone=rk(c(S=10,I=0,S2=0,I2=0,Z=3600,R=K_range[k]),
  #           timestemp,SIZR2,c(w=w_ex,K=K_range[k],e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
  #           ,method="rk45dp7")
  #B_alone=rk(c(S=0,I=0,S2=10,I2=0,Z=3600,R=K_range[k]),
  #           timestemp,SIZR2,c(w=w_ex,K=K_range[k],e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
  #           ,method="rk45dp7")
  AB_together=rk(c(S=5,I=0,S2=5,I2=0,Z=3600,R=K_range[k]),
                 timestemp,SIZR2,c(w=w_ex,K=K_range[k],e=e_ex,th=th_ex,d=d_ex,v=v_ex,s=s_ex,m=m_ex,f=fA,B=betaA,f2=fB,B2=betaB)
                 ,method="rk45dp7")
  
  #change_array[k,1]=max(abs(A_alone[round(length(timestemp)*4/5),2:7]-A_alone[length(timestemp),2:7]))
  #change_array[k,2]=max(abs(B_alone[round(length(timestemp)*4/5),2:7]-B_alone[length(timestemp),2:7]))
  change_array[k,3]=max(abs(AB_together[round(length(timestemp)*4/5),2:7]-AB_together[length(timestemp),2:7]))
  print(k)
  toc()
}

#Confirm that the maximum change observed was very small.
max(change_array)
