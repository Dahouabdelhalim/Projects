setwd("~/Documents/PhD Project/Respirometry/Respirometry")

#only use the 5% dataset for simple method
cobia<-read.csv("MO2_methods10_cobia_20180914.csv")
cobia<-subset(cobia, select=-c(X,NumSMR,Numdat))

#assign NAs to Ccrit and Scrit for 30_28 and 34_28 because the automation method to calculate ccrit from SMR was not accurate for those two fish
cobia$Ccrit[which(cobia$Animal_ID=="C30"&cobia$Test_Temp=="28")]<-NA
cobia$Scrit[which(cobia$Animal_ID=="C30"&cobia$Test_Temp=="28")]<-NA
cobia$Ccrit[which(cobia$Animal_ID=="C34"&cobia$Test_Temp=="28")]<-NA
cobia$Scrit[which(cobia$Animal_ID=="C34"&cobia$Test_Temp=="28")]<-NA

#need to add extra rows for cobia32 and cobia37 for temps that weren't completed on those sharks
cobia[35:75,]<-NA
cobia$Animal_ID[35:75]<-c("C08","C09","C09","C11","C11","C12","C12","C13","C13","C15","C15","C16","C16",
                          "C17","C17","C18","C18","C19","C19","C20","C20","C21","C21","C22","C22",
                          "C23","C24","C25","C26","C27","C27","C28","C30","C31","C31","C32","C32","C33","C34","C35","C35")
cobia$Test_Temp[35:75]<-c(32,28,32,28,32,28,32,28,32,28,32,28,32,28,32,28,32,28,32,28,32,28,32,28,32,24,24,24,24,24,32,24,
                          24,24,32,24,28,24,24,24,28)
cobia$Animal_Mass[35:75]<-c(5.47,4.62,4.62,4.17,4.17,8.02,8.02,8.20,8.20,5.22,5.22,4.17,4.17,
                            4.75,4.75,5.75,5.75,5.45,5.45,4.39,4.39,3.50,3.50,5.88,5.88,
                            4.14,3.86,4.70,9.56,4.16,4.16,6.57,3.45,3.42,3.42,3.14,3.14,5.45,4.6,6.74,6.74)


######
#convert Ccrit to Pcrit
#add 100% saturation in umol/l (it's important to keep order)
cobia$O2_umol.l<-c(234,220,235.7,230,233,232.5,234,234,233,232,232,232,233,233,
                   224.7,206.8,223,206.8,223,206.8,222.2,206,223,222.2,206,223.5,206.8,
                   223.5,206,222.2,205.7,222.2,206,205.7,rep(NA,41))

#convert to mg/l
cobia$O2_mg.l<-cobia$O2_umol.l*32/1000
#add barometric pressure in inches of Hg, from https://www.timeanddate.com/weather/@4794972/historic from Yorktown, VA
cobia$bp<-c(29.81,30.23,29.73,30.07,29.75,30.26,30.01,30.02,29.87,29.95,30.03,30.28,29.96,30.22,
            30.29,29.83,30.24,30.21,30.25,30.16,30.02,30.14,29.94,30.05,30.07,29.99,30.11,29.93,30.24,
            30.04,29.99,29.95,30.12,30.00,rep(NA,41))
#convert barometric pressure to mm Hg
cobia$bp<-cobia$bp*25.4
#add water vapor pressure, from http://intro.chem.okstate.edu/1515sp01/database/vpwater.html
cobia$wvp<-c(22.4,28.3,22.4,22.4,22.4,22.4,22.4,22.4,22.4,22.4,22.4,22.4,22.4,22.4,
             28.3,35.7,28.3,35.7,28.3,35.7,28.3,35.7,28.3,28.3,35.7,28.3,35.7,28.3,35.7,28.3,35.7,28.3,35.7,35.7,rep(NA,41))
#calculate partial pressure of oxygen at 1 atmosphere (should be near 155 mmHg)
cobia$ppO2<-(cobia$bp-cobia$wvp)*0.2095
#convert to pO2 of water with measured oxygen concentrations, solubility
cobia$O2_pO2<-cobia$O2_mg.l/cobia$ppO2
#calculate Pcrit from Ccrit, units are mmHg
cobia$Pcrit<-cobia$Ccrit/cobia$O2_pO2

#get rid of unimportant columns
cobia<-subset(cobia, select=-c(O2_umol.l,O2_mg.l,bp,wvp,ppO2,O2_pO2))


#first need to add TL to dataframe
cobia$TL<-NA
cobia$TL[which(cobia$Animal_ID=="C08")]<-95
cobia$TL[which(cobia$Animal_ID=="C09")]<-88
cobia$TL[which(cobia$Animal_ID=="C11")]<-86
cobia$TL[which(cobia$Animal_ID=="C12")]<-105
cobia$TL[which(cobia$Animal_ID=="C13")]<-108
cobia$TL[which(cobia$Animal_ID=="C15")]<-84
cobia$TL[which(cobia$Animal_ID=="C16")]<-87
cobia$TL[which(cobia$Animal_ID=="C17")]<-88
cobia$TL[which(cobia$Animal_ID=="C18")]<-92
cobia$TL[which(cobia$Animal_ID=="C19")]<-91.5
cobia$TL[which(cobia$Animal_ID=="C20")]<-87.5
cobia$TL[which(cobia$Animal_ID=="C21")]<-82
cobia$TL[which(cobia$Animal_ID=="C22")]<-96
cobia$TL[which(cobia$Animal_ID=="C23")]<-83
cobia$TL[which(cobia$Animal_ID=="C24")]<-81
cobia$TL[which(cobia$Animal_ID=="C25")]<-88
cobia$TL[which(cobia$Animal_ID=="C26")]<-108
cobia$TL[which(cobia$Animal_ID=="C27")]<-88.5
cobia$TL[which(cobia$Animal_ID=="C28")]<-102
cobia$TL[which(cobia$Animal_ID=="C30")]<-81
cobia$TL[which(cobia$Animal_ID=="C31")]<-85
cobia$TL[which(cobia$Animal_ID=="C32")]<-81
cobia$TL[which(cobia$Animal_ID=="C33")]<-96
cobia$TL[which(cobia$Animal_ID=="C34")]<-89
cobia$TL[which(cobia$Animal_ID=="C35")]<-102.5


cobia<-cobia[with(cobia,order(Animal_ID,Test_Temp), decreasing=TRUE),] 


#look to see whether to include individuals that died or not
cobia$Alive_Dead<-ifelse(cobia$Animal_ID=="C31" & cobia$Test_Temp=="28", "dead",
                         ifelse(cobia$Animal_ID=="C24" & cobia$Test_Temp=="32", "dead",
                                ifelse(cobia$Animal_ID=="C32" & cobia$Test_Temp=="32", "dead",
                                       ifelse(cobia$Animal_ID=="C28" & cobia$Test_Temp=="32","dead","alive"))))

par(mfrow=c(2,2))
hist(cobia$MMR[which(cobia$Test_Temp=="32" & cobia$Alive_Dead=="alive")],main="MMR_32_alive",xlab="MMR",xlim=c(200,400))
hist(cobia$MMR[which(cobia$Test_Temp=="28" & cobia$Alive_Dead=="alive")],main="MMR_28_alive",xlab="MMR",xlim=c(150,400))
hist(cobia$MMR[which(cobia$Test_Temp=="32" & cobia$Alive_Dead=="dead")],main="MMR_32_dead",xlab="MMR",xlim=c(200,400))
hist(cobia$MMR[which(cobia$Test_Temp=="28" & cobia$Alive_Dead=="dead")],main="MMR_28_dead",xlab="MMR",xlim=c(150,400))

par(mfrow=c(2,2))
hist(cobia$SMR[which(cobia$Test_Temp=="32" & cobia$Alive_Dead=="alive")],main="SMR_32_alive",xlab="SMR",xlim=c(50,200))
hist(cobia$SMR[which(cobia$Test_Temp=="28" & cobia$Alive_Dead=="alive")],main="SMR_28_alive",xlab="SMR",xlim=c(50,200))
hist(cobia$SMR[which(cobia$Test_Temp=="32" & cobia$Alive_Dead=="dead")],main="SMR_32_dead",xlab="SMR",xlim=c(50,200))
hist(cobia$SMR[which(cobia$Test_Temp=="28" & cobia$Alive_Dead=="dead")],main="SMR_28_dead",xlab="SMR",xlim=c(50,200))

par(mfrow=c(2,2))
hist(cobia$AS[which(cobia$Test_Temp=="32" & cobia$Alive_Dead=="alive")],main="AS_32_alive",xlab="AS",xlim=c(50,240))
hist(cobia$AS[which(cobia$Test_Temp=="28" & cobia$Alive_Dead=="alive")],main="AS_28_alive",xlab="AS",xlim=c(50,240))
hist(cobia$AS[which(cobia$Test_Temp=="32" & cobia$Alive_Dead=="dead")],main="AS_32_dead",xlab="AS",xlim=c(50,240))
hist(cobia$AS[which(cobia$Test_Temp=="28" & cobia$Alive_Dead=="dead")],main="AS_28_dead",xlab="AS",xlim=c(50,240))
#looks like MMR, SMR, and even AS are lower for individuals that died during trial
#try to run stats on just cobia that survived trial


#instead of Zscore, multiple responses by constants so all responses on are similar scale
#for 10%, means: MMR-246, SMR-120, AS-127, Ccrit-2.07, Scrit-29,Pcrit-45, FS-2.1
cobia$MMRZ<-cobia$MMR * 1
cobia$SMRZ<-cobia$SMR * 2.1
cobia$ASZ<-cobia$AS * 1.9
cobia$FSZ<-cobia$FS * 117
cobia$ScritZ<-cobia$Scrit * 8.5
cobia$CcritZ<-cobia$Ccrit * 119
cobia$PcritZ<-cobia$Pcrit * 5.5

cobia1<-cobia[rep(seq_len(nrow(cobia)), 7), ]

cobia1<-cobia1[with(cobia1,order(Animal_ID), decreasing=TRUE),] 

#stacked responses
#scaled responses stacked
Animal_ID<-unique(cobia1$Animal_ID)
#for simp
biglist<-NULL
for(i in Animal_ID){
  a<-which(cobia1$Animal_ID==i)
  a1<-cobia1[a,]
  listR<-NULL
  for(j in 13:19){
    list<-a1[1:3,j]
    listR<-c(listR,list)
  }
  biglist<-c(biglist, listR)
}


#S for scaled
cobia1$responseS<-biglist



vars<-c(rep("MMR",3),rep("SMR",3),rep("AS",3),rep("FS",3),rep("Scrit",3),rep("Ccrit",3), rep("Pcrit",3))
cobia1$VAR<-rep(vars,25)

hist(cobia1$responseS)


#get rid of Scrit, Ccrit and FS
cobia1<-cobia1[-which(cobia1$VAR=="Scrit" | cobia1$VAR=="FS" | cobia1$VAR=="Ccrit"),]


hist(cobia1$responseS)#pretty normal to me, little skewed to right

cobia1$j<-rep(1:3, 100)


write.csv(cobia1, file="FINAL_cobia10_parts_sas.csv")

#turn MMR, SMR, and AS from dead cobia into NAs
cobia_no_dead<-cobia1
cobia_no_dead$responseS<-ifelse(cobia_no_dead$Alive_Dead=="dead",NA,cobia_no_dead$responseS)
hist(cobia_no_dead$responseS)#pretty normal to me, little skewed to right


write.csv(cobia_no_dead, file="FINAL_cobia10_no_dead_parts_sas.csv")




