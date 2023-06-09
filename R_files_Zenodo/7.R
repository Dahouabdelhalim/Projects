#########################################################
##		      	Sharing Carbon Responsibility 		       ##
##												                             ##
##						          By Meng					               ##
##					            SPPM,THU				               ##
##												                             ##
#########################################################
##Part 1.Making Environment
#create folders
setwd( "/Users/limeng/Desktop/CARBONFLOW")
path<- getwd()
pathio<- paste(path,"/","IOTABLE",sep="")
pathresult<- paste(path,"/","RESULT_20190212",sep="")
#library packages
library(readxl)
library(xlsx)
library(xlsxjars)
library(Matrix)
library(XLConnect)
#set dimension
nc<- 41
ni1<- 35
ni2<- 56
ni<- 34
nk<- 5
ny<- 21
nci<- nc*ni
nci1<- nc*ni1
nci2<- nc*ni2
nck<- nc*nk


##Part 2.Import Original Database
#2.1 Import Merging Data
flnm<- "/Users/limeng/Desktop/CARBONFLOW/OTHER/nmerge.xlsx"
ms35<- read_xlsx(flnm,sheet = 4,col_names = FALSE)
ms56<- read_xlsx(flnm,sheet = 5,col_names = FALSE)
ms35<-as.matrix(ms35)
m35<- bdiag(ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35,ms35)
ms56<-as.matrix(ms56)
m56<- bdiag(ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56,ms56)
  
#2.2 Import IO tables
fd<- matrix(0,nc,1)

for (i in 1:ny) {
  year<- 1994+i
  #get wiod
  if(year<=1999){
    #load data
    flnm<- paste(pathio,"/","WIOT13","/","wiot",year,"_row_apr14.csv",sep="")
    wiot<- read.csv(flnm,header = FALSE)
    wiot<- apply(wiot,2,as.numeric )
    wiot<- as.matrix(wiot)
    wiot[is.nan(wiot)]<- 0
    wiot[is.na(wiot)]<- 0
    #merge into 34 sectors
    m1<- nci1; m2<- m1+1; m3<- m1+8; m4<- nci1+nck+1
    wiot<- rbind(t(t(wiot[1:m1,])%*%m35),wiot[m2:m3,])
    wiot<- cbind(wiot[,1:m1]%*%m35,wiot[,m2:m4])
  }else if(year<=2014){
    #load data
    flnm<- paste(pathio,"/","WIOT16","/","WIOT",year,"_October16_ROW.RData",sep="")
    IO<- load(flnm)
    wiot<- wiot[,-(1:5)]
    wiot<- apply(wiot,2,as.numeric )
    wiot<- as.matrix(wiot)
    wiot[is.nan(wiot)]<- 0
    wiot[is.na(wiot)]<- 0
    #merge data into 40+1 countries
    i1<- 7;  a1_1<- (i1-1)*56; a1_2<- (i1-1)*56+1; a1_3<- i1*56; a1_4<- i1*56+1
    i2<- 19; a2_1<- (i2-1)*56; a2_2<- (i2-1)*56+1; a2_3<- i2*56; a2_4<- i2*56+1
    i3<- 33; a3_1<- (i3-1)*56; a3_2<- (i3-1)*56+1; a3_3<- i3*56; a3_4<- i3*56+1
    i4<- 44; a4_1<- (i4-1)*56; a4_2<- (i4-1)*56+1; a4_3<- i4*56; a4_4<- i4*56+1; a4_5<- i4*56+8
    #into less rows
    wiot_row_1<- wiot[a1_2:a1_3,]+wiot[a2_2:a2_3,]+wiot[a3_2:a3_3,]+wiot[a4_2:a4_3,]
    wiot_other<- wiot[a4_4:a4_5,]
    wiot<- rbind(wiot[1:a1_1,],wiot[a1_4:a2_1,],wiot[a2_4:a3_1,],wiot[a3_4:a4_1,],wiot_row_1,wiot_other)
    #into less cols
    #first int into less cols
    wiot_row_2<- wiot[,a1_2:a1_3]+wiot[,a2_2:a2_3]+wiot[,a3_2:a3_3]+wiot[,a4_2:a4_3]
    wiot_int<- cbind(wiot[,1:a1_1],wiot[,a1_4:a2_1],wiot[,a2_4:a3_1],wiot[,a3_4:a4_1],wiot_row_2)
    #then final demand into less cols
    a5_1<- 2464+(i1-1)*5; a5_2<- 2464+(i1-1)*5+1; a5_3<- 2464+i1*5; a5_4<- 2464+i1*5+1
    a6_1<- 2464+(i2-1)*5; a6_2<- 2464+(i2-1)*5+1; a6_3<- 2464+i2*5; a6_4<- 2464+i2*5+1
    a7_1<- 2464+(i3-1)*5; a7_2<- 2464+(i3-1)*5+1; a7_3<- 2464+i3*5; a7_4<- 2464+i3*5+1
    a8_1<- 2464+(i4-1)*5; a8_2<- 2464+(i4-1)*5+1; a8_3<- 2464+i4*5; a8_4<- 2464+i4*5+1
    wiot_row_3<- wiot[,a5_2:a5_3]+wiot[,a6_2:a6_3]+wiot[,a7_2:a7_3]+wiot[,a8_2:a8_3]
    wiot_fd<- cbind(wiot[,2465:a5_1],wiot[,a5_4:a6_1],wiot[,a6_4:a7_1],wiot[,a7_4:a8_1],wiot_row_3,wiot[,a8_4])
    #finally get the 41country wiot
    wiot<- cbind(wiot_int,wiot_fd)
    #merge into 34 sectors
    m1<- nci2; m2<- m1+1; m3<- m1+8; m4<- nci2+nck+1
    wiot<- rbind(t(t(wiot[1:m1,])%*%m56),wiot[m2:m3,])
    wiot<- cbind(wiot[,1:m1]%*%m56,wiot[,m2:m4])
  }else{
    flnm<- paste(pathio,"/","ADB","/","ADB_MRIO_",year,".xlsx",sep="")
    wiot<- read_xlsx(flnm,range = "E8:CSC2220",col_names = FALSE)
    wiot<- apply(wiot,2,as.numeric )
    wiot<- as.matrix(wiot)
    wiot[is.nan(wiot)]<- 0
    wiot[is.na(wiot)]<- 0
    #merge data into 40+1 countries
    i1<- 7;  a1_1<- (i1-1)*35; a1_2<- (i1-1)*35+1; a1_3<- i1*35; a1_4<- i1*35+1
    i2<- 19; a2_1<- (i2-1)*35; a2_2<- (i2-1)*35+1; a2_3<- i2*35; a2_4<- i2*35+1
    i3<- 33; a3_1<- (i3-1)*35; a3_2<- (i3-1)*35+1; a3_3<- i3*35; a3_4<- i3*35+1
    i4<- 44; a4_1<- (i4-1)*35; a4_2<- (i4-1)*35+1; a4_3<- i4*35; a4_4<- i4*35+1
    i5<- 63; a5_1<- i5*35+1; a5_2<- i5*35+8
    #into less rows
    wiot_row_1<- matrix(0,35,2521)
    for (i in c(7,19,33,44:63)) {
      a_1<- (i-1)*35; a_2<- (i-1)*35+1; a_3<- i*35; a_4<- i*35+1
      wiot_row_1<- wiot_row_1+wiot[a_2:a_3,]
    }
    wiot_other<- wiot[a5_1:a5_2,]
    wiot<- rbind(wiot[1:a1_1,],wiot[a1_4:a2_1,],wiot[a2_4:a3_1,],wiot[a3_4:a4_1,],wiot_row_1,wiot_other)
    #into less columns
    #int into less cols
    wiot_row_2<- matrix(0,1443,35)
    for (i in c(7,19,33,44:63)) {
      a_1<- (i-1)*35; a_2<- (i-1)*35+1; a_3<- i*35; a_4<- i*35+1
      wiot_row_2<- wiot_row_2+wiot[,a_2:a_3]
    }
    wiot_int<- cbind(wiot[,1:a1_1],wiot[,a1_4:a2_1],wiot[,a2_4:a3_1],wiot[,a3_4:a4_1],wiot_row_2)
    #final demand into less cols
    a5_1<- 2205+(i1-1)*5; a5_2<- 2205+(i1-1)*5+1; a5_3<- 2205+i1*5; a5_4<- 2205+i1*5+1
    a6_1<- 2205+(i2-1)*5; a6_2<- 2205+(i2-1)*5+1; a6_3<- 2205+i2*5; a6_4<- 2205+i2*5+1
    a7_1<- 2205+(i3-1)*5; a7_2<- 2205+(i3-1)*5+1; a7_3<- 2205+i3*5; a7_4<- 2205+i3*5+1
    a8_1<- 2205+(i4-1)*5; a8_2<- 2205+(i4-1)*5+1; a8_3<- 2205+i4*5; a8_4<- 2205+i4*5+1
    wiot_row_3<- matrix(0,1443,5)
    for (i in c(7,19,33,44:63)) {
      a_1<- 2205+(i-1)*5; a_2<- 2205+(i-1)*5+1; a_3<- 2205+i*5; a_4<- 2205+i*5+1
      wiot_row_3<- wiot_row_3+wiot[,a_2:a_3]
    }
    wiot_fd<- cbind(wiot[,2206:a5_1],wiot[,a5_4:a6_1],wiot[,a6_4:a7_1],wiot[,a7_4:a8_1],wiot_row_3,wiot[,a8_4])
    #finally get the 41country wiot
    wiot<- cbind(wiot_int,wiot_fd)
    #merge into 34 sectors
    m1<- nci1; m2<- m1+1; m3<- m1+8; m4<- nci1+nck+1
    wiot<- rbind(t(t(wiot[1:m1,])%*%m35),wiot[m2:m3,])
    wiot<- cbind(wiot[,1:m1]%*%m35,wiot[,m2:m4])
  }
  dim(wiot) 
  
  #get final demand
  m1<- nci+1
  m2<- nci+nck
  wfd<- wiot[1:nci,m1:m2]
  sumwfd<- colSums(wfd)
  swfd<- matrix(sumwfd,nc,nk,byrow = TRUE)
  dim(swfd)
  nfd<- swfd[,1]+swfd[,2]+swfd[,3]
  fd<- cbind(fd,nfd)
}

fd<- fd[,-1]

  

##Part 6.Save the results
#6.1
flnm<- paste(pathresult,"/","finaldemand.csv",sep="")
write.csv(fd,file = flnm)




