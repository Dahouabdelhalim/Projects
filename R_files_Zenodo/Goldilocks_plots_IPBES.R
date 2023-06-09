###################################################

LostYieldU=function(UvUMSY){
  phi=1.736
  RelYield=UvUMSY * (phi - UvUMSY * (phi - 1))^(1/(phi - 1))
  LostYield=1-RelYield
  return(LostYield)
}


LostYieldB=function(BvBMSY){
  phi=1.736
  phi2=(phi^((phi/(phi-1))))/(phi-1)
  BMSY=phi^(1/(1-phi))
  RelYield=BvBMSY*BMSY * phi2 - phi2*(BvBMSY*BMSY)^phi
  
  
    LostYield=1-RelYield
  LostYield[which(LostYield>1)]=1
  return(LostYield)
}







PlotStatus=function(Data){ #
  
  if (DataType=="B") PLost=LostYieldB(Data)
  if (DataType=="U") PLost=LostYieldU(Data)
 
  FishMore=array(dim=NY,0);FishLess=FishMore;JustRight=FishMore
  
  
  for (is in 1:Nstocks){
    for (iy in 1:NY){
      if (is.na(MSYT[is])==TRUE | is.na(PLost[iy,is])==TRUE) next
      if (PLost[iy,is]>1) {PLost[iy,is]=1}
      CCurrent=(1-PLost[iy,is])* MSYT[is]
      #if (is.na(CCurent)==TRUE)CCurent=0
      JustRight[iy]=JustRight[iy]+CCurrent
      if (Data[iy,is] >1) FishLess[iy]=FishLess[iy]+PLost[iy,is]*MSYT[is]
      if (Data[iy,is] <1) FishMore[iy]=FishMore[iy]+PLost[iy,is]*MSYT[is]
    }
    
  }
  TMSY=sum(MSYT,na.rm=TRUE)
  JustRight=JustRight/TMSY
  FishLess=FishLess/TMSY
  FishMore=FishMore/TMSY
  Tot=JustRight+FishLess+FishMore
  M=array(dim=c(NY,3))
  M[,1]=FishLess/Tot
  M[,2]=JustRight/Tot
  M[,3]=FishMore/Tot
  
  col=c("red","green","blue")
  legend=c("Fish Less","Just Right","Fish More")
  if (DataType=="B"){ #swap order for B
  M[,3]=FishLess/Tot
  M[,1]=FishMore/Tot
    col=c("red","green","blue")
    legend=c("Increase Abundance","Just Right","Decrease Abundance")
      }
  MM=M  #[which(Tot>.5),] plot all
  M=t(MM)
  

  
  names=Years #[which(Tot>0.5)]
  par(yaxt="n",mgp=c(2,1,0))
  hold=barplot(M,beside=FALSE,space=0,col=col,
               xlab="Year", ylab="Fraction of potential yield",names.arg=names)
  
  
  legend(x=40,y=.7,legend,fill=col,bg="white",cex=.5)
  #lines(hold,Tot,lwd=4,col="orange")
  axis(2, at=seq(0, 1, by=0.2), labels = FALSE)
  lablist=c("0.0","0.2","0.4","0.6","0.8","1.0")
    text(y=seq(0, 1, by=0.2), par("usr")[1], labels = lablist, srt = 0, pos = 2, xpd = TRUE)
 write.csv(file="LostYield.csv",data.frame(M))
    
} #end of PlotStatus function

##################################################################################

setwd("C:/Main/Research/IPBES 2022/Goldilocks figure")

load("RAMCore.RData")





#loop over and do all regions
#pdf(file="GoldiloxEurope.pdf")

RegionName="European Union"
RegionName="Atlantic Ocean"
RegionList=c("US Alaska","US East Coast","US Southeast and Gulf","US West Coast")
RegionName="US All Stocks"
RegionList=sort(unique(Region))
RegionName="All Assessed Stocks"
i=which(Region %in% RegionList )


Nstocks=length(StockID[i])

Years=seq(1950,2016)
NY=length(Years)
MSYT=MSY[i]

col=c("red","green","blue")

legend=c("Fish Less","Just Right","Fish More")
Header=paste(RegionName, "based on fishing mortality")
Data=UvU[,i]
DataType="U"
pdf(file="U.pdf",height=4,width=4)
par(mar=c(5,4,4,2))
PlotStatus(Data)
dev.off()

pdf(file="B.pdf",height=4,width=4)
Data=BvT[,i]
legend=c("Reduce Abundance","Just Right","Increase Abundance")
Header=paste(RegionName," based on abundance")
col=c("blue","green","red")
#col=c("red","green","blue")
DataType="B"
PlotStatus(Data)
dev.off()






