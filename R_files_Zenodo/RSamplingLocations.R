library(fields)
read.table("LocationsDec.csv",header=TRUE,sep="\\t")->Locat

#make map of USA
png("SamplingLocations.png",height=500, width=700)
US( xlim=c(-91, -70), ylim = c(37.5, 46), add=FALSE,shift=FALSE,col="grey",lwd=2)
for (i in 1:length(Locat[,1])){
points(-Locat$West[i],Locat$North[i],pch=16,cex=2)
}
for (i in c(1,3,4,5,6,7,9, 11, 12, 13, 14)){
text(-Locat$West[i],Locat$North[i]-0.4,Locat$Pop[i])}
i=2
text(-Locat$West[i]+0.2,Locat$North[i]+0.3,Locat$Pop[i])
i=8
text(-Locat$West[i]+0.2,Locat$North[i]-0.4,Locat$Pop[i])
i=10
text(-Locat$West[i]-0.8,Locat$North[i],Locat$Pop[i])
#text(-80,45.5,"Sampling Locations")
dev.off()
