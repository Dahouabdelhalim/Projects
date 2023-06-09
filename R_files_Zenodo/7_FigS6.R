library(foreign)
library(sandwich)
#library(tidyverse)
#library(caret)

### 1. PREPARE DATA TABLE
dataCore=read.table("MASZK_dataCore.csv", sep=",", header = T)

#not=dataCore[dataCore$Vaccinated==0,c("Id","PfizerRating", "ModernaRating", "AstraRating", "SputRating", "SinoRating")]
not=dataCore[,c("Id","PfizerRating", "ModernaRating", "AstraRating", "SputRating", "SinoRating")]


not[is.na(not)]=0
not$s=not$PfizerRating+ not$ModernaRating + not$AstraRating + not$SputRating + not$SinoRating

not2=not[not$s>0,]

not2=not2[,2:6]



# Acceptance by the number of vaccines taken out

acc=not2
acc[acc<2]=0
acc[acc>0]=1
acc$s=acc$PfizerRating+ acc$ModernaRating + acc$AstraRating + acc$SputRating + acc$SinoRating

table(acc$s)

acc1=acc[acc$s==1,]


library(reshape2)

# co-acceptance matrix
acc_net=acc[acc$s>1 & acc$s<4,1:5]

names(acc_net)=c("Pfizer", "Moderna", "AstraZeneca", "Sputnik", "Sinopharm")

pflinks=melt(acc_net, id.vars=c("Pfizer"))
pflinks=pflinks[pflinks$Pfizer!=0,]
pflinks=aggregate(pflinks$value, by=list(pflinks$variable), FUN=sum)
pflinks$From="Pfizer"
pflinks=pflinks[,c(3,1,2)]
names(pflinks)=c("From", "To", "weight")

molinks=melt(acc_net, id.vars=c("Moderna"))
molinks=molinks[molinks$Moderna!=0,]
molinks=aggregate(molinks$value, by=list(molinks$variable), FUN=sum)
molinks$From="Moderna"
molinks=molinks[,c(3,1,2)]
names(molinks)=c("From", "To", "weight")

aslinks=melt(acc_net, id.vars=c("AstraZeneca"))
aslinks=aslinks[aslinks$AstraZeneca!=0,]
aslinks=aggregate(aslinks$value, by=list(aslinks$variable), FUN=sum)
aslinks$From="AstraZeneca"
aslinks=aslinks[,c(3,1,2)]
names(aslinks)=c("From", "To", "weight")

splinks=melt(acc_net, id.vars=c("Sputnik"))
splinks=splinks[splinks$Sputnik!=0,]
splinks=aggregate(splinks$value, by=list(splinks$variable), FUN=sum)
splinks$From="Sputnik"
splinks=splinks[,c(3,1,2)]
names(splinks)=c("From", "To", "weight")

silinks=melt(acc_net, id.vars=c("Sinopharm"))
silinks=silinks[silinks$Sinopharm!=0,]
silinks=aggregate(silinks$value, by=list(silinks$variable), FUN=sum)
silinks$From="Sinopharm"
silinks=silinks[,c(3,1,2)]
names(silinks)=c("From", "To", "weight")


acc_links=rbind(pflinks, molinks, aslinks, splinks, silinks)

library(igraph)
acc_g=graph_from_data_frame(acc_links, directed = F)

acc_g=simplify(acc_g, remove.multiple = T, edge.attr.comb = "min")

coord=layout_(acc_g, in_circle())

png("Acceptance on non-vaccinated.png", width=500, height=500)
plot(acc_g, layout=coord, edge.width=E(acc_g)$weight/2.5, edge.label=E(acc_g)$weight, edge.color="gray90", 
     edge.label.cex=2, vertex.label.cex=2)
dev.off()
