setwd("C:/Users/emchenoweth/Desktop/New folder/Step 3 Energy Costs")


library("doBy", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
library("plyr", lib.loc="~/R/win-library/3.5")

se <- function(x, na.rm = TRUE) {
  
  sd(x,na.rm = )/sqrt(sum(!is.na(x)))}

EnergyCosts<-read.csv("EnergyCosts2.csv")

AllLunges<-read.csv("AllLunges.csv")
AllDives<-read.csv("AllDives.csv")

Mc<-26942
AllLunges$meta.cost<-(0.5*Mc*(AllLunges$max.speed^2- AllLunges$min.speed^2)+Mc*2.4)/1000/0.16
EnergyCosts$meta.costs.mean<-NA
EnergyCosts$meta.costs.median<-NA
EnergyCosts$tag<-as.character(EnergyCosts$tag) 
AllLunges$tag<-as.character(AllLunges$tag) 
AllLunges<-AllLunges[!is.na(AllLunges$meta.cost),]


for (n in which(is.na(EnergyCosts$dur))){
  vec<-EnergyCosts[n,]
  sub<-AllDives[AllDives$tag == vec$tag &
                  AllDives$dive >= vec$dive.start &
                  AllDives$dive <= vec$dive.stop,]
  EnergyCosts[n,]$dur<- max(sub$sec.max, na.rm = TRUE) - min(sub$sec.min, na.rm = TRUE)
  EnergyCosts[n,]$dive.depth<- median(sub$p.max, na.rm = TRUE)
}


EnergyCosts$meta.costs.se<-NA
EnergyCosts$p.surface<-NA
#lunge costs for above 22.5 m lunges and below for whales without tag or speeds
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$meta.costs.mean<-mean(((0.5*Mc*(AllLunges[AllLunges$p < 22.5,]$max.speed^2-AllLunges[AllLunges$p< 22.5,]$min.speed^2))+Mc*2.4)/1000, na.rm = TRUE)/0.16
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$meta.costs.se<-se(((0.5*Mc*(AllLunges[AllLunges$p < 22.5,]$max.speed^2-AllLunges[AllLunges$p< 22.5,]$min.speed^2))+Mc*2.4)/1000, na.rm = TRUE)/0.16
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$meta.costs.median<-median(((0.5*Mc*(AllLunges[AllLunges$p < 22.5,]$max.speed^2-AllLunges[AllLunges$p< 22.5,]$min.speed^2))+Mc*2.4)/1000, na.rm = TRUE)/0.16

EnergyCosts[EnergyCosts$group == "KrillSS", ]$meta.costs.mean<-mean(0.5*26949*(AllLunges[AllLunges$tag == "mn12_262a",]$max.speed^2-AllLunges[AllLunges$tag == "mn12_262a",]$min.speed^2+2.4)/1000, na.rm = TRUE)/0.16
EnergyCosts[EnergyCosts$group == "KrillSS",]$meta.costs.se<-se(0.5*26949*(AllLunges[AllLunges$tag == "mn12_262a",]$max.speed^2-AllLunges[AllLunges$tag == "mn12_262a",]$min.speed^2+2.4)/1000, na.rm = TRUE)/0.16
EnergyCosts[EnergyCosts$group == "KrillSS",]$meta.costs.median<-median(0.5*26949*(AllLunges[AllLunges$tag == "mn12_262a",]$max.speed^2-AllLunges[AllLunges$tag == "mn12_262a",]$min.speed^2+2.4)/1000, na.rm = TRUE)/0.16

for (n in which(is.na(EnergyCosts$meta.costs.mean))){
  vec<-EnergyCosts[n,]
  sub<-AllLunges[AllLunges$tag == vec$tag &
                   AllLunges$dive >= vec$dive.start &
                   AllLunges$dive <= vec$dive.stop & !is.na(AllLunges$meta.cost),]
  EnergyCosts[n,]$lunge.depth<- mean(sub$p, na.rm = TRUE)
  EnergyCosts[n,]$max.speed<- mean(sub$max.speed, na.rm = TRUE)
  EnergyCosts[n,]$X..speed<-median(sub$delta.speed, na.rm = TRUE)
  #EnergyCosts[n,]$TOD.start<-times(format(min(as.POSIXlt(sub$DateTime, tz ="America/Anchorage"), na.rm = TRUE), "%H:%M:%S"))
  #EnergyCosts[n,]$TOD.end<-times(format(max(as.POSIXlt(sub$DateTime, tz ="America/Anchorage"), na.rm = TRUE), "%H:%M:%S"))
  EnergyCosts[n,]$meta.costs.median<-median(sub$meta.cost,na.rm = TRUE)
  EnergyCosts[n,]$meta.costs.mean<-mean(sub$meta.cost,na.rm = TRUE)
  EnergyCosts[n,]$meta.costs.se<-se(sub$meta.cost,na.rm = TRUE)
  EnergyCosts[n,]$p.surface<-nrow(sub[sub$p<22.5,])/nrow(sub)
}

BMR<-(3.6*Mc^0.74)/1000 #kJ/s
AllDives$mean.meta.cost_lunge<-NA
AllDives$meta.cost.lunging<-NA
AllDives$meta.cost_lunge.se<-NA
AllDives$p.dive.gliding<-NA
AllDives$meta.cost.gliding<-NA
AllDives$meta.cost.swimming <-NA

for (n in which(AllDives$lunges>0)){
  d<-AllDives[n,]$dive
  t<-AllDives[n,]$tag
  
  AllDives[n,]$mean.meta.cost_lunge<- mean(AllLunges[AllLunges$dive == d & AllLunges$tag == t & !is.na(AllLunges$meta.cost),]$meta.cost, na.rm = TRUE)
  AllDives[n,]$meta.cost.lunging<-AllDives[n,]$mean.meta.cost_lunge*AllDives[n,]$lunges
  AllDives[n,]$meta.cost_lunge.se<- se(AllLunges[AllLunges$dive == d & AllLunges$tag == t,]$meta.cost, na.rm = TRUE)
}
AllDives[AllDives$lunges==0 & !is.na(AllDives$lunges),]$meta.cost.lunging<-0
AllDives$p.dive.gliding<-0.859-(28.203/AllDives$p.max)#from williams/100 to make it a proportion
  AllDives[AllDives$p.dive.gliding <0 & !is.na(AllDives$p.dive.gliding),]$p.dive.gliding<-0
  AllDives$meta.cost.gliding<-AllDives$dive.time*AllDives$p.dive.gliding*1.4*BMR
  AllDives$meta.cost.swimming <-3*BMR*(AllDives$surface.time+
                                             (AllDives$dive.time-AllDives$dive.time*AllDives$p.dive.gliding)-
                                             AllDives$lunges*12.2) #time spent lunging

AllDives$total.meta.cost<-AllDives$meta.cost.lunging+AllDives$meta.cost.gliding+AllDives$meta.cost.swimming
AllDives$time.lunging<-AllDives$lunges*12.2
AllDives$time.swimming<-AllDives$surface.time+(AllDives$dive.time-AllDives$dive.time*AllDives$p.dive.gliding)-AllDives$time.lunging
AllDives$time.gliding<-AllDives$dive.time*AllDives$p.dive.gliding
AllDives$p.meta.cost.gliding<-AllDives$meta.cost.gliding/AllDives$total.meta.cost
AllDives$p.meta.cost.swimming<-AllDives$meta.cost.swimming/AllDives$total.meta.cost
AllDives$p.meta.cost.lunging<-AllDives$meta.cost.lunging/AllDives$total.meta.cost
AllDives$lunges_min<-AllDives$lunges/(AllDives$dur/60)

#ggplot()+geom_raster(aes(x = dive.p, y = lunges_min, fill = p.meta.cost.lunging), data = AllDives[!is.na(AllDives$p.max)&
                                                                                          #  !is.na(AllDives$lunges_min)&
                                                                                          #  !is.na(AllDives$p.meta.cost.lunging),])

EnergyCosts$meta.cost.gliding<- NA
EnergyCosts$meta.cost.swimming<- NA
EnergyCosts$meta.cost.lunging<- NA
EnergyCosts$total.meta.cost<- NA
EnergyCosts$meta.costs_min<- NA
EnergyCosts$p.time.lunging<- NA
EnergyCosts$p.time.gliding<- NA
EnergyCosts$p.time.swimming<- NA
EnergyCosts$lunge_min<- NA
#AllDives<-AllDives[,-(18:32)]

for (n in which(EnergyCosts$tag %in% levels(as.factor(AllDives$tag)))){
  vec<-EnergyCosts[n,]
  sub<-AllDives[AllDives$tag == vec$tag &
                  AllDives$dive >= vec$dive.start &
                  AllDives$dive <= vec$dive.stop &
                  !is.na(AllDives$lunges) ,]
  EnergyCosts[n,]$meta.cost.gliding<-sum(sub$meta.cost.gliding, rm.na = TRUE)
  EnergyCosts[n,]$meta.cost.swimming<-sum(sub$meta.cost.swimming, rm.na = TRUE)
  EnergyCosts[n,]$meta.cost.lunging<-sum(sub$meta.cost.lunging, rm.na = TRUE)
  EnergyCosts[n,]$total.meta.cost<-sum(sub$total.meta.cost, rm.na = TRUE)
  EnergyCosts[n,]$meta.costs_min<-EnergyCosts[n,]$total.meta.cost/(EnergyCosts[n,]$dur/60)
  EnergyCosts[n,]$p.time.lunging<-sum(sub$time.lunging, rm.na = TRUE)/(sum(sub$time.lunging, rm.na = TRUE)+sum(sub$time.gliding, rm.na = TRUE)+sum(sub$time.swimming, rm.na = TRUE))
  EnergyCosts[n,]$p.time.gliding<-sum(sub$time.gliding, rm.na = TRUE)/(sum(sub$time.lunging, rm.na = TRUE)+sum(sub$time.gliding, rm.na = TRUE)+sum(sub$time.swimming, rm.na = TRUE))
  EnergyCosts[n,]$p.time.swimming<-sum(sub$time.swimming, rm.na = TRUE)/(sum(sub$time.lunging, rm.na = TRUE)+sum(sub$time.gliding, rm.na = TRUE)+sum(sub$time.swimming, rm.na = TRUE))
}

NoTagSpeed<-which(EnergyCosts$meta.cost.lunging <2| is.na(EnergyCosts$meta.cost.lunging))
EnergyCosts[NoTagSpeed,]$meta.cost.lunging<-EnergyCosts[NoTagSpeed,]$lunges*EnergyCosts[NoTagSpeed,]$meta.costs.mean
EnergyCosts[NoTagSpeed,]$p.time.lunging<-EnergyCosts[NoTagSpeed,]$lunges*12.2/EnergyCosts[NoTagSpeed,]$dur
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$meta.cost.gliding<-0
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$p.time.gliding<-0
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$p.time.swimming<-(EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$dur-EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$lunges*12.2)/EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$dur
EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$meta.cost.swimming<-(EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$dur-EnergyCosts[EnergyCosts$tag %in% c("PJ", "Hobo"),]$lunges*12.2)*BMR*3
EnergyCosts[NoTagSpeed,]$total.meta.cost<-EnergyCosts[NoTagSpeed,]$meta.cost.lunging+EnergyCosts[NoTagSpeed,]$meta.cost.gliding+EnergyCosts[NoTagSpeed,]$meta.cost.swimming
EnergyCosts[NoTagSpeed,]$meta.costs_min<-EnergyCosts[NoTagSpeed,]$total.meta.cost/(EnergyCosts[NoTagSpeed,]$dur/60)

EnergyCosts$breaths_lunge<-EnergyCosts$breaths/EnergyCosts$lunges 
EnergyCosts$lunge_min<-EnergyCosts$lunges/(EnergyCosts$dur/60) 
EnergyCosts$breaths_min<-EnergyCosts$breaths/(EnergyCosts$dur/60)
EnergyCosts$min_breath<-(EnergyCosts$dur/60)/EnergyCosts$breaths
EnergyCosts[is.na(EnergyCosts$dive.depth),]$dive.depth<-EnergyCosts[is.na(EnergyCosts$dive.depth),]$lunge.depth

EnergyCosts$p.meta.cost.gliding<-EnergyCosts$meta.cost.gliding/EnergyCosts$total.meta.cost
EnergyCosts$p.meta.cost.lunging<-EnergyCosts$meta.cost.lunging/EnergyCosts$total.meta.cost
EnergyCosts$p.meta.cost.swimming<-EnergyCosts$meta.cost.swimming/EnergyCosts$total.meta.cost

EnergyCosts

############### Bonus Plots #################


nrow(AllLunges[AllLunges$tag== "mn13_108a" | AllLunges$tag== "mn13_109a" & AllLunges$p<22.5 & !is.na(AllLunges$p) , ])
nrow(AllLunges[AllLunges$tag== "mn13_108a" | AllLunges$tag== "mn13_109a" & AllLunges$p>=22.5, ])
nrow(AllLunges[AllLunges$tag== "mn13_111a" | AllLunges$tag== "mn13_112a"| AllLunges$tag== "mn13_113a"| AllLunges$tag== "mn13_113d" & AllLunges$p<22.5, ])
nrow(AllLunges[AllLunges$tag== "mn13_111a" | AllLunges$tag== "mn13_112a"| AllLunges$tag== "mn13_113a"| AllLunges$tag== "mn13_113d" & AllLunges$p>=22.5, ])

mean(AllLunges[AllLunges$tag== "mn13_108a" | AllLunges$tag== "mn13_109a" & AllLunges$p<22.5, ]$meta.cost)
mean(AllLunges[AllLunges$tag== "mn13_108a" | AllLunges$tag== "mn13_109a" & AllLunges$p>=22.5, ]$meta.cost)
mean(AllLunges[AllLunges$tag== "mn13_111a" | AllLunges$tag== "mn13_112a"| AllLunges$tag== "mn13_113a"| AllLunges$tag== "mn13_113d" & AllLunges$p<22.5, ]$meta.cost)
mean(AllLunges[AllLunges$tag== "mn13_111a" | AllLunges$tag== "mn13_112a"| AllLunges$tag== "mn13_113a"| AllLunges$tag== "mn13_113d" & AllLunges$p>=22.5, ]$meta.cost)

mean(EnergyCosts[EnergyCosts$tag== "mn13_108a" | EnergyCosts$tag== "mn13_109a" & EnergyCosts$lunge.depth<22.5, ]$max.speed)
mean(EnergyCosts[EnergyCosts$tag== "mn13_108a" | EnergyCosts$tag== "mn13_109a" & EnergyCosts$lunge.depth>=22.5, ]$max.speed)
mean(EnergyCosts[EnergyCosts$tag== "mn13_111a" | EnergyCosts$tag== "mn13_112a"| EnergyCosts$tag== "mn13_113a"| EnergyCosts$tag== "mn13_113d" & EnergyCosts$lunge.depth<22.5, ]$max.speed)
mean(EnergyCosts[EnergyCosts$tag== "mn13_111a" | EnergyCosts$tag== "mn13_112a"| EnergyCosts$tag== "mn13_113a"| EnergyCosts$tag== "mn13_113d" & EnergyCosts$lunge.depth>=22.5, ]$max.speed)



ggplot()+geom_point(aes(breaths_min,meta.costs_min, color = prey), data = EnergyCosts[-1,])+
  geom_smooth(aes(breaths_min,meta.costs_min, color = prey), data = EnergyCosts[-1,],method = "lm")
summary(step(lm(breaths_min~ prey+ dive.depth+meta.costs_min, data = EnergyCosts[-16,])))
ggplot()+geom_smooth(aes(breaths_lunge, max.speed), data = EnergyCosts, method = "lm")+geom_point(aes(breaths_lunge, max.speed, color = group.type, shape = prey), data = EnergyCosts)
ggplot()+geom_smooth(aes(breaths_lunge, ?..speed), data = EnergyCosts[-10,], method = "lm")+geom_point(aes(breaths_lunge, ?..speed, color = tag), data = EnergyCosts[-10,])
ggplot()+geom_smooth(aes(breaths_lunge, lunge_min), data = EnergyCosts, method = "lm")+geom_point(aes(breaths_lunge, lunge_min, color = tag), data = EnergyCosts)
ggplot()+geom_smooth(aes(breaths_lunge, meta.costs.median), data = EnergyCosts, method = "lm")+geom_point(aes(breaths_lunge, meta.costs.median, color =tag), data = EnergyCosts)
ggplot()+geom_point(aes(breaths_min, lunge_min, color = prey), size =3, data = EnergyCosts)#+geom_smooth(aes(breaths_min, lunge_min), data = EnergyCosts, method = "lm")
ggplot()+geom_smooth(aes(min_breath, dive.depth), data = EnergyCosts, method = "lm")+geom_point(aes(min_breath, dive.depth, color = prey), data = EnergyCosts[-10,])
ggplot()+geom_point(aes(meta.costs_min, breaths_min, color = tag), data = EnergyCosts[-10,])+
  geom_smooth(aes(meta.costs_min, breaths_min),method = "lm", data = EnergyCosts[-10,])
ggplot()+geom_point(aes(breaths_min, max.speed, color = tag), size= 3,data = EnergyCosts)

mod<-lm(breaths_lunge ~ meta.costs.mean, data = EnergyCosts)
mod<-lm(breaths_min ~ lunge_min+meta.costs.mean, data = EnergyCosts)
mod<-lm(meta.costs_min ~ lunge_min*max.speed , data = EnergyCosts[-10,])
mod<-lm(breaths_min ~ max.speed, data = EnergyCosts[-5,])

ggplot()+geom_raster(aes(dive.depth, lunge_min, fill = p.meta.cost.gliding), data = EnergyCosts)

ggplot()+geom_bar(aes(x=p, y = , value, fill = variable),stat = "identity", data = melt(EnergyCosts[c(12,26,27,25)]))
ggplot()+geom_bar(aes(x=dive.depth, y = , value, fill = variable),stat = "identity",data = melt(EnergyCosts[c(12,32,33,34)]))


summary(mod)
summary(step(mod))
AllLunges$depth.bin<-AllLunges$p < 22.5
summary(step(lm(max.speed~depth.bin+tag, data = AllLunges)))

mod<-lm(breaths_lunge ~?..speed +lunge_min +lunge_min^2 , data = EnergyCosts)
plot(EnergyCosts$?..speed,EnergyCosts$lunge_min)
plot(step(mod))
summary(step(mod))


mod<-lm(breaths_lunge ~max.speed+?..speed  , data = EnergyCosts)
summary(step(mod))



#write.csv(AllLunges,"F:/Analysis/Chapter3/Tags/AllLunges.csv")


ggplot()+geom_boxplot(aes(tag, delta.speed),data = AllLunges)


ggplot()+geom_point(aes(max.speed, delta.speed),data = AllLunges)

ggplot()+geom_point(aes(p, delta.speed, color = tag),data = AllLunges)+geom_smooth(aes(p, delta.speed),data = AllLunges, method = "lm")


