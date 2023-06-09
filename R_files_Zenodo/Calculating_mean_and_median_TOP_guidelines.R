library(ggplot2)
library(reshape)
library(stringr)

# IMPORT
top=read.csv("C:/Users/ipatarc/Desktop/PR_TOP/top-factor -v33.csv")

# # Policies of interest 
# 


policies=which(str_detect(colnames(top),"score"))

filtered.j=top[,policies]
per_score_policy=rbind(apply(filtered.j,2,function(x){length(which(x==1))}),
                       apply(filtered.j,2,function(x){length(which(x==2))}),
                       apply(filtered.j,2,function(x){length(which(x==3))}))

round(apply(filtered.j,2,mean,na.rm=T),1)
round(apply(filtered.j,2,median,na.rm=T),1)
apply(filtered.j,2,table)


mean(unlist(filtered.j),na.rm=T)
median(unlist(filtered.j),na.rm=T)

# what is mean if the journal adopts the guideline?


round(apply(filtered.j,2, function(x){
   return(mean(x[x!=0],na.rm=T))}),1)
      

mean(filtered.j[filtered.j!=0],na.rm=T)

round(apply(filtered.j,2, function(x){
  return(median(x[x!=0],na.rm=T))}),1)


median(filtered.j[filtered.j!=0],na.rm=T)