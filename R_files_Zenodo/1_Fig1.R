# Network to matrices
setwd("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/package")

d=read.table("MASZK_dataCore.csv", sep=",", header = T)

edges=d[d$Vaccinated==1 & d$RejectedAny==1,c(5:10)]
e_c=data.frame(unique(edges$VaccineType),c("Moderna", "Sinopharm", "Pfizer", "AstraZeneca", "Sputnyk", "Janssen"))
names(e_c)=c("VaccineType", "vaccine")
edges=merge(edges, e_c, by="VaccineType")
#edges$VaccineType=edges$vaccine.y

library(reshape2)
e=melt(edges, id.vars = "VaccineType")
e=e[e$value==1,]
e$value=1

e_c2=data.frame(unique(e$variable),c("Pfizer", "Moderna", "AstraZeneca", "Sputnyk", "Sinopharm"))
names(e_c2)=c("variable", "vaccine")
e=merge(e, e_c2, by="variable")
names(e)[c(2,4)]=c("Target", "Source")
e=e[,c(4,2)]
e$val=1
e_to_gephi=aggregate(e$val, by=list(e$Source, e$Target), FUN=sum)
names(e_to_gephi)=c("Source", "Target", "Weight")


library(igraph)
net=graph_from_data_frame(e_to_gephi, directed = T)

plot(net)

write.table(e_to_gephi, file = "e_to_gephi.csv", col.names = T, row.names = F, sep=",")

