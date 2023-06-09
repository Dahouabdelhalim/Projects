
data <- read.csv("C:/Users/elizabeth/Dropbox/specificity field Hilltop/weightsforDIVtrimmedwithone headerrow.csv",header=T, row.names = 1)
names(data)
attach(data)
library(vegan)
library(labdsv)
library(MASS)

head(data)
colnames(data)
summary(data)#grand summary of all data

#totalnative
#data1 <- as.matrix(data[,2:65])

data1 <- as.matrix(data[,27:49])
class(data1)
dim(data1)
rownames(data1)
head(colnames(data1))
#data1[1:5,1:5]
apply(data1, 1, sum)
    data2<-decostand(data1, method="total")
apply(data2, 1, sum)
#Create Forest-stand subsets (AM, ECM, Mixed)
lam<-subset(data2[17:21,])
all<-subset(data2[1:6,])
st<-subset(data2[27:32,])
spi<-subset(data2[22:26,])
cla<-subset(data2[7:11,])
inf<-subset(data2[12:16,])
st


H1=diversity(data2,index="invsimpson")
lamdiv=diversity(lam,index="shannon")
stdiv=diversity(st,index="shannon")
lamdiv=diversity(lam,index="shannon")

H1


plot(H1)
points(stdiv)
