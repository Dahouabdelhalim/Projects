

# This code allows to calculate the averaged values of the climate variables across neighbour 
# cells within a radius of 1º to 8º. E.g., from BIO1 (local value) we obtained from D1BIO1 
# (1º around) to D8BIO1 (8º around).



# ----------- Data loading and preparation

# Open Data_III and copy all


read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\\t",header=header,...)
}

Data=read.excel()


Dist <- function(stanD3,X,Y){
  return (round(sqrt((stanD3$lon - X)^2 + (stanD3$lat - Y)^2)))
}

Neighbours <- function(idparcela,Data,mindist,maxdist){
  lon <- Data$lon
  lat <- Data$lat
  dist <- Dist(Data[Data$Key == idparcela,],lon,lat)
  neigh <- Data[which(dist > mindist & dist < maxdist),]
  return(neigh)
}


# For example, to obtain from D4BIO1 to D4BIO19 (4º around)

D=4

D4matz <- matrix(nrow=nrow(Data),ncol=19)

colnames(D4matz)<-c("D4BIO1","D4BIO2","D4BIO3","D4BIO4","D4BIO5","D4BIO6","D4BIO7",
                    "D4BIO8","D4BIO9","D4BIO10","D4BIO11","D4BIO12","D4BIO13",
                    "D4BIO14","D4BIO15","D4BIO16","D4BIO17","D4BIO18","D4BIO19")

D4matz <- as.data.frame(D4matz)

for(j in 4:22){ # 19 worldclim variables

  for(i in 1:nrow(Data)){

  dists <- Dist(Data[i,],Data$lon[-i],Data$lat[-i])
  a<-Data$Key[i]
  neighbours <- Neighbours(a,Data,0,D+1)
  value<- neighbours[j]
  m<-mean(value[,1])
  m<-ifelse(is.na(m),Data[i,j],m)
  D4matz[i,(j-3)]<-round(m,2)

  }
 print(j)
 
 }

#------------------------------------------------------------------

# To joint from D1BIOs to D8BIOs

Data<-cbind(Data,D1matz,D2matz,D3matz,D4matz,D5matz,D6matz,D7matz,D8matz)

write.csv(Data, 'DBIOsdata.csv')

