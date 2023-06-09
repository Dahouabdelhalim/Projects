sites <- read.csv("Rad_Pop_LATLON_CP.csv", header = T)

#origing GPS
anchor <- c(-89.98607, 38.41056)

#origin1: -84.56253, 38.41419
#origin2: -89.98607, 38.41056
#dist between these origins is 473.0025km

#pop indicator
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
sites <- sites[sites$pop %in% pops,]

#contains disthaversine function, returns distance in meters
library(geosphere)

#empty result file to store distances
results <- data.frame(pop1=character(), dist_origin=numeric(), stringsAsFactors=F)
for(i in 1:nrow(sites)){
  calcdist <- distHaversine(c(sites[i,3],sites[i,2]),anchor)
  results[nrow(results)+1,] <- c(pops[i],calcdist)
}
results$dist_origin <- as.numeric(results$dist_origin)
results$dist_origin <- results$dist_origin/1000

#empty result file to store distances
check <- read.csv("model2.csv")
for(i in 1:nrow(df)){
  if(df[i,3]>0){
    df[i,2] <- (df[i,2]+473.0025)
    }
}

df1 <- read.csv("model1.csv")
df2 <- read.csv("model2.csv")
df3 <- read.csv("model3.csv")
duo <- cbind(df1,df2)
trio <- cbind(duo, df3)