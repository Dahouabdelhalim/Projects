library(geosphere)
library(prevR)
library(rangeBuilder)
library(recluster)
library(rworldmap)
library(rworldxtra)
library(rgeos)
library(plotrix)
source("biodecrypt.R")
map <- getMap(resolution = "low")

# Create an example for a dataset
mat<-rbind(cbind(rnorm(n = 25, mean = 1, sd = 4),rnorm(n = 25, mean = 40, sd = 3)),cbind(rnorm(n = 25, mean = 7, sd = 5),rnorm(n = 25, mean = 45, sd = 2)))

id<-c(rep(1,25),rep(2,25))
id[sample(c(1:50))[1:10]]<-0

#Apply the function to view hulls with custum parameters
# alpha increases to include 95% of default cases and hulls are almost convex 
biodecrypt.view(mat,id, clipToCoast="no")

# With a lower fraction hulls becomes more concave. 
# Excluded dots works as a punctiform sub-hull in the attribution. 
biodecrypt.view(mat,id, alpha=c(1,5), fraction=0.80, clipToCoast="no")


# Make the separation with default parameters 
attribution<-biodecrypt(mat,id, clipToCoast="no")

#plot the results
plot(mat,type="n")
plot.biodecrypt(attribution, attributed="fade")

#Group plots into pies of 5x5 squares
plot(mat,type="n")
plot.biodecrypt(attribution, square=5, minsize=0.8, attributed="fade")

# Make the separation with custom parameters
attribution<-biodecrypt(mat,id, alpha=c(1,5), buffer=20, ratio=2, fraction=0.80, minimum=3, clipToCoast="no")

#plot the results
plot(mat,type="n")
plot.biodecrypt(attribution, attributed="fade")


# Make the cross-validation with default parameters
cross<-biodecrypt.cross(mat,id, clipToCoast="no")
plot(mat,type="n")
plot.biodecrypt(cross)

#Optimize parameters (a few cases and few runs to make it shorter)
wrap_data<-biodecrypt.wrap(mat,id, alpha=c(1,10), ratio=c(2,4), buffer=c(20,50), runs=5, clipToCoast="no")
parameters<-biodecrypt.optimise(wrap_data$table)
#inspect the optimised parameters
parameters

# Use MIR+NUR+NIR instead of MIR^2+NUR+NIR and a higher threshold
parameters<-biodecrypt.optimise(wrap_data$table, coef=c(1,1,1), penalty=20)
#inspect the optimised parameters
parameters

#Run biodecrypt with optimised parameters
attribution<-biodecrypt(mat,id, alpha=rep(parameters$alpha, max(id)), ratio=parameters$ratio,buffer=parameters$buffer,clipToCoast="no")
#plot the results
plot(mat,type="n")
plot.biodecrypt(attribution, attributed="fade")

# Open the dataset for butterflies
data<-read.table("Total_data.txt", sep="\\t",h=T)
head(data)

#Select the genus
genus<-"Melitaea"
mydata<-data[which(data[,4]==genus),]

# Extract the variables
longlat<-mydata[,c(6:5)]
id<-mydata[,2]

# Preview of hulls
biodecrypt.view(longlat, id,map=map)

#Run the analysis with the default parameters
separation_default<-biodecrypt(longlat, id)
plot(map, xlim=c(-8,25),ylim=c(35,70))
plot.biodecrypt(separation_default,fading=60, attributed="fade")

# Or with custom parameters
separation<-biodecrypt(longlat, id,alpha=c(10.4,10.4), buffer=31.7, ratio=2.7)
plot(map, xlim=c(-8,25),ylim=c(35,70))
plot.biodecrypt(separation_default,fading=60, attributed="fade")
#save the results as a txt file
write.table(separation$table,"Results.txt")

# Run the cross validation with default parameters
cross_default<-biodecrypt.cross(longlat, id)
plot(map, xlim=c(-8,25),ylim=c(35,70))
plot.biodecrypt(cross_default)

# Run the wrap analysis to optimize the parameters 
it can take hours, we reduced the number of combinations and imposed runs=5 instead of the default runs=10)
wrap<-biodecrypt.wrap(longlat,id,alpha=c(1,5),ratio=c(2,5),buffer=c(40,160), checkdist=T, xlim=NULL,ylim=NULL,main=NULL,save=F, runs=5)

# Inspect values
wrap

# Compute the optimized parameters
parameters<-biodecrypt.optimise(wrap$table,penalty=10)


# Make the separation with optimized parameters
separation<-biodecrypt(longlat, id, alpha=rep(parameters$alpha, max(id)), ratio=parameters$ratio,buffer=parameters$buffer) 
plot(map, xlim=c(-8,25),ylim=c(35,70))
plot.biodecrypt(separation,fading=60, attributed="fade")


# Make the cross-validation with optimized parameters
cross<-biodecrypt.cross(longlat, id, alpha=rep(parameters$alpha,max(id)), ratio=parameters$ratio,buffer=parameters$buffer) 
plot(map, xlim=c(-8,25),ylim=c(35,70))
plot.biodecrypt(cross)
#Inspect MIR, NIR and NUR
cross$MIR
cross$NIR
cross$NUR


