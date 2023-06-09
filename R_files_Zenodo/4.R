### Temp coral disease analyses ###

data<-read.csv("4._Temp_cordis_8-29-18.csv") # temperature heat stress data + CARIBBEAN coral disease data 1970-2013 
head(data)

# rename
data$temp<-data$decimalFraction_locations_DHWge4_reefs_Had

# check

length(data$decimalFraction_locations_DHWge4_reefs_Had)
length(data$temp)

plot(data$decimalFraction_locations_DHWge4_reefs_Had)
plot(data$temp)

### Spearman's rho ###
cor.test(data$temp,data$proport_3yr,method="spearman",conf.level=0.95) 
# ties warning is because of 0s but OK

##### Plot ### 
# use proport_3yr_10x to be on same scale
head(data)

par(mar = c(5,5,2,5))
with(data, plot(year, temp, type="l", col="red3", 
             ylab="Temp anomaly index with DHWge4_reefs_Had"))

par(new = T)
with(data, plot(year, proport_3yr_10x, type="l", axes=F, xlab=NA, ylab=NA, cex=1.2))
axis(side = 4)
mtext(side = 4, line = 3, 'Disease reports/taxon * 10')
legend("topleft",
       legend=c("Temp anomaly index", "Disease reports/taxon"),
       lty=c(1,1), cex=0.8, col=c("red3", "black"))

