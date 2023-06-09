#Read in your data files

#Assume you term your data file 'df'
#Assume variable of interest is called 'Accel.vectorial.sum'

### find the max and min for our data
max.data = max(df$Accel.vectorial.sum)
min.data = min(df$Accel.vectorial.sum)

#Create bins (from min to max, by 0.1)
breaks = seq(min.data, max.data, by=0.1)

#Cut the data according to bin size
cut.data = cut(df$Accel.vectorial.sum, breaks, right=FALSE)
#Frequency table
frequency = table(cut.data)

#Add up the frequencies in the table
cummul.freq=cumsum(frequency)
cummul.freq

### Calculate the Relative Frequency
relative.frequency=frequency/sum(frequency)
cf=as.data.frame(cummul.freq)
cf
cummul.freq=cf[,1]
cummul.freq
cummul.percentile = cummul.freq/max(cummul.freq)
cbind(frequency,relative.frequency,cummul.freq, cummul.percentile)
graph.cummul.perc = c(0, cummul.percentile)

#Plot
plot(breaks, graph.cummul.perc, ylab="Relative Cumulative Frequency", main="Relative Cumulative Frequency Graph")
lines(breaks, graph.cummul.perc)
abline(v = quantile(df$Accel.vectorial.sum, 0.95), col="red", lwd=3, lty=2) #0.95 quantile
abline(v = quantile(df$Accel.vectorial.sum, 0.99), col="blue", lwd=3, lty=2) #0.99 quantile

#more quantiles of data - user defined
x = c(0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.999, 0.9999)
y=rep(NA, length(x))

#Find values for each quantile as specified above in the vector termed 'x'
#Below code assumes you have 5 separate data files (e.g., individuals) and you want to calculate the above quantiles for each
#Results are saved in files termed L1 to L5
for(i in 1:length(x)){ y[i] = quantile(df$Accel.vectorial.sum, x[i])}
L1 = y
for(i in 1:length(x)){ y[i] = quantile(df2$Accel.vectorial.sum, x[i])}
L2 = y
for(i in 1:length(x)){ y[i] = quantile(df3$Accel.vectorial.sum, x[i])}
L3 = y
for(i in 1:length(x)){ y[i] = quantile(df4$Accel.vectorial.sum, x[i])}
L4 = y
for(i in 1:length(x)){ y[i] = quantile(df5$Accel.vectorial.sum, x[i])}
L5 = y

#Plot results of each file on one plot
plot(x, L1, type = "l", xlab = "Quantile", ylab = "Vector sum of acceleration (g)", ylim = c(1, 6)) #Change specified y-limits when necessary
points(x, L1)
lines(x, L2, col = "red") ; points(x, L2.q, col = "red")
lines(x, L3, col = "blue") ; points(x, L3.q, col = "blue")
lines(x, L4, col = "green") ; points(x, L1.q, col = "green")
lines(x, L5, col = "orange") ; points(x, L5.q, col = "orange")

