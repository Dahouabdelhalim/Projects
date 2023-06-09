data <- read.csv("massChecksChicksCross-foster.csv", stringsAsFactors=FALSE, sep=";")

# remove excess columns
data <- data[,which(!(colnames(data) %in% c("X","X.1","X.2","X.3","X.4","X.5","X.6")))]

# get dates
dates <- colnames(data)
which.cols.dates <- grep("X", dates)

#dates <- dates[which.cols.dates]
dates <- strptime(dates, format="X%d.%m.%y")

# convert hatch date into dates
data$hatchdate <- strptime(data$hatchdate, format="%d.%m.%y")

data2 <- data.frame()

for (i in 1:nrow(data)) {
for (j in which.cols.dates) {
	
	if (!is.na(data[i,j])) {
		age <- as.numeric(dates[j] - data$hatchdate[i] + 1)
		weight <- data[i,j]
		data2 <- rbind(data2, data.frame(ring=data$RING[i], age=age, weight=weight, stringsAsFactors=FALSE))
	}
}
}

data2$age2 <- data2$age^2

plot(data2$age,data2$weight, pch=20, xlab="Age", ylab="Weight", xlim=c(15,50) ,ylim = c(6,18))

abline(h=10.5, lty=2)
abline(v=18, lty=2)

mod <- lm(weight ~ age + age2, data=data2)
newx <- 18:45
lines(newx, predict(mod, newdata=list(age=newx, age2=newx^2)), col="red", lwd=2)


