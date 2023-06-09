#Read data
Hindlimb.data <- read.csv("Metatarsal_estimate_inputdata_FINAL.csv")
View(Hindlimb.data)

#Create vectors for metatarsal length and width data
Meta3<-c(Hindlimb.data[,5])
Meta3width<-c(Hindlimb.data[,6])
plot(Meta3width,Meta3)

#Test normality of these data
shapiro.test(Meta3)$p.value
shapiro.test(Meta3width)$p.value

#Log data to ensure normality
Meta3 <- log(Meta3,10)
Meta3width <- log(Meta3width,10)
plot(Meta3width,Meta3)

#Calculate regression model
hindlimb.model = lm(Meta3 ~ Meta3width)
hindlimb.model
abline(hindlimb.model,col="red")
summary(hindlimb.model)
shapiro.test(hindlimb.model$residuals)
plot(hindlimb.model$fitted.values,hindlimb.model$residuals)

#Use model to estimate
Teleocrater.estimate.logged <- (hindlimb.model$coefficients[2] * log(Hindlimb.data[27,6],10))+hindlimb.model$coefficients[1]
Teleocrater.estimate.nonlogged <- 10^Teleocrater.estimate.logged
Teleocrater.estimate.nonlogged


