library(pwt8)

sazonal <- ts(B_variegatus_JB_comcoordenadas, frequency = 12, start = 2015.5)
View(sazonal)
sazonal

shapiro.test(sazonal)
#nao normal

decomp <-decompose(sazonal)

plot(decomp)

plot(decomp$seasonal)

library(seastests)
summary(wo(sazonal))

plot(sazonal, xlab="Time", ylab="Captures")
text.default(2015.9, 6, labels = "WO-Test = 0
             p-values > 0.5", cex = .8)

births <- scan("http://robjhyndman.com/tsdldata/data/nybirths.dat")

View(births)
births
