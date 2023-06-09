## flower size

ms.field<-c(21,22,18,17,13,15,14,22,18,17,21,19,18,20,15,18,20,17,20,21,18,18,15,15,13,18,18,18,19,23,18,21,21,18,17,20,19,17,20,26)
hh.field<-c(20,24,26,25,21,26,18,25,23,27,19,24,21,27,25,19,28,26,21,25,25,21,26,20,24,27,21,22,24,22,26,25,27,29,25,25,16,27,24,25)

par(mfrow=c(1,2))
hist(hh.field, xlab='diameter (mm)', ylab="Freq", main="Hermaphrodite", col='lightblue', xlim=c(10,30))
hist(ms.field, xlab='diameter (mm)', ylab="", main="Female", col='pink', xlim=c(10,30))

t.test(hh.field, ms.field)