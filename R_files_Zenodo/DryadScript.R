###import data

r2ed=read.csv("\\\\soundbaser2.csv")
library(metafor) #####start metafor package


sounddatamod=r2ed#[which(r2ed$Tot=="Study"),] ####separates summary lines for graphing from raw data
str(sounddatamod)
sounddatamod$STBin=factor(sounddatamod$STBin)
sounddatamod$yi=as.numeric(as.character(sounddatamod$yi))
sounddatamod$vi=as.numeric(as.character(sounddatamod$vi))

soundsplit=split(sounddatamod, sounddatamod$Response) ###split up sound data for each individual response

options(stringsAsFactors = FALSE)


####minimum frequency

#####create object containing moderator variables
modmatmin=as.matrix(cbind(STBin=as.factor(soundsplit$`Min Frequency`$STBin), Taxa=as.factor(soundsplit$`Min Frequency`$Taxa)))


result.min=rma.mv(yi=yi, V=vi, method="REML", measure= "SMD", mods = modmatmin, random=~1|as.factor(AutN),data=soundsplit$`Min Frequency`, sei=soundsplit$`Min Frequency`$se)
result.min
print.rma.mv(result.min)
summary(result.min)


hist(soundsplit$`Min Frequency`$yi)
hist(soundsplit$`Max Frequency`$yi)
hist(soundsplit$Amplitude$yi)
hist(soundsplit$Duration$yi)
hist(soundsplit$F0$yi)
qqnorm(result.min)


confint(result.min)

summary(result.min)

###amplitude
modmatamp=as.matrix(cbind(STBin=soundsplit$Amplitude$STBin, AutN=soundsplit$Amplitude$AutN, Taxa=soundsplit$Amplitude$Taxa))

result.amplitude=rma(yi=yi, vi=vi, method="REML", measure= "SMD", mods = modmatamp, data=soundsplit$Amplitude, sei=se)
result.amplitude
qqnorm(result.amplitude)

###Call Rate
modmatcr=as.matrix(cbind(STBin=soundsplit$Rate$STBin, AutN=soundsplit$Rate$AutN, Taxa=soundsplit$Rate$Taxa))

result.rate=rma(yi=yi, vi=vi, method="REML", measure= "SMD", mods = modmatcr, data=soundsplit$Rate, sei=se)
result.rate

qqnorm(result.rate)

###duration
modmatduration=as.matrix(cbind(STBin=soundsplit$Duration$STBin, AutN=soundsplit$Duration$AutN, Taxa=soundsplit$Duration$Taxa))

result.duration=rma(yi=yi, vi=vi, method="REML", measure= "SMD", mods = modmatduration, data=soundsplit$Duration, sei=se)
result.duration

qqnorm(result.duration)

###f0
modmatf0=as.matrix(cbind(STBin=soundsplit$F0$STBin, AutN=soundsplit$F0$AutN, Taxa=soundsplit$F0$Taxa))

result.f0=rma(yi=yi, vi=vi, method="REML", measure= "SMD", mods = modmatf0, data=soundsplit$F0, sei=se)
result.f0

qqnorm(result.f0)


###maxfreq
modmatmax=as.matrix(cbind(STBin=soundsplit$`Max Frequency`$STBin, AutN=soundsplit$`Max Frequency`$AutN, Taxa=soundsplit$`Max Frequency`$Taxa))

result.max=rma(yi=yi, vi=vi, method="REML", measure= "SMD", mods = modmatmax, data=soundsplit$`Max Frequency`, sei=se)
result.max

qqnorm(result.max)

