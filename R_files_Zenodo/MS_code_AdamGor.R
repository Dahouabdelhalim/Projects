library(multcomp)
library(multcompView)
library(lsmeans)
library(PerformanceAnalytics)

sphragis=read.table(".../MS_data_AdamGor.csv", sep=",", dec=".", header=T, na.strings = "NA")

sphragis$year=as.factor(sphragis$year)
table(sphragis$year, sphragis$app2)

sphragis=sphragis[sphragis$app2=="sphragis",]
sphragis=droplevels(sphragis)

#sl_mean
mod1=lm(sl_mean~year, data=sphragis)
summary(mod1)

#sh_mean
mod2=lm(sh_mean~year, data=sphragis)
summary(mod2)

#sw_mean
mod3=lm(sw_mean~year, data=sphragis)
summary(mod3)

#tukey
sphragis$year=as.factor(sphragis$year)
summary(glht(mod1, linfct=mcp(year="Tukey")))
summary(glht(mod2, linfct=mcp(year="Tukey")))
summary(glht(mod3, linfct=mcp(year="Tukey")))

plot(glht(mod1, linfct=mcp(year="Tukey")))
plot(glht(mod2, linfct=mcp(year="Tukey")))
plot(glht(mod3, linfct=mcp(year="Tukey")))

mod1_t=glht(mod1, linfct=mcp(year="Tukey"))
mod2_t=glht(mod2, linfct=mcp(year="Tukey"))
mod3_t=glht(mod3, linfct=mcp(year="Tukey"))
cld(mod1_t)
cld(mod2_t)
cld(mod3_t)


#plotting shield attributes yearly
par(mfrow=c(3,1))
par(mar = c(0, 4.5, 0.5, 0.5))
boxplot(sl_mean~year, sphragis, xlab="", ylab="Shield length [mm]", las=1, notch=T, xaxt="n", ylim=c(min(sphragis$sl_mean, na.rm=T), 14), cex.lab=1.5, cex.axis=1.2)
mtext("a", 3, at=1, line=-1.7)
mtext("a", 3, at=2, line=-1.7)
mtext("a", 3, at=3, line=-1.7)
mtext("a", 3, at=4, line=-1.7)
mtext("a", 3, at=5, line=-1.7)
mtext("a", 3, at=6, line=-1.7)
par(mar = c(0, 4.5, 0, 0.5))
boxplot(sh_mean~year, sphragis, xlab="", ylab="Shield height [mm]", las=1, notch=T, xaxt="n", ylim=c(min(sphragis$sh_mean, na.rm=T), 5.5), cex.lab=1.5, cex.axis=1.2)
mtext("a", 3, at=1, line=-1.7)
mtext("a", 3, at=2, line=-1.7)
mtext("ab", 3, at=3, line=-1.7)
mtext("b", 3, at=4, line=-1.7)
mtext("b", 3, at=5, line=-1.7)
mtext("ab", 3, at=6, line=-1.7)
par(mar = c(4, 4.5, 0, 0.5))
boxplot(sw_mean~year, sphragis, xlab="Year", ylab="Shield width [mm]", las=1, notch=T, ylim=c(min(sphragis$sw_mean, na.rm=T), 4.2), cex.lab=1.5, cex.axis=1.2)
mtext("bc", 3, at=1, line=-1.7)
mtext("c", 3, at=2, line=-1.7)
mtext("bc", 3, at=3, line=-1.7)
mtext("a", 3, at=4, line=-1.7)
mtext("ab", 3, at=5, line=-1.7)
mtext("a", 3, at=6, line=-1.7)

#correlation
cor_test=cbind(sphragis[["sl_mean"]], sphragis[["sh_mean"]], sphragis[["sw_mean"]])
colnames(cor_test)=c("Shield length [mm]", "Shield height [mm]", "Shield width [mm]")
chart.Correlation(cor_test, histogram = T, pch=16, method = "pearson")
