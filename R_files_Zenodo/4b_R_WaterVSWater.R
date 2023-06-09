rm(list=ls())
library (nlme)

##  Daten einlesen
setwd("D:/ConsumerDemandPaper/AuswertungWasserVsWater")
ww.df <- read.table ('data_WaterVsWater.txt', header= TRUE, sep= '\\t', dec= ',', as.is= FALSE)

str(ww.df)
attach(ww.df)
names(ww.df)


contrasts (ww.df [, 'Corner']) <- contr.sum (2)
ww.df [, 'Day'] <- factor (ww.df [, 'Day'])
contrasts (ww.df [, 'Day']) <- contr.sum (9)

summary (ww.df)


## Model
M4 <- lme (DrinkingEvents ~ Day * Corner, random= ~ 1 | Animal/Day, data= ww.df, method= 'ML')
summary (M4)
anova (M4, type= 'marginal')

qqnorm(residuals(M4)) 
plot(fitted(M4),residuals(M4))


## Annahmen
qqnorm (resid (wvw3))
qqnorm (unlist (ranef (wvw3, level= 1)))
qqnorm (unlist (ranef (wvw3, level= 2)))

scatter.smooth (fitted (wvw3), resid (wvw3))
boxplot (split (resid (wvw3), ww.df [, 'Corner']))
boxplot (split (resid (wvw3), ww.df [, 'Day']))

## SchÃ¤tzungen "Randeffekte"
ww.df <- cbind (ww.df, model.matrix (~ Day + Corner, ww.df) [, -1])
wvw.numContr <- lme (DrinkingEvents ~ Corner1 +
                         Day1 + Day2 + Day3 + Day4 + Day5 + Day6 + Day7 + Day8 +
                         Corner1:Day1 + Corner1:Day2 + Corner1:Day3 + Corner1:Day4 +
                         Corner1:Day5 + Corner1:Day6 + Corner1:Day7 + Corner1:Day8,
                     random= ~ 1 | Animal/Day, data= ww.df, method= 'ML')
summary (wvw.numContr)

#install.packages("contrast")
library (contrast)
Corner.estim <- contrast (wvw.numContr,
                          list (Corner1= c (-1, 1),
                                Day1= 0, Day2= 0, Day3= 0, Day4= 0,
                                Day5= 0, Day6= 0, Day7= 0, Day8= 0))

boxplot (DrinkingEvents ~ Corner, ww.df)
boxplot (DrinkingEvents ~ Day, ww.df)