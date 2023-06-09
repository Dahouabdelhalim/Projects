rm(list=ls())

## Allgemeines
setwd ('')
if (file.exists ('.RData')) load ('.RData')

library (lmerTest)
library (boot)

## Daten einlesen
ICT <- read.table ('data_time spent in IC.txt', header= TRUE, sep= '\\t', dec= ',', as.is= FALSE)
names(ICT)
ICT [, 'Day'] <- ICT [, 'Day'] - 7
contrasts (ICT [, 'Run']) <- contr.sum (3)
ICT <- ICT [!is.na (ICT [, 'sum_IC_time']),]
summary (ICT)

boxplot (sum_IC_time ~ Run + Day, ICT)

## Modell
M1 <- lmer (sum_IC_time ~ Day * Run + (1|Tag/Run), data=ICT)
summary (M1)
anova (M1, type= 'marginal')


qqnorm(residuals(M1)) 
plot(fitted(M1),residuals(M1))

## Bootstrap
extract.ci <- function (x) { ## Hilfsfunktion
    out <- data.frame (numeric (0),
                       numeric (0),
                       numeric (0))
    for (i in 1:length (x [['t0']])) {
        out <- rbind (out, c (x [['t0']] [i],
                              boot.ci (x, index= i,
                                       type= 'perc') [['percent']] [, 4:5]))
    }
    names (out) <- c ('estim', 'lo.ci', 'up.ci')
    out
}

extract.var <- function (mod, part) {
    vars.df <- as.data.frame (VarCorr (mod))
    vars.df [vars.df [, 'grp'] == part, 'vcov']
}

extract.var (M1, 'Tag')
extract.var (M1, 'Run:Tag')
extract.var (M1, 'Residual')

Q.Perc <- function (x) 100 * extract.var (x, 'Tag') / (extract.var (x, 'Tag') +
                                                       extract.var (x, 'Run:Tag') +
                                                       extract.var (x, 'Residual'))
estimQP.raw <- bootMer (M1, Q.Perc, nsim= 1000, .progress= 'win')
estimQP.val <- extract.ci (estimQP.raw)
estimQP.val

Q.Quot <- function (x) extract.var (x, 'Tag') / (extract.var (x, 'Run:Tag') +
                                                 extract.var (x, 'Residual'))
estimQQ.raw <- bootMer (M1, Q.Quot, nsim= 1000, .progress= 'win')
estimQQ.val <- extract.ci (estimQQ.raw)
estimQQ.val