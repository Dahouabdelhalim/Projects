#### SEX RATIOS  ####

setwd("path")
###- 1. sex ratio in the wild fish:  #
sr <- read.delim("sr_wild2019.txt", header = TRUE, dec = ",")
attach(sr)
head(sr)
 library(lme4)
 y<-cbind(male,fem)
 sr1<-glmer(y~species-1+(1|year)+(1|site),family="binomial")
 summary(sr1)
 confint(sr1)
 # Random effects:
 #   Groups Name        Variance Std.Dev.
 # site   (Intercept) 0.173580 0.41663 
 # year   (Intercept) 0.005327 0.07298 
 # Number of obs: 373, groups:  site, 133; year, 7
 # 
 # Fixed effects:
 #           Estimate Std. Error z value Pr(>|z|)    
 # speciesF -0.72738    0.06221 -11.693  < 2e-16 ***
 # speciesK -0.07814    0.11177  -0.699    0.484    
 # speciesO -0.48968    0.06412  -7.637 2.22e-14 ***
 # speciesP -0.80958    0.09658  -8.382  < 2e-16 ***
 
 detach(sr)
 
 ### 2. captive fish sex ratio 
 sr <- read.delim("~/Dropbox/Documents/HALANCICI/MSS KILLI/Sex_Diffs/sexratio_at_birth_Notho.txt")
 names(sr)
 attach(sr)
 m <- lm(Sex_ratio ~ species - 1)
 summary(m)
 plot(m)
 library(effects)
 plot(allEffects(m))
 
 
 #### Lipofuscin #####
 
 lipofus <- read.delim("lipofus.txt")
 attach(lipofus)
 head(lipofus)
 str(lipofus)
 
 # EXPLORATORY PART
 dotchart(tl)
 dotchart(partAll)
 densityplot(partAll)
 densityplot(log(partAll))
 histogram(partAll)
 std.tl<-scale(tl)
 
 # lipo_A <- glmer(partAll ~ sex * age + std.tl + (1 | pop) + (1 | indID), family= poisson, data = lipofus)
 lipo_A <- glmer(partAll ~ sex * age + (1 | pop) + (1 | indID), family= poisson, data = lipofus)
 summary(lipo_A)
 anova(lipo_A)
 plot(allEffects(mod=lipo_A))
 plot(Effect(focal.predictor="sex",mod=lipo_A))
 plot(lipo_A)
 
 #validate model:
 E <-resid(lipo_A)
 F <- fitted(lipo_A)
 plot(F, E, cex = 0.7, pch = 16, xlab= "fitted values", ylab = "Residuals")
 abline(h = 0, lty = 2)
 hist(resid(lipo_A))
 # OK
 
 detach(lipofusc)
 
### kidney ###
 histology <- read.delim("kidTum.txt")
 attach(histology)
 head(histology)
 
 require(lme4)
 require(effects)
 require(lattice)
 require(lmerTest)
 
 # EXPLORATORY PART
 dotchart(kidtumor)
 densityplot(kidtumor)
 histogram(kidtumor)
 std.tl<-scale(TL)
 
 #### ================= Kidney proliferation analysis ================
 # new analysis with tumors and sections:
 kt1 <- glmer(kidtumor ~  sex * age + (1 | pop), family= poisson, data = histology)
 summary(kt1)
 plot(allEffects(mod=kt1))
 #  effect of sex (0.010)
 
 plot(Effect(focal.predictor="sex",mod=kt1))
 plot(Effect(focal.predictor="age",mod=kt1))
 # sex and age approaches singificance
 
 plot(kt1)
 # testing for overdispersion formally:
 E1 <- resid(kt1, type = "pearson")
 N  <- nrow(histology)
 p  <- length(coef(kt1))
 sum(E1^2) / (N - p)  # overdoispersion 1.15 - OK!
 
 
 #### ================= Liver proliferation analysis ================
 histology <- read.delim("histology.txt")
 # new analysis with tumors and sections:
 lt1 <- glmer(liverTumor ~  sex * age + (1 | pop), family= poisson, data = histology)
 summary(lt1)
 plot(allEffects(mod=lt1))
 #  effect of sex (0.010)
 
 plot(Effect(focal.predictor="sex",mod=lt1))
 plot(Effect(focal.predictor="age",mod=lt1))
 # sex and age approaches singificance
 
 plot(kt1)
 
 