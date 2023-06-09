setwd("D:/Florence Jaffrezic/Professionnel/Documents/CATHY/Avril_2022")

data_croissance = read.table("data_croissance_analyses.txt",sep="\\t",header=TRUE)

library(nlme);library(lattice);library(car);library(lsmeans)

options(contrasts = c("contr.sum", "contr.poly")) 

# post-weaning and pubertal

data1 = data_croissance[data_croissance$time<=11,]

xyplot(mesures~time|group,data=data1,groups=as.character(animal),type="l",
	main="Period 1")

model1_1 = lme(mesures ~ group + time + group:time, random= ~1|animal,
data=data1,correlation=corAR1(form=~1|animal),na.action=na.omit)

summary(model1_1)
Anova(model1_1,type=3)

model2_1 = lme(mesures ~ 0 + group:as.factor(time), random= ~1|animal,
data=data1,correlation=corAR1(form=~1|animal),na.action=na.omit)

lsmeans(model2_1,list(pairwise~group|time)) 

# fattening period

data2 = data_croissance[data_croissance$time>=12,]

xyplot(mesures~time|group,data=data2,groups=as.character(animal),type="l",
	main="Period 2")

model1_2 = lme(mesures ~ group + time + group:time, random= ~1|animal,
data=data2,correlation=corAR1(form=~1|animal),na.action=na.omit)

summary(model1_2)
Anova(model1_2,type=3)

model2_2 = lme(mesures ~ 0 + group:as.factor(time), random= ~1|animal,
data=data2,correlation=corAR1(form=~1|animal),na.action=na.omit)

lsmeans(model2_2,list(pairwise~group|time)) 







