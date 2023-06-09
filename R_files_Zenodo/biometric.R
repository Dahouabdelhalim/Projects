attach(morfo)
rm <- morfometria
rm <- variegatusvariegatus

variegatus <- B_variegatus_JB_comcoordenadas
shapiro.test(variegatus[HBL])


sexo
sexo <- factor(c("F", "M", "M", "F"))
levels(sexo)
classe <- factor(c("AD", "JV", "AD", "JV"))
levels(classe)

shapiro.test(HBL)
shapiro.test(peso)

wilcox.test(peso[sexo=="F"], peso[sexo=="M"])
wilcox.test(peso[classe=="AD"], peso[classe=="JV"])

wilcox.test(HBL[sexo=="F"], HBL[sexo=="M"])
wilcox.test(HBL[classe=="AD"], HBL[classe=="JV"])

summary(HBL[sexo=="F"])
summary(HBL[sexo=="M"])
summary(peso[sexo=="F"])
summary(peso[sexo=="M"])

library(Rcmdr)
