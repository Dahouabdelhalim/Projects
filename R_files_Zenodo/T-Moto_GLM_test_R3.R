### Test1 #######################
## Test the difference in Moto between JP and CA subpopulations for
## the larval (S1), early juvenile (S2) and late juvenile (S3) stages
## accounting for the effect of temperature.

All<-read.csv('Moto_T_bin_re_outed_R3.csv',header = T)

hist(All$Moto.mean)#check distribution for determining family function in GLM
shapiro.test(All$Moto.mean)#check distribution
# no significant from normal distribution so using Gaussian family 


glm1<-glm(Moto.mean~Estimated.Temperature*(Region.ID), data=All,family=gaussian(link = "log"))
#considering reviewer's opinion, we use "log" to see if it is better to describe the trend
glm2<-glm(Moto.mean~Estimated.Temperature*(Region.ID), data=All,family=gaussian(link = "identity"))
glm3<-glm(Moto.mean~Estimated.Temperature+(Region.ID), data=All,family=gaussian(link = "identity"))
glm4<-glm(Moto.mean~Estimated.Temperature*(Region.ID)*Stage, data=All,family=gaussian(link = "identity"))
glm5<-glm(Moto.mean~Estimated.Temperature*(Region.ID)+Stage, data=All,family=gaussian(link = "identity"))

AIC(glm1, glm2, glm3, glm4, glm5)# model selection glm4 is the best (Supplementary Table 7)
summary(glm4) #Supplementary Table 8
sink("glm4_summary.txt")
summary(glm4)
sink()
par(mfrow=c(2,2))# model diagnostic plots (Supplementary Figure 9)
plot(glm4) # all look good

### Test1 Ends. ##############################################