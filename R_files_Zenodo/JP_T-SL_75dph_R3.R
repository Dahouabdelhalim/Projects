#### Test2 ##############################################
## Fit quadratic relationship between Standard length and temperature
## for JP sardine at 75 dph (end of early juvenile stage).

library(MASS)

All2<-read.csv('t-SL_moto.csv',header = T,stringsAsFactors=T) # Choose "t-SL_moto.csv"

model1=lm(L4~Tint4+I(Tint4**2),data=All2) 
model2=lm(L4~Tint4,data=All2) 

AIC(model1,model2)
model1.best=stepAIC(model1) # Model selection based on AIC

summary(model1) # Supplementary Table 10
sink("model1_summary.txt")
summary(model1)
sink()
#### Test2 Ends. ##############################################
