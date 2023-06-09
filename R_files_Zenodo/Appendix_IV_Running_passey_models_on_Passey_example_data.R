
source("Appendix II Passey EMeas Model 2021.r")
source("Appendix III Passey Inverse Model 2019.r")

# Example Data from Passey et al., 2005

nsolxns<-100
numtrials1<-100
Length1<-c(11,13,8,10,9,12,11,7,7,8,10,10)
dMeas1<-c(-7.93,-7.67,-6.54,-6.33,-5.41,-4.48,-5.31,-6.41,-7.79,-8.15,-8.34,-8.1)
finit1<-0.25
r1 = .05;
r2 = 1;
r3 = 5;
maxratio = -6.8;
minratio = -7.0;
stdev = .01;
df1 = 0.011;# this is the value that Passey et al. settle on and my code agrees with this
depth1<-rep(20,12)
la1<-50
lm1<-100
avelength1<-round(mean(Length1))
maxlength1<-20
minlength1 <- 2
mindepth1 <- 2

test<-PasseyEMeas1_1(Length1,dMeas1,la1,numtrials1)
hist(test[[1]])
test2<-PasseyInverse(Length1,dMeas1,depth1,finit1,la1,lm1,maxlength1,minlength1,mindepth1,df1,nsolxns)
mean(test2[[1]])

forplot<-data.frame(test2[[2]])

library(ggplot2)
p<-ggplot(forplot, aes(x=totallength,y=mean))+
  theme(panel.grid.major = element_line(size = 0.75, linetype = 'solid',colour = "white"))+
  xlab("Distance")+
  ylab("Delta 18O")+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="grey70")+
  geom_line(size=1,aes(x=totallength,y=mean),colour="black")+
  geom_point(size=2.5,aes(x=totallength,y=mean))
p



