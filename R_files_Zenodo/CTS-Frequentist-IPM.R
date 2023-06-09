### Preliminary frequentist IPM and PVA script for the California Tiger Salamander (Ambystoma californiense)
### Christopher A. Searcy
### 2019

## This script creates the transition matrix and calculates the population growth rate (lambda) when the population is at low density
# Discretize body size continuum
x<-seq(-0.85,5.2,by=0.05)
# Define all of the parameters in the demographic functions
a<-(-0.423)
b<-0.564
c<-2.22
d<-0.588
e<-3.09
f<-2.25
g<-(-2.2)
h<-1.79
i<-0.11
j<-0.257
k<-0.629
l<-(-1.7)
m<-3.6
n<-0.021
o<-0.148
p<-3.3
q<-(-13.1)
r<-4.51
s<-1100
tt<-814
# Define the demographic functions for juvenile/adult and metamorph survival
sa<-1-exp(k*(x-2.47)^2+l*(x-2.47)+m)/((1/n)*(1/(1-o))+exp(k*(x-2.47)^2+l*(x-2.47)+m))
sm<-1-exp(f*(x-2.28)^2+g*(x-2.28)+h)/((1/i)*(1/(1-j))+exp(f*(x-2.28)^2+g*(x-2.28)+h))
# Create empty matrices to receive the metamorph and juvenile/adult growth functions
gm<-array(0,dim=c(length(x),length(x)))
ga<-array(0,dim=c(length(x),length(x)))
# Create vectors that describe the variation around these growth functions
y<-c(0.001,0.001,0.001,0.002,0.004,0.005,0.008,0.011,0.015,0.02,0.026,0.033,0.04,0.048,0.054,0.06,0.065,0.068,0.069,0.068,0.065,0.06,0.054,0.048,0.04,0.033,0.026,0.02,0.015,0.011,0.008,0.005,0.004,0.002,0.001,0.001,0.001)
z<-c(0.001,0.001,0.002,0.003,0.005,0.008,0.011,0.016,0.023,0.03,0.039,0.048,0.058,0.066,0.073,0.077,0.078,0.077,0.073,0.066,0.058,0.048,0.039,0.03,0.023,0.016,0.011,0.008,0.005,0.003,0.002,0.001,0.001)
# Loop filling the metamorph growth matrix
for(v in 32:122)#for a given starting size class, what weight does it correspond to? Around expected size class assign normal distribution ofexpected outcomes
{for(w in 1:37)
{if(u<-(round((a*(v/20-3.29)^2+b*(v/20-3.29)+c)*20))+w-1){gm[u,v]<-y[w]}}}
# Loop filling the juvenile/adult growth matrix
for(v in 1:122)
{for(w in 1:33)
{if(u<-(round((e+d*(v/20-3.85))*20))+w+1){ga[u,v]<-z[w]}}}
# Define the fecundity function, which is the product of maturity, fertility, percent females breeding, percent population that is female, and maximum possible embryonic/juvenile survival
fe<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*(s*(x-3.48)+tt)*0.341*0.5*0.0917
for(w in 1:122)
{if(fe[w]<0){fe[w]<-0}}
# Set the body size distribution of new metamorphs
Phi<-dnorm(x, mean = 2.24, sd = 0.345)
# Create empty matrices that will be the quadrants of the overall transition matrix
uj<-array(0,dim=c(length(x),length(x)))
uk<-array(0,dim=c(length(x),length(x)))
ul<-array(0,dim=c(length(x),length(x)))
up<-array(0,dim=c(length(x),length(x)))
uq<-array(0,dim=c(length(x),length(x)))
ur<-array(0,dim=c(length(x),length(x)))
bl<-array(0,dim=c(length(x),length(x)))
br<-array(0,dim=c(length(x),length(x)))
for(r in 1:122)
{uj[ ,r]<-gm[ ,r]*sm[r]}
for(r in 1:122)
{uk[r, ]<-uj[r, ]*fe[r]}
# Matrix ul will be the upper left quadrant of the transition matrix
for(r in 1:122)
{ul[ ,r]<-(Phi/sum(Phi))*sum(uk[ ,r])}
for(r in 1:122)
{up[ ,r]<-ga[ ,r]*sa[r]}
for(r in 1:122)
{uq[r, ]<-up[r, ]*fe[r]}
# Matrix ur will be the upper right quadrant of the transition matrix
for(r in 1:122)
{ur[ ,r]<-(Phi/sum(Phi))*sum(uq[ ,r])}
# Matrix bl will be the bottom left quadrant of the transition matrix
for(r in 1:122)
{bl[ ,r]<-gm[ ,r]*sm[r]}
# Matrix br will be the bottom right quadrant of the transition matrix
for(r in 1:122)
{br[ ,r]<-ga[ ,r]*sa[r]}
# Create the transition matrix
gt<-array(0,dim=c(2*length(x),2*length(x)))
for(i in 1:122)
{for(j in 1:122)
{gt[i,j]<-ul[i,j]}}
for(i in 1:122)
{for(j in 123:244)
{gt[i,j]<-ur[i,j-122]}}
for(i in 123:244)
{for(j in 1:122)
{gt[i,j]<-bl[i-122,j]}}
for(i in 123:244)
{for(j in 123:244)
{gt[i,j]<-br[i-122,j-122]}}
# Calculate lambda (principal eigenvector of transition matrix)
Output<-eigen(gt)
Output$values[1]
# Define number of years for simulation (T)
T<-100
# Create vector (N) to store population size distribution for each year of simulation
N<-array(NA,dim=c(T,2*length(x)))
# Simulation begins with one individual in each size class
N[1, ]<-1
# Simulation proceeds by applying transition matrix to current population size distribution
for(t in 1:(T-1))
{N[t+1, ]<-gt%*%N[t, ]}

## This script is used to calculate 95% confidence intervals around the population growth rate (lambda) when the population is at low density
# Create a vector to store the lambda values
Record<-rep(NA,1000)
# Create 1000 transition matrices and calculate population growth rate (lambda) for each one
for(G in 1:1000)
{x<-seq(-0.85,5.2,by=0.05)
# Now rather than being defined as a single value, each parameter in the demographic models is a random draw from its parametric distribution
a<-rnorm(1,-0.423,0.186)
b<-rnorm(1,0.564,0.0526)
c<-rnorm(1,2.22,0.126)
d<-rnorm(1,0.588,0.0415)
e<-rnorm(1,3.09,0.124)
f<-rnorm(1,2.25,0.581)
g<-rnorm(1,-2.2,0.225)
h<-rnorm(1,1.79,0.514)
i<-rnorm(1,0.11,0.022)
j<-rnorm(1,0.257,0.00946)
k<-rnorm(1,0.629,0.247)
l<-rnorm(1,-1.7,0.19)
m<-rnorm(1,3.6,0.478)
n<-rnorm(1,0.021,0.0065)
o<-rnorm(1,0.148,0.0114)
p<-rnorm(1,3.3,0.812)
q<-rnorm(1,-13.1,0.746)
r<-rnorm(1,4.51,1.98)
s<-rnorm(1,1100,87.4)
tt<-rnorm(1,814,305)
sa<-1-exp(k*(x-2.47)^2+l*(x-2.47)+m)/((1/n)*(1/(1-o))+exp(k*(x-2.47)^2+l*(x-2.47)+m))
sm<-1-exp(f*(x-2.28)^2+g*(x-2.28)+h)/((1/i)*(1/(1-j))+exp(f*(x-2.28)^2+g*(x-2.28)+h))
gm<-array(0,dim=c(length(x),length(x)))
ga<-array(0,dim=c(length(x),length(x)))
y<-c(0.001,0.001,0.001,0.002,0.004,0.005,0.008,0.011,0.015,0.02,0.026,0.033,0.04,0.048,0.054,0.06,0.065,0.068,0.069,0.068,0.065,0.06,0.054,0.048,0.04,0.033,0.026,0.02,0.015,0.011,0.008,0.005,0.004,0.002,0.001,0.001,0.001)
z<-c(0.001,0.001,0.002,0.003,0.005,0.008,0.011,0.016,0.023,0.03,0.039,0.048,0.058,0.066,0.073,0.077,0.078,0.077,0.073,0.066,0.058,0.048,0.039,0.03,0.023,0.016,0.011,0.008,0.005,0.003,0.002,0.001,0.001)
for(v in 32:106)
{for(w in 1:37)
{if(u<-(round((a*(v/20-3.29)^2+b*(v/20-3.29)+c)*20))+w-1){gm[u,v]<-y[w]}}}
for(v in 32:106)
{for(w in 1:33)
{if(u<-(round((e+d*(v/20-3.85))*20))+w+1){ga[u,v]<-z[w]}}}
fe<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*(s*(x-3.48)+tt)*0.341*0.5*0.0917
for(w in 1:122)
{if(fe[w]<0){fe[w]<-0}}
Phi<-dnorm(x, mean = 2.24, sd = 0.345)
uj<-array(0,dim=c(length(x),length(x)))
uk<-array(0,dim=c(length(x),length(x)))
ul<-array(0,dim=c(length(x),length(x)))
up<-array(0,dim=c(length(x),length(x)))
uq<-array(0,dim=c(length(x),length(x)))
ur<-array(0,dim=c(length(x),length(x)))
bl<-array(0,dim=c(length(x),length(x)))
br<-array(0,dim=c(length(x),length(x)))
for(w in 1:122)
{uj[ ,w]<-gm[ ,w]*sm[w]}
for(w in 1:122)
{uk[w, ]<-uj[w, ]*fe[w]}
for(w in 1:122)
{ul[ ,w]<-(Phi/sum(Phi))*sum(uk[ ,w])}
for(w in 1:122)
{up[ ,w]<-ga[ ,w]*sa[w]}
for(w in 1:122)
{uq[w, ]<-up[w, ]*fe[w]}
for(w in 1:122)
{ur[ ,w]<-(Phi/sum(Phi))*sum(uq[ ,w])}
for(w in 1:122)
{bl[ ,w]<-gm[ ,w]*sm[w]}
for(w in 1:122)
{br[ ,w]<-ga[ ,w]*sa[w]}
gt<-array(0,dim=c(2*length(x),2*length(x)))
for(w in 1:122)
{for(v in 1:122)
{gt[w,v]<-ul[w,v]}}
for(w in 1:122)
{for(v in 123:244)
{gt[w,v]<-ur[w,v-122]}}
for(w in 123:244)
{for(v in 1:122)
{gt[w,v]<-bl[w-122,v]}}
for(w in 123:244)
{for(v in 123:244)
{gt[w,v]<-br[w-122,v-122]}}
Output<-eigen(gt)
Record[G]<-Output$values[1]}
# Sort the list of lambda values in order to identify the 95% confidence intervals (25th and 976th values)
sort(Record)

## This script describes density-dependent population dynamics and is used to calculate carrying capacity and stable size distribution
x<-seq(-0.85,5.2,by=0.05)
# Demographic parameters are again defined as a single value
a<-rnorm(1,-0.423,0)
b<-rnorm(1,0.564,0)
c<-rnorm(1,2.22,0)
d<-rnorm(1,0.588,0)
e<-rnorm(1,3.09,0)
f<-rnorm(1,2.25,0)
g<-rnorm(1,-2.2,0)
h<-rnorm(1,1.79,0)
i<-rnorm(1,0.11,0)
j<-rnorm(1,0.257,0)
k<-rnorm(1,0.629,0)
l<-rnorm(1,-1.7,0)
m<-rnorm(1,3.6,0)
n<-rnorm(1,0.021,0)
o<-rnorm(1,0.148,0)
p<-rnorm(1,3.3,0)
q<-rnorm(1,-13.1,0)
r<-rnorm(1,4.51,0)
s<-rnorm(1,1100,0)
t<-rnorm(1,814,0)
u<-rnorm(1,-4.64,0)#intercept larv survival
v<-rnorm(1,-0.64,0)#slope larv survival (or swith w/ v)
w<-rnorm(1,-0.272,0)#slope mean meta body size (density effect size from Searcy et al. 2015, back-transformed = -0.10 ln(mean metamorph mass) 95% CI = (-0.16,-0.04))
z<-rnorm(1,2.27,0)##intercept mean meta body size
sa<-1-exp(k*(x-2.47)^2+l*(x-2.47)+m)/((1/n)*(1/(1-o))+exp(k*(x-2.47)^2+l*(x-2.47)+m))
sm<-1-exp(f*(x-2.28)^2+g*(x-2.28)+h)/((1/i)*(1/(1-j))+exp(f*(x-2.28)^2+g*(x-2.28)+h))
gm<-array(0,dim=c(length(x),length(x)))
ga<-array(0,dim=c(length(x),length(x)))
v1<-c(0.001,0.001,0.001,0.002,0.004,0.005,0.008,0.011,0.015,0.02,0.026,0.033,0.04,0.048,0.054,0.06,0.065,0.068,0.069,0.068,0.065,0.06,0.054,0.048,0.04,0.033,0.026,0.02,0.015,0.011,0.008,0.005,0.004,0.002,0.001,0.001,0.001)
v2<-c(0.001,0.001,0.002,0.003,0.005,0.008,0.011,0.016,0.023,0.03,0.039,0.048,0.058,0.066,0.073,0.077,0.078,0.077,0.073,0.066,0.058,0.048,0.039,0.03,0.023,0.016,0.011,0.008,0.005,0.003,0.002,0.001,0.001)
for(rep in 32:122)
{for(rep1 in 1:37)
{if(index<-(round((a*(rep/20-3.29)^2+b*(rep/20-3.29)+c)*20))+rep1-1){gm[index,rep]<-v1[rep1]}}}
for(rep in 1:122)
{for(rep1 in 1:33)
{if(index<-(round((d*(rep/20-3.85)+e)*20))+rep1+1){ga[index,rep]<-v2[rep1]}}}
fe<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*(s*(x-3.48)+t)*0.341*0.5
for(rep in 1:122)
{if(fe[rep]<0){fe[rep]<-0}}
# Set the number of years of simulation (T)
T<-100
# Create matrices to store the body size distribution of juveniles/adults (A), metamorphs (M), and the total popuplation (N) in each year
N<-array(NA,dim=c(T,length(x)))
A<-array(NA,dim=c(T,length(x)))
M<-array(NA,dim=c(T,length(x)))
# The simulation begins with one juvenile/adult in each size class
A[1, ]<-1
M[1, ]<-0
N[1, ]<-0
for(time in 1:(T-1))
# The population size distribution of juveniles/adults is equal to the number of metamorphs the survived and grew from the previous time step plus the number of juveniles/adults that survived and grew from the previous time step
{A[time+1, ]<-t(gm%*%(sm*M[time, ]))+t(ga%*%(sa*A[time, ]))
# The mean of the metamorph body size distribution is now density dependent, with density calculated as the ratio of the number of eggs produced (fecundity times the body size distribution of adults) to the pond volume (here set at 101000 cubic meters for Olcott Lake)
Phi<-dnorm(x,mean = z+w*(log(sum(fe*A[time+1, ])/101000)-1.23), sd = 0.227)
# The population size distribution of metamorphs is equal to the metamorph body size distribution (Phi/sum(Phi)) times the number of eggs produced times embryonic/larval survival, which is now density dependent
M[time+1, ]<-(Phi/sum(Phi))*sum(fe*A[time+1, ])*min(exp(u)*exp(v*(log(sum(fe*A[time+1, ])/101000)-3.38)),0.0917)
# The total population size distribution is the sum of the metamorph and juvenile/adult population size distributions
N[time+1, ]<-A[time+1, ]+M[time+1, ]}
# This calculates the carrying capacity
sum(N[100, ])
# Create a vector to store the population size in each year of the simulation
trajectory<-numeric(100)
# Fill the vector to illustrate the trajectory by which the population approahces carrying capacity (exhibits slight overshoot)
for(rep in 1:100)
{trajectory[rep]<-sum(N[rep, ])}
# m1 is the population size distribution of one-year old individuals when the population is at carrying capacity
m1<-t(gm%*%(sm*M[100, ]))
# a1 is the population size distribution of mature one-year olds when the population is at carrying capacity
a1<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m1
# mx is the population size distribution of x-year olds for individuals that did not reach maturity before age x
m2<-t(ga%*%t(sa*(m1-a1)))
# ax is the population size distribution of individuals who reach maturity at age x when the population is at carrying capacity
a2<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m2
m3<-t(ga%*%t(sa*(m2-a2)))
a3<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m3
m4<-t(ga%*%t(sa*(m3-a3)))
a4<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m4
m5<-t(ga%*%t(sa*(m4-a4)))
a5<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m5
m6<-t(ga%*%t(sa*(m5-a5)))
a6<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m6
m7<-t(ga%*%t(sa*(m6-a6)))
a7<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m7
m8<-t(ga%*%t(sa*(m7-a7)))
a8<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m8
m9<-t(ga%*%t(sa*(m8-a8)))
a9<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*m9
# gx is the population size distribution of x-year old individuals when the population is at carrying capacity
g1<-t(gm%*%(sm*M[100, ]))
# ax is the population size distribution of mature x-year old individuals when the population is at carrying capacity
a1<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g1
g2<-t(ga%*%t(sa*g1))
a2<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g2
g3<-t(ga%*%t(sa*g2))
a3<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g3
g4<-t(ga%*%t(sa*g3))
a4<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g4
g5<-t(ga%*%t(sa*g4))
a5<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g5
g6<-t(ga%*%t(sa*g5))
a6<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g6
g7<-t(ga%*%t(sa*g6))
a7<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g7
g8<-t(ga%*%t(sa*g7))
a8<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g8
g9<-t(ga%*%t(sa*g8))
a9<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g9
g10<-t(ga%*%t(sa*g9))
a10<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g10
g11<-t(ga%*%t(sa*g10))
a11<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g11
g12<-t(ga%*%t(sa*g11))
a12<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g12
g13<-t(ga%*%t(sa*g12))
a13<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g13
g14<-t(ga%*%t(sa*g13))
a14<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g14
g15<-t(ga%*%t(sa*g14))
a15<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g15
g16<-t(ga%*%t(sa*g15))
a16<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g16
g17<-t(ga%*%t(sa*g16))
a17<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g17
g18<-t(ga%*%t(sa*g17))
a18<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g18
g19<-t(ga%*%t(sa*g18))
a19<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g19
g20<-t(ga%*%t(sa*g19))
a20<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g20
g21<-t(ga%*%t(sa*g20))
a21<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g21
g22<-t(ga%*%t(sa*g21))
a22<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g22
g23<-t(ga%*%t(sa*g22))
a23<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g23
g24<-t(ga%*%t(sa*g23))
a24<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g24
g25<-t(ga%*%t(sa*g24))
a25<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g25
g26<-t(ga%*%t(sa*g25))
a26<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g26
g27<-t(ga%*%t(sa*g26))
a27<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g27
g28<-t(ga%*%t(sa*g27))
a28<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g28
g29<-t(ga%*%t(sa*g28))
a29<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g29
g30<-t(ga%*%t(sa*g29))
a30<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g30
g31<-t(ga%*%t(sa*g30))
a31<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g31
g32<-t(ga%*%t(sa*g31))
a32<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g32
g33<-t(ga%*%t(sa*g32))
a33<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g33
g34<-t(ga%*%t(sa*g33))
a34<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g34
g35<-t(ga%*%t(sa*g34))
a35<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g35
g36<-t(ga%*%t(sa*g35))
a36<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g36
g37<-t(ga%*%t(sa*g36))
a37<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g37
g38<-t(ga%*%t(sa*g37))
a38<-(1-exp(p*(x-2.45)^2+q*(x-2.45)+r)/(1+exp(p*(x-2.45)^2+q*(x-2.45)+r)))*g38

## This script looks at predicted population sizes through time
x<-seq(-0.85,5.2,by=0.05)
sa<-1-exp(0.629*(x-2.47)^2-1.7*(x-2.47)+3.6)/((1/0.021)*(1/(1-0.148))+exp(0.629*(x-2.47)^2-1.7*(x-2.47)+3.6))
sm<-1-exp(2.25*(x-2.28)^2-2.2*(x-2.28)+1.79)/((1/0.11)*(1/(1-0.257))+exp(2.25*(x-2.28)^2-2.2*(x-2.28)+1.79))
gm<-array(0,dim=c(length(x),length(x)))
ga<-array(0,dim=c(length(x),length(x)))
v<-c(0.001,0.001,0.001,0.002,0.004,0.005,0.008,0.011,0.015,0.02,0.026,0.033,0.04,0.048,0.054,0.06,0.065,0.068,0.069,0.068,0.065,0.06,0.054,0.048,0.04,0.033,0.026,0.02,0.015,0.011,0.008,0.005,0.004,0.002,0.001,0.001,0.001)
z<-c(0.001,0.001,0.002,0.003,0.005,0.008,0.011,0.016,0.023,0.03,0.039,0.048,0.058,0.066,0.073,0.077,0.078,0.077,0.073,0.066,0.058,0.048,0.039,0.03,0.023,0.016,0.011,0.008,0.005,0.003,0.002,0.001,0.001)
for(j in 32:122)
{for(w in 1:37)
{if(i<-(round((-0.423*(j/20-3.29)^2+0.564*(j/20-3.29)+2.22)*20))+w-1){gm[i,j]<-v[w]}}}
for(j in 1:122)
{for(y in 1:33)
{if(i<-(round((3.09+0.588*(j/20-3.85))*20))+y+1){ga[i,j]<-z[y]}}}
f<-(1-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)))*(1100*(x-3.48)+814)*0.5
# Specify functions for the fraction of mature (m) and not mature (nm) individuals
m<-(1-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)))
nm<-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51))
for(w in 1:122)
{if(f[w]<0){f[w]<-0}}
record<-numeric(100)
total<-numeric(10000)
# Read in a data sheet with the 96-year climate record
data<-read.table("C:\\\\Users\\\\CAS383\\\\Desktop\\\\Stochastic Climate Pool.txt",header=T)
# Create vector to represent the 96 years of the climate record
u<-seq(1,96,by=1)
# Set the number of years of simulation plus one (T)
T<-97
# Create matrices to hold the simulated size distributions of metamorphs (M), adults (Ad), and the total population (N)
N<-array(NA,dim=c(T,length(x)))
Ad<-array(NA,dim=c(T,length(x)))
M<-array(NA,dim=c(T,length(x)))
# Simulation begins with juveniles/adults at carrying capacity. Vector A[100, ] comes from the Density-dependent demography script.
Ad[1, ]<-A[100, ]
M[1, ]<-0
N[1, ]<-0
for(t in 1:(T-1))
  # If Oct-Jun rainfall is less than 426 mm, then the number of metamorphs recruited to the population is zero
{if(data[u[t],3]<426)
{Ad[t+1, ]<-t(gm%*%(1*sm*M[t, ]))+t(ga%*%(((sum(1*m*Ad[t, ])+sum(1*nm*Ad[t, ]))/sum(Ad[t, ]))*sa*Ad[t, ]))
Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(f*Ad[t+1, ])/101000)-1.23), sd = 0.227)
M[t+1, ]<-0*(Phi/sum(Phi))*sum(f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(f*Ad[t+1, ])/101000)-3.38)),0.0917)
N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}
  # If Oct-Jun rainfall is between 426-538 mm, then the number of metamorphs recruited to the popualtion is a fraction ((data[u[t],3]-426)/112) of the individuals that survived to late-stage larvae
  if(538>data[u[t],3] && data[u[t],3]>426)
    # Vector fa gives the fraction of breeding females based on Dec-Jan rainfall
  {fa<-min(0.295*(exp(0.00432*(data[u[t],2]-234))),1)
  Ad[t+1, ]<-t(gm%*%(1*sm*M[t, ]))+t(ga%*%(((sum(1*m*Ad[t, ])+sum(1*nm*Ad[t, ]))/sum(Ad[t, ]))*sa*Ad[t, ]))
  # Both metamorph body size distribution and embryonic/larval survival are now influenced by the fraction of breeding females
  Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(fa*f*Ad[t+1, ])/101000)-1.23), sd = 0.227)
  M[t+1, ]<-((data[u[t],3]-426)/112)*(Phi/sum(Phi))*sum(fa*f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(fa*f*Ad[t+1, ])/101000)-3.38)),0.0917)
  N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}
  # If Oct-Jun rainfall is above 538 mm, then the number of metamorphs recruited to the population is equal to the number of individuals that survived to late-stage larvae
  if(data[u[t],3]>538)
  {fa<-min(0.295*(exp(0.00432*(data[u[t],2]-234))),1)
  Ad[t+1, ]<-t(gm%*%(1*sm*M[t, ]))+t(ga%*%(((sum(1*m*Ad[t, ])+sum(1*nm*Ad[t, ]))/sum(Ad[t, ]))*sa*Ad[t, ]))
  Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(fa*f*Ad[t+1, ])/101000)-1.23), sd = 0.227)
  M[t+1, ]<-(Phi/sum(Phi))*sum(fa*f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(fa*f*Ad[t+1, ])/101000)-3.38)),0.0917)
  N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}}
# Create vector to store the simulation output
output<-numeric(96)
# Fill this vector with the population size in each year of the simulation
for(w in 1:96)
{output[w]<-sum(N[w+1, ])}
# Create vector with the number of adults in each year of the simulation
adults<-numeric(96)
for(w in 1:96)
{adults[w]<-sum(Ad[w+1, ])}
# Create vector with the number of metamorphs in each year of the simulation
meta<-numeric(96)
for(w in 1:96)
{meta[w]<-sum(M[w+1, ])}
# Create vector with the expected body size distribution of metamorphs over the course of the mark-recapture study
meta2<-numeric(122)
for(w in 1:122)
{meta2[w]<-sum(M[89:97,w])}
# Create vector (first) with the expected body size distribution of juveniles over the course of the mark-recapture study
juv2<-numeric(122)
for(w in 1:122)
{juv2[w]<-sum(Ad[90:97,w])}
first<-nm*juv2
# Create vector (second) with the expected body size distribution of adults over the course of the mark-recapture study
ad2<-numeric(122)
for(w in 1:122)
{ad2[w]<-sum(Ad[90:97,w])}
second<-m*ad2
# Create histograms of the expected body size distributions
hist(meta2)
hist(first)
hist(second)

## Ths script is used to track population dynamics driven by climatic fluctuations
x<-seq(-0.85,5.2,by=0.05)
sa<-1-exp(0.629*(x-2.47)^2-1.7*(x-2.47)+3.6)/((1/0.021)*(1/(1-0.148))+exp(0.629*(x-2.47)^2-1.7*(x-2.47)+3.6))
sm<-1-exp(2.25*(x-2.28)^2-2.2*(x-2.28)+1.79)/((1/0.11)*(1/(1-0.257))+exp(2.25*(x-2.28)^2-2.2*(x-2.28)+1.79))
gm<-array(0,dim=c(length(x),length(x)))
ga<-array(0,dim=c(length(x),length(x)))
v<-c(0.001,0.001,0.001,0.002,0.004,0.005,0.008,0.011,0.015,0.02,0.026,0.033,0.04,0.048,0.054,0.06,0.065,0.068,0.069,0.068,0.065,0.06,0.054,0.048,0.04,0.033,0.026,0.02,0.015,0.011,0.008,0.005,0.004,0.002,0.001,0.001,0.001)
z<-c(0.001,0.001,0.002,0.003,0.005,0.008,0.011,0.016,0.023,0.03,0.039,0.048,0.058,0.066,0.073,0.077,0.078,0.077,0.073,0.066,0.058,0.048,0.039,0.03,0.023,0.016,0.011,0.008,0.005,0.003,0.002,0.001,0.001)
for(j in 32:122)
{for(w in 1:37)
{if(i<-(round((-0.423*(j/20-3.29)^2+0.564*(j/20-3.29)+2.22)*20))+w-1){gm[i,j]<-v[w]}}}
for(j in 1:122)
{for(y in 1:33)
{if(i<-(round((3.09+0.588*(j/20-3.85))*20))+y+1){ga[i,j]<-z[y]}}}
f<-(1-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)))*(1100*(x-3.48)+814)*0.5
# Specify functions for the fraction of mature (m) and not mature (nm) individuals
m<-(1-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)))
nm<-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51))
for(w in 1:122)
{if(f[w]<0){f[w]<-0}}
# Create vector to hold simulated climate record
u<-numeric(100)
# Create vector to store the minimum population size during the simulation
record<-numeric(100)
# Create vector to store all 100 100-yr population trajectories
total<-numeric(10000)
# Read in a file with the 96-year climate record
data<-read.table("C:\\\\Users\\\\CAS383\\\\Desktop\\\\Stochastic Climate Pool Rev.txt",header=T)
for(a in 1:100)
  # Each entry in vector u is a random draw from the 96-year climate record
{for(w in 1:100)
{u[w]<-sample(1:96,1)}
  # Set number of years of simulation plus one (T)
  T<-101
  # Create matrices to hold the size distribution of metamorphs (M), juveniles/adults (Ad), and the total population (N)
  N<-array(NA,dim=c(T,length(x)))
  Ad<-array(NA,dim=c(T,length(x)))
  M<-array(NA,dim=c(T,length(x)))
  # Simulation begins with juveniles/adults at 63.7% of carrying capacity (the historic average population size). Vector A[100, ] comes from Density-Dependent demography file.
  Ad[1, ]<-0.637*A[100, ]
  M[1, ]<-0
  N[1, ]<-0
  # Enter the percentage of each age class for adults (ap), juveniles (jp), and metamorphs (mp) that are within the preserved area.
  ap<-1
  jp<-1
  mp<-1
  # Enter the volume of the breeding pond in cubic meters (cv)
  cv<-101000
  for(t in 1:(T-1))
  {if(data[u[t],3]<426)
  {Ad[t+1, ]<-t(gm%*%(mp*sm*M[t, ]))+t(ga%*%(jp*nm*sa*Ad[t, ]))+t(ga%*%(ap*m*sa*Ad[t, ]))
  Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(f*Ad[t+1, ])/cv)-1.23), sd = 0.227)
  M[t+1, ]<-0*(Phi/sum(Phi))*sum(f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(f*Ad[t+1, ])/cv)-3.38)),0.0917)
  N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}
    if(538>data[u[t],3] && data[u[t],3]>426)
    {fa<-min(0.295*(exp(0.00432*(data[u[t],2]-234))),1)
    Ad[t+1, ]<-t(gm%*%(mp*sm*M[t, ]))+t(ga%*%(jp*nm*sa*Ad[t, ]))+t(ga%*%(ap*m*sa*Ad[t, ]))
    Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(fa*f*Ad[t+1, ])/cv)-1.23), sd = 0.227)
    M[t+1, ]<-((data[u[t],3]-426)/112)*(Phi/sum(Phi))*sum(fa*f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(fa*f*Ad[t+1, ])/cv)-3.38)),0.0917)
    N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}
    if(data[u[t],3]>538)
    {fa<-min(0.295*(exp(0.00432*(data[u[t],2]-234))),1)
    Ad[t+1, ]<-t(gm%*%(mp*sm*M[t, ]))+t(ga%*%(jp*nm*sa*Ad[t, ]))+t(ga%*%(ap*m*sa*Ad[t, ]))
    Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(fa*f*Ad[t+1, ])/cv)-1.23), sd = 0.227)
    M[t+1, ]<-(Phi/sum(Phi))*sum(fa*f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(fa*f*Ad[t+1, ])/cv)-3.38)),0.0917)
    N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}}
  # Create vector output to store the population size over the 100 years of the simulation
  output<-numeric(100)
  for(w in 1:100)
  {output[w]<-sum(N[w+1, ])}
  # Fill vector that stores minimum population sizes
  record[a]<-min(output)
  # Fill vector that stores all population sizes across all 100 simulations
  total[(100*a-99):(100*a)]<-output}
record
# Create vector to calculate the number of times that simulations dropped below the quasi-extinction threshold of 3
v<-numeric(100)
for(w in 1:100)
{if(record[w]<3|is.nan(record[w])==TRUE)
{v[w]<-1}}
sum(v)

## This script is used to simulate the population dynamics under future climate conditions
x<-seq(-0.85,5.2,by=0.05)
sa<-1-exp(0.629*(x-2.47)^2-1.7*(x-2.47)+3.6)/((1/0.021)*(1/(1-0.148))+exp(0.629*(x-2.47)^2-1.7*(x-2.47)+3.6))
sm<-1-exp(2.25*(x-2.28)^2-2.2*(x-2.28)+1.79)/((1/0.11)*(1/(1-0.257))+exp(2.25*(x-2.28)^2-2.2*(x-2.28)+1.79))
gm<-array(0,dim=c(length(x),length(x)))
ga<-array(0,dim=c(length(x),length(x)))
v<-c(0.001,0.001,0.001,0.002,0.004,0.005,0.008,0.011,0.015,0.02,0.026,0.033,0.04,0.048,0.054,0.06,0.065,0.068,0.069,0.068,0.065,0.06,0.054,0.048,0.04,0.033,0.026,0.02,0.015,0.011,0.008,0.005,0.004,0.002,0.001,0.001,0.001)
z<-c(0.001,0.001,0.002,0.003,0.005,0.008,0.011,0.016,0.023,0.03,0.039,0.048,0.058,0.066,0.073,0.077,0.078,0.077,0.073,0.066,0.058,0.048,0.039,0.03,0.023,0.016,0.011,0.008,0.005,0.003,0.002,0.001,0.001)
for(j in 32:122)
{for(w in 1:37)
{if(i<-(round((-0.423*(j/20-3.29)^2+0.564*(j/20-3.29)+2.22)*20))+w-1){gm[i,j]<-v[w]}}}
for(j in 1:122)
{for(y in 1:33)
{if(i<-(round((3.09+0.588*(j/20-3.85))*20))+y+1){ga[i,j]<-z[y]}}}
f<-(1-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)))*(1100*(x-3.48)+814)*0.5
m<-(1-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)))
nm<-exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51)/(1+exp(3.3*(x-2.45)^2-13.1*(x-2.45)+4.51))
for(w in 1:122)
{if(f[w]<0){f[w]<-0}}
record<-numeric(100)
total<-numeric(10000)
T<-101
# Create matrix to hold simulated future climate record
data<-array(NA,dim=c(100,2))
for(a in 1:100)
{for(w in 1:100)
  # Dec-Jan rainfall in 2100 is simulated as a random draw from a parametric distribution
{data[w,1]<-254+rnorm(1,0,288)}
  for(w in 1:100)
    # Corrects any rainfall totals below zero to zero
  {if(data[w,1]<0)
  {data[w,1]<-0}}
  # Oct-Jun rainfall in 2100 is simulated as a random draw from a parametric distribution around the regression relating Oct-Jun rainfall to Dec-Jan rainfall
  for(w in 1:100)
  {data[w,2]<-1.15*data[w,1]+rnorm(1,0,158)+310}
  # If Oct-Jun rainfall is ever less than Dec-Jan rainfall this corrects them to be equal
  for(w in 1:100)
  {if(data[w,2]<data[w,1])
  {data[w,2]<-data[w,1]}}
  N<-array(NA,dim=c(T,length(x)))
  Ad<-array(NA,dim=c(T,length(x)))
  M<-array(NA,dim=c(T,length(x)))
  # Simulation starts with juveniles/adults at 62.6% of carrying capacity, which is projected to be the future population average
  Ad[1, ]<-0.626*A[100, ]
  M[1, ]<-0
  N[1, ]<-0
  ap<-1
  jp<-1
  mp<-1
  cv<-101000
  for(t in 1:(T-1))
  {if(data[t,2]<426)
  {Ad[t+1, ]<-t(gm%*%(mp*sm*M[t, ]))+t(ga%*%(jp*nm*sa*Ad[t, ]))+t(ga%*%(ap*m*sa*Ad[t, ]))
  Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(f*Ad[t+1, ])/cv)-1.23), sd = 0.227)
  M[t+1, ]<-0*(Phi/sum(Phi))*sum(f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(f*Ad[t+1, ])/cv)-3.38)),0.0917)
  N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}
    if(538>data[t,2] && data[t,2]>426)
    {fa<-min(0.295*(exp(0.00432*(data[t,1]-234))),1)
    Ad[t+1, ]<-t(gm%*%(mp*sm*M[t, ]))+t(ga%*%(jp*nm*sa*Ad[t, ]))+t(ga%*%(ap*m*sa*Ad[t, ]))
    Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(fa*f*Ad[t+1, ])/cv)-1.23), sd = 0.227)
    M[t+1, ]<-((data[t,2]-426)/112)*(Phi/sum(Phi))*sum(fa*f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(fa*f*Ad[t+1, ])/cv)-3.38)),0.0917)
    N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}
    if(data[t,2]>538)
    {fa<-min(0.295*(exp(0.00432*(data[t,1]-234))),1)
    Ad[t+1, ]<-t(gm%*%(mp*sm*M[t, ]))+t(ga%*%(jp*nm*sa*Ad[t, ]))+t(ga%*%(ap*m*sa*Ad[t, ]))
    Phi<-dnorm(x,mean = 2.27-0.272*(log(sum(fa*f*Ad[t+1, ])/cv)-1.23), sd = 0.227)
    M[t+1, ]<-(Phi/sum(Phi))*sum(fa*f*Ad[t+1, ])*min(0.00966*exp(-0.64*(log(sum(fa*f*Ad[t+1, ])/cv)-3.38)),0.0917)
    N[t+1, ]<-Ad[t+1, ]+M[t+1, ]}}
  output<-numeric(100)
  for(w in 1:100)
  {output[w]<-sum(N[w+1, ])}
  record[a]<-min(output)
  total[(100*a-99):(100*a)]<-output}
record
v<-numeric(100)
for(w in 1:100)
{if(record[w]<3|is.nan(record[w])==TRUE)
{v[w]<-1}}
sum(v)